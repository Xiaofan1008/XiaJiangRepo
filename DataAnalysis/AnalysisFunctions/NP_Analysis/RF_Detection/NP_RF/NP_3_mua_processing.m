%% ====================================================================
%  s3_noisegrid_mua.m  —  CJ280 noisegrid, STEP 3 of 4
%
%  MUA threshold crossings for every channel of the AP band. This is the
%  only expensive step (one pass over ~40 GB), so its output is cached and
%  step 4 can then be re-run with different windows and bin sizes in
%  seconds.
%
%  Detection, matching your V1/V2 raster script:
%    - 3rd order Butterworth highpass at 300 Hz, zero-phase (filtfilt)
%    - threshold = -THR_MAD * median(|x|)/0.6745
%    - negative-going crossings, NO refractory period (as requested)
%
%  Two details that matter for correctness:
%
%  1. The file is read in time chunks, not channel by channel. Pulling one
%     channel at a time from a 384-channel interleaved file means a strided
%     read across the whole 40 GB per channel -- hundreds of passes. Time
%     chunks are sequential.
%
%  2. Each chunk is read with PAD_S of context on both sides and the pad is
%     discarded after filtering. filtfilt has edge transients, and without
%     the pad every chunk boundary would produce spurious or missed
%     crossings, roughly 190 artefacts per channel.
%
%  The MAD threshold is computed once from chunks spread across the whole
%  recording and then held fixed, rather than recomputed per chunk, so a
%  noisy period cannot locally raise the threshold and silence real spikes.
%
%  Writes noisegrid_mua.mat into the recording folder.
% ====================================================================

clear; clc; close all;

%% ==================== CONFIG =======================================
REC_ROOT = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/CJ283/CJ283_noisegrid_110713_2026-09-16_11-07-18';
OUT_NAME = 'noisegrid_mua.mat';

AP_DIR   = '';            % empty = auto-detect
CHANNELS = 1:384;         % subset for a quick test, e.g. 1:16

% detection
HP_HZ    = 300;
HP_ORDER = 3;
THR_MAD  = 4;             % threshold = -THR_MAD * MAD

% common median reference across channels before detection.
% Off by default: it suppresses shared noise but also shared signal, and
% on a probe spanning one area it can remove real correlated MUA.
DO_CMR   = false;

% chunking
CHUNK_S  = 10;            % core seconds per chunk
PAD_S    = 0.5;           % context discarded after filtering
CH_BLOCK = 64;            % channels filtered per filtfilt call

% threshold estimation
N_MAD_CHUNKS = 12;        % chunks sampled across the recording

TIME_RANGE = [];          % [] = whole file, or [t0 t1] in seconds

%% ==================== LOCATE THE AP STREAM =========================
ob = dir(fullfile(REC_ROOT,'**','structure.oebin'));
assert(~isempty(ob),'No structure.oebin under %s',REC_ROOT);
rec_dir = ob(1).folder;
info    = jsondecode(fileread(fullfile(rec_dir,'structure.oebin')));
cont    = info.continuous;
if isstruct(cont), cont = num2cell(cont(:)'); end
names   = cellfun(@(s) string(s.folder_name), cont);
if isempty(AP_DIR)
    k = find(contains(names,'ProbeA-AP','IgnoreCase',true),1);
    assert(~isempty(k),'No ProbeA-AP stream'); AP_DIR = char(names(k));
end
iAP  = find(names==string(AP_DIR),1);
FS   = double(cont{iAP}.sample_rate);
NCH  = double(cont{iAP}.num_channels);
BV   = cont{iAP}.channels(1).bit_volts;

apFile = fullfile(rec_dir,'continuous',regexprep(AP_DIR,'[\\/]+$',''),'continuous.dat');
f  = dir(apFile);
nS = f.bytes/(NCH*2);
fprintf('recording : %s\n', rec_dir);
fprintf('AP        : %s\n', AP_DIR);
fprintf('            %g Hz, %d ch, bit_volts %g\n', FS, NCH, BV);
fprintf('            %.2f GB, %d samples, %.1f s\n', f.bytes/2^30, nS, nS/FS);

CHANNELS = CHANNELS(:)';
assert(all(CHANNELS >= 1 & CHANNELS <= NCH), 'CHANNELS out of range 1..%d', NCH);
nCh = numel(CHANNELS);

if isempty(TIME_RANGE)
    sFirst = 1;  sLast = nS;
else
    sFirst = max(1, floor(TIME_RANGE(1)*FS)+1);
    sLast  = min(nS, ceil(TIME_RANGE(2)*FS));
end
fprintf('analysing : %d channels, samples %d..%d (%.1f s)\n\n', ...
        nCh, sFirst, sLast, (sLast-sFirst+1)/FS);

mm = memmapfile(apFile,'Format',{'int16',[NCH nS],'x'});
[bHP,aHP] = butter(HP_ORDER, HP_HZ/(FS/2), 'high');

core = round(CHUNK_S*FS);
pad  = round(PAD_S*FS);
starts = sFirst:core:sLast;
nChunk = numel(starts);
fprintf('chunks: %d x %g s (+%g s pad each side)\n', nChunk, CHUNK_S, PAD_S);

%% ==================== PASS 1: THRESHOLDS ===========================
% MAD from chunks spread across the recording, so the threshold reflects
% the whole session rather than whatever the first chunk happened to hold.
pick = unique(round(linspace(1, nChunk, min(N_MAD_CHUNKS,nChunk))));
fprintf('\nestimating MAD from %d chunks ...\n', numel(pick));

madAcc = zeros(nCh, numel(pick));
tic;
for q = 1:numel(pick)
    [XF, ~, ~] = readFilterChunk(mm, starts(pick(q)), core, pad, sFirst, sLast, ...
                                 CHANNELS, BV, bHP, aHP, CH_BLOCK, DO_CMR);
    madAcc(:,q) = median(abs(XF), 1)' / 0.6745;
    fprintf('  %d/%d (%.1f s)\n', q, numel(pick), toc);
end
madCh = median(madAcc, 2);
thr   = -THR_MAD * madCh;

fprintf('\nMAD (uV): median %.2f, range %.2f .. %.2f\n', ...
        median(madCh), min(madCh), max(madCh));
fprintf('threshold: median %.2f uV\n', median(thr));
madSpread = max(madAcc,[],2)./max(min(madAcc,[],2),eps);
fprintf('MAD stability across chunks: median ratio %.2f, worst %.2f\n', ...
        median(madSpread), max(madSpread));
if max(madSpread) > 3
    warning(['MAD varies more than 3x across the session on at least one ' ...
             'channel — check for a noisy epoch or a moving probe.']);
end

%% ==================== PASS 2: DETECT ===============================
fprintf('\ndetecting ...\n');
buf = cell(nCh, nChunk);
tic;
for k = 1:nChunk
    [XF, r0, coreIx] = readFilterChunk(mm, starts(k), core, pad, sFirst, sLast, ...
                                       CHANNELS, BV, bHP, aHP, CH_BLOCK, DO_CMR);

    a0 = max(coreIx(1), 2);          % need one preceding sample
    i1 = coreIx(2);
    if a0 > i1, continue; end

    S = XF(a0:i1,   :);
    P = XF(a0-1:i1-1, :);
    hit = S < thr' & P >= thr';      % negative-going crossing, no refractory
    [ri, ci] = find(hit);
    if ~isempty(ri)
        samp = r0 + (a0 - 1) + ri - 1;          % absolute sample, 1-based
        tsp  = (samp - 1)/FS;                   % probe clock, same as s1
        for cc = unique(ci)'
            buf{cc,k} = tsp(ci==cc);
        end
    end

    if mod(k,10)==0 || k==nChunk
        el = toc;
        fprintf('  chunk %3d/%d  %.0f s elapsed, ~%.0f s left\n', ...
                k, nChunk, el, el/k*(nChunk-k));
    end
end

spikes = cell(nCh,1);
for cc = 1:nCh
    spikes{cc} = sort(vertcat(buf{cc,:}));
end
clear buf
nSpk = cellfun(@numel, spikes);

fprintf('\nspikes: total %d, per channel median %d, range %d .. %d\n', ...
        sum(nSpk), median(nSpk), min(nSpk), max(nSpk));
rate = nSpk / ((sLast-sFirst+1)/FS);
fprintf('rate  : median %.2f sp/s, range %.2f .. %.2f\n', ...
        median(rate), min(rate), max(rate));
dead = find(nSpk < 10);
if ~isempty(dead)
    fprintf('%d channel(s) with fewer than 10 crossings: %s\n', ...
            numel(dead), mat2str(CHANNELS(dead)));
end

%% ==================== SAVE =========================================
mua = struct();
mua.spikes    = spikes;          % {nCh x 1} spike times, s, probe clock
mua.channels  = CHANNELS(:);
mua.n_spikes  = nSpk(:);
mua.rate      = rate(:);
mua.thr_uv    = thr(:);
mua.mad_uv    = madCh(:);
mua.fs        = FS;
mua.bitvolts  = BV;
mua.hp_hz     = HP_HZ;
mua.hp_order  = HP_ORDER;
mua.thr_mad   = THR_MAD;
mua.refractory_ms = 0;
mua.cmr       = DO_CMR;
mua.chunk_s   = CHUNK_S;
mua.pad_s     = PAD_S;
mua.time_range= [(sFirst-1)/FS (sLast-1)/FS];
mua.ap_dir    = AP_DIR;
mua.rec_dir   = rec_dir;
mua.created   = datetime('now');
mua.script    = mfilename('fullpath');

outFile = fullfile(rec_dir, OUT_NAME);
save(outFile, '-struct', 'mua', '-v7.3');
fprintf('\nsaved -> %s\n', outFile);

%% ==================== DIAGNOSTIC PLOTS =============================
figure('Color','w','Position',[60 60 1150 780],'Name','s3 MUA');

subplot(2,2,1);
plot(CHANNELS, madCh, '.-','MarkerSize',8); grid on;
xlabel('channel'); ylabel('MAD (\muV)');
title(sprintf('noise floor (threshold = -%g x MAD)', THR_MAD));

subplot(2,2,2);
semilogy(CHANNELS, max(nSpk,0.5), '.-','MarkerSize',8); grid on;
xlabel('channel'); ylabel('spikes (log)');
title(sprintf('%d crossings total', sum(nSpk)));

subplot(2,2,3);
edges = 0:1:((sLast-sFirst+1)/FS);
allR = zeros(nCh, numel(edges)-1);
for cc = 1:nCh
    if ~isempty(spikes{cc}), allR(cc,:) = histcounts(spikes{cc}, edges); end
end
imagesc(edges(1:end-1), 1:nCh, allR); colorbar;
xlabel('time (s)'); ylabel('channel index');
title('rate over the session (1 s bins) — look for dropouts or artefact bands');

subplot(2,2,4);
[~,best] = max(nSpk);
if nSpk(best) > 1
    isi = diff(spikes{best})*1000;
    histogram(isi(isi<20), 200);
    xlabel('ISI (ms)'); ylabel('count');
    title(sprintf('ch %d ISI, %d spikes (no refractory imposed)', ...
                  CHANNELS(best), nSpk(best)));
end

%% ==================== LOCAL FUNCTIONS ==============================
function [XF, r0, coreIx] = readFilterChunk(mm, s0, core, pad, sFirst, sLast, ...
                                            CHANNELS, BV, b, a, chBlock, doCMR)
% Read one padded chunk, highpass it, return [nSamples x nChannels] in uV
% together with the absolute index of its first sample (r0) and the range
% of rows that form the core (pad excluded).
    s1 = min(s0 + core - 1, sLast);
    r0 = max(sFirst, s0 - pad);
    r1 = min(sLast,  s1 + pad);

    raw = mm.Data.x(CHANNELS, r0:r1);
    n   = size(raw,2);
    XF  = zeros(n, numel(CHANNELS));

    if doCMR
        % median across channels, per sample, removed before filtering
        med = median(double(raw), 1);
    end

    for c0 = 1:chBlock:numel(CHANNELS)
        c1 = min(c0+chBlock-1, numel(CHANNELS));
        Xb = double(raw(c0:c1,:))';
        if doCMR, Xb = Xb - med'; end
        XF(:,c0:c1) = filtfilt(b, a, Xb * BV);
    end

    coreIx = [s0 - r0 + 1, s1 - r0 + 1];
end