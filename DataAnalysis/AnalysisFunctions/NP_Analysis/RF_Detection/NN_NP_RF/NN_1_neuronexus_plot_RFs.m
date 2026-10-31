%% ====================================================================
%  n1_nnx_rf.m  —  RECEPTIVE FIELDS FROM THE NEURONEXUS (INTAN) PROBE
%
%  Standalone. Reads the Intan recording, recovers the frame times from
%  the photodiode digital line, reuses the noise sequence already
%  reconstructed by s2, and computes ON and OFF receptive fields.
%
%  Two points about the timing, both of which differ from the pipeline
%  used for the Neuropixels data.
%
%  1. TRIGGERS COME FROM THE ANALOG INPUT, NOT loadTrig. cleanTrig keeps
%     digital_in == 1, the stim line, so loadTrig returns one event per
%     trial. The noisegrid needs one per video frame -- ~216k of them at
%     120 Hz, with a noise update every 4th. The photodiode is on a board
%     ADC channel in analogin.dat; the channel carrying it is found by
%     scanning for the one that toggles at the frame rate, and detected
%     with a Schmitt trigger so that ringing on an edge cannot split one
%     flash into several frames.
%
%  2. NO CROSS-CLOCK SYNC IS REQUIRED. Both probes recorded the same
%     presentation through the same photodiode, so this script builds its
%     own frame times on the Intan clock and takes only the noise PATTERN
%     from noisegrid_stimulus.mat (act_idx / act_pol / updates_per_trial).
%     update_t in that file is on the Neuropixels clock and is ignored.
%
%  Spike detection matches s3 -- 300 Hz highpass, -THR_MAD*MAD threshold
%  crossings -- rather than the findpeaks/SD criterion in Spikes_Save.
%  MAD is not inflated by the spikes themselves the way SD is, and n2
%  compares these maps against the Neuropixels ones, so using the same
%  detector on both sides keeps that comparison clean.
%
%  Writes nnx_rf.mat into NNX_DIR.
% ====================================================================

clear; clc; close all;

%% ==================== CONFIG =======================================
NNX_DIR   = 'D:\Data\CJ283\noisegrid_260916_110717';
STIM_FILE = ['D:\Data\CJ283\CJ283_noisegrid_110713_2026-09-16_11-07-18\' ...
             'Record Node 101\experiment1\recording1\noisegrid_stimulus.mat'];

FS        = 30000;
N_CHAN    = 32;
GAIN      = 0.195;          % uV per bit
MAP_COL   = 5;              % ProbeMAP column, as in your scripts
SITE_UM   = 50;             % spacing between neighbouring sites
ANGLE_DEG = 90;
SITE1_IS_DEEPEST = true;    % electrode 1 = deepest. Flip if not.

% ---- photodiode (board ADC, analogin.dat) -----------------------------
ADC_FILE  = '';             % '' = look for analogin.dat, then
                            % board-ANALOG-IN-*.dat, in NNX_DIR
N_ADC     = [];             % [] = derive from the file size against the
                            % amplifier sample count. analogin.dat carries
                            % no channel count of its own.
DIODE_ADC = [];             % [] = auto-detect the channel toggling at the
                            % frame rate. Otherwise the 1-based index
                            % within analogin.dat (AI0 -> 1, AI1 -> 2).
FRAME_HZ_EXPECT = 120;
THR_HI_FRAC  = 0.70;        % Schmitt trigger, fractions of the p1..p99 range
THR_LO_FRAC  = 0.30;
MIN_STATE_US = 100;

% ---- what to plot -----------------------------------------------------
PLOT_CHANNELS = 2:2:32;        % electrode numbers for the detailed figures
PLOT_MOSAIC   = true;
PLOT_DEPTH    = true;

% ---- reverse correlation ---------------------------------------------
LAG_MS   = [-50 250];
BIN_MS   = 10;
PEAK_WIN = 3;
SMOOTH_SD = 1.0;
FLIP_Y   = true;

% ---- spike detection --------------------------------------------------
HP_HZ = 300;  THR_MAD = 4;  REFRAC_MS = 0;
CHUNK_S = 20;  PAD_S = 0.5;

SAVE_OUT = true;

%% ==================== PROBE MAP ====================================
old = pwd;  cleanupDir = onCleanup(@() cd(old));
cd(NNX_DIR);

E_MAP = ProbeMAP;
map   = E_MAP(2:1+N_CHAN, MAP_COL);
map_nums = cellfun(@(l) str2double(l(3:end)), map);
map_nums_plus = map_nums + 1;
if numel(map_nums_plus) > 32, map_nums_plus(33:end) = map_nums_plus(33:end) + 32; end
if numel(map_nums_plus) > 64, map_nums_plus(65:end) = map_nums_plus(65:end) + 32; end
assert(all(map_nums_plus >= 1 & map_nums_plus <= N_CHAN), ...
    'ProbeMAP gave data rows outside 1..%d.', N_CHAN);

depth_um = (0:N_CHAN-1) * SITE_UM * sind(ANGLE_DEG);
if SITE1_IS_DEEPEST, depth_um = flip(depth_um); end   % E1 -> 1550, E32 -> 0
fprintf('probe: %d sites, %g um apart, depths %g to %g um\n', ...
        N_CHAN, SITE_UM*sind(ANGLE_DEG), min(depth_um), max(depth_um));

%% ==================== PHOTODIODE FROM THE ANALOG INPUT =============
% The amplifier file fixes the sample count, and analogin.dat holds the same
% number of samples per channel, so the ADC channel count follows from the
% two file sizes. Nothing is assumed about how many inputs were saved.
afile = fullfile(NNX_DIR,'amplifier.dat');
assert(isfile(afile),'No amplifier.dat in %s', NNX_DIR);
fa = dir(afile);  nS = fa.bytes/(N_CHAN*2);
fprintf('amplifier: %d channels, %d samples (%.1f s)\n', N_CHAN, nS, nS/FS);

perChannelFiles = {};
if isempty(ADC_FILE)
    if isfile(fullfile(NNX_DIR,'analogin.dat'))
        ADC_FILE = fullfile(NNX_DIR,'analogin.dat');
    else
        d = dir(fullfile(NNX_DIR,'board-ANALOG-IN-*.dat'));
        assert(~isempty(d), ...
            ['No analogin.dat or board-ANALOG-IN-*.dat in %s. Set ADC_FILE.'], NNX_DIR);
        [~,o] = sort({d.name});  d = d(o);
        perChannelFiles = arrayfun(@(x) fullfile(x.folder,x.name), d, 'uni', 0);
        fprintf('found %d per-channel ADC files: %s ...\n', numel(d), d(1).name);
    end
end

if isempty(perChannelFiles)
    fi = dir(ADC_FILE);
    if isempty(N_ADC)
        N_ADC = (fi.bytes/2)/nS;
        assert(abs(N_ADC-round(N_ADC)) < 1e-6, ...
            ['analogin.dat holds %.3f channels against the amplifier sample ' ...
             'count, which is not a whole number. Set N_ADC by hand.'], N_ADC);
        N_ADC = round(N_ADC);
    end
    fprintf('analogin.dat: %d channel(s), %d samples\n', N_ADC, nS);
    mmA = memmapfile(ADC_FILE,'Format',{'uint16',[N_ADC nS],'x'});
    getADC = @(c) double(mmA.Data.x(c,:))';
else
    N_ADC = numel(perChannelFiles);
    getADC = @(c) readOneADC(perChannelFiles{c}, nS);
end

% scan each ADC channel for the one toggling at the frame rate
fprintf('\n%-6s %10s %10s %12s %12s\n','ch','p1','p99','pk-pk','period (ms)');
scanDec = 5;
cands = nan(N_ADC,3);
for c = 1:N_ADC
    x = getADC(c);
    xs = x(1:scanDec:end);
    q = prctile(xs,[1 99]);
    hi = q(1) + THR_HI_FRAC*(q(2)-q(1));
    up = find(xs(1:end-1)<=hi & xs(2:end)>hi);
    if numel(up) < 100
        fprintf('%-6d %10.0f %10.0f %12.0f %12s\n', c, q(1), q(2), q(2)-q(1), '-');
        continue;
    end
    per = median(diff(up))*scanDec/FS*1000;
    fprintf('%-6d %10.0f %10.0f %12.0f %12.3f\n', c, q(1), q(2), q(2)-q(1), per);
    cands(c,:) = [c, q(2)-q(1), per];
end
ok = ~isnan(cands(:,1));
assert(any(ok),'No ADC channel shows repeated transitions.');

if isempty(DIODE_ADC)
    sub = cands(ok,:);
    [~,k] = min(abs(sub(:,3) - 1000/FRAME_HZ_EXPECT));
    DIODE_ADC = sub(k,1);
    fprintf('auto-selected ADC channel %d (%.2f Hz, %.0f counts pk-pk)\n', ...
            DIODE_ADC, 1000/sub(k,3), sub(k,2));
else
    fprintf('using configured ADC channel %d\n', DIODE_ADC);
end
assert(DIODE_ADC>=1 && DIODE_ADC<=N_ADC,'DIODE_ADC out of range 1..%d',N_ADC);

d = getADC(DIODE_ADC);
lo = prctile(d,1);  hi = prctile(d,99);
thrHi = lo + THR_HI_FRAC*(hi-lo);
thrLo = lo + THR_LO_FRAC*(hi-lo);
fprintf('p1 %.0f | p99 %.0f | Schmitt hi %.0f / lo %.0f\n', lo, hi, thrHi, thrLo);

% Schmitt trigger with strict alternation, as in s1: a single threshold on
% an analog trace splits one flash into several on any ringing edge.
upC = find(d(1:end-1) <= thrHi & d(2:end) >  thrHi) + 1;
dnC = find(d(1:end-1) >= thrLo & d(2:end) <  thrLo) + 1;
ev  = sortrows([upC(:) ones(numel(upC),1); dnC(:) -ones(numel(dnC),1)],1);
upI = zeros(numel(upC),1);  dnI = zeros(numel(upC),1);
kk = 0;  armed = true;
for e = 1:size(ev,1)
    if armed && ev(e,2) > 0
        kk = kk+1;  upI(kk) = ev(e,1);  armed = false;
    elseif ~armed && ev(e,2) < 0
        dnI(kk) = ev(e,1);  armed = true;
    end
end
if kk > 0 && dnI(kk) == 0, kk = kk-1; end
upI = upI(1:kk);  dnI = dnI(1:kk);
minSamp = max(1, round(MIN_STATE_US/1e6*FS));
keepS = (dnI-upI) >= minSamp;
upI = upI(keepS);  dnI = dnI(keepS);
clear d

frame_t = (upI-1)/FS;
per_ms  = diff(frame_t)*1000;
fprintf('frames: %d, high %.3f ms median, period %.4f ms -> %.2f Hz\n', ...
        numel(frame_t), median((dnI-upI)/FS*1000), median(per_ms), 1000/median(per_ms));

gapThr = max(6*median(per_ms), 100);
brk = find(per_ms > gapThr);
trS = [1; brk+1];  trE = [brk; numel(frame_t)];
nFr = trE - trS + 1;
nTrialD = numel(trS);
inTrial = true(numel(per_ms),1);  inTrial(brk) = false;
frameMs = median(per_ms(inTrial));
nDrop = sum(round(per_ms(inTrial)/frameMs) - 1);
fprintf('trials: %d, %d..%d frames each, %d dropped frames\n', ...
        nTrialD, min(nFr), max(nFr), nDrop);

%% ==================== NOISE SEQUENCE FROM s2 =======================
assert(isfile(STIM_FILE),'Cannot find %s', STIM_FILE);
S = load(STIM_FILE,'act_idx','act_pol','updates_per_trial','frames_per_update', ...
                   'n_updates','n_trials','n_randels','size_h','size_w', ...
                   'randel_row','randel_col','n_active','geometry','validated');
fprintf('\nstimulus file: %d updates over %d trials, %dx%d randels, validated %d\n', ...
        S.n_updates, S.n_trials, S.size_h, S.size_w, S.validated);
assert(S.validated,'The stimulus reconstruction did not validate; do not use it.');
assert(nTrialD == S.n_trials, ...
    ['The photodiode found %d trials but the stimulus file has %d. These ' ...
     'must be the same presentation.'], nTrialD, S.n_trials);

F = S.frames_per_update;
cb = double(S.updates_per_trial(:));
need = (cb-1)*F + 1;
slack = nFr(:) - need;
fprintf('frames per update %d | trailing frames after the last update: min %d, max %d\n', ...
        F, min(slack), max(slack));
assert(all(slack >= 0), ...
    '%d trial(s) have fewer frames than the update count requires.', sum(slack<0));

nUpd = sum(cb);
update_t = zeros(nUpd,1);  q = 0;
for i = 1:nTrialD
    n = cb(i);
    update_t(q+(1:n)) = frame_t(trS(i) + (0:n-1)'*F);
    q = q + n;
end
assert(q == S.n_updates && issorted(update_t), 'Update construction failed.');
fprintf('update onsets: %d, %.4f .. %.4f s on the Intan clock\n', ...
        nUpd, update_t(1), update_t(end));

nPix = double(S.n_randels);  H = double(S.size_h);  W = double(S.size_w);

%% ==================== MUA ==========================================
mm = memmapfile(afile,'Format',{'int16',[N_CHAN nS],'x'});
[bHP,aHP] = butter(3, HP_HZ/(FS/2), 'high');
core = round(CHUNK_S*FS);  pad = round(PAD_S*FS);
starts = 1:core:nS;  nChunk = numel(starts);

pick = unique(round(linspace(1,nChunk,min(12,nChunk))));
madAcc = zeros(N_CHAN,numel(pick));
for q = 1:numel(pick)
    XF = readChunkN(mm,starts(pick(q)),core,pad,nS,N_CHAN,GAIN,bHP,aHP);
    madAcc(:,q) = median(abs(XF),1)'/0.6745;
end
madCh = median(madAcc,2);  thr = -THR_MAD*madCh;
fprintf('MAD %.2f uV median (%.2f .. %.2f), threshold %.2f uV\n', ...
        median(madCh), min(madCh), max(madCh), median(thr));

buf = cell(N_CHAN,nChunk);  refr = round(REFRAC_MS/1000*FS);  tic;
for k = 1:nChunk
    s0 = starts(k);  s1 = min(s0+core-1,nS);  r0 = max(1,s0-pad);
    XF = readChunkN(mm,s0,core,pad,nS,N_CHAN,GAIN,bHP,aHP);
    a0 = max(s0-r0+1,2);  i1 = s1-r0+1;
    if a0 > i1, continue; end
    hit = XF(a0:i1,:) < thr' & XF(a0-1:i1-1,:) >= thr';
    [ri,ci] = find(hit);
    if ~isempty(ri)
        samp = r0 + (a0-1) + ri - 1;
        for cc = unique(ci)'
            s = sort(samp(ci==cc));
            if refr > 0 && numel(s) > 1, s = s([true; diff(s)>refr]); end
            buf{cc,k} = (s-1)/FS;
        end
    end
    if mod(k,10)==0 || k==nChunk
        fprintf('  chunk %3d/%d (%.0f s)\n', k, nChunk, toc);
    end
end
spikesRow = cell(N_CHAN,1);
for cc = 1:N_CHAN, spikesRow{cc} = sort(vertcat(buf{cc,:})); end
clear buf

% reorder from data rows to electrode numbers
spikes = cell(N_CHAN,1);
for e = 1:N_CHAN, spikes{e} = spikesRow{map_nums_plus(e)}; end
nSpk = cellfun(@numel,spikes);
fprintf('spikes per electrode: median %.0f, range %d..%d\n', ...
        median(nSpk), min(nSpk), max(nSpk));

%% ==================== REVERSE CORRELATION ==========================
edges = LAG_MS(1):BIN_MS:LAG_MS(2);
ctrs  = edges(1:end-1)+BIN_MS/2;  nBin = numel(ctrs);

cols = repmat(uint32(1:nUpd), S.n_active, 1);
isOn = S.act_pol > 0;
Mon  = sparse(double(S.act_idx(isOn)),  double(cols(isOn)),  1, nPix, nUpd);
Moff = sparse(double(S.act_idx(~isOn)), double(cols(~isOn)), 1, nPix, nUpd);
clear cols isOn
nOn = full(sum(Mon,2));  nOff = full(sum(Moff,2));
fprintf('\nON %.0f, OFF %.0f presentations per randel\n', mean(nOn), mean(nOff));

zOn = zeros(nPix,nBin,N_CHAN,'single');  zOff = zOn;
psth = zeros(nBin,N_CHAN);
for e = 1:N_CHAN
    cnt = countPerUpdateN(spikes{e}, update_t, edges);
    mAll = mean(cnt,1);  sdAll = std(cnt,0,1);
    mOnP  = full(Mon *cnt)./nOn;
    mOffP = full(Moff*cnt)./nOff;
    zOn(:,:,e)  = single((mOnP  - mAll)./max(sdAll./sqrt(nOn), eps));
    zOff(:,:,e) = single((mOffP - mAll)./max(sdAll./sqrt(nOff),eps));
    psth(:,e)   = mAll/(BIN_MS/1000);
end

sprOn  = squeeze(std(double(zOn),0,1));   sprOff = squeeze(std(double(zOff),0,1));
if N_CHAN==1, sprOn=sprOn(:); sprOff=sprOff(:); end
spr = max(sprOn,sprOff);
posB = find(ctrs>0);  negB = find(ctrs<-BIN_MS);
[pkVal,pr] = max(spr(posB,:),[],1);  pkBin = posB(pr);  pkLag = ctrs(pkBin);
nullLvl = mean(spr(negB,:),1);
snr = pkVal./max(nullLvl,eps);

fprintf('\n%-5s %8s %8s %9s %8s %8s\n','elec','depth','spikes','peak ms','spread','ratio');
for e = 1:N_CHAN
    fprintf('%-5d %8.0f %8d %9.0f %8.2f %8.2f\n', ...
            e, depth_um(e), nSpk(e), pkLag(e), pkVal(e), snr(e));
end

%% ==================== AXES IN DEGREES ==============================
g = S.geometry;
if isstruct(g) && isfield(g,'width') && ~isempty(g.width)
    xC = g.X - g.width/2  + ((1:W)-0.5)*g.width/W;
    yC = g.Y + g.height/2 - ((1:H)-0.5)*g.height/H;
    unitStr = 'deg';
else
    xC = 1:W; yC = 1:H; unitStr = 'randels';
end
if FLIP_Y, yC = fliplr(yC); end
cmap = divMapN(256);

%% ==================== FIGURES ======================================
for e = PLOT_CHANNELS(:)'
    b = pkBin(e);
    wIx = max(1,b-floor(PEAK_WIN/2)):min(nBin,b+floor(PEAK_WIN/2));
    mOnI  = smoothN(mean(double(zOn(:,wIx,e)),2), H,W,SMOOTH_SD);
    mOffI = smoothN(mean(double(zOff(:,wIx,e)),2),H,W,SMOOTH_SD);
    cl = max(prctile(abs([mOnI(:);mOffI(:)]),99.5),1);

    figure('Color','w','Position',[60 60 1150 340], ...
           'Name',sprintf('E%d',e),'NumberTitle','off');
    subplot(1,4,1);
    imagesc(xC,yC,mOnI,[-cl cl]); axis image xy; colormap(gca,cmap); colorbar;
    title(sprintf('E%d ON, %g ms',e,pkLag(e)));
    xlabel(sprintf('azimuth (%s)',unitStr)); ylabel(sprintf('elevation (%s)',unitStr));
    subplot(1,4,2);
    imagesc(xC,yC,mOffI,[-cl cl]); axis image xy; colormap(gca,cmap); colorbar;
    title('OFF');
    subplot(1,4,3);
    imagesc(xC,yC,mOnI-mOffI,[-cl cl]); axis image xy; colormap(gca,cmap); colorbar;
    title('ON - OFF');
    subplot(1,4,4);
    plot(ctrs,sprOn(:,e),'-o','MarkerSize',3); hold on;
    plot(ctrs,sprOff(:,e),'-s','MarkerSize',3);
    yline(nullLvl(e),'k:'); xline(pkLag(e),'r--'); grid on;
    xlabel('lag (ms)'); ylabel('SD of z'); legend('ON','OFF','Location','best');
    title(sprintf('depth %g \\mum, ratio %.2f', depth_um(e), snr(e)));
end

if PLOT_MOSAIC
    for pass = 1:2
        nm = {'ON','OFF'};
        figure('Color','w','Position',[40 40 1250 780], ...
               'Name',sprintf('%s mosaic',nm{pass}),'NumberTitle','off');
        [~,ord] = sort(depth_um);
        for i = 1:N_CHAN
            e = ord(i);
            b = pkBin(e);
            wIx = max(1,b-floor(PEAK_WIN/2)):min(nBin,b+floor(PEAK_WIN/2));
            if pass==1, mp = mean(double(zOn(:,wIx,e)),2);
            else,       mp = mean(double(zOff(:,wIx,e)),2); end
            mp = smoothN(mp,H,W,SMOOTH_SD);
            cl = max(prctile(abs(mp(:)),99.5),1);
            subplot(4,8,i);
            imagesc(mp,[-cl cl]); axis image off xy; colormap(gca,cmap);
            title(sprintf('E%d | %g\\mum | %.1f', e, depth_um(e), snr(e)),'FontSize',7);
        end
        sgtitle(sprintf('%s receptive fields, shallow to deep', nm{pass}));
    end
end

%% ==================== DEPTH PROGRESSION ============================
cx = nan(N_CHAN,1); cy = nan(N_CHAN,1); okFit = false(N_CHAN,1);
for e = 1:N_CHAN
    b = pkBin(e);
    wIx = max(1,b-floor(PEAK_WIN/2)):min(nBin,b+floor(PEAK_WIN/2));
    m = smoothN(mean(double(zOn(:,wIx,e)),2)+mean(double(zOff(:,wIx,e)),2), H,W,SMOOTH_SD)/2;
    [pk,pkIx] = max(m(:));
    if pk < 2 || snr(e) < 1.05, continue; end
    thrm = max(2, 0.5*pk);
    mask = growN(m>=thrm, pkIx, H, W);
    if nnz(mask) < 3, continue; end
    [rr,cc] = find(mask);  w = m(mask)-thrm;
    if sum(w) <= 0, continue; end
    cx(e) = sum(w(:).*xC(cc)')/sum(w);
    cy(e) = sum(w(:).*yC(rr)')/sum(w);
    okFit(e) = true;
end
fprintf('\n%d of %d electrodes gave a fittable receptive field\n', sum(okFit), N_CHAN);

if PLOT_DEPTH && sum(okFit) >= 3
    figure('Color','w','Position',[70 70 1150 400], ...
           'Name','RF vs depth','NumberTitle','off');
    subplot(1,3,1); hold on;
    scatter(cx(okFit), cy(okFit), 40, depth_um(okFit), 'filled');
    [~,o] = sort(depth_um(okFit));  ix = find(okFit);  ix = ix(o);
    plot(cx(ix), cy(ix), 'k-');
    axis equal; grid on; box on; colormap(gca,parula); cb = colorbar;
    cb.Label.String = 'depth (\mum)';
    xlim([min(xC) max(xC)]); ylim([min(yC) max(yC)]);
    xlabel(sprintf('azimuth (%s)',unitStr)); ylabel(sprintf('elevation (%s)',unitStr));
    title('RF centres along the shank');
    subplot(1,3,2);
    plot(depth_um(okFit), cx(okFit), 'o-'); grid on;
    xlabel('depth (\mum)'); ylabel(sprintf('azimuth (%s)',unitStr)); title('azimuth');
    subplot(1,3,3);
    plot(depth_um(okFit), cy(okFit), 'o-'); grid on;
    xlabel('depth (\mum)'); ylabel(sprintf('elevation (%s)',unitStr)); title('elevation');
end

%% ==================== SAVE =========================================
if SAVE_OUT
    rf = struct();
    rf.probe = 'neuronexus';
    rf.channels = (1:N_CHAN)';
    rf.data_row = map_nums_plus(:);
    rf.depth_um = depth_um(:);
    rf.zOn = zOn;  rf.zOff = zOff;  rf.psth = psth;
    rf.lag_ms = ctrs(:);  rf.bin_ms = BIN_MS;
    rf.peak_bin = pkBin(:);  rf.peak_lag_ms = pkLag(:);  rf.snr = snr(:);
    rf.n_spikes = nSpk(:);  rf.mad_uv = madCh(:);
    rf.size_h = H;  rf.size_w = W;
    rf.randel_row = S.randel_row;  rf.randel_col = S.randel_col;
    rf.x_deg = xC;  rf.y_deg = yC;  rf.flip_y = FLIP_Y;
    rf.n_on = nOn;  rf.n_off = nOff;
    rf.cx = cx;  rf.cy = cy;  rf.fitted = okFit;
    rf.stim_file = STIM_FILE;  rf.nnx_dir = NNX_DIR;
    rf.created = datetime('now');
    save(fullfile(NNX_DIR,'nnx_rf.mat'),'-struct','rf','-v7.3');
    fprintf('\nsaved -> %s\n', fullfile(NNX_DIR,'nnx_rf.mat'));
end

%% ==================== LOCAL FUNCTIONS ==============================
function x = readOneADC(fpath, nS)
    fi = dir(fpath);
    n  = fi.bytes/2;
    if n ~= nS
        warning('%s has %d samples, amplifier has %d.', fpath, n, nS);
    end
    fid = fopen(fpath,'r');
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>
    x = double(fread(fid, n, 'uint16'));
end

function XF = readChunkN(mm,s0,core,pad,nS,NCH,GAIN,b,a)
    s1 = min(s0+core-1,nS);
    r0 = max(1,s0-pad);  r1 = min(nS,s1+pad);
    XF = filtfilt(b, a, double(mm.Data.x(:, r0:r1))'*GAIN);
end

function cnt = countPerUpdateN(st,u,edges)
    nU = numel(u); nE = numel(edges);
    if isempty(st), cnt = zeros(nU,nE-1); return; end
    st = st(:);  cum = zeros(nU,nE);
    for e = 1:nE
        q = u + edges(e)/1000;
        h = histcounts(st, [-Inf; q; Inf]);
        cum(:,e) = cumsum(h(1:end-1))';
    end
    cnt = diff(cum,1,2);
end

function m = smoothN(v,H,W,sd)
    m = reshape(v,H,W);
    if sd > 0
        r = ceil(3*sd);
        [gx,gy] = meshgrid(-r:r,-r:r);
        g = exp(-(gx.^2+gy.^2)/(2*sd^2)); g = g/sum(g(:));
        m = conv2(m,g,'same');
    end
end

function mask = growN(bin,seed,H,W)
    mask = false(H,W);
    if ~bin(seed), return; end
    mask(seed) = true;  stack = seed;
    while ~isempty(stack)
        p = stack(end); stack(end) = [];
        r = mod(p-1,H)+1; c = floor((p-1)/H)+1;
        for dr = -1:1
            for dc = -1:1
                rr = r+dr; cc = c+dc;
                if rr>=1 && rr<=H && cc>=1 && cc<=W && bin(rr,cc) && ~mask(rr,cc)
                    mask(rr,cc) = true;  stack(end+1) = (cc-1)*H+rr; %#ok<AGROW>
                end
            end
        end
    end
end

function cm = divMapN(n)
    m = floor(n/2);
    lo = [linspace(0.02,1,m)' linspace(0.19,1,m)' linspace(0.38,1,m)'];
    hi = [linspace(1,0.70,n-m)' linspace(1,0.09,n-m)' linspace(1,0.17,n-m)'];
    cm = [lo; hi];
end