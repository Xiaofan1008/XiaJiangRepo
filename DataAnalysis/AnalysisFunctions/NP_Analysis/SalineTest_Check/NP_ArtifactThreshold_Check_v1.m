%% NP CHECK 4: stimulation artifact in the Neuropixels AP signal (script)
% Read-only. Cuts a window around every trigger on all 384 AP channels (small
% windows only, the 15 GB file is never loaded as a whole), subtracts the
% pre-trigger baseline per channel, and compares trials grouped by amplitude
% (control = visual/blank, 50 uA, 75 uA ...).
%
% "raw"      = the AP stream as recorded (the probe already applies a ~300 Hz
%              hardware high-pass), baseline-subtracted.
% "filtered" = 300-6000 Hz Butterworth, causal (needs Signal Processing Toolbox).
%
% Prints: noise level, the channels with the largest artifact (peak in 0-60 ms),
% how long the artifact lasts after the last pulse (40 ms), and
% plots: Fig 1 mean artifact over channel x time, Fig 2 peak amplitude per channel,
% Fig 3 single trials on the largest / middle / smallest artifact channel.
% Channel numbers are the 1-based index within the 384 recorded channels.
% Nothing is saved.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
npSession = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/2026_10_29/NP/Train_Amp_test1_expTrainAmps_133554_2026-10-29_13-36-51';
condFile  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_Amp_Test1_261029_133651/Train_Amp_test1.expTrainAmps.133554.conds.mat';
niStream  = 'NI-DAQmx-104.PXI-6133';
apStream  = 'Neuropix-PXI-106.ProbeA-AP';
trigAI    = 3;
nCh       = 384;
preMs     = 50;  postMs = 150;      % window around the trigger
peakWinMs = [0 60];                 % where the peak artifact is measured
lastPulseMs = 40;                   % 5 pulses at 100 Hz start at 0,10,20,30,40 ms
tailSigma = 5;                      % artifact "over" when |filtered| < this x noise
bandHz    = [300 6000];
nEx       = 8;                      % single trials drawn in Fig 3
showMs    = [-10 100];              % time range of the plots
%% ================================

npRec = resolveRec_(npSession);
fprintf('Recording: %s\n', npRec);

%% ---- trial labels and amplitude groups ----
cond = trialConds_(npRec);  Q = load(condFile);
label = arrayfun(@(k) char(Q.cond(k).label), cond, 'UniformOutput', false);  label = label(:);
nExp = numel(cond);
amp = nan(nExp,1);
for k = 1:nExp
    if contains(label{k},'spot') || contains(label{k},'blank'), amp(k) = 0;
    else
        tk = regexp(label{k}, '@(\d+)', 'tokens', 'once');
        if isempty(tk), amp(k) = -1; else, amp(k) = str2double(tk{1}); end
    end
end
levels = unique(amp)';
names = arrayfun(@(a) sprintf('%g uA', a), levels, 'UniformOutput', false);
names(levels == 0)  = {'control (visual/blank)'};
names(levels == -1) = {'electrical, amplitude not in label'};
for g = 1:numel(levels), fprintf('  %-34s %3d trials\n', names{g}, sum(amp == levels(g))); end

%% ---- trigger times (NP clock) from the NI analog input ----
niFile = fullfile(npRec,'continuous',niStream,'continuous.dat');
dN = dir(niFile);  NN = dN.bytes/2/8;
mN = memmapfile(niFile,'Format',{'int16',[8 NN],'a'});
x = single(mN.Data.a(trigAI+1,:)) * 0.0003051851;
th = (min(x(1:50:end)) + max(x(1:50:end))) / 2;
on = x > th;
up = find(~on(1:end-1) & on(2:end)) + 1;
clear x on
tsNI = readNpy_(fullfile(npRec,'continuous',niStream,'timestamps.npy'));
tTrig = double(tsNI(up));  clear tsNI mN
nTr = min(numel(tTrig), nExp);
fprintf('Triggers %d | stim-file trials %d -> using %d\n', numel(tTrig), nExp, nTr);
tTrig = tTrig(1:nTr);  ampT = amp(1:nTr);

%% ---- AP stream ----
oe = jsondecode(fileread(fullfile(npRec,'structure.oebin')));
cs = oe.continuous;  if iscell(cs), cs = [cs{:}]; end
iAP = find(arrayfun(@(s) contains(s.folder_name, apStream), cs), 1);
assert(~isempty(iAP), 'AP stream %s not found in structure.oebin', apStream);
fsAP = cs(iAP).sample_rate;
chn = cs(iAP).channels;  if iscell(chn), chn = chn{1}; else, chn = chn(1); end
bV = chn.bit_volts;
fprintf('AP stream: %g Hz, bit_volts %.4g (uV per count; expected 0.195)\n', fsAP, bV);

fA = fullfile(npRec,'continuous',apStream,'continuous.dat');
dA = dir(fA);  N = dA.bytes/2/nCh;
m = memmapfile(fA,'Format',{'int16',[nCh N],'a'});
tsAP = readNpy_(fullfile(npRec,'continuous',apStream,'timestamps.npy'));
idx = round((tTrig - double(tsAP(1))) * fsAP) + 1;
chk = idx >= 1 & idx <= N;
err = nan(nTr,1);  err(chk) = (double(tsAP(idx(chk))) - tTrig(chk)) * 1e3;
fprintf('Trigger -> AP sample index: max |timestamp mismatch| %.3f ms (should be < 0.04)\n', max(abs(err)));
clear tsAP

pre = round(preMs/1e3*fsAP);  post = round(postMs/1e3*fsAP);  nS = pre+post+1;
t = (-pre:post)/fsAP*1e3;
ok = (tTrig >= 0) & (idx-pre >= 1) & (idx+post <= N);
if any(~ok), fprintf('Skipping %d trial(s) (unsynchronised timestamp or window outside the file)\n', sum(~ok)); end

assert(exist('butter','file') > 0, 'butter not found: needs the Signal Processing Toolbox');
[bF, aF] = butter(3, bandHz/(fsAP/2));
bIdx  = find(t >= -preMs+5 & t < -3);                  % baseline / noise window
pkIdx = find(t >= peakWinMs(1) & t <= peakWinMs(2));
tlIdx = find(t >= 0);

%% ---- pass 1: per-trial statistics and group means ----
nG = numel(levels);
sumRaw = zeros(nCh,nS,nG,'single');  sumFlt = zeros(nCh,nS,nG,'single');  nInG = zeros(nG,1);
pkRaw = nan(nCh,nTr);  pkFlt = nan(nCh,nTr);  sigFlt = nan(nCh,nTr);  tailMs = nan(nCh,nTr);
tt = t(tlIdx);
fprintf('\nReading %d trials ...\n', sum(ok));
for k = 1:nTr
    if ~ok(k), continue, end
    [R, F] = cutTrial_(m, idx(k), pre, post, bV, bF, aF, bIdx);
    g = find(levels == ampT(k));
    sumRaw(:,:,g) = sumRaw(:,:,g) + single(R);  sumFlt(:,:,g) = sumFlt(:,:,g) + single(F);  nInG(g) = nInG(g) + 1;
    pkRaw(:,k) = max(abs(R(:,pkIdx)), [], 2);
    pkFlt(:,k) = max(abs(F(:,pkIdx)), [], 2);
    sigFlt(:,k) = 1.4826 * median(abs(F(:,bIdx)), 2);
    above = abs(F(:,tlIdx)) > tailSigma * sigFlt(:,k);
    lastI = max(above .* (1:numel(tlIdx)), [], 2);
    tv = zeros(nCh,1);  nz = lastI > 0;  tv(nz) = tt(lastI(nz));
    tailMs(:,k) = tv;
    if mod(k,50) == 0, fprintf('  %d / %d\n', k, nTr); end
end
meanRaw = sumRaw ./ reshape(max(nInG,1), 1, 1, nG);
meanFlt = sumFlt ./ reshape(max(nInG,1), 1, 1, nG);
clear sumRaw sumFlt

%% ---- per-group medians and printout ----
sigF = median(sigFlt(:,ok), 2, 'omitnan');
mpkR = nan(nCh,nG);  mpkF = nan(nCh,nG);  mtail = nan(nCh,nG);
for g = 1:nG
    s = find(ok & (ampT == levels(g)));
    mpkR(:,g) = median(pkRaw(:,s), 2, 'omitnan');  mpkF(:,g) = median(pkFlt(:,s), 2, 'omitnan');
    mtail(:,g) = median(tailMs(:,s), 2, 'omitnan');
end
gC = find(levels == 0);
fprintf('\nNoise before the trigger (filtered, robust sd): median over channels %.1f uV (min %.1f, max %.1f)\n', median(sigF), min(sigF), max(sigF));
for g = find(levels > 0)
    [~, o] = sort(mpkF(:,g), 'descend');
    fprintf('\n=== %s: %d trials ===\n', names{g}, nInG(g));
    fprintf('Channels with median peak (0-%d ms) above %d x noise: filtered %d of %d\n', peakWinMs(2), tailSigma, sum(mpkF(:,g) > tailSigma*sigF), nCh);
    fprintf('Largest artifact channels:\n  ch | peak raw uV | peak filt uV | control filt peak uV | artifact lasts after last pulse (ms)\n');
    for i = 1:10
        c = o(i);
        if isempty(gC), cp = nan; else, cp = mpkF(c,gC); end
        fprintf('  %3d | %11.0f | %12.0f | %20.0f | %8.1f\n', c, mpkR(c,g), mpkF(c,g), cp, max(0, mtail(c,g) - lastPulseMs));
    end
    c = o(round(nCh/2));  fprintf('  median-rank channel %d: raw %.0f, filt %.0f uV\n', c, mpkR(c,g), mpkF(c,g));
    c = o(end);           fprintf('  smallest channel %d:    raw %.0f, filt %.0f uV\n', c, mpkR(c,g), mpkF(c,g));
    fprintf('  tail over channels above 5 x noise: median %.1f ms, max %.1f ms (after the last pulse)\n', ...
        median(max(0, mtail(mpkF(:,g) > tailSigma*sigF, g) - lastPulseMs)), max(max(0, mtail(:,g) - lastPulseMs)));
end

%% ---- Fig 1: mean artifact over channel x time ----
sh = find(t >= showMs(1) & t <= showMs(2));
figure('Color','w','Position',[40 40 1500 800]);
cm = bwr_(256);
for r = 1:2
    if r == 1, M = meanRaw; rn = 'raw'; else, M = meanFlt; rn = 'filtered 300-6000 Hz'; end
    cl = 0;
    for g = 1:nG
        v = sort(abs(reshape(M(:,sh,g), [], 1)));  cl = max(cl, v(max(1, round(0.995*numel(v)))));
    end
    for g = 1:nG
        subplot(2, nG, (r-1)*nG + g);
        imagesc(t(sh), 1:nCh, M(:,sh,g));  axis xy;  colormap(gca, cm);  caxis([-cl cl] + eps);
        hold on;  plot([0 0], [1 nCh], 'k:');  colorbar;
        title(sprintf('%s | %s | mean, n=%d', names{g}, rn, nInG(g)), 'FontSize', 9);
        xlabel('ms from trigger');  ylabel('channel');
    end
end

%% ---- Fig 2: peak artifact per channel ----
figure('Color','w','Position',[60 60 1300 500]);
cmap = lines(nG);
for p = 1:2
    subplot(1,2,p);  hold on;
    if p == 1, Y = mpkR; yn = 'raw peak |uV|'; else, Y = mpkF; yn = 'filtered peak |uV|'; end
    for g = 1:nG, semilogy(1:nCh, Y(:,g), '.-', 'Color', cmap(g,:), 'DisplayName', names{g}); end
    if p == 2, semilogy(1:nCh, tailSigma*sigF, 'k--', 'DisplayName', sprintf('%d x noise', tailSigma)); end
    set(gca, 'YScale', 'log');  xlabel('channel');  ylabel(yn);  grid on;  legend('Location', 'best');
    title(sprintf('Median over trials, peak in 0-%d ms after the trigger', peakWinMs(2)));
end

%% ---- Fig 3: single trials on the largest / middle / smallest artifact channel ----
gTop = find(levels == max(levels));
[~, o] = sort(mpkF(:,gTop), 'descend');
chSel = [o(1), o(round(nCh/2)), o(end)];  chTxt = {'largest', 'median', 'smallest'};
selE = find(ok & ampT == levels(gTop), nEx);
selC = [];  if ~isempty(gC), selC = find(ok & ampT == 0, nEx); end
ER = nan(3,nS,numel(selE));  EF = ER;  CR = nan(3,nS,numel(selC));  CF = CR;
for i = 1:numel(selE), [R, F] = cutTrial_(m, idx(selE(i)), pre, post, bV, bF, aF, bIdx);  ER(:,:,i) = R(chSel,:);  EF(:,:,i) = F(chSel,:); end
for i = 1:numel(selC), [R, F] = cutTrial_(m, idx(selC(i)), pre, post, bV, bF, aF, bIdx);  CR(:,:,i) = R(chSel,:);  CF(:,:,i) = F(chSel,:); end
figure('Color','w','Position',[80 80 1500 800]);
for r = 1:2
    for c = 1:3
        subplot(2,3,(r-1)*3+c);  hold on;
        for i = 1:numel(selC)
            if r == 1, y = squeeze(CR(c,:,i)); else, y = squeeze(CF(c,:,i)); end
            plot(t, y, 'Color', [0.7 0.7 0.7]);
        end
        for i = 1:numel(selE)
            if r == 1, y = squeeze(ER(c,:,i)); else, y = squeeze(EF(c,:,i)); end
            plot(t, y, 'Color', [0.1 0.4 0.9 0.5]);
        end
        yl = ylim;  plot(0:10:lastPulseMs, yl(2)*ones(1,numel(0:10:lastPulseMs)), 'rv', 'MarkerFaceColor', 'r');
        xlim(showMs);  xlabel('ms from trigger');  ylabel('uV');
        if r == 1, rn = 'raw'; else, rn = 'filtered'; end
        title(sprintf('ch %d (%s artifact) | %s | blue %s, grey control', chSel(c), chTxt{c}, rn, names{gTop}), 'FontSize', 9);
    end
end
sgtitle('Single trials (red triangles = expected pulse starts)');

%% ======================== local functions ==============================
function [R, F] = cutTrial_(m, i0, pre, post, bV, b, a, bIdx)
X = double(m.Data.a(:, i0+(-pre:post))) * bV;      % uV, channels x samples
R = X - median(X(:,bIdx), 2);                      % baseline-subtracted raw
F = filter(b, a, R, [], 2);                        % causal band-pass
end

function cm = bwr_(n)
h = floor(n/2);
cm = [[linspace(0,1,h)', linspace(0,1,h)', ones(h,1)]; [ones(n-h,1), linspace(1,0,n-h)', linspace(1,0,n-h)']];
end

function rec = resolveRec_(p)
if exist(fullfile(p,'structure.oebin'),'file'), rec = p; return, end
c = dir(fullfile(p,'**','structure.oebin'));
assert(~isempty(c), 'No structure.oebin below %s', p);
assert(numel(c) == 1, 'Several recordings below %s:\n%s', p, strjoin(fullfile({c.folder},{c.name}), newline));
rec = c(1).folder;
end

function cond = trialConds_(npRec)
% condition number of every trial, from the 'TRIALID <cond>-<trial>' messages
txt = readNpy_(fullfile(npRec,'events','MessageCenter','text.npy'));
c = nan(0,1);  t = nan(0,1);
for i = 1:size(txt,1)
    tk = regexp(char(txt(i,:)), 'TRIALID (\d+)-(\d+)', 'tokens', 'once');
    if ~isempty(tk) && ~ismember(str2double(tk{2}), t), c(end+1,1) = str2double(tk{1}); t(end+1,1) = str2double(tk{2}); end %#ok<AGROW>
end
[~, o] = sort(t);  cond = c(o);
end

function A = readNpy_(fname)
fid = fopen(fname,'r');  assert(fid > 0, 'Cannot open %s', fname);
cl = onCleanup(@() fclose(fid));
magic = fread(fid,6,'uint8');  assert(isequal(magic', [147 78 85 77 80 89]), 'Not an npy file: %s', fname);
ver = fread(fid,2,'uint8');
if ver(1) == 1, hl = fread(fid,1,'uint16'); else, hl = fread(fid,1,'uint32'); end
hdr = char(fread(fid,hl,'uint8')');
tk = regexp(hdr,'''descr'':\s*''([^'']+)''','tokens','once');  descr = tk{1};
tk = regexp(hdr,'\(([^\)]*)\)','tokens','once');
dims = str2double(strsplit(strtrim(tk{1}),','));  dims = dims(~isnan(dims));
if isempty(dims), dims = 1; end
n = prod(dims);
if contains(descr,'S')                                  % fixed-width byte strings -> n x len uint8
    len = str2double(regexprep(descr,'[^0-9]',''));
    A = fread(fid,[len n],'uint8=>uint8')';  return
end
switch descr(2:end)
    case 'i2', t = 'int16';   case 'u2', t = 'uint16';
    case 'i4', t = 'int32';   case 'u4', t = 'uint32';
    case 'i8', t = 'int64';   case 'u8', t = 'uint64';
    case 'f4', t = 'single';  case 'f8', t = 'double';
    otherwise, error('Unsupported npy type %s', descr);
end
A = fread(fid, n, ['*' t]);
end