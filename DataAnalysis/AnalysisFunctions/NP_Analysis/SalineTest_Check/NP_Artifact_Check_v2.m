%% NP CHECK 6: is every condition delivered as designed? (Neuropixels recording, one session) (script)
% Read-only. Replaces the threshold on the median envelope of np_check5 with:
%   1. the EXPECTED pulse times from the stim PC .conds.mat (amp / ptd / nsp / frequency),
%   2. detection on the channels that carry the artifact (chosen automatically),
%   3. a threshold PER CHANNEL, from that channel's own noise in THIS trial (adaptive),
%   4. a pulse is accepted only when many of those channels step at the same sample
%      (artifact = all channels together; noise spikes and slow drift are not).
% Detection uses the sample-to-sample step, so pulses 2 ms apart are separated and a slow
% baseline drift does not matter.
%
% Prints: run parameters, the channels used, the condition list, one line per condition
% (expected / detected pulses, trials fully correct, missing and extra pulses, worst timing
% error), the failing trials. Fig 1: detected (black) vs expected (red circles) per condition.
% Fig 2+: raw and filtered traces of chosen conditions on chosen channels (one figure per channel).
% Nothing is saved.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
npSession = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/2026_10_29/NP/Train_Amp_test2_expTrainAmps_135452_2026-10-29_13-55-01';
condFile  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_Amp_Test2_261029_135502/Train_Amp_test2.expTrainAmps.135452.conds.mat';
niStream  = 'NI-DAQmx-104.PXI-6133';
apStream  = 'Neuropix-PXI-106.ProbeA-AP';
trigAI    = 3;
nCh       = 384;
preMs     = 20;  postMs = 120;   % window around the trigger (must cover the whole pulse train)
kThr      = 3;      % a channel "steps" when |step| > kThr x its own noise (robust sd of the steps before the trigger)
fracMin   = 0.3;    % a pulse = at least this fraction of the used channels step at the same sample
mergeMs   = 0.7;    % steps closer than this belong to one pulse (the phases of a biphasic pulse)
tolMs     = 0.3;    % allowed timing error of a pulse (after removing the constant trigger offset)
matchMs   = 1.5;    % only for estimating the constant offset: nearest detected pulse within this
nEach     = 30;     % channels used per electrode (the nEach strongest artifact channels of electrode A and of B)
nSelTr    = 15;     % trials (largest amplitude) used to choose the channels
snrMin    = 8;      % a channel must have a pulse step at least this x its noise to be used
% Fig 2 (optional): conditions = numbers of the printed list; channels = AP channel numbers.
% Leave [] to choose automatically (the longest train of each type, the 3 strongest channels).
plotConds = [];
plotChan  = [];
showMs    = [];     % [] = from -5 ms to the end of the longest pulse train + 15 ms
bandHz    = [300 6000];
%% ================================

npRec = resolveRec_(npSession);
fprintf('Recording: %s\n', npRec);

%% ---- design: expected pulse times of every condition ----
Q = load(condFile);  a = Q.args;
fprintf('Run parameters: amp %s uA | freq %g Hz | nPulses %g | order %s | isi_ms %s | nReps %g | pulseTime %g | trialDur %g | visDur %g\n', ...
    mat2str(double(a.amp)), a.freq, a.nPulses, string(a.order), mat2str(double(Q.isi_ms)), a.nReps, a.pulseTime, a.trialDur, a.visDur);
periodMs = 1e3 / double(a.freq);
nCq = numel(Q.cond);
expT = cell(nCq,1);  ampC = zeros(nCq,2);  hasE = false(nCq,2);
for c = 1:nCq
    v = [];
    for e = 1:2
        ampC(c,e) = double(Q.cond(c).amp(e));
        hasE(c,e) = ampC(c,e) > 0 && double(Q.cond(c).nsp(e)) > 0;
        if hasE(c,e)
            v = [v, double(Q.cond(c).ptd(e))/1e3 + (0:double(Q.cond(c).nsp(e))-1)*periodMs]; %#ok<AGROW>
        end
    end
    v = sort(v);  if ~isempty(v), v = v([true, diff(v) > 0.05]); end      % simultaneous pulses = one event
    expT{c} = v;
end
nEcond = cellfun(@numel, expT);
lastExp = max([0; cellfun(@(v) max([0 v]), expT)]);
if lastExp + 10 > postMs
    fprintf('!! the longest pulse train ends at %.0f ms but postMs is %g: increase postMs\n', lastExp, postMs);
end
if isempty(showMs), showMs = [-5, min(postMs, ceil(lastExp) + 15)]; end

%% ---- trial order and list of conditions ----
cond = trialConds_(npRec);  cond = cond(:);
nExpTr = numel(cond);
ucL = unique(cond)';
fprintf('\nConditions in this session:\n');
for c = ucL
    fprintf('  %3d  %-34s %d trials\n', c, char(Q.cond(c).label), sum(cond == c));
end

%% ---- trigger times on the NP clock ----
niFile = fullfile(npRec,'continuous',niStream,'continuous.dat');
dN = dir(niFile);  NN = dN.bytes/2/8;
mN = memmapfile(niFile,'Format',{'int16',[8 NN],'a'});
x = single(mN.Data.a(trigAI+1,:)) * 0.0003051851;
th = (min(x(1:50:end)) + max(x(1:50:end))) / 2;
on = x > th;
up = find(~on(1:end-1) & on(2:end)) + 1;
clear x on
tsNI = readNpy_(fullfile(npRec,'continuous',niStream,'timestamps.npy'));
tTrig = double(tsNI(up));  tTrig = tTrig(:);  clear tsNI mN
nTr = min(numel(tTrig), nExpTr);
fprintf('\nTriggers %d | stim-file trials %d -> using %d\n', numel(tTrig), nExpTr, nTr);
tTrig = tTrig(1:nTr);  cond = cond(1:nTr);

%% ---- AP stream ----
oe = jsondecode(fileread(fullfile(npRec,'structure.oebin')));
cs = oe.continuous;  if iscell(cs), cs = [cs{:}]; end
iAP = find(arrayfun(@(s) contains(s.folder_name, apStream), cs), 1);
assert(~isempty(iAP), 'AP stream %s not found', apStream);
fsAP = cs(iAP).sample_rate;
chn = cs(iAP).channels;  if iscell(chn), chn = chn{1}; else, chn = chn(1); end
bV = chn.bit_volts;
fA = fullfile(npRec,'continuous',apStream,'continuous.dat');
dA = dir(fA);  N = dA.bytes/2/nCh;
m = memmapfile(fA,'Format',{'int16',[nCh N],'a'});
tsAP = readNpy_(fullfile(npRec,'continuous',apStream,'timestamps.npy'));
idx = round((tTrig - double(tsAP(1))) * fsAP) + 1;
clear tsAP

pre = round(preMs/1e3*fsAP);  post = round(postMs/1e3*fsAP);  nS = pre+post+1;
t = (-pre:post)/fsAP*1e3;
ok = (tTrig >= 0) & (idx-pre >= 1) & (idx+post <= N);
bIdx = find(t >= -preMs+3 & t < -2);
mergeS = round(mergeMs/1e3*fsAP);
valid = t(1:end-1) >= -2;
assert(exist('butter','file') > 0, 'butter not found: needs the Signal Processing Toolbox');
[bF, aF] = butter(3, bandHz/(fsAP/2));

%% ---- choose the channels that carry the artifact ----
elec = find(ok & nEcond(cond) > 0);
e1only = find(ok & hasE(cond,1) & ~hasE(cond,2));
e2only = find(ok & hasE(cond,2) & ~hasE(cond,1));
fprintf('Choosing channels from the largest-amplitude trials: electrode A only %d trials, B only %d, any electrical %d\n', numel(e1only), numel(e2only), numel(elec));
selList = [];
for g = {e1only, e2only}
    tr = g{1};
    if numel(tr) < 3, continue, end
    s = chanSnr_(m, idx, pickTr_(tr, cond, ampC, nSelTr), nCh, bV, pre, post, bIdx, valid);
    [v, o] = sort(s, 'descend');  o = o(v >= snrMin);
    selList = [selList; o(1:min(nEach, numel(o)))]; %#ok<AGROW>
end
if numel(unique(selList)) < 8
    s = chanSnr_(m, idx, pickTr_(elec, cond, ampC, nSelTr), nCh, bV, pre, post, bIdx, valid);
    [v, o] = sort(s, 'descend');  o = o(v >= snrMin);
    selList = [selList; o(1:min(2*nEach, numel(o)))];
end
sel = unique(selList, 'stable');
assert(numel(sel) >= 4, 'Fewer than 4 channels with a pulse step above %g x noise: lower snrMin or check the data', snrMin);
fprintf('Using %d channels: %s\n', numel(sel), mat2str(sort(sel(:))'));

%% ---- find the pulses of every trial ----
det = cell(nTr,1);
fprintf('Reading %d trials ...\n', sum(ok));
for k = 1:nTr
    if ~ok(k), continue, end
    X = double(m.Data.a(sel, idx(k)+(-pre:post))) * bV;
    R = X - median(X(:,bIdx), 2);
    db = diff(R(:,bIdx), 1, 2);                                    % steps before the trigger = this trial's noise
    sd = max(1.4826*median(abs(db - median(db,2)), 2), 0.5);       % per channel (uV)
    fr = mean(abs(diff(R,1,2)) > kThr*sd, 1);                      % fraction of channels stepping at each sample
    ia = find(fr >= fracMin & valid);
    if ~isempty(ia)
        br = [1, find(diff(ia) > mergeS) + 1];
        det{k} = t(ia(br)+1);                                      % time of the first step of each pulse
    end
    if mod(k,50) == 0, fprintf('  %d / %d\n', k, nTr); end
end

%% ---- constant offset trigger -> first pulse, then match expected and detected ----
offs = [];
for k = find(ok & nEcond(cond) > 0)'
    e = expT{cond(k)};  d = det{k};
    for i = 1:numel(e)
        if isempty(d), break, end
        [dm, j] = min(abs(d - e(i)));
        if dm <= matchMs, offs(end+1) = d(j) - e(i); end %#ok<AGROW>
    end
end
assert(~isempty(offs), 'No pulse was detected in any electrical trial: lower kThr / fracMin / snrMin');
offset = median(offs);
fprintf('\nPulse start relative to the trigger: median offset %.2f ms (from %d pulses)\n', offset, numel(offs));

nMiss = zeros(nTr,1);  nExtra = zeros(nTr,1);  errMax = nan(nTr,1);  missT = cell(nTr,1);  extraT = cell(nTr,1);
for k = find(ok)'
    e = expT{cond(k)} + offset;  d = det{k}(:)';
    hitD = false(size(d));  ms = false(size(e));  er = [];
    for i = 1:numel(e)
        if isempty(d), break, end
        [dm, j] = min(abs(d - e(i)));
        if dm <= tolMs, ms(i) = true; hitD(j) = true; er(end+1) = dm; end %#ok<AGROW>
    end
    nMiss(k) = sum(~ms);  nExtra(k) = sum(~hitD);
    missT{k} = e(~ms) - offset;  extraT{k} = d(~hitD) - offset;
    if ~isempty(er), errMax(k) = max(er); end
end
good = ok & nMiss == 0 & nExtra == 0;
nDet = cellfun(@numel, det);

%% ---- report per condition ----
uc = unique(cond)';
fprintf('\n%-4s %-34s %3s | %-8s | %-9s | %-10s | %-7s | %-5s | %s\n', 'cond', 'label', 'n', 'expected', 'detected', 'trials ok', 'missing', 'extra', 'worst error ms');
for c = uc
    s = find(ok & cond == c);
    if isempty(s), continue, end
    cnt = nDet(s);
    if all(cnt == cnt(1)), cs_ = sprintf('%d', cnt(1)); else, cs_ = sprintf('%d-%d', min(cnt), max(cnt)); end
    w = max(errMax(s));  if isnan(w), ws = '-'; else, ws = sprintf('%.2f', w); end
    flag = '';  if sum(good(s)) < numel(s), flag = '  <-- CHECK'; end
    fprintf('%-4d %-34s %3d | %-8d | %-9s | %3d of %-3d | %-7d | %-5d | %s%s\n', c, char(Q.cond(c).label), numel(s), nEcond(c), cs_, sum(good(s)), numel(s), sum(nMiss(s)), sum(nExtra(s)), ws, flag);
end
isEl = ok & nEcond(cond) > 0;  isCt = ok & nEcond(cond) == 0;
fprintf('\nElectrical trials with every expected pulse and no extra pulse: %d of %d\n', sum(good & isEl), sum(isEl));
fprintf('Visual / blank trials with NO pulse: %d of %d\n', sum(good & isCt), sum(isCt));
fprintf('Pulses missing %d of %d expected | extra pulses %d\n', sum(nMiss(isEl)), sum(nEcond(cond(isEl))), sum(nExtra(ok)));
bad = find(ok & ~good);
if ~isempty(bad)
    fprintf('\nTrials that do not match (trial, cond, label | missing at ms | extra at ms), first %d of %d:\n', min(40, numel(bad)), numel(bad));
    fm_ = @(v) ['[' strtrim(sprintf('%.1f ', v)) ']'];
    for q = bad(1:min(40,numel(bad)))'
        fprintf('  %3d  %-3d %-30s | %s | %s\n', q, cond(q), char(Q.cond(cond(q)).label), fm_(missT{q}), fm_(extraT{q}));
    end
end

%% ---- Fig 1: detected (black) vs expected (red) per condition ----
figure('Color','w','Position',[40 40 1100 1000]);  hold on;
[~, row] = ismember(cond, uc);
XX = [];  YY = [];
for k = find(ok)'
    nrep = sum(row == row(k) & ok);  rep = sum(row(1:k) == row(k) & ok(1:k));
    y = row(k) + ((rep-0.5)/max(nrep,1) - 0.5) * 0.7;
    XX = [XX, det{k}(:)']; YY = [YY, y*ones(1,numel(det{k}))]; %#ok<AGROW>
end
plot(XX, YY, 'k.', 'MarkerSize', 5);
for i = 1:numel(uc)
    e = expT{uc(i)};
    if ~isempty(e), plot(e + offset, i*ones(size(e)), 'ro', 'MarkerSize', 5); end
end
set(gca, 'YDir', 'reverse', 'YTick', 1:numel(uc), 'YTickLabel', arrayfun(@(c) char(Q.cond(c).label), uc, 'UniformOutput', false), 'TickLabelInterpreter', 'none', 'FontSize', 6);
ylim([0.3 numel(uc)+0.7]);  xlim(showMs);  grid on;  xlabel('ms from trigger');
title('Pulses found on NP (black) vs the design (red circles)');

%% ---- Fig 2+: traces of chosen conditions, one figure per channel ----
if isempty(plotConds)
    plotConds = [];
    for kd = {'sim', 'A->B', 'B->A'}
        best = 0;  bs = -1;
        for c = uc
            if startsWith(char(Q.cond(c).label), kd{1}) && nEcond(c) > 0
                sc = nEcond(c)*1000 + max(ampC(c,:));
                if sc > bs, bs = sc; best = c; end
            end
        end
        if best > 0, plotConds(end+1) = best; end %#ok<AGROW>
    end
    fprintf('\nFig 2 conditions (auto): %s\n', mat2str(plotConds));
end
if isempty(plotChan), plotChan = sel(1:min(3, numel(sel)))'; fprintf('Fig 2 channels (auto): %s\n', mat2str(plotChan)); end
pc = plotConds(ismember(plotConds, uc));
if isempty(pc), fprintf('No valid plotConds: Fig 2 skipped\n'); return, end
nP = numel(pc);  nPC = numel(plotChan);
RRq = cell(nP,1);  FFq = cell(nP,1);
for q = 1:nP
    s = find(ok & cond == pc(q));
    RRq{q} = nan(nS, numel(s), nPC);  FFq{q} = RRq{q};
    for i = 1:numel(s)
        X = double(m.Data.a(plotChan, idx(s(i))+(-pre:post))) * bV;
        R = X - median(X(:,bIdx), 2);
        F = filter(bF, aF, R, [], 2);
        RRq{q}(:,i,:) = R';  FFq{q}(:,i,:) = F';
    end
end
for c = 1:nPC
    figure('Color','w','Position',[60+25*c 60+25*c min(1700, 420*nP) 750]);
    for q = 1:nP
        subplot(2, nP, q);  plot(t, RRq{q}(:,:,c));  xlim(showMs);  xlabel('ms from trigger');  ylabel('uV');
        title(sprintf('%d: %s | raw, ch %d', pc(q), char(Q.cond(pc(q)).label), plotChan(c)), 'Interpreter', 'none', 'FontSize', 8);
        subplot(2, nP, nP+q);  plot(t, FFq{q}(:,:,c));  xlim(showMs);  xlabel('ms from trigger');  ylabel('uV');
        title('filtered 300-6000 Hz', 'FontSize', 9);
    end
    sgtitle(sprintf('All repetitions of each condition overlaid (channel %d)', plotChan(c)));
end

%% ======================== local functions ==============================
function tr = pickTr_(tr, cond, ampC, n)
[~, o] = sort(max(ampC(cond(tr),:), [], 2), 'descend');
tr = tr(o);  tr = tr(1:min(n, numel(tr)));
end

function snr = chanSnr_(m, idx, tr, nCh, bV, pre, post, bIdx, valid)
% median over trials of (largest single step of the channel) / (its own step noise before the trigger)
S = nan(nCh, numel(tr));
for i = 1:numel(tr)
    X = double(m.Data.a(:, idx(tr(i))+(-pre:post))) * bV;
    R = X - median(X(:,bIdx), 2);
    db = diff(R(:,bIdx), 1, 2);
    sd = max(1.4826*median(abs(db - median(db,2)), 2), 0.5);
    d = abs(diff(R,1,2));
    S(:,i) = max(d(:,valid), [], 2) ./ sd;
end
snr = median(S, 2);
end

function rec = resolveRec_(p)
if exist(fullfile(p,'structure.oebin'),'file'), rec = p; return, end
rr = fullfile(p,'Record Node 101','experiment1','recording1');
if exist(fullfile(rr,'structure.oebin'),'file'), rec = rr; return, end
c = dir(fullfile(p,'**','structure.oebin'));
assert(~isempty(c), 'No structure.oebin below %s', p);
assert(numel(c) == 1, 'Several recordings below %s', p);
rec = c(1).folder;
end

function cond = trialConds_(npRec)
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
if contains(descr,'S')
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