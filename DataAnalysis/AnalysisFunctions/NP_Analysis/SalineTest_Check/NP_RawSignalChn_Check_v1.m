%% NP VIEW TRIALS: look at chosen trials of a Neuropixels (AP) recording (script)
% General trial viewer, not tied to the saline tests. Read-only, nothing is saved.
% For every chosen trial (number = order of the NI trigger, as in the other np_check scripts) it shows
%   top left : raw, baseline-subtracted traces of the chosen channels, stacked (clipped at +-stackUV/2),
%              a same-condition reference trial in grey behind them (if the condition is known)
%   top right: the same channels filtered 300-6000 Hz
%   bottom   : heat map of ALL channels (raw, baseline-subtracted, colour clipped at +-heatUV)
%              a noise burst or baseline jump shows as a band across many channels
% and prints, per trial, three numbers that tell whether the trial itself is disturbed, next to the
% same numbers for ~40 typical trials spread over the whole recording:
%   noise    : step noise of the pre-trigger signal (median over channels, uV) -> sets detection thresholds
%   wander   : range of the common-mode (median over channels) before the trigger, uV
%   end off. : common-mode offset at the end of the window relative to the baseline, uV
% condFile is OPTIONAL: with it the figure titles carry the condition label and the expected pulse
% times are marked (red triangles); without it only the trial number is shown.
% To view windows around other times (spikes, events) instead of triggers, set centerSec.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
npSession = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/2026_10_29/NP/Train_Freq_test1_expTrainFreq_200Hz_141342_2026-10-29_14-13-49';
condFile  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_Freq_Test1_261029_141350/Train_Freq_test1.expTrainFreq_200Hz.141342.conds.mat';
niStream  = 'NI-DAQmx-104.PXI-6133';
apStream  = 'Neuropix-PXI-106.ProbeA-AP';
trigAI    = 3;               % NI analog input with the trigger
nCh       = 384; 
trials    = [1:10];   % trial numbers to look at 8 85 135 172 144 93
centerSec = [];              % [] = use the triggers; or NP-clock times (s) to centre the windows on instead
chans     = [321 325 323];              % AP channels for the traces; [] = the nAuto channels with the largest step in these trials
nAuto     = 3;
winMs     = [-20 120];       % window around the trigger (the pre part is used as the baseline: needs at least ~10 ms)
offsetMs  = 0.3;             % trigger -> first pulse delay, only for the red expected-pulse marks
stackUV   = 600;             % spacing of the stacked traces (traces are clipped at +-stackUV/2)
heatUV    = 300;             % colour range of the heat map
showRef   = true;            % grey reference trial of the same condition behind the raw traces
nTyp      = 40;              % trials used for the "typical" numbers
bandHz    = [300 6000];
%% ================================

npRec = resolveRec_(npSession);
fprintf('Recording: %s\n', npRec);

%% ---- optional design / labels ----
hasCond = ~isempty(condFile);
cond = [];  expT = {};  Q = [];
if hasCond
    Q = load(condFile);  a = Q.args;
    periodMs = 1e3 / double(a.freq);
    expT = cell(numel(Q.cond),1);
    for c = 1:numel(Q.cond)
        v = [];
        for e = 1:2
            if double(Q.cond(c).amp(e)) > 0 && double(Q.cond(c).nsp(e)) > 0
                v = [v, double(Q.cond(c).ptd(e))/1e3 + (0:double(Q.cond(c).nsp(e))-1)*periodMs]; %#ok<AGROW>
            end
        end
        v = sort(v);  if ~isempty(v), v = v([true, diff(v) > 0.05]); end
        expT{c} = v; %#ok<SAGROW>
    end
    try, cond = trialConds_(npRec);  cond = cond(:); catch, cond = []; end
end

%% ---- times to centre the windows on (NP clock) ----
if isempty(centerSec)
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
    fprintf('Triggers: %d\n', numel(tTrig));
else
    tTrig = double(centerSec(:));
end
nTr = numel(tTrig);
if ~isempty(cond) && numel(cond) ~= nTr
    fprintf('!! %d triggers but %d trial messages: labels use the first %d\n', nTr, numel(cond), min(nTr, numel(cond)));
end
assert(all(trials >= 1 & trials <= nTr), 'trials must be between 1 and %d', nTr);

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

pre = round(-winMs(1)/1e3*fsAP);  post = round(winMs(2)/1e3*fsAP);  nS = pre+post+1;
t = (-pre:post)/fsAP*1e3;
ok = (tTrig >= 0) & (idx-pre >= 1) & (idx+post <= N);
bIdx = find(t >= winMs(1)+3 & t < -2);
assert(numel(bIdx) > 50, 'Need at least ~2 ms of baseline: make winMs(1) more negative');
assert(all(ok(trials)), 'These trials do not fit in the recording: %s', mat2str(trials(~ok(trials))));
assert(exist('butter','file') > 0, 'butter not found: needs the Signal Processing Toolbox');
[bF, aF] = butter(3, bandHz/(fsAP/2));

%% ---- typical trials (numbers to compare with) ----
pool = find(ok)';
tp = pool(unique(round(linspace(1, numel(pool), min(nTyp, numel(pool))))));
S = nan(numel(tp), 3);
fprintf('Reading %d typical trials ...\n', numel(tp));
for i = 1:numel(tp)
    R = readTrial_(m, idx(tp(i)), pre, post, bV, bIdx);
    [S(i,1), S(i,2), S(i,3)] = trialStats_(R, t, bIdx);
end
typ = [median(S,1); arrayfun(@(j) pct_(S(:,j), 0.05), 1:3); arrayfun(@(j) pct_(S(:,j), 0.95), 1:3)];
fprintf('\nTypical trials (n=%d): median [5%%-95%%]\n', numel(tp));
fprintf('  noise    %.1f uV  [%.1f - %.1f]\n', typ(1,1), typ(2,1), typ(3,1));
fprintf('  wander   %.1f uV  [%.1f - %.1f]\n', typ(1,2), typ(2,2), typ(3,2));
fprintf('  end off. %.1f uV  [%.1f - %.1f]   (absolute value used for the flag)\n', typ(1,3), typ(2,3), typ(3,3));

%% ---- read the chosen trials ----
nT = numel(trials);
Rall = zeros(nCh, nS, nT, 'single');
for i = 1:nT
    Rall(:,:,i) = single(readTrial_(m, idx(trials(i)), pre, post, bV, bIdx));
end
if isempty(chans)
    pk = zeros(nCh, nT);
    for i = 1:nT, pk(:,i) = max(abs(diff(double(Rall(:,:,i)),1,2)), [], 2); end
    [~, o] = sort(median(pk,2), 'descend');
    chans = o(1:min(nAuto, nCh))';
    fprintf('\nChannels (auto): %s\n', mat2str(chans));
end
nC = numel(chans);

%% ---- per-trial report ----
fprintf('\n%-6s %-5s %-30s | %-12s | %-12s | %-12s\n', 'trial', 'cond', 'label', 'noise uV', 'wander uV', 'end off. uV');
for i = 1:nT
    k = trials(i);
    [sd, wd, eo] = trialStats_(double(Rall(:,:,i)), t, bIdx);
    cn = '';  lb = '';
    if ~isempty(cond) && k <= numel(cond), cn = sprintf('%d', cond(k)); if hasCond, lb = char(Q.cond(cond(k)).label); end, end
    fl = '';
    if sd > typ(3,1), fl = [fl ' NOISY']; end %#ok<AGROW>
    if wd > typ(3,2), fl = [fl ' WANDER']; end %#ok<AGROW>
    if abs(eo) > max(abs(typ(2:3,3))), fl = [fl ' OFFSET']; end %#ok<AGROW>
    fprintf('%-6d %-5s %-30s | %-12.1f | %-12.1f | %-12.1f |%s\n', k, cn, lb, sd, wd, eo, fl);
end
fprintf('(flags: above the 95%% value of the typical trials)\n');

%% ---- figures ----
cm = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
half = stackUV/2;
for i = 1:nT
    k = trials(i);
    R = double(Rall(:,:,i));
    Rc = R(chans,:);
    Fc = filter(bF, aF, Rc, [], 2);
    ttl = sprintf('Trial %d', k);
    e = [];  ref = [];
    if ~isempty(cond) && k <= numel(cond)
        ttl = sprintf('Trial %d | cond %d', k, cond(k));
        if hasCond, ttl = sprintf('%s: %s', ttl, char(Q.cond(cond(k)).label)); e = expT{cond(k)}; end
        if showRef
            same = find(cond(1:min(numel(cond),nTr)) == cond(k) & ok(1:min(numel(cond),nTr)));
            same = same(~ismember(same, trials));
            if ~isempty(same)
                ref = readTrial_(m, idx(same(1)), pre, post, bV, bIdx);  ref = ref(chans,:);
                ttl = sprintf('%s | grey = trial %d (same condition)', ttl, same(1));
            end
        end
    end
    figure('Color','w','Position',[40+20*i 40+20*i 1300 850]);
    subplot(2,2,1);  hold on;
    for j = 1:nC
        y0 = (j-1)*stackUV;
        if ~isempty(ref), plot(t, y0 + max(min(ref(j,:), half), -half), 'Color', [0.75 0.75 0.75]); end
        plot(t, y0 + max(min(Rc(j,:), half), -half), 'k');
    end
    if ~isempty(e), plot(e + offsetMs, (nC-0.35)*stackUV*ones(size(e)), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 4); end
    xline(0, ':');  xlim(winMs);  ylim([-stackUV, nC*stackUV]);
    set(gca, 'YTick', (0:nC-1)*stackUV, 'YTickLabel', arrayfun(@(c) sprintf('ch %d', c), chans, 'UniformOutput', false));
    xlabel('ms from trigger');  title(sprintf('raw, baseline-subtracted (clipped at \\pm%g uV)', half));
    subplot(2,2,2);  hold on;
    for j = 1:nC
        y0 = (j-1)*stackUV;
        plot(t, y0 + max(min(Fc(j,:), half), -half), 'k');
    end
    xline(0, ':');  xlim(winMs);  ylim([-stackUV, nC*stackUV]);
    set(gca, 'YTick', (0:nC-1)*stackUV, 'YTickLabel', arrayfun(@(c) sprintf('ch %d', c), chans, 'UniformOutput', false));
    xlabel('ms from trigger');  title(sprintf('filtered %d-%d Hz', bandHz(1), bandHz(2)));
    subplot(2,1,2);
    imagesc(t, 1:nCh, R);  set(gca, 'YDir', 'normal');  caxis([-heatUV heatUV]);  colormap(gca, cm);
    cb = colorbar;  cb.Label.String = 'uV';  hold on;
    plot(t(1)*ones(size(chans)) + 0.01*(t(end)-t(1)), chans, '>k', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
    xline(0, ':');  xlabel('ms from trigger');  ylabel('AP channel');
    title(sprintf('all channels, raw baseline-subtracted (colour \\pm%g uV); yellow = channels above', heatUV));
    sgtitle(ttl, 'Interpreter', 'none');
end

%% ======================== local functions ==============================
function R = readTrial_(m, i0, pre, post, bV, bIdx)
X = double(m.Data.a(:, i0+(-pre:post))) * bV;
R = X - median(X(:,bIdx), 2);
end

function [sdStep, wander, endOff] = trialStats_(R, t, bIdx)
db = diff(R(:,bIdx), 1, 2);
sdStep = median(1.4826*median(abs(db - median(db,2)), 2));
cmn = median(R(:,bIdx), 1);
wander = max(cmn) - min(cmn);
ie = t >= t(end) - 5;
endOff = median(mean(R(:,ie), 2));
end

function v = pct_(x, p)
x = sort(x(~isnan(x)));
v = x(max(1, min(numel(x), round(p*numel(x)))));
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