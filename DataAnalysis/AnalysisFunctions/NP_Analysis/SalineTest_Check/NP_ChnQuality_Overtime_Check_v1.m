%% NP QUALITY OVER TIME: is the NP signal healthy in every trial, and when does it go bad? (script)
% General recording-quality check, not tied to the saline tests. Read-only.
% For EVERY trigger of every listed run it cuts -preMs..+postMs around the trigger on all 384 AP
% channels and measures, from the signal itself (no stimulation design needed):
%   noise    : median over channels of the step noise before the trigger (uV). ~10 = normal, ~0 = no signal
%   flat pre : share of channels with no signal before the trigger (step noise < 1 uV)
%   flat end : share of channels with no signal in the last 20 ms of the window
%              (flat at the end = amplifiers run to their limit during the trial)
%   wander   : range of the common signal (median over channels) before the trigger, uV
%   end off. : common offset at the end of the window relative to the baseline, uV
% Time axis = clock time (start time from the recording folder name + trigger time), so runs can be
% compared and an event in the room shows up in all of them at the same moment.
% Prints per run: trials, share of bad trials, and the bad share in 4 equal parts of the run.
% Fig 1: noise / flat pre / flat end for all trials against clock time, bad trials circled in red.
% Fig 2: bad-trial share per quarter of each run.
% A trial is BAD when flat pre or flat end is above badFlat.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
npRoot  = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/2026_10_29/NP';
runDirs = { ...     % folder names (or full recording1 paths) in recording order
    'Train_Amp_test1_expTrainAmps_133554_2026-10-29_13-36-51', ...
    'Train_Amp_test2_expTrainAmps_135452_2026-10-29_13-55-01', ...
    'Train_Freq_test1_expTrainFreq_200Hz_141342_2026-10-29_14-13-49', ...
    'Train_Freq_test2_expTrainFreq_50Hz_142927_2026-10-29_14-29-33', ...
    'Train_ISI_test1_expTrainISI_144448_2026-10-29_14-44-55', ...
    'Train_ISI_test2_expTrainISI_145939_2026-10-29_14-59-46'};
niStream = 'NI-DAQmx-104.PXI-6133';
apStream = 'Neuropix-PXI-106.ProbeA-AP';
trigAI   = 3;
nCh      = 384;
preMs    = 20;  postMs = 120;
badFlat  = 0.25;      % a trial is BAD when more than this share of channels has no signal (before the trigger or at the end)
saveMat  = '';        % e.g. '/Volumes/.../np_quality_2026_10_29.mat' to keep the numbers (the reading is slow); '' = do not save
%% ================================

nR = numel(runDirs);
Res = cell(nR,1);  lab = cell(nR,1);
for r = 1:nR
    p = runDirs{r};
    if ~startsWith(p, filesep), p = fullfile(npRoot, p); end
    npRec = resolveRec_(p);
    tk = regexp(runDirs{r}, '^(.*?)_exp', 'tokens', 'once');
    if isempty(tk), [~, lab{r}] = fileparts(runDirs{r}); else, lab{r} = tk{1}; end
    fprintf('\n[%d/%d] %s\n', r, nR, lab{r});
    Res{r} = runQuality_(npRec, niStream, apStream, trigAI, nCh, preMs, postMs);
    Res{r}.bad = Res{r}.flatPre > badFlat | Res{r}.flatEnd > badFlat;
end
if ~isempty(saveMat), save(saveMat, 'Res', 'lab', 'runDirs', 'badFlat'); fprintf('Saved %s\n', saveMat); end

%% ---- report ----
fprintf('\n%-26s %6s %8s | bad share in the 4 quarters of the run (first -> last)\n', 'run', 'trials', 'bad');
Q = nan(nR,4);
for r = 1:nR
    b = Res{r}.bad;  v = ~isnan(Res{r}.noise);  n = numel(b);
    q = ceil((1:n)/n*4);
    for j = 1:4, Q(r,j) = mean(b(q == j & v)); end
    fprintf('%-26s %6d %7.0f%% | %s\n', lab{r}, sum(v), 100*mean(b(v)), sprintf('%4.0f%% ', 100*Q(r,:)));
end
fprintf('\nPer run: typical noise %s uV\n', mat2str(round(cellfun(@(s) median(s.noise,'omitnan'), Res)', 1)));

%% ---- Fig 1: against clock time ----
cols = lines(nR);
figure('Color','w','Position',[40 40 1300 850]);
nm = {'noise (uV)', 'flat chns before trigger (share)', 'flat chns at the end of win(share)'};
fld = {'noise', 'flatPre', 'flatEnd'};
for s = 1:3
    ax = subplot(3,1,s);  hold on;
    for r = 1:nR
        x = Res{r}.clock;  y = Res{r}.(fld{s});
        plot(x, y, '.', 'Color', cols(r,:), 'MarkerSize', 7);
        b = Res{r}.bad;
        plot(x(b), y(b), 'o', 'Color', [0.85 0 0], 'MarkerSize', 5);
    end
    ylabel(nm{s});  grid on;
    if s == 1, title('Every trial of every run (red circle = bad trial)'); end
    if s == 3, xlabel('clock time'); end
    if s > 1, ylim([-0.03 1.03]); end
    for r = 1:nR, xline(Res{r}.clock(1), ':', lab{r}, 'LabelVerticalAlignment', 'bottom', 'FontSize', 7, 'Interpreter', 'none'); end
    if s == 1, ax1 = ax; else, linkaxes([ax1 ax], 'x'); end
end

%% ---- Fig 2: bad share per quarter ----
figure('Color','w','Position',[80 80 900 450]);
bar(100*Q);  set(gca, 'XTick', 1:nR, 'XTickLabel', lab, 'TickLabelInterpreter', 'none');  xtickangle(25);
ylabel('bad trials (%)');  legend({'1st quarter', '2nd', '3rd', '4th'}, 'Location', 'northeast');
title('Share of bad trials in each quarter of each run');  grid on;

%% ======================== local functions ==============================
function Rn = runQuality_(npRec, niStream, apStream, trigAI, nCh, preMs, postMs)
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
nTr = numel(tTrig);
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
secFromStart = tTrig - double(tsAP(1));  clear tsAP
tk = regexp(npRec, '(\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})', 'tokens');
if isempty(tk)
    warning('No start time in the folder name: clock time = seconds from the start of this run');
    clk = datetime(2000,1,1) + seconds(secFromStart);
else
    t0 = datetime(tk{end}{1}, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss');
    clk = t0 + seconds(secFromStart);
end
pre = round(preMs/1e3*fsAP);  post = round(postMs/1e3*fsAP);
t = (-pre:post)/fsAP*1e3;
bIdx = find(t >= -preMs+3 & t < -2);
eIdx = find(t >= t(end)-20);
e5 = t >= t(end)-5;
ok = (tTrig >= 0) & (idx-pre >= 1) & (idx+post <= N);
Rn.clock = clk(:)';  Rn.noise = nan(1,nTr);  Rn.flatPre = nan(1,nTr);  Rn.flatEnd = nan(1,nTr);
Rn.wander = nan(1,nTr);  Rn.endOff = nan(1,nTr);
fprintf('  %d triggers, reading ...\n', nTr);
for k = 1:nTr
    if ~ok(k), continue, end
    X = double(m.Data.a(:, idx(k)+(-pre:post))) * bV;
    R = X - median(X(:,bIdx), 2);
    db = diff(R(:,bIdx), 1, 2);
    sd = 1.4826*median(abs(db - median(db,2)), 2);
    Rn.noise(k) = median(sd);  Rn.flatPre(k) = mean(sd < 1);
    cm = median(R(:,bIdx), 1);  Rn.wander(k) = max(cm) - min(cm);
    de = diff(R(:,eIdx), 1, 2);
    se = 1.4826*median(abs(de - median(de,2)), 2);
    Rn.flatEnd(k) = mean(se < 1);
    Rn.endOff(k) = median(mean(R(:,e5), 2));
    if mod(k,50) == 0, fprintf('  %d / %d\n', k, nTr); end
end
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