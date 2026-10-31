%% NP CHECK 5: is every condition delivered correctly? artifact pulses per condition (script)
% Read-only. For every trial it cuts -20..+120 ms around the trigger on all 384 AP
% channels, subtracts the pre-trigger baseline, and finds the stimulation pulses as
% bursts in the median |signal| over channels (the artifact is on all channels).
%
% Prints, for EVERY condition in the stim file: label, number of trials, pulses found
% per trial, pulse start times (ms from trigger, first trial) and the rep-to-rep jitter.
% Fig 1: raster of detected pulse starts, one row per condition (all conditions at once)
% Fig 2+: raw and filtered trace of the chosen conditions, all reps overlaid, one figure per chosen channel
% Nothing is saved.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
npSession = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/2026_10_29/NP/Train_Freq_test1_expTrainFreq_200Hz_141342_2026-10-29_14-13-49';
condFile  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_Freq_Test1_261029_141350/Train_Freq_test1.expTrainFreq_200Hz.141342.conds.mat';
niStream  = 'NI-DAQmx-104.PXI-6133';
apStream  = 'Neuropix-PXI-106.ProbeA-AP';
trigAI    = 3;
nCh       = 384;
preMs     = 20;  postMs = 120;
evThrUV   = 35;        % a pulse = median |signal| over channels above this (uV)
mergeMs   = 1.5;        % bursts closer than this are one pulse (the two phases of a biphasic pulse)
% Fig 2: which conditions to show = the numbers in the "cond" column of the printed list
% (1-5 sim @50, 6-15 A->B @50, 16-25 B->A @50, 26-30 sim @75, 31-40 A->B @75, 41-50 B->A @75,
%  51 blank, 52-55 visual). Default: the longest of each type at 50 uA, then sim 5/5, A->B 10/10, B->A 10/10 at 75 uA.
plotConds = [41:50];
plotChan  = [203];   % channels for Fig 2: ONE FIGURE PER CHANNEL (e.g. [20 100 233] or just 233)
showMs    = [-5 60];
bandHz    = [300 6000];
%% ================================

npRec = resolveRec_(npSession);
fprintf('Recording: %s\n', npRec);

%% ---- labels ----
cond = trialConds_(npRec);  Q = load(condFile);
label = arrayfun(@(k) char(Q.cond(k).label), cond, 'UniformOutput', false);  label = label(:);
nExp = numel(cond);

%% ---- list of conditions (pick numbers for plotConds from here) ----
ucL = unique(cond(:))';
fprintf('\nConditions in this session (use these numbers in plotConds):\n');
for ci = 1:numel(ucL)
    j = find(cond == ucL(ci));
    fprintf('  %3d  %-34s %d trials\n', ucL(ci), label{j(1)}, numel(j));
end
fprintf('Currently plotting: %s\n', mat2str(plotConds));

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
tTrig = double(tsNI(up));  clear tsNI mN
nTr = min(numel(tTrig), nExp);
fprintf('Triggers %d | stim-file trials %d -> using %d\n', numel(tTrig), nExp, nTr);
tTrig = tTrig(1:nTr);  cond = cond(1:nTr);  label = label(1:nTr);

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
assert(exist('butter','file') > 0, 'butter not found: needs the Signal Processing Toolbox');
[bF, aF] = butter(3, bandHz/(fsAP/2));

%% ---- pass 1: find the pulses of every trial ----
evT = cell(nTr,1);  evW = cell(nTr,1);
fprintf('Reading %d trials ...\n', sum(ok));
for k = 1:nTr
    if ~ok(k), continue, end
    X = double(m.Data.a(:, idx(k)+(-pre:post))) * bV;
    R = X - median(X(:,bIdx), 2);
    env = median(abs(R), 1);
    ia = find(env > evThrUV & t >= -2);
    if ~isempty(ia)
        br = [1, find(diff(ia) > mergeS) + 1];
        s = ia(br);  e = ia([br(2:end)-1, numel(ia)]);
        ts_ = t(s);  te_ = t(e);
        evT{k} = ts_(:)';  evW{k} = (te_(:) - ts_(:))';
    end
    if mod(k,50) == 0, fprintf('  %d / %d\n', k, nTr); end
end

%% ---- per-condition report ----
uc = unique(cond(:))';  nC = numel(uc);
fprintf('\n%-4s %-34s %3s | %-14s | %-9s | pulse starts in the first trial (ms from trigger)\n', 'cond', 'label', 'n', 'pulses/trial', 'jitter ms');
condLab = cell(nC,1);  rowOf = zeros(nTr,1);  repOf = zeros(nTr,1);
for ci = 1:nC
    s = find(ok & cond(:) == uc(ci));
    condLab{ci} = label{find(cond == uc(ci), 1)};
    rowOf(s) = ci;  repOf(s) = 1:numel(s);
    if isempty(s), continue, end
    cnt = cellfun(@numel, evT(s));
    if all(cnt == cnt(1)) && cnt(1) > 0
        D = cell2mat(cellfun(@(v) v(:)', evT(s), 'UniformOutput', false));
        jit = max(max(abs(D - D(1,:))));  js = sprintf('%.2f', jit);
    elseif all(cnt == 0), js = '-';
    else, js = 'COUNT VARIES';
    end
    if all(cnt == cnt(1)), cs_ = sprintf('%d', cnt(1)); else, cs_ = sprintf('%d-%d', min(cnt), max(cnt)); end
    v = evT{s(1)};
    fprintf('%-4d %-34s %3d | %-14s | %-9s | %s\n', uc(ci), condLab{ci}, numel(s), cs_, js, strtrim(sprintf('%.1f ', v(1:min(numel(v),14)))));
end

%% ---- Fig 1: raster of pulse starts, one row per condition ----
figure('Color','w','Position',[40 40 1100 1000]);  hold on;
XX = [];  YY = [];
for k = find(ok)'
    if isempty(evT{k}), continue, end
    nrep = sum(rowOf == rowOf(k));
    XX = [XX, evT{k}]; YY = [YY, rowOf(k) + ((repOf(k)-0.5)/max(nrep,1) - 0.5) * 0.7 + zeros(1, numel(evT{k}))]; %#ok<AGROW>
end
plot(XX, YY, 'k.', 'MarkerSize', 5);
set(gca, 'YDir', 'reverse', 'YTick', 1:nC, 'YTickLabel', condLab, 'TickLabelInterpreter', 'none', 'FontSize', 6);
ylim([0.3 nC+0.7]);  xlim(showMs);  grid on;  xlabel('ms from trigger');
title('Detected stimulation pulses per condition (each row = one condition, the 5 reps are the small offsets within the row)');

%% ---- Fig 2: chosen conditions on one channel, all reps overlaid ----
pc = [];
for q = 1:numel(plotConds)
    j = find(uc == plotConds(q), 1);
    if isempty(j), warning('No condition number %d', plotConds(q)); else, pc(end+1) = j; end %#ok<AGROW>
end
assert(~isempty(pc), 'plotConds has no valid condition number (see the "cond" column of the printed list)');
nP = numel(pc);  nPC = numel(plotChan);
RRq = cell(nP,1);  FFq = cell(nP,1);
for q = 1:nP                                           % read every chosen trial once, keep all chosen channels
    s = find(ok & rowOf == pc(q));
    RRq{q} = nan(nS, numel(s), nPC);  FFq{q} = RRq{q};
    for i = 1:numel(s)
        X = double(m.Data.a(:, idx(s(i))+(-pre:post))) * bV;
        R = X - median(X(:,bIdx), 2);
        F = filter(bF, aF, R, [], 2);
        RRq{q}(:,i,:) = R(plotChan,:)';  FFq{q}(:,i,:) = F(plotChan,:)';
    end
end
for c = 1:nPC                                          % one figure per channel
    figure('Color','w','Position',[60+25*c 60+25*c min(1700, 420*nP) 750]);
    for q = 1:nP
        subplot(2, nP, q);  plot(t, RRq{q}(:,:,c));  xlim(showMs);  xlabel('ms from trigger');  ylabel('uV');
        title(sprintf('%d: %s | raw, ch %d', uc(pc(q)), condLab{pc(q)}, plotChan(c)), 'Interpreter', 'none', 'FontSize', 8);
        subplot(2, nP, nP+q);  plot(t, FFq{q}(:,:,c));  xlim(showMs);  xlabel('ms from trigger');  ylabel('uV');
        title('filtered 300-6000 Hz', 'FontSize', 9);
    end
    sgtitle(sprintf('All repetitions of each condition overlaid (channel %d)', plotChan(c)));
end

%% ======================== local functions ==============================
function rec = resolveRec_(p)
if exist(fullfile(p,'structure.oebin'),'file'), rec = p; return, end
c = dir(fullfile(p,'**','structure.oebin'));
assert(~isempty(c), 'No structure.oebin below %s', p);
assert(numel(c) == 1, 'Several recordings below %s:\n%s', p, strjoin(fullfile({c.folder},{c.name}), newline));
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