%% NP CHECK 3: were all trials delivered? trigger + photodiode per trial (script)
% Read-only. For ONE session. Uses only the 8 NI analog inputs (a few hundred MB),
% not the 15 GB AP file.
%
% Per trial it finds: the trigger (AI3, 100 ms pulse) on the NP clock, and the
% photodiode (AI0): number of frame pulses, first / last pulse relative to the
% trigger. Then it compares with the stim-file table:
%   - electrical and blank trials should have NO photodiode pulses
%   - visual-spot trials should have photodiode pulses (~60 at 120 Hz over 500 ms)
% Prints a summary and every disagreement, plots the photodiode, and (optionally)
% saves the per-trial table. Nothing else is written.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
npSession = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/2026_10_29/NP/Train_Amp_test2_expTrainAmps_135452_2026-10-29_13-55-01';
nnFolder  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_Amp_Test2_261029_135502';   % (a) folder with the *_exp_datafile_*.mat built from the stim log
condFile  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_Amp_Test2_261029_135502/Train_Amp_test2.expTrainAmps.135452.conds.mat';
niStream  = 'NI-DAQmx-104.PXI-6133';
fs        = 25000;          % NI analog sample rate
trigAI    = 3;              % AI3: trigger
pdAI      = 0;              % AI0: photodiode
preMs     = 100;  postMs = 800;
pdFrac    = 0.4;            % photodiode threshold = baseline + pdFrac x (max - baseline)
minPdEdges = 10;            % fewer pulses than this = "no photodiode on this trial"
saveTable = '';             % e.g. '/path/NP_trial_table_test2.mat'  ('' = do not save)
%% ================================

npRec = resolveRec_(npSession);
fprintf('Recording: %s\n', npRec);

%% ---- stim-file table ----
if ~isempty(condFile)
    % NP data only: condition per trial from the TRIALID messages, label from the conds file
    cond = trialConds_(npRec);
    Q = load(condFile);
    label = arrayfun(@(k) char(Q.cond(k).label), cond, 'UniformOutput', false);
    label = label(:);  cond = cond(:);
    fprintf('Condition labels: NP TRIALID messages + %s\n', condFile);
else
    f = dir(fullfile(nnFolder,'*_exp_datafile_*.mat'));
    assert(~isempty(f), 'No *_exp_datafile_*.mat in %s (or set condFile instead)', nnFolder);
    S = load(fullfile(nnFolder,f(1).name),'TrialTable');
    cond  = double([S.TrialTable.cond])';
    label = cellstr({S.TrialTable.label}');
    fprintf('Condition labels: %s\n', f(1).name);
end
nExp = numel(cond);
expectPD = contains(label, 'spot');
fprintf('Stim-file trials: %d (visual spot %d, blank %d, electrical %d)\n', nExp, sum(expectPD), sum(contains(label,'blank')), sum(~expectPD & ~contains(label,'blank')));

%% ---- NI analog data ----
fN = fullfile(npRec,'continuous',niStream,'continuous.dat');
d = dir(fN);  N = d.bytes/2/8;
m = memmapfile(fN,'Format',{'int16',[8 N],'a'});
V = 0.0003051851;

x = single(m.Data.a(trigAI+1,:)) * V;
th = (min(x(1:50:end)) + max(x(1:50:end))) / 2;
on = x > th;
up = find(~on(1:end-1) & on(2:end)) + 1;
clear x on
nTrig = numel(up);
fprintf('Triggers on AI%d: %d | stim-file trials: %d\n', trigAI, nTrig, nExp);

tsAll = readNpy_(fullfile(npRec,'continuous',niStream,'timestamps.npy'));
tTrig = double(tsAll(up));  clear tsAll
if any(tTrig < 0), warning('%d trigger(s) have an unsynchronised (-1) timestamp.', sum(tTrig < 0)); end
dT = diff(tTrig);
fprintf('Trigger start-to-start (s): median %.3f, min %.3f, max %.3f\n', median(dT), min(dT), max(dT));
gap = find(dT > 1.5*median(dT));
if isempty(gap), fprintf('  no missing triggers (no gap longer than 1.5 x median)\n');
else, fprintf('  GAPS after trigger(s) %s: possible missing trial(s)\n', mat2str(gap(:)')); end

%% ---- photodiode per trial ----
p = single(m.Data.a(pdAI+1, 1:50:end)) * V;
base = median(p);  hi = max(p);  thPD = base + pdFrac*(hi - base);   clear p
fprintf('Photodiode AI%d: baseline %.3f V, max %.3f V -> threshold %.3f V\n', pdAI, base, hi, thPD);

pre = round(preMs/1e3*fs);  post = round(postMs/1e3*fs);
n = min(nTrig, nExp);
nE = zeros(n,1);  first = nan(n,1);  last = nan(n,1);  ipi = nan(n,1);
W = nan(pre+post+1, 0, 'single');  keep = [];
for k = 1:n
    if up(k)-pre < 1 || up(k)+post > N, continue, end
    y = single(m.Data.a(pdAI+1, up(k)+(-pre:post))) * V;
    o = y > thPD;
    e = find(~o(1:end-1) & o(2:end)) + 1;
    e = e((e-pre-1) > -round(0.05*fs));                % ignore anything long before the trigger
    nE(k) = numel(e);
    if ~isempty(e)
        first(k) = (e(1)-pre-1)/fs*1e3;  last(k) = (e(end)-pre-1)/fs*1e3;
        if numel(e) > 1, ipi(k) = median(diff(e))/fs*1e3; end
        if expectPD(k) && numel(e) >= minPdEdges, W(:,end+1) = y(:); keep(end+1) = k; end %#ok<SAGROW>
    end
end
gotPD = nE >= minPdEdges;

%% ---- compare with the stim file ----
fprintf('\n=== Photodiode vs stim file ===\n');
fprintf('Visual-spot trials with photodiode pulses: %d of %d\n', sum(gotPD(1:n) & expectPD(1:n)), sum(expectPD(1:n)));
fprintf('Other trials (electrical + blank) with photodiode pulses: %d of %d\n', sum(gotPD(1:n) & ~expectPD(1:n)), sum(~expectPD(1:n)));
bad = find(gotPD(1:n) ~= expectPD(1:n));
if isempty(bad), fprintf('All %d trials agree with the stim file.\n', n);
else
    fprintf('DISAGREEMENTS (trial, label, pulses found):\n');
    for q = bad(:)', fprintf('  %3d  %-30s %d\n', q, label{q}, nE(q)); end
end
v = find(expectPD(1:n) & gotPD(1:n));
if ~isempty(v)
    fprintf('\nVisual trials (n=%d):\n', numel(v));
    fprintf('  pulses per trial: %s\n', strjoin(arrayfun(@(a) sprintf('%d x%d', a, sum(nE(v)==a)), unique(nE(v))', 'UniformOutput', false), ', '));
    fprintf('  first pulse after trigger (ms): median %.1f, min %.1f, max %.1f\n', median(first(v)), min(first(v)), max(first(v)));
    fprintf('  last pulse after trigger (ms):  median %.1f, min %.1f, max %.1f\n', median(last(v)), min(last(v)), max(last(v)));
    fprintf('  pulse spacing (ms): median %.2f (120 Hz = 8.33)\n', median(ipi(v)));
    long = v(nE(v) ~= mode(nE(v)));
    if ~isempty(long), fprintf('  trials with a different pulse count than the most common (%d): %s\n', mode(nE(v)), mat2str(long(:)')); end
end

T = table((1:n)', cond(1:n), label(1:n), tTrig(1:n), nE, first, last, ipi, expectPD(1:n), gotPD(1:n), ...
    'VariableNames', {'trial','cond','label','trigTime_s','pdPulses','pdFirst_ms','pdLast_ms','pdSpacing_ms','expectPD','gotPD'});
if ~isempty(saveTable), save(saveTable, 'T'); fprintf('\nTable saved to %s\n', saveTable); end
fprintf('\nThe per-trial table is in the variable T.\n');

%% ---- plots ----
figure('Color','w','Position',[60 60 1400 800]);
subplot(2,2,1);
plot(2:n, dT(1:n-1), '.-'); xlabel('trial'); ylabel('s'); title('Trigger start-to-start interval'); ylim([0 max(3, 1.1*max(dT))]);
subplot(2,2,2); hold on;
tt = (-pre:post)/fs*1e3;
if ~isempty(keep)
    plot(tt, W, 'Color', [0.75 0.75 0.75]);  plot(tt, mean(W,2), 'k', 'LineWidth', 1.2);
    plot([0 0], ylim, 'r:');  xlim([-20 600]);
end
xlabel('ms from trigger'); ylabel('V'); title(sprintf('Photodiode, %d visual-spot trials', numel(keep)));
subplot(2,2,3);
histogram(first(v), 20); xlabel('first photodiode pulse after trigger (ms)'); ylabel('trials'); title('Onset latency (visual trials)');
subplot(2,2,4);
plot(1:n, nE, '.'); hold on; plot(find(expectPD(1:n)), nE(expectPD(1:n)), 'ro', 'MarkerSize', 4);
xlabel('trial'); ylabel('photodiode pulses'); title('Pulses per trial (red = visual-spot trials in stim file)');

%% ======================== local functions ==============================
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