%% ============================================================
% NN_FirstSpikeTime_v1.m
% Step 1 of the PULSE-TRAIN artifact recovery pipeline (marmoset NN data).
%
% For every ELECTRICAL trial, find the first spike that appears after the
% start of that trial's LAST pulse, i.e. where real recording resumes after
% the last artifact. Step 2 (the recovery script) uses it as the end of the
% window that is blanked and filled from the level with one pulse fewer.
%
% It also works out the LEVEL STRUCTURE of the design, from the conds file:
% each condition is a set of pulses (electrode, time, amplitude). Its SOURCE
% condition is the one whose pulse set equals this one minus its latest
% pulse (for a simultaneous level, minus both pulses at the latest time).
% Level 1 of every pattern has no source. The result is printed as a table
% and saved, so step 2 uses exactly the same chain.
%
% READ-ONLY with respect to the spike data: sp_clipped is only inspected.
% Blank and visual conditions are skipped (left NaN), and so are level-1 (single first
% pulse) conditions: they are the starting set that recovery copies from, never recovered.
%
% Inputs  (in data_folder unless entered):
%   times_file   FILTERED spike times (sp_clipped), entered directly
%   *.conds.mat  condition table saved by the stimulus PC (cond, info)
%   *_exp_datafile_*.mat   TrialParams / StimParams for the Intan folder
%   *.trig.dat   trigger times
% Output (next to times_file):
%   <base>.sp_xia_FirstSpikeTimes.mat
%       firstSpikeTimes{ch}(trial)  ms after the trigger (NaN = none found)
%       hasSpike(ch,trial)          1 if a spike was found
%       lastPulse_ms(trial)         start of the last pulse of that trial, ms
%       trialCond(trial)            condition number of every trial
%       condTable                   per condition: kind, level, source, last pulse ...
%       win_ms, FS, n_Trials
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
analysisPath = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/';
                    % added to the MATLAB path (needed for cleanTrig_sabquick only)

data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
                    % folder holding this dataset's trigger (*.trig.dat) and
                    % stim-protocol (*_exp_datafile_*.mat) files.

times_file = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_filtered.mat';
                    % the FILTERED spike-times file (sp_clipped) of this same dataset.

conds_file = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/CJ287_15Pulses3_100Hz.expTrainFreq_100Hz.080922.conds.mat';    % the *.conds.mat saved by the stimulus PC. '' = look for one inside data_folder.
                    % If it is not there, copy it into data_folder or enter its full path here.

win_ms = 10;        % how far after the start of the last pulse to look for the first spike (ms)
FS     = 30000;     % sampling rate (overridden below if saved in the times file)

%% ================= LOAD DATA =================
if ~isempty(analysisPath), addpath(genpath(analysisPath)); end
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

assert(isfile(times_file), 'Cannot find %s.', times_file);
S_t = load(times_file);
assert(isfield(S_t, 'sp_clipped'), 'Variable "sp_clipped" not found in %s.', times_file);
sp_clipped = S_t.sp_clipped;
if isfield(S_t, 'fs'), FS = S_t.fs; end
nChn = numel(sp_clipped);

% triggers (absolute samples). dir() without macOS "._" junk files.
f = listClean_('*.trig.dat');
if isempty(f)
    cleanTrig_sabquick;
    f = listClean_('*.trig.dat');
end
assert(~isempty(f), 'No *.trig.dat in %s.', data_folder);
ft = fopen(f(1).name, 'r');
trig = fread(ft, [1 inf], 'double');
fclose(ft);

% exp file: condition of every trial, pulse parameters
f = listClean_('*_exp_datafile_*.mat');
assert(~isempty(f), 'No *_exp_datafile_*.mat found in %s.', data_folder);
E = load(f(1).name, 'TrialParams', 'StimParams', 'simultaneous_stim', 'n_Trials');
StimParams = E.StimParams;
nEl      = E.simultaneous_stim;
n_Trials = E.n_Trials;
condAll   = cell2mat(E.TrialParams(2:end,2));
trialCond = condAll(1:nEl:end);
trialCond = trialCond(:);
assert(numel(trialCond) == n_Trials, 'exp file: %d trials in TrialParams but n_Trials = %d.', numel(trialCond), n_Trials);
if numel(trig) ~= n_Trials
    warning('FST:ntrig', '%d triggers but %d trials in the exp file. Using the first %d.', ...
        numel(trig), n_Trials, min(numel(trig), n_Trials));
end
nUse = min(numel(trig), n_Trials);

% conds file
if isempty(conds_file)
    f = listClean_('*.conds.mat');
    assert(~isempty(f), ['No *.conds.mat in %s. Copy the conds file saved by the stimulus PC into this ' ...
        'folder, or enter its full path in conds_file.'], data_folder);
    conds_file = fullfile(data_folder, f(1).name);
end
assert(isfile(conds_file), 'Cannot find %s.', conds_file);
K = load(conds_file);
assert(isfield(K,'cond') && isfield(K,'info'), 'conds file needs the variables cond and info.');
nCond = numel(K.cond);
fprintf('Conds file: %s (%d conditions)\n', conds_file, nCond);
assert(max(trialCond) <= nCond, 'exp file refers to condition %d but the conds file has only %d.', max(trialCond), nCond);

%% ================= PULSES OF EVERY CONDITION =================
% pulses{k} = [electrode, time_us, amp_uA], one row per pulse. Electrode i is row i of each trial's block in StimParams.
pulses = cell(nCond,1);
isElec = false(nCond,1);
for k = 1:nCond
    c = K.cond(k);
    amp = double(c.amp(:)');  ptd = double(c.ptd(:)');  nsp = double(c.nsp(:)');
    assert(numel(amp) == nEl, 'cond %d has %d electrodes but the exp file has %d rows per trial.', k, numel(amp), nEl);
    t0 = find(trialCond == k, 1);
    P = zeros(0,3);
    for e = 1:nEl
        if amp(e) > 0 && nsp(e) >= 1
            if isempty(t0)
                error('Condition %d has electrical pulses but no trial in the exp file uses it.', k);
            end
            row    = 1 + (t0-1)*nEl + e;                 % +1 for the header row
            period = double(StimParams{row,9});          % pulse train period, us
            if abs(double(StimParams{row,6}) - ptd(e)) > 1
                warning('FST:ptd', 'cond %d electrode %d: delay in conds file (%g us) differs from exp file (%g us).', ...
                    k, e, ptd(e), double(StimParams{row,6}));
            end
            t = ptd(e) + (0:nsp(e)-1)*period;
            P = [P; repmat(e,numel(t),1), t(:), repmat(amp(e),numel(t),1)]; %#ok<AGROW>
        end
    end
    pulses{k} = sortrows(P);
    isElec(k) = ~isempty(P);
end

%% ================= LEVEL STRUCTURE: SOURCE OF EVERY CONDITION =================
keyOf = @(P) mat2str(sortrows(P));
keyMap = containers.Map('KeyType','char','ValueType','double');
for k = find(isElec)'
    key = keyOf(pulses{k});
    if isKey(keyMap, key)
        warning('FST:dup', 'Conditions %d and %d deliver identical pulses; using %d as the source target.', keyMap(key), k, keyMap(key));
    else
        keyMap(key) = k;
    end
end

lastPulse_ms_cond = nan(nCond,1);
sourceCond        = nan(nCond,1);
for k = find(isElec)'
    P = pulses{k};
    tmax = max(P(:,2));
    lastPulse_ms_cond(k) = tmax/1000;
    rem = P(P(:,2) < tmax, :);                 % this condition minus its latest pulse(s)
    if isempty(rem), continue; end             % level 1: no source
    key = keyOf(rem);
    if isKey(keyMap, key)
        sourceCond(k) = keyMap(key);
    else
        warning('FST:nosrc', 'cond %d (%s): no condition equals it minus the last pulse. No source.', k, K.cond(k).label);
    end
end

%% ================= FIRST SPIKE TIME PER ELECTRICAL TRIAL =================
lastPulse_ms = nan(n_Trials,1);
isElecTrial  = false(n_Trials,1);
for tr = 1:n_Trials
    k = trialCond(tr);
    if isElec(k) && ~isnan(sourceCond(k))      % level 1 (single first pulse, no source) is the starting set: skipped
        lastPulse_ms(tr) = lastPulse_ms_cond(k);
        isElecTrial(tr)  = tr <= nUse && trig(tr) > 0;
    end
end

firstSpikeTimes = cell(nChn,1);
hasSpike = zeros(nChn, n_Trials);
t0_ms = nan(n_Trials,1);
t0_ms(1:nUse) = trig(1:nUse)' / FS * 1000;
idxTr = find(isElecTrial);

fprintf('\nExtracting first spike time after the last pulse (electrical trials only)...\n');
for ch = 1:nChn
    firstSpikeTimes{ch} = nan(n_Trials,1);
    S_ch = sp_clipped{ch};
    if isempty(S_ch), continue; end
    st = sort(S_ch(:,1));                                      % spike times, ms, column
    a  = t0_ms(idxTr) + lastPulse_ms(idxTr);                   % window start (absolute ms)
    b  = a + win_ms;                                           % window end
    ia = lowerBound_(st, a);                                   % first spike >= a
    ok = ia <= numel(st);
    first = nan(size(a));
    first(ok) = st(ia(ok));
    found = ok & (first <= b);
    rel = first - t0_ms(idxTr);
    v = nan(numel(idxTr),1);
    v(found) = rel(found);
    firstSpikeTimes{ch}(idxTr) = v;
    hasSpike(ch, idxTr(found)) = 1;
    fprintf('Ch %2d: first spike found in %d/%d electrical trials.\n', ch, sum(found), numel(idxTr));
end
fprintf('==============================================================\n');

%% ================= CONDITION TABLE =================
condTable = struct('cond',{},'kind',{},'level',{},'nLevels',{},'isFinal',{},'label',{}, ...
                   'sourceCond',{},'lastPulse_ms',{},'nTrials',{},'fractionFound',{});
fprintf('\n cond  kind     level   source  lastPulse(ms)  trials  found(%%)   label\n');
for k = 1:nCond
    inf_k = K.info(k);
    nTr_k = sum(trialCond == k);
    ff = NaN;
    if isElec(k) && ~isnan(sourceCond(k))
        tt = isElecTrial & trialCond == k;
        if any(tt), ff = 100*mean(mean(hasSpike(:, tt), 1)); end
    end
    condTable(k).cond = k;                      condTable(k).kind = inf_k.kind;
    condTable(k).level = inf_k.level;           condTable(k).nLevels = inf_k.nLevels;
    condTable(k).isFinal = inf_k.isFinal;       condTable(k).label = K.cond(k).label;
    condTable(k).sourceCond = sourceCond(k);    condTable(k).lastPulse_ms = lastPulse_ms_cond(k);
    condTable(k).nTrials = nTr_k;               condTable(k).fractionFound = ff;
    fprintf(' %3d   %-7s %2d/%-3d  %5s   %10s      %4d    %6s     %s\n', k, inf_k.kind, inf_k.level, inf_k.nLevels, ...
        numstr_(sourceCond(k)), numstr_(lastPulse_ms_cond(k)), nTr_k, numstr_(ff), K.cond(k).label);
end

%% ================= SAVE =================
out_name = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_FirstSpikeTimes.mat');
if strcmp(out_name, times_file)
    [p,n,e] = fileparts(times_file); out_name = fullfile(p, [n '_FirstSpikeTimes' e]);
end
save(out_name, 'firstSpikeTimes', 'hasSpike', 'lastPulse_ms', 'trialCond', 'condTable', 'win_ms', 'FS', 'n_Trials');
fprintf('\nSaved to: %s\n', out_name);

%% ================= HELPERS =================
function f = listClean_(pat)
% dir() without macOS "._" junk files
f = dir(pat);
if ~isempty(f), f = f(~startsWith({f.name}, '.')); end
end

function idx = lowerBound_(s, x)
% First index i with s(i) >= x (numel(s)+1 if none). s sorted ascending (column), x a column vector.
n  = numel(s);
lo = ones(size(x));
hi = (n+1)*ones(size(x));
while any(lo < hi)
    act = lo < hi;
    mid = floor((lo + hi)/2);
    midc = min(max(mid,1), n);
    right = act & (s(midc) < x);
    lo(right) = mid(right) + 1;
    left = act & ~right;
    hi(left) = mid(left);
end
idx = lo;
end

function s = numstr_(x)
if isnan(x), s = '-'; else, s = num2str(x, '%g'); end
end
