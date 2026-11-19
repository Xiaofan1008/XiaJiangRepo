%% ============================================================
% NN_FirstSpikeTime_v2.m
% Step 1 of the PULSE-TRAIN artifact recovery pipeline (marmoset NN data).
% v2 = v1 WITHOUT the *.conds.mat file: everything comes from the EXPERIMENT FILE (*_exp_datafile_*.mat)
% of the new experiment code (it saves StimMeta, one entry per condition). v1 is not changed.
%
% For every ELECTRICAL trial, find the first spike that appears after the
% start of that trial's LAST pulse, i.e. where real recording resumes after
% the last artifact. Step 2 (the recovery script) uses it as the end of the
% window that is blanked and filled from the level with one pulse fewer.
%
% It also works out the LEVEL STRUCTURE of the design. Each condition is a set of pulses
% (electrode, time, amplitude), read from StimParams (the first trial of the condition). Its SOURCE
% condition is the one whose pulse set equals this one minus its latest pulse (for a simultaneous
% level, minus both pulses at the latest time). Level 1 of every chain has no source.
%
% WHAT IS NEW IN v2
%   - kind / level / nLevels / isFinal of every condition come from exp file.StimMeta
%     (Kind, TrainLevel, TotalTrainLevels, IsFinalLevel); the label is built from them.
%   - The SOURCE is searched only INSIDE THE SAME STIMULUS SET (StimMeta.StimSetIndex = one chain:
%     sim, AB, BA, singleA, singleB). v1 matched on the pulses alone, which would link the level 2 of
%     the 'singleA' chain to level 1 of the AB chain (identical pulse set [A at 0]) and 'singleB' to BA.
%     Every chain, including the single-electrode linearity chains, is self-contained.
%   - kinds 'singleA' / 'singleB' (linearity chains) are handled like any other chain.
%   - checks: the pulses read from StimParams must agree with StimMeta (number of pulses, delay); a
%     warning is printed if the spike file is not inside data_folder (wrong dataset).
%   If the exp file has no StimMeta (older experiment code), use NN_FirstSpikeTime_v1 with the conds file.
%
% READ-ONLY with respect to the spike data: sp_clipped is only inspected.
% Blank and visual conditions are skipped (left NaN), and so are level-1 (single first
% pulse) conditions: they are the starting set that recovery copies from, never recovered.
%
% Inputs  (in data_folder unless entered):
%   times_file   FILTERED spike times (sp_clipped), entered directly
%   *_exp_datafile_*.mat   TrialParams / StimParams / StimMeta of the experiment
%   *.trig.dat   trigger times
% Output (next to times_file):
%   <base>.sp_xia_FirstSpikeTimes.mat
%       firstSpikeTimes{ch}(trial)  ms after the trigger (NaN = none found)
%       hasSpike(ch,trial)          1 if a spike was found
%       lastPulse_ms(trial)         start of the last pulse of that trial, ms
%       trialCond(trial)            condition number of every trial
%       condTable                   per condition: kind, level, source, last pulse ... (+ stimSet)
%       win_ms, FS, n_Trials
% Same variables as v1, so the recovery scripts (step 2) work unchanged. Written without access to
% MATLAB: please report any error.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
analysisPath = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/';
                    % added to the MATLAB path (needed for cleanTrig_sabquick only)

data_folder = '/Volumes/MACData/Data/Data_Xia/CJ288/Xia_15Pulses6_100Hz_30uA';
                    % folder holding this dataset's trigger (*.trig.dat) and
                    % stim-protocol (*_exp_datafile_*.mat) files.

times_file = '/Volumes/MACData/Data/Data_Xia/CJ288/Xia_15Pulses6_100Hz_30uA/Xia_15Pulses6_100Hz_30uA.sp_xia_filtered.mat';
                    % the FILTERED spike-times file (sp_clipped) of this same dataset.


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
E = load(f(1).name, 'TrialParams', 'StimParams', 'StimMeta', 'simultaneous_stim', 'n_Trials');
expFile = fullfile(f(1).folder, f(1).name);
StimParams = E.StimParams;
StimMeta   = E.StimMeta;
if ~startsWith(char(times_file), char(data_folder))
    warning('FST:folder', ['times_file is NOT inside data_folder:\n  times_file : %s\n  data_folder: %s\n' ...
        'The trigger and exp file are read from data_folder: make sure both belong to the SAME dataset.'], times_file, data_folder);
end
nEl      = double(E.simultaneous_stim);
n_Trials = double(E.n_Trials);
condAll   = cell2mat(E.TrialParams(2:end,2));
trialCond = condAll(1:nEl:end);
trialCond = trialCond(:);
assert(numel(trialCond) == n_Trials, 'exp file: %d trials in TrialParams but n_Trials = %d.', numel(trialCond), n_Trials);
if numel(trig) ~= n_Trials
    warning('FST:ntrig', '%d triggers but %d trials in the exp file. Using the first %d.', ...
        numel(trig), n_Trials, min(numel(trig), n_Trials));
end
nUse = min(numel(trig), n_Trials);

% conditions from the experiment file (replaces the conds file of v1)
nCond = numel(StimMeta);
fprintf('Exp file: %s (%d conditions in StimMeta)\n', expFile, nCond);
assert(max(trialCond) <= nCond, 'exp file refers to condition %d but StimMeta has only %d.', max(trialCond), nCond);
mKind = cell(nCond,1);  mLevel = nan(nCond,1);  mNLev = nan(nCond,1);  mFinal = false(nCond,1);
mSet = nan(nCond,1);  mNsp = nan(nCond,nEl);  mPtd = nan(nCond,nEl);
for k = 1:nCond
    sm = StimMeta(k);
    assert(double(sm.ConditionID) == k, 'StimMeta(%d).ConditionID is %g (expected %d).', k, double(sm.ConditionID), k);
    mKind{k}  = char(sm.Kind);
    mLevel(k) = double(sm.TrainLevel(1));   mNLev(k) = double(sm.TotalTrainLevels(1));
    mFinal(k) = double(sm.IsFinalLevel(1)) > 0;
    mSet(k)   = double(sm.StimSetIndex(1));
    mNsp(k,:) = double(sm.NPulsesPerElectrode(:)');   mPtd(k,:) = double(sm.PTD_us(:)');
end

%% ================= PULSES OF EVERY CONDITION =================
% pulses{k} = [electrode, time_us, amp_uA], one row per pulse, read from StimParams (first trial of the condition).
% Electrode i is row i of each trial's block in StimParams (+1 header row). An electrode takes part if its amplitude > 0
% (amplitude -1 = switched off).
pulses = cell(nCond,1);
isElec = false(nCond,1);
for k = 1:nCond
    t0 = find(trialCond == k, 1);
    P = zeros(0,3);
    if isempty(t0)
        warning('FST:notrial', 'Condition %d has no trial in the exp file.', k);
    else
        for e = 1:nEl
            row    = 1 + (t0-1)*nEl + e;
            ampE   = double(StimParams{row,16});         % first phase amplitude, uA
            nE     = double(StimParams{row,8});          % number of pulses in the train
            delay  = double(StimParams{row,6});          % post-trigger delay, us
            period = double(StimParams{row,9});          % pulse train period, us
            if ampE > 0 && nE >= 1
                t = delay + (0:nE-1)*period;
                P = [P; repmat(e,numel(t),1), t(:), repmat(ampE,numel(t),1)]; %#ok<AGROW>
                % cross-check with StimMeta
                if abs(mNsp(k,e) - nE) > 0 || abs(mPtd(k,e) - delay) > 1
                    warning('FST:meta', 'cond %d electrode %d: StimParams has %g pulses, delay %g us; StimMeta has %g pulses, delay %g us.', ...
                        k, e, nE, delay, mNsp(k,e), mPtd(k,e));
                end
            elseif mNsp(k,e) >= 1 && ~strcmp(mKind{k}, 'blank') && ~strcmp(mKind{k}, 'visual') && ampE <= 0 && double(StimMeta(k).Amp_uA(e)) > 0
                warning('FST:meta', 'cond %d electrode %d: StimMeta has an amplitude but StimParams has it switched off.', k, e);
            end
        end
    end
    pulses{k} = sortrows(P);
    isElec(k) = ~isempty(P);
end

%% ================= LEVEL STRUCTURE: SOURCE OF EVERY CONDITION =================
% the source is searched only inside the SAME stimulus set (= the same chain), see header.
keyOf = @(k, P) sprintf('S%g|%s', mSet(k), mat2str(sortrows(P)));
keyMap = containers.Map('KeyType','char','ValueType','double');
for k = find(isElec)'
    key = keyOf(k, pulses{k});
    if isKey(keyMap, key)
        warning('FST:dup', 'Conditions %d and %d deliver identical pulses in the same stimulus set; using %d as the source target.', keyMap(key), k, keyMap(key));
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
    key = keyOf(k, rem);
    if isKey(keyMap, key)
        sourceCond(k) = keyMap(key);
    else
        warning('FST:nosrc', 'cond %d (%s level %g): no condition of the same set equals it minus the last pulse. No source.', k, mKind{k}, mLevel(k));
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
                   'sourceCond',{},'lastPulse_ms',{},'nTrials',{},'fractionFound',{},'stimSet',{});
fprintf('\n cond  set  kind     level   source  lastPulse(ms)  trials  found(%%)   label\n');
for k = 1:nCond
    nTr_k = sum(trialCond == k);
    ff = NaN;
    if isElec(k) && ~isnan(sourceCond(k))
        tt = isElecTrial & trialCond == k;
        if any(tt), ff = 100*mean(mean(hasSpike(:, tt), 1)); end
    end
    if isElec(k)
        lab = sprintf('%s %g/%g [%s]', mKind{k}, mLevel(k), mNLev(k), strtrim(sprintf('%g ', mNsp(k,:))));
    else
        lab = mKind{k};
    end
    condTable(k).cond = k;                      condTable(k).kind = mKind{k};
    condTable(k).level = mLevel(k);             condTable(k).nLevels = mNLev(k);
    condTable(k).isFinal = mFinal(k);           condTable(k).label = lab;
    condTable(k).sourceCond = sourceCond(k);    condTable(k).lastPulse_ms = lastPulse_ms_cond(k);
    condTable(k).nTrials = nTr_k;               condTable(k).fractionFound = ff;
    condTable(k).stimSet = mSet(k);
    fprintf(' %3d   %3s  %-8s %3s/%-3s  %5s   %10s      %4d    %6s     %s\n', k, numstr_(mSet(k)), mKind{k}, numstr_(mLevel(k)), numstr_(mNLev(k)), ...
        numstr_(sourceCond(k)), numstr_(lastPulse_ms_cond(k)), nTr_k, numstr_(ff), lab);
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
