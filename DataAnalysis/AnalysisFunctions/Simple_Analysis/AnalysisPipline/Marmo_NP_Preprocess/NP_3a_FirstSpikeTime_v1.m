%% ============================================================
% NP_3a_FirstSpikeTime_v1.m   (OPTIONAL helper of step 3; only needed for NP_3_SpikeRecovery with boundaryMode = 'firstSpike')
% NP probe | marmoset pulse-TRAIN sessions | FIRST SPIKE AFTER THE LAST PULSE + LEVEL CHAIN | v1
%
% PURPOSE (same as NN_FirstSpikeTime_v1)
%   1. LEVEL CHAIN: the SOURCE of every electrical condition is the condition whose pulse set equals its own
%      minus its latest pulse(s) (for a simultaneous level: minus both pulses at the latest time).
%      Level 1 of every pattern has no source (it is never recovered).
%   2. For every electrical trial that has a source, and every channel: the first spike within win_ms after the
%      START of the trial's last pulse (= where real recording resumes after the last artifact).
%   Read-only for the spike data. Blank / visual trials and level-1 trials are left NaN.
%
% INPUTS (recording1 folder): <short>.NP_TrialTable.mat, <short>.NP_Spikes.mat (or spikeFile)
% OUTPUT : <short>.NP_FirstSpikeTimes.mat   (NN names)
%   firstSpikeTimes{ch}(trial)  ms after the trigger (NaN = none found / not an electrical trial with a source)
%   hasSpike(ch, trial)         1 if a spike was found
%   lastPulse_ms(trial)         start of the last pulse of the trial, ms after the trigger (NaN = not used)
%   trialCond(trial)            condition of every trial (TT.trial.cond)
%   condTable                   per condition: label, kind, level, sourceCond, lastPulse_ms, isElec
%   win_ms, FS, n_Trials, spikeFile, ttFile, version, date
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile  = '';      % '' = <short>.NP_Spikes.mat (the un-recovered spikes of step 2)
win_ms     = 10;      % how far after the START of the last pulse to look for the first spike (ms)
overwrite  = false;   % true = replace an existing file (the old one is first copied to a dated backup)
%% =================================================

codeVersion = 'NP_3a_FirstSpikeTime_v1 (2026-11-12)';
[TT, ttFile, recDir, short] = loadTT_(sessionDir);
if isempty(spikeFile), spikeFile = [short '.NP_Spikes.mat']; end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
assert(isfile(spikeFile), 'Spike file not found: %s', spikeFile);
SP = load(spikeFile, 'sp_clipped', 'sp_trial', 'FS', 'nTrialsProcessed');
FS = SP.FS;  nChn = numel(SP.sp_clipped);  n_Trials = SP.nTrialsProcessed;

%% ---- chain ------------------------------------------------------------------------------------------
condTable = buildChain_(TT);
trialCond = TT.trial.cond(1:n_Trials);
lastPulse_ms = nan(n_Trials, 1);
useTr = false(n_Trials, 1);
for k = 1:n_Trials
    c = trialCond(k);
    if condTable(c).isElec && ~isnan(condTable(c).sourceCond)
        useTr(k) = true;  lastPulse_ms(k) = TT.trial.lastPulseMs(k);      % design + measured delay (from step 1)
    end
end

%% ---- first spike per channel and trial ---------------------------------------------------------------
trigMs = (double(TT.trial.trigAP(1:n_Trials)) - 1)/FS*1e3;
firstSpikeTimes = cell(nChn, 1);  hasSpike = zeros(nChn, n_Trials);
idx = find(useTr);
for ch = 1:nChn
    firstSpikeTimes{ch} = nan(n_Trials, 1);
    st = sort(SP.sp_clipped{ch}(:));
    if isempty(st), continue; end
    a = trigMs(idx) + lastPulse_ms(idx);  b = a + win_ms;     % window starts, in time order
    j = 1;  ns = numel(st);
    for q = 1:numel(idx)                                        % one pass through the sorted spikes
        while j <= ns && st(j) < a(q), j = j + 1; end
        if j > ns, break; end
        if st(j) <= b(q)
            firstSpikeTimes{ch}(idx(q)) = st(j) - trigMs(idx(q));  hasSpike(ch, idx(q)) = 1;
        end
    end
end

%% ---- printout --------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP FIRST SPIKE TIME (%s)\n', codeVersion);
fprintf('============================================================\n');
fprintf('Session : %s\nSpikes  : %s\nWindow  : %g ms after the start of the last pulse\n', short, spikeFile, win_ms);
printChain_(condTable);
fprintf('\n  cond  label                         trials   spike found (%% of channel x trial)   median first spike - last pulse (ms)\n');
for c = find([condTable.isElec] & ~isnan([condTable.sourceCond]))
    kk = find(trialCond == c);
    if isempty(kk), continue; end
    H = hasSpike(:, kk);
    d = cellfun(@(v) median(v(kk) - lastPulse_ms(kk), 'omitnan'), firstSpikeTimes);
    fprintf('  %4d  %-28s  %5d        %5.1f %%                              %5.2f\n', c, condTable(c).label, numel(kk), 100*mean(H(:)), median(d, 'omitnan'));
end

%% ---- save ---------------------------------------------------------------------------------------------------
outFile = fullfile(recDir, [short '.NP_FirstSpikeTimes.mat']);
if isfile(outFile)
    if ~overwrite, error('%s already exists and overwrite = false. Nothing was saved.', outFile); end
    bak = strrep(outFile, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);  copyfile(outFile, bak);
    fprintf('Old file copied to %s\n', bak);
end
version = codeVersion;  date = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<NASGU>
save(outFile, 'firstSpikeTimes', 'hasSpike', 'lastPulse_ms', 'trialCond', 'condTable', 'win_ms', 'FS', 'n_Trials', ...
     'spikeFile', 'ttFile', 'version', 'date');
fprintf('\nSaved %s\n', outFile);

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function [TT, ttFile, recDir, short] = loadTT_(sessionDir)
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
ttFile = fullfile(recDir, tf(1).name);
L = load(ttFile, 'TT');  TT = L.TT;
end

function condTable = buildChain_(TT)
% source of every condition = the condition whose pulse set equals its own minus its latest pulse(s)
nC = numel(TT.cond);
key = cell(nC, 1);  P = cell(nC, 1);  isElec = false(nC, 1);  lastP = nan(nC, 1);
for c = 1:nC
    A = round(TT.cond(c).designMsA(:)*1000);  B = round(TT.cond(c).designMsB(:)*1000);   % us
    Pc = sortrows([[ones(numel(A),1); 2*ones(numel(B),1)], [A; B]]);
    P{c} = Pc;  isElec(c) = ~isempty(Pc);
    if isElec(c), lastP(c) = max(Pc(:,2))/1000; key{c} = mat2str(Pc); end
end
km = containers.Map('KeyType', 'char', 'ValueType', 'double');
for c = find(isElec)'
    if ~isKey(km, key{c}), km(key{c}) = c; end
end
src = nan(nC, 1);
for c = find(isElec)'
    Pc = P{c};  rem = Pc(Pc(:,2) < max(Pc(:,2)), :);
    if isempty(rem), continue; end
    kr = mat2str(sortrows(rem));
    if isKey(km, kr), src(c) = km(kr);
    else, warning('NPCH:nosrc', 'cond %d (%s): no condition equals it minus its last pulse -> no source.', c, TT.cond(c).label);
    end
end
condTable = struct('cond', num2cell((1:nC)'), 'label', {TT.cond.label}', 'kind', {TT.cond.kind}', ...
    'level', {TT.cond.level}', 'sourceCond', num2cell(src), 'lastPulse_ms', num2cell(lastP), 'isElec', num2cell(isElec));
end

function printChain_(condTable)
fprintf('\nLevel chain (source = same pulses minus the latest one(s)):\n  cond  kind   level  source  last pulse (ms, design)  label\n');
for c = 1:numel(condTable)
    C = condTable(c);
    if ~C.isElec, continue; end
    fprintf('  %4d  %-5s  %5d  %6s  %8.2f                 %s\n', c, C.kind, C.level, num2strNaN_(C.sourceCond), C.lastPulse_ms, C.label);
end
end

function s = num2strNaN_(v)
if isnan(v), s = '-'; else, s = sprintf('%d', v); end
end
