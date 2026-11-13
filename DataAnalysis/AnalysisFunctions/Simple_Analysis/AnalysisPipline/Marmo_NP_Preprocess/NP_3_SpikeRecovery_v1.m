%% ============================================================
% NP_3_SpikeRecovery_v1.m   (step 3 of the NP pipeline)
% NP probe | marmoset pulse-TRAIN sessions | CHAINED TRAIN RECOVERY (+ optional RandomEvent simulation) | v1
%
% PURPOSE (same idea as NN_SpikeRecovery_v1)
%   The artifact blank hides the train. Levels are processed in ascending order; in every target trial of level L
%   the replaced stretch is emptied and filled with the spikes of the SAME relative stretch of its paired trial of
%   the source level (already recovered), so level 3 uses recovered level 2, which used recovered level 1 ...
%   Level 1 of every pattern, blank and visual trials, and every spike outside the replaced stretches are untouched.
%   Source of a condition = the condition whose pulses equal its own minus its latest pulse(s) (as NN_FirstSpikeTime).
%   Pairing: i-th target trial <-> i-th source trial (cycling if the source has fewer trials), same for all channels.
%
% boundaryMode  (what is replaced, per channel and trial)
%   'blank'      (DEFAULT) exactly the blanked stretch(es) saved by NP_2_SpikeDetect (blank_ms).
%   'firstSpike' the NN rule: [0, win_end) after the trigger, win_end = first spike after the last pulse
%                (NP_3a_FirstSpikeTime) capped at lastPulse + cap_ms; no spike found: lastPulse + fallback_ms.
%
% simMode
%   'none'        (DEFAULT) normal recovery.
%   'randomEvent' SIMULATION (as NN_SpikeRecovery_RandomEvent_v1): after a level is recovered, ALL spikes in
%                 [t_L, t_L + event_ms) after its last pulse t_L (copied AND genuine ones) of all its trials are
%                 pooled per channel and exactly round(percent/100 x N) of them are kept at random (fixed seed).
%                 Higher levels copy that stretch, so every pulse carries the percentage of ITS OWN level.
%                 percent = 100 everywhere gives exactly the same result as 'none'. Use it only to test whether the
%                 analysis finds a KNOWN reduction; it deliberately removes real spikes.
%
% SETTINGS PER LEVEL (levels start at 2; level = the pattern's own level number, TT.cond.level)
%   sim_default / seq_default   = [fallback_ms  cap_ms  percent  event_ms]
%   sim_override / seq_override = rows [level fallback_ms cap_ms percent event_ms] for levels that differ
%   fallback_ms, cap_ms : only used with boundaryMode = 'firstSpike'
%   percent, event_ms   : only used with simMode = 'randomEvent';  event_ms = NaN -> the pulse spacing of that level
%                         step taken from the design (e.g. 10 ms sim / 5 ms seq at 100 Hz; 5 / 2.5 ms at 200 Hz)
%
% INPUTS (recording1 folder): <short>.NP_TrialTable.mat, <short>.NP_Spikes.mat,
%   <short>.NP_FirstSpikeTimes.mat (only 'firstSpike'), <short>.NP_Spikes_waveforms.mat (only saveWaveforms = true)
%
% OUTPUT  <short>.NP_Spikes_TrainRecovery[_<suffix>].mat   (same layout as NP_Spikes.mat, readable by NP_SP_RasterPSTH)
%   suffix: out_suffix, or automatically 'FirstSpike' / 'RandomEvent' / 'FirstSpike_RandomEvent' ('' for blank + none)
%   sp_clipped{ch}, sp_trial{ch}, nSpk(ch, trial)   recovered spikes (ms from the AP file start), as in NP_Spikes.mat
%   blank_ms{ch, trial}         stretches that are STILL without data after recovery (level 1 keeps its blank)
%   blank_ms_original           the blanks of NP_Spikes.mat
%   rec_ms{ch, trial}           replaced stretches [start end] (ms after the trigger); spikes inside are copied ones
%   srcTrial(trial), isRecovered(trial)   source trial (row of TT.trial) and whether the trial was recovered
%   condTable, depth            the chain (source of every condition) and the chain depth (0 = not recovered)
%   recLog(cond)                per condition, summed over channels: targets, deleted, added, (firstSpike) found /
%                               fallback, (randomEvent) spikes in the window and kept
%   Recovery_params             all settings and the per-condition values used
%   + goodCh, noise_uV, thr_uV, FS, preMs, postMs, nTrialsProcessed, spikeFile, ttFile, version, date
%   <short>.NP_Spikes_TrainRecovery[_<suffix>]_waveforms.mat   ONLY if saveWaveforms = true (separate file)
%
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir   = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile    = '';        % '' = <short>.NP_Spikes.mat
boundaryMode = 'blank';   % 'blank' (default) | 'firstSpike' (needs NP_3a_FirstSpikeTime)
simMode      = 'randomEvent';    % 'none' (default) | 'randomEvent' (simulation)
out_suffix   = 'randomEvent';        % '' = automatic (see header); or your own text, e.g. 'RandomEvent_test2'
saveWaveforms = false;    % true = also recover the waveforms (needs <short>.NP_Spikes_waveforms.mat), separate file
overwrite    = true;     % true = replace existing output (dated backups of the old files are made first)

% ---- settings per level: [fallback_ms  cap_ms  percent  event_ms]   (event_ms NaN = design pulse spacing) ----
sim_default  = [10 10 100 NaN];   % simultaneous pattern, every level not listed in sim_override
sim_override = [
    2,10,10,60,10;
    9,10,10,60,10;
    10,10,10,50,10;
    11,10,10,40,10;
    12,10,10,30,10;
    13,10,10,30,10;
    14,10,10,30,10;
    15,10,10,30,20;
    ];                % rows [level fallback_ms cap_ms percent event_ms], e.g. [10 10 10 50 NaN; 15 10 10 30 20]
seq_default  = [10 10 100 NaN];   % sequential pattern (AB or BA), every level not listed in seq_override
seq_override = [
    25,10,10,40,5;
    26,10,10,40,5;
    27,10,10,15,5;
    28,10,10,15,5;
    29,10,10,5,5;
    30,10,10,5,20;    
];                % rows [level fallback_ms cap_ms percent event_ms]
rand_seed    = 1;                 % seed of the random choice of spikes (randomEvent, percent < 100)
%% =================================================

codeVersion = 'NP_3_SpikeRecovery_v1 (2026-11-13, stretch includes its end sample)';
tStart = tic;
assert(any(strcmp(boundaryMode, {'blank', 'firstSpike'})), 'boundaryMode must be ''blank'' or ''firstSpike''.');
assert(any(strcmp(simMode, {'none', 'randomEvent'})), 'simMode must be ''none'' or ''randomEvent''.');

%% ---- load ---------------------------------------------------------------------------------------------
[TT, ttFile, recDir, short] = loadTT_(sessionDir);
if isempty(spikeFile), spikeFile = [short '.NP_Spikes.mat']; end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
assert(isfile(spikeFile), 'Spike file not found: %s', spikeFile);
SP = load(spikeFile);
FS = SP.FS;  nChn = numel(SP.sp_clipped);  nTr = SP.nTrialsProcessed;
trialCond = TT.trial.cond(1:nTr);
trigMs = (double(TT.trial.trigAP(1:nTr)) - 1)/FS*1e3;
lastP  = TT.trial.lastPulseMs(1:nTr);                     % start of the last pulse, ms after the trigger

condTable = buildChain_(TT);
nCond = numel(condTable);
fstFile = '';
if strcmp(boundaryMode, 'firstSpike')
    fstFile = fullfile(recDir, [short '.NP_FirstSpikeTimes.mat']);
    assert(isfile(fstFile), 'Missing %s: run NP_3a_FirstSpikeTime first (or use boundaryMode = ''blank'').', fstFile);
    F = load(fstFile, 'firstSpikeTimes', 'condTable', 'n_Trials');
    assert(numel(F.firstSpikeTimes) == nChn && F.n_Trials == nTr, 'The first-spike file does not match the spike file (channels / trials).');
    s1 = [condTable.sourceCond];  s2 = [F.condTable.sourceCond];
    assert(isequal(isnan(s1), isnan(s2)) && all(s1(~isnan(s1)) == s2(~isnan(s2))), ...
           'The level chain of the first-spike file differs from the one built here: re-run NP_3a_FirstSpikeTime.');
end
if saveWaveforms
    wfIn = fullfile(recDir, [short '.NP_Spikes_waveforms.mat']);
    assert(isfile(wfIn), 'saveWaveforms = true but %s is missing (run NP_2 with saveWaveforms = true).', wfIn);
    fprintf('loading waveforms ...\n');
    WF = load(wfIn, 'sp_waveforms', 'wfScale_uV', 'wfTime_ms');
end

%% ---- chain depth -------------------------------------------------------------------------------------------
srcOf = [condTable.sourceCond]';  depth = zeros(nCond, 1);
for it = 1:nCond + 1
    changed = false;
    for k = 1:nCond
        if ~isnan(srcOf(k))
            d = depth(srcOf(k)) + 1;
            if d ~= depth(k), depth(k) = d; changed = true; end
        end
    end
    if ~changed, break; end
    assert(it <= nCond, 'Source chain has a loop: check the design.');
end
maxDepth = max(depth);

%% ---- settings per condition ------------------------------------------------------------------------------------
chkSet_('sim', sim_default, sim_override);  chkSet_('seq', seq_default, seq_override);
kindK = {condTable.kind};  lvlK = [condTable.level]';  lpK = [condTable.lastPulse_ms]';
fbK = nan(nCond,1);  capK = nan(nCond,1);  pctK = nan(nCond,1);  evK = nan(nCond,1);  spK = nan(nCond,1);
for k = find(depth > 0)'
    if strcmp(kindK{k}, 'sim'), def = sim_default; ov = sim_override; else, def = seq_default; ov = seq_override; end
    s = def;
    if ~isempty(ov), r = find(ov(:,1) == lvlK(k), 1); if ~isempty(r), s = ov(r, 2:5); end, end
    spK(k) = lpK(k) - lpK(srcOf(k));                       % pulse spacing of this level step (design, ms)
    fbK(k) = s(1);  capK(k) = s(2);  pctK(k) = s(3);  evK(k) = s(4);
    if isnan(evK(k)), evK(k) = spK(k); end
end
if strcmp(simMode, 'none'), pctK(depth > 0) = 100; end
for pp = {'sim', 'seq'}
    if strcmp(pp{1}, 'sim'), ov = sim_override; isP = strcmp(kindK, 'sim'); else, ov = seq_override; isP = ~strcmp(kindK, 'sim'); end
    mx = max([0; lvlK(isP(:) & depth > 0)]);
    for r = 1:size(ov, 1)
        if ov(r,1) > mx, warning('NPSR:level', '%s_override: level %g does not exist (max level %d).', pp{1}, ov(r,1), mx); end
    end
end

%% ---- output names ------------------------------------------------------------------------------------------
suffix = out_suffix;
if isempty(suffix)
    parts = {};
    if strcmp(boundaryMode, 'firstSpike'), parts{end+1} = 'FirstSpike'; end
    if strcmp(simMode, 'randomEvent'), parts{end+1} = 'RandomEvent'; end
    suffix = strjoin(parts, '_');
end
base = [short '.NP_Spikes_TrainRecovery'];  if ~isempty(suffix), base = [base '_' suffix]; end
outFile = fullfile(recDir, [base '.mat']);  wfOut = fullfile(recDir, [base '_waveforms.mat']);
outs = {outFile};  if saveWaveforms, outs{end+1} = wfOut; end
for q = 1:numel(outs)
    if isfile(outs{q}) && ~overwrite
        error('%s already exists and overwrite = false. Nothing was done.', outs{q});
    end
end

%% ---- print settings ---------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP SPIKE RECOVERY (%s)\n', codeVersion);
fprintf('============================================================\n');
fprintf('Session  : %s\nSpikes   : %s\nBoundary : %s | simulation: %s | waveforms: %d\n', short, spikeFile, boundaryMode, simMode, saveWaveforms);
fprintf('%d conditions, chain depth 0..%d (0 = not recovered)\n', nCond, maxDepth);
fprintf('\n  cond  kind   level  source  depth  spacing(ms)');
if strcmp(boundaryMode, 'firstSpike'), fprintf('  fallback(ms)  cap(ms)'); end
if strcmp(simMode, 'randomEvent'), fprintf('  percent  event(ms)'); end
fprintf('  label\n');
for k = find(depth > 0)'
    fprintf('  %4d  %-5s  %5d  %6d  %5d  %10.2f', k, kindK{k}, lvlK(k), srcOf(k), depth(k), spK(k));
    if strcmp(boundaryMode, 'firstSpike'), fprintf('  %11.2f  %7.2f', fbK(k), capK(k)); end
    if strcmp(simMode, 'randomEvent'), fprintf('  %7.1f  %9.2f', pctK(k), evK(k)); end
    fprintf('  %s\n', condTable(k).label);
    if pctK(k) < 100 && evK(k) > spK(k) + 1e-9
        warning('NPSR:event', 'cond %d (level %d): event_ms %g is longer than the pulse spacing %g ms; the part after the next pulse is overwritten by the next level.', k, lvlK(k), evK(k), spK(k));
    end
end
if strcmp(simMode, 'randomEvent'), fprintf('random seed %d\n', rand_seed); end

%% ---- pairs (same for all channels) ------------------------------------------------------------------------------
pairs = cell(nCond, 1);  srcTrial = nan(nTr, 1);  isRecovered = false(nTr, 1);
for k = find(depth > 0)'
    tT = find(trialCond == k);  tS = find(trialCond == srcOf(k));
    if isempty(tT), continue; end
    if isempty(tS), warning('NPSR:nosrc', 'cond %d: its source cond %d has no trial. Not recovered.', k, srcOf(k)); continue; end
    sSel = tS(mod((1:numel(tT)) - 1, numel(tS)) + 1);
    pairs{k} = [tT(:), sSel(:)];  srcTrial(tT) = sSel;  isRecovered(tT) = true;
end

%% ---- recovery, channel by channel -----------------------------------------------------------------------------------
recLog = struct('cond', num2cell((1:nCond)'), 'depth', num2cell(depth), 'source', num2cell(srcOf), ...
    'nTargets', 0, 'nDeleted', 0, 'nAdded', 0, 'nFound', 0, 'nFallback', 0, 'nWindow', 0, 'nKept', 0);
for k = 1:nCond, if ~isempty(pairs{k}), recLog(k).nTargets = size(pairs{k}, 1); end, end
rng(rand_seed);
hS = 0.5/FS*1e3;                                          % half a sample (ms): the stretch includes its first and last sample,
                                                          % as the blank of NP_2 does (also absorbs rounding of the times)
sp_clipped = cell(1, nChn);  sp_trial = cell(1, nChn);  sp_waveforms = cell(1, nChn);
rec_ms = cell(nChn, nTr);  blank_ms = cell(nChn, nTr);
for ch = 1:nChn
    % spikes (and waveforms) of this channel, ms after each trial's trigger, grouped by trial
    kk = SP.sp_trial{ch}(:);  rel = SP.sp_clipped{ch}(:) - trigMs(kk);
    if saveWaveforms, Wch = WF.sp_waveforms{ch}; else, Wch = []; end
    [relBy, wfBy] = groupByTrial_(rel, kk, nTr, saveWaveforms, Wch);
    for d = 1:maxDepth
        for k = find(depth == d)'
            P = pairs{k};
            if isempty(P), continue; end
            for q = 1:size(P, 1)
                tT = P(q,1);  tS = P(q,2);
                if strcmp(boundaryMode, 'blank')
                    R = SP.blank_ms{ch, tT};
                else
                    f = F.firstSpikeTimes{ch}(tT);  lp = lastP(tT);
                    if isfinite(f) && f > 0, we = min(f, lp + capK(k)); recLog(k).nFound = recLog(k).nFound + 1;
                    else,                    we = lp + fbK(k);          recLog(k).nFallback = recLog(k).nFallback + 1;
                    end
                    R = [0 we];
                end
                if isempty(R), continue; end
                inT = inIv_(relBy{tT}, R, hS);  inS = inIv_(relBy{tS}, R, hS);   % start and end samples included
                recLog(k).nDeleted = recLog(k).nDeleted + sum(inT);
                recLog(k).nAdded   = recLog(k).nAdded + sum(inS);
                [relBy{tT}, o] = sort([relBy{tT}(~inT); relBy{tS}(inS)]);
                if saveWaveforms, w = [wfBy{tT}(~inT, :); wfBy{tS}(inS, :)]; wfBy{tT} = w(o, :); end
                rec_ms{ch, tT} = R;
            end
            % ---- RandomEvent thinning of this level ----
            if pctK(k) < 100
                cand = cell(size(P, 1), 1);
                for q = 1:size(P, 1)
                    tT = P(q,1);  r_ = relBy{tT};
                    cand{q} = find(r_ >= lastP(tT) & r_ < lastP(tT) + evK(k));
                end
                nPer = cellfun(@numel, cand);  N = sum(nPer);
                if N > 0
                    keepN = round(pctK(k)/100*N);
                    keepAll = false(N, 1);  keepAll(randperm(N, keepN)) = true;
                    pos = [0; cumsum(nPer(:))];
                    for q = 1:size(P, 1)
                        if nPer(q) == 0, continue; end
                        tT = P(q,1);  drop = cand{q}(~keepAll(pos(q)+1:pos(q+1)));
                        relBy{tT}(drop) = [];
                        if saveWaveforms, wfBy{tT}(drop, :) = []; end
                    end
                    recLog(k).nWindow = recLog(k).nWindow + N;  recLog(k).nKept = recLog(k).nKept + keepN;
                end
            end
        end
    end
    % remaining blank = original blank minus the replaced stretch
    for t = 1:nTr
        blank_ms{ch, t} = subIv_(SP.blank_ms{ch, t}, rec_ms{ch, t});
    end
    % back to the NN layout
    a = cell(nTr, 1);  b = cell(nTr, 1);
    for t = 1:nTr, a{t} = trigMs(t) + relBy{t}; b{t} = t*ones(numel(relBy{t}), 1); end
    sp_clipped{ch} = vertcat(a{:});  sp_trial{ch} = vertcat(b{:});
    if saveWaveforms, sp_waveforms{ch} = vertcat(wfBy{:}); end
    if mod(ch, 64) == 0 || ch == nChn, fprintf('  channel %d/%d (%.0f s)\n', ch, nChn, toc(tStart)); end
end
nSpk = zeros(nChn, nTr);
for ch = 1:nChn, if ~isempty(sp_trial{ch}), nSpk(ch, :) = accumarray(sp_trial{ch}, 1, [nTr 1])'; end, end

%% ---- printout ------------------------------------------------------------------------------------------------------
fprintf('\n  cond  depth  source  targets   deleted     added');
if strcmp(boundaryMode, 'firstSpike'), fprintf('   found  fallback'); end
if strcmp(simMode, 'randomEvent'), fprintf('   window    kept  (kept %%)'); end
fprintf('   (summed over %d channels)\n', nChn);
for k = find(depth > 0)'
    L_ = recLog(k);
    fprintf('  %4d  %5d  %6d  %7d  %8d  %8d', k, L_.depth, L_.source, L_.nTargets, L_.nDeleted, L_.nAdded);
    if strcmp(boundaryMode, 'firstSpike'), fprintf('  %6d  %8d', L_.nFound, L_.nFallback); end
    if strcmp(simMode, 'randomEvent'), fprintf('  %7d  %6d  (%5.1f)', L_.nWindow, L_.nKept, 100*L_.nKept/max(L_.nWindow, 1)); end
    fprintf('\n');
end
fprintf('spikes: %d before, %d after recovery\n', sum(SP.nSpk(:)), sum(nSpk(:)));

%% ---- save ------------------------------------------------------------------------------------------------------------
for q = 1:numel(outs)
    if isfile(outs{q})
        bak = strrep(outs{q}, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
        copyfile(outs{q}, bak);  fprintf('Old file copied to %s\n', bak);
    end
end
blank_ms_original = SP.blank_ms;  goodCh = SP.goodCh;  noise_uV = SP.noise_uV;  thr_uV = SP.thr_uV;
preMs = SP.preMs;  postMs = SP.postMs;  nTrialsProcessed = nTr;
Recovery_params = struct('boundaryMode', boundaryMode, 'simMode', simMode, 'suffix', suffix, ...
    'sim_default', sim_default, 'sim_override', sim_override, 'seq_default', seq_default, 'seq_override', seq_override, ...
    'rand_seed', rand_seed, 'fallback_ms', fbK, 'cap_ms', capK, 'percent', pctK, 'event_ms', evK, 'spacing_ms', spK, ...
    'fstFile', fstFile, 'pairs', {pairs});
version = codeVersion;  date = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<NASGU>
save(outFile, 'sp_clipped', 'sp_trial', 'nSpk', 'blank_ms', 'blank_ms_original', 'rec_ms', 'srcTrial', 'isRecovered', ...
     'condTable', 'depth', 'recLog', 'Recovery_params', 'goodCh', 'noise_uV', 'thr_uV', 'FS', 'preMs', 'postMs', ...
     'nTrialsProcessed', 'spikeFile', 'ttFile', 'version', 'date', '-v7');
fprintf('Saved %s\n', outFile);
if saveWaveforms
    wfScale_uV = WF.wfScale_uV;  wfTime_ms = WF.wfTime_ms; %#ok<NASGU>
    save(wfOut, 'sp_waveforms', 'wfScale_uV', 'wfTime_ms', 'FS', 'version', 'date', '-v7.3');
    fprintf('Saved %s\n', wfOut);
end
fprintf('Done in %.0f s.\n', toc(tStart));

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
% (identical to NP_3a_FirstSpikeTime_v1)
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

function chkSet_(nm, def, ov)
assert(isnumeric(def) && numel(def) == 4, '%s_default must be [fallback_ms cap_ms percent event_ms].', nm);
assert(def(3) >= 0 && def(3) <= 100, '%s_default: percent must be 0..100.', nm);
if isempty(ov), return; end
assert(isnumeric(ov) && size(ov, 2) == 5, '%s_override rows must be [level fallback_ms cap_ms percent event_ms].', nm);
assert(all(ov(:,1) >= 2 & ov(:,1) == round(ov(:,1))), '%s_override: levels must be whole numbers >= 2 (level 1 is never recovered).', nm);
assert(numel(unique(ov(:,1))) == size(ov, 1), '%s_override: a level is listed twice.', nm);
assert(all(ov(:,4) >= 0 & ov(:,4) <= 100), '%s_override: percent must be 0..100.', nm);
end

function [relBy, wfBy] = groupByTrial_(rel, kk, nTr, doWf, W)
relBy = cell(nTr, 1);  wfBy = cell(nTr, 1);
for t = 1:nTr, relBy{t} = zeros(0, 1); end
if doWf
    for t = 1:nTr, wfBy{t} = zeros(0, size(W, 2), 'like', W); end
end
if isempty(kk), return; end
[kkS, o] = sort(kk);  relS = rel(o);
edges = [0; find(diff(kkS)); numel(kkS)];
for i = 1:numel(edges) - 1
    rows = edges(i)+1 : edges(i+1);  t = kkS(rows(1));
    [relBy{t}, oo] = sort(relS(rows));
    if doWf, w = W(o(rows), :); wfBy{t} = w(oo, :); end
end
end

function m = inIv_(t, R, h)
% true for times inside any [start end] row of R, widened by h on both sides (h = half a sample)
m = false(size(t));
for b = 1:size(R, 1), m = m | (t >= R(b,1) - h & t <= R(b,2) + h); end
end

function C = subIv_(A, B)
% intervals A minus intervals B (rows [start end])
C = A;
if isempty(A) || isempty(B), return; end
for b = 1:size(B, 1)
    D = zeros(0, 2);
    for a = 1:size(C, 1)
        s = C(a,1);  e = C(a,2);
        if B(b,2) <= s || B(b,1) >= e, D(end+1, :) = [s e]; continue; end %#ok<AGROW>
        if B(b,1) > s, D(end+1, :) = [s B(b,1)]; end %#ok<AGROW>
        if B(b,2) < e, D(end+1, :) = [B(b,2) e]; end %#ok<AGROW>
    end
    C = D;
end
end
