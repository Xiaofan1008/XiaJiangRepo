%% ============================================================
% NN_SpikeRecovery_Random_v1.m
% Step 2 of the PULSE-TRAIN artifact recovery pipeline (marmoset NN data).
% Variant of NN_SpikeRecovery_v1 with PER-LEVEL settings and a RANDOM SUBSET of trials per level.
% NN_SpikeRecovery_v1 and its output files are NOT touched.
%
% WHAT IS NEW COMPARED WITH v1
%   1. fallback time, cap (max) time and the % of trials to recover can be set PER LEVEL,
%      separately for the SIMULTANEOUS pattern (sim: levels 1..n) and the SEQUENTIAL pattern
%      (seq: levels 1..2n; AB and BA share one block). Settings start at level 2 (level 1 has
%      no source and is never changed).
%   2. At every level only a random X% of the valid trials is recovered. The other trials stay
%      EXACTLY as the original data (with the blank gap) and are flagged isRecovered = false.
%   3. The paired source trial must itself be recovered (level 1 trials always count as
%      recovered), so a higher level never copies a gap. If there are fewer sources than
%      targets, sources are reused cyclically (as in v1).
%   4. Output files carry the suffix _Random so the v1 files are never overwritten.
%
% HOW A LEVEL IS RECOVERED (unchanged logic)
%   Levels are processed in ascending order. For a recovered target trial of level L the window
%   [0, blank_end) after ITS trigger is blanked (every original spike inside is deleted) and
%   filled with the spikes found in the SAME relative window of its paired trial of level L-1,
%   taken from the already RECOVERED data.
%   blank_end (per channel/trial) = first real spike after the last pulse (NN_FirstSpikeTime_v1)
%   capped at lastPulse + cap_ms; if no spike was found: lastPulse + fallback_ms.
%   Pairing: i-th recovered target trial <-> i-th eligible source trial (cycling).
%
% SETTINGS PER LEVEL
%   sim_default / seq_default  = [fallback_ms  cap_ms  percent]      used for every level of that pattern
%   sim_override / seq_override = rows [level fallback_ms cap_ms percent] for levels that differ, e.g.
%       seq_override = [ 2  4  8  50
%                        3  6 10  80 ];
%   Levels not listed use the default. A level is the pattern's own level number (condTable.level).
%   percent = share of the VALID trials of that level that are recovered (rounded, at least 1 if > 0).
%   The trials are chosen at random with the seed rand_seed (same trials for all channels);
%   the same seed always gives the same choice.
%
% OUTPUT (new files next to the input)
%   *.sp_xia_TrainRecovery_<out_suffix>.mat          sp_clipped (times only), blank_end_ms, srcTrial,
%                                                      isRecovered, recLog, Recovery_params, FS
%   *.sp_xia_waveforms_TrainRecovery_<out_suffix>.mat sp_waveforms (only if do_waveforms = true)
%   blank_end_ms is NaN for trials that were not recovered. A spike whose time relative to its
%   trigger is < blank_end_ms{ch}(trial) is a SUBSTITUTED spike.
%
% CHECK AGAINST v1: set compare_v1_file to the v1 *.sp_xia_TrainRecovery.mat. If every percent is
%   100 and the fallback / cap equal the ones of the v1 run, the result must be identical.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
analysisPath = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/';

data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';

times_file = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_filtered.mat';
                   % the FILTERED spike-times file of this dataset (same one used in step 1).
                   % The other files are found next to it by name:
                   %   *.sp_xia_waveforms_filtered.mat  (if do_waveforms)
                   %   *.sp_xia_FirstSpikeTimes.mat     (from step 1)

out_suffix   = 'Random';   % output files: *.sp_xia_TrainRecovery_<out_suffix>.mat
do_waveforms = true;       % true : recover waveform rows too, save the separate waveform file
                           % false: recover spike times only (smaller, faster), no waveform file
FS = 30000;

% ---- settings per level: [fallback_ms  cap_ms  percent]  (levels start at 2) ----
sim_default  = [10 10 100];     % simultaneous pattern, every level not listed in sim_override
sim_override = [5,9,9,90;
                6,8,8,80;
                7,7,7,70;
                9,6,6,60;
               10,5,5,50;
               11,4,4,40;
               12,3,3,30;
               13,2,2,20;
               14,1,1,10;
               15,0.5,0.5,5;
];              % rows [level fallback_ms cap_ms percent]
seq_default  = [10 10 100];     % sequential pattern (AB and BA), every level not listed in seq_override
seq_override = [15,8,8,80;
                16,8,8,80;
                17,7,7,70;
                18,7,7,70;
                19,6,6,60;
                20,6,6,60;
                21,5,5,50;
                22,5,5,50;
                23,4,4,40;
                24,4,4,40;
                25,3,3,30;
                26,3,3,30;
                27,2,2,20;
                28,2,2,20;
                29,1,1,10;
                30,1,1,10;
];              % rows [level fallback_ms cap_ms percent]

rand_seed = 1;                  % seed of the random choice of trials (percent < 100)

compare_v1_file = '';           % '' = no comparison; or the v1 *.sp_xia_TrainRecovery.mat (see header)

%% ================= LOAD =================
if ~isempty(analysisPath), addpath(genpath(analysisPath)); end
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
cd(data_folder);

assert(contains(times_file, '.sp_xia_filtered.mat'), 'times_file must end in .sp_xia_filtered.mat');
wf_file  = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_waveforms_filtered.mat');
fst_file = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_FirstSpikeTimes.mat');
assert(isfile(fst_file), 'Cannot find %s -- run NN_FirstSpikeTime_v1 first.', fst_file);
assert(~isempty(out_suffix) && ischar(out_suffix), 'out_suffix must be a non-empty text.');

if do_waveforms
    assert(isfile(wf_file), 'Cannot find %s.', wf_file);
    W = load(wf_file);
    assert(isfield(W,'sp_waveforms'), 'sp_waveforms not found in %s.', wf_file);
    rows = W.sp_waveforms;  clear W
    if isfile(times_file) && ~isempty(whos('-file', times_file, 'fs')), tmp = load(times_file,'fs'); FS = tmp.fs; end
else
    assert(isfile(times_file), 'Cannot find %s.', times_file);
    T = load(times_file);
    assert(isfield(T,'sp_clipped'), 'sp_clipped not found in %s.', times_file);
    rows = T.sp_clipped;
    if isfield(T,'fs'), FS = T.fs; end
    clear T
end
nChn = numel(rows);

Sf = load(fst_file);
need = {'firstSpikeTimes','lastPulse_ms','trialCond','condTable','n_Trials'};
assert(all(isfield(Sf, need)), '%s is missing expected fields -- re-run NN_FirstSpikeTime_v1.', fst_file);
firstSpikeTimes = Sf.firstSpikeTimes;
lastPulse_ms    = Sf.lastPulse_ms(:);
trialCond       = Sf.trialCond(:);
condTable       = Sf.condTable;
n_Trials        = Sf.n_Trials;
assert(numel(firstSpikeTimes) == nChn, 'FirstSpikeTimes has %d channels, spike file has %d.', numel(firstSpikeTimes), nChn);
if isfield(Sf,'FS') && Sf.FS ~= FS, warning('NSR:fs', 'FS differs between files (%g vs %g).', Sf.FS, FS); end

% triggers (same reading as step 1; skips macOS "._" files)
f = dir('*.trig.dat'); f = f(~startsWith({f.name}, '.'));
assert(~isempty(f), 'No *.trig.dat in %s.', data_folder);
fid = fopen(f(1).name, 'r'); trig = fread(fid, [1 inf], 'double'); fclose(fid);
trig = trig(:);
trigOK = false(n_Trials,1);
nUse = min(numel(trig), n_Trials);
trigOK(1:nUse) = trig(1:nUse) > 0;
t0_ms = nan(n_Trials,1); t0_ms(1:nUse) = trig(1:nUse)/FS*1000;

%% ================= LEVEL ORDER (CHAIN DEPTH) =================
nCond  = numel(condTable);
srcOf  = [condTable.sourceCond]';            % NaN = no source
depth  = zeros(nCond,1);
for it = 1:nCond+1
    changed = false;
    for k = 1:nCond
        if ~isnan(srcOf(k))
            d = depth(srcOf(k)) + 1;
            if d ~= depth(k), depth(k) = d; changed = true; end
        end
    end
    if ~changed, break; end
    if it == nCond+1, error('Source chain has a loop -- check the condTable.'); end
end
maxDepth = max(depth);
fprintf('%d conditions, chain depth 0..%d (0 = not recovered).\n', nCond, maxDepth);

%% ================= SETTINGS PER CONDITION =================
chkSet_('sim', sim_default, sim_override);
chkSet_('seq', seq_default, seq_override);
assert(rand_seed == round(rand_seed), 'rand_seed must be an integer.');

kindK = arrayfun(@(c) char(c.kind), condTable, 'UniformOutput', false);
lvlK  = arrayfun(@(c) double(c.level), condTable);
fbK = nan(nCond,1);  capK = nan(nCond,1);  pctK = nan(nCond,1);
for k = 1:nCond
    if depth(k) == 0, continue; end
    if strcmp(kindK{k}, 'sim')
        def = sim_default;  ov = sim_override;
    else
        if ~any(strcmp(kindK{k}, {'AB','BA'}))
            warning('NSR:kind', 'cond %d has kind ''%s'': the seq settings are used.', k, kindK{k});
        end
        def = seq_default;  ov = seq_override;
    end
    s = def;
    if ~isempty(ov)
        r = find(ov(:,1) == lvlK(k), 1);
        if ~isempty(r), s = ov(r,2:4); end
    end
    fbK(k) = s(1);  capK(k) = s(2);  pctK(k) = s(3);
end
% overrides for levels that do not exist
for pp = 1:2
    if pp == 1, nm = 'sim'; ov = sim_override; isP = strcmp(kindK, 'sim');
    else,       nm = 'seq'; ov = seq_override; isP = ~strcmp(kindK, 'sim') & depth(:)' > 0; end
    if isempty(ov), continue; end
    mx = max([0, lvlK(isP & depth(:)' > 0)]);
    for r = 1:size(ov,1)
        if ov(r,1) > mx, warning('NSR:level', '%s_override: level %g does not exist (max level %d).', nm, ov(r,1), mx); end
    end
end

fprintf('\n--- SETTINGS PER CONDITION ---\n  cond  kind  level  source   fallback(ms)  cap(ms)  percent\n');
for k = find(depth(:)' > 0)
    fprintf('  %3d   %-4s   %3d    %3d      %6.2f      %6.2f    %5.1f\n', k, kindK{k}, lvlK(k), srcOf(k), fbK(k), capK(k), pctK(k));
end
fprintf('random seed %d | output suffix _%s\n', rand_seed, out_suffix);

%% ================= CHAINED RECOVERY =================
blank_end_ms = cell(nChn,1);
for ch = 1:nChn, blank_end_ms{ch} = nan(n_Trials,1); end
srcTrial = nan(n_Trials,1);
isRecovered = false(n_Trials,1);

% a trial may be used as a SOURCE if it is a chain root (depth 0) or was itself recovered
okSrc = false(n_Trials,1);
for k = 1:nCond
    if depth(k) == 0, okSrc(trialCond == k) = trigOK(trialCond == k); end
end

% recLog(k): per condition counts
recLog = repmat(struct('cond',0,'depth',0,'source',NaN,'fallback_ms',NaN,'cap_ms',NaN,'percent',NaN, ...
    'nValid',0,'nTargets',0,'relRecovered',[],'nDetected',0,'nFallback',0,'nAdded',0,'nDeleted',0), nCond, 1);
for k = 1:nCond
    recLog(k).cond = k; recLog(k).depth = depth(k); recLog(k).source = srcOf(k);
    recLog(k).fallback_ms = fbK(k); recLog(k).cap_ms = capK(k); recLog(k).percent = pctK(k);
end

for d = 1:maxDepth
    condsD = find(depth == d)';
    fprintf('\n--- Level depth %d: conditions %s ---\n', d, mat2str(condsD));

    % choose the target trials (random subset) and pair them with eligible source trials (same for all channels)
    pairs = cell(nCond,1);
    for k = condsD
        trAll = find(trialCond == k & trigOK & isfinite(lastPulse_ms));
        recLog(k).nValid = numel(trAll);
        if isempty(trAll), continue; end
        nSel = round(pctK(k)/100*numel(trAll));
        if pctK(k) > 0 && nSel == 0, nSel = 1; end
        if nSel >= numel(trAll)
            trT = trAll(:);
        elseif nSel <= 0
            trT = zeros(0,1);
        else
            rs = RandStream('mt19937ar', 'Seed', rand_seed + k);
            trT = sort(trAll(randperm(rs, numel(trAll), nSel)));
            trT = trT(:);
        end
        if isempty(trT), continue; end
        trS = find(trialCond == srcOf(k) & trigOK & okSrc);
        if isempty(trS)
            warning('NSR:nosrc', 'cond %d: source cond %d has no recovered trial. Left unrecovered.', k, srcOf(k));
            continue;
        end
        sIdx = mod((1:numel(trT))-1, numel(trS)) + 1;
        sSel = trS(sIdx);
        pairs{k} = [trT(:), sSel(:)];
        srcTrial(trT) = sSel;
        isRecovered(trT) = true;
        okSrc(trT) = true;
        recLog(k).nTargets = numel(trT);
        recLog(k).relRecovered = find(ismember(find(trialCond == k), trT))';
    end

    for ch = 1:nChn
        R = rows{ch};
        if isempty(R) && ~any(cellfun(@(p) ~isempty(p), pairs)), continue; end
        tt = R(:,1);                                  % sorted spike times (ms)
        delMask = false(size(R,1),1);
        adds = {};
        for k = condsD
            P = pairs{k};
            if isempty(P), continue; end
            fb = fbK(k);  cp = capK(k);
            for q = 1:size(P,1)
                trT = P(q,1); trS = P(q,2);
                first_ms = firstSpikeTimes{ch}(trT);
                lp = lastPulse_ms(trT);
                if isfinite(first_ms) && first_ms > 0
                    win_end = min(first_ms, lp + cp);
                    recLog(k).nDetected = recLog(k).nDetected + 1;
                else
                    win_end = lp + fb;
                    recLog(k).nFallback = recLog(k).nFallback + 1;
                end
                if win_end <= 0, continue; end
                blank_end_ms{ch}(trT) = win_end;

                % delete original target spikes in [0, win_end)
                a = t0_ms(trT);  b = a + win_end;
                i1 = lb_(tt, a);  i2 = lb_(tt, b);          % rows i1..i2-1 are inside
                if i2 > i1
                    delMask(i1:i2-1) = true;
                    recLog(k).nDeleted = recLog(k).nDeleted + (i2 - i1);
                end

                % spikes of the (already recovered) source trial, same relative window
                sa = t0_ms(trS);  sb = sa + win_end;
                j1 = lb_(tt, sa);  j2 = lb_(tt, sb);
                if j2 > j1
                    add = R(j1:j2-1, :);
                    add(:,1) = add(:,1) - sa + a;
                    adds{end+1} = add; %#ok<AGROW>
                    recLog(k).nAdded = recLog(k).nAdded + size(add,1);
                end
            end
        end
        if isempty(adds) && ~any(delMask), continue; end
        R(delMask,:) = [];
        if ~isempty(adds), R = [R; vertcat(adds{:})]; end %#ok<AGROW>
        [~,o] = sort(R(:,1));
        rows{ch} = R(o,:);
    end

    for k = condsD
        fprintf(' cond %3d %-4s lvl %2d (source %3d) fb %g cap %g pct %g%%: %2d / %2d trials | window end from spike %d, fallback %d | deleted %d, added %d spikes\n', ...
            k, kindK{k}, lvlK(k), srcOf(k), fbK(k), capK(k), pctK(k), recLog(k).nTargets, recLog(k).nValid, ...
            recLog(k).nDetected, recLog(k).nFallback, recLog(k).nDeleted, recLog(k).nAdded);
        if recLog(k).nTargets < recLog(k).nValid
            fprintf('            recovered trials (number within the condition): %s\n', fmtList_(recLog(k).relRecovered));
        end
    end
end

%% ================= SUMMARY =================
fprintf('\n--- RECOVERED TRIALS PER LEVEL (recovered / valid) ---\n');
for pp = 1:2
    if pp == 1, nm = 'sim'; sel = find(depth(:)' > 0 & strcmp(kindK(:)', 'sim'));
    else,       nm = 'seq'; sel = find(depth(:)' > 0 & ~strcmp(kindK(:)', 'sim')); end
    if isempty(sel), continue; end
    fprintf(' %s: ', nm);
    for k = sel, fprintf('L%d %d/%d  ', lvlK(k), recLog(k).nTargets, recLog(k).nValid); end
    fprintf('\n');
end
fprintf('Not recovered trials keep the original data (blank gap) and have isRecovered = false.\n');

%% ================= SAVE =================
Recovery_params = struct('times_file',times_file,'do_waveforms',do_waveforms,'fst_file',fst_file, ...
    'condTable',condTable,'trialCond',trialCond,'depth',depth,'sourceCond',srcOf, ...
    'fallback_ms_k',fbK,'cap_ms_k',capK,'percent_k',pctK,'rand_seed',rand_seed, ...
    'sim_default',sim_default,'sim_override',sim_override,'seq_default',seq_default,'seq_override',seq_override, ...
    'isRecovered',isRecovered,'out_suffix',out_suffix);
if do_waveforms
    sp_clipped = cellfun(@col1_, rows, 'UniformOutput', false);
    sp_waveforms = rows;
else
    sp_clipped = rows;
end
times_out = strrep(times_file, '.sp_xia_filtered.mat', ['.sp_xia_TrainRecovery_' out_suffix '.mat']);
save(times_out, 'sp_clipped', 'blank_end_ms', 'srcTrial', 'isRecovered', 'recLog', 'Recovery_params', 'FS', '-v7.3');
fprintf('\nSaved recovered spike times to %s\n', times_out);
if do_waveforms
    wf_out = strrep(times_file, '.sp_xia_filtered.mat', ['.sp_xia_waveforms_TrainRecovery_' out_suffix '.mat']);
    save(wf_out, 'sp_waveforms', 'blank_end_ms', 'isRecovered', 'FS', '-v7.3');
    fprintf('Saved recovered waveforms to %s (delete after waveform checking)\n', wf_out);
end

%% ================= OPTIONAL: COMPARE WITH v1 =================
if ~isempty(compare_v1_file)
    assert(isfile(compare_v1_file), 'compare_v1_file not found: %s', compare_v1_file);
    C = load(compare_v1_file, 'sp_clipped', 'Recovery_params');
    fprintf('\n--- COMPARISON WITH v1 (%s) ---\n', compare_v1_file);
    if isfield(C.Recovery_params,'fallback_ms') && isfield(C.Recovery_params,'cap_ms')
        fprintf('v1 run used fallback %g ms, cap %g ms.\n', C.Recovery_params.fallback_ms, C.Recovery_params.cap_ms);
    end
    if any(pctK(depth > 0) < 100)
        fprintf('Some percent values are < 100: the result is NOT expected to equal v1.\n');
    end
    nSame = 0;
    for ch = 1:nChn
        a_ = sp_clipped{ch};  b_ = C.sp_clipped{ch};
        if isequal(a_(:), b_(:)), nSame = nSame + 1;
        else, fprintf('  channel %2d differs: %d vs %d spikes\n', ch, numel(a_), numel(b_)); end
    end
    fprintf('identical in %d of %d channels\n', nSame, nChn);
end

%% ================= HELPERS =================
function i = lb_(s, x)
% first index with s(i) >= x (numel(s)+1 if none); s sorted ascending
lo = 1; hi = numel(s) + 1;
while lo < hi
    m = floor((lo + hi)/2);
    if s(m) < x, lo = m + 1; else, hi = m; end
end
i = lo;
end

function c = col1_(x)
if isempty(x), c = x; else, c = x(:,1); end
end

function chkSet_(name, def, ov)
% validate the settings block of one pattern
assert(isnumeric(def) && numel(def) == 3, '%s_default must be [fallback_ms cap_ms percent].', name);
assert(def(1) > 0 && def(2) > 0, '%s_default: fallback_ms and cap_ms must be > 0.', name);
assert(def(3) >= 0 && def(3) <= 100, '%s_default: percent must be between 0 and 100.', name);
if isempty(ov), return; end
assert(isnumeric(ov) && size(ov,2) == 4, '%s_override must have 4 columns: [level fallback_ms cap_ms percent].', name);
assert(all(ov(:,1) >= 2 & ov(:,1) == round(ov(:,1))), '%s_override: levels must be integers >= 2 (level 1 is never recovered).', name);
assert(numel(unique(ov(:,1))) == size(ov,1), '%s_override: a level is listed more than once.', name);
assert(all(ov(:,2) > 0 & ov(:,3) > 0), '%s_override: fallback_ms and cap_ms must be > 0.', name);
assert(all(ov(:,4) >= 0 & ov(:,4) <= 100), '%s_override: percent must be between 0 and 100.', name);
end

function s = fmtList_(v)
% compact MATLAB-style list, e.g. [1 2 4:7]
v = sort(v(:)');
if isempty(v), s = '[]'; return; end
parts = {};  i = 1;
while i <= numel(v)
    j = i;
    while j < numel(v) && v(j+1) == v(j) + 1, j = j + 1; end
    if j == i,        parts{end+1} = sprintf('%d', v(i));          %#ok<AGROW>
    elseif j == i+1,  parts{end+1} = sprintf('%d %d', v(i), v(j)); %#ok<AGROW>
    else,             parts{end+1} = sprintf('%d:%d', v(i), v(j)); %#ok<AGROW>
    end
    i = j + 1;
end
s = ['[' strjoin(parts, ' ') ']'];
end
