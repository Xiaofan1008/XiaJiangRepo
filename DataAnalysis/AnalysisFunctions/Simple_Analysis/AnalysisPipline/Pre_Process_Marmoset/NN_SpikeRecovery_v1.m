 %% ============================================================
% NN_SpikeRecovery_v1.m
% Step 2 of the PULSE-TRAIN artifact recovery pipeline (marmoset NN data).
%
% CHAINED recovery. Levels are processed in ascending order. For a target
% trial of level L, the window [0, blank_end) after ITS trigger is blanked
% (every original spike inside is deleted) and filled with the spikes found
% in the SAME relative window of a paired trial of its source level L-1,
% taken from the already RECOVERED data (so level 3 uses recovered level 2,
% which used recovered level 1 ...). Level 1 (no source) is never changed.
% Blank and visual trials are never changed.
%
% blank_end (per channel/trial) = first real spike after the last pulse
%   (from NN_FirstSpikeTime_v1) capped at lastPulse + cap_ms;
%   if no spike was found: lastPulse + fallback_ms.
% Pairing: i-th trial of the target condition <-> i-th trial of the source
%   condition (cycling if the source has fewer trials).
%
% Every condition is its own level, so every level's recovered result is
% stored in place in the output (one condition = one level). Use the
% condTable / trialCond saved in the output to pick levels afterwards.
%
% Inputs are never modified. New files written next to the input:
%   *.sp_xia_TrainRecovery.mat            sp_clipped (times only), blank_end_ms,
%                                              srcTrial, recLog, Recovery_params, FS
%   *.sp_xia_waveforms_TrainRecovery.mat  sp_waveforms (only if do_waveforms = true;
%                                              delete it after waveform checking)
% A spike whose time relative to its trigger is < blank_end_ms{ch}(trial)
% is a SUBSTITUTED spike; >= is original, untouched data.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
analysisPath = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/';

data_folder = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses';

times_file = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses/Xia_Train_Freq4_100Hz_20Pulses.sp_xia_filtered.mat';
                   % the FILTERED spike-times file of this dataset (same one used in step 1).
                   % The other files are found next to it by name:
                   %   *.sp_xia_waveforms_filtered.mat  (if do_waveforms)
                   %   *.sp_xia_FirstSpikeTimes.mat     (from step 1)

do_waveforms = true;   % true : recover waveform rows too, save the separate waveform file
                       % false: recover spike times only (smaller, faster), no waveform file

fallback_ms = 10;       % no spike found after the last pulse: blank until lastPulse + this (ms)
cap_ms      = 10;      % never blank later than lastPulse + this (ms)
FS = 30000;

%% ================= LOAD =================
if ~isempty(analysisPath), addpath(genpath(analysisPath)); end
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
cd(data_folder);

assert(contains(times_file, '.sp_xia_filtered.mat'), 'times_file must end in .sp_xia_filtered.mat');
wf_file  = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_waveforms_filtered.mat');
fst_file = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_FirstSpikeTimes.mat');
assert(isfile(fst_file), 'Cannot find %s -- run NN_FirstSpikeTime_v1 first.', fst_file);

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

%% ================= CHAINED RECOVERY =================
blank_end_ms = cell(nChn,1);
for ch = 1:nChn, blank_end_ms{ch} = nan(n_Trials,1); end
srcTrial = nan(n_Trials,1);

% recLog(k): per condition counts
recLog = repmat(struct('cond',0,'depth',0,'source',NaN,'nTargets',0,'nDetected',0,'nFallback',0,'nAdded',0,'nDeleted',0), nCond, 1);
for k = 1:nCond
    recLog(k).cond = k; recLog(k).depth = depth(k); recLog(k).source = srcOf(k);
end

for d = 1:maxDepth
    condsD = find(depth == d)';
    fprintf('\n--- Level depth %d: conditions %s ---\n', d, mat2str(condsD));

    % pair every target trial with a source trial (same for all channels)
    pairs = cell(nCond,1);
    for k = condsD
        trT = find(trialCond == k & trigOK & isfinite(lastPulse_ms));
        trS = find(trialCond == srcOf(k) & trigOK);
        if isempty(trT), continue; end
        if isempty(trS)
            warning('NSR:nosrc', 'cond %d: source cond %d has no usable trial. Left unrecovered.', k, srcOf(k));
            continue;
        end
        sIdx = mod((1:numel(trT))-1, numel(trS)) + 1;
        sSel = trS(sIdx);
        pairs{k} = [trT(:), sSel(:)];
        srcTrial(trT) = sSel;
        recLog(k).nTargets = numel(trT);
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
            for q = 1:size(P,1)
                trT = P(q,1); trS = P(q,2);
                first_ms = firstSpikeTimes{ch}(trT);
                lp = lastPulse_ms(trT);
                if isfinite(first_ms) && first_ms > 0
                    win_end = min(first_ms, lp + cap_ms);
                    recLog(k).nDetected = recLog(k).nDetected + 1;
                else
                    win_end = lp + fallback_ms;
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
        fprintf(' cond %3d (source %3d): %3d trials | window end from spike %d, fallback %d | deleted %d, added %d spikes\n', ...
            k, srcOf(k), recLog(k).nTargets, recLog(k).nDetected, recLog(k).nFallback, recLog(k).nDeleted, recLog(k).nAdded);
    end
end

%% ================= SAVE =================
Recovery_params = struct('times_file',times_file,'do_waveforms',do_waveforms,'fallback_ms',fallback_ms,'cap_ms',cap_ms, ...
    'fst_file',fst_file,'condTable',condTable,'trialCond',trialCond,'depth',depth,'sourceCond',srcOf);
if do_waveforms
    sp_clipped = cellfun(@col1_, rows, 'UniformOutput', false);
    sp_waveforms = rows;
else
    sp_clipped = rows;
end
times_out = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_TrainRecovery.mat');
save(times_out, 'sp_clipped', 'blank_end_ms', 'srcTrial', 'recLog', 'Recovery_params', 'FS', '-v7.3');
fprintf('\nSaved recovered spike times to %s\n', times_out);
if do_waveforms
    wf_out = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_waveforms_TrainRecovery.mat');
    save(wf_out, 'sp_waveforms', 'blank_end_ms', 'FS', '-v7.3');
    fprintf('Saved recovered waveforms to %s (delete after waveform checking)\n', wf_out);
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
