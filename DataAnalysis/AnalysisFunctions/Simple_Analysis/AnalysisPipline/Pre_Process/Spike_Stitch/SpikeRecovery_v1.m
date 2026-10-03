%% ============================================================
% SpikeRecovery_SequentialStitch_v1.m
% Step 2 of the sequential-stimulation artifact recovery pipeline.
%
% For each SEQUENTIAL trial (PTD > 0) only, the window [0, blank_end)
% relative to that trial's own trigger is blanked and REPLACED with the
% matching single-pulse (stim ch1 only) recording's own response over
% the same relative window, time-shifted into place. blank_end is the
% first real spike detected after the 2nd pulse (from
% FirstSpikeTime_Sequential_v2.m), falling back to
% PTD + artifact_blank_fallback_ms when no spike was detected for that
% channel/trial.
%
% SIMULTANEOUS trials (PTD == 0) are NEVER selected as a target trial
% here -- the trial-selection mask below requires ~isSimultaneous(tr),
% exactly as in FirstSpikeTime_Sequential_v2.m, so their rows are never
% deleted, never substituted, and pass through to the output completely
% unchanged, byte for byte.
%
% Operates on the FULL waveform files (time + waveform shape together,
% one row per spike) for BOTH datasets, so time and waveform shape can
% never drift apart from each other. The injection logic runs ONCE, on
% these combined rows; the companion times-only output
% (<base>.sp_xia_stitched.mat / sp_clipped) is then simply column 1 of
% the very same result -- not a second, independently-run pass -- so
% the two output files are guaranteed to agree with each other.
%
% Neither input file is modified. Two NEW files are written:
%   <base_name>.sp_xia_stitched.mat           -- sp_clipped only
%   <base_name>.sp_xia_waveforms_stitched.mat -- sp_waveforms only
% Both also carry 'fs', 'Stitch_params', and 'blank_end_ms' (per
% channel/trial, the actual substitution-window end used -- NaN for
% simultaneous trials and any trial/channel never blanked). A spike's
% relative time < blank_end_ms means it's a SUBSTITUTED spike; >= means
% it's original, untouched data. This is what the QC comparison script
% (next step) uses to color-split and mark the stitch boundary, without
% needing a separate per-spike flag.
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
% -------- Single-pulse (source) dataset --------
single_wf_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_waveforms_filtered.mat';
single_folder  = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1';

% -------- Sequential (target) dataset --------
seq_wf_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_waveforms_filtered.mat';
seq_folder  = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1';

% -------- First-spike-time file (from FirstSpikeTime_Sequential_v2.m) --------
fst_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_FirstSpikeTimes.mat';

pulse_offset_ms = 0;     % extra shift applied when injecting into the
                         % sequential timeline (0 = none)
artifact_blank_fallback_ms = 2.5;
                         % NAMED fallback: when no real spike was
                         % detected after the 2nd pulse for a given
                         % channel/trial, the blank window is assumed to
                         % end at PTD + this many ms instead.
FS = 30000;

%% ================= LOAD SINGLE-PULSE (SOURCE) DATASET =================
assert(isfile(single_wf_file), 'Cannot find %s.', single_wf_file);
S1w = load(single_wf_file);
assert(isfield(S1w, 'sp_waveforms'), 'Variable "sp_waveforms" not found in %s.', single_wf_file);
sp_waveforms_single = S1w.sp_waveforms;
if isfield(S1w, 'fs'), FS = S1w.fs; end

cd(single_folder);
trig_single = loadTrig(0);
pfile1 = dir(fullfile(single_folder, '*_exp_datafile_*.mat'));
assert(~isempty(pfile1), 'No *_exp_datafile_*.mat found in %s.', single_folder);
S1p = load(fullfile(single_folder, pfile1(1).name), 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_single  = S1p.StimParams;
E_MAP_single       = S1p.E_MAP;
n_Trials_single    = S1p.n_Trials;
simultaneous_stim1 = S1p.simultaneous_stim;

trialAmps_single = cell2mat(StimParams_single(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);
stimNames_single = StimParams_single(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP_single(2:end));
stimChPerTrial_single = cell(n_Trials_single,1);
for t = 1:n_Trials_single
    rr = (t-1)*simultaneous_stim1 + (1:simultaneous_stim1);
    v  = idx_all_single(rr);
    v  = v(v>0);
    stimChPerTrial_single{t} = v(:).';
end

%% ================= LOAD SEQUENTIAL (TARGET) DATASET =================
assert(isfile(seq_wf_file), 'Cannot find %s.', seq_wf_file);
S2w = load(seq_wf_file);
assert(isfield(S2w, 'sp_waveforms'), 'Variable "sp_waveforms" not found in %s.', seq_wf_file);
sp_waveforms_seq = S2w.sp_waveforms;
nChn = numel(sp_waveforms_seq);
assert(numel(sp_waveforms_single) == nChn, ...
    'Single-pulse (%d channels) and sequential (%d channels) waveform files disagree.', ...
    numel(sp_waveforms_single), nChn);

cd(seq_folder);
trig_seq = loadTrig(0);
pfile2 = dir(fullfile(seq_folder, '*_exp_datafile_*.mat'));
assert(~isempty(pfile2), 'No *_exp_datafile_*.mat found in %s.', seq_folder);
S2p = load(fullfile(seq_folder, pfile2(1).name), 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_seq     = S2p.StimParams;
E_MAP_seq          = S2p.E_MAP;
n_Trials_seq       = S2p.n_Trials;
simultaneous_stim2 = S2p.simultaneous_stim;

trialAmps_seq = cell2mat(StimParams_seq(2:end,16));
trialAmps_seq = trialAmps_seq(1:simultaneous_stim2:end);
stimNames_seq = StimParams_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP_seq(2:end));
stimChPerTrial_seq = cell(n_Trials_seq,1);
for t = 1:n_Trials_seq
    rr = (t-1)*simultaneous_stim2 + (1:simultaneous_stim2);
    v  = idx_all_seq(rr);
    v  = v(v>0);
    stimChPerTrial_seq{t} = v(:).';
end
comb_seq = zeros(n_Trials_seq, simultaneous_stim2);
for t = 1:n_Trials_seq
    v = stimChPerTrial_seq{t};
    comb_seq(t, 1:numel(v)) = v;
end
[uniqueComb_seq, ~, combClass_seq] = unique(comb_seq, 'rows', 'stable');
nSeqSets = size(uniqueComb_seq,1);

%% ================= LOAD FIRST-SPIKE-TIME FILE =================
assert(isfile(fst_file), 'Cannot find %s.', fst_file);
Sf = load(fst_file);
assert(all(isfield(Sf, {'firstSpikeTimes','hasSpike','PTD_ms','isSimultaneous','n_Trials'})), ...
    '%s is missing expected fields -- re-run FirstSpikeTime_Sequential_v2.m.', fst_file);
firstSpikeTimes = Sf.firstSpikeTimes;
hasSpike        = Sf.hasSpike;
PTD_ms          = Sf.PTD_ms;
isSimultaneous  = Sf.isSimultaneous;
assert(Sf.n_Trials == n_Trials_seq, ...
    'FirstSpikeTimes file (%d trials) does not match the sequential dataset (%d trials) -- same run?', ...
    Sf.n_Trials, n_Trials_seq);
assert(numel(firstSpikeTimes) == nChn, ...
    'FirstSpikeTimes file (%d channels) does not match the waveform files (%d channels).', ...
    numel(firstSpikeTimes), nChn);

fprintf('%d/%d trials are simultaneous and will be left completely untouched.\n', ...
    sum(isSimultaneous), n_Trials_seq);

%% ================= STITCH (SEQUENTIAL TRIALS ONLY) =================
sp_waveforms_stitched = sp_waveforms_seq;   % start as an exact copy; only
                                             % sequential-trial windows
                                             % below are ever modified

% blank_end_ms{ch}(trial) records the ACTUAL win_end used for that
% channel/trial (whether from detection or the fallback), i.e. the
% boundary between "substituted" (relative time < blank_end_ms) and
% "original, untouched" (relative time >= blank_end_ms) spikes in the
% stitched output. NaN wherever no blanking was applied (including every
% simultaneous trial, which is never touched at all). This is what lets
% the QC comparison script color-split and mark the stitch boundary
% without needing a separate per-spike flag -- the boundary alone is
% enough, since deletion always removes every original spike below it.
blank_end_ms = cell(nChn,1);
for ch = 1:nChn
    blank_end_ms{ch} = nan(n_Trials_seq,1);
end

unique_amps = unique(trialAmps_seq);
unique_PTDs = unique(PTD_ms(~isSimultaneous));

total_spikes_added_all = 0;
n_windows_detected = 0;   % blank_end came from a real detected spike
n_windows_fallback = 0;   % blank_end came from the PTD+fallback default

for a = 1:numel(unique_amps)
    amp_val = unique_amps(a);
    total_spikes_added_amp = 0;
    mask_single_amp = (trialAmps_single == amp_val);

    for set_id = 1:nSeqSets
        stimVec = uniqueComb_seq(set_id,:);
        stimVec = stimVec(stimVec > 0);
        if isempty(stimVec), continue; end

        % Source = the FIRST stim channel of this set (ch1), matched by
        % NAME across the two datasets' E_MAPs (not by raw index).
        ch_first_idx_seq = stimVec(1);
        ch_name = E_MAP_seq{ch_first_idx_seq + 1};
        ch_idx_in_single_map = find(strcmp(E_MAP_single(2:end), ch_name));
        if isempty(ch_idx_in_single_map)
            fprintf('Warning: Seq Channel %s not found in Single Data Map.\n', ch_name);
            continue;
        end
        mask_single_chan = cellfun(@(x) ismember(ch_idx_in_single_map, x), stimChPerTrial_single);
        single_trials_for_group = find(mask_single_amp & mask_single_chan);
        if isempty(single_trials_for_group)
            fprintf('No single trials for amp %g, Stim %s (SeqID %d)\n', amp_val, ch_name, set_id);
            continue;
        end

        for p = 1:numel(unique_PTDs)
            ptd_val = unique_PTDs(p);

            % Sequential trials in this (amp, set, PTD) group --
            % ~isSimultaneous excludes PTD==0 trials explicitly, exactly
            % as in FirstSpikeTime_Sequential_v2.m, so a simultaneous
            % trial can never be selected as tr_seq below.
            mask_seq_group = (trialAmps_seq == amp_val) & ...
                              (combClass_seq == set_id) & ...
                              (PTD_ms == ptd_val) & ...
                              ~isSimultaneous;
            seq_trials_group = find(mask_seq_group);
            if isempty(seq_trials_group), continue; end

            for gi = 1:numel(seq_trials_group)
                tr_seq = seq_trials_group(gi);

                % deterministic cycling through the matching single-pulse trials
                idx_single = mod(gi-1, numel(single_trials_for_group)) + 1;
                tr_single  = single_trials_for_group(idx_single);

                t0_seq_ms    = trig_seq(tr_seq)/FS*1000;
                t0_single_ms = trig_single(tr_single)/FS*1000;

                for rec_ch = 1:nChn
                    spikes_src = sp_waveforms_single{rec_ch};
                    if isempty(spikes_src), continue; end

                    % -------- Determine this channel/trial's blank window end --------
                    first_ms = NaN;
                    if rec_ch <= numel(firstSpikeTimes) && tr_seq <= numel(firstSpikeTimes{rec_ch})
                        first_ms = firstSpikeTimes{rec_ch}(tr_seq);
                    end

                    win_start = 0;   % always blank from the FIRST pulse onset
                    if isfinite(first_ms) && first_ms > win_start
                        win_end = first_ms;
                        n_windows_detected = n_windows_detected + 1;
                    else
                        win_end = ptd_val + artifact_blank_fallback_ms;
                        n_windows_fallback = n_windows_fallback + 1;
                    end
                    if win_end <= win_start, continue; end

                    blank_end_ms{rec_ch}(tr_seq) = win_end;

                    % -------- Blank: delete existing sequential rows in this window --------
                    seq_rows = sp_waveforms_stitched{rec_ch};
                    if ~isempty(seq_rows)
                        rel_seq = seq_rows(:,1) - t0_seq_ms;
                        del_mask = (rel_seq >= win_start) & (rel_seq < win_end);
                        if any(del_mask)
                            seq_rows(del_mask,:) = [];
                            sp_waveforms_stitched{rec_ch} = seq_rows;
                        end
                    end

                    % -------- Substitute: single-pulse rows over the same relative window --------
                    t_start = t0_single_ms + win_start;
                    t_end   = t0_single_ms + win_end;
                    in_win = (spikes_src(:,1) >= t_start) & (spikes_src(:,1) < t_end);
                    rows_add = spikes_src(in_win,:);
                    if isempty(rows_add), continue; end

                    t_inject = t0_seq_ms + pulse_offset_ms;
                    rows_add(:,1) = rows_add(:,1) - t0_single_ms + t_inject;

                    sp_waveforms_stitched{rec_ch} = sortrows([sp_waveforms_stitched{rec_ch}; rows_add], 1);
                    total_spikes_added_amp = total_spikes_added_amp + size(rows_add,1);
                end
            end
        end
    end

    fprintf('Amplitude %g uA: %d spikes added.\n', amp_val, total_spikes_added_amp);
    total_spikes_added_all = total_spikes_added_all + total_spikes_added_amp;
end

fprintf('\nTotal spikes added across all amplitudes: %d\n', total_spikes_added_all);
fprintf('Blank-window end came from a detected spike in %d cases, from the %g ms fallback in %d cases.\n', ...
    n_windows_detected, artifact_blank_fallback_ms, n_windows_fallback);

%% ================= DERIVE TIMES-ONLY OUTPUT FROM THE SAME RESULT =================
sp_clipped_stitched = cellfun(@(x) ifelse_col1(x), sp_waveforms_stitched, 'UniformOutput', false);

%% ================= SAVE (two files, one computation) =================
Stitch_params = struct();
Stitch_params.single_wf_file  = single_wf_file;
Stitch_params.seq_wf_file     = seq_wf_file;
Stitch_params.fst_file        = fst_file;
Stitch_params.pulse_offset_ms = pulse_offset_ms;
Stitch_params.artifact_blank_fallback_ms = artifact_blank_fallback_ms;
Stitch_params.n_windows_detected = n_windows_detected;
Stitch_params.n_windows_fallback = n_windows_fallback;

wf_out = strrep(seq_wf_file, '.sp_xia_waveforms_filtered.mat', '.sp_xia_waveforms_stitched.mat');
if strcmp(wf_out, seq_wf_file)
    [p,n,e] = fileparts(seq_wf_file); wf_out = fullfile(p, [n '_stitched' e]);
end
sp_waveforms = sp_waveforms_stitched;
save(wf_out, 'sp_waveforms', 'blank_end_ms', 'FS', 'Stitch_params', '-v7.3');
fprintf('\nSaved stitched spike waveforms to %s\n', wf_out);

times_out = strrep(wf_out, '_waveforms_stitched.mat', '_stitched.mat');
sp_clipped = sp_clipped_stitched;
save(times_out, 'sp_clipped', 'blank_end_ms', 'FS', 'Stitch_params', '-v7.3');
fprintf('Saved stitched spike times to %s\n', times_out);

%% ================= LOCAL HELPER =================
function c = ifelse_col1(x)
    if isempty(x)
        c = x;
    else
        c = x(:,1);
    end
end