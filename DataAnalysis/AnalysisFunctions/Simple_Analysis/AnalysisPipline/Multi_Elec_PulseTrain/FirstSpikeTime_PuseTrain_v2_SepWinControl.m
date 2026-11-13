%% ============================================================
%  Extract First Post-Blank Spike Time per Trial
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER INPUT ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

FS = 30000;

%% ====================== WINDOW SETTINGS ======================

% Sequential first-spike search window, relative to last event.
first_spike_search_start_after_last_event_seq_ms = 0;
first_spike_search_end_after_last_event_seq_ms   = 10;

% AutoSim / simultaneous first-spike search window, relative to last event.
first_spike_search_start_after_last_event_autosim_ms = 0;
first_spike_search_end_after_last_event_autosim_ms   = 3;

%% ---------- Recovery fallback window ----------
% The recovery code should use this only when:
%   firstSpikeTimes{ch}(trial) is NaN
% Then:
%   recovery_end = last event time + fallback_after_last_event_*_ms

fallback_after_last_event_seq_ms     = 4;
fallback_after_last_event_autosim_ms = 0;

%% ---------- Recovery maximum cap ----------

max_after_last_event_seq_ms     = 10;
max_after_last_event_autosim_ms = 3;

%% ---------- Trial inclusion ----------
% Include auto simultaneous conditions in first-spike detection?
include_auto_sim = true;

% Zero-current trials usually should not be searched.
include_zero_control = false;

% If true, print example trials to check metadata.
debug_print_trial_content = true;

%% ====================== CHECK WINDOW SETTINGS ======================

if first_spike_search_end_after_last_event_seq_ms <= first_spike_search_start_after_last_event_seq_ms
    error('Sequential first-spike search end must be larger than start.');
end

if first_spike_search_end_after_last_event_autosim_ms <= first_spike_search_start_after_last_event_autosim_ms
    error('AutoSim first-spike search end must be larger than start.');
end

if fallback_after_last_event_seq_ms < 0 || fallback_after_last_event_autosim_ms < 0
    error('Fallback windows must be non-negative.');
end

if max_after_last_event_seq_ms < 0 || max_after_last_event_autosim_ms < 0
    error('Max-after-last-event windows must be non-negative.');
end

if max_after_last_event_seq_ms < fallback_after_last_event_seq_ms
    warning(['Seq max_after_last_event is smaller than fallback_after_last_event. ', ...
             'If firstSpikeTimes is NaN, the recovery code will probably cap the fallback window.']);
end

if max_after_last_event_autosim_ms < fallback_after_last_event_autosim_ms
    warning(['AutoSim max_after_last_event is smaller than fallback_after_last_event. ', ...
             'If firstSpikeTimes is NaN, the recovery code will probably cap the fallback window.']);
end

% Legacy variable kept for compatibility / record keeping.
% The old code used one shared win_ms. The actual detection now uses the
% separate Seq/AutoSim search windows above.
win_ms = max([ ...
    first_spike_search_end_after_last_event_seq_ms - first_spike_search_start_after_last_event_seq_ms, ...
    first_spike_search_end_after_last_event_autosim_ms - first_spike_search_start_after_last_event_autosim_ms]);

%% ====================== LOAD FOLDER ======================

if ~isfolder(data_folder)
    error('Data folder does not exist.');
end

cd(data_folder);
fprintf('\nData folder:\n%s\n', data_folder);

%% ====================== LOAD SPIKES ======================
% This code expects the spike file generated after blanking / spike extraction.
% It should contain sp_clipped.

sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No *sp_xia.mat file found.');

load(sp_file(1).name, 'sp_clipped');
fprintf('Loaded spike file: %s\n', sp_file(1).name);

%% ====================== LOAD TRIGGERS ======================

if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD EXPERIMENT PARAMETERS ======================

param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat file found.');

S = load(param_file(1).name);

StimParams        = S.StimParams;
TrialParams       = S.TrialParams;
simultaneous_stim = S.simultaneous_stim;
n_Trials          = S.n_Trials;

if isfield(S, 'StimMeta')
    StimMeta = S.StimMeta;
else
    error('StimMeta was not found. This code requires the pulse-train StimMeta structure.');
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== REMOVE HEADER ROWS ======================

StimParams_data  = StimParams(2:end,:);
TrialParams_data = TrialParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

if size(TrialParams_data,1) ~= expected_rows
    warning('TrialParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(TrialParams_data,1), expected_rows);
end

%% ====================== TRIAL CONDITION IDS ======================
% TrialParams columns:
%   column 1 = trial number
%   column 2 = condition ID
%   column 3 = internal electrode ID
%
% One trial has simultaneous_stim rows.
% Therefore we use the first row of each trial to get the condition ID.

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));

conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of trial-level condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA FROM StimMeta ======================

stimSet_trial        = NaN(n_Trials,1);
trainLevel_trial     = NaN(n_Trials,1);
totalLevels_trial    = NaN(n_Trials,1);
trialAmps            = NaN(n_Trials,1);
isAutoSim_trial      = false(n_Trials,1);
isZeroControl        = false(n_Trials,1);
eventEnd_ms_trial    = NaN(n_Trials,1);
eventEnd_us_trial    = NaN(n_Trials,1);

eventTimes_ms_trial  = cell(n_Trials,1);
pulseCount_trial     = cell(n_Trials,1);

for tr = 1:n_Trials

    condID = conditionID_trial(tr);

    if condID < 1 || condID > numel(StimMeta)
        warning('Trial %d has invalid condition ID %d.', tr, condID);
        continue;
    end

    meta = StimMeta(condID);

    if isfield(meta, 'StimSetIndex')
        stimSet_trial(tr) = meta.StimSetIndex;
    end

    if isfield(meta, 'TrainLevel')
        trainLevel_trial(tr) = meta.TrainLevel;
    end

    if isfield(meta, 'TotalTrainLevels')
        totalLevels_trial(tr) = meta.TotalTrainLevels;
    end

    if isfield(meta, 'Amp_uA')
        amp_vec = meta.Amp_uA;
        amp_vec = amp_vec(:)';

        if all(amp_vec <= 0)
            trialAmps(tr) = 0;
        else
            trialAmps(tr) = max(amp_vec(amp_vec > 0));
        end
    end

    if isfield(meta, 'IsAutoSimultaneous')
        isAutoSim_trial(tr) = logical(meta.IsAutoSimultaneous);
    end

    if isfield(meta, 'IsZeroCurrentControl')
        isZeroControl(tr) = logical(meta.IsZeroCurrentControl);
    end

    if isfield(meta, 'EventTimesIncluded_ms')
        eventTimes_ms_trial{tr} = meta.EventTimesIncluded_ms;
    else
        eventTimes_ms_trial{tr} = [];
    end

    if isfield(meta, 'EventEndTime_ms')
        eventEnd_ms_trial(tr) = meta.EventEndTime_ms;
        eventEnd_us_trial(tr) = meta.EventEndTime_ms * 1000;
    elseif ~isempty(eventTimes_ms_trial{tr})
        eventEnd_ms_trial(tr) = max(eventTimes_ms_trial{tr});
        eventEnd_us_trial(tr) = eventEnd_ms_trial(tr) * 1000;
    else
        eventEnd_ms_trial(tr) = 0;
        eventEnd_us_trial(tr) = 0;
    end

    if isfield(meta, 'PulseCountPerElectrode')
        pulseCount_trial{tr} = meta.PulseCountPerElectrode;
    else
        pulseCount_trial{tr} = [];
    end
end

% Keep these names because later recovery code expects the old names.
lastActivePTD_ms = eventEnd_ms_trial;
lastActivePTD_us = eventEnd_us_trial;

%% ====================== SELECT TRIALS FOR FIRST-SPIKE DETECTION ======================

% Start with all non-zero pulse-train trials.
isPulseTrainTrial = ~isZeroControl;

% Optionally exclude AutoSim trials.
if ~include_auto_sim
    isPulseTrainTrial = isPulseTrainTrial & ~isAutoSim_trial;
end

% Zero-current control handling.
if include_zero_control
    isPulseTrainTrial = isPulseTrainTrial | isZeroControl;
else
    isPulseTrainTrial = isPulseTrainTrial & ~isZeroControl;
end

% Only use trials that have triggers.
if nTrials_use < n_Trials
    isPulseTrainTrial((nTrials_use+1):end) = false;
end

%% ====================== BUILD TRIAL-LEVEL SEARCH AND RECOVERY WINDOWS ======================
% These arrays record the actual timing windows for each trial.
%
% firstSpikeSearchStart_ms_trial / firstSpikeSearchEnd_ms_trial:
%   Used in this script to find firstSpikeTimes.
%
% recoveryFallbackEnd_ms_trial:
%   Saved for recovery code.
%   If firstSpikeTimes{ch}(trial) is NaN:
%       recovery_end = recoveryFallbackEnd_ms_trial(trial)
%
% recoveryMaxEnd_ms_trial:
%   Saved for recovery code.
%   Recovery code should always apply:
%       recovery_end = min(recovery_end, recoveryMaxEnd_ms_trial(trial))

firstSpikeSearchStart_ms_trial = NaN(n_Trials,1);
firstSpikeSearchEnd_ms_trial   = NaN(n_Trials,1);

recoveryFallbackEnd_ms_trial   = NaN(n_Trials,1);
recoveryMaxEnd_ms_trial        = NaN(n_Trials,1);

for tr = 1:nTrials_use

    if ~isPulseTrainTrial(tr)
        continue;
    end

    if isnan(lastActivePTD_ms(tr))
        continue;
    end

    if isAutoSim_trial(tr)

        % AutoSim / simultaneous condition:
        % use AutoSim-specific settings.
        firstSpikeSearchStart_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            first_spike_search_start_after_last_event_autosim_ms;

        firstSpikeSearchEnd_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            first_spike_search_end_after_last_event_autosim_ms;

        recoveryFallbackEnd_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            fallback_after_last_event_autosim_ms;

        recoveryMaxEnd_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            max_after_last_event_autosim_ms;

    else

        % Sequential condition:
        % use Seq-specific settings.
        firstSpikeSearchStart_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            first_spike_search_start_after_last_event_seq_ms;

        firstSpikeSearchEnd_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            first_spike_search_end_after_last_event_seq_ms;

        recoveryFallbackEnd_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            fallback_after_last_event_seq_ms;

        recoveryMaxEnd_ms_trial(tr) = lastActivePTD_ms(tr) + ...
            max_after_last_event_seq_ms;
    end
end

%% ====================== PRINT DETECTED CONDITIONS ======================

fprintf('\nDetected StimSetIndex values: ');
disp(unique(stimSet_trial(~isnan(stimSet_trial)))');

fprintf('Detected train levels: ');
disp(unique(trainLevel_trial(~isnan(trainLevel_trial)))');

fprintf('Detected amplitudes: ');
disp(unique(trialAmps(~isnan(trialAmps)))');

fprintf('Detected final event times / lastActivePTD_ms: ');
disp(unique(lastActivePTD_ms(~isnan(lastActivePTD_ms)))');

fprintf('Number of trials used for first-spike search: %d\n', ...
    sum(isPulseTrainTrial(1:nTrials_use)));

fprintf('\n================ WINDOW SETTINGS SAVED TO FIRST-SPIKE FILE ================\n');

fprintf('\nFirst-spike search window settings:\n');
fprintf('  Seq:     last event + %.2f ms to last event + %.2f ms\n', ...
    first_spike_search_start_after_last_event_seq_ms, ...
    first_spike_search_end_after_last_event_seq_ms);

fprintf('  AutoSim: last event + %.2f ms to last event + %.2f ms\n', ...
    first_spike_search_start_after_last_event_autosim_ms, ...
    first_spike_search_end_after_last_event_autosim_ms);

fprintf('\nRecovery fallback settings:\n');
fprintf('  Seq:     recovery fallback end = last event + %.2f ms\n', ...
    fallback_after_last_event_seq_ms);
fprintf('  AutoSim: recovery fallback end = last event + %.2f ms\n', ...
    fallback_after_last_event_autosim_ms);

fprintf('\nRecovery max-cap settings:\n');
fprintf('  Seq:     recovery max end = last event + %.2f ms\n', ...
    max_after_last_event_seq_ms);
fprintf('  AutoSim: recovery max end = last event + %.2f ms\n', ...
    max_after_last_event_autosim_ms);

fprintf('\nActual trial-level first-spike search start times:\n');
disp(unique(firstSpikeSearchStart_ms_trial(~isnan(firstSpikeSearchStart_ms_trial)))');

fprintf('Actual trial-level first-spike search end times:\n');
disp(unique(firstSpikeSearchEnd_ms_trial(~isnan(firstSpikeSearchEnd_ms_trial)))');

fprintf('Actual trial-level recovery fallback end times:\n');
disp(unique(recoveryFallbackEnd_ms_trial(~isnan(recoveryFallbackEnd_ms_trial)))');

fprintf('Actual trial-level recovery max end times:\n');
disp(unique(recoveryMaxEnd_ms_trial(~isnan(recoveryMaxEnd_ms_trial)))');

fprintf('==========================================================================\n');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================

if debug_print_trial_content

    fprintf('\n================ DEBUG EVENT-TRAIN FIRST-SPIKE CHECK ================\n');

    set_all = unique(stimSet_trial(~isnan(stimSet_trial)));
    amp_all = unique(trialAmps(~isnan(trialAmps) & trialAmps > 0));
    level_all = unique(trainLevel_trial(~isnan(trainLevel_trial)));

    if ~isempty(set_all) && ~isempty(amp_all) && ~isempty(level_all)

        debug_set = set_all(1);
        debug_amp = amp_all(end);

        for li = 1:numel(level_all)

            level_val = level_all(li);

            tr_debug = find(stimSet_trial == debug_set & ...
                            trialAmps == debug_amp & ...
                            trainLevel_trial == level_val & ...
                            isZeroControl == 0, ...
                            1, 'first');

            if isempty(tr_debug)
                continue;
            end

            rr = (tr_debug-1)*simultaneous_stim + (1:simultaneous_stim);

            stimNames_debug   = StimParams_data(rr,1);
            ptd_debug         = cell2mat(StimParams_data(rr,6));
            pulseNum_debug    = cell2mat(StimParams_data(rr,8));
            pulsePeriod_debug = cell2mat(StimParams_data(rr,9));
            amp_debug         = cell2mat(StimParams_data(rr,16));

            fprintf('\nSet %g | Amp %.1f uA | TrainLevel %g | Trial %d | CondID %d\n', ...
                debug_set, debug_amp, level_val, tr_debug, conditionID_trial(tr_debug));

            fprintf('  IsAutoSimultaneous:      %d\n', isAutoSim_trial(tr_debug));
            fprintf('  StimNames:               ');
            disp(stimNames_debug');

            fprintf('  PTD_us:                  ');
            disp(ptd_debug');

            fprintf('  PulseNum:                ');
            disp(pulseNum_debug');

            fprintf('  PulsePeriod_us:          ');
            disp(pulsePeriod_debug');

            fprintf('  Amp_col16:               ');
            disp(amp_debug');

            fprintf('  EventTimesIncluded_ms:   ');
            disp(eventTimes_ms_trial{tr_debug});

            fprintf('  lastActivePTD_ms:        %.3f\n', lastActivePTD_ms(tr_debug));

            fprintf('  First-spike search win:  %.3f to %.3f ms\n', ...
                firstSpikeSearchStart_ms_trial(tr_debug), ...
                firstSpikeSearchEnd_ms_trial(tr_debug));

            fprintf('  Recovery fallback end:   %.3f ms\n', ...
                recoveryFallbackEnd_ms_trial(tr_debug));

            fprintf('  Recovery max end:        %.3f ms\n', ...
                recoveryMaxEnd_ms_trial(tr_debug));
        end
    end

    fprintf('=====================================================================\n');
end

%% ====================== INITIALIZE OUTPUTS ======================

nChn = numel(sp_clipped);

% Keep the saved variable structure:
%   firstSpikeTimes{ch}(trial)
%   hasSpike(ch,trial)
firstSpikeTimes = cell(nChn,1);
hasSpike = zeros(nChn, n_Trials);

fprintf('\nExtracting first post-blank spike times...\n');

%% ====================== MAIN LOOP ======================

for ch = 1:nChn

    S_ch = sp_clipped{ch};
    firstSpikeTimes{ch} = nan(n_Trials,1);

    if isempty(S_ch)
        continue;
    end

    spike_times_ms = S_ch(:,1);

    for tr = 1:nTrials_use

        if ~isPulseTrainTrial(tr)
            continue;
        end

        win_start = firstSpikeSearchStart_ms_trial(tr);
        win_end   = firstSpikeSearchEnd_ms_trial(tr);

        if isnan(win_start) || isnan(win_end) || win_end <= win_start
            continue;
        end

        % Trigger time in absolute ms.
        t0 = trig(tr) / FS * 1000;

        % Convert absolute spike times to trial-relative spike times.
        rel_t = spike_times_ms - t0;

        % Search in the family-specific window.
        spk_in_win = rel_t(rel_t >= win_start & rel_t <= win_end);

        if ~isempty(spk_in_win)
            firstSpikeTimes{ch}(tr) = spk_in_win(1);
            hasSpike(ch,tr) = 1;
        end
    end

    fprintf('Ch %2d: first spikes found in %d/%d searched trials.\n', ...
        ch, sum(hasSpike(ch,isPulseTrainTrial)), sum(isPulseTrainTrial(1:nTrials_use)));
end

fprintf('==============================================================\n');

%% ====================== SUMMARY BY FAMILY ======================

searched_seq = isPulseTrainTrial(:) & ~isAutoSim_trial(:);
searched_auto = isPulseTrainTrial(:) & isAutoSim_trial(:);

fprintf('\nFirst-spike detection summary:\n');

if any(searched_seq)
    n_seq_search = sum(searched_seq(1:nTrials_use));
    n_seq_found_per_ch = sum(hasSpike(:,searched_seq), 2);
    fprintf('  Seq searched trials: %d\n', n_seq_search);
    fprintf('  Seq found per channel: mean %.1f / %d trials\n', ...
        mean(n_seq_found_per_ch), n_seq_search);
end

if any(searched_auto)
    n_auto_search = sum(searched_auto(1:nTrials_use));
    n_auto_found_per_ch = sum(hasSpike(:,searched_auto), 2);
    fprintf('  AutoSim searched trials: %d\n', n_auto_search);
    fprintf('  AutoSim found per channel: mean %.1f / %d trials\n', ...
        mean(n_auto_found_per_ch), n_auto_search);
end

%% ====================== SAVE ======================
save_name = 'FirstSpikeTimes_Prefix.mat';

save(save_name, ...
     'firstSpikeTimes', ...
     'hasSpike', ...
     'lastActivePTD_ms', ...
     'lastActivePTD_us', ...
     'eventEnd_ms_trial', ...
     'eventEnd_us_trial', ...
     'eventTimes_ms_trial', ...
     'pulseCount_trial', ...
     'stimSet_trial', ...
     'trainLevel_trial', ...
     'totalLevels_trial', ...
     'trialAmps', ...
     'isPulseTrainTrial', ...
     'isAutoSim_trial', ...
     'isZeroControl', ...
     'conditionID_trial', ...
     'win_ms', ...
     'first_spike_search_start_after_last_event_seq_ms', ...
     'first_spike_search_end_after_last_event_seq_ms', ...
     'first_spike_search_start_after_last_event_autosim_ms', ...
     'first_spike_search_end_after_last_event_autosim_ms', ...
     'fallback_after_last_event_seq_ms', ...
     'fallback_after_last_event_autosim_ms', ...
     'max_after_last_event_seq_ms', ...
     'max_after_last_event_autosim_ms', ...
     'firstSpikeSearchStart_ms_trial', ...
     'firstSpikeSearchEnd_ms_trial', ...
     'recoveryFallbackEnd_ms_trial', ...
     'recoveryMaxEnd_ms_trial', ...
     'FS', ...
     'n_Trials', ...
     'nTrials_use', ...
     '-v7.3');

fprintf('Saved to: %s\n', fullfile(data_folder, save_name));