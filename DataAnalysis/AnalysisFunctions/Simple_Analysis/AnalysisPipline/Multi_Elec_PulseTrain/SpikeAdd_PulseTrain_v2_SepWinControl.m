%% ============================================================
%  Recursive Spike Recovery for Pulse-Train / Event-Level Stimulation
%
%  Recovery-end logic:
%    1. Try to use:
%         recovery_end_ms = firstSpikeTimes{ch}(targetTrial)
%    2. If firstSpikeTimes is NaN, use:
%         recovery_end_ms = recoveryFallbackEnd_ms_trial(targetTrial)
%    3. Always apply:
%         recovery_end_ms <= recoveryMaxEnd_ms_trial(targetTrial)
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER INPUT ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Single3';

FS = 30000;

% Include auto-simultaneous conditions in recovery?
include_auto_sim = true;

% Zero-current trials are not recovered.
include_zero_control = false;

% Print detailed recovery summary.
debug_print_recovery_plan = true;

%% ====================== LOAD FOLDER ======================

if ~isfolder(data_folder)
    error('Data folder does not exist.');
end

cd(data_folder);
fprintf('\nData folder:\n%s\n', data_folder);

%% ====================== BASE NAME ======================

parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');

if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ====================== LOAD SPIKES ======================
% Use the blanked / clipped spike file as the starting point.

sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No *sp_xia.mat file found.');

load(sp_file(1).name, 'sp_clipped');
fprintf('Loaded spike file: %s\n', sp_file(1).name);

% Main output variable name.
% sp_seq will be progressively modified.
sp_seq = sp_clipped;

%% ====================== LOAD FIRST SPIKE TIMES AND WINDOW SETTINGS ======================

first_file = dir('FirstSpikeTimes_Prefix.mat');
assert(~isempty(first_file), ...
    'FirstSpikeTimes_Prefix.mat not found. Run the first-spike-time code first.');

F = load(first_file(1).name);

% ----- Required original variables -----
if ~isfield(F, 'firstSpikeTimes')
    error('First-spike file does not contain firstSpikeTimes.');
end

if ~isfield(F, 'hasSpike')
    error('First-spike file does not contain hasSpike.');
end

firstSpikeTimes = F.firstSpikeTimes;
hasSpike        = F.hasSpike;

% ----- Required recovery-window variables -----
if ~isfield(F, 'recoveryFallbackEnd_ms_trial')
    error(['First-spike file does not contain recoveryFallbackEnd_ms_trial. ', ...
           'Please rerun the updated first-spike-time detection code.']);
end

if ~isfield(F, 'recoveryMaxEnd_ms_trial')
    error(['First-spike file does not contain recoveryMaxEnd_ms_trial. ', ...
           'Please rerun the updated first-spike-time detection code.']);
end

recoveryFallbackEnd_ms_trial = F.recoveryFallbackEnd_ms_trial;
recoveryMaxEnd_ms_trial      = F.recoveryMaxEnd_ms_trial;

% ----- Saved parameter values for record keeping -----
% These are not strictly required for the recovery loop because the trial-level
% windows above already contain the actual fallback and cap times.
% However, saving/loading them is useful for checking what settings were used.
fallback_after_last_event_seq_ms     = F.fallback_after_last_event_seq_ms;
fallback_after_last_event_autosim_ms = F.fallback_after_last_event_autosim_ms;

max_after_last_event_seq_ms     = F.max_after_last_event_seq_ms;
max_after_last_event_autosim_ms = F.max_after_last_event_autosim_ms;

fprintf('Loaded first-spike file: %s\n', first_file(1).name);

fprintf('\nRecovery-window settings loaded from first-spike file:\n');
fprintf('  Fallback Seq:     last event + %.2f ms\n', fallback_after_last_event_seq_ms);
fprintf('  Fallback AutoSim: last event + %.2f ms\n', fallback_after_last_event_autosim_ms);
fprintf('  Max cap Seq:      last event + %.2f ms\n', max_after_last_event_seq_ms);
fprintf('  Max cap AutoSim:  last event + %.2f ms\n', max_after_last_event_autosim_ms);

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
    error('StimMeta was not found. This recovery code requires StimMeta.');
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== BASIC CONSISTENCY CHECKS ======================
% These checks help avoid accidentally using a FirstSpikeTimes_Prefix.mat
% file from another folder or another run.

if isfield(F, 'n_Trials')
    if F.n_Trials ~= n_Trials
        warning(['First-spike file n_Trials (%d) does not match current exp file n_Trials (%d). ', ...
                 'Please check that the first-spike file belongs to this folder.'], ...
                 F.n_Trials, n_Trials);
    end
end

if numel(recoveryFallbackEnd_ms_trial) ~= n_Trials
    warning('recoveryFallbackEnd_ms_trial length does not match n_Trials.');
end

if numel(recoveryMaxEnd_ms_trial) ~= n_Trials
    warning('recoveryMaxEnd_ms_trial length does not match n_Trials.');
end

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

% Keep the familiar names.
lastActivePTD_ms = eventEnd_ms_trial;
lastActivePTD_us = eventEnd_us_trial;

%% ====================== CHECK FIRST-SPIKE FILE AGAINST CURRENT METADATA ======================
% The first-spike file also contains lastActivePTD_ms. Compare it against
% the current metadata-derived lastActivePTD_ms. If these differ, you may be
% loading a first-spike file that does not belong to this data folder.

if isfield(F, 'lastActivePTD_ms')

    if numel(F.lastActivePTD_ms) == numel(lastActivePTD_ms)

        diff_last_event = abs(F.lastActivePTD_ms(:) - lastActivePTD_ms(:));
        diff_last_event = diff_last_event(~isnan(diff_last_event));

        if ~isempty(diff_last_event)
            max_diff_last_event = max(diff_last_event);

            if max_diff_last_event > 1e-6
                warning(['lastActivePTD_ms differs between FirstSpikeTimes_Prefix.mat and current StimMeta. ', ...
                         'Max difference = %.6f ms. Check folder/file consistency.'], ...
                         max_diff_last_event);
            end
        end
    else
        warning('F.lastActivePTD_ms length does not match current lastActivePTD_ms length.');
    end
end

%% ====================== RECOVERY GROUPS ======================

all_sets = unique(stimSet_trial(~isnan(stimSet_trial)));
all_sets = all_sets(all_sets > 0);

all_amps = unique(trialAmps(~isnan(trialAmps)));
all_amps = all_amps(all_amps > 0);

all_levels = unique(trainLevel_trial(~isnan(trainLevel_trial)));
all_levels = sort(all_levels(:)');

fprintf('\nDetected sets: ');
disp(all_sets');

fprintf('Detected amplitudes: ');
disp(all_amps');

fprintf('Detected train levels: ');
disp(all_levels);

%% ====================== RECOVERY INFO ======================

nChn = numel(sp_seq);

% usedFallback:
%   1 if firstSpikeTimes{ch}(targetTrial) was NaN and the trial-level
%   fallback end from FirstSpikeTimes_Prefix.mat was used.
usedFallback = zeros(nChn, n_Trials);

% usedCap:
%   1 if recovery_end_ms was shortened by recoveryMaxEnd_ms_trial.
usedCap = zeros(nChn, n_Trials);

% Actual recovery end used for each channel/trial.
recoveryEnd_ms_used = nan(nChn, n_Trials);

% Trial-level cap used for each channel/trial.
% This records the cap that was available, even if the cap was not applied.
recoveryEndCap_ms_used = nan(nChn, n_Trials);

nAddedSpikes = zeros(nChn, n_Trials);
nDeletedSpikes = zeros(nChn, n_Trials);

fprintf('\nStarting recursive pulse-train spike recovery...\n');

%% ====================== MAIN RECOVERY LOOP ======================
% Important:
%   We process levels in ascending order.
%   This means when Level 3 is being recovered, Level 2 has already been
%   recovered and is available inside sp_seq.
%
% Recovery chain:
%   Level 1 -> Level 2 -> Level 3 -> ...

for li = 2:numel(all_levels)

    targetLevel = all_levels(li);
    sourceLevel = all_levels(li-1);

    fprintf('\nRecovering Level %.0f using recovered Level %.0f...\n', ...
        targetLevel, sourceLevel);

    for si = 1:numel(all_sets)

        setID = all_sets(si);

        for ai = 1:numel(all_amps)

            ampVal = all_amps(ai);

            % We recover auto-sim and non-auto-sim separately.
            for autoFlag = [0 1]

                if autoFlag == 1 && ~include_auto_sim
                    continue;
                end

                if autoFlag == 1
                    familyLabel_this = 'AutoSim';
                else
                    familyLabel_this = 'Seq';
                end

                % Source trials from previous level.
                sourceTrials = find(stimSet_trial == setID & ...
                                    trialAmps == ampVal & ...
                                    trainLevel_trial == sourceLevel & ...
                                    isAutoSim_trial == logical(autoFlag));

                % Target trials from current level.
                targetTrials = find(stimSet_trial == setID & ...
                                    trialAmps == ampVal & ...
                                    trainLevel_trial == targetLevel & ...
                                    isAutoSim_trial == logical(autoFlag));

                if ~include_zero_control
                    sourceTrials = sourceTrials(isZeroControl(sourceTrials) == 0);
                    targetTrials = targetTrials(isZeroControl(targetTrials) == 0);
                end

                sourceTrials = sourceTrials(sourceTrials <= nTrials_use);
                targetTrials = targetTrials(targetTrials <= nTrials_use);

                if isempty(sourceTrials) || isempty(targetTrials)
                    continue;
                end

                sourceTrials = sort(sourceTrials(:));
                targetTrials = sort(targetTrials(:));

                if debug_print_recovery_plan
                    fprintf('  Set %g | Amp %.1f | %s | Source N=%d | Target N=%d\n', ...
                        setID, ampVal, familyLabel_this, ...
                        numel(sourceTrials), numel(targetTrials));
                end

                %% ----- Target trial loop -----
                for tt = 1:numel(targetTrials)

                    targetTrial = targetTrials(tt);

                    % Match source trial by repeat index.
                    % If the number of source and target trials differs,
                    % cycle through source trials.
                    sourceTrial = sourceTrials(mod(tt-1, numel(sourceTrials)) + 1);

                    target_t0 = trig(targetTrial) / FS * 1000;
                    source_t0 = trig(sourceTrial) / FS * 1000;

                    %% ----- Channel loop -----
                    for ch = 1:nChn

                        if isempty(sp_seq{ch})
                            continue;
                        end

                        % ====================================================
                        % Recovery end is channel/trial specific.
                        %
                        % Step 1:
                        %   Try firstSpikeTimes{ch}(targetTrial).
                        %
                        % Step 2:
                        %   If NaN, use recoveryFallbackEnd_ms_trial(targetTrial).
                        %   This value was saved by the first-spike-time code.
                        %
                        % Step 3:
                        %   Always cap using recoveryMaxEnd_ms_trial(targetTrial).
                        %   This value was also saved by the first-spike-time code.
                        %
                        % This makes the recovery code use exactly the same
                        % window settings that were saved with the first-spike
                        % detection file.
                        % ====================================================

                        recovery_end_ms = firstSpikeTimes{ch}(targetTrial);

                        fallback_end_ms = recoveryFallbackEnd_ms_trial(targetTrial);
                        cap_end_ms      = recoveryMaxEnd_ms_trial(targetTrial);

                        recoveryEndCap_ms_used(ch,targetTrial) = cap_end_ms;

                        if isnan(recovery_end_ms)

                            % Fallback when no first spike was detected.
                            recovery_end_ms = fallback_end_ms;
                            usedFallback(ch,targetTrial) = 1;
                        end

                        if isnan(recovery_end_ms) || recovery_end_ms <= 0
                            continue;
                        end

                        % Apply the saved trial-level maximum cap.
                        if ~isnan(cap_end_ms) && recovery_end_ms > cap_end_ms
                            recovery_end_ms = cap_end_ms;
                            usedCap(ch,targetTrial) = 1;
                        end

                        if isnan(recovery_end_ms) || recovery_end_ms <= 0
                            continue;
                        end

                        recoveryEnd_ms_used(ch,targetTrial) = recovery_end_ms;

                        % Recovery window is from 0 to recovery_end_ms.
                        target_abs_start = target_t0;
                        target_abs_end   = target_t0 + recovery_end_ms;

                        source_abs_start = source_t0;
                        source_abs_end   = source_t0 + recovery_end_ms;

                        %% ----- Delete target spikes in recovery window -----
                        targetRowsToDelete = sp_seq{ch}(:,1) >= target_abs_start & ...
                                             sp_seq{ch}(:,1) <= target_abs_end;

                        nDeletedSpikes(ch,targetTrial) = nDeletedSpikes(ch,targetTrial) + ...
                            sum(targetRowsToDelete);

                        sp_seq{ch}(targetRowsToDelete,:) = [];

                        %% ----- Copy source spikes from recovered source level -----
                        % Because sp_seq has already been updated level by level,
                        % sourceTrial contains recovered spikes if sourceLevel
                        % has already been recovered.
                        sourceRowsToCopy = sp_seq{ch}(:,1) >= source_abs_start & ...
                                           sp_seq{ch}(:,1) <= source_abs_end;

                        sourceRows = sp_seq{ch}(sourceRowsToCopy,:);

                        if isempty(sourceRows)
                            continue;
                        end

                        % Convert source spike times to relative times.
                        sourceRelTimes = sourceRows(:,1) - source_t0;

                        % Move source spikes to target trial time.
                        sourceRows(:,1) = target_t0 + sourceRelTimes;

                        %% ----- Add recovered spikes to target trial -----
                        sp_seq{ch} = [sp_seq{ch}; sourceRows];

                        nAddedSpikes(ch,targetTrial) = nAddedSpikes(ch,targetTrial) + ...
                            size(sourceRows,1);

                        % Sort by spike time.
                        [~, sortIdx] = sort(sp_seq{ch}(:,1));
                        sp_seq{ch} = sp_seq{ch}(sortIdx,:);
                    end
                end
            end
        end
    end
end

fprintf('\nRecursive recovery finished.\n');

fprintf('Total added spikes: %d\n', sum(nAddedSpikes(:)));
fprintf('Total deleted spikes before recovery: %d\n', sum(nDeletedSpikes(:)));
fprintf('Total fallback windows used: %d\n', sum(usedFallback(:)));
fprintf('Total capped recovery windows: %d\n', sum(usedCap(:)));

%% ====================== SAVE ======================

save_name = [base_name '.sp_xia_PrefixRecovery.mat'];

save(save_name, ...
     'sp_seq', ...
     'usedFallback', ...
     'usedCap', ...
     'recoveryEnd_ms_used', ...
     'recoveryEndCap_ms_used', ...
     'nAddedSpikes', ...
     'nDeletedSpikes', ...
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
     'isAutoSim_trial', ...
     'isZeroControl', ...
     'conditionID_trial', ...
     'recoveryFallbackEnd_ms_trial', ...
     'recoveryMaxEnd_ms_trial', ...
     'fallback_after_last_event_seq_ms', ...
     'fallback_after_last_event_autosim_ms', ...
     'max_after_last_event_seq_ms', ...
     'max_after_last_event_autosim_ms', ...
     'FS', ...
     'n_Trials', ...
     'nTrials_use', ...
     '-v7.3');

fprintf('Saved recovered spike file to:\n%s\n', fullfile(data_folder, save_name));