%% ========================================================================
% BAD-TRIAL REFERENCE CRITERIA (PRINT ONLY)
%
% PURPOSE
%   Load a saved <base_name>_TrialSpikeCounts.mat file (from
%   TrialSpikeCounts_v2.m) and flag trials whose baseline-corrected
%   analysis-window spike count total is < lowThreshold OR > highThreshold,
%   then print the flagged trials - grouped by condition (Set x Amp x PTD)
%   - using the RELATIVE trial ID (ConditionTrialIndex, i.e. the trial's
%   position within that specific condition).
%
%   This is a GLOBAL bad-trial reference: a trial is flagged for the whole
%   dataset, not per channel. Nothing is saved - this is reference only,
%   printed to the command window, for you to use when running
%   BadTrial_ManualConfirm afterwards.
%
% METRIC
%   Uses the analysis window baked into the saved file at the time
%   TrialSpikeCounts_v2.m was run (TotalBaselineCorrectedAnalysisWindowSpikes
%   / Channel_Results.BaselineCorrected_AnalysisWindow_Spikes). To use a
%   different window, re-run TrialSpikeCounts_v2.m with a different
%   analysis_win_ms first.
% ========================================================================

clear; clc;

%% ============================ USER SETTINGS ===========================

% Full path to the saved <base_name>_TrialSpikeCounts.mat file
matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Linearity_S4_100_400um_SimSeq1_TrialSpikeCounts.mat';

% Flag a trial if its (baseline-corrected, analysis-window) spike count
% total is < lowThreshold OR > highThreshold
lowThreshold  = 0;
highThreshold = 200;

% Which channels to sum for the per-trial total - double choice:
%   Leave empty [] to use the SAME channels the file was originally saved
%   with (count_channels), and reuse the already-computed TrialSummary
%   totals directly (fastest).
%   Or specify a custom channel list, e.g. InspectChannels = 1:32, to
%   recompute totals over a different channel subset - this works because
%   Channel_Results (saved per trial) covers every map-position channel,
%   not just the channels originally used for totals.
InspectChannels = [];

%% =========================== LOAD THE FILE ============================

if ~isfile(matFilePath)
    error('File does not exist:\n%s',matFilePath);
end

fprintf('\n============================================================\n');
fprintf('BAD-TRIAL REFERENCE CRITERIA\n');
fprintf('============================================================\n');
fprintf('File: %s\n',matFilePath);

loadedData = load(matFilePath,'SpikeCounts','TrialSummary','count_channels', ...
    'selected_channels','baseline_win_ms','analysis_win_ms','Electrode_Type','sim_stim');

SpikeCounts       = loadedData.SpikeCounts;
TrialSummary      = loadedData.TrialSummary;
saved_count_chns  = loadedData.count_channels;
selected_channels = loadedData.selected_channels;
baseline_win_ms   = loadedData.baseline_win_ms;
analysis_win_ms   = loadedData.analysis_win_ms;

nTrials = height(TrialSummary);

fprintf('Baseline window used in this file: [%g,%g) ms\n',baseline_win_ms);
fprintf('Analysis window used in this file: [%g,%g) ms\n',analysis_win_ms);
fprintf('Trials in file: %d\n',nTrials);
fprintf('Threshold: flag if value < %g OR > %g\n',lowThreshold,highThreshold);

%% ===================== DETERMINE CHANNELS TO USE ======================

useCustomChannels = ~isempty(InspectChannels);

if useCustomChannels

    inspectChannels = unique(double(InspectChannels(:).'),'stable');

    validChannels = isfinite(inspectChannels) & ...
        inspectChannels >= 1 & ...
        inspectChannels <= numel(selected_channels) & ...
        fix(inspectChannels) == inspectChannels;

    if any(~validChannels)
        warning('Ignoring invalid InspectChannels: %s', ...
            num2str(inspectChannels(~validChannels)));
        inspectChannels = inspectChannels(validChannels);
    end

    if isempty(inspectChannels)
        error('No valid InspectChannels remain.');
    end

    fprintf('Using CUSTOM channel subset (%d channels): %s\n', ...
        numel(inspectChannels),num2str(inspectChannels));
else
    inspectChannels = saved_count_chns;
    fprintf('Using the ORIGINAL saved count_channels (%d channels): %s\n', ...
        numel(inspectChannels),num2str(inspectChannels));
end

%% ===================== COMPUTE PER-TRIAL METRIC =======================

metricValue = nan(nTrials,1);

if ~useCustomChannels
    % Fast path: reuse the already-computed trial totals directly.
    metricValue = TrialSummary.TotalBaselineCorrectedAnalysisWindowSpikes;
else
    % Recompute totals over the custom channel subset, using the
    % per-channel baseline-corrected analysis-window counts already
    % saved in each trial's Channel_Results table.
    for trial_id = 1:nTrials

        si  = TrialSummary.SetIndex(trial_id);
        ai  = TrialSummary.AmplitudeIndex(trial_id);
        pi  = TrialSummary.PTDIndex(trial_id);
        cti = TrialSummary.ConditionTrialIndex(trial_id);

        T = SpikeCounts.set(si).amp(ai).ptd(pi).trial(cti);
        CR = T.Channel_Results;

        [foundChn,chnPos] = ismember(inspectChannels,CR.Channel_Index);

        if any(~foundChn)
            warning(['Trial %d: some InspectChannels not found in ' ...
                'Channel_Results, skipping those for this trial.'],trial_id);
        end

        chnPos = chnPos(foundChn);
        metricValue(trial_id) = sum(CR.BaselineCorrected_AnalysisWindow_Spikes(chnPos));
    end
end

%% ===================== FLAG TRIALS BY THRESHOLD ========================

isLow  = metricValue < lowThreshold;
isHigh = metricValue > highThreshold;

fprintf('\nTotal trials flagged LOW  (<%g): %d\n',lowThreshold,sum(isLow));
fprintf('Total trials flagged HIGH (>%g): %d\n',highThreshold,sum(isHigh));

%% ===================== PRINT PER CONDITION ============================

fprintf('\n------------------------------------------------------------\n');
fprintf('Flagged trials by condition (relative trial ID = position\n');
fprintf('within that condition, e.g. 1-30):\n');
fprintf('------------------------------------------------------------\n');

nSets = numel(SpikeCounts.set);

for si = 1:nSets

    set_label = SpikeCounts.set(si).set_name;

    dist_um = SpikeCounts.set(si).distance_um;
    if isnan(dist_um)
        dist_text = sprintf('Dist n/a (%s)',SpikeCounts.set(si).distance_note);
    else
        dist_text = sprintf('Dist %g um',dist_um);
    end

    nAMP = numel(SpikeCounts.set(si).amp);

    for ai = 1:nAMP

        amp_label = SpikeCounts.set(si).amp(ai).amp_label;
        nPTD = numel(SpikeCounts.set(si).amp(ai).ptd);

        for pi = 1:nPTD

            PTD_ms_val = SpikeCounts.set(si).amp(ai).ptd(pi).PTD_ms;

            conditionMask = ...
                TrialSummary.SetIndex == si & ...
                TrialSummary.AmplitudeIndex == ai & ...
                TrialSummary.PTDIndex == pi;

            if ~any(conditionMask)
                continue;
            end

            relativeIDs_low  = TrialSummary.ConditionTrialIndex(conditionMask & isLow);
            relativeIDs_high = TrialSummary.ConditionTrialIndex(conditionMask & isHigh);

            if isempty(relativeIDs_low) && isempty(relativeIDs_high)
                continue;
            end

            fprintf('\nSet %d (%s) | %s | PTD %g ms | %s:\n', ...
                si,set_label,amp_label,PTD_ms_val,dist_text);

            if ~isempty(relativeIDs_low)
                fprintf('  LOW  (<%g):  trials %s\n', ...
                    lowThreshold,num2str(sort(relativeIDs_low(:).')));
            end

            if ~isempty(relativeIDs_high)
                fprintf('  HIGH (>%g):  trials %s\n', ...
                    highThreshold,num2str(sort(relativeIDs_high(:).')));
            end
        end
    end
end

fprintf('\n============================================================\n');
fprintf('INSPECTION COMPLETE (no files modified or saved)\n');
fprintf('============================================================\n');