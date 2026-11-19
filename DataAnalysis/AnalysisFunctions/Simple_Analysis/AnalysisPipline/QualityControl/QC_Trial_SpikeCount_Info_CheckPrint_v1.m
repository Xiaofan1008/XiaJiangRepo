%% ========================================================================
% INSPECT_MULTIISI_TRIAL_SPIKE_COUNTS
%
% PURPOSE
%   Load a saved *_MultiISI_TrialSpikeCounts_PerTrial.mat file and flag
%   trials whose 2-40ms baseline-corrected analysis-window spike count is
%   below a low threshold or above a high threshold (e.g. <0 or >80),
%   then print the flagged trials -- grouped by condition (Set x Amp x
%   PTD) -- using the RELATIVE trial ID (i.e. the trial's position within
%   that specific condition, matching ConditionTrialIndex / the
%   "trials 1-30 for this condition" numbering).
%
%   No files are modified or saved -- this is inspection only, printed to
%   the command window.
%
% METRIC
%   Uses the analysis window baked into the saved file (2-40ms baseline-
%   corrected spike count, i.e. TotalBaselineCorrectedAnalysisWindowSpikes
%   / Channel_Results.BaselineCorrected_AnalysisWindow_Spikes). If you
%   need a different window (e.g. true 0-40ms), you'd need to re-run the
%   original spike-count script with a different analysis_win_ms first --
%   that window is fixed at save time, not something this script can
%   change after the fact.
%
%   Just edit the "USER SETTINGS" block below and press Run.
% ========================================================================

clear; clc;

%% ============================ USER SETTINGS ===========================

% Full path to the saved *_MultiISI_TrialSpikeCounts_PerTrial.mat file
matFilePath = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1/Xia_ISI_SimSeq1_MultiISI_TrialSpikeCounts_PerTrial.mat';

% Flag a trial if its (baseline-corrected, 2-40ms) spike count total is
% < lowThreshold OR > highThreshold
lowThreshold  = 0;
highThreshold = 50;

% Which channels to sum for the per-trial total.
%   Leave empty [] to use the SAME channels the file was originally saved
%   with (count_channels), and reuse the already-computed
%   TrialSummary totals directly (fastest).
%   Or specify a custom channel list, e.g. InspectChannels = 35:64, to
%   recompute totals over a different channel subset -- this works
%   because Channel_Results (saved per trial) contains the baseline-
%   corrected analysis-window count for EVERY recording channel, not just
%   the channels originally used for totals.
InspectChannels = [1:16,33:63];

%% =========================== LOAD THE FILE ============================

if ~isfile(matFilePath)
    error('File does not exist:\n%s', matFilePath);
end

fprintf('\n============================================================\n');
fprintf('INSPECT MULTI-ISI TRIAL SPIKE COUNTS\n');
fprintf('============================================================\n');
fprintf('File: %s\n', matFilePath);

loadedData = load(matFilePath, 'SpikeCounts', 'TrialSummary', ...
    'count_channels', 'selected_channels', 'analysis_win_ms');

SpikeCounts       = loadedData.SpikeCounts;
TrialSummary      = loadedData.TrialSummary;
saved_count_chns  = loadedData.count_channels;
selected_channels = loadedData.selected_channels;
analysis_win_ms   = loadedData.analysis_win_ms;

nTrials = height(TrialSummary);

fprintf('Analysis window used in this file: [%g,%g) ms\n', analysis_win_ms);
fprintf('Trials in file: %d\n', nTrials);
fprintf('Threshold: flag if value < %g OR > %g\n', lowThreshold, highThreshold);

%% ===================== DETERMINE CHANNELS TO USE ======================

useCustomChannels = ~isempty(InspectChannels);

if useCustomChannels

    inspectChannels = unique(double(InspectChannels(:).'), 'stable');

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
        numel(inspectChannels), num2str(inspectChannels));
else
    inspectChannels = saved_count_chns;
    fprintf('Using the ORIGINAL saved count_channels (%d channels): %s\n', ...
        numel(inspectChannels), num2str(inspectChannels));
end

%% ===================== COMPUTE PER-TRIAL METRIC =======================

metricValue = nan(nTrials, 1);

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

        [foundChn, chnPos] = ismember(inspectChannels, CR.Channel_Index);

        if any(~foundChn)
            warning(['Trial %d: some InspectChannels not found in ' ...
                'Channel_Results, skipping those for this trial.'], trial_id);
        end

        chnPos = chnPos(foundChn);
        metricValue(trial_id) = sum(CR.BaselineCorrected_AnalysisWindow_Spikes(chnPos));
    end
end

%% ===================== FLAG TRIALS BY THRESHOLD ========================

isLow  = metricValue < lowThreshold;
isHigh = metricValue > highThreshold;

fprintf('\nTotal trials flagged LOW  (<%g): %d\n', lowThreshold, sum(isLow));
fprintf('Total trials flagged HIGH (>%g): %d\n', highThreshold, sum(isHigh));

%% ===================== PRINT PER CONDITION ============================

fprintf('\n------------------------------------------------------------\n');
fprintf('Flagged trials by condition (relative trial ID = position\n');
fprintf('within that condition, e.g. 1-30):\n');
fprintf('------------------------------------------------------------\n');

nSets = numel(SpikeCounts.set);

for si = 1:nSets

    set_label = SpikeCounts.set(si).set_name;
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

            fprintf('\nSet %d (%s) | %s | PTD %g ms:\n', ...
                si, set_label, amp_label, PTD_ms_val);

            if ~isempty(relativeIDs_low)
                fprintf('  LOW  (<%g):  trials %s\n', ...
                    lowThreshold, num2str(sort(relativeIDs_low(:).')));
            end

            if ~isempty(relativeIDs_high)
                fprintf('  HIGH (>%g):  trials %s\n', ...
                    highThreshold, num2str(sort(relativeIDs_high(:).')));
            end
        end
    end
end

fprintf('\n============================================================\n');
fprintf('INSPECTION COMPLETE (no files modified or saved)\n');
fprintf('============================================================\n');