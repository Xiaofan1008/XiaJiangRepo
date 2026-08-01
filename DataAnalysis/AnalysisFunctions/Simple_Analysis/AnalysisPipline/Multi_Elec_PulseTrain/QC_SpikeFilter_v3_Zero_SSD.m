%% ============================================================
% SpikeFiltering_PulseTrainRecovery_SSD.m
%
% Purpose:
%   Clean spike waveforms after pulse-train spike recovery using:
%
%     1) waveform range threshold
%     2) optional zero-crossing / positive-phase check
%     3) sum-of-squared-difference (SSD) waveform outlier filter
%
% Input:
%   *.sp_xia_PrefixRecovery.mat
%   variable: sp_seq
%
% Output:
%   *.sp_xia_PrefixRecovery_SSD.mat
%   final cleaned variable: sp_corr
%
% Notes:
%   - Correlation filtering is intentionally removed.
%   - SSD filtering is applied globally to all spikes in each channel.
%   - The original recovered spike file is not overwritten.
% ============================================================

clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

FS = 30000;

%% ================= FILTERING PARAMETERS =================

% 1) Waveform range filter
% Remove spikes whose waveform range exceeds this value:
%   range = max(waveform) - min(waveform)
do_range_filter = true;
waveform_range_thresh_uV = 500;

% 2) Zero-crossing / positive-phase check
% Remove waveforms that do not contain any positive phase.
do_zero_crossing_check = true;

% 3) SSD filter
% For each channel:
%   template = mean waveform after range + zero-crossing filtering
%   SSD_i = sum((waveform_i - template).^2)
%   remove if SSD_i > ssd_multiplier * mean(SSD)
do_ssd_filter = true;
ssd_multiplier = 16;

% Minimum number of candidate-clean spikes needed to build SSD template.
min_template_spikes = 10;

%% ================= SUMMARY WINDOW =================
% Filtering is global.
% This window is only used for the RemovalSummary table.
summary_window_ms = [0 60];

%% ================= FOLDER & BASE NAME =================

if ~isfolder(data_folder)
    error('Folder does not exist.');
end

cd(data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');

if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Data folder:\n%s\n', data_folder);
fprintf('Base name: %s\n', base_name);

%% ================= LOAD PULSE-TRAIN RECOVERED SPIKES =================
% Input should be the recovered spike file generated from the recursive
% pulse-train recovery code.
%
% Required variable:
%   sp_seq

rec_file = dir('*sp_xia_PrefixRecovery.mat');
rec_file = rec_file(~contains({rec_file.name}, 'SSD'));

assert(~isempty(rec_file), ...
    'Cannot find non-SSD *sp_xia_PrefixRecovery.mat. Please run recovery first.');

fname_sp = rec_file(1).name;
S_in = load(fname_sp);

if isfield(S_in, 'sp_seq')
    sp_in = S_in.sp_seq;
else
    error('Variable "sp_seq" not found in %s.', fname_sp);
end

nCh = numel(sp_in);

fprintf('\nLoaded recovered spike file: %s\n', fname_sp);
fprintf('Number of spike channels: %d\n', nCh);

%% ================= LOAD TRIGGERS =================

if isempty(dir('*.trig.dat'))
    cur = pwd;
    cleanTrig_sabquick;
    cd(cur);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ================= LOAD EXPERIMENT PARAMETERS =================

param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat found.');

S_exp = load(param_file(1).name);

StimParams        = S_exp.StimParams;
TrialParams       = S_exp.TrialParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;

if isfield(S_exp, 'StimMeta')
    StimMeta = S_exp.StimMeta;
else
    error('StimMeta was not found. This pulse-train filter code requires StimMeta.');
end

if n_Trials ~= nTrig
    warning('n_Trials in exp file (%d) does not match trigger number (%d). Using nTrials_use = min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

%% ================= REMOVE HEADER ROWS =================

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

%% ================= TRIAL CONDITION IDS FROM TrialParams =================

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of condition IDs does not match n_Trials.');
end

%% ================= BUILD PULSE-TRAIN TRIAL METADATA FROM StimMeta =================

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

    % Use actual randomized StimParams rows to get amplitude.
    % This avoids problems if metadata and randomized order differ.
    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);

    if max(rr) <= size(StimParams_data,1)
        amp_vec = cell2mat(StimParams_data(rr,16));
        amp_vec = amp_vec(:)';

        if all(amp_vec <= 0)
            trialAmps(tr) = 0;
        else
            trialAmps(tr) = max(amp_vec(amp_vec > 0));
        end
    end
end

% Keep old variable names for compatibility with later checking scripts.
% In pulse-train files, lastActivePTD is simply the last event time.
lastActivePTD_ms = eventEnd_ms_trial;
lastActivePTD_us = eventEnd_us_trial;

fprintf('\nDetected train levels: ');
disp(unique(trainLevel_trial(~isnan(trainLevel_trial)))');

fprintf('Detected stimulation sets: ');
disp(unique(stimSet_trial(~isnan(stimSet_trial)))');

fprintf('Detected amplitudes: ');
disp(unique(trialAmps(~isnan(trialAmps)))');

fprintf('Detected AutoSim flags: ');
disp(unique(isAutoSim_trial)');

fprintf('Detected final event times / lastActivePTD_ms: ');
disp(unique(lastActivePTD_ms(~isnan(lastActivePTD_ms)))');

%% ================= PRE-COUNT =================

n_spikes_original = zeros(nCh,1);

for ch = 1:nCh
    if isempty(sp_in{ch})
        n_spikes_original(ch) = 0;
    else
        n_spikes_original(ch) = size(sp_in{ch},1);
    end
end

%% ================= 1) GLOBAL WAVEFORM RANGE FILTER =================

sp_range = sp_in;
removed_range = zeros(nCh,1);

if do_range_filter

    fprintf('\nGlobal waveform range filter...\n');
    fprintf('Removing waveforms with range > %.1f uV.\n', waveform_range_thresh_uV);

    for ch = 1:nCh

        if isempty(sp_range{ch})
            continue;
        end

        wfs = sp_range{ch}(:,2:end);

        wf_range = max(wfs, [], 2) - min(wfs, [], 2);
        keep_mask = wf_range <= waveform_range_thresh_uV;

        removed_range(ch) = sum(~keep_mask);

        if removed_range(ch) > 0
            sp_range{ch} = sp_range{ch}(keep_mask, :);
            fprintf('Ch %2d: Removed %d spikes by waveform range filter.\n', ...
                ch, removed_range(ch));
        end
    end

else
    fprintf('\nWaveform range filter skipped.\n');
end

%% ================= 2) GLOBAL ZERO-CROSSING / POSITIVE-PHASE CHECK =================

sp_zcross = sp_range;
removed_zcross = zeros(nCh,1);

if do_zero_crossing_check

    fprintf('\nGlobal zero-crossing / positive-phase check...\n');
    fprintf('Removing waveforms without positive phase.\n');

    for ch = 1:nCh

        if isempty(sp_zcross{ch})
            continue;
        end

        wfs = sp_zcross{ch}(:,2:end);

        % Keep spikes that contain at least one positive sample.
        has_positive_phase = max(wfs, [], 2) > 0;

        removed_zcross(ch) = sum(~has_positive_phase);

        if removed_zcross(ch) > 0
            sp_zcross{ch} = sp_zcross{ch}(has_positive_phase, :);
            fprintf('Ch %2d: Removed %d spikes by zero-crossing check.\n', ...
                ch, removed_zcross(ch));
        end
    end

else
    fprintf('\nZero-crossing / positive-phase check skipped.\n');
end

%% ================= 3) GLOBAL SSD FILTER =================

sp_corr = sp_zcross;

ssd_template = cell(nCh,1);
ssd_source_n = zeros(nCh,1);
ssd_mean = nan(nCh,1);
ssd_thresh = nan(nCh,1);
removed_ssd = zeros(nCh,1);

if do_ssd_filter

    fprintf('\nGlobal SSD waveform outlier filter...\n');
    fprintf('Template source: all spikes after range + zero-crossing filtering.\n');
    fprintf('SSD threshold: SSD > %.1f x mean SSD.\n', ssd_multiplier);

    for ch = 1:nCh

        if isempty(sp_corr{ch})
            continue;
        end

        wfs = sp_corr{ch}(:,2:end);
        nSpikes = size(wfs,1);

        if nSpikes < min_template_spikes
            fprintf('Ch %2d: Only %d candidate-clean spikes. Skipping SSD filter.\n', ...
                ch, nSpikes);
            continue;
        end

        % Build channel-level mean waveform template after range + zero-crossing.
        template = mean(wfs, 1);
        ssd_template{ch} = template;
        ssd_source_n(ch) = nSpikes;

        % Compute sum-of-squared difference from template.
        diff_wfs = wfs - template;
        ssd_vals = sum(diff_wfs.^2, 2);

        mean_ssd_this = mean(ssd_vals);

        ssd_mean(ch) = mean_ssd_this;
        ssd_thresh(ch) = ssd_multiplier * mean_ssd_this;

        % Protect against perfectly identical/flat waveforms.
        if isnan(mean_ssd_this) || mean_ssd_this <= 0
            fprintf('Ch %2d: mean SSD is %.3g. Skipping SSD filter.\n', ...
                ch, mean_ssd_this);
            continue;
        end

        keep_mask = ssd_vals <= ssd_thresh(ch);

        removed_ssd(ch) = sum(~keep_mask);

        if removed_ssd(ch) > 0
            sp_corr{ch} = sp_corr{ch}(keep_mask, :);
        end

        fprintf('Ch %2d: Template N = %d | mean SSD = %.3g | threshold = %.3g | Removed %d spikes by SSD.\n', ...
            ch, nSpikes, ssd_mean(ch), ssd_thresh(ch), removed_ssd(ch));
    end

else
    fprintf('\nSSD filter skipped.\n');
end

%% ================= FINAL PER-CHANNEL COUNT SUMMARY =================

n_after_range = zeros(nCh,1);
n_after_zcross = zeros(nCh,1);
n_after_ssd = zeros(nCh,1);

for ch = 1:nCh

    if isempty(sp_range{ch})
        n_after_range(ch) = 0;
    else
        n_after_range(ch) = size(sp_range{ch},1);
    end

    if isempty(sp_zcross{ch})
        n_after_zcross(ch) = 0;
    else
        n_after_zcross(ch) = size(sp_zcross{ch},1);
    end

    if isempty(sp_corr{ch})
        n_after_ssd(ch) = 0;
    else
        n_after_ssd(ch) = size(sp_corr{ch},1);
    end
end

ChannelRemovalSummary = table( ...
    (1:nCh)', ...
    n_spikes_original, ...
    removed_range, ...
    n_after_range, ...
    removed_zcross, ...
    n_after_zcross, ...
    removed_ssd, ...
    n_after_ssd, ...
    ssd_source_n, ...
    ssd_mean, ...
    ssd_thresh, ...
    'VariableNames', { ...
        'SpikeChannel', ...
        'OriginalSpikeCount', ...
        'RemovedByRange', ...
        'AfterRangeSpikeCount', ...
        'RemovedByZeroCrossing', ...
        'AfterZeroCrossingSpikeCount', ...
        'RemovedBySSD', ...
        'FinalSpikeCount', ...
        'SSDTemplateSourceN', ...
        'MeanSSD', ...
        'SSDThreshold'});

fprintf('\n================ Channel Removal Summary ================\n');
disp(ChannelRemovalSummary);
fprintf('=========================================================\n');

%% ================= REMOVAL SUMMARY BY PULSE-TRAIN CONDITION =================
% This summary is only for quick QC in the stimulation-related window.
% Filtering itself was applied globally.

SummaryRows = {};
row_i = 0;

%% ----- Zero-control summary -----
zero_trials = find(isZeroControl);
zero_trials = zero_trials(zero_trials <= nTrials_use);

if ~isempty(zero_trials)

    [original_count, range_count, zcross_count, final_count] = countSpikesInTrials_4Stages( ...
        sp_in, sp_range, sp_zcross, sp_corr, trig, FS, zero_trials, summary_window_ms);

    row_i = row_i + 1;
    SummaryRows(row_i,1:11) = { ...
        'ZeroControl', ...
        'ZeroControl', ...
        NaN, ...
        NaN, ...
        NaN, ...
        original_count, ...
        range_count, ...
        zcross_count, ...
        final_count, ...
        original_count - final_count, ...
        percentRemoved(original_count, final_count)};
end

%% ----- Seq and AutoSim level summary -----
for autoFlag = [0 1]

    if autoFlag == 0
        familyName = 'Seq';
    else
        familyName = 'AutoSim';
    end

    level_list = unique(trainLevel_trial(isAutoSim_trial == logical(autoFlag) & ...
                                         ~isZeroControl & ...
                                         ~isnan(trainLevel_trial)));
    level_list = sort(level_list(:))';

    for li = 1:numel(level_list)

        level_val = level_list(li);

        trial_list = find(isAutoSim_trial == logical(autoFlag) & ...
                          isZeroControl == 0 & ...
                          trainLevel_trial == level_val);

        trial_list = trial_list(trial_list <= nTrials_use);

        if isempty(trial_list)
            continue;
        end

        [original_count, range_count, zcross_count, final_count] = countSpikesInTrials_4Stages( ...
            sp_in, sp_range, sp_zcross, sp_corr, trig, FS, trial_list, summary_window_ms);

        row_i = row_i + 1;
        SummaryRows(row_i,1:11) = { ...
            sprintf('%s_Level_%g', familyName, level_val), ...
            familyName, ...
            level_val, ...
            NaN, ...
            NaN, ...
            original_count, ...
            range_count, ...
            zcross_count, ...
            final_count, ...
            original_count - final_count, ...
            percentRemoved(original_count, final_count)};
    end
end

RemovalSummary = cell2table(SummaryRows, ...
    'VariableNames', { ...
        'ConditionLabel', ...
        'ConditionFamily', ...
        'TrainLevel', ...
        'StimSet', ...
        'Amplitude_uA', ...
        'OriginalCountInSummaryWindow', ...
        'AfterRangeCountInSummaryWindow', ...
        'AfterZeroCrossingCountInSummaryWindow', ...
        'FinalCountInSummaryWindow', ...
        'TotalRemovedInSummaryWindow', ...
        'RemovedPercentInSummaryWindow'});

fprintf('\n================ Removal Summary inside summary window ================\n');
fprintf('Summary window: %.1f to %.1f ms relative to trigger.\n', ...
    summary_window_ms(1), summary_window_ms(2));
disp(RemovalSummary);
fprintf('======================================================================\n');

%% ================= SAVE OUTPUT =================

QC_params = struct();

QC_params.input_file = fname_sp;
QC_params.output_type = 'PulseTrainPrefixRecovery_SSD';

QC_params.do_range_filter = do_range_filter;
QC_params.waveform_range_thresh_uV = waveform_range_thresh_uV;

QC_params.do_zero_crossing_check = do_zero_crossing_check;

QC_params.do_ssd_filter = do_ssd_filter;
QC_params.ssd_multiplier = ssd_multiplier;
QC_params.min_template_spikes = min_template_spikes;
QC_params.ssd_template_strategy = 'channel_level_mean_after_range_and_zero_crossing';
QC_params.ssd_filter_scope = 'global_all_spikes';

QC_params.do_corr_filter = false;
QC_params.corr_filter_note = 'Correlation filter was intentionally removed to avoid over-filtering. SSD filter is used instead.';

QC_params.filter_pipeline = 'sp_in -> global_range -> global_zero_crossing -> global_SSD -> sp_corr';
QC_params.note = ['Waveform range filter and SSD filter follow the lab-style waveform outlier rejection logic. ', ...
                  'SSD is computed relative to the channel mean waveform after range and zero-crossing filtering.'];

QC_params.summary_window_ms = summary_window_ms;
QC_params.FS = FS;

% Output name:
%   Input:
%     xxx.sp_xia_PrefixRecovery.mat
%   Output:
%     xxx.sp_xia_PrefixRecovery_SSD.mat
if contains(fname_sp, '.sp_xia_PrefixRecovery.mat')
    out_name = strrep(fname_sp, '.sp_xia_PrefixRecovery.mat', '.sp_xia_PrefixRecovery_SSD.mat');
else
    out_name = [base_name '.sp_xia_PrefixRecovery_SSD.mat'];
end

if isfile(out_name)
    delete(out_name);
end

save(out_name, ...
    'sp_corr', ...
    'QC_params', ...
    'ChannelRemovalSummary', ...
    'RemovalSummary', ...
    'ssd_template', ...
    'ssd_source_n', ...
    'ssd_mean', ...
    'ssd_thresh', ...
    'removed_range', ...
    'removed_zcross', ...
    'removed_ssd', ...
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
    'n_Trials', ...
    'nTrials_use', ...
    '-v7.3');

fprintf('\nSaved cleaned pulse-train recovered spikes to:\n%s\n', fullfile(data_folder, out_name));

fprintf('\nTotal original spikes:       %d\n', sum(n_spikes_original));
fprintf('Total removed by range:      %d\n', sum(removed_range));
fprintf('Total removed by zcross:     %d\n', sum(removed_zcross));
fprintf('Total removed by SSD:        %d\n', sum(removed_ssd));
fprintf('Total final spikes:          %d\n', sum(n_after_ssd));

fprintf('\nFinished pulse-train recovered spike filtering.\n');

%% ================= LOCAL HELPER FUNCTIONS =================

function [count_original, count_range, count_zcross, count_final] = countSpikesInTrials_4Stages( ...
    sp_original, sp_range, sp_zcross, sp_final, trig, FS, trial_list, summary_window_ms)
    % Count spikes inside summary_window_ms across a list of trials.
    %
    % This helper is used only for the removal summary.
    % Filtering itself is global.

    count_original = 0;
    count_range    = 0;
    count_zcross   = 0;
    count_final    = 0;

    nCh = numel(sp_original);

    for ch = 1:nCh

        if ~isempty(sp_original{ch})
            spt_original = sp_original{ch}(:,1);
        else
            spt_original = [];
        end

        if ~isempty(sp_range{ch})
            spt_range = sp_range{ch}(:,1);
        else
            spt_range = [];
        end

        if ~isempty(sp_zcross{ch})
            spt_zcross = sp_zcross{ch}(:,1);
        else
            spt_zcross = [];
        end

        if ~isempty(sp_final{ch})
            spt_final = sp_final{ch}(:,1);
        else
            spt_final = [];
        end

        for tt = 1:numel(trial_list)

            tr = trial_list(tt);

            t0 = trig(tr) / FS * 1000;

            win_start = t0 + summary_window_ms(1);
            win_end   = t0 + summary_window_ms(2);

            count_original = count_original + sum(spt_original >= win_start & spt_original < win_end);
            count_range    = count_range    + sum(spt_range    >= win_start & spt_range    < win_end);
            count_zcross   = count_zcross   + sum(spt_zcross   >= win_start & spt_zcross   < win_end);
            count_final    = count_final    + sum(spt_final    >= win_start & spt_final    < win_end);
        end
    end
end

function p = percentRemoved(before_count, after_count)

    if before_count > 0
        p = (before_count - after_count) / before_count * 100;
    else
        p = NaN;
    end
end