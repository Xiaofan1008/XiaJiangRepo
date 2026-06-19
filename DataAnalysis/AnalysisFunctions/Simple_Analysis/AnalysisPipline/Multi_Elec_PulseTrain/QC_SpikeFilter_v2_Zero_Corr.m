%% ============================================================
% SpikeFiltering_PulseTrainRecovery_Cleanup.m
%
% Purpose:
%   Clean spike waveforms after pulse-train spike recovery.
%
% Input:
%   *.sp_xia_PrefixRecovery.mat
%   variable: sp_seq
%
% Output:
%   *.sp_xia_PrefixRecovery_SSD.mat
%   final cleaned variable: sp_corr
% ============================================================

clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Single3';

FS = 30000;

%% ================= FILTERING PARAMETERS =================
do_zero_crossing_check = true;

corr_thresh = 0.6;               % lower = more permissive, higher = stricter
do_corr_filter = 1;

% Minimum number of candidate-clean spikes needed to build channel template.
min_template_spikes = 10;

%% ================= WINDOWS =================

% Cleanup window:
% Correlation filter is applied only inside this time range relative to each trigger.
%
% For pulse-train recovery, [0 60] ms is still reasonable because the longest
% 2-electrode 3-pulse interleaved train ends around 25 ms, and the evoked
% response can continue after that.
cleanup_window_ms = [0 60];

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
assert(~isempty(rec_file), 'Cannot find *sp_xia_PrefixRecovery.mat. Please run recovery first.');

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
simultaneous_stim = S_exp.simultaneous_stim;   % rows/slots per trial
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
% TrialParams columns:
%   column 1 = trial number
%   column 2 = condition ID
%   column 3 = internal electrode ID
%
% Each trial has simultaneous_stim rows.
% We use the first row of each trial to get the trial-level condition ID.

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

%% ================= 1) ZERO-CROSSING CHECK =================
sp_zcross = sp_in;
removed_zcross = zeros(nCh,1);

if do_zero_crossing_check

    fprintf('\nZero-Crossing Check (removing waveforms without positive phase)...\n');

    for ch = 1:nCh

        if isempty(sp_zcross{ch})
            continue;
        end

        wfs = sp_zcross{ch}(:,2:end);

        % Keep spikes that contain at least one positive sample.
        has_positive_phase = max(wfs, [], 2) > 0;

        removed = sum(~has_positive_phase);
        removed_zcross(ch) = removed;

        if removed > 0
            sp_zcross{ch} = sp_zcross{ch}(has_positive_phase, :);
            fprintf('Ch %2d: Removed %d spikes by zero-crossing check.\n', ch, removed);
        end
    end
else
    fprintf('\nZero-Crossing Check skipped.\n');
end

%% ================= 2) CHANNEL-LEVEL TEMPLATE CORRELATION FILTER =================

sp_corr = sp_zcross;

corr_template = cell(nCh,1);
corr_source_n = zeros(nCh,1);
removed_corr = zeros(nCh,1);

if do_corr_filter

    fprintf('\nChannel-Level Correlation Filtering...\n');
    fprintf('Template source: all candidate-clean spikes after zero-crossing.\n');
    fprintf('Cleanup window: %.1f to %.1f ms relative to each trigger.\n', ...
        cleanup_window_ms(1), cleanup_window_ms(2));

    for ch = 1:nCh

        if isempty(sp_corr{ch})
            continue;
        end

        spt = sp_corr{ch}(:,1);      % spike times in ms
        wfs = sp_corr{ch}(:,2:end);  % spike waveforms
        nSpikes = size(wfs,1);

        %% ----- A. Build one channel-level template -----
        if nSpikes < min_template_spikes
            fprintf('Ch %2d: Only %d candidate-clean spikes. Skipping correlation filter.\n', ...
                ch, nSpikes);
            continue;
        end

        template = mean(wfs, 1);
        corr_template{ch} = template;
        corr_source_n(ch) = nSpikes;

        %% ----- B. Correlate each spike with the template -----
        corr_vals = zeros(nSpikes,1);

        for i = 1:nSpikes
            corr_vals(i) = corr(template(:), wfs(i,:)');
        end

        % Protect against flat waveforms or zero-variance template/waveform.
        corr_vals(isnan(corr_vals)) = -1;

        %% ----- C. Apply correlation cleanup only inside cleanup_window_ms -----
        bad_mask = false(nSpikes,1);

        for tr = 1:nTrials_use

            t0 = trig(tr) / FS * 1000;

            in_cleanup_window = (spt >= t0 + cleanup_window_ms(1)) & ...
                                (spt <= t0 + cleanup_window_ms(2));

            if any(in_cleanup_window)
                bad_mask(in_cleanup_window & (corr_vals < corr_thresh)) = true;
            end
        end

        removed_count = sum(bad_mask);
        removed_corr(ch) = removed_count;

        sp_corr{ch} = sp_corr{ch}(~bad_mask, :);

        fprintf('Ch %2d: Template N = %d | Removed %d spikes by correlation filter.\n', ...
            ch, nSpikes, removed_count);
    end
else
    fprintf('\nChannel-Level Correlation Filtering skipped.\n');
end

%% ================= REMOVAL SUMMARY BY PULSE-TRAIN CONDITION =================

SummaryRows = {};
row_i = 0;

%% ----- Zero-control summary -----
zero_trials = find(isZeroControl);
zero_trials = zero_trials(zero_trials <= nTrials_use);

if ~isempty(zero_trials)

    [before_count, after_count] = countSpikesInTrials( ...
        sp_zcross, sp_corr, trig, FS, zero_trials, cleanup_window_ms);

    removed_count = before_count - after_count;

    if before_count > 0
        removed_percent = removed_count / before_count * 100;
    else
        removed_percent = NaN;
    end

    row_i = row_i + 1;
    SummaryRows(row_i,1:8) = { ...
        'ZeroControl', ...
        'ZeroControl', ...
        NaN, ...
        NaN, ...
        NaN, ...
        before_count, ...
        after_count, ...
        removed_percent};
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

        [before_count, after_count] = countSpikesInTrials( ...
            sp_zcross, sp_corr, trig, FS, trial_list, cleanup_window_ms);

        removed_count = before_count - after_count;

        if before_count > 0
            removed_percent = removed_count / before_count * 100;
        else
            removed_percent = NaN;
        end

        row_i = row_i + 1;
        SummaryRows(row_i,1:8) = { ...
            sprintf('%s_Level_%g', familyName, level_val), ...
            familyName, ...
            level_val, ...
            NaN, ...
            NaN, ...
            before_count, ...
            after_count, ...
            removed_percent};
    end
end

RemovalSummary = cell2table(SummaryRows, ...
    'VariableNames', {'ConditionLabel', 'ConditionFamily', 'TrainLevel', ...
                      'StimSet', 'Amplitude_uA', ...
                      'BeforeCorr', 'AfterCorr', 'RemovedPercent'});

fprintf('\n================ Removal Summary inside cleanup window ================\n');
disp(RemovalSummary);
fprintf('======================================================================\n');

%% ================= SAVE OUTPUT =================

QC_params = struct();

QC_params.input_file = fname_sp;
QC_params.output_type = 'PulseTrainPrefixRecovery_SSD';

% Record which filters were actually used.
QC_params.do_start_slope_filter = false;
QC_params.do_morphology_filter  = false;
QC_params.do_zero_crossing_check = do_zero_crossing_check;
QC_params.do_corr_filter = do_corr_filter;

% Record the current filter philosophy.
QC_params.filter_pipeline = 'sp_in -> zero_crossing -> correlation -> sp_corr';
QC_params.note = ['Start-slope and morphology filters were intentionally skipped ', ...
                  'to avoid deleting potentially real recovered pulse-train spikes.'];

QC_params.cleanup_window_ms = cleanup_window_ms;

QC_params.corr_thresh = corr_thresh;
QC_params.min_template_spikes = min_template_spikes;
QC_params.template_strategy = 'channel_level_all_zero_crossing_passed_spikes';

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

% Important:
%   Only sp_corr is saved as the cleaned spike variable.
%
% Not saved:
%   sp_in
%   sp_zcross
%   sp_slope
%   sp_morph
%
% Reason:
%   The final downstream analysis should use sp_corr only.

save(out_name, ...
    'sp_corr', ...
    'QC_params', ...
    'corr_template', ...
    'corr_source_n', ...
    'removed_zcross', ...
    'removed_corr', ...
    'RemovalSummary', ...
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

%% ================= LOCAL HELPER FUNCTION =================

function [before_count, after_count] = countSpikesInTrials( ...
    sp_before, sp_after, trig, FS, trial_list, cleanup_window_ms)
    % Count spikes inside cleanup_window_ms across a list of trials.
    %
    % This helper is used only for the removal summary.
    %
    % before_count:
    %   spike count before correlation filtering, using sp_zcross
    %
    % after_count:
    %   spike count after correlation filtering, using sp_corr

    before_count = 0;
    after_count  = 0;

    nCh = numel(sp_before);

    for ch = 1:nCh

        if ~isempty(sp_before{ch})
            spt_before = sp_before{ch}(:,1);
        else
            spt_before = [];
        end

        if ~isempty(sp_after{ch})
            spt_after = sp_after{ch}(:,1);
        else
            spt_after = [];
        end

        for tt = 1:numel(trial_list)

            tr = trial_list(tt);
            t0 = trig(tr) / FS * 1000;

            before_count = before_count + sum(spt_before >= t0 + cleanup_window_ms(1) & ...
                                              spt_before <= t0 + cleanup_window_ms(2));

            after_count = after_count + sum(spt_after >= t0 + cleanup_window_ms(1) & ...
                                            spt_after <= t0 + cleanup_window_ms(2));
        end
    end
end