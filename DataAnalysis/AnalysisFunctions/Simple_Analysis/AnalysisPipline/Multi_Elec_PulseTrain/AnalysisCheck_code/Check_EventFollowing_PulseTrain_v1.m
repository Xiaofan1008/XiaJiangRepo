%% ============================================================
%  Quick Epoch-wise Response Analysis for Pulse-Train Data
%
%  Purpose:
%    Compare how response is distributed across non-overlapping 5 ms epochs
%    for sequential vs AutoSim pulse trains.
%
%  Why this analysis:
%    Sequential events are spaced every 5 ms:
%       0, 5, 10, 15, 20, 25 ms
%
%    Therefore, event-locked windows like 2–8 ms would overlap.
%    Instead, this code uses fixed non-overlapping 5 ms epochs:
%       0–5, 5–10, 10–15, 15–20, 20–25, 25–30 ms
%
%  Interpretation:
%    AutoSim may show stronger responses in epochs after simultaneous events:
%       0–5, 10–15, 20–25 ms
%
%    Sequential may show more continuous / evenly distributed response
%    across all six epochs.
%
%  Input:
%    Default: *.sp_xia_PrefixRecovery_SSD.mat
%             variable: sp_corr
%
%  Important:
%    This is a quick checker.
%    It does NOT exclude bad trials.
%    It does NOT use responding-channel files.
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train1';

% Recording channels to analyse.
channels_to_analyse = [1:7 17:28];

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Spike source:
%   'recovery_ssd' = load filtered recovered spikes, sp_corr
%   'recovery'     = load recovered but unfiltered spikes, sp_seq
spike_source = 'recovery_ssd';

% Non-overlapping epoch edges.
%
% For sequential Level 6 with PTD = 5 ms:
%   Seq events = 0, 5, 10, 15, 20, 25 ms
%
% So these epochs correspond to:
%   0–5, 5–10, 10–15, 15–20, 20–25, 25–30 ms
epoch_edges_ms = 0:5:30;

% Baseline window.
baseline_window_ms = [-60 -10];

% Apply baseline correction to each epoch?
do_baseline_correction = true;

% Remove negative baseline-corrected values before calculating distribution metrics?
%
% Entropy and fraction metrics require non-negative values.
% This option is only used for distribution metrics, not for raw epoch counts.
set_negative_bins_to_zero_for_metrics = true;

% Final train levels to compare.
seq_train_level = 6;
autosim_train_level = 3;

% Plot these families.
plot_sequential = true;
plot_auto_sim = true;

% Stimulation set IDs to include.
% [] means all detected sets.
SetIDs_to_plot = [];

% Amplitudes to include.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% Pool sequential sets?
%
% Since the current question is Sequential vs AutoSim, not order difference,
% default is true.
pool_seq_sets = true;

% Pool AutoSim sets?
pool_autosim_sets = true;

% Save result?
save_epoch_result = false;

%% ====================== CHECK FOLDER ======================

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end

cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

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

switch lower(spike_source)

    case 'recovery_ssd'
        sp_files = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));
        assert(~isempty(sp_files), 'No *sp_xia_PrefixRecovery_SSD.mat file found.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading filtered recovered spike file:\n%s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp_corr')
            sp_use = S_sp.sp_corr;
        else
            error('Variable sp_corr not found in %s.', sp_filename);
        end

    case 'recovery'
        sp_files = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
        sp_files = sp_files(~contains({sp_files.name}, 'SSD'));
        assert(~isempty(sp_files), 'No *sp_xia_PrefixRecovery.mat file found.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading recovered spike file:\n%s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp_seq')
            sp_use = S_sp.sp_seq;
        else
            error('Variable sp_seq not found in %s.', sp_filename);
        end

    otherwise
        error('Unknown spike_source: %s. Use recovery_ssd or recovery.', spike_source);
end

nCh = numel(sp_use);
fprintf('Number of spike channels: %d\n', nCh);

%% ====================== LOAD TRIGGERS ======================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD EXPERIMENT DATAFILE ======================

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
TrialParams       = S_exp.TrialParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;
E_MAP             = S_exp.E_MAP;

if isfield(S_exp, 'StimMeta')
    StimMeta = S_exp.StimMeta;
else
    error('StimMeta was not found. This code requires StimMeta.');
end

if isfield(S_exp, 'trainMode')
    trainMode = S_exp.trainMode;
else
    trainMode = NaN;
end

fprintf('Loaded exp datafile: %s\n', fileDIR(1).name);
fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);
fprintf('trainMode: %g\n', trainMode);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match trigger number (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== SAMPLING RATE ======================

FS = 30000;

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

%% ====================== TRIAL-LEVEL CONDITION IDS ======================

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA FROM StimMeta ======================

stimSet_trial      = NaN(n_Trials,1);
trainLevel_trial   = NaN(n_Trials,1);
totalLevels_trial  = NaN(n_Trials,1);
amp_trial          = NaN(n_Trials,1);
isAutoSim_trial    = false(n_Trials,1);
isZero_trial       = false(n_Trials,1);
eventEnd_ms_trial  = NaN(n_Trials,1);

eventTimes_ms_trial = cell(n_Trials,1);
pulseCount_trial    = cell(n_Trials,1);

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
        isZero_trial(tr) = logical(meta.IsZeroCurrentControl);
    end

    if isfield(meta, 'EventTimesIncluded_ms')
        eventTimes_ms_trial{tr} = meta.EventTimesIncluded_ms;
    else
        eventTimes_ms_trial{tr} = [];
    end

    if isfield(meta, 'EventEndTime_ms')
        eventEnd_ms_trial(tr) = meta.EventEndTime_ms;
    elseif ~isempty(eventTimes_ms_trial{tr})
        eventEnd_ms_trial(tr) = max(eventTimes_ms_trial{tr});
    else
        eventEnd_ms_trial(tr) = 0;
    end

    if isfield(meta, 'PulseCountPerElectrode')
        pulseCount_trial{tr} = meta.PulseCountPerElectrode;
    else
        pulseCount_trial{tr} = [];
    end

    % Use actual randomized StimParams rows to get amplitude.
    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);

    if max(rr) <= size(StimParams_data,1)

        amp_vec = cell2mat(StimParams_data(rr,16));
        amp_vec = amp_vec(:)';

        if all(amp_vec <= 0)
            amp_trial(tr) = 0;
        else
            amp_trial(tr) = max(amp_vec(amp_vec > 0));
        end
    end
end

%% ====================== APPLY USER FILTERS ======================

% Amplitudes.
all_amps = unique(amp_trial(~isnan(amp_trial)));
all_amps = all_amps(all_amps > 0);

if isempty(Amps_to_plot)
    Amps_selected = all_amps;
else
    Amps_selected = intersect(all_amps, Amps_to_plot);
end

% Set IDs.
all_sets = unique(stimSet_trial(~isnan(stimSet_trial)));
all_sets = all_sets(all_sets > 0);

if isempty(SetIDs_to_plot)
    SetIDs_selected = all_sets;
else
    SetIDs_selected = intersect(all_sets, SetIDs_to_plot);
end

fprintf('\nDetected amplitudes: ');
disp(all_amps');

fprintf('Selected amplitudes: ');
disp(Amps_selected');

fprintf('\nDetected set IDs: ');
disp(all_sets');

fprintf('Selected set IDs: ');
disp(SetIDs_selected');

fprintf('\nDetected Seq levels: ');
disp(unique(trainLevel_trial(isAutoSim_trial == 0 & isZero_trial == 0))');

fprintf('Detected AutoSim levels: ');
disp(unique(trainLevel_trial(isAutoSim_trial == 1 & isZero_trial == 0))');

fprintf('\nEpoch edges: ');
disp(epoch_edges_ms);

fprintf('Baseline window: %.1f to %.1f ms\n', baseline_window_ms(1), baseline_window_ms(2));
fprintf('Baseline correction: %d\n', do_baseline_correction);
fprintf('Pool Seq sets: %d\n', pool_seq_sets);
fprintf('Pool AutoSim sets: %d\n', pool_autosim_sets);

%% ====================== BUILD CONDITION LIST ======================
% With pooling:
%   Sequential pooled L6
%   AutoSim pooled L3
%
% Without pooling:
%   Seq Set2 L6
%   Seq Set3 L6
%   AutoSim Set1 L3

CondList = struct( ...
    'Family', {}, ...
    'SetID', {}, ...
    'Level', {}, ...
    'Label', {}, ...
    'ShortLabel', {}, ...
    'AutoFlag', {}, ...
    'TrialList', {});

cond_i = 0;

%% ----- Sequential conditions -----
if plot_sequential

    if pool_seq_sets

        trial_pool = find(isAutoSim_trial == 0 & ...
                          isZero_trial == 0 & ...
                          trainLevel_trial == seq_train_level & ...
                          ismember(amp_trial, Amps_selected) & ...
                          ismember(stimSet_trial, SetIDs_selected));

        trial_pool = trial_pool(trial_pool <= nTrials_use);

        if ~isempty(trial_pool)
            cond_i = cond_i + 1;

            CondList(cond_i).Family = 'Seq';
            CondList(cond_i).SetID = NaN;
            CondList(cond_i).Level = seq_train_level;
            CondList(cond_i).Label = sprintf('Sequential pooled L%d', seq_train_level);
            CondList(cond_i).ShortLabel = 'Sequential';
            CondList(cond_i).AutoFlag = false;
            CondList(cond_i).TrialList = trial_pool;
        end

    else

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trial_this = find(isAutoSim_trial == 0 & ...
                              isZero_trial == 0 & ...
                              trainLevel_trial == seq_train_level & ...
                              stimSet_trial == setID & ...
                              ismember(amp_trial, Amps_selected));

            trial_this = trial_this(trial_this <= nTrials_use);

            if isempty(trial_this)
                continue;
            end

            setLabel = buildSetLabelFromTrials(trial_this, StimParams_data, ...
                simultaneous_stim, E_MAP, false);

            cond_i = cond_i + 1;

            CondList(cond_i).Family = 'Seq';
            CondList(cond_i).SetID = setID;
            CondList(cond_i).Level = seq_train_level;
            CondList(cond_i).Label = sprintf('Seq Set%d %s L%d', setID, setLabel, seq_train_level);
            CondList(cond_i).ShortLabel = sprintf('Seq S%d', setID);
            CondList(cond_i).AutoFlag = false;
            CondList(cond_i).TrialList = trial_this;
        end
    end
end

%% ----- AutoSim conditions -----
if plot_auto_sim

    if pool_autosim_sets

        trial_pool = find(isAutoSim_trial == 1 & ...
                          isZero_trial == 0 & ...
                          trainLevel_trial == autosim_train_level & ...
                          ismember(amp_trial, Amps_selected) & ...
                          ismember(stimSet_trial, SetIDs_selected));

        trial_pool = trial_pool(trial_pool <= nTrials_use);

        if ~isempty(trial_pool)
            cond_i = cond_i + 1;

            CondList(cond_i).Family = 'AutoSim';
            CondList(cond_i).SetID = NaN;
            CondList(cond_i).Level = autosim_train_level;
            CondList(cond_i).Label = sprintf('AutoSim pooled L%d', autosim_train_level);
            CondList(cond_i).ShortLabel = 'AutoSim';
            CondList(cond_i).AutoFlag = true;
            CondList(cond_i).TrialList = trial_pool;
        end

    else

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trial_this = find(isAutoSim_trial == 1 & ...
                              isZero_trial == 0 & ...
                              trainLevel_trial == autosim_train_level & ...
                              stimSet_trial == setID & ...
                              ismember(amp_trial, Amps_selected));

            trial_this = trial_this(trial_this <= nTrials_use);

            if isempty(trial_this)
                continue;
            end

            setLabel = buildSetLabelFromTrials(trial_this, StimParams_data, ...
                simultaneous_stim, E_MAP, true);

            cond_i = cond_i + 1;

            CondList(cond_i).Family = 'AutoSim';
            CondList(cond_i).SetID = setID;
            CondList(cond_i).Level = autosim_train_level;
            CondList(cond_i).Label = sprintf('AutoSim Set%d %s L%d', setID, setLabel, autosim_train_level);
            CondList(cond_i).ShortLabel = sprintf('AutoSim S%d', setID);
            CondList(cond_i).AutoFlag = true;
            CondList(cond_i).TrialList = trial_this;
        end
    end
end

if isempty(CondList)
    error('No matching conditions found. Check set IDs, amplitude, and train level settings.');
end

fprintf('\n================ Conditions for epoch-wise analysis ================\n');
for c = 1:numel(CondList)
    fprintf('%2d) %s | total trials across selected amps = %d\n', ...
        c, CondList(c).Label, numel(CondList(c).TrialList));
end
fprintf('===================================================================\n');

ShortLabels = {CondList.ShortLabel};

%% ====================== CHANNEL MAP AND EPOCH SETUP ======================

d = Depth_s(Electrode_Type);

channels_to_analyse = channels_to_analyse(:)';
nAnalyseCh = numel(channels_to_analyse);

epoch_edges_ms = epoch_edges_ms(:)';
epoch_centres_ms = epoch_edges_ms(1:end-1) + diff(epoch_edges_ms)/2;
nEpochs = numel(epoch_edges_ms) - 1;
epoch_width_ms = diff(epoch_edges_ms);

base_dur_ms = baseline_window_ms(2) - baseline_window_ms(1);

if base_dur_ms <= 0
    error('baseline_window_ms must have positive duration.');
end

%% ====================== BUILD EPOCH RESPONSE MATRICES ======================
% EpochMat dimensions:
%   selected channel x epoch x condition x amplitude
%
% Each value:
%   mean baseline-corrected spike count / trial

nCond = numel(CondList);
nAmp = numel(Amps_selected);

EpochMat = NaN(nAnalyseCh, nEpochs, nCond, nAmp);
EpochSEMMat = NaN(nAnalyseCh, nEpochs, nCond, nAmp);

for c = 1:nCond

    for aa = 1:nAmp

        amp_val = Amps_selected(aa);

        trial_list = CondList(c).TrialList;
        trial_list = trial_list(amp_trial(trial_list) == amp_val);
        trial_list = trial_list(trial_list <= nTrials_use);

        if isempty(trial_list)
            continue;
        end

        for ich_i = 1:nAnalyseCh

            depth_idx = channels_to_analyse(ich_i);
            spike_ch = d(depth_idx);

            if spike_ch > nCh || isempty(sp_use{spike_ch})
                continue;
            end

            spike_times = sp_use{spike_ch}(:,1);

            trial_epoch_counts = NaN(numel(trial_list), nEpochs);

            for tt = 1:numel(trial_list)

                tr = trial_list(tt);
                t0_ms = trig(tr) / FS * 1000;

                rel_t = spike_times - t0_ms;

                baseline_count = sum(rel_t >= baseline_window_ms(1) & ...
                                     rel_t <  baseline_window_ms(2));

                baseline_rate_per_ms = baseline_count / base_dur_ms;

                raw_epoch_counts = histcounts(rel_t, epoch_edges_ms);

                expected_baseline_per_epoch = baseline_rate_per_ms .* epoch_width_ms;

                if do_baseline_correction
                    corrected_epoch_counts = raw_epoch_counts - expected_baseline_per_epoch;
                else
                    corrected_epoch_counts = raw_epoch_counts;
                end

                trial_epoch_counts(tt,:) = corrected_epoch_counts;
            end

            EpochMat(ich_i,:,c,aa) = mean(trial_epoch_counts, 1, 'omitnan');
            EpochSEMMat(ich_i,:,c,aa) = std(trial_epoch_counts, 0, 1, 'omitnan') ./ ...
                                        sqrt(sum(~isnan(trial_epoch_counts),1));
        end
    end
end

%% ====================== CALCULATE EPOCH METRICS ======================
% For each condition and amplitude:
%
% Epoch temporal vector:
%   summed across selected channels for each epoch
%
% We calculate:
%   total response
%   epoch fraction
%   peak epoch fraction
%   epoch entropy
%   effective active epochs
%   evenness index

MetricRows = {};
row_i = 0;

EpochProfile = NaN(nCond, nAmp, nEpochs);
EpochFraction = NaN(nCond, nAmp, nEpochs);

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    for c = 1:nCond

        M_raw = squeeze(EpochMat(:,:,c,aa));

        if all(isnan(M_raw(:)))
            continue;
        end

        M = M_raw;
        M(isnan(M)) = 0;

        M_metric = M;

        if set_negative_bins_to_zero_for_metrics
            M_metric(M_metric < 0) = 0;
        end

        epoch_vec = sum(M_metric, 1);
        total_response = sum(epoch_vec);

        EpochProfile(c,aa,:) = epoch_vec;

        if total_response <= 0
            epoch_fraction = NaN(1,nEpochs);
            peak_epoch_fraction = NaN;
            epoch_entropy = NaN;
            epoch_entropy_norm = NaN;
            effective_epochs = NaN;
            evenness_index = NaN;
        else
            epoch_fraction = epoch_vec ./ total_response;

            p = epoch_fraction(epoch_fraction > 0);

            epoch_entropy = -sum(p .* log(p));
            epoch_entropy_norm = epoch_entropy / log(nEpochs);
            effective_epochs = exp(epoch_entropy);

            % Same as normalized entropy, but named intuitively.
            % 1 means response is evenly distributed across epochs.
            % 0 means response is concentrated in one epoch.
            evenness_index = epoch_entropy_norm;

            peak_epoch_fraction = max(epoch_fraction);
        end

        EpochFraction(c,aa,:) = epoch_fraction;

        row_i = row_i + 1;

        MetricRows(row_i,1:13) = { ...
            CondList(c).Family, ...
            CondList(c).SetID, ...
            CondList(c).Level, ...
            CondList(c).Label, ...
            amp_val, ...
            total_response, ...
            peak_epoch_fraction, ...
            epoch_entropy, ...
            epoch_entropy_norm, ...
            effective_epochs, ...
            evenness_index, ...
            nEpochs, ...
            nAnalyseCh};
    end
end

EpochMetricTable = cell2table(MetricRows, ...
    'VariableNames', { ...
    'Family', ...
    'SetID', ...
    'TrainLevel', ...
    'ConditionLabel', ...
    'Amplitude_uA', ...
    'TotalResponse', ...
    'PeakEpochFraction', ...
    'EpochEntropy', ...
    'EpochEntropyNorm', ...
    'EffectiveEpochs', ...
    'EvennessIndex', ...
    'NEpochs', ...
    'NSelectedChannels'});

fprintf('\n================ Epoch-wise Response Metrics ================\n');
disp(EpochMetricTable);
fprintf('=============================================================\n');

%% ====================== FIGURE 1: EPOCH RESPONSE PROFILE ======================
% One figure per amplitude.
%
% Y-axis:
%   summed baseline-corrected spike count across selected channels.
%
% This shows absolute response in each non-overlapping epoch.

cond_colors = lines(nCond);

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    figure('Color','w', ...
           'Name', sprintf('EpochProfile_Response_Amp%g', amp_val), ...
           'Position', [100 100 900 520]);

    hold on; box off;

    for c = 1:nCond

        y = squeeze(EpochProfile(c,aa,:));

        if all(isnan(y))
            continue;
        end

        plot(epoch_centres_ms, y, '-o', ...
            'Color', cond_colors(c,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 7, ...
            'DisplayName', CondList(c).ShortLabel);
    end

    xlabel('Epoch centre time after trigger (ms)');
    ylabel('Summed baseline-corrected spike count');

    title(sprintf('Epoch-wise response profile | %.1f \\muA', amp_val), ...
        'Interpreter','none');

    legend('Box','off','Location','best');
    grid on;

    xlim([epoch_edges_ms(1) epoch_edges_ms(end)]);

    xticks(epoch_centres_ms);
    xticklabels(arrayfun(@(i) sprintf('%.0f–%.0f', epoch_edges_ms(i), epoch_edges_ms(i+1)), ...
        1:nEpochs, 'UniformOutput', false));
    xtickangle(30);
end

%% ====================== FIGURE 2: EPOCH RESPONSE FRACTION ======================
% One figure per amplitude.
%
% Y-axis:
%   fraction of total response in each epoch.
%
% This removes total spike-count differences and focuses on distribution shape.

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    figure('Color','w', ...
           'Name', sprintf('EpochProfile_Fraction_Amp%g', amp_val), ...
           'Position', [130 130 900 520]);

    hold on; box off;

    for c = 1:nCond

        y = squeeze(EpochFraction(c,aa,:));

        if all(isnan(y))
            continue;
        end

        plot(epoch_centres_ms, y, '-o', ...
            'Color', cond_colors(c,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 7, ...
            'DisplayName', CondList(c).ShortLabel);
    end

    xlabel('Epoch centre time after trigger (ms)');
    ylabel('Fraction of total response');

    title(sprintf('Epoch-wise response fraction | %.1f \\muA', amp_val), ...
        'Interpreter','none');

    legend('Box','off','Location','best');
    grid on;

    xlim([epoch_edges_ms(1) epoch_edges_ms(end)]);

    xticks(epoch_centres_ms);
    xticklabels(arrayfun(@(i) sprintf('%.0f–%.0f', epoch_edges_ms(i), epoch_edges_ms(i+1)), ...
        1:nEpochs, 'UniformOutput', false));
    xtickangle(30);
end

%% ====================== FIGURE 3: EPOCH METRIC SUMMARY ======================
% One figure per amplitude.
%
% Metrics:
%   Total response
%   Peak epoch fraction
%   Normalized epoch entropy
%   Effective epochs
%   Evenness index

metricNames = { ...
    'TotalResponse', ...
    'PeakEpochFraction', ...
    'EpochEntropyNorm', ...
    'EffectiveEpochs', ...
    'EvennessIndex'};

metricTitles = { ...
    'Total response', ...
    'Peak epoch fraction', ...
    'Normalized epoch entropy', ...
    'Effective active epochs', ...
    'Evenness index'};

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    figure('Color','w', ...
           'Name', sprintf('EpochMetrics_Amp%g', amp_val), ...
           'Position', [120 80 1100 650]);

    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for mm = 1:numel(metricNames)

        ax = nexttile;
        hold(ax,'on'); box(ax,'off');

        vals = NaN(1,nCond);

        for c = 1:nCond

            rowMatch = strcmp(EpochMetricTable.ConditionLabel, CondList(c).Label) & ...
                       EpochMetricTable.Amplitude_uA == amp_val;

            if any(rowMatch)
                vals(c) = EpochMetricTable{find(rowMatch,1,'first'), metricNames{mm}};
            end
        end

        bar(ax, vals);

        xticks(ax, 1:nCond);
        xticklabels(ax, ShortLabels);
        xtickangle(ax, 30);

        ylabel(ax, metricTitles{mm});
        title(ax, metricTitles{mm}, 'Interpreter','none');
        grid(ax,'on');

        for c = 1:nCond
            if ~isnan(vals(c))
                text(ax, c, vals(c), sprintf('%.2f', vals(c)), ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom');
            end
        end
    end

    sgtitle(sprintf('Epoch-wise response metrics | %.1f \\muA', amp_val), ...
        'Interpreter','none');
end

%% ====================== FIGURE 4: SEQ VS AUTOSIM ACROSS AMPLITUDES ======================
% Most useful when there are only two conditions:
%   Sequential pooled
%   AutoSim pooled

if nCond == 2

    figure('Color','w', ...
           'Name', 'EpochMetrics_SeqVsAutoSim_AcrossAmps', ...
           'Position', [150 100 1100 650]);

    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for mm = 1:numel(metricNames)

        ax = nexttile;
        hold(ax,'on'); box(ax,'off');

        for c = 1:nCond

            vals = NaN(1,nAmp);

            for aa = 1:nAmp

                amp_val = Amps_selected(aa);

                rowMatch = strcmp(EpochMetricTable.ConditionLabel, CondList(c).Label) & ...
                           EpochMetricTable.Amplitude_uA == amp_val;

                if any(rowMatch)
                    vals(aa) = EpochMetricTable{find(rowMatch,1,'first'), metricNames{mm}};
                end
            end

            plot(ax, Amps_selected, vals, '-o', ...
                'Color', cond_colors(c,:), ...
                'LineWidth', 2, ...
                'MarkerSize', 7, ...
                'DisplayName', CondList(c).ShortLabel);
        end

        xlabel(ax, 'Amplitude (\muA)');
        ylabel(ax, metricTitles{mm});
        title(ax, metricTitles{mm}, 'Interpreter','none');
        grid(ax,'on');

        if mm == 1
            legend(ax, 'Box','off','Location','best');
        end
    end

    sgtitle('Epoch-wise response metrics: Sequential vs AutoSim', ...
        'Interpreter','none');
end

%% ====================== SAVE RESULT ======================

if save_epoch_result

    out_name = sprintf('%s_EpochWiseResponse_%s_Epoch%.0fto%.0fms.mat', ...
        base_name, spike_source, epoch_edges_ms(1), epoch_edges_ms(end));

    save(out_name, ...
        'EpochMetricTable', ...
        'EpochMat', ...
        'EpochSEMMat', ...
        'EpochProfile', ...
        'EpochFraction', ...
        'CondList', ...
        'ShortLabels', ...
        'channels_to_analyse', ...
        'Amps_selected', ...
        'epoch_edges_ms', ...
        'epoch_centres_ms', ...
        'baseline_window_ms', ...
        'do_baseline_correction', ...
        'set_negative_bins_to_zero_for_metrics', ...
        'spike_source', ...
        'pool_seq_sets', ...
        'pool_autosim_sets', ...
        '-v7.3');

    fprintf('\nSaved epoch-wise response result to:\n%s\n', fullfile(data_folder, out_name));
end

fprintf('\nFinished epoch-wise response analysis.\n');

%% ====================== LOCAL HELPER FUNCTIONS ======================

function setLabel = buildSetLabelFromTrials(trial_list, StimParams_data, simultaneous_stim, E_MAP, isAutoSim)
    % Build stimulation label from a list of trials.
    %
    % Output examples:
    %   Seq:     Ch18→Ch22
    %   AutoSim: Ch18+Ch22

    stimNames_all = {};

    for ii = 1:numel(trial_list)

        tr = trial_list(ii);

        rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);

        if max(rr) > size(StimParams_data,1)
            continue;
        end

        amp_vec = cell2mat(StimParams_data(rr,16));
        activeRows = amp_vec > 0;

        stimNames_this = StimParams_data(rr(activeRows),1);

        stimNames_all = [stimNames_all; stimNames_this(:)]; %#ok<AGROW>
    end

    if isempty(stimNames_all)
        setLabel = 'NoActiveStim';
        return;
    end

    stimNames_all = unique(stimNames_all, 'stable');

    labelParts = cell(1, numel(stimNames_all));

    for i = 1:numel(stimNames_all)

        stimName = stimNames_all{i};

        chNum = convertStimNameUsingEMap(stimName, E_MAP);

        if isnan(chNum)
            labelParts{i} = sprintf('%s', stimName);
        else
            labelParts{i} = sprintf('Ch%d', chNum);
        end
    end

    if isAutoSim
        setLabel = strjoin(labelParts, '+');
    else
        setLabel = strjoin(labelParts, '→');
    end
end

function chNum = convertStimNameUsingEMap(stimName, E_MAP)
    % Convert Intan stimulation label such as 'A-017' to channel number.
    %
    % In your E_MAP format:
    %   E_MAP{1,1} is the array/map name.
    %   E_MAP{2,1} is channel 1.
    %   E_MAP{3,1} is channel 2.
    %
    % Therefore:
    %   if E_MAP{row,1} matches stimName,
    %   channel number = row - 1.

    chNum = NaN;

    if isempty(stimName)
        return;
    end

    if isstring(stimName)
        stimName = char(stimName);
    end

    stimName = strtrim(stimName);

    if iscell(E_MAP)

        for r = 1:size(E_MAP,1)

            thisName = E_MAP{r,1};

            if isstring(thisName)
                thisName = char(thisName);
            end

            if ischar(thisName)

                thisName = strtrim(thisName);

                if strcmp(thisName, stimName)
                    chNum = r - 1;
                    return;
                end
            end
        end
    end

    if isstring(E_MAP)

        hit = find(strcmp(strtrim(E_MAP), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    if ischar(E_MAP)

        nameList = cellstr(E_MAP);
        hit = find(strcmp(strtrim(nameList), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
        warning('Stim channel %s was not found in E_MAP. Falling back to parsed number %d.', ...
            stimName, chNum);
    else
        warning('Stim channel %s was not found in E_MAP and could not be parsed.', stimName);
    end
end