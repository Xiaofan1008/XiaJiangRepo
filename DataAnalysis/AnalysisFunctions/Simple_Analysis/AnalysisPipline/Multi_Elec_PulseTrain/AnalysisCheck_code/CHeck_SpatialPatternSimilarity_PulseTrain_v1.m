%% ============================================================
%  Quick Spatial Pattern Similarity Checker for Pulse-Train Data
%
%  Purpose:
%    Test whether sequential pulse-train stimulation produces more distinct
%    spatial response patterns across recording channels.
%
%  Main idea:
%    For each condition, build a spatial response vector:
%
%        [Ch1 count, Ch2 count, Ch3 count, ..., ChN count]
%
%    Then compare these vectors between conditions using:
%       1) Pearson correlation
%       2) Cosine similarity
%       3) Euclidean distance
%
%  Input:
%    Default: *.sp_xia_PrefixRecovery_SSD.mat
%             variable: sp_corr
%
%  Conditions compared:
%    Seq final level     = seq_train_level
%    AutoSim final level = autosim_train_level
%
%  Figures:
%    Figure 1:
%      Spatial response vector across selected recording channels
%
%    Figure 2:
%      Pearson similarity matrix
%
%    Figure 3:
%      Pairwise spatial distance bar plot
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train1';

% Recording channels to analyse.
% These are Depth_s indices, not raw Intan channel numbers.
%
% Examples:
%   channels_to_analyse = 1:16;
%   channels_to_analyse = [1:7 17:28];
%   channels_to_analyse = [1:7 9:12 17:28]; % manually remove noisy channels
channels_to_analyse = [1:7 17:28];

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Spike source:
%   'recovery_ssd' = load filtered recovered spikes, sp_corr
%   'recovery'     = load recovered but unfiltered spikes, sp_seq
spike_source = 'recovery_ssd';

% Response window for spike counting, relative to trigger.
%
% For spatial pattern, I would start with [0 60].
% This includes the train and early post-train response.
response_window_ms = [0 60];

% Baseline window, relative to trigger.
baseline_window_ms = [-60 -10];

% Apply baseline correction?
do_baseline_correction = true;

% Final train levels to compare.
seq_train_level = 6;
autosim_train_level = 3;

% Plot these families.
plot_sequential = true;
plot_auto_sim = true;

% Stimulation set IDs to include.
% [] means all detected sets that match selected family/level.
%
% For your current data, this will likely detect:
%   AutoSim Set1
%   Seq Set2
%   Seq Set3
SetIDs_to_plot = [];

% Amplitudes to include.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% Pool sequential sets?
%
% For this spatial-pattern test, I recommend false.
% We want to know whether Seq Set2 and Seq Set3 produce different patterns.
pool_seq_sets = false;

% Pool AutoSim sets?
%
% Usually false. If there is only one AutoSim set, this does not matter.
pool_autosim_sets = false;

% Normalise each spatial vector before similarity analysis?
%
% false:
%   similarity includes response magnitude and spatial shape.
%
% true:
%   each condition vector is divided by its vector norm, so the comparison
%   focuses more on spatial shape than overall response magnitude.
%
% I recommend false first, then try true as a sensitivity check.
normalise_vectors_before_distance = false;

% Save summary result?
save_similarity_result = false;

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

fprintf('\nResponse window: %.1f to %.1f ms\n', response_window_ms(1), response_window_ms(2));
fprintf('Baseline window: %.1f to %.1f ms\n', baseline_window_ms(1), baseline_window_ms(2));
fprintf('Baseline correction: %d\n', do_baseline_correction);

%% ====================== BUILD CONDITION LIST ======================
% Each condition is family + set + level.
%
% Example:
%   AutoSim Set1 Ch18+Ch22 L3
%   Seq Set2 Ch18→Ch22 L6
%   Seq Set3 Ch22→Ch18 L6

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
            CondList(cond_i).Label = sprintf('Seq L%d pooled sets', seq_train_level);
            CondList(cond_i).ShortLabel = sprintf('Seq pooled L%d', seq_train_level);
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
            CondList(cond_i).ShortLabel = sprintf('Seq S%d L%d', setID, seq_train_level);
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
            CondList(cond_i).Label = sprintf('AutoSim L%d pooled sets', autosim_train_level);
            CondList(cond_i).ShortLabel = sprintf('AutoSim pooled L%d', autosim_train_level);
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
            CondList(cond_i).ShortLabel = sprintf('AutoSim S%d L%d', setID, autosim_train_level);
            CondList(cond_i).AutoFlag = true;
            CondList(cond_i).TrialList = trial_this;
        end
    end
end

if isempty(CondList)
    error('No matching conditions found. Check set IDs, amplitude, and train level settings.');
end

fprintf('\n================ Conditions for spatial pattern comparison ================\n');
for c = 1:numel(CondList)
    fprintf('%2d) %s | total trials across selected amps = %d\n', ...
        c, CondList(c).Label, numel(CondList(c).TrialList));
end
fprintf('=========================================================================\n');

CondLabels = {CondList.Label};
ShortLabels = {CondList.ShortLabel};

%% ====================== CHANNEL MAP ======================

d = Depth_s(Electrode_Type);

channels_to_analyse = channels_to_analyse(:)';
nAnalyseCh = numel(channels_to_analyse);

%% ====================== COUNT BASELINE-CORRECTED SPIKES ======================
% For each:
%   channel × condition × amplitude
%
% This creates the spatial response vector for each condition:
%   MeanMat(:, condition, amplitude)

resp_dur_ms = response_window_ms(2) - response_window_ms(1);
base_dur_ms = baseline_window_ms(2) - baseline_window_ms(1);

if resp_dur_ms <= 0
    error('response_window_ms must have positive duration.');
end

if base_dur_ms <= 0
    error('baseline_window_ms must have positive duration.');
end

MeanMat = NaN(nAnalyseCh, numel(CondList), numel(Amps_selected));
SEMMat  = NaN(nAnalyseCh, numel(CondList), numel(Amps_selected));
NMat    = NaN(nAnalyseCh, numel(CondList), numel(Amps_selected));

SummaryRows = {};
row_i = 0;

for ich_i = 1:nAnalyseCh

    depth_idx = channels_to_analyse(ich_i);
    spike_ch = d(depth_idx);

    if spike_ch > nCh
        warning('Depth index %d maps to spike channel %d, but file has only %d channels. Skipped.', ...
            depth_idx, spike_ch, nCh);
        continue;
    end

    if isempty(sp_use{spike_ch})
        warning('Spike channel %d is empty. Skipped.', spike_ch);
        continue;
    end

    spike_times = sp_use{spike_ch}(:,1);

    for c = 1:numel(CondList)

        for aa = 1:numel(Amps_selected)

            amp_val = Amps_selected(aa);

            trial_list = CondList(c).TrialList;
            trial_list = trial_list(amp_trial(trial_list) == amp_val);
            trial_list = trial_list(trial_list <= nTrials_use);

            if isempty(trial_list)
                continue;
            end

            corr_counts = NaN(numel(trial_list),1);

            for tt = 1:numel(trial_list)

                tr = trial_list(tt);
                t0_ms = trig(tr) / FS * 1000;

                rel_t = spike_times - t0_ms;

                response_count = sum(rel_t >= response_window_ms(1) & ...
                                     rel_t <  response_window_ms(2));

                baseline_count = sum(rel_t >= baseline_window_ms(1) & ...
                                     rel_t <  baseline_window_ms(2));

                expected_baseline_spikes = baseline_count * (resp_dur_ms / base_dur_ms);

                if do_baseline_correction
                    corrected_count = response_count - expected_baseline_spikes;
                else
                    corrected_count = response_count;
                end

                corr_counts(tt) = corrected_count;
            end

            mean_corr = mean(corr_counts, 'omitnan');
            sem_corr  = std(corr_counts, 'omitnan') ./ sqrt(sum(~isnan(corr_counts)));

            MeanMat(ich_i,c,aa) = mean_corr;
            SEMMat(ich_i,c,aa)  = sem_corr;
            NMat(ich_i,c,aa)    = numel(trial_list);

            row_i = row_i + 1;

            SummaryRows(row_i,1:10) = { ...
                depth_idx, ...
                spike_ch, ...
                CondList(c).Family, ...
                CondList(c).SetID, ...
                CondList(c).Level, ...
                CondList(c).Label, ...
                amp_val, ...
                numel(trial_list), ...
                mean_corr, ...
                sem_corr};
        end
    end
end

SummaryTable = cell2table(SummaryRows, ...
    'VariableNames', { ...
    'DepthIndex', ...
    'SpikeChannel', ...
    'Family', ...
    'SetID', ...
    'TrainLevel', ...
    'ConditionLabel', ...
    'Amplitude_uA', ...
    'NTrials', ...
    'MeanBaselineCorrectedCount', ...
    'SEMBaselineCorrectedCount'});

fprintf('\n================ Spatial Pattern Count Summary ================\n');
disp(SummaryTable);
fprintf('===============================================================\n');

%% ====================== SIMILARITY ANALYSIS ======================
% For each amplitude:
%   response vector = MeanMat(:, condition, amplitude)
%
% Similarity matrices:
%   PearsonCorrMat
%   CosineSimMat
%
% Distance matrices:
%   PearsonDistMat = 1 - PearsonCorrMat
%   CosineDistMat  = 1 - CosineSimMat
%   EuclidDistMat

nCond = numel(CondList);
nAmp = numel(Amps_selected);

PearsonCorrMat = NaN(nCond, nCond, nAmp);
CosineSimMat   = NaN(nCond, nCond, nAmp);
PearsonDistMat = NaN(nCond, nCond, nAmp);
CosineDistMat  = NaN(nCond, nCond, nAmp);
EuclidDistMat  = NaN(nCond, nCond, nAmp);

PairwiseRows = {};
pair_i = 0;

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    for i = 1:nCond

        v1 = MeanMat(:,i,aa);

        for j = 1:nCond

            v2 = MeanMat(:,j,aa);

            valid = ~isnan(v1) & ~isnan(v2);

            if sum(valid) < 2
                continue;
            end

            x = v1(valid);
            y = v2(valid);

            if normalise_vectors_before_distance
                x_norm = norm(x);
                y_norm = norm(y);

                if x_norm > 0
                    x = x ./ x_norm;
                end

                if y_norm > 0
                    y = y ./ y_norm;
                end
            end

            % Pearson correlation.
            if std(x) == 0 || std(y) == 0
                r = NaN;
            else
                r = corr(x, y, 'type', 'Pearson', 'rows', 'complete');
            end

            % Cosine similarity.
            if norm(x) == 0 || norm(y) == 0
                cos_sim = NaN;
            else
                cos_sim = dot(x,y) / (norm(x) * norm(y));
            end

            % Euclidean distance.
            euc_dist = sqrt(sum((x - y).^2));

            PearsonCorrMat(i,j,aa) = r;
            CosineSimMat(i,j,aa)   = cos_sim;
            PearsonDistMat(i,j,aa) = 1 - r;
            CosineDistMat(i,j,aa)  = 1 - cos_sim;
            EuclidDistMat(i,j,aa)  = euc_dist;
        end
    end

    % Pairwise table, only upper triangle.
    for i = 1:nCond-1
        for j = i+1:nCond

            pair_i = pair_i + 1;

            PairwiseRows(pair_i,1:9) = { ...
                amp_val, ...
                CondList(i).Label, ...
                CondList(j).Label, ...
                PearsonCorrMat(i,j,aa), ...
                PearsonDistMat(i,j,aa), ...
                CosineSimMat(i,j,aa), ...
                CosineDistMat(i,j,aa), ...
                EuclidDistMat(i,j,aa), ...
                sum(~isnan(MeanMat(:,i,aa)) & ~isnan(MeanMat(:,j,aa)))};
        end
    end
end

PairwiseTable = cell2table(PairwiseRows, ...
    'VariableNames', { ...
    'Amplitude_uA', ...
    'Condition1', ...
    'Condition2', ...
    'PearsonCorrelation', ...
    'PearsonDistance_1minusR', ...
    'CosineSimilarity', ...
    'CosineDistance_1minusCosine', ...
    'EuclideanDistance', ...
    'NChannelsCompared'});

fprintf('\n================ Pairwise Spatial Pattern Similarity ================\n');
disp(PairwiseTable);
fprintf('====================================================================\n');

%% ====================== FIGURE 1: SPATIAL RESPONSE VECTORS ======================
% One figure per amplitude.
%
% X-axis:
%   selected recording channels, compact position
%
% Tick labels:
%   real Depth_s channel index
%
% Y-axis:
%   baseline-corrected spike count / trial
%
% This shows the actual spatial response pattern for each condition.
% cond_colors = lines(nCond);
% xpos = 1:nAnalyseCh;
% 
% for aa = 1:nAmp
% 
%     amp_val = Amps_selected(aa);
% 
%     figure('Color','w', ...
%            'Name', sprintf('SpatialPattern_ResponseVector_Amp%g', amp_val), ...
%            'Position', [100 100 950 560]);
% 
%     hold on; box off;
% 
%     for c = 1:nCond
% 
%         y = squeeze(MeanMat(:,c,aa));
%         e = squeeze(SEMMat(:,c,aa));
% 
%         if all(isnan(y))
%             continue;
%         end
% 
%         errorbar(xpos, y, e, '-o', ...
%             'Color', cond_colors(c,:), ...
%             'LineWidth', 1.8, ...
%             'MarkerSize', 6, ...
%             'DisplayName', CondList(c).Label);
%     end
% 
%     yline(0, 'k:', 'LineWidth', 1);
% 
%     xticks(xpos);
%     xticklabels(arrayfun(@num2str, channels_to_analyse, 'UniformOutput', false));
% 
%     if numel(channels_to_analyse) > 12
%         xtickangle(45);
%     end
% 
%     xlabel('Recording channel index (Depth_s)');
%     ylabel('Baseline-corrected spike count / trial');
% 
%     title(sprintf('Spatial response vectors | %.1f \\muA | Response %.0f–%.0f ms', ...
%         amp_val, response_window_ms(1), response_window_ms(2)), ...
%         'Interpreter','none');
% 
%     legend('Box','off', 'Location','best');
%     xlim([0.5 nAnalyseCh+0.5]);
% end

%% ====================== FIGURE 2: PEARSON SIMILARITY MATRIX ======================
% One figure per amplitude.
%
% Value:
%   Pearson correlation between spatial response vectors.
%
% Interpretation:
%   Higher value  = more similar spatial pattern.
%   Lower value   = more distinct spatial pattern.

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    M = PearsonCorrMat(:,:,aa);

    figure('Color','w', ...
           'Name', sprintf('SpatialPattern_PearsonSimilarity_Amp%g', amp_val), ...
           'Position', [150 150 650 560]);

    imagesc(M);
    axis square;
    colormap(parula);
    colorbar;

    caxis([-1 1]);

    xticks(1:nCond);
    yticks(1:nCond);
    xticklabels(ShortLabels);
    yticklabels(ShortLabels);
    xtickangle(45);

    title(sprintf('Spatial pattern Pearson similarity | %.1f \\muA', amp_val), ...
        'Interpreter','none');

    % Add values inside matrix.
    for i = 1:nCond
        for j = 1:nCond
            if ~isnan(M(i,j))
                text(j, i, sprintf('%.2f', M(i,j)), ...
                    'HorizontalAlignment','center', ...
                    'FontWeight','bold', ...
                    'Color','w');
            end
        end
    end
end

%% ====================== FIGURE 3: PAIRWISE SPATIAL DISTANCE BAR PLOT ======================
% One figure per amplitude.
%
% Metric:
%   Pearson distance = 1 - Pearson correlation
%
% Interpretation:
%   Larger distance = more different spatial pattern.

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    pair_labels = {};
    pair_values = [];

    for i = 1:nCond-1
        for j = i+1:nCond

            pair_labels{end+1} = sprintf('%s vs %s', ShortLabels{i}, ShortLabels{j}); %#ok<SAGROW>
            pair_values(end+1) = PearsonDistMat(i,j,aa); %#ok<SAGROW>
        end
    end

    if isempty(pair_values)
        continue;
    end

    figure('Color','w', ...
           'Name', sprintf('SpatialPattern_PairwiseDistance_Amp%g', amp_val), ...
           'Position', [200 200 850 500]);

    bar(pair_values);
    box off;
    grid on;

    xticks(1:numel(pair_values));
    xticklabels(pair_labels);
    xtickangle(35);

    ylabel('Spatial pattern distance (1 - Pearson r)');
    xlabel('Condition pair');

    title(sprintf('Pairwise spatial pattern distance | %.1f \\muA | Larger = more distinct', ...
        amp_val), ...
        'Interpreter','none');

    yline(0, 'k:', 'LineWidth', 1);

    % Add values above bars.
    for k = 1:numel(pair_values)
        if ~isnan(pair_values(k))
            text(k, pair_values(k), sprintf('%.2f', pair_values(k)), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom');
        end
    end
end

%% ====================== SAVE RESULT ======================

if save_similarity_result

    out_name = sprintf('%s_SpatialPatternSimilarity_%s_Resp%.0fto%.0fms.mat', ...
        base_name, spike_source, response_window_ms(1), response_window_ms(2));

    save(out_name, ...
        'SummaryTable', ...
        'PairwiseTable', ...
        'MeanMat', ...
        'SEMMat', ...
        'NMat', ...
        'PearsonCorrMat', ...
        'CosineSimMat', ...
        'PearsonDistMat', ...
        'CosineDistMat', ...
        'EuclidDistMat', ...
        'CondList', ...
        'CondLabels', ...
        'ShortLabels', ...
        'channels_to_analyse', ...
        'Amps_selected', ...
        'response_window_ms', ...
        'baseline_window_ms', ...
        'do_baseline_correction', ...
        'normalise_vectors_before_distance', ...
        'spike_source', ...
        '-v7.3');

    fprintf('\nSaved spatial-pattern similarity result to:\n%s\n', fullfile(data_folder, out_name));
end

fprintf('\nFinished spatial pattern similarity analysis.\n');

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

    % Fallback only if E_MAP matching fails.
    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
        warning('Stim channel %s was not found in E_MAP. Falling back to parsed number %d.', ...
            stimName, chNum);
    else
        warning('Stim channel %s was not found in E_MAP and could not be parsed.', stimName);
    end
end