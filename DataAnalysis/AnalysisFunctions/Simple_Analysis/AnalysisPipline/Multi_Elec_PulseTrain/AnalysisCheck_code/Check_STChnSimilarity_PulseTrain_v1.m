%% ============================================================
%  Quick Response Distribution / Pattern Richness Checker
%  for Pulse-Train Recovery Data
%
%  Purpose:
%    Test whether sequential pulse-train stimulation produces a more
%    distributed / less synchronized / richer spatiotemporal response
%    than AutoSim/simultaneous pulse train.
%
%  Main question:
%    AutoSim may produce more spikes, but is it more temporally concentrated?
%    Sequential may produce fewer spikes, but is it more distributed across
%    time and channel-time bins?
%
%  Metrics:
%    1) Total response
%       Sum of baseline-corrected response across channel x time bins.
%
%    2) Peak-to-total ratio
%       max temporal response / total temporal response.
%       Higher value = more synchronized / concentrated response.
%
%    3) Temporal entropy
%       Entropy of response distribution over time.
%       Higher value = response distributed across more time bins.
%
%    4) Temporal 80% width
%       Time range containing the middle 80% of response.
%       Higher value = broader temporal response.
%
%    5) Spatiotemporal entropy
%       Entropy of response distribution across channel x time bins.
%       Higher value = response distributed across more channel-time bins.
%
%    6) Effective number of active bins
%       exp(entropy), interpretable as number of meaningfully active bins.
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

% Analysis window.
%
% For pulse trains, [0 60] is a good first choice.
analysis_window_ms = [0 60];

% Time bin size.
% I recommend checking both:
%   bin_ms = 2;
%   bin_ms = 5;
bin_ms = 5;

% Baseline window.
baseline_window_ms = [-60 -10];

% Apply baseline correction to each channel x time bin?
do_baseline_correction = true;

% Remove negative values before entropy calculation?
set_negative_bins_to_zero_for_metrics = true;

% Optional smoothing along time dimension before metrics.
%
% 0 = no smoothing.
% 3 = smooth each channel's temporal profile with 3 bins.
smooth_bins = 0;

% Final train levels to compare.
seq_train_level = 6;
autosim_train_level = 3;

% Plot these families.
plot_sequential = true;
plot_auto_sim = true;

% Stimulation set IDs to include.
% [] means all detected sets that match selected family/level.
SetIDs_to_plot = [];

% Amplitudes to include.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% Pool sequential sets?
%
% Since your current question is Seq vs AutoSim rather than Seq order,
% default is true.
%
% true:
%   Seq Set2 and Seq Set3 are pooled into one Sequential condition.
%
% false:
%   Seq Set2 and Seq Set3 are shown separately.
pool_seq_sets = true;

% Pool AutoSim sets?
pool_autosim_sets = true;

% Save result?
save_distribution_result = false;

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

fprintf('\nAnalysis window: %.1f to %.1f ms\n', analysis_window_ms(1), analysis_window_ms(2));
fprintf('Bin size: %.1f ms\n', bin_ms);
fprintf('Baseline window: %.1f to %.1f ms\n', baseline_window_ms(1), baseline_window_ms(2));
fprintf('Baseline correction: %d\n', do_baseline_correction);
fprintf('Pool Seq sets: %d\n', pool_seq_sets);
fprintf('Pool AutoSim sets: %d\n', pool_autosim_sets);

%% ====================== BUILD CONDITION LIST ======================
% Each condition is family + set + level.
%
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

fprintf('\n================ Conditions for distribution analysis ================\n');
for c = 1:numel(CondList)
    fprintf('%2d) %s | total trials across selected amps = %d\n', ...
        c, CondList(c).Label, numel(CondList(c).TrialList));
end
fprintf('=====================================================================\n');

ShortLabels = {CondList.ShortLabel};

%% ====================== CHANNEL MAP AND BINNING ======================

d = Depth_s(Electrode_Type);

channels_to_analyse = channels_to_analyse(:)';
nAnalyseCh = numel(channels_to_analyse);

bin_edges = analysis_window_ms(1):bin_ms:analysis_window_ms(2);
bin_centres = bin_edges(1:end-1) + bin_ms/2;
nBins = numel(bin_edges) - 1;

base_dur_ms = baseline_window_ms(2) - baseline_window_ms(1);

if base_dur_ms <= 0
    error('baseline_window_ms must have positive duration.');
end

%% ====================== BUILD CHANNEL x TIME MATRICES ======================
% STMat dimensions:
%   selected channel x time bin x condition x amplitude
%
% Each value:
%   mean baseline-corrected spike count / trial

nCond = numel(CondList);
nAmp = numel(Amps_selected);

STMat = NaN(nAnalyseCh, nBins, nCond, nAmp);

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

            trial_bin_counts = NaN(numel(trial_list), nBins);

            for tt = 1:numel(trial_list)

                tr = trial_list(tt);
                t0_ms = trig(tr) / FS * 1000;

                rel_t = spike_times - t0_ms;

                baseline_count = sum(rel_t >= baseline_window_ms(1) & ...
                                     rel_t <  baseline_window_ms(2));

                baseline_rate_per_ms = baseline_count / base_dur_ms;
                expected_baseline_per_bin = baseline_rate_per_ms * bin_ms;

                raw_bin_counts = histcounts(rel_t, bin_edges);

                if do_baseline_correction
                    corrected_bin_counts = raw_bin_counts - expected_baseline_per_bin;
                else
                    corrected_bin_counts = raw_bin_counts;
                end

                trial_bin_counts(tt,:) = corrected_bin_counts;
            end

            channel_time_vector = mean(trial_bin_counts, 1, 'omitnan');

            if smooth_bins > 1
                channel_time_vector = movmean(channel_time_vector, smooth_bins, 'omitnan');
            end

            STMat(ich_i,:,c,aa) = channel_time_vector;
        end
    end
end

%% ====================== CALCULATE DISTRIBUTION METRICS ======================
% For each condition and amplitude:
%
% Temporal response:
%   sum over channels -> response over time
%
% Spatiotemporal response:
%   full channel x time matrix
%
% Entropy calculations use non-negative response values.

MetricRows = {};
row_i = 0;

MetricStruct = struct();

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    for c = 1:nCond

        M_raw = squeeze(STMat(:,:,c,aa));

        if all(isnan(M_raw(:)))
            continue;
        end

        % Replace NaN with 0 for metric calculation.
        M = M_raw;
        M(isnan(M)) = 0;

        % For entropy/distribution metrics, response weights must be non-negative.
        M_metric = M;

        if set_negative_bins_to_zero_for_metrics
            M_metric(M_metric < 0) = 0;
        end

        % Temporal response vector:
        % sum across channels for each time bin.
        temporal_vec = sum(M_metric, 1);

        % Spatiotemporal vector:
        % all channel x time bins.
        st_vec = M_metric(:);

        total_response = sum(st_vec);

        if total_response <= 0
            peak_to_total = NaN;
            temporal_entropy = NaN;
            temporal_entropy_norm = NaN;
            temporal_effective_bins = NaN;
            temporal_width80_ms = NaN;
            st_entropy = NaN;
            st_entropy_norm = NaN;
            st_effective_bins = NaN;
            st_participation_ratio = NaN;
        else

            %% ----- Peak-to-total ratio -----
            peak_to_total = max(temporal_vec) / total_response;

            %% ----- Temporal entropy -----
            p_t = temporal_vec ./ sum(temporal_vec);
            p_t = p_t(p_t > 0);

            temporal_entropy = -sum(p_t .* log(p_t));
            temporal_entropy_norm = temporal_entropy / log(nBins);
            temporal_effective_bins = exp(temporal_entropy);

            %% ----- Temporal 80% width -----
            cdf_t = cumsum(temporal_vec) ./ sum(temporal_vec);

            idx10 = find(cdf_t >= 0.10, 1, 'first');
            idx90 = find(cdf_t >= 0.90, 1, 'first');

            if isempty(idx10) || isempty(idx90)
                temporal_width80_ms = NaN;
            else
                temporal_width80_ms = bin_centres(idx90) - bin_centres(idx10);
            end

            %% ----- Spatiotemporal entropy -----
            p_st = st_vec ./ sum(st_vec);
            p_st = p_st(p_st > 0);

            nSTBins = numel(st_vec);

            st_entropy = -sum(p_st .* log(p_st));
            st_entropy_norm = st_entropy / log(nSTBins);
            st_effective_bins = exp(st_entropy);

            %% ----- Participation ratio -----
            % Another effective-bin metric:
            %   1 / sum(p^2)
            % Larger value = more distributed response.
            st_participation_ratio = 1 / sum(p_st.^2);
        end

        row_i = row_i + 1;

        MetricRows(row_i,1:15) = { ...
            CondList(c).Family, ...
            CondList(c).SetID, ...
            CondList(c).Level, ...
            CondList(c).Label, ...
            amp_val, ...
            total_response, ...
            peak_to_total, ...
            temporal_entropy, ...
            temporal_entropy_norm, ...
            temporal_effective_bins, ...
            temporal_width80_ms, ...
            st_entropy, ...
            st_entropy_norm, ...
            st_effective_bins, ...
            st_participation_ratio};

        MetricStruct(aa,c).Amp = amp_val;
        MetricStruct(aa,c).ConditionLabel = CondList(c).Label;
        MetricStruct(aa,c).TemporalVector = temporal_vec;
        MetricStruct(aa,c).STMatrix = M_metric;
    end
end

MetricTable = cell2table(MetricRows, ...
    'VariableNames', { ...
    'Family', ...
    'SetID', ...
    'TrainLevel', ...
    'ConditionLabel', ...
    'Amplitude_uA', ...
    'TotalResponse', ...
    'PeakToTotalRatio', ...
    'TemporalEntropy', ...
    'TemporalEntropyNorm', ...
    'TemporalEffectiveBins', ...
    'TemporalWidth80_ms', ...
    'SpatiotemporalEntropy', ...
    'SpatiotemporalEntropyNorm', ...
    'SpatiotemporalEffectiveBins', ...
    'SpatiotemporalParticipationRatio'});

fprintf('\n================ Response Distribution / Pattern Richness Metrics ================\n');
disp(MetricTable);
fprintf('================================================================================\n');

%% ====================== FIGURE 1: TEMPORAL DISTRIBUTION ======================
% Shows whether AutoSim is more peak-concentrated and Seq is more distributed.

cond_colors = lines(nCond);

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    figure('Color','w', ...
           'Name', sprintf('Distribution_TemporalProfile_Amp%g', amp_val), ...
           'Position', [100 100 900 520]);

    hold on; box off;

    for c = 1:nCond

        M_raw = squeeze(STMat(:,:,c,aa));
        M_raw(isnan(M_raw)) = 0;

        if set_negative_bins_to_zero_for_metrics
            M_raw(M_raw < 0) = 0;
        end

        temporal_vec = sum(M_raw, 1);

        if all(temporal_vec == 0)
            continue;
        end

        plot(bin_centres, temporal_vec, '-o', ...
            'Color', cond_colors(c,:), ...
            'LineWidth', 2, ...
            'MarkerSize', 5, ...
            'DisplayName', CondList(c).ShortLabel);
    end

    xlabel('Time after trigger (ms)');
    ylabel('Summed baseline-corrected spike count');

    title(sprintf('Temporal distribution of response | %.1f \\muA | bin %.1f ms', ...
        amp_val, bin_ms), ...
        'Interpreter','none');

    legend('Box','off','Location','best');
    grid on;
    xlim(analysis_window_ms);
end

%% ====================== FIGURE 2: SUMMARY METRIC BAR PLOTS ======================
% One figure per amplitude.
%
% Metrics:
%   Total response
%   Peak-to-total ratio
%   Temporal entropy
%   Temporal 80% width
%   Spatiotemporal entropy
%   Effective channel-time bins

metricNames = { ...
    'TotalResponse', ...
    'PeakToTotalRatio', ...
    'TemporalEntropyNorm', ...
    'TemporalWidth80_ms', ...
    'SpatiotemporalEntropyNorm', ...
    'SpatiotemporalEffectiveBins'};

metricTitles = { ...
    'Total response', ...
    'Peak-to-total ratio', ...
    'Normalized temporal entropy', ...
    'Temporal 80% width (ms)', ...
    'Normalized spatiotemporal entropy', ...
    'Effective channel-time bins'};

for aa = 1:nAmp

    amp_val = Amps_selected(aa);

    figure('Color','w', ...
           'Name', sprintf('Distribution_MetricSummary_Amp%g', amp_val), ...
           'Position', [120 80 1200 700]);

    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for mm = 1:numel(metricNames)

        ax = nexttile;
        hold(ax,'on'); box(ax,'off');

        vals = NaN(1,nCond);

        for c = 1:nCond

            rowMatch = strcmp(MetricTable.ConditionLabel, CondList(c).Label) & ...
                       MetricTable.Amplitude_uA == amp_val;

            if any(rowMatch)
                vals(c) = MetricTable{find(rowMatch,1,'first'), metricNames{mm}};
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

    sgtitle(sprintf('Response distribution metrics | %.1f \\muA', amp_val), ...
        'Interpreter','none');
end

%% ====================== FIGURE 3: SEQ VS AUTOSIM DIRECT COMPARISON ======================
% This figure is most useful when pool_seq_sets = true and pool_autosim_sets = true.
% It directly compares Sequential and AutoSim across amplitudes.

if nCond == 2

    figure('Color','w', ...
           'Name', 'Distribution_SeqVsAutoSim_AcrossAmps', ...
           'Position', [150 100 1100 700]);

    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for mm = 1:numel(metricNames)

        ax = nexttile;
        hold(ax,'on'); box(ax,'off');

        for c = 1:nCond

            vals = NaN(1,nAmp);

            for aa = 1:nAmp

                amp_val = Amps_selected(aa);

                rowMatch = strcmp(MetricTable.ConditionLabel, CondList(c).Label) & ...
                           MetricTable.Amplitude_uA == amp_val;

                if any(rowMatch)
                    vals(aa) = MetricTable{find(rowMatch,1,'first'), metricNames{mm}};
                end
            end

            plot(ax, Amps_selected, vals, '-o', ...
                'Color', cond_colors(c,:), ...
                'LineWidth', 2, ...
                'MarkerSize', 6, ...
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

    sgtitle('Sequential vs AutoSim response distribution across amplitudes', ...
        'Interpreter','none');
end

%% ====================== SAVE RESULT ======================

if save_distribution_result

    out_name = sprintf('%s_ResponseDistribution_%s_Bin%.0fms_Resp%.0fto%.0fms.mat', ...
        base_name, spike_source, bin_ms, analysis_window_ms(1), analysis_window_ms(2));

    save(out_name, ...
        'MetricTable', ...
        'MetricStruct', ...
        'STMat', ...
        'CondList', ...
        'ShortLabels', ...
        'channels_to_analyse', ...
        'Amps_selected', ...
        'analysis_window_ms', ...
        'bin_edges', ...
        'bin_centres', ...
        'bin_ms', ...
        'baseline_window_ms', ...
        'do_baseline_correction', ...
        'set_negative_bins_to_zero_for_metrics', ...
        'smooth_bins', ...
        'spike_source', ...
        'pool_seq_sets', ...
        'pool_autosim_sets', ...
        '-v7.3');

    fprintf('\nSaved response-distribution result to:\n%s\n', fullfile(data_folder, out_name));
end

fprintf('\nFinished response distribution / pattern richness analysis.\n');

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