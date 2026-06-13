%% ============================================================
%  Step 2: Build and Plot Shifted Linear Predictions
%  for Single- vs Multi-electrode Pulse-train Data
%
%  Purpose:
%    Visually test whether the measured two-electrode response can be
%    predicted from the sum of single-electrode responses.
%
%  Prediction logic:
%
%    Predicted AutoSim Ch35+Ch39 =
%        Single Ch35 train
%      + Single Ch39 train
%
%    Predicted Seq Ch35→Ch39 =
%        Single Ch35 train
%      + Single Ch39 train shifted by +5 ms
%
%    Predicted Seq Ch39→Ch35 =
%        Single Ch39 train
%      + Single Ch35 train shifted by +5 ms
%
%  Output:
%    For each amplitude:
%      predicted matrix
%      actual matrix
%      residual matrix = actual - predicted
%
%    Raster plots:
%      synthetic predicted raster
%      actual multi-electrode raster
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

single_data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Single2';
multi_data_folder  = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

% Recording channels to analyse.
% These are Depth_s indices, not raw Intan channel numbers.
channels_to_analyse = [33:37,40:42,44:48,50:52,54:55,57,61:64];

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Spike source.
% Use the same processed spike type for both datasets.
spike_source = 'recovery_ssd';

% Time window and binning.
analysis_window_ms = [0 60];
baseline_window_ms = [-60 -10];

% Use 1 ms bins so the 5 ms shift is exact.
bin_ms = 1;

% Optional smoothing after binning.
% 0 = no smoothing.
% 3 = smooth over 3 ms.
smooth_bins = 3;

% Sequential pulse temporal delay.
shift_ms = 5;

% Amplitudes to analyse.
Amps_to_analyse = [5 10];

%% ====================== RASTER SETTINGS ======================
% ADDED SECTION
%
% These raster plots are for visual checking.
% Start with only a few representative channels, otherwise you will get
% many figures.

plot_rasters = true;

% Depth_s channel indices to plot.
% You can change these after checking which channels show clear responses.
raster_channels_to_plot = [50:52];

% Raster plotting window.
raster_window_ms = [0 60];

% Number of trials to show in predicted and actual raster.
n_raster_trials_to_plot = 30;

% Marker size for raster dots.
raster_marker_size = 6;

%% ----- Condition IDs from your inspection table -----

CondID.Single_Ch35 = containers.Map({'5','10'}, [5 6]);
CondID.Single_Ch39 = containers.Map({'5','10'}, [11 12]);

CondID.AutoSim     = containers.Map({'5','10'}, [5 6]);
CondID.Seq_35to39  = containers.Map({'5','10'}, [17 18]);
CondID.Seq_39to35  = containers.Map({'5','10'}, [29 30]);

% Save result?
save_result = false;

%% ====================== PREPARE TIME BINS ======================

bin_edges = analysis_window_ms(1):bin_ms:analysis_window_ms(2);
bin_centres = bin_edges(1:end-1) + bin_ms/2;
nBins = numel(bin_centres);

shift_bins = round(shift_ms / bin_ms);

if abs(shift_bins * bin_ms - shift_ms) > 1e-6
    error('shift_ms must be an integer multiple of bin_ms.');
end

fprintf('\nAnalysis window: %.1f to %.1f ms\n', analysis_window_ms(1), analysis_window_ms(2));
fprintf('Bin size: %.1f ms\n', bin_ms);
fprintf('Shift: %.1f ms = %d bins\n', shift_ms, shift_bins);
fprintf('Smoothing bins: %d\n', smooth_bins);

%% ====================== LOAD BOTH DATASETS ======================

SingleData = loadPulseTrainDataset(single_data_folder, spike_source);
MultiData  = loadPulseTrainDataset(multi_data_folder, spike_source);

%% ====================== CHANNEL MAP ======================
% Depth_s requires read_Intan_RHS2000_file, which reads info.rhs from pwd.
% Therefore, we need to temporarily cd into the dataset folder containing
% info.rhs before calling Depth_s.

channels_to_analyse = channels_to_analyse(:)';
nAnalyseCh = numel(channels_to_analyse);

cur_dir = pwd;

try
    cd(multi_data_folder);
    d = Depth_s(Electrode_Type);
    cd(cur_dir);

    fprintf('\nUsing Depth_s mapping from info.rhs in:\n%s\n', multi_data_folder);

catch ME
    cd(cur_dir);

    warning('Depth_s failed even after changing directory.');
    warning('%s', ME.message);

    % Fallback for quick check only.
    % This assumes sp_use{1} corresponds to recording channel 1, etc.
    d = 1:max(channels_to_analyse);

    fprintf('\nFalling back to direct spike-channel mapping: d(i) = i\n');
end

fprintf('Number of selected recording channels: %d\n', nAnalyseCh);

%% ====================== MAIN ANALYSIS ======================

PredictionResult = struct();

for aa = 1:numel(Amps_to_analyse)

    amp_val = Amps_to_analyse(aa);
    amp_key = sprintf('%g', amp_val);

    fprintf('\n============================================================\n');
    fprintf('Amplitude: %.1f uA\n', amp_val);
    fprintf('============================================================\n');

    %% ----- Get condition IDs -----
    cond_single_ch35 = CondID.Single_Ch35(amp_key);
    cond_single_ch39 = CondID.Single_Ch39(amp_key);

    cond_autosim     = CondID.AutoSim(amp_key);
    cond_seq_35to39  = CondID.Seq_35to39(amp_key);
    cond_seq_39to35  = CondID.Seq_39to35(amp_key);

    fprintf('Single Ch35 condition: %d\n', cond_single_ch35);
    fprintf('Single Ch39 condition: %d\n', cond_single_ch39);
    fprintf('AutoSim condition    : %d\n', cond_autosim);
    fprintf('Seq Ch35→Ch39 cond   : %d\n', cond_seq_35to39);
    fprintf('Seq Ch39→Ch35 cond   : %d\n', cond_seq_39to35);

    %% ----- Build single-electrode response matrices -----
    Single_Ch35 = buildConditionMatrix( ...
        SingleData, cond_single_ch35, channels_to_analyse, d, ...
        bin_edges, baseline_window_ms, smooth_bins);

    Single_Ch39 = buildConditionMatrix( ...
        SingleData, cond_single_ch39, channels_to_analyse, d, ...
        bin_edges, baseline_window_ms, smooth_bins);

    %% ----- Build actual multi-electrode response matrices -----
    Actual_AutoSim = buildConditionMatrix( ...
        MultiData, cond_autosim, channels_to_analyse, d, ...
        bin_edges, baseline_window_ms, smooth_bins);

    Actual_Seq_35to39 = buildConditionMatrix( ...
        MultiData, cond_seq_35to39, channels_to_analyse, d, ...
        bin_edges, baseline_window_ms, smooth_bins);

    Actual_Seq_39to35 = buildConditionMatrix( ...
        MultiData, cond_seq_39to35, channels_to_analyse, d, ...
        bin_edges, baseline_window_ms, smooth_bins);

    %% ----- Build predictions -----
    Pred_AutoSim = Single_Ch35 + Single_Ch39;

    Pred_Seq_35to39 = Single_Ch35 + shiftMatrixRight(Single_Ch39, shift_bins);

    Pred_Seq_39to35 = Single_Ch39 + shiftMatrixRight(Single_Ch35, shift_bins);

    %% ----- Residuals -----
    Resid_AutoSim    = Actual_AutoSim    - Pred_AutoSim;
    Resid_Seq_35to39 = Actual_Seq_35to39 - Pred_Seq_35to39;
    Resid_Seq_39to35 = Actual_Seq_39to35 - Pred_Seq_39to35;

    %% ----- Store -----
    PredictionResult(aa).Amplitude_uA = amp_val;

    PredictionResult(aa).Single_Ch35 = Single_Ch35;
    PredictionResult(aa).Single_Ch39 = Single_Ch39;

    PredictionResult(aa).Pred_AutoSim = Pred_AutoSim;
    PredictionResult(aa).Actual_AutoSim = Actual_AutoSim;
    PredictionResult(aa).Resid_AutoSim = Resid_AutoSim;

    PredictionResult(aa).Pred_Seq_35to39 = Pred_Seq_35to39;
    PredictionResult(aa).Actual_Seq_35to39 = Actual_Seq_35to39;
    PredictionResult(aa).Resid_Seq_35to39 = Resid_Seq_35to39;

    PredictionResult(aa).Pred_Seq_39to35 = Pred_Seq_39to35;
    PredictionResult(aa).Actual_Seq_39to35 = Actual_Seq_39to35;
    PredictionResult(aa).Resid_Seq_39to35 = Resid_Seq_39to35;

    %% ----- Print simple total response check -----
    fprintf('\nTotal response check, summed across selected channels and time:\n');
    fprintf('Single Ch35       : %.3f\n', sum(Single_Ch35(:), 'omitnan'));
    fprintf('Single Ch39       : %.3f\n', sum(Single_Ch39(:), 'omitnan'));
    fprintf('Pred AutoSim      : %.3f\n', sum(Pred_AutoSim(:), 'omitnan'));
    fprintf('Actual AutoSim    : %.3f\n', sum(Actual_AutoSim(:), 'omitnan'));
    fprintf('Pred Seq 35→39    : %.3f\n', sum(Pred_Seq_35to39(:), 'omitnan'));
    fprintf('Actual Seq 35→39  : %.3f\n', sum(Actual_Seq_35to39(:), 'omitnan'));
    fprintf('Pred Seq 39→35    : %.3f\n', sum(Pred_Seq_39to35(:), 'omitnan'));
    fprintf('Actual Seq 39→35  : %.3f\n', sum(Actual_Seq_39to35(:), 'omitnan'));

    %% ----- Plot spatiotemporal matrices -----
    plotPredActualResidual( ...
        Pred_AutoSim, Actual_AutoSim, Resid_AutoSim, ...
        channels_to_analyse, bin_centres, ...
        sprintf('AutoSim Ch35+Ch39 | %.1f uA', amp_val));

    plotPredActualResidual( ...
        Pred_Seq_35to39, Actual_Seq_35to39, Resid_Seq_35to39, ...
        channels_to_analyse, bin_centres, ...
        sprintf('Seq Ch35→Ch39 | %.1f uA', amp_val));

    plotPredActualResidual( ...
        Pred_Seq_39to35, Actual_Seq_39to35, Resid_Seq_39to35, ...
        channels_to_analyse, bin_centres, ...
        sprintf('Seq Ch39→Ch35 | %.1f uA', amp_val));

    %% ----- Plot summed temporal profiles -----
    plotTemporalPrediction( ...
        Pred_AutoSim, Actual_AutoSim, ...
        bin_centres, ...
        sprintf('AutoSim Ch35+Ch39 | %.1f uA', amp_val));

    plotTemporalPrediction( ...
        Pred_Seq_35to39, Actual_Seq_35to39, ...
        bin_centres, ...
        sprintf('Seq Ch35→Ch39 | %.1f uA', amp_val));

    plotTemporalPrediction( ...
        Pred_Seq_39to35, Actual_Seq_39to35, ...
        bin_centres, ...
        sprintf('Seq Ch39→Ch35 | %.1f uA', amp_val));

    %% ----- ADDED: Plot predicted vs actual rasters -----
    if plot_rasters

        for rr_ch = 1:numel(raster_channels_to_plot)

            depth_ch = raster_channels_to_plot(rr_ch);

            % AutoSim prediction:
            % Single Ch35 + Single Ch39, no temporal shift.
            plotPredictedVsActualRaster( ...
                SingleData, MultiData, ...
                cond_single_ch35, cond_single_ch39, cond_autosim, ...
                depth_ch, d, ...
                0, ...
                raster_window_ms, ...
                n_raster_trials_to_plot, ...
                raster_marker_size, ...
                sprintf('AutoSim Ch35+Ch39 | %.1f uA | Depth_s Ch%d', amp_val, depth_ch));

            % Sequential Ch35→Ch39 prediction:
            % Single Ch35 + Single Ch39 shifted by +5 ms.
            plotPredictedVsActualRaster( ...
                SingleData, MultiData, ...
                cond_single_ch35, cond_single_ch39, cond_seq_35to39, ...
                depth_ch, d, ...
                shift_ms, ...
                raster_window_ms, ...
                n_raster_trials_to_plot, ...
                raster_marker_size, ...
                sprintf('Seq Ch35→Ch39 | %.1f uA | Depth_s Ch%d', amp_val, depth_ch));

            % Sequential Ch39→Ch35 prediction:
            % Single Ch39 + Single Ch35 shifted by +5 ms.
            plotPredictedVsActualRaster( ...
                SingleData, MultiData, ...
                cond_single_ch39, cond_single_ch35, cond_seq_39to35, ...
                depth_ch, d, ...
                shift_ms, ...
                raster_window_ms, ...
                n_raster_trials_to_plot, ...
                raster_marker_size, ...
                sprintf('Seq Ch39→Ch35 | %.1f uA | Depth_s Ch%d', amp_val, depth_ch));
        end
    end
end

%% ====================== SAVE RESULT ======================

if save_result

    out_name = sprintf('Step2_ShiftedLinearPrediction_Bin%.0fms_Shift%.0fms.mat', ...
        bin_ms, shift_ms);

    save(fullfile(multi_data_folder, out_name), ...
        'PredictionResult', ...
        'single_data_folder', ...
        'multi_data_folder', ...
        'channels_to_analyse', ...
        'Electrode_Type', ...
        'spike_source', ...
        'analysis_window_ms', ...
        'baseline_window_ms', ...
        'bin_ms', ...
        'bin_edges', ...
        'bin_centres', ...
        'smooth_bins', ...
        'shift_ms', ...
        'shift_bins', ...
        'CondID', ...
        'plot_rasters', ...
        'raster_channels_to_plot', ...
        'raster_window_ms', ...
        'n_raster_trials_to_plot', ...
        '-v7.3');

    fprintf('\nSaved prediction result to:\n%s\n', fullfile(multi_data_folder, out_name));
end

fprintf('\nFinished Step 2 shifted linear prediction plotting.\n');

%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================

function Data = loadPulseTrainDataset(data_folder, spike_source)

    fprintf('\nLoading dataset:\n%s\n', data_folder);

    if ~isfolder(data_folder)
        error('Folder does not exist:\n%s', data_folder);
    end

    %% ----- Load spikes -----
    switch lower(spike_source)

        case 'recovery_ssd'
            sp_files = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));
            assert(~isempty(sp_files), 'No *sp_xia_PrefixRecovery_SSD.mat found in:\n%s', data_folder);

            sp_filename = fullfile(data_folder, sp_files(1).name);
            S_sp = load(sp_filename);

            if isfield(S_sp, 'sp_corr')
                sp_use = S_sp.sp_corr;
            else
                error('sp_corr not found in:\n%s', sp_filename);
            end

        case 'recovery'
            sp_files = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
            sp_files = sp_files(~contains({sp_files.name}, 'SSD'));
            assert(~isempty(sp_files), 'No *sp_xia_PrefixRecovery.mat found in:\n%s', data_folder);

            sp_filename = fullfile(data_folder, sp_files(1).name);
            S_sp = load(sp_filename);

            if isfield(S_sp, 'sp_seq')
                sp_use = S_sp.sp_seq;
            else
                error('sp_seq not found in:\n%s', sp_filename);
            end

        otherwise
            error('Unknown spike_source: %s', spike_source);
    end

    fprintf('Loaded spike file:\n%s\n', sp_filename);

    %% ----- Load triggers -----
    if isempty(dir(fullfile(data_folder, '*.trig.dat')))
        cur_dir = pwd;
        cd(data_folder);
        cleanTrig_sabquick;
        cd(cur_dir);
    end

    cur_dir = pwd;
    cd(data_folder);
    trig = loadTrig(0);
    cd(cur_dir);

    fprintf('Loaded %d triggers.\n', numel(trig));

    %% ----- Load exp datafile -----
    fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
    assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found in:\n%s', data_folder);

    exp_file = fullfile(data_folder, fileDIR(1).name);
    S_exp = load(exp_file);

    fprintf('Loaded exp datafile:\n%s\n', exp_file);

    StimParams        = S_exp.StimParams;
    TrialParams       = S_exp.TrialParams;
    simultaneous_stim = S_exp.simultaneous_stim;
    n_Trials          = S_exp.n_Trials;

    StimParams_data  = StimParams(2:end,:);
    TrialParams_data = TrialParams(2:end,:);

    firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);
    conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
    conditionID_trial = conditionID_trial(:);

    nTrials_use = min(n_Trials, numel(trig));

    Data.folder = data_folder;
    Data.sp_use = sp_use;
    Data.trig = trig;
    Data.FS = 30000;
    Data.n_Trials = n_Trials;
    Data.nTrials_use = nTrials_use;
    Data.StimParams = StimParams;
    Data.TrialParams = TrialParams;
    Data.StimParams_data = StimParams_data;
    Data.TrialParams_data = TrialParams_data;
    Data.simultaneous_stim = simultaneous_stim;
    Data.conditionID_trial = conditionID_trial;

    fprintf('n_Trials: %d | nTriggers: %d | using: %d\n', ...
        n_Trials, numel(trig), nTrials_use);
end

function M = buildConditionMatrix(Data, condID, channels_to_analyse, d, ...
                                  bin_edges, baseline_window_ms, smooth_bins)
    % Build channel x time matrix for one condition.
    %
    % Each value:
    %   mean baseline-corrected spike count per trial per bin

    FS = Data.FS;
    trig = Data.trig;
    sp_use = Data.sp_use;

    nCh = numel(sp_use);
    nAnalyseCh = numel(channels_to_analyse);
    nBins = numel(bin_edges) - 1;

    bin_ms = mean(diff(bin_edges));
    base_dur_ms = baseline_window_ms(2) - baseline_window_ms(1);

    if base_dur_ms <= 0
        error('Invalid baseline window.');
    end

    trial_list = find(Data.conditionID_trial == condID);
    trial_list = trial_list(trial_list <= Data.nTrials_use);

    if isempty(trial_list)
        warning('No trials found for condition ID %d in dataset:\n%s', condID, Data.folder);
        M = NaN(nAnalyseCh, nBins);
        return;
    end

    M = NaN(nAnalyseCh, nBins);

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

            raw_counts = histcounts(rel_t, bin_edges);

            corrected_counts = raw_counts - expected_baseline_per_bin;

            trial_bin_counts(tt,:) = corrected_counts;
        end

        y = mean(trial_bin_counts, 1, 'omitnan');

        if smooth_bins > 1
            y = movmean(y, smooth_bins, 'omitnan');
        end

        M(ich_i,:) = y;
    end
end

function M_shift = shiftMatrixRight(M, shift_bins)
    % Shift channel x time matrix right by shift_bins.
    %
    % Right shift by +5 ms means response appears later.
    %
    % New empty early bins are filled with 0.
    % Late bins shifted beyond the window are discarded.

    M_shift = zeros(size(M));

    if shift_bins == 0
        M_shift = M;
        return;
    end

    if shift_bins > 0
        M_shift(:, shift_bins+1:end) = M(:, 1:end-shift_bins);
    else
        shift_bins = abs(shift_bins);
        M_shift(:, 1:end-shift_bins) = M(:, shift_bins+1:end);
    end
end

function plotPredActualResidual(Pred, Actual, Resid, channels_to_analyse, bin_centres, figTitle)

    figure('Color','w', ...
           'Name', ['PredActualResidual_' figTitle], ...
           'Position', [100 100 1200 420]);

    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    allVals = [Pred(:); Actual(:)];
    climMax = max(abs(allVals), [], 'omitnan');

    if isempty(climMax) || isnan(climMax) || climMax == 0
        climMax = 1;
    end

    residMax = max(abs(Resid(:)), [], 'omitnan');
    if isempty(residMax) || isnan(residMax) || residMax == 0
        residMax = 1;
    end

    %% ----- Predicted -----
    nexttile;
    imagesc(bin_centres, 1:numel(channels_to_analyse), Pred);
    axis xy;
    colorbar;
    caxis([0 climMax]);
    xlabel('Time after trigger (ms)');
    ylabel('Depth_s channel index');
    title('Predicted', 'Interpreter','none');
    yticks(1:numel(channels_to_analyse));
    yticklabels(arrayfun(@num2str, channels_to_analyse, 'UniformOutput', false));

    %% ----- Actual -----
    nexttile;
    imagesc(bin_centres, 1:numel(channels_to_analyse), Actual);
    axis xy;
    colorbar;
    caxis([0 climMax]);
    xlabel('Time after trigger (ms)');
    ylabel('Depth_s channel index');
    title('Actual', 'Interpreter','none');
    yticks(1:numel(channels_to_analyse));
    yticklabels(arrayfun(@num2str, channels_to_analyse, 'UniformOutput', false));

    %% ----- Residual -----
    nexttile;
    imagesc(bin_centres, 1:numel(channels_to_analyse), Resid);
    axis xy;
    colorbar;
    caxis([-residMax residMax]);
    xlabel('Time after trigger (ms)');
    ylabel('Depth_s channel index');
    title('Residual: Actual - Predicted', 'Interpreter','none');
    yticks(1:numel(channels_to_analyse));
    yticklabels(arrayfun(@num2str, channels_to_analyse, 'UniformOutput', false));

    sgtitle(figTitle, 'Interpreter','none');
end

function plotTemporalPrediction(Pred, Actual, bin_centres, figTitle)

    pred_t = sum(Pred, 1, 'omitnan');
    actual_t = sum(Actual, 1, 'omitnan');
    resid_t = actual_t - pred_t;

    figure('Color','w', ...
           'Name', ['TemporalPrediction_' figTitle], ...
           'Position', [150 150 850 480]);

    hold on; box off;

    plot(bin_centres, pred_t, '-', ...
        'LineWidth', 2, ...
        'DisplayName', 'Predicted');

    plot(bin_centres, actual_t, '-', ...
        'LineWidth', 2, ...
        'DisplayName', 'Actual');

    plot(bin_centres, resid_t, '--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Residual');

    yline(0, 'k:');

    xlabel('Time after trigger (ms)');
    ylabel('Summed baseline-corrected spike count');
    title(figTitle, 'Interpreter','none');
    legend('Box','off','Location','best');
    grid on;
end

%% ====================== ADDED LOCAL FUNCTIONS FOR RASTER PLOTS ======================

function plotPredictedVsActualRaster( ...
    SingleData, MultiData, ...
    cond_single_first, cond_single_second, cond_actual_multi, ...
    depth_ch, d, ...
    second_shift_ms, ...
    raster_window_ms, ...
    n_raster_trials_to_plot, ...
    raster_marker_size, ...
    figTitle)
% Plot predicted vs actual raster for one recording channel.
%
% Predicted raster:
%   single first-electrode spikes
%   +
%   single second-electrode spikes shifted by second_shift_ms
%
% For AutoSim:
%   second_shift_ms = 0
%
% For sequential:
%   second_shift_ms = 5
%
% Important:
%   The predicted raster is synthetic. It combines spikes from separate
%   single-electrode trials. It is not a real recorded trial.

    spike_ch = d(depth_ch);

    %% ----- Get single-electrode spike trials -----
    first_trials = getConditionTrialSpikeTimes( ...
        SingleData, cond_single_first, spike_ch, raster_window_ms);

    second_trials = getConditionTrialSpikeTimes( ...
        SingleData, cond_single_second, spike_ch, raster_window_ms);

    %% ----- Get actual multi-electrode spike trials -----
    actual_trials = getConditionTrialSpikeTimes( ...
        MultiData, cond_actual_multi, spike_ch, raster_window_ms);

    if isempty(first_trials) || isempty(second_trials) || isempty(actual_trials)
        warning('Skipping raster because one condition has no trials: %s', figTitle);
        return;
    end

    %% ----- Build predicted trials -----
    nPred = min([numel(first_trials), numel(second_trials), n_raster_trials_to_plot]);
    nAct  = min(numel(actual_trials), n_raster_trials_to_plot);

    predicted_trials = cell(nPred, 1);

    for tr = 1:nPred

        sp1 = first_trials{tr};
        sp2 = second_trials{tr} + second_shift_ms;

        sp_pred = [sp1(:); sp2(:)];
        sp_pred = sp_pred(sp_pred >= raster_window_ms(1) & ...
                          sp_pred <= raster_window_ms(2));
        sp_pred = sort(sp_pred);

        predicted_trials{tr} = sp_pred;
    end

    actual_trials_plot = actual_trials(1:nAct);

    %% ----- Plot -----
    figure('Color','w', ...
           'Name', ['Raster_' figTitle], ...
           'Position', [150 150 950 650]);

    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    %% Predicted raster
    ax1 = nexttile;
    hold(ax1, 'on'); box(ax1, 'off');

    for tr = 1:nPred
        sp = predicted_trials{tr};
        if isempty(sp)
            continue;
        end

        y = tr * ones(size(sp));
        plot(ax1, sp, y, 'k.', 'MarkerSize', raster_marker_size);
    end

    xlim(ax1, raster_window_ms);
    ylim(ax1, [0 nPred+1]);
    ylabel(ax1, 'Predicted trial');
    title(ax1, 'Predicted raster from single-electrode responses', 'Interpreter','none');
    grid(ax1, 'on');

    %% Actual raster
    ax2 = nexttile;
    hold(ax2, 'on'); box(ax2, 'off');

    for tr = 1:nAct
        sp = actual_trials_plot{tr};
        if isempty(sp)
            continue;
        end

        y = tr * ones(size(sp));
        plot(ax2, sp, y, 'k.', 'MarkerSize', raster_marker_size);
    end

    xlim(ax2, raster_window_ms);
    ylim(ax2, [0 nAct+1]);
    xlabel(ax2, 'Time after trigger (ms)');
    ylabel(ax2, 'Actual trial');
    title(ax2, 'Actual multi-electrode raster', 'Interpreter','none');
    grid(ax2, 'on');

    sgtitle(figTitle, 'Interpreter','none');

    %% Add event timing guide lines
    addRasterEventLines(ax1, figTitle);
    addRasterEventLines(ax2, figTitle);
end

function trial_spike_cells = getConditionTrialSpikeTimes(Data, condID, spike_ch, raster_window_ms)
% Return cell array of relative spike times for one condition and one channel.
%
% Output:
%   trial_spike_cells{trial_i} = spike times relative to trigger, in ms

    trial_list = find(Data.conditionID_trial == condID);
    trial_list = trial_list(trial_list <= Data.nTrials_use);

    if isempty(trial_list)
        trial_spike_cells = {};
        return;
    end

    if spike_ch > numel(Data.sp_use) || isempty(Data.sp_use{spike_ch})
        trial_spike_cells = {};
        return;
    end

    spike_times = Data.sp_use{spike_ch}(:,1);
    trig = Data.trig;
    FS = Data.FS;

    trial_spike_cells = cell(numel(trial_list), 1);

    for tt = 1:numel(trial_list)

        tr = trial_list(tt);
        t0_ms = trig(tr) / FS * 1000;

        rel_t = spike_times - t0_ms;

        rel_t = rel_t(rel_t >= raster_window_ms(1) & ...
                      rel_t <= raster_window_ms(2));

        trial_spike_cells{tt} = rel_t(:);
    end
end

function addRasterEventLines(ax, figTitle)
% Add stimulation timing guide lines based on condition title.

    yl = ylim(ax);

    if contains(figTitle, 'AutoSim')

        eventTimes = [0 10 20];

        for i = 1:numel(eventTimes)
            xline(ax, eventTimes(i), '--', ...
                'Color', [0.85 0 0], ...
                'LineWidth', 1.2);
        end

    elseif contains(figTitle, 'Seq')

        eventTimes = [0 5 10 15 20 25];

        for i = 1:numel(eventTimes)

            if mod(i,2) == 1
                lineColor = [0.85 0 0];      % first electrode events
            else
                lineColor = [0 0.2 0.85];    % second electrode events
            end

            xline(ax, eventTimes(i), '--', ...
                'Color', lineColor, ...
                'LineWidth', 1.2);
        end
    end

    ylim(ax, yl);
end