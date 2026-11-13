%% ============================================================
%  Step 2A: Generic Raster + PSTH Visual Check
%  for Shifted Linear Prediction
%
%  Purpose:
%    Visually check whether two-electrode responses can be predicted from
%    single-electrode responses.
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

single_data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Single2';
multi_data_folder  = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

% Recording channels to plot.
% These are Depth_s indices, not raw Intan channel numbers.
% channels_to_plot = [17:27 33:48 50:52 54:55 57 61:63];
channels_to_plot = [61:63];

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Spike source.
spike_source = 'recovery_ssd';

% Time window and binning.
analysis_window_ms = [-10 60];
baseline_window_ms = [-60 -10];

% Use 1 ms bins so the ISI shift is exact.
bin_ms = 1;

% Smooth PSTH after binning.
% 0 = no smoothing.
% 3 = 3 ms smoothing.
smooth_bins = 3;

% Baseline correction.
do_baseline_correction = true;

% Sequential pulse temporal delay.
shift_ms = 5;

% Amplitudes to plot.
Amps_to_plot = [10];

% Plot SEM shaded region?
%
% SEM = standard deviation across trials / sqrt(number of trials)
%
% The line is mean spike count per trial per bin.
% The shaded region is SEM across trials.
%
% For visual checking, false is often cleaner.
plot_sem = false;

% Raster settings.
n_raster_trials_to_plot = 30;
raster_marker_size = 6;

% Save figures?
save_figures = false;
fig_out_folder = fullfile(multi_data_folder, 'Step2A_Generic_RasterPSTH_LinearPrediction_Check');

%% ====================== CONDITION IDS ======================
CondID.Single_First     = containers.Map({'5','10'}, [5 6]);
CondID.Single_Second    = containers.Map({'5','10'}, [11 12]);

CondID.AutoSim          = containers.Map({'5','10'}, [5 6]);
CondID.Seq_FirstToSecond = containers.Map({'5','10'}, [17 18]);
CondID.Seq_SecondToFirst = containers.Map({'5','10'}, [29 30]);

%% ====================== COLOUR SETTINGS ======================

Color.First    = [0 0.4470 0.7410];        % blue
Color.Second   = [0.8500 0.3250 0.0980];   % orange
Color.Pred     = [0.4940 0.1840 0.5560];   % purple
Color.Actual   = [0.4660 0.6740 0.1880];   % green

Color.EventAuto   = [0.85 0 0];            % red
Color.EventFirst  = [0.85 0 0];            % red
Color.EventSecond = [0 0.2 0.85];          % dark blue
Color.EventGrey   = [0.45 0.45 0.45];      % grey

%% ====================== PREPARE ======================

if save_figures && ~isfolder(fig_out_folder)
    mkdir(fig_out_folder);
end

bin_edges = analysis_window_ms(1):bin_ms:analysis_window_ms(2);
bin_centres = bin_edges(1:end-1) + bin_ms/2;

shift_bins = round(shift_ms / bin_ms);

if abs(shift_bins * bin_ms - shift_ms) > 1e-6
    error('shift_ms must be an integer multiple of bin_ms.');
end

fprintf('\nGeneric raster + PSTH visual check settings:\n');
fprintf('Analysis window: %.1f to %.1f ms\n', analysis_window_ms(1), analysis_window_ms(2));
fprintf('Bin size: %.1f ms\n', bin_ms);
fprintf('Shift: %.1f ms = %d bins\n', shift_ms, shift_bins);
fprintf('Smoothing bins: %d\n', smooth_bins);
fprintf('Baseline correction: %d\n', do_baseline_correction);
fprintf('Plot SEM: %d\n', plot_sem);

%% ====================== LOAD DATASETS ======================

SingleData = loadPulseTrainDataset(single_data_folder, spike_source);
MultiData  = loadPulseTrainDataset(multi_data_folder, spike_source);

%% ====================== CHANNEL MAP ======================

channels_to_plot = channels_to_plot(:)';

cur_dir = pwd;

try
    cd(multi_data_folder);
    d = Depth_s(Electrode_Type);
    cd(cur_dir);

    fprintf('\nUsing Depth_s mapping from info.rhs in:\n%s\n', multi_data_folder);

catch ME
    cd(cur_dir);

    warning('Depth_s failed. Falling back to direct mapping.');
    warning('%s', ME.message);

    d = 1:max(channels_to_plot);
end

%% ====================== MAIN PLOTTING ======================

for aa = 1:numel(Amps_to_plot)

    amp_val = Amps_to_plot(aa);
    amp_key = sprintf('%g', amp_val);

    %% ----- Condition IDs -----
    cond_single_first  = CondID.Single_First(amp_key);
    cond_single_second = CondID.Single_Second(amp_key);

    cond_autosim       = CondID.AutoSim(amp_key);
    cond_seq_first_second = CondID.Seq_FirstToSecond(amp_key);
    cond_seq_second_first = CondID.Seq_SecondToFirst(amp_key);

    %% ----- Automatically extract stimulation labels -----
    firstLabel = getConditionStimLabel(SingleData, cond_single_first);
    secondLabel = getConditionStimLabel(SingleData, cond_single_second);

    autoLabel = getConditionStimLabel(MultiData, cond_autosim);
    seqFirstSecondLabel = sprintf('%s→%s', firstLabel, secondLabel);
    seqSecondFirstLabel = sprintf('%s→%s', secondLabel, firstLabel);

    fprintf('\n============================================================\n');
    fprintf('Amplitude %.1f uA\n', amp_val);
    fprintf('Single first cond : %d | label: %s\n', cond_single_first, firstLabel);
    fprintf('Single second cond: %d | label: %s\n', cond_single_second, secondLabel);
    fprintf('AutoSim cond      : %d | label: %s\n', cond_autosim, autoLabel);
    fprintf('Seq first→second  : %d | label: %s\n', cond_seq_first_second, seqFirstSecondLabel);
    fprintf('Seq second→first  : %d | label: %s\n', cond_seq_second_first, seqSecondFirstLabel);

    for cc = 1:numel(channels_to_plot)

        depth_ch = channels_to_plot(cc);
        spike_ch = d(depth_ch);

        fprintf('Plotting Depth_s Ch%d, spike channel %d\n', depth_ch, spike_ch);

        %% ----- Single-electrode PSTHs -----
        [SingleFirst_mean, SingleFirst_sem] = buildConditionPSTH( ...
            SingleData, cond_single_first, spike_ch, ...
            bin_edges, baseline_window_ms, smooth_bins, do_baseline_correction);

        [SingleSecond_mean, SingleSecond_sem] = buildConditionPSTH( ...
            SingleData, cond_single_second, spike_ch, ...
            bin_edges, baseline_window_ms, smooth_bins, do_baseline_correction);

        %% ----- Actual multi-electrode PSTHs -----
        [ActualAuto_mean, ActualAuto_sem] = buildConditionPSTH( ...
            MultiData, cond_autosim, spike_ch, ...
            bin_edges, baseline_window_ms, smooth_bins, do_baseline_correction);

        [ActualSeqFirstSecond_mean, ActualSeqFirstSecond_sem] = buildConditionPSTH( ...
            MultiData, cond_seq_first_second, spike_ch, ...
            bin_edges, baseline_window_ms, smooth_bins, do_baseline_correction);

        [ActualSeqSecondFirst_mean, ActualSeqSecondFirst_sem] = buildConditionPSTH( ...
            MultiData, cond_seq_second_first, spike_ch, ...
            bin_edges, baseline_window_ms, smooth_bins, do_baseline_correction);

        %% ----- Predicted PSTHs -----
        PredAuto_mean = SingleFirst_mean + SingleSecond_mean;
        PredAuto_sem  = sqrt(SingleFirst_sem.^2 + SingleSecond_sem.^2);

        SingleSecond_shift_mean = shiftVectorRight(SingleSecond_mean, shift_bins);
        SingleSecond_shift_sem  = shiftVectorRight(SingleSecond_sem, shift_bins);

        SingleFirst_shift_mean = shiftVectorRight(SingleFirst_mean, shift_bins);
        SingleFirst_shift_sem  = shiftVectorRight(SingleFirst_sem, shift_bins);

        PredSeqFirstSecond_mean = SingleFirst_mean + SingleSecond_shift_mean;
        PredSeqFirstSecond_sem  = sqrt(SingleFirst_sem.^2 + SingleSecond_shift_sem.^2);

        PredSeqSecondFirst_mean = SingleSecond_mean + SingleFirst_shift_mean;
        PredSeqSecondFirst_sem  = sqrt(SingleSecond_sem.^2 + SingleFirst_shift_sem.^2);

        %% ----- Figure 1: AutoSim raster + PSTH check -----
        figTitleAuto = sprintf('AutoSim prediction check | %.1f uA | Depth_s Ch%d | %s', ...
            amp_val, depth_ch, autoLabel);

        plotAutoSimRasterPSTHCheck( ...
            SingleData, MultiData, ...
            cond_single_first, cond_single_second, cond_autosim, ...
            spike_ch, depth_ch, ...
            bin_centres, ...
            SingleFirst_mean, SingleFirst_sem, ...
            SingleSecond_mean, SingleSecond_sem, ...
            PredAuto_mean, PredAuto_sem, ...
            ActualAuto_mean, ActualAuto_sem, ...
            plot_sem, ...
            analysis_window_ms, ...
            n_raster_trials_to_plot, ...
            raster_marker_size, ...
            Color, ...
            firstLabel, secondLabel, ...
            autoLabel, ...
            figTitleAuto);

        if save_figures
            saveas(gcf, fullfile(fig_out_folder, safeFigName(figTitleAuto, 'png')));
        end

        %% ----- Figure 2: Seq first→second raster + PSTH check -----
        figTitleSeq1 = sprintf('Seq %s prediction check | %.1f uA | Depth_s Ch%d', ...
            seqFirstSecondLabel, amp_val, depth_ch);

        plotSeqRasterPSTHCheck( ...
            SingleData, MultiData, ...
            cond_single_first, cond_single_second, cond_seq_first_second, ...
            spike_ch, depth_ch, ...
            shift_ms, ...
            bin_centres, ...
            SingleFirst_mean, SingleFirst_sem, ...
            SingleSecond_mean, SingleSecond_sem, ...
            PredSeqFirstSecond_mean, PredSeqFirstSecond_sem, ...
            ActualSeqFirstSecond_mean, ActualSeqFirstSecond_sem, ...
            plot_sem, ...
            analysis_window_ms, ...
            n_raster_trials_to_plot, ...
            raster_marker_size, ...
            Color, ...
            firstLabel, secondLabel, ...
            figTitleSeq1);

        if save_figures
            saveas(gcf, fullfile(fig_out_folder, safeFigName(figTitleSeq1, 'png')));
        end

        %% ----- Figure 3: Seq second→first raster + PSTH check -----
        figTitleSeq2 = sprintf('Seq %s prediction check | %.1f uA | Depth_s Ch%d', ...
            seqSecondFirstLabel, amp_val, depth_ch);

        plotSeqRasterPSTHCheck( ...
            SingleData, MultiData, ...
            cond_single_second, cond_single_first, cond_seq_second_first, ...
            spike_ch, depth_ch, ...
            shift_ms, ...
            bin_centres, ...
            SingleSecond_mean, SingleSecond_sem, ...
            SingleFirst_mean, SingleFirst_sem, ...
            PredSeqSecondFirst_mean, PredSeqSecondFirst_sem, ...
            ActualSeqSecondFirst_mean, ActualSeqSecondFirst_sem, ...
            plot_sem, ...
            analysis_window_ms, ...
            n_raster_trials_to_plot, ...
            raster_marker_size, ...
            Color, ...
            secondLabel, firstLabel, ...
            figTitleSeq2);

        if save_figures
            saveas(gcf, fullfile(fig_out_folder, safeFigName(figTitleSeq2, 'png')));
        end
    end
end

fprintf('\nFinished generic Step 2A raster + PSTH prediction check.\n');

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

    %% ----- Load E_MAP if available -----
    if isfield(S_exp, 'E_MAP')
        E_MAP = S_exp.E_MAP;
    else
        E_MAP = [];
        warning('E_MAP not found in exp datafile. Stim labels may be less clean.');
    end

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
    Data.E_MAP = E_MAP;

    fprintf('n_Trials: %d | nTriggers: %d | using: %d\n', ...
        n_Trials, numel(trig), nTrials_use);
end

function label = getConditionStimLabel(Data, condID)
% Automatically build stimulation label for a condition.
%
% Uses active rows in StimParams for the first trial of the condition.
% Converts stimulation names to Ch numbers using E_MAP if available.

    trial_list = find(Data.conditionID_trial == condID);
    trial_list = trial_list(trial_list <= Data.nTrials_use);

    if isempty(trial_list)
        label = sprintf('Cond%d', condID);
        return;
    end

    tr = trial_list(1);

    rr = (tr-1)*Data.simultaneous_stim + (1:Data.simultaneous_stim);

    rr = rr(rr <= size(Data.StimParams_data,1));

    if isempty(rr)
        label = sprintf('Cond%d', condID);
        return;
    end

    stimNames = Data.StimParams_data(rr,1);

    try
        ampVec = cell2mat(Data.StimParams_data(rr,16));
    catch
        ampVec = ones(numel(rr),1);
    end

    activeRows = ampVec > 0;

    stimNames = stimNames(activeRows);

    if isempty(stimNames)
        label = sprintf('Cond%d_NoActiveStim', condID);
        return;
    end

    stimNames = unique(stimNames, 'stable');

    labelParts = cell(1, numel(stimNames));

    for i = 1:numel(stimNames)

        stimName = stimNames{i};

        chNum = convertStimNameUsingEMap(stimName, Data.E_MAP);

        if isnan(chNum)
            labelParts{i} = char(string(stimName));
        else
            labelParts{i} = sprintf('Ch%d', chNum);
        end
    end

    if numel(labelParts) == 1
        label = labelParts{1};
    else
        label = strjoin(labelParts, '+');
    end
end

function chNum = convertStimNameUsingEMap(stimName, E_MAP)
% Convert stimulation name to channel number using E_MAP.
%
% Important:
%   Your previous E_MAP logic used:
%     row 2 corresponds to channel 1
%     row n corresponds to channel n-1
%
% Therefore, channel number = row index - 1.

    chNum = NaN;

    if isempty(stimName)
        return;
    end

    if isstring(stimName)
        stimName = char(stimName);
    end

    stimName = strtrim(stimName);

    if isempty(E_MAP)
        chNum = parseNumberFromStimName(stimName);
        return;
    end

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

    elseif isstring(E_MAP)

        E_list = cellstr(E_MAP);

        for r = 1:numel(E_list)
            if strcmp(strtrim(E_list{r}), stimName)
                chNum = r - 1;
                return;
            end
        end

    elseif ischar(E_MAP)

        E_list = cellstr(E_MAP);

        for r = 1:numel(E_list)
            if strcmp(strtrim(E_list{r}), stimName)
                chNum = r - 1;
                return;
            end
        end
    end

    chNum = parseNumberFromStimName(stimName);
end

function chNum = parseNumberFromStimName(stimName)
% Fallback if E_MAP is unavailable.
% Extracts number from string, if possible.

    chNum = NaN;

    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
    end
end

function [psth_mean, psth_sem] = buildConditionPSTH( ...
    Data, condID, spike_ch, bin_edges, baseline_window_ms, smooth_bins, do_baseline_correction)
% Build PSTH for one condition and one recording channel.
%
% Output:
%   psth_mean = mean spike count per trial per bin
%   psth_sem  = SEM across trials per bin

    FS = Data.FS;
    trig = Data.trig;
    sp_use = Data.sp_use;

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
        psth_mean = NaN(1,nBins);
        psth_sem = NaN(1,nBins);
        return;
    end

    if spike_ch > numel(sp_use) || isempty(sp_use{spike_ch})
        warning('Spike channel %d missing or empty in dataset:\n%s', spike_ch, Data.folder);
        psth_mean = NaN(1,nBins);
        psth_sem = NaN(1,nBins);
        return;
    end

    spike_times = sp_use{spike_ch}(:,1);

    trial_counts = NaN(numel(trial_list), nBins);

    for tt = 1:numel(trial_list)

        tr = trial_list(tt);
        t0_ms = trig(tr) / FS * 1000;

        rel_t = spike_times - t0_ms;

        raw_counts = histcounts(rel_t, bin_edges);

        if do_baseline_correction

            baseline_count = sum(rel_t >= baseline_window_ms(1) & ...
                                 rel_t <  baseline_window_ms(2));

            baseline_rate_per_ms = baseline_count / base_dur_ms;
            expected_baseline_per_bin = baseline_rate_per_ms * bin_ms;

            trial_counts(tt,:) = raw_counts - expected_baseline_per_bin;

        else
            trial_counts(tt,:) = raw_counts;
        end
    end

    psth_mean = mean(trial_counts, 1, 'omitnan');
    psth_sem  = std(trial_counts, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(trial_counts),1));

    if smooth_bins > 1
        psth_mean = movmean(psth_mean, smooth_bins, 'omitnan');
        psth_sem  = movmean(psth_sem, smooth_bins, 'omitnan');
    end
end

function trial_spike_cells = getConditionTrialSpikeTimes(Data, condID, spike_ch, raster_window_ms)
% Return relative spike times for one condition and one recording channel.

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

function y_shift = shiftVectorRight(y, shift_bins)
% Shift vector right by shift_bins.
% Empty early bins are filled with 0.
% Late bins shifted beyond the window are discarded.

    y_shift = zeros(size(y));

    if shift_bins == 0
        y_shift = y;
        return;
    end

    if shift_bins > 0
        y_shift(shift_bins+1:end) = y(1:end-shift_bins);
    else
        shift_bins = abs(shift_bins);
        y_shift(1:end-shift_bins) = y(shift_bins+1:end);
    end
end

function shifted_trials = shiftRasterTrials(trial_spike_cells, shift_ms, raster_window_ms)
% Shift all spike times in raster trials by shift_ms.

    shifted_trials = cell(size(trial_spike_cells));

    for tt = 1:numel(trial_spike_cells)

        sp = trial_spike_cells{tt} + shift_ms;

        sp = sp(sp >= raster_window_ms(1) & ...
                sp <= raster_window_ms(2));

        shifted_trials{tt} = sp(:);
    end
end

function combined_trials = combineRasterTrials(first_trials, second_trials, raster_window_ms)
% Combine two raster trial sets trial-by-trial.
%
% predicted trial i =
%   first_trials{i} + second_trials{i}

    nComb = min(numel(first_trials), numel(second_trials));
    combined_trials = cell(nComb,1);

    for tt = 1:nComb

        sp = [first_trials{tt}(:); second_trials{tt}(:)];

        sp = sp(sp >= raster_window_ms(1) & ...
                sp <= raster_window_ms(2));

        combined_trials{tt} = sort(sp);
    end
end

function plotAutoSimRasterPSTHCheck( ...
    SingleData, MultiData, ...
    cond_single_first, cond_single_second, cond_autosim, ...
    spike_ch, depth_ch, ...
    bin_centres, ...
    SingleFirst_mean, SingleFirst_sem, ...
    SingleSecond_mean, SingleSecond_sem, ...
    PredAuto_mean, PredAuto_sem, ...
    ActualAuto_mean, ActualAuto_sem, ...
    plot_sem, ...
    plot_window_ms, ...
    n_raster_trials_to_plot, ...
    raster_marker_size, ...
    Color, ...
    firstLabel, secondLabel, ...
    autoLabel, ...
    figTitle)

    %% ----- Extract raster trials -----
    first_trials = getConditionTrialSpikeTimes( ...
        SingleData, cond_single_first, spike_ch, plot_window_ms);

    second_trials = getConditionTrialSpikeTimes( ...
        SingleData, cond_single_second, spike_ch, plot_window_ms);

    actualAuto_trials = getConditionTrialSpikeTimes( ...
        MultiData, cond_autosim, spike_ch, plot_window_ms);

    %% ----- Create figure -----
    figure('Color','w', ...
           'Name', figTitle, ...
           'Position', [100 80 700 850]);

    tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

    %% Raster: Single first
    ax1 = nexttile;
    plotRasterOnly(ax1, first_trials, n_raster_trials_to_plot, ...
        raster_marker_size, plot_window_ms, ...
        sprintf('Single %s raster', firstLabel), Color.First);
    addEventLines(ax1, [0 10 20], Color.EventGrey);

    %% Raster: Single second
    ax2 = nexttile;
    plotRasterOnly(ax2, second_trials, n_raster_trials_to_plot, ...
        raster_marker_size, plot_window_ms, ...
        sprintf('Single %s raster', secondLabel), Color.Second);
    addEventLines(ax2, [0 10 20], Color.EventGrey);

    %% Raster: Actual AutoSim
    ax3 = nexttile;
    plotRasterOnly(ax3, actualAuto_trials, n_raster_trials_to_plot, ...
        raster_marker_size, plot_window_ms, ...
        sprintf('Actual AutoSim %s raster', autoLabel), Color.Actual);
    addEventLines(ax3, [0 10 20], Color.EventAuto);

    %% PSTH overlay
    ax4 = nexttile;
    hold(ax4,'on'); box(ax4,'off');

    plotPSTHWithSEM(ax4, bin_centres, SingleFirst_mean, SingleFirst_sem, ...
        plot_sem, Color.First, sprintf('Single %s', firstLabel));

    plotPSTHWithSEM(ax4, bin_centres, SingleSecond_mean, SingleSecond_sem, ...
        plot_sem, Color.Second, sprintf('Single %s', secondLabel));

    plotPSTHWithSEM(ax4, bin_centres, PredAuto_mean, PredAuto_sem, ...
        plot_sem, Color.Pred, 'Pred AutoSim');

    plotPSTHWithSEM(ax4, bin_centres, ActualAuto_mean, ActualAuto_sem, ...
        plot_sem, Color.Actual, 'Actual AutoSim');

    addEventLines(ax4, [0 10 20], Color.EventAuto);

    xlabel(ax4, 'Time after trigger (ms)');
    ylabel(ax4, 'Spike count / trial / bin');
    title(ax4, 'PSTH overlay', 'Interpreter','none');
    legend(ax4, 'Box','off','Location','best');
    xlim(ax4, plot_window_ms);

    improveAxes([ax1 ax2 ax3 ax4]);

    sgtitle(sprintf('%s | spike channel %d', figTitle, spike_ch), ...
        'Interpreter','none');
end

function plotSeqRasterPSTHCheck( ...
    SingleData, MultiData, ...
    cond_single_first, cond_single_second, cond_actual_seq, ...
    spike_ch, depth_ch, ...
    shift_ms, ...
    bin_centres, ...
    SingleFirst_mean, SingleFirst_sem, ...
    SingleSecond_mean, SingleSecond_sem, ...
    PredSeq_mean, PredSeq_sem, ...
    ActualSeq_mean, ActualSeq_sem, ...
    plot_sem, ...
    plot_window_ms, ...
    n_raster_trials_to_plot, ...
    raster_marker_size, ...
    Color, ...
    firstLabel, secondLabel, ...
    figTitle)

    %% ----- Extract raster trials -----
    first_trials = getConditionTrialSpikeTimes( ...
        SingleData, cond_single_first, spike_ch, plot_window_ms);

    second_trials_raw = getConditionTrialSpikeTimes( ...
        SingleData, cond_single_second, spike_ch, plot_window_ms);

    actualSeq_trials = getConditionTrialSpikeTimes( ...
        MultiData, cond_actual_seq, spike_ch, plot_window_ms);

    %% ----- Shift second single-electrode raster -----
    second_trials_shifted = shiftRasterTrials(second_trials_raw, shift_ms, plot_window_ms);

    %% ----- Build predicted sequential synthetic raster -----
    predictedSeq_trials = combineRasterTrials( ...
        first_trials, second_trials_shifted, plot_window_ms);

    %% ----- Shift second PSTH for component display -----
    if numel(bin_centres) > 1
        bin_ms_local = mean(diff(bin_centres));
    else
        bin_ms_local = 1;
    end

    shift_bins = round(shift_ms / bin_ms_local);

    SingleSecond_shift_mean = shiftVectorRight(SingleSecond_mean, shift_bins);
    SingleSecond_shift_sem  = shiftVectorRight(SingleSecond_sem, shift_bins);

    %% ----- Create figure -----
    figure('Color','w', ...
           'Name', figTitle, ...
           'Position', [100 50 700 850]);

    tiledlayout(5,1,'TileSpacing','compact','Padding','compact');

    %% Raster: Single first electrode
    ax1 = nexttile;
    plotRasterOnly(ax1, first_trials, n_raster_trials_to_plot, ...
        raster_marker_size, plot_window_ms, ...
        sprintf('Single %s raster', firstLabel), Color.First);
    addEventLines(ax1, [0 10 20], Color.EventGrey);

    %% Raster: Single second electrode shifted
    ax2 = nexttile;
    plotRasterOnly(ax2, second_trials_shifted, n_raster_trials_to_plot, ...
        raster_marker_size, plot_window_ms, ...
        sprintf('Single %s raster shifted +%.1f ms', secondLabel, shift_ms), Color.Second);
    addEventLines(ax2, [5 15 25], Color.EventGrey);

    %% Raster: predicted sequential
    ax3 = nexttile;
    plotRasterOnly(ax3, predictedSeq_trials, n_raster_trials_to_plot, ...
        raster_marker_size, plot_window_ms, ...
        sprintf('Predicted Seq %s→%s raster', firstLabel, secondLabel), Color.Pred);
    addSeqEventLines(ax3, [0 5 10 15 20 25], Color);

    %% Raster: actual sequential
    ax4 = nexttile;
    plotRasterOnly(ax4, actualSeq_trials, n_raster_trials_to_plot, ...
        raster_marker_size, plot_window_ms, ...
        sprintf('Actual Seq %s→%s raster', firstLabel, secondLabel), Color.Actual);
    addSeqEventLines(ax4, [0 5 10 15 20 25], Color);

    %% PSTH overlay
    ax5 = nexttile;
    hold(ax5,'on'); box(ax5,'off');

    % Components as thinner lines.
    plotPSTHWithSEM(ax5, bin_centres, SingleFirst_mean, SingleFirst_sem, ...
        false, Color.First, sprintf('Single %s', firstLabel));

    plotPSTHWithSEM(ax5, bin_centres, SingleSecond_shift_mean, SingleSecond_shift_sem, ...
        false, Color.Second, sprintf('Single %s shifted', secondLabel));

    % Main predicted vs actual comparison.
    plotPSTHWithSEM(ax5, bin_centres, PredSeq_mean, PredSeq_sem, ...
        plot_sem, Color.Pred, 'Pred Seq');

    plotPSTHWithSEM(ax5, bin_centres, ActualSeq_mean, ActualSeq_sem, ...
        plot_sem, Color.Actual, 'Actual Seq');

    addSeqEventLines(ax5, [0 5 10 15 20 25], Color);

    xlabel(ax5, 'Time after trigger (ms)');
    ylabel(ax5, 'Spike count / trial / bin');
    title(ax5, 'PSTH overlay', 'Interpreter','none');
    legend(ax5, 'Box','off','Location','best');
    xlim(ax5, plot_window_ms);

    improveAxes([ax1 ax2 ax3 ax4 ax5]);

    sgtitle(sprintf('%s | spike channel %d', figTitle, spike_ch), ...
        'Interpreter','none');
end

function plotRasterOnly(ax, trial_spike_cells, n_trials_to_plot, marker_size, x_window, titleText, rasterColor)

    hold(ax,'on'); box(ax,'off');

    if isempty(trial_spike_cells)
        title(ax, [titleText ' | no trials'], 'Interpreter','none');
        xlim(ax, x_window);
        return;
    end

    nPlot = min(numel(trial_spike_cells), n_trials_to_plot);

    for tt = 1:nPlot

        sp = trial_spike_cells{tt};

        if isempty(sp)
            continue;
        end

        y = tt .* ones(size(sp));

        plot(ax, sp, y, '.', ...
            'Color', rasterColor, ...
            'MarkerSize', marker_size);
    end

    xlim(ax, x_window);
    ylim(ax, [0 nPlot+1]);
    ylabel(ax, 'Trial');
    title(ax, titleText, 'Interpreter','none');
end

function plotPSTHWithSEM(ax, x, y, sem, plot_sem, lineColor, labelText)

    if plot_sem && ~all(isnan(sem))

        upper = y + sem;
        lower = y - sem;

        fill(ax, [x fliplr(x)], [upper fliplr(lower)], ...
            lineColor, ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end

    plot(ax, x, y, ...
        'Color', lineColor, ...
        'LineWidth', 2, ...
        'DisplayName', labelText);
end

function addEventLines(ax, eventTimes, lineColor)

    yl = ylim(ax);

    for i = 1:numel(eventTimes)
        xline(ax, eventTimes(i), '--', ...
            'Color', lineColor, ...
            'LineWidth', 1.0, ...
            'Alpha', 0.55, ...
            'HandleVisibility', 'off');
    end

    ylim(ax, yl);
end

function addSeqEventLines(ax, eventTimes, Color)

    yl = ylim(ax);

    for i = 1:numel(eventTimes)

        if mod(i,2) == 1
            lineColor = Color.EventFirst;
        else
            lineColor = Color.EventSecond;
        end

        xline(ax, eventTimes(i), '--', ...
            'Color', lineColor, ...
            'LineWidth', 1.0, ...
            'Alpha', 0.55, ...
            'HandleVisibility', 'off');
    end

    ylim(ax, yl);
end

function improveAxes(axList)

    for ii = 1:numel(axList)

        ax = axList(ii);

        ax.FontSize = 10;
        ax.LineWidth = 0.8;
        ax.TickDir = 'out';
        ax.Box = 'off';

        ax.GridAlpha = 0.15;
        ax.MinorGridAlpha = 0.08;
    end
end

function outName = safeFigName(figTitle, ext)
% Make safe file name from figure title.

    outName = regexprep(figTitle, '[^\w\d\-_. ]', '_');
    outName = strrep(outName, ' ', '_');
    outName = [outName '.' ext];
end