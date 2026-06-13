%% ============================================================
%  Quick Spike Count Checker for Pulse-Train Recovery Result
%
%  Purpose:
%    Quickly compare total spike count between final sequential pulse train
%    and final AutoSim / simultaneous pulse train.
%
%  Main question:
%    Does AutoSim really produce higher total spike count than sequential,
%    or was the apparent difference mainly due to PSTH peak shape?
%
%  Input spike file options:
%    recovery_ssd = *.sp_xia_PrefixRecovery_SSD.mat, variable sp_corr
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

% Recording channels to analyse.
% These are Depth_s indices, not raw Intan channel numbers.
channels_to_analyse = [33:37,40:42,44:48,50:52,54:55,57,61:64];

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Spike source:
%   'recovery_ssd' = load filtered recovered spikes, sp_corr
%   'recovery'     = load recovered but unfiltered spikes, sp_seq
spike_source = 'recovery_ssd';

% Response window for spike counting, relative to trigger.
%
% Good windows to test:
%   [0 40]   train + early response
%   [0 50]   train + slightly longer response
%   [0 60]   train + delayed response
%   [25 60]  post-train response only
response_window_ms = [0 40];

% Baseline window, relative to trigger.
baseline_window_ms = [-60 -10];

% Apply baseline correction?
% Corrected count =
%   response_count - expected_baseline_spikes_in_response_window
do_baseline_correction = true;

% Final train levels to compare.
%
% Sequential final level:
%   Level 6 = [0 5 10 15 20 25] ms for 2 electrodes × 3 pulses
%
% AutoSim final level:
%   Level 3 = [0 10 20] ms for 3 simultaneous events
seq_train_level = 6;
autosim_train_level = 3;

% Plot these families.
plot_sequential = true;
plot_auto_sim = true;

% Stimulation set IDs to include.
% [] means all detected sets that match the selected family/level.
SetIDs_to_plot = [];

% Amplitudes to include.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% Pool sequential sets?
%
% false:
%   Seq Set 2 and Seq Set 3 are plotted separately.
%
% true:
%   All sequential sets are pooled into one Seq curve.
%
% For first quick check, false is better because order effects may matter.
pool_seq_sets = false;

% Pool AutoSim sets?
%
% Usually AutoSim may only be Set 1, but this option keeps the code flexible.
pool_autosim_sets = false;

% Save quick summary table?
save_summary_table = false;

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
    error('StimMeta was not found. This quick-check code requires StimMeta.');
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

%% ====================== BUILD TRIAL METADATA FROM StimMeta ======================
% This replaces old prefix metadata columns 26–30.
%
% Important trial-level variables:
%   stimSet_trial
%   trainLevel_trial
%   amp_trial
%   isAutoSim_trial
%   isZero_trial
%   eventTimes_ms_trial
%   eventEnd_ms_trial

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
    % This is safer than relying only on StimMeta.
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
% Each condition is a family + set + level.
%
% For your current comparison:
%   Seq final level     = level 6
%   AutoSim final level = level 3
%
% AutoSim may have a different set ID from Seq.

CondList = struct( ...
    'Family', {}, ...
    'SetID', {}, ...
    'Level', {}, ...
    'Label', {}, ...
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
            CondList(cond_i).AutoFlag = true;
            CondList(cond_i).TrialList = trial_this;
        end
    end
end

if isempty(CondList)
    error('No matching conditions found. Check SetIDs_to_plot, amplitudes, and level settings.');
end

fprintf('\n================ Conditions to compare ================\n');
for c = 1:numel(CondList)
    fprintf('%2d) %s | N total trials across amps = %d\n', ...
        c, CondList(c).Label, numel(CondList(c).TrialList));
end
fprintf('=======================================================\n');

%% ====================== CHANNEL MAP ======================

d = Depth_s(Electrode_Type);

channels_to_analyse = channels_to_analyse(:)';
nAnalyseCh = numel(channels_to_analyse);

%% ====================== COUNT SPIKES ======================
% For each:
%   channel × condition × amplitude
%
% We calculate trial-by-trial:
%   response_count
%   baseline_count
%   expected baseline spikes during response window
%   baseline-corrected spike count

resp_dur_ms = response_window_ms(2) - response_window_ms(1);
base_dur_ms = baseline_window_ms(2) - baseline_window_ms(1);

if resp_dur_ms <= 0
    error('response_window_ms must have positive duration.');
end

if base_dur_ms <= 0
    error('baseline_window_ms must have positive duration.');
end

SummaryRows = {};
row_i = 0;

% Store mean and SEM for plotting.
MeanMat = NaN(nAnalyseCh, numel(CondList), numel(Amps_selected));
SEMMat  = NaN(nAnalyseCh, numel(CondList), numel(Amps_selected));
NMat    = NaN(nAnalyseCh, numel(CondList), numel(Amps_selected));

for ich_i = 1:nAnalyseCh

    ch_depth_index = channels_to_analyse(ich_i);
    ch_spike = d(ch_depth_index);

    if ch_spike > nCh
        warning('Depth index %d maps to spike channel %d, but file has only %d channels. Skipped.', ...
            ch_depth_index, ch_spike, nCh);
        continue;
    end

    if isempty(sp_use{ch_spike})
        warning('Spike channel %d is empty. Skipped.', ch_spike);
        continue;
    end

    spike_times = sp_use{ch_spike}(:,1);

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
            raw_resp_counts = NaN(numel(trial_list),1);
            expected_base_counts = NaN(numel(trial_list),1);

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
                raw_resp_counts(tt) = response_count;
                expected_base_counts(tt) = expected_baseline_spikes;
            end

            mean_corr = mean(corr_counts, 'omitnan');
            sem_corr  = std(corr_counts, 'omitnan') ./ sqrt(sum(~isnan(corr_counts)));

            mean_raw_resp = mean(raw_resp_counts, 'omitnan');
            mean_expected_base = mean(expected_base_counts, 'omitnan');

            MeanMat(ich_i, c, aa) = mean_corr;
            SEMMat(ich_i, c, aa)  = sem_corr;
            NMat(ich_i, c, aa)    = numel(trial_list);

            row_i = row_i + 1;

            SummaryRows(row_i,1:15) = { ...
                ch_depth_index, ...
                ch_spike, ...
                CondList(c).Family, ...
                CondList(c).SetID, ...
                CondList(c).Level, ...
                CondList(c).Label, ...
                amp_val, ...
                numel(trial_list), ...
                response_window_ms(1), ...
                response_window_ms(2), ...
                baseline_window_ms(1), ...
                baseline_window_ms(2), ...
                mean_raw_resp, ...
                mean_expected_base, ...
                mean_corr};
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
    'RespWinStart_ms', ...
    'RespWinEnd_ms', ...
    'BaseWinStart_ms', ...
    'BaseWinEnd_ms', ...
    'MeanRawResponseCount', ...
    'MeanExpectedBaselineCount', ...
    'MeanBaselineCorrectedCount'});

fprintf('\n================ Quick Spike Count Summary ================\n');
disp(SummaryTable);
fprintf('===========================================================\n');

if save_summary_table
    out_summary_name = sprintf('%s_QuickSpikeCount_%s_Resp%.0fto%.0fms.mat', ...
        base_name, spike_source, response_window_ms(1), response_window_ms(2));

    save(out_summary_name, ...
        'SummaryTable', ...
        'MeanMat', ...
        'SEMMat', ...
        'NMat', ...
        'CondList', ...
        'channels_to_analyse', ...
        'Amps_selected', ...
        'response_window_ms', ...
        'baseline_window_ms', ...
        'do_baseline_correction', ...
        'seq_train_level', ...
        'autosim_train_level', ...
        'spike_source', ...
        '-v7.3');

    fprintf('Saved quick spike-count summary to:\n%s\n', fullfile(data_folder, out_summary_name));
end

%% ====================== FIGURE 1: AVERAGED SPIKE COUNT VS AMPLITUDE ======================
% One figure only.
%
% This figure averages across all selected recording channels.
%
% X-axis:
%   amplitude
%
% Y-axis:
%   mean baseline-corrected spike count / trial across selected channels
%
% Error bar:
%   SEM across selected recording channels
%
% Important:
%   This gives a quick population-level view.
%   However, it may hide channel-specific differences, so Figure 2 below
%   still shows each channel separately.

cond_colors = lines(numel(CondList));

MeanAcrossCh = NaN(numel(CondList), numel(Amps_selected));
SEMAcrossCh  = NaN(numel(CondList), numel(Amps_selected));
NAcrossCh    = NaN(numel(CondList), numel(Amps_selected));

for c = 1:numel(CondList)

    for aa = 1:numel(Amps_selected)

        % MeanMat dimensions:
        %   selected channel × condition × amplitude
        y_ch = squeeze(MeanMat(:,c,aa));

        % Remove channels with no valid value.
        y_ch = y_ch(~isnan(y_ch));

        if isempty(y_ch)
            continue;
        end

        MeanAcrossCh(c,aa) = mean(y_ch, 'omitnan');
        SEMAcrossCh(c,aa)  = std(y_ch, 'omitnan') ./ sqrt(numel(y_ch));
        NAcrossCh(c,aa)    = numel(y_ch);
    end
end

figure('Color','w', ...
       'Name', 'QuickCount_AverageAcrossSelectedChannels_vsAmp', ...
       'Position', [100 100 800 560]);

hold on; box off;

for c = 1:numel(CondList)

    y = MeanAcrossCh(c,:);
    e = SEMAcrossCh(c,:);

    if all(isnan(y))
        continue;
    end

    errorbar(Amps_selected, y, e, '-o', ...
        'Color', cond_colors(c,:), ...
        'LineWidth', 2.0, ...
        'MarkerSize', 7, ...
        'DisplayName', CondList(c).Label);
end

yline(0, 'k:', 'LineWidth', 1);

xlabel('Amplitude (\muA)');
ylabel('Baseline-corrected spike count / trial');

title(sprintf(['Average across selected channels | Response %.0f–%.0f ms | ' ...
               'Baseline %.0f–%.0f ms'], ...
      response_window_ms(1), response_window_ms(2), ...
      baseline_window_ms(1), baseline_window_ms(2)), ...
      'Interpreter','none');

legend('Box','off', 'Location','best');

% Print selected channel list in command window for record keeping.
fprintf('\nFigure 1 averaged across these Depth_s channel indices:\n');
disp(channels_to_analyse);

%% ====================== FIGURE 2: CHANNEL-WISE COMPARISON ======================
% One figure per amplitude.
%
% X-axis:
%   compact channel position
%
% Tick labels:
%   real recording channel index from Depth_s
%
% This is useful when channels_to_analyse is not continuous, for example:
%   channels_to_analyse = [1:7 17:28]
%
% In that case, the plot will not have a large visual gap between 7 and 17,
% but the x-axis labels will still show the real channel numbers.

xpos = 1:numel(channels_to_analyse);

for aa = 1:numel(Amps_selected)

    amp_val = Amps_selected(aa);

    figure('Color','w', ...
           'Name', sprintf('QuickCount_ChannelCompare_Amp%g', amp_val), ...
           'Position', [150 150 900 550]);

    hold on; box off;

    for c = 1:numel(CondList)

        y = squeeze(MeanMat(:,c,aa));
        e = squeeze(SEMMat(:,c,aa));

        if all(isnan(y))
            continue;
        end

        errorbar(xpos, y, e, '-o', ...
            'Color', cond_colors(c,:), ...
            'LineWidth', 1.8, ...
            'MarkerSize', 6, ...
            'DisplayName', CondList(c).Label);
    end

    yline(0, 'k:', 'LineWidth', 1);

    % Compact x-axis positions, but real channel labels.
    xticks(xpos);
    xticklabels(arrayfun(@num2str, channels_to_analyse, 'UniformOutput', false));

    % Rotate labels slightly if many channels are plotted.
    if numel(channels_to_analyse) > 12
        xtickangle(45);
    end

    xlabel('Recording channel index (Depth_s)');
    ylabel('Baseline-corrected spike count / trial');

    title(sprintf('Channel-wise quick spike count | %.1f \\muA | Response %.0f–%.0f ms', ...
          amp_val, response_window_ms(1), response_window_ms(2)), ...
          'Interpreter','none');

    legend('Box','off', 'Location','best');
    % Keep a little space at both ends of the compact x-axis.
    xlim([0.5 numel(channels_to_analyse)+0.5]);
end

fprintf('\nFinished quick spike-count check.\n');

%% ====================== FIGURE 3: CHANNEL-WISE BAR PLOT ======================

xpos = 1:numel(channels_to_analyse);

for aa = 1:numel(Amps_selected)

    amp_val = Amps_selected(aa);

    % Build channel × condition matrix for this amplitude.
    % Rows = selected channels
    % Columns = conditions
    Y = squeeze(MeanMat(:,:,aa));

    % If only one channel or one condition, squeeze can sometimes change
    % dimensions, so force it back to the expected shape.
    if numel(channels_to_analyse) == 1
        Y = reshape(Y, 1, []);
    end

    if numel(CondList) == 1
        Y = reshape(Y, [], 1);
    end

    if all(isnan(Y(:)))
        continue;
    end

    figure('Color','w', ...
           'Name', sprintf('QuickCount_ChannelBar_Amp%g', amp_val), ...
           'Position', [150 150 950 560]);

    hold on; box off;

    % Grouped bar plot.
    % Each channel group contains bars for AutoSim / Seq conditions.
    b = bar(xpos, Y, 'grouped');

    for c = 1:numel(CondList)
        b(c).FaceColor = cond_colors(c,:);
        b(c).DisplayName = CondList(c).Label;
    end

    yline(0, 'k:', 'LineWidth', 1);

    % Compact x-axis positions, but real channel labels.
    xticks(xpos);
    xticklabels(arrayfun(@num2str, channels_to_analyse, 'UniformOutput', false));

    if numel(channels_to_analyse) > 12
        xtickangle(45);
    end

    xlabel('Recording channel index (Depth_s)');
    ylabel('Baseline-corrected spike count / trial');

    title(sprintf(['Channel-wise baseline-corrected spike count | %.1f \\muA | ' ...
                   'Response %.0f–%.0f ms'], ...
          amp_val, response_window_ms(1), response_window_ms(2)), ...
          'Interpreter','none');

    legend('Box','off', 'Location','best');

    xlim([0.5 numel(channels_to_analyse)+0.5]);
end

%% ====================== LOCAL HELPER FUNCTIONS ======================

function setLabel = buildSetLabelFromTrials(trial_list, StimParams_data, simultaneous_stim, E_MAP, isAutoSim)
    % Build stimulation label from a list of trials.
    %
    % This avoids the problem where sequential Level 1 may only show one
    % active electrode. We collect active stimulation channels from the
    % available trials and then convert them through E_MAP.
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
    % Convert Intan stimulation label such as 'A-017' to the corresponding
    % channel number using E_MAP.
    %
    % In your E_MAP format:
    %   E_MAP{1,1} is the array/map name.
    %   E_MAP{2,1} is channel 1.
    %   E_MAP{3,1} is channel 2.
    %   ...
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

    %% ----- Main case: E_MAP is a cell array -----
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

    %% ----- If E_MAP is a string array -----
    if isstring(E_MAP)

        hit = find(strcmp(strtrim(E_MAP), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    %% ----- If E_MAP is a char matrix -----
    if ischar(E_MAP)

        nameList = cellstr(E_MAP);
        hit = find(strcmp(strtrim(nameList), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    %% ----- Fallback only if E_MAP matching fails -----
    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
        warning('Stim channel %s was not found in E_MAP. Falling back to parsed number %d.', ...
            stimName, chNum);
    else
        warning('Stim channel %s was not found in E_MAP and could not be parsed.', stimName);
    end
end