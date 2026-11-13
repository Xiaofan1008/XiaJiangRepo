%% =============================================================
%  Spike Time Bin Summary for Final Pulse-Train Conditions
%
%  Purpose:
%    Summarize spike timing distribution in small time bins for the
%    final recovered pulse-train conditions only.
%
%  This code DOES NOT delete any spikes.
%
%  Main output:
%    SpikeTimeSummary_FinalOnly_PulseTrain.mat
% =============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ========================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train3';

Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Depth_s channels to summarize.
% These are the channel numbers you normally use before Depth_s mapping.
depth_channels_to_summarize = 1:64;

% Spike source:
%   'auto'         = prefer PrefixRecovery_SSD, then PrefixRecovery, then base
%   'recovery_ssd' = *_sp_xia_PrefixRecovery_SSD.mat, variable sp_corr
%   'recovery'     = *_sp_xia_PrefixRecovery.mat, variable sp_seq
%   'base'         = *_sp_xia.mat, variable sp_clipped or sp
spike_source = 'auto';

% Amplitudes to summarize.
% [] means all non-zero amplitudes.
Plot_Amps = [];

% Stimulation set IDs to summarize.
% [] means all detected sets.
Plot_SetIDs = [];

% Which final condition families to summarize.
plot_auto_sim   = true;
plot_sequential = true;

% Include zero-current control?
% Usually false for artifact-bin checking.
include_zero_control = false;

% Detection window for spike-bin summary, relative to trigger.
% Example:
%   [0 40] with 0.5 ms bin gives:
%   0-0.5, 0.5-1.0, ..., 39.5-40 ms
summary_window_ms = [0 40];

% Bin size for summary.
% Smaller bins are better for detecting very aligned artifacts.
summary_bin_ms = 0.5;

% Top-bin reporting thresholds.
% These only affect TopBins inside each condition.
% AllBins still includes every bin, including zero-spike bins.
min_trial_fraction_for_top_bins = 0.20;
min_total_spikes_for_top_bins   = 1;

% Save output name.
summary_output_name = 'SpikeTimeSummary_PulseTrain.mat';

% Print progress?
verbose = true;

%% ====================== CHECK FOLDER =========================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end

cd(data_folder);
fprintf('Spike-time summary in folder:\n%s\n\n', data_folder);

%% ====================== LOAD SPIKES ==========================

sp = [];
spike_file_used = '';
spike_variable_used = '';

switch lower(spike_source)

    case 'auto'

        ssd_recovery_file = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));
        prefix_file       = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
        prefix_file       = prefix_file(~contains({prefix_file.name}, 'SSD'));
        base_file         = dir(fullfile(data_folder, '*sp_xia.mat'));

        if ~isempty(ssd_recovery_file)

            spike_file_used = ssd_recovery_file(1).name;
            S_sp = load(fullfile(data_folder, spike_file_used));

            if isfield(S_sp, 'sp_corr')
                sp = S_sp.sp_corr;
                spike_variable_used = 'sp_corr';
            else
                error('PrefixRecovery_SSD file found but variable sp_corr was not found.');
            end

        elseif ~isempty(prefix_file)

            spike_file_used = prefix_file(1).name;
            S_sp = load(fullfile(data_folder, spike_file_used));

            if isfield(S_sp, 'sp_seq')
                sp = S_sp.sp_seq;
                spike_variable_used = 'sp_seq';
            else
                error('PrefixRecovery file found but variable sp_seq was not found.');
            end

        elseif ~isempty(base_file)

            spike_file_used = base_file(1).name;
            S_sp = load(fullfile(data_folder, spike_file_used));

            if isfield(S_sp, 'sp_clipped')
                sp = S_sp.sp_clipped;
                spike_variable_used = 'sp_clipped';
            elseif isfield(S_sp, 'sp')
                sp = S_sp.sp;
                spike_variable_used = 'sp';
            else
                error('Base spike file found but no usable spike variable was found.');
            end

        else
            error('No usable spike file found.');
        end

    case 'recovery_ssd'

        sp_file = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));
        assert(~isempty(sp_file), 'No *sp_xia_PrefixRecovery_SSD.mat file found.');

        spike_file_used = sp_file(1).name;
        S_sp = load(fullfile(data_folder, spike_file_used));

        if isfield(S_sp, 'sp_corr')
            sp = S_sp.sp_corr;
            spike_variable_used = 'sp_corr';
        else
            error('sp_corr not found in %s.', spike_file_used);
        end

    case 'recovery'

        sp_file = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
        sp_file = sp_file(~contains({sp_file.name}, 'SSD'));
        assert(~isempty(sp_file), 'No *sp_xia_PrefixRecovery.mat file found.');

        spike_file_used = sp_file(1).name;
        S_sp = load(fullfile(data_folder, spike_file_used));

        if isfield(S_sp, 'sp_seq')
            sp = S_sp.sp_seq;
            spike_variable_used = 'sp_seq';
        else
            error('sp_seq not found in %s.', spike_file_used);
        end

    case 'base'

        sp_file = dir(fullfile(data_folder, '*sp_xia.mat'));
        assert(~isempty(sp_file), 'No *sp_xia.mat file found.');

        spike_file_used = sp_file(1).name;
        S_sp = load(fullfile(data_folder, spike_file_used));

        if isfield(S_sp, 'sp_clipped')
            sp = S_sp.sp_clipped;
            spike_variable_used = 'sp_clipped';
        elseif isfield(S_sp, 'sp')
            sp = S_sp.sp;
            spike_variable_used = 'sp';
        else
            error('No usable spike variable found in %s.', spike_file_used);
        end

    otherwise
        error('Unknown spike_source: %s', spike_source);
end

nCh = numel(sp);

fprintf('Loaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);
fprintf('Loaded %d spike channels.\n', nCh);

%% ====================== LOAD TRIGGERS ========================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD EXPERIMENT PARAMETERS ===========

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
TrialParams       = S_exp.TrialParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;

if isfield(S_exp, 'E_MAP')
    E_MAP = S_exp.E_MAP;
else
    E_MAP = [];
    warning('E_MAP not found. Stim labels may be less clean.');
end

if isfield(S_exp, 'StimMeta')
    StimMeta = S_exp.StimMeta;
else
    error('StimMeta was not found. This pulse-train summary requires StimMeta.');
end

fprintf('Loaded exp datafile: %s\n', fileDIR(1).name);
fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match nTrig (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== SAMPLING RATE ========================

try
    [~, freq_params] = read_Intan_RHS2000_file;
    FS = freq_params.amplifier_sample_rate;
catch
    FS = 30000;
    warning('Could not read info.rhs. Using FS = 30000 Hz.');
end

fprintf('Sampling rate: %.1f Hz\n', FS);

%% ====================== REMOVE HEADER ROWS ===================

StimParams_data  = StimParams(2:end,:);
TrialParams_data = TrialParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(StimParams_data,1), expected_rows);
end

if size(TrialParams_data,1) ~= expected_rows
    warning('TrialParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(TrialParams_data,1), expected_rows);
end

%% ====================== TRIAL CONDITION IDS ==================

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of trial-level condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA =================
% This section uses StimMeta, not old mixed-prefix columns.

stimSet_trial        = NaN(n_Trials,1);
trainLevel_trial     = NaN(n_Trials,1);
totalLevels_trial    = NaN(n_Trials,1);
trialAmps            = NaN(n_Trials,1);
isAutoSim_trial      = false(n_Trials,1);
isZeroControl        = false(n_Trials,1);
eventEnd_ms_trial    = NaN(n_Trials,1);

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
            trialAmps(tr) = 0;
        else
            trialAmps(tr) = max(amp_vec(amp_vec > 0));
        end
    end
end

%% ====================== APPLY USER FILTERS ===================

% Amplitudes.
all_amps = unique(trialAmps(~isnan(trialAmps)));
all_amps = all_amps(all_amps > 0);
all_amps = sort(all_amps(:))';

if isempty(Plot_Amps)
    Plot_Amps_selected = all_amps;
else
    Plot_Amps_selected = intersect(all_amps, Plot_Amps);
end

% Sets.
all_sets = unique(stimSet_trial(~isnan(stimSet_trial)));
all_sets = all_sets(all_sets > 0);
all_sets = sort(all_sets(:))';

if isempty(Plot_SetIDs)
    SetIDs_selected = all_sets;
else
    SetIDs_selected = intersect(all_sets, Plot_SetIDs);
end

fprintf('\nSelected amplitudes: ');
disp(Plot_Amps_selected);

fprintf('Selected set IDs: ');
disp(SetIDs_selected);

fprintf('Summarize AutoSim:    %d\n', plot_auto_sim);
fprintf('Summarize sequential: %d\n', plot_sequential);

%% ====================== ELECTRODE MAP ========================

cur_dir = pwd;

try
    cd(data_folder);
    d = Depth_s(Electrode_Type);
    cd(cur_dir);
catch ME
    cd(cur_dir);
    warning('Depth_s failed. Falling back to direct mapping.');
    warning('%s', ME.message);
    d = 1:max(depth_channels_to_summarize);
end

depth_channels_to_summarize = depth_channels_to_summarize(:)';
depth_channels_to_summarize = depth_channels_to_summarize( ...
    depth_channels_to_summarize >= 1 & ...
    depth_channels_to_summarize <= numel(d));

if isempty(depth_channels_to_summarize)
    error('No valid Depth_s channels selected.');
end

fprintf('Number of Depth_s channels to summarize: %d\n', numel(depth_channels_to_summarize));

%% ====================== SUMMARY TIME BINS ====================

if numel(summary_window_ms) ~= 2 || summary_window_ms(2) <= summary_window_ms(1)
    error('summary_window_ms must be [start end] with end > start.');
end

if summary_bin_ms <= 0
    error('summary_bin_ms must be positive.');
end

bin_edges = summary_window_ms(1):summary_bin_ms:summary_window_ms(2);

% Make sure the final edge is exactly the requested end.
if bin_edges(end) < summary_window_ms(2)
    bin_edges = [bin_edges summary_window_ms(2)];
end

bin_starts = bin_edges(1:end-1);
bin_ends   = bin_edges(2:end);
nBins      = numel(bin_starts);

fprintf('\nSummary window: %.3f to %.3f ms\n', summary_window_ms(1), summary_window_ms(2));
fprintf('Summary bin size: %.3f ms\n', summary_bin_ms);
fprintf('Number of bins: %d\n', nBins);

%% ====================== BUILD FINAL CONDITION LIST ============
% Each FinalCond is one final recovered condition:
%   AutoSim: highest TrainLevel for that amp/set/family
%   Seq:     highest TrainLevel for that amp/set/family

FinalConds = struct( ...
    'ConditionName', {}, ...
    'ConditionFamily', {}, ...
    'StimSetID', {}, ...
    'Amplitude_uA', {}, ...
    'FinalTrainLevel', {}, ...
    'TrialList', {}, ...
    'ConditionLabel', {}, ...
    'EventTimes_ms', {}, ...
    'EventEnd_ms', {});

for ai = 1:numel(Plot_Amps_selected)

    ampVal = Plot_Amps_selected(ai);

    %% ---------------- AutoSim final conditions ----------------
    if plot_auto_sim

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trials_base = find(abs(trialAmps - ampVal) < 1e-9 & ...
                               stimSet_trial == setID & ...
                               isAutoSim_trial == true);

            if ~include_zero_control
                trials_base = trials_base(isZeroControl(trials_base) == false);
            end

            trials_base = trials_base(trials_base <= nTrials_use);

            if isempty(trials_base)
                continue;
            end

            finalLevel = max(trainLevel_trial(trials_base));
            trials_final = trials_base(trainLevel_trial(trials_base) == finalLevel);

            if isempty(trials_final)
                continue;
            end

            tr_rep = trials_final(1);

            conditionLabel = getStimLabelForTrial( ...
                tr_rep, StimParams_data, simultaneous_stim, E_MAP, true);

            conditionName = sprintf('AutoSim | %s | Set %g | %.1f uA | final level %g', ...
                conditionLabel, setID, ampVal, finalLevel);

            C.ConditionName = conditionName;
            C.ConditionFamily = 'AutoSim';
            C.StimSetID = setID;
            C.Amplitude_uA = ampVal;
            C.FinalTrainLevel = finalLevel;
            C.TrialList = trials_final(:)';
            C.ConditionLabel = conditionLabel;
            C.EventTimes_ms = eventTimes_ms_trial{tr_rep};
            C.EventEnd_ms = eventEnd_ms_trial(tr_rep);

            FinalConds(end+1) = C; %#ok<SAGROW>
        end
    end

    %% ---------------- Sequential final conditions --------------
    if plot_sequential

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trials_base = find(abs(trialAmps - ampVal) < 1e-9 & ...
                               stimSet_trial == setID & ...
                               isAutoSim_trial == false);

            if ~include_zero_control
                trials_base = trials_base(isZeroControl(trials_base) == false);
            end

            trials_base = trials_base(trials_base <= nTrials_use);

            if isempty(trials_base)
                continue;
            end

            finalLevel = max(trainLevel_trial(trials_base));
            trials_final = trials_base(trainLevel_trial(trials_base) == finalLevel);

            if isempty(trials_final)
                continue;
            end

            tr_rep = trials_final(1);

            conditionLabel = getStimLabelForTrial( ...
                tr_rep, StimParams_data, simultaneous_stim, E_MAP, false);

            conditionName = sprintf('Seq | %s | Set %g | %.1f uA | final level %g', ...
                conditionLabel, setID, ampVal, finalLevel);

            C.ConditionName = conditionName;
            C.ConditionFamily = 'Seq';
            C.StimSetID = setID;
            C.Amplitude_uA = ampVal;
            C.FinalTrainLevel = finalLevel;
            C.TrialList = trials_final(:)';
            C.ConditionLabel = conditionLabel;
            C.EventTimes_ms = eventTimes_ms_trial{tr_rep};
            C.EventEnd_ms = eventEnd_ms_trial(tr_rep);

            FinalConds(end+1) = C; %#ok<SAGROW>
        end
    end
end

if isempty(FinalConds)
    error('No final conditions found. Check Plot_Amps, Plot_SetIDs, and metadata.');
end

fprintf('\n================ FINAL CONDITIONS SUMMARIZED ================\n');
for ci = 1:numel(FinalConds)
    fprintf('%2d | %s | Trials %d\n', ...
        ci, FinalConds(ci).ConditionName, numel(FinalConds(ci).TrialList));

    fprintf('     Events: ');
    disp(FinalConds(ci).EventTimes_ms);
end
fprintf('=============================================================\n');

%% ====================== INITIALIZE OUTPUT STRUCTURE ===========

SummarySettings = struct();

SummarySettings.DataFolder = data_folder;
SummarySettings.SpikeFileUsed = spike_file_used;
SummarySettings.SpikeVariableUsed = spike_variable_used;
SummarySettings.ExpDataFileUsed = fileDIR(1).name;

SummarySettings.Electrode_Type = Electrode_Type;
SummarySettings.DepthChannelsSummarized = depth_channels_to_summarize;
SummarySettings.DepthToSpikeChannelMap = d;

SummarySettings.Plot_Amps_Input = Plot_Amps;
SummarySettings.Plot_Amps_Selected = Plot_Amps_selected;
SummarySettings.Plot_SetIDs_Input = Plot_SetIDs;
SummarySettings.SetIDs_Selected = SetIDs_selected;

SummarySettings.PlotAutoSim = plot_auto_sim;
SummarySettings.PlotSequential = plot_sequential;
SummarySettings.IncludeZeroControl = include_zero_control;
SummarySettings.OnlyFinalTrainLevel = true;

SummarySettings.SummaryWindow_ms = summary_window_ms;
SummarySettings.SummaryBin_ms = summary_bin_ms;
SummarySettings.BinEdges_ms = bin_edges;

SummarySettings.MinTrialFractionForTopBins = min_trial_fraction_for_top_bins;
SummarySettings.MinTotalSpikesForTopBins = min_total_spikes_for_top_bins;

SummarySettings.FS = FS;
SummarySettings.nTrials_use = nTrials_use;
SummarySettings.DateGenerated = datestr(now);

SpikeSummary = struct();
SpikeSummary.Settings = SummarySettings;
SpikeSummary.FinalConditionInfo = rmfield(FinalConds, 'TrialList');

maxDepthCh = max(depth_channels_to_summarize);

emptyConditionStruct = struct( ...
    'ConditionName', '', ...
    'ConditionFamily', '', ...
    'ConditionLabel', '', ...
    'StimSetID', NaN, ...
    'Amplitude_uA', NaN, ...
    'FinalTrainLevel', NaN, ...
    'NumberOfTrials', NaN, ...
    'EventTimes_ms', [], ...
    'EventEnd_ms', NaN, ...
    'AllBins', table(), ...
    'TopBins', table());

emptyChannelStruct = struct( ...
    'DepthChannel', [], ...
    'Overview', table(), ...
    'Condition', emptyConditionStruct);

SpikeSummary.ByChannel = repmat(emptyChannelStruct, maxDepthCh, 1);

%% ====================== MAIN SUMMARY LOOP =====================

fprintf('\nBuilding spike-time bin summary...\n');

for di = 1:numel(depth_channels_to_summarize)

    depthCh = depth_channels_to_summarize(di);
    spikeCh = d(depthCh);

    if verbose
        fprintf('  Depth channel %d (%d/%d)\n', ...
            depthCh, di, numel(depth_channels_to_summarize));
    end

    SpikeSummary.ByChannel(depthCh).DepthChannel = depthCh;

    if spikeCh < 1 || spikeCh > nCh
        warning('DepthChannel %d maps to spike channel %d, outside spike file range. Skipped.', ...
            depthCh, spikeCh);
        continue;
    end

    if isempty(sp{spikeCh})
        sp_times = [];
    else
        sp_times = sp{spikeCh}(:,1);
    end

    % Prepare overview rows for this channel.
    Overview_ConditionIndex = [];
    Overview_ConditionName = {};
    Overview_ConditionFamily = {};
    Overview_ConditionLabel = {};
    Overview_Amplitude_uA = [];
    Overview_NumberOfTrials = [];
    Overview_TotalSpikeCountInDetectionWindow = [];
    Overview_MaxTrialFractionWithSpike = [];
    Overview_MostSuspiciousBinStart_ms = [];
    Overview_MostSuspiciousBinEnd_ms = [];
    Overview_MostSuspiciousBinTotalSpikeCount = [];

    % Preallocate condition structure for this channel.
    ChannelConditions = repmat(emptyConditionStruct, numel(FinalConds), 1);

    for ci = 1:numel(FinalConds)

        C = FinalConds(ci);

        trials_this = C.TrialList;
        trials_this = trials_this(trials_this <= nTrials_use);
        nTrials_cond = numel(trials_this);

        if nTrials_cond == 0
            continue;
        end

        total_counts = zeros(1, nBins);
        trials_with_spike = zeros(1, nBins);

        %% ----- Count spikes per time bin -----
        for ti = 1:nTrials_cond

            tr = trials_this(ti);

            if tr > numel(trig)
                continue;
            end

            t0_ms = trig(tr) / FS * 1000;

            if isempty(sp_times)
                rel_t = [];
            else
                rel_t = sp_times - t0_ms;
                rel_t = rel_t(rel_t >= summary_window_ms(1) & ...
                              rel_t <  summary_window_ms(2));
            end

            if isempty(rel_t)
                trial_counts = zeros(1, nBins);
            else
                trial_counts = histcounts(rel_t, bin_edges);
            end

            total_counts = total_counts + trial_counts;
            trials_with_spike = trials_with_spike + double(trial_counts > 0);
        end

        trial_fraction = trials_with_spike ./ nTrials_cond;

        %% ----- Clean bin table for this channel-condition -----
        AllBins = table( ...
            bin_starts(:), ...
            bin_ends(:), ...
            total_counts(:), ...
            trials_with_spike(:), ...
            trial_fraction(:), ...
            'VariableNames', { ...
                'TimeBinStart_ms', ...
                'TimeBinEnd_ms', ...
                'TotalSpikeCount', ...
                'TrialsWithAtLeastOneSpike', ...
                'TrialFractionWithSpike'});

        top_mask = AllBins.TrialFractionWithSpike >= min_trial_fraction_for_top_bins & ...
                   AllBins.TotalSpikeCount >= min_total_spikes_for_top_bins;

        TopBins = AllBins(top_mask, :);

        if ~isempty(TopBins)
            TopBins = sortrows(TopBins, ...
                {'TrialFractionWithSpike', 'TotalSpikeCount'}, ...
                {'descend', 'descend'});
        end

        %% ----- Overview values -----
        total_spikes_window = sum(total_counts);

        if total_spikes_window == 0
            max_trial_fraction = 0;
            most_bin_start = NaN;
            most_bin_end = NaN;
            most_bin_total_spikes = 0;
        else
            max_trial_fraction = max(trial_fraction);

            % If multiple bins share the same trial fraction, choose the one
            % with the largest total spike count.
            candidate_bins = find(abs(trial_fraction - max_trial_fraction) < 1e-12);

            if numel(candidate_bins) > 1
                [~, local_best] = max(total_counts(candidate_bins));
                best_bin = candidate_bins(local_best);
            else
                best_bin = candidate_bins(1);
            end

            most_bin_start = bin_starts(best_bin);
            most_bin_end = bin_ends(best_bin);
            most_bin_total_spikes = total_counts(best_bin);
        end

        %% ----- Store condition-level result -----
        ChannelConditions(ci).ConditionName = C.ConditionName;
        ChannelConditions(ci).ConditionFamily = C.ConditionFamily;
        ChannelConditions(ci).ConditionLabel = C.ConditionLabel;
        ChannelConditions(ci).StimSetID = C.StimSetID;
        ChannelConditions(ci).Amplitude_uA = C.Amplitude_uA;
        ChannelConditions(ci).FinalTrainLevel = C.FinalTrainLevel;
        ChannelConditions(ci).NumberOfTrials = nTrials_cond;
        ChannelConditions(ci).EventTimes_ms = C.EventTimes_ms;
        ChannelConditions(ci).EventEnd_ms = C.EventEnd_ms;
        ChannelConditions(ci).AllBins = AllBins;
        ChannelConditions(ci).TopBins = TopBins;

        %% ----- Add overview row -----
        Overview_ConditionIndex = [Overview_ConditionIndex; ci]; %#ok<AGROW>
        Overview_ConditionName = [Overview_ConditionName; {C.ConditionName}]; %#ok<AGROW>
        Overview_ConditionFamily = [Overview_ConditionFamily; {C.ConditionFamily}]; %#ok<AGROW>
        Overview_ConditionLabel = [Overview_ConditionLabel; {C.ConditionLabel}]; %#ok<AGROW>
        Overview_Amplitude_uA = [Overview_Amplitude_uA; C.Amplitude_uA]; %#ok<AGROW>
        Overview_NumberOfTrials = [Overview_NumberOfTrials; nTrials_cond]; %#ok<AGROW>
        Overview_TotalSpikeCountInDetectionWindow = [Overview_TotalSpikeCountInDetectionWindow; total_spikes_window]; %#ok<AGROW>
        Overview_MaxTrialFractionWithSpike = [Overview_MaxTrialFractionWithSpike; max_trial_fraction]; %#ok<AGROW>
        Overview_MostSuspiciousBinStart_ms = [Overview_MostSuspiciousBinStart_ms; most_bin_start]; %#ok<AGROW>
        Overview_MostSuspiciousBinEnd_ms = [Overview_MostSuspiciousBinEnd_ms; most_bin_end]; %#ok<AGROW>
        Overview_MostSuspiciousBinTotalSpikeCount = [Overview_MostSuspiciousBinTotalSpikeCount; most_bin_total_spikes]; %#ok<AGROW>
    end

    %% ----- Channel overview table -----
    ChannelOverview = table( ...
        Overview_ConditionIndex, ...
        Overview_ConditionName, ...
        Overview_ConditionFamily, ...
        Overview_ConditionLabel, ...
        Overview_Amplitude_uA, ...
        Overview_NumberOfTrials, ...
        Overview_TotalSpikeCountInDetectionWindow, ...
        Overview_MaxTrialFractionWithSpike, ...
        Overview_MostSuspiciousBinStart_ms, ...
        Overview_MostSuspiciousBinEnd_ms, ...
        Overview_MostSuspiciousBinTotalSpikeCount, ...
        'VariableNames', { ...
            'ConditionIndex', ...
            'ConditionName', ...
            'ConditionFamily', ...
            'ConditionLabel', ...
            'Amplitude_uA', ...
            'NumberOfTrials', ...
            'TotalSpikeCountInDetectionWindow', ...
            'MaxTrialFractionWithSpike', ...
            'MostSuspiciousBinStart_ms', ...
            'MostSuspiciousBinEnd_ms', ...
            'MostSuspiciousBinTotalSpikeCount'});

    SpikeSummary.ByChannel(depthCh).Overview = ChannelOverview;
    SpikeSummary.ByChannel(depthCh).Condition = ChannelConditions;
end

%% ====================== SAVE OUTPUT ===========================

summary_output_path = fullfile(data_folder, summary_output_name);

save(summary_output_path, ...
    'SpikeSummary', ...
    'SummarySettings', ...
    '-v7.3');

fprintf('\nSaved spike-time summary to:\n%s\n', summary_output_path);

fprintf('\nSaved variables:\n');
fprintf('  SpikeSummary\n');
fprintf('  SummarySettings\n');

fprintf('\nExample commands after loading:\n');
fprintf('  load(''%s'')\n', summary_output_name);
fprintf('  SpikeSummary.ByChannel(25).Overview\n');
fprintf('  SpikeSummary.ByChannel(25).Condition(1).ConditionName\n');
fprintf('  SpikeSummary.ByChannel(25).Condition(1).AllBins\n');
fprintf('  SpikeSummary.ByChannel(25).Condition(1).TopBins\n');

fprintf('\nFinished spike-time summary. No spikes were deleted.\n');

%% =============================================================
%  LOCAL HELPER FUNCTIONS
%% =============================================================

function label = getStimLabelForTrial(tr, StimParams_data, simultaneous_stim, E_MAP, isAutoSim)
    % Build stimulation label for one representative trial.

    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);
    rr = rr(rr <= size(StimParams_data,1));

    if isempty(rr)
        label = sprintf('Trial%d', tr);
        return;
    end

    stimNames = StimParams_data(rr,1);

    try
        ampVec = cell2mat(StimParams_data(rr,16));
    catch
        ampVec = ones(numel(rr),1);
    end

    activeRows = ampVec > 0;
    stimNames_active = stimNames(activeRows);

    label = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim);
end

function setLabel = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim)
    % Build stimulation label for condition identification.
    %
    % AutoSim:
    %   Ch35+Ch39
    %
    % Sequential:
    %   Ch35→Ch39

    if isempty(stimNames_active)
        setLabel = 'NoActiveStim';
        return;
    end

    stimNames_active = unique(stimNames_active, 'stable');

    labelParts = cell(1, numel(stimNames_active));

    for i = 1:numel(stimNames_active)

        stimName = stimNames_active{i};

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
    % Convert Intan stimulation label to channel number using E_MAP.
    %
    % E_MAP convention:
    %   E_MAP{1,1} = header / map name
    %   E_MAP{2,1} = channel 1
    %   E_MAP{3,1} = channel 2
    %
    % Therefore:
    %   channel number = row index - 1

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
    % Fallback if E_MAP lookup fails.

    chNum = NaN;

    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
    end
end