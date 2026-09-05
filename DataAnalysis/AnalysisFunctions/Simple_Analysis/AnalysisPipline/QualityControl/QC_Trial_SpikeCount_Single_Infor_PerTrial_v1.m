%% ========================================================================
% SINGLE-STIMULATION TRIAL SPIKE-COUNT TABLE
%
% PURPOSE
%   Count baseline and post-stimulation spikes for every recording channel
%   and every original trial in a single-stimulation dataset.
%
% OUTPUTS
%   1. SpikeCounts.set(si).amp(ai).ptd(1).trial(k)
%      Detailed per-channel counts for each condition-specific raster row.
%
%   2. TrialSummary
%      One row per original experiment trial. Rows remain in their original
%      trial order and are never sorted by this script.
%
% IMPORTANT TRIAL IDENTIFIERS
%   AbsoluteTrialID:
%       Original experiment/trigger trial number. Use this number later
%       when entering bad trials.
%
%   ConditionTrialIndex:
%       Trial row within the corresponding stimulation set and amplitude.
%       This matches the row position in a condition-specific raster.
%
% DATA RULES
%   - Requires sp_corr from *.sp_xia_SSD.mat.
%   - Automatically includes every recording channel.
%   - Does not use responding-channel results.
%   - Does not remove bad channels or bad trials.
%   - Does not modify the original experimental data.
% ========================================================================

clear;
close all;

addpath(genpath( ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ========================= USER SETTINGS ==============================

data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_Single1';

% Electrode type:
%   0 = rigid single shank
%   1 = flexible single shank
%   2 = four-shank flexible probe
Electrode_Type = 2;

% Sampling rate in Hz
FS = 30000;

% Baseline window
% The upper boundary is excluded: [-50,-5) ms
baseline_win_ms = [-50 -5];

% Short response window used for responding-channel analysis
% The upper boundary is excluded: [2,20) ms
response_win_ms = [2 20];

% Whole-response window intended for later multi-ISI analysis
% The upper boundary is excluded: [2,40) ms
analysis_win_ms = [2 40];

% If an output file already exists, create a timestamped backup
Create_Backup_If_Output_Exists = true;

%% ======================== INITIAL CHECKS ==============================

if ~isfolder(data_folder)
    error('Dataset folder does not exist:\n%s',data_folder);
end

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(response_win_ms,'response_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');

% Remember the current folder and restore it when the script finishes
starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

% Change automatically to the selected dataset folder
cd(data_folder);

fprintf('\n============================================================\n');
fprintf('SINGLE-STIMULATION TRIAL SPIKE-COUNT TABLE\n');
fprintf('============================================================\n');
fprintf('Dataset: %s\n',data_folder);
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Short response window: [%g,%g) ms\n',response_win_ms);
fprintf('Whole-analysis window: [%g,%g) ms\n',analysis_win_ms);
fprintf('Responding-channel file used: NO\n');
fprintf('Bad channels removed: NO\n');
fprintf('Bad trials removed: NO\n');

%% ========================== LOAD sp_corr ==============================

ssd_files = dir('*.sp_xia_SSD.mat');

if isempty(ssd_files)
    error('No *.sp_xia_SSD.mat file was found.');
end

if numel(ssd_files) > 1
    warning('Multiple SSD files found. Using: %s',ssd_files(1).name);
end

ssd_file = ssd_files(1).name;
base_name = erase(ssd_file,'.sp_xia_SSD.mat');

% This script intentionally requires the final filtered spike variable
if ~ismember('sp_corr',who('-file',ssd_file))
    error('The SSD file does not contain sp_corr:\n%s',ssd_file);
end

SpikeLoad = load(ssd_file,'sp_corr');
sp = SpikeLoad.sp_corr;
nSpChannels = numel(sp);

fprintf('Spike file: %s\n',ssd_file);
fprintf('Spike-data channels: %d\n',nSpChannels);

%% ========================== LOAD TRIGGERS =============================

% Create the cleaned trigger file if necessary
if isempty(dir('*.trig.dat'))
    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);
end

if isempty(dir('*.trig.dat'))
    error('No *.trig.dat file could be found or created.');
end

trig = double(loadTrig(0));
trig = trig(:);
nTrig = numel(trig);

%% ================= LOAD EXPERIMENT PARAMETERS =========================

experiment_files = dir('*_exp_datafile_*.mat');

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found.');
end

if numel(experiment_files) > 1
    warning('Multiple experiment files found. Using: %s', ...
        experiment_files(1).name);
end

ExpLoad = load(experiment_files(1).name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');

required_fields = {'StimParams','simultaneous_stim','E_MAP','n_Trials'};

for field_id = 1:numel(required_fields)
    if ~isfield(ExpLoad,required_fields{field_id})
        error('Experiment file is missing variable: %s', ...
            required_fields{field_id});
    end
end

StimParams = ExpLoad.StimParams;
sim_stim   = ExpLoad.simultaneous_stim;
E_MAP      = ExpLoad.E_MAP;
n_Trials   = ExpLoad.n_Trials;

% This version is specifically for single-electrode stimulation
if sim_stim ~= 1
    error(['This script requires single stimulation, but ' ...
        'simultaneous_stim = %d.'],sim_stim);
end

% Check agreement between trial and trigger counts
if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.', ...
        nTrig,n_Trials);

elseif nTrig > n_Trials
    warning('%d triggers were loaded for %d trials; extras are ignored.', ...
        nTrig,n_Trials);
end

fprintf('Experiment file: %s\n',experiment_files(1).name);
fprintf('Trials/triggers: %d/%d\n',n_Trials,nTrig);

%% ===================== AMPLITUDES AND SETS ============================

% Stimulation amplitude is stored in column 16
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:sim_stim:end);

if numel(trialAmps) ~= n_Trials
    error('Decoded %d amplitudes for %d trials.', ...
        numel(trialAmps),n_Trials);
end

% Convert the experimental zero-amplitude representation to zero
trialAmps(trialAmps == -1) = 0;

[Amps,~,ampIdx] = unique(trialAmps(:));
nAMP = numel(Amps);

% Map stimulation electrode names to electrode-map indices
stimNames = StimParams(2:end,1);
[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));

if any(~isMapped)
    error('%d stimulation entries could not be mapped through E_MAP.', ...
        sum(~isMapped));
end

% Decode the stimulation electrode used in every trial
stimSeq = zeros(n_Trials,sim_stim);

for trial_id = 1:n_Trials
    rows_this_trial = (trial_id-1)*sim_stim+(1:sim_stim);

    stimulation_indices = idx_all(rows_this_trial);
    stimulation_indices = stimulation_indices(stimulation_indices > 0);

    stimSeq(trial_id,1:numel(stimulation_indices)) = ...
        stimulation_indices;
end

% Each unique stimulation electrode becomes one stimulation set
[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

% Single stimulation has no pulse-time delay
% A PTD=0 layer is retained for structural compatibility
PTDs = 0;

fprintf('Amplitudes: %s uA\n',num2str(Amps(:).'));
fprintf('Stimulation sets: %d\n',nSets);

for si = 1:nSets
    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    fprintf('  Set %d: %s\n', ...
        si,channel_list(stim_channels));
end

%% ========================= CHANNEL SELECTION ==========================

% Depth_s returns the mapping between displayed channel indices and the
% corresponding cells in sp_corr.
d = double(Depth_s(Electrode_Type));
d = d(:);

if isempty(d)
    error('Depth_s returned an empty channel map.');
end

% Automatically include every recording channel represented by Depth_s.
%
% Examples:
%   A 32-channel dataset normally produces selected_channels = 1:32.
%   A 64-channel dataset normally produces selected_channels = 1:64.
%
% No manual channel_start or channel_end setting is required.
selected_channels = 1:numel(d);

nSelectedChannels = numel(selected_channels);

% Convert displayed channel indices into the corresponding sp_corr cells
spike_data_channels = d(selected_channels);

% Cache the spike-time columns once so the full waveform matrices do not
% need to be repeatedly accessed for every trial.
selected_spike_times = cell(nSelectedChannels,1);
channel_has_spike_data = false(nSelectedChannels,1);

for channel_position = 1:nSelectedChannels

    spike_channel = spike_data_channels(channel_position);

    if spike_channel >= 1 && ...
            spike_channel <= nSpChannels && ...
            ~isempty(sp{spike_channel})

        % Column 1 contains absolute spike times in milliseconds
        selected_spike_times{channel_position} = ...
            double(sp{spike_channel}(:,1));

        channel_has_spike_data(channel_position) = true;

    else
        selected_spike_times{channel_position} = [];
    end
end

fprintf('Automatically detected channel indices: 1:%d\n', ...
    nSelectedChannels);

fprintf('Channels with spike data: %d/%d\n', ...
    sum(channel_has_spike_data),nSelectedChannels);

%% ================= INITIALIZE NESTED STRUCTURE ========================

SpikeCounts = struct();

SpikeCounts.metadata.analysis_name = ...
    'Single-stimulation trial spike-count table';

SpikeCounts.metadata.created = char(datetime('now', ...
    'Format','yyyy-MM-dd HH:mm:ss'));

SpikeCounts.metadata.data_folder = data_folder;
SpikeCounts.metadata.spike_file = ssd_file;
SpikeCounts.metadata.spike_variable = 'sp_corr';
SpikeCounts.metadata.experiment_file = experiment_files(1).name;
SpikeCounts.metadata.sampling_rate_hz = FS;
SpikeCounts.metadata.electrode_type = Electrode_Type;

SpikeCounts.metadata.baseline_window_ms = baseline_win_ms;
SpikeCounts.metadata.response_window_ms = response_win_ms;
SpikeCounts.metadata.analysis_window_ms = analysis_win_ms;

SpikeCounts.metadata.responding_channels_used = false;
SpikeCounts.metadata.bad_channels_removed = false;
SpikeCounts.metadata.bad_trials_removed = false;
SpikeCounts.metadata.trial_summary_sorted = false;

SpikeCounts.metadata.channel_selection = ...
    'Automatically included all channel indices returned by Depth_s';

SpikeCounts.metadata.selected_channels = selected_channels;

% Create the condition structure
for si = 1:nSets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    SpikeCounts.set(si).set_index = si;
    SpikeCounts.set(si).set_name = channel_list(stim_channels);
    SpikeCounts.set(si).stim_channels = stim_channels;

    for ai = 1:nAMP

        condition_trials = find( ...
            combClass == si & ampIdx == ai);

        SpikeCounts.set(si).amp(ai).amp_index = ai;
        SpikeCounts.set(si).amp(ai).amp_value = Amps(ai);
        SpikeCounts.set(si).amp(ai).amp_label = ...
            sprintf('%g uA',Amps(ai));

        SpikeCounts.set(si).amp(ai).ptd(1).PTD_ms = 0;

        SpikeCounts.set(si).amp(ai).ptd(1).condition_summary = ...
            sprintf('Set %s | Amp %g uA | Single', ...
            channel_list(stim_channels),Amps(ai));

        SpikeCounts.set(si).amp(ai).ptd(1).trial_ids = ...
            condition_trials(:).';

        SpikeCounts.set(si).amp(ai).ptd(1).total_trials = ...
            numel(condition_trials);
    end
end

%% ================= INITIALIZE FLAT TRIAL SUMMARY ======================

% Absolute trial IDs are the original experiment/trigger trial numbers
AbsoluteTrialID = (1:n_Trials).';
OriginalTrialSequence = AbsoluteTrialID;

% ConditionTrialIndex is the trial number within its stimulation condition
ConditionTrialIndex = zeros(n_Trials,1);

SetIndex = combClass(:);
StimChannelIndex = zeros(n_Trials,1);
AmplitudeIndex = ampIdx(:);
Amplitude_uA = trialAmps(:);
PTD_ms = zeros(n_Trials,1);

N_SelectedChannels = repmat(nSelectedChannels,n_Trials,1);

% Total counts across all recording channels
TotalBaselineSpikes = zeros(n_Trials,1);
TotalResponseWindowSpikes = zeros(n_Trials,1);
TotalAnalysisWindowSpikes = zeros(n_Trials,1);

% Signed baseline-corrected totals
TotalBaselineCorrectedResponseWindowSpikes = zeros(n_Trials,1);
TotalBaselineCorrectedAnalysisWindowSpikes = zeros(n_Trials,1);

% Largest raw count observed in one recording channel
MaxBaselineSpikes_OneChannel = zeros(n_Trials,1);
MaxResponseSpikes_OneChannel = zeros(n_Trials,1);
MaxAnalysisSpikes_OneChannel = zeros(n_Trials,1);

condition_counters = zeros(nSets,nAMP);

% Calculate window durations once
baseline_duration_ms = diff(baseline_win_ms);
response_duration_ms = diff(response_win_ms);
analysis_duration_ms = diff(analysis_win_ms);

%% ================= COUNT IN ORIGINAL TRIAL ORDER ======================

fprintf('\nCounting spikes in original trial order...\n');

for trial_id = 1:n_Trials

    si = combClass(trial_id);
    ai = ampIdx(trial_id);

    condition_counters(si,ai) = ...
        condition_counters(si,ai)+1;

    condition_trial_index = condition_counters(si,ai);

    ConditionTrialIndex(trial_id) = condition_trial_index;

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    if ~isempty(stim_channels)
        StimChannelIndex(trial_id) = stim_channels(1);
    end

    % Convert the trigger sample to absolute time in milliseconds
    trigger_ms = trig(trial_id)/FS*1000;

    % Create absolute counting-window boundaries
    baseline_absolute = trigger_ms+baseline_win_ms;
    response_absolute = trigger_ms+response_win_ms;
    analysis_absolute = trigger_ms+analysis_win_ms;

    baseline_counts = zeros(nSelectedChannels,1);
    response_counts = zeros(nSelectedChannels,1);
    analysis_counts = zeros(nSelectedChannels,1);

    % Count spikes independently for every recording channel
    for channel_position = 1:nSelectedChannels

        if ~channel_has_spike_data(channel_position)
            continue;
        end

        spike_times = selected_spike_times{channel_position};

        % Baseline count in [-50,-5) ms
        baseline_counts(channel_position) = sum( ...
            spike_times >= baseline_absolute(1) & ...
            spike_times <  baseline_absolute(2));

        % Short response count in [2,20) ms
        response_counts(channel_position) = sum( ...
            spike_times >= response_absolute(1) & ...
            spike_times <  response_absolute(2));

        % Whole-analysis count in [2,40) ms
        analysis_counts(channel_position) = sum( ...
            spike_times >= analysis_absolute(1) & ...
            spike_times <  analysis_absolute(2));
    end

    %% ---------------- SIGNED BASELINE CORRECTION ----------------------

    % Estimate the expected number of baseline spikes for a window with
    % the same duration as the short response window.
    expected_baseline_response = baseline_counts * ...
        (response_duration_ms/baseline_duration_ms);

    % Estimate the expected number of baseline spikes for a window with
    % the same duration as the whole-analysis window.
    expected_baseline_analysis = baseline_counts * ...
        (analysis_duration_ms/baseline_duration_ms);

    % Preserve signed values. Negative and fractional values are valid and
    % useful when identifying unusually silent or noisy trials.
    corrected_response_counts = ...
        response_counts-expected_baseline_response;

    corrected_analysis_counts = ...
        analysis_counts-expected_baseline_analysis;

    %% ---------------- PER-CHANNEL RESULT TABLE ------------------------

    Channel_Results = table( ...
        selected_channels(:), ...
        baseline_counts, ...
        response_counts, ...
        analysis_counts, ...
        corrected_response_counts, ...
        corrected_analysis_counts, ...
        'VariableNames',{ ...
        'Channel_Index', ...
        'Baseline_Spikes', ...
        'ResponseWindow_Spikes', ...
        'AnalysisWindow_Spikes', ...
        'BaselineCorrected_ResponseWindow_Spikes', ...
        'BaselineCorrected_AnalysisWindow_Spikes'});

    %% ---------------- STORE THIS INDIVIDUAL TRIAL ---------------------

    T = struct();

    % Trial number within this stimulation condition
    T.relative_trial_ID = condition_trial_index;
    T.ConditionTrialIndex = condition_trial_index;

    % Original experiment trial ID
    T.absolute_trial_ID = trial_id;
    T.AbsoluteTrialID = trial_id;
    T.OriginalTrialSequence = trial_id;

    T.PTD_ms = 0;
    T.Channel_Results = Channel_Results;

    % Raw totals across all automatically detected recording channels
    T.total_baseline_spikes = sum(baseline_counts);

    % Retained for compatibility with older scripts
    T.total_post_spikes = sum(response_counts);

    T.total_response_window_spikes = sum(response_counts);
    T.total_analysis_window_spikes = sum(analysis_counts);

    % Signed baseline-corrected totals
    T.total_baseline_corrected_response_window_spikes = ...
        sum(corrected_response_counts);

    T.total_baseline_corrected_analysis_window_spikes = ...
        sum(corrected_analysis_counts);

    % Store in its stimulation-set and amplitude condition
    SpikeCounts.set(si).amp(ai).ptd(1).trial( ...
        condition_trial_index) = T;

    %% ---------------- UPDATE FLAT TRIAL SUMMARY -----------------------

    TotalBaselineSpikes(trial_id) = ...
        sum(baseline_counts);

    TotalResponseWindowSpikes(trial_id) = ...
        sum(response_counts);

    TotalAnalysisWindowSpikes(trial_id) = ...
        sum(analysis_counts);

    TotalBaselineCorrectedResponseWindowSpikes(trial_id) = ...
        sum(corrected_response_counts);

    TotalBaselineCorrectedAnalysisWindowSpikes(trial_id) = ...
        sum(corrected_analysis_counts);

    MaxBaselineSpikes_OneChannel(trial_id) = ...
        max(baseline_counts);

    MaxResponseSpikes_OneChannel(trial_id) = ...
        max(response_counts);

    MaxAnalysisSpikes_OneChannel(trial_id) = ...
        max(analysis_counts);
end

%% ======================== BUILD SUMMARY TABLE =========================

% Do not sort this table.
% Row n corresponds directly to original experiment trial n.
TrialSummary = table( ...
    AbsoluteTrialID, ...
    OriginalTrialSequence, ...
    ConditionTrialIndex, ...
    SetIndex, ...
    StimChannelIndex, ...
    AmplitudeIndex, ...
    Amplitude_uA, ...
    PTD_ms, ...
    N_SelectedChannels, ...
    TotalBaselineSpikes, ...
    TotalResponseWindowSpikes, ...
    TotalAnalysisWindowSpikes, ...
    TotalBaselineCorrectedResponseWindowSpikes, ...
    TotalBaselineCorrectedAnalysisWindowSpikes, ...
    MaxBaselineSpikes_OneChannel, ...
    MaxResponseSpikes_OneChannel, ...
    MaxAnalysisSpikes_OneChannel);

%% =========================== SAVE RESULT ==============================

save_name = sprintf( ...
    '%s_Single_TrialSpikeCounts_PerTrial.mat',base_name);

full_save_path = fullfile(data_folder,save_name);

% Protect an existing result before overwriting it
if isfile(full_save_path) && Create_Backup_If_Output_Exists

    timestamp = datestr(now,'yyyymmdd_HHMMSS');

    backup_name = sprintf( ...
        '%s_Single_TrialSpikeCounts_PerTrial_BACKUP_%s.mat', ...
        base_name,timestamp);

    backup_path = fullfile(data_folder,backup_name);

    copyfile(full_save_path,backup_path);

    fprintf('Existing output backed up to:\n%s\n', ...
        backup_path);
end

% post_win_ms is retained as an alias for compatibility with older scripts
post_win_ms = response_win_ms;

save(full_save_path, ...
    'SpikeCounts', ...
    'TrialSummary', ...
    'baseline_win_ms', ...
    'response_win_ms', ...
    'analysis_win_ms', ...
    'post_win_ms', ...
    'Amps', ...
    'PTDs', ...
    'uniqueComb', ...
    'selected_channels');

fprintf('\n============================================================\n');
fprintf('SPIKE-COUNT EXTRACTION COMPLETE\n');
fprintf('Trials saved: %d\n',height(TrialSummary));
fprintf('Recording channels detected: %d\n',nSelectedChannels);
fprintf('TrialSummary remains in original trial order: YES\n');
fprintf('Use AbsoluteTrialID when entering bad trials.\n');
fprintf('Output: %s\n',full_save_path);
fprintf('============================================================\n');

%% ========================== LOCAL FUNCTIONS ===========================

function validate_window(value,name)
% Confirm that an analysis window has the form [start end] with end > start.

if ~isnumeric(value) || ...
        numel(value) ~= 2 || ...
        any(~isfinite(value)) || ...
        value(2) <= value(1)

    error('%s must be [start end], with end > start.',name);
end
end

function text_value = channel_list(values)
% Convert channel indices into a readable label.

if isempty(values)
    text_value = '(none)';
else
    labels = arrayfun(@(x) sprintf('Ch%d',x), ...
        values,'UniformOutput',false);

    text_value = strjoin(labels,' + ');
end
end