%% ========================================================================
% MULTI-ISI TRIAL SPIKE-COUNT TABLE
%
% PURPOSE
%   Extract per-channel spike counts for every original trial in a paired
%   simultaneous/sequential stimulation dataset.
%
% CHANNEL BEHAVIOR
%   - Channel_Results stores every recording channel.
%   - Trial-level totals and maximum values use only Count_Channels.
%
% OUTPUT
%   <base_name>_MultiISI_TrialSpikeCounts_PerTrial.mat
% ========================================================================

clear;
% close all;

addpath(genpath( ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';

% Electrode type:
%   0 = rigid single-shank
%   1 = flexible single-shank
%   2 = four-shank flexible
Electrode_Type = 2;

FS = 30000;

% Counting windows relative to the first stimulation trigger
baseline_win_ms = [-50 -5];
response_win_ms = [2 20];
analysis_win_ms = [2 40];

% Channels used for trial-level totals and maximum values.
%
% Empty means use all recording channels.
%
% Examples:
%   Count_Channels = 35:64;
%   Count_Channels = [35:40 42:48 50:64];
Count_Channels = [];

Create_Backup_If_Output_Exists = true;

%% =========================== INITIAL CHECKS ===========================

if ~isfolder(data_folder)
    error('Dataset folder does not exist:\n%s',data_folder);
end

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(response_win_ms,'response_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');

if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('FS must be a positive sampling rate.');
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

cd(data_folder);

baseline_duration_ms = diff(baseline_win_ms);
response_duration_ms = diff(response_win_ms);
analysis_duration_ms = diff(analysis_win_ms);

fprintf('\n============================================================\n');
fprintf('MULTI-ISI TRIAL SPIKE-COUNT TABLE\n');
fprintf('============================================================\n');
fprintf('Dataset: %s\n',data_folder);
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Response window: [%g,%g) ms\n',response_win_ms);
fprintf('Analysis window: [%g,%g) ms\n',analysis_win_ms);
fprintf('Responding-channel file used: NO\n');
fprintf('Bad-channel file used: NO\n');
fprintf('Bad-trial file used: NO\n');

%% ============================ LOAD sp_corr ============================

ssd_files = dir('*.sp_xia_SSD.mat');
ssd_files = remove_metadata_and_backups(ssd_files);

if isempty(ssd_files)
    error('No *.sp_xia_SSD.mat file was found.');
end

if numel(ssd_files) > 1
    error('Multiple current *.sp_xia_SSD.mat files were found.');
end

ssd_file = ssd_files(1).name;
base_name = erase(ssd_file,'.sp_xia_SSD.mat');

if ~ismember('sp_corr',who('-file',ssd_file))
    error('The SSD file does not contain sp_corr:\n%s',ssd_file);
end

SpikeLoad = load(ssd_file,'sp_corr');
sp = SpikeLoad.sp_corr;
nSpChannels = numel(sp);

if ~iscell(sp)
    error('sp_corr must be a cell array.');
end

fprintf('\nSpike file: %s\n',ssd_file);
fprintf('Spike variable: sp_corr\n');
fprintf('Spike-data cells: %d\n',nSpChannels);

%% =========================== LOAD TRIGGERS ============================

if isempty(dir('*.trig.dat'))
    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);
end

trigger_files = dir('*.trig.dat');

if isempty(trigger_files)
    error('No *.trig.dat file could be found or created.');
end

trig = double(loadTrig(0));
trig = trig(:);
nTrig = numel(trig);

fprintf('Trigger file: %s\n',trigger_files(1).name);
fprintf('Triggers loaded: %d\n',nTrig);

%% ==================== LOAD EXPERIMENT PARAMETERS =====================

experiment_files = dir('*_exp_datafile_*.mat');
experiment_files = remove_metadata_and_backups(experiment_files);

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found.');
end

if numel(experiment_files) > 1
    error('Multiple current experiment files were found.');
end

experiment_file = experiment_files(1).name;

ExpLoad = load(experiment_file, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');

required_variables = { ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials'};

for variable_index = 1:numel(required_variables)

    variable_name = required_variables{variable_index};

    if ~isfield(ExpLoad,variable_name)
        error('Experiment file is missing variable: %s', ...
            variable_name);
    end
end

StimParams = ExpLoad.StimParams;
sim_stim = double(ExpLoad.simultaneous_stim);
E_MAP = ExpLoad.E_MAP;
n_Trials = double(ExpLoad.n_Trials);

if sim_stim ~= 2
    error(['This script requires paired stimulation, but ' ...
        'simultaneous_stim = %d.'],sim_stim);
end

if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.', ...
        nTrig,n_Trials);

elseif nTrig > n_Trials
    warning('%d triggers loaded for %d trials; extras are ignored.', ...
        nTrig,n_Trials);

    trig = trig(1:n_Trials);
end

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Trials used: %d\n',n_Trials);

%% =========================== DECODE AMPLITUDES ========================

amplitudes_all = cell2mat(StimParams(2:end,16));
amplitudes_all = double(amplitudes_all(:));

first_event_rows = 1:sim_stim:numel(amplitudes_all);
second_event_rows = 2:sim_stim:numel(amplitudes_all);

trialAmps = amplitudes_all(first_event_rows);
secondPulseAmps = amplitudes_all(second_event_rows);

trialAmps = trialAmps(1:n_Trials);
secondPulseAmps = secondPulseAmps(1:n_Trials);

trialAmps(trialAmps == -1) = 0;
secondPulseAmps(secondPulseAmps == -1) = 0;

if any(abs(trialAmps-secondPulseAmps) > 1e-6)
    warning(['Some trials contain different first- and second-pulse ' ...
        'amplitudes. Conditions use the first-pulse amplitude.']);
end

[Amps,~,ampIdx] = unique(trialAmps(:));
nAMP = numel(Amps);

%% ============================== DECODE PTDs ===========================

PTD_all_us = cell2mat(StimParams(2:end,6));
PTD_all_us = double(PTD_all_us(:));

trialPTD_us = PTD_all_us(second_event_rows);
trialPTD_us = trialPTD_us(1:n_Trials);

trialPTD_ms = trialPTD_us/1000;

[PTDs,~,ptdIdx] = unique(trialPTD_us(:));
PTDs_ms = PTDs/1000;
nPTD = numel(PTDs);

%% ==================== DECODE STIMULATION ORDERS ======================

stimNames = StimParams(2:end,1);
stimNames = stimNames(1:n_Trials*sim_stim);

[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));

if any(~isMapped)
    error('%d stimulation entries could not be mapped.', ...
        sum(~isMapped));
end

stimSeq = zeros(n_Trials,sim_stim);

for trial_id = 1:n_Trials

    rows_this_trial = ...
        (trial_id-1)*sim_stim+(1:sim_stim);

    mapped_channels = idx_all(rows_this_trial);
    mapped_channels = mapped_channels(mapped_channels > 0);

    stimSeq(trial_id,1:numel(mapped_channels)) = ...
        mapped_channels(:).';
end

% Preserve stimulation order
[uniqueComb,~,combClass] = unique( ...
    stimSeq,'rows','stable');

nSets = size(uniqueComb,1);

fprintf('\nAmplitudes: %s uA\n',num2str(Amps(:).'));
fprintf('PTDs: %s ms\n',num2str(PTDs_ms(:).'));
fprintf('Ordered stimulation sets: %d\n',nSets);

for si = 1:nSets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    fprintf('  Set %d: %s\n', ...
        si,format_stimulation_order(stim_channels));
end

%% ======================= AUTOMATIC CHANNEL MAPPING ====================

d = double(Depth_s(Electrode_Type));
d = d(:);

if isempty(d)
    error('Depth_s returned an empty channel map.');
end

% All recording channels remain available in Channel_Results.
selected_channels = 1:numel(d);
nSelectedChannels = numel(selected_channels);

spike_data_channels = d(selected_channels);

valid_channel_mapping = ...
    isfinite(spike_data_channels) & ...
    spike_data_channels >= 1 & ...
    spike_data_channels <= nSpChannels & ...
    fix(spike_data_channels) == spike_data_channels;

if any(~valid_channel_mapping)
    warning(['%d channel indices map outside the available sp_corr ' ...
        'cells. Their counts will remain zero.'], ...
        sum(~valid_channel_mapping));
end

%% ====================== SELECT CHANNELS FOR TOTALS ====================

if isempty(Count_Channels)

    count_channels = selected_channels;

else
    count_channels = unique( ...
        double(Count_Channels(:).'),'stable');

    valid_count_channels = ...
        isfinite(count_channels) & ...
        count_channels >= 1 & ...
        count_channels <= numel(d) & ...
        fix(count_channels) == count_channels;

    if any(~valid_count_channels)

        warning(['Invalid Count_Channels were ignored: %s'], ...
            num2str(count_channels(~valid_count_channels)));

        count_channels = ...
            count_channels(valid_count_channels);
    end
end

if isempty(count_channels)
    error('No valid Count_Channels remain.');
end

% Convert displayed channel indices into positions in the count arrays
[found_count_channels,count_channel_positions] = ...
    ismember(count_channels,selected_channels);

if any(~found_count_channels)
    error('Some Count_Channels could not be mapped.');
end

nCountChannels = numel(count_channels);

fprintf('\nAll channels stored in Channel_Results: 1:%d\n', ...
    nSelectedChannels);

fprintf('Channels used for trial totals: %s\n', ...
    number_list(count_channels));

fprintf('Number of channels used for totals: %d\n', ...
    nCountChannels);

%% ====================== CACHE ALL SPIKE TIMES =========================

selected_spike_times = cell(nSelectedChannels,1);
channel_has_spike_data = false(nSelectedChannels,1);

for channel_position = 1:nSelectedChannels

    if ~valid_channel_mapping(channel_position)
        continue;
    end

    spike_channel = spike_data_channels(channel_position);

    if isempty(sp{spike_channel})
        continue;
    end

    selected_spike_times{channel_position} = ...
        double(sp{spike_channel}(:,1));

    channel_has_spike_data(channel_position) = true;
end

fprintf('Channels containing spike data: %d/%d\n', ...
    sum(channel_has_spike_data),nSelectedChannels);

%% ======================== INITIALIZE OUTPUT ===========================

SpikeCounts = struct();

SpikeCounts.metadata.analysis_name = ...
    'Multi-ISI trial spike-count table';

SpikeCounts.metadata.created = char(datetime( ...
    'now','Format','yyyy-MM-dd HH:mm:ss'));

SpikeCounts.metadata.data_folder = data_folder;
SpikeCounts.metadata.spike_file = ssd_file;
SpikeCounts.metadata.spike_variable = 'sp_corr';
SpikeCounts.metadata.experiment_file = experiment_file;
SpikeCounts.metadata.sampling_rate_hz = FS;
SpikeCounts.metadata.electrode_type = Electrode_Type;

SpikeCounts.metadata.baseline_window_ms = baseline_win_ms;
SpikeCounts.metadata.response_window_ms = response_win_ms;
SpikeCounts.metadata.analysis_window_ms = analysis_win_ms;

SpikeCounts.metadata.responding_channels_used = false;
SpikeCounts.metadata.bad_channels_removed = false;
SpikeCounts.metadata.bad_trials_removed = false;
SpikeCounts.metadata.trial_summary_sorted = false;
SpikeCounts.metadata.order_sensitive_sets = true;

SpikeCounts.metadata.amplitudes_uA = Amps;
SpikeCounts.metadata.PTDs_us = PTDs;
SpikeCounts.metadata.PTDs_ms = PTDs_ms;
SpikeCounts.metadata.unique_stimulation_orders = uniqueComb;

% selected_channels contains every stored channel.
SpikeCounts.metadata.selected_channels = selected_channels;

% count_channels contains only channels used for trial-level totals.
SpikeCounts.metadata.count_channels = count_channels;
SpikeCounts.metadata.n_count_channels = nCountChannels;

%% ==================== BUILD CONDITION STRUCTURE =======================

for si = 1:nSets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);
    set_label = format_stimulation_order(stim_channels);

    SpikeCounts.set(si).set_index = si;
    SpikeCounts.set(si).set_name = set_label;
    SpikeCounts.set(si).stim_channels = stim_channels;

    for ai = 1:nAMP

        SpikeCounts.set(si).amp(ai).amp_index = ai;
        SpikeCounts.set(si).amp(ai).amp_value = Amps(ai);
        SpikeCounts.set(si).amp(ai).amp_label = ...
            sprintf('%g uA',Amps(ai));

        for pi = 1:nPTD

            condition_trials = find( ...
                combClass == si & ...
                ampIdx == ai & ...
                ptdIdx == pi);

            SpikeCounts.set(si).amp(ai).ptd(pi).ptd_index = pi;
            SpikeCounts.set(si).amp(ai).ptd(pi).PTD_us = PTDs(pi);
            SpikeCounts.set(si).amp(ai).ptd(pi).PTD_ms = PTDs_ms(pi);

            SpikeCounts.set(si).amp(ai).ptd(pi).condition_summary = ...
                sprintf('Set %s | Amp %g uA | PTD %g ms', ...
                set_label,Amps(ai),PTDs_ms(pi));

            SpikeCounts.set(si).amp(ai).ptd(pi).trial_ids = ...
                condition_trials(:).';

            SpikeCounts.set(si).amp(ai).ptd(pi).total_trials = ...
                numel(condition_trials);
        end
    end
end

%% ================= INITIALIZE FLAT TRIAL SUMMARY ======================

AbsoluteTrialID = (1:n_Trials).';
OriginalTrialSequence = AbsoluteTrialID;
ConditionTrialIndex = zeros(n_Trials,1);

SetIndex = combClass(:);
FirstStimChannel = stimSeq(:,1);
SecondStimChannel = stimSeq(:,2);

AmplitudeIndex = ampIdx(:);
Amplitude_uA = trialAmps(:);

PTDIndex = ptdIdx(:);
PTD_us = trialPTD_us(:);
PTD_ms = trialPTD_ms(:);

% Number of channels actually used for trial-level totals
N_CountChannels = repmat(nCountChannels,n_Trials,1);

TotalBaselineSpikes = zeros(n_Trials,1);
TotalResponseWindowSpikes = zeros(n_Trials,1);
TotalAnalysisWindowSpikes = zeros(n_Trials,1);

TotalBaselineCorrectedResponseWindowSpikes = zeros(n_Trials,1);
TotalBaselineCorrectedAnalysisWindowSpikes = zeros(n_Trials,1);

MaxBaselineSpikes_OneChannel = zeros(n_Trials,1);
MaxResponseSpikes_OneChannel = zeros(n_Trials,1);
MaxAnalysisSpikes_OneChannel = zeros(n_Trials,1);

condition_counters = zeros(nSets,nAMP,nPTD);

%% ================= COUNT IN ORIGINAL TRIAL ORDER ======================

fprintf('\nCounting spikes in original trial order...\n');

for trial_id = 1:n_Trials

    si = combClass(trial_id);
    ai = ampIdx(trial_id);
    pi = ptdIdx(trial_id);

    condition_counters(si,ai,pi) = ...
        condition_counters(si,ai,pi)+1;

    condition_trial_index = ...
        condition_counters(si,ai,pi);

    ConditionTrialIndex(trial_id) = ...
        condition_trial_index;

    trigger_ms = trig(trial_id)/FS*1000;

    baseline_absolute = trigger_ms+baseline_win_ms;
    response_absolute = trigger_ms+response_win_ms;
    analysis_absolute = trigger_ms+analysis_win_ms;

    % Counts remain available for every recording channel
    baseline_counts = zeros(nSelectedChannels,1);
    response_counts = zeros(nSelectedChannels,1);
    analysis_counts = zeros(nSelectedChannels,1);

    %% ---------------- COUNT EVERY RECORDING CHANNEL -------------------

    for channel_position = 1:nSelectedChannels

        if ~channel_has_spike_data(channel_position)
            continue;
        end

        spike_times = ...
            selected_spike_times{channel_position};

        baseline_counts(channel_position) = sum( ...
            spike_times >= baseline_absolute(1) & ...
            spike_times <  baseline_absolute(2));

        response_counts(channel_position) = sum( ...
            spike_times >= response_absolute(1) & ...
            spike_times <  response_absolute(2));

        analysis_counts(channel_position) = sum( ...
            spike_times >= analysis_absolute(1) & ...
            spike_times <  analysis_absolute(2));
    end

    %% ---------------- SIGNED BASELINE CORRECTION ----------------------

    expected_baseline_response = baseline_counts * ...
        (response_duration_ms/baseline_duration_ms);

    expected_baseline_analysis = baseline_counts * ...
        (analysis_duration_ms/baseline_duration_ms);

    corrected_response_counts = ...
        response_counts-expected_baseline_response;

    corrected_analysis_counts = ...
        analysis_counts-expected_baseline_analysis;

    %% ---------------- PER-CHANNEL TABLE -------------------------------

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

    %% ---------------- SELECT COUNTS FOR TOTALS ------------------------

    summary_baseline_counts = ...
        baseline_counts(count_channel_positions);

    summary_response_counts = ...
        response_counts(count_channel_positions);

    summary_analysis_counts = ...
        analysis_counts(count_channel_positions);

    summary_corrected_response_counts = ...
        corrected_response_counts(count_channel_positions);

    summary_corrected_analysis_counts = ...
        corrected_analysis_counts(count_channel_positions);

    %% ---------------- STORE INDIVIDUAL TRIAL --------------------------

    T = struct();

    T.relative_trial_ID = condition_trial_index;
    T.ConditionTrialIndex = condition_trial_index;

    T.absolute_trial_ID = trial_id;
    T.AbsoluteTrialID = trial_id;
    T.OriginalTrialSequence = trial_id;

    T.SetIndex = si;
    T.FirstStimChannel = stimSeq(trial_id,1);
    T.SecondStimChannel = stimSeq(trial_id,2);
    T.Amplitude_uA = trialAmps(trial_id);
    T.PTD_us = trialPTD_us(trial_id);
    T.PTD_ms = trialPTD_ms(trial_id);

    T.Channel_Results = Channel_Results;

    % Record which channels were used for these totals
    T.count_channels = count_channels;
    T.n_count_channels = nCountChannels;

    % Trial totals use Count_Channels only
    T.total_baseline_spikes = ...
        sum(summary_baseline_counts);

    T.total_response_window_spikes = ...
        sum(summary_response_counts);

    T.total_analysis_window_spikes = ...
        sum(summary_analysis_counts);

    T.total_baseline_corrected_response_window_spikes = ...
        sum(summary_corrected_response_counts);

    T.total_baseline_corrected_analysis_window_spikes = ...
        sum(summary_corrected_analysis_counts);

    % The duplicate T.total_post_spikes field has been removed.

    SpikeCounts.set(si).amp(ai).ptd(pi).trial( ...
        condition_trial_index) = T;

    %% ---------------- UPDATE FLAT TRIAL SUMMARY -----------------------

    TotalBaselineSpikes(trial_id) = ...
        sum(summary_baseline_counts);

    TotalResponseWindowSpikes(trial_id) = ...
        sum(summary_response_counts);

    TotalAnalysisWindowSpikes(trial_id) = ...
        sum(summary_analysis_counts);

    TotalBaselineCorrectedResponseWindowSpikes(trial_id) = ...
        sum(summary_corrected_response_counts);

    TotalBaselineCorrectedAnalysisWindowSpikes(trial_id) = ...
        sum(summary_corrected_analysis_counts);

    MaxBaselineSpikes_OneChannel(trial_id) = ...
        max(summary_baseline_counts);

    MaxResponseSpikes_OneChannel(trial_id) = ...
        max(summary_response_counts);

    MaxAnalysisSpikes_OneChannel(trial_id) = ...
        max(summary_analysis_counts);
end

%% ======================== BUILD SUMMARY TABLE =========================

% Row n remains original experiment trial n.
TrialSummary = table( ...
    AbsoluteTrialID, ...
    OriginalTrialSequence, ...
    ConditionTrialIndex, ...
    SetIndex, ...
    FirstStimChannel, ...
    SecondStimChannel, ...
    AmplitudeIndex, ...
    Amplitude_uA, ...
    PTDIndex, ...
    PTD_us, ...
    PTD_ms, ...
    N_CountChannels, ...
    TotalBaselineSpikes, ...
    TotalResponseWindowSpikes, ...
    TotalAnalysisWindowSpikes, ...
    TotalBaselineCorrectedResponseWindowSpikes, ...
    TotalBaselineCorrectedAnalysisWindowSpikes, ...
    MaxBaselineSpikes_OneChannel, ...
    MaxResponseSpikes_OneChannel, ...
    MaxAnalysisSpikes_OneChannel);

%% ======================== PRINT CONDITION SUMMARY =====================

fprintf('\nCondition summary:\n');

for si = 1:nSets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);
    set_label = format_stimulation_order(stim_channels);

    for ai = 1:nAMP
        for pi = 1:nPTD

            nConditionTrials = ...
                SpikeCounts.set(si).amp(ai).ptd(pi).total_trials;

            if nConditionTrials == 0
                continue;
            end

            fprintf(['  Set %d | %s | Amp %g uA | PTD %g ms | ' ...
                'Trials %d\n'], ...
                si, ...
                set_label, ...
                Amps(ai), ...
                PTDs_ms(pi), ...
                nConditionTrials);
        end
    end
end

%% =========================== SAVE RESULT ==============================

save_name = sprintf( ...
    '%s_MultiISI_TrialSpikeCounts_PerTrial.mat',base_name);

full_save_path = fullfile(data_folder,save_name);

if isfile(full_save_path) && ...
        Create_Backup_If_Output_Exists

    timestamp = datestr(now,'yyyymmdd_HHMMSS');

    backup_name = sprintf( ...
        '%s_MultiISI_TrialSpikeCounts_PerTrial_BACKUP_%s.mat', ...
        base_name,timestamp);

    backup_path = fullfile(data_folder,backup_name);

    copyfile(full_save_path,backup_path);

    fprintf('\nExisting output backed up to:\n%s\n', ...
        backup_path);
end

% Retained only as a window-name compatibility alias
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
    'PTDs_ms', ...
    'uniqueComb', ...
    'selected_channels', ...
    'count_channels', ...
    'FS', ...
    'Electrode_Type', ...
    '-v7.3');

fprintf('\n============================================================\n');
fprintf('MULTI-ISI SPIKE-COUNT EXTRACTION COMPLETE\n');
fprintf('Trials saved: %d\n',height(TrialSummary));
fprintf('Channels stored per trial: %d\n',nSelectedChannels);
fprintf('Channels used for totals: %d\n',nCountChannels);
fprintf('Count channels: %s\n',number_list(count_channels));
fprintf('TrialSummary remains in original trial order: YES\n');
fprintf('Output file:\n%s\n',full_save_path);
fprintf('Original experiment files modified: NO\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function validate_window(window_ms,variable_name)

if ~isnumeric(window_ms) || ...
        numel(window_ms) ~= 2 || ...
        any(~isfinite(window_ms)) || ...
        window_ms(2) <= window_ms(1)

    error('%s must be [start end], with end > start.', ...
        variable_name);
end
end

function files = remove_metadata_and_backups(files)

if isempty(files)
    return;
end

keep = true(size(files));

for file_index = 1:numel(files)

    file_name = files(file_index).name;

    if startsWith(file_name,'._') || ...
            contains(file_name,'BACKUP','IgnoreCase',true)

        keep(file_index) = false;
    end
end

files = files(keep);
end

function label = format_stimulation_order(stim_channels)

stim_channels = stim_channels(stim_channels > 0);

if isempty(stim_channels)

    label = '(none)';

elseif numel(stim_channels) == 1

    label = sprintf('Ch%d',stim_channels(1));

else
    labels = arrayfun(@(channel_number) ...
        sprintf('Ch%d',channel_number), ...
        stim_channels, ...
        'UniformOutput',false);

    label = strjoin(labels,' -> ');
end
end

function output = number_list(values)

if isempty(values)
    output = '(none)';
else
    output = strtrim(sprintf('%g ',values));
end
end