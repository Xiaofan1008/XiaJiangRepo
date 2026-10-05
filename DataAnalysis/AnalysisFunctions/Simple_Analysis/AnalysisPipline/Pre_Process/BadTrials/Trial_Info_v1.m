%% ========================================================================
% TRIAL SPIKE-COUNT TABLE
%
% PURPOSE
%   Extract per-trial, per-channel spike counts for every original trial
%   in this dataset (single-pulse, simultaneous, or sequential paired
%   stimulation). This is the input the bad-trial reference/manual-confirm
%   tools read from - it does not decide anything about bad trials
%   itself.
%
% CHANNEL BEHAVIOR
%   - Channel_Results (stored per trial) covers every map-position
%     channel.
%   - Trial-level totals (used later for bad-trial flagging) use only
%     Count_Channels - leave Count_Channels empty to use every channel, or
%     specify a subset.
%
% WINDOWS
%   Just two: a baseline window (for the signed baseline correction) and
%   an analysis window (the post-stimulation count that bad-trial
%   flagging actually reads). No separate "response" window.
%
% INPUT CONVENTION
%   Same as RespondingChn_Raster_AllChn_v2.m / RespondingChn_ReferenceCriteria_v1.m:
%   sp_clipped / sp_waveforms (physical-channel indexed), ChnMap(Electrode_Type, nCh)
%   for the map-position -> physical-channel lookup, StimParams/E_MAP/n_Trials
%   from the dataset's own *_exp_datafile_*.mat for condition decoding.
%
% OUTPUT
%   <base_name>_TrialSpikeCounts.mat, in dataset_folder
% ========================================================================

clear;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_stitched.mat';
dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1';

% Electrode type (must match this dataset - probe differs by animal):
%   0 = rigid single-shank probe
%   1 = flexible single-shank probe
%   2 = four-shank flexible probe
%   3 = 64chn+32chn hybrid
Electrode_Type = 3;

FS_fallback = 30000;   % used only if fs/FS is absent from spike_file

% Counting windows relative to the first stimulation trigger (t=0).
% Upper boundary excluded.
baseline_win_ms = [-50 -5];
analysis_win_ms = [2 20];

% Channels used for TRIAL-LEVEL TOTALS (what the bad-trial tools read).
% Empty = use every map-position channel. Or specify a subset, e.g.
% Count_Channels = [1:16 33:48];
% Channel_Results (stored per trial, for every channel) is unaffected
% either way.
Count_Channels = [];

Create_Backup_If_Output_Exists = true;

%% =========================== INITIAL CHECKS ===========================

if ~isfile(spike_file)
    error('Spike file does not exist:\n%s',spike_file);
end
if ~isfolder(dataset_folder)
    error('Dataset folder does not exist:\n%s',dataset_folder);
end

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

baseline_duration_ms = diff(baseline_win_ms);
analysis_duration_ms = diff(analysis_win_ms);

fprintf('\n============================================================\n');
fprintf('TRIAL SPIKE-COUNT TABLE\n');
fprintf('============================================================\n');
fprintf('Spike file: %s\n',spike_file);
fprintf('Dataset folder: %s\n',dataset_folder);
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Analysis window: [%g,%g) ms\n',analysis_win_ms);

%% ========================= LOAD SPIKE TIMES ============================

spike_vars_present = who('-file',spike_file);

if ismember('sp_clipped',spike_vars_present)
    SpikeLoad = load(spike_file,'sp_clipped');
    sp_raw = SpikeLoad.sp_clipped;
    spike_source_var = 'sp_clipped';
elseif ismember('sp_waveforms',spike_vars_present)
    SpikeLoad = load(spike_file,'sp_waveforms');
    sp_raw = SpikeLoad.sp_waveforms;
    spike_source_var = 'sp_waveforms';
else
    error('Spike file contains neither sp_clipped nor sp_waveforms:\n%s',spike_file);
end

nCh = numel(sp_raw);

sp = cell(nCh,1);
for ch = 1:nCh
    cellData = sp_raw{ch};
    if isempty(cellData)
        sp{ch} = [];
    else
        spike_times = double(cellData(:,1));
        spike_times = spike_times(isfinite(spike_times));
        if ~issorted(spike_times)
            spike_times = sort(spike_times);
        end
        sp{ch} = spike_times;
    end
end
clear sp_raw SpikeLoad;

if ismember('fs',spike_vars_present)
    FsLoad = load(spike_file,'fs');
    FS = double(FsLoad.fs);
elseif ismember('FS',spike_vars_present)
    FsLoad = load(spike_file,'FS');
    FS = double(FsLoad.FS);
else
    warning('fs/FS not found in spike file - using fallback FS = %g Hz.',FS_fallback);
    FS = FS_fallback;
end

fprintf('\nSpike variable used: %s\n',spike_source_var);
fprintf('Physical channels in spike file: %d\n',nCh);
fprintf('Sampling rate: %g Hz\n',FS);

%% ========================= ELECTRODE MAPPING ===========================

map_nums_plus = double(ChnMap(Electrode_Type,nCh));
map_nums_plus = map_nums_plus(:);

if isempty(map_nums_plus)
    error('ChnMap(%d,%d) returned an empty channel map.',Electrode_Type,nCh);
end

nMapPos = numel(map_nums_plus);
fprintf('Map positions: %d\n',nMapPos);

selected_channels = 1:nMapPos;
nSelectedChannels = nMapPos;

valid_channel_mapping = isfinite(map_nums_plus) & map_nums_plus >= 1 & ...
    map_nums_plus <= nCh & fix(map_nums_plus) == map_nums_plus;

if any(~valid_channel_mapping)
    warning(['%d map positions map outside the available spike-data channels. ' ...
        'Their counts will remain zero.'],sum(~valid_channel_mapping));
end

%% =========================== LOAD TRIGGERS ============================

cd(dataset_folder);

if isempty(dir('*.trig.dat'))
    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);
end

trigger_files = dir('*.trig.dat');
if isempty(trigger_files)
    error('No *.trig.dat file could be found or created in:\n%s',dataset_folder);
end

trig = double(loadTrig(0));
trig = trig(:);
nTrig = numel(trig);

fprintf('Triggers loaded: %d\n',nTrig);

%% ==================== LOAD EXPERIMENT PARAMETERS =====================

experiment_files = dir('*_exp_datafile_*.mat');
experiment_files = remove_metadata_and_backups(experiment_files);

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found in:\n%s',dataset_folder);
end
if numel(experiment_files) > 1
    error('Multiple current experiment files were found in:\n%s',dataset_folder);
end

experiment_file = experiment_files(1).name;

ExpLoad = load(experiment_file,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required_variables = {'StimParams','simultaneous_stim','E_MAP','n_Trials'};
for variable_index = 1:numel(required_variables)
    if ~isfield(ExpLoad,required_variables{variable_index})
        error('Experiment file is missing variable: %s',required_variables{variable_index});
    end
end

StimParams = ExpLoad.StimParams;
sim_stim = double(ExpLoad.simultaneous_stim);
E_MAP = ExpLoad.E_MAP;
n_Trials = double(ExpLoad.n_Trials);

if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.',nTrig,n_Trials);
elseif nTrig > n_Trials
    warning('%d triggers loaded for %d trials; extras are ignored.',nTrig,n_Trials);
    trig = trig(1:n_Trials);
end

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Stimulation events per trial (simultaneous_stim): %d\n',sim_stim);
fprintf('Trials used: %d\n',n_Trials);

%% =========================== DECODE AMPLITUDES ========================

amplitudes_all = cell2mat(StimParams(2:end,16));
amplitudes_all = double(amplitudes_all(:));

first_event_rows = 1:sim_stim:numel(amplitudes_all);
trialAmps = amplitudes_all(first_event_rows);
trialAmps = trialAmps(1:n_Trials);
trialAmps(trialAmps == -1) = 0;

if sim_stim >= 2
    second_event_rows = 2:sim_stim:numel(amplitudes_all);
    secondPulseAmps = amplitudes_all(second_event_rows);
    secondPulseAmps = secondPulseAmps(1:n_Trials);
    secondPulseAmps(secondPulseAmps == -1) = 0;

    if any(abs(trialAmps-secondPulseAmps) > 1e-6)
        warning(['Some trials contain different first- and second-pulse ' ...
            'amplitudes. Conditions use the first-pulse amplitude.']);
    end
end

[Amps,~,ampIdx] = unique(trialAmps(:));
nAMP = numel(Amps);

%% ============================== DECODE PTDs =============================

if sim_stim >= 2
    PTD_all_us = cell2mat(StimParams(2:end,6));
    PTD_all_us = double(PTD_all_us(:));

    second_event_rows = 2:sim_stim:numel(PTD_all_us);
    trialPTD_us = PTD_all_us(second_event_rows);
    trialPTD_us = trialPTD_us(1:n_Trials);

    [PTDs,~,ptdIdx] = unique(trialPTD_us(:));
    PTDs_ms = PTDs/1000;
    nPTD = numel(PTDs);
else
    trialPTD_us = zeros(n_Trials,1);
    PTDs = 0;
    PTDs_ms = 0;
    nPTD = 1;
    ptdIdx = ones(n_Trials,1);
end

trialPTD_ms = trialPTD_us/1000;

%% ==================== DECODE STIMULATION CHANNEL ORDERS =================

stimNames = StimParams(2:end,1);
stimNames = stimNames(1:n_Trials*sim_stim);

[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));
if any(~isMapped)
    error('%d stimulation entries could not be mapped through E_MAP.',sum(~isMapped));
end

stimSeq = zeros(n_Trials,sim_stim);
for trial_id = 1:n_Trials
    rows_this_trial = (trial_id-1)*sim_stim+(1:sim_stim);
    mapped_channels = idx_all(rows_this_trial);
    mapped_channels = mapped_channels(mapped_channels > 0);
    stimSeq(trial_id,1:numel(mapped_channels)) = mapped_channels(:).';
end

[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

fprintf('\nAmplitudes: %s uA\n',num2str(Amps(:).'));
if sim_stim >= 2
    fprintf('PTDs: %s ms\n',num2str(PTDs_ms(:).'));
end
fprintf('Ordered stimulation sets: %d\n',nSets);
for si = 1:nSets
    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);
    fprintf('  Set %d: %s\n',si,format_order(stim_channels));
end

FirstStimChannel = stimSeq(:,1);
if sim_stim >= 2
    SecondStimChannel = stimSeq(:,2);
else
    SecondStimChannel = zeros(n_Trials,1);
end

%% ====================== SELECT CHANNELS FOR TOTALS ====================

if isempty(Count_Channels)
    count_channels = selected_channels;
else
    count_channels = unique(double(Count_Channels(:).'),'stable');
    valid_count_channels = isfinite(count_channels) & count_channels >= 1 & ...
        count_channels <= nMapPos & fix(count_channels) == count_channels;

    if any(~valid_count_channels)
        warning('Invalid Count_Channels were ignored: %s', ...
            num2str(count_channels(~valid_count_channels)));
        count_channels = count_channels(valid_count_channels);
    end
end

if isempty(count_channels)
    error('No valid Count_Channels remain.');
end

[found_count_channels,count_channel_positions] = ismember(count_channels,selected_channels);
if any(~found_count_channels)
    error('Some Count_Channels could not be mapped.');
end

nCountChannels = numel(count_channels);

fprintf('\nAll channels stored in Channel_Results: 1:%d\n',nSelectedChannels);
fprintf('Channels used for trial totals: %s\n',number_list(count_channels));
fprintf('Number of channels used for totals: %d\n',nCountChannels);

%% ====================== CACHE ALL SPIKE TIMES =========================

selected_spike_times = cell(nSelectedChannels,1);
channel_has_spike_data = false(nSelectedChannels,1);

for channel_position = 1:nSelectedChannels

    if ~valid_channel_mapping(channel_position)
        continue;
    end

    phys_ch = map_nums_plus(channel_position);

    if isempty(sp{phys_ch})
        continue;
    end

    selected_spike_times{channel_position} = sp{phys_ch};
    channel_has_spike_data(channel_position) = true;
end

fprintf('Channels containing spike data: %d/%d\n', ...
    sum(channel_has_spike_data),nSelectedChannels);

clear sp;

%% ======================== INITIALIZE OUTPUT ===========================

SpikeCounts = struct();

SpikeCounts.metadata.analysis_name = 'Trial spike-count table';
SpikeCounts.metadata.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
SpikeCounts.metadata.spike_file = spike_file;
SpikeCounts.metadata.spike_source_var = spike_source_var;
SpikeCounts.metadata.dataset_folder = dataset_folder;
SpikeCounts.metadata.experiment_file = experiment_file;
SpikeCounts.metadata.sampling_rate_hz = FS;
SpikeCounts.metadata.Electrode_Type = Electrode_Type;

SpikeCounts.metadata.baseline_window_ms = baseline_win_ms;
SpikeCounts.metadata.analysis_window_ms = analysis_win_ms;

SpikeCounts.metadata.amplitudes_uA = Amps;
SpikeCounts.metadata.PTDs_us = PTDs;
SpikeCounts.metadata.PTDs_ms = PTDs_ms;
SpikeCounts.metadata.unique_stimulation_orders = uniqueComb;

SpikeCounts.metadata.selected_channels = selected_channels;
SpikeCounts.metadata.count_channels = count_channels;
SpikeCounts.metadata.n_count_channels = nCountChannels;

%% ==================== BUILD CONDITION STRUCTURE =======================

for si = 1:nSets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);
    set_label = format_order(stim_channels);

    SpikeCounts.set(si).set_index = si;
    SpikeCounts.set(si).set_name = set_label;
    SpikeCounts.set(si).stim_channels = stim_channels;

    dist_um = NaN;
    dist_note = '';
    if numel(stim_channels) == 2
        [dist_um,dist_note] = ChnPairDistance(stim_channels(1),stim_channels(2),Electrode_Type);
    end
    SpikeCounts.set(si).distance_um = dist_um;
    SpikeCounts.set(si).distance_note = dist_note;

    for ai = 1:nAMP

        SpikeCounts.set(si).amp(ai).amp_index = ai;
        SpikeCounts.set(si).amp(ai).amp_value = Amps(ai);
        SpikeCounts.set(si).amp(ai).amp_label = sprintf('%g uA',Amps(ai));

        for pi = 1:nPTD

            condition_trials = find(combClass == si & ampIdx == ai & ptdIdx == pi);

            SpikeCounts.set(si).amp(ai).ptd(pi).ptd_index = pi;
            SpikeCounts.set(si).amp(ai).ptd(pi).PTD_us = PTDs(pi);
            SpikeCounts.set(si).amp(ai).ptd(pi).PTD_ms = PTDs_ms(pi);

            SpikeCounts.set(si).amp(ai).ptd(pi).condition_summary = ...
                sprintf('Set %s | Amp %g uA | PTD %g ms',set_label,Amps(ai),PTDs_ms(pi));

            SpikeCounts.set(si).amp(ai).ptd(pi).trial_ids = condition_trials(:).';
            SpikeCounts.set(si).amp(ai).ptd(pi).total_trials = numel(condition_trials);
        end
    end
end

%% ================= INITIALIZE FLAT TRIAL SUMMARY ======================

AbsoluteTrialID = (1:n_Trials).';
OriginalTrialSequence = AbsoluteTrialID;
ConditionTrialIndex = zeros(n_Trials,1);

SetIndex = combClass(:);

AmplitudeIndex = ampIdx(:);
Amplitude_uA = trialAmps(:);

PTDIndex = ptdIdx(:);
PTD_us = trialPTD_us(:);
PTD_ms = trialPTD_ms(:);

N_CountChannels = repmat(nCountChannels,n_Trials,1);

TotalBaselineSpikes = zeros(n_Trials,1);
TotalAnalysisWindowSpikes = zeros(n_Trials,1);
TotalBaselineCorrectedAnalysisWindowSpikes = zeros(n_Trials,1);
MaxBaselineSpikes_OneChannel = zeros(n_Trials,1);
MaxAnalysisSpikes_OneChannel = zeros(n_Trials,1);

condition_counters = zeros(nSets,nAMP,nPTD);

%% ================= COUNT IN ORIGINAL TRIAL ORDER ======================

fprintf('\nCounting spikes in original trial order...\n');

for trial_id = 1:n_Trials

    si = combClass(trial_id);
    ai = ampIdx(trial_id);
    pi = ptdIdx(trial_id);

    condition_counters(si,ai,pi) = condition_counters(si,ai,pi)+1;
    condition_trial_index = condition_counters(si,ai,pi);
    ConditionTrialIndex(trial_id) = condition_trial_index;

    trigger_ms = trig(trial_id)/FS*1000;

    baseline_absolute = trigger_ms+baseline_win_ms;
    analysis_absolute = trigger_ms+analysis_win_ms;

    baseline_counts = zeros(nSelectedChannels,1);
    analysis_counts = zeros(nSelectedChannels,1);

    for channel_position = 1:nSelectedChannels

        if ~channel_has_spike_data(channel_position)
            continue;
        end

        spike_times = selected_spike_times{channel_position};

        baseline_counts(channel_position) = sum( ...
            spike_times >= baseline_absolute(1) & spike_times < baseline_absolute(2));

        analysis_counts(channel_position) = sum( ...
            spike_times >= analysis_absolute(1) & spike_times < analysis_absolute(2));
    end

    %% ---------------- SIGNED BASELINE CORRECTION ----------------------

    expected_baseline_analysis = baseline_counts*(analysis_duration_ms/baseline_duration_ms);
    corrected_analysis_counts = analysis_counts-expected_baseline_analysis;

    %% ---------------- PER-CHANNEL TABLE -------------------------------

    Channel_Results = table(selected_channels(:),baseline_counts,analysis_counts, ...
        corrected_analysis_counts, ...
        'VariableNames',{'Channel_Index','Baseline_Spikes','AnalysisWindow_Spikes', ...
        'BaselineCorrected_AnalysisWindow_Spikes'});

    %% ---------------- SELECT COUNTS FOR TOTALS ------------------------

    summary_baseline_counts = baseline_counts(count_channel_positions);
    summary_analysis_counts = analysis_counts(count_channel_positions);
    summary_corrected_analysis_counts = corrected_analysis_counts(count_channel_positions);

    %% ---------------- STORE INDIVIDUAL TRIAL --------------------------

    T = struct();
    T.relative_trial_ID = condition_trial_index;
    T.ConditionTrialIndex = condition_trial_index;
    T.absolute_trial_ID = trial_id;
    T.AbsoluteTrialID = trial_id;
    T.OriginalTrialSequence = trial_id;

    T.SetIndex = si;
    T.FirstStimChannel = FirstStimChannel(trial_id);
    T.SecondStimChannel = SecondStimChannel(trial_id);
    T.Amplitude_uA = trialAmps(trial_id);
    T.PTD_us = trialPTD_us(trial_id);
    T.PTD_ms = trialPTD_ms(trial_id);

    T.Channel_Results = Channel_Results;
    T.count_channels = count_channels;
    T.n_count_channels = nCountChannels;

    T.total_baseline_spikes = sum(summary_baseline_counts);
    T.total_analysis_window_spikes = sum(summary_analysis_counts);
    T.total_baseline_corrected_analysis_window_spikes = sum(summary_corrected_analysis_counts);

    SpikeCounts.set(si).amp(ai).ptd(pi).trial(condition_trial_index) = T;

    %% ---------------- UPDATE FLAT TRIAL SUMMARY -----------------------

    TotalBaselineSpikes(trial_id) = sum(summary_baseline_counts);
    TotalAnalysisWindowSpikes(trial_id) = sum(summary_analysis_counts);
    TotalBaselineCorrectedAnalysisWindowSpikes(trial_id) = sum(summary_corrected_analysis_counts);

    MaxBaselineSpikes_OneChannel(trial_id) = max(summary_baseline_counts);
    MaxAnalysisSpikes_OneChannel(trial_id) = max(summary_analysis_counts);
end

%% ======================== BUILD SUMMARY TABLE =========================

TrialSummary = table(AbsoluteTrialID,OriginalTrialSequence,ConditionTrialIndex, ...
    SetIndex,FirstStimChannel,SecondStimChannel,AmplitudeIndex,Amplitude_uA, ...
    PTDIndex,PTD_us,PTD_ms,N_CountChannels,TotalBaselineSpikes, ...
    TotalAnalysisWindowSpikes,TotalBaselineCorrectedAnalysisWindowSpikes, ...
    MaxBaselineSpikes_OneChannel,MaxAnalysisSpikes_OneChannel);

%% ======================== PRINT CONDITION SUMMARY =====================

fprintf('\nCondition summary:\n');

for si = 1:nSets

    set_label = SpikeCounts.set(si).set_name;
    dist_um = SpikeCounts.set(si).distance_um;
    if isnan(dist_um)
        dist_text = sprintf('Dist n/a (%s)',SpikeCounts.set(si).distance_note);
    else
        dist_text = sprintf('Dist %g um',dist_um);
    end

    for ai = 1:nAMP
        for pi = 1:nPTD

            nConditionTrials = SpikeCounts.set(si).amp(ai).ptd(pi).total_trials;
            if nConditionTrials == 0
                continue;
            end

            fprintf('  Set %d | %s | Amp %g uA | PTD %g ms | Trials %d | %s\n', ...
                si,set_label,Amps(ai),PTDs_ms(pi),nConditionTrials,dist_text);
        end
    end
end

%% =========================== SAVE RESULT ==============================

output_name = sprintf('%s_TrialSpikeCounts.mat',regexprep(experiment_file,'_exp_datafile_.*$',''));
full_save_path = fullfile(dataset_folder,output_name);

if isfile(full_save_path) && Create_Backup_If_Output_Exists
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    backup_name = sprintf('%s_TrialSpikeCounts_BACKUP_%s.mat', ...
        regexprep(experiment_file,'_exp_datafile_.*$',''),timestamp);
    backup_path = fullfile(dataset_folder,backup_name);
    copyfile(full_save_path,backup_path);
    fprintf('\nExisting output backed up to:\n%s\n',backup_path);
end

save(full_save_path,'SpikeCounts','TrialSummary','baseline_win_ms', ...
    'analysis_win_ms','Amps','PTDs','PTDs_ms','uniqueComb','selected_channels', ...
    'count_channels','FS','Electrode_Type','sim_stim','-v7.3');

fprintf('\n============================================================\n');
fprintf('TRIAL SPIKE-COUNT EXTRACTION COMPLETE\n');
fprintf('Trials saved: %d\n',height(TrialSummary));
fprintf('Channels stored per trial: %d\n',nSelectedChannels);
fprintf('Channels used for totals: %d\n',nCountChannels);
fprintf('Count channels: %s\n',number_list(count_channels));
fprintf('Output file:\n%s\n',full_save_path);
fprintf('Original experiment/spike files modified: NO\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function validate_window(window_ms,variable_name)
if ~isnumeric(window_ms) || numel(window_ms) ~= 2 || any(~isfinite(window_ms)) || window_ms(2) <= window_ms(1)
    error('%s must be [start end], with end > start.',variable_name);
end
end

function files = remove_metadata_and_backups(files)
if isempty(files), return; end
keep = true(size(files));
for file_index = 1:numel(files)
    file_name = files(file_index).name;
    if startsWith(file_name,'._') || contains(file_name,'BACKUP','IgnoreCase',true)
        keep(file_index) = false;
    end
end
files = files(keep);
end

function label = format_order(stim_channels)
stim_channels = stim_channels(stim_channels > 0);
if isempty(stim_channels)
    label = '(none)';
elseif numel(stim_channels) == 1
    label = sprintf('Ch%d',stim_channels(1));
else
    labels = arrayfun(@(c) sprintf('Ch%d',c),stim_channels,'UniformOutput',false);
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