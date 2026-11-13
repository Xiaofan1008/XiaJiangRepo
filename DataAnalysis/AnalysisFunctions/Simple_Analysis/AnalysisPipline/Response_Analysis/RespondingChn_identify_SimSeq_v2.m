%% ========================================================================
% MULTI-ISI RESPONDING-CHANNEL DETECTION
%
% PURPOSE
%   Detect responding recording channels for paired simultaneous and
%   sequential stimulation using:
%
%       Mean post-stimulation firing rate
%           >=
%       Mean baseline firing rate + 3 * baseline firing-rate SD
%
% CONDITIONS
%   The analysis is performed independently for every:
%
%       stimulation order/set × amplitude × recorded PTD × channel
%
%   Therefore, A->B and B->A sequential stimulation remain separate and
%   can have different responding-channel populations.
%
% WINDOWS
%   Baseline: [-50,-5) ms relative to the first trigger
%   Response: [2,40) ms relative to the first trigger
%
%   The [2,40) ms response window is used because the multi-ISI experiment
%   contains PTDs up to approximately 20 ms. This window can include the
%   responses following both stimulation pulses.
%
% SPIKE DATA
%   This script strictly requires:
%
%       *.sp_xia_SSD.mat
%       variable: sp_corr
%
% CHANNELS
%   All channel indices returned by Depth_s(Electrode_Type) are analysed
%   automatically. No manual channel range is required.
%
% BAD TRIALS AND BAD CHANNELS
%   This initial version does not remove bad trials or bad channels.
%   Responding channels can be reviewed and manually corrected afterward.
%
% OUTPUT
%   <base_name>_MultiISI_RespondingChannels.mat
%
%   The main saved variable is named:
%
%       Responding
%
%   This preserves compatibility with the existing Responding structure.
% ========================================================================

clear;
close all;

addpath(genpath( ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';

% Electrode type:
%   0 = rigid single-shank probe
%   1 = flexible single-shank probe
%   2 = four-shank flexible probe
Electrode_Type = 2;

% Sampling rate in Hz
FS = 30000;

% Baseline window relative to the first stimulation trigger
% The upper boundary is excluded.
baseline_win_ms = [-50 -5];

% Whole paired-response window relative to the first stimulation trigger
% The upper boundary is excluded.
post_win_ms = [2 40];

% Responding-channel threshold:
% mean post FR >= mean baseline FR + k_SD * baseline FR SD
k_SD = 3;

% Additional safeguards
%
% These prevent nearly silent channels from being classified as responsive
% only because their baseline mean and SD are both zero.
Use_Additional_Safeguards = true;

% Minimum total post-stimulation spikes across all trials
min_total_post_spikes = 2;

% Minimum fraction of trials containing at least one post-stimulus spike
min_frac_trials_with_spikes = 0.13;

% Minimum mean post-stimulation firing rate
min_abs_post_FR = 5;

% Create a timestamped backup if the output already exists
Create_Backup_If_Output_Exists = true;

%% =========================== INITIAL CHECKS ===========================

if ~isfolder(data_folder)
    error('Dataset folder does not exist:\n%s',data_folder);
end

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(post_win_ms,'post_win_ms');

if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('FS must be a positive sampling rate in Hz.');
end

if ~isscalar(k_SD) || ~isfinite(k_SD) || k_SD < 0
    error('k_SD must be a finite nonnegative number.');
end

% Return to the original MATLAB folder when the script finishes
starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

% Automatically move to the selected dataset folder
cd(data_folder);

baseline_duration_s = diff(baseline_win_ms)/1000;
post_duration_s = diff(post_win_ms)/1000;

fprintf('\n============================================================\n');
fprintf('MULTI-ISI RESPONDING-CHANNEL DETECTION\n');
fprintf('============================================================\n');
fprintf('Dataset: %s\n',data_folder);
fprintf('Electrode type: %d\n',Electrode_Type);
fprintf('Sampling rate: %g Hz\n',FS);
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Response window: [%g,%g) ms\n',post_win_ms);
fprintf('Detection rule: mean baseline FR + %g SD\n',k_SD);

if Use_Additional_Safeguards
    fprintf('Additional response safeguards: ON\n');
    fprintf('  Minimum total post spikes: %g\n', ...
        min_total_post_spikes);
    fprintf('  Minimum active-trial fraction: %.3f\n', ...
        min_frac_trials_with_spikes);
    fprintf('  Minimum mean post FR: %g spikes/s\n', ...
        min_abs_post_FR);
else
    fprintf('Additional response safeguards: OFF\n');
end

fprintf('Bad-channel file used: NO\n');
fprintf('Bad-trial file used: NO\n');

%% ============================ LOAD sp_corr ============================

ssd_files = dir('*.sp_xia_SSD.mat');

if isempty(ssd_files)
    error('No *.sp_xia_SSD.mat file was found.');
end

if numel(ssd_files) > 1
    warning('Multiple SSD files found. Using: %s',ssd_files(1).name);
end

ssd_file = ssd_files(1).name;
base_name = erase(ssd_file,'.sp_xia_SSD.mat');

if ~ismember('sp_corr',who('-file',ssd_file))
    error(['The SSD file does not contain sp_corr.\n' ...
        'Required file: %s'],ssd_file);
end

SpikeLoad = load(ssd_file,'sp_corr');
sp = SpikeLoad.sp_corr;
nSpChannels = numel(sp);

if ~iscell(sp)
    error('sp_corr must be a cell array containing one cell per channel.');
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
trig_ms = trig/FS*1000;
nTrig = numel(trig);

fprintf('Trigger file: %s\n',trigger_files(1).name);
fprintf('Triggers loaded: %d\n',nTrig);

%% ==================== LOAD EXPERIMENT PARAMETERS =====================

experiment_files = dir('*_exp_datafile_*.mat');

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found.');
end

if numel(experiment_files) > 1
    warning('Multiple experiment files found. Using: %s', ...
        experiment_files(1).name);
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
sim_stim   = double(ExpLoad.simultaneous_stim);
E_MAP      = ExpLoad.E_MAP;
n_Trials   = double(ExpLoad.n_Trials);

% This script is designed for paired stimulation
if sim_stim ~= 2
    error(['This script requires two stimulation events per trial, but ' ...
        'simultaneous_stim = %d.'],sim_stim);
end

if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.', ...
        nTrig,n_Trials);

elseif nTrig > n_Trials
    warning('%d triggers were loaded for %d trials; extras are ignored.', ...
        nTrig,n_Trials);

    trig = trig(1:n_Trials);
    trig_ms = trig_ms(1:n_Trials);
end

expected_parameter_rows = 1+n_Trials*sim_stim;

if size(StimParams,1) < expected_parameter_rows
    error(['StimParams contains insufficient rows. Expected at least %d ' ...
        'rows but found %d.'], ...
        expected_parameter_rows,size(StimParams,1));
end

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Trials/triggers used: %d/%d\n',n_Trials,n_Trials);
fprintf('Stimulation events per trial: %d\n',sim_stim);

%% =========================== DECODE AMPLITUDES ========================

% The first stimulation event in every trial is used as the trial-level
% amplitude. The code also checks whether the two pulses have equal
% amplitudes.
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

amplitude_mismatch = ...
    abs(trialAmps-secondPulseAmps) > 1e-6;

if any(amplitude_mismatch)
    warning(['The two pulses have different amplitudes in %d trials. ' ...
        'Conditions will be grouped using the first-pulse amplitude.'], ...
        sum(amplitude_mismatch));
end

[Amps,~,ampIdx] = unique(trialAmps(:));
nAMP = numel(Amps);

%% ============================== DECODE PTDs ===========================

% PTD is taken from the second stimulation event in every trial.
% StimParams column 6 stores PTD in microseconds.
PTD_all_us = cell2mat(StimParams(2:end,6));
PTD_all_us = double(PTD_all_us(:));

trialPTD_us = PTD_all_us(second_event_rows);
trialPTD_us = trialPTD_us(1:n_Trials);

[PTDs,~,ptdIdx] = unique(trialPTD_us(:));
PTDs_ms = PTDs/1000;
nPTD = numel(PTDs);

%% ==================== DECODE STIMULATION ORDERS ======================

stimNames = StimParams(2:end,1);
stimNames = stimNames(1:n_Trials*sim_stim);

[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));

if any(~isMapped)
    error('%d stimulation entries could not be mapped through E_MAP.', ...
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

% Order is deliberately preserved:
% [A B] and [B A] become different sets.
[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

fprintf('\nAmplitudes: %s uA\n',num2str(Amps(:).'));
fprintf('Observed PTDs: %s ms\n',num2str(PTDs_ms(:).'));
fprintf('Ordered stimulation sets: %d\n',nSets);

for si = 1:nSets
    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    fprintf('  Set %d: %s\n', ...
        si,format_stimulation_order(stim_channels));
end

%% ======================= AUTOMATIC CHANNEL MAPPING ====================

% Depth_s reads the recording configuration from the selected dataset and
% maps the displayed channel indices to the corresponding sp_corr cells.
d = double(Depth_s(Electrode_Type));
d = d(:);

if isempty(d)
    error('Depth_s returned an empty channel map.');
end

channel_indices = 1:numel(d);
nChannelIndices = numel(channel_indices);

fprintf('\nAutomatically detected channel indices: 1:%d\n', ...
    nChannelIndices);

% Validate the mapped sp_corr indices
invalid_mapped_channels = ...
    ~isfinite(d) | ...
    d < 1 | ...
    d > nSpChannels | ...
    fix(d) ~= d;

if any(invalid_mapped_channels)
    warning(['%d channel indices map outside the available sp_corr ' ...
        'cells. Those channels will be marked nonresponsive.'], ...
        sum(invalid_mapped_channels));
end

%% ======================== INITIALIZE OUTPUT ===========================

Responding = struct();

Responding.metadata.analysis_name = ...
    'Multi-ISI responding-channel detection';

Responding.metadata.created = char(datetime('now', ...
    'Format','yyyy-MM-dd HH:mm:ss'));

Responding.metadata.data_folder = data_folder;
Responding.metadata.spike_file = ssd_file;
Responding.metadata.spike_variable = 'sp_corr';
Responding.metadata.experiment_file = experiment_file;
Responding.metadata.sampling_rate_hz = FS;
Responding.metadata.electrode_type = Electrode_Type;

Responding.metadata.baseline_window_ms = baseline_win_ms;
Responding.metadata.post_window_ms = post_win_ms;
Responding.metadata.k_SD = k_SD;

Responding.metadata.detection_rule = ...
    'mean post FR >= mean baseline FR + k_SD * baseline FR SD';

Responding.metadata.use_additional_safeguards = ...
    Use_Additional_Safeguards;

Responding.metadata.min_total_post_spikes = ...
    min_total_post_spikes;

Responding.metadata.min_frac_trials_with_spikes = ...
    min_frac_trials_with_spikes;

Responding.metadata.min_abs_post_FR = ...
    min_abs_post_FR;

Responding.metadata.bad_channels_removed = false;
Responding.metadata.bad_trials_removed = false;
Responding.metadata.order_sensitive_sets = true;

Responding.metadata.amplitudes_uA = Amps;
Responding.metadata.PTDs_us = PTDs;
Responding.metadata.PTDs_ms = PTDs_ms;
Responding.metadata.unique_stimulation_orders = uniqueComb;
Responding.metadata.channel_indices = channel_indices;
Responding.metadata.depth_map = d;

%% =====================================================================
% MAIN LOOP: SET × AMPLITUDE × PTD × RECORDING CHANNEL
% ======================================================================

fprintf('\nRunning responding-channel detection...\n');

for si = 1:nSets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    Responding.set(si).set_index = si;
    Responding.set(si).stimChannels = stim_channels;
    Responding.set(si).stimOrderLabel = ...
        format_stimulation_order(stim_channels);

    for ai = 1:nAMP

        Responding.set(si).amp(ai).amp_index = ai;
        Responding.set(si).amp(ai).amp_value = Amps(ai);

        for pi = 1:nPTD

            current_ptd_us = PTDs(pi);
            current_ptd_ms = PTDs_ms(pi);

            % Find the original trials belonging to this exact condition
            trials_this = find( ...
                combClass == si & ...
                ampIdx == ai & ...
                ptdIdx == pi);

            condition = struct();

            condition.PTD_us = current_ptd_us;
            condition.PTD_ms = current_ptd_ms;
            condition.amp_value = Amps(ai);
            condition.set_index = si;
            condition.amp_index = ai;
            condition.ptd_index = pi;
            condition.trial_ids = trials_this(:).';
            condition.n_trials = numel(trials_this);

            % Preallocate one consistent result structure for every channel
            channel_template = empty_channel_result();

            condition.channel = repmat( ...
                channel_template,1,nChannelIndices);

            for ich = 1:nChannelIndices

                R = channel_template;

                R.channel_index = ich;
                R.spike_data_channel = d(ich);
                R.trial_ids = trials_this(:).';
                R.n_trials = numel(trials_this);

                % Empty experimental condition
                if isempty(trials_this)
                    R.status = 'no_trials';
                    condition.channel(ich) = R;
                    continue;
                end

                % Invalid channel mapping
                if invalid_mapped_channels(ich)
                    R.status = 'invalid_channel_mapping';
                    condition.channel(ich) = R;
                    continue;
                end

                spike_channel = d(ich);

                % Channel contains no accepted spikes
                if isempty(sp{spike_channel})
                    R.status = 'no_spike_data';
                    condition.channel(ich) = R;
                    continue;
                end

                spike_times = double(sp{spike_channel}(:,1));

                nTr = numel(trials_this);

                baseline_counts = zeros(1,nTr);
                post_counts = zeros(1,nTr);

                FR_baseline = zeros(1,nTr);
                FR_post = zeros(1,nTr);

                %% ---------------- COUNT EACH TRIAL --------------------

                for trial_position = 1:nTr

                    trial_id = trials_this(trial_position);
                    trigger_time_ms = trig_ms(trial_id);

                    baseline_start = ...
                        trigger_time_ms+baseline_win_ms(1);

                    baseline_end = ...
                        trigger_time_ms+baseline_win_ms(2);

                    post_start = ...
                        trigger_time_ms+post_win_ms(1);

                    post_end = ...
                        trigger_time_ms+post_win_ms(2);

                    baseline_mask = ...
                        spike_times >= baseline_start & ...
                        spike_times < baseline_end;

                    post_mask = ...
                        spike_times >= post_start & ...
                        spike_times < post_end;

                    baseline_counts(trial_position) = ...
                        sum(baseline_mask);

                    post_counts(trial_position) = ...
                        sum(post_mask);

                    FR_baseline(trial_position) = ...
                        baseline_counts(trial_position) / ...
                        baseline_duration_s;

                    FR_post(trial_position) = ...
                        post_counts(trial_position) / ...
                        post_duration_s;
                end

                %% ---------------- BASELINE + 3 SD RULE ----------------

                mean_baseline_FR = mean(FR_baseline);
                sd_baseline_FR = std(FR_baseline,0);

                response_threshold_FR = ...
                    mean_baseline_FR+k_SD*sd_baseline_FR;

                mean_post_FR = mean(FR_post);

                total_baseline_spikes = sum(baseline_counts);
                total_post_spikes = sum(post_counts);

                frac_post_trials = mean(post_counts > 0);

                pass_FR_rule = ...
                    mean_post_FR >= response_threshold_FR;

                pass_total_post_spikes = ...
                    total_post_spikes >= min_total_post_spikes;

                pass_frac_trials = ...
                    frac_post_trials >= ...
                    min_frac_trials_with_spikes;

                pass_abs_post_FR = ...
                    mean_post_FR >= min_abs_post_FR;

                if Use_Additional_Safeguards
                    isResp = ...
                        pass_FR_rule && ...
                        pass_total_post_spikes && ...
                        pass_frac_trials && ...
                        pass_abs_post_FR;
                else
                    isResp = pass_FR_rule;
                end

                %% ---------------- SAVE CHANNEL RESULT -----------------

                R.status = 'analysed';
                R.is_responsive = logical(isResp);

                R.mean_baseline_FR = mean_baseline_FR;
                R.sd_baseline_FR = sd_baseline_FR;
                R.response_threshold_FR = response_threshold_FR;
                R.mean_post_FR = mean_post_FR;

                R.FR_baseline_all = FR_baseline;
                R.FR_post_all = FR_post;

                R.baseline_spike_counts_all = baseline_counts;
                R.post_spike_counts_all = post_counts;

                R.total_baseline_spikes = ...
                    total_baseline_spikes;

                R.total_post_spikes = ...
                    total_post_spikes;

                R.frac_post_trials = frac_post_trials;

                R.pass_FR_rule = logical(pass_FR_rule);
                R.pass_total_post_spikes = ...
                    logical(pass_total_post_spikes);

                R.pass_frac_trials = ...
                    logical(pass_frac_trials);

                R.pass_abs_post_FR = ...
                    logical(pass_abs_post_FR);

                condition.channel(ich) = R;
            end

            % Store this set × amplitude × PTD result
            Responding.set(si).amp(ai).ptd(pi) = condition;

            % Print a compact condition summary
            responsive_mask = false(1,nChannelIndices);

            for ich = 1:nChannelIndices
                responsive_mask(ich) = ...
                    condition.channel(ich).is_responsive;
            end

            fprintf(['  Set %d | %s | Amp %g uA | PTD %g ms | ' ...
                'Trials %d | Responding %d/%d\n'], ...
                si, ...
                format_stimulation_order(stim_channels), ...
                Amps(ai), ...
                current_ptd_ms, ...
                numel(trials_this), ...
                sum(responsive_mask), ...
                nChannelIndices);
        end
    end
end

%% =========================== SAVE RESULT ==============================

output_name = sprintf( ...
    '%s_MultiISI_RespondingChannels.mat',base_name);

full_output_path = fullfile(data_folder,output_name);

% Protect an existing file before overwriting it
if isfile(full_output_path) && Create_Backup_If_Output_Exists

    timestamp = datestr(now,'yyyymmdd_HHMMSS');

    backup_name = sprintf( ...
        '%s_MultiISI_RespondingChannels_BACKUP_%s.mat', ...
        base_name,timestamp);

    backup_path = fullfile(data_folder,backup_name);

    copyfile(full_output_path,backup_path);

    fprintf('\nExisting output backed up to:\n%s\n', ...
        backup_path);
end

% Retain this variable for compatibility with older analysis scripts
Detection_Mode = 1;

save(full_output_path, ...
    'Responding', ...
    'Detection_Mode', ...
    'baseline_win_ms', ...
    'post_win_ms', ...
    'k_SD', ...
    'Use_Additional_Safeguards', ...
    'min_total_post_spikes', ...
    'min_frac_trials_with_spikes', ...
    'min_abs_post_FR', ...
    'Amps', ...
    'PTDs', ...
    'PTDs_ms', ...
    'uniqueComb', ...
    'FS', ...
    'Electrode_Type', ...
    '-v7.3');

fprintf('\n============================================================\n');
fprintf('MULTI-ISI RESPONDING-CHANNEL DETECTION COMPLETE\n');
fprintf('Output file:\n%s\n',full_output_path);
fprintf('Main saved variable: Responding\n');
fprintf('Bad trials removed: NO\n');
fprintf('Bad channels removed: NO\n');
fprintf('Original experiment files modified: NO\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function validate_window(window_ms,variable_name)
% Confirm that a window has the form [start end], where end > start.

if ~isnumeric(window_ms) || ...
        numel(window_ms) ~= 2 || ...
        any(~isfinite(window_ms)) || ...
        window_ms(2) <= window_ms(1)

    error('%s must be [start end], with end > start.', ...
        variable_name);
end
end

function label = format_stimulation_order(stim_channels)
% Produce a readable stimulation-order label.

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

function R = empty_channel_result()
% Create a consistent empty result structure for one recording channel.

R = struct();

R.channel_index = NaN;
R.spike_data_channel = NaN;

R.status = 'not_analysed';
R.is_responsive = false;

R.trial_ids = [];
R.n_trials = 0;

R.mean_baseline_FR = NaN;
R.sd_baseline_FR = NaN;
R.response_threshold_FR = NaN;
R.mean_post_FR = NaN;

R.FR_baseline_all = [];
R.FR_post_all = [];

R.baseline_spike_counts_all = [];
R.post_spike_counts_all = [];

R.total_baseline_spikes = NaN;
R.total_post_spikes = NaN;
R.frac_post_trials = NaN;

R.pass_FR_rule = false;
R.pass_total_post_spikes = false;
R.pass_frac_trials = false;
R.pass_abs_post_FR = false;
end