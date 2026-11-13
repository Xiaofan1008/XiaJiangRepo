%% SINGLE-STIMULATION RESPONDING-CHANNEL DETECTION
% Primary rule: mean post FR >= mean baseline FR + 3*baseline SD.
% Additional safeguards require a minimum spike count, active-trial
% fraction, and absolute post-stimulation firing rate.
%
% This script:
%   - requires sp_corr from *.sp_xia_SSD.mat;
%   - accepts manually identified artifact channels per stimulation set;
%   - does not remove bad trials at this stage;
%   - preserves the familiar Responding.set(si).amp(ai).ptd(1).channel(ich)
%     structure;
%   - saves *_RespondingChannels.mat and backs up an existing output.

clear;
close all;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ========================= USER SETTINGS ==============================

data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_Single1';

% 0 = rigid single shank; 1 = flexible single shank; 2 = four shank
Electrode_Type = 2;
FS = 30000;

baseline_win_ms = [-50 -10];
post_win_ms     = [2 20];
k_SD            = 3;

% Safeguards in addition to the baseline + 3 SD rule
min_total_post_spikes      = 2;
min_fraction_active_trials = 0.1;
min_mean_post_FR           = 5;

% Artifact-biased recording channels identified from waveform review.
% Channel numbers are the Depth_s/channel indices displayed in figures.
% Each cell corresponds to one stimulation set. Examples:
% BadChannels_PerSet = {
%     [12 18 25]   % Set 1
%     [12 19 26]   % Set 2
% };
BadChannels_PerSet = {};

Create_Backup_If_Output_Exists = true;

%% ======================== INITIAL CHECKS ==============================

if ~isfolder(data_folder)
    error('Dataset folder does not exist:\n%s',data_folder);
end

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(post_win_ms,'post_win_ms');

if ~isscalar(k_SD) || ~isfinite(k_SD) || k_SD < 0
    error('k_SD must be a finite nonnegative scalar.');
end
if ~isscalar(min_total_post_spikes) || min_total_post_spikes < 0
    error('min_total_post_spikes must be nonnegative.');
end
if ~isscalar(min_fraction_active_trials) || ...
        min_fraction_active_trials < 0 || min_fraction_active_trials > 1
    error('min_fraction_active_trials must be between 0 and 1.');
end
if ~isscalar(min_mean_post_FR) || min_mean_post_FR < 0
    error('min_mean_post_FR must be nonnegative.');
end

baseline_duration_s = diff(baseline_win_ms)/1000;
post_duration_s     = diff(post_win_ms)/1000;

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>
cd(data_folder);

fprintf('\n============================================================\n');
fprintf('SINGLE-STIMULATION RESPONDING-CHANNEL DETECTION\n');
fprintf('============================================================\n');
fprintf('Dataset: %s\n',data_folder);
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Post window: [%g,%g) ms\n',post_win_ms);
fprintf('Primary rule: mean baseline + %g SD\n',k_SD);
fprintf('Bad trials excluded: NO\n');

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

if ~ismember('sp_corr',who('-file',ssd_file))
    error('The SSD file does not contain sp_corr:\n%s',ssd_file);
end

SpikeLoad = load(ssd_file,'sp_corr');
sp = SpikeLoad.sp_corr;
nSpChannels = numel(sp);

fprintf('Spike file: %s\n',ssd_file);
fprintf('Spike channels: %d\n',nSpChannels);

%% ========================== LOAD TRIGGERS =============================

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
trig_ms = trig/FS*1000;
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

if sim_stim ~= 1
    error(['This script requires single stimulation, but ' ...
        'simultaneous_stim = %d.'],sim_stim);
end
if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.',nTrig,n_Trials);
elseif nTrig > n_Trials
    warning('%d triggers were loaded for %d trials; extras are ignored.', ...
        nTrig,n_Trials);
end

fprintf('Experiment file: %s\n',experiment_files(1).name);
fprintf('Trials/triggers: %d/%d\n',n_Trials,nTrig);

%% ===================== AMPLITUDES AND SETS ============================

trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:sim_stim:end);
if numel(trialAmps) ~= n_Trials
    error('Decoded %d amplitudes for %d trials.',numel(trialAmps),n_Trials);
end
trialAmps(trialAmps == -1) = 0;
[Amps,~,ampIdx] = unique(trialAmps(:));
nAMP = numel(Amps);

% Keep one PTD=0 layer for compatibility with existing scripts.
PTDs = 0;
ptdIdx = ones(n_Trials,1);

stimNames = StimParams(2:end,1);
[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));
if any(~isMapped)
    error('%d stimulation entries could not be mapped through E_MAP.', ...
        sum(~isMapped));
end

stimSeq = zeros(n_Trials,sim_stim);
for trial_id = 1:n_Trials
    rows_this_trial = (trial_id-1)*sim_stim+(1:sim_stim);
    stimulation_indices = idx_all(rows_this_trial);
    stimulation_indices = stimulation_indices(stimulation_indices > 0);
    stimSeq(trial_id,1:numel(stimulation_indices)) = stimulation_indices;
end
[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

fprintf('Amplitudes: %s uA\n',num2str(Amps(:).'));
fprintf('Stimulation sets: %d\n',nSets);
for si = 1:nSets
    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);
    fprintf('  Set %d: %s\n',si,channel_list(stim_channels));
end

%% ======================== DEPTH CHANNEL MAP ===========================

d = double(Depth_s(Electrode_Type));
d = d(:);
if isempty(d)
    error('Depth_s returned an empty channel map.');
end
nDepthChannels = numel(d);

%% ================= MANUAL BAD-CHANNEL SETTINGS =======================

if isempty(BadChannels_PerSet)
    BadChannels_PerSet = cell(nSets,1);
elseif ~iscell(BadChannels_PerSet)
    error('BadChannels_PerSet must be a cell array.');
elseif numel(BadChannels_PerSet) < nSets
    BadChannels_PerSet(end+1:nSets) = {[]};
elseif numel(BadChannels_PerSet) > nSets
    warning('Extra BadChannels_PerSet entries will be ignored.');
    BadChannels_PerSet = BadChannels_PerSet(1:nSets);
end

for si = 1:nSets
    bad_channels = unique(double(BadChannels_PerSet{si}(:).'));
    valid = bad_channels >= 1 & bad_channels <= nDepthChannels & ...
        fix(bad_channels) == bad_channels;
    if any(~valid)
        warning('Invalid bad-channel indices for Set %d were ignored.',si);
    end
    BadChannels_PerSet{si} = bad_channels(valid);
    fprintf('Set %d artifact channels: %s\n',si, ...
        channel_list(BadChannels_PerSet{si}));
end

%% ====================== INITIALIZE OUTPUT =============================

Responding = struct();
Responding.metadata.analysis_name = ...
    'Single-stimulation responding-channel detection';
Responding.metadata.created = char(datetime('now', ...
    'Format','yyyy-MM-dd HH:mm:ss'));
Responding.metadata.data_folder = data_folder;
Responding.metadata.spike_file = ssd_file;
Responding.metadata.spike_variable = 'sp_corr';
Responding.metadata.experiment_file = experiment_files(1).name;
Responding.metadata.sampling_rate_hz = FS;
Responding.metadata.electrode_type = Electrode_Type;
Responding.metadata.baseline_window_ms = baseline_win_ms;
Responding.metadata.post_window_ms = post_win_ms;
Responding.metadata.k_SD = k_SD;
Responding.metadata.min_total_post_spikes = min_total_post_spikes;
Responding.metadata.min_fraction_active_trials = ...
    min_fraction_active_trials;
Responding.metadata.min_mean_post_FR = min_mean_post_FR;
Responding.metadata.bad_trials_excluded = false;
Responding.metadata.provisional_until_bad_trial_review = true;
Responding.metadata.response_rule = ...
    ['mean_post_FR >= mean_baseline_FR + k_SD*sd_baseline_FR, ' ...
    'plus minimum spike/fraction/FR safeguards'];

%% ========================= RESPONSE DETECTION =========================

fprintf('\nRunning response detection...\n');

for si = 1:nSets
    bad_channels_this_set = BadChannels_PerSet{si};
    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    Responding.set(si).set_index = si;
    Responding.set(si).stimChannels = stim_channels;
    Responding.set(si).badChannels = bad_channels_this_set;

    for ai = 1:nAMP
        amp_value = Amps(ai);
        pi = 1;
        condition_trials = find(combClass == si & ampIdx == ai & ...
            ptdIdx == pi);

        Responding.set(si).amp(ai).amp_value = amp_value;

        P = struct();
        P.PTD_us = 0;
        P.PTD_ms = 0;
        P.amp_value = amp_value;
        P.set_index = si;
        P.amp_index = ai;
        P.ptd_index = pi;
        P.trial_ids = condition_trials(:).';
        P.n_trials = numel(condition_trials);
        P.bad_trials_excluded = false;

        responsive_channels = [];

        for ich = 1:nDepthChannels
            R = blank_channel_result();
            R.channel_index = ich;
            R.spike_data_channel = d(ich);
            R.trial_ids = condition_trials(:).';
            R.n_trials = numel(condition_trials);
            R.is_excluded_bad_channel = ...
                ismember(ich,bad_channels_this_set);

            spike_channel = d(ich);
            if spike_channel < 1 || spike_channel > nSpChannels
                R.exclusion_reason = 'Depth map outside spike-data range';
                P.channel(ich) = R;
                continue;
            end
            if R.is_excluded_bad_channel
                R.exclusion_reason = ...
                    'Manually excluded artifact-biased channel';
                P.channel(ich) = R;
                continue;
            end
            if isempty(sp{spike_channel})
                R.exclusion_reason = 'No spikes in sp_corr';
                P.channel(ich) = R;
                continue;
            end
            if isempty(condition_trials)
                R.exclusion_reason = 'No trials for this condition';
                P.channel(ich) = R;
                continue;
            end

            R.has_spike_data = true;
            sp_times = double(sp{spike_channel}(:,1));
            nConditionTrials = numel(condition_trials);

            baseline_counts = zeros(1,nConditionTrials);
            post_counts = zeros(1,nConditionTrials);

            for trial_position = 1:nConditionTrials
                trial_id = condition_trials(trial_position);
                t0 = trig_ms(trial_id);
                baseline_counts(trial_position) = sum( ...
                    sp_times >= t0+baseline_win_ms(1) & ...
                    sp_times <  t0+baseline_win_ms(2));
                post_counts(trial_position) = sum( ...
                    sp_times >= t0+post_win_ms(1) & ...
                    sp_times <  t0+post_win_ms(2));
            end

            FR_baseline = baseline_counts/baseline_duration_s;
            FR_post = post_counts/post_duration_s;

            mean_baseline_FR = mean(FR_baseline);
            sd_baseline_FR = std(FR_baseline,0);
            mean_post_FR = mean(FR_post);
            response_threshold_FR = ...
                mean_baseline_FR+k_SD*sd_baseline_FR;

            total_baseline_spikes = sum(baseline_counts);
            total_post_spikes = sum(post_counts);
            fraction_active_trials = mean(post_counts > 0);

            pass_3SD_rule = mean_post_FR >= response_threshold_FR;
            pass_minimum_spikes = ...
                total_post_spikes >= min_total_post_spikes;
            pass_active_fraction = ...
                fraction_active_trials >= min_fraction_active_trials;
            pass_minimum_FR = mean_post_FR >= min_mean_post_FR;

            R.exclusion_reason = '';
            R.baseline_counts_all = baseline_counts;
            R.post_counts_all = post_counts;
            R.FR_baseline_all = FR_baseline;
            R.FR_post_all = FR_post;
            R.mean_baseline_FR = mean_baseline_FR;
            R.sd_baseline_FR = sd_baseline_FR;
            R.response_threshold_FR = response_threshold_FR;
            R.mean_post_FR = mean_post_FR;
            R.total_baseline_spikes = total_baseline_spikes;
            R.total_post_spikes = total_post_spikes;
            R.fraction_active_trials = fraction_active_trials;
            R.pass_3SD_rule = pass_3SD_rule;
            R.pass_minimum_spikes = pass_minimum_spikes;
            R.pass_active_fraction = pass_active_fraction;
            R.pass_minimum_FR = pass_minimum_FR;
            R.is_responsive = pass_3SD_rule && pass_minimum_spikes && ...
                pass_active_fraction && pass_minimum_FR;

            P.channel(ich) = R;
            if R.is_responsive
                responsive_channels(end+1) = ich; %#ok<SAGROW>
            end
        end

        P.responsive_channels = responsive_channels;
        P.n_responsive_channels = numel(responsive_channels);
        Responding.set(si).amp(ai).ptd(pi) = P;

        fprintf(['Set %d | Amp %g uA | trials %d | ' ...
            'responding (%d): %s\n'],si,amp_value, ...
            numel(condition_trials),numel(responsive_channels), ...
            channel_list(responsive_channels));
    end
end

%% =========================== SAVE RESULT ==============================

outfile = sprintf('%s_RespondingChannels.mat',base_name);
full_output_path = fullfile(data_folder,outfile);

if isfile(full_output_path) && Create_Backup_If_Output_Exists
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    backup_file = sprintf('%s_RespondingChannels_BACKUP_%s.mat', ...
        base_name,timestamp);
    backup_path = fullfile(data_folder,backup_file);
    copyfile(full_output_path,backup_path);
    fprintf('\nExisting output backed up to:\n%s\n',backup_path);
end

% Familiar top-level variables retained for existing plotting/manual code.
Detection_Mode = 1;
save(full_output_path,'Responding','Detection_Mode','baseline_win_ms', ...
    'post_win_ms','k_SD','min_total_post_spikes', ...
    'min_fraction_active_trials','min_mean_post_FR','Amps','PTDs', ...
    'BadChannels_PerSet');

fprintf('\n============================================================\n');
fprintf('RESPONDING-CHANNEL DETECTION COMPLETE\n');
fprintf('Output: %s\n',full_output_path);
fprintf('Bad trials excluded: NO\n');
fprintf('Result remains provisional until bad-trial review.\n');
fprintf('============================================================\n');

%% ========================== LOCAL FUNCTIONS ===========================

function validate_window(value,name)
if ~isnumeric(value) || numel(value) ~= 2 || any(~isfinite(value)) || ...
        value(2) <= value(1)
    error('%s must be [start end], with end > start.',name);
end
end

function R = blank_channel_result()
R = struct('channel_index',NaN,'spike_data_channel',NaN, ...
    'trial_ids',[],'n_trials',0,'has_spike_data',false, ...
    'is_excluded_bad_channel',false,'exclusion_reason','', ...
    'baseline_counts_all',[],'post_counts_all',[], ...
    'FR_baseline_all',[],'FR_post_all',[], ...
    'mean_baseline_FR',NaN,'sd_baseline_FR',NaN, ...
    'response_threshold_FR',NaN,'mean_post_FR',NaN, ...
    'total_baseline_spikes',NaN,'total_post_spikes',NaN, ...
    'fraction_active_trials',NaN,'pass_3SD_rule',false, ...
    'pass_minimum_spikes',false,'pass_active_fraction',false, ...
    'pass_minimum_FR',false,'is_responsive',false);
end

function text_value = channel_list(values)
if isempty(values)
    text_value = '(none)';
else
    labels = arrayfun(@(x) sprintf('Ch%d',x),values, ...
        'UniformOutput',false);
    text_value = strjoin(labels,' ');
end
end
