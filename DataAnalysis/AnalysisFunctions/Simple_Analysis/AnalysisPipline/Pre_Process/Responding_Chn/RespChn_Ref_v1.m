%% ========================================================================
% RESPONDING-CHANNEL REFERENCE CRITERIA (PRINTED TO COMMAND WINDOW ONLY)
%
% PURPOSE
%   Automated reference for the responding-channel list: for every
%   stimulation condition (set x amplitude x PTD, or set x amplitude for
%   single-pulse files), flag channels whose mean post-stimulation firing
%   rate clears baseline by k_SD, same rule as
%   RespondingChn_identify_SimSeq_v2.m:
%
%       mean post FR  >=  mean baseline FR + k_SD * baseline FR SD
%
%   This is a REFERENCE only - meant to be read alongside
%   RespondingChn_Raster_AllChn_v2.m's figures, not applied automatically.
%   Nothing is saved to disk. The final, manually-confirmed list is
%   produced by a separate tool (not built yet).
%
% INPUT CONVENTION
%   Same as RespondingChn_Raster_AllChn_v2.m: sp_clipped / sp_waveforms
%   (physical-channel indexed), ChnMap(Electrode_Type, nCh) for the
%   map-position -> physical-channel lookup, StimParams/E_MAP/n_Trials
%   from the dataset's own *_exp_datafile_*.mat for condition decoding.
%
% WINDOWS
%   Baseline: [-50,-5) ms relative to the first stimulation pulse (t=0)
%   Response: [2,40) ms relative to the first stimulation pulse
%   These are independent of the fixed 0-20 ms window used later for the
%   spike-count linearity metric - this tool is about detecting ANY
%   evoked response for QC, not measuring its exact count.
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

% Empty means all available sets/amplitudes/PTDs
Plot_Sets = [];
Plot_Amps = [];
Plot_PTDs = [];

Condition_Tolerance = 1e-4;

%% ----------------------- RESPONSE-DETECTION SETTINGS -------------------

% Baseline/response windows relative to the first stimulation trigger.
% Upper boundary excluded, same convention as RespondingChn_identify_SimSeq_v2.m.
baseline_win_ms = [-50 -5];
post_win_ms     = [2 20];

% mean post FR >= mean baseline FR + k_SD * baseline FR SD
k_SD = 3;

% Additional safeguards (prevent near-silent channels passing only
% because baseline mean and SD are both ~0)
Use_Additional_Safeguards = true;
min_total_post_spikes       = 2;
min_frac_trials_with_spikes = 0.13;
min_abs_post_FR              = 5;

%% =========================== INITIAL CHECKS ===========================

if ~isfile(spike_file)
    error('Spike file does not exist:\n%s',spike_file);
end
if ~isfolder(dataset_folder)
    error('Dataset folder does not exist:\n%s',dataset_folder);
end

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(post_win_ms,'post_win_ms');

if ~isscalar(k_SD) || ~isfinite(k_SD) || k_SD < 0
    error('k_SD must be a finite nonnegative number.');
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

baseline_duration_s = diff(baseline_win_ms)/1000;
post_duration_s = diff(post_win_ms)/1000;

fprintf('\n============================================================\n');
fprintf('RESPONDING-CHANNEL REFERENCE CRITERIA\n');
fprintf('============================================================\n');
fprintf('Spike file: %s\n',spike_file);
fprintf('Dataset folder: %s\n',dataset_folder);
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Response window: [%g,%g) ms\n',post_win_ms);
fprintf('Detection rule: mean baseline FR + %g SD\n',k_SD);

if Use_Additional_Safeguards
    fprintf('Additional response safeguards: ON\n');
    fprintf('  Minimum total post spikes: %g\n',min_total_post_spikes);
    fprintf('  Minimum active-trial fraction: %.3f\n',min_frac_trials_with_spikes);
    fprintf('  Minimum mean post FR: %g spikes/s\n',min_abs_post_FR);
else
    fprintf('Additional response safeguards: OFF\n');
end

fprintf('This is a REFERENCE ONLY. Nothing is saved to disk.\n');

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

invalid_mapped_channels = ~isfinite(map_nums_plus) | map_nums_plus < 1 | ...
    map_nums_plus > nCh | fix(map_nums_plus) ~= map_nums_plus;

if any(invalid_mapped_channels)
    warning(['%d map positions map outside the available spike-data channels. ' ...
        'Those channels will be marked nonresponsive.'],sum(invalid_mapped_channels));
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
trig_ms = trig/FS*1000;

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
    trig_ms = trig_ms(1:n_Trials);
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
    PTDs = 0;
    PTDs_ms = 0;
    nPTD = 1;
    ptdIdx = ones(n_Trials,1);
end

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

%% ======================= SELECT CONDITIONS TO ANALYSE ===================

if isempty(Plot_Sets)
    selected_sets = 1:nSets;
else
    selected_sets = unique(double(Plot_Sets(:).'),'stable');
    valid_sets = selected_sets >= 1 & selected_sets <= nSets & fix(selected_sets) == selected_sets;
    if any(~valid_sets)
        warning('Invalid Plot_Sets were ignored.');
        selected_sets = selected_sets(valid_sets);
    end
end
if isempty(selected_sets)
    error('No valid stimulation sets were selected.');
end

if isempty(Plot_Amps)
    selected_amps = Amps(:).';
else
    selected_amps = double(Plot_Amps(:).');
end

if isempty(Plot_PTDs)
    selected_ptds = PTDs_ms(:).';
else
    selected_ptds = double(Plot_PTDs(:).');
end

%% =========================== INITIALIZE OUTPUT ===========================

Responding = struct();
Responding.metadata.analysis_name = 'Responding-channel reference criteria (not saved)';
Responding.metadata.spike_file = spike_file;
Responding.metadata.spike_source_var = spike_source_var;
Responding.metadata.dataset_folder = dataset_folder;
Responding.metadata.experiment_file = experiment_file;
Responding.metadata.Electrode_Type = Electrode_Type;
Responding.metadata.sampling_rate_hz = FS;
Responding.metadata.baseline_window_ms = baseline_win_ms;
Responding.metadata.post_window_ms = post_win_ms;
Responding.metadata.k_SD = k_SD;
Responding.metadata.use_additional_safeguards = Use_Additional_Safeguards;
Responding.metadata.amplitudes_uA = Amps;
Responding.metadata.PTDs_ms = PTDs_ms;
Responding.metadata.unique_stimulation_orders = uniqueComb;
Responding.metadata.map_positions = 1:nMapPos;

%% =====================================================================
% MAIN LOOP: SET x AMPLITUDE x PTD x MAP-POSITION CHANNEL
% ======================================================================

fprintf('\nRunning responding-channel reference detection...\n');
fprintf('(Reference only - compare against the raster/PSTH figures before deciding.)\n\n');

for si = selected_sets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    Responding.set(si).set_index = si;
    Responding.set(si).stimChannels = stim_channels;
    Responding.set(si).stimOrderLabel = format_order(stim_channels);

    % Physical distance between the stim channels of this set (NaN, with
    % a note, if undefined - e.g. channels on different probes)
    dist_label = '';
    dist_um = NaN;
    dist_note = '';
    if numel(stim_channels) == 2
        [dist_um,dist_note] = ChnPairDistance(stim_channels(1),stim_channels(2),Electrode_Type);
        if isnan(dist_um)
            dist_label = sprintf(' | Dist n/a (%s)',dist_note);
        else
            dist_label = sprintf(' | Dist %g um',dist_um);
        end
    end
    Responding.set(si).distance_um = dist_um;
    Responding.set(si).distance_note = dist_note;

    for amp_value = selected_amps

        ai = find(abs(Amps-amp_value) < Condition_Tolerance,1);
        if isempty(ai), continue; end

        Responding.set(si).amp(ai).amp_index = ai;
        Responding.set(si).amp(ai).amp_value = Amps(ai);

        for ptd_value = selected_ptds

            pi = find(abs(PTDs_ms-ptd_value) < Condition_Tolerance,1);
            if isempty(pi), continue; end

            current_ptd_ms = PTDs_ms(pi);

            trials_this = find(combClass == si & ampIdx == ai & ptdIdx == pi);

            condition = struct();
            condition.PTD_ms = current_ptd_ms;
            condition.amp_value = Amps(ai);
            condition.set_index = si;
            condition.amp_index = ai;
            condition.ptd_index = pi;
            condition.trial_ids = trials_this(:).';
            condition.n_trials = numel(trials_this);

            channel_template = empty_channel_result();
            condition.channel = repmat(channel_template,1,nMapPos);

            if isempty(trials_this)
                for ich = 1:nMapPos
                    R = channel_template;
                    R.channel_index = ich;
                    R.status = 'no_trials';
                    condition.channel(ich) = R;
                end
                Responding.set(si).amp(ai).ptd(pi) = condition;
                continue;
            end

            responsive_mask = false(1,nMapPos);

            for ich = 1:nMapPos

                R = channel_template;
                R.channel_index = ich;
                R.trial_ids = trials_this(:).';
                R.n_trials = numel(trials_this);

                if invalid_mapped_channels(ich)
                    R.status = 'invalid_channel_mapping';
                    condition.channel(ich) = R;
                    continue;
                end

                phys_ch = map_nums_plus(ich);
                R.spike_data_channel = phys_ch;

                if isempty(sp{phys_ch})
                    R.status = 'no_spike_data';
                    condition.channel(ich) = R;
                    continue;
                end

                spike_times = sp{phys_ch};
                nTr = numel(trials_this);

                baseline_counts = zeros(1,nTr);
                post_counts = zeros(1,nTr);
                FR_baseline = zeros(1,nTr);
                FR_post = zeros(1,nTr);

                for trial_position = 1:nTr

                    trial_id = trials_this(trial_position);
                    trigger_time_ms = trig_ms(trial_id);

                    baseline_start = trigger_time_ms+baseline_win_ms(1);
                    baseline_end   = trigger_time_ms+baseline_win_ms(2);
                    post_start     = trigger_time_ms+post_win_ms(1);
                    post_end       = trigger_time_ms+post_win_ms(2);

                    baseline_mask = spike_times >= baseline_start & spike_times < baseline_end;
                    post_mask     = spike_times >= post_start & spike_times < post_end;

                    baseline_counts(trial_position) = sum(baseline_mask);
                    post_counts(trial_position) = sum(post_mask);

                    FR_baseline(trial_position) = baseline_counts(trial_position)/baseline_duration_s;
                    FR_post(trial_position) = post_counts(trial_position)/post_duration_s;
                end

                mean_baseline_FR = mean(FR_baseline);
                sd_baseline_FR = std(FR_baseline,0);
                response_threshold_FR = mean_baseline_FR+k_SD*sd_baseline_FR;
                mean_post_FR = mean(FR_post);

                total_post_spikes = sum(post_counts);
                frac_post_trials = mean(post_counts > 0);

                pass_FR_rule = mean_post_FR >= response_threshold_FR;
                pass_total_post_spikes = total_post_spikes >= min_total_post_spikes;
                pass_frac_trials = frac_post_trials >= min_frac_trials_with_spikes;
                pass_abs_post_FR = mean_post_FR >= min_abs_post_FR;

                if Use_Additional_Safeguards
                    isResp = pass_FR_rule && pass_total_post_spikes && pass_frac_trials && pass_abs_post_FR;
                else
                    isResp = pass_FR_rule;
                end

                R.status = 'analysed';
                R.is_responsive = logical(isResp);
                R.mean_baseline_FR = mean_baseline_FR;
                R.sd_baseline_FR = sd_baseline_FR;
                R.response_threshold_FR = response_threshold_FR;
                R.mean_post_FR = mean_post_FR;
                R.total_post_spikes = total_post_spikes;
                R.frac_post_trials = frac_post_trials;
                R.pass_FR_rule = logical(pass_FR_rule);
                R.pass_total_post_spikes = logical(pass_total_post_spikes);
                R.pass_frac_trials = logical(pass_frac_trials);
                R.pass_abs_post_FR = logical(pass_abs_post_FR);

                condition.channel(ich) = R;
                responsive_mask(ich) = isResp;
            end

            Responding.set(si).amp(ai).ptd(pi) = condition;

            %% ---------------- PRINT CONDITION SUMMARY ------------------

            if sim_stim == 1
                stim_mode = 'Single';
                cond_label = sprintf('Set %d | %s | Amp %g uA | nTrials %d | %s%s', ...
                    si,format_order(stim_channels),Amps(ai),numel(trials_this),stim_mode,dist_label);
            elseif abs(current_ptd_ms) < Condition_Tolerance
                stim_mode = 'Simultaneous';
                cond_label = sprintf('Set %d | %s | Amp %g uA | PTD %g ms | nTrials %d | %s%s', ...
                    si,format_simultaneous(stim_channels),Amps(ai),current_ptd_ms,numel(trials_this),stim_mode,dist_label);
            else
                stim_mode = 'Sequential';
                cond_label = sprintf('Set %d | %s | Amp %g uA | PTD %g ms | nTrials %d | %s%s', ...
                    si,format_order(stim_channels),Amps(ai),current_ptd_ms,numel(trials_this),stim_mode,dist_label);
            end

            responding_channels = find(responsive_mask);

            fprintf('%s\n',cond_label);
            fprintf('  Responding (reference): %d/%d\n',numel(responding_channels),nMapPos);
            fprintf('  Channels: %s\n\n',channel_list_text(responding_channels));
        end
    end
end

fprintf('============================================================\n');
fprintf('RESPONDING-CHANNEL REFERENCE CRITERIA COMPLETE\n');
fprintf('Nothing saved to disk. Full detail available in workspace variable: Responding\n');
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

function label = format_simultaneous(stim_channels)
stim_channels = stim_channels(stim_channels > 0);
if isempty(stim_channels)
    label = '(none)';
elseif numel(stim_channels) == 1
    label = sprintf('Ch%d',stim_channels(1));
else
    labels = arrayfun(@(c) sprintf('Ch%d',c),stim_channels,'UniformOutput',false);
    label = strjoin(labels,' + ');
end
end

function output = channel_list_text(values)
if isempty(values)
    output = '(none)';
else
    output = strtrim(sprintf('%d ',values));
end
end

function R = empty_channel_result()
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
R.total_post_spikes = NaN;
R.frac_post_trials = NaN;
R.pass_FR_rule = false;
R.pass_total_post_spikes = false;
R.pass_frac_trials = false;
R.pass_abs_post_FR = false;
end