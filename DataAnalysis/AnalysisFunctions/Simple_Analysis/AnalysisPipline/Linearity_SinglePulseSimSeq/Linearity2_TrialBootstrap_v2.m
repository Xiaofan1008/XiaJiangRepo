%% ========================================================================
%  LINEARITY ANALYSIS 2: TRIAL BOOTSTRAP
%
%  Independently resamples the cleaned A-alone, B-alone, and observed
%  combined-condition trials to estimate uncertainty in:
%
%      Difference = Observed - (A alone + B alone)
%
%  This script is independent of Analysis 1 and the Luke export workflow.
%  It loads the original experiment folders directly and creates no files.
%
%  Required external functions: Depth_s and loadTrig
%  Window convention: [start, end), excluding the right endpoint.
% ========================================================================

clear;

%% ============================= USER SETTINGS ============================

single_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Single4_new';
sim_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim4';
seq_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim4';

% Enter the same path for sim_folder and seq_folder when both conditions
% are stored together. The dataset will then be loaded only once.

analysis_functions_folder = ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE';

Electrode_Type = 2;
FS = 30000;

Pairs_To_Analyze = 'all';       % 'all' or detected pair numbers
Channels_To_Analyze = 'responding'; % 'responding', 'all', or numeric
Amplitudes_To_Analyze = [];
Include_Zero_Amplitude = false;

baseline_win_ms = [-50 -5];
post_win_ms = [2 20];
Sequential_PTD_ms = 5;
PTD_Tolerance_ms = 0.01;

Trial_Mode = 'all';             % 'all' or 'balanced'
Random_Seed = 1;

N_Bootstrap = 10000;
Confidence_Level = 95;          % percentile bootstrap interval

Sequential_Mode = 'separate';   % 'separate', 'merged', or 'both'
Curve_Groups_To_Plot = 'all';   % 'all', 'simultaneous', 'sequential'

Plot_Bootstrap_Difference = true;
Channels_Per_Figure = 12;
Print_Result_Table = true;

% Optional selected-channel average. Trial IDs are resampled jointly across
% channels, preserving the within-trial relationship between channels.
Run_Selected_Channel_Average = true;
Exclude_Stim_Contacts_From_Average = true;
Plot_Selected_Channel_Average = true;
Print_Average_Result_Table = true;

% Optional empirical single-trial distribution plots. The predicted
% distribution contains every A-trial + B-trial combination. This is
% separate from the bootstrap confidence interval for the mean difference.
Plot_Trial_Distributions = false;
Distribution_Channels = [];     % depth channels, e.g. [2 3]
Distribution_Amplitudes = [];   % e.g. [5 10]; empty uses all selected amps
Distribution_Conditions = {'Sim'}; % Sim, A_to_B, B_to_A

%% ============================== VALIDATION ===============================

validate_settings(single_folder,sim_folder,seq_folder, ...
    analysis_functions_folder,baseline_win_ms,post_win_ms, ...
    Sequential_PTD_ms,Trial_Mode,Sequential_Mode,N_Bootstrap, ...
    Confidence_Level,Channels_Per_Figure,FS);

addpath(genpath(analysis_functions_folder));
if exist('Depth_s','file') ~= 2
    error('Linearity2:MissingDepthFunction', ...
        'Depth_s was not found under analysis_functions_folder.');
end
if exist('loadTrig','file') ~= 2
    error('Linearity2:MissingTriggerFunction', ...
        'loadTrig was not found under analysis_functions_folder.');
end

depth_to_spike_channel = double(Depth_s(Electrode_Type));
depth_to_spike_channel = depth_to_spike_channel(:);
nDepthChannels = numel(depth_to_spike_channel);

fprintf('\nTrial-bootstrap linearity analysis\n');
fprintf('Baseline window: [%g, %g) ms\n',baseline_win_ms);
fprintf('Response window: [%g, %g) ms\n',post_win_ms);
fprintf('Bootstrap repetitions: %d\n',N_Bootstrap);
fprintf('Confidence level: %g%%\n',Confidence_Level);
fprintf('Trial mode: %s\n\n',upper(char(string(Trial_Mode))));

%% ============================== LOAD DATA ================================

Single = load_dataset(single_folder,'Single',FS,nDepthChannels);
Sim = load_dataset(sim_folder,'Simultaneous',FS,nDepthChannels);
if same_folder(sim_folder,seq_folder)
    Seq = Sim;
    Seq.role = 'Sequential';
    fprintf('Sim and Seq paths are identical; the paired dataset was loaded once.\n\n');
else
    Seq = load_dataset(seq_folder,'Sequential',FS,nDepthChannels);
end

%% ========================== DETECT VALID PAIRS ===========================

DetectedPairs = detect_complete_pairs(Single,Sim,Seq, ...
    Sequential_PTD_ms,PTD_Tolerance_ms);
if isempty(DetectedPairs)
    error('Linearity2:NoCompletePairs', ...
        ['No complete pair has A alone, B alone, Sim, A->B, and B->A ' ...
         'at the selected PTD.']);
end

fprintf('Detected complete stimulation pairs\n');
disp(pair_table(DetectedPairs));
pair_indices = resolve_pair_selection(Pairs_To_Analyze,numel(DetectedPairs));
[plot_codes,plot_labels,plot_colors] = resolve_display_conditions( ...
    Curve_Groups_To_Plot,Sequential_Mode);

all_amplitudes = unique([Single.trial_amp_uA(:);Sim.trial_amp_uA(:); ...
    Seq.trial_amp_uA(:)]).';
all_amplitudes = all_amplitudes(isfinite(all_amplitudes));
if ~Include_Zero_Amplitude
    all_amplitudes(abs(all_amplitudes) < 1e-9) = [];
end
if isempty(Amplitudes_To_Analyze)
    selected_amplitudes = all_amplitudes;
else
    selected_amplitudes = unique(double(Amplitudes_To_Analyze(:).'),'stable');
    if ~Include_Zero_Amplitude
        selected_amplitudes(abs(selected_amplitudes) < 1e-9) = [];
    end
end
if isempty(selected_amplitudes)
    error('Linearity2:NoAmplitudes','No amplitudes remain for analysis.');
end

%% ============================= MAIN ANALYSIS =============================

base_condition_codes = {'A','B','Sim','A_to_B','B_to_A'};
display_condition_codes = {'Sim','A_to_B','B_to_A','SeqMerged'};
nBaseConditions = numel(base_condition_codes);

BootstrapResults = struct([]);
ResultTables = cell(1,numel(pair_indices));
random_stream = RandStream('mt19937ar','Seed',Random_Seed);

for ipOut = 1:numel(pair_indices)
    detected_index = pair_indices(ipOut);
    Pair = DetectedPairs(detected_index);
    channels = resolve_channels(Channels_To_Analyze,Pair,Sim,Seq, ...
        nDepthChannels,Sequential_PTD_ms,PTD_Tolerance_ms);
    if isempty(channels)
        warning('Linearity2:NoChannels','No requested channels for pair %d.',detected_index);
        continue;
    end

    nChannels = numel(channels);
    nAmp = numel(selected_amplitudes);
    nDisplayConditions = numel(display_condition_codes);

    observed_mean = nan(nChannels,nAmp,nDisplayConditions);
    linear_mean = nan(nChannels,nAmp);
    difference = nan(nChannels,nAmp,nDisplayConditions);
    ci_low = nan(nChannels,nAmp,nDisplayConditions);
    ci_high = nan(nChannels,nAmp,nDisplayConditions);
    probability_below_zero = nan(nChannels,nAmp,nDisplayConditions);
    probability_above_zero = nan(nChannels,nAmp,nDisplayConditions);
    observed_n = nan(nChannels,nAmp,nDisplayConditions);
    single_n_A = nan(nChannels,nAmp);
    single_n_B = nan(nChannels,nAmp);
    interpretation = strings(nChannels,nAmp,nDisplayConditions);
    bootstrap_difference = cell(nChannels,nAmp,nDisplayConditions);
    trial_values = cell(nChannels,nAmp,nBaseConditions);
    all_clean_trial_ids = cell(nChannels,nAmp,nBaseConditions);
    all_clean_trial_values = cell(nChannels,nAmp,nBaseConditions);

    stimulation_depth_A = find_electrode_depth_index(Sim,Pair.A);
    stimulation_depth_B = find_electrode_depth_index(Sim,Pair.B);

    for jc = 1:nChannels
        depth_channel = channels(jc);
        spike_channel = depth_to_spike_channel(depth_channel);

        for ia = 1:nAmp
            amp = selected_amplitudes(ia);
            trial_ids = cell(1,nBaseConditions);
            trial_ids{1} = select_clean_trials(Single,Pair.A,Pair.key, ...
                'single',amp,0,depth_channel,PTD_Tolerance_ms);
            trial_ids{2} = select_clean_trials(Single,Pair.B,Pair.key, ...
                'single',amp,0,depth_channel,PTD_Tolerance_ms);
            trial_ids{3} = select_clean_trials(Sim,'',Pair.key, ...
                'sim',amp,0,depth_channel,PTD_Tolerance_ms);
            trial_ids{4} = select_clean_trials(Seq,Pair.A_to_B_key,Pair.key, ...
                'seq',amp,Sequential_PTD_ms,depth_channel,PTD_Tolerance_ms);
            trial_ids{5} = select_clean_trials(Seq,Pair.B_to_A_key,Pair.key, ...
                'seq',amp,Sequential_PTD_ms,depth_channel,PTD_Tolerance_ms);

            % Retain every clean trial for the joint across-channel
            % bootstrap. Individual-channel balancing is applied below to
            % a separate copy and therefore cannot break trial alignment.
            all_trial_ids_for_average = trial_ids;

            if is_balanced_mode(Trial_Mode)
                nBalanced = min(cellfun(@numel,trial_ids));
                if nBalanced < 1
                    trial_ids(:) = {[]};
                else
                    for ic = 1:nBaseConditions
                        chosen = randperm(random_stream,numel(trial_ids{ic}),nBalanced);
                        trial_ids{ic} = sort(trial_ids{ic}(chosen));
                    end
                end
            end

            data_sources = {Single,Single,Sim,Seq,Seq};
            for ic = 1:nBaseConditions
                trial_values{jc,ia,ic} = calculate_evoked_counts( ...
                    data_sources{ic},spike_channel,trial_ids{ic}, ...
                    baseline_win_ms,post_win_ms);
                all_clean_trial_ids{jc,ia,ic} = all_trial_ids_for_average{ic};
                all_clean_trial_values{jc,ia,ic} = calculate_evoked_counts( ...
                    data_sources{ic},spike_channel,all_trial_ids_for_average{ic}, ...
                    baseline_win_ms,post_win_ms);
            end

            A = trial_values{jc,ia,1};
            B = trial_values{jc,ia,2};
            SimValues = trial_values{jc,ia,3};
            AtoB = trial_values{jc,ia,4};
            BtoA = trial_values{jc,ia,5};
            single_n_A(jc,ia) = numel(A);
            single_n_B(jc,ia) = numel(B);

            [actual_linear,actual_observed,boot_delta,n_observed] = ...
                bootstrap_all_conditions(A,B,SimValues,AtoB,BtoA, ...
                N_Bootstrap,random_stream);

            linear_mean(jc,ia) = actual_linear;
            observed_mean(jc,ia,:) = reshape(actual_observed,1,1,[]);
            observed_n(jc,ia,:) = reshape(n_observed,1,1,[]);

            for ic = 1:nDisplayConditions
                samples = boot_delta(:,ic);
                samples = samples(isfinite(samples));
                if isempty(samples) || ~isfinite(actual_observed(ic)) || ...
                        ~isfinite(actual_linear)
                    interpretation(jc,ia,ic) = "Unavailable";
                    continue;
                end
                difference(jc,ia,ic) = actual_observed(ic)-actual_linear;
                alpha = (100-Confidence_Level)/2;
                ci_low(jc,ia,ic) = simple_percentile(samples,alpha);
                ci_high(jc,ia,ic) = simple_percentile(samples,100-alpha);
                probability_below_zero(jc,ia,ic) = mean(samples < 0);
                probability_above_zero(jc,ia,ic) = mean(samples > 0);
                bootstrap_difference{jc,ia,ic} = samples;
                interpretation(jc,ia,ic) = classify_interval( ...
                    ci_low(jc,ia,ic),ci_high(jc,ia,ic));
            end
        end
    end

    P = struct();
    P.detected_pair_index = detected_index;
    P.electrode_A = Pair.A;
    P.electrode_B = Pair.B;
    P.pair_key = Pair.key;
    P.stimulation_depth_A = stimulation_depth_A;
    P.stimulation_depth_B = stimulation_depth_B;
    P.depth_channels = channels;
    P.spike_channels = depth_to_spike_channel(channels);
    P.amplitudes_uA = selected_amplitudes;
    P.base_condition_codes = base_condition_codes;
    P.display_condition_codes = display_condition_codes;
    P.observed_mean = observed_mean;
    P.linear_mean = linear_mean;
    P.difference = difference;
    P.ci_low = ci_low;
    P.ci_high = ci_high;
    P.probability_below_zero = probability_below_zero;
    P.probability_above_zero = probability_above_zero;
    P.observed_n = observed_n;
    P.single_n_A = single_n_A;
    P.single_n_B = single_n_B;
    P.interpretation = interpretation;
    P.bootstrap_difference = bootstrap_difference;
    P.trial_values = trial_values;
    P.all_clean_trial_ids = all_clean_trial_ids;
    P.all_clean_trial_values = all_clean_trial_values;
    P.settings = struct('baseline_win_ms',baseline_win_ms, ...
        'post_win_ms',post_win_ms,'trial_mode',char(string(Trial_Mode)), ...
        'random_seed',Random_Seed,'n_bootstrap',N_Bootstrap, ...
        'confidence_level',Confidence_Level, ...
        'sequential_ptd_ms',Sequential_PTD_ms);

    AverageTable = table();
    if Run_Selected_Channel_Average
        AverageResult = calculate_selected_channel_average(P,Trial_Mode, ...
            N_Bootstrap,Confidence_Level,random_stream, ...
            Exclude_Stim_Contacts_From_Average);
        P.selected_channel_average = AverageResult;
        AverageTable = build_average_result_table(P,AverageResult,plot_codes);
    else
        P.selected_channel_average = struct([]);
    end

    if isempty(BootstrapResults)
        BootstrapResults = P;
    else
        BootstrapResults(end+1) = P;
    end

    ResultTables{ipOut} = build_result_table(P,plot_codes);
    fprintf('\nPair %d: %s and %s\n',detected_index,Pair.A,Pair.B);
    fprintf('Selected depth channels: %s\n',num2str(channels));
    if Print_Result_Table
        disp(ResultTables{ipOut});
    end
    if Run_Selected_Channel_Average && Print_Average_Result_Table
        fprintf('Selected-channel average bootstrap\n');
        disp(AverageTable);
    end

    if Plot_Bootstrap_Difference
        plot_bootstrap_differences(P,plot_codes,plot_labels,plot_colors, ...
            Channels_Per_Figure,Confidence_Level);
    end
    if Run_Selected_Channel_Average && Plot_Selected_Channel_Average
        plot_selected_channel_average(P,P.selected_channel_average, ...
            plot_codes,plot_labels,plot_colors,Confidence_Level, ...
            Exclude_Stim_Contacts_From_Average);
    end
    if Plot_Trial_Distributions
        plot_trial_distributions(P,Distribution_Channels, ...
            Distribution_Amplitudes,Distribution_Conditions);
    end
end

fprintf('\nAnalysis complete. No files were created or modified.\n');
fprintf('Results remain in the MATLAB workspace as BootstrapResults.\n');
fprintf(['Interpretations describe uncertainty across trials in these datasets; ' ...
    'they do not represent replication across animals.\n']);

%% ============================== FUNCTIONS ================================

function validate_settings(single_folder,sim_folder,seq_folder,functions_folder, ...
        baseline_win,post_win,seq_ptd,trial_mode,seq_mode,n_boot, ...
        confidence,page_size,FS)
folders = {single_folder,sim_folder,seq_folder,functions_folder};
names = {'single_folder','sim_folder','seq_folder','analysis_functions_folder'};
for k = 1:numel(folders)
    if ~isfolder(folders{k})
        error('Linearity2:FolderNotFound','%s not found: %s',names{k},folders{k});
    end
end
validate_window(baseline_win,'baseline_win_ms');
validate_window(post_win,'post_win_ms');
if baseline_win(2) > post_win(1)
    error('Linearity2:WindowOverlap','Baseline and response windows overlap.');
end
if ~isscalar(seq_ptd) || ~isfinite(seq_ptd) || seq_ptd <= 0
    error('Linearity2:InvalidPTD','Sequential_PTD_ms must be positive.');
end
if ~any(strcmpi(string(trial_mode),["all","balanced"]))
    error('Linearity2:InvalidTrialMode','Trial_Mode must be all or balanced.');
end
if ~any(strcmpi(string(seq_mode),["separate","merged","both"]))
    error('Linearity2:InvalidSeqMode', ...
        'Sequential_Mode must be separate, merged, or both.');
end
if ~isscalar(n_boot) || n_boot < 100 || fix(n_boot) ~= n_boot
    error('Linearity2:InvalidBootstrapCount','N_Bootstrap must be an integer >= 100.');
end
if ~isscalar(confidence) || confidence <= 0 || confidence >= 100
    error('Linearity2:InvalidConfidence','Confidence_Level must be between 0 and 100.');
end
if ~isscalar(page_size) || page_size < 1 || fix(page_size) ~= page_size
    error('Linearity2:InvalidPageSize','Channels_Per_Figure must be a positive integer.');
end
if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('Linearity2:InvalidFS','FS must be positive.');
end
end

function validate_window(window,name)
if ~isnumeric(window) || numel(window) ~= 2 || any(~isfinite(window)) || ...
        window(2) <= window(1)
    error('Linearity2:InvalidWindow','%s must be [start end], with end > start.',name);
end
end

function D = load_dataset(folder,role,FS,nDepthChannels)
folder = char(string(folder));
fprintf('Loading %s dataset\nFolder: %s\n',role,folder);

[spike_file,spike_variable] = find_spike_source(folder);
loaded_spikes = load(spike_file,spike_variable);
sp = loaded_spikes.(spike_variable);
if ~iscell(sp)
    error('Linearity2:SpikeNotCell','%s in %s is not a cell array.', ...
        spike_variable,spike_file);
end
fprintf('Spike file: %s\nSpike variable: %s (%d channels)\n', ...
    spike_file,spike_variable,numel(sp));

trig_files = clean_file_list(dir(fullfile(folder,'*.trig.dat')));
if isempty(trig_files)
    error('Linearity2:NoTrigger','No existing *.trig.dat file in %s.',folder);
elseif numel(trig_files) > 1
    error('Linearity2:AmbiguousTrigger','Multiple trigger files in %s.',folder);
end
trig = read_triggers(folder);
trig = double(trig(:));
fprintf('Trigger file: %s\nTriggers: %d\n', ...
    fullfile(trig_files(1).folder,trig_files(1).name),numel(trig));

exp_files = clean_file_list(dir(fullfile(folder,'*_exp_datafile_*.mat')));
if isempty(exp_files)
    error('Linearity2:NoExperimentFile','No experiment data file in %s.',folder);
elseif numel(exp_files) > 1
    error('Linearity2:AmbiguousExperimentFile','Multiple experiment files in %s.',folder);
end
experiment_file = fullfile(exp_files(1).folder,exp_files(1).name);
E = load(experiment_file,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required = {'StimParams','simultaneous_stim','E_MAP'};
for k = 1:numel(required)
    if ~isfield(E,required{k})
        error('Linearity2:MissingVariable','%s missing from %s.', ...
            required{k},experiment_file);
    end
end
fprintf('Experiment file: %s\n',experiment_file);

StimParams = E.StimParams;
simN = double(E.simultaneous_stim);
nRows = size(StimParams,1)-1;
if ~isscalar(simN) || simN < 1 || fix(simN) ~= simN || mod(nRows,simN) ~= 0
    error('Linearity2:InvalidStimParams','Invalid simultaneous_stim/StimParams in %s.', ...
        experiment_file);
end
nTrials = nRows/simN;
if isfield(E,'n_Trials') && ~isempty(E.n_Trials) && double(E.n_Trials) ~= nTrials
    error('Linearity2:TrialMismatch','n_Trials disagrees with StimParams in %s.', ...
        experiment_file);
end
if numel(trig) < nTrials
    warning('Linearity2:FewerTriggers', ...
        '%s has %d trials but %d triggers; unmatched trials will be excluded.', ...
        role,nTrials,numel(trig));
end

emap_labels = normalize_emap_labels(E.E_MAP);
trial_labels = cell(nTrials,1);
order_key = strings(nTrials,1);
pair_key = strings(nTrials,1);
trial_amp_uA = nan(nTrials,1);
trial_ptd_ms = zeros(nTrials,1);
for t = 1:nTrials
    first_row = 2+(t-1)*simN;
    labels = strings(1,simN);
    for pulse = 1:simN
        labels(pulse) = normalize_label(StimParams{first_row+pulse-1,1});
        if ~isempty(emap_labels) && ~any(labels(pulse) == emap_labels)
            error('Linearity2:LabelNotInMap','%s was not found in E_MAP for %s.', ...
                labels(pulse),experiment_file);
        end
    end
    trial_labels{t} = labels;
    order_key(t) = join(labels,">");
    pair_key(t) = join(sort(labels),"|");
    trial_amp_uA(t) = numeric_value(StimParams{first_row,16},'amplitude',t);
    if trial_amp_uA(t) == -1
        trial_amp_uA(t) = 0;
    end
    if simN > 1
        trial_ptd_ms(t) = numeric_value(StimParams{first_row+1,6},'PTD',t)/1000;
    end
end
[set_keys,~,set_index] = unique(order_key,'stable');

[BadCh,bad_ch_file] = load_bad_channels(folder);
[BadTrials,bad_trial_file] = load_bad_trials(folder,nDepthChannels);
[Responding,responding_file] = load_responding(folder);
print_optional_file('Bad-channel file',bad_ch_file);
if isempty(bad_trial_file)
    fprintf('Bad-trial file: none; no exclusions applied\n');
else
    fprintf('Bad-trial file: %s\n',bad_trial_file);
end
print_optional_file('Responding-channel file',responding_file);
fprintf('\n');

D=struct();
D.role=role; D.folder=folder; D.FS=FS; D.sp=sp;
D.spike_file=spike_file; D.spike_variable=spike_variable; D.trig=trig;
D.experiment_file=experiment_file; D.simN=simN; D.E_MAP=E.E_MAP;
D.emap_labels=emap_labels; D.nTrials=nTrials; D.trial_labels=trial_labels;
D.order_key=order_key; D.pair_key=pair_key; D.trial_amp_uA=trial_amp_uA;
D.trial_ptd_ms=trial_ptd_ms; D.set_keys=set_keys; D.set_index=set_index;
D.BadCh=BadCh; D.BadTrials=BadTrials; D.Responding=Responding;
end

function [file_path,variable_name] = find_spike_source(folder)
files = clean_file_list(dir(fullfile(folder,'*sp*.mat')));
keep = true(size(files));
for k = 1:numel(files)
    name = lower(files(k).name);
    keep(k) = ~contains(name,'spikecount') && ~contains(name,'result');
end
files = files(keep);
priority = {'sp_corr','sp_SSD','sp_clipped','sp'};
for ip = 1:numel(priority)
    matches = strings(0,1);
    for k = 1:numel(files)
        candidate = fullfile(files(k).folder,files(k).name);
        info = whos('-file',candidate);
        if any(strcmp({info.name},priority{ip}))
            matches(end+1,1) = string(candidate); %#ok<AGROW>
        end
    end
    if numel(matches) > 1
        error('Linearity2:AmbiguousSpikeFile', ...
            'Multiple files contain %s in %s.',priority{ip},folder);
    elseif isscalar(matches)
        file_path = char(matches);
        variable_name = priority{ip};
        return;
    end
end
error('Linearity2:NoSpikeVariable', ...
    'No sp_corr, sp_SSD, sp_clipped, or sp variable in %s.',folder);
end

function trig = read_triggers(folder)
original = pwd;
cleanup = onCleanup(@() cd(original));
cd(folder);
trig = loadTrig(0);
clear cleanup;
cd(original);
end

function labels = normalize_emap_labels(E_MAP)
labels = strings(0,1);
values = E_MAP(:);
if numel(values) >= 2
    values = values(2:end);
end
for k = 1:numel(values)
    if iscell(values), value = values{k}; else, value = values(k); end
    label = normalize_label(value);
    if strlength(label) > 0
        labels(end+1,1) = label; %#ok<AGROW>
    end
end
labels = unique(labels,'stable');
end

function label = normalize_label(value)
if iscell(value) && isscalar(value), value = value{1}; end
if ischar(value) || (isstring(value) && isscalar(value))
    label = upper(strtrim(string(value)));
elseif isnumeric(value) && isscalar(value) && isfinite(value)
    label = string(value);
else
    error('Linearity2:InvalidLabel','A stimulation label is invalid.');
end
label = regexprep(label,'\s+','');
end

function value = numeric_value(raw,name,trial)
if iscell(raw) && isscalar(raw), raw = raw{1}; end
if isnumeric(raw) && isscalar(raw)
    value = double(raw);
elseif ischar(raw) || (isstring(raw) && isscalar(raw))
    value = str2double(raw);
else
    value = NaN;
end
if ~isfinite(value)
    error('Linearity2:InvalidStimValue','Invalid %s in trial %d.',name,trial);
end
end

function [BadCh,file_path] = load_bad_channels(folder)
files = find_named_mat_files(folder,'badchannels');
if isempty(files)
    BadCh = []; file_path = ''; return;
elseif numel(files) > 1
    error('Linearity2:AmbiguousBadChannels','Multiple bad-channel files in %s.',folder);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if isfield(S,'BadCh_perSet'), BadCh = S.BadCh_perSet;
elseif isfield(S,'BadCh'), BadCh = S.BadCh;
else, error('Linearity2:BadChannelVariable','No bad-channel variable in %s.',file_path);
end
end

function [BadTrials,file_path] = load_bad_trials(folder,nDepth)
preferred = find_named_mat_files(folder,'simseqbadtrials');
if numel(preferred) > 1
    error('Linearity2:AmbiguousBadTrials','Multiple SimSeq bad-trial files in %s.',folder);
elseif isscalar(preferred)
    files = preferred;
else
    files = find_named_mat_files(folder,'badtrials');
    keep = true(size(files));
    for k = 1:numel(files)
        normalized = regexprep(lower(files(k).name),'[^a-z0-9]','');
        keep(k) = ~contains(normalized,'simseqbadtrials');
    end
    files = files(keep);
end
if isempty(files)
    BadTrials = cell(nDepth,1); file_path = ''; return;
elseif numel(files) > 1
    error('Linearity2:AmbiguousBadTrials','Multiple bad-trial files in %s.',folder);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if ~isfield(S,'BadTrials') || ~iscell(S.BadTrials)
    error('Linearity2:BadTrialVariable','BadTrials cell array missing from %s.',file_path);
end
BadTrials = S.BadTrials(:);
if numel(BadTrials) < nDepth
    BadTrials(end+1:nDepth) = {[]};
end
end

function [Responding,file_path] = load_responding(folder)
files = find_named_mat_files(folder,'respondingchannels');
if isempty(files)
    Responding = []; file_path = ''; return;
elseif numel(files) > 1
    error('Linearity2:AmbiguousResponding','Multiple responding files in %s.',folder);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if ~isfield(S,'Responding')
    error('Linearity2:RespondingVariable','Responding missing from %s.',file_path);
end
Responding = S.Responding;
end

function files = find_named_mat_files(folder,token)
files = clean_file_list(dir(fullfile(folder,'*.mat')));
keep = false(size(files));
for k = 1:numel(files)
    normalized = regexprep(lower(files(k).name),'[^a-z0-9]','');
    keep(k) = contains(normalized,lower(token));
end
files = files(keep);
end

function files = clean_file_list(files)
if isempty(files), return; end
keep = true(size(files));
for k = 1:numel(files), keep(k) = ~startsWith(files(k).name,'._'); end
files = files(keep);
end

function print_optional_file(label,path)
if isempty(path), fprintf('%s: none\n',label);
else, fprintf('%s: %s\n',label,path);
end
end

function tf = same_folder(a,b)
tf = strcmp(strip_separator(char(string(a))),strip_separator(char(string(b))));
end

function out = strip_separator(in)
out = in;
while numel(out) > 1 && any(out(end) == ['/' '\'])
    out(end) = [];
end
end

function Pairs = detect_complete_pairs(Single,Sim,Seq,target_ptd,tol)
sim_mask = abs(Sim.trial_ptd_ms) < tol & cellfun(@numel,Sim.trial_labels) == 2;
seq_mask = abs(Seq.trial_ptd_ms-target_ptd) < tol & ...
    cellfun(@numel,Seq.trial_labels) == 2;
candidate = intersect(unique(Sim.pair_key(sim_mask),'stable'), ...
    unique(Seq.pair_key(seq_mask),'stable'),'stable');
single_labels = unique(Single.order_key(cellfun(@numel,Single.trial_labels)==1));
Pairs = repmat(struct('key','', 'A','', 'B','', 'A_to_B_key','', ...
    'B_to_A_key','', 'sim_set_indices',[], 'a_to_b_set_indices',[], ...
    'b_to_a_set_indices',[]),1,0);
for k = 1:numel(candidate)
    labels = sort(split(candidate(k),'|'));
    if numel(labels) ~= 2 || labels(1) == labels(2), continue; end
    A = labels(1); B = labels(2);
    AtoB = A+">"+B; BtoA = B+">"+A;
    if ~any(single_labels==A) || ~any(single_labels==B) || ...
            ~any(seq_mask & Seq.order_key==AtoB) || ...
            ~any(seq_mask & Seq.order_key==BtoA)
        continue;
    end
    P = struct('key',char(candidate(k)),'A',char(A),'B',char(B), ...
        'A_to_B_key',char(AtoB),'B_to_A_key',char(BtoA), ...
        'sim_set_indices',unique(Sim.set_index(sim_mask & Sim.pair_key==candidate(k))).', ...
        'a_to_b_set_indices',unique(Seq.set_index(seq_mask & Seq.order_key==AtoB)).', ...
        'b_to_a_set_indices',unique(Seq.set_index(seq_mask & Seq.order_key==BtoA)).');
    Pairs(end+1) = P; %#ok<AGROW>
end
end

function T = pair_table(Pairs)
n = numel(Pairs); Pair = (1:n).';
ElectrodeA = strings(n,1); ElectrodeB = strings(n,1);
SimSets = strings(n,1); AtoBSets = strings(n,1); BtoASets = strings(n,1);
for k = 1:n
    ElectrodeA(k)=string(Pairs(k).A); ElectrodeB(k)=string(Pairs(k).B);
    SimSets(k)=vector_text(Pairs(k).sim_set_indices);
    AtoBSets(k)=vector_text(Pairs(k).a_to_b_set_indices);
    BtoASets(k)=vector_text(Pairs(k).b_to_a_set_indices);
end
T = table(Pair,ElectrodeA,ElectrodeB,SimSets,AtoBSets,BtoASets);
end

function text = vector_text(values)
if isempty(values), text="none"; else, text=string(strtrim(num2str(values))); end
end

function indices = resolve_pair_selection(request,nPairs)
if ischar(request) || (isstring(request) && isscalar(request))
    if ~strcmpi(strtrim(string(request)),'all')
        error('Linearity2:PairSelection','Text pair selection must be all.');
    end
    indices = 1:nPairs;
elseif isnumeric(request)
    indices = unique(double(request(:).'),'stable');
    if isempty(indices) || any(indices<1 | indices>nPairs | fix(indices)~=indices)
        error('Linearity2:PairSelection','Invalid detected pair number.');
    end
else
    error('Linearity2:PairSelection','Pairs_To_Analyze must be all or numeric.');
end
end

function channels = resolve_channels(request,Pair,Sim,Seq,nDepth,target_ptd,tol)
if isnumeric(request)
    channels = unique(double(request(:).'),'stable');
    if isempty(channels) || any(channels<1 | channels>nDepth | fix(channels)~=channels)
        error('Linearity2:ChannelSelection','Invalid depth channel selection.');
    end
    return;
end
switch lower(strtrim(char(string(request))))
    case 'all'
        channels = 1:nDepth;
    case 'responding'
        sim_mask = responding_mask(Sim,Pair.sim_set_indices,0,tol,nDepth);
        seq_sets = unique([Pair.a_to_b_set_indices Pair.b_to_a_set_indices]);
        seq_mask = responding_mask(Seq,seq_sets,target_ptd,tol,nDepth);
        channels = find(sim_mask | seq_mask).';
        if isempty(channels)
            error('Linearity2:NoRespondingChannels', ...
                'No responding labels found for %s.',Pair.key);
        end
    otherwise
        error('Linearity2:ChannelSelection','Use responding, all, or numeric channels.');
end
end

function mask = responding_mask(D,set_indices,target_ptd,tol,nDepth)
mask = false(nDepth,1);
if isempty(D.Responding) || ~isfield(D.Responding,'set'), return; end
for set_idx = set_indices
    if set_idx<1 || set_idx>numel(D.Responding.set), continue; end
    ptd_values = unique(D.trial_ptd_ms(D.set_index==set_idx));
    pi = find(abs(ptd_values-target_ptd)<tol,1);
    if isempty(pi), continue; end
    amps = D.Responding.set(set_idx).amp;
    for ia = 1:numel(amps)
        if ~isfield(amps(ia),'ptd') || numel(amps(ia).ptd)<pi || ...
                ~isfield(amps(ia).ptd(pi),'channel'), continue; end
        C = amps(ia).ptd(pi).channel;
        for ch = 1:min(numel(C),nDepth)
            if isfield(C(ch),'is_responsive') && logical(C(ch).is_responsive)
                mask(ch)=true;
            end
        end
    end
end
end

function ids = select_clean_trials(D,ordered_key,pair_key,mode,amp,ptd,depth,tol)
mask = abs(D.trial_amp_uA-amp)<1e-6 & abs(D.trial_ptd_ms-ptd)<tol;
if strcmpi(mode,'sim')
    mask = mask & D.pair_key==string(pair_key);
else
    mask = mask & D.order_key==string(ordered_key);
end
trigger_available = (1:D.nTrials).' <= numel(D.trig);
ids = find(mask & trigger_available);
if isempty(ids), return; end
keep = true(size(ids));
for k = 1:numel(ids)
    keep(k)=~is_bad_channel(D.BadCh,D.set_index(ids(k)),depth);
end
ids = ids(keep);
if depth<=numel(D.BadTrials) && ~isempty(D.BadTrials{depth})
    ids = ids(~ismember(ids,double(D.BadTrials{depth}(:))));
end
ids = ids(:).';
end

function tf = is_bad_channel(BadCh,set_idx,depth)
tf=false;
if isempty(BadCh), return; end
if iscell(BadCh)
    if set_idx<=numel(BadCh) && ~isempty(BadCh{set_idx})
        tf=ismember(depth,double(BadCh{set_idx}(:)));
    end
elseif isnumeric(BadCh) || islogical(BadCh)
    tf=ismember(depth,double(BadCh(:)));
end
end

function values = calculate_evoked_counts(D,spike_channel,ids,baseline_win,post_win)
if isempty(ids) || spike_channel<1 || spike_channel>numel(D.sp) || ...
        isempty(D.sp{spike_channel})
    values=nan(0,1); return;
end
spike_data=D.sp{spike_channel};
if ~isnumeric(spike_data) || size(spike_data,2)<1
    values=nan(0,1); return;
end
times=double(spike_data(:,1));
base_duration=diff(baseline_win); post_duration=diff(post_win);
values=nan(numel(ids),1);
for k=1:numel(ids)
    tr=ids(k); t0=D.trig(tr)/D.FS*1000;
    base=sum(times>=t0+baseline_win(1) & times<t0+baseline_win(2));
    post=sum(times>=t0+post_win(1) & times<t0+post_win(2));
    values(k)=post-base*post_duration/base_duration;
end
values=values(isfinite(values));
end

function [linear,observed,boot_delta,n_observed] = bootstrap_all_conditions( ...
        A,B,Sim,AtoB,BtoA,nBoot,stream)
linear=NaN; observed=nan(1,4); boot_delta=nan(nBoot,4);
n_observed=[numel(Sim) numel(AtoB) numel(BtoA) NaN];
if isempty(A) || isempty(B), return; end
linear=mean(A)+mean(B);
Aboot=bootstrap_means(A,nBoot,stream);
Bboot=bootstrap_means(B,nBoot,stream);
linear_boot=Aboot+Bboot;
observed_sets={Sim,AtoB,BtoA};
observed_boot=nan(nBoot,3);
for k=1:3
    x=observed_sets{k};
    if isempty(x), continue; end
    observed(k)=mean(x);
    observed_boot(:,k)=bootstrap_means(x,nBoot,stream);
    boot_delta(:,k)=observed_boot(:,k)-linear_boot;
end
if ~isempty(AtoB) && ~isempty(BtoA)
    observed(4)=(observed(2)+observed(3))/2;
    boot_delta(:,4)=(observed_boot(:,2)+observed_boot(:,3))/2-linear_boot;
end
end

function means = bootstrap_means(values,nBoot,stream)
values=double(values(:)); n=numel(values);
indices=randi(stream,n,[n nBoot]);
means=mean(values(indices),1).';
end

function value = simple_percentile(samples,percent)
x=sort(double(samples(:)));
n=numel(x);
if n==0, value=NaN; return; end
position=1+(n-1)*percent/100;
lo=floor(position); hi=ceil(position);
if lo==hi, value=x(lo);
else, value=x(lo)+(position-lo)*(x(hi)-x(lo));
end
end

function label = classify_interval(low,high)
if ~isfinite(low) || ~isfinite(high), label="Unavailable";
elseif high<0, label="Sublinear";
elseif low>0, label="Supralinear";
else, label="Consistent with linearity";
end
end

function depth = find_electrode_depth_index(D,label)
depth=NaN; target=normalize_label(label); values=D.E_MAP(:);
if numel(values)<2, return; end
for k=2:numel(values)
    if iscell(values), value=values{k}; else, value=values(k); end
    try
        current=normalize_label(value);
    catch
        continue;
    end
    if current==target
        depth=k-1;
        return;
    end
end
end

function R = calculate_selected_channel_average(P,trial_mode,nBoot,confidence, ...
        stream,exclude_stim_contacts)
nAmp=numel(P.amplitudes_uA); nDisplay=numel(P.display_condition_codes);
R.linear_mean=nan(1,nAmp);
R.observed_mean=nan(nAmp,nDisplay);
R.difference=nan(nAmp,nDisplay);
R.ci_low=nan(nAmp,nDisplay);
R.ci_high=nan(nAmp,nDisplay);
R.probability_below_zero=nan(nAmp,nDisplay);
R.probability_above_zero=nan(nAmp,nDisplay);
R.interpretation=strings(nAmp,nDisplay);
R.bootstrap_difference=cell(nAmp,nDisplay);
R.n_channels=zeros(1,nAmp);
R.depth_channels=cell(1,nAmp);
R.n_trials=nan(nAmp,5);
R.exclude_stimulation_contacts=logical(exclude_stim_contacts);

for ia=1:nAmp
    eligible=true(numel(P.depth_channels),1);
    if exclude_stim_contacts
        stimulation_depths=[P.stimulation_depth_A P.stimulation_depth_B];
        stimulation_depths=stimulation_depths(isfinite(stimulation_depths));
        eligible=eligible & ~ismember(P.depth_channels(:),stimulation_depths);
    end
    for jc=1:numel(P.depth_channels)
        if ~eligible(jc), continue; end
        for ic=1:5
            ids=P.all_clean_trial_ids{jc,ia,ic};
            values=P.all_clean_trial_values{jc,ia,ic};
            if isempty(ids) || numel(ids)~=numel(values)
                eligible(jc)=false;
                break;
            end
        end
    end
    channel_indices=find(eligible);
    if isempty(channel_indices)
        R.interpretation(ia,:)="Unavailable";
        continue;
    end

    common_ids=cell(1,5);
    for ic=1:5
        ids=P.all_clean_trial_ids{channel_indices(1),ia,ic};
        common_ids{ic}=double(ids(:).');
        for j=2:numel(channel_indices)
            ids=P.all_clean_trial_ids{channel_indices(j),ia,ic};
            common_ids{ic}=intersect(common_ids{ic},double(ids(:).'),'stable');
        end
    end
    if any(cellfun(@isempty,common_ids))
        R.interpretation(ia,:)="Unavailable";
        continue;
    end

    if is_balanced_mode(trial_mode)
        nBalanced=min(cellfun(@numel,common_ids));
        for ic=1:5
            chosen=randperm(stream,numel(common_ids{ic}),nBalanced);
            common_ids{ic}=sort(common_ids{ic}(chosen));
        end
    end

    condition_values=cell(1,5);
    matrices_are_complete=true;
    for ic=1:5
        matrix=nan(numel(common_ids{ic}),numel(channel_indices));
        for j=1:numel(channel_indices)
            jc=channel_indices(j);
            ids=double(P.all_clean_trial_ids{jc,ia,ic}(:));
            values=double(P.all_clean_trial_values{jc,ia,ic}(:));
            [found,positions]=ismember(common_ids{ic}(:),ids);
            if any(~found)
                matrices_are_complete=false;
                break;
            end
            matrix(:,j)=values(positions);
        end
        if ~matrices_are_complete, break; end
        % All entries are present, so this gives every channel equal weight
        % and preserves the channel pattern within each resampled trial.
        condition_values{ic}=mean(matrix,2);
        R.n_trials(ia,ic)=size(matrix,1);
    end
    if ~matrices_are_complete
        R.interpretation(ia,:)="Unavailable";
        continue;
    end

    [linear,observed,boot_delta,~]=bootstrap_all_conditions( ...
        condition_values{1},condition_values{2},condition_values{3}, ...
        condition_values{4},condition_values{5},nBoot,stream);
    R.linear_mean(ia)=linear;
    R.observed_mean(ia,:)=observed;
    R.n_channels(ia)=numel(channel_indices);
    R.depth_channels{ia}=P.depth_channels(channel_indices);
    alpha=(100-confidence)/2;
    for ic=1:nDisplay
        samples=boot_delta(:,ic);
        samples=samples(isfinite(samples));
        if isempty(samples) || ~isfinite(observed(ic)) || ~isfinite(linear)
            R.interpretation(ia,ic)="Unavailable";
            continue;
        end
        R.difference(ia,ic)=observed(ic)-linear;
        R.ci_low(ia,ic)=simple_percentile(samples,alpha);
        R.ci_high(ia,ic)=simple_percentile(samples,100-alpha);
        R.probability_below_zero(ia,ic)=mean(samples<0);
        R.probability_above_zero(ia,ic)=mean(samples>0);
        R.interpretation(ia,ic)=classify_interval(R.ci_low(ia,ic),R.ci_high(ia,ic));
        R.bootstrap_difference{ia,ic}=samples;
    end
end
end

function tf = is_balanced_mode(request)
tf=strcmpi(strtrim(char(string(request))),'balanced');
end

function [codes,labels,colors] = resolve_display_conditions(groups_request,seq_mode)
groups=normalize_groups(groups_request);
include_sim=any(groups=="all" | groups=="simultaneous" | groups=="sim");
include_seq=any(groups=="all" | groups=="sequential" | groups=="seq");
codes=strings(0,1); labels=strings(0,1); colors=zeros(0,3);
if include_sim
    codes(end+1)="Sim"; labels(end+1)="A+B simultaneous";
    colors(end+1,:)=[0.00 0.35 0.80];
end
if include_seq
    switch lower(char(string(seq_mode)))
        case 'separate'
            codes=[codes;"A_to_B";"B_to_A"];
            labels=[labels;"A then B";"B then A"];
            colors=[colors;0.85 0.33 0.10;0.55 0.20 0.65];
        case 'merged'
            codes(end+1)="SeqMerged"; labels(end+1)="Sequential order average";
            colors(end+1,:)=[0.20 0.60 0.25];
        case 'both'
            codes=[codes;"A_to_B";"B_to_A";"SeqMerged"];
            labels=[labels;"A then B";"B then A";"Sequential order average"];
            colors=[colors;0.85 0.33 0.10;0.55 0.20 0.65;0.20 0.60 0.25];
    end
end
if isempty(codes), error('Linearity2:NoCurves','No curve group selected.'); end
end

function groups = normalize_groups(request)
if ischar(request) || (isstring(request) && isscalar(request))
    groups=lower(strtrim(string(request)));
elseif iscell(request) || isstring(request)
    groups=lower(strtrim(string(request(:))));
else
    error('Linearity2:CurveGroups','Invalid Curve_Groups_To_Plot.');
end
valid=["all","simultaneous","sim","sequential","seq"];
if isempty(groups) || any(~ismember(groups,valid))
    error('Linearity2:CurveGroups','Use all, simultaneous, or sequential.');
end
end

function idx = display_index(P,code)
idx=find(string(P.display_condition_codes)==string(code),1);
if isempty(idx), error('Linearity2:DisplayCode','Unknown condition %s.',code); end
end

function T = build_result_table(P,codes)
nRows=numel(P.depth_channels)*numel(P.amplitudes_uA)*numel(codes);
Pair=repmat(string(P.electrode_A)+" + "+string(P.electrode_B),nRows,1);
DepthChannel=nan(nRows,1); Amplitude_uA=nan(nRows,1);
Condition=strings(nRows,1); N_A=nan(nRows,1); N_B=nan(nRows,1);
N_Observed=nan(nRows,1);
Observed=nan(nRows,1); Linear=nan(nRows,1); Difference=nan(nRows,1);
CI_Low=nan(nRows,1); CI_High=nan(nRows,1);
P_Below_Zero=nan(nRows,1); P_Above_Zero=nan(nRows,1);
Interpretation=strings(nRows,1); row=0;
for jc=1:numel(P.depth_channels)
    for ia=1:numel(P.amplitudes_uA)
        for k=1:numel(codes)
            row=row+1; idx=display_index(P,codes(k));
            DepthChannel(row)=P.depth_channels(jc);
            Amplitude_uA(row)=P.amplitudes_uA(ia);
            Condition(row)=codes(k);
            N_A(row)=P.single_n_A(jc,ia); N_B(row)=P.single_n_B(jc,ia);
            N_Observed(row)=P.observed_n(jc,ia,idx);
            Observed(row)=P.observed_mean(jc,ia,idx);
            Linear(row)=P.linear_mean(jc,ia);
            Difference(row)=P.difference(jc,ia,idx);
            CI_Low(row)=P.ci_low(jc,ia,idx); CI_High(row)=P.ci_high(jc,ia,idx);
            P_Below_Zero(row)=P.probability_below_zero(jc,ia,idx);
            P_Above_Zero(row)=P.probability_above_zero(jc,ia,idx);
            Interpretation(row)=P.interpretation(jc,ia,idx);
        end
    end
end
T=table(Pair,DepthChannel,Amplitude_uA,Condition,N_A,N_B,N_Observed,Observed,Linear, ...
    Difference,CI_Low,CI_High,P_Below_Zero,P_Above_Zero,Interpretation);
end

function T = build_average_result_table(P,R,codes)
nRows=numel(P.amplitudes_uA)*numel(codes);
Pair=repmat(string(P.electrode_A)+" + "+string(P.electrode_B),nRows,1);
Amplitude_uA=nan(nRows,1); Condition=strings(nRows,1);
N_Channels=nan(nRows,1); N_A=nan(nRows,1); N_B=nan(nRows,1);
N_Observed=nan(nRows,1); Observed=nan(nRows,1); Linear=nan(nRows,1);
Difference=nan(nRows,1); CI_Low=nan(nRows,1); CI_High=nan(nRows,1);
P_Below_Zero=nan(nRows,1); P_Above_Zero=nan(nRows,1);
Interpretation=strings(nRows,1); row=0;
for ia=1:numel(P.amplitudes_uA)
    for k=1:numel(codes)
        row=row+1; idx=display_index(P,codes(k));
        Amplitude_uA(row)=P.amplitudes_uA(ia); Condition(row)=codes(k);
        N_Channels(row)=R.n_channels(ia); N_A(row)=R.n_trials(ia,1);
        N_B(row)=R.n_trials(ia,2);
        if idx<=3, N_Observed(row)=R.n_trials(ia,idx+2); end
        Observed(row)=R.observed_mean(ia,idx); Linear(row)=R.linear_mean(ia);
        Difference(row)=R.difference(ia,idx); CI_Low(row)=R.ci_low(ia,idx);
        CI_High(row)=R.ci_high(ia,idx);
        P_Below_Zero(row)=R.probability_below_zero(ia,idx);
        P_Above_Zero(row)=R.probability_above_zero(ia,idx);
        Interpretation(row)=R.interpretation(ia,idx);
    end
end
T=table(Pair,Amplitude_uA,Condition,N_Channels,N_A,N_B,N_Observed, ...
    Observed,Linear,Difference,CI_Low,CI_High,P_Below_Zero, ...
    P_Above_Zero,Interpretation);
end

function plot_bootstrap_differences(P,codes,labels,colors,page_size,confidence)
nChannels=numel(P.depth_channels); nPages=ceil(nChannels/page_size);
for page=1:nPages
    first=(page-1)*page_size+1; last=min(page*page_size,nChannels);
    selected=first:last; nTiles=numel(selected);
    nColumns=ceil(sqrt(nTiles)); nRows=ceil(nTiles/nColumns);
    figure('Color','w','Position',[50 50 1500 850]);
    tiledlayout(nRows,nColumns,'TileSpacing','compact','Padding','compact');
    sgtitle(sprintf('%s + %s: bootstrap nonlinear difference (%g%% CI), page %d of %d', ...
        P.electrode_A,P.electrode_B,confidence,page,nPages),'FontWeight','bold');
    for jc=selected
        ax=nexttile; hold(ax,'on');
        yline(ax,0,'--','Color',[0.25 0.25 0.25],'HandleVisibility','off');
        for k=1:numel(codes)
            idx=display_index(P,codes(k));
            x=P.amplitudes_uA(:);
            y=squeeze(P.difference(jc,:,idx)); y=y(:);
            low=squeeze(P.ci_low(jc,:,idx)); low=low(:);
            high=squeeze(P.ci_high(jc,:,idx)); high=high(:);
            lower_error=y-low; upper_error=high-y;
            valid=isfinite(x) & isfinite(y) & ...
                isfinite(lower_error) & isfinite(upper_error);
            if any(valid)
                errorbar(ax,x(valid),y(valid), ...
                    lower_error(valid),upper_error(valid),'-o','Color',colors(k,:), ...
                    'MarkerFaceColor',colors(k,:),'LineWidth',1.5,'CapSize',6, ...
                    'DisplayName',char(labels(k)));
            end
        end
        title(ax,sprintf('Depth channel %d',P.depth_channels(jc)));
        xlabel(ax,'Amplitude (uA)'); ylabel(ax,'Observed - linear'); box(ax,'off');
        if jc==selected(1), legend(ax,'Location','best','Box','off'); end
    end
end
end

function plot_selected_channel_average(P,R,codes,labels,colors,confidence,excluded_stim)
figure('Color','w','Position',[100 100 800 560]);
ax=axes(); hold(ax,'on');
yline(ax,0,'--','Color',[0.25 0.25 0.25],'HandleVisibility','off');
x=P.amplitudes_uA(:);
for k=1:numel(codes)
    idx=display_index(P,codes(k));
    y=R.difference(:,idx); low=R.ci_low(:,idx); high=R.ci_high(:,idx);
    lower_error=y-low; upper_error=high-y;
    valid=isfinite(x) & isfinite(y) & isfinite(lower_error) & isfinite(upper_error);
    if any(valid)
        errorbar(ax,x(valid),y(valid),lower_error(valid),upper_error(valid), ...
            '-o','Color',colors(k,:),'MarkerFaceColor',colors(k,:), ...
            'LineWidth',1.7,'MarkerSize',6,'CapSize',7, ...
            'DisplayName',char(labels(k)));
    end
end
if excluded_stim
    contact_text='stimulation contacts excluded';
else
    contact_text='stimulation contacts included';
end
title(ax,sprintf('%s + %s: selected-channel mean difference (%g%% CI; %s)', ...
    P.electrode_A,P.electrode_B,confidence,contact_text), ...
    'Interpreter','none','FontWeight','bold');
xlabel(ax,'Amplitude (uA)'); ylabel(ax,'Mean observed - linear');
legend(ax,'Location','best','Box','off'); box(ax,'off');
yl=ylim(ax);
for ia=1:numel(x)
    if R.n_channels(ia)>0
        text(ax,x(ia),yl(1)+0.03*(yl(2)-yl(1)),sprintf('nCh=%d',R.n_channels(ia)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',8,'Color',[0.3 0.3 0.3]);
    end
end
end

function plot_trial_distributions(P,requested_channels,requested_amps,requested_conditions)
if isempty(requested_channels)
    warning('Linearity2:NoDistributionChannels', ...
        'Plot_Trial_Distributions is true but Distribution_Channels is empty.');
    return;
end
channels=intersect(P.depth_channels,double(requested_channels(:).'),'stable');
if isempty(channels)
    warning('Linearity2:DistributionChannelsMissing', ...
        'None of the requested distribution channels are in this pair.');
    return;
end
if isempty(requested_amps), amps=P.amplitudes_uA;
else, amps=intersect(P.amplitudes_uA,double(requested_amps(:).'),'stable');
end
conditions=string(requested_conditions(:));
valid_conditions=["Sim","A_to_B","B_to_A"];
if any(~ismember(conditions,valid_conditions))
    error('Linearity2:DistributionConditions', ...
        'Distribution conditions must be Sim, A_to_B, or B_to_A.');
end
base_indices=[3 4 5];
for depth=channels
    jc=find(P.depth_channels==depth,1);
    for amp=amps
        ia=find(abs(P.amplitudes_uA-amp)<1e-9,1);
        A=P.trial_values{jc,ia,1}; B=P.trial_values{jc,ia,2};
        if isempty(A) || isempty(B), continue; end
        predicted=reshape(A(:)+B(:).',[],1);
        for k=1:numel(conditions)
            base_idx=base_indices(valid_conditions==conditions(k));
            observed=P.trial_values{jc,ia,base_idx};
            if isempty(observed), continue; end
            figure('Color','w','Position',[100 100 760 500]); hold on;
            histogram(predicted,'Normalization','probability','DisplayStyle','stairs', ...
                'LineWidth',1.8,'EdgeColor',[0 0 0],'DisplayName','A+B combinations');
            histogram(observed,'Normalization','probability','FaceAlpha',0.35, ...
                'FaceColor',[0.00 0.35 0.80],'EdgeColor',[0.00 0.35 0.80], ...
                'DisplayName',char(conditions(k)));
            xline(mean(predicted),'--k','LineWidth',1.5,'DisplayName','Predicted mean');
            xline(mean(observed),'--','Color',[0.00 0.35 0.80], ...
                'LineWidth',1.5,'DisplayName','Observed mean');
            title(sprintf('%s + %s, depth %d, %g uA, %s', ...
                P.electrode_A,P.electrode_B,depth,amp,conditions(k)), ...
                'Interpreter','none');
            xlabel('Baseline-corrected spikes/trial'); ylabel('Probability');
            legend('Location','best','Box','off'); box off;
        end
    end
end
end
