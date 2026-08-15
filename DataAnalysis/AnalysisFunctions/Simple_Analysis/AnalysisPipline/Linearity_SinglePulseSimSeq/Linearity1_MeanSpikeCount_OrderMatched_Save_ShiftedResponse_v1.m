%% ========================================================================
%  LINEARITY ANALYSIS 1: MEAN SPIKE COUNT VERSUS AMPLITUDE
%
%  Purpose
%    Compare observed dual-electrode responses with condition-specific
%    additive references obtained from the two single-electrode responses:
%
%       Sim linear = unshifted A alone + unshifted B alone
%       Seq linear = first single component + second single component
%                    after shifting its spike times by the PTD
%
%  Main properties
%    - Independent of the Luke export workflow.
%    - Loads the original Single, Simultaneous, and Sequential datasets.
%    - Sim and Seq may be in separate folders or in the same folder.
%    - Detects valid electrode pairs automatically from StimParams/E_MAP.
%    - Excludes bad channels and channel-specific bad trials.
%    - Uses baseline-corrected evoked spike count by default.
%    - Uses all clean trials by default; deterministic balancing is optional.
%    - Analyses two order-matched branches separately:
%         Sim trials from the A->B set versus sequential A->B
%         Sim trials from the B->A set versus sequential B->A
%    - Uses a separate responding-channel union for each matched branch.
%    - Optionally saves one compact E1 MAT result for later cross-dataset
%      analysis; saving can be disabled in USER SETTINGS.
%
%  Required external analysis functions
%    - Depth_s
%    - loadTrig
%
%  Window convention
%    All windows are half-open: [start, end). A spike exactly at the end
%    of a window is not included.
% ========================================================================

clear;

%% ============================= USER SETTINGS ============================

% Original experiment folders.
single_folder = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_Exp1_Single4';
sim_folder    = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_Exp1_Sim4';
seq_folder    = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_Exp1_Seq4';

% If Sim and Seq are stored together, enter exactly the same path for
% sim_folder and seq_folder. The folder will be loaded only once.

analysis_functions_folder = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE';

Electrode_Type = 2;  % 0 rigid, 1 single-shank flex, 2 four-shank flex
FS = 30000;

% Empty or 'all' analyses every automatically detected complete pair.
% A numeric vector selects detected pair numbers printed by this script.
Pairs_To_Analyze = 'all';

% Channel choices:
%   'responding' = a separate union for each matched branch:
%                  Sim(A->B set) union Seq(A->B), and
%                  Sim(B->A set) union Seq(B->A)
%   'all'        = all depth channels returned by Depth_s
%   numeric      = selected depth-channel indices, for example [1 3 5]
Channels_To_Analyze = 'responding';

% Empty uses every amplitude found across the three datasets. Missing
% comparisons remain NaN; amplitudes are never interpolated.
Amplitudes_To_Analyze = [];

% false ignores detected 0 uA/sham entries in all calculations, tables,
% and figures. Set true only when a complete 0 uA comparison is desired.
Include_Zero_Amplitude = false;

baseline_win_ms = [-50 -5];
post_win_ms     = [2 20];
Sequential_PTD_ms = 5;
PTD_Tolerance_ms = 0.01;

% 'all' or 'balanced'. Balanced mode uses the minimum clean-trial count
% across A, B, Sim, A->B, and B->A separately for each channel/amplitude.
Trial_Mode = 'all';
Random_Seed = 1;

% Retained for compatibility with the original version. In this
% order-matched preview, each branch is always plotted separately.
Sequential_Mode = 'separate';

% Broad groups shown in figures and tables:
%   'all', 'simultaneous', 'sequential', or a cell combination such as
%   {'simultaneous','sequential'}.
Curve_Groups_To_Plot = 'all';

% 'individual', 'average', or 'both'. Across-channel averages use only
% channels valid for every displayed curve at each amplitude.
Plot_Figure_Mode = 'average';
Channels_Per_Figure = 12;

Plot_Count_Curves = true;
Plot_Difference_Curves = true;
Show_Error_Bars = true;

% Keep command-window output compact. Detailed channel and late-window
% results remain available in the saved MAT file.
Print_Average_Results = true;
Print_Channel_Results = false;
Print_Late_Window_Check = false;

% ---- Result saving ----
% One MAT file is created for this dataset/run. It contains the underlying
% channel table plus the exact across-channel values used for the
% spike-count-versus-amplitude and difference-versus-amplitude plots.
Save_Results = true;
Results_Save_Folder = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse';
Overwrite_Existing_Result = true;

%% ============================== VALIDATION ===============================

validate_settings(single_folder,sim_folder,seq_folder, ...
    analysis_functions_folder,baseline_win_ms,post_win_ms, ...
    Sequential_PTD_ms,Trial_Mode,Sequential_Mode,Plot_Figure_Mode, ...
    Channels_Per_Figure,FS);

addpath(genpath(analysis_functions_folder));
if exist('Depth_s','file') ~= 2
    error('Linearity1:MissingDepthFunction', ...
        'Depth_s was not found under analysis_functions_folder.');
end
if exist('loadTrig','file') ~= 2
    error('Linearity1:MissingTriggerFunction', ...
        'loadTrig was not found under analysis_functions_folder.');
end

% depth_to_spike_channel = Depth_s(Electrode_Type);

% Depth_s must run inside a dataset folder containing the .rhs file
original_folder = pwd;
restore_folder = onCleanup(@() cd(original_folder));

cd(single_folder);
depth_to_spike_channel = Depth_s(Electrode_Type);

clear restore_folder;
cd(original_folder);

depth_to_spike_channel = double(depth_to_spike_channel(:));
nDepthChannels = numel(depth_to_spike_channel);

fprintf('\nE1 order-matched mean spike-count linearity\n');
fprintf('Dataset: %s\n',folder_leaf_name(sim_folder));
fprintf('Windows: baseline [%g,%g) ms; response [%g,%g) ms; PTD %g ms\n', ...
    baseline_win_ms(1),baseline_win_ms(2), ...
    post_win_ms(1),post_win_ms(2),Sequential_PTD_ms);
fprintf(['Linear references: Sim = A+B unshifted; Seq = first component ' ...
    '+ PTD-shifted second component\n']);

%% ============================== LOAD DATA ================================

Single = load_linearity_dataset(single_folder,'Single',FS,nDepthChannels);
Sim = load_linearity_dataset(sim_folder,'Simultaneous',FS,nDepthChannels);

if same_folder(sim_folder,seq_folder)
    Seq = Sim;
    Seq.role = 'Sequential';
else
    Seq = load_linearity_dataset(seq_folder,'Sequential',FS,nDepthChannels);
end

%% ========================== DETECT VALID PAIRS ===========================

DetectedPairs = detect_complete_pairs(Single,Sim,Seq, ...
    Sequential_PTD_ms,PTD_Tolerance_ms);

if isempty(DetectedPairs)
    error('Linearity1:NoCompletePairs', ...
         ['No analyzable electrode pair was found across Single, Sim, and Seq. ' ...
         'Each pair requires A alone and B alone plus at least one complete ' ...
         'order-matched branch: simultaneous and sequential A->B, or ' ...
         'simultaneous and sequential B->A.']);
end

DetectedPairTable = pair_table(DetectedPairs);
fprintf('\nDetected matched pair/set information\n');
disp(DetectedPairTable);

pair_indices = resolve_pair_selection(Pairs_To_Analyze,numel(DetectedPairs));
[plot_individual,plot_average] = resolve_figure_mode(Plot_Figure_Mode);

%% ============================= MAIN ANALYSIS =============================

condition_codes = {'A','B','SimMatched','Sequential'};
nConditions = numel(condition_codes);
post_duration_ms = diff(post_win_ms);
baseline_duration_ms = diff(baseline_win_ms);
late_win_ms = [post_win_ms(2)-Sequential_PTD_ms, post_win_ms(2)];

all_amplitudes = unique([Single.trial_amp_uA(:); Sim.trial_amp_uA(:); ...
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
    missing_requested = selected_amplitudes(~ismembertol( ...
        selected_amplitudes,all_amplitudes,1e-9));
    if ~isempty(missing_requested)
        warning('Linearity1:RequestedAmplitudesMissing', ...
            'Requested amplitudes not detected and retained as NaN: %s uA.', ...
            num2str(missing_requested));
    end
end
if isempty(selected_amplitudes)
    error('Linearity1:NoSelectedAmplitudes', ...
        ['No amplitudes remain after applying Amplitudes_To_Analyze and ' ...
         'Include_Zero_Amplitude.']);
end

% Start with a genuinely empty structure. A preallocated no-field structure
% cannot accept a populated result and causes "dissimilar structures".
LinearityResults = struct([]);
ResultTables = {};
LateWindowTables = {};

random_stream = RandStream('mt19937ar','Seed',Random_Seed);

for ip = pair_indices
    Pair = DetectedPairs(ip);
    Branches = make_order_matched_branches(Pair);

    for ib = 1:numel(Branches)
        Branch = Branches(ib);
        branch_channels = resolve_channels_for_branch(Channels_To_Analyze, ...
            Branch,Sim,Seq,nDepthChannels,Sequential_PTD_ms,PTD_Tolerance_ms);
        if isempty(branch_channels)
            warning('Linearity1:NoSelectedChannels', ...
                'Pair %d, %s branch has no channels under the requested setting.', ...
                ip,Branch.code);
            continue;
        end

        nSelectedChannels = numel(branch_channels);
        nAmp = numel(selected_amplitudes);
        trial_mean = nan(nSelectedChannels,nAmp,nConditions);
        trial_sem = nan(nSelectedChannels,nAmp,nConditions);
        trial_sd = nan(nSelectedChannels,nAmp,nConditions);
        trial_n = zeros(nSelectedChannels,nAmp,nConditions);
        seq_component_mean = nan(nSelectedChannels,nAmp,2);
        seq_component_sem = nan(nSelectedChannels,nAmp,2);
        seq_component_sd = nan(nSelectedChannels,nAmp,2);
        seq_component_n = zeros(nSelectedChannels,nAmp,2);
        late_mean = nan(nSelectedChannels,nAmp,2);

        for jc = 1:nSelectedChannels
            depth_channel = branch_channels(jc);
            spike_channel = depth_to_spike_channel(depth_channel);
            for ia = 1:nAmp
                amp = selected_amplitudes(ia);
                trial_ids = cell(1,nConditions);
                trial_ids{1} = select_clean_trials(Single,Pair.A,Pair.key, ...
                    'single',amp,0,depth_channel,PTD_Tolerance_ms);
                trial_ids{2} = select_clean_trials(Single,Pair.B,Pair.key, ...
                    'single',amp,0,depth_channel,PTD_Tolerance_ms);
                trial_ids{3} = select_clean_trials(Sim,Branch.sim_order_key,Pair.key, ...
                    'sim_order',amp,0,depth_channel,PTD_Tolerance_ms);
                trial_ids{4} = select_clean_trials(Seq,Branch.seq_order_key,Pair.key, ...
                    'seq',amp,Sequential_PTD_ms,depth_channel,PTD_Tolerance_ms);

                if is_balanced_trial_mode(Trial_Mode)
                    n_balanced = min(cellfun(@numel,trial_ids));
                    if n_balanced > 0
                        for ic = 1:nConditions
                            chosen = randperm(random_stream,numel(trial_ids{ic}),n_balanced);
                            trial_ids{ic} = sort(trial_ids{ic}(chosen));
                        end
                    else
                        trial_ids(:) = {[]};
                    end
                end

                datasets = {Single,Single,Sim,Seq};
                for ic = 1:nConditions
                    values = calculate_evoked_counts(datasets{ic},spike_channel, ...
                        trial_ids{ic},baseline_win_ms,post_win_ms);
                    [trial_mean(jc,ia,ic),trial_sd(jc,ia,ic), ...
                        trial_sem(jc,ia,ic),trial_n(jc,ia,ic)] = ...
                        summarize_values(values);
                end

                % Construct the sequential linear reference from spike
                % times. The first single component remains at time zero;
                % the second single component is shifted later by the PTD
                % before both are counted in the common [2,20) window.
                if strcmp(Branch.code,'AB')
                    seq_single_indices = [1 2]; % A first, B second
                else
                    seq_single_indices = [2 1]; % B first, A second
                end
                seq_delays_ms = [0 Sequential_PTD_ms];
                for iseq = 1:2
                    source_idx = seq_single_indices(iseq);
                    shifted_values = calculate_shifted_evoked_counts( ...
                        Single,spike_channel,trial_ids{source_idx}, ...
                        baseline_win_ms,post_win_ms,seq_delays_ms(iseq));
                    [seq_component_mean(jc,ia,iseq), ...
                        seq_component_sd(jc,ia,iseq), ...
                        seq_component_sem(jc,ia,iseq), ...
                        seq_component_n(jc,ia,iseq)] = ...
                        summarize_values(shifted_values);
                end

                late_A = calculate_evoked_counts(Single,spike_channel, ...
                    trial_ids{1},baseline_win_ms,late_win_ms);
                late_B = calculate_evoked_counts(Single,spike_channel, ...
                    trial_ids{2},baseline_win_ms,late_win_ms);
                late_mean(jc,ia,1) = mean(late_A,'omitnan');
                late_mean(jc,ia,2) = mean(late_B,'omitnan');
            end
        end

        sim_linear_mean = trial_mean(:,:,1)+trial_mean(:,:,2);
        sim_linear_sem = sqrt(trial_sem(:,:,1).^2+trial_sem(:,:,2).^2);
        seq_linear_mean = seq_component_mean(:,:,1)+seq_component_mean(:,:,2);
        seq_linear_sem = sqrt(seq_component_sem(:,:,1).^2+ ...
            seq_component_sem(:,:,2).^2);
        displayed_mean = cat(3,trial_mean(:,:,3),trial_mean(:,:,4));
        displayed_sem = cat(3,trial_sem(:,:,3),trial_sem(:,:,4));
        displayed_n = cat(3,trial_n(:,:,3),trial_n(:,:,4));
        linear_mean_by_display = cat(3,sim_linear_mean,seq_linear_mean);
        linear_sem_by_display = cat(3,sim_linear_sem,seq_linear_sem);
        difference_mean = displayed_mean-linear_mean_by_display;
        difference_sem = sqrt(displayed_sem.^2+linear_sem_by_display.^2);

        P = struct();
        P.detected_pair_index = ip;
        P.electrode_A = Pair.A;
        P.electrode_B = Pair.B;
        P.pair_key = Pair.key;
        P.branch_code = Branch.code;
        P.branch_label = Branch.label;
        P.sim_order_key = Branch.sim_order_key;
        P.seq_order_key = Branch.seq_order_key;
        P.sim_set_indices = Branch.sim_set_indices;
        P.seq_set_indices = Branch.seq_set_indices;
        P.stimulation_depth_A = find_electrode_depth_index(Sim,Pair.A);
        P.stimulation_depth_B = find_electrode_depth_index(Sim,Pair.B);
        P.depth_channels = branch_channels;
        P.spike_channels = depth_to_spike_channel(branch_channels);
        P.amplitudes_uA = selected_amplitudes;
        P.condition_codes = condition_codes;
        P.trial_mean = trial_mean;
        P.trial_sd = trial_sd;
        P.trial_sem = trial_sem;
        P.trial_n = trial_n;
        P.sim_linear_mean = sim_linear_mean;
        P.sim_linear_sem = sim_linear_sem;
        P.seq_linear_mean = seq_linear_mean;
        P.seq_linear_sem = seq_linear_sem;
        P.linear_mean_by_display = linear_mean_by_display;
        P.linear_sem_by_display = linear_sem_by_display;
        % Legacy aliases explicitly refer to the simultaneous reference.
        P.linear_mean = sim_linear_mean;
        P.linear_sem = sim_linear_sem;
        P.seq_component_mean = seq_component_mean;
        P.seq_component_sd = seq_component_sd;
        P.seq_component_sem = seq_component_sem;
        P.seq_component_n = seq_component_n;
        P.seq_component_roles = {'First_single','Second_single_shifted'};
        P.seq_component_delay_ms = [0 Sequential_PTD_ms];
        P.seq_component_source_windows_ms = [post_win_ms; ...
            post_win_ms-Sequential_PTD_ms];
        if strcmp(Branch.code,'AB')
            P.seq_component_electrodes = {Pair.A,Pair.B};
        else
            P.seq_component_electrodes = {Pair.B,Pair.A};
        end
        P.display_condition_codes = {Branch.sim_display_code,Branch.seq_display_code};
        P.displayed_mean = displayed_mean;
        P.displayed_sem = displayed_sem;
        P.displayed_n = displayed_n;
        P.difference_mean = difference_mean;
        P.difference_sem = difference_sem;
        P.late_window_ms = late_win_ms;
        P.late_mean_A_B = late_mean;
        P.settings = struct('baseline_win_ms',baseline_win_ms, ...
            'post_win_ms',post_win_ms,'trial_mode',char(string(Trial_Mode)), ...
            'random_seed',Random_Seed,'sequential_ptd_ms',Sequential_PTD_ms, ...
            'channel_population_mode','order_matched', ...
            'linear_prediction_mode', ...
            'sim_unshifted_seq_second_component_shifted', ...
            'baseline_prediction_mode','each_single_component_own_baseline');

        % Save the exact dataset-average values underlying both average
        % figures. The same common channel mask is used for both linear
        % references and both observed curves at each amplitude.
        [P.spike_count_vs_amplitude,P.difference_vs_amplitude, ...
            P.channels_used_by_amplitude] = build_average_result_tables(P);
        P.n_responding_channels_detected = numel(P.depth_channels);
        P.n_channels_used_by_amplitude = cellfun( ...
            @numel,P.channels_used_by_amplitude);
        P.average_summary = build_compact_average_table(P);

        if isempty(LinearityResults)
            LinearityResults = P;
        else
            LinearityResults(end+1) = P;
        end

        result_index = numel(ResultTables)+1;
        plot_codes = string(P.display_condition_codes);
        plot_labels = string({Branch.sim_plot_label,Branch.seq_plot_label});
        plot_colors = [0 0.35 0.85; Branch.seq_color];
        ResultTables{result_index} = build_result_table(P,plot_codes);
        LateWindowTables{result_index} = build_late_table(P);

        fprintf('\n%s + %s | %s branch\n',Pair.A,Pair.B,Branch.label);
        fprintf('Sets: Sim %s; Sequential %s\n', ...
            num2str(Branch.sim_set_indices),num2str(Branch.seq_set_indices));
        fprintf('Responding channels: %d\n',numel(branch_channels));
        if Print_Average_Results
            disp(P.average_summary);
        end
        if Print_Channel_Results
            disp(ResultTables{result_index});
        end
        if Print_Late_Window_Check
            fprintf('Single-response late-window check [%g, %g) ms\n',late_win_ms);
            disp(LateWindowTables{result_index});
        end

        if Plot_Count_Curves
            if plot_individual
                plot_individual_curves(P,plot_codes,plot_labels,plot_colors, ...
                    'count',Show_Error_Bars,Channels_Per_Figure);
            end
            if plot_average
                plot_average_curves(P,plot_codes,plot_labels,plot_colors, ...
                    'count',Show_Error_Bars);
            end
        end
        if Plot_Difference_Curves
            if plot_individual
                plot_individual_curves(P,plot_codes,plot_labels,plot_colors, ...
                    'difference',Show_Error_Bars,Channels_Per_Figure);
            end
            if plot_average
                plot_average_curves(P,plot_codes,plot_labels,plot_colors, ...
                    'difference',Show_Error_Bars);
            end
        end
    end
end

if isempty(LinearityResults)
    error('Linearity1:NoCompletedResults', ...
        'No order-matched branch produced a complete result.');
end

if Save_Results
    Linearity1Saved = build_saved_result(LinearityResults,ResultTables, ...
        LateWindowTables,DetectedPairTable,Single,Sim,Seq, ...
        single_folder,sim_folder,seq_folder,Electrode_Type,FS, ...
        baseline_win_ms,post_win_ms,Sequential_PTD_ms, ...
        PTD_Tolerance_ms,Trial_Mode,Random_Seed,Include_Zero_Amplitude, ...
        selected_amplitudes,Channels_To_Analyze);

    output_path = save_linearity_result(Linearity1Saved, ...
        Results_Save_Folder,Overwrite_Existing_Result);
    fprintf('\nSaved result: %s\n',output_path);
else
    fprintf('\nAnalysis complete; saving disabled.\n');
end

%% ============================== FUNCTIONS ================================

function validate_settings(single_folder,sim_folder,seq_folder,functions_folder, ...
        baseline_win,post_win,seq_ptd,trial_mode,seq_mode,figure_mode, ...
        channels_per_figure,FS)
folders = {single_folder,sim_folder,seq_folder,functions_folder};
folder_names = {'single_folder','sim_folder','seq_folder', ...
    'analysis_functions_folder'};
for k = 1:numel(folders)
    if ~isfolder(folders{k})
        error('Linearity1:FolderNotFound','%s not found: %s', ...
            folder_names{k},folders{k});
    end
end
validate_window(baseline_win,'baseline_win_ms');
validate_window(post_win,'post_win_ms');
if baseline_win(2) > post_win(1)
    error('Linearity1:OverlappingWindows', ...
        'Baseline and response windows must not overlap.');
end
if ~isscalar(seq_ptd) || ~isfinite(seq_ptd) || seq_ptd <= 0
    error('Linearity1:InvalidPTD','Sequential_PTD_ms must be positive.');
end
if ~any(strcmpi(string(trial_mode),["all","balanced"]))
    error('Linearity1:InvalidTrialMode','Trial_Mode must be all or balanced.');
end
if ~any(strcmpi(string(seq_mode),["separate","merged","both"]))
    error('Linearity1:InvalidSequentialMode', ...
        'Sequential_Mode must be separate, merged, or both.');
end
if ~any(strcmpi(string(figure_mode),["individual","average","both","all"]))
    error('Linearity1:InvalidFigureMode', ...
        'Plot_Figure_Mode must be individual, average, or both.');
end
if ~isscalar(channels_per_figure) || channels_per_figure < 1 || ...
        fix(channels_per_figure) ~= channels_per_figure
    error('Linearity1:InvalidPageSize','Channels_Per_Figure must be a positive integer.');
end
if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('Linearity1:InvalidFS','FS must be positive.');
end
end

function validate_window(window,name)
if ~isnumeric(window) || numel(window) ~= 2 || ...
        any(~isfinite(window)) || window(2) <= window(1)
    error('Linearity1:InvalidWindow','%s must be [start end] with end > start.',name);
end
end

function D = load_linearity_dataset(folder,role,FS,nDepthChannels)
folder = char(string(folder));

[spike_file,spike_variable] = find_spike_source(folder);
spike_loaded = load(spike_file,spike_variable);
sp = spike_loaded.(spike_variable);
if ~iscell(sp)
    error('Linearity1:SpikeVariableNotCell', ...
        '%s in %s is not a cell array.',spike_variable,spike_file);
end

trigger_files = clean_file_list(dir(fullfile(folder,'*.trig.dat')));
if isempty(trigger_files)
    error('Linearity1:MissingTrigger', ...
        'No existing *.trig.dat file was found in %s. No trigger file was generated.',folder);
end
if numel(trigger_files) > 1
    error('Linearity1:AmbiguousTrigger', ...
        'Multiple *.trig.dat files were found in %s. Resolve the ambiguity first.',folder);
end
trig = read_triggers_without_writing(folder);
trig = double(trig(:));

experiment_files = clean_file_list(dir(fullfile(folder,'*_exp_datafile_*.mat')));
if isempty(experiment_files)
    error('Linearity1:MissingExperimentFile', ...
        'No *_exp_datafile_*.mat file was found in %s.',folder);
end
if numel(experiment_files) > 1
    error('Linearity1:AmbiguousExperimentFile', ...
        'Multiple experiment files were found in %s. Resolve the ambiguity first.',folder);
end
experiment_file = fullfile(experiment_files(1).folder,experiment_files(1).name);
E = load(experiment_file,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required = {'StimParams','simultaneous_stim','E_MAP'};
for k = 1:numel(required)
    if ~isfield(E,required{k})
        error('Linearity1:MissingExperimentVariable', ...
            '%s is missing from %s.',required{k},experiment_file);
    end
end

StimParams = E.StimParams;
simN = double(E.simultaneous_stim);
if ~isscalar(simN) || simN < 1 || fix(simN) ~= simN
    error('Linearity1:InvalidSimultaneousStim', ...
        'simultaneous_stim must be a positive integer in %s.',experiment_file);
end
nRows = size(StimParams,1)-1;
if mod(nRows,simN) ~= 0
    error('Linearity1:StimParamRowMismatch', ...
        'StimParams rows are not divisible by simultaneous_stim in %s.',experiment_file);
end
nTrialsCalculated = nRows/simN;
if isfield(E,'n_Trials') && ~isempty(E.n_Trials)
    nTrials = double(E.n_Trials);
    if nTrials ~= nTrialsCalculated
        error('Linearity1:TrialCountMismatch', ...
            'n_Trials disagrees with StimParams in %s.',experiment_file);
    end
else
    nTrials = nTrialsCalculated;
end
if numel(trig) < nTrials
    warning('Linearity1:FewerTriggers', ...
        '%s has %d trials but only %d triggers. Untriggered trials will be excluded.', ...
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
        if strlength(labels(pulse)) == 0
            error('Linearity1:EmptyStimLabel', ...
                'Empty stimulation label in trial %d of %s.',t,experiment_file);
        end
        if ~isempty(emap_labels) && ~any(labels(pulse) == emap_labels)
            error('Linearity1:StimLabelNotInMap', ...
                'Stimulus %s in trial %d was not found in E_MAP for %s.', ...
                labels(pulse),t,experiment_file);
        end
    end
    trial_labels{t} = labels;
    order_key(t) = join(labels,">");
    pair_key(t) = join(sort(labels),"|");
    trial_amp_uA(t) = numeric_cell_value(StimParams{first_row,16}, ...
        'amplitude',t,experiment_file);
    if trial_amp_uA(t) == -1
        trial_amp_uA(t) = 0;
    end
    if simN > 1
        ptd_us = numeric_cell_value(StimParams{first_row+1,6}, ...
            'PTD',t,experiment_file);
        trial_ptd_ms(t) = ptd_us/1000;
    end
end

[set_keys,~,set_index] = unique(order_key,'stable');

[BadCh,bad_channel_file] = load_bad_channels(folder);
[BadTrials,bad_trial_file] = load_bad_trials(folder,nDepthChannels);
[Responding,responding_file] = load_responding(folder);

D = struct();
D.role = role;
D.folder = folder;
D.FS = FS;
D.sp = sp;
D.spike_file = spike_file;
D.spike_variable = spike_variable;
D.trig = trig;
D.trigger_file = fullfile(trigger_files(1).folder,trigger_files(1).name);
D.experiment_file = experiment_file;
D.StimParams = StimParams;
D.simN = simN;
D.E_MAP = E.E_MAP;
D.emap_labels = emap_labels;
D.nTrials = nTrials;
D.trial_labels = trial_labels;
D.order_key = order_key;
D.pair_key = pair_key;
D.trial_amp_uA = trial_amp_uA;
D.trial_ptd_ms = trial_ptd_ms;
D.set_keys = set_keys;
D.set_index = set_index;
D.BadCh = BadCh;
D.BadTrials = BadTrials;
D.Responding = Responding;
D.bad_channel_file = bad_channel_file;
D.bad_trial_file = bad_trial_file;
D.responding_file = responding_file;
end

function [file_path,variable_name] = find_spike_source(folder)
files = clean_file_list(dir(fullfile(folder,'*sp*.mat')));
if isempty(files)
    error('Linearity1:MissingSpikeFile','No spike MAT file was found in %s.',folder);
end

% Ignore known derived count/result files that happen to contain "sp".
keep = true(size(files));
for k = 1:numel(files)
    lower_name = lower(files(k).name);
    if contains(lower_name,'spikecount') || contains(lower_name,'result')
        keep(k) = false;
    end
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
        error('Linearity1:AmbiguousSpikeFile', ...
            'Multiple files in %s contain the priority variable %s:\n%s', ...
            folder,priority{ip},strjoin(cellstr(matches),newline));
    elseif isscalar(matches)
        file_path = char(matches(1));
        variable_name = priority{ip};
        return;
    end
end
error('Linearity1:NoUsableSpikeVariable', ...
    'No sp_corr, sp_SSD, sp_clipped, or sp variable was found in %s.',folder);
end

function trig = read_triggers_without_writing(folder)
original_folder = pwd;
cleanup = onCleanup(@() cd(original_folder));
cd(folder);
trig = loadTrig(0);
clear cleanup;
cd(original_folder);
end

function labels = normalize_emap_labels(E_MAP)
labels = strings(0,1);
if isempty(E_MAP)
    return;
end
values = E_MAP(:);
if numel(values) >= 2
    values = values(2:end);
end
for k = 1:numel(values)
    if iscell(values)
        value = values{k};
    else
        value = values(k);
    end
    label = normalize_label(value);
    if strlength(label) > 0
        labels(end+1,1) = label; %#ok<AGROW>
    end
end
labels = unique(labels,'stable');
end

function label = normalize_label(value)
if iscell(value) && isscalar(value)
    value = value{1};
end
if ischar(value) || (isstring(value) && isscalar(value))
    label = upper(strtrim(string(value)));
elseif isnumeric(value) && isscalar(value) && isfinite(value)
    label = string(value);
else
    error('Linearity1:InvalidStimLabel','A stimulation label could not be converted to text.');
end
label = regexprep(label,'\s+','');
end

function value = numeric_cell_value(raw,name,trial,file)
if iscell(raw) && isscalar(raw)
    raw = raw{1};
end
if isnumeric(raw) && isscalar(raw)
    value = double(raw);
elseif ischar(raw) || (isstring(raw) && isscalar(raw))
    value = str2double(raw);
else
    value = NaN;
end
if ~isfinite(value)
    error('Linearity1:InvalidNumericStimParam', ...
        'Invalid %s value in trial %d of %s.',name,trial,file);
end
end

function [BadCh,file_path] = load_bad_channels(folder)
files = find_named_mat_files(folder,'badchannels');
if isempty(files)
    BadCh = [];
    file_path = '';
    return;
end
if numel(files) > 1
    error('Linearity1:AmbiguousBadChannelFile', ...
        'Multiple bad-channel files were found in %s.',folder);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if isfield(S,'BadCh_perSet')
    BadCh = S.BadCh_perSet;
elseif isfield(S,'BadCh')
    BadCh = S.BadCh;
else
    error('Linearity1:InvalidBadChannelFile', ...
        'No BadCh_perSet or BadCh variable was found in %s.',file_path);
end
end

function [BadTrials,file_path] = load_bad_trials(folder,nDepthChannels)
% Only find the standard *.BadTrials.mat file.
files = find_named_mat_files(folder,'badtrials');
if isempty(files)
    BadTrials = cell(nDepthChannels,1);
    file_path = '';
    return;
end
if numel(files) > 1
    matched_names = strjoin(string({files.name}),newline);

    error('Linearity1:AmbiguousBadTrialFile', ...
        ['Multiple standard BadTrials files were found in:\n%s\n' ...
         'Matched files:\n%s'], ...
        folder,matched_names);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if ~isfield(S,'BadTrials')
    error('Linearity1:InvalidBadTrialFile', ...
        'BadTrials is missing from %s.',file_path);
end
BadTrials = S.BadTrials;
if ~iscell(BadTrials)
    error('Linearity1:InvalidBadTrials', ...
        'BadTrials in %s must be a channel-indexed cell array.',file_path);
end
if numel(BadTrials) < nDepthChannels
    BadTrials(end+1:nDepthChannels) = {[]};
end
BadTrials = BadTrials(:);
end

function [Responding,file_path] = load_responding(folder)
files = find_named_mat_files(folder,'respondingchannels');
if isempty(files)
    Responding = [];
    file_path = '';
    return;
end
if numel(files) > 1
    error('Linearity1:AmbiguousRespondingFile', ...
        'Multiple RespondingChannels files were found in %s.',folder);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if ~isfield(S,'Responding')
    error('Linearity1:InvalidRespondingFile', ...
        'Responding is missing from %s.',file_path);
end
Responding = S.Responding;
end

function files = find_named_mat_files(folder,file_label)
files = clean_file_list(dir(fullfile(folder,'*.mat')));
keep = false(size(files));
% Require:
%   beginning of filename OR ".", "_" or "-"
%   followed by the exact requested label
%   followed immediately by ".mat" at the end
pattern = ['(^|[._-])' ...
           regexptranslate('escape',lower(file_label)) ...
           '\.mat$'];
for k = 1:numel(files)
    keep(k) = ~isempty(regexpi(files(k).name,pattern,'once'));
end
files = files(keep);
end

function files = clean_file_list(files)
if isempty(files)
    return;
end
keep = true(size(files));
for k = 1:numel(files)
    keep(k) = ~startsWith(files(k).name,'._');
end
files = files(keep);
end

function tf = same_folder(a,b)
a = strip_trailing_separator(char(string(a)));
b = strip_trailing_separator(char(string(b)));
tf = strcmp(a,b);
end

function path_out = strip_trailing_separator(path_in)
path_out = path_in;
while numel(path_out) > 1 && any(path_out(end) == ['/' '\'])
    path_out(end) = [];
end
end

function Pairs = detect_complete_pairs(Single,Sim,Seq,target_ptd,tol)
sim_mask = abs(Sim.trial_ptd_ms) < tol & cellfun(@numel,Sim.trial_labels) == 2;
seq_mask = abs(Seq.trial_ptd_ms-target_ptd) < tol & cellfun(@numel,Seq.trial_labels) == 2;
sim_keys = unique(Sim.pair_key(sim_mask),'stable');
seq_keys = unique(Seq.pair_key(seq_mask),'stable');
candidate_keys = intersect(sim_keys,seq_keys,'stable');
single_labels = unique(Single.order_key(cellfun(@numel,Single.trial_labels) == 1));

Pairs = repmat(struct('key','', 'A','', 'B','', 'A_to_B_key','', ...
    'B_to_A_key','', 'sim_set_indices',[], ...
    'sim_a_to_b_set_indices',[], 'sim_b_to_a_set_indices',[], ...
    'a_to_b_set_indices',[], 'b_to_a_set_indices',[]),1,0);

for k = 1:numel(candidate_keys)
    labels = split(candidate_keys(k),'|');
    if numel(labels) ~= 2 || labels(1) == labels(2)
        continue;
    end
    labels = sort(labels);
    A = labels(1);
    B = labels(2);
    AtoB = A+">"+B;
    BtoA = B+">"+A;
    has_singles = any(single_labels == A) && any(single_labels == B);
    has_sim_AtoB = any(sim_mask & Sim.order_key == AtoB);
    has_sim_BtoA = any(sim_mask & Sim.order_key == BtoA);
    has_AtoB = any(seq_mask & Seq.order_key == AtoB);
    has_BtoA = any(seq_mask & Seq.order_key == BtoA);
    % Treat the two sequential orders as independent optional branches.
    % Retain the pair when at least one order has both its nominally matched
    % simultaneous condition and its sequential condition. Never substitute
    % the opposite simultaneous set for a missing matched set.
    has_complete_AtoB = has_sim_AtoB && has_AtoB;
    has_complete_BtoA = has_sim_BtoA && has_BtoA;
    if ~(has_singles && (has_complete_AtoB || has_complete_BtoA))
        continue;
    end
    P = struct();
    P.key = char(candidate_keys(k));
    P.A = char(A);
    P.B = char(B);
    P.A_to_B_key = char(AtoB);
    P.B_to_A_key = char(BtoA);
    P.sim_set_indices = unique(Sim.set_index(sim_mask & Sim.pair_key == candidate_keys(k))).';
    P.sim_a_to_b_set_indices = unique(Sim.set_index( ...
        sim_mask & Sim.order_key == AtoB)).';
    P.sim_b_to_a_set_indices = unique(Sim.set_index( ...
        sim_mask & Sim.order_key == BtoA)).';
    P.a_to_b_set_indices = unique(Seq.set_index(seq_mask & Seq.order_key == AtoB)).';
    P.b_to_a_set_indices = unique(Seq.set_index(seq_mask & Seq.order_key == BtoA)).';
    Pairs(end+1) = P; %#ok<AGROW>
end
end

function T = pair_table(Pairs)
n = numel(Pairs);
Pair = (1:n).';
ElectrodeA = strings(n,1);
ElectrodeB = strings(n,1);
SimAtoBSets = strings(n,1);
SimBtoASets = strings(n,1);
AtoBSets = strings(n,1);
BtoASets = strings(n,1);
for k = 1:n
    ElectrodeA(k) = string(Pairs(k).A);
    ElectrodeB(k) = string(Pairs(k).B);
    SimAtoBSets(k) = vector_text(Pairs(k).sim_a_to_b_set_indices);
    SimBtoASets(k) = vector_text(Pairs(k).sim_b_to_a_set_indices);
    AtoBSets(k) = vector_text(Pairs(k).a_to_b_set_indices);
    BtoASets(k) = vector_text(Pairs(k).b_to_a_set_indices);
end
T = table(Pair,ElectrodeA,ElectrodeB,SimAtoBSets,SimBtoASets, ...
    AtoBSets,BtoASets);
end

function Branches = make_order_matched_branches(Pair)
Branches = repmat(struct(),1,2);

Branches(1).code = 'AB';
Branches(1).label = 'A->B matched';
Branches(1).sim_order_key = Pair.A_to_B_key;
Branches(1).seq_order_key = Pair.A_to_B_key;
Branches(1).sim_set_indices = Pair.sim_a_to_b_set_indices;
Branches(1).seq_set_indices = Pair.a_to_b_set_indices;
Branches(1).sim_display_code = 'Sim_AB_Matched';
Branches(1).seq_display_code = 'A_to_B';
Branches(1).sim_plot_label = 'Sim matched to A->B set';
Branches(1).seq_plot_label = 'A then B';
Branches(1).seq_color = [0.85 0.33 0.10];

Branches(2).code = 'BA';
Branches(2).label = 'B->A matched';
Branches(2).sim_order_key = Pair.B_to_A_key;
Branches(2).seq_order_key = Pair.B_to_A_key;
Branches(2).sim_set_indices = Pair.sim_b_to_a_set_indices;
Branches(2).seq_set_indices = Pair.b_to_a_set_indices;
Branches(2).sim_display_code = 'Sim_BA_Matched';
Branches(2).seq_display_code = 'B_to_A';
Branches(2).sim_plot_label = 'Sim matched to B->A set';
Branches(2).seq_plot_label = 'B then A';
Branches(2).seq_color = [0.55 0.15 0.70];

% A dataset may contain only one sequential order. Analyse only branches
% for which both the order-matched simultaneous set and sequential set are
% present. The absent order remains absent from plots and saved results.
keep_branch = [ ...
    ~isempty(Pair.sim_a_to_b_set_indices) && ~isempty(Pair.a_to_b_set_indices), ...
    ~isempty(Pair.sim_b_to_a_set_indices) && ~isempty(Pair.b_to_a_set_indices)];
Branches = Branches(keep_branch);
end

function out = vector_text(values)
if isempty(values)
    out = "none";
else
    out = string(strtrim(num2str(values)));
end
end

function indices = resolve_pair_selection(request,nPairs)
if ischar(request) || (isstring(request) && isscalar(request))
    if ~strcmpi(strtrim(string(request)),'all')
        error('Linearity1:InvalidPairSelection', ...
            'Pairs_To_Analyze text must be all.');
    end
    indices = 1:nPairs;
elseif isnumeric(request)
    indices = unique(double(request(:).'),'stable');
    if isempty(indices) || any(indices < 1 | indices > nPairs | fix(indices) ~= indices)
        error('Linearity1:InvalidPairSelection', ...
            'Pairs_To_Analyze contains an invalid detected pair number.');
    end
else
    error('Linearity1:InvalidPairSelection', ...
        'Pairs_To_Analyze must be all or a numeric vector.');
end
end

function channels = resolve_channels_for_branch(request,Branch,Sim,Seq,nDepth,target_ptd,tol)
if isnumeric(request)
    channels = unique(double(request(:).'),'stable');
    if isempty(channels) || any(channels < 1 | channels > nDepth | fix(channels) ~= channels)
        error('Linearity1:InvalidChannels', ...
            'Numeric Channels_To_Analyze contains an invalid depth channel.');
    end
    return;
end
mode = lower(strtrim(char(string(request))));
switch mode
    case 'all'
        channels = 1:nDepth;
    case 'responding'
        sim_mask = responding_mask_for_sets( ...
            Sim,Branch.sim_set_indices,0,tol,nDepth);
        seq_mask = responding_mask_for_sets( ...
            Seq,Branch.seq_set_indices,target_ptd,tol,nDepth);
        channels = find(sim_mask | seq_mask).';
        if isempty(channels)
            error('Linearity1:NoRespondingLabels', ...
                ['No responding-channel labels were found for the %s branch. ' ...
                 'Use Channels_To_Analyze=''all'' or verify RespondingChannels files.'], ...
                Branch.label);
        end
    otherwise
        error('Linearity1:InvalidChannels', ...
            'Channels_To_Analyze must be responding, all, or a numeric vector.');
end
end

function mask = responding_mask_for_sets(D,set_indices,target_ptd,tol,nDepth)
mask = false(nDepth,1);
if isempty(D.Responding) || ~isfield(D.Responding,'set')
    return;
end
for set_idx = set_indices
    if set_idx > numel(D.Responding.set) || set_idx < 1
        continue;
    end
    ptd_values = unique(D.trial_ptd_ms(D.set_index == set_idx));
    ptd_position = find(abs(ptd_values-target_ptd) < tol,1);
    if isempty(ptd_position)
        continue;
    end
    amp_struct = D.Responding.set(set_idx).amp;
    for ia = 1:numel(amp_struct)
        if ~isfield(amp_struct(ia),'ptd') || numel(amp_struct(ia).ptd) < ptd_position
            continue;
        end
        ptd_struct = amp_struct(ia).ptd(ptd_position);
        if ~isfield(ptd_struct,'channel')
            continue;
        end
        channel_struct = ptd_struct.channel;
        for ch = 1:min(numel(channel_struct),nDepth)
            if isfield(channel_struct(ch),'is_responsive') && ...
                    logical(channel_struct(ch).is_responsive)
                mask(ch) = true;
            end
        end
    end
end
end

function trial_ids = select_clean_trials(D,ordered_key,pair_key,mode,amp,ptd,depth,tol)
amp_mask = abs(D.trial_amp_uA-amp) < 1e-6;
ptd_mask = abs(D.trial_ptd_ms-ptd) < tol;
switch lower(mode)
    case 'single'
        condition_mask = D.order_key == string(ordered_key);
    case 'sim'
        condition_mask = D.pair_key == string(pair_key);
    case 'sim_order'
        condition_mask = D.pair_key == string(pair_key) & ...
            D.order_key == string(ordered_key);
    case 'seq'
        condition_mask = D.order_key == string(ordered_key);
    otherwise
        error('Linearity1:InternalConditionMode','Unknown internal condition mode.');
end
trial_ids = find(amp_mask & ptd_mask & condition_mask);
if isempty(trial_ids)
    return;
end

keep = trial_ids <= numel(D.trig);
trial_ids = trial_ids(keep);
if isempty(trial_ids)
    return;
end

good_set = true(size(trial_ids));
for k = 1:numel(trial_ids)
    good_set(k) = ~is_bad_channel(D.BadCh,D.set_index(trial_ids(k)),depth);
end
trial_ids = trial_ids(good_set);
if isempty(trial_ids)
    return;
end

if depth <= numel(D.BadTrials) && ~isempty(D.BadTrials{depth})
    bad_trials = double(D.BadTrials{depth}(:));
    trial_ids = trial_ids(~ismember(trial_ids,bad_trials));
end
trial_ids = trial_ids(:).';
end

function tf = is_bad_channel(BadCh,set_index,depth)
tf = false;
if isempty(BadCh)
    return;
end
if iscell(BadCh)
    if set_index <= numel(BadCh) && ~isempty(BadCh{set_index})
        tf = ismember(depth,double(BadCh{set_index}(:)));
    end
elseif isnumeric(BadCh) || islogical(BadCh)
    tf = ismember(depth,double(BadCh(:)));
end
end

function values = calculate_evoked_counts(D,spike_channel,trial_ids,baseline_win,post_win)
if isempty(trial_ids) || spike_channel < 1 || ...
        spike_channel > numel(D.sp) || isempty(D.sp{spike_channel})
    values = nan(0,1);
    return;
end
spike_data = D.sp{spike_channel};
if ~isnumeric(spike_data) || size(spike_data,2) < 1
    values = nan(0,1);
    return;
end
spike_times_ms = double(spike_data(:,1));
baseline_duration_ms = diff(baseline_win);
post_duration_ms = diff(post_win);
values = nan(numel(trial_ids),1);
for k = 1:numel(trial_ids)
    tr = trial_ids(k);
    if tr < 1 || tr > numel(D.trig)
        continue;
    end
    t0_ms = D.trig(tr)/D.FS*1000;
    baseline_count = sum(spike_times_ms >= t0_ms+baseline_win(1) & ...
        spike_times_ms < t0_ms+baseline_win(2));
    post_count = sum(spike_times_ms >= t0_ms+post_win(1) & ...
        spike_times_ms < t0_ms+post_win(2));
    expected_baseline = baseline_count*post_duration_ms/baseline_duration_ms;
    values(k) = post_count-expected_baseline;
end
end

function values = calculate_shifted_evoked_counts(D,spike_channel,trial_ids, ...
        baseline_win,post_win,delay_ms)
% Shift each spike time relative to the single-stimulation trigger, then
% count the shifted spikes that remain inside the common analysis window.
% The baseline estimate is taken from the original single trial and is
% scaled to the duration of the common response window.
if isempty(trial_ids) || spike_channel < 1 || ...
        spike_channel > numel(D.sp) || isempty(D.sp{spike_channel})
    values = nan(0,1);
    return;
end
spike_data = D.sp{spike_channel};
if ~isnumeric(spike_data) || size(spike_data,2) < 1
    values = nan(0,1);
    return;
end
if ~isscalar(delay_ms) || ~isfinite(delay_ms) || delay_ms < 0
    error('Linearity1:InvalidShiftDelay', ...
        'The sequential component delay must be a finite nonnegative scalar.');
end

spike_times_ms = double(spike_data(:,1));
baseline_duration_ms = diff(baseline_win);
post_duration_ms = diff(post_win);
values = nan(numel(trial_ids),1);
for k = 1:numel(trial_ids)
    tr = trial_ids(k);
    if tr < 1 || tr > numel(D.trig)
        continue;
    end
    t0_ms = D.trig(tr)/D.FS*1000;
    relative_times_ms = spike_times_ms-t0_ms;
    baseline_count = sum(relative_times_ms >= baseline_win(1) & ...
        relative_times_ms < baseline_win(2));
    shifted_times_ms = relative_times_ms+delay_ms;
    post_count = sum(shifted_times_ms >= post_win(1) & ...
        shifted_times_ms < post_win(2));
    expected_baseline = baseline_count*post_duration_ms/baseline_duration_ms;
    values(k) = post_count-expected_baseline;
end
end

function [mean_value,sd_value,sem_value,n] = summarize_values(values)
values = values(isfinite(values));
n = numel(values);
if n == 0
    mean_value = NaN;
    sd_value = NaN;
    sem_value = NaN;
elseif n == 1
    mean_value = values(1);
    sd_value = NaN;
    sem_value = NaN;
else
    mean_value = mean(values);
    sd_value = std(values,0);
    sem_value = sd_value/sqrt(n);
end
end

function depth = find_electrode_depth_index(D,label)
depth = NaN;
target = normalize_label(label);
raw = D.E_MAP(:);
if numel(raw) < 2
    return;
end
for k = 2:numel(raw)
    if iscell(raw)
        value = raw{k};
    else
        value = raw(k);
    end
    try
        current = normalize_label(value);
    catch
        continue;
    end
    if current == target
        depth = k-1;
        return;
    end
end
end

function [codes,labels,colors] = resolve_display_conditions(groups_request,seq_mode)
groups = normalize_group_request(groups_request);
include_sim = any(groups == "all" | groups == "simultaneous" | groups == "sim");
include_seq = any(groups == "all" | groups == "sequential" | groups == "seq");
if ~include_sim && ~include_seq
    error('Linearity1:EmptyCurveGroups','No recognized curve group was selected.');
end
codes = strings(0,1);
labels = strings(0,1);
colors = zeros(0,3);
if include_sim
    codes(end+1,1) = "Sim";
    labels(end+1,1) = "A+B simultaneous";
    colors(end+1,:) = [0.00 0.35 0.80];
end
if include_seq
    switch lower(char(string(seq_mode)))
        case 'separate'
            codes = [codes; "A_to_B"; "B_to_A"];
            labels = [labels; "A then B"; "B then A"];
            colors = [colors; 0.85 0.33 0.10; 0.55 0.20 0.65];
        case 'merged'
            codes(end+1,1) = "SeqMerged";
            labels(end+1,1) = "Sequential order average";
            colors(end+1,:) = [0.20 0.60 0.25];
        case 'both'
            codes = [codes; "A_to_B"; "B_to_A"; "SeqMerged"];
            labels = [labels; "A then B"; "B then A"; ...
                "Sequential order average"];
            colors = [colors; 0.85 0.33 0.10; 0.55 0.20 0.65; ...
                0.20 0.60 0.25];
    end
end
end

function groups = normalize_group_request(request)
if ischar(request) || (isstring(request) && isscalar(request))
    groups = lower(strtrim(string(request)));
elseif isstring(request)
    groups = lower(strtrim(request(:)));
elseif iscell(request)
    groups = lower(strtrim(string(request(:))));
else
    error('Linearity1:InvalidCurveGroups', ...
        'Curve_Groups_To_Plot must be text or a cell array of text.');
end
valid = ["all","simultaneous","sim","sequential","seq"];
if isempty(groups) || any(~ismember(groups,valid))
    error('Linearity1:InvalidCurveGroups', ...
        'Use all, simultaneous, sequential, or a combination.');
end
end

function [individual,average] = resolve_figure_mode(request)
switch lower(strtrim(char(string(request))))
    case 'individual'
        individual = true;
        average = false;
    case 'average'
        individual = false;
        average = true;
    case {'both','all'}
        individual = true;
        average = true;
    otherwise
        error('Linearity1:InvalidFigureMode','Invalid Plot_Figure_Mode.');
end
end

function tf = is_balanced_trial_mode(request)
tf = strcmpi(strtrim(char(string(request))),'balanced');
end

function [SpikeTable,DifferenceTable,channels_used] = ...
        build_average_result_tables(P)
nAmp = numel(P.amplitudes_uA);
nDisplayed = numel(P.display_condition_codes);
channels_used = cell(1,nAmp);

% Spike-count table: both singles, the condition-specific linear
% references, matched Sim, and matched Seq.
nSpikeRows = nAmp*(4+nDisplayed);
Pair = repmat(string(P.electrode_A)+" + "+string(P.electrode_B),nSpikeRows,1);
Branch = repmat(string(P.branch_code),nSpikeRows,1);
Amplitude_uA = nan(nSpikeRows,1);
Condition = strings(nSpikeRows,1);
NRespondingChannels = repmat(numel(P.depth_channels),nSpikeRows,1);
NChannelsUsed = zeros(nSpikeRows,1);
MeanSpikeCountPerTrial = nan(nSpikeRows,1);
SEM_AcrossChannels = nan(nSpikeRows,1);

nDifferenceRows = nAmp*nDisplayed;
DPair = repmat(string(P.electrode_A)+" + "+string(P.electrode_B),nDifferenceRows,1);
DBranch = repmat(string(P.branch_code),nDifferenceRows,1);
DAmplitude_uA = nan(nDifferenceRows,1);
DCondition = strings(nDifferenceRows,1);
DNRespondingChannels = repmat(numel(P.depth_channels),nDifferenceRows,1);
DNChannelsUsed = zeros(nDifferenceRows,1);
MeanDifference = nan(nDifferenceRows,1);
DSEM_AcrossChannels = nan(nDifferenceRows,1);

spike_row = 0;
difference_row = 0;
for ia = 1:nAmp
    valid = isfinite(P.trial_mean(:,ia,1)) & ...
        isfinite(P.trial_mean(:,ia,2)) & ...
        isfinite(P.sim_linear_mean(:,ia)) & ...
        isfinite(P.seq_linear_mean(:,ia));
    for id = 1:nDisplayed
        valid = valid & isfinite(P.displayed_mean(:,ia,id));
    end
    channels_used{ia} = P.depth_channels(valid);

    spike_codes = ["Single_A","Single_B","SimLinear","SeqLinear", ...
        string(P.display_condition_codes)];
    spike_values = cell(1,4+nDisplayed);
    spike_values{1} = P.trial_mean(:,ia,1);
    spike_values{2} = P.trial_mean(:,ia,2);
    spike_values{3} = P.sim_linear_mean(:,ia);
    spike_values{4} = P.seq_linear_mean(:,ia);
    for id = 1:nDisplayed
        spike_values{4+id} = P.displayed_mean(:,ia,id);
    end

    for ic = 1:numel(spike_codes)
        spike_row = spike_row+1;
        x = spike_values{ic}(valid);
        [avg,se,n] = summarize_channel_vector(x);
        Amplitude_uA(spike_row) = P.amplitudes_uA(ia);
        Condition(spike_row) = spike_codes(ic);
        NChannelsUsed(spike_row) = n;
        MeanSpikeCountPerTrial(spike_row) = avg;
        SEM_AcrossChannels(spike_row) = se;
    end

    for id = 1:nDisplayed
        difference_row = difference_row+1;
        x = P.difference_mean(valid,ia,id);
        [avg,se,n] = summarize_channel_vector(x);
        DAmplitude_uA(difference_row) = P.amplitudes_uA(ia);
        DCondition(difference_row) = string(P.display_condition_codes{id});
        DNChannelsUsed(difference_row) = n;
        MeanDifference(difference_row) = avg;
        DSEM_AcrossChannels(difference_row) = se;
    end
end

SpikeTable = table(Pair,Branch,Amplitude_uA,Condition, ...
    NRespondingChannels,NChannelsUsed,MeanSpikeCountPerTrial, ...
    SEM_AcrossChannels);
DifferenceTable = table(DPair,DBranch,DAmplitude_uA,DCondition, ...
    DNRespondingChannels,DNChannelsUsed,MeanDifference, ...
    DSEM_AcrossChannels, ...
    'VariableNames',{'Pair','Branch','Amplitude_uA','Condition', ...
    'NRespondingChannels','NChannelsUsed','MeanDifference', ...
    'SEM_AcrossChannels'});
end

function T = build_compact_average_table(P)
Amplitude_uA = P.amplitudes_uA(:);
nAmp = numel(Amplitude_uA);
NResponding = repmat(numel(P.depth_channels),nAmp,1);
NUsed = P.n_channels_used_by_amplitude(:);
SimLinear = nan(nAmp,1);
SeqLinear = nan(nAmp,1);
SimMatched = nan(nAmp,1);
Sequential = nan(nAmp,1);
SimMinusLinear = nan(nAmp,1);
SeqMinusLinear = nan(nAmp,1);

S = P.spike_count_vs_amplitude;
D = P.difference_vs_amplitude;
sim_code = string(P.display_condition_codes{1});
seq_code = string(P.display_condition_codes{2});
for ia = 1:nAmp
    s_amp = abs(S.Amplitude_uA-Amplitude_uA(ia)) < 1e-9;
    d_amp = abs(D.Amplitude_uA-Amplitude_uA(ia)) < 1e-9;
    SimLinear(ia) = compact_table_value(S,s_amp,"SimLinear", ...
        'MeanSpikeCountPerTrial');
    SeqLinear(ia) = compact_table_value(S,s_amp,"SeqLinear", ...
        'MeanSpikeCountPerTrial');
    SimMatched(ia) = compact_table_value(S,s_amp,sim_code, ...
        'MeanSpikeCountPerTrial');
    Sequential(ia) = compact_table_value(S,s_amp,seq_code, ...
        'MeanSpikeCountPerTrial');
    SimMinusLinear(ia) = compact_table_value(D,d_amp,sim_code, ...
        'MeanDifference');
    SeqMinusLinear(ia) = compact_table_value(D,d_amp,seq_code, ...
        'MeanDifference');
end

T = table(Amplitude_uA,NResponding,NUsed,SimLinear,SeqLinear, ...
    SimMatched,Sequential,SimMinusLinear,SeqMinusLinear);
end

function value = compact_table_value(T,amp_mask,condition,variable_name)
idx = find(amp_mask & T.Condition == condition,1);
if isempty(idx)
    value = NaN;
else
    value = T.(variable_name)(idx);
end
end

function [avg,se,n] = summarize_channel_vector(x)
x = double(x(:));
x = x(isfinite(x));
n = numel(x);
if n == 0
    avg = NaN;
    se = NaN;
    return;
end
avg = mean(x);
if n > 1
    se = std(x,0)/sqrt(n);
else
    se = NaN;
end
end

function Saved = build_saved_result(LinearityResults,ResultTables, ...
        LateWindowTables,DetectedPairTable,Single,Sim,Seq, ...
        single_folder,sim_folder,seq_folder,Electrode_Type,FS, ...
        baseline_win,post_win,seq_ptd,ptd_tol,trial_mode,random_seed, ...
        include_zero,amplitudes,channel_request)

Saved = struct();
Saved.format_version = '2.0';
Saved.analysis_name = 'Linearity1_MeanSpikeCount_ConditionSpecificLinear';
Saved.created_at = datestr(now,31);
Saved.animal_id = detect_animal_id(sim_folder);
Saved.dataset_id = folder_leaf_name(sim_folder);
Saved.folders = struct('single',char(string(single_folder)), ...
    'simultaneous',char(string(sim_folder)), ...
    'sequential',char(string(seq_folder)));
Saved.settings = struct('electrode_type',Electrode_Type,'FS',FS, ...
    'baseline_win_ms',baseline_win,'post_win_ms',post_win, ...
    'sequential_ptd_ms',seq_ptd,'ptd_tolerance_ms',ptd_tol, ...
    'trial_mode',char(string(trial_mode)),'random_seed',random_seed, ...
    'include_zero_amplitude',logical(include_zero), ...
    'amplitudes_uA',amplitudes, ...
    'channels_to_analyze',channel_request, ...
    'channel_population_mode','order_matched', ...
    'linear_prediction_mode','sim_unshifted_seq_second_component_shifted', ...
    'baseline_prediction_mode','each_single_component_own_baseline');
Saved.detected_pairs = DetectedPairTable;

Saved.source_files = struct();
Saved.source_files.single = source_file_metadata(Single);
Saved.source_files.simultaneous = source_file_metadata(Sim);
Saved.source_files.sequential = source_file_metadata(Seq);

Saved.branch_results = LinearityResults;
Saved.responding_channel_summary = ...
    build_responding_channel_summary(LinearityResults);
spike_tables = arrayfun(@(x) x.spike_count_vs_amplitude, ...
    LinearityResults,'UniformOutput',false);
difference_tables = arrayfun(@(x) x.difference_vs_amplitude, ...
    LinearityResults,'UniformOutput',false);
Saved.spike_count_vs_amplitude = vertcat(spike_tables{:});
Saved.difference_vs_amplitude = vertcat(difference_tables{:});
Saved.channel_results = vertcat(ResultTables{:});
Saved.late_window_results = vertcat(LateWindowTables{:});
end

function T = build_responding_channel_summary(Results)
pair_indices = unique([Results.detected_pair_index],'stable');
nRows = 3*numel(pair_indices);
DetectedPairIndex = nan(nRows,1);
Pair = strings(nRows,1);
Population = strings(nRows,1);
NRespondingChannels = zeros(nRows,1);
RespondingChannels = cell(nRows,1);
row = 0;
for ip = pair_indices
    pair_results = Results([Results.detected_pair_index] == ip);
    ab = pair_results(strcmp({pair_results.branch_code},'AB'));
    ba = pair_results(strcmp({pair_results.branch_code},'BA'));

    % Save every available branch independently. Older datasets may have
    % only one sequential order, so do not require both AB and BA here.
    available_channel_sets = {};
    available_names = strings(0,1);
    if ~isempty(ab)
        available_channel_sets{end+1} = ab(1).depth_channels; %#ok<AGROW>
        available_names(end+1,1) = "AB"; %#ok<AGROW>
    end
    if ~isempty(ba)
        available_channel_sets{end+1} = ba(1).depth_channels; %#ok<AGROW>
        available_names(end+1,1) = "BA"; %#ok<AGROW>
    end
    if isempty(available_channel_sets)
        continue;
    end

    across_union = available_channel_sets{1};
    for k = 2:numel(available_channel_sets)
        across_union = union(across_union,available_channel_sets{k},'stable');
    end
    channel_sets = [available_channel_sets,{across_union}];
    population_names = [available_names;"AcrossBranchesUnion"];

    for k = 1:numel(channel_sets)
        row = row+1;
        DetectedPairIndex(row) = ip;
        Pair(row) = string(pair_results(1).electrode_A)+" + "+ ...
            string(pair_results(1).electrode_B);
        Population(row) = population_names(k);
        RespondingChannels{row} = channel_sets{k};
        NRespondingChannels(row) = numel(channel_sets{k});
    end
end
keep = 1:row;
T = table(DetectedPairIndex(keep),Pair(keep),Population(keep), ...
    NRespondingChannels(keep),RespondingChannels(keep), ...
    'VariableNames',{'DetectedPairIndex','Pair','Population', ...
    'NRespondingChannels','RespondingChannels'});
end

function M = source_file_metadata(D)
M = struct('folder',D.folder,'spike_file',D.spike_file, ...
    'spike_variable',D.spike_variable, ...
    'trigger_file',D.trigger_file, ...
    'experiment_file',D.experiment_file, ...
    'bad_channel_file',D.bad_channel_file, ...
    'bad_trial_file',D.bad_trial_file, ...
    'responding_file',D.responding_file);
end

function output_path = save_linearity_result(Linearity1Saved,save_folder,overwrite)
save_folder = char(string(save_folder));
if isempty(strtrim(save_folder))
    error('Linearity1:EmptySaveFolder', ...
        'Results_Save_Folder cannot be empty when Save_Results is true.');
end
if ~isfolder(save_folder)
    [ok,message] = mkdir(save_folder);
    if ~ok
        error('Linearity1:CannotCreateSaveFolder', ...
            'Could not create %s: %s',save_folder,message);
    end
end

animal = sanitize_filename(Linearity1Saved.animal_id);
dataset = sanitize_filename(Linearity1Saved.dataset_id);
filename = sprintf('%s_%s_Linearity1_ShiftedSeq_v2.mat',animal,dataset);
output_path = fullfile(save_folder,filename);
if isfile(output_path) && ~overwrite
    error('Linearity1:ResultAlreadyExists', ...
        ['Result already exists:\n%s\nSet Overwrite_Existing_Result=true ' ...
         'to replace it.'],output_path);
end
save(output_path,'Linearity1Saved','-v7.3');
end

function animal_id = detect_animal_id(folder)
animal_id = regexp(upper(char(string(folder))),'DX\d+', ...
    'match','once');
if isempty(animal_id)
    animal_id = 'AnimalUnknown';
end
end

function name = folder_leaf_name(folder)
[~,name] = fileparts(strip_trailing_separator(char(string(folder))));
if isempty(name)
    name = 'DatasetUnknown';
end
end

function out = sanitize_filename(value)
out = regexprep(char(string(value)),'[^A-Za-z0-9_-]','_');
if isempty(out)
    out = 'Unknown';
end
end

function T = build_result_table(P,plot_codes)
nRows = numel(P.depth_channels)*numel(P.amplitudes_uA)*numel(plot_codes);
Pair = repmat(string(P.electrode_A)+" + "+string(P.electrode_B),nRows,1);
DepthChannel = nan(nRows,1);
Amplitude_uA = nan(nRows,1);
Condition = strings(nRows,1);
NTrials = nan(nRows,1);
Observed = nan(nRows,1);
Linear = nan(nRows,1);
Difference = nan(nRows,1);
row = 0;
for jc = 1:numel(P.depth_channels)
    for ia = 1:numel(P.amplitudes_uA)
        for ic = 1:numel(plot_codes)
            row = row+1;
            idx = display_code_index(P,plot_codes(ic));
            DepthChannel(row) = P.depth_channels(jc);
            Amplitude_uA(row) = P.amplitudes_uA(ia);
            Condition(row) = plot_codes(ic);
            NTrials(row) = P.displayed_n(jc,ia,idx);
            Observed(row) = P.displayed_mean(jc,ia,idx);
            Linear(row) = P.linear_mean_by_display(jc,ia,idx);
            Difference(row) = P.difference_mean(jc,ia,idx);
        end
    end
end
T = table(Pair,DepthChannel,Amplitude_uA,Condition,NTrials, ...
    Observed,Linear,Difference);
end

function T = build_late_table(P)
nRows = numel(P.depth_channels)*numel(P.amplitudes_uA);
Pair = repmat(string(P.electrode_A)+" + "+string(P.electrode_B),nRows,1);
Branch = repmat(string(P.branch_code),nRows,1);
DepthChannel = nan(nRows,1);
Amplitude_uA = nan(nRows,1);
A_Late = nan(nRows,1);
B_Late = nan(nRows,1);
row = 0;
for jc = 1:numel(P.depth_channels)
    for ia = 1:numel(P.amplitudes_uA)
        row = row+1;
        DepthChannel(row) = P.depth_channels(jc);
        Amplitude_uA(row) = P.amplitudes_uA(ia);
        A_Late(row) = P.late_mean_A_B(jc,ia,1);
        B_Late(row) = P.late_mean_A_B(jc,ia,2);
    end
end
T = table(Pair,Branch,DepthChannel,Amplitude_uA,A_Late,B_Late);
end

function idx = display_code_index(P,code)
idx = find(string(P.display_condition_codes) == string(code),1);
if isempty(idx)
    error('Linearity1:InternalDisplayCode','Unknown displayed condition code: %s',code);
end
end

function plot_individual_curves(P,codes,labels,colors,mode,show_errors,page_size)
nChannels = numel(P.depth_channels);
nPages = ceil(nChannels/page_size);
for page = 1:nPages
    first = (page-1)*page_size+1;
    last = min(page*page_size,nChannels);
    selected = first:last;
    nTiles = numel(selected);
    nColumns = ceil(sqrt(nTiles));
    nRows = ceil(nTiles/nColumns);
    figure('Color','w','Position',[50 50 1500 850]);
    tiledlayout(nRows,nColumns,'TileSpacing','compact','Padding','compact');
    if strcmp(mode,'count')
        figure_title = sprintf('%s + %s: spike count versus amplitude', ...
            P.electrode_A,P.electrode_B);
    else
        figure_title = sprintf('%s + %s: observed minus linear reference', ...
            P.electrode_A,P.electrode_B);
    end
    figure_title = sprintf('%s | %s',figure_title,P.branch_label);
    if nPages > 1
        figure_title = sprintf('%s (page %d of %d)',figure_title,page,nPages);
    end
    sgtitle(figure_title,'FontWeight','bold');
    for jc = selected
        ax = nexttile;
        hold(ax,'on');
        if strcmp(mode,'count')
            plot_one_curve(ax,P.amplitudes_uA,P.sim_linear_mean(jc,:), ...
                P.sim_linear_sem(jc,:),[0 0 0],'--o', ...
                'Sim linear (A+B)',show_errors);
            plot_one_curve(ax,P.amplitudes_uA,P.seq_linear_mean(jc,:), ...
                P.seq_linear_sem(jc,:),[0.4 0.4 0.4],'--s', ...
                'Seq linear (shifted)',show_errors);
        else
            yline(ax,0,'--','Color',[0.25 0.25 0.25], ...
                'HandleVisibility','off');
        end
        for k = 1:numel(codes)
            idx = display_code_index(P,codes(k));
            if strcmp(mode,'count')
                y = P.displayed_mean(jc,:,idx);
                se = P.displayed_sem(jc,:,idx);
            else
                y = P.difference_mean(jc,:,idx);
                se = P.difference_sem(jc,:,idx);
            end
            plot_one_curve(ax,P.amplitudes_uA,y,se,colors(k,:),'-o', ...
                char(labels(k)),show_errors);
        end
        title_text = sprintf('Depth channel %d',P.depth_channels(jc));
        stim_tags = strings(0,1);
        if P.depth_channels(jc) == P.stimulation_depth_A
            stim_tags(end+1) = "Stim A"; %#ok<AGROW>
        end
        if P.depth_channels(jc) == P.stimulation_depth_B
            stim_tags(end+1) = "Stim B"; %#ok<AGROW>
        end
        if ~isempty(stim_tags)
            title_text = sprintf('%s (%s)',title_text,strjoin(cellstr(stim_tags),', '));
        end
        title(ax,title_text,'Interpreter','none');
        xlabel(ax,'Amplitude (uA)');
        if strcmp(mode,'count')
            ylabel(ax,'Evoked spikes/trial');
        else
            ylabel(ax,'Observed - linear');
        end
        box(ax,'off');
        if jc == selected(1)
            legend(ax,'Location','best','Box','off');
        end
    end
end
end

function plot_average_curves(P,codes,labels,colors,mode,show_errors)
figure('Color','w','Position',[100 100 760 540]);
ax = axes();
hold(ax,'on');
nAmp = numel(P.amplitudes_uA);
nCodes = numel(codes);
common_valid = false(numel(P.depth_channels),nAmp);
for ia = 1:nAmp
    valid = isfinite(P.sim_linear_mean(:,ia)) & ...
        isfinite(P.seq_linear_mean(:,ia));
    for k = 1:nCodes
        idx = display_code_index(P,codes(k));
        valid = valid & isfinite(P.displayed_mean(:,ia,idx));
    end
    common_valid(:,ia) = valid;
end

if strcmp(mode,'count')
    [sim_linear_avg,sim_linear_sem] = summarize_across_channels( ...
        P.sim_linear_mean,common_valid);
    [seq_linear_avg,seq_linear_sem] = summarize_across_channels( ...
        P.seq_linear_mean,common_valid);
    plot_one_curve(ax,P.amplitudes_uA,sim_linear_avg,sim_linear_sem, ...
        [0 0 0],'--o','Sim linear (A+B)',show_errors);
    plot_one_curve(ax,P.amplitudes_uA,seq_linear_avg,seq_linear_sem, ...
        [0.4 0.4 0.4],'--s','Seq linear (shifted)',show_errors);
else
    yline(ax,0,'--','Color',[0.25 0.25 0.25],'HandleVisibility','off');
end

for k = 1:nCodes
    idx = display_code_index(P,codes(k));
    if strcmp(mode,'count')
        values = P.displayed_mean(:,:,idx);
    else
        values = P.difference_mean(:,:,idx);
    end
    [avg,se] = summarize_across_channels(values,common_valid);
    plot_one_curve(ax,P.amplitudes_uA,avg,se,colors(k,:),'-o', ...
        char(labels(k)),show_errors);
end

if strcmp(mode,'count')
    title_text = sprintf('%s + %s: selected-channel average', ...
        P.electrode_A,P.electrode_B);
    ylabel_text = 'Mean evoked spikes/trial';
else
    title_text = sprintf('%s + %s: average nonlinear difference', ...
        P.electrode_A,P.electrode_B);
    ylabel_text = 'Mean observed - linear';
end
title_text = sprintf('%s | %s',title_text,P.branch_label);
title(ax,title_text,'FontWeight','bold','Interpreter','none');
xlabel(ax,'Amplitude (uA)');
ylabel(ax,ylabel_text);
legend(ax,'Location','best','Box','off');
box(ax,'off');

n_common = sum(common_valid,1);
for ia = 1:nAmp
    if n_common(ia) > 0
        yl = ylim(ax);
        text(ax,P.amplitudes_uA(ia),yl(1)+0.03*(yl(2)-yl(1)), ...
            sprintf('n=%d',n_common(ia)),'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontSize',8,'Color',[0.3 0.3 0.3]);
    end
end
end

function [avg,se] = summarize_across_channels(values,mask)
nAmp = size(values,2);
avg = nan(1,nAmp);
se = nan(1,nAmp);
for ia = 1:nAmp
    x = values(mask(:,ia),ia);
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end
    avg(ia) = mean(x);
    if numel(x) > 1
        se(ia) = std(x,0)/sqrt(numel(x));
    end
end
end

function plot_one_curve(ax,x,y,se,color,line_style,label,show_errors)
x = double(x(:).');
y = double(y(:).');
se = double(se(:).');
valid = isfinite(x) & isfinite(y);
if ~any(valid)
    return;
end
if show_errors
    se_plot = se;
    se_plot(~isfinite(se_plot)) = 0;
    errorbar(ax,x(valid),y(valid),se_plot(valid),line_style, ...
        'Color',color,'MarkerFaceColor',color,'LineWidth',1.5, ...
        'MarkerSize',5,'CapSize',6,'DisplayName',label);
else
    plot(ax,x(valid),y(valid),line_style,'Color',color, ...
        'MarkerFaceColor',color,'LineWidth',1.5,'MarkerSize',5, ...
        'DisplayName',label);
end
end
