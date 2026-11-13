%% ========================================================================
%  LINEARITY ANALYSIS 3: SHORT-LATENCY PSTH TEMPORAL LINEARITY
%
%  Purpose
%    Compare observed simultaneous and sequential dual-electrode PSTHs
%    with predictions constructed from the two single-electrode PSTHs.
%
%       Sim prediction       = A(t) + B(t)
%       A->B Seq prediction  = A(t) + B(t-PTD)
%       B->A Seq prediction  = B(t) + A(t-PTD)
%
%  Main properties
%    - Independent of the Luke workflow and of E1/E2 execution.
%    - Sim and Seq may be stored in separate folders or the same folder.
%    - Detects electrode labels and analyzable pairs automatically.
%    - A pair may contain A->B only, B->A only, or both orders.
%    - Uses a separate responding-channel union for each available branch.
%    - Excludes standard bad channels and channel-specific bad trials.
%    - Ignores MultiISI QC/Responding files by exact filename matching.
%    - Uses sp_corr first; then sp_SSD, sp_clipped, or sp.
%    - Reads existing triggers only and never creates or modifies them.
%    - Uses baseline-corrected PSTHs and half-open windows [start,end).
%    - Symmetric Gaussian smoothing is used for plotting only.
%    - Numerical residuals use unsmoothed, baseline-corrected spike counts.
%    - Does not save results in this exploratory version.
%
%  Required external functions
%    - Depth_s
%    - loadTrig
% ========================================================================

clear;

%% ============================= USER SETTINGS ============================

single_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Single4_new';
sim_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim4';
seq_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim4';

% Enter the same path for sim_folder and seq_folder when both conditions
% are stored together. That folder is loaded only once.
analysis_functions_folder = ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE';

Electrode_Type = 2; % 0 rigid; 1 single-shank flex; 2 four-shank flex
FS = 30000;

% Empty or 'all' analyses every detected pair. A numeric vector selects
% pair numbers from the detected-pair table printed by this script.
Pairs_To_Analyze = 'all';

% 'responding' uses a separate Sim-union-Seq population for each branch.
% 'all' uses all depth channels. A numeric vector selects depth indices.
Channels_To_Analyze = 'responding';

% Empty uses all detected amplitudes. Missing condition combinations are
% retained as unavailable rather than interpolated.
Amplitudes_To_Analyze = [];
Include_Zero_Amplitude = false;

baseline_win_ms = [-50 -5];
response_win_ms = [2 20];

% The response summary is divided into early and late components. These
% two half-open windows must exactly partition response_win_ms.
early_win_ms = [2 7];
late_win_ms  = [7 20];

Sequential_PTD_ms = 5;
PTD_Tolerance_ms = 0.01;

% PSTH settings. The compute window is automatically extended beyond the
% display window to avoid smoothing and shifting edge effects.
display_win_ms = [-10 30];
bin_ms = 1;
smooth_sigma_ms = 1.5;

% 'all' uses all clean trials. 'balanced' deterministically selects the
% same minimum number among A, B, matched Sim, and matched Seq separately
% for every channel/amplitude.
Trial_Mode = 'all';
Random_Seed = 1;

% One page contains this many amplitude rows. Each row has prediction
% components, observed-versus-predicted PSTHs, and residual PSTHs.
Amplitudes_Per_Figure = 4;
Show_SEM = true;

% Optional diagnostic figures show clean-trial rasters and PSTHs for only
% the requested depth channels. Keep false for the population analysis.
Plot_Channel_Diagnostics = false;
Diagnostic_Channels = [];
Diagnostic_Amplitudes = [];

%% ============================== VALIDATION ===============================

validate_settings(single_folder,sim_folder,seq_folder, ...
    analysis_functions_folder,baseline_win_ms,response_win_ms, ...
    early_win_ms,late_win_ms,display_win_ms,Sequential_PTD_ms, ...
    PTD_Tolerance_ms,bin_ms,smooth_sigma_ms,Trial_Mode, ...
    Amplitudes_Per_Figure,FS);

addpath(genpath(analysis_functions_folder));
if exist('Depth_s','file') ~= 2
    error('LinearityPSTH:MissingDepthFunction', ...
        'Depth_s was not found under analysis_functions_folder.');
end
if exist('loadTrig','file') ~= 2
    error('LinearityPSTH:MissingTriggerFunction', ...
        'loadTrig was not found under analysis_functions_folder.');
end

% Depth_s reads the Intan header in the current dataset folder.
original_folder = pwd;
restore_folder = onCleanup(@() cd(original_folder));
cd(single_folder);
depth_to_spike_channel = Depth_s(Electrode_Type);
clear restore_folder;
cd(original_folder);

depth_to_spike_channel = double(depth_to_spike_channel(:));
nDepthChannels = numel(depth_to_spike_channel);

fprintf('\nPSTH temporal linearity analysis\n');
fprintf('Dataset: %s\n',folder_leaf_name(sim_folder));
fprintf(['Windows: baseline [%g,%g) ms; response [%g,%g) ms; ' ...
    'PTD %g ms\n'],baseline_win_ms(1),baseline_win_ms(2), ...
    response_win_ms(1),response_win_ms(2),Sequential_PTD_ms);

%% ============================== LOAD DATA ================================

Single = load_linearity_dataset(single_folder,'Single',FS,nDepthChannels);
Sim = load_linearity_dataset(sim_folder,'Simultaneous',FS,nDepthChannels);
if same_folder(sim_folder,seq_folder)
    Seq = Sim;
    Seq.role = 'Sequential';
else
    Seq = load_linearity_dataset(seq_folder,'Sequential',FS,nDepthChannels);
end

%% ======================= DETECT ANALYZABLE PAIRS =========================

DetectedPairs = detect_analyzable_pairs(Single,Sim,Seq, ...
    Sequential_PTD_ms,PTD_Tolerance_ms);
if isempty(DetectedPairs)
    error('LinearityPSTH:NoAnalyzablePairs', ...
        ['No analyzable electrode pair was found. Each pair requires A ' ...
         'alone and B alone plus at least one complete matched branch: ' ...
         'Sim A->B with Seq A->B, or Sim B->A with Seq B->A.']);
end

DetectedPairTable = pair_table(DetectedPairs);
fprintf('\nDetected matched pair/set information\n');
disp(DetectedPairTable);
pair_indices = resolve_pair_selection(Pairs_To_Analyze,numel(DetectedPairs));

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
end
if isempty(selected_amplitudes)
    error('LinearityPSTH:NoAmplitudes','No amplitudes remain for analysis.');
end

% An extended compute window protects the displayed curve from shift and
% symmetric-convolution edge effects.
margin_ms = max(Sequential_PTD_ms,ceil(4*smooth_sigma_ms))+bin_ms;
compute_win_ms = [display_win_ms(1)-margin_ms, ...
                  display_win_ms(2)+margin_ms];
edges_compute = compute_win_ms(1):bin_ms:compute_win_ms(2);
ctrs_compute = edges_compute(1:end-1)+bin_ms/2;
display_mask = ctrs_compute >= display_win_ms(1) & ...
               ctrs_compute < display_win_ms(2);
ctrs_display = ctrs_compute(display_mask);
kernel = gaussian_kernel(bin_ms,smooth_sigma_ms);

%% ============================= MAIN ANALYSIS =============================

PSTHResults = struct([]);
random_stream = RandStream('mt19937ar','Seed',Random_Seed);

for ip = pair_indices
    Pair = DetectedPairs(ip);
    Branches = make_available_branches(Pair);

    for ib = 1:numel(Branches)
        Branch = Branches(ib);
        branch_channels = resolve_channels_for_branch(Channels_To_Analyze, ...
            Branch,Sim,Seq,nDepthChannels,Sequential_PTD_ms,PTD_Tolerance_ms);
        if isempty(branch_channels)
            warning('LinearityPSTH:NoSelectedChannels', ...
                'Pair %d, %s has no selected channels and was skipped.', ...
                ip,Branch.label);
            continue;
        end

        nCh = numel(branch_channels);
        nAmp = numel(selected_amplitudes);
        nBins = numel(ctrs_compute);

        A_raw = nan(nCh,nAmp,nBins);
        B_raw = nan(nCh,nAmp,nBins);
        Sim_observed_raw = nan(nCh,nAmp,nBins);
        Seq_observed_raw = nan(nCh,nAmp,nBins);
        Sim_linear_raw = nan(nCh,nAmp,nBins);
        Seq_linear_raw = nan(nCh,nAmp,nBins);
        First_component_raw = nan(nCh,nAmp,nBins);
        Second_shifted_raw = nan(nCh,nAmp,nBins);
        trial_n = zeros(nCh,nAmp,4);
        trial_ids_used = cell(nCh,nAmp,4);

        % Residual summaries: third dimension is Sim, Seq.
        full_difference = nan(nCh,nAmp,2);
        early_difference = nan(nCh,nAmp,2);
        late_difference = nan(nCh,nAmp,2);

        for jc = 1:nCh
            depth_channel = branch_channels(jc);
            spike_channel = depth_to_spike_channel(depth_channel);

            for ia = 1:nAmp
                amp = selected_amplitudes(ia);
                trial_ids = cell(1,4);
                trial_ids{1} = select_clean_trials(Single,Pair.A,Pair.key, ...
                    'single',amp,0,depth_channel,PTD_Tolerance_ms);
                trial_ids{2} = select_clean_trials(Single,Pair.B,Pair.key, ...
                    'single',amp,0,depth_channel,PTD_Tolerance_ms);
                trial_ids{3} = select_clean_trials(Sim,Branch.sim_order_key, ...
                    Pair.key,'sim_order',amp,0,depth_channel,PTD_Tolerance_ms);
                trial_ids{4} = select_clean_trials(Seq,Branch.seq_order_key, ...
                    Pair.key,'seq',amp,Sequential_PTD_ms,depth_channel, ...
                    PTD_Tolerance_ms);

                if is_balanced_trial_mode(Trial_Mode)
                    nBalanced = min(cellfun(@numel,trial_ids));
                    if nBalanced > 0
                        for ic = 1:4
                            pick = randperm(random_stream,numel(trial_ids{ic}),nBalanced);
                            trial_ids{ic} = sort(trial_ids{ic}(pick));
                        end
                    else
                        trial_ids(:) = {[]};
                    end
                end

                for ic = 1:4
                    trial_n(jc,ia,ic) = numel(trial_ids{ic});
                    trial_ids_used{jc,ia,ic} = trial_ids{ic};
                end
                if any(cellfun(@isempty,trial_ids))
                    continue;
                end

                A = calculate_evoked_psth(Single,spike_channel,trial_ids{1}, ...
                    baseline_win_ms,edges_compute);
                B = calculate_evoked_psth(Single,spike_channel,trial_ids{2}, ...
                    baseline_win_ms,edges_compute);
                ObsSim = calculate_evoked_psth(Sim,spike_channel,trial_ids{3}, ...
                    baseline_win_ms,edges_compute);
                ObsSeq = calculate_evoked_psth(Seq,spike_channel,trial_ids{4}, ...
                    baseline_win_ms,edges_compute);
                if any(~isfinite([A B ObsSim ObsSeq]))
                    continue;
                end

                % Only the prespecified reliable response interval from a
                % single pulse contributes to the short-latency prediction.
                A_component = gate_response(A,ctrs_compute,response_win_ms);
                B_component = gate_response(B,ctrs_compute,response_win_ms);
                SimLinear = A_component+B_component;

                if strcmp(Branch.code,'AB')
                    First = A_component;
                    SecondShifted = shift_curve(B_component,ctrs_compute, ...
                        Sequential_PTD_ms);
                else
                    First = B_component;
                    SecondShifted = shift_curve(A_component,ctrs_compute, ...
                        Sequential_PTD_ms);
                end
                SeqLinear = First+SecondShifted;

                A_raw(jc,ia,:) = reshape(A,1,1,[]);
                B_raw(jc,ia,:) = reshape(B,1,1,[]);
                Sim_observed_raw(jc,ia,:) = reshape(ObsSim,1,1,[]);
                Seq_observed_raw(jc,ia,:) = reshape(ObsSeq,1,1,[]);
                Sim_linear_raw(jc,ia,:) = reshape(SimLinear,1,1,[]);
                Seq_linear_raw(jc,ia,:) = reshape(SeqLinear,1,1,[]);
                First_component_raw(jc,ia,:) = reshape(First,1,1,[]);
                Second_shifted_raw(jc,ia,:) = reshape(SecondShifted,1,1,[]);

                windows = {response_win_ms,early_win_ms,late_win_ms};
                for iw = 1:3
                    W = windows{iw};
                    meanA = mean_evoked_count(Single,spike_channel,trial_ids{1}, ...
                        baseline_win_ms,W);
                    meanB = mean_evoked_count(Single,spike_channel,trial_ids{2}, ...
                        baseline_win_ms,W);
                    meanSim = mean_evoked_count(Sim,spike_channel,trial_ids{3}, ...
                        baseline_win_ms,W);
                    meanSeq = mean_evoked_count(Seq,spike_channel,trial_ids{4}, ...
                        baseline_win_ms,W);

                    simPrediction = meanA+meanB;
                    if strcmp(Branch.code,'AB')
                        firstPrediction = meanA;
                        secondPrediction = mean_evoked_count(Single,spike_channel, ...
                            trial_ids{2},baseline_win_ms,shifted_reliable_window( ...
                            W,Sequential_PTD_ms,response_win_ms));
                    else
                        firstPrediction = meanB;
                        secondPrediction = mean_evoked_count(Single,spike_channel, ...
                            trial_ids{1},baseline_win_ms,shifted_reliable_window( ...
                            W,Sequential_PTD_ms,response_win_ms));
                    end
                    seqPrediction = firstPrediction+secondPrediction;
                    values = [meanSim-simPrediction,meanSeq-seqPrediction];
                    if iw == 1
                        full_difference(jc,ia,:) = reshape(values,1,1,2);
                    elseif iw == 2
                        early_difference(jc,ia,:) = reshape(values,1,1,2);
                    else
                        late_difference(jc,ia,:) = reshape(values,1,1,2);
                    end
                end
            end
        end

        valid_by_amplitude = false(nCh,nAmp);
        for ia = 1:nAmp
            valid_by_amplitude(:,ia) = all(trial_n(:,ia,:) > 0,3) & ...
                all(isfinite(full_difference(:,ia,:)),3);
        end

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
        P.depth_channels = branch_channels;
        P.spike_channels = depth_to_spike_channel(branch_channels);
        P.amplitudes_uA = selected_amplitudes;
        P.ctrs_compute_ms = ctrs_compute;
        P.display_mask = display_mask;
        P.ctrs_display_ms = ctrs_display;
        P.A_raw = A_raw;
        P.B_raw = B_raw;
        P.Sim_observed_raw = Sim_observed_raw;
        P.Seq_observed_raw = Seq_observed_raw;
        P.Sim_linear_raw = Sim_linear_raw;
        P.Seq_linear_raw = Seq_linear_raw;
        P.First_component_raw = First_component_raw;
        P.Second_shifted_raw = Second_shifted_raw;
        P.trial_n = trial_n;
        P.trial_ids_used = trial_ids_used;
        P.valid_by_amplitude = valid_by_amplitude;
        P.full_difference = full_difference;
        P.early_difference = early_difference;
        P.late_difference = late_difference;
        P.summary = build_summary_table(P);

        if isempty(PSTHResults)
            PSTHResults = P;
        else
            PSTHResults(end+1) = P; %#ok<SAGROW>
        end

        fprintf('\n%s + %s | %s branch\n',Pair.A,Pair.B,Branch.label);
        fprintf('Sets: Sim %s; Sequential %s\n', ...
            vector_text(Branch.sim_set_indices), ...
            vector_text(Branch.seq_set_indices));
        fprintf('Responding channels: %d\n',numel(branch_channels));
        disp(P.summary);

        plot_branch_population(P,kernel,display_win_ms, ...
            response_win_ms,Sequential_PTD_ms,Show_SEM, ...
            Amplitudes_Per_Figure);

        if Plot_Channel_Diagnostics
            plot_selected_channel_diagnostics(P,Single,Sim,Seq,kernel, ...
                display_win_ms,response_win_ms,Sequential_PTD_ms, ...
                Diagnostic_Channels,Diagnostic_Amplitudes, ...
                depth_to_spike_channel);
        end
    end
end

if isempty(PSTHResults)
    error('LinearityPSTH:NoCompletedResults', ...
        'No available branch produced a complete PSTH result.');
end

fprintf('\nPSTH temporal linearity analysis complete; saving disabled.\n');

%% ============================== FUNCTIONS ================================

function validate_settings(single_folder,sim_folder,seq_folder,functions_folder, ...
        baseline_win,response_win,early_win,late_win,display_win,seq_ptd, ...
        ptd_tol,bin_ms,smooth_sigma,trial_mode,amps_per_figure,FS)
folders = {single_folder,sim_folder,seq_folder,functions_folder};
names = {'single_folder','sim_folder','seq_folder','analysis_functions_folder'};
for k = 1:numel(folders)
    if ~isfolder(folders{k})
        error('LinearityPSTH:FolderNotFound','%s not found: %s', ...
            names{k},folders{k});
    end
end
validate_window(baseline_win,'baseline_win_ms');
validate_window(response_win,'response_win_ms');
validate_window(early_win,'early_win_ms');
validate_window(late_win,'late_win_ms');
validate_window(display_win,'display_win_ms');
if abs(early_win(1)-response_win(1)) > 1e-9 || ...
        abs(early_win(2)-late_win(1)) > 1e-9 || ...
        abs(late_win(2)-response_win(2)) > 1e-9
    error('LinearityPSTH:WindowPartition', ...
        'early_win_ms and late_win_ms must partition response_win_ms.');
end
if baseline_win(2) > response_win(1)
    error('LinearityPSTH:OverlappingWindows', ...
        'Baseline and response windows must not overlap.');
end
if ~isscalar(seq_ptd) || ~isfinite(seq_ptd) || seq_ptd <= 0
    error('LinearityPSTH:InvalidPTD','Sequential_PTD_ms must be positive.');
end
if ~isscalar(ptd_tol) || ~isfinite(ptd_tol) || ptd_tol <= 0
    error('LinearityPSTH:InvalidPTDTolerance','PTD tolerance must be positive.');
end
if ~isscalar(bin_ms) || ~isfinite(bin_ms) || bin_ms <= 0
    error('LinearityPSTH:InvalidBin','bin_ms must be positive.');
end
if ~isscalar(smooth_sigma) || ~isfinite(smooth_sigma) || smooth_sigma <= 0
    error('LinearityPSTH:InvalidSmoothing','smooth_sigma_ms must be positive.');
end
if ~any(strcmpi(string(trial_mode),["all","balanced"]))
    error('LinearityPSTH:InvalidTrialMode','Trial_Mode must be all or balanced.');
end
if ~isscalar(amps_per_figure) || amps_per_figure < 1 || ...
        fix(amps_per_figure) ~= amps_per_figure
    error('LinearityPSTH:InvalidPageSize', ...
        'Amplitudes_Per_Figure must be a positive integer.');
end
if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('LinearityPSTH:InvalidFS','FS must be positive.');
end
end

function validate_window(window,name)
if ~isnumeric(window) || numel(window) ~= 2 || any(~isfinite(window)) || ...
        window(2) <= window(1)
    error('LinearityPSTH:InvalidWindow', ...
        '%s must be [start end] with end > start.',name);
end
end

function D = load_linearity_dataset(folder,role,FS,nDepthChannels)
folder = char(string(folder));
[spike_file,spike_variable] = find_spike_source(folder);
Ssp = load(spike_file,spike_variable);
sp = Ssp.(spike_variable);
if ~iscell(sp)
    error('LinearityPSTH:SpikeVariableNotCell', ...
        '%s in %s is not a cell array.',spike_variable,spike_file);
end

trigger_files = clean_file_list(dir(fullfile(folder,'*.trig.dat')));
if isempty(trigger_files)
    error('LinearityPSTH:MissingTrigger', ...
        'No existing *.trig.dat file was found in %s.',folder);
end
if numel(trigger_files) > 1
    error('LinearityPSTH:AmbiguousTrigger', ...
        'Multiple trigger files were found in %s.',folder);
end
trig = read_triggers_without_writing(folder);
trig = double(trig(:));

experiment_files = clean_file_list(dir(fullfile(folder,'*_exp_datafile_*.mat')));
if isempty(experiment_files)
    error('LinearityPSTH:MissingExperimentFile', ...
        'No *_exp_datafile_*.mat file was found in %s.',folder);
end
if numel(experiment_files) > 1
    error('LinearityPSTH:AmbiguousExperimentFile', ...
        'Multiple experiment files were found in %s.',folder);
end
experiment_file = fullfile(experiment_files(1).folder,experiment_files(1).name);
E = load(experiment_file,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required = {'StimParams','simultaneous_stim','E_MAP'};
for k = 1:numel(required)
    if ~isfield(E,required{k})
        error('LinearityPSTH:MissingExperimentVariable', ...
            '%s is missing from %s.',required{k},experiment_file);
    end
end

StimParams = E.StimParams;
simN = double(E.simultaneous_stim);
if ~isscalar(simN) || simN < 1 || fix(simN) ~= simN
    error('LinearityPSTH:InvalidSimultaneousStim', ...
        'simultaneous_stim is invalid in %s.',experiment_file);
end
nRows = size(StimParams,1)-1;
if mod(nRows,simN) ~= 0
    error('LinearityPSTH:StimParamRowMismatch', ...
        'StimParams rows do not match simultaneous_stim in %s.',experiment_file);
end
nCalculated = nRows/simN;
if isfield(E,'n_Trials') && ~isempty(E.n_Trials)
    nTrials = double(E.n_Trials);
    if nTrials ~= nCalculated
        error('LinearityPSTH:TrialCountMismatch', ...
            'n_Trials disagrees with StimParams in %s.',experiment_file);
    end
else
    nTrials = nCalculated;
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
            error('LinearityPSTH:EmptyStimLabel', ...
                'Empty stimulation label in trial %d of %s.',t,experiment_file);
        end
        if ~isempty(emap_labels) && ~any(labels(pulse) == emap_labels)
            error('LinearityPSTH:StimLabelNotInMap', ...
                'Stimulus %s in trial %d was not found in E_MAP.',labels(pulse),t);
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
        error('LinearityPSTH:AmbiguousSpikeFile', ...
            'Multiple files contain %s in %s:\n%s',priority{ip},folder, ...
            strjoin(cellstr(matches),newline));
    elseif isscalar(matches)
        file_path = char(matches(1));
        variable_name = priority{ip};
        return;
    end
end
error('LinearityPSTH:NoUsableSpikeVariable', ...
    'No sp_corr, sp_SSD, sp_clipped, or sp variable was found in %s.',folder);
end

function trig = read_triggers_without_writing(folder)
original = pwd;
cleanup = onCleanup(@() cd(original));
cd(folder);
trig = loadTrig(0);
clear cleanup;
cd(original);
end

function labels = normalize_emap_labels(E_MAP)
labels = strings(0,1);
if isempty(E_MAP), return; end
values = E_MAP(:);
if numel(values) >= 2, values = values(2:end); end
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
    error('LinearityPSTH:InvalidStimLabel', ...
        'A stimulation label could not be converted to text.');
end
label = regexprep(label,'\s+','');
end

function value = numeric_cell_value(raw,name,trial,file)
if iscell(raw) && isscalar(raw), raw = raw{1}; end
if isnumeric(raw) && isscalar(raw)
    value = double(raw);
elseif ischar(raw) || (isstring(raw) && isscalar(raw))
    value = str2double(raw);
else
    value = NaN;
end
if ~isfinite(value)
    error('LinearityPSTH:InvalidNumericStimParam', ...
        'Invalid %s in trial %d of %s.',name,trial,file);
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
    error('LinearityPSTH:AmbiguousBadChannelFile', ...
        'Multiple standard BadChannels files were found in %s.',folder);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if isfield(S,'BadCh_perSet')
    BadCh = S.BadCh_perSet;
elseif isfield(S,'BadCh')
    BadCh = S.BadCh;
else
    error('LinearityPSTH:InvalidBadChannelFile', ...
        'No BadCh_perSet or BadCh variable was found in %s.',file_path);
end
end

function [BadTrials,file_path] = load_bad_trials(folder,nDepthChannels)
files = find_named_mat_files(folder,'badtrials');
if isempty(files)
    BadTrials = cell(nDepthChannels,1);
    file_path = '';
    return;
end
if numel(files) > 1
    error('LinearityPSTH:AmbiguousBadTrialFile', ...
        'Multiple standard BadTrials files were found in %s:\n%s',folder, ...
        strjoin(string({files.name}),newline));
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if ~isfield(S,'BadTrials') || ~iscell(S.BadTrials)
    error('LinearityPSTH:InvalidBadTrials', ...
        'BadTrials must be a channel-indexed cell array in %s.',file_path);
end
BadTrials = S.BadTrials;
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
    error('LinearityPSTH:AmbiguousRespondingFile', ...
        'Multiple standard RespondingChannels files were found in %s.',folder);
end
file_path = fullfile(files(1).folder,files(1).name);
S = load(file_path);
if ~isfield(S,'Responding')
    error('LinearityPSTH:InvalidRespondingFile', ...
        'Responding is missing from %s.',file_path);
end
Responding = S.Responding;
end

function files = find_named_mat_files(folder,file_label)
files = clean_file_list(dir(fullfile(folder,'*.mat')));
keep = false(size(files));
pattern = ['(^|[._-])' regexptranslate('escape',lower(file_label)) '\.mat$'];
for k = 1:numel(files)
    keep(k) = ~isempty(regexpi(files(k).name,pattern,'once'));
end
files = files(keep);
end

function files = clean_file_list(files)
if isempty(files), return; end
keep = true(size(files));
for k = 1:numel(files)
    keep(k) = ~startsWith(files(k).name,'._');
end
files = files(keep);
end

function tf = same_folder(a,b)
tf = strcmp(strip_trailing_separator(char(string(a))), ...
    strip_trailing_separator(char(string(b))));
end

function out = strip_trailing_separator(in)
out = in;
while numel(out) > 1 && any(out(end) == ['/' '\'])
    out(end) = [];
end
end

function Pairs = detect_analyzable_pairs(Single,Sim,Seq,target_ptd,tol)
sim_mask = abs(Sim.trial_ptd_ms) < tol & ...
    cellfun(@numel,Sim.trial_labels) == 2;
seq_mask = abs(Seq.trial_ptd_ms-target_ptd) < tol & ...
    cellfun(@numel,Seq.trial_labels) == 2;
candidate_keys = intersect(unique(Sim.pair_key(sim_mask),'stable'), ...
    unique(Seq.pair_key(seq_mask),'stable'),'stable');
single_labels = unique(Single.order_key( ...
    cellfun(@numel,Single.trial_labels) == 1));

Pairs = repmat(struct('key','','A','','B','','A_to_B_key','', ...
    'B_to_A_key','','sim_a_to_b_set_indices',[], ...
    'sim_b_to_a_set_indices',[],'a_to_b_set_indices',[], ...
    'b_to_a_set_indices',[],'has_AB',false,'has_BA',false),1,0);

for k = 1:numel(candidate_keys)
    labels = sort(split(candidate_keys(k),'|'));
    if numel(labels) ~= 2 || labels(1) == labels(2), continue; end
    A = labels(1); B = labels(2);
    AtoB = A+">"+B;
    BtoA = B+">"+A;
    if ~(any(single_labels == A) && any(single_labels == B)), continue; end

    simAB = unique(Sim.set_index(sim_mask & Sim.order_key == AtoB)).';
    simBA = unique(Sim.set_index(sim_mask & Sim.order_key == BtoA)).';
    seqAB = unique(Seq.set_index(seq_mask & Seq.order_key == AtoB)).';
    seqBA = unique(Seq.set_index(seq_mask & Seq.order_key == BtoA)).';
    hasAB = ~isempty(simAB) && ~isempty(seqAB);
    hasBA = ~isempty(simBA) && ~isempty(seqBA);
    if ~(hasAB || hasBA), continue; end

    P = struct('key',char(candidate_keys(k)),'A',char(A),'B',char(B), ...
        'A_to_B_key',char(AtoB),'B_to_A_key',char(BtoA), ...
        'sim_a_to_b_set_indices',simAB, ...
        'sim_b_to_a_set_indices',simBA, ...
        'a_to_b_set_indices',seqAB,'b_to_a_set_indices',seqBA, ...
        'has_AB',hasAB,'has_BA',hasBA);
    Pairs(end+1) = P; %#ok<AGROW>
end
end

function T = pair_table(Pairs)
n = numel(Pairs);
Pair = (1:n).';
ElectrodeA = strings(n,1); ElectrodeB = strings(n,1);
SimAtoBSets = strings(n,1); SimBtoASets = strings(n,1);
AtoBSets = strings(n,1); BtoASets = strings(n,1);
AvailableBranches = strings(n,1);
for k = 1:n
    ElectrodeA(k) = string(Pairs(k).A);
    ElectrodeB(k) = string(Pairs(k).B);
    SimAtoBSets(k) = vector_text(Pairs(k).sim_a_to_b_set_indices);
    SimBtoASets(k) = vector_text(Pairs(k).sim_b_to_a_set_indices);
    AtoBSets(k) = vector_text(Pairs(k).a_to_b_set_indices);
    BtoASets(k) = vector_text(Pairs(k).b_to_a_set_indices);
    names = strings(0,1);
    if Pairs(k).has_AB, names(end+1) = "A->B"; end %#ok<AGROW>
    if Pairs(k).has_BA, names(end+1) = "B->A"; end %#ok<AGROW>
    AvailableBranches(k) = join(names,", ");
end
T = table(Pair,ElectrodeA,ElectrodeB,SimAtoBSets,SimBtoASets, ...
    AtoBSets,BtoASets,AvailableBranches);
end

function Branches = make_available_branches(Pair)
Branches = repmat(struct('code','','label','','sim_order_key','', ...
    'seq_order_key','','sim_set_indices',[],'seq_set_indices',[], ...
    'seq_plot_label','','seq_color',[]),1,0);
if Pair.has_AB
    B = struct('code','AB','label','A->B matched', ...
        'sim_order_key',Pair.A_to_B_key,'seq_order_key',Pair.A_to_B_key, ...
        'sim_set_indices',Pair.sim_a_to_b_set_indices, ...
        'seq_set_indices',Pair.a_to_b_set_indices, ...
        'seq_plot_label','A then B','seq_color',[0.85 0.33 0.10]);
    Branches(end+1) = B;
end
if Pair.has_BA
    B = struct('code','BA','label','B->A matched', ...
        'sim_order_key',Pair.B_to_A_key,'seq_order_key',Pair.B_to_A_key, ...
        'sim_set_indices',Pair.sim_b_to_a_set_indices, ...
        'seq_set_indices',Pair.b_to_a_set_indices, ...
        'seq_plot_label','B then A','seq_color',[0.55 0.20 0.65]);
    Branches(end+1) = B;
end
end

function indices = resolve_pair_selection(request,nPairs)
if ischar(request) || (isstring(request) && isscalar(request))
    if ~strcmpi(strtrim(string(request)),'all')
        error('LinearityPSTH:InvalidPairSelection','Text selection must be all.');
    end
    indices = 1:nPairs;
elseif isnumeric(request)
    indices = unique(double(request(:).'),'stable');
    if isempty(indices) || any(indices < 1 | indices > nPairs | ...
            fix(indices) ~= indices)
        error('LinearityPSTH:InvalidPairSelection','Invalid pair number.');
    end
else
    error('LinearityPSTH:InvalidPairSelection', ...
        'Pairs_To_Analyze must be all or a numeric vector.');
end
end

function channels = resolve_channels_for_branch(request,Branch,Sim,Seq, ...
        nDepth,target_ptd,tol)
if isnumeric(request)
    channels = unique(double(request(:).'),'stable');
    if isempty(channels) || any(channels < 1 | channels > nDepth | ...
            fix(channels) ~= channels)
        error('LinearityPSTH:InvalidChannels','Invalid depth channel selection.');
    end
    return;
end
switch lower(strtrim(char(string(request))))
    case 'all'
        channels = 1:nDepth;
    case 'responding'
        sim_mask = responding_mask_for_sets(Sim,Branch.sim_set_indices,0,tol,nDepth);
        seq_mask = responding_mask_for_sets(Seq,Branch.seq_set_indices, ...
            target_ptd,tol,nDepth);
        channels = find(sim_mask | seq_mask).';
        if isempty(channels)
            error('LinearityPSTH:NoRespondingLabels', ...
                ['No responding channels were found for %s. Use ' ...
                 'Channels_To_Analyze=''all'' or verify the standard ' ...
                 'RespondingChannels file.'],Branch.label);
        end
    otherwise
        error('LinearityPSTH:InvalidChannels', ...
            'Channels_To_Analyze must be responding, all, or numeric.');
end
end

function mask = responding_mask_for_sets(D,set_indices,target_ptd,tol,nDepth)
mask = false(nDepth,1);
if isempty(D.Responding) || ~isfield(D.Responding,'set'), return; end
for set_idx = set_indices
    if set_idx < 1 || set_idx > numel(D.Responding.set), continue; end
    ptd_values = unique(D.trial_ptd_ms(D.set_index == set_idx));
    ptd_position = find(abs(ptd_values-target_ptd) < tol,1);
    if isempty(ptd_position), continue; end
    amp_struct = D.Responding.set(set_idx).amp;
    for ia = 1:numel(amp_struct)
        if ~isfield(amp_struct(ia),'ptd') || ...
                numel(amp_struct(ia).ptd) < ptd_position
            continue;
        end
        ptd_struct = amp_struct(ia).ptd(ptd_position);
        if ~isfield(ptd_struct,'channel'), continue; end
        C = ptd_struct.channel;
        for ch = 1:min(numel(C),nDepth)
            if isfield(C(ch),'is_responsive') && logical(C(ch).is_responsive)
                mask(ch) = true;
            end
        end
    end
end
end

function trial_ids = select_clean_trials(D,ordered_key,pair_key,mode,amp, ...
        ptd,depth,tol)
amp_mask = abs(D.trial_amp_uA-amp) < 1e-6;
ptd_mask = abs(D.trial_ptd_ms-ptd) < tol;
switch lower(mode)
    case 'single'
        condition_mask = D.order_key == string(ordered_key);
    case 'sim_order'
        condition_mask = D.pair_key == string(pair_key) & ...
            D.order_key == string(ordered_key);
    case 'seq'
        condition_mask = D.pair_key == string(pair_key) & ...
            D.order_key == string(ordered_key);
    otherwise
        error('LinearityPSTH:InternalConditionMode','Unknown condition mode.');
end
trial_ids = find(amp_mask & ptd_mask & condition_mask);
trial_ids = trial_ids(trial_ids <= numel(D.trig));
if isempty(trial_ids), return; end

good = true(size(trial_ids));
for k = 1:numel(trial_ids)
    good(k) = ~is_bad_channel(D.BadCh,D.set_index(trial_ids(k)),depth);
end
trial_ids = trial_ids(good);
if isempty(trial_ids), return; end
if depth <= numel(D.BadTrials) && ~isempty(D.BadTrials{depth})
    bad = double(D.BadTrials{depth}(:));
    trial_ids = trial_ids(~ismember(trial_ids,bad));
end
trial_ids = trial_ids(:).';
end

function tf = is_bad_channel(BadCh,set_index,depth)
tf = false;
if isempty(BadCh), return; end
if iscell(BadCh)
    if set_index <= numel(BadCh) && ~isempty(BadCh{set_index})
        tf = ismember(depth,double(BadCh{set_index}(:)));
    end
elseif isnumeric(BadCh) || islogical(BadCh)
    tf = ismember(depth,double(BadCh(:)));
end
end

function psth = calculate_evoked_psth(D,spike_channel,trial_ids, ...
        baseline_win,edges)
nBins = numel(edges)-1;
psth = nan(1,nBins);
if isempty(trial_ids) || spike_channel < 1 || ...
        spike_channel > numel(D.sp) || isempty(D.sp{spike_channel})
    return;
end
spike_data = D.sp{spike_channel};
if ~isnumeric(spike_data) || size(spike_data,2) < 1, return; end
spike_times = double(spike_data(:,1));
bin_s = diff(edges(1:2))/1000;
baseline_s = diff(baseline_win)/1000;
trial_rates = nan(numel(trial_ids),nBins);
for k = 1:numel(trial_ids)
    tr = trial_ids(k);
    if tr < 1 || tr > numel(D.trig), continue; end
    t0 = D.trig(tr)/D.FS*1000;
    relative = spike_times-t0;
    counts = histcounts(relative,edges);
    baseline_count = sum(relative >= baseline_win(1) & ...
        relative < baseline_win(2));
    trial_rates(k,:) = counts/bin_s-baseline_count/baseline_s;
end
if any(all(isfinite(trial_rates),2))
    psth = mean(trial_rates,1,'omitnan');
end
end

function value = mean_evoked_count(D,spike_channel,trial_ids,baseline_win,count_win)
% An empty or zero-width reliable component contributes zero to a linear
% prediction. Missing trials/channels remain NaN.
if count_win(2) <= count_win(1)
    value = 0;
    return;
end
values = calculate_evoked_counts(D,spike_channel,trial_ids, ...
    baseline_win,count_win);
if isempty(values) || ~any(isfinite(values))
    value = NaN;
else
    value = mean(values,'omitnan');
end
end

function values = calculate_evoked_counts(D,spike_channel,trial_ids, ...
        baseline_win,count_win)
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
spike_times = double(spike_data(:,1));
base_duration = diff(baseline_win);
count_duration = diff(count_win);
values = nan(numel(trial_ids),1);
for k = 1:numel(trial_ids)
    tr = trial_ids(k);
    if tr < 1 || tr > numel(D.trig), continue; end
    t0 = D.trig(tr)/D.FS*1000;
    base_count = sum(spike_times >= t0+baseline_win(1) & ...
        spike_times < t0+baseline_win(2));
    response_count = sum(spike_times >= t0+count_win(1) & ...
        spike_times < t0+count_win(2));
    values(k) = response_count-base_count*count_duration/base_duration;
end
end

function Wsingle = shifted_reliable_window(Wglobal,delay,reliable_window)
% Convert a global window to the second pulse's time coordinates and retain
% only the reliable single-pulse response interval.
Wsingle = [max(Wglobal(1)-delay,reliable_window(1)), ...
           min(Wglobal(2)-delay,reliable_window(2))];
if Wsingle(2) < Wsingle(1)
    Wsingle(2) = Wsingle(1);
end
end

function y = gate_response(x,t,response_win)
y = zeros(size(x));
mask = t >= response_win(1) & t < response_win(2);
y(mask) = x(mask);
end

function shifted = shift_curve(x,t,delay)
shifted = interp1(t,x,t-delay,'linear',0);
end

function kernel = gaussian_kernel(bin_ms,sigma_ms)
half_width = ceil(4*sigma_ms/bin_ms);
x = (-half_width:half_width)*bin_ms;
kernel = exp(-0.5*(x/sigma_ms).^2);
kernel = kernel/sum(kernel);
end

function tf = is_balanced_trial_mode(request)
tf = strcmpi(strtrim(char(string(request))),'balanced');
end

function T = build_summary_table(P)
nAmp = numel(P.amplitudes_uA);
Amplitude_uA = P.amplitudes_uA(:);
NResponding = repmat(numel(P.depth_channels),nAmp,1);
NUsed = zeros(nAmp,1);
SimFull = nan(nAmp,1); SeqFull = nan(nAmp,1);
SimEarly = nan(nAmp,1); SeqEarly = nan(nAmp,1);
SimLate = nan(nAmp,1); SeqLate = nan(nAmp,1);
for ia = 1:nAmp
    valid = P.valid_by_amplitude(:,ia);
    NUsed(ia) = sum(valid);
    if ~any(valid), continue; end
    SimFull(ia) = mean(P.full_difference(valid,ia,1),'omitnan');
    SeqFull(ia) = mean(P.full_difference(valid,ia,2),'omitnan');
    SimEarly(ia) = mean(P.early_difference(valid,ia,1),'omitnan');
    SeqEarly(ia) = mean(P.early_difference(valid,ia,2),'omitnan');
    SimLate(ia) = mean(P.late_difference(valid,ia,1),'omitnan');
    SeqLate(ia) = mean(P.late_difference(valid,ia,2),'omitnan');
end
T = table(Amplitude_uA,NResponding,NUsed,SimFull,SeqFull, ...
    SimEarly,SeqEarly,SimLate,SeqLate);
end

function plot_branch_population(P,kernel,display_win,response_win,ptd, ...
        show_sem,amps_per_figure)
nAmp = numel(P.amplitudes_uA);
nPages = ceil(nAmp/amps_per_figure);
blue = [0.00 0.35 0.85];
orange = [0.85 0.33 0.10];
grey = [0.35 0.35 0.35];
green = [0.20 0.60 0.30];
t = P.ctrs_display_ms;

for page = 1:nPages
    amp_indices = (page-1)*amps_per_figure+1: ...
        min(page*amps_per_figure,nAmp);
    nRows = numel(amp_indices);
    figure('Color','w','Position',[60 60 1550 max(420,310*nRows)]);
    tl = tiledlayout(nRows,3,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s + %s | %s PSTH temporal linearity (page %d/%d)', ...
        P.electrode_A,P.electrode_B,P.branch_label,page,nPages), ...
        'FontWeight','bold');

    for row = 1:nRows
        ia = amp_indices(row);
        valid = P.valid_by_amplitude(:,ia);
        amp = P.amplitudes_uA(ia);
        nUsed = sum(valid);

        % Panel 1: sequential prediction components.
        ax1 = nexttile(tl,(row-1)*3+1); hold(ax1,'on'); box(ax1,'off');
        if nUsed > 0
            [first,~,~] = curve_summary(P.First_component_raw(valid,ia,:), ...
                kernel,P.display_mask);
            [second,~,~] = curve_summary(P.Second_shifted_raw(valid,ia,:), ...
                kernel,P.display_mask);
            [seqLinear,~,~] = curve_summary(P.Seq_linear_raw(valid,ia,:), ...
                kernel,P.display_mask);
            plot(ax1,t,first,'Color',green,'LineWidth',1.4, ...
                'DisplayName','First single component');
            plot(ax1,t,second,'Color',orange,'LineWidth',1.4, ...
                'DisplayName','Shifted second component');
            plot(ax1,t,seqLinear,'k--','LineWidth',1.8, ...
                'DisplayName','Sequential linear');
        end
        add_pulse_lines(ax1,ptd);
        xlim(ax1,display_win);
        ylabel(ax1,'Evoked rate (sp/s)');
        title(ax1,sprintf('%.3g uA | components | n=%d',amp,nUsed));
        if row == 1, legend(ax1,'Location','best','Box','off'); end

        % Panel 2: observed versus condition-specific predictions.
        ax2 = nexttile(tl,(row-1)*3+2); hold(ax2,'on'); box(ax2,'off');
        if nUsed > 0
            [simObs,simObsSEM,~] = curve_summary( ...
                P.Sim_observed_raw(valid,ia,:),kernel,P.display_mask);
            [seqObs,seqObsSEM,~] = curve_summary( ...
                P.Seq_observed_raw(valid,ia,:),kernel,P.display_mask);
            [simLin,simLinSEM,~] = curve_summary( ...
                P.Sim_linear_raw(valid,ia,:),kernel,P.display_mask);
            [seqLin,seqLinSEM,~] = curve_summary( ...
                P.Seq_linear_raw(valid,ia,:),kernel,P.display_mask);
            if show_sem
                shaded_sem(ax2,t,simObs,simObsSEM,blue);
                shaded_sem(ax2,t,seqObs,seqObsSEM,orange);
                shaded_sem(ax2,t,simLin,simLinSEM,[0 0 0]);
                shaded_sem(ax2,t,seqLin,seqLinSEM,grey);
            end
            plot(ax2,t,simLin,'k--','LineWidth',1.7, ...
                'DisplayName','Sim linear');
            plot(ax2,t,simObs,'Color',blue,'LineWidth',1.8, ...
                'DisplayName','Sim observed');
            plot(ax2,t,seqLin,'--','Color',grey,'LineWidth',1.7, ...
                'DisplayName','Seq linear');
            plot(ax2,t,seqObs,'Color',orange,'LineWidth',1.8, ...
                'DisplayName','Seq observed');
        end
        add_pulse_lines(ax2,ptd);
        xlim(ax2,display_win);
        title(ax2,'Observed versus prediction');
        if row == 1, legend(ax2,'Location','best','Box','off'); end

        % Panel 3: time-resolved observed-minus-linear residual.
        ax3 = nexttile(tl,(row-1)*3+3); hold(ax3,'on'); box(ax3,'off');
        if nUsed > 0
            simResidual = P.Sim_observed_raw(valid,ia,:)- ...
                P.Sim_linear_raw(valid,ia,:);
            seqResidual = P.Seq_observed_raw(valid,ia,:)- ...
                P.Seq_linear_raw(valid,ia,:);
            [simDiff,simDiffSEM,~] = curve_summary(simResidual,kernel, ...
                P.display_mask);
            [seqDiff,seqDiffSEM,~] = curve_summary(seqResidual,kernel, ...
                P.display_mask);
            analysis_mask = t >= response_win(1) & t < response_win(2);
            simDiff(~analysis_mask) = NaN;
            seqDiff(~analysis_mask) = NaN;
            simDiffSEM(~analysis_mask) = NaN;
            seqDiffSEM(~analysis_mask) = NaN;
            if show_sem
                shaded_sem(ax3,t,simDiff,simDiffSEM,blue);
                shaded_sem(ax3,t,seqDiff,seqDiffSEM,orange);
            end
            plot(ax3,t,simDiff,'Color',blue,'LineWidth',1.8, ...
                'DisplayName','Sim - linear');
            plot(ax3,t,seqDiff,'Color',orange,'LineWidth',1.8, ...
                'DisplayName','Seq - linear');
        end
        yline(ax3,0,'--','Color',[0.4 0.4 0.4], ...
            'HandleVisibility','off');
        add_pulse_lines(ax3,ptd);
        xlim(ax3,display_win);
        title(ax3,'Observed - condition-specific linear');
        if row == 1, legend(ax3,'Location','best','Box','off'); end

        if row == nRows
            xlabel(ax1,'Time from first pulse (ms)');
            xlabel(ax2,'Time from first pulse (ms)');
            xlabel(ax3,'Time from first pulse (ms)');
        end
    end
    unify_row_limits(tl,nRows);
end
end

function [avg,sem,n] = curve_summary(data,kernel,display_mask)
M = squeeze(data);
if isvector(M), M = reshape(M,1,[]); end
valid = all(isfinite(M),2);
M = M(valid,:);
n = size(M,1);
if n == 0
    avg = nan(1,sum(display_mask));
    sem = avg;
    return;
end
for k = 1:n
    M(k,:) = conv(M(k,:),kernel,'same');
end
M = M(:,display_mask);
avg = mean(M,1,'omitnan');
if n > 1
    sem = std(M,0,1,'omitnan')/sqrt(n);
else
    sem = nan(size(avg));
end
end

function shaded_sem(ax,x,y,se,color)
valid = isfinite(x) & isfinite(y) & isfinite(se);
if sum(valid) < 2, return; end
x = x(valid); y = y(valid); se = se(valid);
fill(ax,[x fliplr(x)],[y+se fliplr(y-se)],color, ...
    'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
end

function add_pulse_lines(ax,ptd)
xline(ax,0,'r--','LineWidth',1,'HandleVisibility','off');
xline(ax,ptd,':','Color',[0.2 0.2 0.2],'LineWidth',1, ...
    'HandleVisibility','off');
end

function unify_row_limits(tl,nRows)
% Use one common scale for every amplitude row on the current page.
% The component/observed panels share an evoked-rate scale, while all
% residual panels share a separate observed-minus-linear scale.
evoked_axes = gobjects(2*nRows,1);
residual_axes = gobjects(nRows,1);
evoked_limits = nan(2*nRows,2);
residual_limits = nan(nRows,2);

for row = 1:nRows
    ax1 = nexttile(tl,(row-1)*3+1);
    ax2 = nexttile(tl,(row-1)*3+2);
    ax3 = nexttile(tl,(row-1)*3+3);

    evoked_axes(2*row-1:2*row) = [ax1; ax2];
    residual_axes(row) = ax3;
    evoked_limits(2*row-1,:) = ylim(ax1);
    evoked_limits(2*row,:) = ylim(ax2);
    residual_limits(row,:) = ylim(ax3);
end

evoked_range = [min(evoked_limits(:,1),[],'omitnan'), ...
                 max(evoked_limits(:,2),[],'omitnan')];
residual_range = [min(residual_limits(:,1),[],'omitnan'), ...
                   max(residual_limits(:,2),[],'omitnan')];

if all(isfinite(evoked_range)) && evoked_range(2) > evoked_range(1)
    set(evoked_axes,'YLim',evoked_range);
end
if all(isfinite(residual_range)) && residual_range(2) > residual_range(1)
    set(residual_axes,'YLim',residual_range);
end
end

function plot_selected_channel_diagnostics(P,Single,Sim,Seq,kernel, ...
        display_win,response_win,ptd,requested_channels,requested_amps, ...
        depth_to_spike_channel)
if isempty(requested_channels)
    warning('LinearityPSTH:NoDiagnosticChannels', ...
        'Plot_Channel_Diagnostics is true but Diagnostic_Channels is empty.');
    return;
end
channels = intersect(P.depth_channels,double(requested_channels(:).'),'stable');
if isempty(requested_amps)
    amp_indices = 1:numel(P.amplitudes_uA);
else
    amp_indices = find(ismembertol(P.amplitudes_uA, ...
        double(requested_amps(:).'),1e-9));
end
for depth = channels
    jc = find(P.depth_channels == depth,1);
    spike_channel = depth_to_spike_channel(depth);
    for ia = amp_indices
        if ~P.valid_by_amplitude(jc,ia), continue; end
        figure('Color','w','Position',[120 80 950 1000]);
        tl = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('%s + %s | %s | depth %d | %.3g uA', ...
            P.electrode_A,P.electrode_B,P.branch_label,depth, ...
            P.amplitudes_uA(ia)));
        datasets = {Single,Single,Sim,Seq};
        labels = {'A alone','B alone','Matched simultaneous','Sequential'};
        for ic = 1:4
            ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
            draw_raster(ax,datasets{ic},spike_channel, ...
                P.trial_ids_used{jc,ia,ic},display_win);
            xline(ax,0,'r--','HandleVisibility','off');
            if ic == 4
                xline(ax,ptd,':','Color',[0.2 0.2 0.2], ...
                    'HandleVisibility','off');
            end
            title(ax,labels{ic}); xlim(ax,display_win);
            ylabel(ax,'Trial');
        end
        ax = nexttile(tl); hold(ax,'on'); box(ax,'off');
        t = P.ctrs_display_ms;
        fields = {'Sim_linear_raw','Sim_observed_raw', ...
            'Seq_linear_raw','Seq_observed_raw'};
        colors = [0 0 0; 0 0.35 0.85; 0.35 0.35 0.35; 0.85 0.33 0.10];
        styles = {'--','-','--','-'};
        names = {'Sim linear','Sim observed','Seq linear','Seq observed'};
        for k = 1:4
            x = squeeze(P.(fields{k})(jc,ia,:)).';
            x = conv(x,kernel,'same');
            plot(ax,t,x(P.display_mask),styles{k},'Color',colors(k,:), ...
                'LineWidth',1.6,'DisplayName',names{k});
        end
        add_pulse_lines(ax,ptd);
        xlim(ax,display_win); xline(ax,response_win(2),':', ...
            'HandleVisibility','off');
        xlabel(ax,'Time from first pulse (ms)');
        ylabel(ax,'Evoked rate (sp/s)');
        legend(ax,'Location','best','Box','off');
    end
end
end

function draw_raster(ax,D,spike_channel,trial_ids,display_win)
if isempty(trial_ids) || spike_channel < 1 || ...
        spike_channel > numel(D.sp) || isempty(D.sp{spike_channel})
    return;
end
spike_times = double(D.sp{spike_channel}(:,1));
for k = 1:numel(trial_ids)
    tr = trial_ids(k);
    t0 = D.trig(tr)/D.FS*1000;
    tt = spike_times(spike_times >= t0+display_win(1) & ...
        spike_times < t0+display_win(2))-t0;
    for j = 1:numel(tt)
        plot(ax,[tt(j) tt(j)],[k-0.35 k+0.35],'k-','LineWidth',0.7);
    end
end
ylim(ax,[0.5 max(1.5,numel(trial_ids)+0.5)]);
end

function out = vector_text(values)
if isempty(values), out = "none"; else, out = string(strtrim(num2str(values))); end
end

function name = folder_leaf_name(folder)
[~,name] = fileparts(strip_trailing_separator(char(string(folder))));
end
