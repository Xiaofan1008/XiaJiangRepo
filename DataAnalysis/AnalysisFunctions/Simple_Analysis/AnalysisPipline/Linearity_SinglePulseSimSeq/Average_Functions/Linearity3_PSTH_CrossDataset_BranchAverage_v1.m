%% ========================================================================
%  Cross-dataset PSTH temporal linearity (branch-weighted)
%
%  Input:
%    PSTHLinearitySaved files created by
%    Linearity3_PSTH_TemporalLinearity_Save_v1.m
%
%  Averaging unit:
%    One saved order-matched branch (A->B or B->A) = one observation.
%    Branches are NOT merged within a dataset and are NOT weighted by the
%    number of responding channels.
%
%  Main outputs:
%    1. Absolute observed and linear-prediction PSTHs (spikes/s)
%    2. Absolute residual PSTHs: observed - linear (spikes/s)
%    3. Optional normalized PSTHs and residuals
%    4. Full, early, and late residual summaries (spikes/trial)
%
%  Normalization:
%    A single scale is calculated for each Branch x Amplitude from the
%    largest absolute value of its two smoothed condition-specific linear
%    predictions during the response window. The SAME scale is applied to
%    Sim and Seq observed, prediction, and residual curves. Weak branches
%    are excluded only from normalized plots, never from absolute plots.
%
%  IMPORTANT:
%    This is a descriptive branch-weighted analysis. A dataset containing
%    both orders contributes two branches, while a one-order dataset
%    contributes one branch.
% ========================================================================
clear;

%% ========================== USER SETTINGS ===============================
Results_Folder = ...
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq';
File_Pattern = 'PSTHLinearity_*.mat';

% File-selection mode:
%   false = automatically load every file matching File_Pattern above.
%   true  = load only the explicitly listed files below.
Use_Manual_Result_Files = true;
Manual_Result_Files = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX005_Xia_Exp1_Sim.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX009_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX009_Xia_Exp1_Sim5.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX010_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX010_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX010_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX010_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX010_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX010_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX010_Xia_Exp1_Sim8.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX011_Xia_Exp1_Sim9.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX012_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX012_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX012_Xia_Exp1_Sim6.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX013_Xia_Exp1_Seq_Sim8.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX014_Xia_Seq_Sim6.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX015_Xia_Seq_Sim7.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX016_Xia_Exp1_Seq_Full_4.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX018_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX018_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX018_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_PSTH_SinglePulse_SimSeq/PSTHLinearity_DX018_Xia_Exp1_Sim4.mat';

    
    };

% Empty = use the union of all non-zero amplitudes found in the files.
Target_Amplitudes_uA = [];
Include_Zero_Amplitude = false;
Amplitude_Tolerance_uA = 1e-6;

% -------------------------------------------------------------------------
% Optional exclusions. Each row is:
%   Animal ID, Dataset ID, Amplitude (uA), Branch label
%
% Use '*' to match every value in a text column. Branch examples are
% 'A->B matched' and 'B->A matched'. Exclusion is amplitude-specific.
%
% Examples:
% ExcludeRules = {
%     'DX015', 'Xia_Seq_Sim1', 6,   '*';
%     'DX014', 'Xia_Seq_Sim4', 10,  'A->B matched';
% };
% -------------------------------------------------------------------------
ExcludeRules = {
    % AnimalID, DatasetID, Amplitude_uA, BranchLabel
};

% Plot smoothing only. Saved numerical window summaries remain unchanged.
Plot_Smooth_Sigma_ms = 1.5;
Plot_Window_ms = [-10 30];

% Display location of the second pulse. This is deliberately independent
% of saved PTD metadata because some experiments label the same practical
% timing as either 5 or 5.5 ms. Set [] to hide the second-pulse marker.
Second_Pulse_Marker_ms = 5;

% Normalized temporal-shape analysis.
Plot_Normalized_Results = true;
Min_Normalization_Scale_sp_s = 5;

% Display options.
% Keep this false for a lighter, faster figure. The window-summary figure
% always includes the individual branch points.
Plot_Individual_Branch_Traces = false;
Individual_Trace_Max_Number = inf;  % inf = show every contributing branch
Individual_Trace_Color_Strength = 0.18;
Max_Amplitudes_Per_Figure = 4;
Figure_Position = [60 60 1550 900];

% This prevents accidentally counting two saved versions of the same
% dataset. Set true only if duplicate files are intentional.
Allow_Duplicate_Dataset_Files = false;

% Optional saving of the compact cross-dataset result.
Save_CrossDataset_Result = false;
CrossDataset_Save_File = fullfile(Results_Folder, ...
    'CrossDataset_PSTHLinearity_BranchAverage_v1.mat');

%% ============================== CHECKS ==================================
if ~Use_Manual_Result_Files && ~isfolder(Results_Folder)
    error('CrossPSTH:FolderNotFound', ...
        'Results folder was not found:\n%s',Results_Folder);
end
if Plot_Smooth_Sigma_ms < 0
    error('CrossPSTH:InvalidSmoothing','Plot_Smooth_Sigma_ms must be >= 0.');
end
if Min_Normalization_Scale_sp_s <= 0
    error('CrossPSTH:InvalidNormalizationThreshold', ...
        'Min_Normalization_Scale_sp_s must be greater than zero.');
end
if Max_Amplitudes_Per_Figure < 1
    error('CrossPSTH:InvalidPageSize', ...
        'Max_Amplitudes_Per_Figure must be at least 1.');
end
if ~isempty(Second_Pulse_Marker_ms) && ...
        (~isscalar(Second_Pulse_Marker_ms) || ...
         ~isfinite(Second_Pulse_Marker_ms))
    error('CrossPSTH:InvalidSecondPulseMarker', ...
        'Second_Pulse_Marker_ms must be empty or one finite number.');
end
if Save_CrossDataset_Result
    save_parent = fileparts(CrossDataset_Save_File);
    if ~isfolder(save_parent)
        error('CrossPSTH:SaveFolderNotFound', ...
            'Cross-dataset save folder was not found:\n%s',save_parent);
    end
end

if Use_Manual_Result_Files
    result_file_paths = string(Manual_Result_Files(:));
    result_file_paths = result_file_paths(strlength(strtrim(result_file_paths))>0);
    if isempty(result_file_paths)
        error('CrossPSTH:NoManualFiles', ...
            ['Use_Manual_Result_Files is true, but no paths were entered ' ...
             'in Manual_Result_Files.']);
    end
    for k = 1:numel(result_file_paths)
        if ~isfile(result_file_paths(k))
            error('CrossPSTH:ManualFileNotFound', ...
                'Manual result file was not found:\n%s',result_file_paths(k));
        end
    end
else
    files = dir(fullfile(Results_Folder,File_Pattern));
    files = files(~startsWith({files.name},'._'));
    if isempty(files)
        error('CrossPSTH:NoFiles', ...
            'No files matching %s were found in %s.', ...
            File_Pattern,Results_Folder);
    end
    result_file_paths = strings(numel(files),1);
    for k = 1:numel(files)
        result_file_paths(k) = fullfile(files(k).folder,files(k).name);
    end
end

%% ======================== LOAD BRANCH OBSERVATIONS ======================
Observations = repmat(empty_observation(),0,1);
reference_time_ms = [];
reference_settings = struct();
loaded_file_paths = strings(0,1);
dataset_keys = strings(0,1);
saved_ptd_values_ms = nan(0,1);

for ifile = 1:numel(result_file_paths)
    file_path = char(result_file_paths(ifile));
    X = load(file_path,'PSTHLinearitySaved');
    if ~isfield(X,'PSTHLinearitySaved')
        error('CrossPSTH:MissingSavedVariable', ...
            'PSTHLinearitySaved is missing from:\n%s',file_path);
    end
    R = X.PSTHLinearitySaved;
    validate_saved_result(R,file_path);
    saved_ptd_values_ms(end+1,1) = ...
        double(R.settings.sequential_ptd_ms); %#ok<SAGROW>

    animal_id = char(string(R.animal_id));
    dataset_id = char(string(R.dataset_id));
    single_bad_trial_file = '';
    if isfield(R,'source_files') && isfield(R.source_files,'single') && ...
            isfield(R.source_files.single,'bad_trial_file')
        single_bad_trial_file = char(string( ...
            R.source_files.single.bad_trial_file));
    end
    has_single_bad_trial_file = ~isempty(strtrim(single_bad_trial_file));
    this_dataset_key = string(animal_id) + " | " + string(dataset_id);
    if ~Allow_Duplicate_Dataset_Files && any(dataset_keys == this_dataset_key)
        error('CrossPSTH:DuplicateDataset', ...
            ['More than one result file was found for %s. Remove the old ' ...
             'version or set Allow_Duplicate_Dataset_Files = true.'], ...
            this_dataset_key);
    end
    dataset_keys(end+1,1) = this_dataset_key; %#ok<SAGROW>
    loaded_file_paths(end+1,1) = string(file_path); %#ok<SAGROW>

    if isempty(reference_time_ms)
        reference_time_ms = double(R.branch_averages(1).time_ms(:).');
        reference_settings = R.settings;
    else
        assert_compatible_settings(reference_settings,R.settings,file_path);
    end

    for ib = 1:numel(R.branch_averages)
        B = R.branch_averages(ib);
        this_time_ms = double(B.time_ms(:).');
        assert_matching_time_grid(reference_time_ms,this_time_ms,file_path,ib);

        amps = double(B.amplitudes_uA(:).');
        nAmp = numel(amps);
        for ia = 1:nAmp
            amp = amps(ia);
            if ~Include_Zero_Amplitude && abs(amp) <= Amplitude_Tolerance_uA
                continue;
            end
            if ~isempty(Target_Amplitudes_uA) && ...
                    ~any(abs(Target_Amplitudes_uA - amp) <= Amplitude_Tolerance_uA)
                continue;
            end
            if is_excluded(animal_id,dataset_id,amp,B.branch_label, ...
                    ExcludeRules,Amplitude_Tolerance_uA)
                continue;
            end

            O = empty_observation();
            O.file_path = file_path;
            O.animal_id = animal_id;
            O.dataset_id = dataset_id;
            O.dataset_key = char(this_dataset_key);
            O.single_bad_trial_file = single_bad_trial_file;
            O.has_single_bad_trial_file = has_single_bad_trial_file;
            O.pair_key = char(string(B.pair_key));
            O.electrode_A = char(string(B.electrode_A));
            O.electrode_B = char(string(B.electrode_B));
            O.branch_code = char(string(B.branch_code));
            O.branch_label = char(string(B.branch_label));
            O.amplitude_uA = amp;
            O.n_responding_channels = double(B.n_responding_channels_detected);
            O.n_channels_used = double(B.n_channels_used_by_amplitude(ia));
            O.time_ms = this_time_ms;

            O.sim_observed = row_curve(B.curves.sim_observed.mean,ia);
            O.sim_linear = row_curve(B.curves.sim_linear.mean,ia);
            O.sim_residual = row_curve(B.curves.sim_residual.mean,ia);
            O.seq_observed = row_curve(B.curves.seq_observed.mean,ia);
            O.seq_linear = row_curve(B.curves.seq_linear.mean,ia);
            O.seq_residual = row_curve(B.curves.seq_residual.mean,ia);

            O.full_difference = double(B.full_difference.mean(ia,:));
            O.early_difference = double(B.early_difference.mean(ia,:));
            O.late_difference = double(B.late_difference.mean(ia,:));

            required_curves = [O.sim_observed; O.sim_linear; O.sim_residual; ...
                O.seq_observed; O.seq_linear; O.seq_residual];
            O.is_absolute_valid = O.n_channels_used > 0 && ...
                all(isfinite(required_curves),'all') && ...
                all(isfinite([O.full_difference O.early_difference ...
                              O.late_difference]));

            if O.is_absolute_valid
                sim_linear_s = smooth_curve(O.sim_linear,O.time_ms, ...
                    Plot_Smooth_Sigma_ms);
                seq_linear_s = smooth_curve(O.seq_linear,O.time_ms, ...
                    Plot_Smooth_Sigma_ms);
                response_mask = O.time_ms >= R.settings.response_win_ms(1) & ...
                    O.time_ms < R.settings.response_win_ms(2);
                scale_values = abs([sim_linear_s(response_mask), ...
                                    seq_linear_s(response_mask)]);
                O.normalization_scale_sp_s = max(scale_values,[],'omitnan');
                O.is_normalized_valid = isfinite(O.normalization_scale_sp_s) && ...
                    O.normalization_scale_sp_s >= Min_Normalization_Scale_sp_s;
            end

            Observations(end+1,1) = O; %#ok<SAGROW>
        end
    end
end

if isempty(Observations)
    error('CrossPSTH:NoObservations', ...
        'No branch observations remained after amplitude selection/exclusion.');
end

valid_absolute = [Observations.is_absolute_valid].';
if ~any(valid_absolute)
    error('CrossPSTH:NoValidObservations', ...
        'No branch observation contained a complete valid result.');
end

Observations = Observations(valid_absolute);
all_amps = unique([Observations.amplitude_uA]);
if ~isempty(Target_Amplitudes_uA)
    target_order = Target_Amplitudes_uA(:).';
    Amplitudes_uA = target_order(arrayfun(@(a) ...
        any(abs(all_amps-a)<=Amplitude_Tolerance_uA),target_order));
else
    Amplitudes_uA = sort(all_amps);
end

%% ======================= BUILD CROSS-DATASET RESULT =====================
nAmp = numel(Amplitudes_uA);
nTime = numel(reference_time_ms);
CrossDatasetPSTH = struct();
CrossDatasetPSTH.analysis_name = ...
    'Cross-dataset PSTH temporal linearity: direct branch average';
CrossDatasetPSTH.created_at = datestr(now,31);
CrossDatasetPSTH.results_folder = Results_Folder;
CrossDatasetPSTH.loaded_files = loaded_file_paths;
CrossDatasetPSTH.time_ms = reference_time_ms;
CrossDatasetPSTH.amplitudes_uA = Amplitudes_uA;
CrossDatasetPSTH.averaging_unit = ...
    'One order-matched branch; branches are not merged or channel-weighted';
CrossDatasetPSTH.normalization_description = ...
    ['Per Branch x Amplitude maximum absolute smoothed Sim/Seq linear ' ...
     'prediction within the response window'];
CrossDatasetPSTH.min_normalization_scale_sp_s = ...
    Min_Normalization_Scale_sp_s;
CrossDatasetPSTH.reference_settings = reference_settings;
CrossDatasetPSTH.saved_sequential_ptd_values_ms = ...
    unique(saved_ptd_values_ms(isfinite(saved_ptd_values_ms))).';
CrossDatasetPSTH.second_pulse_plot_marker_ms = Second_Pulse_Marker_ms;
CrossDatasetPSTH.observations = Observations;

curve_names = {'sim_observed','sim_linear','sim_residual', ...
    'seq_observed','seq_linear','seq_residual'};
for k = 1:numel(curve_names)
    name = curve_names{k};
    CrossDatasetPSTH.absolute.(name).mean = nan(nAmp,nTime);
    CrossDatasetPSTH.absolute.(name).sem = nan(nAmp,nTime);
    CrossDatasetPSTH.normalized.(name).mean = nan(nAmp,nTime);
    CrossDatasetPSTH.normalized.(name).sem = nan(nAmp,nTime);
end

CrossDatasetPSTH.absolute.full.mean = nan(nAmp,2);
CrossDatasetPSTH.absolute.full.sem = nan(nAmp,2);
CrossDatasetPSTH.absolute.early.mean = nan(nAmp,2);
CrossDatasetPSTH.absolute.early.sem = nan(nAmp,2);
CrossDatasetPSTH.absolute.late.mean = nan(nAmp,2);
CrossDatasetPSTH.absolute.late.sem = nan(nAmp,2);
CrossDatasetPSTH.branch_indices_by_amplitude = cell(nAmp,1);
CrossDatasetPSTH.normalized_branch_indices_by_amplitude = cell(nAmp,1);

NBranches = zeros(nAmp,1);
NDatasets = zeros(nAmp,1);
NAnimals = zeros(nAmp,1);
NNormalized = zeros(nAmp,1);
NNormalizedDatasets = zeros(nAmp,1);
TotalBranchChannels = zeros(nAmp,1);
MedianChannelsPerBranch = nan(nAmp,1);

for ia = 1:nAmp
    amp = Amplitudes_uA(ia);
    idx = find(abs([Observations.amplitude_uA] - amp) <= ...
        Amplitude_Tolerance_uA);
    idx_norm = idx([Observations(idx).is_normalized_valid]);
    CrossDatasetPSTH.branch_indices_by_amplitude{ia} = idx;
    CrossDatasetPSTH.normalized_branch_indices_by_amplitude{ia} = idx_norm;

    NBranches(ia) = numel(idx);
    NDatasets(ia) = numel(unique(string({Observations(idx).dataset_key})));
    NAnimals(ia) = numel(unique(string({Observations(idx).animal_id})));
    NNormalized(ia) = numel(idx_norm);
    if ~isempty(idx_norm)
        NNormalizedDatasets(ia) = numel(unique( ...
            string({Observations(idx_norm).dataset_key})));
    end
    channel_counts = [Observations(idx).n_channels_used];
    TotalBranchChannels(ia) = sum(channel_counts);
    MedianChannelsPerBranch(ia) = median(channel_counts);

    for k = 1:numel(curve_names)
        name = curve_names{k};
        M = vertcat(Observations(idx).(name));
        [mu,se] = mean_and_sem(M);
        CrossDatasetPSTH.absolute.(name).mean(ia,:) = mu;
        CrossDatasetPSTH.absolute.(name).sem(ia,:) = se;

        if ~isempty(idx_norm)
            Mnorm = nan(numel(idx_norm),nTime);
            for j = 1:numel(idx_norm)
                io = idx_norm(j);
                Mnorm(j,:) = Observations(io).(name) ./ ...
                    Observations(io).normalization_scale_sp_s;
            end
            [mu_norm,se_norm] = mean_and_sem(Mnorm);
            CrossDatasetPSTH.normalized.(name).mean(ia,:) = mu_norm;
            CrossDatasetPSTH.normalized.(name).sem(ia,:) = se_norm;
        end
    end

    [CrossDatasetPSTH.absolute.full.mean(ia,:), ...
        CrossDatasetPSTH.absolute.full.sem(ia,:)] = ...
        mean_and_sem(vertcat(Observations(idx).full_difference));
    [CrossDatasetPSTH.absolute.early.mean(ia,:), ...
        CrossDatasetPSTH.absolute.early.sem(ia,:)] = ...
        mean_and_sem(vertcat(Observations(idx).early_difference));
    [CrossDatasetPSTH.absolute.late.mean(ia,:), ...
        CrossDatasetPSTH.absolute.late.sem(ia,:)] = ...
        mean_and_sem(vertcat(Observations(idx).late_difference));
end

Inventory = table(Amplitudes_uA(:),NBranches,NDatasets,NAnimals, ...
    NNormalized,NNormalizedDatasets,TotalBranchChannels, ...
    MedianChannelsPerBranch, ...
    'VariableNames',{'Amplitude_uA','NBranches','NDatasets','NAnimals', ...
    'NNormalizedBranches','NNormalizedDatasets','TotalBranchChannels', ...
    'MedianChannelsPerBranch'});
CrossDatasetPSTH.inventory = Inventory;

%% =========================== SIMPLE PRINTOUT =============================
fprintf('\nCross-dataset PSTH temporal linearity\n');
fprintf('Files: %d | Datasets: %d | Animals: %d | Branch observations: %d\n', ...
    numel(result_file_paths),numel(unique(string({Observations.dataset_key}))), ...
    numel(unique(string({Observations.animal_id}))),numel(Observations));
fprintf('Averaging: all available branches directly; no channel weighting\n');
fprintf('Saved PTD labels (ms): ');
fprintf('%.3g ',CrossDatasetPSTH.saved_sequential_ptd_values_ms);
fprintf('| plot marker: ');
if isempty(Second_Pulse_Marker_ms)
    fprintf('hidden\n');
else
    fprintf('%.3g ms\n',Second_Pulse_Marker_ms);
end
unique_dataset_keys = unique(string({Observations.dataset_key}));
single_qc_dataset_keys = unique(string( ...
    {Observations([Observations.has_single_bad_trial_file]).dataset_key}));
fprintf('Single-stimulation BadTrials files: %d/%d datasets\n', ...
    numel(single_qc_dataset_keys),numel(unique_dataset_keys));
fprintf('Normalization threshold: %.3g spikes/s\n\n', ...
    Min_Normalization_Scale_sp_s);
disp(Inventory);

WindowSummary = table(Amplitudes_uA(:),NBranches, ...
    CrossDatasetPSTH.absolute.full.mean(:,1), ...
    CrossDatasetPSTH.absolute.full.mean(:,2), ...
    CrossDatasetPSTH.absolute.early.mean(:,1), ...
    CrossDatasetPSTH.absolute.early.mean(:,2), ...
    CrossDatasetPSTH.absolute.late.mean(:,1), ...
    CrossDatasetPSTH.absolute.late.mean(:,2), ...
    'VariableNames',{'Amplitude_uA','NBranches','SimFull','SeqFull', ...
    'SimEarly','SeqEarly','SimLate','SeqLate'});
CrossDatasetPSTH.window_summary = WindowSummary;
fprintf('\nMean observed-minus-linear result (spikes/trial)\n');
disp(WindowSummary);

%% =============================== PLOTS ==================================
time_ms = CrossDatasetPSTH.time_ms;
display_mask = time_ms >= Plot_Window_ms(1) & time_ms <= Plot_Window_ms(2);
response_win = reference_settings.response_win_ms;
ptd_ms = Second_Pulse_Marker_ms;

% Smooth every branch first, then calculate the displayed mean and SEM.
% This gives the correct variability of the smoothed branch curves.
AbsPlot = build_plot_result(Observations, ...
    CrossDatasetPSTH.branch_indices_by_amplitude,curve_names,time_ms, ...
    Plot_Smooth_Sigma_ms,false);
NormPlot = build_plot_result(Observations, ...
    CrossDatasetPSTH.normalized_branch_indices_by_amplitude,curve_names, ...
    time_ms,Plot_Smooth_Sigma_ms,true);

sim_col = [0.00 0.35 0.85];
seq_col = [0.90 0.25 0.05];
linear_col = [0.10 0.10 0.10];
zero_col = [0.45 0.45 0.45];

nPages = ceil(nAmp/Max_Amplitudes_Per_Figure);
for ipage = 1:nPages
    row_start = (ipage-1)*Max_Amplitudes_Per_Figure + 1;
    row_end = min(ipage*Max_Amplitudes_Per_Figure,nAmp);
    page_amp_indices = row_start:row_end;
    nRows = numel(page_amp_indices);

    % -------- Figure 1: absolute observed versus prediction --------
    figure('Color','w','Position',Figure_Position, ...
        'Name',sprintf('Absolute PSTH page %d',ipage));
    tiledlayout(nRows,2,'TileSpacing','compact','Padding','compact');
    sgtitle(sprintf(['Cross-dataset absolute PSTHs | branch-weighted ' ...
        '(page %d/%d)'],ipage,nPages),'FontWeight','bold');

    for ir = 1:nRows
        ia = page_amp_indices(ir);
        amp = Amplitudes_uA(ia);
        idx = CrossDatasetPSTH.branch_indices_by_amplitude{ia};

        ax = nexttile; hold(ax,'on');
        if Plot_Individual_Branch_Traces
            plot_individual_condition(ax,Observations,idx,'sim_observed', ...
                'sim_linear',time_ms,display_mask,Plot_Smooth_Sigma_ms, ...
                sim_col,linear_col,Individual_Trace_Color_Strength, ...
                Individual_Trace_Max_Number,false);
        end
        shaded_line(ax,time_ms(display_mask), ...
            AbsPlot.sim_linear.mean(ia,display_mask), ...
            AbsPlot.sim_linear.sem(ia,display_mask),linear_col,'--',1.8);
        shaded_line(ax,time_ms(display_mask), ...
            AbsPlot.sim_observed.mean(ia,display_mask), ...
            AbsPlot.sim_observed.sem(ia,display_mask),sim_col,'-',2.0);
        decorate_time_axis(ax,Plot_Window_ms,ptd_ms);
        title(ax,sprintf('%.3g uA | Simultaneous | n=%d',amp,NBranches(ia)));
        ylabel(ax,'Baseline-corrected rate (spikes/s)');
        if ir == nRows, xlabel(ax,'Time from first pulse (ms)'); end
        if ir == 1
            legend(ax,{'Linear prediction','Observed'}, ...
                'Location','best','Box','off');
        end

        ax = nexttile; hold(ax,'on');
        if Plot_Individual_Branch_Traces
            plot_individual_condition(ax,Observations,idx,'seq_observed', ...
                'seq_linear',time_ms,display_mask,Plot_Smooth_Sigma_ms, ...
                seq_col,linear_col,Individual_Trace_Color_Strength, ...
                Individual_Trace_Max_Number,false);
        end
        shaded_line(ax,time_ms(display_mask), ...
            AbsPlot.seq_linear.mean(ia,display_mask), ...
            AbsPlot.seq_linear.sem(ia,display_mask),linear_col,'--',1.8);
        shaded_line(ax,time_ms(display_mask), ...
            AbsPlot.seq_observed.mean(ia,display_mask), ...
            AbsPlot.seq_observed.sem(ia,display_mask),seq_col,'-',2.0);
        decorate_time_axis(ax,Plot_Window_ms,ptd_ms);
        title(ax,sprintf('%.3g uA | Sequential | n=%d',amp,NBranches(ia)));
        if ir == nRows, xlabel(ax,'Time from first pulse (ms)'); end
        if ir == 1
            legend(ax,{'Linear prediction','Observed'}, ...
                'Location','best','Box','off');
        end
    end
    apply_common_ylim(gcf);

    % -------- Figure 2: absolute and normalized residual PSTHs --------
    figure('Color','w','Position',Figure_Position, ...
        'Name',sprintf('Residual PSTH page %d',ipage));
    tiledlayout(nRows,2,'TileSpacing','compact','Padding','compact');
    sgtitle(sprintf(['Cross-dataset temporal residuals | branch-weighted ' ...
        '(page %d/%d)'],ipage,nPages),'FontWeight','bold');

    for ir = 1:nRows
        ia = page_amp_indices(ir);
        amp = Amplitudes_uA(ia);
        idx = CrossDatasetPSTH.branch_indices_by_amplitude{ia};
        idx_norm = CrossDatasetPSTH.normalized_branch_indices_by_amplitude{ia};

        ax = nexttile; hold(ax,'on');
        if Plot_Individual_Branch_Traces
            plot_individual_residuals(ax,Observations,idx,time_ms, ...
                display_mask,Plot_Smooth_Sigma_ms,sim_col,seq_col, ...
                Individual_Trace_Color_Strength, ...
                Individual_Trace_Max_Number,false);
        end
        shaded_line(ax,time_ms(display_mask), ...
            AbsPlot.sim_residual.mean(ia,display_mask), ...
            AbsPlot.sim_residual.sem(ia,display_mask),sim_col,'-',2.0);
        shaded_line(ax,time_ms(display_mask), ...
            AbsPlot.seq_residual.mean(ia,display_mask), ...
            AbsPlot.seq_residual.sem(ia,display_mask),seq_col,'-',2.0);
        yline(ax,0,'--','Color',zero_col,'HandleVisibility','off');
        decorate_time_axis(ax,Plot_Window_ms,ptd_ms);
        title(ax,sprintf('%.3g uA | absolute | n=%d',amp,numel(idx)));
        ylabel(ax,'Observed - linear (spikes/s)');
        if ir == nRows, xlabel(ax,'Time from first pulse (ms)'); end
        if ir == 1
            legend(ax,{'Sim - linear','Seq - linear'}, ...
                'Location','best','Box','off');
        end

        ax = nexttile; hold(ax,'on');
        if Plot_Normalized_Results && ~isempty(idx_norm)
            if Plot_Individual_Branch_Traces
                plot_individual_residuals(ax,Observations,idx_norm,time_ms, ...
                    display_mask,Plot_Smooth_Sigma_ms,sim_col,seq_col, ...
                    Individual_Trace_Color_Strength, ...
                    Individual_Trace_Max_Number,true);
            end
            shaded_line(ax,time_ms(display_mask), ...
                NormPlot.sim_residual.mean(ia,display_mask), ...
                NormPlot.sim_residual.sem(ia,display_mask),sim_col,'-',2.0);
            shaded_line(ax,time_ms(display_mask), ...
                NormPlot.seq_residual.mean(ia,display_mask), ...
                NormPlot.seq_residual.sem(ia,display_mask),seq_col,'-',2.0);
        else
            text(ax,0.5,0.5,'Normalized plot disabled or no valid branches', ...
                'Units','normalized','HorizontalAlignment','center');
        end
        yline(ax,0,'--','Color',zero_col,'HandleVisibility','off');
        decorate_time_axis(ax,Plot_Window_ms,ptd_ms);
        title(ax,sprintf('%.3g uA | normalized | n=%d',amp,numel(idx_norm)));
        ylabel(ax,'Normalized observed - linear');
        if ir == nRows, xlabel(ax,'Time from first pulse (ms)'); end
        if ir == 1 && Plot_Normalized_Results && ~isempty(idx_norm)
            legend(ax,{'Sim - linear','Seq - linear'}, ...
                'Location','best','Box','off');
        end
    end
    apply_common_ylim_by_column(gcf,2);

    % -------- Figure 3: optional normalized observed/predicted PSTHs -----
    if Plot_Normalized_Results
        figure('Color','w','Position',Figure_Position, ...
            'Name',sprintf('Normalized PSTH page %d',ipage));
        tiledlayout(nRows,2,'TileSpacing','compact','Padding','compact');
        sgtitle(sprintf(['Cross-dataset normalized PSTHs | common Sim/Seq ' ...
            'scale (page %d/%d)'],ipage,nPages),'FontWeight','bold');

        for ir = 1:nRows
            ia = page_amp_indices(ir);
            amp = Amplitudes_uA(ia);
            idx_norm = ...
                CrossDatasetPSTH.normalized_branch_indices_by_amplitude{ia};

            ax = nexttile; hold(ax,'on');
            if ~isempty(idx_norm)
                if Plot_Individual_Branch_Traces
                    plot_individual_condition(ax,Observations,idx_norm, ...
                        'sim_observed','sim_linear',time_ms,display_mask, ...
                        Plot_Smooth_Sigma_ms,sim_col,linear_col, ...
                        Individual_Trace_Color_Strength, ...
                        Individual_Trace_Max_Number,true);
                end
                shaded_line(ax,time_ms(display_mask), ...
                    NormPlot.sim_linear.mean(ia,display_mask), ...
                    NormPlot.sim_linear.sem(ia,display_mask), ...
                    linear_col,'--',1.8);
                shaded_line(ax,time_ms(display_mask), ...
                    NormPlot.sim_observed.mean(ia,display_mask), ...
                    NormPlot.sim_observed.sem(ia,display_mask), ...
                    sim_col,'-',2.0);
            end
            decorate_time_axis(ax,Plot_Window_ms,ptd_ms);
            title(ax,sprintf('%.3g uA | Sim | normalized n=%d', ...
                amp,numel(idx_norm)));
            ylabel(ax,'Normalized rate');
            if ir == nRows, xlabel(ax,'Time from first pulse (ms)'); end

            ax = nexttile; hold(ax,'on');
            if ~isempty(idx_norm)
                if Plot_Individual_Branch_Traces
                    plot_individual_condition(ax,Observations,idx_norm, ...
                        'seq_observed','seq_linear',time_ms,display_mask, ...
                        Plot_Smooth_Sigma_ms,seq_col,linear_col, ...
                        Individual_Trace_Color_Strength, ...
                        Individual_Trace_Max_Number,true);
                end
                shaded_line(ax,time_ms(display_mask), ...
                    NormPlot.seq_linear.mean(ia,display_mask), ...
                    NormPlot.seq_linear.sem(ia,display_mask), ...
                    linear_col,'--',1.8);
                shaded_line(ax,time_ms(display_mask), ...
                    NormPlot.seq_observed.mean(ia,display_mask), ...
                    NormPlot.seq_observed.sem(ia,display_mask), ...
                    seq_col,'-',2.0);
            end
            decorate_time_axis(ax,Plot_Window_ms,ptd_ms);
            title(ax,sprintf('%.3g uA | Seq | normalized n=%d', ...
                amp,numel(idx_norm)));
            if ir == nRows, xlabel(ax,'Time from first pulse (ms)'); end
        end
        apply_common_ylim(gcf);
    end
end

% -------- Figure 4: exact unsmoothed window summaries --------------------
window_names = {'Full','Early','Late'};
window_fields = {'full','early','late'};
figure('Color','w','Position',[80 100 1550 520], ...
    'Name','Window residual summary');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
sgtitle('Cross-dataset observed-minus-linear summary | branch-weighted', ...
    'FontWeight','bold');

for iw = 1:3
    field_name = window_fields{iw};
    ax = nexttile; hold(ax,'on');
    for ia = 1:nAmp
        idx = CrossDatasetPSTH.branch_indices_by_amplitude{ia};
        values = observation_matrix(Observations,idx, ...
            [field_name '_difference']);
        jitter = deterministic_jitter(size(values,1),0.10);
        scatter(ax,Amplitudes_uA(ia)-0.025+jitter,values(:,1),18, ...
            lighten_color(sim_col,0.72),'filled','HandleVisibility','off');
        scatter(ax,Amplitudes_uA(ia)+0.025+jitter,values(:,2),18, ...
            lighten_color(seq_col,0.72),'filled','HandleVisibility','off');
    end
    sim_sem_plot = sem_for_plot( ...
        CrossDatasetPSTH.absolute.(field_name).mean(:,1), ...
        CrossDatasetPSTH.absolute.(field_name).sem(:,1));
    seq_sem_plot = sem_for_plot( ...
        CrossDatasetPSTH.absolute.(field_name).mean(:,2), ...
        CrossDatasetPSTH.absolute.(field_name).sem(:,2));
    errorbar(ax,Amplitudes_uA, ...
        CrossDatasetPSTH.absolute.(field_name).mean(:,1), ...
        sim_sem_plot,'-o', ...
        'Color',sim_col,'MarkerFaceColor',sim_col,'LineWidth',1.8, ...
        'CapSize',8,'DisplayName','Sim - linear');
    errorbar(ax,Amplitudes_uA, ...
        CrossDatasetPSTH.absolute.(field_name).mean(:,2), ...
        seq_sem_plot,'-o', ...
        'Color',seq_col,'MarkerFaceColor',seq_col,'LineWidth',1.8, ...
        'CapSize',8,'DisplayName','Seq - linear');
    yline(ax,0,'--','Color',zero_col,'HandleVisibility','off');
    xlabel(ax,'Amplitude (uA)');
    ylabel(ax,'Observed - linear (spikes/trial)');
    title(ax,sprintf('%s window',window_names{iw}));
    box(ax,'off');
    if iw == 1
        legend(ax,'Location','best','Box','off');
    end
end
apply_common_ylim(gcf);

%% ============================= OPTIONAL SAVE =============================
if Save_CrossDataset_Result
    save(CrossDataset_Save_File,'CrossDatasetPSTH','-v7.3');
    fprintf('\nSaved cross-dataset result:\n%s\n',CrossDataset_Save_File);
end

%% ============================ LOCAL FUNCTIONS ============================
function O = empty_observation()
O = struct( ...
    'file_path','', ...
    'animal_id','', ...
    'dataset_id','', ...
    'dataset_key','', ...
    'single_bad_trial_file','', ...
    'has_single_bad_trial_file',false, ...
    'pair_key','', ...
    'electrode_A','', ...
    'electrode_B','', ...
    'branch_code','', ...
    'branch_label','', ...
    'amplitude_uA',NaN, ...
    'n_responding_channels',NaN, ...
    'n_channels_used',NaN, ...
    'time_ms',[], ...
    'sim_observed',[], ...
    'sim_linear',[], ...
    'sim_residual',[], ...
    'seq_observed',[], ...
    'seq_linear',[], ...
    'seq_residual',[], ...
    'full_difference',[NaN NaN], ...
    'early_difference',[NaN NaN], ...
    'late_difference',[NaN NaN], ...
    'normalization_scale_sp_s',NaN, ...
    'is_absolute_valid',false, ...
    'is_normalized_valid',false);
end

function validate_saved_result(R,file_path)
required_top = {'animal_id','dataset_id','settings','branch_averages', ...
    'n_saved_branches'};
for k = 1:numel(required_top)
    if ~isfield(R,required_top{k})
        error('CrossPSTH:InvalidSavedResult', ...
            '%s is missing from %s.',required_top{k},file_path);
    end
end
if isempty(R.branch_averages)
    error('CrossPSTH:NoSavedBranches','No branch average exists in %s.',file_path);
end
required_settings = {'response_win_ms','early_win_ms','late_win_ms', ...
    'sequential_ptd_ms','bin_ms'};
for k = 1:numel(required_settings)
    if ~isfield(R.settings,required_settings{k})
        error('CrossPSTH:InvalidSettings', ...
            'settings.%s is missing from %s.',required_settings{k},file_path);
    end
end
required_branch = {'amplitudes_uA','time_ms','branch_label','branch_code', ...
    'pair_key','n_responding_channels_detected', ...
    'n_channels_used_by_amplitude','curves','full_difference', ...
    'early_difference','late_difference'};
for k = 1:numel(required_branch)
    if ~isfield(R.branch_averages(1),required_branch{k})
        error('CrossPSTH:InvalidBranch', ...
            'branch_averages.%s is missing from %s.', ...
            required_branch{k},file_path);
    end
end
end

function assert_matching_time_grid(reference_time,this_time,file_path,branch_idx)
if numel(reference_time) ~= numel(this_time) || ...
        any(abs(reference_time-this_time) > 1e-9)
    error('CrossPSTH:TimeGridMismatch', ...
        ['The time grid in branch %d of %s differs from previous files. ' ...
         'Cross-dataset interpolation is intentionally not applied.'], ...
        branch_idx,file_path);
end
end

function assert_compatible_settings(reference_settings,this_settings,file_path)
% PTD is intentionally not included here. Some datasets encode the same
% practical sequential timing as 5 or 5.5 ms. The stored curves are used as
% calculated, and the display marker is controlled by the user setting.
numeric_fields = {'response_win_ms','early_win_ms','late_win_ms','bin_ms'};
for k = 1:numel(numeric_fields)
    field_name = numeric_fields{k};
    a = double(reference_settings.(field_name));
    b = double(this_settings.(field_name));
    if ~isequal(size(a),size(b)) || any(abs(a(:)-b(:)) > 1e-9)
        error('CrossPSTH:IncompatibleSettings', ...
            ['settings.%s in %s differs from the first result file. ' ...
             'Results with different analysis windows or bins must not be ' ...
             'averaged together.'],field_name,file_path);
    end
end
end

function y = row_curve(M,row_idx)
y = double(M(row_idx,:));
y = y(:).';
end

function tf = is_excluded(animal_id,dataset_id,amp,branch_label,rules,tol)
tf = false;
if isempty(rules), return; end
if size(rules,2) ~= 4
    error('CrossPSTH:InvalidExcludeRules', ...
        'ExcludeRules must contain four columns.');
end
for k = 1:size(rules,1)
    animal_match = text_rule_match(animal_id,rules{k,1});
    dataset_match = text_rule_match(dataset_id,rules{k,2});
    amp_rule = double(rules{k,3});
    amp_match = isfinite(amp_rule) && abs(amp-amp_rule) <= tol;
    branch_match = text_rule_match(branch_label,rules{k,4});
    if animal_match && dataset_match && amp_match && branch_match
        tf = true;
        return;
    end
end
end

function tf = text_rule_match(value,rule)
rule = char(string(rule));
tf = strcmp(rule,'*') || strcmpi(char(string(value)),rule);
end

function y = smooth_curve(x,time_ms,sigma_ms)
x = double(x(:).');
if sigma_ms <= 0 || numel(x) < 2
    y = x;
    return;
end
dt = median(diff(time_ms));
sigma_bins = sigma_ms/dt;
half_width = max(1,ceil(4*sigma_bins));
q = -half_width:half_width;
g = exp(-0.5*(q/sigma_bins).^2);
g = g/sum(g);
y = conv(x,g,'same');
end

function [mu,se] = mean_and_sem(M)
if isempty(M)
    mu = NaN(1,size(M,2));
    se = NaN(1,size(M,2));
    return;
end
mu = mean(M,1,'omitnan');
n = sum(isfinite(M),1);
se = nan(size(mu));
valid = n > 1;
sd = std(M,0,1,'omitnan');
se(valid) = sd(valid)./sqrt(n(valid));
end

function se_plot = sem_for_plot(mu,se)
se_plot = se;
se_plot(~isfinite(se_plot) & isfinite(mu)) = 0;
end

function Out = build_plot_result(Obs,indices_by_amp,curve_names,time_ms, ...
        sigma_ms,normalize)
nAmp = numel(indices_by_amp);
nTime = numel(time_ms);
Out = struct();
for k = 1:numel(curve_names)
    name = curve_names{k};
    Out.(name).mean = nan(nAmp,nTime);
    Out.(name).sem = nan(nAmp,nTime);
    for ia = 1:nAmp
        idx = indices_by_amp{ia};
        if isempty(idx), continue; end
        M = nan(numel(idx),nTime);
        for j = 1:numel(idx)
            io = idx(j);
            curve = Obs(io).(name);
            if normalize
                curve = curve/Obs(io).normalization_scale_sp_s;
            end
            M(j,:) = smooth_curve(curve,time_ms,sigma_ms);
        end
        [Out.(name).mean(ia,:),Out.(name).sem(ia,:)] = mean_and_sem(M);
    end
end
end

function M = observation_matrix(Obs,idx,field_name)
if isempty(idx)
    M = zeros(0,0);
    return;
end
first_value = Obs(idx(1)).(field_name);
M = nan(numel(idx),numel(first_value));
for k = 1:numel(idx)
    M(k,:) = double(Obs(idx(k)).(field_name));
end
end

function shaded_line(ax,x,y,se,col,line_style,line_width)
x = x(:).'; y = y(:).'; se = se(:).';
valid_fill = isfinite(x) & isfinite(y) & isfinite(se);
if sum(valid_fill) >= 2
    xf = x(valid_fill); yf = y(valid_fill); sf = se(valid_fill);
    fill(ax,[xf fliplr(xf)],[yf+sf fliplr(yf-sf)],col, ...
        'FaceAlpha',0.14,'EdgeColor','none','HandleVisibility','off');
end
plot(ax,x,y,'Color',col,'LineStyle',line_style, ...
    'LineWidth',line_width);
end

function plot_individual_condition(ax,Obs,idx,observed_field,linear_field, ...
        time_ms,mask,sigma_ms,observed_col,linear_col,strength,max_number, ...
        normalize)
idx = choose_trace_indices(idx,max_number);
for j = 1:numel(idx)
    io = idx(j);
    obs = Obs(io).(observed_field);
    lin = Obs(io).(linear_field);
    if normalize
        obs = obs/Obs(io).normalization_scale_sp_s;
        lin = lin/Obs(io).normalization_scale_sp_s;
    end
    obs = smooth_curve(obs,time_ms,sigma_ms);
    lin = smooth_curve(lin,time_ms,sigma_ms);
    plot(ax,time_ms(mask),lin(mask),'-', ...
        'Color',lighten_color(linear_col,1-strength), ...
        'LineWidth',0.5,'HandleVisibility','off');
    plot(ax,time_ms(mask),obs(mask),'-', ...
        'Color',lighten_color(observed_col,1-strength), ...
        'LineWidth',0.5,'HandleVisibility','off');
end
end

function plot_individual_residuals(ax,Obs,idx,time_ms,mask,sigma_ms, ...
        sim_col,seq_col,strength,max_number,normalize)
idx = choose_trace_indices(idx,max_number);
for j = 1:numel(idx)
    io = idx(j);
    sim = Obs(io).sim_residual;
    seq = Obs(io).seq_residual;
    if normalize
        sim = sim/Obs(io).normalization_scale_sp_s;
        seq = seq/Obs(io).normalization_scale_sp_s;
    end
    sim = smooth_curve(sim,time_ms,sigma_ms);
    seq = smooth_curve(seq,time_ms,sigma_ms);
    plot(ax,time_ms(mask),sim(mask),'-', ...
        'Color',lighten_color(sim_col,1-strength), ...
        'LineWidth',0.5,'HandleVisibility','off');
    plot(ax,time_ms(mask),seq(mask),'-', ...
        'Color',lighten_color(seq_col,1-strength), ...
        'LineWidth',0.5,'HandleVisibility','off');
end
end

function idx_out = choose_trace_indices(idx,max_number)
if isinf(max_number) || numel(idx) <= max_number
    idx_out = idx;
else
    positions = round(linspace(1,numel(idx),max_number));
    idx_out = idx(unique(positions));
end
end

function col_out = lighten_color(col,amount)
amount = max(0,min(1,amount));
col_out = col + (1-col)*amount;
end

function decorate_time_axis(ax,plot_window,ptd_ms)
xline(ax,0,'r--','LineWidth',0.8,'HandleVisibility','off');
if ~isempty(ptd_ms) && isfinite(ptd_ms(1)) && ptd_ms(1) > 0
    xline(ax,ptd_ms,':','Color',[0.35 0.35 0.35], ...
        'LineWidth',0.9,'HandleVisibility','off');
end
xlim(ax,plot_window);
box(ax,'off');
end

function apply_common_ylim(fig)
axes_list = findall(fig,'Type','axes');
if isempty(axes_list), return; end
lo = inf; hi = -inf;
for k = 1:numel(axes_list)
    yl = ylim(axes_list(k));
    lo = min(lo,yl(1)); hi = max(hi,yl(2));
end
if isfinite(lo) && isfinite(hi) && hi > lo
    pad = 0.03*(hi-lo);
    for k = 1:numel(axes_list)
        ylim(axes_list(k),[lo-pad hi+pad]);
    end
end
end

function apply_common_ylim_by_column(fig,nColumns)
axes_list = flipud(findall(fig,'Type','axes'));
if isempty(axes_list), return; end
for col = 1:nColumns
    group = axes_list(col:nColumns:end);
    lo = inf; hi = -inf;
    for k = 1:numel(group)
        yl = ylim(group(k));
        lo = min(lo,yl(1)); hi = max(hi,yl(2));
    end
    if isfinite(lo) && isfinite(hi) && hi > lo
        pad = 0.03*(hi-lo);
        for k = 1:numel(group)
            ylim(group(k),[lo-pad hi+pad]);
        end
    end
end
end

function jitter = deterministic_jitter(n,width)
if n <= 1
    jitter = zeros(n,1);
else
    jitter = linspace(-width,width,n).';
end
end
