%% ========================================================================
%  CROSS-DATASET E1: CONDITION-SPECIFIC LINEAR PREDICTIONS
%
%  Figure 1
%    Cross-dataset mean +/- SEM of:
%       Sim linear (A+B unshifted), Simultaneous observed,
%       Seq linear (second component shifted), Sequential observed.
%
%  Figure 2
%    Cross-dataset mean +/- SEM, with dataset scatter, of:
%       Simultaneous - Sim linear
%       Sequential   - shifted Seq linear
%
%  Averaging rule
%    - Complete pair/order branches are averaged within each dataset first.
%    - Every dataset then receives equal weight.
%    - Missing amplitudes remain missing; no interpolation is performed.
%
%  Dataset-amplitude exclusion
%    A rule can exclude one dataset at one amplitude without excluding that
%    dataset at its other amplitudes. Saved E1 MAT files are never changed.
% ========================================================================

clear;

%% ============================= USER SETTINGS ============================

Results_Folder = ...
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse';

% false: automatically load every new ShiftedSeq_v2 result in the folder.
% true: load only the full paths listed in Result_Files below.
Use_Manual_Result_Files = true;

% Empty loads every file matching Result_File_Pattern.
% Alternatively, enter selected full MAT-file paths in this cell array.
Result_Files = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX005_Xia_Exp1_Sim_Linearity1_ShiftedSeq_v2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX009_Xia_Exp1_Sim5_Linearity1_ShiftedSeq_v2.mat';
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX010_Xia_Exp1_Sim1_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX010_Xia_Exp1_Sim2_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX010_Xia_Exp1_Sim4_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX010_Xia_Exp1_Sim5_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX010_Xia_Exp1_Sim6_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX010_Xia_Exp1_Sim7_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX010_Xia_Exp1_Sim8_Linearity1_ShiftedSeq_v2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim1_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim2_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim3_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim4_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim5_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim6_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim7_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim8_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX011_Xia_Exp1_Sim9_Linearity1_ShiftedSeq_v2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX012_Xia_Exp1_Sim1_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX012_Xia_Exp1_Sim4_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX012_Xia_Exp1_Sim6_Linearity1_ShiftedSeq_v2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX013_Xia_Exp1_Seq_Sim2_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX013_Xia_Exp1_Seq_Sim3_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX013_Xia_Exp1_Seq_Sim4_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX013_Xia_Exp1_Seq_Sim5_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX013_Xia_Exp1_Seq_Sim6_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX013_Xia_Exp1_Seq_Sim7_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX013_Xia_Exp1_Seq_Sim8_Linearity1_ShiftedSeq_v2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX014_Xia_Seq_Sim1_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX014_Xia_Seq_Sim2_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX014_Xia_Seq_Sim3_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX014_Xia_Seq_Sim4_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX014_Xia_Seq_Sim5_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX014_Xia_Seq_Sim6_Linearity1_ShiftedSeq_v2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX015_Xia_Seq_Sim1_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX015_Xia_Seq_Sim2_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX015_Xia_Seq_Sim3_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX015_Xia_Seq_Sim4_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX015_Xia_Seq_Sim5_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX015_Xia_Seq_Sim6_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX015_Xia_Seq_Sim7_Linearity1_ShiftedSeq_v2.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX016_Xia_Exp1_Seq_Full_1_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX016_Xia_Exp1_Seq_Full_2_Linearity1_ShiftedSeq_v2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX016_Xia_Exp1_Seq_Full_3_Linearity1_ShiftedSeq_v2.mat';

    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX018_Xia_Exp1_Sim1_Linearity1_ShiftedSeq_v2.mat'; 
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX018_Xia_Exp1_Sim2_Linearity1_ShiftedSeq_v2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX018_Xia_Exp1_Sim3_Linearity1_ShiftedSeq_v2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_ShiftedResponse/DX018_Xia_Exp1_Sim4_Linearity1_ShiftedSeq_v2.mat';
    
};
Result_File_Pattern = '*_Linearity1_ShiftedSeq_v2.mat';

% Empty uses the union of every positive amplitude in the saved files.
Amplitudes_To_Plot = [];
Include_Zero_Amplitude = false;
Amplitude_Tolerance_uA = 1e-6;

% ------------------------------------------------------------------------
% Optional dataset-amplitude exclusions
%
% Format: Dataset label, Amplitude in uA
% Copy the exact DatasetLabel shown in DatasetAmplitudeResults.
%
% Example:
% Exclude_Dataset_Amplitude = {
%     'DX018 | Xia_Exp1_Sim3', 6;
%     'DX014 | Xia_Seq_Sim2', 10;
% };
%
% Each rule affects only the specified amplitude. Empty means no exclusion.
% ------------------------------------------------------------------------
Exclude_Dataset_Amplitude = {
    % 'DX018 | Xia_Exp1_Sim3', 6;
    % 'DX018 | Xia_Exp1_Sim3', 1;
};

Show_Error_Bars = true;
Figure_Position = [100 100 900 600];

% Light dataset-level scatter points are drawn behind the summary curves.
% Excluded and missing dataset-amplitude values are not plotted.
Show_Dataset_Scatter = true;
Scatter_Marker_Size = 18;
Scatter_Alpha = 0.18;
Scatter_Jitter_Width = 0.10;
Scatter_Condition_Offset = 0.10;
Scatter_Random_Seed = 1;

% Optional plotting-only origin anchor. When true, (0 uA, 0) is prepended
% to every summary curve in both figures so the line connects to the
% origin. This is an assumed visual anchor, not a measured data point, and
% it is not added to calculations, tables, scatter points, or sample sizes.
Add_Zero_Plot_Anchor = true;

%% ============================== FIND FILES ===============================

if ~Use_Manual_Result_Files
    Result_Files = {};
end

if isempty(Result_Files)
    if ~isfolder(Results_Folder)
        error('CrossDataset:FolderNotFound', ...
            'Results folder not found: %s',Results_Folder);
    end
    file_info = dir(fullfile(Results_Folder,Result_File_Pattern));
    file_info = file_info(~startsWith({file_info.name},'._'));
    Result_Files = arrayfun(@(x) fullfile(x.folder,x.name), ...
        file_info,'UniformOutput',false);
else
    Result_Files = cellstr(string(Result_Files(:)));
end

Result_Files = unique(Result_Files,'stable');
if isempty(Result_Files)
    error('CrossDataset:NoResultFiles', ...
        'No result files matched %s in %s.', ...
        Result_File_Pattern,Results_Folder);
end

if ~iscell(Exclude_Dataset_Amplitude) || ...
        (~isempty(Exclude_Dataset_Amplitude) && size(Exclude_Dataset_Amplitude,2) ~= 2)
    error('CrossDataset:InvalidExclusionRules', 'Exclude_Dataset_Amplitude must be an empty or two-column cell array.');
end

%% ============================== LOAD FILES ===============================

nDatasets = numel(Result_Files);
SavedData = cell(nDatasets,1);
AnimalID = strings(nDatasets,1);
DatasetID = strings(nDatasets,1);
DatasetLabel = strings(nDatasets,1);
SavedPTD_ms = nan(nDatasets,1);
all_amplitudes = [];
reference_baseline_win = [];
reference_post_win = [];

for id = 1:nDatasets
    file_path = Result_Files{id};
    if ~isfile(file_path)
        error('CrossDataset:FileNotFound','Result file not found: %s',file_path);
    end

    L = load(file_path,'Linearity1Saved');
    if ~isfield(L,'Linearity1Saved')
        error('CrossDataset:MissingSavedVariable', 'Linearity1Saved is missing from %s.',file_path);
    end
    Saved = L.Linearity1Saved;
    if ~isfield(Saved,'format_version') || ...
            ~strcmp(char(string(Saved.format_version)),'2.0')
        error('CrossDataset:WrongResultVersion', ...
            ['%s is not a ShiftedSeq v2 result. Do not mix the old ' ...
             'common-linear files with this analysis.'],file_path);
    end
    if ~isfield(Saved,'settings') || ...
            ~isfield(Saved.settings,'linear_prediction_mode') || ...
            ~strcmp(char(string(Saved.settings.linear_prediction_mode)), ...
            'sim_unshifted_seq_second_component_shifted')
        error('CrossDataset:WrongLinearPredictionMode', ...
            'The shifted sequential prediction metadata is missing from %s.',file_path);
    end
    required_settings = {'baseline_win_ms','post_win_ms','sequential_ptd_ms'};
    if ~all(isfield(Saved.settings,required_settings))
        error('CrossDataset:MissingSettings', ...
            'Required window or PTD settings are missing from %s.',file_path);
    end
    this_baseline_win = double(Saved.settings.baseline_win_ms(:).');
    this_post_win = double(Saved.settings.post_win_ms(:).');
    SavedPTD_ms(id) = double(Saved.settings.sequential_ptd_ms);
    if id == 1
        reference_baseline_win = this_baseline_win;
        reference_post_win = this_post_win;
    elseif ~isequal(this_baseline_win,reference_baseline_win) || ...
            ~isequal(this_post_win,reference_post_win)
        error('CrossDataset:IncompatibleWindows', ...
            ['Baseline or response windows in %s differ from the first ' ...
             'result file.'],file_path);
    end
    if ~isfield(Saved,'branch_results') || isempty(Saved.branch_results)
        error('CrossDataset:MissingBranchResults', 'No branch_results were found in %s.',file_path);
    end

    [AnimalID(id),DatasetID(id),DatasetLabel(id)] = identify_dataset(Saved,file_path);

    for ib = 1:numel(Saved.branch_results)
        P = Saved.branch_results(ib);
        if ~isfield(P,'average_summary') || ~istable(P.average_summary)
            error('CrossDataset:MissingAverageSummary', 'A branch in %s does not contain average_summary.',file_path);
        end
        required_columns = {'Amplitude_uA','NUsed','SimLinear', ...
            'SeqLinear','SimMatched','Sequential','SimMinusLinear', ...
            'SeqMinusLinear'};
        if ~all(ismember(required_columns, ...
                P.average_summary.Properties.VariableNames))
            error('CrossDataset:InvalidAverageSummary', 'average_summary has missing columns in %s.',file_path);
        end
        all_amplitudes = [all_amplitudes; ...
            double(P.average_summary.Amplitude_uA(:))]; %#ok<AGROW>
    end
    SavedData{id} = Saved;
end

all_amplitudes = unique(all_amplitudes(isfinite(all_amplitudes))).';
if ~Include_Zero_Amplitude
    all_amplitudes(abs(all_amplitudes) < Amplitude_Tolerance_uA) = [];
end

if isempty(Amplitudes_To_Plot)
    Amplitudes = sort(all_amplitudes);
else
    Amplitudes = unique(double(Amplitudes_To_Plot(:).'),'stable');
    if ~Include_Zero_Amplitude
        Amplitudes(abs(Amplitudes) < Amplitude_Tolerance_uA) = [];
    end
end
if isempty(Amplitudes)
    error('CrossDataset:NoAmplitudes','No amplitudes remain for plotting.');
end

%% ==================== CALCULATE ONE VALUE PER DATASET ====================

nAmp = numel(Amplitudes);
SimLinearByDataset = nan(nDatasets,nAmp);
SeqLinearByDataset = nan(nDatasets,nAmp);
SimByDataset = nan(nDatasets,nAmp);
SeqByDataset = nan(nDatasets,nAmp);
SimDifferenceByDataset = nan(nDatasets,nAmp);
SeqDifferenceByDataset = nan(nDatasets,nAmp);
NBranchesUsed = zeros(nDatasets,nAmp);

for id = 1:nDatasets
    BranchResults = SavedData{id}.branch_results;

    for ia = 1:nAmp
        amp = Amplitudes(ia);
        branch_sim_linear = [];
        branch_seq_linear = [];
        branch_sim = [];
        branch_seq = [];

        for ib = 1:numel(BranchResults)
            T = BranchResults(ib).average_summary;
            source_row = find(abs(double(T.Amplitude_uA)-amp) < Amplitude_Tolerance_uA,1);
            if isempty(source_row)
                continue;
            end

            sim_linear_value = double(T.SimLinear(source_row));
            seq_linear_value = double(T.SeqLinear(source_row));
            sim_value = double(T.SimMatched(source_row));
            seq_value = double(T.Sequential(source_row));
            n_used = double(T.NUsed(source_row));

            % Use the same complete branch for all four response curves.
            if n_used <= 0 || ~isfinite(sim_linear_value) || ...
                    ~isfinite(seq_linear_value) || ~isfinite(sim_value) || ...
                    ~isfinite(seq_value)
                continue;
            end

            branch_sim_linear(end+1) = sim_linear_value; %#ok<SAGROW>
            branch_seq_linear(end+1) = seq_linear_value; %#ok<SAGROW>
            branch_sim(end+1) = sim_value; %#ok<SAGROW>
            branch_seq(end+1) = seq_value; %#ok<SAGROW>
        end

        if isempty(branch_sim_linear)
            continue;
        end

        % Equal pair/order-branch weight within the dataset.
        SimLinearByDataset(id,ia) = mean(branch_sim_linear);
        SeqLinearByDataset(id,ia) = mean(branch_seq_linear);
        SimByDataset(id,ia) = mean(branch_sim);
        SeqByDataset(id,ia) = mean(branch_seq);
        SimDifferenceByDataset(id,ia) = mean(branch_sim-branch_sim_linear);
        SeqDifferenceByDataset(id,ia) = mean(branch_seq-branch_seq_linear);
        NBranchesUsed(id,ia) = numel(branch_sim_linear);
    end
end

ValidComparison = isfinite(SimLinearByDataset) & ...
    isfinite(SeqLinearByDataset) & isfinite(SimByDataset) & ...
    isfinite(SeqByDataset);

%% ====================== APPLY OPTIONAL EXCLUSIONS ========================

Excluded = false(nDatasets,nAmp);
AppliedDataset = strings(0,1);
AppliedAmplitude = zeros(0,1);

for ir = 1:size(Exclude_Dataset_Amplitude,1)
    requested_label = strtrim(string(Exclude_Dataset_Amplitude{ir,1}));
    requested_amp = double(Exclude_Dataset_Amplitude{ir,2});
    if ~isscalar(requested_label) || strlength(requested_label) == 0 || ~isscalar(requested_amp) || ~isfinite(requested_amp)
        error('CrossDataset:InvalidExclusionRule', 'Exclusion rule %d has an invalid dataset label or amplitude.',ir);
    end

    id_match = find(DatasetLabel == requested_label);
    ia_match = find(abs(Amplitudes-requested_amp) < Amplitude_Tolerance_uA);
    if isempty(id_match)
        warning('CrossDataset:UnmatchedExclusionDataset', 'Exclusion dataset was not found: %s',requested_label);
        continue;
    end
    if numel(id_match) > 1
        error('CrossDataset:DuplicateDatasetLabel', 'Dataset label is not unique: %s',requested_label);
    end
    if isempty(ia_match)
        warning('CrossDataset:UnmatchedExclusionAmplitude', 'Exclusion amplitude %.3g uA was not selected.',requested_amp);
        continue;
    end
    if ~ValidComparison(id_match,ia_match)
        warning('CrossDataset:ExclusionAlreadyMissing', '%s already has no valid result at %.3g uA.', ...
            requested_label,requested_amp);
    end

    Excluded(id_match,ia_match) = true;
    AppliedDataset(end+1,1) = requested_label; %#ok<AGROW>
    AppliedAmplitude(end+1,1) = requested_amp; %#ok<AGROW>
end

UseForSummary = ValidComparison & ~Excluded;

%% ===================== DATASET-AMPLITUDE TABLE ===========================

nRows = nDatasets*nAmp;
RowAnimal = strings(nRows,1);
RowDataset = strings(nRows,1);
RowLabel = strings(nRows,1);
Amplitude_uA = nan(nRows,1);
Valid = false(nRows,1);
IsExcluded = false(nRows,1);
UseInSummary = false(nRows,1);
RowNBranches = zeros(nRows,1);
SimLinear = nan(nRows,1);
SeqLinear = nan(nRows,1);
SimMatched = nan(nRows,1);
Sequential = nan(nRows,1);
SimMinusLinear = nan(nRows,1);
SeqMinusLinear = nan(nRows,1);
SimSublinear = nan(nRows,1);
SeqSublinear = nan(nRows,1);
SeqGreaterSim = nan(nRows,1);

row = 0;
for id = 1:nDatasets
    for ia = 1:nAmp
        row = row+1;
        RowAnimal(row) = AnimalID(id);
        RowDataset(row) = DatasetID(id);
        RowLabel(row) = DatasetLabel(id);
        Amplitude_uA(row) = Amplitudes(ia);
        Valid(row) = ValidComparison(id,ia);
        IsExcluded(row) = Excluded(id,ia);
        UseInSummary(row) = UseForSummary(id,ia);
        RowNBranches(row) = NBranchesUsed(id,ia);
        SimLinear(row) = SimLinearByDataset(id,ia);
        SeqLinear(row) = SeqLinearByDataset(id,ia);
        SimMatched(row) = SimByDataset(id,ia);
        Sequential(row) = SeqByDataset(id,ia);
        SimMinusLinear(row) = SimDifferenceByDataset(id,ia);
        SeqMinusLinear(row) = SeqDifferenceByDataset(id,ia);
        if Valid(row)
            SimSublinear(row) = double(SimMatched(row) < SimLinear(row));
            SeqSublinear(row) = double(Sequential(row) < SeqLinear(row));
            SeqGreaterSim(row) = double(Sequential(row) > SimMatched(row));
        end
    end
end

DatasetAmplitudeResults = table(RowAnimal,RowDataset,RowLabel, ...
    Amplitude_uA,Valid,IsExcluded,UseInSummary,RowNBranches, ...
    SimLinear,SeqLinear,SimMatched,Sequential,SimMinusLinear,SeqMinusLinear, ...
    SimSublinear,SeqSublinear,SeqGreaterSim, ...
    'VariableNames',{'Animal','Dataset','DatasetLabel','Amplitude_uA', ...
    'Valid','Excluded','Used','NBranches','SimLinear','SeqLinear', ...
    'SimMatched','Sequential','SimMinusLinear','SeqMinusLinear','SimSublinear', ...
    'SeqSublinear','SeqGreaterSim'});

%% ============================= SUMMARIES ================================

NDatasets = zeros(nAmp,1);
NAnimals = zeros(nAmp,1);
NBranches = zeros(nAmp,1);
MeanSimLinear = nan(nAmp,1); SEMSimLinear = nan(nAmp,1);
MeanSeqLinear = nan(nAmp,1); SEMSeqLinear = nan(nAmp,1);
MeanSim = nan(nAmp,1); SEMSim = nan(nAmp,1);
MeanSeq = nan(nAmp,1); SEMSeq = nan(nAmp,1);
MeanSimDifference = nan(nAmp,1); SEMSimDifference = nan(nAmp,1);
MeanSeqDifference = nan(nAmp,1); SEMSeqDifference = nan(nAmp,1);
MedianSimDifference = nan(nAmp,1);
Q25SimDifference = nan(nAmp,1); Q75SimDifference = nan(nAmp,1);
MedianSeqDifference = nan(nAmp,1);
Q25SeqDifference = nan(nAmp,1); Q75SeqDifference = nan(nAmp,1);
PctSimSublinear = nan(nAmp,1);
PctSeqSublinear = nan(nAmp,1);
PctSeqGreaterSim = nan(nAmp,1);

for ia = 1:nAmp
    use = UseForSummary(:,ia);
    [MeanSimLinear(ia),SEMSimLinear(ia),n_sim_linear] = ...
        mean_sem(SimLinearByDataset(use,ia));
    [MeanSeqLinear(ia),SEMSeqLinear(ia),n_seq_linear] = ...
        mean_sem(SeqLinearByDataset(use,ia));
    [MeanSim(ia),SEMSim(ia),n_sim] = mean_sem(SimByDataset(use,ia));
    [MeanSeq(ia),SEMSeq(ia),n_seq] = mean_sem(SeqByDataset(use,ia));
    [MeanSimDifference(ia),SEMSimDifference(ia),n_sim_diff] = mean_sem(SimDifferenceByDataset(use,ia));
    [MeanSeqDifference(ia),SEMSeqDifference(ia),n_seq_diff] = mean_sem(SeqDifferenceByDataset(use,ia));
    [MedianSimDifference(ia),Q25SimDifference(ia),Q75SimDifference(ia)] = median_iqr(SimDifferenceByDataset(use,ia));
    [MedianSeqDifference(ia),Q25SeqDifference(ia),Q75SeqDifference(ia)] = median_iqr(SeqDifferenceByDataset(use,ia));

    counts = [n_sim_linear,n_seq_linear,n_sim,n_seq,n_sim_diff,n_seq_diff];
    if numel(unique(counts)) ~= 1
        error('CrossDataset:InternalPairingError', 'Dataset counts became inconsistent at %.3g uA.',Amplitudes(ia));
    end
    NDatasets(ia) = n_sim_linear;
    NAnimals(ia) = numel(unique(AnimalID(use)));
    NBranches(ia) = sum(NBranchesUsed(use,ia));
    if n_sim_linear > 0
        PctSimSublinear(ia) = 100*mean( ...
            SimByDataset(use,ia) < SimLinearByDataset(use,ia));
        PctSeqSublinear(ia) = 100*mean( ...
            SeqByDataset(use,ia) < SeqLinearByDataset(use,ia));
        PctSeqGreaterSim(ia) = 100*mean( SeqByDataset(use,ia) > SimByDataset(use,ia));
    end
end

ResponseSummary = table(Amplitudes(:),NDatasets,NAnimals,NBranches, ...
    MeanSimLinear,SEMSimLinear,MeanSim,SEMSim, ...
    MeanSeqLinear,SEMSeqLinear,MeanSeq,SEMSeq, ...
    'VariableNames',{'Amplitude_uA','NDatasets','NAnimals','NBranches', ...
    'MeanSimLinear','SEMSimLinear','MeanSim','SEMSim', ...
    'MeanSeqLinear','SEMSeqLinear','MeanSeq','SEMSeq'});

DifferenceSummary = table(Amplitudes(:),NDatasets,NAnimals,NBranches, ...
    MeanSimDifference, ...
    SEMSimDifference,MedianSimDifference,Q25SimDifference,Q75SimDifference, ...
    MeanSeqDifference,SEMSeqDifference,MedianSeqDifference, ...
    Q25SeqDifference,Q75SeqDifference,PctSimSublinear,PctSeqSublinear, ...
    PctSeqGreaterSim, ...
    'VariableNames',{'Amplitude_uA','NDatasets','NAnimals','NBranches', ...
    'MeanSimDifference', ...
    'SEMSimDifference','MedianSimDifference','Q25SimDifference', ...
    'Q75SimDifference','MeanSeqDifference','SEMSeqDifference', ...
    'MedianSeqDifference','Q25SeqDifference','Q75SeqDifference', ...
    'PctSimSublinear','PctSeqSublinear','PctSeqGreaterSim'});


%% ===================== RESPONDING CHANNEL COUNTS =========================
% TotalUniqueChannels:
%   A depth channel is counted once per dataset and amplitude, even if it
%   appears in multiple pair/order branches.
%
% TotalBranchChannelObservations:
%   Sums NUsed across every included branch. The same physical channel can
%   therefore be counted more than once if used in multiple branches.

TotalUniqueChannels = zeros(nAmp,1);
TotalBranchChannelObservations = zeros(nAmp,1);

for ia = 1:nAmp

    amp = Amplitudes(ia);

    for id = find(UseForSummary(:,ia)).'

        BranchResults = SavedData{id}.branch_results;
        dataset_unique_channels = [];

        for ib = 1:numel(BranchResults)

            P = BranchResults(ib);
            T = P.average_summary;

            source_row = find(abs(double(T.Amplitude_uA)-amp) < Amplitude_Tolerance_uA,1);

            if isempty(source_row)
                continue;
            end

            % Require the same complete branch used by the main analysis.
            if T.NUsed(source_row) <= 0 || ...
                    ~isfinite(T.SimLinear(source_row)) || ...
                    ~isfinite(T.SeqLinear(source_row)) || ...
                    ~isfinite(T.SimMatched(source_row)) || ...
                    ~isfinite(T.Sequential(source_row))
                continue;
            end

            % Count channel observations across branches.
            TotalBranchChannelObservations(ia) = TotalBranchChannelObservations(ia) + T.NUsed(source_row);

            % Collect unique physical depth channels within this dataset.
            p_amp = find(abs(double(P.amplitudes_uA)-amp) < ...
                Amplitude_Tolerance_uA,1);

            if ~isempty(p_amp) && ...
                    isfield(P,'channels_used_by_amplitude') && numel(P.channels_used_by_amplitude) >= p_amp

                channels_this_branch = ...
                    P.channels_used_by_amplitude{p_amp};

                dataset_unique_channels = union(dataset_unique_channels, double(channels_this_branch(:)).', 'stable');
            end
        end

        TotalUniqueChannels(ia) = TotalUniqueChannels(ia) + ...
            numel(dataset_unique_channels);
    end
end

RespondingChannelSummary = table( ...
    Amplitudes(:), NDatasets, TotalUniqueChannels, TotalBranchChannelObservations, 'VariableNames', {'Amplitude_uA', 'NDatasets', 'TotalUniqueRespondingChannels', 'TotalBranchChannelObservations'});

fprintf('\nResponding channels included at each amplitude\n');
disp(RespondingChannelSummary);

%% ================================ PRINT =================================

fprintf('\nCross-dataset shifted-response E1 summary\n');
fprintf('Datasets loaded: %d\n',nDatasets);
fprintf('Animals loaded: %d\n',numel(unique(AnimalID)));
fprintf('Windows: baseline [%g,%g) ms; response [%g,%g) ms\n', ...
    reference_baseline_win(1),reference_baseline_win(2), ...
    reference_post_win(1),reference_post_win(2));
fprintf('Saved PTD labels (ms): ');
fprintf('%g ',unique(SavedPTD_ms(isfinite(SavedPTD_ms))));
fprintf('\n');
fprintf(['Averaging: branches averaged within dataset; datasets receive ' ...
    'equal weight\n']);
fprintf('Dataset-amplitude exclusions applied: %d\n',numel(AppliedAmplitude));
if ~isempty(AppliedAmplitude)
    disp(table(AppliedDataset,AppliedAmplitude, ...
        'VariableNames',{'DatasetLabel','Amplitude_uA'}));
end

fprintf('\nRaw response mean +/- SEM\n');
disp(ResponseSummary);
fprintf('\nNonlinear difference and consistency summary\n');
disp(DifferenceSummary);

%% ============================== FIGURE 1 ================================

sim_linear_color = [0 0 0];
seq_linear_color = [0.45 0.45 0.45];
sim_color = [0 0.35 0.85];
seq_color = [0.85 0.33 0.10];

% Plotting-only response arrays. Calculated results remain unchanged.
ResponsePlotX = Amplitudes;
PlotMeanSimLinear = MeanSimLinear;
PlotSEMSimLinear = SEMSimLinear;
PlotMeanSeqLinear = MeanSeqLinear;
PlotSEMSeqLinear = SEMSeqLinear;
PlotMeanSim = MeanSim;
PlotSEMSim = SEMSim;
PlotMeanSeq = MeanSeq;
PlotSEMSeq = SEMSeq;
if Add_Zero_Plot_Anchor && ...
        ~any(abs(ResponsePlotX) < Amplitude_Tolerance_uA)
    ResponsePlotX = [0 ResponsePlotX];
    PlotMeanSimLinear = [0; PlotMeanSimLinear];
    PlotSEMSimLinear = [0; PlotSEMSimLinear];
    PlotMeanSeqLinear = [0; PlotMeanSeqLinear];
    PlotSEMSeqLinear = [0; PlotSEMSeqLinear];
    PlotMeanSim = [0; PlotMeanSim];
    PlotSEMSim = [0; PlotSEMSim];
    PlotMeanSeq = [0; PlotMeanSeq];
    PlotSEMSeq = [0; PlotSEMSeq];
end

figure('Color','w','Position',Figure_Position);
ax = axes();
hold(ax,'on');

% Plot dataset distributions first so the summary curves remain on top.
if Show_Dataset_Scatter
    rng(Scatter_Random_Seed,'twister');
    plot_dataset_scatter(ax,Amplitudes,SimLinearByDataset,UseForSummary, ...
        sim_linear_color,-1.5*Scatter_Condition_Offset,Scatter_Jitter_Width, ...
        Scatter_Marker_Size,Scatter_Alpha);
    plot_dataset_scatter(ax,Amplitudes,SimByDataset,UseForSummary, ...
        sim_color,-0.5*Scatter_Condition_Offset,Scatter_Jitter_Width, ...
        Scatter_Marker_Size,Scatter_Alpha);
    plot_dataset_scatter(ax,Amplitudes,SeqLinearByDataset,UseForSummary, ...
        seq_linear_color,0.5*Scatter_Condition_Offset,Scatter_Jitter_Width, ...
        Scatter_Marker_Size,Scatter_Alpha);
    plot_dataset_scatter(ax,Amplitudes,SeqByDataset,UseForSummary, ...
        seq_color,1.5*Scatter_Condition_Offset,Scatter_Jitter_Width, ...
        Scatter_Marker_Size,Scatter_Alpha);
end

plot_mean_curve(ax,ResponsePlotX,PlotMeanSimLinear,PlotSEMSimLinear, ...
    sim_linear_color,'--o','Sim linear (A+B)',Show_Error_Bars);
plot_mean_curve(ax,ResponsePlotX,PlotMeanSim,PlotSEMSim,sim_color, ...
    '-o','Simultaneous',Show_Error_Bars);
plot_mean_curve(ax,ResponsePlotX,PlotMeanSeqLinear,PlotSEMSeqLinear, ...
    seq_linear_color,'--s','Seq linear (shifted)',Show_Error_Bars);
plot_mean_curve(ax,ResponsePlotX,PlotMeanSeq,PlotSEMSeq,seq_color, ...
    '-o','Sequential',Show_Error_Bars);
xlabel(ax,'Amplitude (uA)');
ylabel(ax,'Baseline-corrected evoked spikes/trial');
title(ax,'Condition-specific response magnitude [2,20) ms', ...
    'FontWeight','bold');
legend(ax,'Location','northeast','Box','off');
box(ax,'off');
set(ax,'XTick',ResponsePlotX);
% add_n_labels(ax,Amplitudes,NDatasets);

%% ============================== FIGURE 2 ================================

figure('Color','w','Position',Figure_Position+[60 40 0 0]);
ax = axes();
hold(ax,'on');
yline(ax,0,'--','Color',[0.25 0.25 0.25], ...
    'LineWidth',1.2,'HandleVisibility','off');

if Show_Dataset_Scatter
    rng(Scatter_Random_Seed+1,'twister');
    plot_dataset_scatter(ax,Amplitudes,SimDifferenceByDataset, ...
        UseForSummary,sim_color,-0.8*Scatter_Condition_Offset, ...
        Scatter_Jitter_Width,Scatter_Marker_Size,Scatter_Alpha);
    plot_dataset_scatter(ax,Amplitudes,SeqDifferenceByDataset, ...
        UseForSummary,seq_color,0.8*Scatter_Condition_Offset, ...
        Scatter_Jitter_Width,Scatter_Marker_Size,Scatter_Alpha);
end

% Plotting-only nonlinear-difference arrays.
DifferencePlotX = Amplitudes;
PlotMeanSimDifference = MeanSimDifference;
PlotSEMSimDifference = SEMSimDifference;
PlotMeanSeqDifference = MeanSeqDifference;
PlotSEMSeqDifference = SEMSeqDifference;
if Add_Zero_Plot_Anchor && ...
        ~any(abs(DifferencePlotX) < Amplitude_Tolerance_uA)
    DifferencePlotX = [0 DifferencePlotX];
    PlotMeanSimDifference = [0; PlotMeanSimDifference];
    PlotSEMSimDifference = [0; PlotSEMSimDifference];
    PlotMeanSeqDifference = [0; PlotMeanSeqDifference];
    PlotSEMSeqDifference = [0; PlotSEMSeqDifference];
end
% plot_median_iqr_curve(ax,Amplitudes,MedianSimDifference, ...
%     Q25SimDifference,Q75SimDifference,sim_color, ...
%     'Simultaneous - linear',Show_Error_Bars);
% plot_median_iqr_curve(ax,Amplitudes,MedianSeqDifference, ...
%     Q25SeqDifference,Q75SeqDifference,seq_color, ...
%     'Sequential - linear',Show_Error_Bars);

plot_mean_curve(ax,DifferencePlotX,PlotMeanSimDifference, ...
    PlotSEMSimDifference,sim_color,'-o', ...
    'Simultaneous - Sim linear',Show_Error_Bars);

plot_mean_curve(ax,DifferencePlotX,PlotMeanSeqDifference, ...
    PlotSEMSeqDifference,seq_color,'-o', ...
    'Sequential - shifted Seq linear',Show_Error_Bars);

xlabel(ax,'Amplitude (uA)');
ylabel(ax,'Observed - condition-specific linear (spikes/trial)');
title(ax,'Cross-dataset nonlinear difference','FontWeight','bold');
legend(ax,'Location','northeast','Box','off');
box(ax,'off');
set(ax,'XTick',DifferencePlotX);
% add_n_labels(ax,Amplitudes,NDatasets);

%% ============================== FUNCTIONS ================================

function [animal,dataset,label] = identify_dataset(Saved,file_path)
[~,file_name] = fileparts(file_path);
animal = "";
dataset = "";
if isfield(Saved,'animal_id') && ~isempty(Saved.animal_id)
    animal = string(Saved.animal_id);
end
if isfield(Saved,'dataset_id') && ~isempty(Saved.dataset_id)
    dataset = string(Saved.dataset_id);
end
if strlength(animal) == 0
    animal = "unknown";
end
if strlength(dataset) == 0
    dataset = string(file_name);
end
label = animal+" | "+dataset;
end

function [avg,se,n] = mean_sem(values)
values = double(values(:));
values = values(isfinite(values));
n = numel(values);
if n == 0
    avg = NaN;
    se = NaN;
    return;
end
avg = mean(values);
if n > 1
    se = std(values,0)/sqrt(n);
else
    se = NaN;
end
end

function [med,q25,q75] = median_iqr(values)
values = sort(double(values(:)));
values = values(isfinite(values));
if isempty(values)
    med = NaN;
    q25 = NaN;
    q75 = NaN;
    return;
end
med = percentile_linear(values,50);
q25 = percentile_linear(values,25);
q75 = percentile_linear(values,75);
end

function value = percentile_linear(sorted_values,percent)
n = numel(sorted_values);
if n == 1
    value = sorted_values(1);
    return;
end
position = 1+(n-1)*(percent/100);
lower_index = floor(position);
upper_index = ceil(position);
fraction = position-lower_index;
value = sorted_values(lower_index)+fraction* ...
    (sorted_values(upper_index)-sorted_values(lower_index));
end

function plot_mean_curve(ax,x,y,se,color,line_style,label,show_errors)
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
        'Color',color,'MarkerFaceColor',color,'LineWidth',2, ...
        'MarkerSize',6,'CapSize',7,'DisplayName',label);
else
    plot(ax,x(valid),y(valid),line_style,'Color',color, ...
        'MarkerFaceColor',color,'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',label);
end
end

function plot_dataset_scatter(ax,x,values,use_mask,color,x_offset, ...
        jitter_width,marker_size,marker_alpha)
% Plot one point per included dataset and amplitude. Horizontal offsets
% separate conditions; random jitter prevents identical points overlapping.
x = double(x(:).');
values = double(values);
use_mask = logical(use_mask);
for ia = 1:numel(x)
    valid = use_mask(:,ia) & isfinite(values(:,ia));
    n = sum(valid);
    if n == 0
        continue;
    end
    x_points = x(ia)+x_offset+ ...
        (rand(n,1)-0.5)*jitter_width;
    scatter(ax,x_points,values(valid,ia),marker_size, ...
        'MarkerFaceColor',color,'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',marker_alpha, ...
        'MarkerEdgeAlpha',marker_alpha, ...
        'HandleVisibility','off');
end
end

function plot_median_iqr_curve(ax,x,med,q25,q75,color,label,show_errors)
x = double(x(:).');
med = double(med(:).');
q25 = double(q25(:).');
q75 = double(q75(:).');
valid = isfinite(x) & isfinite(med);
if ~any(valid)
    return;
end
if show_errors
    lower = med-q25;
    upper = q75-med;
    lower(~isfinite(lower)) = 0;
    upper(~isfinite(upper)) = 0;
    errorbar(ax,x(valid),med(valid),lower(valid),upper(valid),'-o', ...
        'Color',color,'MarkerFaceColor',color,'LineWidth',2, ...
        'MarkerSize',6,'CapSize',7,'DisplayName',label);
else
    plot(ax,x(valid),med(valid),'-o','Color',color, ...
        'MarkerFaceColor',color,'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',label);
end
end

function add_n_labels(ax,x,n_values)
valid_n = n_values > 0;
if ~any(valid_n)
    return;
end
yl = ylim(ax);
span = yl(2)-yl(1);
if ~isfinite(span) || span <= 0
    return;
end
for k = find(valid_n(:).')
    text(ax,x(k),yl(1)+0.025*span,sprintf('n=%d',n_values(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',8,'Color',[0.3 0.3 0.3]);
end
end
