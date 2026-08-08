%% ========================================================================
%  CROSS-DATASET E1: INSPECT RAW DATASET-LEVEL RESULTS
%
%  Purpose
%    Calculate one result per dataset and amplitude before deciding how to
%    present the cross-dataset average.
%
%  For datasets containing multiple pairs/order branches:
%    1. Keep only branches with complete Linear, Sim, and Seq results.
%    2. Average those branches within the dataset.
%    3. Give the dataset one value at that amplitude.
%
%  Calculated values
%    Linear
%    SimMatched
%    Sequential
%    SimMinusLinear = SimMatched - Linear
%    SeqMinusLinear = Sequential - Linear
%    SimSublinear   = SimMatched < Linear
%    SeqSublinear   = Sequential < Linear
%    SeqGreaterSim  = Sequential > SimMatched
%
%  Missing amplitudes remain NaN. Nothing is interpolated or replaced by 0.
%  This inspection script does not create figures or save a new MAT file.
% ========================================================================

clear;

%% ============================= USER SETTINGS ============================

Results_Folder = ...
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq';

% Empty loads all files matching Result_File_Pattern.
% Alternatively, enter selected full MAT-file paths in this cell array.
Result_Files = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX005_Xia_Exp1_Sim_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX009_Xia_Exp1_Sim5_Linearity1_v1.mat';
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim7_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim8_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim7_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim8_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim9_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX012_Xia_Exp1_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX012_Xia_Exp1_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX012_Xia_Exp1_Sim6_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim7_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim8_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim6_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim1_Linearity1_v1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim7_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX016_Xia_Exp1_Seq_Full_1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX016_Xia_Exp1_Seq_Full_2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX016_Xia_Exp1_Seq_Full_3_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim4_Linearity1_v1.mat';
    


};
% Result_File_Pattern = '*_Linearity1_v1.mat';

% Empty uses the union of all positive amplitudes in the saved files.
Amplitudes_To_Analyze = [];
Include_Zero_Amplitude = false;
Amplitude_Tolerance_uA = 1e-6;

% true prints one readable table for every amplitude.
Print_Results = true;

% true hides missing dataset rows in the printed tables only. Missing rows
% are always retained in DatasetAmplitudeResults as NaN.
Print_Only_Valid_Rows = true;

%% ============================== FIND FILES ===============================

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

%% ============================== LOAD FILES ===============================

nDatasets = numel(Result_Files);
SavedData = cell(nDatasets,1);
AnimalID = strings(nDatasets,1);
DatasetID = strings(nDatasets,1);
DatasetLabel = strings(nDatasets,1);
all_amplitudes = [];

for id = 1:nDatasets
    file_path = Result_Files{id};
    if ~isfile(file_path)
        error('CrossDataset:FileNotFound','Result file not found: %s',file_path);
    end

    L = load(file_path,'Linearity1Saved');
    if ~isfield(L,'Linearity1Saved')
        error('CrossDataset:MissingSavedVariable', ...
            'Linearity1Saved is missing from %s.',file_path);
    end
    Saved = L.Linearity1Saved;
    if ~isfield(Saved,'branch_results') || isempty(Saved.branch_results)
        error('CrossDataset:MissingBranchResults', ...
            'No branch_results were found in %s.',file_path);
    end

    [AnimalID(id),DatasetID(id),DatasetLabel(id)] = ...
        identify_dataset(Saved,file_path);

    for ib = 1:numel(Saved.branch_results)
        P = Saved.branch_results(ib);
        if ~isfield(P,'average_summary') || ~istable(P.average_summary)
            error('CrossDataset:MissingAverageSummary', ...
                'A branch in %s does not contain average_summary.',file_path);
        end
        required_columns = {'Amplitude_uA','NUsed','Linear', ...
            'SimMatched','Sequential'};
        if ~all(ismember(required_columns, ...
                P.average_summary.Properties.VariableNames))
            error('CrossDataset:InvalidAverageSummary', ...
                'average_summary has missing columns in %s.',file_path);
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

if isempty(Amplitudes_To_Analyze)
    Amplitudes = sort(all_amplitudes);
else
    Amplitudes = unique(double(Amplitudes_To_Analyze(:).'),'stable');
    if ~Include_Zero_Amplitude
        Amplitudes(abs(Amplitudes) < Amplitude_Tolerance_uA) = [];
    end
end
if isempty(Amplitudes)
    error('CrossDataset:NoAmplitudes','No amplitudes remain for inspection.');
end

%% ==================== CALCULATE ONE VALUE PER DATASET ====================

nAmp = numel(Amplitudes);
nRows = nDatasets*nAmp;

RowAnimal = strings(nRows,1);
RowDataset = strings(nRows,1);
RowLabel = strings(nRows,1);
RowFile = strings(nRows,1);
Amplitude_uA = nan(nRows,1);
ValidComparison = false(nRows,1);
NBranchesUsed = zeros(nRows,1);
Linear = nan(nRows,1);
SimMatched = nan(nRows,1);
Sequential = nan(nRows,1);
SimMinusLinear = nan(nRows,1);
SeqMinusLinear = nan(nRows,1);
SimSublinear = nan(nRows,1);
SeqSublinear = nan(nRows,1);
SeqGreaterSim = nan(nRows,1);

row_out = 0;
for id = 1:nDatasets
    BranchResults = SavedData{id}.branch_results;

    for ia = 1:nAmp
        row_out = row_out+1;
        amp = Amplitudes(ia);

        RowAnimal(row_out) = AnimalID(id);
        RowDataset(row_out) = DatasetID(id);
        RowLabel(row_out) = DatasetLabel(id);
        RowFile(row_out) = string(Result_Files{id});
        Amplitude_uA(row_out) = amp;

        branch_linear = [];
        branch_sim = [];
        branch_seq = [];

        for ib = 1:numel(BranchResults)
            T = BranchResults(ib).average_summary;
            source_row = find(abs(double(T.Amplitude_uA)-amp) < ...
                Amplitude_Tolerance_uA,1);
            if isempty(source_row)
                continue;
            end

            linear_value = double(T.Linear(source_row));
            sim_value = double(T.SimMatched(source_row));
            seq_value = double(T.Sequential(source_row));
            n_used = double(T.NUsed(source_row));

            % Sim and Seq must come from the same complete saved branch.
            if n_used <= 0 || ~isfinite(linear_value) || ...
                    ~isfinite(sim_value) || ~isfinite(seq_value)
                continue;
            end

            branch_linear(end+1) = linear_value; %#ok<SAGROW>
            branch_sim(end+1) = sim_value; %#ok<SAGROW>
            branch_seq(end+1) = seq_value; %#ok<SAGROW>
        end

        if isempty(branch_linear)
            continue;
        end

        % Equal branch/pair weight within this dataset.
        Linear(row_out) = mean(branch_linear);
        SimMatched(row_out) = mean(branch_sim);
        Sequential(row_out) = mean(branch_seq);
        SimMinusLinear(row_out) = mean(branch_sim-branch_linear);
        SeqMinusLinear(row_out) = mean(branch_seq-branch_linear);
        NBranchesUsed(row_out) = numel(branch_linear);
        ValidComparison(row_out) = true;

        SimSublinear(row_out) = double(SimMatched(row_out) < Linear(row_out));
        SeqSublinear(row_out) = double(Sequential(row_out) < Linear(row_out));
        SeqGreaterSim(row_out) = double(Sequential(row_out) > SimMatched(row_out));
    end
end

DatasetAmplitudeResults = table(RowAnimal,RowDataset,RowLabel,Amplitude_uA, ...
    ValidComparison,NBranchesUsed,Linear,SimMatched,Sequential, ...
    SimMinusLinear,SeqMinusLinear,SimSublinear,SeqSublinear, ...
    SeqGreaterSim,RowFile, ...
    'VariableNames',{'Animal','Dataset','DatasetLabel','Amplitude_uA', ...
    'Valid','NBranches','Linear','SimMatched','Sequential', ...
    'SimMinusLinear','SeqMinusLinear','SimSublinear','SeqSublinear', ...
    'SeqGreaterSim','ResultFile'});

%% ======================= AMPLITUDE CHECK SUMMARY =========================

NDatasets = zeros(nAmp,1);
MeanLinear = nan(nAmp,1);
MeanSim = nan(nAmp,1);
MeanSeq = nan(nAmp,1);
MeanSimDifference = nan(nAmp,1);
MeanSeqDifference = nan(nAmp,1);
MedianSimDifference = nan(nAmp,1);
MedianSeqDifference = nan(nAmp,1);
PctSimSublinear = nan(nAmp,1);
PctSeqSublinear = nan(nAmp,1);
PctSeqGreaterSim = nan(nAmp,1);

for ia = 1:nAmp
    rows = DatasetAmplitudeResults.Valid & ...
        abs(DatasetAmplitudeResults.Amplitude_uA-Amplitudes(ia)) < ...
        Amplitude_Tolerance_uA;

    NDatasets(ia) = sum(rows);
    if NDatasets(ia) == 0
        continue;
    end

    MeanLinear(ia) = mean(DatasetAmplitudeResults.Linear(rows));
    MeanSim(ia) = mean(DatasetAmplitudeResults.SimMatched(rows));
    MeanSeq(ia) = mean(DatasetAmplitudeResults.Sequential(rows));
    MeanSimDifference(ia) = mean(DatasetAmplitudeResults.SimMinusLinear(rows));
    MeanSeqDifference(ia) = mean(DatasetAmplitudeResults.SeqMinusLinear(rows));
    MedianSimDifference(ia) = median( ...
        DatasetAmplitudeResults.SimMinusLinear(rows));
    MedianSeqDifference(ia) = median( ...
        DatasetAmplitudeResults.SeqMinusLinear(rows));
    PctSimSublinear(ia) = 100*mean( ...
        DatasetAmplitudeResults.SimSublinear(rows) == 1);
    PctSeqSublinear(ia) = 100*mean( ...
        DatasetAmplitudeResults.SeqSublinear(rows) == 1);
    PctSeqGreaterSim(ia) = 100*mean( ...
        DatasetAmplitudeResults.SeqGreaterSim(rows) == 1);
end

AmplitudeCheckSummary = table(Amplitudes(:),NDatasets,MeanLinear, ...
    MeanSim,MeanSeq,MeanSimDifference,MeanSeqDifference, ...
    MedianSimDifference,MedianSeqDifference,PctSimSublinear, ...
    PctSeqSublinear,PctSeqGreaterSim, ...
    'VariableNames',{'Amplitude_uA','NDatasets','MeanLinear','MeanSim', ...
    'MeanSeq','MeanSimDifference','MeanSeqDifference', ...
    'MedianSimDifference','MedianSeqDifference','PctSimSublinear', ...
    'PctSeqSublinear','PctSeqGreaterSim'});

%% ================================ PRINT =================================

fprintf('\nCross-dataset raw-result inspection\n');
fprintf('Datasets loaded: %d\n',nDatasets);
fprintf('No ratio or denominator threshold was used.\n');
fprintf('Sublinear columns: 1 = yes, 0 = no, NaN = missing.\n');

if Print_Results
    display_variables = {'DatasetLabel','Amplitude_uA','NBranches', ...
        'Linear','SimMatched','Sequential','SimMinusLinear', ...
        'SeqMinusLinear','SimSublinear','SeqSublinear','SeqGreaterSim'};

    for ia = 1:nAmp
        rows = abs(DatasetAmplitudeResults.Amplitude_uA-Amplitudes(ia)) < ...
            Amplitude_Tolerance_uA;
        if Print_Only_Valid_Rows
            rows = rows & DatasetAmplitudeResults.Valid;
        end
        fprintf('\nAmplitude %.3g uA\n',Amplitudes(ia));
        disp(DatasetAmplitudeResults(rows,display_variables));
    end
end

fprintf('\nAmplitude check summary\n');
disp(AmplitudeCheckSummary);

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
