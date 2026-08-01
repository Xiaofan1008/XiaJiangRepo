%% Stage A: validate and index the three experiment datasets
% This is the only Stage A script that you need to run.
% Enter the settings below, then run the complete script. No pop-ups are used.
% It does not modify any experiment data.

clear;
clc;

%% ========================= USER SETTINGS ==============================
% Enter the full paths to the three datasets.
% If Sim and Seq are stored together, use the same path for both.
single_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Single4_new';
sim_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim4';
seq_folder    = sim_folder;

% Acquisition/electrode settings.
Electrode_Type = 2;       % 0=rigid, 1=single-shank flex, 2=four-shank flex
FS             = 30000;   % sampling rate in Hz

% Conditions and analysis windows used by later checking stages.
sim_PTD_ms       = 0;
seq_PTD_ms       = 5;
baseline_win_ms  = [-50 -5];
count_win_ms     = [2 20];
raster_win_ms    = [-50 80];

%% ======================= INITIAL SETUP ================================
script_folder = fileparts(mfilename('fullpath'));
helpers_folder = fullfile(script_folder, 'helpers');
if ~isfolder(helpers_folder)
    error('StageA:MissingHelpers', ...
        'The helpers folder must remain beside this script:\n%s', helpers_folder);
end
addpath(helpers_folder, '-begin');

if ~ismember(Electrode_Type, [0 1 2])
    error('StageA:InvalidElectrodeType', ...
        'Electrode type must be 0, 1, or 2.');
end
if ~isfinite(FS) || FS <= 0
    error('StageA:InvalidSamplingRate', ...
        'Sampling rate must be a positive number.');
end

%% ====================== AUTOMATIC SETTINGS ============================
CFG = struct();
CFG.single_folder = single_folder;
CFG.sim_folder = sim_folder;
CFG.seq_folder = seq_folder;
CFG.FS = FS;
CFG.Electrode_Type = Electrode_Type;

% Derive Data_Xia/AnalysisFunctions from Data_Xia/DXxxx/ExperimentFolder.
data_xia_folder = fileparts(fileparts(single_folder));
CFG.analysis_functions_folder = fullfile(data_xia_folder, 'AnalysisFunctions');
if ~isfolder(CFG.analysis_functions_folder)
    error('StageA:MissingAnalysisFunctions', ...
        ['AnalysisFunctions folder was not found automatically:\n%s\n' ...
         'Check that the dataset follows Data_Xia/DXxxx/ExperimentFolder.'], ...
        CFG.analysis_functions_folder);
end

% AnalysisFunctions may contain older copies of this checking package.
% Put that broad tree at the END of the path, then force the helpers shipped
% beside this script back to the BEGINNING to prevent function shadowing.
addpath(genpath(CFG.analysis_functions_folder), '-end');
addpath(helpers_folder, '-begin');
clear qc_validate_dataset qc_resolve_dataset_files qc_decode_stim_params ...
    qc_load_bad_trials qc_load_bad_channels qc_load_responding ...
    qc_validate_cross_dataset qc_write_stage_a_report;
rehash;

resolved_validator = which('qc_validate_dataset');
expected_validator = fullfile(helpers_folder, 'qc_validate_dataset.m');
if ~strcmp(resolved_validator, expected_validator)
    error('StageA:HelperShadowing', ...
        ['MATLAB resolved the wrong qc_validate_dataset function.\n' ...
         'Expected: %s\nResolved: %s'], ...
        expected_validator, resolved_validator);
end

CFG.sim_PTD_ms = sim_PTD_ms;
CFG.seq_PTD_ms = seq_PTD_ms;
CFG.PTD_tolerance_ms = 0.01;
CFG.baseline_win_ms = baseline_win_ms;
CFG.count_win_ms = count_win_ms;
CFG.raster_win_ms = raster_win_ms;
CFG.require_paired_qc_files = true;
CFG.output_folder = fullfile(script_folder, 'check_results');

if ~isfolder(CFG.output_folder)
    mkdir(CFG.output_folder);
end

fprintf('\nSelected Stage A inputs:\n');
fprintf('  Single: %s\n', CFG.single_folder);
fprintf('  Sim:    %s\n', CFG.sim_folder);
fprintf('  Seq:    %s\n', CFG.seq_folder);
fprintf('  Electrode type: %d\n', CFG.Electrode_Type);
fprintf('  Sampling rate: %.0f Hz\n', CFG.FS);
fprintf('  AnalysisFunctions: %s\n', CFG.analysis_functions_folder);

%% ====================== VALIDATE DATASETS =============================
fprintf('\n============================================================\n');
fprintf('STAGE A: DATASET VALIDATION\n');
fprintf('============================================================\n');

roles = {'single','sim','seq'};
folders = {CFG.single_folder, CFG.sim_folder, CFG.seq_folder};
require_qc = [false, CFG.require_paired_qc_files, CFG.require_paired_qc_files];

DatasetInfo = struct();
all_messages = {};
has_error = false;

for iRole = 1:numel(roles)
    role = roles{iRole};
    fprintf('\n--- Validating %s dataset ---\n', upper(role));

    try
        info = qc_validate_dataset(role, folders{iRole}, CFG, require_qc(iRole));
        DatasetInfo.(role) = info;
        all_messages = [all_messages; info.messages(:)]; %#ok<AGROW>

        fprintf('Folder: %s\n', info.folder);
        fprintf('Spike file: %s (variable: %s)\n', ...
            info.files.spike.name, info.spike.variable);
        fprintf('Experiment file: %s\n', info.files.experiment.name);
        fprintf('Trials/triggers: %d / %d\n', ...
            info.experiment.n_trials, info.trigger.n_triggers);
        fprintf('Recording channels: %d\n', info.spike.n_channels);
        fprintf('Stimulation events per trial: %d\n', ...
            info.experiment.simultaneous_stim);
        fprintf('Amplitudes: %s uA\n', num2str(info.experiment.amplitudes(:).'));
        fprintf('PTDs: %s ms\n', num2str(info.experiment.ptds_ms(:).'));

        if info.bad_trials.file_found
            fprintf('Bad trials: %d global exclusions; global consistency = %s\n', ...
                numel(info.bad_trials.global_trials), ...
                char(string(info.bad_trials.is_global)));
        else
            fprintf('Bad trials: no file found\n');
        end

        if info.responding.file_found
            fprintf('Responding file: %s\n', info.files.responding.name);
        else
            fprintf('Responding file: not present (allowed for SINGLE only)\n');
        end

        if ~info.validation_passed
            has_error = true;
        end
    catch ME
        has_error = true;
        DatasetInfo.(role) = struct('role',role, 'folder',folders{iRole}, ...
            'validation_passed',false, 'fatal_error',ME.message);
        msg = sprintf('[ERROR][%s] %s', upper(role), ME.message);
        all_messages{end+1,1} = msg; %#ok<SAGROW>
        fprintf('%s\n', msg);
        fprintf('Detailed MATLAB error report:\n%s\n', ...
            getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

CrossDataset = struct();
have_all_metadata = all(isfield(DatasetInfo, roles)) && ...
    all(arrayfun(@(k) isfield(DatasetInfo.(roles{k}), 'experiment'), ...
    1:numel(roles)));

if have_all_metadata
    CrossDataset = qc_validate_cross_dataset(DatasetInfo, CFG);
    all_messages = [all_messages; CrossDataset.messages(:)];
    if ~CrossDataset.validation_passed
        has_error = true;
    end

    fprintf('\n--- Cross-dataset checks ---\n');
    fprintf('Sampling rates consistent: %s\n', ...
        char(string(CrossDataset.fs_consistent)));
    fprintf('Recording-channel counts consistent: %s\n', ...
        char(string(CrossDataset.channel_count_consistent)));
    fprintf('Single-to-Sim electrode-name coverage: %.1f%%\n', ...
        100 * CrossDataset.single_to_sim_name_coverage);
    fprintf('Single-to-Seq electrode-name coverage: %.1f%%\n', ...
        100 * CrossDataset.single_to_seq_name_coverage);
    fprintf('Sim PTD %.3f ms available: %s\n', CFG.sim_PTD_ms, ...
        char(string(CrossDataset.sim_target_ptd_available)));
    fprintf('Seq PTD %.3f ms available: %s\n', CFG.seq_PTD_ms, ...
        char(string(CrossDataset.seq_target_ptd_available)));
end

%% ========================= SAVE & REPORT ==============================
StageA = struct();
StageA.created_on = datetime('now');
StageA.config = CFG;
StageA.datasets = DatasetInfo;
StageA.cross_dataset = CrossDataset;
StageA.validation_passed = ~has_error;
StageA.messages = all_messages;

animal_match = regexp(single_folder, 'DX\d+', 'match', 'once');
if isempty(animal_match)
    animal_match = 'UnknownDataset';
end
[~, pair_folder_name] = fileparts(sim_folder);
dataset_number = regexp(pair_folder_name, '\d+', 'match');
if isempty(dataset_number)
    dataset_tag = animal_match;
else
    dataset_tag = sprintf('%s_D%s', animal_match, dataset_number{end});
end

output_mat = fullfile(CFG.output_folder, ...
    [dataset_tag '_StageA_DatasetInfo.mat']);
output_txt = fullfile(CFG.output_folder, ...
    [dataset_tag '_StageA_Report.txt']);

save(output_mat, 'StageA', '-v7.3');
qc_write_stage_a_report(output_txt, StageA, true);

fprintf('\n============================================================\n');
if StageA.validation_passed
    fprintf('STAGE A PASSED\n');
else
    fprintf('STAGE A FOUND ITEMS THAT REQUIRE REVIEW\n');
end
fprintf('MAT result: %s\n', output_mat);
fprintf('Text report: %s\n', output_txt);
fprintf('============================================================\n\n');
