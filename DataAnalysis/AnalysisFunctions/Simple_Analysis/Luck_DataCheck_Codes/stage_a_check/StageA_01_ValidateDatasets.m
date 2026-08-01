%% Stage A: validate and index Single, Simultaneous, and Sequential datasets
% Run 00_CheckSettings.m first by editing its paths, then run this script.
% This stage is read-only with respect to the experiment folders.

clear;
clc;

stage_a_dir = fileparts(mfilename('fullpath'));
run(fullfile(stage_a_dir, 'StageA_00_CheckSettings.m'));
addpath(fullfile(stage_a_dir, 'helpers'));

if isfolder(CFG.analysis_functions_folder)
    addpath(genpath(CFG.analysis_functions_folder));
else
    error('StageA:MissingAnalysisFunctions', ...
        'Analysis-functions folder not found: %s', ...
        CFG.analysis_functions_folder);
end

if ~isfolder(CFG.output_folder)
    mkdir(CFG.output_folder);
end

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
    end
end

CrossDataset = struct();
if all(isfield(DatasetInfo, roles)) && ...
        all(arrayfun(@(k) isfield(DatasetInfo.(roles{k}), 'experiment'), 1:numel(roles)))
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

StageA = struct();
StageA.created_on = datetime('now');
StageA.config = CFG;
StageA.datasets = DatasetInfo;
StageA.cross_dataset = CrossDataset;
StageA.validation_passed = ~has_error;
StageA.messages = all_messages;

output_mat = fullfile(CFG.output_folder, 'DX014_D4_StageA_DatasetInfo.mat');
save(output_mat, 'StageA', '-v7.3');

output_txt = fullfile(CFG.output_folder, 'DX014_D4_StageA_Report.txt');
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
