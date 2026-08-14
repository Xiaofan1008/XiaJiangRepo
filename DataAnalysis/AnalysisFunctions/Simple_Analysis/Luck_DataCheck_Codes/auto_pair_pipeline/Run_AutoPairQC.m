%% ONE-RUN AUTOMATIC STIMULATION-PAIR QC PIPELINE
%
% This is the only file that you need to edit and run.
%
% The pipeline automatically:
%   1. validates the Single, Sim and Seq datasets;
%   2. discovers every two-electrode stimulation pair;
%   3. requires A, B, AB and at least one sequential order;
%   4. processes every eligible responsive pair;
%   5. aligns sp_corr spikes to the first pulse;
%   6. saves one QC data file per complete pair; and
%   7. optionally displays QC figures and exports Luke's model-data file.
%
% It never chooses a "best" pair and never modifies experiment files.
% Intermediate Stage A/B files are retained in the result folder as an
% audit trail, but they do not need to be opened or selected manually.

clear;
clc;

%% ========================= USER SETTINGS ==============================

% Enter the three dataset folders. If Sim and Seq are stored together,
% use the same folder for both.
single_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Single3';
sim_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim3';
seq_folder    = sim_folder;

% Recording and probe settings.
Electrode_Type = 2;       % 0=rigid, 1=single-shank flex, 2=four-shank flex
FS = 30000;               % sampling rate in Hz

% Leave empty to derive Data_Xia/AnalysisFunctions automatically from the
% Single folder. Enter a full path only for a nonstandard folder layout.
analysis_functions_folder = '';

% All results are placed here. Leave empty to create auto_pair_results
% beside this runner. The result folder must not be inside experiment data.
output_folder = '';

% Conditions used to define a complete pair.
sim_PTD_ms = 0;
seq_PTD_ms = 5;

% Spike windows in milliseconds relative to the first pulse.
raster_win_ms   = [-50 80];
baseline_win_ms = [-50 -5];
count_win_ms    = [2 20];

% Channel selection saved in each final pair QC file:
%   'responding' = manually corrected paired-condition response union
%   'all_common' = all channels not excluded for this pair
%   'manual'     = manual_depth_channels below
channel_mode = 'responding';
manual_depth_channels = [];

% Empty means all structurally complete amplitudes for each pair.
amplitudes_to_check = [];

% Optional global exclusions for the Single dataset. These can remain
% empty until Single raster review has been completed.
single_bad_trials = [];
single_review_complete = false;

% Export one self-contained ModelData MAT file for every accepted pair.
% All clean trials are retained. The balanced indices are an optional,
% reproducible subset and do not remove or alter the all-trial data.
Export_ModelData = true;
balancing_random_seed = 20260802;

% Existing modelling exports are protected by default. Change this to true
% only when you deliberately want a repeated run to replace them.
allow_modeldata_overwrite = false;

% Expected physical geometry of every stimulation pair.
require_same_stim_shank = true;
expected_stim_spacing_um = 200;
geometry_tolerance_um = 0.01;

%% =========================== QC FIGURES ===============================

Plot_QC_Figures = true;

% C2 figure choices:
%   1 = trial-count overview (one figure per pair)
%   2 = raster + PSTH (one figure per selected channel)
%   3 = per-trial spike-count dots (one figure per selected channel)
%   4 = amplitude curves (four channels per figure)
%   5 = average amplitude curves (one figure per pair)
%
% The default [1 2 5] retains the most useful pair-selection QC and avoids
% duplicating the heavier analyses available in the later D2/D3 scripts.
Plot_Figures = [1 5];

% 'responding', 'all', or a numeric depth-channel list.
Channels_To_Plot = 'responding';

% Empty means every amplitude saved in the pair QC file.
Amplitudes_To_Plot = [];

psth_bin_ms = 1;
psth_smooth_ms = 5;
baseline_correct_amplitude_curves = true;
channels_per_curve_figure = 4;
figure_position_large = [30 40 1800 950];
figure_position_medium = [100 80 1250 800];
Close_Existing_Figures = true;

% Continue to other automatically discovered pairs if one pair fails.
% The failure and its message will be saved in the run summary.
Continue_After_Pair_Error = true;

%% ====================== INITIALIZE PIPELINE ===========================
script_folder = fileparts(mfilename('fullpath'));
stage_a_script = fullfile(script_folder,'Pipeline_StageA.m');
stage_b_script = fullfile(script_folder,'Pipeline_StageB.m');
stage_c1_script = fullfile(script_folder,'Pipeline_StageC1.m');
stage_c2_script = fullfile(script_folder,'Pipeline_StageC2.m');
stage_d1_script = fullfile(script_folder,'Pipeline_StageD1.m');
required_scripts = {stage_a_script,stage_b_script,stage_c1_script, ...
    stage_c2_script,stage_d1_script};
for k = 1:numel(required_scripts)
    if ~isfile(required_scripts{k})
        error('AutoPairQC:MissingPipelineFile', ...
            'Required pipeline file is missing:\n%s',required_scripts{k});
    end
end
if ~isfolder(fullfile(script_folder,'helpers'))
    error('AutoPairQC:MissingHelpers', ...
        'The helpers folder must remain beside Run_AutoPairQC.m.');
end

validate_folder(single_folder,'single_folder');
validate_folder(sim_folder,'sim_folder');
validate_folder(seq_folder,'seq_folder');
if ~ismember(Electrode_Type,[0 1 2])
    error('AutoPairQC:ElectrodeType','Electrode_Type must be 0, 1, or 2.');
end
if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('AutoPairQC:SamplingRate','FS must be one positive number.');
end

if isempty(analysis_functions_folder)
    data_xia_folder = fileparts(fileparts(single_folder));
    analysis_functions_folder = fullfile(data_xia_folder,'AnalysisFunctions');
end
validate_folder(analysis_functions_folder,'analysis_functions_folder');

if isempty(output_folder)
    output_folder = fullfile(script_folder,'auto_pair_results');
end
output_folder = char(string(output_folder));
experiment_folders = {single_folder,sim_folder,seq_folder};
for k = 1:numel(experiment_folders)
    if path_is_inside(output_folder,experiment_folders{k})
        error('AutoPairQC:UnsafeOutput', ...
            'output_folder cannot be inside an experiment folder.');
    end
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end
if Close_Existing_Figures
    close all;
end

PipelineConfig = struct();
PipelineConfig.single_folder = char(string(single_folder));
PipelineConfig.sim_folder = char(string(sim_folder));
PipelineConfig.seq_folder = char(string(seq_folder));
PipelineConfig.Electrode_Type = Electrode_Type;
PipelineConfig.FS = FS;
PipelineConfig.analysis_functions_folder = char(string(analysis_functions_folder));
PipelineConfig.output_folder = output_folder;
PipelineConfig.sim_PTD_ms = sim_PTD_ms;
PipelineConfig.seq_PTD_ms = seq_PTD_ms;
PipelineConfig.baseline_win_ms = baseline_win_ms;
PipelineConfig.count_win_ms = count_win_ms;
PipelineConfig.raster_win_ms = raster_win_ms;
PipelineConfig.channel_mode = channel_mode;
PipelineConfig.manual_depth_channels = manual_depth_channels;
PipelineConfig.amplitudes_to_check = amplitudes_to_check;
PipelineConfig.single_bad_trials = single_bad_trials;
PipelineConfig.single_review_complete = single_review_complete;
PipelineConfig.Export_ModelData = Export_ModelData;
PipelineConfig.balancing_random_seed = balancing_random_seed;
PipelineConfig.allow_modeldata_overwrite = allow_modeldata_overwrite;
PipelineConfig.model_export_folder = fullfile(output_folder,'Luke_ModelPackage');
PipelineConfig.require_same_stim_shank = require_same_stim_shank;
PipelineConfig.expected_stim_spacing_um = expected_stim_spacing_um;
PipelineConfig.geometry_tolerance_um = geometry_tolerance_um;
PipelineConfig.Plot_Figures = Plot_Figures;
PipelineConfig.Channels_To_Plot = Channels_To_Plot;
PipelineConfig.Amplitudes_To_Plot = Amplitudes_To_Plot;
PipelineConfig.psth_bin_ms = psth_bin_ms;
PipelineConfig.psth_smooth_ms = psth_smooth_ms;
PipelineConfig.baseline_correct_amplitude_curves = ...
    baseline_correct_amplitude_curves;
PipelineConfig.channels_per_curve_figure = channels_per_curve_figure;
PipelineConfig.figure_position_large = figure_position_large;
PipelineConfig.figure_position_medium = figure_position_medium;

dataset_tag = derive_dataset_tag(single_folder,sim_folder);
PipelineConfig.dataset_tag = dataset_tag;
config_file = fullfile(output_folder,[dataset_tag '_AutoPairQC_Config.mat']);
save(config_file,'PipelineConfig');

old_config_environment = getenv('AUTO_PAIR_QC_CONFIG');
environment_cleanup = onCleanup(@() setenv( ...
    'AUTO_PAIR_QC_CONFIG',old_config_environment));
setenv('AUTO_PAIR_QC_CONFIG',config_file);

fprintf('\n============================================================\n');
fprintf('ONE-RUN AUTOMATIC PAIR QC PIPELINE\n');
fprintf('============================================================\n');
fprintf('Single: %s\n',single_folder);
fprintf('Sim:    %s\n',sim_folder);
fprintf('Seq:    %s\n',seq_folder);
fprintf('Output: %s\n',output_folder);
fprintf('Stimulation pairs will be discovered automatically.\n');
fprintf('Experiment files will not be modified.\n');

%% ===================== INTERNAL STAGE A ===============================
fprintf('\n>>> Internal step 1/5: validating datasets...\n');
run_stage_isolated(stage_a_script);
stage_a_file = fullfile(output_folder, ...
    [dataset_tag '_StageA_DatasetInfo.mat']);
if ~isfile(stage_a_file)
    error('AutoPairQC:StageAOutputMissing', ...
        'Expected Stage A output was not created:\n%s',stage_a_file);
end
tmp_a = load(stage_a_file,'StageA');
if ~tmp_a.StageA.validation_passed
    error('AutoPairQC:StageAFailed', ...
        'Dataset validation failed. Review the Stage A report before continuing.');
end

%% ===================== INTERNAL STAGE B ===============================
fprintf('\n>>> Internal step 2/5: discovering stimulation pairs...\n');
PipelineConfig.stage_a_file = stage_a_file;
save(config_file,'PipelineConfig');
run_stage_isolated(stage_b_script);
stage_b_file = fullfile(output_folder, ...
    [dataset_tag '_StageB_PairSummary.mat']);
if ~isfile(stage_b_file)
    error('AutoPairQC:StageBOutputMissing', ...
        'Expected pair-summary output was not created:\n%s',stage_b_file);
end
tmp_b = load(stage_b_file,'StageB');
StageB = tmp_b.StageB;

complete_pair_indices = find([StageB.pairs.is_exportable]);
fprintf('\n================ AUTOMATIC PAIR DISCOVERY ================\n');
disp(StageB.PairSummary);
fprintf('Candidate pairs: %d\n',StageB.n_candidate_pairs);
fprintf('Complete responsive pairs: %d\n',numel(complete_pair_indices));
if isempty(complete_pair_indices)
    error('AutoPairQC:NoCompletePairs', ...
        ['No pair with A, B, simultaneous A+B, at least one common ' ...
         'sequential order, and responding channels was found.']);
end

%% ===================== INTERNAL STAGE C1/C2 ===========================
fprintf('\n>>> Internal step 3/5: preparing every accepted pair...\n');
nComplete = numel(complete_pair_indices);
result_rows = cell(nComplete,9);

for iComplete = 1:nComplete
    pair_index = complete_pair_indices(iComplete);
    pair_info = StageB.pairs(pair_index);
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Processing complete pair %d/%d: Pair %d, %s\n', ...
        iComplete,nComplete,pair_index,pair_info.key);

    PipelineConfig.stage_b_file = stage_b_file;
    PipelineConfig.current_pair_index = pair_index;
    save(config_file,'PipelineConfig');

    pair_status = 'FAILED';
    error_message = '';
    model_file = '';
    c1_file = fullfile(output_folder,sprintf( ...
        '%s_StageC1_Pair%d_QCData.mat',dataset_tag,pair_index));
    try
        run_stage_isolated(stage_c1_script);
        if ~isfile(c1_file)
            error('AutoPairQC:C1OutputMissing', ...
                'Expected pair QC file was not created: %s',c1_file);
        end
        pair_status = 'QC_DATA_SAVED';

        if Plot_QC_Figures
            PipelineConfig.current_stage_c1_file = c1_file;
            save(config_file,'PipelineConfig');
            tmp_c1 = load(c1_file,'StageC1');
            estimated_figures = estimate_c2_figures(Plot_Figures, ...
                tmp_c1.StageC1,Channels_To_Plot,channels_per_curve_figure);
            fprintf('Expected QC figures for this pair: %d\n',estimated_figures);
            run_stage_isolated(stage_c2_script);
            pair_status = 'QC_DATA_AND_FIGURES_COMPLETE';
        end

        if Export_ModelData
            fprintf('\n>>> Internal step 4/5: exporting model data for this pair...\n');
            PipelineConfig.current_stage_c1_file = c1_file;
            save(config_file,'PipelineConfig');
            run_stage_isolated(stage_d1_script);
            model_file = fullfile(PipelineConfig.model_export_folder,sprintf( ...
                '%s_Pair_%s_%s_ModelData.mat',dataset_tag, ...
                safe_name(pair_info.electrode_A), ...
                safe_name(pair_info.electrode_B)));
            if ~isfile(model_file)
                error('AutoPairQC:D1OutputMissing', ...
                    'Expected model-data file was not created: %s',model_file);
            end
            if Plot_QC_Figures
                pair_status = 'QC_FIGURES_AND_MODEL_DATA_COMPLETE';
            else
                pair_status = 'MODEL_DATA_EXPORTED';
            end
        end
    catch ME
        error_message = ME.message;
        fprintf(2,'PAIR FAILED: %s\n',ME.message);
        fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off'));
        if ~Continue_After_Pair_Error
            rethrow(ME);
        end
    end

    result_rows(iComplete,:) = {pair_index,string(pair_info.key), ...
        string(pair_info.electrode_A),string(pair_info.electrode_B), ...
        string(strtrim(num2str(pair_info.complete_amplitudes_uA))), ...
        pair_info.n_responding_union,string(pair_status),string(model_file), ...
        string(error_message)};
end

fprintf('\n>>> Internal step 5/5: saving the combined run summary...\n');
PairRunSummary = cell2table(result_rows,'VariableNames', ...
    {'PairIndex','PairKey','ElectrodeA','ElectrodeB', ...
     'CompleteAmplitudes_uA','RespondingChannels','RunStatus', ...
     'ModelDataFile','ErrorMessage'});

PipelineRun = struct();
PipelineRun.created_on = datetime('now');
PipelineRun.dataset_tag = dataset_tag;
PipelineRun.config = PipelineConfig;
PipelineRun.stage_a_file = stage_a_file;
PipelineRun.stage_b_file = stage_b_file;
PipelineRun.complete_pair_indices = complete_pair_indices;
PipelineRun.PairSummary = StageB.PairSummary;
PipelineRun.PairRunSummary = PairRunSummary;
PipelineRun.experiment_files_modified = false;

summary_file = fullfile(output_folder, ...
    [dataset_tag '_AutoPairQC_RunSummary.mat']);
save(summary_file,'PipelineRun');

fprintf('\n==================== FINAL RUN SUMMARY ====================\n');
disp(PairRunSummary);
fprintf('Summary file: %s\n',summary_file);
fprintf('Output folder: %s\n',output_folder);
fprintf('Experiment files modified: NO\n');
fprintf('============================================================\n\n');

clear environment_cleanup;
setenv('AUTO_PAIR_QC_CONFIG',old_config_environment);

%% =========================== FUNCTIONS ================================
function run_stage_isolated(stage_file)
% Scripts are executed inside this function so their clear commands and
% working variables cannot erase or overwrite the master-runner workspace.
run(stage_file);
end

function validate_folder(folder_value,setting_name)
if ~isfolder(folder_value)
    error('AutoPairQC:MissingFolder','%s does not exist:\n%s', ...
        setting_name,folder_value);
end
end

function tag = derive_dataset_tag(single_folder,sim_folder)
animal_match = regexp(single_folder,'DX\d+','match','once');
if isempty(animal_match)
    animal_match = 'UnknownDataset';
end
[~,pair_folder_name] = fileparts(sim_folder);
dataset_number = regexp(pair_folder_name,'\d+','match');
if isempty(dataset_number)
    tag = animal_match;
else
    tag = sprintf('%s_D%s',animal_match,dataset_number{end});
end
end

function n = estimate_c2_figures(ids,C,channel_choice,channels_per_curve)
if isnumeric(channel_choice)
    nChannels = numel(unique(channel_choice));
elseif strcmpi(char(string(channel_choice)),'responding')
    nChannels = sum(C.ChannelInfo.IsRespondingUnion);
else
    nChannels = height(C.ChannelInfo);
end
n = 0;
if ismember(1,ids), n = n+1; end
if ismember(2,ids), n = n+nChannels; end
if ismember(3,ids), n = n+nChannels; end
if ismember(4,ids), n = n+ceil(nChannels/channels_per_curve); end
if ismember(5,ids), n = n+1; end
end

function tf = path_is_inside(candidate,parent)
candidate = normalize_path(candidate);
parent = normalize_path(parent);
tf = strcmp(candidate,parent) || ...
    startsWith([candidate filesep],[parent filesep]);
end

function value = normalize_path(value)
value = char(string(value));
while numel(value) > 1 && value(end) == filesep
    value(end) = [];
end
end

function value = safe_name(value)
value = regexprep(char(string(value)),'[^A-Za-z0-9]+','');
if isempty(value)
    value = 'Unknown';
end
end
