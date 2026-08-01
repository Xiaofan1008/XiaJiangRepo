%% Stage A settings: edit this file before running Stage A
% This script only defines CFG. It does not load or modify experiment data.

CFG = struct();

% Dataset folders. sim_folder and seq_folder may point to the same folder.
CFG.single_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Single4_new';
CFG.sim_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim4';
CFG.seq_folder    = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim4';

% Existing analysis functions used only by the internal checking workflow.
% This is derived automatically from the usual folder structure:
%   Data_Xia/DXxxx/ExperimentFolder
%   Data_Xia/AnalysisFunctions
% You normally do not need to edit this path.
data_xia_folder = fileparts(fileparts(CFG.single_folder));
CFG.analysis_functions_folder = fullfile(data_xia_folder, 'AnalysisFunctions');

% Acquisition/electrode settings.
CFG.FS             = 30000;
CFG.Electrode_Type = 2;  % 0=rigid, 1=single-shank flex, 2=four-shank flex

% Conditions that later stages will use. Stage A only validates availability.
CFG.sim_PTD_ms = 0;
CFG.seq_PTD_ms = 5;
CFG.PTD_tolerance_ms = 0.01;

% Default checking windows. Later scripts may override these.
CFG.baseline_win_ms = [-50 -5];
CFG.count_win_ms    = [2 20];
CFG.raster_win_ms   = [-50 80];

% Stage A output location.
this_settings_dir = fileparts(mfilename('fullpath'));
CFG.output_folder = fullfile(this_settings_dir, 'check_results');

% If true, Stage A treats missing paired-data QC/response files as errors.
% The single dataset is allowed to have no QC or RespondingChannels files.
CFG.require_paired_qc_files = true;
