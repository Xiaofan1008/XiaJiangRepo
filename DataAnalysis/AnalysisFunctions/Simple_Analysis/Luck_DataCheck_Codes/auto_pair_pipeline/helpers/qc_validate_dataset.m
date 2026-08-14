function info = qc_validate_dataset(role, folder, CFG, require_qc)
%QC_VALIDATE_DATASET Resolve files, decode metadata, and perform validation.

arguments
    role (1,:) char
    folder (1,:) char
    CFG struct
    require_qc (1,1) logical
end

if ~isfolder(folder)
    error('StageA:FolderNotFound', 'Dataset folder not found: %s', folder);
end

info = struct();
info.role = role;
info.folder = folder;
info.messages = {};
info.validation_passed = true;

files = qc_resolve_dataset_files(folder);
info.files = files;

if ~files.spike.found
    error('StageA:MissingSpikeFile', ...
        'No unique *sp_xia_SSD.mat file found in %s', folder);
end
if ~files.experiment.found
    error('StageA:MissingExperimentFile', ...
        'No unique *_exp_datafile_*.mat file found in %s', folder);
end
if ~files.trigger.found
    error('StageA:MissingTriggerFile', ...
        'No unique *.trig.dat file found in %s. Stage A will not create one.', folder);
end

% Inspect sp_corr without loading the large cell array.
sp_vars = whos('-file', files.spike.path);
sp_idx = find(strcmp({sp_vars.name}, 'sp_corr'), 1);
if isempty(sp_idx)
    error('StageA:MissingSpCorr', 'sp_corr is missing from %s', files.spike.path);
end

sp_meta = sp_vars(sp_idx);
if ~strcmp(sp_meta.class, 'cell')
    error('StageA:InvalidSpCorr', 'sp_corr exists but is not a cell array in %s', ...
        files.spike.path);
end

info.spike = struct();
info.spike.variable = 'sp_corr';
info.spike.class = sp_meta.class;
info.spike.size = sp_meta.size;
info.spike.bytes = sp_meta.bytes;
info.spike.n_channels = numel_from_size(sp_meta.size);

% Load only experiment metadata fields.
E = load(files.experiment.path, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
required_fields = {'StimParams','simultaneous_stim','E_MAP','n_Trials'};
for k = 1:numel(required_fields)
    if ~isfield(E, required_fields{k})
        error('StageA:MissingExperimentField', ...
            'Missing %s in %s', required_fields{k}, files.experiment.path);
    end
end
info.experiment = qc_decode_stim_params(E);

% Load triggers using the established project function.
if exist('loadTrig', 'file') ~= 2
    error('StageA:MissingLoadTrig', ...
        'loadTrig was not found after adding %s', CFG.analysis_functions_folder);
end
old_folder = pwd;
folder_cleanup = onCleanup(@() cd(old_folder));
cd(folder);
trig = loadTrig(0);
clear folder_cleanup;
cd(old_folder);

info.trigger = struct();
info.trigger.file = files.trigger.path;
info.trigger.n_triggers = numel(trig);
info.trigger.first_sample = first_or_nan(trig);
info.trigger.last_sample = last_or_nan(trig);
info.trigger.FS = CFG.FS;

if info.trigger.n_triggers ~= info.experiment.n_trials
    info = add_problem(info, 'ERROR', sprintf( ...
        'Trigger count (%d) does not equal n_Trials (%d).', ...
        info.trigger.n_triggers, info.experiment.n_trials));
end

expected_rows = 1 + info.experiment.n_trials * ...
    info.experiment.simultaneous_stim;
if size(info.experiment.StimParams,1) ~= expected_rows
    info = add_problem(info, 'ERROR', sprintf( ...
        'StimParams has %d rows; expected %d.', ...
        size(info.experiment.StimParams,1), expected_rows));
end

if exist('Depth_s', 'file') ~= 2
    error('StageA:MissingDepthMap', ...
        'Depth_s was not found after adding %s', CFG.analysis_functions_folder);
end

% Depth_s calls read_Intan_RHS2000_file, which expects info.rhs in the
% current dataset folder. Change folder only for this call and always restore
% the user's original MATLAB folder, including when Depth_s throws an error.
old_folder = pwd;
folder_cleanup = onCleanup(@() cd(old_folder));
cd(folder);
depth_map = Depth_s(CFG.Electrode_Type);
clear folder_cleanup;
cd(old_folder);

depth_map = depth_map(:).';
info.channel_map = struct('depth_to_hardware', depth_map, ...
    'n_depth_channels',numel(depth_map));

if any(depth_map < 1) || any(depth_map > info.spike.n_channels)
    info = add_problem(info, 'ERROR', ...
        'Depth_s contains indices outside the sp_corr channel range.');
end

% Load small QC structures and final response structure.
info.bad_trials = qc_load_bad_trials(files.bad_trials);
info.bad_channels = qc_load_bad_channels(files.bad_channels);
info.responding = qc_load_responding(files.responding);

if info.bad_trials.file_found && ~info.bad_trials.is_global
    info = add_problem(info, 'ERROR', ...
        'BadTrials entries are not identical across recording channels.');
end

if require_qc
    if ~info.bad_trials.file_found
        info = add_problem(info, 'ERROR', ...
            'Required paired-data bad-trial file was not found.');
    end
    if ~info.bad_channels.file_found
        info = add_problem(info, 'ERROR', ...
            'Required paired-data bad-channel file was not found.');
    end
    if ~info.responding.file_found
        info = add_problem(info, 'ERROR', ...
            'Required paired-data RespondingChannels file was not found.');
    end
end

if strcmp(role, 'single')
    if ~info.bad_trials.file_found
        info.messages{end+1,1} = ...
            '[INFO][SINGLE] No bad-trial file; single trials remain unreviewed.';
    end
    if ~info.responding.file_found
        info.messages{end+1,1} = ...
            '[INFO][SINGLE] No responding file; paired-data labels will be used.';
    end
end
end

function n = numel_from_size(sz)
n = prod(double(sz));
end

function v = first_or_nan(x)
if isempty(x), v = NaN; else, v = x(1); end
end

function v = last_or_nan(x)
if isempty(x), v = NaN; else, v = x(end); end
end

function info = add_problem(info, severity, message)
info.messages{end+1,1} = sprintf('[%s][%s] %s', ...
    severity, upper(info.role), message);
if strcmp(severity, 'ERROR')
    info.validation_passed = false;
end
end
