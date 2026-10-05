%% ========================================================================
% RESPONDING-CHANNEL MANUAL CONFIRMATION
%
%   Edit Overrides below (one row per condition you're confirming this
%   run - you don't have to do every condition at once), then run.
%
% OVERRIDES TABLE
%   Overrides = {
%     % Set  Amp  PTD(ms)  ResetAllFirst  ForceRespond            ForceSilent
%       1,   5,   0,       true,          [3,4,6:11,13,15:16],    [];
%       1,   5,   3,       true,          [1:11,13,15,32:37],     [];
%   };
%
%   ResetAllFirst = true  : ignore whatever is already saved for this
%                           condition. The new final list = ForceRespond.
%   ResetAllFirst = false : start from whatever is already saved for this
%                           condition (empty, with a warning, if nothing
%                           was saved yet), then add ForceRespond and
%                           remove ForceSilent - a patch, so you can fix
%                           one or two channels without retyping the
%                           whole list.
%   ForceRespond and ForceSilent must not share a channel on the same
%   row - that is treated as a contradiction and raises an error.
%
%   PTD is ignored (leave as 0) for a single-pulse dataset
%   (simultaneous_stim == 1).
%
% OUTPUT
%   <base_name>_RespondingChannels.mat, in dataset_folder
%   Main saved variable: RespondingConfirmed (struct array, one row per
%   confirmed condition). A matching condition from a previous run is
%   REPLACED (not duplicated); conditions not mentioned in Overrides this
%   run are left untouched. A timestamped backup of the previous file is
%   kept.
%
% NOTE
%   This script does not load spike data - it only needs the experiment
%   file (to decode/validate conditions) and the channel count (to bound
%   channel numbers and compute distance).
% ========================================================================

clear;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_500_700um_SimSeq1';

% Electrode type (must match this dataset - probe differs by animal):
%   0 = rigid single-shank probe
%   1 = flexible single-shank probe
%   2 = four-shank flexible probe
%   3 = 64chn+32chn hybrid
Electrode_Type = 3;

% Physical channel count in this dataset's spike file (printed by
% RespondingChn_Raster_AllChn_v2.m / RespondingChn_ReferenceCriteria_v1.m
% as "Physical channels in spike file: N" when you ran them on the same
% spike_file). Needed here only to bound channel numbers via ChnMap -
% spike data itself is not loaded.
nPhysicalChannels = 96;

Condition_Tolerance = 1e-4;

%% --------------------- CONDITIONS BEING CONFIRMED NOW ------------------

Overrides = {
  % Set  Amp  PTD(ms)  ResetAllFirst  ForceRespond                        ForceSilent
    1,   5,   0,       true,          [1:50,52:64,65:72,86,90], [];
    1,   10,   0,      true,          [1:50,52:57,60:64,65,67,71,72,74:76,86,88,89,96], [];
    1,   5,   5,       true,          [1:50,52:64,66:68,70,72:74,80:82,89,91], [];
    1,   10,   5,      true,          [1:64,65:75,81:82,84,89,91:94], [];

    2,   5,   5,       true,          [1:50,52:58,60:64,65:67,69:77,82:94], [];
    2,   10,   5,      true,          [1:50,52:58,60:64,65:76,78:87,89:94], [];

    3,   5,   5,       true,          [1:50,52:58,60:64,65:69,71,74,75,77:93], [];
    3,   10,   5,      true,          [1:50,52:64,65:94], [];
    
    4,   5,   0,       true,          [1:50,52:58,60:64,65:66,69:74,86,90,92], [];
    4,   10,   0,      true,          [1:50,52:58,60:64,65:68,71:74,86,88,91:92], [];
    4,   5,   5,       true,          [1:50,52:64,65:67,70,72:75,79:86,90:92], [];
    4,   10,   5,      true,          [1:50,52:64,65:95], [];

    5,   5,   5,       true,          [1:50,52:64,65:68,70:82,86,88:92], [];
    5,   10,   5,      true,          [1:50,52:64,65:95], [];

    6,   5,   0,       true,          [1:50,52:55,57:64,65:67,71,83:84,92], [];
    6,   10,   0,      true,          [1:50,52:58,61:64,65:68,71:74,81:82,86,88:89,91:93], [];
    6,   5,   5,       true,          [1:50,52:64,65:67,70:75,77,80:82,86,91:92], [];
    6,   10,   5,      true,          [1:64,65:94], [];

    % 7,   5,   5,       true,          [1:50,52:61,63:64,65:67,71,74,83,85], [];
    % 7,   10,   5,      true,          [1:50,52:64,65:75,78:92], [];
    % 
    % 8,   5,   5,       true,          [2:50,52:58,60:61,63:64,65:66,74,83:91], [];
    % 8,   10,   5,      true,          [2:50,52:64,65:86,88:95], [];

};

Create_Backup_If_Output_Exists = true;

%% =========================== INITIAL CHECKS ===========================

if ~isfolder(dataset_folder)
    error('Dataset folder does not exist:\n%s',dataset_folder);
end
if isempty(Overrides)
    error('Overrides is empty - nothing to confirm.');
end
if size(Overrides,2) ~= 6
    error('Overrides must have 6 columns: Set, Amp, PTD(ms), ResetAllFirst, ForceRespond, ForceSilent.');
end
if ~isscalar(nPhysicalChannels) || nPhysicalChannels < 1
    error('nPhysicalChannels must be a positive integer.');
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('RESPONDING-CHANNEL MANUAL CONFIRMATION\n');
fprintf('============================================================\n');
fprintf('Dataset folder: %s\n',dataset_folder);
fprintf('Rows in Overrides: %d\n',size(Overrides,1));

%% ========================= ELECTRODE MAPPING ===========================

map_nums_plus = double(ChnMap(Electrode_Type,nPhysicalChannels));
map_nums_plus = map_nums_plus(:);
nMapPos = numel(map_nums_plus);

if nMapPos < 1
    error('ChnMap(%d,%d) returned an empty channel map.',Electrode_Type,nPhysicalChannels);
end
fprintf('Map positions: %d\n',nMapPos);

%% ==================== LOAD EXPERIMENT PARAMETERS =====================

cd(dataset_folder);

experiment_files = dir('*_exp_datafile_*.mat');
experiment_files = remove_metadata_and_backups(experiment_files);

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found in:\n%s',dataset_folder);
end
if numel(experiment_files) > 1
    error('Multiple current experiment files were found in:\n%s',dataset_folder);
end

experiment_file = experiment_files(1).name;

ExpLoad = load(experiment_file,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required_variables = {'StimParams','simultaneous_stim','E_MAP','n_Trials'};
for variable_index = 1:numel(required_variables)
    if ~isfield(ExpLoad,required_variables{variable_index})
        error('Experiment file is missing variable: %s',required_variables{variable_index});
    end
end

StimParams = ExpLoad.StimParams;
sim_stim = double(ExpLoad.simultaneous_stim);
E_MAP = ExpLoad.E_MAP;
n_Trials = double(ExpLoad.n_Trials);

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Stimulation events per trial (simultaneous_stim): %d\n',sim_stim);

%% =========================== DECODE AMPLITUDES ========================

amplitudes_all = cell2mat(StimParams(2:end,16));
amplitudes_all = double(amplitudes_all(:));

first_event_rows = 1:sim_stim:numel(amplitudes_all);
trialAmps = amplitudes_all(first_event_rows);
trialAmps = trialAmps(1:n_Trials);
trialAmps(trialAmps == -1) = 0;

[Amps,~,ampIdx] = unique(trialAmps(:));

%% ============================== DECODE PTDs =============================

if sim_stim >= 2
    PTD_all_us = cell2mat(StimParams(2:end,6));
    PTD_all_us = double(PTD_all_us(:));

    second_event_rows = 2:sim_stim:numel(PTD_all_us);
    trialPTD_us = PTD_all_us(second_event_rows);
    trialPTD_us = trialPTD_us(1:n_Trials);

    [PTDs,~,ptdIdx] = unique(trialPTD_us(:));
    PTDs_ms = PTDs/1000;
else
    PTDs_ms = 0;
    ptdIdx = ones(n_Trials,1);
end

%% ==================== DECODE STIMULATION CHANNEL ORDERS =================

stimNames = StimParams(2:end,1);
stimNames = stimNames(1:n_Trials*sim_stim);

[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));
if any(~isMapped)
    error('%d stimulation entries could not be mapped through E_MAP.',sum(~isMapped));
end

stimSeq = zeros(n_Trials,sim_stim);
for trial_id = 1:n_Trials
    rows_this_trial = (trial_id-1)*sim_stim+(1:sim_stim);
    mapped_channels = idx_all(rows_this_trial);
    mapped_channels = mapped_channels(mapped_channels > 0);
    stimSeq(trial_id,1:numel(mapped_channels)) = mapped_channels(:).';
end

[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

fprintf('Ordered stimulation sets: %d\n',nSets);
for si = 1:nSets
    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);
    fprintf('  Set %d: %s\n',si,format_order(stim_channels));
end

%% ================ LOAD EXISTING SAVED FILE (IF ANY) ====================

base_name = regexprep(experiment_file,'_exp_datafile_.*$','');
output_name = sprintf('%s_RespondingChannels.mat',base_name);
output_path = fullfile(dataset_folder,output_name);

if isfile(output_path)
    ExistingLoad = load(output_path,'RespondingConfirmed');
    if isfield(ExistingLoad,'RespondingConfirmed')
        RespondingConfirmed = ExistingLoad.RespondingConfirmed;
        fprintf('\nExisting saved file found (%d condition(s) already confirmed):\n%s\n', ...
            numel(RespondingConfirmed),output_path);
    else
        warning('Existing output file does not contain RespondingConfirmed - starting fresh.');
        RespondingConfirmed = empty_confirmed_struct();
    end
else
    fprintf('\nNo existing saved file - this run creates it.\n');
    RespondingConfirmed = empty_confirmed_struct();
end

%% ===================== VALIDATE + APPLY EACH OVERRIDE ROW ==============

nRows = size(Overrides,1);

for r = 1:nRows

    si           = Overrides{r,1};
    amp_value    = Overrides{r,2};
    ptd_value    = Overrides{r,3};
    resetFirst   = Overrides{r,4};
    forceRespond = Overrides{r,5};
    forceSilent  = Overrides{r,6};

    if isempty(si) || si < 1 || si > nSets || fix(si) ~= si
        error('Row %d: set %s is not a valid stimulation set (1-%d).',r,mat2str(si),nSets);
    end

    ai = find(abs(Amps-amp_value) < Condition_Tolerance,1);
    if isempty(ai)
        error('Row %d: amp %g uA does not match any decoded amplitude (%s).', ...
            r,amp_value,num2str(Amps(:).'));
    end

    if sim_stim >= 2
        pi = find(abs(PTDs_ms-ptd_value) < Condition_Tolerance,1);
        if isempty(pi)
            error('Row %d: ptd %g ms does not match any decoded PTD (%s).', ...
                r,ptd_value,num2str(PTDs_ms(:).'));
        end
        current_ptd_ms = PTDs_ms(pi);
    else
        pi = 1;
        current_ptd_ms = 0;
    end

    trials_this = find(combClass == si & ampIdx == ai & ptdIdx == pi);

    forceRespond = validate_channel_list(forceRespond,nMapPos,r,'ForceRespond');
    forceSilent  = validate_channel_list(forceSilent,nMapPos,r,'ForceSilent');

    overlap = intersect(forceRespond,forceSilent);
    if ~isempty(overlap)
        error(['Row %d: channel(s) %s appear in BOTH ForceRespond and ForceSilent - ' ...
            'contradictory input.'],r,num2str(overlap));
    end

    %% ---------------- FIND EXISTING SAVED ENTRY (IF ANY) ---------------

    match_idx = [];
    for k = 1:numel(RespondingConfirmed)
        if RespondingConfirmed(k).set == si && ...
                abs(RespondingConfirmed(k).amp-Amps(ai)) < Condition_Tolerance && ...
                abs(RespondingConfirmed(k).ptd_ms-current_ptd_ms) < Condition_Tolerance
            match_idx = k;
            break;
        end
    end

    if resetFirst || isempty(match_idx)
        if ~resetFirst && isempty(match_idx)
            warning(['Row %d: ResetAllFirst is false but no existing saved entry was found ' ...
                'for this condition - starting from an empty list.'],r);
        end
        baseline_channels = [];
    else
        baseline_channels = RespondingConfirmed(match_idx).channels;
    end

    final_channels = union(baseline_channels(:).',forceRespond);
    final_channels = setdiff(final_channels,forceSilent);
    final_channels = sort(unique(final_channels),'ascend');

    %% ---------------- BUILD RECORD --------------------------------------

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    dist_um = NaN;
    dist_note = '';
    if numel(stim_channels) == 2
        [dist_um,dist_note] = ChnPairDistance(stim_channels(1),stim_channels(2),Electrode_Type);
    end

    if sim_stim == 1
        stim_type = 'single';
    elseif abs(current_ptd_ms) < Condition_Tolerance
        stim_type = 'simultaneous';
    else
        stim_type = 'sequential';
    end

    R = empty_confirmed_struct();
    R(1).set = si;
    R(1).stimChannels = stim_channels;
    R(1).stimOrderLabel = format_order(stim_channels);
    R(1).stimType = stim_type;
    R(1).amp = Amps(ai);
    R(1).ptd_ms = current_ptd_ms;
    R(1).distance_um = dist_um;
    R(1).distance_note = dist_note;
    R(1).n_trials = numel(trials_this);
    R(1).n_total_channels = nMapPos;
    R(1).channels = final_channels;
    R(1).n_responding = numel(final_channels);
    R(1).confirmed_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

    if isempty(match_idx)
        RespondingConfirmed(end+1) = R; %#ok<SAGROW>
    else
        RespondingConfirmed(match_idx) = R;
    end

    fprintf('\nRow %d: Set %d (%s) | Amp %g uA | PTD %g ms | %s | ResetAllFirst=%d\n', ...
        r,si,R(1).stimOrderLabel,R(1).amp,R(1).ptd_ms,stim_type,resetFirst);
    fprintf('  Responding channels (%d/%d): %s\n', ...
        R(1).n_responding,nMapPos,channel_list_text(final_channels));
end

%% =========================== SAVE RESULT ==============================

if isfile(output_path) && Create_Backup_If_Output_Exists
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    backup_name = sprintf('%s_RespondingChannels_BACKUP_%s.mat',base_name,timestamp);
    backup_path = fullfile(dataset_folder,backup_name);
    copyfile(output_path,backup_path);
    fprintf('\nExisting output backed up to:\n%s\n',backup_path);
end

save(output_path,'RespondingConfirmed','-v7.3');

fprintf('\n============================================================\n');
fprintf('RESPONDING-CHANNEL MANUAL CONFIRMATION COMPLETE\n');
fprintf('Saved to: %s\n',output_path);
fprintf('Total conditions confirmed (all-time, this file): %d\n',numel(RespondingConfirmed));
fprintf('============================================================\n');

fprintf('\nAll confirmed conditions in this file:\n');
for k = 1:numel(RespondingConfirmed)
    R = RespondingConfirmed(k);
    fprintf('  Set %d (%s) | Amp %g uA | PTD %g ms | %s | Responding %d/%d\n', ...
        R.set,R.stimOrderLabel,R.amp,R.ptd_ms,R.stimType,R.n_responding,R.n_total_channels);
end

%% =========================== LOCAL FUNCTIONS ==========================

function files = remove_metadata_and_backups(files)
if isempty(files), return; end
keep = true(size(files));
for file_index = 1:numel(files)
    file_name = files(file_index).name;
    if startsWith(file_name,'._') || contains(file_name,'BACKUP','IgnoreCase',true)
        keep(file_index) = false;
    end
end
files = files(keep);
end

function label = format_order(stim_channels)
stim_channels = stim_channels(stim_channels > 0);
if isempty(stim_channels)
    label = '(none)';
elseif numel(stim_channels) == 1
    label = sprintf('Ch%d',stim_channels(1));
else
    labels = arrayfun(@(c) sprintf('Ch%d',c),stim_channels,'UniformOutput',false);
    label = strjoin(labels,' -> ');
end
end

function output = channel_list_text(values)
if isempty(values)
    output = '(none)';
else
    output = strtrim(sprintf('%d ',values));
end
end

function channels = validate_channel_list(channels,nMapPos,row_index,field_name)
if isempty(channels)
    channels = [];
    return;
end
channels = unique(double(channels(:).'));
invalid = ~isfinite(channels) | channels < 1 | channels > nMapPos | fix(channels) ~= channels;
if any(invalid)
    error('Row %d: %s contains invalid channel numbers (must be 1-%d): %s', ...
        row_index,field_name,nMapPos,num2str(channels(invalid)));
end
end

function R = empty_confirmed_struct()
R = struct('set',{},'stimChannels',{},'stimOrderLabel',{},'stimType',{}, ...
    'amp',{},'ptd_ms',{},'distance_um',{},'distance_note',{},'n_trials',{}, ...
    'n_total_channels',{},'channels',{},'n_responding',{},'confirmed_at',{});
end