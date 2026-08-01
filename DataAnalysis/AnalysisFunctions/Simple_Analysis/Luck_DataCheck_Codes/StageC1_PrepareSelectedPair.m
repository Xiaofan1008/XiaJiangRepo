%% Stage C1: prepare trigger-aligned QC data for one selected pair
% Reads confirmed Stage A/B results and sp_corr. Saves relative spike times
% and numerical QC results only. Does not modify experiment files.

clear;
clc;

%% ========================= USER SETTINGS ==============================
stage_b_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/stage_a_check_v5/check_results/DX014_D4_StageB_PairSummary.mat'];

pair_index = 1;

% Channel selection: 'responding', 'all_common', or 'manual'.
channel_mode = 'responding';
manual_depth_channels = [];

% Empty means every structurally complete amplitude saved by Stage B.
amplitudes_to_check = [];

% Time windows in milliseconds. Count windows use [start, end).
raster_win_ms   = [-50 80];
baseline_win_ms = [-50 -5];
count_win_ms    = [2 20];

% Single data have not yet been reviewed. After inspecting Stage C2 rasters,
% enter global Single trial IDs here and set review_complete to true.
single_bad_trials = [];
single_review_complete = false;

% Expected physical geometry of the selected stimulation pair.
require_same_stim_shank = true;
expected_stim_spacing_um = 200;
geometry_tolerance_um = 0.01;

%% ======================= LOAD STAGE B/A ===============================
if ~isfile(stage_b_file)
    error('StageC1:MissingStageBFile', ...
        'Stage B result not found:\n%s', stage_b_file);
end

SB = load(stage_b_file, 'StageB');
if ~isfield(SB,'StageB')
    error('StageC1:InvalidStageBFile', ...
        'The selected file does not contain StageB: %s', stage_b_file);
end
StageB = SB.StageB;

if ~isfield(StageB,'source_stage_a_file') || ...
        ~isfile(StageB.source_stage_a_file)
    error('StageC1:MissingStageAFile', ...
        'The Stage A source linked by Stage B was not found.');
end

SA = load(StageB.source_stage_a_file, 'StageA');
if ~isfield(SA,'StageA') || ~SA.StageA.validation_passed
    error('StageC1:InvalidStageAResult', ...
        'The linked Stage A result is missing or did not pass.');
end
StageA = SA.StageA;

if pair_index < 1 || pair_index > numel(StageB.pairs) || ...
        mod(pair_index,1) ~= 0
    error('StageC1:InvalidPairIndex', ...
        'pair_index must be an integer from 1 to %d.', numel(StageB.pairs));
end
PairB = StageB.pairs(pair_index);
if ~PairB.is_structurally_complete
    error('StageC1:IncompletePair', ...
        'Pair %d (%s) is not structurally complete.', pair_index, PairB.key);
end

Dsingle = StageA.datasets.single;
Dsim = StageA.datasets.sim;
Dseq = StageA.datasets.seq;
Esingle = Dsingle.experiment;
Esim = Dsim.experiment;
Eseq = Dseq.experiment;
FS = StageA.config.FS;
Electrode_Type = StageA.config.Electrode_Type;

% Add only the established analysis tree needed for loadTrig and ProbeMAP.
if ~isfolder(StageA.config.analysis_functions_folder)
    error('StageC1:MissingAnalysisFunctions', ...
        'AnalysisFunctions folder not found: %s', ...
        StageA.config.analysis_functions_folder);
end
addpath(genpath(StageA.config.analysis_functions_folder), '-end');

%% ========================= WRITE SAFETY ===============================
check_results_folder = fileparts(stage_b_file);
dataset_folders = {StageA.config.single_folder, StageA.config.sim_folder, ...
    StageA.config.seq_folder};
for k = 1:numel(dataset_folders)
    if path_is_inside(check_results_folder, dataset_folders{k})
        error('StageC1:UnsafeOutputFolder', ...
            ['Refusing to write because check_results is inside an ' ...
             'experiment folder.\nOutput: %s\nExperiment: %s'], ...
            check_results_folder, dataset_folders{k});
    end
end

%% ==================== VALIDATE USER SETTINGS ==========================
validate_window(raster_win_ms, 'raster_win_ms');
validate_window(baseline_win_ms, 'baseline_win_ms');
validate_window(count_win_ms, 'count_win_ms');
if baseline_win_ms(1) < raster_win_ms(1) || ...
        baseline_win_ms(2) > raster_win_ms(2) || ...
        count_win_ms(1) < raster_win_ms(1) || ...
        count_win_ms(2) > raster_win_ms(2)
    error('StageC1:WindowOutsideRaster', ...
        'Baseline and count windows must lie inside raster_win_ms.');
end

single_bad_trials = unique(double(single_bad_trials(:).'));
if ~isscalar(single_review_complete) || ...
        ~(islogical(single_review_complete) || isnumeric(single_review_complete))
    error('StageC1:InvalidSingleReviewFlag', ...
        'single_review_complete must be one true/false value.');
end
single_review_complete = logical(single_review_complete);
if any(single_bad_trials < 1) || ...
        any(single_bad_trials > Esingle.n_trials) || ...
        any(mod(single_bad_trials,1) ~= 0)
    error('StageC1:InvalidSingleBadTrials', ...
        'single_bad_trials contains an invalid Single trial ID.');
end

complete_amps = PairB.complete_amplitudes_uA(:).';
if isempty(amplitudes_to_check)
    selected_amps = complete_amps;
else
    requested_amps = unique(double(amplitudes_to_check(:).'));
    selected_amps = zeros(size(requested_amps));
    for ia = 1:numel(requested_amps)
        idx = find(abs(complete_amps-requested_amps(ia)) <= 1e-6, 1);
        if isempty(idx)
            error('StageC1:UnavailableAmplitude', ...
                'Amplitude %.6g uA is not complete for pair %s.', ...
                requested_amps(ia), PairB.key);
        end
        selected_amps(ia) = complete_amps(idx);
    end
end
if isempty(selected_amps)
    error('StageC1:NoAmplitudes', 'No complete amplitudes were selected.');
end

switch lower(strtrim(channel_mode))
    case 'responding'
        selected_depth_channels = PairB.responding_union_depth_channels;
    case 'all_common'
        selected_depth_channels = PairB.common_good_depth_channels;
    case 'manual'
        selected_depth_channels = unique(double(manual_depth_channels(:).'), ...
            'stable');
        if isempty(selected_depth_channels)
            error('StageC1:EmptyManualChannels', ...
                'manual_depth_channels is empty.');
        end
        if ~all(ismember(selected_depth_channels, ...
                PairB.common_good_depth_channels))
            error('StageC1:InvalidManualChannel', ...
                'Every manual channel must be a common good depth channel.');
        end
    otherwise
        error('StageC1:InvalidChannelMode', ...
            'channel_mode must be responding, all_common, or manual.');
end
selected_depth_channels = selected_depth_channels(:).';
if isempty(selected_depth_channels)
    error('StageC1:NoSelectedChannels', ...
        'The selected channel mode produced no recording channels.');
end

%% ==================== STIMULATION/CHANNEL MAPPING =====================
nDepth = Dsim.channel_map.n_depth_channels;
stimA = map_stimulation_contact(PairB.electrode_A, Electrode_Type, ...
    nDepth, Dsim.channel_map.depth_to_hardware);
stimB = map_stimulation_contact(PairB.electrode_B, Electrode_Type, ...
    nDepth, Dsim.channel_map.depth_to_hardware);

% The experimental E_MAP positions and ProbeMAP positions should agree.
if isempty(PairB.sim_electrode_indices.A) || ...
        PairB.sim_electrode_indices.A ~= stimA.depth_channel
    error('StageC1:StimAMappingMismatch', ...
        ['Stim A mapping mismatch for %s. Stage B E_MAP index=%s; ' ...
         'ProbeMAP depth channel=%d.'], PairB.electrode_A, ...
        mat2str(PairB.sim_electrode_indices.A), stimA.depth_channel);
end
if isempty(PairB.sim_electrode_indices.B) || ...
        PairB.sim_electrode_indices.B ~= stimB.depth_channel
    error('StageC1:StimBMappingMismatch', ...
        ['Stim B mapping mismatch for %s. Stage B E_MAP index=%s; ' ...
         'ProbeMAP depth channel=%d.'], PairB.electrode_B, ...
        mat2str(PairB.sim_electrode_indices.B), stimB.depth_channel);
end

stim_pair_distance_um = euclidean_distance(stimA.x_um, stimA.y_um, ...
    stimB.x_um, stimB.y_um);
stim_pair_same_shank = stimA.shank == stimB.shank;
if require_same_stim_shank && ~stim_pair_same_shank
    error('StageC1:StimPairShankMismatch', ...
        'Mapped stimulation contacts are not on the same shank.');
end
if isfinite(expected_stim_spacing_um) && ...
        abs(stim_pair_distance_um-expected_stim_spacing_um) > ...
        geometry_tolerance_um
    error('StageC1:StimPairSpacingMismatch', ...
        ['Mapped stimulation spacing is %.3f um; expected %.3f um. ' ...
         'Check the electrode map before continuing.'], ...
        stim_pair_distance_um, expected_stim_spacing_um);
end

stim_midpoint = struct('x_um',(stimA.x_um+stimB.x_um)/2, ...
    'y_um',(stimA.y_um+stimB.y_um)/2);

ChannelInfo = build_channel_table(selected_depth_channels, ...
    Dsingle.channel_map.depth_to_hardware, ...
    Dsim.channel_map.depth_to_hardware, ...
    Dseq.channel_map.depth_to_hardware, Electrode_Type, ...
    stimA, stimB, stim_midpoint, PairB.responding_union_depth_channels);

%% ==================== LOAD TRIGGERS & SPIKES ==========================
fprintf('\n============================================================\n');
fprintf('STAGE C1: PREPARE SELECTED-PAIR QC DATA\n');
fprintf('============================================================\n');
fprintf('Pair %d: %s\n', pair_index, PairB.key);
fprintf('Amplitudes: %s uA\n', num2str(selected_amps));
fprintf('Selected depth channels (%s): %s\n', ...
    channel_mode, num2str(selected_depth_channels));
fprintf('Loading relative spike data only; waveform columns are not retained.\n');

trig_single = load_triggers_read_only(Dsingle.folder, Esingle.n_trials);
if strcmp(Dsim.folder,Dseq.folder)
    trig_sim = load_triggers_read_only(Dsim.folder, Esim.n_trials);
    trig_seq = trig_sim;
else
    trig_sim = load_triggers_read_only(Dsim.folder, Esim.n_trials);
    trig_seq = load_triggers_read_only(Dseq.folder, Eseq.n_trials);
end

hw_single = Dsingle.channel_map.depth_to_hardware(selected_depth_channels);
hw_sim = Dsim.channel_map.depth_to_hardware(selected_depth_channels);
hw_seq = Dseq.channel_map.depth_to_hardware(selected_depth_channels);

fprintf('Loading %d selected channels from Single sp_corr...\n', ...
    numel(selected_depth_channels));
spike_times_single = load_selected_spike_times( ...
    Dsingle.files.spike.path, hw_single);

fprintf('Loading %d selected channels from Sim sp_corr...\n', ...
    numel(selected_depth_channels));
spike_times_sim = load_selected_spike_times(Dsim.files.spike.path, hw_sim);

same_pair_spike_source = strcmp(Dsim.files.spike.path,Dseq.files.spike.path) && ...
    isequal(hw_sim,hw_seq);
if same_pair_spike_source
    spike_times_seq = spike_times_sim;
    fprintf('Seq uses the same sp_corr source and mapping as Sim; data reused.\n');
else
    fprintf('Loading %d selected channels from Seq sp_corr...\n', ...
        numel(selected_depth_channels));
    spike_times_seq = load_selected_spike_times(Dseq.files.spike.path, hw_seq);
end

%% ==================== EXTRACT & CALCULATE =============================
condition_order = {'A','B','AB','A_to_B','B_to_A'};
condition_labels = {PairB.electrode_A, PairB.electrode_B, ...
    [PairB.electrode_A '+' PairB.electrode_B], ...
    [PairB.electrode_A '->' PairB.electrode_B], ...
    [PairB.electrode_B '->' PairB.electrode_A]};
condition_sources = {'single','single','sim','seq','seq'};

nSelected = numel(selected_depth_channels);
nAmp = numel(selected_amps);
nCondition = numel(condition_order);
Data = repmat(empty_condition_data(), nSelected, nAmp, nCondition);
summary_rows = cell(nSelected*nAmp*nCondition, 18);
summary_row = 0;
ConditionTrialCounts = zeros(nAmp,nCondition);

for ia = 1:nAmp
    amp = selected_amps(ia);
    amp_idx_B = find(abs([PairB.amplitudes.amplitude_uA]-amp) <= 1e-6, 1);
    if isempty(amp_idx_B) || ~PairB.amplitudes(amp_idx_B).is_complete
        error('StageC1:StageBAmplitudeMismatch', ...
            'Stage B no longer reports %.6g uA as complete.', amp);
    end
    AmpB = PairB.amplitudes(amp_idx_B);

    trial_sets = {AmpB.clean_trial_ids.A, AmpB.clean_trial_ids.B, ...
        AmpB.clean_trial_ids.AB, AmpB.clean_trial_ids.A_to_B, ...
        AmpB.clean_trial_ids.B_to_A};
    trial_sets{1} = setdiff(trial_sets{1},single_bad_trials,'stable');
    trial_sets{2} = setdiff(trial_sets{2},single_bad_trials,'stable');

    for icond = 1:nCondition
        trial_ids = trial_sets{icond}(:).';
        if isempty(trial_ids)
            error('StageC1:EmptyConditionAfterExclusion', ...
                'No trials remain for %.6g uA, condition %s.', ...
                amp, condition_order{icond});
        end
        ConditionTrialCounts(ia,icond) = numel(trial_ids);

        switch condition_sources{icond}
            case 'single'
                trig = trig_single;
                spike_source = spike_times_single;
            case 'sim'
                trig = trig_sim;
                spike_source = spike_times_sim;
            case 'seq'
                trig = trig_seq;
                spike_source = spike_times_seq;
        end

        if any(trial_ids > numel(trig))
            error('StageC1:TrialExceedsTriggerCount', ...
                'A Stage B trial ID exceeds the trigger count.');
        end

        for ichSel = 1:nSelected
            R = extract_condition(spike_source{ichSel}, trig, trial_ids, FS, ...
                raster_win_ms, baseline_win_ms, count_win_ms);
            R.condition = condition_order{icond};
            R.condition_label = condition_labels{icond};
            R.source_dataset = condition_sources{icond};
            R.amplitude_uA = amp;
            R.depth_channel = selected_depth_channels(ichSel);
            R.second_pulse_ms = second_pulse_time(condition_order{icond}, ...
                StageB.config.seq_PTD_ms);
            Data(ichSel,ia,icond) = R;

            summary_row = summary_row + 1;
            CI = ChannelInfo(ichSel,:);
            summary_rows(summary_row,:) = { ...
                selected_depth_channels(ichSel), CI.SpCellIndex_Sim, ...
                CI.Shank, CI.LocalSite, CI.MinimumDistance_um, amp, ...
                condition_order{icond}, numel(trial_ids), ...
                R.summary.mean_post_count, R.summary.sd_post_count, ...
                R.summary.sem_post_count, R.summary.mean_post_rate_hz, ...
                R.summary.response_probability, ...
                R.summary.mean_baseline_rate_hz, ...
                R.summary.mean_evoked_count, ...
                R.summary.mean_first_spike_latency_ms, ...
                min(trial_ids), max(trial_ids)};
        end
    end
end

summary_names = {'DepthChannel','SpCellIndex','Shank','LocalSite', ...
    'MinimumDistance_um','Amplitude_uA','Condition','NTrials', ...
    'MeanPostCount','SDPostCount','SEMPostCount','MeanPostRate_Hz', ...
    'ResponseProbability','MeanBaselineRate_Hz','MeanEvokedCount', ...
    'MeanFirstSpikeLatency_ms','MinSourceTrialID','MaxSourceTrialID'};
TrialSummary = cell2table(summary_rows(1:summary_row,:), ...
    'VariableNames',summary_names);

ConditionTrialCountTable = array2table(ConditionTrialCounts, ...
    'VariableNames',condition_order, 'RowNames', ...
    cellstr(compose('%g_uA',selected_amps(:))));

%% ============================= SAVE ==================================
StageC1 = struct();
StageC1.created_on = datetime('now');
StageC1.source_stage_b_file = stage_b_file;
StageC1.source_stage_a_file = StageB.source_stage_a_file;
StageC1.pair_index = pair_index;
StageC1.pair_key = PairB.key;
StageC1.FS = FS;
StageC1.Electrode_Type = Electrode_Type;
StageC1.electrode_A = PairB.electrode_A;
StageC1.electrode_B = PairB.electrode_B;
StageC1.amplitudes_uA = selected_amps;
StageC1.condition_order = condition_order;
StageC1.condition_labels = condition_labels;
StageC1.selected_depth_channels = selected_depth_channels;
StageC1.channel_mode = channel_mode;
StageC1.windows = struct('raster_ms',raster_win_ms, ...
    'baseline_ms',baseline_win_ms, 'count_ms',count_win_ms, ...
    'interval_convention','[start,end), except raster includes end');
StageC1.single_qc = struct('bad_trial_ids',single_bad_trials, ...
    'review_complete',logical(single_review_complete));
StageC1.stimulation = struct('A',stimA, 'B',stimB, ...
    'same_shank',stim_pair_same_shank, ...
    'pair_distance_um',stim_pair_distance_um, ...
    'midpoint',stim_midpoint, ...
    'expected_pair_distance_um',expected_stim_spacing_um);
StageC1.ChannelInfo = ChannelInfo;
StageC1.ConditionTrialCounts = ConditionTrialCountTable;
StageC1.TrialSummary = TrialSummary;
StageC1.Data = Data;
StageC1.waveforms_retained = false;
StageC1.absolute_spike_times_retained = false;
StageC1.trigger_times_retained = false;
StageC1.source_spike_time_unit = 'milliseconds';
StageC1.experiment_files_modified = false;

[~, stage_b_stem] = fileparts(stage_b_file);
dataset_tag = strrep(stage_b_stem,'_StageB_PairSummary','');
output_mat = fullfile(check_results_folder, sprintf( ...
    '%s_StageC1_Pair%d_QCData.mat',dataset_tag,pair_index));
output_txt = fullfile(check_results_folder, sprintf( ...
    '%s_StageC1_Pair%d_Report.txt',dataset_tag,pair_index));

for k = 1:numel(dataset_folders)
    if path_is_inside(fileparts(output_mat),dataset_folders{k}) || ...
            path_is_inside(fileparts(output_txt),dataset_folders{k})
        error('StageC1:UnsafeWriteBlocked', ...
            'Stage C1 output resolved inside an experiment folder.');
    end
end

% Remove large absolute-time source arrays from the workspace before saving.
clear spike_times_single spike_times_sim spike_times_seq ...
    trig_single trig_sim trig_seq;

save(output_mat,'StageC1','-v7.3');
write_stage_c1_report(output_txt,StageC1);

fprintf('\n============================================================\n');
fprintf('STAGE C1 COMPLETE\n');
fprintf('MAT result: %s\n',output_mat);
fprintf('Text report: %s\n',output_txt);
fprintf('Experiment files modified: NO\n');
fprintf('============================================================\n\n');

%% =========================== FUNCTIONS ================================
function validate_window(value,name)
if ~isnumeric(value) || numel(value) ~= 2 || ...
        any(~isfinite(value)) || value(2) <= value(1)
    error('StageC1:InvalidWindow','%s must be [start end], with end>start.',name);
end
end

function stim = map_stimulation_contact(native_name,electrode_type,nDepth,depth_map)
if exist('ProbeMAP','file') ~= 2
    error('StageC1:MissingProbeMAP','ProbeMAP was not found on the MATLAB path.');
end
probe_map = ProbeMAP;
switch electrode_type
    case 0
        map_column = 3;
    case 1
        map_column = 5;
    case 2
        map_column = 6;
    otherwise
        error('StageC1:InvalidElectrodeType','Unknown Electrode_Type.');
end
if size(probe_map,1) < nDepth+1 || size(probe_map,2) < map_column
    error('StageC1:ProbeMapTooSmall','ProbeMAP does not cover all depth channels.');
end
map_names = probe_map(2:nDepth+1,map_column);
matches = find(strcmp(map_names,native_name));
if numel(matches) ~= 1
    error('StageC1:StimNameMappingFailed', ...
        'Expected one ProbeMAP match for %s; found %d.',native_name,numel(matches));
end
depth_channel = matches(1);
sp_cell_index = depth_map(depth_channel);

expected_hardware = native_name_to_matlab_index(native_name);
if sp_cell_index ~= expected_hardware
    error('StageC1:HardwareIndexMismatch', ...
        ['%s mapped to sp_corr cell %d through Depth_s, but its native ' ...
         'name implies MATLAB index %d.'],native_name,sp_cell_index, ...
        expected_hardware);
end

[shank,local_site,x_um,y_um] = depth_channel_coordinates( ...
    depth_channel,electrode_type,nDepth);
stim = struct('native_name',native_name, ...
    'depth_channel',depth_channel, ...
    'sp_cell_index',sp_cell_index, ...
    'shank',shank, 'local_site',local_site, ...
    'x_um',x_um, 'y_um',y_um);
end

function index = native_name_to_matlab_index(native_name)
tokens = regexp(native_name,'^([A-D])-([0-9]+)$','tokens','once');
if isempty(tokens)
    error('StageC1:InvalidNativeName','Cannot parse native name %s.',native_name);
end
bank = double(tokens{1})-double('A');
native_number = str2double(tokens{2});
index = bank*32 + native_number + 1;
end

function [shank,local_site,x_um,y_um] = depth_channel_coordinates(ich,type,nDepth)
vertical_pitch = 50;
inter_shank_pitch = 200;
if type == 0 || type == 1
    if ich < 1 || ich > 32
        error('StageC1:SingleShankChannelRange', ...
            'Single-shank depth channel must be from 1 to 32.');
    end
    shank = 1;
    local_site = ich;
    x_um = 0;
    y_um = (local_site-1)*vertical_pitch;
elseif type == 2
    if nDepth > 64 || ich < 1 || ich > 64
        error('StageC1:FourShankChannelRange', ...
            'The supplied four-shank geometry supports depth channels 1 to 64.');
    end
    if ich <= 16
        shank = 1; local_site = ich;
    elseif ich <= 32
        shank = 4; local_site = ich-16;
    elseif ich <= 48
        shank = 2; local_site = ich-32;
    else
        shank = 3; local_site = ich-48;
    end
    x_um = (shank-1)*inter_shank_pitch;
    y_um = (local_site-1)*vertical_pitch;
else
    error('StageC1:InvalidElectrodeType','Unknown Electrode_Type.');
end
end

function T = build_channel_table(depth_channels,hw_single_map,hw_sim_map, ...
        hw_seq_map,electrode_type,stimA,stimB,midpoint,responding_union)
n = numel(depth_channels);
DepthChannel = depth_channels(:);
SpCellIndex_Single = hw_single_map(depth_channels);
SpCellIndex_Sim = hw_sim_map(depth_channels);
SpCellIndex_Seq = hw_seq_map(depth_channels);
SpCellIndex_Single = SpCellIndex_Single(:);
SpCellIndex_Sim = SpCellIndex_Sim(:);
SpCellIndex_Seq = SpCellIndex_Seq(:);
Shank = zeros(n,1);
LocalSite = zeros(n,1);
X_um = zeros(n,1);
Y_um = zeros(n,1);
DistanceToA_um = zeros(n,1);
DistanceToB_um = zeros(n,1);
MinimumDistance_um = zeros(n,1);
DistanceToMidpoint_um = zeros(n,1);
IsRespondingUnion = ismember(DepthChannel,responding_union(:));

for k = 1:n
    [Shank(k),LocalSite(k),X_um(k),Y_um(k)] = ...
        depth_channel_coordinates(DepthChannel(k),electrode_type, ...
        numel(hw_sim_map));
    DistanceToA_um(k) = euclidean_distance(X_um(k),Y_um(k),stimA.x_um,stimA.y_um);
    DistanceToB_um(k) = euclidean_distance(X_um(k),Y_um(k),stimB.x_um,stimB.y_um);
    MinimumDistance_um(k) = min(DistanceToA_um(k),DistanceToB_um(k));
    DistanceToMidpoint_um(k) = euclidean_distance( ...
        X_um(k),Y_um(k),midpoint.x_um,midpoint.y_um);
end

T = table(DepthChannel,SpCellIndex_Single,SpCellIndex_Sim, ...
    SpCellIndex_Seq,Shank,LocalSite,X_um,Y_um,DistanceToA_um, ...
    DistanceToB_um,MinimumDistance_um,DistanceToMidpoint_um, ...
    IsRespondingUnion);
end

function distance = euclidean_distance(x1,y1,x2,y2)
distance = sqrt((x1-x2).^2 + (y1-y2).^2);
end

function trig = load_triggers_read_only(folder,expected_count)
if exist('loadTrig','file') ~= 2
    error('StageC1:MissingLoadTrig','loadTrig was not found on the MATLAB path.');
end
old_folder = pwd;
restore_folder = onCleanup(@() cd(old_folder));
cd(folder);
if isempty(dir('*.trig.dat'))
    error('StageC1:MissingTriggerFile', ...
        'No trig.dat file exists in %s. Stage C1 will not create one.',folder);
end
trig = loadTrig(0);
clear restore_folder;
cd(old_folder);
if numel(trig) ~= expected_count
    error('StageC1:TriggerCountMismatch', ...
        'Trigger count %d does not match expected trial count %d in %s.', ...
        numel(trig),expected_count,folder);
end
end

function selected = load_selected_spike_times(spike_file,hardware_indices)
meta = whos('-file',spike_file,'sp_corr');
if isempty(meta)
    error('StageC1:MissingSpCorr','sp_corr is missing from %s.',spike_file);
end
hardware_indices = hardware_indices(:).';
if any(hardware_indices < 1) || any(hardware_indices > prod(meta.size))
    error('StageC1:HardwareIndexOutOfRange', ...
        'Requested sp_corr channel is outside the saved cell array.');
end

raw_cells = [];
partial_ok = false;
try
    M = matfile(spike_file);
    if meta.size(1) == 1
        raw_cells = M.sp_corr(1,hardware_indices);
    elseif meta.size(2) == 1
        raw_cells = M.sp_corr(hardware_indices,1);
    else
        error('StageC1:UnexpectedSpCorrShape', ...
            'sp_corr must be a row or column cell array.');
    end
    partial_ok = iscell(raw_cells);
catch ME
    warning('StageC1:PartialLoadFailed', ...
        ['Selective matfile loading failed (%s). Stage C1 will load the ' ...
         'complete sp_corr variable for this dataset.'],ME.message);
end

if ~partial_ok
    S = load(spike_file,'sp_corr');
    raw_cells = S.sp_corr(hardware_indices);
end
raw_cells = raw_cells(:).';

selected = cell(size(raw_cells));
for k = 1:numel(raw_cells)
    matrix = raw_cells{k};
    if isempty(matrix)
        selected{k} = zeros(0,1);
    elseif ~isnumeric(matrix) || size(matrix,2) < 1
        error('StageC1:InvalidSpikeMatrix', ...
            'sp_corr{%d} is not a valid numeric spike matrix.',hardware_indices(k));
    else
        selected{k} = double(matrix(:,1));
    end
end
end

function R = extract_condition(abs_spike_times,trig,trial_ids,FS, ...
        raster_win,baseline_win,count_win)
nTrials = numel(trial_ids);
spike_times_ms = cell(nTrials,1);
baseline_count = zeros(nTrials,1);
baseline_rate_hz = zeros(nTrials,1);
post_count = zeros(nTrials,1);
post_rate_hz = zeros(nTrials,1);
has_response = false(nTrials,1);
first_spike_latency_ms = nan(nTrials,1);
evoked_count = zeros(nTrials,1);
base_duration_s = diff(baseline_win)/1000;
post_duration_s = diff(count_win)/1000;

for k = 1:nTrials
    tr = trial_ids(k);
    t0_ms = double(trig(tr))/FS*1000;
    keep = abs_spike_times >= t0_ms+raster_win(1) & ...
        abs_spike_times <= t0_ms+raster_win(2);
    relative = abs_spike_times(keep)-t0_ms;
    spike_times_ms{k} = relative(:);

    base_mask = relative >= baseline_win(1) & relative < baseline_win(2);
    post_mask = relative >= count_win(1) & relative < count_win(2);
    baseline_count(k) = sum(base_mask);
    baseline_rate_hz(k) = baseline_count(k)/base_duration_s;
    post_count(k) = sum(post_mask);
    post_rate_hz(k) = post_count(k)/post_duration_s;
    has_response(k) = post_count(k) > 0;
    if has_response(k)
        post_times = relative(post_mask);
        first_spike_latency_ms(k) = min(post_times);
    end
    evoked_count(k) = post_count(k)-baseline_rate_hz(k)*post_duration_s;
end

R = empty_condition_data();
R.source_trial_ids = trial_ids(:);
R.spike_times_ms = spike_times_ms;
R.trial_metrics = struct('baseline_count',baseline_count, ...
    'baseline_rate_hz',baseline_rate_hz, 'post_count',post_count, ...
    'post_rate_hz',post_rate_hz, 'has_response',has_response, ...
    'first_spike_latency_ms',first_spike_latency_ms, ...
    'evoked_count',evoked_count);
R.summary = struct( ...
    'n_trials',nTrials, ...
    'mean_post_count',mean(post_count), ...
    'sd_post_count',std(post_count,0), ...
    'sem_post_count',std(post_count,0)/sqrt(nTrials), ...
    'mean_post_rate_hz',mean(post_rate_hz), ...
    'response_probability',mean(has_response), ...
    'mean_baseline_rate_hz',mean(baseline_rate_hz), ...
    'median_baseline_rate_hz',median(baseline_rate_hz), ...
    'mean_evoked_count',mean(evoked_count), ...
    'mean_first_spike_latency_ms',mean(first_spike_latency_ms,'omitnan'));
end

function value = second_pulse_time(condition,seq_ptd)
if strcmp(condition,'A_to_B') || strcmp(condition,'B_to_A')
    value = seq_ptd;
else
    value = NaN;
end
end

function R = empty_condition_data()
R = struct('condition','', 'condition_label','', 'source_dataset','', ...
    'amplitude_uA',NaN, 'depth_channel',NaN, 'second_pulse_ms',NaN, ...
    'source_trial_ids',[], 'spike_times_ms',{{}}, ...
    'trial_metrics',struct(), 'summary',struct());
end

function write_stage_c1_report(filename,C)
fid = fopen(filename,'w');
if fid < 0
    error('StageC1:ReportWriteFailed','Could not write report: %s',filename);
end
cleanup = onCleanup(@() fclose(fid));

dual_print(fid,'\n################ COMPLETE STAGE C1 REPORT ################\n');
dual_print(fid,'Created: %s\n',char(string(C.created_on)));
dual_print(fid,'Pair: %s\n',C.pair_key);
dual_print(fid,'Amplitudes: %s uA\n',num2str(C.amplitudes_uA));
dual_print(fid,'Selected depth channels: %s\n',num2str(C.selected_depth_channels));
dual_print(fid,'Single review complete: %s\n', ...
    char(string(C.single_qc.review_complete)));
dual_print(fid,'Single excluded trials: %s\n',vector_text(C.single_qc.bad_trial_ids));
dual_print(fid,'Waveforms retained: NO\n');
dual_print(fid,'Absolute spike/trigger times retained: NO\n');
dual_print(fid,'Experiment files modified: NO\n\n');

dual_print(fid,'[STIMULATION MAPPING]\n');
dual_print(fid,'%s -> depth %d, sp_corr cell %d, shank %d, site %d, (%.1f, %.1f) um\n', ...
    C.stimulation.A.native_name,C.stimulation.A.depth_channel, ...
    C.stimulation.A.sp_cell_index,C.stimulation.A.shank, ...
    C.stimulation.A.local_site,C.stimulation.A.x_um,C.stimulation.A.y_um);
dual_print(fid,'%s -> depth %d, sp_corr cell %d, shank %d, site %d, (%.1f, %.1f) um\n', ...
    C.stimulation.B.native_name,C.stimulation.B.depth_channel, ...
    C.stimulation.B.sp_cell_index,C.stimulation.B.shank, ...
    C.stimulation.B.local_site,C.stimulation.B.x_um,C.stimulation.B.y_um);
dual_print(fid,'Same shank: %s\n',char(string(C.stimulation.same_shank)));
dual_print(fid,'Stimulation-pair distance: %.1f um\n\n', ...
    C.stimulation.pair_distance_um);

dual_print(fid,'[CONDITION TRIAL COUNTS]\n%s\n', ...
    evalc('disp(C.ConditionTrialCounts)'));
dual_print(fid,'[SELECTED CHANNELS AND DISTANCES]\n%s\n', ...
    evalc('disp(C.ChannelInfo)'));

dual_print(fid,'############## END COMPLETE STAGE C1 REPORT ##############\n');
end

function text_value = vector_text(values)
if isempty(values)
    text_value = '(none)';
else
    text_value = strtrim(num2str(values(:).'));
end
end

function dual_print(fid,varargin)
fprintf(1,varargin{:});
fprintf(fid,varargin{:});
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
