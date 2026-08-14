%% Stage D1: export a simple modelling package for one selected pair
% Reads the confirmed Stage C1 result and its linked Stage A/B metadata.
% Writes one self-contained MAT file. Original experiment and QC files are
% never modified.

clear;
clc;

%% ========================= USER SETTINGS ==============================
stage_c1_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/stage_a_check_v5/check_results/DX014_D4_StageC1_Pair1_QCData.mat'];

% Leave empty to create a folder named Luke_ModelPackage beside C1.
export_folder = '';

% Leave empty to generate a name from the C1 filename and pair identity.
output_filename = '';

% Balanced subsets are optional views of the complete exported trials.
% The original trials are never removed. A fixed seed makes the subset
% identical every time this script is run on the same input.
balancing_random_seed = 20260802;

% Existing exports are protected unless this is deliberately changed.
allow_overwrite = false;

%% ============================ LOAD ====================================
if ~isfile(stage_c1_file)
    error('StageD1:MissingC1','Stage C1 file not found:\n%s',stage_c1_file);
end

S = load(stage_c1_file,'StageC1');
if ~isfield(S,'StageC1')
    error('StageD1:InvalidC1','The selected file does not contain StageC1.');
end
C = S.StageC1;

required_c1 = {'Data','ChannelInfo','amplitudes_uA','condition_order', ...
    'windows','stimulation','source_stage_b_file','pair_index','pair_key', ...
    'electrode_A','electrode_B','FS','Electrode_Type','single_qc'};
require_fields(C,required_c1,'StageC1');

if ~isfile(C.source_stage_b_file)
    error('StageD1:MissingStageB', ...
        'The Stage B file linked by C1 was not found:\n%s',C.source_stage_b_file);
end
SB = load(C.source_stage_b_file,'StageB');
if ~isfield(SB,'StageB')
    error('StageD1:InvalidStageB','The linked file does not contain StageB.');
end
StageB = SB.StageB;

if ~isfield(StageB,'source_stage_a_file') || ...
        ~isfile(StageB.source_stage_a_file)
    error('StageD1:MissingStageA','The Stage A file linked by Stage B was not found.');
end
SA = load(StageB.source_stage_a_file,'StageA');
if ~isfield(SA,'StageA') || ~SA.StageA.validation_passed
    error('StageD1:InvalidStageA','The linked Stage A result is missing or invalid.');
end
StageA = SA.StageA;

if C.pair_index < 1 || C.pair_index > numel(StageB.pairs)
    error('StageD1:PairIndexMismatch','C1 pair index is not present in Stage B.');
end
PairB = StageB.pairs(C.pair_index);
if ~strcmp(PairB.key,C.pair_key)
    error('StageD1:PairKeyMismatch','C1 and Stage B pair identities disagree.');
end

expected_conditions = {'A','B','AB','A_to_B','B_to_A'};
if ~isequal(C.condition_order,expected_conditions) || size(C.Data,3) ~= 5
    error('StageD1:ConditionMismatch', ...
        'C1 must contain A, B, AB, A_to_B, and B_to_A in the expected order.');
end

%% ======================= OUTPUT LOCATION ==============================
if isempty(export_folder)
    export_folder = fullfile(fileparts(stage_c1_file),'Luke_ModelPackage');
end
export_folder = char(string(export_folder));

experiment_folders = {StageA.config.single_folder,StageA.config.sim_folder, ...
    StageA.config.seq_folder};
for k = 1:numel(experiment_folders)
    if path_is_inside(export_folder,experiment_folders{k})
        error('StageD1:UnsafeOutputFolder', ...
            ['The export folder cannot be inside an experiment folder.\n' ...
             'Export: %s\nExperiment: %s'], ...
            export_folder,experiment_folders{k});
    end
end

[~,c1_stem] = fileparts(stage_c1_file);
dataset_tag = regexprep(c1_stem,'_StageC1_.*$','');
if isempty(output_filename)
    pair_tag = sprintf('%s_%s',safe_name(C.electrode_A),safe_name(C.electrode_B));
    output_filename = sprintf('%s_Pair_%s_ModelData.mat',dataset_tag,pair_tag);
else
    output_filename = char(string(output_filename));
    [~,~,ext] = fileparts(output_filename);
    if isempty(ext)
        output_filename = [output_filename '.mat'];
    elseif ~strcmpi(ext,'.mat')
        error('StageD1:OutputExtension','output_filename must have a .mat extension.');
    end
end
output_file = fullfile(export_folder,output_filename);
if isfile(output_file) && ~allow_overwrite
    error('StageD1:ExistingOutput', ...
        ['The output already exists and was not overwritten:\n%s\n' ...
         'Set allow_overwrite=true only when replacement is intended.'],output_file);
end

%% ====================== RESPONDING CHANNELS ===========================
CI = C.ChannelInfo;
required_channel_variables = {'DepthChannel','Shank','LocalSite','X_um','Y_um', ...
    'DistanceToA_um','DistanceToB_um','MinimumDistance_um', ...
    'DistanceToMidpoint_um','IsRespondingUnion'};
for k = 1:numel(required_channel_variables)
    if ~ismember(required_channel_variables{k},CI.Properties.VariableNames)
        error('StageD1:MissingChannelVariable', ...
            'C1 ChannelInfo is missing %s.',required_channel_variables{k});
    end
end

export_c1_rows = find(CI.IsRespondingUnion(:)).';
if isempty(export_c1_rows)
    error('StageD1:NoRespondingChannels', ...
        'C1 contains no responding-union recording channels.');
end

DepthChannel = CI.DepthChannel(export_c1_rows);
ChannelIndex = (1:numel(export_c1_rows)).';
Shank = CI.Shank(export_c1_rows);
LocalSite = CI.LocalSite(export_c1_rows);
X_um = CI.X_um(export_c1_rows);
Y_um = CI.Y_um(export_c1_rows);
DistanceToA_um = CI.DistanceToA_um(export_c1_rows);
DistanceToB_um = CI.DistanceToB_um(export_c1_rows);
MinimumStimDistance_um = CI.MinimumDistance_um(export_c1_rows);
DistanceToMidpoint_um = CI.DistanceToMidpoint_um(export_c1_rows);
IsStimA = DepthChannel == C.stimulation.A.depth_channel;
IsStimB = DepthChannel == C.stimulation.B.depth_channel;
IsResponding = true(numel(export_c1_rows),1);

ChannelTable = table(ChannelIndex,DepthChannel,Shank,LocalSite,X_um,Y_um, ...
    DistanceToA_um,DistanceToB_um,MinimumStimDistance_um, ...
    DistanceToMidpoint_um,IsStimA,IsStimB,IsResponding);

%% ==================== CONDITION DEFINITIONS ===========================
condition_labels = {'A alone','B alone','A+B simultaneous', ...
    'A then B, 5 ms','B then A, 5 ms'};
condition_sources = {'single','single','sim','seq','seq'};
electrode_orders = {{C.electrode_A},{C.electrode_B}, ...
    {C.electrode_A,C.electrode_B},{C.electrode_A,C.electrode_B}, ...
    {C.electrode_B,C.electrode_A}};
pulse_times = {[0],[0],[0 0],[0 StageB.config.seq_PTD_ms], ...
    [0 StageB.config.seq_PTD_ms]};

ConditionDefinitions = repmat(struct('index',[], 'code','', 'label','', ...
    'source_dataset','', 'electrode_order',{{}}, 'pulse_times_ms',[]),1,5);
for ic = 1:5
    ConditionDefinitions(ic).index = ic;
    ConditionDefinitions(ic).code = expected_conditions{ic};
    ConditionDefinitions(ic).label = condition_labels{ic};
    ConditionDefinitions(ic).source_dataset = condition_sources{ic};
    ConditionDefinitions(ic).electrode_order = electrode_orders{ic};
    ConditionDefinitions(ic).pulse_times_ms = pulse_times{ic};
end

%% ================= STIMULATION PARAMETER HEADERS ======================
E_single = StageA.datasets.single.experiment;
E_sim = StageA.datasets.sim.experiment;
E_seq = StageA.datasets.seq.experiment;
header_single = normalize_parameter_names(E_single.StimParams(1,:));
header_sim = normalize_parameter_names(E_sim.StimParams(1,:));
header_seq = normalize_parameter_names(E_seq.StimParams(1,:));
if ~isequal(header_single,header_sim) || ~isequal(header_single,header_seq)
    error('StageD1:ParameterHeaderMismatch', ...
        'StimParams column headers differ between Single, Sim, and Seq datasets.');
end
parameter_names = header_single;

%% =================== BALANCED SUBSET INDICES ==========================
nAmp = numel(C.amplitudes_uA);
nCondition = 5;
all_counts = zeros(nAmp,nCondition);
balanced_counts = zeros(nAmp,nCondition);
balanced_indices = cell(nAmp,nCondition);
random_stream = RandStream('mt19937ar','Seed',balancing_random_seed);

for ia = 1:nAmp
    for ic = 1:nCondition
        all_counts(ia,ic) = numel(C.Data(1,ia,ic).source_trial_ids);
    end
    nBalanced = min(all_counts(ia,:));
    if nBalanced < 1
        error('StageD1:EmptyCondition', ...
            'At least one condition has no trials at %g uA.',C.amplitudes_uA(ia));
    end
    for ic = 1:nCondition
        chosen = randperm(random_stream,all_counts(ia,ic),nBalanced);
        balanced_indices{ia,ic} = sort(chosen(:));
        balanced_counts(ia,ic) = nBalanced;
    end
end

%% ====================== BUILD EXPORTED DATA ===========================
nChannels = numel(export_c1_rows);
Conditions = repmat(struct('index',[], 'code','', 'label','', ...
    'source_dataset','', 'electrode_order',{{}}, 'pulse_times_ms',[], ...
    'amplitude',[]),1,nCondition);
all_parameter_blocks_constant = true;

for ic = 1:nCondition
    Conditions(ic).index = ic;
    Conditions(ic).code = expected_conditions{ic};
    Conditions(ic).label = condition_labels{ic};
    Conditions(ic).source_dataset = condition_sources{ic};
    Conditions(ic).electrode_order = electrode_orders{ic};
    Conditions(ic).pulse_times_ms = pulse_times{ic};

    AmpBlocks = repmat(empty_amplitude_block(),1,nAmp);
    for ia = 1:nAmp
        amp = C.amplitudes_uA(ia);
        reference_data = C.Data(export_c1_rows(1),ia,ic);
        trial_ids = double(reference_data.source_trial_ids(:));
        nTrials = numel(trial_ids);
        if nTrials ~= all_counts(ia,ic)
            error('StageD1:TrialCountMismatch','Unexpected C1 trial count mismatch.');
        end

        spike_times_ms = cell(nTrials,nChannels);
        baseline_count = nan(nTrials,nChannels);
        baseline_rate_hz = nan(nTrials,nChannels);
        post_count = nan(nTrials,nChannels);
        post_rate_hz = nan(nTrials,nChannels);
        expected_baseline_post_count = nan(nTrials,nChannels);
        evoked_spike_count = nan(nTrials,nChannels);
        has_response_spike = false(nTrials,nChannels);
        first_spike_latency_ms = nan(nTrials,nChannels);

        for jc = 1:nChannels
            D = C.Data(export_c1_rows(jc),ia,ic);
            if ~isequal(double(D.source_trial_ids(:)),trial_ids)
                error('StageD1:ChannelTrialMismatch', ...
                    'Source trial IDs differ across channels for %g uA, %s.', ...
                    amp,expected_conditions{ic});
            end
            if numel(D.spike_times_ms) ~= nTrials
                error('StageD1:SpikeTrialMismatch', ...
                    'Spike cell count differs from source trial count.');
            end
            for it = 1:nTrials
                spike_times_ms{it,jc} = double(D.spike_times_ms{it}(:));
            end

            M = D.trial_metrics;
            required_metrics = {'baseline_count','baseline_rate_hz','post_count', ...
                'post_rate_hz','has_response','first_spike_latency_ms','evoked_count'};
            require_fields(M,required_metrics,'C1 trial_metrics');
            baseline_count(:,jc) = double(M.baseline_count(:));
            baseline_rate_hz(:,jc) = double(M.baseline_rate_hz(:));
            post_count(:,jc) = double(M.post_count(:));
            post_rate_hz(:,jc) = double(M.post_rate_hz(:));
            has_response_spike(:,jc) = logical(M.has_response(:));
            first_spike_latency_ms(:,jc) = double(M.first_spike_latency_ms(:));
            evoked_spike_count(:,jc) = double(M.evoked_count(:));
        end

        post_duration_s = diff(C.windows.count_ms)/1000;
        expected_baseline_post_count = baseline_rate_hz*post_duration_s;
        if any(abs((post_count-expected_baseline_post_count)- ...
                evoked_spike_count) > 1e-9,'all')
            error('StageD1:EvokedCountMismatch', ...
                'Stored C1 evoked counts do not match the documented calculation.');
        end

        ChannelIndexSummary = ChannelIndex;
        DepthChannelSummary = DepthChannel;
        NTrials = repmat(nTrials,nChannels,1);
        MeanPostCount = mean(post_count,1).';
        SDPostCount = std(post_count,0,1).';
        SEMPostCount = SDPostCount/sqrt(nTrials);
        MeanPostRate_Hz = mean(post_rate_hz,1).';
        ResponseProbability = mean(has_response_spike,1).';
        MeanBaselineRate_Hz = mean(baseline_rate_hz,1).';
        MeanEvokedCount = mean(evoked_spike_count,1).';
        MeanFirstSpikeLatency_ms = mean(first_spike_latency_ms,1,'omitnan').';
        Summary = table(ChannelIndexSummary,DepthChannelSummary,NTrials, ...
            MeanPostCount,SDPostCount,SEMPostCount,MeanPostRate_Hz, ...
            ResponseProbability,MeanBaselineRate_Hz,MeanEvokedCount, ...
            MeanFirstSpikeLatency_ms);

        amp_idx_B = find(abs([PairB.amplitudes.amplitude_uA]-amp) <= 1e-6,1);
        if isempty(amp_idx_B)
            error('StageD1:StageBAmplitudeMissing', ...
                'Amplitude %g uA was not found in Stage B.',amp);
        end
        AmpB = PairB.amplitudes(amp_idx_B);
        switch expected_conditions{ic}
            case {'A','B'}
                response_label_available = false;
                is_responsive = nan(nChannels,1);
            case 'AB'
                response_label_available = true;
                is_responsive = ismember(DepthChannel, ...
                    AmpB.responding_depth_channels.AB);
            case 'A_to_B'
                response_label_available = true;
                is_responsive = ismember(DepthChannel, ...
                    AmpB.responding_depth_channels.A_to_B);
            case 'B_to_A'
                response_label_available = true;
                is_responsive = ismember(DepthChannel, ...
                    AmpB.responding_depth_channels.B_to_A);
        end

        switch condition_sources{ic}
            case 'single', E = E_single;
            case 'sim',    E = E_sim;
            case 'seq',    E = E_seq;
        end
        StimEvents = collect_stimulation_events(E,trial_ids, ...
            electrode_orders{ic},pulse_times{ic});
        event_constant = [StimEvents.constant_across_trials];
        all_parameter_blocks_constant = all_parameter_blocks_constant && ...
            all(event_constant);

        B = empty_amplitude_block();
        B.amplitude_uA = amp;
        B.n_trials_all = nTrials;
        B.n_trials_balanced = numel(balanced_indices{ia,ic});
        B.source_trial_ids = trial_ids;
        B.all_trial_indices = (1:nTrials).';
        B.balanced_trial_indices = balanced_indices{ia,ic};
        B.balanced_source_trial_ids = trial_ids(B.balanced_trial_indices);
        B.spike_times_ms = spike_times_ms;
        B.trial_metrics = struct( ...
            'baseline_count',baseline_count, ...
            'baseline_rate_hz',baseline_rate_hz, ...
            'post_count',post_count, ...
            'post_rate_hz',post_rate_hz, ...
            'expected_baseline_post_count',expected_baseline_post_count, ...
            'evoked_spike_count',evoked_spike_count, ...
            'has_response_spike',has_response_spike, ...
            'first_spike_latency_ms',first_spike_latency_ms);
        B.summary_all_trials = Summary;
        B.response_label_available = response_label_available;
        B.is_responsive = is_responsive(:);
        B.stimulation_events = StimEvents;
        AmpBlocks(ia) = B;
    end
    Conditions(ic).amplitude = AmpBlocks;
end

trial_names = {'A','B','AB','A_to_B','B_to_A'};
AllTrialCounts = array2table(all_counts,'VariableNames',trial_names, ...
    'RowNames',cellstr(compose('%g_uA',C.amplitudes_uA(:))));
BalancedTrialCounts = array2table(balanced_counts,'VariableNames',trial_names, ...
    'RowNames',cellstr(compose('%g_uA',C.amplitudes_uA(:))));

%% ======================== TOP-LEVEL PACKAGE ===========================
ModelData = struct();
ModelData.format_name = 'SinglePairStimulationModelData';
ModelData.format_version = '1.0';
ModelData.created_on = datetime('now');
ModelData.dataset_id = dataset_tag;
ModelData.pair_key = C.pair_key;

ModelData.metadata = struct( ...
    'sampling_rate_hz',C.FS, ...
    'time_unit','ms', ...
    'spike_time_reference','first stimulation pulse at 0 ms', ...
    'stored_spike_window_ms',C.windows.raster_ms, ...
    'default_baseline_window_ms',C.windows.baseline_ms, ...
    'default_count_window_ms',C.windows.count_ms, ...
    'window_interval_convention','[start,end), except stored raster includes end', ...
    'sequential_PTD_ms',StageB.config.seq_PTD_ms);

ModelData.qc = struct('single_trial_review_complete', ...
    logical(C.single_qc.review_complete));

ModelData.trial_selection = struct( ...
    'available_modes',{{'all','balanced'}}, ...
    'default_mode','all', ...
    'balancing_method','fixed-seed random subset within amplitude and condition', ...
    'balancing_random_seed',balancing_random_seed, ...
    'all_trial_counts',AllTrialCounts, ...
    'balanced_trial_counts',BalancedTrialCounts);

ModelData.stimulation = struct( ...
    'electrode_A',C.electrode_A, ...
    'electrode_B',C.electrode_B, ...
    'pair_distance_um',C.stimulation.pair_distance_um, ...
    'same_shank',C.stimulation.same_shank, ...
    'parameter_names',{parameter_names}, ...
    'all_parameter_blocks_constant_across_trials', ...
        all_parameter_blocks_constant);

ModelData.amplitudes_uA = double(C.amplitudes_uA(:).');
ModelData.condition_definitions = ConditionDefinitions;
ModelData.channels = ChannelTable;
ModelData.conditions = Conditions;
ModelData.source = struct( ...
    'stage_c1_filename',file_name_only(stage_c1_file), ...
    'stage_b_filename',file_name_only(C.source_stage_b_file), ...
    'stage_a_filename',file_name_only(StageB.source_stage_a_file), ...
    'spike_variable','sp_corr');

ModelData.experiment_files_modified = false;

%% ====================== VALIDATE FINAL PACKAGE ========================
validate_exported_package(ModelData);

%% ======================= COMMAND REPORT ===============================
fprintf('\n============================================================\n');
fprintf('STAGE D1: EXPORT MODELLING DATA\n');
fprintf('============================================================\n');
fprintf('C1 source: %s\n',stage_c1_file);
fprintf('Pair: %s\n',C.pair_key);
fprintf('Amplitudes: %s uA\n',num2str(C.amplitudes_uA));
fprintf('Responding channels exported: %d\n',height(ChannelTable));
fprintf('Depth channels: %s\n',num2str(ChannelTable.DepthChannel.'));
fprintf('Stored spike window: [%g %g] ms\n',C.windows.raster_ms);
fprintf('Default baseline window: [%g %g) ms\n',C.windows.baseline_ms);
fprintf('Default count window: [%g %g) ms\n',C.windows.count_ms);
fprintf('Single-trial review complete: %s\n', ...
    logical_text(C.single_qc.review_complete));
if ~C.single_qc.review_complete
    fprintf(['NOTICE: Single trials are currently retained without a ' ...
        'completed manual review.\n']);
end

fprintf('\n[ALL INCLUDED TRIALS]\n%s\n',evalc('disp(AllTrialCounts)'));
fprintf('[OPTIONAL BALANCED SUBSETS]\n%s\n', ...
    evalc('disp(BalancedTrialCounts)'));
fprintf('Balancing random seed: %d\n',balancing_random_seed);
fprintf('Default trial mode: all\n');

fprintf('\n[STIMULATION]\n');
fprintf('Electrode A: %s\n',C.electrode_A);
fprintf('Electrode B: %s\n',C.electrode_B);
fprintf('Pair distance: %.1f um\n',C.stimulation.pair_distance_um);
fprintf('Sequential PTD: %.3f ms\n',StageB.config.seq_PTD_ms);
fprintf('Recorded parameter columns: %d\n',numel(parameter_names));
fprintf('Recorded parameter names:\n  %s\n',strjoin(parameter_names,' | '));
fprintf('All parameter blocks constant across included trials: %s\n', ...
    logical_text(all_parameter_blocks_constant));
if ~all_parameter_blocks_constant
    fprintf(['NOTICE: At least one recorded stimulation parameter differs ' ...
        'across included trials. All unique rows are retained in the package.\n']);
end

fprintf('\n[OUTPUT]\n');
fprintf('MAT file: %s\n',output_file);
fprintf('Original experiment files modified: NO\n');

%% ============================= SAVE ==================================
if ~isfolder(export_folder)
    mkdir(export_folder);
end
save(output_file,'ModelData','-v7.3');

fprintf('\n============================================================\n');
fprintf('STAGE D1 COMPLETE\n');
fprintf('Exported file: %s\n',output_file);
fprintf('Original experiment files modified: NO\n');
fprintf('============================================================\n\n');

%% =========================== FUNCTIONS ================================
function B = empty_amplitude_block()
B = struct('amplitude_uA',NaN, 'n_trials_all',0, ...
    'n_trials_balanced',0, 'source_trial_ids',[], ...
    'all_trial_indices',[], 'balanced_trial_indices',[], ...
    'balanced_source_trial_ids',[], 'spike_times_ms',{{}}, ...
    'trial_metrics',struct(), 'summary_all_trials',table(), ...
    'response_label_available',false, 'is_responsive',[], ...
    'stimulation_events',struct([]));
end

function Events = collect_stimulation_events(E,trial_ids,electrode_order,pulse_times)
nEvents = numel(electrode_order);
nParams = size(E.StimParams,2);
Events = repmat(struct('pulse_number',[], 'electrode','', ...
    'pulse_time_ms',NaN, 'representative_parameter_row',{{}}, ...
    'unique_parameter_rows',{{}}, 'n_unique_parameter_rows',0, ...
    'parameter_row_index_per_trial',[], ...
    'constant_across_trials',false),1,nEvents);

for ie = 1:nEvents
    rows_for_event = cell(numel(trial_ids),nParams);
    for it = 1:numel(trial_ids)
        tr = trial_ids(it);
        if tr < 1 || tr > E.n_trials
            error('StageD1:StimTrialOutOfRange', ...
                'A source trial ID is outside the stimulation parameter data.');
        end
        names = E.stim_names_per_trial(tr,:);
        match = find(strcmp(names,electrode_order{ie}),1);
        if isempty(match)
            error('StageD1:StimElectrodeMissing', ...
                'Trial %d does not contain expected electrode %s.', ...
                tr,electrode_order{ie});
        end
        source_rows = E.trial_stimparam_rows{tr};
        rows_for_event(it,:) = E.StimParams(source_rows(match),:);
    end
    [unique_rows,row_index_per_trial] = unique_cell_rows(rows_for_event);
    Events(ie).pulse_number = ie;
    Events(ie).electrode = electrode_order{ie};
    Events(ie).pulse_time_ms = pulse_times(ie);
    Events(ie).representative_parameter_row = rows_for_event(1,:);
    Events(ie).unique_parameter_rows = unique_rows;
    Events(ie).n_unique_parameter_rows = size(unique_rows,1);
    Events(ie).parameter_row_index_per_trial = row_index_per_trial;
    Events(ie).constant_across_trials = size(unique_rows,1) == 1;
end
end

function [rows_out,row_index] = unique_cell_rows(rows_in)
if isempty(rows_in)
    rows_out = rows_in;
    row_index = zeros(0,1);
    return;
end
keep = false(size(rows_in,1),1);
representatives = zeros(0,1);
row_index = zeros(size(rows_in,1),1);
for r = 1:size(rows_in,1)
    is_duplicate = false;
    for j = 1:numel(representatives)
        if isequaln(rows_in(r,:),rows_in(representatives(j),:))
            is_duplicate = true;
            row_index(r) = j;
            break;
        end
    end
    if ~is_duplicate
        keep(r) = true;
        representatives(end+1,1) = r; %#ok<AGROW>
        row_index(r) = numel(representatives);
    end
end
rows_out = rows_in(keep,:);
end

function names = normalize_parameter_names(header)
names = cell(size(header));
for k = 1:numel(header)
    value = header{k};
    if isempty(value)
        names{k} = sprintf('Parameter_%d',k);
    else
        names{k} = strtrim(char(string(value)));
        if isempty(names{k})
            names{k} = sprintf('Parameter_%d',k);
        end
    end
end
names = names(:).';
end

function validate_exported_package(M)
if height(M.channels) < 1
    error('StageD1:PackageNoChannels','Exported package has no channels.');
end
if numel(M.conditions) ~= 5
    error('StageD1:PackageConditionCount','Exported package must have five conditions.');
end
nChannels = height(M.channels);
nAmp = numel(M.amplitudes_uA);
for ic = 1:5
    if numel(M.conditions(ic).amplitude) ~= nAmp
        error('StageD1:PackageAmplitudeCount','Exported amplitude count mismatch.');
    end
    for ia = 1:nAmp
        B = M.conditions(ic).amplitude(ia);
        if size(B.spike_times_ms,1) ~= B.n_trials_all || ...
                size(B.spike_times_ms,2) ~= nChannels
            error('StageD1:PackageSpikeShape','Exported spike cell shape is invalid.');
        end
        metric_names = fieldnames(B.trial_metrics);
        for k = 1:numel(metric_names)
            values = B.trial_metrics.(metric_names{k});
            if ~isequal(size(values),[B.n_trials_all nChannels])
                error('StageD1:PackageMetricShape', ...
                    'Exported metric %s has an invalid shape.',metric_names{k});
            end
        end
        if any(B.balanced_trial_indices < 1) || ...
                any(B.balanced_trial_indices > B.n_trials_all) || ...
                numel(unique(B.balanced_trial_indices)) ~= B.n_trials_balanced
            error('StageD1:PackageBalancedIndices', ...
                'Balanced trial indices are invalid.');
        end
    end
end
end

function require_fields(S,names,label)
for k = 1:numel(names)
    if ~isfield(S,names{k})
        error('StageD1:MissingField','%s.%s is missing.',label,names{k});
    end
end
end

function name = file_name_only(path_value)
[~,stem,ext] = fileparts(path_value);
name = [stem ext];
end

function value = safe_name(value)
value = regexprep(char(string(value)),'[^A-Za-z0-9]+','');
if isempty(value)
    value = 'Unknown';
end
end

function value = logical_text(tf)
if tf, value = 'true'; else, value = 'false'; end
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
