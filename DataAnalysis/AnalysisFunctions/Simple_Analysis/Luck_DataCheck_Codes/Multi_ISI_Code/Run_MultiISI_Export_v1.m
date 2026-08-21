%% MULTI-ISI MODELLING DATA EXPORT
%
% Run this single script to create one clean ModelData MAT file for every
% stimulation pair found in a multi-ISI simultaneous/sequential dataset.
%
% IMPORTANT DATA RULES
%   - Uses sp_corr from *sp_xia_SSD.mat only.
%   - Aligns every spike time to the first stimulation pulse (0 ms).
%   - Removes bad trials before saving. Bad-trial lists are not exported.
%   - Supports global or recording-channel-specific bad trials.
%   - Removes bad recording channels for the stimulation pair.
%   - Exports only the manually corrected responding-channel union.
%   - Saves every observed amplitude, PTD, and sequential order.
%   - Does not balance or downsample trials.
%   - Does not modify any original experiment file.
%
% CLEAN SPIKE LOCATION
%   ModelData.conditions(condition_index) ...
%            .amplitude(amplitude_index) ...
%            .channel(channel_index).spike_times_ms{trial_index}
%
% A different number of clean trials is allowed for each recording channel.

clear;
clc;

%% ========================= USER SETTINGS ==============================

% Single-electrode stimulation dataset containing A-alone and B-alone.
single_folder = '/Volumes/MACData/Data/Data_Xia/DX020/Xia_ISI_Single1';

% Combined multi-ISI dataset containing simultaneous (PTD=0) and
% sequential (PTD>0) stimulation.
multiisi_folder = '/Volumes/MACData/Data/Data_Xia/DX020/Xia_ISI_SimSeq1';

Electrode_Type = 2;       % 0=rigid, 1=single-shank flex, 2=four-shank flex
FS = 30000;               % sampling rate in Hz

% Leave empty to derive Data_Xia/AnalysisFunctions from single_folder.
analysis_functions_folder = '';

% Leave empty to create multiisi_model_exports beside this script.
output_folder = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/Multi_ISI_Output';

% Empty means export every observed positive PTD. Simultaneous PTD=0 is
% always included when it exists.
Positive_PTDs_To_Export = [];

% Empty means export every amplitude observed within each condition.
Amplitudes_To_Export = [];

% Spike and analysis windows in ms relative to the first pulse.
stored_spike_window_ms = [-50 80];
baseline_window_ms     = [-50 -5];
whole_response_ms      = [2 45];
early_window_ms        = [2 20];

% For a sequential trial, this window is shifted by its PTD. For example,
% at PTD=10 ms, [2 20] becomes [12 30] ms relative to pulse 1.
second_window_relative_ms = [2 20];

% Optional additional exclusions for the Single dataset. These trial IDs
% are removed from every recording channel in addition to any saved file.
single_bad_trials_manual = [];

% Protect an existing export unless replacement is explicitly requested.
allow_overwrite = false;

% Optional geometry checks. Use NaN to disable the expected-spacing check.
require_same_stim_shank = true;
expected_stim_spacing_um = 200;
geometry_tolerance_um = 0.01;

%% ====================== INITIAL VALIDATION ============================
validate_folder(single_folder,'single_folder');
validate_folder(multiisi_folder,'multiisi_folder');
if ~ismember(Electrode_Type,[0 1 2])
    error('MultiISI:ElectrodeType','Electrode_Type must be 0, 1, or 2.');
end
if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('MultiISI:SamplingRate','FS must be one positive number.');
end
validate_window(stored_spike_window_ms,'stored_spike_window_ms');
validate_window(baseline_window_ms,'baseline_window_ms');
validate_window(whole_response_ms,'whole_response_ms');
validate_window(early_window_ms,'early_window_ms');
validate_window(second_window_relative_ms,'second_window_relative_ms');
if baseline_window_ms(1) < stored_spike_window_ms(1) || ...
        baseline_window_ms(2) > stored_spike_window_ms(2) || ...
        whole_response_ms(1) < stored_spike_window_ms(1) || ...
        whole_response_ms(2) > stored_spike_window_ms(2) || ...
        early_window_ms(1) < stored_spike_window_ms(1) || ...
        early_window_ms(2) > stored_spike_window_ms(2)
    error('MultiISI:WindowOutsideStoredRange', ...
        'Baseline, whole-response and early windows must lie inside the stored range.');
end

single_bad_trials_manual = unique(double(single_bad_trials_manual(:).'));
Positive_PTDs_To_Export = unique(double(Positive_PTDs_To_Export(:).'));
Amplitudes_To_Export = unique(double(Amplitudes_To_Export(:).'));
if any(Positive_PTDs_To_Export <= 0)
    error('MultiISI:PTDSelection', ...
        'Positive_PTDs_To_Export must contain positive PTDs only.');
end

if isempty(analysis_functions_folder)
    data_xia_folder = fileparts(fileparts(single_folder));
    analysis_functions_folder = fullfile(data_xia_folder,'AnalysisFunctions');
end
validate_folder(analysis_functions_folder,'analysis_functions_folder');
addpath(genpath(analysis_functions_folder),'-end');

script_folder = fileparts(mfilename('fullpath'));
if isempty(output_folder)
    output_folder = fullfile(script_folder,'multiisi_model_exports');
end
output_folder = char(string(output_folder));
if path_is_inside(output_folder,single_folder) || ...
        path_is_inside(output_folder,multiisi_folder)
    error('MultiISI:UnsafeOutput', ...
        'output_folder cannot be inside an experiment folder.');
end
if ~isfolder(output_folder)
    mkdir(output_folder);
end

%% ============================ LOAD ====================================
fprintf('\n============================================================\n');
fprintf('MULTI-ISI MODELLING DATA EXPORT\n');
fprintf('============================================================\n');
fprintf('Single dataset:   %s\n',single_folder);
fprintf('Multi-ISI dataset: %s\n',multiisi_folder);
fprintf('Output folder:    %s\n',output_folder);
fprintf('Original experiment files will not be modified.\n');

Single = load_dataset(single_folder,Electrode_Type,FS,false);
Multi = load_dataset(multiisi_folder,Electrode_Type,FS,true);

if Single.experiment.simultaneous_stim ~= 1
    error('MultiISI:SingleStimCount', ...
        'The Single dataset must contain one stimulation event per trial.');
end
if Multi.experiment.simultaneous_stim ~= 2
    error('MultiISI:MultiStimCount', ...
        'The Multi-ISI dataset must contain two stimulation events per trial.');
end
if Single.spike.n_channels ~= Multi.spike.n_channels
    error('MultiISI:SpikeChannelCount', ...
        'Single and Multi-ISI sp_corr channel counts differ.');
end
if numel(Single.depth_to_hardware) ~= numel(Multi.depth_to_hardware)
    error('MultiISI:DepthChannelCount', ...
        'Single and Multi-ISI depth-channel counts differ.');
end
if ~Multi.responding.file_found
    error('MultiISI:MissingResponding', ...
        ['A manually corrected MultiISI RespondingChannels file was not ' ...
         'found in the Multi-ISI folder.']);
end
single_headers = normalize_parameter_names(Single.experiment.StimParams(1,:));
multi_headers = normalize_parameter_names(Multi.experiment.StimParams(1,:));
if ~isequal(single_headers,multi_headers)
    error('MultiISI:ParameterHeaders', ...
        'Single and Multi-ISI StimParams column headers differ.');
end

available_positive_ptds = unique(Multi.experiment.trial_ptd_ms( ...
    Multi.experiment.trial_ptd_ms > 0));
if isempty(Positive_PTDs_To_Export)
    selected_positive_ptds = available_positive_ptds(:).';
else
    selected_positive_ptds = available_positive_ptds(ismembertol( ...
        available_positive_ptds,Positive_PTDs_To_Export,1e-6)).';
end
if ~isempty(selected_positive_ptds) && ...
        max(selected_positive_ptds)+second_window_relative_ms(2) > ...
        stored_spike_window_ms(2)
    error('MultiISI:SecondWindowOutsideStoredRange', ...
        ['The largest selected PTD plus the second-window endpoint exceeds ' ...
         'the stored spike window.']);
end

pair_defs = discover_pairs(Multi.experiment,selected_positive_ptds);
if isempty(pair_defs)
    error('MultiISI:NoPairs','No two-electrode stimulation pairs were found.');
end

fprintf('\nObserved positive PTDs selected for export: %s ms\n', ...
    vector_text(selected_positive_ptds));
fprintf('Candidate stimulation pairs: %d\n',numel(pair_defs));

%% ======================= EXPORT EACH PAIR =============================
exported_files = strings(0,1);
skipped_pairs = strings(0,1);
nDepth = numel(Multi.depth_to_hardware);
amp_tol = 1e-6;
ptd_tol = 1e-6;

for ipair = 1:numel(pair_defs)
    A = pair_defs(ipair).A;
    B = pair_defs(ipair).B;
    pair_key = pair_defs(ipair).key;

    fprintf('\n------------------------------------------------------------\n');
    fprintf('PAIR %d/%d: %s\n',ipair,numel(pair_defs),pair_key);

    descriptors = build_condition_descriptors(Single.experiment, ...
        Multi.experiment,A,B,selected_positive_ptds,ptd_tol);
    descriptor_keep = false(size(descriptors));
    for id = 1:numel(descriptors)
        if strcmp(descriptors(id).source_dataset,'single')
            Echeck = Single.experiment;
        else
            Echeck = Multi.experiment;
        end
        descriptor_keep(id) = ~isempty(selected_amplitudes( ...
            Echeck.trial_amplitude(descriptors(id).trial_ids_all), ...
            Amplitudes_To_Export,amp_tol));
    end
    descriptors = descriptors(descriptor_keep);
    if isempty(descriptors)
        skipped_pairs(end+1,1) = string(pair_key); %#ok<SAGROW>
        fprintf('No recorded conditions selected; pair skipped.\n');
        continue;
    end

    % Remove channels marked bad for any stimulation set belonging to this
    % pair. This preserves one simple recording-channel list for the package.
    paired_trials = pair_trials_for_descriptors(descriptors);
    paired_sets = unique(Multi.experiment.set_index(paired_trials)).';
    bad_depth_union = bad_channels_for_sets(Multi.bad_channels,paired_sets);
    bad_depth_union = bad_depth_union( ...
        bad_depth_union >= 1 & bad_depth_union <= nDepth);

    % Build the manually corrected response union across every observed
    % paired condition and amplitude that will be exported.
    response_union = false(1,nDepth);
    for id = 1:numel(descriptors)
        if strcmp(descriptors(id).source_dataset,'single')
            continue;
        end
        trial_ids = descriptors(id).trial_ids_all;
        amps = selected_amplitudes(Multi.experiment.trial_amplitude(trial_ids), ...
            Amplitudes_To_Export,amp_tol);
        for amp = amps
            amp_trials = trial_ids(abs( ...
                Multi.experiment.trial_amplitude(trial_ids)-amp) <= amp_tol);
            set_indices = unique(Multi.experiment.set_index(amp_trials)).';
            response_union = response_union | responding_mask(Multi, ...
                set_indices,amp,descriptors(id).PTD_ms,amp_tol,ptd_tol,nDepth);
        end
    end
    response_union(bad_depth_union) = false;
    selected_depth_channels = find(response_union);
    if isempty(selected_depth_channels)
        skipped_pairs(end+1,1) = string(pair_key); %#ok<SAGROW>
        fprintf('No manually corrected responding channels; pair skipped.\n');
        continue;
    end

    stimA = map_stimulation_contact(A,Electrode_Type,nDepth, ...
        Multi.depth_to_hardware);
    stimB = map_stimulation_contact(B,Electrode_Type,nDepth, ...
        Multi.depth_to_hardware);
    pair_distance_um = euclidean_distance(stimA.x_um,stimA.y_um, ...
        stimB.x_um,stimB.y_um);
    same_shank = stimA.shank == stimB.shank;
    if require_same_stim_shank && ~same_shank
        error('MultiISI:StimShankMismatch', ...
            'Pair %s is not mapped to one shank.',pair_key);
    end
    if isfinite(expected_stim_spacing_um) && ...
            abs(pair_distance_um-expected_stim_spacing_um) > geometry_tolerance_um
        error('MultiISI:StimSpacingMismatch', ...
            'Pair %s spacing is %.3f um; expected %.3f um.', ...
            pair_key,pair_distance_um,expected_stim_spacing_um);
    end

    ChannelTable = build_channel_table(selected_depth_channels,Single,Multi, ...
        Electrode_Type,stimA,stimB);

    fprintf('Responding channels exported: %d\n',height(ChannelTable));
    fprintf('Depth channels: %s\n',num2str(ChannelTable.DepthChannel.'));
    fprintf('Loading selected sp_corr channels...\n');
    spike_single = load_selected_spike_times(Single.files.spike.path, ...
        ChannelTable.SpCellIndex_Single.');
    spike_multi = load_selected_spike_times(Multi.files.spike.path, ...
        ChannelTable.SpCellIndex_MultiISI.');

    Conditions = repmat(empty_condition(),1,numel(descriptors));
    condition_rows = cell(numel(descriptors),6);
    all_parameter_blocks_constant = true;

    for id = 1:numel(descriptors)
        desc = descriptors(id);
        if strcmp(desc.source_dataset,'single')
            Dataset = Single;
            spike_source = spike_single;
        else
            Dataset = Multi;
            spike_source = spike_multi;
        end

        raw_amps = Dataset.experiment.trial_amplitude(desc.trial_ids_all);
        amps = selected_amplitudes(raw_amps,Amplitudes_To_Export,amp_tol);
        AmpBlocks = repmat(empty_amplitude(),1,numel(amps));

        for ia = 1:numel(amps)
            amp = amps(ia);
            raw_trials = desc.trial_ids_all(abs( ...
                Dataset.experiment.trial_amplitude(desc.trial_ids_all)-amp) <= amp_tol);
            ChannelBlocks = repmat(empty_channel_data(),1,height(ChannelTable));

            % Stimulation parameters do not depend on recording channel.
            StimEvents = collect_stimulation_events(Dataset.experiment, ...
                raw_trials,desc.electrode_order,desc.pulse_times_ms);
            if ~isempty(StimEvents)
                all_parameter_blocks_constant = all_parameter_blocks_constant && ...
                    all([StimEvents.constant_across_trials]);
            end

            for jc = 1:height(ChannelTable)
                depth_ch = ChannelTable.DepthChannel(jc);
                saved_bad = bad_trials_for_channel(Dataset.bad_trials,depth_ch);
                if strcmp(desc.source_dataset,'single')
                    saved_bad = union(saved_bad,single_bad_trials_manual);
                end
                clean_trials = setdiff(raw_trials(:).',saved_bad,'stable');

                ChannelBlocks(jc) = extract_clean_channel( ...
                    spike_source{jc},Dataset.trig,clean_trials,FS, ...
                    stored_spike_window_ms,baseline_window_ms, ...
                    whole_response_ms,early_window_ms, ...
                    second_window_relative_ms,desc.PTD_ms, ...
                    strcmp(desc.stimulation_type,'sequential'), ...
                    ChannelTable.ChannelIndex(jc),depth_ch);
            end

            Bamp = empty_amplitude();
            Bamp.amplitude_uA = amp;
            Bamp.channel = ChannelBlocks;
            Bamp.stimulation_events = StimEvents;
            Bamp.n_trials_per_channel = [ChannelBlocks.n_trials].';
            Bamp.minimum_clean_trials = min(Bamp.n_trials_per_channel);
            Bamp.maximum_clean_trials = max(Bamp.n_trials_per_channel);
            AmpBlocks(ia) = Bamp;
        end

        C = empty_condition();
        C.index = id;
        C.condition_id = desc.condition_id;
        C.code = desc.code;
        C.label = desc.label;
        C.stimulation_type = desc.stimulation_type;
        C.source_dataset = desc.source_dataset;
        C.electrode_order = desc.electrode_order;
        C.PTD_ms = desc.PTD_ms;
        C.pulse_times_ms = desc.pulse_times_ms;
        C.amplitude = AmpBlocks;
        Conditions(id) = C;

        condition_rows(id,:) = {id,string(desc.condition_id), ...
            string(desc.stimulation_type),string(strjoin(desc.electrode_order,' -> ')), ...
            desc.PTD_ms,string(vector_text(amps))};

        fprintf('  %-24s | PTD %g ms | amplitudes %s uA\n', ...
            desc.condition_id,desc.PTD_ms,vector_text(amps));
    end

    ConditionTable = cell2table(condition_rows,'VariableNames', ...
        {'ConditionIndex','ConditionID','StimulationType', ...
         'ElectrodeOrder','PTD_ms','Amplitudes_uA'});

    parameter_names = normalize_parameter_names( ...
        Multi.experiment.StimParams(1,:));
    dataset_id = derive_dataset_id(single_folder,multiisi_folder);

    ModelData = struct();
    ModelData.format_name = 'MultiISIStimulationModelData';
    ModelData.format_version = '1.0';
    ModelData.created_on = datetime('now');
    ModelData.dataset_id = dataset_id;
    ModelData.pair_key = pair_key;
    ModelData.metadata = struct( ...
        'sampling_rate_hz',FS, ...
        'time_unit','ms', ...
        'spike_time_reference','first stimulation pulse at 0 ms', ...
        'stored_spike_window_ms',stored_spike_window_ms, ...
        'baseline_window_ms',baseline_window_ms, ...
        'whole_response_window_ms',whole_response_ms, ...
        'early_window_ms',early_window_ms, ...
        'second_window_relative_to_second_pulse_ms', ...
            second_window_relative_ms, ...
        'window_interval_convention','[start,end), except stored spike window includes end', ...
        'trials_balanced',false);
    ModelData.stimulation = struct( ...
        'electrode_A',A, ...
        'electrode_B',B, ...
        'pair_distance_um',pair_distance_um, ...
        'same_shank',same_shank, ...
        'parameter_names',{parameter_names}, ...
        'all_parameter_blocks_constant_across_trials', ...
            all_parameter_blocks_constant);
    ModelData.channels = ChannelTable;
    ModelData.ConditionTable = ConditionTable;
    ModelData.conditions = Conditions;
    ModelData.source = struct( ...
        'single_dataset_name',file_name_only(single_folder), ...
        'multiisi_dataset_name',file_name_only(multiisi_folder), ...
        'spike_variable','sp_corr');
    ModelData.experiment_files_modified = false;

    validate_model_data(ModelData);

    pair_tag = sprintf('%s_%s',safe_name(A),safe_name(B));
    output_name = sprintf('%s_Pair_%s_MultiISI_ModelData.mat', ...
        dataset_id,pair_tag);
    output_file = fullfile(output_folder,output_name);
    if isfile(output_file) && ~allow_overwrite
        error('MultiISI:ExistingOutput', ...
            ['Output already exists and was not overwritten:\n%s\n' ...
             'Set allow_overwrite=true only when replacement is intended.'], ...
            output_file);
    end
    save(output_file,'ModelData','-v7.3');
    exported_files(end+1,1) = string(output_file); %#ok<SAGROW>
    fprintf('Saved: %s\n',output_file);
end

%% ============================ SUMMARY =================================
fprintf('\n============================================================\n');
fprintf('MULTI-ISI EXPORT COMPLETE\n');
fprintf('Files exported: %d\n',numel(exported_files));
for k = 1:numel(exported_files)
    fprintf('  %s\n',exported_files(k));
end
if ~isempty(skipped_pairs)
    fprintf('Pairs skipped because no responding data were available: %s\n', ...
        strjoin(cellstr(skipped_pairs),', '));
end
fprintf('Trials were not balanced or downsampled.\n');
fprintf('Bad-trial lists were not exported.\n');
fprintf('Original experiment files modified: NO\n');
fprintf('============================================================\n\n');

%% =========================== FUNCTIONS ================================

function D = load_dataset(folder,electrode_type,FS,require_multi_qc)
files = resolve_dataset_files(folder,require_multi_qc);
if ~files.spike.found
    error('MultiISI:MissingSpike','No unique *sp_xia_SSD.mat in %s.',folder);
end
meta = whos('-file',files.spike.path,'sp_corr');
if isempty(meta) || ~strcmp(meta.class,'cell')
    error('MultiISI:MissingSpCorr','sp_corr cell array is missing from %s.', ...
        files.spike.path);
end
spike = struct('n_channels',prod(double(meta.size)),'size',meta.size);

Eraw = load(files.experiment.path,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required = {'StimParams','simultaneous_stim','E_MAP','n_Trials'};
require_fields(Eraw,required,'experiment file');
experiment = decode_stim_params(Eraw);
trig = load_triggers(folder,experiment.n_trials);
depth_to_hardware = load_depth_map(folder,electrode_type,spike.n_channels);
bad_trials = load_bad_trials(files.bad_trials,experiment.n_trials);
bad_channels = load_bad_channels(files.bad_channels);
responding = load_responding(files.responding);

if require_multi_qc && ~files.responding.found
    error('MultiISI:MissingResponding','No MultiISI responding file in %s.',folder);
end
D = struct('folder',folder,'files',files,'spike',spike, ...
    'experiment',experiment,'trig',trig,'FS',FS, ...
    'depth_to_hardware',depth_to_hardware,'bad_trials',bad_trials, ...
    'bad_channels',bad_channels,'responding',responding);
end

function files = resolve_dataset_files(folder,prefer_multi)
files = struct();
files.spike = unique_file(folder,'*sp_xia_SSD.mat',true,{});
files.experiment = unique_file(folder,'*_exp_datafile_*.mat',true,{});
files.trigger = unique_file(folder,'*.trig.dat',true,{});
if prefer_multi
    files.responding = priority_file(folder,{ ...
        '*_MultiISIRespondingChannels.mat', ...
        '*MultiISIRespondingChannels.mat', ...
        '*_RespondingChannels.mat'},{'BACKUP'});
    files.bad_trials = priority_file(folder,{ ...
        '*.MultiISIsBadTrials.mat','*.MultiSIsBadTrials.mat', ...
        '*._MultiISIsBadTrials.mat','*.SimSeqBadTrials.mat', ...
        '*.BadTrials.mat'},{'BACKUP'});
    files.bad_channels = priority_file(folder,{ ...
        '*.MultiISIsBadChannels.mat','*.MultiSIsBadChannels.mat', ...
        '*.BadChannels.mat'},{'BACKUP'});
else
    files.responding = empty_file_record();
    files.bad_trials = priority_file(folder,{ ...
        '*.BadTrials.mat','*.MultiISIsBadTrials.mat'},{'BACKUP'});
    files.bad_channels = empty_file_record();
end
end

function out = unique_file(folder,pattern,required,exclusions)
listing = filtered_listing(dir(fullfile(folder,pattern)),exclusions);
if isempty(listing)
    if required
        error('MultiISI:MissingFile','No %s file in %s.',pattern,folder);
    end
    out = empty_file_record();
    return;
end
if numel(listing) > 1
    error('MultiISI:AmbiguousFile','Multiple %s files in %s: %s', ...
        pattern,folder,strjoin({listing.name},', '));
end
out = file_record(listing(1),folder);
end

function out = priority_file(folder,patterns,exclusions)
out = empty_file_record();
for k = 1:numel(patterns)
    listing = filtered_listing(dir(fullfile(folder,patterns{k})),exclusions);
    if isempty(listing), continue; end
    if numel(listing) > 1
        error('MultiISI:AmbiguousFile','Multiple %s files in %s: %s', ...
            patterns{k},folder,strjoin({listing.name},', '));
    end
    out = file_record(listing(1),folder);
    return;
end
end

function listing = filtered_listing(listing,exclusions)
if isempty(listing), return; end
keep = ~[listing.isdir];
for k = 1:numel(exclusions)
    keep = keep & ~contains({listing.name},exclusions{k},'IgnoreCase',true);
end
listing = listing(keep);
end

function out = empty_file_record()
out = struct('found',false,'name','','path','');
end

function out = file_record(item,folder)
out = struct('found',true,'name',item.name,'path',fullfile(folder,item.name));
end

function E = decode_stim_params(S)
StimParams = S.StimParams;
simN = double(S.simultaneous_stim);
nTrials = double(S.n_Trials);
if ~iscell(StimParams) || size(StimParams,2) < 16 || ...
        size(StimParams,1) < 1+nTrials*simN
    error('MultiISI:InvalidStimParams','StimParams has an invalid shape.');
end
names = cell(nTrials,simN);
indices = zeros(nTrials,simN);
amps = nan(nTrials,simN);
ptd_us = nan(nTrials,simN);
trial_rows = cell(nTrials,1);
map_names = normalize_names(S.E_MAP(2:end));
for tr = 1:nTrials
    rows = (tr-1)*simN+(2:simN+1);
    trial_rows{tr} = rows;
    names(tr,:) = normalize_names(StimParams(rows,1));
    [tf,idx] = ismember(names(tr,:),map_names);
    idx(~tf) = 0;
    indices(tr,:) = idx;
    amps(tr,:) = numeric_cells(StimParams(rows,16),'amplitude',tr);
    ptd_us(tr,:) = numeric_cells(StimParams(rows,6),'PTD',tr);
end
trial_amp = amps(:,1);
trial_amp(trial_amp == -1) = 0;
if simN > 1
    trial_ptd = ptd_us(:,2)/1000;
else
    trial_ptd = zeros(nTrials,1);
end
[unique_sets,~,set_index] = unique(indices,'rows','stable');
E = struct('n_trials',nTrials,'simultaneous_stim',simN, ...
    'E_MAP',{S.E_MAP},'electrode_names',{map_names}, ...
    'StimParams',{StimParams},'trial_stimparam_rows',{trial_rows}, ...
    'stim_names_per_trial',{names},'stim_indices_per_trial',indices, ...
    'amp_per_event',amps,'trial_amplitude',trial_amp, ...
    'trial_ptd_ms',trial_ptd,'unique_sets',unique_sets, ...
    'set_index',set_index,'amplitudes',unique(trial_amp), ...
    'ptds_ms',unique(trial_ptd));
end

function names = normalize_names(values)
names = cell(size(values));
for k = 1:numel(values)
    if isempty(values{k})
        names{k} = '';
    else
        names{k} = strtrim(char(string(values{k})));
    end
end
end

function values = numeric_cells(cells_in,label,trial)
values = nan(1,numel(cells_in));
for k = 1:numel(cells_in)
    if ~isnumeric(cells_in{k}) || ~isscalar(cells_in{k})
        error('MultiISI:NonNumericParameter', ...
            'Non-numeric %s at trial %d, event %d.',label,trial,k);
    end
    values(k) = double(cells_in{k});
end
end

function trig = load_triggers(folder,nTrials)
old = pwd;
cleanup = onCleanup(@() cd(old));
cd(folder);
trig = loadTrig(0);
clear cleanup;
cd(old);
if numel(trig) ~= nTrials
    error('MultiISI:TriggerCount', ...
        'Trigger count %d does not match n_Trials %d in %s.', ...
        numel(trig),nTrials,folder);
end
trig = double(trig(:));
end

function depth = load_depth_map(folder,electrode_type,nSpikeChannels)
old = pwd;
cleanup = onCleanup(@() cd(old));
cd(folder);
depth = Depth_s(electrode_type);
clear cleanup;
cd(old);
depth = double(depth(:).');
if any(depth < 1) || any(depth > nSpikeChannels)
    error('MultiISI:DepthMap','Depth_s contains an invalid sp_corr index.');
end
end

function B = load_bad_trials(file,nTrials)
B = struct('file_found',false,'mode','none','per_channel',{{}},'global_trials',[]);
if ~file.found, return; end
S = load(file.path);
if ~isfield(S,'BadTrials')
    error('MultiISI:BadTrialsVariable','BadTrials is missing from %s.',file.path);
end
value = S.BadTrials;
if isnumeric(value)
    global_ids = normalize_trial_ids(value,nTrials,file.path);
    B = struct('file_found',true,'mode','global','per_channel',{{}}, ...
        'global_trials',global_ids);
    return;
end
if ~iscell(value)
    error('MultiISI:BadTrialsType','BadTrials must be numeric or a cell array.');
end
normalized = cell(size(value));
for k = 1:numel(value)
    normalized{k} = normalize_trial_ids(value{k},nTrials,file.path);
end
if isempty(normalized)
    mode = 'none'; global_ids = [];
elseif all(cellfun(@(x) isequal(x,normalized{1}),normalized))
    mode = 'global'; global_ids = normalized{1};
else
    mode = 'channel_specific'; global_ids = [];
end
B = struct('file_found',true,'mode',mode, ...
    'per_channel',{normalized},'global_trials',global_ids);
end

function ids = normalize_trial_ids(value,nTrials,source)
if isempty(value)
    ids = [];
elseif isnumeric(value)
    ids = unique(double(value(:).'));
else
    error('MultiISI:BadTrialEntry','A BadTrials entry is not numeric in %s.',source);
end
if any(ids < 1) || any(ids > nTrials) || any(mod(ids,1) ~= 0)
    error('MultiISI:BadTrialRange','BadTrials contains an invalid trial ID in %s.',source);
end
end

function ids = bad_trials_for_channel(B,depth_channel)
switch B.mode
    case 'global'
        ids = B.global_trials;
    case 'channel_specific'
        if depth_channel <= numel(B.per_channel)
            ids = B.per_channel{depth_channel};
        else
            ids = [];
        end
    otherwise
        ids = [];
end
end

function B = load_bad_channels(file)
B = struct('file_found',false,'per_set',{{}});
if ~file.found, return; end
S = load(file.path);
if isfield(S,'BadCh_perSet')
    value = S.BadCh_perSet;
elseif isfield(S,'BadCh')
    value = S.BadCh;
else
    error('MultiISI:BadChannelVariable', ...
        'BadCh_perSet or BadCh is missing from %s.',file.path);
end
if ~iscell(value), value = {value}; end
for k = 1:numel(value)
    if isempty(value{k})
        value{k} = [];
    elseif isnumeric(value{k})
        value{k} = unique(double(value{k}(:).'));
    else
        error('MultiISI:BadChannelType','Bad-channel entry %d is invalid.',k);
    end
end
B = struct('file_found',true,'per_set',{value});
end

function R = load_responding(file)
R = struct('file_found',false,'Responding',[]);
if ~file.found, return; end
S = load(file.path);
if ~isfield(S,'Responding')
    error('MultiISI:RespondingVariable','Responding is missing from %s.',file.path);
end
R = struct('file_found',true,'Responding',S.Responding);
end

function pairs = discover_pairs(E,positive_ptds)
keep_ptd = abs(E.trial_ptd_ms) <= 1e-6;
for p = positive_ptds
    keep_ptd = keep_ptd | abs(E.trial_ptd_ms-p) <= 1e-6;
end
keys = {}; As = {}; Bs = {};
for tr = find(keep_ptd).'
    names = E.stim_names_per_trial(tr,1:2);
    if any(cellfun(@isempty,names)) || strcmp(names{1},names{2}), continue; end
    sorted = sort(string(names));
    A = char(sorted(1)); B = char(sorted(2));
    key = [A ' + ' B];
    if ~ismember(key,keys)
        keys{end+1} = key; As{end+1} = A; Bs{end+1} = B; %#ok<AGROW>
    end
end
pairs = repmat(struct('key','','A','','B',''),1,numel(keys));
for k = 1:numel(keys)
    pairs(k) = struct('key',keys{k},'A',As{k},'B',Bs{k});
end
end

function D = build_condition_descriptors(Esingle,Emulti,A,B,positive_ptds,tol)
D = repmat(empty_descriptor(),1,0);
trA = find(strcmp(Esingle.stim_names_per_trial(:,1),A));
trB = find(strcmp(Esingle.stim_names_per_trial(:,1),B));
if ~isempty(trA)
    D(end+1) = make_descriptor('A','A','single','single',{A},0,0,trA); %#ok<AGROW>
end
if ~isempty(trB)
    D(end+1) = make_descriptor('B','B','single','single',{B},0,0,trB); %#ok<AGROW>
end
trAB = find_pair_trials(Emulti,A,B,0,tol,'unordered');
if ~isempty(trAB)
    D(end+1) = make_descriptor('AB','AB','simultaneous','multiisi', ...
        {A,B},0,[0 0],trAB); %#ok<AGROW>
end
for p = positive_ptds
    tr = find_pair_trials(Emulti,A,B,p,tol,'A_to_B');
    if ~isempty(tr)
        id = ['A_to_B_PTD_' number_tag(p) 'ms'];
        D(end+1) = make_descriptor(id,'A_to_B','sequential','multiisi', ...
            {A,B},p,[0 p],tr); %#ok<AGROW>
    end
    tr = find_pair_trials(Emulti,A,B,p,tol,'B_to_A');
    if ~isempty(tr)
        id = ['B_to_A_PTD_' number_tag(p) 'ms'];
        D(end+1) = make_descriptor(id,'B_to_A','sequential','multiisi', ...
            {B,A},p,[0 p],tr); %#ok<AGROW>
    end
end
end

function D = make_descriptor(id,code,type,source,order,ptd,pulses,trials)
D = empty_descriptor();
D.condition_id = id; D.code = code; D.stimulation_type = type;
D.source_dataset = source; D.electrode_order = order;
D.PTD_ms = ptd; D.pulse_times_ms = pulses; D.trial_ids_all = trials(:).';
if strcmp(code,'A'), D.label = 'A alone';
elseif strcmp(code,'B'), D.label = 'B alone';
elseif strcmp(code,'AB'), D.label = 'A+B simultaneous';
elseif strcmp(code,'A_to_B'), D.label = sprintf('A then B, %g ms',ptd);
else, D.label = sprintf('B then A, %g ms',ptd);
end
end

function D = empty_descriptor()
D = struct('condition_id','','code','','label','', ...
    'stimulation_type','','source_dataset','','electrode_order',{{}}, ...
    'PTD_ms',NaN,'pulse_times_ms',[],'trial_ids_all',[]);
end

function trials = find_pair_trials(E,A,B,ptd,tol,mode)
n1 = E.stim_names_per_trial(:,1); n2 = E.stim_names_per_trial(:,2);
ptd_ok = abs(E.trial_ptd_ms-ptd) <= tol;
switch mode
    case 'unordered'
        name_ok = (strcmp(n1,A)&strcmp(n2,B)) | (strcmp(n1,B)&strcmp(n2,A));
    case 'A_to_B'
        name_ok = strcmp(n1,A)&strcmp(n2,B);
    case 'B_to_A'
        name_ok = strcmp(n1,B)&strcmp(n2,A);
end
trials = find(ptd_ok & name_ok);
end

function trials = pair_trials_for_descriptors(D)
trials = [];
for k = 1:numel(D)
    if ~strcmp(D(k).source_dataset,'single')
        trials = [trials D(k).trial_ids_all]; %#ok<AGROW>
    end
end
trials = unique(trials);
end

function amps = selected_amplitudes(observed,requested,tol)
amps = unique(double(observed(:).'));
if ~isempty(requested)
    keep = false(size(amps));
    for k = 1:numel(requested)
        keep = keep | abs(amps-requested(k)) <= tol;
    end
    amps = amps(keep);
end
end

function bad = bad_channels_for_sets(B,set_indices)
bad = [];
if ~B.file_found, return; end
for si = set_indices(:).'
    if si >= 1 && si <= numel(B.per_set)
        bad = [bad B.per_set{si}]; %#ok<AGROW>
    end
end
bad = unique(bad);
end

function mask = responding_mask(D,set_indices,amp,ptd,amp_tol,ptd_tol,nDepth)
mask = false(1,nDepth);
E = D.experiment; Resp = D.responding.Responding;
ai = find(abs(E.amplitudes-amp) <= amp_tol,1);
pi = find(abs(E.ptds_ms-ptd) <= ptd_tol,1);
if isempty(ai) || isempty(pi), return; end
for si = set_indices(:).'
    if si > numel(Resp.set) || ai > numel(Resp.set(si).amp) || ...
            ~isfield(Resp.set(si).amp(ai),'ptd') || ...
            pi > numel(Resp.set(si).amp(ai).ptd) || ...
            ~isfield(Resp.set(si).amp(ai).ptd(pi),'channel')
        continue;
    end
    channels = Resp.set(si).amp(ai).ptd(pi).channel;
    for ch = 1:min(nDepth,numel(channels))
        if isfield(channels(ch),'is_responsive') && channels(ch).is_responsive
            mask(ch) = true;
        end
    end
end
end

function stim = map_stimulation_contact(name,type,nDepth,depth_map)
P = ProbeMAP;
if type == 0, col = 3; elseif type == 1, col = 5; else, col = 6; end
map_names = P(2:nDepth+1,col);
matches = find(strcmp(map_names,name));
if numel(matches) ~= 1
    error('MultiISI:StimMap','Expected one ProbeMAP match for %s.',name);
end
depth_ch = matches(1);
sp_index = depth_map(depth_ch);
expected = native_name_to_index(name);
if sp_index ~= expected
    error('MultiISI:StimHardwareMap', ...
        '%s maps to sp_corr cell %d, but its name implies %d.', ...
        name,sp_index,expected);
end
[shank,site,x,y] = depth_coordinates(depth_ch,type,nDepth);
stim = struct('native_name',name,'depth_channel',depth_ch, ...
    'sp_cell_index',sp_index,'shank',shank,'local_site',site, ...
    'x_um',x,'y_um',y);
end

function index = native_name_to_index(name)
t = regexp(name,'^([A-D])-([0-9]+)$','tokens','once');
if isempty(t), error('MultiISI:NativeName','Cannot parse %s.',name); end
index = (double(t{1})-double('A'))*32+str2double(t{2})+1;
end

function [shank,site,x,y] = depth_coordinates(ch,type,nDepth)
if type == 0 || type == 1
    if ch < 1 || ch > 32, error('MultiISI:Geometry','Invalid single-shank channel.'); end
    shank = 1; site = ch; x = 0; y = (site-1)*50;
else
    if nDepth > 64 || ch < 1 || ch > 64
        error('MultiISI:Geometry','Four-shank geometry expects channels 1:64.');
    end
    if ch <= 16, shank=1; site=ch;
    elseif ch <= 32, shank=4; site=ch-16;
    elseif ch <= 48, shank=2; site=ch-32;
    else, shank=3; site=ch-48;
    end
    x = (shank-1)*200; y = (site-1)*50;
end
end

function T = build_channel_table(depth_channels,Single,Multi,type,stimA,stimB)
n = numel(depth_channels);
ChannelIndex = (1:n).';
DepthChannel = depth_channels(:);
SpCellIndex_Single = Single.depth_to_hardware(depth_channels).';
SpCellIndex_MultiISI = Multi.depth_to_hardware(depth_channels).';
Shank=zeros(n,1); LocalSite=zeros(n,1); X_um=zeros(n,1); Y_um=zeros(n,1);
DistanceToA_um=zeros(n,1); DistanceToB_um=zeros(n,1); MinimumStimDistance_um=zeros(n,1);
for k = 1:n
    [Shank(k),LocalSite(k),X_um(k),Y_um(k)] = ...
        depth_coordinates(DepthChannel(k),type,numel(Multi.depth_to_hardware));
    DistanceToA_um(k) = euclidean_distance(X_um(k),Y_um(k),stimA.x_um,stimA.y_um);
    DistanceToB_um(k) = euclidean_distance(X_um(k),Y_um(k),stimB.x_um,stimB.y_um);
    MinimumStimDistance_um(k) = min(DistanceToA_um(k),DistanceToB_um(k));
end
T = table(ChannelIndex,DepthChannel,SpCellIndex_Single,SpCellIndex_MultiISI, ...
    Shank,LocalSite,X_um,Y_um,DistanceToA_um,DistanceToB_um,MinimumStimDistance_um);
end

function selected = load_selected_spike_times(file,indices)
meta = whos('-file',file,'sp_corr');
indices = double(indices(:).');
selected = cell(1,numel(indices));
partial_ok = false;
try
    M = matfile(file);
    for k = 1:numel(indices)
        if meta.size(1) == 1
            one = M.sp_corr(1,indices(k));
        else
            one = M.sp_corr(indices(k),1);
        end
        selected{k} = one{1};
    end
    partial_ok = true;
catch ME
    warning('MultiISI:PartialSpikeLoad', ...
        ['Per-channel loading failed (%s). Loading the complete sp_corr ' ...
         'variable for this dataset.'],ME.message);
end
if ~partial_ok
    S = load(file,'sp_corr');
    selected = S.sp_corr(indices);
    selected = selected(:).';
end
for k = 1:numel(indices)
    matrix = selected{k};
    if isempty(matrix), selected{k} = zeros(0,1);
    elseif ~isnumeric(matrix) || size(matrix,2)<1
        error('MultiISI:SpikeMatrix','sp_corr{%d} is invalid.',indices(k));
    else, selected{k} = double(matrix(:,1));
    end
end
end

function C = extract_clean_channel(abs_times,trig,trials,FS,stored_win, ...
        baseline_win,whole_win,early_win,second_relative,ptd,is_seq, ...
        channel_index,depth_channel)
n = numel(trials);
spikes = cell(n,1);
baseline_count=zeros(n,1); baseline_rate=zeros(n,1);
whole_count=zeros(n,1); whole_rate=zeros(n,1);
early_count=zeros(n,1); early_rate=zeros(n,1);
if is_seq, second_count=zeros(n,1); second_rate=zeros(n,1);
else, second_count=nan(n,1); second_rate=nan(n,1); end
for it = 1:n
    t0 = trig(trials(it))/FS*1000;
    keep = abs_times >= t0+stored_win(1) & abs_times <= t0+stored_win(2);
    rel = abs_times(keep)-t0;
    spikes{it} = rel(:);
    baseline_count(it) = count_window(rel,baseline_win);
    whole_count(it) = count_window(rel,whole_win);
    early_count(it) = count_window(rel,early_win);
    if is_seq
        second_count(it) = count_window(rel,ptd+second_relative);
    end
end
baseline_rate = baseline_count/(diff(baseline_win)/1000);
whole_rate = whole_count/(diff(whole_win)/1000);
early_rate = early_count/(diff(early_win)/1000);
if is_seq, second_rate = second_count/(diff(second_relative)/1000); end
expected_baseline_whole = baseline_rate*(diff(whole_win)/1000);
baseline_corrected_whole = whole_count-expected_baseline_whole;
C = empty_channel_data();
C.channel_index = channel_index; C.depth_channel = depth_channel;
C.n_trials = n; C.source_trial_ids = trials(:); C.spike_times_ms = spikes;
C.trial_metrics = struct('baseline_count',baseline_count, ...
    'baseline_rate_hz',baseline_rate,'whole_response_count',whole_count, ...
    'whole_response_rate_hz',whole_rate,'early_window_count',early_count, ...
    'early_window_rate_hz',early_rate,'second_aligned_count',second_count, ...
    'second_aligned_rate_hz',second_rate, ...
    'expected_baseline_whole_count',expected_baseline_whole, ...
    'baseline_corrected_whole_count',baseline_corrected_whole);
C.summary = struct('mean_whole_response_count',mean_or_nan(whole_count), ...
    'mean_baseline_corrected_whole_count',mean_or_nan(baseline_corrected_whole), ...
    'mean_early_window_count',mean_or_nan(early_count), ...
    'mean_second_aligned_count',mean_or_nan(second_count));
end

function n = count_window(times,win)
n = sum(times >= win(1) & times < win(2));
end

function value = mean_or_nan(values)
if isempty(values) || all(isnan(values)), value = NaN;
else, value = mean(values,'omitnan'); end
end

function Events = collect_stimulation_events(E,trials,order,pulse_times)
if isempty(trials), Events = struct([]); return; end
nEvents = numel(order); nParams = size(E.StimParams,2);
Events = repmat(struct('pulse_number',[],'electrode','', ...
    'pulse_time_ms',NaN,'representative_parameter_row',{{}}, ...
    'unique_parameter_rows',{{}},'n_unique_parameter_rows',0, ...
    'constant_across_trials',false),1,nEvents);
for ie = 1:nEvents
    rows = cell(numel(trials),nParams);
    for it = 1:numel(trials)
        tr = trials(it);
        match = find(strcmp(E.stim_names_per_trial(tr,:),order{ie}),1);
        if isempty(match)
            error('MultiISI:StimEvent','Trial %d lacks electrode %s.',tr,order{ie});
        end
        source_rows = E.trial_stimparam_rows{tr};
        rows(it,:) = E.StimParams(source_rows(match),:);
    end
    unique_rows = unique_cell_rows(rows);
    Events(ie).pulse_number = ie; Events(ie).electrode = order{ie};
    Events(ie).pulse_time_ms = pulse_times(ie);
    Events(ie).representative_parameter_row = rows(1,:);
    Events(ie).unique_parameter_rows = unique_rows;
    Events(ie).n_unique_parameter_rows = size(unique_rows,1);
    Events(ie).constant_across_trials = size(unique_rows,1)==1;
end
end

function out = unique_cell_rows(in)
keep = false(size(in,1),1); representatives = zeros(0,1);
for r = 1:size(in,1)
    duplicate = false;
    for j = 1:numel(representatives)
        if isequaln(in(r,:),in(representatives(j),:)), duplicate=true; break; end
    end
    if ~duplicate, keep(r)=true; representatives(end+1,1)=r; end %#ok<AGROW>
end
out = in(keep,:);
end

function C = empty_condition()
C = struct('index',[],'condition_id','','code','','label','', ...
    'stimulation_type','','source_dataset','','electrode_order',{{}}, ...
    'PTD_ms',NaN,'pulse_times_ms',[],'amplitude',[]);
end

function A = empty_amplitude()
A = struct('amplitude_uA',NaN,'channel',[],'stimulation_events',struct([]), ...
    'n_trials_per_channel',[],'minimum_clean_trials',0,'maximum_clean_trials',0);
end

function C = empty_channel_data()
C = struct('channel_index',NaN,'depth_channel',NaN,'n_trials',0, ...
    'source_trial_ids',[],'spike_times_ms',{{}}, ...
    'trial_metrics',struct(),'summary',struct());
end

function validate_model_data(M)
if ~strcmp(M.format_name,'MultiISIStimulationModelData') || ...
        isempty(M.conditions) || height(M.channels)<1
    error('MultiISI:InvalidPackage','The final ModelData structure is incomplete.');
end
nChannels = height(M.channels);
for ic = 1:numel(M.conditions)
    for ia = 1:numel(M.conditions(ic).amplitude)
        B = M.conditions(ic).amplitude(ia);
        if numel(B.channel) ~= nChannels
            error('MultiISI:ChannelBlockCount','Condition channel count mismatch.');
        end
        for ch = 1:nChannels
            C = B.channel(ch);
            if numel(C.spike_times_ms) ~= C.n_trials || ...
                    numel(C.source_trial_ids) ~= C.n_trials
                error('MultiISI:TrialShape','Clean trial data shape mismatch.');
            end
            names = fieldnames(C.trial_metrics);
            for k = 1:numel(names)
                if numel(C.trial_metrics.(names{k})) ~= C.n_trials
                    error('MultiISI:MetricShape','Metric %s has the wrong size.',names{k});
                end
            end
        end
    end
end
end

function names = normalize_parameter_names(header)
names = cell(size(header));
for k = 1:numel(header)
    if isempty(header{k}), names{k}=sprintf('Parameter_%d',k);
    else, names{k}=strtrim(char(string(header{k}))); end
end
names = names(:).';
end

function tag = number_tag(value)
tag = regexprep(sprintf('%.12g',value),'\.','p');
tag = regexprep(tag,'-','m');
end

function id = derive_dataset_id(single_folder,multi_folder)
animal = regexp(single_folder,'DX\d+','match','once');
if isempty(animal), animal='Dataset'; end
[~,multi_name] = fileparts(multi_folder);
id = [safe_name(animal) '_' safe_name(multi_name)];
end

function value = safe_name(value)
value = regexprep(char(string(value)),'[^A-Za-z0-9]+','');
if isempty(value), value='Unknown'; end
end

function name = file_name_only(path_value)
[~,name,ext] = fileparts(path_value);
name = [name ext];
end

function distance = euclidean_distance(x1,y1,x2,y2)
distance = sqrt((x1-x2).^2+(y1-y2).^2);
end

function validate_folder(folder,name)
if ~isfolder(folder), error('MultiISI:Folder','%s does not exist:\n%s',name,folder); end
end

function validate_window(value,name)
if ~isnumeric(value) || numel(value)~=2 || any(~isfinite(value)) || value(2)<=value(1)
    error('MultiISI:Window','%s must be [start end], with end>start.',name);
end
end

function require_fields(S,names,label)
for k = 1:numel(names)
    if ~isfield(S,names{k}), error('MultiISI:MissingField','%s.%s is missing.',label,names{k}); end
end
end

function text = vector_text(values)
if isempty(values), text='(none)'; else, text=strtrim(num2str(values(:).')); end
end

function tf = path_is_inside(candidate,parent)
candidate=normalize_path(candidate); parent=normalize_path(parent);
tf=strcmp(candidate,parent)||startsWith([candidate filesep],[parent filesep]);
end

function value = normalize_path(value)
value=char(string(value));
while numel(value)>1 && value(end)==filesep, value(end)=[]; end
end
