%% Stage B: summarize structurally complete stimulation pairs
% Reads one confirmed Stage A result. Does not load sp_corr and does not
% write to any experiment folder.

clear;
clc;

%% ========================= USER SETTING ===============================
% Enter the full path to the confirmed Stage A result.
stage_a_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/stage_a_check_v5/check_results/DX014_D1_StageA_DatasetInfo.mat'];

%% ======================= INITIAL SETUP ================================
check_results_folder = fileparts(stage_a_file);

if ~isfile(stage_a_file)
    error('StageB:MissingStageAFile', ...
        'Stage A result not found:\n%s', stage_a_file);
end

S = load(stage_a_file, 'StageA');
if ~isfield(S, 'StageA')
    error('StageB:InvalidStageAFile', ...
        'The selected file does not contain StageA: %s', stage_a_file);
end
StageA = S.StageA;

if ~isfield(StageA, 'validation_passed') || ~StageA.validation_passed
    error('StageB:StageAFailed', ...
        'Stage A did not pass. Stage B will not continue.');
end

required_roles = {'single','sim','seq'};
for k = 1:numel(required_roles)
    if ~isfield(StageA.datasets, required_roles{k})
        error('StageB:MissingDataset', ...
            'Stage A is missing the %s dataset.', required_roles{k});
    end
end

% Stage B must write only beside the Stage A result, never in source data.
dataset_folders = {StageA.config.single_folder, StageA.config.sim_folder, ...
    StageA.config.seq_folder};
for k = 1:numel(dataset_folders)
    if path_is_inside(check_results_folder, dataset_folders{k})
        error('StageB:UnsafeOutputFolder', ...
            ['Refusing to write Stage B results because the output folder is ' ...
             'inside an experiment folder.\nOutput: %s\nExperiment: %s'], ...
            check_results_folder, dataset_folders{k});
    end
end

Dsingle = StageA.datasets.single;
Dsim = StageA.datasets.sim;
Dseq = StageA.datasets.seq;
Esingle = Dsingle.experiment;
Esim = Dsim.experiment;
Eseq = Dseq.experiment;

sim_ptd_ms = StageA.config.sim_PTD_ms;
seq_ptd_ms = StageA.config.seq_PTD_ms;
ptd_tol = StageA.config.PTD_tolerance_ms;
amp_tol = 1e-6;

fprintf('\n============================================================\n');
fprintf('STAGE B: STIMULATION-PAIR SUMMARY\n');
fprintf('============================================================\n');
fprintf('Stage A input: %s\n', stage_a_file);
fprintf('Output folder: %s\n', check_results_folder);
fprintf('No experiment files will be modified.\n');

%% ===================== BUILD CANDIDATE PAIRS ==========================
pair_defs = collect_candidate_pairs(Esim, Eseq, sim_ptd_ms, seq_ptd_ms, ptd_tol);
if isempty(pair_defs)
    error('StageB:NoCandidatePairs', ...
        'No two-electrode pairs were found at PTD 0 or %.3f ms.', seq_ptd_ms);
end

nDepth = min([Dsingle.channel_map.n_depth_channels, ...
    Dsim.channel_map.n_depth_channels, Dseq.channel_map.n_depth_channels]);

Pairs = repmat(empty_pair_result(), numel(pair_defs), 1);
trial_count_rows = {};

for ip = 1:numel(pair_defs)
    A = pair_defs(ip).A;
    B = pair_defs(ip).B;
    pair_key = pair_defs(ip).key;

    P = empty_pair_result();
    P.pair_index = ip;
    P.key = pair_key;
    P.electrode_A = A;
    P.electrode_B = B;

    P.single_electrode_indices = struct( ...
        'A', find(strcmp(Esingle.electrode_names, A), 1), ...
        'B', find(strcmp(Esingle.electrode_names, B), 1));
    P.sim_electrode_indices = struct( ...
        'A', find(strcmp(Esim.electrode_names, A), 1), ...
        'B', find(strcmp(Esim.electrode_names, B), 1));
    P.seq_electrode_indices = struct( ...
        'A', find(strcmp(Eseq.electrode_names, A), 1), ...
        'B', find(strcmp(Eseq.electrode_names, B), 1));

    all_amps = collect_pair_amplitudes(Esingle, Esim, Eseq, ...
        A, B, sim_ptd_ms, seq_ptd_ms, ptd_tol);

    % Determine all paired-data set indices associated with this pair.
    sim_pair_trials_all = find_pair_trials(Esim, A, B, sim_ptd_ms, ptd_tol, 'unordered');
    seq_ab_trials_all = find_pair_trials(Eseq, A, B, seq_ptd_ms, ptd_tol, 'A_to_B');
    seq_ba_trials_all = find_pair_trials(Eseq, A, B, seq_ptd_ms, ptd_tol, 'B_to_A');
    sim_set_indices = unique(Esim.set_index(sim_pair_trials_all)).';
    seq_ab_set_indices = unique(Eseq.set_index(seq_ab_trials_all)).';
    seq_ba_set_indices = unique(Eseq.set_index(seq_ba_trials_all)).';

    P.set_indices = struct('sim_0ms',sim_set_indices, ...
        'A_to_B_5ms',seq_ab_set_indices, ...
        'B_to_A_5ms',seq_ba_set_indices);

    bad_depth = unique([ ...
        bad_channels_for_sets(Dsim.bad_channels, sim_set_indices), ...
        bad_channels_for_sets(Dseq.bad_channels, seq_ab_set_indices), ...
        bad_channels_for_sets(Dseq.bad_channels, seq_ba_set_indices)]);
    bad_depth = bad_depth(bad_depth >= 1 & bad_depth <= nDepth);
    common_good_depth = setdiff(1:nDepth, bad_depth, 'stable');

    P.bad_depth_channels_union = bad_depth;
    P.common_good_depth_channels = common_good_depth;
    P.common_good_hardware_channels = ...
        Dsim.channel_map.depth_to_hardware(common_good_depth);

    amp_results = repmat(empty_amplitude_result(), numel(all_amps), 1);
    pair_resp_union = false(1,nDepth);
    condition_exists_any = false(1,5);
    has_parameter_problem = false;

    for ia = 1:numel(all_amps)
        amp = all_amps(ia);
        R = empty_amplitude_result();
        R.amplitude_uA = amp;

        trA_raw = find_single_trials(Esingle, A, amp, amp_tol);
        trB_raw = find_single_trials(Esingle, B, amp, amp_tol);
        trAB_raw = filter_trials_by_amp(Esim, sim_pair_trials_all, amp, amp_tol);
        trAtoB_raw = filter_trials_by_amp(Eseq, seq_ab_trials_all, amp, amp_tol);
        trBtoA_raw = filter_trials_by_amp(Eseq, seq_ba_trials_all, amp, amp_tol);

        R.original_trial_ids = struct('A',trA_raw, 'B',trB_raw, ...
            'AB',trAB_raw, 'A_to_B',trAtoB_raw, 'B_to_A',trBtoA_raw);

        % Single trials are currently unreviewed; paired trials use their
        % respective confirmed global bad-trial lists.
        trA = remove_global_bad_trials(trA_raw, Dsingle.bad_trials);
        trB = remove_global_bad_trials(trB_raw, Dsingle.bad_trials);
        trAB = remove_global_bad_trials(trAB_raw, Dsim.bad_trials);
        trAtoB = remove_global_bad_trials(trAtoB_raw, Dseq.bad_trials);
        trBtoA = remove_global_bad_trials(trBtoA_raw, Dseq.bad_trials);

        R.clean_trial_ids = struct('A',trA, 'B',trB, 'AB',trAB, ...
            'A_to_B',trAtoB, 'B_to_A',trBtoA);
        R.clean_trial_count = [numel(trA), numel(trB), numel(trAB), ...
            numel(trAtoB), numel(trBtoA)];
        R.condition_available = R.clean_trial_count > 0;
        condition_exists_any = condition_exists_any | R.condition_available;

        R.sim_parameter_consistent = simultaneous_order_parameters_agree( ...
            Esim, trAB_raw, A, B);
        if all(R.condition_available) && ~R.sim_parameter_consistent
            has_parameter_problem = true;
        end

        R.is_complete = all(R.condition_available) && ...
            R.sim_parameter_consistent;
        if all(R.condition_available)
            R.potential_balanced_count = min(R.clean_trial_count);
        else
            R.potential_balanced_count = 0;
        end

        sim_sets_amp = unique(Esim.set_index(trAB_raw)).';
        seq_ab_sets_amp = unique(Eseq.set_index(trAtoB_raw)).';
        seq_ba_sets_amp = unique(Eseq.set_index(trBtoA_raw)).';
        R.set_indices = struct('sim_0ms',sim_sets_amp, ...
            'A_to_B_5ms',seq_ab_sets_amp, ...
            'B_to_A_5ms',seq_ba_sets_amp);

        resp_sim = responding_mask(Dsim, sim_sets_amp, amp, sim_ptd_ms, ...
            amp_tol, ptd_tol, nDepth);
        resp_ab = responding_mask(Dseq, seq_ab_sets_amp, amp, seq_ptd_ms, ...
            amp_tol, ptd_tol, nDepth);
        resp_ba = responding_mask(Dseq, seq_ba_sets_amp, amp, seq_ptd_ms, ...
            amp_tol, ptd_tol, nDepth);

        % Apply the common-good-channel population to every response mask.
        common_good_mask = false(1,nDepth);
        common_good_mask(common_good_depth) = true;
        resp_sim = resp_sim & common_good_mask;
        resp_ab = resp_ab & common_good_mask;
        resp_ba = resp_ba & common_good_mask;
        resp_union = resp_sim | resp_ab | resp_ba;

        R.responding_depth_channels = struct( ...
            'AB',find(resp_sim), ...
            'A_to_B',find(resp_ab), ...
            'B_to_A',find(resp_ba), ...
            'union',find(resp_union));
        R.responding_count = [sum(resp_sim), sum(resp_ab), ...
            sum(resp_ba), sum(resp_union)];

        if R.is_complete
            pair_resp_union = pair_resp_union | resp_union;
        end

        trial_count_rows(end+1,:) = {ip, pair_key, A, B, amp, ... %#ok<SAGROW>
            R.clean_trial_count(1), R.clean_trial_count(2), ...
            R.clean_trial_count(3), R.clean_trial_count(4), ...
            R.clean_trial_count(5), R.potential_balanced_count, ...
            R.sim_parameter_consistent, R.is_complete, ...
            R.responding_count(1), R.responding_count(2), ...
            R.responding_count(3), R.responding_count(4)};

        amp_results(ia) = R;
    end

    P.amplitudes = amp_results;
    if isempty(amp_results)
        complete_mask = false(0,1);
    else
        complete_mask = [amp_results.is_complete];
    end
    P.complete_amplitudes_uA = [amp_results(complete_mask).amplitude_uA];
    P.n_common_good_channels = numel(common_good_depth);
    P.responding_union_depth_channels = find(pair_resp_union);
    P.n_responding_union = sum(pair_resp_union);
    P.condition_type_exists = condition_exists_any;

    if has_parameter_problem
        P.status = 'QC_PROBLEM';
    elseif ~all(condition_exists_any)
        P.status = 'INCOMPLETE';
    elseif isempty(P.complete_amplitudes_uA)
        P.status = 'NO_COMMON_AMPLITUDE';
    elseif isempty(common_good_depth)
        P.status = 'NO_COMMON_GOOD_CHANNELS';
    elseif P.n_responding_union == 0
        P.status = 'NO_RESPONDING_CHANNELS';
    else
        P.status = 'COMPLETE';
    end

    P.is_structurally_complete = strcmp(P.status, 'COMPLETE') || ...
        strcmp(P.status, 'NO_RESPONDING_CHANNELS');
    P.single_trials_reviewed = Dsingle.bad_trials.file_found;
    Pairs(ip) = P;
end

%% ========================== SUMMARY TABLES ============================
PairIndex = (1:numel(Pairs)).';
PairKey = string({Pairs.key}).';
ElectrodeA = string({Pairs.electrode_A}).';
ElectrodeB = string({Pairs.electrode_B}).';
Status = string({Pairs.status}).';
CompleteAmplitudes_uA = strings(numel(Pairs),1);
for ip = 1:numel(Pairs)
    CompleteAmplitudes_uA(ip) = string( ...
        vector_text(Pairs(ip).complete_amplitudes_uA));
end
CommonGoodChannels = [Pairs.n_common_good_channels].';
RespondingUnion = [Pairs.n_responding_union].';
SingleTrialsReviewed = [Pairs.single_trials_reviewed].';

PairSummary = table(PairIndex, PairKey, ElectrodeA, ElectrodeB, Status, ...
    CompleteAmplitudes_uA, CommonGoodChannels, RespondingUnion, ...
    SingleTrialsReviewed);

trial_variable_names = {'PairIndex','PairKey','ElectrodeA','ElectrodeB', ...
    'Amplitude_uA','N_A','N_B','N_AB','N_AtoB','N_BtoA', ...
    'PotentialBalancedN','SimParametersAgree','IsComplete', ...
    'NResp_AB','NResp_AtoB','NResp_BtoA','NResp_Union'};
TrialCountSummary = cell2table(trial_count_rows, ...
    'VariableNames', trial_variable_names);

%% ============================= SAVE ==================================
StageB = struct();
StageB.created_on = datetime('now');
StageB.source_stage_a_file = stage_a_file;
StageB.config = struct('sim_PTD_ms',sim_ptd_ms, ...
    'seq_PTD_ms',seq_ptd_ms, 'PTD_tolerance_ms',ptd_tol, ...
    'amplitude_tolerance_uA',amp_tol, ...
    'single_trials_reviewed',Dsingle.bad_trials.file_found);
StageB.condition_order = {'A','B','AB','A_to_B','B_to_A'};
StageB.responding_condition_order = {'AB','A_to_B','B_to_A','union'};
StageB.pairs = Pairs;
StageB.PairSummary = PairSummary;
StageB.TrialCountSummary = TrialCountSummary;
StageB.n_candidate_pairs = numel(Pairs);
StageB.n_complete_pairs = sum(strcmp({Pairs.status}, 'COMPLETE'));
StageB.experiment_files_modified = false;

[~, stage_a_stem] = fileparts(stage_a_file);
dataset_tag = strrep(stage_a_stem, '_StageA_DatasetInfo', '');
output_mat = fullfile(check_results_folder, ...
    [dataset_tag '_StageB_PairSummary.mat']);
output_txt = fullfile(check_results_folder, ...
    [dataset_tag '_StageB_Report.txt']);

% Final safety check immediately before the only two write operations.
for k = 1:numel(dataset_folders)
    if path_is_inside(fileparts(output_mat), dataset_folders{k}) || ...
            path_is_inside(fileparts(output_txt), dataset_folders{k})
        error('StageB:UnsafeWriteBlocked', ...
            'Stage B output resolved inside an experiment folder. Save blocked.');
    end
end

save(output_mat, 'StageB', '-v7.3');
write_stage_b_report(output_txt, StageB);

fprintf('\n============================================================\n');
fprintf('STAGE B COMPLETE\n');
fprintf('Candidate pairs: %d\n', StageB.n_candidate_pairs);
fprintf('Complete responsive pairs: %d\n', StageB.n_complete_pairs);
fprintf('MAT result: %s\n', output_mat);
fprintf('Text report: %s\n', output_txt);
fprintf('Experiment files modified: NO\n');
fprintf('============================================================\n\n');

%% =========================== FUNCTIONS ================================
function pair_defs = collect_candidate_pairs(Esim, Eseq, sim_ptd, seq_ptd, tol)
keys = {};
As = {};
Bs = {};
datasets = {Esim, Eseq};
target_ptds = [sim_ptd, seq_ptd];

for id = 1:2
    E = datasets{id};
    if E.simultaneous_stim ~= 2
        error('StageB:UnsupportedStimCount', ...
            'Stage B currently requires exactly two stimulation events/trial.');
    end
    trials = find(abs(E.trial_ptd_ms - target_ptds(id)) <= tol);
    for tr = trials(:).'
        names = E.stim_names_per_trial(tr,1:2);
        if any(cellfun(@isempty,names)) || strcmp(names{1},names{2})
            continue;
        end
        sorted_names = sort(string(names));
        name_A = char(sorted_names(1));
        name_B = char(sorted_names(2));
        key = [name_A ' + ' name_B];
        if ~ismember(key, keys)
            keys{end+1} = key; %#ok<AGROW>
            As{end+1} = name_A; %#ok<AGROW>
            Bs{end+1} = name_B; %#ok<AGROW>
        end
    end
end

pair_defs = repmat(struct('key','','A','','B',''), numel(keys), 1);
for k = 1:numel(keys)
    pair_defs(k).key = keys{k};
    pair_defs(k).A = As{k};
    pair_defs(k).B = Bs{k};
end
end

function amps = collect_pair_amplitudes(Esingle, Esim, Eseq, A, B, sim_ptd, seq_ptd, tol)
trA = find(strcmp(Esingle.stim_names_per_trial(:,1), A));
trB = find(strcmp(Esingle.stim_names_per_trial(:,1), B));
trAB = find_pair_trials(Esim, A, B, sim_ptd, tol, 'unordered');
trAtoB = find_pair_trials(Eseq, A, B, seq_ptd, tol, 'A_to_B');
trBtoA = find_pair_trials(Eseq, A, B, seq_ptd, tol, 'B_to_A');
amps = unique([Esingle.trial_amplitude([trA;trB]); ...
    Esim.trial_amplitude(trAB); Eseq.trial_amplitude([trAtoB;trBtoA])]);
amps = amps(isfinite(amps)).';
end

function trials = find_pair_trials(E, A, B, ptd, tol, mode)
if E.simultaneous_stim < 2
    trials = [];
    return;
end
n1 = E.stim_names_per_trial(:,1);
n2 = E.stim_names_per_trial(:,2);
ptd_mask = abs(E.trial_ptd_ms - ptd) <= tol;

switch mode
    case 'unordered'
        name_mask = (strcmp(n1,A) & strcmp(n2,B)) | ...
            (strcmp(n1,B) & strcmp(n2,A));
    case 'A_to_B'
        name_mask = strcmp(n1,A) & strcmp(n2,B);
    case 'B_to_A'
        name_mask = strcmp(n1,B) & strcmp(n2,A);
    otherwise
        error('StageB:UnknownPairMode', 'Unknown pair mode: %s', mode);
end
trials = find(ptd_mask & name_mask);
end

function trials = find_single_trials(E, electrode_name, amp, tol)
name_mask = strcmp(E.stim_names_per_trial(:,1), electrode_name);
amp_mask = abs(E.trial_amplitude - amp) <= tol;
trials = find(name_mask & amp_mask);
end

function trials = filter_trials_by_amp(E, candidates, amp, tol)
if isempty(candidates)
    trials = [];
else
    keep = abs(E.trial_amplitude(candidates) - amp) <= tol;
    trials = candidates(keep);
end
end

function clean = remove_global_bad_trials(trials, bad_trial_info)
if bad_trial_info.file_found && bad_trial_info.is_global
    clean = setdiff(trials(:).', bad_trial_info.global_trials, 'stable');
else
    clean = trials(:).';
end
end

function bad = bad_channels_for_sets(bad_info, set_indices)
bad = [];
if ~bad_info.file_found || isempty(set_indices)
    return;
end
for si = set_indices(:).'
    if si >= 1 && si <= numel(bad_info.per_set)
        bad = [bad, bad_info.per_set{si}]; %#ok<AGROW>
    end
end
bad = unique(bad);
end

function mask = responding_mask(D, set_indices, amp, ptd, amp_tol, ptd_tol, nDepth)
mask = false(1,nDepth);
if ~D.responding.file_found || isempty(set_indices)
    return;
end

E = D.experiment;
ai = find(abs(E.amplitudes - amp) <= amp_tol, 1);
pi = find(abs(E.ptds_ms - ptd) <= ptd_tol, 1);
if isempty(ai) || isempty(pi)
    return;
end

Resp = D.responding.Responding;
for si = set_indices(:).'
    if si > numel(Resp.set) || ai > numel(Resp.set(si).amp) || ...
            ~isfield(Resp.set(si).amp(ai),'ptd') || ...
            pi > numel(Resp.set(si).amp(ai).ptd) || ...
            ~isfield(Resp.set(si).amp(ai).ptd(pi),'channel')
        continue;
    end
    channels = Resp.set(si).amp(ai).ptd(pi).channel;
    for ich = 1:min(nDepth,numel(channels))
        if isfield(channels(ich),'is_responsive') && ...
                channels(ich).is_responsive
            mask(ich) = true;
        end
    end
end
end

function agrees = simultaneous_order_parameters_agree(E, trials, A, B)
% Compare representative canonical parameter rows when both order labels exist.
if isempty(trials)
    agrees = true;
    return;
end

names1 = E.stim_names_per_trial(trials,1);
names2 = E.stim_names_per_trial(trials,2);
ab = trials(strcmp(names1,A) & strcmp(names2,B));
ba = trials(strcmp(names1,B) & strcmp(names2,A));

% Every simultaneous event should use the amplitude represented by the trial.
trial_amp = E.trial_amplitude(trials);
event_amps = E.amp_per_event(trials,1:2);
amp_agrees = all(abs(event_amps(:,1)-trial_amp) <= 1e-6) && ...
    all(abs(event_amps(:,2)-trial_amp) <= 1e-6);
if ~amp_agrees
    agrees = false;
    return;
end

if isempty(ab) || isempty(ba)
    agrees = true;
    return;
end

params_ab = canonical_parameter_rows(E, ab(1), A, B);
params_ba = canonical_parameter_rows(E, ba(1), A, B);
agrees = isequaln(params_ab, params_ba);
end

function rows_out = canonical_parameter_rows(E, trial_id, A, B)
source_rows = E.trial_stimparam_rows{trial_id};
names = E.stim_names_per_trial(trial_id,:);
idxA = find(strcmp(names,A),1);
idxB = find(strcmp(names,B),1);
if isempty(idxA) || isempty(idxB)
    rows_out = {};
    return;
end
rows_out = E.StimParams(source_rows([idxA idxB]),:);
end

function tf = path_is_inside(candidate, parent)
candidate = normalize_path(candidate);
parent = normalize_path(parent);
tf = strcmp(candidate,parent) || ...
    startsWith([candidate filesep], [parent filesep]);
end

function value = normalize_path(value)
value = char(string(value));
while numel(value) > 1 && value(end) == filesep
    value(end) = [];
end
end

function P = empty_pair_result()
P = struct('pair_index',[], 'key','', 'electrode_A','', 'electrode_B','', ...
    'single_electrode_indices',struct(), 'sim_electrode_indices',struct(), ...
    'seq_electrode_indices',struct(), 'set_indices',struct(), ...
    'bad_depth_channels_union',[], 'common_good_depth_channels',[], ...
    'common_good_hardware_channels',[], 'amplitudes',[], ...
    'complete_amplitudes_uA',[], 'n_common_good_channels',0, ...
    'responding_union_depth_channels',[], 'n_responding_union',0, ...
    'condition_type_exists',false(1,5), 'status','', ...
    'is_structurally_complete',false, 'single_trials_reviewed',false);
end

function R = empty_amplitude_result()
R = struct('amplitude_uA',NaN, 'original_trial_ids',struct(), ...
    'clean_trial_ids',struct(), 'clean_trial_count',zeros(1,5), ...
    'condition_available',false(1,5), ...
    'sim_parameter_consistent',true, 'is_complete',false, ...
    'potential_balanced_count',0, 'set_indices',struct(), ...
    'responding_depth_channels',struct(), 'responding_count',zeros(1,4));
end

function text_value = vector_text(values)
if isempty(values)
    text_value = '(none)';
else
    text_value = strtrim(num2str(values(:).'));
end
end

function write_stage_b_report(filename, StageB)
fid = fopen(filename,'w');
if fid < 0
    error('StageB:ReportWriteFailed', 'Could not write: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

dual_print(fid, '\n################ COMPLETE STAGE B REPORT ################\n');
dual_print(fid, 'Created: %s\n', char(string(StageB.created_on)));
dual_print(fid, 'Stage A source: %s\n', StageB.source_stage_a_file);
dual_print(fid, 'Candidate pairs: %d\n', StageB.n_candidate_pairs);
dual_print(fid, 'Complete responsive pairs: %d\n', StageB.n_complete_pairs);
dual_print(fid, 'Experiment files modified: NO\n\n');

dual_print(fid, '[PAIR OVERVIEW]\n%s\n', evalc('disp(StageB.PairSummary)'));

for ip = 1:numel(StageB.pairs)
    P = StageB.pairs(ip);
    dual_print(fid, '\n------------------------------------------------------------\n');
    dual_print(fid, 'PAIR %d: %s\n', P.pair_index, P.key);
    dual_print(fid, 'Status: %s\n', P.status);
    dual_print(fid, 'Complete amplitudes: %s uA\n', ...
        vector_text(P.complete_amplitudes_uA));
    dual_print(fid, 'Common good channels: %d\n', P.n_common_good_channels);
    dual_print(fid, 'Responding-channel union: %d\n', P.n_responding_union);
    dual_print(fid, 'Single trials reviewed: %s\n', ...
        char(string(P.single_trials_reviewed)));
    dual_print(fid, ['Amp      A      B     A+B    A->B   B->A   ' ...
        'Balance  SimParams  RespUnion  Complete\n']);

    for ia = 1:numel(P.amplitudes)
        R = P.amplitudes(ia);
        dual_print(fid, '%4g %6d %6d %6d %6d %6d %8d %10s %10d %9s\n', ...
            R.amplitude_uA, R.clean_trial_count(1), R.clean_trial_count(2), ...
            R.clean_trial_count(3), R.clean_trial_count(4), ...
            R.clean_trial_count(5), R.potential_balanced_count, ...
            char(string(R.sim_parameter_consistent)), ...
            R.responding_count(4), char(string(R.is_complete)));
    end
end

dual_print(fid, '\n############## END COMPLETE STAGE B REPORT ##############\n');
end

function dual_print(fid, varargin)
fprintf(1, varargin{:});
fprintf(fid, varargin{:});
end
