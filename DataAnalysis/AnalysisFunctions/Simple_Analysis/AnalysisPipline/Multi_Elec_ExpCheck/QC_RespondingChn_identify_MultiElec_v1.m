%% ========================================================================
%  Responding Channel Detection for Mixed-Prefix Multi-Electrode Stimulation
%
%  Purpose:
%    Identify responding recording channels for the new mixed-prefix
%    multi-electrode stimulation experiment.
%  ConditionType:
%       0 = zero-control
%       1 = sequential prefix/recovery
%       2 = full simultaneous
%
%  Output structure:
%       Responding.set(si).amp(ai).prefix(pi).channel(ich)
%       Responding.set(si).amp(ai).sim.channel(ich)   % only created if sim trials exist
%
%  ResponseScope:
%       1 = detect only full sequential prefix + simultaneous
%       2 = detect all prefixes + simultaneous
% ========================================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER INPUT =================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

Electrode_Type = 2;    % 0 = rigid, 1 = single-shank flex, 2 = four-shank flex
FS = 30000;            % sampling rate

%% ================= RESPONSE SCOPE =================
% ResponseScope = 1:
%   Detect responding channels only for:
%       full sequential prefix, e.g. Prefix 5
%       simultaneous condition
%
% ResponseScope = 2:
%   Detect responding channels for:
%       Prefix 1, Prefix 2, Prefix 3, Prefix 4, Prefix 5
%       simultaneous condition

ResponseScope = 1;

% Full sequential prefix length.
% For your 5-electrode experiment, this should usually be 5.
full_prefix_length = 5;

%% ================= CONDITION FILTERS =================

% ISIs to analyze, in ms.
% [] means all detected ISIs.
ISI_to_analyze_ms = [5];

% Stimulation set/order IDs to analyze.
% [] means all detected sets.
SetIDs_to_analyze = [];

% Amplitudes to analyze.
% [] means all non-zero amplitudes.
Amps_to_analyze = [];

% Whether to pool simultaneous trials across all set IDs.
% If true, the same simultaneous result is attached to each set.
% If false, simultaneous trials must match the current setID.
pool_simultaneous_across_sets = false;

%% ================= WINDOW DEFINITIONS =================

% Baseline window relative to trigger.
baseline_win_ms = [-50 -5];

% ResponseWindowMode:
%   'train':
%       response window = [response_start_ms, final artifact + response_after_final_pulse_ms]
%
%   'after_final':
%       response window = [final artifact + response_after_final_start_ms,
%                          final artifact + response_after_final_end_ms]
%
%   'fixed':
%       response window = fixed_post_win_ms

ResponseWindowMode = 'fixed';

% Used when ResponseWindowMode = 'train'
response_start_ms = 2;
response_after_final_pulse_ms = 20;

% Used when ResponseWindowMode = 'after_final'
response_after_final_start_ms = 2;
response_after_final_end_ms   = 20;

% Used when ResponseWindowMode = 'fixed'
fixed_post_win_ms = [2 40];

%% ================= FR RESPONDING RULE PARAMETERS =================

% Robust baseline threshold:
%   responsive if mean_post_FR >= median_baseline_FR + k_SD * robust_SD
k_SD = 3;

% Noise floor prevents silent baseline channels from having threshold = 0.
noise_floor_hz = 2.0;

% Extra robustness thresholds.
min_total_baseline_spikes   = 0;      % across all trials
min_total_post_spikes       = 2;      % across all trials
min_frac_trials_with_spikes = 0.13;   % fraction of trials with >=1 post spike
min_abs_post_FR             = 5;      % mean post FR must be >= this value

%% ================= SAVE DETAIL OPTION =================
% If false, only compact summary values are saved for each channel.
% If true, trial-by-trial FR arrays and rule-pass details are also saved.
%
% Recommendation:
%   false for normal analysis
%   true only when debugging the response detection rule

save_detailed_diagnostics = false;

%% ================= PREPARE FOLDER =================

if ~isfolder(data_folder)
    error('Data folder does not exist: %s', data_folder);
end

cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');

if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ================= LOAD SPIKES =================

sp = [];

ssd_file = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
prefix_file = dir(fullfile(data_folder, '*.sp_xia_PrefixRecovery.mat'));
base_file = dir(fullfile(data_folder, '*.sp_xia.mat'));

if ~isempty(ssd_file)

    spike_file_used = ssd_file(1).name;
    S_sp = load(fullfile(data_folder, spike_file_used));

    if isfield(S_sp, 'sp_corr')
        sp = S_sp.sp_corr;
        spike_variable_used = 'sp_corr';
    elseif isfield(S_sp, 'sp_SSD')
        sp = S_sp.sp_SSD;
        spike_variable_used = 'sp_SSD';
    elseif isfield(S_sp, 'sp_in')
        sp = S_sp.sp_in;
        spike_variable_used = 'sp_in';
    else
        error('SSD file found but no usable spike variable was found.');
    end

elseif ~isempty(prefix_file)

    spike_file_used = prefix_file(1).name;
    S_sp = load(fullfile(data_folder, spike_file_used));

    if isfield(S_sp, 'sp_seq')
        sp = S_sp.sp_seq;
        spike_variable_used = 'sp_seq';
    else
        error('PrefixRecovery file found but variable sp_seq was not found.');
    end

elseif ~isempty(base_file)

    spike_file_used = base_file(1).name;
    S_sp = load(fullfile(data_folder, spike_file_used));

    if isfield(S_sp, 'sp_clipped')
        sp = S_sp.sp_clipped;
        spike_variable_used = 'sp_clipped';
    elseif isfield(S_sp, 'sp')
        sp = S_sp.sp;
        spike_variable_used = 'sp';
    else
        error('Base spike file found but no usable spike variable was found.');
    end

else
    error('No spike file found. Expected *.sp_xia_SSD.mat, *.sp_xia_PrefixRecovery.mat, or *.sp_xia.mat.');
end

nCh = numel(sp);

fprintf('\nLoaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);
fprintf('Number of spike channels: %d\n', nCh);

%% ================= LOAD TRIGGERS =================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
trig_ms = trig / FS * 1000;
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ================= LOAD EXPERIMENT PARAMETERS =================

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;
E_MAP             = S_exp.E_MAP;

if isfield(S_exp, 'CHN')
    CHN = S_exp.CHN;
else
    CHN = [];
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ================= REMOVE HEADER ROW =================

StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;
if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

%% ================= TRIAL METADATA FROM STIMPARAMS =================

if size(StimParams_data,2) < 30
    error('StimParams does not contain columns 26–30. Cannot use mixed-prefix metadata.');
end

firstRow_eachTrial = 1:simultaneous_stim:size(StimParams_data,1);

activeCount_trial    = cell2mat(StimParams_data(firstRow_eachTrial,26));
prefixLength_trial   = cell2mat(StimParams_data(firstRow_eachTrial,27));
isi_ms_trial         = cell2mat(StimParams_data(firstRow_eachTrial,28));
conditionType_trial  = cell2mat(StimParams_data(firstRow_eachTrial,29));
conditionSetID_trial = cell2mat(StimParams_data(firstRow_eachTrial,30));

activeCount_trial    = activeCount_trial(:);
prefixLength_trial   = prefixLength_trial(:);
isi_ms_trial         = isi_ms_trial(:);
conditionType_trial  = conditionType_trial(:);
conditionSetID_trial = conditionSetID_trial(:);

fprintf('\nUsing trial metadata from StimParams columns 26–30.\n');

%% ================= AMPLITUDE PER TRIAL =================

trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

[Amps, ~, ampIdx_all] = unique(trialAmps(:)); %#ok<ASGLU>

%% ================= FINAL ACTIVE ARTIFACT TIME =================

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials

    if conditionType_trial(tr) == 0
        lastActivePTD_us(tr) = 0;
        continue;
    end

    if conditionType_trial(tr) == 1
        nActive_this = prefixLength_trial(tr);
    elseif conditionType_trial(tr) == 2
        nActive_this = activeCount_trial(tr);
    else
        nActive_this = activeCount_trial(tr);
    end

    if isnan(nActive_this) || nActive_this < 1
        lastActivePTD_us(tr) = 0;
        continue;
    end

    nActive_this = min(round(nActive_this), simultaneous_stim);

    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);
    activeRows = rr(1:nActive_this);

    ptd_us = cell2mat(StimParams_data(activeRows,6));
    pulseNum = cell2mat(StimParams_data(activeRows,8));
    pulsePeriod_us = cell2mat(StimParams_data(activeRows,9));

    pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
    pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

    rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowFinalArtifact_us);
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

fprintf('\nDetected final active artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ================= APPLY USER FILTERS =================

% Set/order IDs.
SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));

if isempty(SetIDs_to_analyze)
    SetIDs_selected = SetIDs_all;
else
    SetIDs_selected = intersect(SetIDs_all, SetIDs_to_analyze);
end

% ISIs.
ISI_all = unique(isi_ms_trial(conditionType_trial == 1));

if isempty(ISI_to_analyze_ms)
    ISIs_selected = ISI_all;
else
    ISIs_selected = intersect(ISI_all, ISI_to_analyze_ms);
end

% Amplitudes.
all_nonzero_amps = Amps(Amps > 0);

if isempty(Amps_to_analyze)
    Amps_selected = all_nonzero_amps;
else
    Amps_selected = intersect(all_nonzero_amps, Amps_to_analyze);
end

% Prefixes.
Prefix_all = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));
Prefix_all = sort(Prefix_all(:))';

if ResponseScope == 1

    if ~ismember(full_prefix_length, Prefix_all)
        warning('full_prefix_length = %d was not detected. Using max detected prefix instead.', ...
            full_prefix_length);

        if isempty(Prefix_all)
            error('No sequential prefix trials were detected.');
        end

        Prefixes_selected = max(Prefix_all);
        full_prefix_length = Prefixes_selected;
    else
        Prefixes_selected = full_prefix_length;
    end

    scope_label = 'FullSeqAndSim';

elseif ResponseScope == 2

    Prefixes_selected = Prefix_all;
    scope_label = 'AllPrefixesAndSim';

else
    error('ResponseScope must be 1 or 2.');
end

fprintf('\nSelected set IDs: ');
disp(SetIDs_selected');

fprintf('Selected amplitudes: ');
disp(Amps_selected');

fprintf('Selected ISIs (ms): ');
disp(ISIs_selected');

fprintf('Detected prefixes: ');
disp(Prefix_all);

fprintf('Prefixes to analyze: ');
disp(Prefixes_selected);

fprintf('Response scope: %d (%s)\n', ResponseScope, scope_label);
fprintf('Pool simultaneous across sets: %d\n', pool_simultaneous_across_sets);
fprintf('Save detailed diagnostics: %d\n', save_detailed_diagnostics);

%% ================= DEPTH MAP =================

d = Depth_s(Electrode_Type);

%% ================= INITIALIZE OUTPUT =================

Responding = struct();

RespondingInfo = struct();
RespondingInfo.data_folder = data_folder;
RespondingInfo.base_name = base_name;
RespondingInfo.spike_file_used = spike_file_used;
RespondingInfo.spike_variable_used = spike_variable_used;
RespondingInfo.ResponseScope = ResponseScope;
RespondingInfo.scope_label = scope_label;
RespondingInfo.full_prefix_length = full_prefix_length;
RespondingInfo.pool_simultaneous_across_sets = pool_simultaneous_across_sets;
RespondingInfo.ResponseWindowMode = ResponseWindowMode;
RespondingInfo.baseline_win_ms = baseline_win_ms;
RespondingInfo.response_start_ms = response_start_ms;
RespondingInfo.response_after_final_pulse_ms = response_after_final_pulse_ms;
RespondingInfo.response_after_final_start_ms = response_after_final_start_ms;
RespondingInfo.response_after_final_end_ms = response_after_final_end_ms;
RespondingInfo.fixed_post_win_ms = fixed_post_win_ms;
RespondingInfo.k_SD = k_SD;
RespondingInfo.noise_floor_hz = noise_floor_hz;
RespondingInfo.min_total_baseline_spikes = min_total_baseline_spikes;
RespondingInfo.min_total_post_spikes = min_total_post_spikes;
RespondingInfo.min_frac_trials_with_spikes = min_frac_trials_with_spikes;
RespondingInfo.min_abs_post_FR = min_abs_post_FR;
RespondingInfo.save_detailed_diagnostics = save_detailed_diagnostics;
RespondingInfo.FS = FS;

%% ================= HELPER SETTINGS =================

baseDur_s = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;

fprintf('\nRunning responding-channel detection using FR rule only...\n');

%% ======================================================================
%                   MAIN LOOP: SET × AMP × PREFIX/SIM × CHANNEL
%% ======================================================================

for si = 1:numel(SetIDs_selected)

    setID = SetIDs_selected(si);

    %% ----- Store set metadata -----
    Responding.set(si).setID = setID;

    if ~isempty(CHN) && setID <= size(CHN,1)
        stimVec = CHN(setID,:);
        stimVec = stimVec(stimVec > 0);
        Responding.set(si).stimChannels = stimVec;
        Responding.set(si).stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), ...
            stimVec, 'UniformOutput', false), '→');
    else
        Responding.set(si).stimChannels = [];
        Responding.set(si).stimLabel = sprintf('Set%d', setID);
    end

    fprintf('\nProcessing Set %d/%d: setID = %d (%s)\n', ...
        si, numel(SetIDs_selected), setID, Responding.set(si).stimLabel);

    for ai = 1:numel(Amps_selected)

        ampVal = Amps_selected(ai);

        Responding.set(si).amp(ai).amp_value = ampVal;
        Responding.set(si).amp(ai).amp_index = ai;

        %% ==============================================================
        %  Sequential prefix conditions
        %% ==============================================================

        for pi = 1:numel(Prefixes_selected)

            prefixVal = Prefixes_selected(pi);

            Responding.set(si).amp(ai).prefix(pi).prefix_length = prefixVal;

            trials_prefix = find(conditionType_trial == 1 & ...
                                 conditionSetID_trial == setID & ...
                                 prefixLength_trial == prefixVal & ...
                                 ismember(isi_ms_trial, ISIs_selected) & ...
                                 trialAmps == ampVal);

            trials_prefix = trials_prefix(trials_prefix <= nTrials_use);

            Responding.set(si).amp(ai).prefix(pi).isi_ms_selected = ISIs_selected;
            Responding.set(si).amp(ai).prefix(pi).trial_ids = trials_prefix(:)';

            fprintf('  Amp %.1f uA | Prefix %d | Trials = %d\n', ...
                ampVal, prefixVal, numel(trials_prefix));

            for ich = 1:length(d)

                ch = d(ich);

                result = initialize_channel_result();

                if ch > nCh || isempty(sp{ch}) || isempty(trials_prefix)
                    result.is_responsive = false;
                    Responding.set(si).amp(ai).prefix(pi).channel(ich) = result;
                    continue;
                end

                sp_times = sp{ch}(:,1);

                result = compute_channel_response( ...
                    result, sp_times, trials_prefix, trig_ms, ...
                    baseline_win_ms, baseDur_s, ...
                    ResponseWindowMode, fixed_post_win_ms, ...
                    response_start_ms, response_after_final_pulse_ms, ...
                    response_after_final_start_ms, response_after_final_end_ms, ...
                    lastActivePTD_ms, ...
                    k_SD, noise_floor_hz, ...
                    min_total_baseline_spikes, min_total_post_spikes, ...
                    min_frac_trials_with_spikes, min_abs_post_FR, ...
                    save_detailed_diagnostics);

                Responding.set(si).amp(ai).prefix(pi).channel(ich) = result;
            end
        end

        %% ==============================================================
        %  Simultaneous condition
        %  XIA MODIFICATION:
        %    Only create Responding.set(si).amp(ai).sim if sim trials exist.
        %% ==============================================================

        if pool_simultaneous_across_sets

            trials_sim = find(conditionType_trial == 2 & ...
                              activeCount_trial == full_prefix_length & ...
                              trialAmps == ampVal);

        else

            trials_sim = find(conditionType_trial == 2 & ...
                              conditionSetID_trial == setID & ...
                              activeCount_trial == full_prefix_length & ...
                              trialAmps == ampVal);
        end

        trials_sim = trials_sim(trials_sim <= nTrials_use);

        fprintf('  Amp %.1f uA | Simultaneous | Trials = %d\n', ...
            ampVal, numel(trials_sim));

        % Do not create the sim field if no simultaneous trials exist.
        if ~isempty(trials_sim)

            Responding.set(si).amp(ai).sim.trial_ids = trials_sim(:)';
            Responding.set(si).amp(ai).sim.pooled_across_sets = pool_simultaneous_across_sets;

            for ich = 1:length(d)

                ch = d(ich);

                result = initialize_channel_result();

                if ch > nCh || isempty(sp{ch})
                    result.is_responsive = false;
                    Responding.set(si).amp(ai).sim.channel(ich) = result;
                    continue;
                end

                sp_times = sp{ch}(:,1);

                result = compute_channel_response( ...
                    result, sp_times, trials_sim, trig_ms, ...
                    baseline_win_ms, baseDur_s, ...
                    ResponseWindowMode, fixed_post_win_ms, ...
                    response_start_ms, response_after_final_pulse_ms, ...
                    response_after_final_start_ms, response_after_final_end_ms, ...
                    lastActivePTD_ms, ...
                    k_SD, noise_floor_hz, ...
                    min_total_baseline_spikes, min_total_post_spikes, ...
                    min_frac_trials_with_spikes, min_abs_post_FR, ...
                    save_detailed_diagnostics);

                Responding.set(si).amp(ai).sim.channel(ich) = result;
            end
        end

    end
end

%% ================= SAVE RESULT =================

outfile = sprintf('%s_RespondingChannels_%s.mat', base_name, scope_label);
full_out_path = fullfile(data_folder, outfile);

save(full_out_path, ...
    'Responding', ...
    'RespondingInfo', ...
    'ResponseScope', ...
    'scope_label', ...
    'ResponseWindowMode', ...
    'baseline_win_ms', ...
    'fixed_post_win_ms', ...
    'response_start_ms', ...
    'response_after_final_pulse_ms', ...
    'response_after_final_start_ms', ...
    'response_after_final_end_ms', ...
    'k_SD', ...
    'noise_floor_hz', ...
    'min_total_baseline_spikes', ...
    'min_total_post_spikes', ...
    'min_frac_trials_with_spikes', ...
    'min_abs_post_FR', ...
    'save_detailed_diagnostics', ...
    'Amps', ...
    'Amps_selected', ...
    'SetIDs_selected', ...
    'ISIs_selected', ...
    'Prefixes_selected', ...
    'full_prefix_length', ...
    'pool_simultaneous_across_sets', ...
    'lastActivePTD_ms', ...
    'lastActivePTD_us', ...
    'activeCount_trial', ...
    'prefixLength_trial', ...
    'isi_ms_trial', ...
    'conditionType_trial', ...
    'conditionSetID_trial', ...
    'trialAmps', ...
    '-v7.3');

fprintf('\nSaved responding-channel results:\n%s\n', full_out_path);

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function result = initialize_channel_result()

    result = struct();

    result.mean_baseline_FR = NaN;
    result.sd_baseline_FR_robust = NaN;

    result.mean_post_FR = NaN;
    result.total_post_spikes = NaN;
    result.frac_post_trials = NaN;

    % Keep this as the final field.
    result.is_responsive = false;
end

function result = compute_channel_response( ...
    result, sp_times, trials_this, trig_ms, ...
    baseline_win_ms, baseDur_s, ...
    ResponseWindowMode, fixed_post_win_ms, ...
    response_start_ms, response_after_final_pulse_ms, ...
    response_after_final_start_ms, response_after_final_end_ms, ...
    lastActivePTD_ms, ...
    k_SD, noise_floor_hz, ...
    min_total_baseline_spikes, min_total_post_spikes, ...
    min_frac_trials_with_spikes, min_abs_post_FR, ...
    save_detailed_diagnostics)

    nTr = numel(trials_this);

    FR_baseline = zeros(1,nTr);
    FR_post = zeros(1,nTr);

    post_win_start_ms_all = zeros(1,nTr);
    post_win_end_ms_all = zeros(1,nTr);
    post_win_duration_s_all = zeros(1,nTr);

    %% ================= BASELINE FR =================
    for k = 1:nTr

        tr = trials_this(k);
        t0 = trig_ms(tr);

        base_start = t0 + baseline_win_ms(1);
        base_end   = t0 + baseline_win_ms(2);

        mask_base = sp_times >= base_start & sp_times < base_end;

        FR_baseline(k) = sum(mask_base) / baseDur_s;
    end

    %% ================= POST FR =================
    for k = 1:nTr

        tr = trials_this(k);
        t0 = trig_ms(tr);

        finalArtifact_ms = lastActivePTD_ms(tr);

        switch lower(ResponseWindowMode)

            case 'train'
                post_start_rel = response_start_ms;
                post_end_rel   = finalArtifact_ms + response_after_final_pulse_ms;

            case 'after_final'
                post_start_rel = finalArtifact_ms + response_after_final_start_ms;
                post_end_rel   = finalArtifact_ms + response_after_final_end_ms;

            case 'fixed'
                post_start_rel = fixed_post_win_ms(1);
                post_end_rel   = fixed_post_win_ms(2);

            otherwise
                error('Unknown ResponseWindowMode: %s', ResponseWindowMode);
        end

        if post_end_rel <= post_start_rel
            warning('Invalid post window for trial %d: %.2f to %.2f ms. Setting FR_post to 0.', ...
                tr, post_start_rel, post_end_rel);

            post_win_start_ms_all(k) = post_start_rel;
            post_win_end_ms_all(k) = post_end_rel;
            post_win_duration_s_all(k) = NaN;
            FR_post(k) = 0;
            continue;
        end

        postDur_s = (post_end_rel - post_start_rel) / 1000;

        post_start = t0 + post_start_rel;
        post_end   = t0 + post_end_rel;

        mask_post = sp_times >= post_start & sp_times < post_end;

        FR_post(k) = sum(mask_post) / postDur_s;

        post_win_start_ms_all(k) = post_start_rel;
        post_win_end_ms_all(k) = post_end_rel;
        post_win_duration_s_all(k) = postDur_s;
    end

    %% ================= BASELINE / POST SUMMARY =================
    median_b = median(FR_baseline);
    MAD_b = median(abs(FR_baseline - median_b));
    sd_b = MAD_b / 0.6745;
    sd_b_robust = max(sd_b, noise_floor_hz);

    mean_b = mean(FR_baseline);
    mean_p = mean(FR_post);

    total_baseline_spikes = sum(FR_baseline) * baseDur_s;

    % Since post window duration may differ across trials, calculate total
    % post spikes by summing FR * individual duration.
    total_post_spikes = sum(FR_post .* post_win_duration_s_all, 'omitnan');

    frac_post_trials = mean(FR_post > 0);

    %% ================= FR RESPONDING RULE =================
    pass_FR_rule = (mean_p >= median_b + k_SD * sd_b_robust);
    pass_baseline_spikes = (total_baseline_spikes >= min_total_baseline_spikes);
    pass_post_spikes     = (total_post_spikes     >= min_total_post_spikes);
    pass_frac_trials     = (frac_post_trials      >= min_frac_trials_with_spikes);
    pass_abs_FR          = (mean_p                >= min_abs_post_FR);

    isResp = pass_FR_rule & ...
             pass_baseline_spikes & ...
             pass_post_spikes & ...
             pass_frac_trials & ...
             pass_abs_FR;

    %% ================= SAVE COMPACT SUMMARY =================
    result.mean_baseline_FR = mean_b;
    result.sd_baseline_FR_robust = sd_b_robust;

    result.mean_post_FR = mean_p;
    result.total_post_spikes = total_post_spikes;
    result.frac_post_trials = frac_post_trials;

    %% ================= OPTIONAL DETAILED DIAGNOSTICS =================
    if save_detailed_diagnostics

        result.FR_baseline_all = FR_baseline;
        result.FR_post_all = FR_post;

        result.total_baseline_spikes = total_baseline_spikes;

        result.post_win_start_ms_all = post_win_start_ms_all;
        result.post_win_end_ms_all = post_win_end_ms_all;
        result.post_win_duration_s_all = post_win_duration_s_all;

        result.pass_FR_rule = pass_FR_rule;
        result.pass_baseline_spikes = pass_baseline_spikes;
        result.pass_post_spikes = pass_post_spikes;
        result.pass_frac_trials = pass_frac_trials;
        result.pass_abs_FR = pass_abs_FR;
    end

    %% ================= FINAL FIELD =================
    result.is_responsive = isResp;
end