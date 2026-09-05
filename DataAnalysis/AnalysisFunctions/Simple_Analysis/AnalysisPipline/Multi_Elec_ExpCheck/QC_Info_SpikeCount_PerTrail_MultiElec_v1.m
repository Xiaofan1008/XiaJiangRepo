%% ========================================================================
%  Spike Count Summary Extractor for Mixed-Prefix Multi-Electrode Stimulation
%
%  Purpose:
%    Extract baseline and post-stimulation spike counts per trial and
%    per recording channel.
%  Input spike file priority:
%    1) *.sp_xia_SSD.mat              variable: sp_corr
%    2) *.sp_xia_PrefixRecovery.mat   variable: sp_seq
%    3) *.sp_xia.mat                  variable: sp_clipped
%  ConditionType:
%       0 = zero-control
%       1 = sequential prefix/recovery
%       2 = full simultaneous
%  Output structure:
%       SpikeCounts.set(si).amp(ai).prefix(pi).trial(k)
%       SpikeCounts.set(si).amp(ai).sim.trial(k)
% ========================================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ========================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

Electrode_Type = 2;       % 0 rigid, 1 single-shank flex, 2 four-shank flex

raster_chn_start = 1;     % Depth_s index start
raster_chn_end   = 64;    % Depth_s index end

FS = 30000;

%% ====================== COUNTING WINDOWS =====================

% Baseline window relative to trigger.
baseline_win_ms = [-50 -5];

% ResponseWindowMode:
%   'fixed':
%       post window = fixed_post_win_ms
%
%   'train':
%       post window = [response_start_ms, final artifact + response_after_final_pulse_ms]
%
%   'after_final':
%       post window = [final artifact + response_after_final_start_ms,
%                      final artifact + response_after_final_end_ms]
%
% Recommendation for bad-trial/burst checking:
%   Use 'fixed' with [2 40] so all conditions use the same counting window.

ResponseWindowMode = 'fixed';

% Used when ResponseWindowMode = 'fixed'
fixed_post_win_ms = [2 40];

% Used when ResponseWindowMode = 'train'
response_start_ms = 2;
response_after_final_pulse_ms = 20;

% Used when ResponseWindowMode = 'after_final'
response_after_final_start_ms = 2;
response_after_final_end_ms   = 20;

%% ====================== CONDITION FILTERS ====================

% Which condition types to process:
%   1 = sequential prefix/recovery
%   2 = simultaneous
ConditionTypes_to_process = [1 2];

% Which prefix lengths to process for sequential trials.
% [] means all detected prefixes.
% For full-sequence only, use [5].
PrefixLengths_to_process = [1 2 3 4 5];

% Which ISIs to process.
% [] means all detected ISIs.
ISI_to_process_ms = [5];

% Which set/order IDs to process.
% [] means all detected sets.
SetIDs_to_process = [];

% Which amplitudes to process.
% [] means all non-zero amplitudes.
Amps_to_process = [];

% Full simultaneous active count.
% For 5-electrode experiment, usually 5.
full_prefix_length = 5;

% Whether to pool simultaneous trials across all set/order IDs.
% For bad-trial checking, false is cleaner.
pool_simultaneous_across_sets = false;

%% ====================== RESPONDING FILE SETTINGS ==============

% Use responding-channel file to mark responsive channels in Channel_Results.
use_responding_file = true;

% Which responding-channel file to load:
%   'FullSeqAndSim'
%   'AllPrefixesAndSim'
%
% This should match the responding-channel file you generated.
RespondingFileMode = 'AllPrefixesAndSim';

%% ===================== INITIAL SETUP =========================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end

cd(data_folder);
fprintf('Extracting trial spike counts in folder:\n%s\n\n', data_folder);

% ---- Extract base_name ----
parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');

if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ===================== LOAD SPIKES ===========================

sp = [];

ssd_file    = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
prefix_file = dir(fullfile(data_folder, '*.sp_xia_PrefixRecovery.mat'));
base_file   = dir(fullfile(data_folder, '*.sp_xia.mat'));

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
    elseif isfield(S_sp, 'sp_pca')
        sp = S_sp.sp_pca;
        spike_variable_used = 'sp_pca';
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
    error('No spike file found.');
end

nCh = numel(sp);

fprintf('Loaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);
fprintf('Number of spike channels: %d\n', nCh);

%% ===================== LOAD RESPONDING CHANNELS ==============

Resp = [];
hasResp = false;

if use_responding_file

    resp_file = sprintf('%s_RespondingChannels_%s.mat', base_name, RespondingFileMode);

    if isfile(resp_file)

        tmp = load(resp_file);

        if isfield(tmp, 'Responding')
            Resp = tmp.Responding;
            hasResp = true;
            fprintf('Loaded responding-channel file: %s\n', resp_file);
        else
            warning('Responding file found but Responding struct missing.');
        end

    else
        fprintf('No responding-channel file found: %s\n', resp_file);
        fprintf('Responsive-channel summary will be saved as NaN.\n');
    end
end

%% ===================== LOAD TRIGGERS =========================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ===================== LOAD STIM PARAMS ======================

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;
E_MAP             = S_exp.E_MAP; %#ok<NASGU>

if isfield(S_exp, 'CHN')
    CHN = S_exp.CHN;
else
    CHN = [];
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match nTrig (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ===================== REMOVE HEADER ROW =====================

StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(StimParams_data,1), expected_rows);
end

%% ===================== TRIAL METADATA FROM STIMPARAMS =========

if size(StimParams_data,2) < 30
    error('StimParams does not contain columns 26–30.');
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

%% ===================== AMPLITUDE PER TRIAL ====================

trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

Amps = unique(trialAmps(:));
Amps = sort(Amps(:))';

%% ===================== FINAL ACTIVE ARTIFACT TIME =============

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

fprintf('Detected final active artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ===================== APPLY USER FILTERS ====================

% Set/order IDs.
SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));

if isempty(SetIDs_to_process)
    SetIDs_selected = SetIDs_all;
else
    SetIDs_selected = intersect(SetIDs_all, SetIDs_to_process);
end

% Amplitudes.
all_nonzero_amps = Amps(Amps > 0);

if isempty(Amps_to_process)
    Amps_selected = all_nonzero_amps;
else
    Amps_selected = intersect(all_nonzero_amps, Amps_to_process);
end

% Condition types.
ConditionTypes_all = unique(conditionType_trial(conditionType_trial > 0));

if isempty(ConditionTypes_to_process)
    ConditionTypes_selected = ConditionTypes_all;
else
    ConditionTypes_selected = intersect(ConditionTypes_all, ConditionTypes_to_process);
end

% Prefixes.
Prefix_all = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));
Prefix_all = sort(Prefix_all(:))';

if isempty(PrefixLengths_to_process)
    Prefixes_selected = Prefix_all;
else
    Prefixes_selected = intersect(Prefix_all, PrefixLengths_to_process);
end

% ISIs.
ISI_all = unique(isi_ms_trial(conditionType_trial == 1));

if isempty(ISI_to_process_ms)
    ISIs_selected = ISI_all;
else
    ISIs_selected = intersect(ISI_all, ISI_to_process_ms);
end

fprintf('\nSelected set IDs: ');
disp(SetIDs_selected');

fprintf('Selected amplitudes: ');
disp(Amps_selected');

fprintf('Selected condition types: ');
disp(ConditionTypes_selected');

fprintf('Selected prefixes: ');
disp(Prefixes_selected);

fprintf('Selected ISIs (ms): ');
disp(ISIs_selected');

fprintf('Pool simultaneous across sets: %d\n', pool_simultaneous_across_sets);

%% ===================== ELECTRODE MAP =========================

d = Depth_s(Electrode_Type);
depth_range = raster_chn_start : min(raster_chn_end, numel(d));
nChns = numel(depth_range);

%% ===================== INITIALIZE OUTPUT =====================

SpikeCounts = struct();

SummaryRows = {};
rowCounter = 0;

fprintf('\nStarting spike counting loop...\n');

%% ========================================================================
%  MAIN COUNTING LOOP
%% ========================================================================

for si = 1:numel(SetIDs_selected)

    setID = SetIDs_selected(si);

    %% ----- Set metadata -----
    SpikeCounts.set(si).setID = setID;

    if ~isempty(CHN) && setID <= size(CHN,1)

        stim_channels = CHN(setID,:);
        stim_channels = stim_channels(stim_channels > 0);

        set_label = strjoin(arrayfun(@(x) sprintf('Ch%d', x), ...
            stim_channels, 'UniformOutput', false), ' + ');

        SpikeCounts.set(si).set_name = set_label;
        SpikeCounts.set(si).stim_channels = stim_channels;

    else
        set_label = sprintf('Set%d', setID);

        SpikeCounts.set(si).set_name = set_label;
        SpikeCounts.set(si).stim_channels = [];
    end

    fprintf('\nProcessing Set %d/%d: setID = %d (%s)\n', ...
        si, numel(SetIDs_selected), setID, set_label);

    for ai = 1:numel(Amps_selected)

        amp_val = Amps_selected(ai);

        SpikeCounts.set(si).amp(ai).amp_value = amp_val;
        SpikeCounts.set(si).amp(ai).amp_label = sprintf('%.1f uA', amp_val);

        %% ==============================================================
        %  Sequential prefix/recovery conditions
        %% ==============================================================

        if ismember(1, ConditionTypes_selected)

            for pi = 1:numel(Prefixes_selected)

                prefixVal = Prefixes_selected(pi);

                for ii = 1:numel(ISIs_selected)

                    isiVal = ISIs_selected(ii);

                    trials_this_condition = find(conditionType_trial == 1 & ...
                                                 conditionSetID_trial == setID & ...
                                                 prefixLength_trial == prefixVal & ...
                                                 isi_ms_trial == isiVal & ...
                                                 trialAmps == amp_val);

                    trials_this_condition = trials_this_condition(trials_this_condition <= nTrials_use);

                    cond_str = sprintf('Set %d (%s) | Prefix %d | ISI %.1f ms | Amp %.1f uA', ...
                        setID, set_label, prefixVal, isiVal, amp_val);

                    SpikeCounts.set(si).amp(ai).prefix(pi).condition_summary = cond_str;
                    SpikeCounts.set(si).amp(ai).prefix(pi).prefix_length = prefixVal;
                    SpikeCounts.set(si).amp(ai).prefix(pi).isi_ms = isiVal;
                    SpikeCounts.set(si).amp(ai).prefix(pi).total_trials = numel(trials_this_condition);
                    SpikeCounts.set(si).amp(ai).prefix(pi).trial_ids = trials_this_condition(:)';

                    if isempty(trials_this_condition)
                        continue;
                    end

                    fprintf('  Amp %.1f uA | Prefix %d | ISI %.1f ms | Trials = %d\n', ...
                        amp_val, prefixVal, isiVal, numel(trials_this_condition));

                    for k = 1:numel(trials_this_condition)

                        tr = trials_this_condition(k);

                        if tr > nTrig
                            fprintf('WARNING: Trial %d exists in StimParams but not trigger file. Skipping.\n', tr);
                            continue;
                        end

                        % Count spikes for this trial.
                        [Channel_Results, trialSummary] = count_one_trial( ...
                            sp, nCh, d, depth_range, nChns, ...
                            trig, FS, tr, ...
                            baseline_win_ms, ...
                            ResponseWindowMode, fixed_post_win_ms, ...
                            response_start_ms, response_after_final_pulse_ms, ...
                            response_after_final_start_ms, response_after_final_end_ms, ...
                            lastActivePTD_ms, ...
                            hasResp, Resp, si, ai, pi, [], ...
                            'prefix');

                        % Save nested trial result.
                        SpikeCounts.set(si).amp(ai).prefix(pi).trial(k).relative_trial_ID = k;
                        SpikeCounts.set(si).amp(ai).prefix(pi).trial(k).absolute_trial_ID = tr;
                        SpikeCounts.set(si).amp(ai).prefix(pi).trial(k).finalArtifact_ms = lastActivePTD_ms(tr);
                        SpikeCounts.set(si).amp(ai).prefix(pi).trial(k).Channel_Results = Channel_Results;

                        % Add row to flat summary.
                        rowCounter = rowCounter + 1;

                        SummaryRows(rowCounter,:) = make_summary_row( ...
                            setID, set_label, amp_val, ...
                            1, 'Prefix', prefixVal, isiVal, ...
                            k, tr, lastActivePTD_ms(tr), ...
                            trialSummary);
                    end
                end
            end
        end

        %% ==============================================================
        %  Simultaneous condition
        %% ==============================================================

        if ismember(2, ConditionTypes_selected)

            if pool_simultaneous_across_sets

                trials_this_condition = find(conditionType_trial == 2 & ...
                                             activeCount_trial == full_prefix_length & ...
                                             trialAmps == amp_val);

                simLabel = 'Simultaneous pooled';

            else

                trials_this_condition = find(conditionType_trial == 2 & ...
                                             conditionSetID_trial == setID & ...
                                             activeCount_trial == full_prefix_length & ...
                                             trialAmps == amp_val);

                simLabel = 'Simultaneous';
            end

            trials_this_condition = trials_this_condition(trials_this_condition <= nTrials_use);

            fprintf('  Amp %.1f uA | %s | Trials = %d\n', ...
                amp_val, simLabel, numel(trials_this_condition));

            % Only create sim field if simultaneous trials exist.
            if ~isempty(trials_this_condition)

                cond_str = sprintf('Set %d (%s) | %s | Amp %.1f uA', ...
                    setID, set_label, simLabel, amp_val);

                SpikeCounts.set(si).amp(ai).sim.condition_summary = cond_str;
                SpikeCounts.set(si).amp(ai).sim.total_trials = numel(trials_this_condition);
                SpikeCounts.set(si).amp(ai).sim.trial_ids = trials_this_condition(:)';
                SpikeCounts.set(si).amp(ai).sim.pooled_across_sets = pool_simultaneous_across_sets;

                for k = 1:numel(trials_this_condition)

                    tr = trials_this_condition(k);

                    if tr > nTrig
                        fprintf('WARNING: Trial %d exists in StimParams but not trigger file. Skipping.\n', tr);
                        continue;
                    end

                    [Channel_Results, trialSummary] = count_one_trial( ...
                        sp, nCh, d, depth_range, nChns, ...
                        trig, FS, tr, ...
                        baseline_win_ms, ...
                        ResponseWindowMode, fixed_post_win_ms, ...
                        response_start_ms, response_after_final_pulse_ms, ...
                        response_after_final_start_ms, response_after_final_end_ms, ...
                        lastActivePTD_ms, ...
                        hasResp, Resp, si, ai, [], [], ...
                        'sim');

                    SpikeCounts.set(si).amp(ai).sim.trial(k).relative_trial_ID = k;
                    SpikeCounts.set(si).amp(ai).sim.trial(k).absolute_trial_ID = tr;
                    SpikeCounts.set(si).amp(ai).sim.trial(k).finalArtifact_ms = lastActivePTD_ms(tr);
                    SpikeCounts.set(si).amp(ai).sim.trial(k).Channel_Results = Channel_Results;

                    rowCounter = rowCounter + 1;

                    SummaryRows(rowCounter,:) = make_summary_row( ...
                        setID, set_label, amp_val, ...
                        2, simLabel, NaN, NaN, ...
                        k, tr, lastActivePTD_ms(tr), ...
                        trialSummary);
                end
            end
        end
    end
end

%% ===================== CREATE SUMMARY TABLE ==================

if isempty(SummaryRows)

    TrialSummaryTable = table();

else

    TrialSummaryTable = cell2table(SummaryRows, ...
        'VariableNames', { ...
        'SetID', ...
        'SetLabel', ...
        'Amp_uA', ...
        'ConditionType', ...
        'ConditionLabel', ...
        'PrefixLength', ...
        'ISI_ms', ...
        'RelativeTrialID', ...
        'AbsoluteTrialID', ...
        'FinalArtifact_ms', ...
        'TotalBaseline_AllChannels', ...
        'TotalPost_AllChannels', ...
        'TotalPostMinusBaseline_AllChannels', ...
        'MeanPost_AllChannels', ...
        'MaxPost_OneChannel', ...
        'NumChannelsWithPostSpikes', ...
        'NumResponsiveChannels', ...
        'TotalBaseline_RespChannels', ...
        'TotalPost_RespChannels', ...
        'TotalPostMinusBaseline_RespChannels', ...
        'MeanPost_RespChannels'});

end

%% ===================== SAVE RESULTS ==========================

save_name = sprintf('%s_MixedPrefix_TrialSpikeCounts_PerTrial.mat', base_name);
full_save_path = fullfile(data_folder, save_name);

SpikeCountInfo = struct();
SpikeCountInfo.data_folder = data_folder;
SpikeCountInfo.base_name = base_name;
SpikeCountInfo.spike_file_used = spike_file_used;
SpikeCountInfo.spike_variable_used = spike_variable_used;
SpikeCountInfo.use_responding_file = use_responding_file;
SpikeCountInfo.hasResp = hasResp;
SpikeCountInfo.RespondingFileMode = RespondingFileMode;
SpikeCountInfo.baseline_win_ms = baseline_win_ms;
SpikeCountInfo.ResponseWindowMode = ResponseWindowMode;
SpikeCountInfo.fixed_post_win_ms = fixed_post_win_ms;
SpikeCountInfo.response_start_ms = response_start_ms;
SpikeCountInfo.response_after_final_pulse_ms = response_after_final_pulse_ms;
SpikeCountInfo.response_after_final_start_ms = response_after_final_start_ms;
SpikeCountInfo.response_after_final_end_ms = response_after_final_end_ms;
SpikeCountInfo.ConditionTypes_selected = ConditionTypes_selected;
SpikeCountInfo.Prefixes_selected = Prefixes_selected;
SpikeCountInfo.ISIs_selected = ISIs_selected;
SpikeCountInfo.SetIDs_selected = SetIDs_selected;
SpikeCountInfo.Amps_selected = Amps_selected;
SpikeCountInfo.full_prefix_length = full_prefix_length;
SpikeCountInfo.pool_simultaneous_across_sets = pool_simultaneous_across_sets;
SpikeCountInfo.depth_range = depth_range;
SpikeCountInfo.FS = FS;

save(full_save_path, ...
     'SpikeCounts', ...
     'TrialSummaryTable', ...
     'SpikeCountInfo', ...
     'baseline_win_ms', ...
     'ResponseWindowMode', ...
     'fixed_post_win_ms', ...
     'lastActivePTD_ms', ...
     'lastActivePTD_us', ...
     'Amps_selected', ...
     'SetIDs_selected', ...
     'ConditionTypes_selected', ...
     'Prefixes_selected', ...
     'ISIs_selected', ...
     '-v7.3');

fprintf('\nSuccess! Extracted trial spike counts.\n');
fprintf('Data saved to:\n%s\n', full_save_path);

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function [Channel_Results, trialSummary] = count_one_trial( ...
    sp, nCh, d, depth_range, nChns, ...
    trig, FS, tr, ...
    baseline_win_ms, ...
    ResponseWindowMode, fixed_post_win_ms, ...
    response_start_ms, response_after_final_pulse_ms, ...
    response_after_final_start_ms, response_after_final_end_ms, ...
    lastActivePTD_ms, ...
    hasResp, Resp, si, ai, pi, dummy, conditionKind) %#ok<INUSD>

    t0 = trig(tr) / FS * 1000;

    %% ----- Baseline window -----
    win_base_start = t0 + baseline_win_ms(1);
    win_base_end   = t0 + baseline_win_ms(2);

    %% ----- Post window -----
    finalArtifact_ms = lastActivePTD_ms(tr);

    switch lower(ResponseWindowMode)

        case 'fixed'
            post_start_rel = fixed_post_win_ms(1);
            post_end_rel   = fixed_post_win_ms(2);

        case 'train'
            post_start_rel = response_start_ms;
            post_end_rel   = finalArtifact_ms + response_after_final_pulse_ms;

        case 'after_final'
            post_start_rel = finalArtifact_ms + response_after_final_start_ms;
            post_end_rel   = finalArtifact_ms + response_after_final_end_ms;

        otherwise
            error('Unknown ResponseWindowMode: %s', ResponseWindowMode);
    end

    win_post_start = t0 + post_start_rel;
    win_post_end   = t0 + post_end_rel;

    %% ----- Pre-allocate arrays -----
    depth_idx_arr = zeros(nChns, 1);
    intan_ch_arr  = zeros(nChns, 1);

    base_spikes = zeros(nChns, 1);
    post_spikes = zeros(nChns, 1);
    post_minus_baseline = zeros(nChns, 1);

    is_responsive_arr = false(nChns, 1);

    %% ----- Channel loop -----
    for idxDepth = 1:nChns

        ich = depth_range(idxDepth);
        ch  = d(ich);

        depth_idx_arr(idxDepth) = ich;
        intan_ch_arr(idxDepth)  = ch;

        %% ----- Is responsive? -----
        if hasResp

            if strcmpi(conditionKind, 'prefix')

                is_responsive_arr(idxDepth) = get_resp_status_prefix(Resp, si, ai, pi, ich);

            elseif strcmpi(conditionKind, 'sim')

                is_responsive_arr(idxDepth) = get_resp_status_sim(Resp, si, ai, ich);

            end
        end

        %% ----- Spike count -----
        if ch < 1 || ch > nCh || isempty(sp{ch})
            continue;
        end

        sp_times = sp{ch}(:,1);

        base_spikes(idxDepth) = sum(sp_times >= win_base_start & sp_times < win_base_end);
        post_spikes(idxDepth) = sum(sp_times >= win_post_start & sp_times < win_post_end);

        post_minus_baseline(idxDepth) = post_spikes(idxDepth) - base_spikes(idxDepth);
    end

    %% ----- Channel table -----
    Channel_Results = table( ...
        depth_idx_arr, ...
        intan_ch_arr, ...
        base_spikes, ...
        post_spikes, ...
        post_minus_baseline, ...
        is_responsive_arr, ...
        'VariableNames', { ...
        'Depth_Index', ...
        'Intan_Channel', ...
        'Baseline_Spikes', ...
        'Post_Spikes', ...
        'Post_minus_Baseline', ...
        'Is_Responsive'});

    %% ----- All-channel summary -----
    total_base_all = sum(base_spikes);
    total_post_all = sum(post_spikes);
    total_post_minus_base_all = sum(post_minus_baseline);

    mean_post_all = mean(post_spikes);
    max_post_one_channel = max(post_spikes);
    num_channels_with_post = sum(post_spikes > 0);

    %% ----- Responsive-channel summary -----
    if hasResp && any(is_responsive_arr)

        resp_idx = is_responsive_arr;

        num_resp_channels = sum(resp_idx);

        total_base_resp = sum(base_spikes(resp_idx));
        total_post_resp = sum(post_spikes(resp_idx));
        total_post_minus_base_resp = sum(post_minus_baseline(resp_idx));
        mean_post_resp = mean(post_spikes(resp_idx));

    else

        num_resp_channels = 0;

        total_base_resp = NaN;
        total_post_resp = NaN;
        total_post_minus_base_resp = NaN;
        mean_post_resp = NaN;
    end

    %% ----- Save trial summary -----
    trialSummary = struct();

    trialSummary.TotalBaseline_AllChannels = total_base_all;
    trialSummary.TotalPost_AllChannels = total_post_all;
    trialSummary.TotalPostMinusBaseline_AllChannels = total_post_minus_base_all;
    trialSummary.MeanPost_AllChannels = mean_post_all;
    trialSummary.MaxPost_OneChannel = max_post_one_channel;
    trialSummary.NumChannelsWithPostSpikes = num_channels_with_post;

    trialSummary.NumResponsiveChannels = num_resp_channels;
    trialSummary.TotalBaseline_RespChannels = total_base_resp;
    trialSummary.TotalPost_RespChannels = total_post_resp;
    trialSummary.TotalPostMinusBaseline_RespChannels = total_post_minus_base_resp;
    trialSummary.MeanPost_RespChannels = mean_post_resp;
end

function row = make_summary_row( ...
    setID, set_label, amp_val, ...
    conditionType, conditionLabel, prefixVal, isiVal, ...
    relativeTrialID, absoluteTrialID, finalArtifact_ms, ...
    trialSummary)

    row = { ...
        setID, ...
        set_label, ...
        amp_val, ...
        conditionType, ...
        conditionLabel, ...
        prefixVal, ...
        isiVal, ...
        relativeTrialID, ...
        absoluteTrialID, ...
        finalArtifact_ms, ...
        trialSummary.TotalBaseline_AllChannels, ...
        trialSummary.TotalPost_AllChannels, ...
        trialSummary.TotalPostMinusBaseline_AllChannels, ...
        trialSummary.MeanPost_AllChannels, ...
        trialSummary.MaxPost_OneChannel, ...
        trialSummary.NumChannelsWithPostSpikes, ...
        trialSummary.NumResponsiveChannels, ...
        trialSummary.TotalBaseline_RespChannels, ...
        trialSummary.TotalPost_RespChannels, ...
        trialSummary.TotalPostMinusBaseline_RespChannels, ...
        trialSummary.MeanPost_RespChannels};
end

function isResp = get_resp_status_prefix(Resp, si, ai, pi, ich)

    isResp = false;

    try
        if si <= numel(Resp.set) && ...
           ai <= numel(Resp.set(si).amp) && ...
           isfield(Resp.set(si).amp(ai), 'prefix') && ...
           pi <= numel(Resp.set(si).amp(ai).prefix) && ...
           isfield(Resp.set(si).amp(ai).prefix(pi), 'channel') && ...
           ich <= numel(Resp.set(si).amp(ai).prefix(pi).channel)

            R = Resp.set(si).amp(ai).prefix(pi).channel(ich);

            if isfield(R, 'is_responsive') && R.is_responsive
                isResp = true;
            end
        end
    catch
        isResp = false;
    end
end

function isResp = get_resp_status_sim(Resp, si, ai, ich)

    isResp = false;

    try
        if si <= numel(Resp.set) && ...
           ai <= numel(Resp.set(si).amp) && ...
           isfield(Resp.set(si).amp(ai), 'sim') && ...
           isfield(Resp.set(si).amp(ai).sim, 'channel') && ...
           ich <= numel(Resp.set(si).amp(ai).sim.channel)

            R = Resp.set(si).amp(ai).sim.channel(ich);

            if isfield(R, 'is_responsive') && R.is_responsive
                isResp = true;
            end
        end
    catch
        isResp = false;
    end
end