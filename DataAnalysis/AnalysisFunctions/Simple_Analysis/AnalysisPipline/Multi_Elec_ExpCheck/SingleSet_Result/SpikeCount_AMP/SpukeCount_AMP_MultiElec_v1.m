%% ============================================================
%  Spike Count Analysis for Mixed-Prefix Multi-Electrode Stimulation
%
%  Main comparison:
%       Full sequential stimulation
%       vs
%       Full simultaneous stimulation
%  Sequential condition:
%       conditionType_trial == 1
%       prefixLength_trial  == TargetSeqPrefix
%       isi_ms_trial        == TargetSeqISI_ms
%  Simultaneous condition:
%       conditionType_trial == 2
%       activeCount_trial   == TargetSimActiveCount
%  Output:
%       Result.Info
%       Result.Data
%       Result.Summary
%       Result.Stats
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex
FS = 30000;

%% ================= ANALYSIS TARGET ===========================

% Full sequential condition.
TargetSeqPrefix = 5;
TargetSeqISI_ms = 5;

% Full simultaneous condition.
TargetSimActiveCount = 5;

% Which responding-channel file to use:
%   'FullSeqAndSim'
%   'AllPrefixesAndSim'
RespondingFileMode = 'FullSeqAndSim';

% Include zero amplitude?
IncludeZeroAmp = false;

% Specific amplitudes to analyze.
% [] means all detected non-zero amplitudes.
Plot_Amps = [];

% Specific set IDs to analyze.
% [] means all detected sets.
Plot_SetIDs = [];

%% ================= BAD TRIAL SETTINGS ========================

remove_bad_trials = true;

%% ================= SPIKE COUNT SETTINGS ======================

% Spike counting window, relative to trigger.
post_win_ms = [2 40];

% PSTH window.
psth_win_ms = [-50 100];

bin_ms     = 1;
sigma_bins = 3;

%% ================= PLOT SETTINGS =============================

jitter_width  = 0.2;
scatter_alpha = 0.4;
dot_size      = 20;

fig_position = [100 100 900 650];

%% ================= SAVE SETTINGS =============================

save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Electrode5_Results/SpikeCount_AMP/DX027';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% ================= INITIAL SETUP =============================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end

cd(data_folder);

fprintf('\nRunning mixed-prefix spike count analysis in:\n%s\n\n', data_folder);

%% ================= BASE NAME ================================

parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');

if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ================= LOAD SPIKES ===============================

[sp, spike_file_used, spike_variable_used] = load_spike_file(data_folder);
nCh = numel(sp);

fprintf('Loaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);
fprintf('Number of spike channels: %d\n', nCh);

%% ================= LOAD RESPONDING CHANNELS ==================

resp_file = sprintf('%s_RespondingChannels_%s.mat', base_name, RespondingFileMode);
resp_path = fullfile(data_folder, resp_file);

if ~isfile(resp_path)
    error('Responding-channel file not found:\n%s', resp_path);
end

tmp = load(resp_path);

if ~isfield(tmp, 'Responding')
    error('Responding variable not found in %s', resp_file);
end

Responding = tmp.Responding;

fprintf('Loaded responding-channel file: %s\n', resp_file);

%% ================= LOAD BAD TRIALS ===========================

BadTrials = {};
BadTrialFile = '';

if remove_bad_trials

    bad_file = sprintf('%s.MixedPrefixBadTrials.mat', base_name);
    bad_path = fullfile(data_folder, bad_file);

    if isfile(bad_path)
        tmp = load(bad_path);

        if isfield(tmp, 'BadTrials')
            BadTrials = tmp.BadTrials;
            BadTrialFile = bad_file;
            fprintf('Loaded bad-trial file: %s\n', bad_file);
        else
            warning('Bad-trial file found but BadTrials variable missing. No bad trials removed.');
        end
    else
        fprintf('No bad-trial file found: %s\n', bad_file);
        fprintf('No bad trials will be removed.\n');
    end
else
    fprintf('remove_bad_trials = false. No bad trials removed.\n');
end

%% ================= LOAD TRIGGERS =============================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ================= LOAD EXPERIMENT FILE ======================

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;

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

%% ================= REMOVE HEADER ROW =========================

StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(StimParams_data,1), expected_rows);
end

%% ================= TRIAL METADATA ============================

% StimParams columns:
%   26 = ActiveElectrodeCount
%   27 = PrefixLength
%   28 = ISI_ms
%   29 = ConditionType
%   30 = ConditionSetID

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

%% ================= AMPLITUDES ================================

trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

Amps_all = unique(trialAmps(:));
Amps_all = sort(Amps_all(:))';

if IncludeZeroAmp
    Amps_available = Amps_all;
else
    Amps_available = Amps_all(Amps_all > 0);
end

if isempty(Plot_Amps)
    Amps = Amps_available;
else
    Amps = intersect(Amps_available, Plot_Amps);
end

nAmp = numel(Amps);

fprintf('Selected amplitudes: ');
disp(Amps);

%% ================= SET IDs ==================================

SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));
SetIDs_all = sort(SetIDs_all(:))';

if isempty(Plot_SetIDs)
    SetIDs = SetIDs_all;
else
    SetIDs = intersect(SetIDs_all, Plot_SetIDs);
end

nSet = numel(SetIDs);

fprintf('Selected set IDs: ');
disp(SetIDs);

%% ================= ELECTRODE MAP =============================

d = Depth_s(Electrode_Type);
nCh_Total = length(d);

%% ================= PSTH KERNEL ===============================

edges_psth = psth_win_ms(1):bin_ms:psth_win_ms(2);
ctrs_psth  = edges_psth(1:end-1) + diff(edges_psth)/2;
nBins      = numel(edges_psth) - 1;
bin_s      = bin_ms / 1000;

kernel_size = 2 * ceil(2*sigma_bins) + 1;
g_sym = gausswin(kernel_size);
g_sym = g_sym / sum(g_sym);

%% =============================================================
%  BUILD UNION RESPONDING CHANNELS PER SET
%% =============================================================

resp_channels_per_set = cell(nSet, 1);

fprintf('\nBuilding union responding channels per set:\n');

for ss = 1:nSet

    setID = SetIDs(ss);
    si = find_set_index(Responding, setID);

    local_resp_mask = false(nCh_Total, 1);

    if si == 0
        fprintf('  Set %d not found in Responding. No channels selected.\n', setID);
        resp_channels_per_set{ss} = [];
        continue;
    end

    for ai = 1:numel(Responding.set(si).amp)

        ampVal = Responding.set(si).amp(ai).amp_value;

        if ~ismembertol(ampVal, Amps, 1e-6)
            continue;
        end

        %% ----- Sequential target condition -----
        if isfield(Responding.set(si).amp(ai), 'prefix')

            for pi = 1:numel(Responding.set(si).amp(ai).prefix)

                P = Responding.set(si).amp(ai).prefix(pi);

                if ~is_target_prefix(P, TargetSeqPrefix, TargetSeqISI_ms)
                    continue;
                end

                if isfield(P, 'channel')
                    for ch = 1:min(numel(P.channel), nCh_Total)
                        if isfield(P.channel(ch), 'is_responsive') && P.channel(ch).is_responsive
                            local_resp_mask(ch) = true;
                        end
                    end
                end
            end
        end

        %% ----- Simultaneous target condition, if available -----
        if isfield(Responding.set(si).amp(ai), 'sim') && ...
           isfield(Responding.set(si).amp(ai).sim, 'channel')

            S = Responding.set(si).amp(ai).sim;

            % If active_count field exists, check it. If not, accept sim field.
            if isfield(S, 'active_count')
                if abs(S.active_count - TargetSimActiveCount) > 1e-6
                    continue;
                end
            elseif isfield(S, 'activeCount')
                if abs(S.activeCount - TargetSimActiveCount) > 1e-6
                    continue;
                end
            end

            for ch = 1:min(numel(S.channel), nCh_Total)
                if isfield(S.channel(ch), 'is_responsive') && S.channel(ch).is_responsive
                    local_resp_mask(ch) = true;
                end
            end
        end
    end

    resp_channels_per_set{ss} = find(local_resp_mask);

    setLabel = get_set_label(CHN, setID);
    fprintf('  Set %d (%s): %d responding channels\n', ...
        setID, setLabel, numel(resp_channels_per_set{ss}));
end

%% =============================================================
%  INITIALIZE OUTPUT ARRAYS
%% =============================================================

SC_seq = nan(nCh_Total, nAmp, nSet);
SC_sim = nan(nCh_Total, nAmp, nSet);

PSTH_seq = nan(nCh_Total, nBins, nAmp, nSet);
PSTH_sim = nan(nCh_Total, nBins, nAmp, nSet);

nTrialsRaw_seq     = zeros(nAmp, nSet);
nTrialsClean_seq   = zeros(nAmp, nSet);
nTrialsRemoved_seq = zeros(nAmp, nSet);

nTrialsRaw_sim     = zeros(nAmp, nSet);
nTrialsClean_sim   = zeros(nAmp, nSet);
nTrialsRemoved_sim = zeros(nAmp, nSet);

SummaryRows = {};
rowCounter = 0;

%% =============================================================
%  COMPUTE SPIKE COUNTS
%% =============================================================

fprintf('\nComputing spike counts...\n');

for ss = 1:nSet

    setID = SetIDs(ss);
    current_channels = resp_channels_per_set{ss};

    setLabel = get_set_label(CHN, setID);

    if isempty(current_channels)
        fprintf('  Set %d (%s): no responding channels. Skipping.\n', setID, setLabel);
        continue;
    end

    for ai = 1:nAmp

        ampVal = Amps(ai);

        %% =====================================================
        %  Sequential target condition
        %% =====================================================

        tr_seq_raw = find(conditionType_trial == 1 & ...
                          conditionSetID_trial == setID & ...
                          prefixLength_trial == TargetSeqPrefix & ...
                          abs(isi_ms_trial - TargetSeqISI_ms) < 1e-6 & ...
                          abs(trialAmps - ampVal) < 1e-6);

        tr_seq_raw = tr_seq_raw(tr_seq_raw <= nTrials_use);

        nTrialsRaw_seq(ai, ss) = numel(tr_seq_raw);

        %% =====================================================
        %  Simultaneous target condition
        %% =====================================================

        tr_sim_raw = find(conditionType_trial == 2 & ...
                          conditionSetID_trial == setID & ...
                          activeCount_trial == TargetSimActiveCount & ...
                          abs(trialAmps - ampVal) < 1e-6);

        tr_sim_raw = tr_sim_raw(tr_sim_raw <= nTrials_use);

        nTrialsRaw_sim(ai, ss) = numel(tr_sim_raw);

        %% =====================================================
        %  Channel loop
        %% =====================================================

        seq_channel_values = [];
        sim_channel_values = [];

        for cc = 1:numel(current_channels)

            ch_idx = current_channels(cc);   % depth index
            recCh  = d(ch_idx);              % spike channel index

            if recCh < 1 || recCh > nCh || isempty(sp{recCh})
                continue;
            end

            S_ch = sp{recCh};

            %% ----- Bad trials for this channel -----
            bad_trs = [];

            if remove_bad_trials && ~isempty(BadTrials) && ch_idx <= numel(BadTrials)
                bad_trs = BadTrials{ch_idx};
            end

            %% ----- Sequential clean trials -----
            tr_seq_clean = setdiff(tr_seq_raw(:), bad_trs(:), 'stable');

            %% ----- Sim clean trials -----
            tr_sim_clean = setdiff(tr_sim_raw(:), bad_trs(:), 'stable');

            %% ----- Sequential spike count -----
            if ~isempty(tr_seq_clean)

                [count_val, psth_curve] = get_spike_count( ...
                    tr_seq_clean, trig, S_ch, post_win_ms, ...
                    edges_psth, g_sym, bin_s, FS);

                SC_seq(ch_idx, ai, ss) = count_val;
                PSTH_seq(ch_idx, :, ai, ss) = psth_curve;

                seq_channel_values = [seq_channel_values; count_val]; %#ok<AGROW>
            end

            %% ----- Sim spike count -----
            if ~isempty(tr_sim_clean)

                [count_val, psth_curve] = get_spike_count( ...
                    tr_sim_clean, trig, S_ch, post_win_ms, ...
                    edges_psth, g_sym, bin_s, FS);

                SC_sim(ch_idx, ai, ss) = count_val;
                PSTH_sim(ch_idx, :, ai, ss) = psth_curve;

                sim_channel_values = [sim_channel_values; count_val]; %#ok<AGROW>
            end
        end

        %% ----- Trial counts after bad-trial removal -----
        tr_seq_clean_for_count = remove_bad_trials_global_for_count(tr_seq_raw, BadTrials);
        tr_sim_clean_for_count = remove_bad_trials_global_for_count(tr_sim_raw, BadTrials);

        nTrialsClean_seq(ai, ss) = numel(tr_seq_clean_for_count);
        nTrialsRemoved_seq(ai, ss) = nTrialsRaw_seq(ai, ss) - nTrialsClean_seq(ai, ss);

        nTrialsClean_sim(ai, ss) = numel(tr_sim_clean_for_count);
        nTrialsRemoved_sim(ai, ss) = nTrialsRaw_sim(ai, ss) - nTrialsClean_sim(ai, ss);

        %% ----- Summary rows -----
        if ~isempty(seq_channel_values)

            rowCounter = rowCounter + 1;

            SummaryRows(rowCounter,:) = { ...
                setID, ...
                setLabel, ...
                'Sequential', ...
                ampVal, ...
                TargetSeqPrefix, ...
                NaN, ...
                TargetSeqISI_ms, ...
                numel(current_channels), ...
                nTrialsRaw_seq(ai, ss), ...
                nTrialsClean_seq(ai, ss), ...
                nTrialsRemoved_seq(ai, ss), ...
                mean(seq_channel_values, 'omitnan'), ...
                std(seq_channel_values, 0, 'omitnan') ./ sqrt(sum(~isnan(seq_channel_values)))};
        end

        if ~isempty(sim_channel_values)

            rowCounter = rowCounter + 1;

            SummaryRows(rowCounter,:) = { ...
                setID, ...
                setLabel, ...
                'Simultaneous', ...
                ampVal, ...
                NaN, ...
                TargetSimActiveCount, ...
                NaN, ...
                numel(current_channels), ...
                nTrialsRaw_sim(ai, ss), ...
                nTrialsClean_sim(ai, ss), ...
                nTrialsRemoved_sim(ai, ss), ...
                mean(sim_channel_values, 'omitnan'), ...
                std(sim_channel_values, 0, 'omitnan') ./ sqrt(sum(~isnan(sim_channel_values)))};
        end
    end
end

%% =============================================================
%  CREATE SUMMARY TABLE
%% =============================================================

if isempty(SummaryRows)

    SummaryTable = table();

else

    SummaryTable = cell2table(SummaryRows, ...
        'VariableNames', { ...
        'SetID', ...
        'SetLabel', ...
        'Condition', ...
        'Amp_uA', ...
        'PrefixLength', ...
        'ActiveCount', ...
        'ISI_ms', ...
        'NChannels', ...
        'NTrialsRaw', ...
        'NTrialsClean', ...
        'NTrialsRemoved', ...
        'MeanSpikeCount', ...
        'SEMSpikeCount'});
end

% fprintf('\nSummary table:\n');
% disp(SummaryTable);

%% =============================================================
%  PLOT: ONE LINE PER SEQUENTIAL SET AND SIM SET
%% =============================================================

figure('Color','w', 'Position', fig_position);
hold on;

seq_colors = lines(max(nSet, 1));
sim_colors = gray(max(nSet + 2, 3));

%% ----- Plot sequential lines -----
for ss = 1:nSet

    setID = SetIDs(ss);
    setLabel = get_set_label(CHN, setID);

    data_set = squeeze(SC_seq(:, :, ss));

    if all(isnan(data_set), 'all')
        continue;
    end

    AvgSeq = mean(data_set, 1, 'omitnan');
    SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));

    col = seq_colors(ss, :);

    %% Scatter channel values.
    for ai = 1:nAmp
        valid = ~isnan(data_set(:, ai));

        if any(valid)
            x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(ai);

            scatter(x_jit, data_set(valid, ai), dot_size, col, ...
                'filled', ...
                'MarkerFaceAlpha', scatter_alpha, ...
                'MarkerEdgeAlpha', scatter_alpha, ...
                'HandleVisibility','off');
        end
    end

    plot_shaded_error(Amps, AvgSeq, SEMSeq, col);

    plot(Amps, AvgSeq, '-s', ...
        'Color', col, ...
        'LineWidth', 2, ...
        'MarkerFaceColor', 'w', ...
        'DisplayName', sprintf('Seq Set %d (%s)', setID, setLabel));
end

%% ----- Plot simultaneous lines -----
for ss = 1:nSet

    setID = SetIDs(ss);
    setLabel = get_set_label(CHN, setID);

    data_set = squeeze(SC_sim(:, :, ss));

    if all(isnan(data_set), 'all')
        continue;
    end

    AvgSim = mean(data_set, 1, 'omitnan');
    SEMSim = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));

    col = sim_colors(min(ss + 1, size(sim_colors,1)), :);

    %% Scatter channel values.
    for ai = 1:nAmp
        valid = ~isnan(data_set(:, ai));

        if any(valid)
            x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(ai);

            scatter(x_jit, data_set(valid, ai), dot_size, col, ...
                'filled', ...
                'MarkerFaceAlpha', scatter_alpha, ...
                'MarkerEdgeAlpha', scatter_alpha, ...
                'HandleVisibility','off');
        end
    end

    plot_shaded_error(Amps, AvgSim, SEMSim, col);

    plot(Amps, AvgSim, '-o', ...
        'Color', col, ...
        'LineWidth', 2.2, ...
        'MarkerFaceColor', col, ...
        'DisplayName', sprintf('Sim Set %d (%s)', setID, setLabel));
end

xlabel('Amplitude (\muA)', 'FontWeight','bold');
ylabel(sprintf('Mean spike count / trial (%d-%d ms)', post_win_ms(1), post_win_ms(2)), ...
    'FontWeight','bold');

title(sprintf('Spike Count: Seq Prefix %d, ISI %.1f ms vs Sim ActiveCount %d', ...
    TargetSeqPrefix, TargetSeqISI_ms, TargetSimActiveCount), ...
    'FontWeight','bold');

legend('Location','best', 'Box','off');
box off;
%% =============================================================
%  EXPLORATORY STATISTICS
%% =============================================================

Stats = struct();

%% ---------- A. Paired-set exploratory comparison ----------
% Only sets with both seq and sim are included.

paired_y = [];
paired_group = {};
paired_amp = [];
paired_set = [];

for ss = 1:nSet

    has_seq = ~all(isnan(SC_seq(:, :, ss)), 'all');
    has_sim = ~all(isnan(SC_sim(:, :, ss)), 'all');

    if ~(has_seq && has_sim)
        continue;
    end

    for ai = 1:nAmp

        d_seq = squeeze(SC_seq(:, ai, ss));
        d_sim = squeeze(SC_sim(:, ai, ss));

        d_seq = d_seq(~isnan(d_seq));
        d_sim = d_sim(~isnan(d_sim));

        if ~isempty(d_seq)
            paired_y = [paired_y; d_seq]; %#ok<AGROW>
            paired_group = [paired_group; repmat({'Sequential'}, numel(d_seq), 1)]; %#ok<AGROW>
            paired_amp = [paired_amp; repmat(Amps(ai), numel(d_seq), 1)]; %#ok<AGROW>
            paired_set = [paired_set; repmat(SetIDs(ss), numel(d_seq), 1)]; %#ok<AGROW>
        end

        if ~isempty(d_sim)
            paired_y = [paired_y; d_sim]; %#ok<AGROW>
            paired_group = [paired_group; repmat({'Simultaneous'}, numel(d_sim), 1)]; %#ok<AGROW>
            paired_amp = [paired_amp; repmat(Amps(ai), numel(d_sim), 1)]; %#ok<AGROW>
            paired_set = [paired_set; repmat(SetIDs(ss), numel(d_sim), 1)]; %#ok<AGROW>
        end
    end
end

if numel(unique(paired_group)) == 2 && ~isempty(paired_y)

    fprintf('\n=== Exploratory ANOVA: Paired available sets only ===\n');

    [p_paired, tbl_paired, stats_paired] = anovan( ...
        paired_y, ...
        {paired_group, paired_amp}, ...
        'model', 'interaction', ...
        'varnames', {'StimType', 'Amplitude'}, ...
        'display', 'on');

    Stats.PairedAvailableSets.p = p_paired;
    Stats.PairedAvailableSets.table = tbl_paired;
    Stats.PairedAvailableSets.stats = stats_paired;

else

    fprintf('\nNo paired seq/sim sets available for paired-set ANOVA.\n');
    Stats.PairedAvailableSets = [];
end

%% ---------- B. Common-reference exploratory comparison ----------
% Uses all sequential sets and any available simultaneous set.
% This is exploratory if simultaneous exists only for one set.

common_y = [];
common_group = {};
common_amp = [];
common_set = [];

for ss = 1:nSet

    for ai = 1:nAmp

        d_seq = squeeze(SC_seq(:, ai, ss));
        d_seq = d_seq(~isnan(d_seq));

        if ~isempty(d_seq)
            common_y = [common_y; d_seq]; %#ok<AGROW>
            common_group = [common_group; repmat({'Sequential'}, numel(d_seq), 1)]; %#ok<AGROW>
            common_amp = [common_amp; repmat(Amps(ai), numel(d_seq), 1)]; %#ok<AGROW>
            common_set = [common_set; repmat(SetIDs(ss), numel(d_seq), 1)]; %#ok<AGROW>
        end

        d_sim = squeeze(SC_sim(:, ai, ss));
        d_sim = d_sim(~isnan(d_sim));

        if ~isempty(d_sim)
            common_y = [common_y; d_sim]; %#ok<AGROW>
            common_group = [common_group; repmat({'Simultaneous'}, numel(d_sim), 1)]; %#ok<AGROW>
            common_amp = [common_amp; repmat(Amps(ai), numel(d_sim), 1)]; %#ok<AGROW>
            common_set = [common_set; repmat(SetIDs(ss), numel(d_sim), 1)]; %#ok<AGROW>
        end
    end
end

sim_sets_available = [];

for ss = 1:nSet
    if ~all(isnan(SC_sim(:, :, ss)), 'all')
        sim_sets_available = [sim_sets_available, SetIDs(ss)]; %#ok<AGROW>
    end
end

if numel(sim_sets_available) == 1
    fprintf('\nWARNING: Only Set %d has simultaneous data.\n', sim_sets_available);
    fprintf('Common-reference comparison is exploratory and may include set/location effects.\n');
end

if numel(unique(common_group)) == 2 && ~isempty(common_y)

    fprintf('\n=== Exploratory ANOVA: All Seq Sets vs Available Sim Reference ===\n');

    [p_common, tbl_common, stats_common] = anovan( ...
        common_y, ...
        {common_group, common_amp}, ...
        'model', 'interaction', ...
        'varnames', {'StimType', 'Amplitude'}, ...
        'display', 'on');

    Stats.CommonReference.p = p_common;
    Stats.CommonReference.table = tbl_common;
    Stats.CommonReference.stats = stats_common;
    Stats.CommonReference.sim_sets_available = sim_sets_available;

else

    fprintf('\nNot enough seq/sim data for common-reference ANOVA.\n');
    Stats.CommonReference = [];
end

%% =============================================================
%  SAVE RESULTS
%% =============================================================

Result = struct();

Result.Info = struct();
Result.Info.data_folder = data_folder;
Result.Info.base_name = base_name;
Result.Info.created = datestr(now);
Result.Info.metric = 'Mean spike count per trial per responding channel';
Result.Info.post_win_ms = post_win_ms;
Result.Info.psth_win_ms = psth_win_ms;
Result.Info.TargetSeqPrefix = TargetSeqPrefix;
Result.Info.TargetSeqISI_ms = TargetSeqISI_ms;
Result.Info.TargetSimActiveCount = TargetSimActiveCount;
Result.Info.RespondingFileMode = RespondingFileMode;
Result.Info.RespondingFile = resp_file;
Result.Info.BadTrialFile = BadTrialFile;
Result.Info.remove_bad_trials = remove_bad_trials;
Result.Info.Amps = Amps;
Result.Info.SetIDs = SetIDs;
Result.Info.spike_file_used = spike_file_used;
Result.Info.spike_variable_used = spike_variable_used;
Result.Info.FS = FS;

Result.Data = struct();
Result.Data.SC_seq = SC_seq;
Result.Data.SC_sim = SC_sim;
Result.Data.PSTH_seq = PSTH_seq;
Result.Data.PSTH_sim = PSTH_sim;
Result.Data.PSTH_time_ms = ctrs_psth;
Result.Data.RespChannelsPerSet = resp_channels_per_set;
Result.Data.nTrialsRaw_seq = nTrialsRaw_seq;
Result.Data.nTrialsClean_seq = nTrialsClean_seq;
Result.Data.nTrialsRemoved_seq = nTrialsRemoved_seq;
Result.Data.nTrialsRaw_sim = nTrialsRaw_sim;
Result.Data.nTrialsClean_sim = nTrialsClean_sim;
Result.Data.nTrialsRemoved_sim = nTrialsRemoved_sim;

Result.Summary = struct();
Result.Summary.Table = SummaryTable;

Result.Stats = Stats;

out_filename = fullfile(save_dir, ...
    sprintf('Result_SpikeCount_MixedPrefix_%s.mat', base_name));

save(out_filename, 'Result', '-v7.3');

fprintf('\n>>> Results saved to:\n%s\n', out_filename);

%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function [sp, spike_file_used, spike_variable_used] = load_spike_file(data_folder)

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
        error('No spike file found in %s.', data_folder);
    end
end

function si = find_set_index(Responding, target_set)

    si = 0;

    if ~isfield(Responding, 'set') || isempty(Responding.set)
        return;
    end

    for k = 1:numel(Responding.set)

        if isfield(Responding.set(k), 'setID')
            if Responding.set(k).setID == target_set
                si = k;
                return;
            end
        else
            if k == target_set
                si = k;
                return;
            end
        end
    end
end

function tf = is_target_prefix(P, TargetSeqPrefix, TargetSeqISI_ms)

    tf = false;

    %% ----- Prefix check -----
    if isfield(P, 'prefix_length')
        prefix_val = P.prefix_length;
    elseif isfield(P, 'prefixLength')
        prefix_val = P.prefixLength;
    else
        prefix_val = NaN;
    end

    if isnan(prefix_val) || abs(prefix_val - TargetSeqPrefix) > 1e-6
        return;
    end

    %% ----- ISI check -----
    if isfield(P, 'isi_ms')
        isi_val = P.isi_ms;
    elseif isfield(P, 'ISI_ms')
        isi_val = P.ISI_ms;
    else
        % If ISI is not stored in responding structure, accept it.
        % The trial matching later still uses StimParams metadata.
        isi_val = TargetSeqISI_ms;
    end

    if abs(isi_val - TargetSeqISI_ms) > 1e-6
        return;
    end

    tf = true;
end

function label = get_set_label(CHN, setID)

    if ~isempty(CHN) && setID <= size(CHN,1)

        stimVec = CHN(setID,:);
        stimVec = stimVec(stimVec > 0);

        if isempty(stimVec)
            label = sprintf('Set%d', setID);
        else
            label = strjoin(arrayfun(@(x) sprintf('Ch%d', x), ...
                stimVec, 'UniformOutput', false), '→');
        end
    else
        label = sprintf('Set%d', setID);
    end
end

function tr_clean = remove_bad_trials_global_for_count(tr_raw, BadTrials)

    if isempty(BadTrials)
        tr_clean = tr_raw(:);
        return;
    end

    bad_all = [];

    for ii = 1:numel(BadTrials)
        if ~isempty(BadTrials{ii})
            bad_all = [bad_all; BadTrials{ii}(:)]; %#ok<AGROW>
        end
    end

    bad_all = unique(bad_all);

    tr_clean = setdiff(tr_raw(:), bad_all(:), 'stable');
end

function [count_val, psth_trace] = get_spike_count( ...
    tr_ids, trig, sp_data, count_win, psth_edges, g_sym, bin_s, FS)

    nTr = numel(tr_ids);
    all_psth_spikes = [];
    total_spikes_in_window = 0;

    for k = 1:nTr

        tr = tr_ids(k);
        t0 = trig(tr) / FS * 1000;

        tt = sp_data(:,1) - t0;

        %% ----- Count spikes in response window -----
        mask_count = tt >= count_win(1) & tt <= count_win(2);
        total_spikes_in_window = total_spikes_in_window + sum(mask_count);

        %% ----- Collect spikes for PSTH -----
        all_psth_spikes = [all_psth_spikes; ...
            tt(tt >= psth_edges(1) & tt <= psth_edges(end))]; %#ok<AGROW>
    end

    if nTr > 0
        count_val = total_spikes_in_window / nTr;
    else
        count_val = NaN;
    end

    h_psth = histcounts(all_psth_spikes, psth_edges);
    rate_psth = h_psth / (nTr * bin_s);
    psth_trace = conv(rate_psth, g_sym, 'same');
end

function plot_shaded_error(x, y, se, col)

    if numel(x) < 2
        return;
    end

    x = x(:)';
    y = y(:)';
    se = se(:)';

    valid = ~isnan(y);

    x = x(valid);
    y = y(valid);
    se = se(valid);

    if isempty(x)
        return;
    end

    fill([x fliplr(x)], ...
         [y + se fliplr(y - se)], ...
         col, ...
         'FaceAlpha', 0.15, ...
         'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
end