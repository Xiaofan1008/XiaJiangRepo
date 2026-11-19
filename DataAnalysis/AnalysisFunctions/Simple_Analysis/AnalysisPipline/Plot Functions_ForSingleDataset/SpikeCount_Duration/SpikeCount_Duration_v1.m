%% ============================================================
%   COMBINED SINGLE-DATASET ANALYSIS: SPIKE COUNT + DURATION
%   - Logic:
%       1. Load one dataset.
%       2. Parse amplitudes, PTDs, and stimulation sets.
%       3. Define ONE common responsive channel pool per set
%          using the existing RespondingChannels result.
%       4. For the same channel pool, calculate:
%            a) Mean Spike Count (2-20 ms)
%            b) Response Duration (3*SD threshold method)
%       5. Save channel-level raw values and summary values
%          for each Set x Amp x Condition.
%
%   - Purpose:
%       This file is intended for later population-level analysis,
%       such as:
%           * same spike count -> compare duration
%           * same duration -> compare spike count
%           * other cross-metric analyses
%
%   - Output:
%       CombinedData structure saved to .mat
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim8';
Electrode_Type = 2;   % 0: single shank rigid; 1: single shank flex; 2: four shank flex

% --- Spike Count Settings ---
post_win_ms = [2 20];     % Spike count window

% --- Duration Settings ---
baseline_win_ms = [-50 -10];
analysis_win_ms = [2 20];
psth_bin_ms     = 1;
std_threshold   = 3;      % Duration threshold = baseline mean + N*SD

% --- Condition Settings ---
Seq_PTD_Target  = 5;      % Sequential PTD target in ms

% --- Plotting (optional diagnostic plot) ---
make_plot = true;
fig_pos = [100 100 900 500];
FS = 30000;

%% =================== 1. LOAD DATA ====================
[R, sp, trig, S, QC] = load_experiment_data(data_folder);

%% =================== 2. PARSE STIM PARAMS ====================
Stim = S.StimParams;
simN = S.simultaneous_stim;
E_MAP = S.E_MAP;

if isfield(S, 'n_Trials')
    nTr = S.n_Trials;
else
    nTr = (size(Stim, 1) - 1) / simN;
end

% --- Amplitudes ---
amps_all = cell2mat(Stim(2:end,16));
trialAmps = amps_all(1:simN:end);
[Amps,~,ampIdx] = unique(trialAmps);
Amps(Amps==-1) = 0;

% --- PTDs ---
if simN > 1
    PTD_all_us = cell2mat(Stim(3:simN:end,6));
else
    PTD_all_us = zeros(nTr,1);
end
PTD_all_ms = PTD_all_us / 1000;
[PTDs_ms,~,ptdIdx] = unique(PTD_all_ms);

% --- Parse Sets ---
stimNames = Stim(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

comb = zeros(nTr, simN);
for t = 1:nTr
    rr = (t-1)*simN + (1:simN);
    v = idx_all(rr);
    v = v(v>0);
    comb(t,1:numel(v)) = v(:).';
end
[uniqueComb,~,combClass] = unique(comb, 'rows', 'stable');
nSets = size(uniqueComb,1);

% --- PTD indices ---
ptd_sim_idx = find(abs(PTDs_ms - 0) < 0.001);
ptd_seq_idx = find(abs(PTDs_ms - Seq_PTD_Target) < 0.001);

% --- Channel depth mapping ---
d = Depth_s(Electrode_Type);
nCh_Total = length(d);

%% =================== 3. PREPARE OUTPUT STRUCT ====================
CombinedData = struct();

CombinedData.Settings.data_folder      = data_folder;
CombinedData.Settings.Electrode_Type   = Electrode_Type;
CombinedData.Settings.post_win_ms      = post_win_ms;
CombinedData.Settings.baseline_win_ms  = baseline_win_ms;
CombinedData.Settings.analysis_win_ms  = analysis_win_ms;
CombinedData.Settings.psth_bin_ms      = psth_bin_ms;
CombinedData.Settings.std_threshold    = std_threshold;
CombinedData.Settings.Seq_PTD_Target   = Seq_PTD_Target;
CombinedData.Settings.FS               = FS;

CombinedData.Amps     = Amps;
CombinedData.PTDs_ms  = PTDs_ms;
CombinedData.uniqueComb = uniqueComb;

fprintf('Analyzing combined spike count + duration...\n');
fprintf('Dataset contains %d Sets and %d Amplitudes.\n', nSets, length(Amps));

% Define edges once for duration
full_edges = (baseline_win_ms(1)-10) : psth_bin_ms : (analysis_win_ms(2)+10);
full_centers = full_edges(1:end-1) + psth_bin_ms/2;

%% =================== 4. PROCESS SET BY SET ====================
for ss = 1:nSets

    % --- Set label ---
    stimCh = uniqueComb(ss,:);
    stimCh = stimCh(stimCh>0);

    set_label = ['Set ' num2str(ss) ' (Ch ' num2str(stimCh) ')'];
    CombinedData.Set(ss).Label = set_label;
    CombinedData.Set(ss).StimCh = stimCh;

    % ============================================================
    % A. DEFINE COMMON RESPONSIVE CHANNEL POOL FOR THIS SET
    %    Use union of responsive channels from RespondingChannels file
    %    across:
    %       - all amplitudes
    %       - Sim (PTD = 0)
    %       - Seq (PTD = Seq_PTD_Target)
    % ============================================================
    local_resp_mask = false(nCh_Total, 1);

    % --- Check Sim channels ---
    if ~isempty(ptd_sim_idx)
        for ai = 1:length(Amps)
            try
                this = R.set(ss).amp(ai).ptd(ptd_sim_idx).channel;
                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch) = true;
                    end
                end
            catch
            end
        end
    end

    % --- Check Seq channels ---
    if ~isempty(ptd_seq_idx)
        for ai = 1:length(Amps)
            try
                this = R.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch) = true;
                    end
                end
            catch
            end
        end
    end

    local_resp_indices = find(local_resp_mask);

    CombinedData.Set(ss).ResponsiveChannels = local_resp_indices;
    CombinedData.Set(ss).N_ResponsiveChannels = length(local_resp_indices);

    fprintf('  Set %d (Ch %s): %d responsive channels in common pool\n', ...
        ss, num2str(stimCh), length(local_resp_indices));

    if isempty(local_resp_indices)
        continue;
    end

    % ============================================================
    % B. LOOP OVER AMPLITUDES
    % ============================================================
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);

        CombinedData.Set(ss).Amp(ai).Val = curr_amp;

        % --- Preallocate per-channel results for this amplitude ---
        spike_sim_ch = nan(length(local_resp_indices), 1);
        spike_seq_ch = nan(length(local_resp_indices), 1);
        dur_sim_ch   = nan(length(local_resp_indices), 1);
        dur_seq_ch   = nan(length(local_resp_indices), 1);

        % --- Trial IDs for this set x amplitude x condition ---
        tr_sim = [];
        tr_seq = [];

        if ~isempty(ptd_sim_idx)
            tr_sim = find(ampIdx == ai & combClass == ss & ptdIdx == ptd_sim_idx);
        end

        if ~isempty(ptd_seq_idx)
            tr_seq = find(ampIdx == ai & combClass == ss & ptdIdx == ptd_seq_idx);
        end

        % ========================================================
        % C. LOOP OVER CHANNELS IN COMMON RESPONSIVE POOL
        % ========================================================
        for kk = 1:length(local_resp_indices)
            ch_idx = local_resp_indices(kk);
            recCh  = d(ch_idx);
            sp_ch  = sp{recCh};

            % --- Bad trial handling ---
            bad_trs = [];
            if ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
                bad_trs = QC.BadTrials{ch_idx};
            end

            % --- Bad channel handling ---
            if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
                continue;
            end

            % ---------------------------
            % Simultaneous condition
            % ---------------------------
            if ~isempty(tr_sim)
                tr_sim_use = setdiff(tr_sim, bad_trs);

                if ~isempty(tr_sim_use)
                    spike_sim_ch(kk) = get_spike_count(tr_sim_use, trig, sp_ch, post_win_ms, FS);

                    dur_sim_ch(kk) = calc_duration( ...
                        tr_sim_use, trig, sp_ch, ...
                        full_edges, full_centers, ...
                        baseline_win_ms, analysis_win_ms, ...
                        std_threshold, FS);
                end
            end

            % ---------------------------
            % Sequential condition
            % ---------------------------
            if ~isempty(tr_seq)
                tr_seq_use = setdiff(tr_seq, bad_trs);

                if ~isempty(tr_seq_use)
                    spike_seq_ch(kk) = get_spike_count(tr_seq_use, trig, sp_ch, post_win_ms, FS);

                    dur_seq_ch(kk) = calc_duration( ...
                        tr_seq_use, trig, sp_ch, ...
                        full_edges, full_centers, ...
                        baseline_win_ms, analysis_win_ms, ...
                        std_threshold, FS);
                end
            end
        end

        % ========================================================
        % D. SAVE CHANNEL-LEVEL RAW RESULTS
        % ========================================================
        CombinedData.Set(ss).Amp(ai).ChannelIdx = local_resp_indices(:);

        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_All = spike_sim_ch;
        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_All = spike_seq_ch;

        CombinedData.Set(ss).Amp(ai).Sim.Duration_All   = dur_sim_ch;
        CombinedData.Set(ss).Amp(ai).Seq.Duration_All   = dur_seq_ch;

        % ========================================================
        % E. SAVE SUMMARY RESULTS
        % ========================================================
        % --- Spike Count summaries ---
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_Mean   = mean(spike_sim_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_Median = median(spike_sim_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_SEM    = std(spike_sim_ch, 'omitnan') / sqrt(sum(~isnan(spike_sim_ch)));
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_N      = sum(~isnan(spike_sim_ch));

        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_Mean   = mean(spike_seq_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_Median = median(spike_seq_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_SEM    = std(spike_seq_ch, 'omitnan') / sqrt(sum(~isnan(spike_seq_ch)));
        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_N      = sum(~isnan(spike_seq_ch));

        % --- Duration summaries ---
        CombinedData.Set(ss).Amp(ai).Sim.Duration_Mean   = mean(dur_sim_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Sim.Duration_Median = median(dur_sim_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Sim.Duration_SEM    = std(dur_sim_ch, 'omitnan') / sqrt(sum(~isnan(dur_sim_ch)));
        CombinedData.Set(ss).Amp(ai).Sim.Duration_N      = sum(~isnan(dur_sim_ch));

        CombinedData.Set(ss).Amp(ai).Seq.Duration_Mean   = mean(dur_seq_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Seq.Duration_Median = median(dur_seq_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Seq.Duration_SEM    = std(dur_seq_ch, 'omitnan') / sqrt(sum(~isnan(dur_seq_ch)));
        CombinedData.Set(ss).Amp(ai).Seq.Duration_N      = sum(~isnan(dur_seq_ch));

        % ========================================================
        % F. PRINT SUMMARY TO COMMAND WINDOW
        % ========================================================
        fprintf(['    Set %d | Amp %4.1f uA | ' ...
                 'Sim Spike = %5.2f (N=%2d), Seq Spike = %5.2f (N=%2d) | ' ...
                 'Sim Dur = %5.2f ms (N=%2d), Seq Dur = %5.2f ms (N=%2d)\n'], ...
                 ss, curr_amp, ...
                 CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_Mean, CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_N, ...
                 CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_Mean, CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_N, ...
                 CombinedData.Set(ss).Amp(ai).Sim.Duration_Mean,   CombinedData.Set(ss).Amp(ai).Sim.Duration_N, ...
                 CombinedData.Set(ss).Amp(ai).Seq.Duration_Mean,   CombinedData.Set(ss).Amp(ai).Seq.Duration_N);
    end
end

%% =================== 5. OPTIONAL DIAGNOSTIC PLOT ====================
% This plot is only for quick checking.
% It plots mean Spike Count and mean Duration for one dataset.
if make_plot
    fprintf('Generating quick diagnostic plots...\n');

    figure('Color','w', 'Position', fig_pos);

    % ------------------------------------------------------------
    % Plot 1: Spike Count
    % ------------------------------------------------------------
    subplot(1,2,1); hold on;
    title('Spike Count', 'FontWeight', 'bold');

    for ss = 1:nSets
        sim_mean = nan(1, length(Amps));
        seq_mean = nan(1, length(Amps));

        for ai = 1:length(Amps)
            if length(CombinedData.Set) >= ss && isfield(CombinedData.Set(ss), 'Amp') && length(CombinedData.Set(ss).Amp) >= ai
                sim_mean(ai) = CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_Mean;
                seq_mean(ai) = CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_Mean;
            end
        end

        if any(~isnan(sim_mean))
            plot(Amps, sim_mean, '--o', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', ...
                'DisplayName', ['Sim Set ' num2str(ss)]);
        end
        if any(~isnan(seq_mean))
            plot(Amps, seq_mean, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', ...
                'DisplayName', ['Seq Set ' num2str(ss)]);
        end
    end

    xlabel('Amplitude (uA)');
    ylabel('Spike Count');
    box off;

    % ------------------------------------------------------------
    % Plot 2: Duration
    % ------------------------------------------------------------
    subplot(1,2,2); hold on;
    title('Response Duration', 'FontWeight', 'bold');

    for ss = 1:nSets
        sim_mean = nan(1, length(Amps));
        seq_mean = nan(1, length(Amps));

        for ai = 1:length(Amps)
            if length(CombinedData.Set) >= ss && isfield(CombinedData.Set(ss), 'Amp') && length(CombinedData.Set(ss).Amp) >= ai
                sim_mean(ai) = CombinedData.Set(ss).Amp(ai).Sim.Duration_Mean;
                seq_mean(ai) = CombinedData.Set(ss).Amp(ai).Seq.Duration_Mean;
            end
        end

        if any(~isnan(sim_mean))
            plot(Amps, sim_mean, '--o', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', ...
                'DisplayName', ['Sim Set ' num2str(ss)]);
        end
        if any(~isnan(seq_mean))
            plot(Amps, seq_mean, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', ...
                'DisplayName', ['Seq Set ' num2str(ss)]);
        end
    end

    xlabel('Amplitude (uA)');
    ylabel('Duration (ms)');
    box off;
end

%% =================== 6. FEASIBILITY CHECK: SPIKE COUNT vs DURATION ====================
% This section is only for checking whether a cross-metric analysis is worth doing.
% It pools all valid channel-level values within this dataset and tests:
%   1) Is duration related to spike count?
%   2) After accounting for spike count, does condition (Sim vs Seq) still matter?

make_relationship_plot = true;

% ------------------------------------------------------------
% A. Collect pooled channel-level values across all Sets x Amps
% ------------------------------------------------------------
Spike_All = [];
Duration_All = [];
Condition_All = [];   % 0 = Sim, 1 = Seq
Set_All = [];
Amp_All = [];
Channel_All = [];

for ss = 1:nSets
    if ~isfield(CombinedData.Set(ss), 'Amp'), continue; end

    for ai = 1:length(Amps)
        if length(CombinedData.Set(ss).Amp) < ai, continue; end

        % --- Sim ---
        if isfield(CombinedData.Set(ss).Amp(ai), 'Sim')
            spike_sim = CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_All;
            dur_sim   = CombinedData.Set(ss).Amp(ai).Sim.Duration_All;

            valid_sim = ~isnan(spike_sim) & ~isnan(dur_sim);
            if any(valid_sim)
                Spike_All     = [Spike_All; spike_sim(valid_sim)];
                Duration_All  = [Duration_All; dur_sim(valid_sim)];
                Condition_All = [Condition_All; zeros(sum(valid_sim),1)];
                Set_All       = [Set_All; ss * ones(sum(valid_sim),1)];
                Amp_All       = [Amp_All; Amps(ai) * ones(sum(valid_sim),1)];
                Channel_All   = [Channel_All; CombinedData.Set(ss).Amp(ai).ChannelIdx(valid_sim)];
            end
        end

        % --- Seq ---
        if isfield(CombinedData.Set(ss).Amp(ai), 'Seq')
            spike_seq = CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_All;
            dur_seq   = CombinedData.Set(ss).Amp(ai).Seq.Duration_All;

            valid_seq = ~isnan(spike_seq) & ~isnan(dur_seq);
            if any(valid_seq)
                Spike_All     = [Spike_All; spike_seq(valid_seq)];
                Duration_All  = [Duration_All; dur_seq(valid_seq)];
                Condition_All = [Condition_All; ones(sum(valid_seq),1)];
                Set_All       = [Set_All; ss * ones(sum(valid_seq),1)];
                Amp_All       = [Amp_All; Amps(ai) * ones(sum(valid_seq),1)];
                Channel_All   = [Channel_All; CombinedData.Set(ss).Amp(ai).ChannelIdx(valid_seq)];
            end
        end
    end
end

% Save pooled diagnostic data into output structure as well
CombinedData.RelationshipCheck.Spike_All = Spike_All;
CombinedData.RelationshipCheck.Duration_All = Duration_All;
CombinedData.RelationshipCheck.Condition_All = Condition_All;
CombinedData.RelationshipCheck.Set_All = Set_All;
CombinedData.RelationshipCheck.Amp_All = Amp_All;
CombinedData.RelationshipCheck.Channel_All = Channel_All;

fprintf('\n============================================================\n');
fprintf('FEASIBILITY CHECK: SPIKE COUNT vs DURATION\n');
fprintf('============================================================\n');
fprintf('Total pooled valid points: %d\n', length(Spike_All));
fprintf('Sim points: %d | Seq points: %d\n', sum(Condition_All==0), sum(Condition_All==1));

% ------------------------------------------------------------
% B. Spearman correlation (more robust than Pearson here)
% ------------------------------------------------------------
sim_mask = (Condition_All == 0);
seq_mask = (Condition_All == 1);

rho_sim = NaN; p_sim = NaN;
rho_seq = NaN; p_seq = NaN;

if sum(sim_mask) >= 3
    [rho_sim, p_sim] = corr(Spike_All(sim_mask), Duration_All(sim_mask), ...
        'Type', 'Spearman', 'Rows', 'complete');
end

if sum(seq_mask) >= 3
    [rho_seq, p_seq] = corr(Spike_All(seq_mask), Duration_All(seq_mask), ...
        'Type', 'Spearman', 'Rows', 'complete');
end

fprintf('\nSpearman correlation:\n');
fprintf('  Sim: rho = %.4f, p = %s\n', rho_sim, format_p_short(p_sim));
fprintf('  Seq: rho = %.4f, p = %s\n', rho_seq, format_p_short(p_seq));

CombinedData.RelationshipCheck.Spearman.Sim.rho = rho_sim;
CombinedData.RelationshipCheck.Spearman.Sim.p   = p_sim;
CombinedData.RelationshipCheck.Spearman.Seq.rho = rho_seq;
CombinedData.RelationshipCheck.Spearman.Seq.p   = p_seq;

% ------------------------------------------------------------
% C. Simple linear model:
%    Duration = b0 + b1*SpikeCount + b2*Condition
%    where Condition: 0 = Sim, 1 = Seq
% ------------------------------------------------------------
mdl = [];
tbl = table(Spike_All, Duration_All, Condition_All, ...
    'VariableNames', {'SpikeCount', 'Duration', 'Condition'});

if height(tbl) >= 5
    mdl = fitlm(tbl, 'Duration ~ SpikeCount + Condition');
    
    fprintf('\nLinear model: Duration ~ SpikeCount + Condition\n');
    disp(mdl);

    CombinedData.RelationshipCheck.LinearModel.Coefficients = mdl.Coefficients;
    CombinedData.RelationshipCheck.LinearModel.Rsquared     = mdl.Rsquared;
    CombinedData.RelationshipCheck.LinearModel.Formula      = 'Duration ~ SpikeCount + Condition';
else
    fprintf('\nLinear model skipped: not enough pooled points.\n');
end

% ------------------------------------------------------------
% D. Optional scatter plot for visual screening
% ------------------------------------------------------------
if make_relationship_plot && ~isempty(Spike_All)
    figure('Color','w', 'Position', [150 150 650 550]);
    hold on;

    % Raw pooled scatter
    scatter(Spike_All(sim_mask), Duration_All(sim_mask), 22, [0.7 0.7 0.7], 'o', ...
        'LineWidth', 0.8, 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Simultaneous');
    scatter(Spike_All(seq_mask), Duration_All(seq_mask), 22, [0.3 0.3 0.3], 's', ...
        'filled', 'MarkerFaceAlpha', 0.4, 'DisplayName', 'Sequential');

    % Add simple fitted lines separately for visual help only
    if sum(sim_mask) >= 2
        x_sim = Spike_All(sim_mask);
        y_sim = Duration_All(sim_mask);
        pfit_sim = polyfit(x_sim, y_sim, 1);
        xq_sim = linspace(min(x_sim), max(x_sim), 100);
        yq_sim = polyval(pfit_sim, xq_sim);
        plot(xq_sim, yq_sim, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    if sum(seq_mask) >= 2
        x_seq = Spike_All(seq_mask);
        y_seq = Duration_All(seq_mask);
        pfit_seq = polyfit(x_seq, y_seq, 1);
        xq_seq = linspace(min(x_seq), max(x_seq), 100);
        yq_seq = polyval(pfit_seq, xq_seq);
        plot(xq_seq, yq_seq, '-k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    xlabel('Spike Count', 'FontWeight', 'bold');
    ylabel('Duration (ms)', 'FontWeight', 'bold');
    title('Single-Dataset Feasibility Check: Spike Count vs Duration', 'FontWeight', 'bold');
    legend('Location', 'best', 'Box', 'off');
    box off;
end

%% =================== 7. SAVE RESULTS ====================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

parts = split(data_folder, filesep);
exp_id = parts{end};

out_filename = fullfile(save_dir, ['Result_Combined_SpikeDuration_DX013_' exp_id '.mat']);

save(out_filename, ...
    'CombinedData', 'Amps', 'PTDs_ms', ...
    'post_win_ms', 'baseline_win_ms', 'analysis_win_ms', ...
    'psth_bin_ms', 'std_threshold', 'Seq_PTD_Target');

fprintf('\n>>> Combined results saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function count_val = get_spike_count(tr_ids, trig, sp_data, count_win, FS)
    % ------------------------------------------------------------
    % Calculate mean spike count per trial in a fixed post-stim window
    % ------------------------------------------------------------
    nTr = numel(tr_ids);
    total_spikes = 0;

    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr) / FS * 1000;   % trigger time in ms
        tt = sp_data(:,1) - t0;      % spike times aligned to trigger

        mask = tt >= count_win(1) & tt <= count_win(2);
        total_spikes = total_spikes + sum(mask);
    end

    if nTr > 0
        count_val = total_spikes / nTr;
    else
        count_val = NaN;
    end
end

function duration = calc_duration(tr_ids, trig, sp_data, edges, centers, ...
                                  base_win, anal_win, nSD, FS)
    % ------------------------------------------------------------
    % Calculate response duration using PSTH threshold crossing
    %
    % Steps:
    %   1. Collect all spikes across valid trials
    %   2. Build PSTH
    %   3. Smooth PSTH
    %   4. Compute baseline mean and SD
    %   5. Threshold = mean + nSD*SD (minimum 10 Hz)
    %   6. Duration = total time above threshold in analysis window
    % ------------------------------------------------------------
    if isempty(tr_ids)
        duration = NaN;
        return;
    end

    all_spikes = [];
    for k = 1:length(tr_ids)
        tr = tr_ids(k);
        t0 = trig(tr) / FS * 1000;

        win_pad = [edges(1), edges(end)];
        tt = sp_data(:,1) - t0;
        mask = tt >= win_pad(1) & tt <= win_pad(2);

        all_spikes = [all_spikes; tt(mask)]; %#ok<AGROW>
    end

    if isempty(all_spikes)
        duration = NaN;
        return;
    end

    bin_size_ms = edges(2) - edges(1);

    counts = histcounts(all_spikes, edges);
    rate_hz = (counts / length(tr_ids)) * (1000 / bin_size_ms);

    % Smooth PSTH
    smooth_rate = smoothdata(rate_hz, 'gaussian', 10);

    % Baseline threshold
    base_mask = centers >= base_win(1) & centers <= base_win(2);
    base_data = smooth_rate(base_mask);

    mu = mean(base_data);
    sigma = std(base_data);
    if sigma == 0
        sigma = 1;
    end

    threshold = max(mu + nSD * sigma, 10);

    % Duration inside analysis window
    anal_mask = centers >= anal_win(1) & centers <= anal_win(2);
    anal_rate = smooth_rate(anal_mask);

    if ~any(anal_rate > threshold)
        duration = NaN;
    else
        valid_mask = anal_rate > threshold;
        duration = sum(valid_mask) * bin_size_ms;
    end
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    % ------------------------------------------------------------
    % Load all required files from one dataset folder
    % ------------------------------------------------------------
    cd(folder);

    f = dir('*_RespondingChannels.mat');
    if isempty(f)
        error('No RespondingChannels file in %s', folder);
    end
    R = load(f(1).name).Responding;

    f = dir('*sp_xia_SSD.mat');
    if isempty(f)
        f = dir('*sp_xia.mat');
    end
    if isempty(f)
        error('No spike file in %s', folder);
    end

    S_sp = load(f(1).name);
    if isfield(S_sp, 'sp_corr')
        sp = S_sp.sp_corr;
    elseif isfield(S_sp, 'sp_SSD')
        sp = S_sp.sp_SSD;
    else
        sp = S_sp.sp_in;
    end

    if isempty(dir('*.trig.dat'))
        cleanTrig_sabquick;
    end
    trig = loadTrig(0);

    f_exp = dir('*_exp_datafile_*.mat');
    if isempty(f_exp)
        error('No exp_datafile found in %s', folder);
    end
    S = load(f_exp(1).name);

    QC.BadCh = [];
    QC.BadTrials = [];

    f_bc = dir('*.BadChannels.mat');
    if ~isempty(f_bc)
        tmp = load(f_bc(1).name);
        QC.BadCh = tmp.BadCh_perSet;
    end

    f_bt = dir('*.BadTrials.mat');
    if ~isempty(f_bt)
        tmp = load(f_bt(1).name);
        QC.BadTrials = tmp.BadTrials;
    end
end

function p_str = format_p_short(p)
    if isnan(p)
        p_str = 'NaN';
    elseif p < 1e-4
        p_str = sprintf('%.2e', p);
    else
        p_str = sprintf('%.4f', p);
    end
end