%% ============================================================
%   COMBINED SINGLE-DATASET ANALYSIS: SPIKE COUNT + DURATION
%   (SEPARATE SIM / SEQ FOLDERS VERSION + FEASIBILITY CHECK)
%
%   - Logic:
%       1. Load Sim and Seq folders separately.
%       2. Parse amplitudes, PTDs, and stimulation sets separately.
%       3. Check that Sim and Seq share the same amplitude axis and set structure.
%       4. Define ONE common responsive channel pool per set
%          using the union of responsive channels from Rsim and Rseq.
%       5. For the same channel pool, calculate:
%            a) Mean Spike Count (2-20 ms)
%            b) Response Duration (3*SD threshold method)
%       6. Save channel-level raw values and summary values
%          for each Set x Amp x Condition.
%       7. NEW: Feasibility check for Spike Count vs Duration relationship.
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
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Sim';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Seq';
Electrode_Type = 2;   % 0: single shank rigid; 1: single shank flex; 2: four shank flex

% --- Spike Count Settings ---
post_win_ms = [2 20];     % Spike count window

% --- Duration Settings ---
baseline_win_ms = [-50 -10];
analysis_win_ms = [2 20];
psth_bin_ms     = 1;
std_threshold   = 3;      % Duration threshold = baseline mean + N*SD

% --- Condition Settings ---
Seq_PTD_Target  = 5.5;    % Sequential PTD target in ms

% --- Plotting (optional diagnostic plot) ---
make_plot = true;
make_feasibility_plot = true;
fig_pos = [100 100 900 500];
FS = 30000;

%% =================== 1. LOAD DATA ====================
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

%% =================== 2. PARSE STIM PARAMS ====================
% ------------------------------------------------------------
% Parse SIM folder
% ------------------------------------------------------------
Stim_sim = Ssim.StimParams;
simN_sim = Ssim.simultaneous_stim;
E_MAP_sim = Ssim.E_MAP;

if isfield(Ssim, 'n_Trials')
    nTr_sim = Ssim.n_Trials;
else
    nTr_sim = (size(Stim_sim, 1) - 1) / simN_sim;
end

amps_all_sim = cell2mat(Stim_sim(2:end,16));
trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim);
Amps_sim(Amps_sim == -1) = 0;

stimNames_sim = Stim_sim(2:end,1);
[~, idx_all_sim] = ismember(stimNames_sim, E_MAP_sim(2:end));

comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim
    rr = (t-1)*simN_sim + (1:simN_sim);
    v = idx_all_sim(rr);
    v = v(v > 0);
    comb_sim(t,1:numel(v)) = v(:).';
end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim, 'rows', 'stable');
nSets_sim = size(uniqueComb_sim,1);

% ------------------------------------------------------------
% Parse SEQ folder
% ------------------------------------------------------------
Stim_seq = Sseq.StimParams;
simN_seq = Sseq.simultaneous_stim;
E_MAP_seq = Sseq.E_MAP;

if isfield(Sseq, 'n_Trials')
    nTr_seq = Sseq.n_Trials;
else
    nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq;
end

amps_all_seq = cell2mat(Stim_seq(2:end,16));
trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq);
Amps_seq(Amps_seq == -1) = 0;

PTD_all_us_seq = cell2mat(Stim_seq(3:simN_seq:end,6));
PTD_all_ms_seq = PTD_all_us_seq / 1000;
[PTDs_ms_seq,~,ptdIdx_seq] = unique(PTD_all_ms_seq);

stimNames_seq = Stim_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP_seq(2:end));

comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq
    rr = (t-1)*simN_seq + (1:simN_seq);
    v = idx_all_seq(rr);
    v = v(v > 0);
    comb_seq(t,1:numel(v)) = v(:).';
end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq, 'rows', 'stable');
nSets_seq = size(uniqueComb_seq,1);

% ------------------------------------------------------------
% Safety checks
% ------------------------------------------------------------
% if length(Amps_sim) ~= length(Amps_seq) || any(abs(Amps_sim - Amps_seq) > 0.001)
%     error('Sim and Seq amplitude axes do not match.');
% end
% 
% if nSets_sim ~= nSets_seq || ~isequal(uniqueComb_sim, uniqueComb_seq)
%     error('Sim and Seq stimulation sets do not match.');
% end

% Use shared axis / structure
Amps = Amps_sim;
nSets = nSets_sim;
uniqueComb = uniqueComb_sim;

% PTD index for sequential target
ptd_seq_idx = find(abs(PTDs_ms_seq - Seq_PTD_Target) < 0.001);
if isempty(ptd_seq_idx)
    error('Could not find Seq PTD target %.3f ms in sequential folder.', Seq_PTD_Target);
end

% Channel mapping
d = Depth_s(Electrode_Type);
nCh_Total = length(d);

%% =================== 3. PREPARE OUTPUT STRUCT ====================
CombinedData = struct();

CombinedData.Settings.folder_sim         = folder_sim;
CombinedData.Settings.folder_seq         = folder_seq;
CombinedData.Settings.Electrode_Type     = Electrode_Type;
CombinedData.Settings.post_win_ms        = post_win_ms;
CombinedData.Settings.baseline_win_ms    = baseline_win_ms;
CombinedData.Settings.analysis_win_ms    = analysis_win_ms;
CombinedData.Settings.psth_bin_ms        = psth_bin_ms;
CombinedData.Settings.std_threshold      = std_threshold;
CombinedData.Settings.Seq_PTD_Target     = Seq_PTD_Target;
CombinedData.Settings.FS                 = FS;

CombinedData.Amps         = Amps;
CombinedData.PTDs_ms_Seq  = PTDs_ms_seq;
CombinedData.uniqueComb   = uniqueComb;

fprintf('Analyzing combined spike count + duration (separate folders)...\n');
fprintf('Dataset contains %d Sets and %d Amplitudes.\n', nSets, length(Amps));

% Define edges once for duration
full_edges = (baseline_win_ms(1)-10) : psth_bin_ms : (analysis_win_ms(2)+10);
full_centers = full_edges(1:end-1) + psth_bin_ms/2;

%% =================== 4. PROCESS SET BY SET ====================
for ss = 1:nSets

    % --- Set label ---
    stimCh = uniqueComb(ss,:);
    stimCh = stimCh(stimCh > 0);

    set_label = ['Set ' num2str(ss) ' (Ch ' num2str(stimCh) ')'];
    CombinedData.Set(ss).Label = set_label;
    CombinedData.Set(ss).StimCh = stimCh;

    % ============================================================
    % A. DEFINE COMMON RESPONSIVE CHANNEL POOL FOR THIS SET
    %    Use union of responsive channels from:
    %       - Rsim
    %       - Rseq
    % ============================================================
    local_resp_mask = false(nCh_Total, 1);

    % --- Check Sim responsive channels ---
    try
        for ai = 1:length(Amps)
            this = Rsim.set(ss).amp(ai).ptd(1).channel;   % Sim folder usually PTD=0 at index 1
            for ch = 1:min(length(this), nCh_Total)
                if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                    local_resp_mask(ch) = true;
                end
            end
        end
    catch
    end

    % --- Check Seq responsive channels ---
    try
        for ai = 1:length(Amps)
            this = Rseq.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
            for ch = 1:min(length(this), nCh_Total)
                if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                    local_resp_mask(ch) = true;
                end
            end
        end
    catch
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
        tr_sim = find(combClass_sim == ss & ampIdx_sim == ai);
        tr_seq = find(combClass_seq == ss & ampIdx_seq == ai & ptdIdx_seq == ptd_seq_idx);

        % ========================================================
        % C. LOOP OVER CHANNELS IN COMMON RESPONSIVE POOL
        % ========================================================
        for kk = 1:length(local_resp_indices)
            ch_idx = local_resp_indices(kk);
            recCh  = d(ch_idx);

            % --- Bad channel handling ---
            is_bad = false;

            if ~isempty(QC_Sim.BadCh)
                if iscell(QC_Sim.BadCh)
                    if ss <= length(QC_Sim.BadCh) && ismember(ch_idx, QC_Sim.BadCh{ss})
                        is_bad = true;
                    end
                elseif ismember(ch_idx, QC_Sim.BadCh)
                    is_bad = true;
                end
            end

            if ~isempty(QC_Seq.BadCh)
                if iscell(QC_Seq.BadCh)
                    if ss <= length(QC_Seq.BadCh) && ismember(ch_idx, QC_Seq.BadCh{ss})
                        is_bad = true;
                    end
                elseif ismember(ch_idx, QC_Seq.BadCh)
                    is_bad = true;
                end
            end

            if is_bad
                continue;
            end

            % --- Bad trial handling ---
            bt_sim = [];
            if ~isempty(QC_Sim.BadTrials) && ch_idx <= length(QC_Sim.BadTrials)
                bt_sim = QC_Sim.BadTrials{ch_idx};
            end

            bt_seq = [];
            if ~isempty(QC_Seq.BadTrials) && ch_idx <= length(QC_Seq.BadTrials)
                bt_seq = QC_Seq.BadTrials{ch_idx};
            end

            sp_ch_sim = sp_sim{recCh};
            sp_ch_seq = sp_seq{recCh};

            % ---------------------------
            % Simultaneous condition
            % ---------------------------
            tr_sim_use = setdiff(tr_sim, bt_sim);
            if ~isempty(tr_sim_use)
                spike_sim_ch(kk) = get_spike_count(tr_sim_use, trig_sim, sp_ch_sim, post_win_ms, FS);

                dur_sim_ch(kk) = calc_duration( ...
                    tr_sim_use, trig_sim, sp_ch_sim, ...
                    full_edges, full_centers, ...
                    baseline_win_ms, analysis_win_ms, ...
                    std_threshold, FS);
            end

            % ---------------------------
            % Sequential condition
            % ---------------------------
            tr_seq_use = setdiff(tr_seq, bt_seq);
            if ~isempty(tr_seq_use)
                spike_seq_ch(kk) = get_spike_count(tr_seq_use, trig_seq, sp_ch_seq, post_win_ms, FS);

                dur_seq_ch(kk) = calc_duration( ...
                    tr_seq_use, trig_seq, sp_ch_seq, ...
                    full_edges, full_centers, ...
                    baseline_win_ms, analysis_win_ms, ...
                    std_threshold, FS);
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
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_Mean   = mean(spike_sim_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_Median = median(spike_sim_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_SEM    = std(spike_sim_ch, 'omitnan') / sqrt(sum(~isnan(spike_sim_ch)));
        CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_N      = sum(~isnan(spike_sim_ch));

        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_Mean   = mean(spike_seq_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_Median = median(spike_seq_ch, 'omitnan');
        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_SEM    = std(spike_seq_ch, 'omitnan') / sqrt(sum(~isnan(spike_seq_ch)));
        CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_N      = sum(~isnan(spike_seq_ch));

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

%% =================== 6. FEASIBILITY CHECK: SPIKE COUNT VS DURATION ====================
fprintf('\n====================================================================\n');
fprintf('FEASIBILITY CHECK: SPIKE COUNT VS DURATION\n');
fprintf('====================================================================\n');

all_spike_sim = [];
all_dur_sim   = [];

all_spike_seq = [];
all_dur_seq   = [];

for ss = 1:length(CombinedData.Set)
    if ~isfield(CombinedData.Set(ss), 'Amp')
        continue;
    end

    for ai = 1:length(CombinedData.Set(ss).Amp)
        if ~isfield(CombinedData.Set(ss).Amp(ai), 'Sim') || ~isfield(CombinedData.Set(ss).Amp(ai), 'Seq')
            continue;
        end

        sim_sp = CombinedData.Set(ss).Amp(ai).Sim.SpikeCount_All;
        sim_du = CombinedData.Set(ss).Amp(ai).Sim.Duration_All;

        seq_sp = CombinedData.Set(ss).Amp(ai).Seq.SpikeCount_All;
        seq_du = CombinedData.Set(ss).Amp(ai).Seq.Duration_All;

        valid_sim = ~isnan(sim_sp) & ~isnan(sim_du);
        valid_seq = ~isnan(seq_sp) & ~isnan(seq_du);

        all_spike_sim = [all_spike_sim; sim_sp(valid_sim)]; %#ok<AGROW>
        all_dur_sim   = [all_dur_sim;   sim_du(valid_sim)]; %#ok<AGROW>

        all_spike_seq = [all_spike_seq; seq_sp(valid_seq)]; %#ok<AGROW>
        all_dur_seq   = [all_dur_seq;   seq_du(valid_seq)]; %#ok<AGROW>
    end
end

fprintf('Total pooled valid points: %d\n', length(all_spike_sim) + length(all_spike_seq));
fprintf('Sim points: %d | Seq points: %d\n', length(all_spike_sim), length(all_spike_seq));

if length(all_spike_sim) >= 3 && length(all_spike_seq) >= 3

    % ------------------------------------------------------------
    % 1. Spearman correlation (non-parametric trend check)
    % ------------------------------------------------------------
    [rho_sim, p_sim_corr] = corr(all_spike_sim, all_dur_sim, 'Type', 'Spearman', 'Rows', 'complete');
    [rho_seq, p_seq_corr] = corr(all_spike_seq, all_dur_seq, 'Type', 'Spearman', 'Rows', 'complete');

    fprintf('\nSpearman correlation:\n');
    fprintf('  Sim: rho = %.4f, p = %.4e\n', rho_sim, p_sim_corr);
    fprintf('  Seq: rho = %.4f, p = %.4e\n', rho_seq, p_seq_corr);

    % ------------------------------------------------------------
    % 2. Linear model:
    %    Duration ~ SpikeCount + Condition
    %    Condition: 0 = Sim, 1 = Seq
    % ------------------------------------------------------------
    spike_all = [all_spike_sim; all_spike_seq];
    dur_all   = [all_dur_sim;   all_dur_seq];
    cond_all  = [zeros(length(all_spike_sim),1); ones(length(all_spike_seq),1)];

    T = table(spike_all, dur_all, cond_all, ...
        'VariableNames', {'SpikeCount', 'Duration', 'Condition'});

    mdl = fitlm(T, 'Duration ~ SpikeCount + Condition');

    fprintf('\nLinear model: Duration ~ SpikeCount + Condition\n');
    disp(mdl);

    % ------------------------------------------------------------
    % 3. Optional scatter plot
    % ------------------------------------------------------------
    if make_feasibility_plot
        figure('Color','w', 'Position', [200 200 650 550]); hold on;

        scatter(all_spike_sim, all_dur_sim, 12, [0.8 0.8 0.8], 'o', ...
            'LineWidth', 0.6, 'MarkerEdgeAlpha', 0.8, 'DisplayName', 'Simultaneous');

        scatter(all_spike_seq, all_dur_seq, 12, [0.6 0.6 0.6], 's', ...
            'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Sequential');

        % Predicted lines from the linear model
        x_line = linspace(min(spike_all), max(spike_all), 100)';

        T_sim_line = table(x_line, zeros(size(x_line)), ...
            'VariableNames', {'SpikeCount', 'Condition'});
        y_sim_line = predict(mdl, T_sim_line);

        T_seq_line = table(x_line, ones(size(x_line)), ...
            'VariableNames', {'SpikeCount', 'Condition'});
        y_seq_line = predict(mdl, T_seq_line);

        plot(x_line, y_sim_line, '--k', 'LineWidth', 1.2, 'DisplayName', 'Fit Sim');
        plot(x_line, y_seq_line, '-k',  'LineWidth', 1.2, 'DisplayName', 'Fit Seq');

        xlabel('Spike Count', 'FontWeight', 'bold');
        ylabel('Duration (ms)', 'FontWeight', 'bold');
        title('Single-Dataset Feasibility Check: Spike Count vs Duration', 'FontWeight', 'bold');
        legend('Location', 'southeast', 'Box', 'off');
        box off;
    end

else
    fprintf('\nNot enough valid points for feasibility check.\n');
end

%% =================== 7. SAVE RESULTS ====================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX005';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

parts = split(folder_sim, filesep);
exp_id = parts{end};

out_filename = fullfile(save_dir, ['Result_Combined_SpikeDuration_Separate_DX005_' exp_id '.mat']);

save(out_filename, ...
    'CombinedData', 'Amps', 'PTDs_ms_seq', ...
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

    smooth_rate = smoothdata(rate_hz, 'gaussian', 10);

    base_mask = centers >= base_win(1) & centers <= base_win(2);
    base_data = smooth_rate(base_mask);

    mu = mean(base_data);
    sigma = std(base_data);
    if sigma == 0
        sigma = 1;
    end

    threshold = max(mu + nSD * sigma, 10);

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

    f = dir('*RespondingChannels.mat');
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