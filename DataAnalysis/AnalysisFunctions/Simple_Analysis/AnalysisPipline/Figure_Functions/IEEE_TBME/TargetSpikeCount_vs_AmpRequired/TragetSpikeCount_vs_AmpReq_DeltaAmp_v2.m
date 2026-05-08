%% ============================================================
%   GRAND AVERAGE: REQUIRED AMPLITUDE FOR TARGET NORMALIZED SPIKE COUNT
%   - Logic: Inverts each DATASET recruitment curve first, then aggregates.
%   - Grouping: Bins data by target spike count.
%   - Interpolation: Linear interpolation within each dataset only.
%   - Stage 1: Uses a slightly broader target range for exploration.
%   - Output 1: Required amplitude vs target spike count
%   - Output 2: Delta amplitude (Sim - Seq) vs target spike count
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% --- INVERSE CURVE SETTINGS ---
% Stage 1 exploratory target range (slightly broader)
target_spike_counts = 0.20 : 0.1 : 1.60;

% Minimum number of paired datasets contributing at a target
stats_min_n_threshold = 2;

% Whether to force monotonic increase before inverse interpolation
% (Recommended for stage 1 to reduce instability)
enforce_monotonic = true;

% Plot Settings
save_figure = true;
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure2/AmpRequired_vs_SpikeCount';
fig_name_inverse = 'RequiredAmplitude_vs_TargetSpikeCount_v1.tiff';
fig_name_delta   = 'DeltaAmplitude_vs_TargetSpikeCount_v1.tiff';

% List all result files to include in the Grand Average
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX009/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX009/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX005/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
};

%% ================= 2. AGGREGATE DATA =================
fprintf('Processing %d datasets...\n', length(file_paths));

Pool_Sim = [];         % [target_spike_count, required_amp]
Pool_Seq = [];         % [target_spike_count, required_amp]
Pool_Paired = [];      % [target_spike_count, required_amp_sim, required_amp_seq]

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    D = load(file_paths{i});
    if ~isfield(D, 'ResultNorm'), continue; end

    R = D.ResultNorm;
    Amps = R.Amps(:);

    sim_ds_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));

    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_ds_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_ds_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end

    sim_ds_mean = sim_ds_mean(:);
    seq_ds_mean = seq_ds_mean(:);

    % Force 0 amplitude to 0 response
    zero_idx = abs(Amps) < 0.001;
    sim_ds_mean(zero_idx) = 0;
    seq_ds_mean(zero_idx) = 0;

    % Remove NaNs separately
    valid_sim = ~isnan(Amps) & ~isnan(sim_ds_mean);
    valid_seq = ~isnan(Amps) & ~isnan(seq_ds_mean);

    amps_sim = Amps(valid_sim);
    resp_sim = sim_ds_mean(valid_sim);

    amps_seq = Amps(valid_seq);
    resp_seq = seq_ds_mean(valid_seq);

    % Sort by amplitude
    [amps_sim, idx_sim] = sort(amps_sim);
    resp_sim = resp_sim(idx_sim);

    [amps_seq, idx_seq] = sort(amps_seq);
    resp_seq = resp_seq(idx_seq);

    % Optional monotonic enforcement for inverse interpolation
    if enforce_monotonic
        resp_sim = cummax(resp_sim);
        resp_seq = cummax(resp_seq);
    end

    % Remove duplicate response values before interp1
    [resp_sim_u, ia_sim] = unique(resp_sim, 'stable');
    amps_sim_u = amps_sim(ia_sim);

    [resp_seq_u, ia_seq] = unique(resp_seq, 'stable');
    amps_seq_u = amps_seq(ia_seq);

    % Need at least 2 points to interpolate
    if length(resp_sim_u) < 2 || length(resp_seq_u) < 2
        continue;
    end

    % Invert each dataset first: target spike count -> required amplitude
    for t = 1:length(target_spike_counts)
        target_val = target_spike_counts(t);

        req_amp_sim = nan;
        req_amp_seq = nan;

        if target_val >= min(resp_sim_u) && target_val <= max(resp_sim_u)
            req_amp_sim = interp1(resp_sim_u, amps_sim_u, target_val, 'linear');
        end

        if target_val >= min(resp_seq_u) && target_val <= max(resp_seq_u)
            req_amp_seq = interp1(resp_seq_u, amps_seq_u, target_val, 'linear');
        end

        % Unpaired pools
        if ~isnan(req_amp_sim)
            Pool_Sim = [Pool_Sim; target_val, req_amp_sim];
        end
        if ~isnan(req_amp_seq)
            Pool_Seq = [Pool_Seq; target_val, req_amp_seq];
        end

        % Paired pool
        if ~isnan(req_amp_sim) && ~isnan(req_amp_seq)
            Pool_Paired = [Pool_Paired; target_val, req_amp_sim, req_amp_seq];
        end
    end
end

%% ================= 3. STATS & VECTORS =================
Unique_Targets = unique([Pool_Sim(:,1); Pool_Seq(:,1)]);
Unique_Targets = sort(Unique_Targets);

Grand_Sim_Mean = []; Grand_Sim_SEM = []; Grand_Sim_N = [];
Grand_Seq_Mean = []; Grand_Seq_SEM = []; Grand_Seq_N = [];

Grand_Delta_Mean = []; Grand_Delta_SEM = []; Grand_Delta_N = [];

for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);

    vs = Pool_Sim(abs(Pool_Sim(:,1)-target_val)<0.001, 2);
    Grand_Sim_Mean(k) = mean(vs, 'omitnan');
    Grand_Sim_SEM(k)  = std(vs, 'omitnan') / sqrt(length(vs));
    Grand_Sim_N(k)    = length(vs);

    vq = Pool_Seq(abs(Pool_Seq(:,1)-target_val)<0.001, 2);
    Grand_Seq_Mean(k) = mean(vq, 'omitnan');
    Grand_Seq_SEM(k)  = std(vq, 'omitnan') / sqrt(length(vq));
    Grand_Seq_N(k)    = length(vq);

    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);
    delta_vals = current_pairs(:,2) - current_pairs(:,3);   % Sim - Seq

    Grand_Delta_Mean(k) = mean(delta_vals, 'omitnan');
    Grand_Delta_SEM(k)  = std(delta_vals, 'omitnan') / sqrt(length(delta_vals));
    Grand_Delta_N(k)    = length(delta_vals);
end

% --- Add origin point for plotting only ---
Plot_Targets    = [0; Unique_Targets];

Plot_Sim_Mean   = [0, Grand_Sim_Mean];
Plot_Sim_SEM    = [0, Grand_Sim_SEM];

Plot_Seq_Mean   = [0, Grand_Seq_Mean];
Plot_Seq_SEM    = [0, Grand_Seq_SEM];

Plot_Delta_Mean = [0, Grand_Delta_Mean];
Plot_Delta_SEM  = [0, Grand_Delta_SEM];

%% ================= 4. PLOT 1: INVERSE CURVES =================
figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w', 'PaperPositionMode', 'auto'); hold on;

% A. Scatter (Background)
jitter_w = 0.02;
for i = 1:size(Pool_Seq, 1)
    scatter(Pool_Seq(i,1)+(rand-0.5)*jitter_w, Pool_Seq(i,2), 12, [0.7 0.7 0.7], 's', ...
        'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
for i = 1:size(Pool_Sim, 1)
    scatter(Pool_Sim(i,1)+(rand-0.5)*jitter_w, Pool_Sim(i,2), 12, [0.8 0.8 0.8], 'o', ...
        'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end

% B. Main Curves
% errorbar(Unique_Targets, Grand_Sim_Mean, Grand_Sim_SEM, '.', 'Color', 'k', ...
%     'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
% p1 = plot(Unique_Targets, Grand_Sim_Mean, '--o', 'Color', 'k', 'LineWidth', 1.5, ...
%     'MarkerSize', 6, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');
% 
% errorbar(Unique_Targets, Grand_Seq_Mean, Grand_Seq_SEM, '.', 'Color', 'k', ...
%     'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
% p2 = plot(Unique_Targets, Grand_Seq_Mean, '-s', 'Color', 'k', 'LineWidth', 1.5, ...
%     'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

errorbar(Plot_Targets, Plot_Sim_Mean, Plot_Sim_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p1 = plot(Plot_Targets, Plot_Sim_Mean, '--o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');

errorbar(Plot_Targets, Plot_Seq_Mean, Plot_Seq_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p2 = plot(Plot_Targets, Plot_Seq_Mean, '-s', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

% C. Paired Signed Rank Test
fprintf('\n=== INVERSE CURVE STATS (Paired Signed Rank) ===\n');
fprintf('Min N = %d\n', stats_min_n_threshold);

for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);

    data_s = current_pairs(:,2);   % required amp for Sim
    data_q = current_pairs(:,3);   % required amp for Seq
    n_pairs = length(data_s);

    if n_pairs < stats_min_n_threshold
        fprintf('Target %.2f: Skipped (Low N=%d)\n', target_val, n_pairs);
        continue;
    end

    p = signrank(data_s, data_q);

    txt = '';
    if p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    end

    if ~isempty(txt)
        y_top = max(Grand_Sim_Mean(k) + Grand_Sim_SEM(k), Grand_Seq_Mean(k) + Grand_Seq_SEM(k));
        text(target_val, y_top + 0.30, txt, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontName', 'Arial');
        fprintf('Target %.2f: MARKED with %s (p=%.5f)\n', target_val, txt, p);
    else
        fprintf('Target %.2f: Not Significant (p=%.5f)\n', target_val, p);
    end
end

% --- Formatting ---
box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.5);
axis square;

xlabel('Matched Response Level (a.u.)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Required Amplitude (\muA)', 'FontSize', 9, 'FontName', 'Arial');

legend([p1, p2], 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

% xlim([min(Unique_Targets), max(Unique_Targets)]);
% ylim([0, max([Grand_Sim_Mean + Grand_Sim_SEM, Grand_Seq_Mean + Grand_Seq_SEM]) + 1.4]);

xlim([0, max(Unique_Targets)]);
% ylim([0, max([Plot_Sim_Mean + Plot_Sim_SEM, Plot_Seq_Mean + Plot_Seq_SEM]) + 1.4]);
ylim([0, 12]);
set(gca, 'XTick', 0 : 0.2 : max(Unique_Targets));

if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    save_path = fullfile(save_dir, fig_name_inverse);
    print(gcf, save_path, '-dtiff', '-r300');
    fprintf('Inverse figure saved to: %s\n', save_path);
end

%% ================= 5. PLOT 2: DELTA AMPLITUDE =================
figure('Units', 'centimeters', 'Position', [12, 2, 8.89, 8.89], 'Color', 'w', 'PaperPositionMode', 'auto'); hold on;

% A. Scatter (Background)
jitter_w = 0.02;
for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);
    delta_vals = current_pairs(:,2) - current_pairs(:,3);   % Sim - Seq

    for i = 1:length(delta_vals)
        scatter(target_val+(rand-0.5)*jitter_w, delta_vals(i), 12, [0.7 0.7 0.7], 'o', ...
            'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
    end
end

% B. Main Curve
% errorbar(Unique_Targets, Grand_Delta_Mean, Grand_Delta_SEM, '.', 'Color', 'k', ...
%     'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
% p3 = plot(Unique_Targets, Grand_Delta_Mean, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
%     'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Sim - Seq');

errorbar(Plot_Targets, Plot_Delta_Mean, Plot_Delta_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p3 = plot(Plot_Targets, Plot_Delta_Mean, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Sim - Seq');

yline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');

% C. Stats
fprintf('\n=== DELTA AMPLITUDE STATS (Paired Signed Rank) ===\n');
for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);

    data_s = current_pairs(:,2);
    data_q = current_pairs(:,3);
    n_pairs = length(data_s);

    if n_pairs < stats_min_n_threshold
        fprintf('Target %.2f: Skipped (Low N=%d)\n', target_val, n_pairs);
        continue;
    end

    p = signrank(data_s, data_q);

    txt = '';
    if p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    end

    if ~isempty(txt)
        y_top = Grand_Delta_Mean(k) + Grand_Delta_SEM(k);
        text(target_val, y_top + 0.15, txt, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontName', 'Arial');
    end
end

% --- Formatting ---
box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.5);
axis square;

xlabel('Matched Response Level (a.u.)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('\Delta Required Amplitude (Sim - Seq, \muA)', 'FontSize', 9, 'FontName', 'Arial');

legend([p3], 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

% xlim([min(Unique_Targets), max(Unique_Targets)]);
xlim([0, max(Unique_Targets)]);
ylim([-2,6]);
set(gca, 'XTick', 0 : 0.2 : max(Unique_Targets));

if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    save_path = fullfile(save_dir, fig_name_delta);
    print(gcf, save_path, '-dtiff', '-r300');
    fprintf('Delta figure saved to: %s\n', save_path);
end