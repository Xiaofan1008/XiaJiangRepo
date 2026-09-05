%% ============================================================
%   POPULATION SPATIAL ANALYSIS: AUC + MATCHED EFFECTIVE SPATIAL SPREAD
%   - Input: Result_Spatial_RawRadius_*.mat files
%   - Original Logic:
%       1. Spatial Integration: Calculates the Area Under the Curve (AUC)
%          for the spatial activation profile of every stimulation set.
%       2. Trapezoidal Rule: Uses trapz() to integrate activation over distance.
%       3. Statistics: Paired Wilcoxon signed-rank test (Sim vs. Seq).
%       4. Plotting: Amplitude vs. AUC (Mean ± SEM).
%
%   - NEW Logic:
%       5. Swapped-Axis Test: For each stimulation set, invert the AUC curve
%          to estimate the required amplitude for a matched effective spatial spread.
%       6. Plot matched-spread required amplitude and delta required amplitude.
% ============================================================
clear;

%% ================= USER SETTINGS =================
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX005_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX009_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX009_Xia_Exp1_Seq5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq2_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX012_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/TargetSpread_AmpReq';
save_figures = false;
tiff_dpi     = 600;

% --- NEW: Swapped-axis test settings (AUC-based matched spread) ---
target_auc_values = 50:50:400;   % Adjust later if needed
stats_min_n_threshold = 2;
enforce_monotonic = true;

%% ================= 1. DYNAMIC AMPLITUDE HARVESTING =================
fprintf('Scanning explicitly provided .mat files...\n');
num_files = length(file_paths);

all_amps_collected = [];
for f = 1:num_files
    tmp = load(file_paths{f}, 'Amps');
    all_amps_collected = [all_amps_collected; tmp.Amps(:)];
end
MasterAmps = unique(all_amps_collected);
MasterAmps(MasterAmps == 0) = [];

fprintf('Loaded %d files. Identified Amplitudes: ', num_files);
fprintf('%.1f ', MasterAmps); fprintf('uA\n');

%% ================= 2. MASTER BUCKET POOLING (N = Sets) =================
Master_AUC_Sim = cell(length(MasterAmps), 1);
Master_AUC_Seq = cell(length(MasterAmps), 1);
total_sets_pooled = 0;

for f = 1:num_files
    data = load(file_paths{f});
    nSets = length(data.SpatialResults.Set);

    for ss = 1:nSets
        total_sets_pooled = total_sets_pooled + 1;

        for ai = 1:length(data.SpatialResults.Set(ss).Amp)
            D = data.SpatialResults.Set(ss).Amp(ai);
            if isempty(D.Val) || D.Val == 0, continue; end

            m_idx = find(abs(MasterAmps - D.Val) < 0.001);
            if isempty(m_idx), continue; end

            % --- AUC CALCULATION (Numerical Integration) ---
            dists = D.Dist_to_Boundary;
            [sorted_dists, sort_idx] = sort(dists);

            auc_sim = trapz(sorted_dists, D.Clean_Sim(sort_idx));
            auc_seq = trapz(sorted_dists, D.Clean_Seq(sort_idx));

            if D.TotalPerc_Sim == 0 && D.TotalPerc_Seq == 0
                % Keep as 0 for now
            end

            Master_AUC_Sim{m_idx} = [Master_AUC_Sim{m_idx}; auc_sim];
            Master_AUC_Seq{m_idx} = [Master_AUC_Seq{m_idx}; auc_seq];
        end
    end
end
fprintf('Total unique stimulation sets pooled (N): %d\n', total_sets_pooled);

%% ================= 3. POPULATION MATH & STATISTICS =================
PopResults = struct();
for a = 1:length(MasterAmps)
    PopResults.Amp(a).Val = MasterAmps(a);

    auc_sim = Master_AUC_Sim{a};
    auc_seq = Master_AUC_Seq{a};

    if isempty(auc_sim), continue; end

    N_count = sum(~isnan(auc_sim) & ~isnan(auc_seq));
    PopResults.Amp(a).N = N_count;

    PopResults.Amp(a).Sim_Mean = mean(auc_sim, 'omitnan');
    PopResults.Amp(a).Seq_Mean = mean(auc_seq, 'omitnan');
    PopResults.Amp(a).Sim_SEM  = std(auc_sim, 'omitnan') / sqrt(sum(~isnan(auc_sim)));
    PopResults.Amp(a).Seq_SEM  = std(auc_seq, 'omitnan') / sqrt(sum(~isnan(auc_seq)));

    try
        [pval, ~] = signrank(auc_sim, auc_seq);
    catch
        pval = NaN;
    end
    PopResults.Amp(a).pval = pval;
end

%% ================= 4. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     POPULATION STATISTICS: SPATIAL LEAKAGE MAGNITUDE (AUC)         \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-4s | %-15s | %-15s | %-10s | %-6s\n', 'Amp', 'N', 'AUC Sim', 'AUC Seq', 'p-value', 'Sig');
fprintf('--------------------------------------------------------------------\n');
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N == 0, continue; end

    sig_str = '';
    if P.pval < 0.001, sig_str = '***';
    elseif P.pval < 0.01, sig_str = '**';
    elseif P.pval < 0.05, sig_str = '*'; end

    fprintf('%4.1fuA | %-4d | %5.1f ± %-5.1f | %5.1f ± %-5.1f | %-10.10f | %s\n', ...
        P.Val, P.N, P.Sim_Mean, P.Sim_SEM, P.Seq_Mean, P.Seq_SEM, P.pval, sig_str);
end

%% ================= 5. PLOTTING: AUC TUNING CURVE =================
fig = figure('Units', 'centimeters', 'Position', [5, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Population Spatial AUC');
hold on;

v_amps = []; sim_m = []; sim_s = []; seq_m = []; seq_s = []; p_list = [];
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N < 2, continue; end
    v_amps = [v_amps; P.Val];
    sim_m = [sim_m; P.Sim_Mean]; sim_s = [sim_s; P.Sim_SEM];
    seq_m = [seq_m; P.Seq_Mean]; seq_s = [seq_s; P.Seq_SEM];
    p_list = [p_list; P.pval];
end

errorbar(v_amps, sim_m, sim_s, '--ok', 'LineWidth', 1, 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous', 'CapSize', 8);
errorbar(v_amps, seq_m, seq_s, '-sk', 'LineWidth', 1, 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential', 'CapSize', 8);

for i = 1:length(v_amps)
    sig_str = '';
    if p_list(i) < 0.001, sig_str = '***';
    elseif p_list(i) < 0.01, sig_str = '**';
    elseif p_list(i) < 0.05, sig_str = '*';
    end

    if ~isempty(sig_str)
        y_pos = max([sim_m(i)+sim_s(i), seq_m(i)+seq_s(i)]) + (max(sim_m)*0.08);
        text(v_amps(i), y_pos, sig_str, 'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', 'Helvetica');
    end
end

xlabel('Amplitude (µA)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Effective Spatial Spread Index (a.u.)', 'FontSize', 9, 'FontName', 'Arial');
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1, 'FontSize', 9, 'FontName', 'Helvetica');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
xlim([0, 10]);
ylim([0, 400]);
axis square;

if save_figures
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    filename = fullfile(save_dir, 'Population_AUC_TuningCurve.tiff');
    exportgraphics(fig, filename, 'Resolution', tiff_dpi);
end

%% ================= 6. SMALL SWAPPED-AXIS TEST: MATCHED EFFECTIVE SPATIAL SPREAD =================
fprintf('\nGenerating swapped-axis test: required amplitude for matched effective spatial spread...\n');

Pool_Target_Sim = [];     % [target_auc, required_amp]
Pool_Target_Seq = [];     % [target_auc, required_amp]
Pool_Paired_Target = [];  % [target_auc, required_amp_sim, required_amp_seq]

for f = 1:num_files
    data = load(file_paths{f});
    nSets = length(data.SpatialResults.Set);

    for ss = 1:nSets
        amps_set = [];
        auc_sim_set = [];
        auc_seq_set = [];

        for ai = 1:length(data.SpatialResults.Set(ss).Amp)
            D = data.SpatialResults.Set(ss).Amp(ai);
            if isempty(D.Val) || D.Val == 0, continue; end

            % --- AUC CALCULATION (same logic as original code) ---
            dists = D.Dist_to_Boundary;
            [sorted_dists, sort_idx] = sort(dists);

            auc_sim = trapz(sorted_dists, D.Clean_Sim(sort_idx));
            auc_seq = trapz(sorted_dists, D.Clean_Seq(sort_idx));

            amps_set(end+1,1) = D.Val; %#ok<SAGROW>
            auc_sim_set(end+1,1) = auc_sim; %#ok<SAGROW>
            auc_seq_set(end+1,1) = auc_seq; %#ok<SAGROW>
        end

        if length(amps_set) < 2
            continue;
        end

        [amps_set, idx_sort] = sort(amps_set);
        auc_sim_set = auc_sim_set(idx_sort);
        auc_seq_set = auc_seq_set(idx_sort);

        valid_sim = ~isnan(amps_set) & ~isnan(auc_sim_set);
        valid_seq = ~isnan(amps_set) & ~isnan(auc_seq_set);

        amps_sim = amps_set(valid_sim);
        amps_seq = amps_set(valid_seq);
        auc_sim  = auc_sim_set(valid_sim);
        auc_seq  = auc_seq_set(valid_seq);

        if length(amps_sim) < 2 || length(amps_seq) < 2
            continue;
        end

        if enforce_monotonic
            auc_sim = cummax(auc_sim);
            auc_seq = cummax(auc_seq);
        end

        [auc_sim_u, ia_sim] = unique(auc_sim, 'stable');
        amps_sim_u = amps_sim(ia_sim);

        [auc_seq_u, ia_seq] = unique(auc_seq, 'stable');
        amps_seq_u = amps_seq(ia_seq);

        if length(auc_sim_u) < 2 || length(auc_seq_u) < 2
            continue;
        end

        for t = 1:length(target_auc_values)
            target_val = target_auc_values(t);

            req_amp_sim = NaN;
            req_amp_seq = NaN;

            if target_val >= min(auc_sim_u) && target_val <= max(auc_sim_u)
                req_amp_sim = interp1(auc_sim_u, amps_sim_u, target_val, 'linear');
            end

            if target_val >= min(auc_seq_u) && target_val <= max(auc_seq_u)
                req_amp_seq = interp1(auc_seq_u, amps_seq_u, target_val, 'linear');
            end

            if ~isnan(req_amp_sim)
                Pool_Target_Sim = [Pool_Target_Sim; target_val, req_amp_sim];
            end
            if ~isnan(req_amp_seq)
                Pool_Target_Seq = [Pool_Target_Seq; target_val, req_amp_seq];
            end
            if ~isnan(req_amp_sim) && ~isnan(req_amp_seq)
                Pool_Paired_Target = [Pool_Paired_Target; target_val, req_amp_sim, req_amp_seq];
            end
        end
    end
end

%% ================= 7. POPULATION MATH FOR MATCHED-SPREAD ANALYSIS =================
Unique_Targets = unique([Pool_Target_Sim(:,1); Pool_Target_Seq(:,1)]);
Unique_Targets = sort(Unique_Targets);

Grand_Target_Sim_Mean = []; Grand_Target_Sim_SEM = []; Grand_Target_Sim_N = [];
Grand_Target_Seq_Mean = []; Grand_Target_Seq_SEM = []; Grand_Target_Seq_N = [];
Grand_Target_Del_Mean = []; Grand_Target_Del_SEM = []; Grand_Target_Del_N = [];
Grand_Target_P = [];

for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);

    vs = Pool_Target_Sim(abs(Pool_Target_Sim(:,1)-target_val)<0.001, 2);
    vq = Pool_Target_Seq(abs(Pool_Target_Seq(:,1)-target_val)<0.001, 2);

    Grand_Target_Sim_Mean(k) = mean(vs, 'omitnan');
    Grand_Target_Sim_SEM(k)  = std(vs, 'omitnan') / sqrt(length(vs));
    Grand_Target_Sim_N(k)    = length(vs);

    Grand_Target_Seq_Mean(k) = mean(vq, 'omitnan');
    Grand_Target_Seq_SEM(k)  = std(vq, 'omitnan') / sqrt(length(vq));
    Grand_Target_Seq_N(k)    = length(vq);

    current_pairs = Pool_Paired_Target(abs(Pool_Paired_Target(:,1)-target_val)<0.001, :);
    delta_vals = current_pairs(:,2) - current_pairs(:,3);   % Sim - Seq

    Grand_Target_Del_Mean(k) = mean(delta_vals, 'omitnan');
    Grand_Target_Del_SEM(k)  = std(delta_vals, 'omitnan') / sqrt(length(delta_vals));
    Grand_Target_Del_N(k)    = length(delta_vals);

    if length(delta_vals) >= stats_min_n_threshold
        try
            Grand_Target_P(k) = signrank(current_pairs(:,2), current_pairs(:,3));
        catch
            Grand_Target_P(k) = NaN;
        end
    else
        Grand_Target_P(k) = NaN;
    end
end

%% ================= 8. PLOT: REQUIRED AMPLITUDE FOR MATCHED EFFECTIVE SPATIAL SPREAD =================
% --- NEW: Add origin point for plotting only ---
Plot_Targets = [0; Unique_Targets(:)];

Plot_Target_Sim_Mean = [0, Grand_Target_Sim_Mean];
Plot_Target_Sim_SEM  = [0, Grand_Target_Sim_SEM];

Plot_Target_Seq_Mean = [0, Grand_Target_Seq_Mean];
Plot_Target_Seq_SEM  = [0, Grand_Target_Seq_SEM];

Plot_Target_Del_Mean = [0, Grand_Target_Del_Mean];
Plot_Target_Del_SEM  = [0, Grand_Target_Del_SEM];

fig2 = figure('Units', 'centimeters', 'Position', [15, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Matched_AUC_RequiredAmplitude');
hold on;

% jitter_w = 8;
% for i = 1:size(Pool_Target_Seq, 1)
%     scatter(Pool_Target_Seq(i,1)+(rand-0.5)*jitter_w, Pool_Target_Seq(i,2), 12, [0.7 0.7 0.7], 's', ...
%         'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
% end
% for i = 1:size(Pool_Target_Sim, 1)
%     scatter(Pool_Target_Sim(i,1)+(rand-0.5)*jitter_w, Pool_Target_Sim(i,2), 12, [0.8 0.8 0.8], 'o', ...
%         'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
% end

% errorbar(Unique_Targets, Grand_Target_Sim_Mean, Grand_Target_Sim_SEM, '.', 'Color', 'k', ...
%     'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
% p1 = plot(Unique_Targets, Grand_Target_Sim_Mean, '--o', 'Color', 'k', 'LineWidth', 1.5, ...
%     'MarkerFaceColor', 'w', 'MarkerSize', 5, 'DisplayName', 'Simultaneous');
% 
% errorbar(Unique_Targets, Grand_Target_Seq_Mean, Grand_Target_Seq_SEM, '.', 'Color', 'k', ...
%     'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
% p2 = plot(Unique_Targets, Grand_Target_Seq_Mean, '-s', 'Color', 'k', 'LineWidth', 1.5, ...
%     'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sequential');

errorbar(Plot_Targets, Plot_Target_Sim_Mean, Plot_Target_Sim_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p1 = plot(Plot_Targets, Plot_Target_Sim_Mean, '--o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5, 'DisplayName', 'Simultaneous');

errorbar(Plot_Targets, Plot_Target_Seq_Mean, Plot_Target_Seq_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p2 = plot(Plot_Targets, Plot_Target_Seq_Mean, '-s', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sequential');

for k = 1:length(Unique_Targets)
    if isnan(Grand_Target_P(k)), continue; end

    sig_str = '';
    if Grand_Target_P(k) < 0.001, sig_str = '***';
    elseif Grand_Target_P(k) < 0.01, sig_str = '**';
    elseif Grand_Target_P(k) < 0.05, sig_str = '*';
    end

    if ~isempty(sig_str)
        y_pos = max([Grand_Target_Sim_Mean(k)+Grand_Target_Sim_SEM(k), Grand_Target_Seq_Mean(k)+Grand_Target_Seq_SEM(k)]) + 0.35;
        text(Unique_Targets(k), y_pos, sig_str, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontName', 'Arial');
    end
end

xlabel('Matched Effective Spread Index (a.u.)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Required Amplitude (µA)', 'FontSize', 9, 'FontName', 'Arial');
legend([p1, p2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1, 'FontSize', 9, 'FontName', 'Arial');
% xlim([min(Unique_Targets)-10, max(Unique_Targets)+10]);
xlim([0, max(Unique_Targets)]);
% set(gca, 'XTick', Unique_Targets);
set(gca, 'XTick', Plot_Targets);
% ylim([0, ceil(max([Grand_Target_Sim_Mean + Grand_Target_Sim_SEM, Grand_Target_Seq_Mean + Grand_Target_Seq_SEM]) + 1)]);
ylim([0, ceil(max([Grand_Target_Sim_Mean + Grand_Target_Sim_SEM, Grand_Target_Seq_Mean + Grand_Target_Seq_SEM]))]);

axis square;

if save_figures
    filename = fullfile(save_dir, 'Population_AUC_MatchedSpread_RequiredAmplitude.tiff');
    exportgraphics(fig2, filename, 'Resolution', tiff_dpi);
end

%% ================= 9. PLOT: DELTA REQUIRED AMPLITUDE FOR MATCHED EFFECTIVE SPATIAL SPREAD =================
fig3 = figure('Units', 'centimeters', 'Position', [25, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Matched_AUC_DeltaAmplitude');
hold on;

% jitter_w = 8;
% for k = 1:length(Unique_Targets)
%     target_val = Unique_Targets(k);
%     current_pairs = Pool_Paired_Target(abs(Pool_Paired_Target(:,1)-target_val)<0.001, :);
%     delta_vals = current_pairs(:,2) - current_pairs(:,3);   % Sim - Seq
% 
%     for i = 1:length(delta_vals)
%         scatter(target_val+(rand-0.5)*jitter_w, delta_vals(i), 12, [0.7 0.7 0.7], 'o', ...
%             'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
%     end
% end

% errorbar(Unique_Targets, Grand_Target_Del_Mean, Grand_Target_Del_SEM, '.', 'Color', 'k', ...
%     'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
% p3 = plot(Unique_Targets, Grand_Target_Del_Mean, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
%     'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sim - Seq');

errorbar(Plot_Targets, Plot_Target_Del_Mean, Plot_Target_Del_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p3 = plot(Plot_Targets, Plot_Target_Del_Mean, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sim - Seq');

yline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Matched Effective Spread Index (a.u.)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Required Amplitude Difference (Sim - Seq, µA)', 'FontSize', 9, 'FontName', 'Arial');
legend([p3], {'Sim - Seq'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1, 'FontSize', 9, 'FontName', 'Arial');
% xlim([min(Unique_Targets)-10, max(Unique_Targets)+10]);
xlim([0, max(Unique_Targets)]);
% set(gca, 'XTick', Unique_Targets);
set(gca, 'XTick', Plot_Targets);
ylim([floor(min([Grand_Target_Del_Mean - Grand_Target_Del_SEM, 0]) - 0.5), ceil(max([Grand_Target_Del_Mean + Grand_Target_Del_SEM]) + 0.8)]);
axis square;

if save_figures
    filename = fullfile(save_dir, 'Population_AUC_MatchedSpread_DeltaAmplitude.tiff');
    exportgraphics(fig3, filename, 'Resolution', tiff_dpi);
end

%% ================= 10. COMMAND WINDOW SUMMARY: MATCHED-SPREAD ANALYSIS =================
fprintf('\n====================================================================\n');
fprintf('   MATCHED EFFECTIVE SPATIAL SPREAD: REQUIRED AMPLITUDE SUMMARY      \n');
fprintf('====================================================================\n');
fprintf('%-12s | %-6s | %-15s | %-15s | %-15s | %-10s\n', ...
    'Target AUC', 'N', 'Req Amp Sim', 'Req Amp Seq', 'Delta (Sim-Seq)', 'p-value');
fprintf('----------------------------------------------------------------------------------------\n');

for k = 1:length(Unique_Targets)
    n_pairs = Grand_Target_Del_N(k);

    if n_pairs < stats_min_n_threshold
        fprintf('%-12.1f | %-6d | %-15.2f | %-15.2f | %-15.2f | Skipped\n', ...
            Unique_Targets(k), n_pairs, Grand_Target_Sim_Mean(k), Grand_Target_Seq_Mean(k), Grand_Target_Del_Mean(k));
        continue;
    end

    fprintf('%-12.1f | %-6d | %5.2f ± %-8.2f | %5.2f ± %-8.2f | %5.2f ± %-8.2f | %s\n', ...
        Unique_Targets(k), n_pairs, ...
        Grand_Target_Sim_Mean(k), Grand_Target_Sim_SEM(k), ...
        Grand_Target_Seq_Mean(k), Grand_Target_Seq_SEM(k), ...
        Grand_Target_Del_Mean(k), Grand_Target_Del_SEM(k), ...
        format_p(Grand_Target_P(k)));
end

%% ================= HELPER FUNCTIONS =================
function p_str = format_p(p)
    if isnan(p)
        p_str = 'NaN';
    else
        p_str = sprintf('%.20f', p);
    end
end