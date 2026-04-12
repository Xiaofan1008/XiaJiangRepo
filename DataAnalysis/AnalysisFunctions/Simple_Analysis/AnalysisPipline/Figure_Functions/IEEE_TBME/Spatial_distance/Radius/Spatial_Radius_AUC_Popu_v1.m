%% ============================================================
%   POPULATION SPATIAL ANALYSIS: AUC (Area Under the Curve)
%   - Input: Result_Spatial_RawRadius_*.mat files
%   - Logic: 
%       1. Spatial Integration: Calculates the Area Under the Curve (AUC)
%          for the spatial activation profile of every stimulation set.
%       2. Trapezoidal Rule: Uses trapz() to integrate activation over distance.
%       3. Statistics: Paired Wilcoxon signed-rank test (Sim vs. Seq).
%       4. Plotting: Amplitude vs. AUC (Mean ± SEM).
% ============================================================
clear;

%% ================= USER SETTINGS =================
% Paste your specific list of saved .mat files here
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

save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Spatial_Radius';
save_figures = true; 
tiff_dpi     = 600; 

%% ================= 1. DYNAMIC AMPLITUDE HARVESTING =================
fprintf('Scanning explicitly provided .mat files...\n');
num_files = length(file_paths);

% Find all unique amplitudes across the explicitly listed files
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
            % Logic: Sort electrodes by distance and integrate the binary response
            dists = D.Dist_to_Boundary;
            [sorted_dists, sort_idx] = sort(dists);
            
            % Calculate area using the Trapezoidal Rule
            auc_sim = trapz(sorted_dists, D.Clean_Sim(sort_idx));
            auc_seq = trapz(sorted_dists, D.Clean_Seq(sort_idx));
            
            % SAFETY: If 0% activation, AUC is 0 (Valid for dose-response)
            if D.TotalPerc_Sim == 0 && D.TotalPerc_Seq == 0
                % Optional: Set to NaN if you only want to see AUC for "active" sets
                % auc_sim = NaN; auc_seq = NaN;
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
    
    % MEANS & SEM
    PopResults.Amp(a).Sim_Mean = mean(auc_sim, 'omitnan');
    PopResults.Amp(a).Seq_Mean = mean(auc_seq, 'omitnan');
    PopResults.Amp(a).Sim_SEM  = std(auc_sim, 'omitnan') / sqrt(sum(~isnan(auc_sim)));
    PopResults.Amp(a).Seq_SEM  = std(auc_seq, 'omitnan') / sqrt(sum(~isnan(auc_seq)));
    
    % Statistics
    try [pval, ~] = signrank(auc_sim, auc_seq); catch, pval = NaN; end
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
% fig = figure('Units', 'centimeters', 'Position', [5, 5, 12, 11], 'Color', 'w', 'Name', 'Population Spatial AUC');
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

% Plotting (IEEE Style)
errorbar(v_amps, sim_m, sim_s, '--ok', 'LineWidth', 1, 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous', 'CapSize', 8);
errorbar(v_amps, seq_m, seq_s, '-sk', 'LineWidth', 1, 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential', 'CapSize', 8);

% % Asterisks
% for i = 1:length(v_amps)
%     if p_list(i) < 0.05
%         y_pos = max([sim_m(i)+sim_s(i), seq_m(i)+seq_s(i)]) + 10; 
%         text(v_amps(i), y_pos, '*', 'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
%     end
% end

% --- MULTI-LEVEL SIGNIFICANCE STARS ---
for i = 1:length(v_amps)
    sig_str = '';
    if p_list(i) < 0.001, sig_str = '***';
    elseif p_list(i) < 0.01, sig_str = '**';
    elseif p_list(i) < 0.05, sig_str = '*';
    end
    
    if ~isempty(sig_str)
        % Place asterisk string above the highest error bar
        y_pos = max([sim_m(i)+sim_s(i), seq_m(i)+seq_s(i)]) + (max(sim_m)*0.08); 
        text(v_amps(i), y_pos, sig_str, 'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', 'Helvetica');
    end
end

% xlabel('Current Amplitude (\muA)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Integrated Spatial Leakage (AUC)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Total Response Footprint (Trapezoidal AUC)', 'FontSize', 13);
% set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 11);
% legend('Location', 'northwest', 'Box', 'off');
% xlim([0, max(v_amps)+1]); axis square;

% --- IEEE FORMATTING ---
xlabel('Amplitude (µA)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Effective Activation Spread (µm)', 'FontSize', 9, 'FontName', 'Arial');
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1, 'FontSize', 9, 'FontName', 'Helvetica');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
% xlim([0, max(v_amps)+1]); 
% ylim([0, max([sim_m; seq_m])*1.3]); % Extra headroom for stars
xlim([0, 10]); 
ylim([0, 400]); 
axis square;

if save_figures
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    filename = fullfile(save_dir, 'Population_AUC_TuningCurve.tiff');
    exportgraphics(fig, filename, 'Resolution', tiff_dpi);
end