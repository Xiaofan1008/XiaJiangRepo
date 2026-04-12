%% ============================================================
%   POPULATION SPATIAL ANALYSIS: Simultaneous vs. Sequential
%   - Logic: 
%       1. Dynamic Re-binning: Switch between 50um and 100um via 'bin_size'.
%       2. Master Pooling: Correctly averages/sums raw indices to ensure 0-100um visibility.
%       3. Statistics: Paired Wilcoxon signed-rank tests.
%       4. Plotting: Separated IEEE figures, Arial 9pt, Grey/Black theme.
%       5. Export: High-resolution .tiff files (600 DPI).
% ============================================================
clear; close all;

%% ================= USER SETTINGS =================
% 1. Toggle Bin Size Here (Set to 50 or 100)
bin_size = 100; 

% 2. File Paths
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX006_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX006_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX006_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX006_Xia_Exp1_Seq4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX009_Xia_Exp1_Seq5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq2_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX012_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Outer_Spatial_Analysis_100um_to600um';
save_figures = true; % Switch to true to export .tiff
tiff_dpi     = 600; 

%% ================= 1. DYNAMIC DATA HARVESTING =================
fprintf('Scanning files...\n');
num_files = length(file_paths);
raw_res   = 50; % Resolution of the raw data vectors in the .mat files

% A. Find all unique amplitudes
all_amps = [];
for f = 1:num_files
    tmp = load(file_paths{f}, 'Amps');
    all_amps = [all_amps; tmp.Amps(:)];
end
MasterAmps = unique(all_amps);
MasterAmps(MasterAmps == 0) = []; 

% B. Establish New Dynamic Bins
dist_bin_edges = 0 : bin_size : 600; 
bin_centers    = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;
num_bins       = length(bin_centers);
stride         = bin_size / raw_res; 

fprintf('Resolution: %dum. Bins: %d. Amps: ', bin_size, num_bins);
fprintf('%.1f ', MasterAmps); fprintf('uA\n');

%% ================= 2. MASTER BUCKET POOLING =================
Master_Prob_Sim    = cell(length(MasterAmps), 1);
Master_Prob_Seq    = cell(length(MasterAmps), 1);
Master_Density_Sim = cell(length(MasterAmps), 1);
Master_Density_Seq = cell(length(MasterAmps), 1);
Master_D90_Sim     = cell(length(MasterAmps), 1);
Master_D90_Seq     = cell(length(MasterAmps), 1);

for f = 1:num_files
    data = load(file_paths{f});
    for ss = 1:length(data.SpatialResults.Set)
        for ai = 1:length(data.SpatialResults.Set(ss).Amp)
            D = data.SpatialResults.Set(ss).Amp(ai);
            if isempty(D.Val) || D.Val == 0, continue; end
            
            m_idx = find(abs(MasterAmps - D.Val) < 0.001);
            if isempty(m_idx), continue; end
            
            % --- DYNAMIC RE-BINNING LOOP ---
            raw_p_sim = D.Prob_Outer_Sim(:)'; raw_p_seq = D.Prob_Outer_Seq(:)';
            raw_d_sim = D.Density_Outer_Sim(:)'; raw_d_seq = D.Density_Outer_Seq(:)';
            
            p_sim = nan(1, num_bins); p_seq = nan(1, num_bins);
            d_sim = nan(1, num_bins); d_seq = nan(1, num_bins);
            
            for b = 1:num_bins
                s_idx = ((b-1) * stride) + 1;
                e_idx = b * stride;
                
                if e_idx <= length(raw_p_sim)
                    p_sim(b) = mean(raw_p_sim(s_idx:e_idx), 'omitnan');
                    p_seq(b) = mean(raw_p_seq(s_idx:e_idx), 'omitnan');
                    d_sim(b) = sum(raw_d_sim(s_idx:e_idx), 'omitnan');
                    d_seq(b) = sum(raw_d_seq(s_idx:e_idx), 'omitnan');
                end
            end
            
            Master_Prob_Sim{m_idx}    = [Master_Prob_Sim{m_idx}; p_sim];
            Master_Prob_Seq{m_idx}    = [Master_Prob_Seq{m_idx}; p_seq];
            Master_Density_Sim{m_idx} = [Master_Density_Sim{m_idx}; d_sim];
            Master_Density_Seq{m_idx} = [Master_Density_Seq{m_idx}; d_seq];
            
            % D90 logic (independent of binning)
            d90_sim_v = D.D90_Sim; if D.TotalPerc_Sim == 0, d90_sim_v = NaN; end
            d90_seq_v = D.D90_Seq; if D.TotalPerc_Seq == 0, d90_seq_v = NaN; end
            Master_D90_Sim{m_idx} = [Master_D90_Sim{m_idx}; d90_sim_v];
            Master_D90_Seq{m_idx} = [Master_D90_Seq{m_idx}; d90_seq_v];
        end
    end
end

%% ================= 3. POPULATION MATH & STATISTICS =================
PopResults = struct();
for a = 1:length(MasterAmps)
    PopResults.Amp(a).Val = MasterAmps(a);
    p_sim = Master_Prob_Sim{a}; p_seq = Master_Prob_Seq{a};
    d_sim = Master_Density_Sim{a}; d_seq = Master_Density_Seq{a};
    d90_sim = Master_D90_Sim{a}; d90_seq = Master_D90_Seq{a};
    
    if isempty(p_sim), continue; end
    PopResults.Amp(a).N = sum(~isnan(d90_sim) | ~isnan(d90_seq));
    
    % Means & SEM
    PopResults.Amp(a).Prob_Sim_Mean = mean(p_sim, 1, 'omitnan');
    PopResults.Amp(a).Prob_Seq_Mean = mean(p_seq, 1, 'omitnan');
    PopResults.Amp(a).Dens_Sim_Mean = mean(d_sim, 1, 'omitnan');
    PopResults.Amp(a).Dens_Seq_Mean = mean(d_seq, 1, 'omitnan');
    PopResults.Amp(a).D90_Sim_Mean  = mean(d90_sim, 'omitnan');
    PopResults.Amp(a).D90_Seq_Mean  = mean(d90_seq, 'omitnan');
    
    PopResults.Amp(a).Prob_Sim_SEM = std(p_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(p_sim), 1));
    PopResults.Amp(a).Prob_Seq_SEM = std(p_seq, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(p_seq), 1));
    PopResults.Amp(a).Dens_Sim_SEM = std(d_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d_sim), 1));
    PopResults.Amp(a).Dens_Seq_SEM = std(d_seq, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d_seq), 1));
    PopResults.Amp(a).D90_Sim_SEM  = std(d90_sim, 'omitnan') / sqrt(sum(~isnan(d90_sim)));
    PopResults.Amp(a).D90_Seq_SEM  = std(d90_seq, 'omitnan') / sqrt(sum(~isnan(d90_seq)));
    
    % Stats
    try [pval_d90, ~] = signrank(d90_sim, d90_seq); catch; pval_d90 = NaN; end
    PopResults.Amp(a).D90_pval = pval_d90;
    
    PopResults.Amp(a).Prob_pvals = nan(1, num_bins);
    PopResults.Amp(a).Dens_pvals = nan(1, num_bins);
    for b = 1:num_bins
        try [p1, ~] = signrank(p_sim(:,b), p_seq(:,b)); catch; p1 = NaN; end
        try [p2, ~] = signrank(d_sim(:,b), d_seq(:,b)); catch; p2 = NaN; end
        PopResults.Amp(a).Prob_pvals(b) = p1;
        PopResults.Amp(a).Dens_pvals(b) = p2;
    end
end

%% ================= 4. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     PART A: D90 CONTAINMENT BOUNDARY SUMMARY                       \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-4s | %-15s | %-15s | %-10s | %-6s\n', 'Amp', 'N', 'D90 Sim (um)', 'D90 Seq (um)', 'p-value', 'Sig');
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a); if isempty(P.N) || P.N == 0, continue; end
    sig = ''; if P.D90_pval < 0.001, sig='***'; elseif P.D90_pval < 0.01, sig='**'; elseif P.D90_pval < 0.05, sig='*'; end
    fprintf('%4.1fuA | %-4d | %5.1f ± %-5.1f | %5.1f ± %-5.1f | %-10.6f | %s\n', P.Val, P.N, P.D90_Sim_Mean, P.D90_Sim_SEM, P.D90_Seq_Mean, P.D90_Seq_SEM, P.D90_pval, sig);
end

fprintf('\n====================================================================\n');
fprintf('     PART B: SPATIAL BREAKDOWN (%dum Bins)                       \n', bin_size);
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a); if isempty(P.N) || P.N == 0, continue; end
    fprintf('\n--- AMPLITUDE: %.1f uA (N = %d) ---\n', P.Val, P.N);
    fprintf('%-10s | %-18s | %-18s | %-8s\n', 'Bin (um)', 'Prob Sim (%)', 'Prob Seq (%)', 'p-val');
    for b = 1:num_bins
        fprintf('%4.0f-%-4.0f | %5.1f ± %-5.1f     | %5.1f ± %-5.1f     | %-8.4f\n', dist_bin_edges(b), dist_bin_edges(b+1), P.Prob_Sim_Mean(b)*100, P.Prob_Sim_SEM(b)*100, P.Prob_Seq_Mean(b)*100, P.Prob_Seq_SEM(b)*100, P.Prob_pvals(b));
    end
end

%% ================= 5. PLOTTING: SEPARATED FIGURES =================
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a); if isempty(P.N) || P.N < 2, continue; end 
    
    % --- FIGURE 1: PROBABILITY (CAPPED ERROR BARS) ---
    fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 8.8, 8], 'Color', 'w','Name', sprintf('Prob_%.1fuA', P.Val)); hold on;
    y1 = P.Prob_Sim_Mean * 100; e1 = P.Prob_Sim_SEM * 100;
    y2 = P.Prob_Seq_Mean * 100; e2 = P.Prob_Seq_SEM * 100;
    errorbar(bin_centers, y1, e1, '--ok', 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'CapSize', 3, 'DisplayName', 'Simultaneous');
    errorbar(bin_centers, y2, e2, '-sk', 'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'CapSize', 3, 'DisplayName', 'Sequential');
    for b = 1:num_bins
        p = P.Prob_pvals(b); s = ''; if p < 0.001, s='***'; elseif p < 0.01, s='**'; elseif p < 0.05, s='*'; end
        if ~isempty(s), text(bin_centers(b), max([y1(b)+e1(b), y2(b)+e2(b)])+8, s, 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end
    end
    xlabel('Distance (µm)', 'FontSize', 9, 'FontName', 'Arial'); ylabel('Activation Probability (%)', 'FontSize', 9, 'FontName', 'Arial');
    set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XTick', 0:100:700);
    xlim([0 600]); ylim([0 100]); axis square; legend('Location', 'northeast', 'Box', 'off');

    % --- FIGURE 2: DENSITY (WIDE BARS, NO CAPS) ---
    fig2 = figure('Units', 'centimeters', 'Position', [11, 2, 8.8, 8], 'Color', 'w','Name', sprintf('Dens_%.1fuA', P.Val)); hold on;
    b_plot = bar(bin_centers, [P.Dens_Sim_Mean(:), P.Dens_Seq_Mean(:)], 'grouped', 'BarWidth', 1.0);
    b_plot(1).FaceColor = [0.8 0.8 0.8]; 
    b_plot(1).DisplayName = 'Simultaneous';
    b_plot(2).FaceColor = [0 0 0];
    b_plot(2).DisplayName = 'Sequential';
    errorbar(b_plot(1).XEndPoints, P.Dens_Sim_Mean, zeros(size(P.Dens_Sim_SEM)), P.Dens_Sim_SEM, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 0, 'HandleVisibility','off');
    errorbar(b_plot(2).XEndPoints, P.Dens_Seq_Mean, zeros(size(P.Dens_Seq_SEM)), P.Dens_Seq_SEM, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 0, 'HandleVisibility','off');
    % for b = 1:num_bins
    %     p = P.Dens_pvals(b); s = ''; if p < 0.001, s='***'; elseif p < 0.01, s='**'; elseif p < 0.05, s='*'; end
    %     if ~isempty(s), text(bin_centers(b), max([P.Dens_Sim_Mean(b)+P.Dens_Sim_SEM(b), P.Dens_Seq_Mean(b)+P.Dens_Seq_SEM(b)])+0.02, s, 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end
    % end
    xlabel('Distance (µm)', 'FontSize', 9, 'FontName', 'Arial'); ylabel('Fraction of Total Response', 'FontSize', 9, 'FontName', 'Arial');
    set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off', 'XTick', 0:100:700);
    xlim([0 600]); ylim([0 0.25]); axis square; legend('Location', 'northeast', 'Box', 'off');

    if save_figures
        % Ensure the directory exists
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        
        % 1. Save the Probability Figure
        % Filename example: Prob_Profile_3.0uA_100um.tif
        fname1 = sprintf('Prob_Profile_%.1fuA_%dum.tif', P.Val, bin_size);
        fullPath1 = fullfile(save_dir, fname1);
        exportgraphics(fig1, fullPath1, 'Resolution', tiff_dpi, 'BackgroundColor', 'w');
        
        % 2. Save the Density Figure
        % Filename example: Dens_Profile_3.0uA_100um.tif
        fname2 = sprintf('Dens_Profile_%.1fuA_%dum.tif', P.Val, bin_size);
        fullPath2 = fullfile(save_dir, fname2);
        exportgraphics(fig2, fullPath2, 'Resolution', tiff_dpi, 'BackgroundColor', 'w');
        
        fprintf('   >>> Figures saved for %.1f uA\n', P.Val);
        
        % Optional: Close figures after saving to prevent MATLAB from slowing down
        % close(fig1); close(fig2); 
    end

end