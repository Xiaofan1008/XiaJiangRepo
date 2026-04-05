%% ============================================================
%   POPULATION SPATIAL ANALYSIS: Simultaneous vs. Sequential
%   - Logic: 
%       1. Explicit File Selection: Uses a hardcoded list of specific .mat files.
%       2. Dynamic Harvesting: Scans the selected files and aligns tested amplitudes.
%       3. Master Pooling: Pools every individual Set as an independent sample (N = Sets).
%       4. Statistics: Paired Wilcoxon signed-rank tests with NaN-omission for dead trials.
%       5. Plotting: Generates 1x2 spatial profiles per amplitude and 1 master D90 curve.
%       6. Export: Saves high-resolution .tiff files for publication.
% ============================================================
clear; close all;

%% ================= USER SETTINGS =================
% Paste your specific list of .mat files here
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
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim7.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim8.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX014_Xia_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim2.mat';
    % % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim5.mat';
    % % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim7.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Outer_Spatial_Analysis';
save_figures = false; % Set to true to export .tiff files
tiff_dpi     = 300;  % High resolution for IEEE publication

%% ================= 1. DYNAMIC DATA HARVESTING =================
fprintf('Scanning explicitly provided .mat files...\n');
num_files = length(file_paths);

if num_files == 0
    error('The file_paths cell array is empty. Please provide file paths.');
end

% A. Find all unique amplitudes across the explicitly listed files
all_amps_collected = [];
for f = 1:num_files
    % Modified to load directly from the cell array
    tmp = load(file_paths{f}, 'Amps');
    all_amps_collected = [all_amps_collected; tmp.Amps(:)];
end
MasterAmps = unique(all_amps_collected);
MasterAmps(MasterAmps == 0) = []; % Remove 0uA baseline if it exists

% B. Load the spatial bins from the first file to establish the X-axes
tmp = load(file_paths{1}, 'dist_bin_edges');
dist_bin_edges = tmp.dist_bin_edges;
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;
num_bins = length(bin_centers);

fprintf('Loaded %d files. Identified Master Amplitudes: ', num_files);
fprintf('%.1f ', MasterAmps); fprintf('uA\n');

%% ================= 2. MASTER BUCKET POOLING (N = Sets) =================
% Preallocate cell arrays to hold stacked data for each amplitude
% Structure: Master_Data{amp_idx} = Matrix of [Total Sets x Num Bins]
Master_Prob_Sim    = cell(length(MasterAmps), 1);
Master_Prob_Seq    = cell(length(MasterAmps), 1);
Master_Density_Sim = cell(length(MasterAmps), 1);
Master_Density_Seq = cell(length(MasterAmps), 1);

% For D90, it's just a 1D value per set, so Matrix of [Total Sets x 1]
Master_D90_Sim = cell(length(MasterAmps), 1);
Master_D90_Seq = cell(length(MasterAmps), 1);

total_sets_pooled = 0;

for f = 1:num_files
    % Modified to load directly from the cell array
    data = load(file_paths{f});
    nSets = length(data.SpatialResults.Set);
    
    for ss = 1:nSets
        total_sets_pooled = total_sets_pooled + 1;
        
        for ai = 1:length(data.SpatialResults.Set(ss).Amp)
            D = data.SpatialResults.Set(ss).Amp(ai);
            if isempty(D.Val) || D.Val == 0, continue; end
            
            % Find where this amplitude belongs in the Master list
            master_idx = find(abs(MasterAmps - D.Val) < 0.001);
            if isempty(master_idx), continue; end
            
            % This prevents the vertcat crash if older files have different bin lengths
            p_sim = D.Prob_Outer_Sim(:)'; p_seq = D.Prob_Outer_Seq(:)';
            d_sim = D.Density_Outer_Sim(:)'; d_seq = D.Density_Outer_Seq(:)';
            
            % Pad with NaNs if the loaded array is shorter than master num_bins
            if length(p_sim) < num_bins, p_sim = [p_sim, nan(1, num_bins - length(p_sim))]; end
            if length(p_seq) < num_bins, p_seq = [p_seq, nan(1, num_bins - length(p_seq))]; end
            if length(d_sim) < num_bins, d_sim = [d_sim, nan(1, num_bins - length(d_sim))]; end
            if length(d_seq) < num_bins, d_seq = [d_seq, nan(1, num_bins - length(d_seq))]; end
            
            % Truncate if the loaded array is longer than master num_bins
            p_sim = p_sim(1:num_bins); p_seq = p_seq(1:num_bins);
            d_sim = d_sim(1:num_bins); d_seq = d_seq(1:num_bins);
            
            % Append spatial profile data safely
            Master_Prob_Sim{master_idx}    = [Master_Prob_Sim{master_idx}; p_sim];
            Master_Prob_Seq{master_idx}    = [Master_Prob_Seq{master_idx}; p_seq];
            Master_Density_Sim{master_idx} = [Master_Density_Sim{master_idx}; d_sim];
            Master_Density_Seq{master_idx} = [Master_Density_Seq{master_idx}; d_seq];
            
            % Append D90 data (Magnitude). 
            % SAFETY FIX: If potency was exactly 0%, D90 should be NaN so it doesn't average to 0um.
            d90_sim_val = D.D90_Sim; if D.TotalPerc_Sim == 0, d90_sim_val = NaN; end
            d90_seq_val = D.D90_Seq; if D.TotalPerc_Seq == 0, d90_seq_val = NaN; end
            
            Master_D90_Sim{master_idx} = [Master_D90_Sim{master_idx}; d90_sim_val];
            Master_D90_Seq{master_idx} = [Master_D90_Seq{master_idx}; d90_seq_val];
        end
    end
end
fprintf('Total unique stimulation sets pooled (N): %d\n', total_sets_pooled);

%% ================= 3. POPULATION MATH & STATISTICS =================
PopResults = struct();

for a = 1:length(MasterAmps)
    PopResults.Amp(a).Val = MasterAmps(a);
    
    % Get data matrices for this amplitude
    p_sim = Master_Prob_Sim{a}; d_sim = Master_Density_Sim{a}; d90_sim = Master_D90_Sim{a};
    p_seq = Master_Prob_Seq{a}; d_seq = Master_Density_Seq{a}; d90_seq = Master_D90_Seq{a};
    
    % If no data exists for this amplitude, skip math
    if isempty(p_sim), continue; end
    
    % Valid N-count (ignoring NaNs)
    N_count = sum(~isnan(d90_sim) | ~isnan(d90_seq));
    PopResults.Amp(a).N = N_count;
    
    % --- MEANS (OmitNaN) ---
    PopResults.Amp(a).Prob_Sim_Mean = mean(p_sim, 1, 'omitnan');
    PopResults.Amp(a).Prob_Seq_Mean = mean(p_seq, 1, 'omitnan');
    PopResults.Amp(a).Dens_Sim_Mean = mean(d_sim, 1, 'omitnan');
    PopResults.Amp(a).Dens_Seq_Mean = mean(d_seq, 1, 'omitnan');
    PopResults.Amp(a).D90_Sim_Mean  = mean(d90_sim, 'omitnan');
    PopResults.Amp(a).D90_Seq_Mean  = mean(d90_seq, 'omitnan');
    
    % --- STANDARD ERROR OF THE MEAN (SEM) ---
    PopResults.Amp(a).Prob_Sim_SEM = std(p_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(p_sim), 1));
    PopResults.Amp(a).Prob_Seq_SEM = std(p_seq, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(p_seq), 1));
    PopResults.Amp(a).Dens_Sim_SEM = std(d_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d_sim), 1));
    PopResults.Amp(a).Dens_Seq_SEM = std(d_seq, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(d_seq), 1));
    PopResults.Amp(a).D90_Sim_SEM  = std(d90_sim, 'omitnan') / sqrt(sum(~isnan(d90_sim)));
    PopResults.Amp(a).D90_Seq_SEM  = std(d90_seq, 'omitnan') / sqrt(sum(~isnan(d90_seq)));
    
    % --- WILCOXON SIGNED-RANK TEST (Paired Non-Parametric) ---
    % D90 Stats
    try [pval_d90, ~] = signrank(d90_sim, d90_seq); catch; pval_d90 = NaN; end
    PopResults.Amp(a).D90_pval = pval_d90;
    
    % Spatial Array Stats (Calculated per distance bin)
    PopResults.Amp(a).Prob_pvals = nan(1, num_bins);
    PopResults.Amp(a).Dens_pvals = nan(1, num_bins);
    for b = 1:num_bins
        try [p_prob, ~] = signrank(p_sim(:,b), p_seq(:,b)); catch; p_prob = NaN; end
        try [p_dens, ~] = signrank(d_sim(:,b), d_seq(:,b)); catch; p_dens = NaN; end
        PopResults.Amp(a).Prob_pvals(b) = p_prob;
        PopResults.Amp(a).Dens_pvals(b) = p_dens;
    end
end

%% ================= 4. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     POPULATION STATISTICS: D90 CONTAINMENT BOUNDARY (Sim vs Seq)   \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-4s | %-15s | %-15s | %-10s | %-6s\n', 'Amp', 'N', 'D90 Sim (um)', 'D90 Seq (um)', 'p-value', 'Sig');
fprintf('--------------------------------------------------------------------\n');
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N == 0, continue; end
    
    sig_str = ''; if P.D90_pval < 0.05, sig_str = '*'; end
    if P.D90_pval < 0.01, sig_str = '**'; end
    if P.D90_pval < 0.001, sig_str = '***'; end
    
    fprintf('%4.1fuA | %-4d | %5.1f ± %-5.1f | %5.1f ± %-5.1f | %-10.4f | %s\n', ...
        P.Val, P.N, P.D90_Sim_Mean, P.D90_Sim_SEM, P.D90_Seq_Mean, P.D90_Seq_SEM, P.D90_pval, sig_str);
end
fprintf('\n');

%% ================= 5. PLOTTING: 1x2 PER AMPLITUDE (TIFF EXPORT) =================
if save_figures && ~exist(save_dir, 'dir'), mkdir(save_dir); end

for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N < 2, continue; end % Need at least N=2 for error bars
    
    fig = figure('Units', 'centimeters', 'Position', [2, 2, 20, 9], 'Color', 'w', ...
                 'Name', sprintf('Pop Spatial - %.1fuA', P.Val));
    t = tiledlayout(1, 2, 'TileSpacing', 'compact');
    
    % --- PANEL A: PROBABILITY DECAY ---
    nexttile; hold on;
    % Plot lines and standard vertical error bars
    errorbar(bin_centers, P.Prob_Sim_Mean, P.Prob_Sim_SEM, '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
    errorbar(bin_centers, P.Prob_Seq_Mean, P.Prob_Seq_SEM, '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
    
    % Add Significance Asterisks
    for b = 1:num_bins
        if P.Prob_pvals(b) < 0.05
            y_pos = max([P.Prob_Sim_Mean(b), P.Prob_Seq_Mean(b)]) + 0.1;
            text(bin_centers(b), y_pos, '*', 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
    
    title(sprintf('A. Leakage Probability (%.1f \\muA, N=%d)', P.Val, P.N)); 
    xlabel('Distance from Boundary (\mum)'); ylabel('Probability of Activation');
    ylim([0 1.2]); xlim([0 800]); grid off; box off;
    legend('Location', 'northeast', 'Box', 'off'); 

    % --- PANEL B: DENSITY FOOTPRINT ---
    nexttile; hold on;
    bar_data = [P.Dens_Sim_Mean(:), P.Dens_Seq_Mean(:)];
    b_plot = bar(bin_centers, bar_data, 'grouped');
    b_plot(1).FaceColor = 'w'; b_plot(1).EdgeColor = 'k'; b_plot(1).LineWidth = 1.2;
    b_plot(2).FaceColor = 'k'; b_plot(2).EdgeColor = 'k'; b_plot(2).LineWidth = 1.2;
    
    % Advanced Trick: Add error bars to grouped bars using XEndPoints
    x1 = b_plot(1).XEndPoints; x2 = b_plot(2).XEndPoints;
    errorbar(x1, P.Dens_Sim_Mean, P.Dens_Sim_SEM, 'k', 'LineStyle', 'none', 'LineWidth', 1);
    errorbar(x2, P.Dens_Seq_Mean, P.Dens_Seq_SEM, 'k', 'LineStyle', 'none', 'LineWidth', 1);
    
    % Add Significance Asterisks
    for b = 1:num_bins
        if P.Dens_pvals(b) < 0.05
            y_pos = max([P.Dens_Sim_Mean(b)+P.Dens_Sim_SEM(b), P.Dens_Seq_Mean(b)+P.Dens_Seq_SEM(b)]) + 0.05;
            text(bin_centers(b), y_pos, '*', 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
    
    title(sprintf('B. Spatial Density (%.1f \\muA, N=%d)', P.Val, P.N)); 
    xlabel('Distance from Boundary (\mum)'); ylabel('Fraction of Total Response');
    ylim([0 0.3]); xlim([0 800]); grid off; box off;
    
    % Export to TIFF
    if save_figures
        filename = fullfile(save_dir, sprintf('Pop_Spatial_Profile_%.1fuA.tiff', P.Val));
        exportgraphics(fig, filename, 'Resolution', tiff_dpi);
    end
end

%% ================= 6. PLOTTING: D90 TUNING CURVE (TIFF EXPORT) =================
fig_d90 = figure('Units', 'centimeters', 'Position', [5, 5, 12, 10], 'Color', 'w', 'Name', 'Population D90 Curve');
hold on; 

% Extract valid data for plotting
valid_amps = []; sim_means = []; sim_sems = []; seq_means = []; seq_sems = []; pvals = [];
for a = 1:length(MasterAmps)
    if isempty(PopResults.Amp(a).N) || PopResults.Amp(a).N < 2, continue; end
    valid_amps = [valid_amps; PopResults.Amp(a).Val];
    sim_means = [sim_means; PopResults.Amp(a).D90_Sim_Mean];
    sim_sems  = [sim_sems; PopResults.Amp(a).D90_Sim_SEM];
    seq_means = [seq_means; PopResults.Amp(a).D90_Seq_Mean];
    seq_sems  = [seq_sems; PopResults.Amp(a).D90_Seq_SEM];
    pvals     = [pvals; PopResults.Amp(a).D90_pval];
end

% Plot Lines with Error Bars
errorbar(valid_amps, sim_means, sim_sems, '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');
errorbar(valid_amps, seq_means, seq_sems, '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

% Add Significance Asterisks above the highest line at that amplitude
for i = 1:length(valid_amps)
    if pvals(i) < 0.05
        y_pos = max([sim_means(i)+sim_sems(i), seq_means(i)+seq_sems(i)]) + 20; % Offset visually
        text(valid_amps(i), y_pos, '*', 'FontSize', 16, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

xlabel('Current Amplitude (\muA)'); ylabel('Containment Boundary D_{90} (\mum)'); 
title('Population Maximum Spatial Extent'); axis square;
set(gca, 'TickDir', 'out', 'Box', 'off');
legend('Location', 'northwest', 'Box', 'off');

% Export to TIFF
if save_figures
    filename = fullfile(save_dir, 'Pop_D90_TuningCurve.tiff');
    exportgraphics(fig_d90, filename, 'Resolution', tiff_dpi);
    fprintf('\n>>> All High-Resolution .tiff Figures successfully saved to:\n%s\n', save_dir);
end