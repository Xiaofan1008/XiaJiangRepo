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
clear; 

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
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX013_Xia_Exp1_Seq_Sim3.mat';
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
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage/Result_Spatial_D90_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Outer_Spatial_Analysis';
save_figures = false; % Set to true to export .tiff files
tiff_dpi     = 600;  % High resolution for IEEE publication

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
% tmp = load(file_paths{1}, 'dist_bin_edges');
% dist_bin_edges = tmp.dist_bin_edges;
% bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;
% num_bins = length(bin_centers);

dist_bin_edges = 0 : 100 : 700; % Changed from 50 to 100
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;
num_bins = length(bin_centers);

fprintf('Loaded %d files. Bins set to 100um. Master Amplitudes: ', num_files);
fprintf('%.1f ', MasterAmps); fprintf('uA\n');

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
            % p_sim = D.Prob_Outer_Sim(:)'; p_seq = D.Prob_Outer_Seq(:)';
            % d_sim = D.Density_Outer_Sim(:)'; d_seq = D.Density_Outer_Seq(:)';
            % 
            % % Pad with NaNs if the loaded array is shorter than master num_bins
            % if length(p_sim) < num_bins, p_sim = [p_sim, nan(1, num_bins - length(p_sim))]; end
            % if length(p_seq) < num_bins, p_seq = [p_seq, nan(1, num_bins - length(p_seq))]; end
            % if length(d_sim) < num_bins, d_sim = [d_sim, nan(1, num_bins - length(d_sim))]; end
            % if length(d_seq) < num_bins, d_seq = [d_seq, nan(1, num_bins - length(d_seq))]; end
            % 
            % % Truncate if the loaded array is longer than master num_bins
            % p_sim = p_sim(1:num_bins); p_seq = p_seq(1:num_bins);
            % d_sim = d_sim(1:num_bins); d_seq = d_seq(1:num_bins);

            % Grab the raw 50um vectors saved in your .mat files
            raw_p_sim = D.Prob_Outer_Sim(:)'; raw_p_seq = D.Prob_Outer_Seq(:)';
            raw_d_sim = D.Density_Outer_Sim(:)'; raw_d_seq = D.Density_Outer_Seq(:)';

            % Preallocate the new 100um vectors
            p_sim = nan(1, num_bins); p_seq = nan(1, num_bins);
            d_sim = nan(1, num_bins); d_seq = nan(1, num_bins);

            for b = 1:num_bins
                % Map the 100um bin back to the two 50um indices in the file
                idx1 = (b*2) - 1; 
                idx2 = b*2;
                
                % Check if these indices exist in the raw file
                if idx2 <= length(raw_p_sim)
                    % Probability: Average the two 50um windows
                    p_sim(b) = mean([raw_p_sim(idx1), raw_p_sim(idx2)], 'omitnan');
                    p_seq(b) = mean([raw_p_seq(idx1), raw_p_seq(idx2)], 'omitnan');
                    
                    % Density: Sum the two 50um windows (since it's a fraction of total)
                    d_sim(b) = sum([raw_d_sim(idx1), raw_d_sim(idx2)], 'omitnan');
                    d_seq(b) = sum([raw_d_seq(idx1), raw_d_seq(idx2)], 'omitnan');
                end
            end
            
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
fprintf('     PART A: D90 CONTAINMENT BOUNDARY SUMMARY (Sim vs Seq)          \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-4s | %-15s | %-15s | %-10s | %-6s\n', 'Amp', 'N', 'D90 Sim (um)', 'D90 Seq (um)', 'p-value', 'Sig');
fprintf('--------------------------------------------------------------------\n');

for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N == 0, continue; end
    
    % Tiered Star Logic
    sig_str = ''; 
    if P.D90_pval < 0.001, sig_str = '***';
    elseif P.D90_pval < 0.01, sig_str = '**';
    elseif P.D90_pval < 0.05, sig_str = '*'; end
    
    fprintf('%4.1fuA | %-4d | %5.1f ± %-5.1f | %5.1f ± %-5.1f | %-10.6f | %s\n', ...
        P.Val, P.N, P.D90_Sim_Mean, P.D90_Sim_SEM, P.D90_Seq_Mean, P.D90_Seq_SEM, P.D90_pval, sig_str);
end

fprintf('\n====================================================================\n');
fprintf('     PART B: DETAILED SPATIAL PROFILES (100um Bins)                 \n');
fprintf('====================================================================\n');

for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N == 0, continue; end
    
    fprintf('\n--- AMPLITUDE: %.1f uA (N = %d) ---\n', P.Val, P.N);
    fprintf('%-10s | %-20s | %-20s | %-10s | %-4s\n', 'Bin (um)', 'Prob Sim (%)', 'Prob Seq (%)', 'p-val', 'Sig');
    fprintf('----------------------------------------------------------------------------\n');
    
    for b = 1:num_bins
        % Probability Stars
        p_prob = P.Prob_pvals(b);
        s_prob = '';
        % if p_prob < 0.001, s_prob = '***'; elseif p_prob < 0.01, s_prob = '**'; elseif p_prob < 0.05, s_prob = '*'; end
        
        % Note: We convert probability to % here to match your Figure 1
        fprintf('%4.0f-%-4.0f | %5.1f ± %-5.1f       | %5.1f ± %-5.1f       | %-10.4f | %s\n', ...
            dist_bin_edges(b), dist_bin_edges(b+1), ...
            P.Prob_Sim_Mean(b)*100, P.Prob_Sim_SEM(b)*100, ...
            P.Prob_Seq_Mean(b)*100, P.Prob_Seq_SEM(b)*100, ...
            p_prob, s_prob);
    end
    
    % Optional: Add a small gap or Density summary here if needed
    fprintf('   [Spatial Density Mean (Sim vs Seq)]: ');
    fprintf('%.2f vs %.2f at first bin.\n', P.Dens_Sim_Mean(1), P.Dens_Seq_Mean(1));
end
fprintf('\n======================= END OF STATISTICS ==========================\n');

%% ================= 5. PLOTTING: SEPARATED IEEE FIGURES (Arial 9pt) =================
% if save_figures && ~exist(save_dir, 'dir'), mkdir(save_dir); end
% 
% for a = 1:length(MasterAmps)
%     P = PopResults.Amp(a);
%     if isempty(P.N) || P.N < 2, continue; end 
% 
%     % --- FIGURE 1: LEAKAGE PROBABILITY (SHADED RIBBONS) ---
%     fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 8.8, 8], 'Color', 'w', 'Name', sprintf('Prob_%.1fuA', P.Val));
%     hold on;
% 
%     % Data Prep (convert to %)
%     y_sim = P.Prob_Sim_Mean * 100; err_sim = P.Prob_Sim_SEM * 100;
%     y_seq = P.Prob_Seq_Mean * 100; err_seq = P.Prob_Seq_SEM * 100;
% 
%     % Use standard errorbar with Caps
%     errorbar(bin_centers, y_sim, err_sim, '--ok', 'LineWidth', 1.2, 'MarkerSize', 4, ...
%         'MarkerFaceColor', 'w', 'CapSize', 4, 'DisplayName', 'Sim');
% 
%     errorbar(bin_centers, y_seq, err_seq, '-sk', 'LineWidth', 1.2, 'MarkerSize', 4, ...
%         'MarkerFaceColor', 'k', 'CapSize', 4, 'DisplayName', 'Seq');
% 
%     % Tiered Significance Stars
%     for b = 1:num_bins
%         pval = P.Prob_pvals(b);
%         sig_str = '';
%         if pval < 0.001, sig_str = '***'; elseif pval < 0.01, sig_str = '**'; elseif pval < 0.05, sig_str = '*'; end
% 
%         if ~isempty(sig_str)
%             y_pos = max([y_sim(b), y_seq(b)]) + 10;
%             text(bin_centers(b), y_pos, sig_str, 'FontSize', 9, 'FontName', 'Arial', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
%         end
%     end
% 
%     % Formatting
%     xlabel('Distance (µm)', 'FontSize', 9, 'FontName', 'Arial');
%     ylabel('Activation Probability (%)', 'FontSize', 9, 'FontName', 'Arial');
%     set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off');
%     xlim([0 700]); ylim([0 100]); axis square;
%     legend('Location', 'northeast', 'Box', 'off');
% 
%     if save_figures
%         exportgraphics(fig1, fullfile(save_dir, sprintf('Prob_Profile_%.1fuA.tiff', P.Val)), 'Resolution', tiff_dpi);
%     end
% 
%     % --- FIGURE 2: SPATIAL DENSITY (CLEAR LEGEND & WIDE BARS) ---
%     fig2 = figure('Units', 'centimeters', 'Position', [11, 2, 8.8, 8], 'Color', 'w', ...
%                   'Name', sprintf('Dens_%.1fuA', P.Val));
%     hold on;
% 
%     % 1. Increase BarWidth to 1.0 so bars fill the space better
%     % 2. Explicitly set DisplayNames for a clear legend
%     b_plot = bar(bin_centers, [P.Dens_Sim_Mean(:), P.Dens_Seq_Mean(:)], 'grouped', 'BarWidth', 1.0);
% 
%     % Simultaneous Style (Grey)
%     b_plot(1).FaceColor = [0.8 0.8 0.8]; 
%     b_plot(1).EdgeColor = 'k'; 
%     b_plot(1).LineWidth = 0.8;
%     b_plot(1).DisplayName = 'Simultaneous'; % <--- Change this for legend
% 
%     % Sequential Style (Black)
%     b_plot(2).FaceColor = [0 0 0];       
%     b_plot(2).EdgeColor = 'k'; 
%     b_plot(2).LineWidth = 0.8;
%     b_plot(2).DisplayName = 'Sequential';   % <--- Change this for legend
% 
%     % --- REFINED ERROR BARS: No Caps, Slightly Thicker ---
%     errorbar(b_plot(1).XEndPoints, P.Dens_Sim_Mean, zeros(size(P.Dens_Sim_SEM)), P.Dens_Sim_SEM, ...
%         'k', 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 0, 'HandleVisibility','off');
% 
%     errorbar(b_plot(2).XEndPoints, P.Dens_Seq_Mean, zeros(size(P.Dens_Seq_SEM)), P.Dens_Seq_SEM, ...
%         'k', 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 0, 'HandleVisibility','off');
% 
%     % (Significance stars logic here remains the same)
% 
%     % --- FORMATTING (Arial 9pt) ---
%     xlabel('Distance (\mum)', 'FontSize', 9, 'FontName', 'Arial');
%     ylabel('Fraction of Total Response', 'FontSize', 9, 'FontName', 'Arial');
%     set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.5, 'Box', 'off');
%     set(gca, 'XTick', 0:100:700); % Finer X-axis ticks for wider bars
%     xlim([0 700]); ylim([0 0.2]); % Extra headroom for stars
%     axis square;
% 
%     % Legend with clear labels
%     legend('Location', 'northeast', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
% 
%     if save_figures
%         exportgraphics(fig2, fullfile(save_dir, sprintf('Dens_Profile_%.1fuA.tiff', P.Val)), 'Resolution', tiff_dpi);
%     end
% end


if save_figures && ~exist(save_dir, 'dir'), mkdir(save_dir); end

for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N < 2, continue; end 
    
    % --- FIGURE 1: LEAKAGE PROBABILITY (CAPPED ERROR BARS) ---
    fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 8.8, 8], 'Color', 'w', ...
                  'Name', sprintf('Prob_%.1fuA', P.Val));
    hold on;
    
    y_sim = P.Prob_Sim_Mean * 100; err_sim = P.Prob_Sim_SEM * 100;
    y_seq = P.Prob_Seq_Mean * 100; err_seq = P.Prob_Seq_SEM * 100;
    
    % Use standard errorbar with Caps
    errorbar(bin_centers, y_sim, err_sim, '--ok', 'LineWidth', 1.2, 'MarkerSize', 4, ...
        'MarkerFaceColor', 'w', 'CapSize', 4, 'DisplayName', 'Simultaneous');
    errorbar(bin_centers, y_seq, err_seq, '-sk', 'LineWidth', 1.2, 'MarkerSize', 4, ...
        'MarkerFaceColor', 'k', 'CapSize', 4, 'DisplayName', 'Sequential');
        
    % Tiered Significance Stars
    for b = 1:num_bins
        pval = P.Prob_pvals(b);
        sig_str = '';
        if pval < 0.001, sig_str = '***'; elseif pval < 0.01, sig_str = '**'; elseif pval < 0.05, sig_str = '*'; end
        if ~isempty(sig_str)
            y_pos = max([y_sim(b)+err_sim(b), y_seq(b)+err_seq(b)]) + 8;
            text(bin_centers(b), y_pos, sig_str, 'FontSize', 9, 'FontName', 'Arial', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
    
    xlabel('Distance (µm)', 'FontSize', 9, 'FontName', 'Arial');
    ylabel('Activation Probability (%)', 'FontSize', 9, 'FontName', 'Arial');
    set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off');
    set(gca, 'XTick', 0:100:700);
    xlim([0 700]); ylim([0 100]); % Extra headroom for stars
    axis square; legend('Location', 'northeast', 'Box', 'off');
   
    if save_figures
        exportgraphics(fig1, fullfile(save_dir, sprintf('Prob_Profile_%.1fuA.tiff', P.Val)), 'Resolution', tiff_dpi);
    end

    % --- FIGURE 2: SPATIAL DENSITY (WIDE BARS, NO CAPS) ---
    fig2 = figure('Units', 'centimeters', 'Position', [11, 2, 8.8, 8], 'Color', 'w', ...
                  'Name', sprintf('Dens_%.1fuA', P.Val));
    hold on;
    
    b_plot = bar(bin_centers, [P.Dens_Sim_Mean(:), P.Dens_Seq_Mean(:)], 'grouped', 'BarWidth', 1.0);
    b_plot(1).FaceColor = [0.8 0.8 0.8]; b_plot(1).EdgeColor = 'k'; b_plot(1).LineWidth = 0.8; b_plot(1).DisplayName = 'Simultaneous';
    b_plot(2).FaceColor = [0 0 0];       b_plot(2).EdgeColor = 'k'; b_plot(2).LineWidth = 0.8; b_plot(2).DisplayName = 'Sequential';
    
    % Error Bars: No Caps, clean vertical lines
    errorbar(b_plot(1).XEndPoints, P.Dens_Sim_Mean, zeros(size(P.Dens_Sim_SEM)), P.Dens_Sim_SEM, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 0, 'HandleVisibility','off');
    errorbar(b_plot(2).XEndPoints, P.Dens_Seq_Mean, zeros(size(P.Dens_Seq_SEM)), P.Dens_Seq_SEM, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 0, 'HandleVisibility','off');
    
    % Tiered Significance Stars
    % for b = 1:num_bins
    %     pval = P.Dens_pvals(b);
    %     sig_str = '';
    %     if pval < 0.001, sig_str = '***'; elseif pval < 0.01, sig_str = '**'; elseif pval < 0.05, sig_str = '*'; end
    %     if ~isempty(sig_str)
    %         y_top = max([P.Dens_Sim_Mean(b)+P.Dens_Sim_SEM(b), P.Dens_Seq_Mean(b)+P.Dens_Seq_SEM(b)]) + 0.02;
    %         text(bin_centers(b), y_top, sig_str, 'FontSize', 9, 'FontName', 'Arial', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    %     end
    % end
    
    xlabel('Distance (µm)', 'FontSize', 9, 'FontName', 'Arial');
    ylabel('Fraction of Total Response', 'FontSize', 9, 'FontName', 'Arial');
    set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1, 'Box', 'off');
    set(gca, 'XTick', 0:100:700);
    xlim([0 700]); ylim([0 0.25]); % Increased limit for 100um binning
    axis square; legend('Location', 'northeast', 'Box', 'off');
    
    if save_figures
        exportgraphics(fig2, fullfile(save_dir, sprintf('Dens_Profile_%.1fuA.tiff', P.Val)), 'Resolution', tiff_dpi);
    end
end