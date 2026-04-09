%% ============================================================
%   GRAND AVERAGE: SPATIAL CROSSTALK (ALL AMPLITUDES)
%   - Metric: Baseline-Subtracted Spikes, Normalized to "Top 3" Sim Peak.
%   - Visual: 50 um bins (Grouped Bars + Individual Scatter).
%   - Stats: Macro-Zone Integration (200 um wide bins).
%   - Loop: Automatically detects all amplitudes and generates 1 figure per amp.
%   - Standard: IEEE TBME (8.89x8.89 cm, Arial 9pt).
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% --- A. File Paths ---
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX005_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX006_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX006_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX006_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX006_Xia_Exp1_Seq4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX009_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX009_Xia_Exp1_Seq5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq2_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX012_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance/SpikeCount_vs_Distance_DX016_Xia_Exp1_Seq_Full_4.mat';
};

% --- B. Spatial & Normalization Settings ---
visual_bin_width = 50;   % Resolution of the bar chart (um)
stats_zone_width = 200;  % Resolution of the statistical testing (um)
top_n_peak       = 3;    % Use mean of top N Simultaneous bins for 1.0 Normalization
min_n_threshold  = 5;    % Minimum number of valid datasets needed to run a stat test

% --- C. Plot Settings ---
save_figure = false;
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Spatial_Spread';

%% ================= PRE-SCAN: FIND ALL AMPLITUDES =================
fprintf('Scanning datasets to dynamically find all unique amplitudes...\n');
All_Amps = [];
for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    tmp = load(file_paths{i}, 'Results');
    if isfield(tmp, 'Results')
        All_Amps = [All_Amps, tmp.Results.Amps(:)'];
    end
end
Unique_Amps = unique(All_Amps);
Unique_Amps(Unique_Amps < 0.001) = []; 
fprintf('Found %d valid amplitudes to process.\n', length(Unique_Amps));
disp(Unique_Amps);

%% ================= MASTER AMPLITUDE LOOP =================
for aa = 1:length(Unique_Amps)
    target_amp = Unique_Amps(aa);
    
    fprintf('\n============================================================\n');
    fprintf('>>> PROCESSING AMPLITUDE: %.1f uA (%d of %d) <<<\n', target_amp, aa, length(Unique_Amps));
    fprintf('============================================================\n');
    
    bin_edges = -825:visual_bin_width:825; 
    bin_centers = -800:visual_bin_width:800;
    num_bins = length(bin_centers);
    
    Binned_Sim_All = []; 
    Binned_Seq_All = [];
    dataset_count = 0; 
    
    %% ---------------- 2. THE HARVESTER ----------------
    for i = 1:length(file_paths)
        if ~exist(file_paths{i}, 'file'), continue; end
        D = load(file_paths{i}); 
        if ~isfield(D, 'Results'), continue; end
        R = D.Results;
        
        amp_idx = find(abs(R.Amps - target_amp) < 0.001, 1);
        if isempty(amp_idx), continue; end
        
        for ss = 1:length(R.SpatialSets)
            if isempty(R.SpatialSets(ss).stimCh), continue; end
            SetData = R.SpatialSets(ss).Amp(amp_idx);
            
            dist = SetData.distances;
            sim_raw = SetData.Sim_Spikes;
            seq_raw = SetData.Seq_Spikes;
            if isempty(dist), continue; end
            
            b_sim = nan(1, num_bins); 
            b_seq = nan(1, num_bins);
            idx_bin = discretize(dist, bin_edges);
            for b = 1:num_bins
                hit = (idx_bin == b);
                if any(hit)
                    b_sim(b) = mean(sim_raw(hit), 'omitnan');
                    b_seq(b) = mean(seq_raw(hit), 'omitnan');
                end
            end
            
            % 3. ROBUST PEAK SELF-NORMALIZATION
            valid_sim = b_sim(~isnan(b_sim));
            valid_sim_sorted = sort(valid_sim, 'descend');
            
            valid_seq = b_seq(~isnan(b_seq));
            valid_seq_sorted = sort(valid_seq, 'descend');
            
            if isempty(valid_sim_sorted) || isempty(valid_seq_sorted), continue; end
            
            n_take_sim = min(top_n_peak, length(valid_sim_sorted));
            norm_factor_sim = mean(valid_sim_sorted(1:n_take_sim));
            
            n_take_seq = min(top_n_peak, length(valid_seq_sorted));
            norm_factor_seq = mean(valid_seq_sorted(1:n_take_seq));
            
            if norm_factor_sim <= 0 || norm_factor_seq <= 0, continue; end 
            
            % 4. Normalize and Store
            Binned_Sim_All(end+1, :) = b_sim / norm_factor_sim;
            Binned_Seq_All(end+1, :) = b_seq / norm_factor_seq;
            dataset_count = dataset_count + 1;
        end
    end
    fprintf('Successfully extracted %d valid sets for %.1f uA.\n', dataset_count, target_amp);
    
    if dataset_count == 0
        fprintf('No valid data found. Skipping figure generation.\n');
        continue; 
    end
    
    %% ---------------- 3. GRAND MEAN & STATS ENGINE ----------------
    Grand_Sim_Mean = mean(Binned_Sim_All, 1, 'omitnan');
    Grand_Sim_SEM  = std(Binned_Sim_All, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Binned_Sim_All), 1));
    
    Grand_Seq_Mean = mean(Binned_Seq_All, 1, 'omitnan');
    Grand_Seq_SEM  = std(Binned_Seq_All, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Binned_Seq_All), 1));
    
    max_extent = 800 - (stats_zone_width/2);
    zone_centers = -max_extent : stats_zone_width : max_extent;
    n_zones = length(zone_centers);
    Stat_Results = struct();
    
    for z = 1:n_zones
        ZC = zone_centers(z);
        z_min = ZC - (stats_zone_width/2);
        z_max = ZC + (stats_zone_width/2);
        bin_hits = find(bin_centers > z_min & bin_centers < z_max);
        if isempty(bin_hits), continue; end
        
        zone_sim_vals = mean(Binned_Sim_All(:, bin_hits), 2, 'omitnan');
        zone_seq_vals = mean(Binned_Seq_All(:, bin_hits), 2, 'omitnan');
        
        valid_pairs = ~isnan(zone_sim_vals) & ~isnan(zone_seq_vals);
        data_s = zone_sim_vals(valid_pairs);
        data_q = zone_seq_vals(valid_pairs);
        n_pairs = length(data_s);
        
        if n_pairs < min_n_threshold, continue; end
        
        p_val = signrank(data_s, data_q);
        mean_s = mean(data_s);
        mean_q = mean(data_q);
        perc_diff = ((mean_q - mean_s) / mean_s) * 100;
        
        txt = '';
        if p_val < 0.001, txt = '***';
        elseif p_val < 0.01, txt = '**';
        elseif p_val < 0.05, txt = '*';
        else, txt = 'n.s.';
        end
        
        fprintf('Zone: %+4.0f to %+4.0f um | N: %2d | Sim: %.2f | Seq: %.2f | Diff: %+5.1f%% | p = %.4f (%s)\n', ...
                z_min, z_max, n_pairs, mean_s, mean_q, perc_diff, p_val, txt);
                
        Stat_Results(z).z_min = z_min;
        Stat_Results(z).z_max = z_max;
        Stat_Results(z).txt = txt;
        Stat_Results(z).p_val = p_val;
        Stat_Results(z).bin_hits = bin_hits;
    end

    % ---> ADDED: 3.5 SHAPE SIMILARITY METRICS <---
    % Convert spatial means to Probability Distributions (Area Under Curve = 1)
    prob_sim = Grand_Sim_Mean / sum(Grand_Sim_Mean, 'omitnan');
    prob_seq = Grand_Seq_Mean / sum(Grand_Seq_Mean, 'omitnan');
    
    % Calculate Cumulative Distribution Functions (CDFs)
    cdf_sim = cumsum(prob_sim, 'omitnan');
    cdf_seq = cumsum(prob_seq, 'omitnan');
    
    % 1. K-S Statistic (Maximum distance between the two CDFs)
    ks_stat = max(abs(cdf_sim - cdf_seq));
    
    % 2. Cosine Similarity (Vector angle match: 1.0 is a perfect match)
    cos_sim = dot(prob_sim, prob_seq) / (norm(prob_sim) * norm(prob_seq));
    
    fprintf('\n--- SHAPE SIMILARITY (Drop-off Comparison) ---\n');
    fprintf('Cosine Similarity Score : %.3f (1.0 = mathematically identical shape)\n', cos_sim);
    fprintf('K-S Distance Statistic  : %.3f (Closer to 0 = identical spatial spread)\n', ks_stat);
    fprintf('----------------------------------------------\n\n');
    % -----------------------------------------------------------

    %% ---------------- 4. TBME-STYLE PLOTTER ----------------
    fig_name = sprintf('Spatial_Crosstalk_%.1fuA.tiff', target_amp);
    
    figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w', 'PaperPositionMode', 'auto', 'Name', fig_name); 
    hold on;
    
    % A. Background Scatter 
    jitter = 12; 
    for c = 1:num_bins
        vals_s = Binned_Sim_All(:, c); vals_s = vals_s(~isnan(vals_s));
        if ~isempty(vals_s)
            x_s = bin_centers(c) - 10 + (rand(size(vals_s))-0.5)*jitter; 
            scatter(x_s, vals_s, 6, [0.4 0.6 0.9], 'filled', 'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');
        end
        vals_q = Binned_Seq_All(:, c); vals_q = vals_q(~isnan(vals_q));
        if ~isempty(vals_q)
            x_q = bin_centers(c) + 10 + (rand(size(vals_q))-0.5)*jitter; 
            scatter(x_q, vals_q, 6, [0.9 0.6 0.4], 'filled', 'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');
        end
    end
    
    % B. Grouped Bar Chart
    b = bar(bin_centers, [Grand_Sim_Mean', Grand_Seq_Mean'], 'grouped', 'BarWidth', 0.85);
    b(1).FaceColor = [0 0.3 0.8]; b(1).EdgeColor = 'none'; b(1).DisplayName = 'Simultaneous';
    b(2).FaceColor = [0.85 0.33 0.10]; b(2).EdgeColor = 'none'; b(2).DisplayName = 'Sequential';
    
    % C. Add SEM Error Bars
    x1 = b(1).XEndPoints; x2 = b(2).XEndPoints;
    errorbar(x1, Grand_Sim_Mean, Grand_Sim_SEM, 'k', 'linestyle', 'none', 'LineWidth', 0.8, 'CapSize', 2, 'HandleVisibility', 'off');
    errorbar(x2, Grand_Seq_Mean, Grand_Seq_SEM, 'k', 'linestyle', 'none', 'LineWidth', 0.8, 'CapSize', 2, 'HandleVisibility', 'off');
    
    % D. Plot Statistical Brackets
    for z = 1:length(Stat_Results)
        if isempty(Stat_Results(z).txt) || strcmp(Stat_Results(z).txt, 'n.s.'), continue; end
        hits = Stat_Results(z).bin_hits;
        max_err_sim = max(Grand_Sim_Mean(hits) + Grand_Sim_SEM(hits));
        max_err_seq = max(Grand_Seq_Mean(hits) + Grand_Seq_SEM(hits));
        y_top = max(max_err_sim, max_err_seq) + 0.1; 
        
        x_start = min(bin_centers(hits)) - 20;
        x_end   = max(bin_centers(hits)) + 20;
        
        plot([x_start, x_end], [y_top, y_top], '-k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        plot([x_start, x_start], [y_top-0.05, y_top], '-k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        plot([x_end, x_end], [y_top-0.05, y_top], '-k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        text(mean([x_start, x_end]), y_top + 0.02, Stat_Results(z).txt, ...
            'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', 'Arial');
    end
    
    % E. Add Zero Lines
    xline(-100, '--k', 'HandleVisibility', 'off');
    xline(100, '--k', 'HandleVisibility', 'off');
    
    % F. Formatting
    box off; 
    set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.0);
    xlabel('Distance from Centroid (\mum)', 'FontSize', 9, 'FontName', 'Arial');
    ylabel('Normalized Spike Count', 'FontSize', 9,  'FontName', 'Arial');
    xlim([-800 800]);
    ylim([0 1.6]); 
    set(gca, 'YTick', 0:0.4:2.0);
    legend('Location','northeast', 'Box','off', 'FontSize', 8, 'FontName', 'Arial');
    title(sprintf('Spatial Recruitment (%.1f \\muA)', target_amp), 'FontSize', 10, 'FontName', 'Arial');
    
    if save_figure
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        save_path = fullfile(save_dir, fig_name);
        print(gcf, save_path, '-dtiff', '-r300');
        fprintf('>>> Figure saved to: %s\n', save_path);
    end
end % END OF MASTER AMPLITUDE LOOP

fprintf('\n>>> All Amplitudes Processed <<<\n');