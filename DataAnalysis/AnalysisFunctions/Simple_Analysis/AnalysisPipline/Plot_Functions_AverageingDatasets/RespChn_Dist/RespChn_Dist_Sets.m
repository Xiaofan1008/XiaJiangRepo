%% ============================================================
%   GROUP SPATIAL ANALYSIS: Pooled Focality & Decay
%   - Input: Result_SpatialRaw_*.mat files (from Single Analysis)
%   - Logic:
%       1. Loop through all files -> Sets -> Amplitudes.
%       2. For each Set: Calculate Probability Profile & CDF.
%       3. Extract Metric: R90 (Radius containing 90% of response).
%       4. Pool data by Amplitude.
%       5. Plot:
%          - Fig 1: Avg Probability Histogram (Per Amp)
%          - Fig 2: Avg CDF Curve (Per Amp)
%          - Fig 3: Summary R90 vs Amplitude
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
% List your single-dataset result files here
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX012/Result_SpatialRaw_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX012/Result_SpatialRaw_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX012/Result_SpatialRaw_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX011/Result_SpatialRaw_Xia_Exp1_Sim9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX010/Result_SpatialRaw_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX009/Result_SpatialRaw_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX009/Result_SpatialRaw_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX006/Result_SpatialRaw_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX006/Result_SpatialRaw_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX006/Result_SpatialRaw_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX006/Result_SpatialRaw_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX005/Result_SpatialRaw_Xia_Exp1_Sim.mat';
};

% Plot Settings
save_figures = true;
save_dir     = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Group_Analysis/RespChn_Dist/';

% Re-define bins to ensure consistency across all files
dist_bin_edges = 0 : 50 : 1300; 

%% ================= 1. POOL DATA =================
fprintf('Pooling data from %d files...\n', length(file_paths));

% Structure to hold pooled data
% Pooled.Amp_X.Sim_Probs = [N_Sets x N_Bins]
% Pooled.Amp_X.Sim_CDFs  = [N_Sets x N_Bins]
% Pooled.Amp_X.Sim_R90   = [N_Sets x 1]
Pooled = struct();
All_Amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), warning('Missing: %s', file_paths{f}); continue; end
    D = load(file_paths{f});
    
    if ~isfield(D, 'SpatialData'), continue; end
    SD = D.SpatialData;
    
    % Loop Sets
    for s = 1:length(SD.Set)
        if isempty(SD.Set(s).Amp), continue; end
        
        % Loop Amplitudes
        for a = 1:length(SD.Set(s).Amp)
            Data = SD.Set(s).Amp(a);
            val = Data.Val;
            if val == 0, continue; end
            
            % Create Field Name (e.g., A_5p0)
            fName = sprintf('A_%.1f', val); fName = strrep(fName, '.', 'p');
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_Probs = []; Pooled.(fName).Seq_Probs = [];
                Pooled.(fName).Sim_CDFs  = []; Pooled.(fName).Seq_CDFs  = [];
                Pooled.(fName).Sim_R90   = []; Pooled.(fName).Seq_R90   = [];
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            % --- CALC METRICS FOR THIS SET ---
            [probs_sim, cdf_sim, r90_sim] = calc_metrics(Data.Dist, Data.Sim_Resp, Data.Sim_Valid, dist_bin_edges);
            [probs_seq, cdf_seq, r90_seq] = calc_metrics(Data.Dist, Data.Seq_Resp, Data.Seq_Valid, dist_bin_edges);
            
            % --- STORE ---
            Pooled.(fName).Sim_Probs = [Pooled.(fName).Sim_Probs; probs_sim];
            Pooled.(fName).Seq_Probs = [Pooled.(fName).Seq_Probs; probs_seq];
            
            Pooled.(fName).Sim_CDFs  = [Pooled.(fName).Sim_CDFs; cdf_sim];
            Pooled.(fName).Seq_CDFs  = [Pooled.(fName).Seq_CDFs; cdf_seq];
            
            Pooled.(fName).Sim_R90   = [Pooled.(fName).Sim_R90; r90_sim];
            Pooled.(fName).Seq_R90   = [Pooled.(fName).Seq_R90; r90_seq];
        end
    end
end

All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;

if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ================= 2. PLOT PROFILES (Per Amplitude) =================
fprintf('Generating Profiles for %d Amplitudes...\n', length(All_Amps));

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    % Extract Data
    Sim_P = Pooled.(fName).Sim_Probs; Seq_P = Pooled.(fName).Seq_Probs;
    Sim_C = Pooled.(fName).Sim_CDFs;  Seq_C = Pooled.(fName).Seq_CDFs;
    
    if isempty(Sim_P), continue; end
    
    % --- Calc Statistics (Mean & SEM) ---
    % Probs
    Avg_P_Sim = mean(Sim_P, 1, 'omitnan'); SEM_P_Sim = std(Sim_P, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Sim_P), 1));
    Avg_P_Seq = mean(Seq_P, 1, 'omitnan'); SEM_P_Seq = std(Seq_P, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Seq_P), 1));
    
    % CDFs
    Avg_C_Sim = mean(Sim_C, 1, 'omitnan'); SEM_C_Sim = std(Sim_C, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Sim_C), 1));
    Avg_C_Seq = mean(Seq_C, 1, 'omitnan'); SEM_C_Seq = std(Seq_C, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Seq_C), 1));
    
    % === FIGURE A: SPATIAL HISTOGRAM (Probability) ===
    figNameA = sprintf('SpatialHist_Group_%.1fuA', curr_amp);
    figure('Color','w', 'Position', [100 100 600 400], 'Name', figNameA); hold on;
    
    % Plot Grouped Bars (Sim vs Seq)
    b = bar(bin_centers, [Avg_P_Sim; Avg_P_Seq]', 'grouped');
    b(1).FaceColor = 'b'; b(1).EdgeColor = 'none'; b(1).DisplayName = 'Simultaneous';
    b(2).FaceColor = 'r'; b(2).EdgeColor = 'none'; b(2).DisplayName = 'Sequential';
    
    % Add Error Bars
    % Calculate X positions of the bars
    x_offset = b(1).XEndPoints - b(1).XData; 
    errorbar(bin_centers + x_offset(1), Avg_P_Sim, SEM_P_Sim, 'k.', 'LineWidth', 1, 'HandleVisibility','off');
    errorbar(bin_centers + x_offset(2), Avg_P_Seq, SEM_P_Seq, 'k.', 'LineWidth', 1, 'HandleVisibility','off');
    
    xlabel('Distance (\mum)'); ylabel('Response Probability');
    ylim([0 1.05]); title(sprintf('Avg Spatial Probability @ %.1f \\muA', curr_amp));
    legend('Location','best'); box off;
    
    % if save_figures, saveas(gcf, fullfile(save_dir, [figNameA '.png'])); end
    if save_figures, saveas(gcf, fullfile(save_dir, [figNameA '.fig'])); end
    
    % === FIGURE B: CDF CURVE (Cumulative) ===
    figNameB = sprintf('SpatialCDF_Group_%.1fuA', curr_amp);
    figure('Color','w', 'Position', [750 100 600 400], 'Name', figNameB); hold on;
    
    % Plot Shaded Error
    plot_shaded_error(bin_centers, Avg_C_Sim, SEM_C_Sim, 'b');
    plot_shaded_error(bin_centers, Avg_C_Seq, SEM_C_Seq, 'r');
    
    % Plot Mean Lines
    plot(bin_centers, Avg_C_Sim, '-o', 'Color', 'b', 'LineWidth', 2, 'MarkerFaceColor','b', 'DisplayName','Sim CDF');
    plot(bin_centers, Avg_C_Seq, '-s', 'Color', 'r', 'LineWidth', 2, 'MarkerFaceColor','r', 'DisplayName','Seq CDF');
    
    yline(0.9, '--k', '90% Recruitment', 'HandleVisibility','off');
    
    xlabel('Distance (\mum)'); ylabel('Cumulative Fraction');
    ylim([0 1.05]); title(sprintf('Avg Spatial Spread @ %.1f \\muA', curr_amp));
    legend('Location','southeast'); box off;
    
    % if save_figures, saveas(gcf, fullfile(save_dir, [figNameB '.png'])); end
    if save_figures, saveas(gcf, fullfile(save_dir, [figNameB '.fig'])); end
end

%% ================= 3. SUMMARY PLOT (R90 vs Amplitude) =================
fprintf('Generating Summary R90 Plot...\n');

Sim_R90_Mean = []; Sim_R90_SEM = [];
Seq_R90_Mean = []; Seq_R90_SEM = [];

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    r90_sim = Pooled.(fName).Sim_R90;
    r90_seq = Pooled.(fName).Seq_R90;
    
    % Filter out NaNs (e.g. if response was too weak to reach 90%)
    r90_sim(isnan(r90_sim)) = []; 
    r90_seq(isnan(r90_seq)) = [];
    
    Sim_R90_Mean(i) = mean(r90_sim); 
    Sim_R90_SEM(i)  = std(r90_sim) / sqrt(length(r90_sim));
    
    Seq_R90_Mean(i) = mean(r90_seq); 
    Seq_R90_SEM(i)  = std(r90_seq) / sqrt(length(r90_seq));
end

figNameSum = 'Spatial_Summary_R90';
figure('Color','w', 'Position', [400 400 600 500], 'Name', figNameSum); hold on;

errorbar(All_Amps, Sim_R90_Mean, Sim_R90_SEM, '-o', 'Color', 'b', 'LineWidth', 2, 'MarkerFaceColor','b', 'DisplayName','Simultaneous');
errorbar(All_Amps, Seq_R90_Mean, Seq_R90_SEM, '-s', 'Color', 'r', 'LineWidth', 2, 'MarkerFaceColor','r', 'DisplayName','Sequential');

xlabel('Amplitude (\muA)'); ylabel('Radius of 90% Recruitment (\mum)');
title('Spatial Focality Summary (R_{90})');
legend('Location','best'); box off; grid on;

% if save_figures, saveas(gcf, fullfile(save_dir, [figNameSum '.png'])); end
if save_figures, saveas(gcf, fullfile(save_dir, [figNameSum '.fig'])); end
fprintf('>>> Group Analysis Complete.\n');

%% ================= HELPER FUNCTIONS =================
function [prob_vec, cdf_vec, r90_val] = calc_metrics(dist, resp, valid, edges)
    % Calculates Probability Profile, CDF, and R90 for ONE Set
    
    prob_vec = nan(1, length(edges)-1);
    counts   = zeros(1, length(edges)-1);
    
    for b = 1:length(edges)-1
        mask = dist >= edges(b) & dist < edges(b+1);
        mask_valid = mask & valid;
        
        % Probability: # Resp / # Valid (Geometry Normalized)
        if sum(mask_valid) > 0
            n_resp = sum(resp(mask_valid));
            n_valid = sum(mask_valid);
            prob_vec(b) = n_resp / n_valid;
        else
            prob_vec(b) = NaN; % No electrode at this distance
        end
        
        % Counts (for CDF)
        if sum(mask_valid) > 0
            counts(b) = sum(resp(mask_valid));
        end
    end
    
    % CDF Calculation (Normalized 0-1)
    cdf_vec = cumsum(counts);
    total_resp = max(cdf_vec);
    if total_resp > 0
        cdf_vec = cdf_vec / total_resp;
    else
        cdf_vec = zeros(size(cdf_vec));
    end
    
    % R90 Calculation
    % Find first distance bin where CDF >= 0.90
    idx_90 = find(cdf_vec >= 0.90, 1, 'first');
    if ~isempty(idx_90)
        centers = edges(1:end-1) + diff(edges)/2;
        r90_val = centers(idx_90);
    else
        r90_val = NaN; % Didn't reach 90%
    end
end

function plot_shaded_error(x, y, se, col)
    % Helper to plot shaded error cloud
    valid = ~isnan(y) & ~isnan(se);
    x = x(valid); y = y(valid); se = se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end