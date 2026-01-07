%% ============================================================
%   GROUP SPATIAL ANALYSIS: Pooled Focality & Decay
%   - Input: Result_SpatialRaw_*.mat files (from Single Analysis)
%   - Outputs (Per Amplitude):
%       1. Spatial Density Hist (Count / Valid Channels in Bin)
%       2. Global Occupancy Hist (Count / Total Valid Channels in Array)
%   - Summary Output:
%       3. Absolute Range (R_ext) vs Amplitude
%   - STYLE: Black & White, 50um Ticks
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
% List your single-dataset result files here
file_paths = {
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX015/Result_SpatialRaw_Xia_Seq_Sim7.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX014/Result_SpatialRaw_Xia_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim7.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/RespChn_Dist/DX013/Result_SpatialRaw_Xia_Exp1_Seq_Sim8.mat';
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
save_figures = false;
save_dir     = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Group_Analysis/RespChn_Dist/';
dist_bin_edges = 0 : 50 : 1300; 

%% ================= 1. POOL DATA =================
fprintf('Pooling data from %d files...\n', length(file_paths));

Pooled = struct();
All_Amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), warning('Missing: %s', file_paths{f}); continue; end
    D = load(file_paths{f});
    if ~isfield(D, 'SpatialData'), continue; end
    SD = D.SpatialData;
    
    for s = 1:length(SD.Set)
        if isempty(SD.Set(s).Amp), continue; end
        
        for a = 1:length(SD.Set(s).Amp)
            Data = SD.Set(s).Amp(a);
            val = Data.Val;
            if val == 0, continue; end
            
            fName = sprintf('A_%.1f', val); fName = strrep(fName, '.', 'p');
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_Probs   = []; Pooled.(fName).Seq_Probs   = [];
                Pooled.(fName).Sim_Global  = []; Pooled.(fName).Seq_Global  = [];
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            % --- Calc Metrics ---
            [prob_sim] = calc_metrics(Data.Dist, Data.Sim_Resp, Data.Sim_Valid, dist_bin_edges, 'prob');
            [prob_seq] = calc_metrics(Data.Dist, Data.Seq_Resp, Data.Seq_Valid, dist_bin_edges, 'prob');
            
            [glob_sim] = calc_metrics(Data.Dist, Data.Sim_Resp, Data.Sim_Valid, dist_bin_edges, 'global');
            [glob_seq] = calc_metrics(Data.Dist, Data.Seq_Resp, Data.Seq_Valid, dist_bin_edges, 'global');
            
            % --- Store ---
            Pooled.(fName).Sim_Probs   = [Pooled.(fName).Sim_Probs; prob_sim];
            Pooled.(fName).Seq_Probs   = [Pooled.(fName).Seq_Probs; prob_seq];
            
            Pooled.(fName).Sim_Global  = [Pooled.(fName).Sim_Global; glob_sim];
            Pooled.(fName).Seq_Global  = [Pooled.(fName).Seq_Global; glob_seq];
        end
    end
end

All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;

if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ================= 2. PLOT PROFILES (Per Amplitude) =================
fprintf('Generating Profiles for %d Amplitudes...\n', length(All_Amps));

% Arrays for Summary Plot
Sim_Rext = []; Seq_Rext = [];

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    % Extract Data
    Sim_P = Pooled.(fName).Sim_Probs;   Seq_P = Pooled.(fName).Seq_Probs;
    Sim_G = Pooled.(fName).Sim_Global;  Seq_G = Pooled.(fName).Seq_Global;
    
    if isempty(Sim_P), Sim_Rext(i)=NaN; Seq_Rext(i)=NaN; continue; end
    
    % --- Statistics (Mean & SEM) ---
    Avg_P_Sim = mean(Sim_P, 1, 'omitnan'); SEM_P_Sim = std(Sim_P, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Sim_P), 1));
    Avg_P_Seq = mean(Seq_P, 1, 'omitnan'); SEM_P_Seq = std(Seq_P, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Seq_P), 1));
    
    Avg_G_Sim = mean(Sim_G, 1, 'omitnan'); SEM_G_Sim = std(Sim_G, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Sim_G), 1));
    Avg_G_Seq = mean(Seq_G, 1, 'omitnan'); SEM_G_Seq = std(Seq_G, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Seq_G), 1));
    
    % --- Calc R_ext (Occupancy > 2%) ---
    threshold = 0.0035; 
    idx_sim = find(Avg_G_Sim > threshold, 1, 'last');
    if ~isempty(idx_sim), Sim_Rext(i) = bin_centers(idx_sim); else, Sim_Rext(i) = 0; end
    
    idx_seq = find(Avg_G_Seq > threshold, 1, 'last');
    if ~isempty(idx_seq), Seq_Rext(i) = bin_centers(idx_seq); else, Seq_Rext(i) = 0; end
    
    % === FIG 1: SPATIAL DENSITY (Prob within Bin) [B&W] ===
    figNameA = sprintf('SpatialDensity_Group_%.1fuA', curr_amp);
    figure('Color','w', 'Position', [100 100 600 400], 'Name', figNameA); hold on;
    
    b = bar(bin_centers, [Avg_P_Sim; Avg_P_Seq]', 'grouped');
    
    % --- STYLE: Black & White ---
    % Sim: White Face, Black Edge
    b(1).FaceColor = 'w'; b(1).EdgeColor = 'k'; b(1).LineWidth = 1.0; b(1).DisplayName = 'Simultaneous';
    % Seq: Black Face, No Edge
    b(2).FaceColor = 'k'; b(2).EdgeColor = 'none'; b(2).DisplayName = 'Sequential';
    
    % Error Bars (Black)
    % Error Bars (Black) - Use XEndPoints for exact centering
    errorbar(b(1).XEndPoints, Avg_P_Sim, SEM_P_Sim, 'k.', 'LineWidth', 1, 'HandleVisibility','off');
    errorbar(b(2).XEndPoints, Avg_P_Seq, SEM_P_Seq, 'k.', 'LineWidth', 1, 'HandleVisibility','off');

    xlabel('Distance (\mum)'); ylabel('Response Probability');
    ylim([0 1.05]); title(sprintf('Spatial Density @ %.1f \\muA', curr_amp));
    
    % --- [UPDATED] X-Axis Ticks Every 50um ---
    set(gca, 'XTick', 0 : 50 : max(dist_bin_edges));
    xlim([0, max(dist_bin_edges)]);
    
    legend('Location','best'); box off;
    % if save_figures, saveas(gcf, fullfile(save_dir, [figNameA '.png'])); end
    
    % === FIG 2: GLOBAL OCCUPANCY (Count / Total Array) [B&W] ===
    figNameB = sprintf('GlobalOccupancy_Group_%.1fuA', curr_amp);
    figure('Color','w', 'Position', [750 100 600 400], 'Name', figNameB); hold on;
    
    b = bar(bin_centers, [Avg_G_Sim; Avg_G_Seq]', 'grouped');
    
    % --- STYLE: Black & White ---
    b(1).FaceColor = 'w'; b(1).EdgeColor = 'k'; b(1).LineWidth = 1.0; b(1).DisplayName = 'Simultaneous';
    b(2).FaceColor = 'k'; b(2).EdgeColor = 'none'; b(2).DisplayName = 'Sequential';
    
    % Error Bars (Black) - Use XEndPoints for exact centering
    errorbar(b(1).XEndPoints, Avg_G_Sim, SEM_G_Sim, 'k.', 'LineWidth', 1, 'HandleVisibility','off');
    errorbar(b(2).XEndPoints, Avg_G_Seq, SEM_G_Seq, 'k.', 'LineWidth', 1, 'HandleVisibility','off');

    yline(threshold, '--k', 'Range Threshold', 'HandleVisibility','off');
    
    xlabel('Distance (\mum)'); ylabel('Fraction of Array Recruited');
    title(sprintf('Global Occupancy @ %.1f \\muA', curr_amp));
    
    % --- [UPDATED] X-Axis Ticks Every 50um ---
    set(gca, 'XTick', 0 : 50 : max(dist_bin_edges));
    xlim([0, max(dist_bin_edges)]);
    
    legend('Location','best'); box off;
    % if save_figures, saveas(gcf, fullfile(save_dir, [figNameB '.png'])); end
end

%% ================= 3. SUMMARY PLOT (Absolute Range) [B&W] =================
fprintf('Generating Summary R_ext Plot...\n');

figNameSum = 'Spatial_Summary_R_ext';
figure('Color','w', 'Position', [400 400 600 500], 'Name', figNameSum); hold on;

% --- STYLE: B&W Lines ---
% Sim: Dashed Line, Open Circle
plot(All_Amps, Sim_Rext, '--o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Simultaneous');

% Seq: Solid Line, Filled Square
plot(All_Amps, Seq_Rext, '-s', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', 'Sequential');

xlabel('Amplitude (µA)'); ylabel('Absolute Range (R_{ext}) [µm]');
title('Absolute Spatial Range');
legend('Location','best'); box off;

if save_figures, saveas(gcf, fullfile(save_dir, [figNameSum '.fig'])); end
fprintf('>>> Group Analysis Complete.\n');

%% ================= HELPER FUNCTIONS =================
function [metric_vec] = calc_metrics(dist, resp, valid, edges, mode)
    metric_vec = nan(1, length(edges)-1);
    counts     = zeros(1, length(edges)-1);
    
    for b = 1:length(edges)-1
        mask = dist >= edges(b) & dist < edges(b+1);
        mask_valid = mask & valid;
        if sum(mask_valid) > 0
            counts(b) = sum(resp(mask_valid));
        end
    end
    
    if strcmp(mode, 'prob')
        for b = 1:length(edges)-1
            mask = dist >= edges(b) & dist < edges(b+1);
            mask_valid = mask & valid;
            n_valid = sum(mask_valid);
            if n_valid > 0, metric_vec(b) = counts(b) / n_valid; else, metric_vec(b) = NaN; end
        end
    elseif strcmp(mode, 'global')
        n_total_valid = sum(valid);
        if n_total_valid > 0, metric_vec = counts / n_total_valid; else, metric_vec = zeros(size(counts)); end
    end
end