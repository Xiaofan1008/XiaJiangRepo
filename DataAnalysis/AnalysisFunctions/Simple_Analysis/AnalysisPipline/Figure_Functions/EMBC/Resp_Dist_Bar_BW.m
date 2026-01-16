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
save_dir     = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE_EMBC/Figures/6_Resp_Dist';
dist_bin_edges = 0 : 100 : 1000; 

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
    threshold = 0.004; 
    idx_sim = find(Avg_G_Sim > threshold, 1, 'last');
    if ~isempty(idx_sim), Sim_Rext(i) = bin_centers(idx_sim); else, Sim_Rext(i) = 0; end
    
    idx_seq = find(Avg_G_Seq > threshold, 1, 'last');
    if ~isempty(idx_seq), Seq_Rext(i) = bin_centers(idx_seq); else, Seq_Rext(i) = 0; end
    
    % === FIG 1: SPATIAL DENSITY (Prob within Bin) [B&W] ===
    figNameA = sprintf('Resp_Prob_%.1fuA', curr_amp);
    % figure('Color','w', 'Position', [100 100 600 400], 'Name', figNameA); hold on;
    figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5],'Name', figNameA, 'Color', 'w'); hold on;

    b = bar(bin_centers, [Avg_P_Sim; Avg_P_Seq]', 'grouped');
    
    % --- STYLE: Black & White ---
    % Sim: White Face, Black Edge
    % b(1).FaceColor = 'w'; 
    b(1).FaceColor = [0.8 0.8 0.8];
    b(1).EdgeColor = 'k'; b(1).LineWidth = 0.1; b(1).DisplayName = 'Simultaneous';
    % Seq: Black Face, No Edge
    b(2).FaceColor = 'k'; b(2).EdgeColor = 'none'; b(2).DisplayName = 'Sequential';
    
    % Error Bars (Black)
    % Error Bars (Black) - Use XEndPoints for exact centering
    errorbar(b(1).XEndPoints, Avg_P_Sim, SEM_P_Sim, 'k.', 'LineWidth', 0.5, 'HandleVisibility','off');
    errorbar(b(2).XEndPoints, Avg_P_Seq, SEM_P_Seq, 'k.', 'LineWidth', 0.5, 'HandleVisibility','off');

    xlabel('Distance (\mum)'); ylabel('Response Probability');
    ylim([0 1.05]); title(sprintf('Spatial Density @ %.1f \\muA', curr_amp));
    
    % --- X-Axis Ticks Every 50um ---
    set(gca, 'XTick', 0 : 50 : max(dist_bin_edges));
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 2);

    xlim([0, max(dist_bin_edges)]);
    
    legend('Location','best'); box off;
    % if save_figures, saveas(gcf, fullfile(save_dir, [figNameA '.fig'])); end
    
    % === FIG 2: GLOBAL OCCUPANCY (Count / Total Array) [B&W] ===
    figNameB = sprintf('Global_Resp_%.1fuA', curr_amp);
    % figure('Color','w', 'Position', [750 100 600 400], 'Name', figNameB); hold on;
    figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'Name', figNameB,'Color', 'w'); hold on;


    b = bar(bin_centers, [Avg_G_Sim; Avg_G_Seq]' *100, 'grouped');

    % --- STYLE: Black & White ---
    % b(1).FaceColor = 'w';
    b(1).FaceColor = [0.8 0.8 0.8];
    b(1).EdgeColor = 'k'; b(1).LineWidth = 0.1; b(1).DisplayName = 'Simultaneous';

    b(2).FaceColor = 'k'; 
    % b(2).EdgeColor = 'none'; 
    b(2).LineWidth = 0.1;b(2).DisplayName = 'Sequential';

    % Double side - Error Bars (Black) - Use XEndPoints for exact centering
    % errorbar(b(1).XEndPoints, Avg_G_Sim, SEM_G_Sim, 'k.', 'LineWidth', 1, 'HandleVisibility','off');
    % errorbar(b(2).XEndPoints, Avg_G_Seq, SEM_G_Seq, 'k.', 'LineWidth', 1, 'HandleVisibility','off');

    % Threshold line
    % yline(threshold, '--k', 'Range Threshold', 'HandleVisibility','off');
    % yline(threshold, '--k', 'HandleVisibility','off');

    % Single side error bar
    errorbar(b(1).XEndPoints, Avg_G_Sim * 100, zeros(size(SEM_G_Sim)), SEM_G_Sim * 100, ...
        'k.', 'LineWidth', 0.5, 'CapSize', 6, 'HandleVisibility','off');
        
    errorbar(b(2).XEndPoints, Avg_G_Seq * 100, zeros(size(SEM_G_Seq)), SEM_G_Seq * 100, ...
        'k.', 'LineWidth', 0.5, 'CapSize', 6, 'HandleVisibility','off');

    xlabel('Distance (\mum)','FontName', 'Times New Roman','FontSize', 14); 
    ylabel('Percentage of Total Channels (%)','FontName', 'Times New Roman','FontSize', 14);
    % title(sprintf('Global Occupancy @ %.1f \\muA', curr_amp));

    % --- X-Axis Ticks ---
    % set(gca, 'XTick', 0 : 50 : max(dist_bin_edges));
    set(gca, 'XTick', 0 : 100 : max(dist_bin_edges));
    set(gca, 'YTick', 0 : 4 : 20);
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.5);
    xlim([0, max(dist_bin_edges)]);
    ylim([0,20]);
    axis square;

    legend('Location','best','Box','off'); box off;
    if save_figures 
        % saveas(gcf, fullfile(save_dir, [figNameB '.fig']));
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
            exportgraphics(gcf, fullfile(save_dir, [figNameB '.tiff']), 'ContentType', 'vector');
    end
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

xlabel('Amplitude (\muA)','FontName', 'Times New Roman','FontSize', 14); 
ylabel('Absolute Range (R_{ext}) [\mum]','FontName', 'Times New Roman','FontSize', 14);
title('Absolute Spatial Range (Occupancy > 2%)');
legend('Location','best'); box off; 

if save_figures 
    % saveas(gcf, fullfile(save_dir, [figNameB '.fig']));
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        exportgraphics(gcf, fullfile(save_dir, [figNameSum '.tiff']), 'ContentType', 'vector');
end
fprintf('>>> Group Analysis Complete.\n');

%% ================= 4. SUMMARY PLOT: TOTAL RECRUITMENT vs. AMPLITUDE =================
fprintf('Generating Summary Plot: Total Recruitment (Sum of Occupancy)...\n');

% 1. Pre-allocate
Recruit_Sim_Mean = nan(size(All_Amps)); Recruit_Sim_SEM = nan(size(All_Amps));
Recruit_Seq_Mean = nan(size(All_Amps)); Recruit_Seq_SEM = nan(size(All_Amps));
Recruit_Pvals    = nan(size(All_Amps));

% 2. Calculate Stats
for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');

    if abs(curr_amp - 8.0) < 0.1  % Use tolerance for float comparison
        continue; 
    end
    
    % Retrieve the matrix of histograms: [N_Datasets x N_Bins]
    % Each row is one dataset's spatial histogram
    Sim_Hists = Pooled.(fName).Sim_Global;
    Seq_Hists = Pooled.(fName).Seq_Global;
    
    if isempty(Sim_Hists) || isempty(Seq_Hists), continue; end
    
    % --- Sum across Distance Bins to get TOTAL RECRUITMENT per dataset ---
    % Result is a vector: [N_Datasets x 1] (e.g., 0.15 means 15% of array active)
    Total_Sim_Vec = sum(Sim_Hists, 2); 
    Total_Seq_Vec = sum(Seq_Hists, 2);
    
    % --- Simultaneous Stats ---
    Recruit_Sim_Mean(i) = mean(Total_Sim_Vec, 'omitnan');
    Recruit_Sim_SEM(i)  = std(Total_Sim_Vec, 0, 'omitnan') / sqrt(sum(~isnan(Total_Sim_Vec)));
    
    % --- Sequential Stats ---
    Recruit_Seq_Mean(i) = mean(Total_Seq_Vec, 'omitnan');
    Recruit_Seq_SEM(i)  = std(Total_Seq_Vec, 0, 'omitnan') / sqrt(sum(~isnan(Total_Seq_Vec)));
    
    % --- Statistical Test (Rank Sum) ---
    p = ranksum(Total_Sim_Vec, Total_Seq_Vec);
    Recruit_Pvals(i) = p;
end

% 3. Plotting
% figNameRec = 'Summary_Total_Recruitment';
% figure('Color','w', 'Position', [500 400 600 500], 'Name', figNameRec); hold on;
% 
% % --- Plot Lines (B&W Style) ---
% % Sim: Dashed, Open Circle
% errorbar(All_Amps, Recruit_Sim_Mean, Recruit_Sim_SEM, '--o', ...
%     'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
%     'MarkerFaceColor', 'w', 'CapSize', 10, 'DisplayName', 'Simultaneous');
% 
% % Seq: Solid, Filled Square
% errorbar(All_Amps, Recruit_Seq_Mean, Recruit_Seq_SEM, '-s', ...
%     'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
%     'MarkerFaceColor', 'k', 'CapSize', 10, 'DisplayName', 'Sequential');
% 
% % --- Add Significance Stars ---
% for i = 1:length(All_Amps)
%     p = Recruit_Pvals(i);
%     if isnan(p), continue; end
% 
%     txt = ''; 
%     if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
% 
%     if ~isempty(txt)
%         % Place star slightly above the higher error bar
%         y_top = max(Recruit_Sim_Mean(i)+Recruit_Sim_SEM(i), Recruit_Seq_Mean(i)+Recruit_Seq_SEM(i));
%         text(All_Amps(i), y_top * 1.05, txt, 'FontSize', 20, ...
%             'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
%     end
% end
% 
% % --- Formatting ---
% xlabel('Amplitude (µA)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
% ylabel('Total Active Fraction (Area)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
% % title('Total Recruitment vs. Amplitude', 'FontSize', 14, 'FontName', 'Times New Roman');
% 
% legend('Location','northwest', 'Box','off', 'FontSize', 16, 'FontName', 'Times New Roman');
% box off; 
% set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 2);
% 
% % xlim([1-0.5, max(All_Amps)+0.5]);
% % xticks(1:1:max(All_Amps));
% xlim([1-0.5, 10+0.5]);
% xticks(1:1:10);
% 
% axis square;
% % Adjust Y-Limit (0 to Max + padding)
% ymax = max([Recruit_Sim_Mean + Recruit_Sim_SEM, Recruit_Seq_Mean + Recruit_Seq_SEM], [], 'all');
% % ylim([0, ymax * 1.2]); 
% ylim([0, 1]);
% 
% % --- Save ---
% if save_figures
%     exportgraphics(gcf, fullfile(save_dir, [figNameRec '.tiff']), 'ContentType', 'vector');
%     fprintf('Summary Recruitment plot saved: %s\n', figNameRec);
% end

% ================= 3. PLOT 3: TOTAL RECRUITMENT (Line Plot) =================
figNameRec = 'Summary_Total_Recruitment_v2';
% figure('Color','w', 'Position', [500 400 600 500], 'Name', figNameRec); hold on;
figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'Name', figNameRec,'Color', 'w'); hold on;


% --- FILTER OUT NaNs BEFORE PLOTTING ---
% Find indices where we actually calculated data (excludes the skipped 8.0uA)
valid_mask = ~isnan(Recruit_Sim_Mean);

% Extract only the valid data points for plotting
Plot_Amps     = All_Amps(valid_mask);
Plot_Sim_Mean = Recruit_Sim_Mean(valid_mask); 
Plot_Sim_SEM  = Recruit_Sim_SEM(valid_mask);
Plot_Seq_Mean = Recruit_Seq_Mean(valid_mask); 
Plot_Seq_SEM  = Recruit_Seq_SEM(valid_mask);

% Amp = 0 case
Plot_Amps     = [0, Plot_Amps];
Plot_Sim_Mean = [0, Plot_Sim_Mean]; Plot_Sim_SEM = [0, Plot_Sim_SEM];
Plot_Seq_Mean = [0, Plot_Seq_Mean]; Plot_Seq_SEM = [0, Plot_Seq_SEM];

% --- Plot Lines (Using the Filtered Data) ---
% Sim: Dashed, Open Circle
errorbar(Plot_Amps, Plot_Sim_Mean, Plot_Sim_SEM, '--o', ...
    'Color', 'k', 'LineWidth', 1, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'w', 'CapSize', 8, 'DisplayName', 'Simultaneous');

% Seq: Solid, Filled Square
errorbar(Plot_Amps, Plot_Seq_Mean, Plot_Seq_SEM, '-s', ...
    'Color', 'k', 'LineWidth', 1, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k', 'CapSize', 8, 'DisplayName', 'Sequential');

% --- Add Significance Stars ---
% for i = 1:length(All_Amps)
%     p = Recruit_Pvals(i);
%     % Check if p is valid (not NaN) before plotting
%     if isnan(p), continue; end
% 
%     txt = ''; 
%     if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
% 
%     if ~isempty(txt)
%         % Use the original full vectors here because 'i' corresponds to All_Amps
%         y_top = max(Recruit_Sim_Mean(i)+Recruit_Sim_SEM(i), Recruit_Seq_Mean(i)+Recruit_Seq_SEM(i));
%         text(All_Amps(i), y_top * 1.05, txt, 'FontSize', 20, ...
%             'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
%     end
% end

% --- Formatting ---
xlabel('Amplitude (µA)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Fraction of Active Channels', 'FontSize', 14, 'FontName', 'Times New Roman');

legend('Location','northwest', 'Box','off', 'FontSize', 14, 'FontName', 'Times New Roman');
box off; 
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.5);
set(gca, 'YTick', 0 : 0.2 : 1);
% X-Limits: Set to cover the full range (e.g., 1 to 10), preserving the gap visually on the axis
xlim([0, 10]);
xticks(0:2:10);
axis square;
ylim([0, 1]);

% --- Save ---
if save_figures
    exportgraphics(gcf, fullfile(save_dir, [figNameRec '.tiff']), 'ContentType', 'vector');
    fprintf('Summary Recruitment plot saved: %s\n', figNameRec);
end


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