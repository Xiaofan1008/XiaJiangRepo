%% ============================================================
%   GRAND AVERAGE: NORMALIZED RECRUITMENT CURVE (Discrete Bins)
%   - Logic: Aggregates pre-calculated normalized results.
%   - Grouping: Bins data by exact amplitude (No interpolation).
%   - Statistics: Mean +/- SEM across DATASETS.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));
%% ================= 1. USER SETTINGS (PASTE FILES HERE) =================
% List all result files to include in the Grand Average
file_paths = {
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
% Plot Settings
save_figure = false;
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE EMBC/Figures/3_Total_Spike_Count';
fig_name    = 'GrandAverage_Recruitment_Curve_BarSEM.tiff';
%% ================= 2. AGGREGATE DATA =================
fprintf('Processing %d datasets...\n', length(file_paths));
% Storage: [Amplitude, Mean_Value]
Pool_Sim = []; 
Pool_Seq = [];
for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file')
        warning('File not found: %s', file_paths{i}); continue; 
    end
    
    D = load(file_paths{i});
    if ~isfield(D, 'ResultNorm'), continue; end
    R = D.ResultNorm;
    
    Amps = R.Amps; % Amplitudes in this dataset
    
    % --- 1. Collapse Channels & Sets to get Dataset Mean ---
    % Input: [Channels x Amps x Sets]
    
    % Simultaneous
    % Mean across Sets (Dim 3), then Mean across Channels (Dim 1)
    sim_dataset_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));
    
    % Sequential (Force PTD=5ms check if needed, usually stored correctly)
    % Mean across Sets (Dim 3), then Mean across Channels (Dim 1)
    % Note: If dim 4 exists (PTDs), squeeze handles it if we filtered correctly before.
    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        % If 4D, average PTD dimension (Dim 4)
        seq_dataset_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_dataset_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end
    
    % --- 2. Store in Master Pool ---
    % We store simply as pairs: [Amp, Value]
    for a = 1:length(Amps)
        % Sim
        if ~isnan(sim_dataset_mean(a))
            Pool_Sim = [Pool_Sim; Amps(a), sim_dataset_mean(a)]; %#ok<*AGROW>
        end
        % Seq
        if ~isnan(seq_dataset_mean(a))
            Pool_Seq = [Pool_Seq; Amps(a), seq_dataset_mean(a)];
        end
    end
end
%% ================= 3. CALCULATE GROUP STATISTICS =================
% Find all unique amplitudes that exist in any dataset
Unique_Amps = unique([Pool_Sim(:,1); Pool_Seq(:,1)]);
Unique_Amps = sort(Unique_Amps);
% Init Output Vectors
Grand_Sim_Mean = nan(size(Unique_Amps)); Grand_Sim_SEM = nan(size(Unique_Amps));
Grand_Seq_Mean = nan(size(Unique_Amps)); Grand_Seq_SEM = nan(size(Unique_Amps));
N_Sim = zeros(size(Unique_Amps)); N_Seq = zeros(size(Unique_Amps));
fprintf('\n=== GRAND AVERAGE STATISTICS (N = Datasets) ===\n');
fprintf('Amp(uA)\tN_Sim\tMean_Sim\tSEM_Sim\t\tN_Seq\tMean_Seq\tSEM_Seq\n');
for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    
    % --- Sim ---
    vals = Pool_Sim(Pool_Sim(:,1) == amp, 2);
    % n = length(vals);
    n = 9;
    if n > 0
        Grand_Sim_Mean(k) = mean(vals);
        Grand_Sim_SEM(k)  = std(vals) / sqrt(n);
        % Grand_Sim_SEM(k)  = std(vals);
        N_Sim(k) = n;
    end
    
    % --- Seq ---
    vals = Pool_Seq(Pool_Seq(:,1) == amp, 2);
    % n = length(vals);
    n = 9;
    if n > 0
        Grand_Seq_Mean(k) = mean(vals);
        Grand_Seq_SEM(k)  = std(vals) / sqrt(n);
        % Grand_Seq_SEM(k)  = std(vals);
        N_Seq(k) = n;
    end
    
    fprintf('%.1f\t%d\t%.4f\t%.4f\t\t%d\t%.4f\t%.4f\n', ...
        amp, N_Sim(k), Grand_Sim_Mean(k), Grand_Sim_SEM(k), N_Seq(k), Grand_Seq_Mean(k), Grand_Seq_SEM(k));
end
%% ================= 4. PLOT GRAND AVERAGE =================
figure('Color','w', 'Position',[100 100 700 500]); hold on;
% % --- Plot Sim with Error Bars ---
% valid = ~isnan(Grand_Sim_Mean);
% errorbar(Unique_Amps(valid), Grand_Sim_Mean(valid), Grand_Sim_SEM(valid), '-o', ...
%     'Color', [0 0.3 0.8], 'LineWidth', 2, 'MarkerFaceColor', [0 0.3 0.8], ...
%     'DisplayName', 'Simultaneous', 'CapSize', 10);
% 
% % --- Plot Seq with Error Bars ---
% valid = ~isnan(Grand_Seq_Mean);
% errorbar(Unique_Amps(valid), Grand_Seq_Mean(valid), Grand_Seq_SEM(valid), '-s', ...
%     'Color', [0.7 0 0], 'LineWidth', 2, 'MarkerFaceColor', 'w', ...
%     'DisplayName', 'Sequential', 'CapSize', 10);
% % --- Plot Sim Thick & Thin Bar ---
% valid = ~isnan(Grand_Sim_Mean);
% % 1. Thin Error Bars (LineWidth = 1)
% errorbar(Unique_Amps(valid), Grand_Sim_Mean(valid), Grand_Sim_SEM(valid), '.', ...
%     'Color', [0 0.3 0.8], 'LineWidth', 1, 'CapSize', 10, 'HandleVisibility', 'off');
% % 2. Thick Main Line (LineWidth = 2)
% plot(Unique_Amps(valid), Grand_Sim_Mean(valid), '-o', ...
%     'Color', [0 0.3 0.8], 'LineWidth', 2, 'MarkerFaceColor', [0 0.3 0.8], ...
%     'DisplayName', 'Simultaneous');
% 
% % --- Plot Seq ---
% valid = ~isnan(Grand_Seq_Mean);
% % 1. Thin Error Bars (LineWidth = 1)
% errorbar(Unique_Amps(valid), Grand_Seq_Mean(valid), Grand_Seq_SEM(valid), '.', ...
%     'Color', [0.7 0 0], 'LineWidth', 1, 'CapSize', 10, 'HandleVisibility', 'off');
% % 2. Thick Main Line (LineWidth = 2)
% plot(Unique_Amps(valid), Grand_Seq_Mean(valid), '-s', ...
%     'Color', [0.7 0 0], 'LineWidth', 2, 'MarkerFaceColor', 'w', ...
%     'DisplayName', 'Sequential');
%% ================= 4. PLOT GRAND AVERAGE (B&W) =================
% figure('Color','w', 'Position',[100 100 700 500]); hold on;
% % --- Plot Sim (Dashed Line, Open Markers) ---
% valid = ~isnan(Grand_Sim_Mean);
% % 1. Error Bars (Black)
% errorbar(Unique_Amps(valid), Grand_Sim_Mean(valid), Grand_Sim_SEM(valid), '.', ...
%     'Color', 'k', 'LineWidth', 1, 'CapSize', 10, 'HandleVisibility', 'off');
% % 2. Main Line (Black, Dashed '--', White Face 'w')
% plot(Unique_Amps(valid), Grand_Sim_Mean(valid), '--o', ...
%     'Color', 'k', 'LineWidth', 3, 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
%     'DisplayName', 'Simultaneous');
% % --- Plot Seq (Solid Line, Filled Markers) ---
% valid = ~isnan(Grand_Seq_Mean);
% % 1. Error Bars (Black)
% errorbar(Unique_Amps(valid), Grand_Seq_Mean(valid), Grand_Seq_SEM(valid), '.', ...
%     'Color', 'k', 'LineWidth', 1, 'CapSize', 10, 'HandleVisibility', 'off');
% % 2. Main Line (Black, Solid '-', Black Face 'k')
% plot(Unique_Amps(valid), Grand_Seq_Mean(valid), '-s', ...
%     'Color', 'k', 'LineWidth', 3, 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
%     'DisplayName', 'Sequential');
% % --- Formatting ---
% % yline(1.0, '--k', 'Ref (5uA)', 'HandleVisibility','off');
% % yline(1.0, '--k', 'HandleVisibility','off');
% 
% xlabel('Amplitude (\muA)', 'FontSize', 14, 'FontWeight','bold', 'FontName', 'Times New Roman');
% ylabel('Normalized Spike Count (a.u.)', 'FontSize', 14, 'FontWeight','bold', 'FontName', 'Times New Roman');
% % title(sprintf('Spike Count vs. AMP (N=%d)', length(file_paths)), 'FontSize', 12, 'FontName', 'Times New Roman');
% 
% % legend('Location','best','Box','off', 'FontName', 'Times New Roman');
% legend('Location','northwest', 'Box','off', 'FontSize', 14, 'FontName', 'Times New Roman');
% box off; 
% 
% set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
% % xlim([min(Unique_Amps)-0.5, max(Unique_Amps)+0.5]);
% xlim([1, max(Unique_Amps)]);
% ylim([0 2]);
% % --- Save ---
% if save_figure
%     if ~exist(save_dir, 'dir'), mkdir(save_dir); end
%     saveas(gcf, fullfile(save_dir, fig_name));
%     fprintf('\nFigure saved to: %s\n', fullfile(save_dir, fig_name));
% end

%% ================= 4. PLOT GRAND AVERAGE (B&W) =================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% --- A. Plot Individual Data Points (Background Scatter) ---
% Jitter width helps separate dots at the same amplitude
jitter_w = 0.4; 

% 1. Sequential Scatter (Dark Grey Dots)
for i = 1:size(Pool_Seq, 1)
    amp = Pool_Seq(i,1); val = Pool_Seq(i,2);
    x_jit = amp + (rand()-0.5) * jitter_w;
    % scatter(x_jit, val, 30, [0.7 0.7 0.7], 's','filled', 'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');
    scatter(x_jit, val, 30, [0.7 0.7 0.7], 's', 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');

end

% 2. Simultaneous Scatter (Light Grey Dots)
for i = 1:size(Pool_Sim, 1)
    amp = Pool_Sim(i,1); val = Pool_Sim(i,2);
    x_jit = amp + (rand()-0.5) * jitter_w;
    % scatter(x_jit, val, 30, [0.8 0.8 0.8], 'o','filled', 'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');
    scatter(x_jit, val, 30, [0.8 0.8 0.8], 'o', 'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');

end

% --- B. Plot Main Curves (Foreground) ---
% --- Plot Sim (Dashed Line, Open Markers) ---
valid = ~isnan(Grand_Sim_Mean);
% 1. Error Bars (Black)
errorbar(Unique_Amps(valid), Grand_Sim_Mean(valid), Grand_Sim_SEM(valid), '.', ...
    'Color', 'k', 'LineWidth', 2, 'CapSize', 12, 'HandleVisibility', 'off');
% 2. Main Line (Black, Dashed '--', White Face 'w')
plot(Unique_Amps(valid), Grand_Sim_Mean(valid), '--o', ...
    'Color', 'k', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerFaceColor', 'w', ...
    'DisplayName', 'Simultaneous');

% --- Plot Seq (Solid Line, Filled Markers) ---
valid = ~isnan(Grand_Seq_Mean);
% 1. Error Bars (Black)
errorbar(Unique_Amps(valid), Grand_Seq_Mean(valid), Grand_Seq_SEM(valid), '.', ...
    'Color', 'k', 'LineWidth', 2, 'CapSize', 12, 'HandleVisibility', 'off');
% 2. Main Line (Black, Solid '-', Black Face 'k')
plot(Unique_Amps(valid), Grand_Seq_Mean(valid), '-s', ...
    'Color', 'k', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Sequential');

% --- Formatting ---
% yline(1.0, '--k', 'Ref (5uA)', 'HandleVisibility','off');
xlabel('Amplitude (\muA)', 'FontSize', 20, 'FontWeight','bold', 'FontName', 'Times New Roman');
ylabel('Normalized Spike Count (a.u.)', 'FontSize', 20, 'FontWeight','bold', 'FontName', 'Times New Roman');

legend('Location','northwest', 'Box','off', 'FontSize', 18, 'FontName', 'Times New Roman');

box off; 
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 3);
axis square;
xlim([1, max(Unique_Amps)]);
ylim([0 2]); % Fixed Y-limit

% --- Save ---
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    exportgraphics(gcf, fullfile(save_dir, fig_name),'ContentType', 'vector');
    fprintf('\nFigure saved to: %s\n', fullfile(save_dir, fig_name));
end



%% ================= 5. PLOT NAKA-RUSHTON FIT =================
figure('Color','w', 'Position',[150 150 700 500]); hold on;
% High-res X axis for smooth curve
x_smooth = linspace(min(Unique_Amps), max(Unique_Amps), 100);

% --- A. Fit & Plot Simultaneous ---
valid_sim = ~isnan(Grand_Sim_Mean);
if any(valid_sim)
    x_dat = Unique_Amps(valid_sim); y_dat = Grand_Sim_Mean(valid_sim);
    [p_sim, y_fit_sim] = fit_naka_rushton(x_dat, y_dat, x_smooth);
    
    % Plot Data Points (Light Grey)
    errorbar(x_dat, y_dat, Grand_Sim_SEM(valid_sim), 'o', ...
        'Color', [0.6 0.6 0.6], 'MarkerFaceColor', 'w', ...
        'LineWidth', 1, 'CapSize', 0, 'HandleVisibility', 'off');
    
    % Plot Fitted Curve (Black Dashed)
    plot(x_smooth, y_fit_sim, '--k', 'LineWidth', 2, 'DisplayName', 'Simultaneous (Fit)');
end

% --- B. Fit & Plot Sequential ---
valid_seq = ~isnan(Grand_Seq_Mean);
if any(valid_seq)
    x_dat = Unique_Amps(valid_seq); y_dat = Grand_Seq_Mean(valid_seq);
    [p_seq, y_fit_seq] = fit_naka_rushton(x_dat, y_dat, x_smooth);
    
    % Plot Data Points (Dark Grey)
    errorbar(x_dat, y_dat, Grand_Seq_SEM(valid_seq), 's', ...
        'Color', [0.4 0.4 0.4], 'MarkerFaceColor', [0.4 0.4 0.4], ...
        'LineWidth', 1, 'CapSize', 0, 'HandleVisibility', 'off');
    
    % Plot Fitted Curve (Black Solid)
    plot(x_smooth, y_fit_seq, '-k', 'LineWidth', 2, 'DisplayName', 'Sequential (Fit)');
end

% --- C. Formatting (Times New Roman Style) ---
yline(1.0, ':k', 'HandleVisibility','off');
xlabel('Amplitude (\muA)', 'FontSize', 12, 'FontWeight','bold', 'FontName', 'Times New Roman');
ylabel('Normalized Spike Count', 'FontSize', 12, 'FontWeight','bold', 'FontName', 'Times New Roman');
title('Naka-Rushton Fit', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location','northwest','Box','off', 'FontName', 'Times New Roman');
box off; set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
xlim([1, max(Unique_Amps)]);
%% ================= LOCAL FUNCTIONS =================
function plot_shaded(x, y, err, col)
    % Filters NaNs and plots filled polygon
    valid = ~isnan(y) & ~isnan(err);
    if sum(valid) < 2, return; end
    
    x = x(valid); y = y(valid); err = err(valid);
    
    % Create polygon points
    x_poly = [x; flipud(x)];
    y_poly = [y + err; flipud(y - err)];
    
    fill(x_poly, y_poly, col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function [params, y_fit] = fit_naka_rushton(x, y, x_query)
    % Fits Naka-Rushton Equation: R(c) = Rmax * (c^n / (c^n + c50^n)) + M
    % Params = [Rmax, c50, n, M]
    
    x = x(:); y = y(:); % Ensure column vectors
    
    % Initial Guesses
    M0 = min(y);
    Rmax0 = max(y) - M0;
    c500 = mean(x);
    n0 = 2; 
    
    start_params = [Rmax0, c500, n0, M0];
    
    % Objective Function: Minimize Sum of Squared Errors
    % Force params to be real via abs() if unconstrained optimization pushes them negative
    cost_func = @(p) sum((y - (p(4) + p(1) .* (x.^p(3) ./ (x.^p(3) + p(2).^p(3))))) .^ 2);
    
    % Optimization
    options = optimset('Display', 'off', 'MaxFunEvals', 1000);
    params = fminsearch(cost_func, start_params, options);
    
    % Generate smooth curve for plotting
    y_fit = params(4) + params(1) .* (x_query.^params(3) ./ (x_query.^params(3) + params(2).^params(3)));
end