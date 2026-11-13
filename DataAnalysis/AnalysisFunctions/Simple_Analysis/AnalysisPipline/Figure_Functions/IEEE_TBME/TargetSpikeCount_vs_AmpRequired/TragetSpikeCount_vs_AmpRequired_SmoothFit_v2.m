%% ============================================================
%   ISO-RESPONSE PLOT: REQUIRED AMPLITUDE FOR TARGET OUTPUT
%   - Logic: Cubic Smoothing Spline (csaps) per dataset.
%   - Style: IEEE TBME Publication Style (Grayscale, Percentage Axis)
%   - Statistics: Paired T-Test with Command Window Printout
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% --- ISO-RESPONSE SETTINGS ---
% Internal math targets (0.2 to 2.0). These will be converted to % for plotting.
Target_Spikes = 0.2 : 0.2 : 2.0; 

% [NEW] SMOOTHING SPLINE PARAMETER
% 1.0 = Connect the dots exactly (Jagged)
% 0.0 = Perfectly straight line (Linear Regression)
% 0.8 = The "Sweet Spot" (Smooth curve that follows the data trend)
spline_smoothness = 0.5; 

% --- STATISTICAL FILTERS ---
stats_min_n_threshold = 1; 

% List all result files
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_4.mat';
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
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure2/IsoResponse';
fig_name    = 'IsoResponse_Required_Amplitude_v1.tiff';

%% ================= 2. AGGREGATE DATA (Harvesting) =================
fprintf('Harvesting recruitment data from %d datasets...\n', length(file_paths));
Pool_Flipped_Raw = []; % Stores [SpikeCount, Amplitude, StrategyID (1=Sim, 2=Seq)]

% We need to keep tracks of individual datasets for paired fitting
Dataset_Amps = {};
Dataset_SimVals = {};
Dataset_SeqVals = {};

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    D = load(file_paths{i}); if ~isfield(D, 'ResultNorm'), continue; end
    R = D.ResultNorm; Amps = R.Amps; 
    
    sim_ds = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));
    
    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_ds = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_ds = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end
    
    Dataset_Amps{end+1} = Amps;
    Dataset_SimVals{end+1} = sim_ds;
    Dataset_SeqVals{end+1} = seq_ds;
    
    % Store for the background scatter plot (Flipped)
    for a = 1:length(Amps)
        if ~isnan(sim_ds(a)), Pool_Flipped_Raw = [Pool_Flipped_Raw; sim_ds(a), Amps(a), 1]; end
        if ~isnan(seq_ds(a)), Pool_Flipped_Raw = [Pool_Flipped_Raw; seq_ds(a), Amps(a), 2]; end
    end
end

%% ================= 3. SMOOTHING SPLINE ENGINE (csaps) =================
fprintf('Fitting Cubic Smoothing Splines (p=%.2f) for target spike levels...\n', spline_smoothness);

% Results holders: [TargetIndex, DatasetIndex]
Cost_Sim = nan(length(Target_Spikes), length(Dataset_Amps));
Cost_Seq = nan(length(Target_Spikes), length(Dataset_Amps));

for d = 1:length(Dataset_Amps)
    curr_amps = Dataset_Amps{d}(:).';
    curr_sim  = Dataset_SimVals{d}(:).';
    curr_seq  = Dataset_SeqVals{d}(:).';
    
    [u_sim, idx_s] = unique(curr_sim);
    [u_seq, idx_q] = unique(curr_seq);
    
    % Fit Simultaneous Curve (Requires at least 2 points for a spline)
    if length(u_sim) >= 2
        % [MODIFIED] csaps evaluates the smoothing spline at Target_Spikes
        Cost_Sim(:, d) = csaps(u_sim, curr_amps(idx_s), spline_smoothness, Target_Spikes);
        
        % Safety: Do not extrapolate beyond the highest recorded spike count
        Cost_Sim(Target_Spikes > max(u_sim), d) = NaN;
    end
    
    % Fit Sequential Curve
    if length(u_seq) >= 2
        Cost_Seq(:, d) = csaps(u_seq, curr_amps(idx_q), spline_smoothness, Target_Spikes);
        
        % Safety: Do not extrapolate beyond the highest recorded spike count
        Cost_Seq(Target_Spikes > max(u_seq), d) = NaN;
    end
end

% Calculate Grand Summary Statistics across Datasets
Grand_Sim_Cost_Mean = mean(Cost_Sim, 2, 'omitnan');
Grand_Sim_Cost_SEM  = std(Cost_Sim, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(Cost_Sim), 2));

Grand_Seq_Cost_Mean = mean(Cost_Seq, 2, 'omitnan');
Grand_Seq_Cost_SEM  = std(Cost_Seq, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(Cost_Seq), 2));

%% ================= 4. PLOT (Flipped Iso-Response with Percentages) =================
figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w'); hold on;

% Convert Math targets to Percentages for Plotting
Target_Percents = Target_Spikes * 100;

% A. Background Scatter (Flipped Axes & Converted to Percentage)
jitter_x = 2; % Increased jitter slightly because the scale is now 0-200
raw_sim = Pool_Flipped_Raw(Pool_Flipped_Raw(:,3)==1, :);
raw_seq = Pool_Flipped_Raw(Pool_Flipped_Raw(:,3)==2, :);

scatter((raw_seq(:,1)*100) + (rand(size(raw_seq,1),1)-0.5)*jitter_x, raw_seq(:,2), 10, [0.7 0.7 0.7], 's', 'filled', 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
scatter((raw_sim(:,1)*100) + (rand(size(raw_sim,1),1)-0.5)*jitter_x, raw_sim(:,2), 10, [0.8 0.8 0.8], 'o', 'MarkerEdgeAlpha', 0.3, 'HandleVisibility', 'off');

% B. Main Cost Curves (Origin connection included)
Plot_X = [0, Target_Percents];
Plot_Sim_Y = [0, Grand_Sim_Cost_Mean'];
Plot_Seq_Y = [0, Grand_Seq_Cost_Mean'];

errorbar(Target_Percents, Grand_Sim_Cost_Mean, Grand_Sim_Cost_SEM, '.', 'Color', 'k', 'LineWidth', 1, 'CapSize', 4, 'HandleVisibility', 'off');
p1 = plot(Plot_X, Plot_Sim_Y, '--o', 'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');

errorbar(Target_Percents, Grand_Seq_Cost_Mean, Grand_Seq_Cost_SEM, '.', 'Color', 'k', 'LineWidth', 1, 'CapSize', 4, 'HandleVisibility', 'off');
p2 = plot(Plot_X, Plot_Seq_Y, '-s', 'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

% C. Statistics (Paired T-Test at each Target Level)
fprintf('\n=== ISO-RESPONSE STATS (Paired T-Test) ===\n');
for t = 1:length(Target_Spikes)
    data_s = Cost_Sim(t, :);
    data_q = Cost_Seq(t, :);
    
    % Clean NaNs for this specific level
    valid = ~isnan(data_s) & ~isnan(data_q);
    n_pairs = sum(valid);
    
    if n_pairs >= stats_min_n_threshold
        [~, p] = ttest(data_s(valid), data_q(valid));
        
        txt = '';
        if p < (0.001), txt = '***';
        elseif p < (0.01), txt = '**';
        elseif p < (0.05), txt = '*';
        end
        
        % Print results directly to the Command Window
        fprintf('Target %3.0f%% | N = %2d pairs | p = %.5f %s\n', Target_Percents(t), n_pairs, p, txt);
        
        if ~isempty(txt)
            y_mid = (Grand_Sim_Cost_Mean(t) + Grand_Seq_Cost_Mean(t)) / 2;
            text(Target_Percents(t), y_mid, txt, 'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', 'Arial');
        end
    end
end

% --- Formatting (IEEE TBME Standards) ---
box off; set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.2);
axis square;
xlabel('Relative Neural Recruitment (%)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Required Amplitude (\muA)', 'FontSize', 9, 'FontName', 'Arial');
legend([p1, p2], 'Location', 'northwest', 'Box', 'off', 'FontSize', 8, 'FontName', 'Arial');

xlim([0 200]); 
ylim([0 12]); 
set(gca, 'YTick', 0:2:12);
set(gca, 'XTick', 0:20:200);

if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    print(gcf, fullfile(save_dir, fig_name), '-dtiff', '-r300');
end