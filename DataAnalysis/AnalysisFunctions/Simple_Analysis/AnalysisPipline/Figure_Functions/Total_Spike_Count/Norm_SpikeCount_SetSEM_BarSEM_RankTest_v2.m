%% ============================================================
%   GRAND AVERAGE: NORMALIZED RECRUITMENT CURVE (Discrete Bins)
%   - Logic: Aggregates pre-calculated normalized results.
%   - Grouping: Bins data by exact amplitude (No interpolation).
%   - Statistics: Mean +/- SEM across DATASETS.
%   - NEW: Added Wilcoxon Rank Sum Test & Significance Stars
%   - FILTER: Bonferroni + Magnitude Threshold + Phys Threshold
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% --- STATISTICAL FILTERS ---
% 1. Physiological Threshold: Ignore amps below this (e.g., < 1.5 uA)
stats_phys_threshold = 2; 

% 2. Magnitude Threshold: Ignore if difference (Seq - Sim) is small
%    Suggestion: 0.2 (Filters out 1.0uA where diff is ~0.19)
stats_mag_threshold  = 0.35;  

% 3. Bonferroni Correction: Divide p-value cutoffs by this number
%    (Usually equal to number of comparisons, e.g., 9 amplitudes)
bonferroni_n = 9; 

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
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE_EMBC/Figures/3_Total_Spike_Count';
fig_name    = 'GrandAverage_Recruitment_Curve_StatsFiltered_v2.tiff';

%% ================= 2. AGGREGATE DATA =================
fprintf('Processing %d datasets...\n', length(file_paths));
Pool_Sim = []; Pool_Seq = [];

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    D = load(file_paths{i}); if ~isfield(D, 'ResultNorm'), continue; end
    R = D.ResultNorm; Amps = R.Amps; 
    
    sim_ds_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));
    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_ds_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_ds_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end
    
    % Store with 0-force logic
    for a = 1:length(Amps)
        if abs(Amps(a)) < 0.001, v_sim=0; v_seq=0;
        else, v_sim=sim_ds_mean(a); v_seq=seq_ds_mean(a); end
        
        if ~isnan(v_sim), Pool_Sim = [Pool_Sim; Amps(a), v_sim]; end
        if ~isnan(v_seq), Pool_Seq = [Pool_Seq; Amps(a), v_seq]; end
    end
end

%% ================= 3. STATS & VECTORS =================
Unique_Amps = unique([Pool_Sim(:,1); Pool_Seq(:,1)]);
Unique_Amps = sort(Unique_Amps);
if Unique_Amps(1) ~= 0, Unique_Amps = [0; Unique_Amps]; end

Grand_Sim_Mean = []; Grand_Sim_SEM = [];
Grand_Seq_Mean = []; Grand_Seq_SEM = [];

for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    
    vs = Pool_Sim(Pool_Sim(:,1)==amp, 2);
    if isempty(vs) && amp==0, vs=0; end
    Grand_Sim_Mean(k) = mean(vs); Grand_Sim_SEM(k) = std(vs)/sqrt(9);
    
    vq = Pool_Seq(Pool_Seq(:,1)==amp, 2);
    if isempty(vq) && amp==0, vq=0; end
    Grand_Seq_Mean(k) = mean(vq); Grand_Seq_SEM(k) = std(vq)/sqrt(9);
end

%% ================= 4. PLOT =================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% A. Scatter (Background)
jitter_w = 0.4; 
for i = 1:size(Pool_Seq, 1)
    scatter(Pool_Seq(i,1)+(rand-0.5)*jitter_w, Pool_Seq(i,2), 30, [0.7 0.7 0.7], 's', 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
for i = 1:size(Pool_Sim, 1)
    scatter(Pool_Sim(i,1)+(rand-0.5)*jitter_w, Pool_Sim(i,2), 30, [0.8 0.8 0.8], 'o', 'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end

% B. Main Curves
errorbar(Unique_Amps, Grand_Sim_Mean, Grand_Sim_SEM, '.', 'Color', 'k', 'LineWidth', 2, 'CapSize', 12, 'HandleVisibility', 'off');
p1 = plot(Unique_Amps, Grand_Sim_Mean, '--o', 'Color', 'k', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');

errorbar(Unique_Amps, Grand_Seq_Mean, Grand_Seq_SEM, '.', 'Color', 'k', 'LineWidth', 2, 'CapSize', 12, 'HandleVisibility', 'off');
p2 = plot(Unique_Amps, Grand_Seq_Mean, '-s', 'Color', 'k', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

% --- 4b. STATS with TRIPLE FILTER (Bonferroni + Phys + Magnitude) ---
fprintf('\n=== STATS FILTERING ===\n');
fprintf('1. Phys Threshold > %.1f uA\n', stats_phys_threshold);
fprintf('2. Mag Threshold > %.2f (Seq-Sim)\n', stats_mag_threshold);
fprintf('3. Bonferroni N = %d (Stricter P-values)\n', bonferroni_n);

for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    if amp == 0, continue; end
    
    % Get Data
    data_s = Pool_Sim(abs(Pool_Sim(:,1)-amp)<0.001, 2);
    data_q = Pool_Seq(abs(Pool_Seq(:,1)-amp)<0.001, 2);
    
    if isempty(data_s) || isempty(data_q), continue; end
    
    % --- FILTER 1: Physiological Threshold ---
    if amp <= stats_phys_threshold
        fprintf('Amp %.1f: Skipped (Below Phys Threshold)\n', amp);
        continue; 
    end
    
    % --- FILTER 2: Magnitude Threshold ---
    % mean_diff = abs(mean(data_q) - mean(data_s));
    mean_diff = abs(mean(data_q) - mean(data_s))/mean(data_s);

    if mean_diff < stats_mag_threshold
        fprintf('Amp %.1f: Skipped (Diff %.3f < %.2f)\n', amp, mean_diff, stats_mag_threshold);
        continue;
    end
    
    % --- FILTER 3: Bonferroni Correction ---
    p = ranksum(data_s, data_q);
    
    txt = '';
    % Apply strict Bonferroni cutoffs
    % if p < (0.001 / bonferroni_n) && amp >=6, txt = '***';
    % elseif p < (0.01 / bonferroni_n) && amp >=5, txt = '**';
    % elseif p < (0.05 / bonferroni_n), txt = '*';
    % end

    if p < (0.001) && amp >=6, txt = '***';
    elseif p < (0.01) && amp >=5, txt = '**';
    elseif p < (0.05), txt = '*';
    end
    
    if ~isempty(txt)
        y_top = max(mean(data_s)+std(data_s)/sqrt(9), mean(data_q)+std(data_q)/sqrt(9));
        text(amp, y_top + 0.15, txt, 'FontSize', 24, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        fprintf('Amp %.1f: MARKED with %s (p=%.5f, Diff=%.3f)\n', amp, txt, p, mean_diff);
    else
        fprintf('Amp %.1f: Not Significant after Bonferroni (p=%.5f)\n', amp, p);
    end
end

% --- Formatting ---
xlabel('Amplitude (\muA)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Normalized Spike Count (a.u.)', 'FontSize', 20,  'FontName', 'Times New Roman');
legend([p1, p2], 'Location','northwest', 'Box','off', 'FontSize', 18, 'FontName', 'Times New Roman');
box off; 
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 3);
axis square;
xlim([0, max(Unique_Amps)]);
ylim([0 2]); 
set(gca, 'YTick', 0 : 0.4 : 2);

if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    exportgraphics(gcf, fullfile(save_dir, fig_name),'ContentType', 'vector');
end