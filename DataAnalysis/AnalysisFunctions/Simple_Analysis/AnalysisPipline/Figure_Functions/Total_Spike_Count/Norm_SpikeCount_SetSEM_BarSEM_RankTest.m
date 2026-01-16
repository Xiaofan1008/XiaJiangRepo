%% ============================================================
%   GRAND AVERAGE: NORMALIZED RECRUITMENT CURVE (Discrete Bins)
%   - Logic: Aggregates pre-calculated normalized results.
%   - Grouping: Bins data by exact amplitude (No interpolation).
%   - Statistics: Mean +/- SEM across DATASETS.
%   - NEW: Added Wilcoxon Rank Sum Test & Significance Stars
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
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE_EMBC/Figures/3_Total_Spike_Count';
fig_name    = 'GrandAverage_Recruitment_Curve_BarSEM_BW_Stats.tiff';

%% ================= 2. AGGREGATE DATA =================
fprintf('Processing %d datasets...\n', length(file_paths));
Pool_Sim = []; 
Pool_Seq = [];

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file')
        warning('File not found: %s', file_paths{i}); continue; 
    end
    
    D = load(file_paths{i});
    if ~isfield(D, 'ResultNorm'), continue; end
    R = D.ResultNorm;
    
    Amps = R.Amps; 
    
    % --- 1. Calculate Dataset Mean ---
    sim_dataset_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));
    
    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_dataset_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_dataset_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end
    
    % --- 2. Store in Master Pool (With Zero Force at 0uA) ---
    for a = 1:length(Amps)
        if abs(Amps(a)) < 0.001
            val_sim = 0; 
            val_seq = 0;
        else
            val_sim = sim_dataset_mean(a);
            val_seq = seq_dataset_mean(a);
        end
        
        if ~isnan(val_sim), Pool_Sim = [Pool_Sim; Amps(a), val_sim]; end
        if ~isnan(val_seq), Pool_Seq = [Pool_Seq; Amps(a), val_seq]; end
    end
end

%% ================= 3. CALCULATE GROUP STATISTICS =================
Unique_Amps = unique([Pool_Sim(:,1); Pool_Seq(:,1)]);
Unique_Amps = sort(Unique_Amps);

Grand_Sim_Mean = nan(size(Unique_Amps)); Grand_Sim_SEM = nan(size(Unique_Amps));
Grand_Seq_Mean = nan(size(Unique_Amps)); Grand_Seq_SEM = nan(size(Unique_Amps));
N_Sim = zeros(size(Unique_Amps)); N_Seq = zeros(size(Unique_Amps));

fprintf('\n=== GRAND AVERAGE STATISTICS ===\n');
for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    
    % Sim
    vals_sim = Pool_Sim(Pool_Sim(:,1) == amp, 2);
    n = 9; % Fixed N as per your previous code, or use length(vals_sim)
    if ~isempty(vals_sim)
        Grand_Sim_Mean(k) = mean(vals_sim);
        Grand_Sim_SEM(k)  = std(vals_sim) / sqrt(n);
        N_Sim(k) = n;
    end
    
    % Seq
    vals_seq = Pool_Seq(Pool_Seq(:,1) == amp, 2);
    n = 9;
    if ~isempty(vals_seq)
        Grand_Seq_Mean(k) = mean(vals_seq);
        Grand_Seq_SEM(k)  = std(vals_seq) / sqrt(n);
        N_Seq(k) = n;
    end
end

% Force 0 in Arrays if needed (usually handled by pool, but safe to keep)
if Unique_Amps(1) ~= 0
    Unique_Amps     = [0; Unique_Amps];
    Grand_Sim_Mean  = [0; Grand_Sim_Mean]; Grand_Sim_SEM = [0; Grand_Sim_SEM];
    Grand_Seq_Mean  = [0; Grand_Seq_Mean]; Grand_Seq_SEM = [0; Grand_Seq_SEM];
end

%% ================= 4. PLOT GRAND AVERAGE (B&W) =================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% --- A. Scatter ---
jitter_w = 0.4; 
% Seq Scatter
for i = 1:size(Pool_Seq, 1)
    amp = Pool_Seq(i,1); val = Pool_Seq(i,2);
    x_jit = amp + (rand()-0.5) * jitter_w;
    scatter(x_jit, val, 30, [0.7 0.7 0.7], 's', 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
% Sim Scatter
for i = 1:size(Pool_Sim, 1)
    amp = Pool_Sim(i,1); val = Pool_Sim(i,2);
    x_jit = amp + (rand()-0.5) * jitter_w;
    scatter(x_jit, val, 30, [0.8 0.8 0.8], 'o', 'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end

% --- B. Main Curves ---
valid = ~isnan(Grand_Sim_Mean);
errorbar(Unique_Amps(valid), Grand_Sim_Mean(valid), Grand_Sim_SEM(valid), '.', 'Color', 'k', 'LineWidth', 2, 'CapSize', 12, 'HandleVisibility', 'off');
p1 = plot(Unique_Amps(valid), Grand_Sim_Mean(valid), '--o', 'Color', 'k', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');

valid = ~isnan(Grand_Seq_Mean);
errorbar(Unique_Amps(valid), Grand_Seq_Mean(valid), Grand_Seq_SEM(valid), '.', 'Color', 'k', 'LineWidth', 2, 'CapSize', 12, 'HandleVisibility', 'off');
p2 = plot(Unique_Amps(valid), Grand_Seq_Mean(valid), '-s', 'Color', 'k', 'LineWidth', 4, 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

% --- 4b. ADD STATISTICS (Wilcoxon Rank Sum) ---
fprintf('\n=== STATISTICAL TESTING (Wilcoxon Rank Sum) ===\n');
for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    if amp == 0, continue; end % Skip 0uA
    
    % Get Raw Data for this amp
    data_sim = Pool_Sim(abs(Pool_Sim(:,1)-amp)<0.001, 2);
    data_seq = Pool_Seq(abs(Pool_Seq(:,1)-amp)<0.001, 2);
    
    if ~isempty(data_sim) && ~isempty(data_seq)
        % Perform Test
        p = ranksum(data_sim, data_seq);
        
        % Determine Star
        txt = '';
        if p < 0.001, txt = '***';
        elseif p < 0.01, txt = '**';
        elseif p < 0.05, txt = '*';
        end
        
        if ~isempty(txt)
            % Calculate position (above the highest error bar)
            mu_s = mean(data_sim); sem_s = std(data_sim)/sqrt(9);
            mu_q = mean(data_seq); sem_q = std(data_seq)/sqrt(9);
            y_top = max(mu_s+sem_s, mu_q+sem_q);
            
            % Draw Star
            text(amp, y_top + 0.1, txt, 'FontSize', 24, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            fprintf('Amp %.1f uA: p = %.4f (%s)\n', amp, p, txt);
        else
            fprintf('Amp %.1f uA: p = %.4f (n.s.)\n', amp, p);
        end
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

% --- Save ---
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    exportgraphics(gcf, fullfile(save_dir, fig_name),'ContentType', 'vector');
    fprintf('\nFigure saved to: %s\n', fullfile(save_dir, fig_name));
end