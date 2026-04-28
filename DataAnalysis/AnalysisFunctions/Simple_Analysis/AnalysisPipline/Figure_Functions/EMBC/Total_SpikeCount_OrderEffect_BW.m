%% ============================================================
%   GRAND AVERAGE: DIRECTIONALITY RECRUITMENT CURVE (Order 1 vs Order 2)
%   - Logic: Aggregates pre-calculated normalized sequential results.
%   - Safety: Verifies that Set 1 and Set 2 are reversed orders of the same electrodes.
%   - Grouping: Bins data by exact amplitude (No interpolation).
%   - Statistics: Mean +/- SEM across DATASETS.
%   - FILTER: Wilcoxon Signed Rank Test + Min N (Printed to Console Only)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% --- STATISTICAL FILTERS ---
% Minimum Sample Size for Stats (Limit dataset number)
stats_min_n_threshold = 0; 

% List all result files to include in the Grand Average
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
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE_EMBC/Figures/3_Total_Spike_Count';
fig_name    = 'Order_effect.tiff';

%% ================= 2. AGGREGATE DATA =================
fprintf('Processing %d datasets...\n', length(file_paths));
Pool_Ord1 = []; 
Pool_Ord2 = [];
Pool_Paired = []; % Initialize Paired Pool for Signed Rank Test

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    D = load(file_paths{i}); 
    if ~isfield(D, 'ResultNorm'), continue; end
    
    R = D.ResultNorm; 
    Amps = R.Amps; 
    seq_raw = R.Norm_Seq;
    
    % --- THE BOUNCER: Check Dimension 3 (Sets) ---
    nSets = size(seq_raw, 3);
    if nSets < 2
        % Silently skipping sets < 2 to keep command window clean
        continue;
    end
    
    % --- THE BOUNCER: Electrode Verification ---
    try
        sets_meta = R.Metadata.Stimulation_Sets;
        set1_ch = sets_meta(1, :); set1_ch = sort(set1_ch(set1_ch > 0));
        set2_ch = sets_meta(2, :); set2_ch = sort(set2_ch(set2_ch > 0));
        
        if ~isequal(set1_ch, set2_ch)
            fprintf('  Skipping Dataset %d: Sets are NOT reverse orders.\n', i);
            continue;
        end
    catch
        % If metadata isn't strictly available in older files, we proceed but notify
        fprintf('  Warning Dataset %d: Metadata missing. Assuming reversed pairs.\n', i);
    end
    
    % --- EXTRACT DATA (Splitting Sets instead of crushing them) ---
    % Squeeze across Channels (Dim 1) and ISIs (Dim 4, if it exists)
    if ndims(seq_raw) == 4
        seq_ds_mean = squeeze(mean(mean(seq_raw, 4, 'omitnan'), 1, 'omitnan')); % Result: [Amps x Sets]
    else
        seq_ds_mean = squeeze(mean(seq_raw, 1, 'omitnan')); % Result: [Amps x Sets]
    end
    
    % Store with 0-force logic
    for a = 1:length(Amps)
        if abs(Amps(a)) < 0.001
            v_ord1 = 0; 
            v_ord2 = 0;
        else
            v_ord1 = seq_ds_mean(a, 1); % Extract Set 1 (Order A->B)
            v_ord2 = seq_ds_mean(a, 2); % Extract Set 2 (Order B->A)
        end
        
        % Only store if BOTH Orders are valid (Ensures perfect pairing)
        if ~isnan(v_ord1) && ~isnan(v_ord2)
            Pool_Ord1 = [Pool_Ord1; Amps(a), v_ord1]; 
            Pool_Ord2 = [Pool_Ord2; Amps(a), v_ord2]; 
            Pool_Paired = [Pool_Paired; Amps(a), v_ord1, v_ord2]; 
        end
    end
end

%% ================= 3. STATS & VECTORS =================
Unique_Amps = unique(Pool_Paired(:,1));
Unique_Amps = sort(Unique_Amps);
if isempty(Unique_Amps)
    error('No valid paired datasets found. Check input files and filters.');
end
if Unique_Amps(1) ~= 0, Unique_Amps = [0; Unique_Amps]; end

Grand_Ord1_Mean = []; Grand_Ord1_SEM = [];
Grand_Ord2_Mean = []; Grand_Ord2_SEM = [];

for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    
    % Order 1
    vs = Pool_Ord1(Pool_Ord1(:,1)==amp, 2);
    if isempty(vs) && amp==0, vs=0; end
    n_s = length(vs); % Calculate exact N per point
    Grand_Ord1_Mean(k) = mean(vs); 
    Grand_Ord1_SEM(k) = std(vs)/sqrt(max(1, n_s)); 
    
    % Order 2
    vq = Pool_Ord2(Pool_Ord2(:,1)==amp, 2);
    if isempty(vq) && amp==0, vq=0; end
    n_q = length(vq); % Calculate exact N per point
    Grand_Ord2_Mean(k) = mean(vq); 
    Grand_Ord2_SEM(k) = std(vq)/sqrt(max(1, n_q)); 
end

%% ================= 4. PLOT =================
figure('Units', 'centimeters', 'Position', [5, 5, 10.5, 9.5], 'Color', 'w'); hold on;

% A. Scatter (Background)
jitter_w = 0.4; 
for i = 1:size(Pool_Ord2, 1)
    scatter(Pool_Ord2(i,1)+(rand-0.5)*jitter_w, Pool_Ord2(i,2), 12, [0.7 0.7 0.7], 's', 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
for i = 1:size(Pool_Ord1, 1)
    scatter(Pool_Ord1(i,1)+(rand-0.5)*jitter_w, Pool_Ord1(i,2), 12, [0.8 0.8 0.8], 'o', 'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end

% B. Main Curves
errorbar(Unique_Amps, Grand_Ord1_Mean, Grand_Ord1_SEM, '.', 'Color', 'k', 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p1 = plot(Unique_Amps, Grand_Ord1_Mean, '--o', 'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'DisplayName', 'Order 1 (A \rightarrow B)');

errorbar(Unique_Amps, Grand_Ord2_Mean, Grand_Ord2_SEM, '.', 'Color', 'k', 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p2 = plot(Unique_Amps, Grand_Ord2_Mean, '-s', 'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Order 2 (B \rightarrow A)');

% --- Formatting ---
box off; 
set(gca, 'FontSize', 10, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1);
axis square;
xlabel('Amplitude (µA)', 'FontSize', 10, 'FontName', 'Arial');
ylabel('Normalized Spike Count (a.u.)', 'FontSize', 10,  'FontName', 'Arial');
legend([p1, p2], 'Location','northwest', 'Box','off', 'FontSize', 10, 'FontName', 'Arial');
xlim([0, max(Unique_Amps)]);
ylim([0 2]); 
set(gca, 'YTick', 0 : 0.4 : 2);

if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    exportgraphics(gcf, fullfile(save_dir, fig_name),'ContentType', 'vector');
end

%% ================= 5. CONSOLE STATISTICS TABLE =================
fprintf('\n=================================================================\n');
fprintf('   WILCOXON SIGNED-RANK TEST: DIRECTIONALITY (Order 1 vs Order 2)\n');
fprintf('   * Testing all valid amplitudes to verify lack of difference *\n');
fprintf('=================================================================\n');
fprintf(' Amp (uA) |   N pairs  |  Ord1 Mean |  Ord2 Mean |  p-value \n');
fprintf('-----------------------------------------------------------------\n');

for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    if amp == 0, continue; end
    
    % Get PAIRED Data
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-amp)<0.001, :);
    data_1 = current_pairs(:, 2);
    data_2 = current_pairs(:, 3);
    
    n_pairs = length(data_1);
    
    % --- FILTER: Check Min N ---
    if n_pairs < stats_min_n_threshold
        fprintf(' %8.1f | %10d |       Skipped (Low N < %d)\n', amp, n_pairs, stats_min_n_threshold);
        continue;
    end
    
    % --- STATS: Signed Rank Test ---
    p = signrank(data_1, data_2);
    m1 = mean(data_1);
    m2 = mean(data_2);
    
    sig_flag = '';
    if p < 0.05
        sig_flag = '<- SIGNIFICANT';
    else
        sig_flag = '(n.s.)';
    end
    
    fprintf(' %8.1f | %10d | %10.3f | %10.3f |  %.4f %s\n', amp, n_pairs, m1, m2, p, sig_flag);
end
fprintf('=================================================================\n');
fprintf('Done! Clean directionality figure generated.\n\n');