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
save_figure = true;
save_dir    = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Group_Analysis/';
fig_name    = 'GrandAverage_Recruitment_Curve.fig';

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

% --- Helper to Plot Shaded Area ---
% We define anonymous functions or just loop for simplicity
plot_shaded(Unique_Amps, Grand_Sim_Mean, Grand_Sim_SEM, [0 0.4 0.8]); % Blue
plot_shaded(Unique_Amps, Grand_Seq_Mean, Grand_Seq_SEM, [0.8 0 0.2]); % Red

% --- Plot Lines ---
% Sim (Blue)
valid = ~isnan(Grand_Sim_Mean);
plot(Unique_Amps(valid), Grand_Sim_Mean(valid), '-o', ...
    'Color', [0 0.3 0.8], 'LineWidth', 2, 'MarkerFaceColor', [0 0.3 0.8], 'DisplayName', 'Simultaneous');

% Seq (Red)
valid = ~isnan(Grand_Seq_Mean);
plot(Unique_Amps(valid), Grand_Seq_Mean(valid), '-s', ...
    'Color', [0.7 0 0], 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

% --- Formatting ---
yline(1.0, '--k', 'Ref (5uA)', 'HandleVisibility','off');
xlabel('Amplitude (\muA)', 'FontSize', 10, 'FontWeight','bold');
ylabel('Normalized Spike Count', 'FontSize', 10, 'FontWeight','bold');
title(sprintf('Spike Count vs. AMP (N=%d)', length(file_paths)), 'FontSize', 10);
legend('Location','best','Box','off');
box off; set(gca, 'FontSize', 12);
% xlim([min(Unique_Amps)-0.5, max(Unique_Amps)+0.5]);
xlim([1, max(Unique_Amps)]);

% --- Save ---
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    saveas(gcf, fullfile(save_dir, fig_name));
    fprintf('\nFigure saved to: %s\n', fullfile(save_dir, fig_name));
end

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