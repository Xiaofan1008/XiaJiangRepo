%% ============================================================
%   GRAND AVERAGE: NORMALIZED RECRUITMENT CURVE (Sigmoid Fit)
%   - Logic: Aggregates pre-calculated normalized results.
%   - Grouping: Bins data by exact amplitude.
%   - Statistics: Mean +/- SEM.
%   - FITTING: Fits a Sigmoid curve to the pooled data.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
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

save_figure = true;
save_dir    = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Group_Analysis/';
fig_name    = 'GrandAverage_Recruitment_Curve_Fitted.fig';

%% ================= 2. AGGREGATE DATA =================
fprintf('Processing %d datasets...\n', length(file_paths));

Pool_Sim = []; 
Pool_Seq = [];

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    
    D = load(file_paths{i});
    if ~isfield(D, 'ResultNorm'), continue; end
    R = D.ResultNorm;
    
    Amps = R.Amps;
    
    % Collapse Channels & Sets
    sim_dataset_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));
    
    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_dataset_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_dataset_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end
    
    % Store
    for a = 1:length(Amps)
        if ~isnan(sim_dataset_mean(a))
            Pool_Sim = [Pool_Sim; Amps(a), sim_dataset_mean(a)]; %#ok<*AGROW>
        end
        if ~isnan(seq_dataset_mean(a))
            Pool_Seq = [Pool_Seq; Amps(a), seq_dataset_mean(a)];
        end
    end
end

%% ================= 3. STATISTICS & FITTING =================
Unique_Amps = unique([Pool_Sim(:,1); Pool_Seq(:,1)]);
Unique_Amps = sort(Unique_Amps);

Grand_Sim_Mean = nan(size(Unique_Amps)); Grand_Sim_SEM = nan(size(Unique_Amps));
Grand_Seq_Mean = nan(size(Unique_Amps)); Grand_Seq_SEM = nan(size(Unique_Amps));

% A. Calculate Discrete Stats (For Error Bars)
for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    
    % Sim
    vals = Pool_Sim(Pool_Sim(:,1) == amp, 2);
    n = 9; % User forced N=9 (Animals)
    if ~isempty(vals)
        Grand_Sim_Mean(k) = mean(vals);
        Grand_Sim_SEM(k)  = std(vals) / sqrt(n);
    end
    
    % Seq
    vals = Pool_Seq(Pool_Seq(:,1) == amp, 2);
    if ~isempty(vals)
        Grand_Seq_Mean(k) = mean(vals);
        Grand_Seq_SEM(k)  = std(vals) / sqrt(n);
    end
end

% B. Curve Fitting (Sigmoid)
fprintf('Fitting Curves...\n');
x_smooth = linspace(min(Unique_Amps), max(Unique_Amps), 100);

[fit_sim_y, ~] = fit_sigmoid(Pool_Sim(:,1), Pool_Sim(:,2), x_smooth);
[fit_seq_y, ~] = fit_sigmoid(Pool_Seq(:,1), Pool_Seq(:,2), x_smooth);

%% ================= 4. PLOT =================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% 1. Plot Shaded Error Bars (Discrete Data)
plot_shaded(Unique_Amps, Grand_Sim_Mean, Grand_Sim_SEM, [0 0.4 0.8]); % Blue
plot_shaded(Unique_Amps, Grand_Seq_Mean, Grand_Seq_SEM, [0.8 0 0.2]); % Red

% 2. Plot Fitted Lines (Sigmoid)
plot(x_smooth, fit_sim_y, '-', 'Color', [0 0.3 0.8], 'LineWidth', 2.5, 'DisplayName', 'Simultaneous (Fit)');
plot(x_smooth, fit_seq_y, '-', 'Color', [0.7 0 0], 'LineWidth', 2.5, 'DisplayName', 'Sequential (Fit)');

% 3. Plot Discrete Mean Points (Markers Only)
plot(Unique_Amps, Grand_Sim_Mean, 'o', 'Color', [0 0.3 0.8], 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'HandleVisibility','off');
plot(Unique_Amps, Grand_Seq_Mean, 's', 'Color', [0.7 0 0], 'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'HandleVisibility','off');

% --- Formatting ---
yline(1.0, '--k', 'Ref (5uA)', 'HandleVisibility','off');
xlabel('Amplitude (\muA)', 'FontSize', 10, 'FontWeight','bold');
ylabel('Normalized Spike Count', 'FontSize', 10, 'FontWeight','bold');
title(sprintf('Grand Average Recruitment', length(file_paths)), 'FontSize', 10);
legend('Location','best','Box','off');
box off; set(gca, 'FontSize', 12);
xlim([1, max(Unique_Amps)]);

% --- Save ---
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    saveas(gcf, fullfile(save_dir, fig_name));
    fprintf('\nFigure saved to: %s\n', fullfile(save_dir, fig_name));
end

%% ================= LOCAL FUNCTIONS =================
function plot_shaded(x, y, err, col)
    valid = ~isnan(y) & ~isnan(err);
    if sum(valid) < 2, return; end
    x = x(valid); y = y(valid); err = err(valid);
    fill([x; flipud(x)], [y+err; flipud(y-err)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function [y_fit, coeffs] = fit_sigmoid(x_data, y_data, x_query)
    % Fits a 4-parameter Sigmoid (Boltzmann/Logistic)
    % y = bottom + (top - bottom) / (1 + exp(-(x - v50) / slope))
    
    % Initial Guesses
    y_min = min(y_data);
    y_max = max(y_data);
    x_half = median(x_data);
    slope = 1;
    
    p0 = [y_min, y_max, x_half, slope];
    
    % Sigmoid Function Definition
    sigmoid_fun = @(p, x) p(1) + (p(2) - p(1)) ./ (1 + exp(-(x - p(3)) ./ p(4)));
    
    % Fit (Robust using fminsearch to avoid Toolbox dependency)
    loss_fun = @(p) sum((y_data - sigmoid_fun(p, x_data)).^2);
    opts = optimset('Display','off');
    coeffs = fminsearch(loss_fun, p0, opts);
    
    % Evaluate
    y_fit = sigmoid_fun(coeffs, x_query);
end