%% ============================================================
%   POOLED GRAND AVERAGE: PEAK LATENCY (Final Refinement)
%   - Metric: Time of Peak Firing Rate (ms)
%   - Fixes:
%       1. Aligns Curve Y-Axis (%) to match Bar Plot
%       2. Uses "Slate" colors (lower saturation)
%       3. Tighter smoothing bandwidth to show real peaks
% ============================================================
clear; clc;

%% ================= 1. SETTINGS & FILES =================
target_amp = 10; 

% --- Histogram Settings ---
hist_bin_width = 1; % 2ms bins (0-2, 2-4, 4-6...)
hist_edges = 0 : hist_bin_width : 20; 

% --- Smoothing Settings ---
% Lower bandwidth = Sharper peaks. Higher = Smoother curve.
kde_bandwidth = 1; 

% --- Color Palette (Low Saturation) ---
col_sim = [0.27 0.45 0.70]; % Slate Blue
col_seq = [0.85 0.50 0.30]; % Burnt Orange

file_list = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set1_Latency_SepPop_5ms_Xia_Exp1_Sim1_251125_112055.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set4_Latency_SepPop_5ms_Xia_Exp1_Sim4_251125_152849.mat';
};

%% ================= 2. DATA POOLING =================
All_Amps_Found = [];
Data_Holder = {}; 
for f = 1:length(file_list)
    if ~exist(file_list{f}, 'file'), error('File not found: %s', file_list{f}); end
    loaded = load(file_list{f});
    if isfield(loaded, 'ResultLat'), D = loaded.ResultLat; else, error('No ResultLat'); end
    Data_Holder{f} = D;
    All_Amps_Found = [All_Amps_Found; D.Metadata.Amps];
end
Master_Amps = unique(All_Amps_Found);
nAmps = length(Master_Amps);

Pooled_Sim = cell(1, nAmps);
Pooled_Seq = cell(1, nAmps);

for f = 1:length(file_list)
    D = Data_Holder{f};
    [~, loc] = ismember(D.Metadata.Amps, Master_Amps);
    for i = 1:length(loc)
        idx = loc(i);
        % Sim
        raw = D.Latency.Sim(:, i); 
        Pooled_Sim{idx} = [Pooled_Sim{idx}; raw(~isnan(raw))];
        % Seq
        raw = D.Latency.Seq(:, i, :); 
        Pooled_Seq{idx} = [Pooled_Seq{idx}; raw(~isnan(raw))];
    end
end

%% ================= 3. COMPUTE STATS =================
Grand_Mean_Sim = nan(1, nAmps); Grand_SEM_Sim = nan(1, nAmps);
Grand_Mean_Seq = nan(1, nAmps); Grand_SEM_Seq = nan(1, nAmps);

for i = 1:nAmps
    d = Pooled_Sim{i}; 
    if ~isempty(d), Grand_Mean_Sim(i)=mean(d); Grand_SEM_Sim(i)=std(d)/sqrt(length(d)); end
    d = Pooled_Seq{i}; 
    if ~isempty(d), Grand_Mean_Seq(i)=mean(d); Grand_SEM_Seq(i)=std(d)/sqrt(length(d)); end
end

%% ================= 4. FIG 1: TUNING CURVE (Scatter + Line) =================
figure('Color','w', 'Position',[100 500 600 500]); hold on;

% 1. Jittered Scatter (Raw Data)
jitter = 0.2;
for i = 1:nAmps
    v = Pooled_Sim{i};
    if ~isempty(v), scatter(Master_Amps(i)+(rand(size(v))-.5)*jitter, v, 20, col_sim, 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off'); end
    v = Pooled_Seq{i};
    if ~isempty(v), scatter(Master_Amps(i)+(rand(size(v))-.5)*jitter, v, 20, col_seq, 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off'); end
end

% 2. Shaded Error Bars
plot_shaded_error(Master_Amps, Grand_Mean_Sim, Grand_SEM_Sim, col_sim);
plot_shaded_error(Master_Amps, Grand_Mean_Seq, Grand_SEM_Seq, col_seq);

% 3. Mean Lines
plot(Master_Amps, Grand_Mean_Sim, '-o', 'Color', col_sim, 'LineWidth', 2, 'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');
plot(Master_Amps, Grand_Mean_Seq, '-s', 'Color', col_seq, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

xlabel('Amplitude (\muA)', 'FontWeight','bold'); ylabel('Peak Latency (ms)', 'FontWeight','bold');
title('Latency vs Amplitude', 'FontWeight','bold'); legend('Location','best','Box','off'); box off;

%% ================= 5. FIG 2: HISTOGRAM (Bar Plot) =================
idx_box = find(Master_Amps == target_amp);
data_sim = Pooled_Sim{idx_box};
data_seq = Pooled_Seq{idx_box};

figure('Color','w', 'Position',[750 500 600 500]); hold on;

if ~isempty(data_sim) || ~isempty(data_seq)
    [c1, ~] = histcounts(data_sim, hist_edges);
    [c2, ~] = histcounts(data_seq, hist_edges);
    
    % Convert to Percentage
    p1 = c1 / length(data_sim) * 100;
    p2 = c2 / length(data_seq) * 100;
    
    % Plot
    centers = hist_edges(1:end-1) + hist_bin_width/2;
    b = bar(centers, [p1; p2]', 'grouped');
    
    % Styling
    b(1).FaceColor = col_sim; b(1).FaceAlpha = 0.6; b(1).EdgeColor = 'none'; b(1).DisplayName = 'Simultaneous';
    b(2).FaceColor = col_seq; b(2).FaceAlpha = 0.6; b(2).EdgeColor = 'none'; b(2).DisplayName = 'Sequential';
    
    xlabel('Peak Latency Bin (ms)', 'FontWeight','bold'); ylabel('% of Channels', 'FontWeight','bold');
    title(sprintf('Latency Distribution at %.0f µA', target_amp), 'FontWeight','bold');
    legend('Location','northeast','Box','off'); 
    xticks(hist_edges); xlim([0 16]); box off;
end

%% ================= 6. FIG 3: SMOOTHED DISTRIBUTION (Aligned) =================
figure('Color','w', 'Position',[750 100 600 500]); hold on;

% Helper to compute Scaled KDE (Converts Density -> Percentage)
% Formula: Density * BinWidth * 100 = Percentage (roughly matches bar height)
[xi_sim, f_sim_pct] = get_scaled_kde(data_sim, kde_bandwidth, hist_bin_width);
[xi_seq, f_seq_pct] = get_scaled_kde(data_seq, kde_bandwidth, hist_bin_width);

% Plot Curves
if ~isempty(xi_sim)
    fill([xi_sim fliplr(xi_sim)], [f_sim_pct zeros(size(f_sim_pct))], col_sim, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(xi_sim, f_sim_pct, 'Color', col_sim, 'LineWidth', 3, 'DisplayName', 'Simultaneous');
end
if ~isempty(xi_seq)
    fill([xi_seq fliplr(xi_seq)], [f_seq_pct zeros(size(f_seq_pct))], col_seq, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(xi_seq, f_seq_pct, 'Color', col_seq, 'LineWidth', 3, 'DisplayName', 'Sequential');
end

xlabel('Peak Latency (ms)', 'FontWeight','bold'); ylabel('Estimated % of Population', 'FontWeight','bold');
title(sprintf('Smoothed Distribution (%.0f µA)', target_amp), 'FontWeight','bold');
legend('Location','northeast','Box','off'); box off; xlim([0 16]);

%% ================= HELPER FUNCTIONS =================
function plot_shaded_error(x, y, se, col)
    if numel(x)<2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid=~isnan(y)&~isnan(se); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function [xi, f_scaled] = get_scaled_kde(data, bw, bin_width)
    if length(data) < 2
        xi = []; f_scaled = []; return;
    end
    % Compute standard KDE (Density)
    [f, xi] = ksdensity(data, 'Bandwidth', bw, 'BoundaryCorrection', 'reflection');
    
    % Convert Density to Percentage to match the Histogram scale
    % Probability = Density * Bin_Width
    % Percentage  = Probability * 100
    f_scaled = f * bin_width * 100;
end