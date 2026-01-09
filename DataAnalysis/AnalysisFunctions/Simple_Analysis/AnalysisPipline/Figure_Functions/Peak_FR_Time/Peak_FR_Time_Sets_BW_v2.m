%% ============================================================
%   GRAND AVERAGE: PEAK LATENCY DISTRIBUTION (Pooled)
%   - Logic: 
%       1. Loads all single-dataset PeakLatency result files.
%       2. Pools ALL peak times together for each amplitude.
%       3. Plots Grouped Histograms (Sim vs Seq) per Amplitude.
%   - STYLE: Black & White (Open vs Filled bars).
%   - UPDATES: Full X-ticks, Dynamic Y-limit.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% List all your Result_PeakLatency_*.mat files here
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX012/Result_PeakLatency_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX012/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX012/Result_PeakLatency_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX009/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX009/Result_PeakLatency_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX005/Result_PeakLatency_Separated_Xia_Exp1_Sim.mat';
};

% Plot Settings
save_figure = false;
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE EMBC/Figures/4_Peak_FR_Time';
hist_bin_width = 1; % ms (Width of visual bars)

%% ================= 2. POOL DATA ACROSS DATASETS =================
fprintf('Pooling data from %d datasets...\n', length(file_paths));

Pooled = struct();
All_Amps = [];

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    
    D = load(file_paths{i});
    if ~isfield(D, 'PeakData'), continue; end
    
    P = D.PeakData;
    Sets = fieldnames(P.Set); 
    
    for s = 1:length(P.Set)
        if isempty(P.Set(s).Amp), continue; end
        for a = 1:length(P.Set(s).Amp)
            ThisAmpData = P.Set(s).Amp(a);
            val = ThisAmpData.Val;
            fName = sprintf('A_%.1f', val);
            fName = strrep(fName, '.', 'p');
            
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim = [];
                Pooled.(fName).Seq = [];
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            Pooled.(fName).Sim = [Pooled.(fName).Sim; ThisAmpData.Sim(:)];
            Pooled.(fName).Seq = [Pooled.(fName).Seq; ThisAmpData.Seq(:)];
        end
    end
end

All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);

%% ================= 3. PLOT HISTOGRAMS (Black & White) =================
% if ~exist(save_dir, 'dir'), mkdir(save_dir); end
% 
% fprintf('Generating plots for %d Amplitudes...\n', length(All_Amps));
% 
% for i = 1:length(All_Amps)
%     curr_amp = All_Amps(i);
%     fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
% 
%     sim_data = Pooled.(fName).Sim;
%     seq_data = Pooled.(fName).Seq;
% 
%     if isempty(sim_data) && isempty(seq_data), continue; end
% 
%     % --- Create Figure ---
%     fig_name = sprintf('PeakLatency_Dist_%.1fuA_BW', curr_amp);
%     figure('Color','w', 'Position', [100 100 600 400], 'Name', fig_name);
%     hold on;
% 
%     % 1. Define Bins
%     all_vals = [sim_data; seq_data];
%     min_t = floor(min(all_vals)); 
%     max_t = ceil(max(all_vals));
%     if min_t < 2, min_t = 2; end 
% 
%     edges = min_t : hist_bin_width : max_t;
%     centers = edges(1:end-1) + hist_bin_width/2;
% 
%     % 2. Calculate Counts (Probability)
%     N_sim = histcounts(sim_data, edges, 'Normalization', 'probability');
%     N_seq = histcounts(seq_data, edges, 'Normalization', 'probability');
% 
%     % 3. Plot Grouped Bars (B&W Style)
%     b = bar(centers, [N_sim; N_seq]', 'grouped');
% 
%     % --- STYLE SETTINGS (Open vs Filled) ---
%     b(1).FaceColor = 'w'; 
%     b(1).EdgeColor = 'k';
%     b(1).LineWidth = 1.0;
%     b(1).DisplayName = sprintf('Simultaneous (N=%d)', length(sim_data));
% 
%     b(2).FaceColor = 'k'; 
%     b(2).EdgeColor = 'none';
%     b(2).DisplayName = sprintf('Sequential (N=%d)', length(seq_data));
% 
%     % 4. Add Median Lines (Distinct Styles)
%     if ~isempty(sim_data)
%         xline(median(sim_data), '--k', 'LineWidth', 1.5, ...
%             'DisplayName', sprintf('Sim Median: %.1fms', median(sim_data)));
%     end
%     if ~isempty(seq_data)
%         xline(median(seq_data), '-k', 'LineWidth', 1.5, ...
%             'DisplayName', sprintf('Seq Median: %.1fms', median(seq_data)));
%     end
% 
%     % Format
%     xlabel('Time (ms)', 'FontSize', 10, 'FontWeight','bold');
%     ylabel('Fraction', 'FontSize', 10, 'FontWeight','bold');
%     title(sprintf('Peak Latency Distribution @ %.1f \\muA', curr_amp), 'FontSize', 10);
%     legend('Location','best', 'Box','off');
%     box off; set(gca, 'FontSize', 10);    
%     xlim([min_t-0.5, max_t+0.5]);
% 
%     % Force X-axis to show every integer
%     set(gca, 'XTick', min_t : 1 : max_t);
% 
%     % Dynamic Y-limit to remove blank space (Max bar height + 10%)
%     max_bar_height = max([N_sim, N_seq]);
%     if max_bar_height > 0
%         % ylim([0, max_bar_height * 1.1]); 
%         ylim([0, 0.35]); 
%     else
%         ylim([0 1]);
%     end
% 
%     % --- Save ---
%     if save_figure
%         saveas(gcf, fullfile(save_dir, [fig_name '.tiff']));
%         fprintf('  Saved: %s\n', fig_name);
%     end
% end
% 
% fprintf('\n>>> All plots saved to: %s\n', save_dir);

if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Generating plots for %d Amplitudes...\n', length(All_Amps));

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    if curr_amp == 0, continue; end

    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    sim_data = Pooled.(fName).Sim;
    seq_data = Pooled.(fName).Seq;

    % --- FILTER: REMOVE LATE PEAKS ---
    max_cutoff_seq = 18; % Defined Window of Interest
    max_cutoff_sim = 14;
    sim_data = sim_data(sim_data <= max_cutoff_sim);
    seq_data = seq_data(seq_data <= max_cutoff_seq);
    
    if isempty(sim_data) && isempty(seq_data), continue; end
    
    % --- Create Figure ---
    fig_name = sprintf('PeakLatency_Dist_%.1fuA_BW', curr_amp);
    figure('Color','w', 'Position', [100 100 600 500], 'Name', fig_name);
    hold on;
    
    % 1. Define Bins
    all_vals = [sim_data; seq_data];
    min_t = floor(min(all_vals)); 
    max_t = ceil(max(all_vals));
    if min_t < 2, min_t = 2; end 
    
    edges = min_t : hist_bin_width : max_t;
    centers = edges(1:end-1) + hist_bin_width/2;
    
    % 2. Calculate Counts (Probability)
    N_sim = histcounts(sim_data, edges, 'Normalization', 'probability');
    N_seq = histcounts(seq_data, edges, 'Normalization', 'probability');
    
    % 3. Plot Grouped Bars (BarWidth = 1 to remove gaps)
    b = bar(centers, [N_sim; N_seq]', 'grouped', 'BarWidth', 0.9);
    
    % --- STYLE SETTINGS ---
    % Sim: Light Grey (Visible against white bg) with Black Edge
    b(1).FaceColor = [0.9 0.9 0.9]; 
    b(1).EdgeColor = 'k';
    b(1).LineWidth = 1.0;
    % RESTORED: Legend with N
    % b(1).DisplayName = sprintf('Simultaneous (N=%d)', length(sim_data));
    b(1).DisplayName = sprintf('Simultaneous');
    
    % Seq: Black Filled
    b(2).FaceColor = 'k'; 
    b(2).EdgeColor = 'none';
    % RESTORED: Legend with N
    b(2).DisplayName = sprintf('Sequential (N=%d)', length(seq_data));
    
    % 4. Add Median Lines (Thicker Lines = 2.0)
    % RESTORED: Legend with Median Values
    if ~isempty(sim_data)
        xline(median(sim_data), '--k', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Sim Median: %.1fms', median(sim_data)));
    end
    if ~isempty(seq_data)
        xline(median(seq_data), '-k', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Seq Median: %.1fms', median(seq_data)));
    end
    
    % --- Formatting ---
    % Labels (Size 16)
    xlabel('Time (ms)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
    ylabel('Fraction', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
    
    % Title (Size 14)
    title(sprintf('Peak Latency @ %.1f \\muA', curr_amp), 'FontSize', 14, 'FontName', 'Times New Roman');
    
    % Legend (Original Style with Data)
    % Note: You might need to adjust 'FontSize' to 12 if the text is too long for the box
    legend('Location','northeast', 'Box','off', 'FontSize', 16, 'FontName', 'Times New Roman');
    
    % Axis Appearance (Thicker Lines, Square Shape)
    box off; 
    set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.5);
    axis square;
    xticks(1:1:max(All_Amps));
    xlim([min_t-0.5, max_t+0.5]);
    set(gca, 'XTick', min_t : 1 : max_t);
    
    % Y-limit (Fixed or Dynamic)
    ylim([0, 0.35]); 
    
    % --- Save ---
    if save_figure
        % Export as PDF for best quality
        exportgraphics(gcf, fullfile(save_dir, [fig_name '.tiff']), 'ContentType', 'vector');
        % saveas(gcf, fullfile(save_dir, [fig_name '.fig']));
        fprintf('  Saved: %s\n', fig_name);
    end
end


%% ================= 4. SUMMARY PLOT: MEDIAN LATENCY vs. AMPLITUDE =================
fprintf('Generating Summary Plot (Median Latency vs. Amp)...\n');

% 1. Pre-allocate storage
Med_Sim = nan(size(All_Amps)); Err_Sim = nan(size(All_Amps));
Med_Seq = nan(size(All_Amps)); Err_Seq = nan(size(All_Amps));
% jitter_w = 0.4; % Width of scatter cloud

% Plotting
fig_name_sum = 'Summary_MedianLatency_vs_Amp';
figure('Color','w', 'Position', [150 150 600 500], 'Name', fig_name_sum); 
hold on;
% 2. Calculate Medians & Errors (Interquartile Range / sqrt(N) as proxy for SE)
for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    if curr_amp == 0, continue; end

    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');

    sim_data = Pooled.(fName).Sim;
    seq_data = Pooled.(fName).Seq;

    sim_data = sim_data(sim_data <= max_cutoff_sim);
    seq_data = seq_data(seq_data <= max_cutoff_seq);

    % --- Simultaneous ---
    if ~isempty(sim_data)
        Med_Sim(i) = median(sim_data);
        % Robust Standard Error estimate: IQR / (1.35 * sqrt(N))
        % Or just use IQR/2 for variability visibility. Let's use simple SEM estimate.
        % Actually, for non-normal distributions (latency), displaying the 
        % 25th-75th percentile range divided by sqrt(N) is a safe "error bar" for the median.
        % Here we use a standard bootstrapping proxy: 1.57 * IQR / sqrt(N)
        iqr_val = iqr(sim_data);
        n_val = length(sim_data);
        % Err_Sim(i) = (1.57 * iqr_val) / sqrt(n_val); 
        Err_Sim(i) = iqr_val / 2;
        % x_jit = curr_amp + (rand(size(sim_data))-0.5) * jitter_w;
        % scatter(x_jit, sim_data, 30, [0.8 0.8 0.8], 'o', 'filled', ...
        %     'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
    end

    % --- Sequential ---
    if ~isempty(seq_data)
        Med_Seq(i) = median(seq_data);
        iqr_val = iqr(seq_data);
        n_val = length(seq_data);
        % Err_Seq(i) = (1.57 * iqr_val) / sqrt(n_val);
        Err_Seq(i) = iqr_val / 2;
        % x_jit = curr_amp + (rand(size(seq_data))-0.5) * jitter_w;
        % scatter(x_jit, seq_data, 30, [0.7 0.7 0.7], 's', 'filled', ...
        %     'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
    end
end

% --- Plot Simultaneous (Dashed, Open Circles) ---
errorbar(All_Amps, Med_Sim, Err_Sim, '--o', ...
    'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'w', 'CapSize', 10, ...
    'DisplayName', 'Simultaneous');

% --- Plot Sequential (Solid, Filled Squares) ---
errorbar(All_Amps, Med_Seq, Err_Seq, '-s', ...
    'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'CapSize', 10, ...
    'DisplayName', 'Sequential');

% --- Formatting ---
xlabel('Amplitude (\muA)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
ylabel('Median Latency (ms)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
% title('Latency vs. Amplitude', 'FontSize', 14, 'FontName', 'Times New Roman');

legend('Location','best', 'Box','off', 'FontSize', 16, 'FontName', 'Times New Roman');

box off; 
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 1.5);
axis square;

% xlim([min(All_Amps)-0.5, max(All_Amps)+0.5]);
xlim([1-0.5, max(All_Amps)+0.5]);
ylim([0 16]); % Adjust based on your data range if needed

% --- Save ---
if save_figure
    exportgraphics(gcf, fullfile(save_dir, [fig_name_sum '.tiff']), 'ContentType', 'vector');
    saveas(gcf, fullfile(save_dir, [fig_name_sum '.fig']));
    fprintf('Summary plot saved: %s\n', fig_name_sum);
    fprintf('\n>>> All plots saved to: %s\n', save_dir);
end


