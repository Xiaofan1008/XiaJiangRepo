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
% List Result_PeakLatency_*.mat files
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

% --- DEFINE CUTOFFS GLOBALLY HERE ---
% (Applied to both Channel calculations and Histograms)
max_cutoff_sim = 14; 
max_cutoff_seq = 18;

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
            
            % Initialize
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim = [];     % For Histogram (All Spikes)
                Pooled.(fName).Seq = [];
                Pooled.(fName).Sim_Ch = [];  % NEW: For Summary (Channel Medians)
                Pooled.(fName).Seq_Ch = [];
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            % 1. Get Raw Data
            raw_sim = ThisAmpData.Sim(:);
            raw_seq = ThisAmpData.Seq(:);
            
            % 2. Filter Raw Data (Apply Cutoffs)
            raw_sim = raw_sim(raw_sim <= max_cutoff_sim);
            raw_seq = raw_seq(raw_seq <= max_cutoff_seq);
            
            % 3. Store ALL spikes (for Histogram)
            Pooled.(fName).Sim = [Pooled.(fName).Sim; raw_sim];
            Pooled.(fName).Seq = [Pooled.(fName).Seq; raw_seq];
            
            % 4. Store CHANNEL MEDIAN (for Summary Plot)
            if ~isempty(raw_sim)
                Pooled.(fName).Sim_Ch = [Pooled.(fName).Sim_Ch; median(raw_sim)];
            end
            if ~isempty(raw_seq)
                Pooled.(fName).Seq_Ch = [Pooled.(fName).Seq_Ch; median(raw_seq)];
            end
        end
    end
end
All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);

%% ================= 3. PLOT HISTOGRAMS (Black & White) =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Generating plots for %d Amplitudes...\n', length(All_Amps));
for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    if curr_amp == 0, continue; end
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    sim_data = Pooled.(fName).Sim;
    seq_data = Pooled.(fName).Seq;
    
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
    
    % 3. Plot Grouped Bars
    b = bar(centers, [N_sim; N_seq]', 'grouped', 'BarWidth', 0.9);
    
    % --- STYLE SETTINGS ---
    b(1).FaceColor = [0.9 0.9 0.9]; b(1).EdgeColor = 'k'; b(1).LineWidth = 1.0;
    b(1).DisplayName = 'Simultaneous';
    
    b(2).FaceColor = 'k'; b(2).EdgeColor = 'none';
    b(2).DisplayName = sprintf('Sequential (N=%d)', length(seq_data)); % Keeping N=spikes here is fine for histograms
    
    % 4. Add Median Lines
    if ~isempty(sim_data)
        xline(median(sim_data), '--k', 'LineWidth', 1.5, 'DisplayName', sprintf('Simultaneous Median: %.1fms', median(sim_data)));
    end
    if ~isempty(seq_data)
        xline(median(seq_data), '-k', 'LineWidth', 1.5, 'DisplayName', sprintf('Sequential Median: %.1fms', median(seq_data)));
    end
    
    % --- Formatting ---
    xlabel('Time (ms)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
    ylabel('Fraction', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
    % title(sprintf('Peak Latency @ %.1f \\muA', curr_amp), 'FontSize', 14, 'FontName', 'Times New Roman');
    legend('Location','northeast', 'Box','off', 'FontSize', 14, 'FontName', 'Times New Roman');
    
    box off; 
    set(gca, 'FontSize', 13, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 2);
    % axis square;
    % xlim([min_t-0.5, max_t+0.5]);
    xlim([0-0.5, 18+0.5]);
    % set(gca, 'XTick', min_t : 1 : max_t);
    set(gca, 'XTick', 0 : 1 : 18);
    ylim([0, 0.35]); 
    
    if save_figure
        exportgraphics(gcf, fullfile(save_dir, [fig_name '.tiff']), 'ContentType', 'vector');
    end
end

%% ================= 4. SUMMARY PLOT: MEDIAN LATENCY (Per Channel) =================
% fprintf('Generating Summary Plot (Mean of Channel Medians)...\n');
% % 1. Pre-allocate
% Med_Sim = nan(size(All_Amps)); Err_Sim = nan(size(All_Amps));
% Med_Seq = nan(size(All_Amps)); Err_Seq = nan(size(All_Amps));
% 
% % 2. Calculate Stats (Using Channel Medians)
% for i = 1:length(All_Amps)
%     curr_amp = All_Amps(i);
%     if curr_amp == 0, continue; end
%     fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
% 
%     % Get List of Medians (One per Channel)
%     sim_ch = Pooled.(fName).Sim_Ch;
%     seq_ch = Pooled.(fName).Seq_Ch;
% 
%     % --- Simultaneous ---
%     if ~isempty(sim_ch)
%         Med_Sim(i) = mean(sim_ch); % Mean of the channel medians
%         % Standard Error of the Mean (SEM) = STD / sqrt(N_channels)
%         Err_Sim(i) = std(sim_ch) / sqrt(length(sim_ch)); 
%     end
% 
%     % --- Sequential ---
%     if ~isempty(seq_ch)
%         Med_Seq(i) = mean(seq_ch);
%         Err_Seq(i) = std(seq_ch) / sqrt(length(seq_ch));
%     end
% end
% 
% % 3. Plotting
% fig_name_sum = 'Summary_MedianLatency_vs_Amp';
% figure('Color','w', 'Position', [150 150 600 500], 'Name', fig_name_sum); 
% hold on;
% 
% % --- Plot Simultaneous ---
% errorbar(All_Amps, Med_Sim, Err_Sim, '--o', ...
%     'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
%     'MarkerFaceColor', 'w', 'CapSize', 10, ...
%     'DisplayName', 'Simultaneous');
% 
% % --- Plot Sequential ---
% errorbar(All_Amps, Med_Seq, Err_Seq, '-s', ...
%     'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
%     'MarkerFaceColor', 'k', 'CapSize', 10, ...
%     'DisplayName', 'Sequential');
% 
% % --- Formatting ---
% xlabel('Amplitude (\muA)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
% ylabel('Peak Firing Rate Time (ms)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
% legend('Location','best', 'Box','off', 'FontSize', 16, 'FontName', 'Times New Roman');
% 
% box off; 
% set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 2);
% axis square;
% 
% xlim([1-0.5, max(All_Amps)+0.5]);
% xticks(1:1:max(All_Amps));
% ylim([0 16]); % Adjusted for cleaner view
% 
% % --- Save ---
% if save_figure
%     exportgraphics(gcf, fullfile(save_dir, [fig_name_sum '.tiff']), 'ContentType', 'vector');
%     % saveas(gcf, fullfile(save_dir, [fig_name_sum '.fig']));
%     fprintf('Summary plot saved: %s\n', fig_name_sum);
%     fprintf('\n>>> All plots saved to: %s\n', save_dir);
% end


%% ================= 4. SUMMARY PLOT: CORRELATION TEST (Lat vs Amp) =================
% fprintf('Generating Summary Plot & Testing Correlations...\n');
% 
% % 1. Pre-allocate for Plotting
% Med_Sim = nan(size(All_Amps)); Err_Sim = nan(size(All_Amps));
% Med_Seq = nan(size(All_Amps)); Err_Seq = nan(size(All_Amps));
% 
% % 2. Pre-allocate for Correlation Test (Grouping all data points)
% % We need paired lists: [Amp1, Lat1; Amp1, Lat2; Amp2, Lat3...]
% CorrData_Sim_Amp = []; CorrData_Sim_Lat = [];
% CorrData_Seq_Amp = []; CorrData_Seq_Lat = [];
% 
% % 3. Calculate Stats & Accumulate Data
% for i = 1:length(All_Amps)
%     curr_amp = All_Amps(i);
%     if curr_amp == 0, continue; end
%     fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
% 
%     % Get List of Medians (One per Channel)
%     sim_ch = Pooled.(fName).Sim_Ch;
%     seq_ch = Pooled.(fName).Seq_Ch;
% 
%     % --- Simultaneous ---
%     if ~isempty(sim_ch)
%         Med_Sim(i) = mean(sim_ch); 
%         Err_Sim(i) = std(sim_ch) / sqrt(length(sim_ch));
% 
%         % Collect for Correlation
%         CorrData_Sim_Amp = [CorrData_Sim_Amp; repmat(curr_amp, length(sim_ch), 1)];
%         CorrData_Sim_Lat = [CorrData_Sim_Lat; sim_ch];
%     end
% 
%     % --- Sequential ---
%     if ~isempty(seq_ch)
%         Med_Seq(i) = mean(seq_ch);
%         Err_Seq(i) = std(seq_ch) / sqrt(length(seq_ch));
% 
%         % Collect for Correlation
%         CorrData_Seq_Amp = [CorrData_Seq_Amp; repmat(curr_amp, length(seq_ch), 1)];
%         CorrData_Seq_Lat = [CorrData_Seq_Lat; seq_ch];
%     end
% end
% 
% % 4. PERFORM CORRELATION TEST (Spearman)
% % Null Hypothesis: No correlation (rho = 0). 
% % If p > 0.05, we fail to reject Null -> No significant relationship (Result you want).
% 
% [rho_sim, p_sim] = corr(CorrData_Sim_Amp, CorrData_Sim_Lat, 'Type', 'Spearman', 'Rows','complete');
% [rho_seq, p_seq] = corr(CorrData_Seq_Amp, CorrData_Seq_Lat, 'Type', 'Spearman', 'Rows','complete');
% 
% fprintf('\n=== Correlation Results (Latency vs. Amplitude) ===\n');
% fprintf('Simultaneous: rho = %.3f, p = %.4f\n', rho_sim, p_sim);
% fprintf('Sequential:   rho = %.3f, p = %.4f\n', rho_seq, p_seq);
% if p_sim > 0.05, fprintf('>> Sim Result: NO significant correlation (Supported)\n'); end
% if p_seq > 0.05, fprintf('>> Seq Result: NO significant correlation (Supported)\n'); end
% 
% 
% % 5. Plotting
% fig_name_sum = 'Summary_MedianLatency_vs_Amp_Corr';
% figure('Color','w', 'Position', [150 150 600 500], 'Name', fig_name_sum); 
% hold on;
% 
% % --- Plot Simultaneous ---
% errorbar(All_Amps, Med_Sim, Err_Sim, '--o', ...
%     'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
%     'MarkerFaceColor', 'w', 'CapSize', 10, ...
%     'DisplayName', 'Simultaneous');
% 
% % --- Plot Sequential ---
% errorbar(All_Amps, Med_Seq, Err_Seq, '-s', ...
%     'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 10, ...
%     'MarkerFaceColor', 'k', 'CapSize', 10, ...
%     'DisplayName', 'Sequential');
% 
% % --- Formatting ---
% xlabel('Amplitude (\muA)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
% ylabel('Peak Firing Rate Time (ms)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
% 
% % Add Stats to Title or Legend
% title_str = sprintf('Latency vs. Amp (Sim: p=%.2f, Seq: p=%.2f)', p_sim, p_seq);
% title(title_str, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight','normal');
% 
% legend('Location','best', 'Box','off', 'FontSize', 16, 'FontName', 'Times New Roman');
% 
% box off; 
% set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 2);
% axis square;
% xlim([1-0.5, max(All_Amps)+0.5]);
% xticks(1:1:max(All_Amps));
% ylim([0 16]); 
% 
% % --- Save ---
% if save_figure
%     exportgraphics(gcf, fullfile(save_dir, [fig_name_sum '.tiff']), 'ContentType', 'vector');
%     fprintf('Summary plot saved: %s\n', fig_name_sum);
%     fprintf('\n>>> All plots saved to: %s\n', save_dir);
% end


%% ================= 4. SUMMARY PLOT: BINNED AMPLITUDES (Low/Mid/High) =================
fprintf('Generating Summary Plot (Binned Groups) & Correlation Test...\n');

% --- 1. Define Groups ---
% Structure to hold pooled channel medians for each group
Groups = struct();
Groups.Low.Sim = []; Groups.Low.Seq = [];   % Range: <= 3 uA
Groups.Mid.Sim = []; Groups.Mid.Seq = [];   % Range: 4 - 6 uA
Groups.High.Sim = []; Groups.High.Seq = []; % Range: >= 8 uA

% Arrays for Correlation Test (Continuous Data)
CorrData_Sim_Amp = []; CorrData_Sim_Lat = [];
CorrData_Seq_Amp = []; CorrData_Seq_Lat = [];

% --- 2. Loop through Amplitudes & Sort Data ---
for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    if curr_amp == 0, continue; end
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    % Get Channel Medians
    sim_ch = Pooled.(fName).Sim_Ch;
    seq_ch = Pooled.(fName).Seq_Ch;
    
    % A. Accumulate for Correlation (Continuous)
    if ~isempty(sim_ch)
        CorrData_Sim_Amp = [CorrData_Sim_Amp; repmat(curr_amp, length(sim_ch), 1)];
        CorrData_Sim_Lat = [CorrData_Sim_Lat; sim_ch];
    end
    if ~isempty(seq_ch)
        CorrData_Seq_Amp = [CorrData_Seq_Amp; repmat(curr_amp, length(seq_ch), 1)];
        CorrData_Seq_Lat = [CorrData_Seq_Lat; seq_ch];
    end
    
    % B. Sort into Bins (Low/Mid/High)
    if curr_amp <= 3.1
        % Low (1-3 uA)
        Groups.Low.Sim = [Groups.Low.Sim; sim_ch];
        Groups.Low.Seq = [Groups.Low.Seq; seq_ch];
    elseif curr_amp >= 3.9 && curr_amp <= 6.1
        % Mid (4-6 uA)
        Groups.Mid.Sim = [Groups.Mid.Sim; sim_ch];
        Groups.Mid.Seq = [Groups.Mid.Seq; seq_ch];
    elseif curr_amp >= 7.9
        % High (8-10 uA)
        Groups.High.Sim = [Groups.High.Sim; sim_ch];
        Groups.High.Seq = [Groups.High.Seq; seq_ch];
    end
end

% --- 3. Calculate Stats for Plotting (Mean +/- SEM) ---
PlotVals.Sim = []; PlotErr.Sim = [];
PlotVals.Seq = []; PlotErr.Seq = [];

categories = {'Low', 'Mid', 'High'};
fields = {'Low', 'Mid', 'High'};

for k = 1:length(fields)
    f = fields{k};
    
    % Simultaneous
    data_sim = Groups.(f).Sim;
    if ~isempty(data_sim)
        PlotVals.Sim(k) = mean(data_sim);
        PlotErr.Sim(k)  = std(data_sim) / sqrt(length(data_sim));
    else
        PlotVals.Sim(k) = NaN; PlotErr.Sim(k) = NaN;
    end
    
    % Sequential
    data_seq = Groups.(f).Seq;
    if ~isempty(data_seq)
        PlotVals.Seq(k) = mean(data_seq);
        PlotErr.Seq(k)  = std(data_seq) / sqrt(length(data_seq));
    else
        PlotVals.Seq(k) = NaN; PlotErr.Seq(k) = NaN;
    end
end

% --- 4. Perform Correlation Test (Spearman on Raw Data) ---
[rho_sim, p_sim] = corr(CorrData_Sim_Amp, CorrData_Sim_Lat, 'Type', 'Spearman', 'Rows','complete');
[rho_seq, p_seq] = corr(CorrData_Seq_Amp, CorrData_Seq_Lat, 'Type', 'Spearman', 'Rows','complete');

fprintf('\n=== Correlation Results (Latency vs. Amplitude) ===\n');
fprintf('Simultaneous: rho = %.3f, p = %.4f\n', rho_sim, p_sim);
fprintf('Sequential:   rho = %.3f, p = %.4f\n', rho_seq, p_seq);


% --- 5. Plotting ---
fig_name_sum = 'Summary_Binned_Latency_vs_Amp';
figure('Color','w', 'Position', [150 150 600 500], 'Name', fig_name_sum); 
hold on;

x_pts = 1:3; % Positions for Low, Mid, High

% Plot Simultaneous
errorbar(x_pts, PlotVals.Sim, PlotErr.Sim, '--o', ...
    'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'w', 'CapSize', 15, ...
    'DisplayName', 'Simultaneous');

% Plot Sequential
errorbar(x_pts, PlotVals.Seq, PlotErr.Seq, '-s', ...
    'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 12, ...
    'MarkerFaceColor', 'k', 'CapSize', 15, ...
    'DisplayName', 'Sequential');

% Formatting
ylabel('Peak Firing Rate Time (ms)', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');
xlabel('Amplitude Range', 'FontSize', 18, 'FontWeight','bold', 'FontName', 'Times New Roman');

% Set X-Axis Labels
set(gca, 'XTick', 1:3, 'XTickLabel', {'Low (1-3\muA)', 'Mid (4-6\muA)', 'High (8-10\muA)'});

% Title with Stats
% title_str = sprintf('Trend Analysis (Sim: p=%.2f, Seq: p=%.3f)', p_sim, p_seq);
% title(title_str, 'FontSize', 14, 'FontName', 'Times New Roman', 'FontWeight','normal');

legend('Location','best', 'Box','off', 'FontSize', 16, 'FontName', 'Times New Roman');

box off; 
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 2);
xlim([0.5, 3.5]);
ylim([0 16]); % Adjust Y-limit as needed

% Save
if save_figure
    exportgraphics(gcf, fullfile(save_dir, [fig_name_sum '.tiff']), 'ContentType', 'vector');
    fprintf('Summary plot saved: %s\n', fig_name_sum);
end