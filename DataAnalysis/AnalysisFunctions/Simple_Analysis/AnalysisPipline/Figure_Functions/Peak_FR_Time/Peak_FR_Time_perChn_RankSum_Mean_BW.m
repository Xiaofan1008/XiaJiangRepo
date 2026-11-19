%% ============================================================
%   GRAND AVERAGE: PEAK LATENCY DISTRIBUTION (Pooled)
%   - Logic: 
%       1. Loads all single-dataset PeakLatency result files.
%       2. Pools ALL peak times together for each amplitude.
%       3. Plots Grouped Histograms (Sim vs Seq) per Amplitude.
%   - STYLE: Black & White (Open vs Filled bars).
%   - UPDATES: Full X-ticks, Dynamic Y-limit, Stats Bracket.
%   - STATS: Per-Amplitude RankSum + Global RankSum Test.
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
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE_EMBC/Figures/4_Peak_FR_Time';
hist_bin_width = 1; % ms (Width of visual bars)

%% ================= 2. POOL DATA ACROSS DATASETS =================
fprintf('Pooling data from %d datasets...\n', length(file_paths));
Pooled = struct();
All_Amps = [];
% --- DEFINE CUTOFFS GLOBALLY HERE ---
max_cutoff_sim = 14; 
max_cutoff_seq = 16;
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

%% ================= 3. PLOT HISTOGRAMS (With Bracket Stats) =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('\nGenerating plots and Calculating Stats (RankSum on Channels)...\n');
% [NEW] Initialize Global Arrays for All-Data Test
Global_Sim_Ch = [];
Global_Seq_Ch = [];
for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    if curr_amp == 0, continue; end
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    % Get Data
    sim_data = Pooled.(fName).Sim; % Spikes for Plotting
    seq_data = Pooled.(fName).Seq;
    
    % Get Channel Medians for STATISTICS (Important for "795 channels")
    sim_ch = Pooled.(fName).Sim_Ch;
    seq_ch = Pooled.(fName).Seq_Ch;
    
    if isempty(sim_data) && isempty(seq_data), continue; end
    % [NEW] Accumulate for Global Test
    Global_Sim_Ch = [Global_Sim_Ch; sim_ch];
    Global_Seq_Ch = [Global_Seq_Ch; seq_ch];
    
    % --- PERFORM STATISTIC TEST (Rank Sum) & CALC MEAN/SEM ---
    p_val = ranksum(sim_ch, seq_ch);
    
    % Calculate Mean and SEM for this amplitude
    mean_sim = mean(sim_ch); sem_sim = std(sim_ch)/sqrt(length(sim_ch));
    mean_seq = mean(seq_ch); sem_seq = std(seq_ch)/sqrt(length(seq_ch));
    
    fprintf('Amp %.1f uA: Sim=%.2f+/-%.2f ms, Seq=%.2f+/-%.2f ms | Diff=%.2f, p=%.5f\n', ...
        curr_amp, mean_sim, sem_sim, mean_seq, sem_seq, median(seq_ch)-median(sim_ch), p_val);
    
    % --- Create Figure ---
    fig_name = sprintf('PeakLatency_Dist_%.1fuA_RankSum_BW', curr_amp);
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
    b = bar(centers, [N_sim; N_seq]'*100, 'grouped', 'BarWidth', 1);
    
    % --- STYLE SETTINGS ---
    b(1).FaceColor = [0.9 0.9 0.9]; b(1).EdgeColor = 'k'; b(1).LineWidth = 1.0;
    b(1).DisplayName = 'Simultaneous';
    
    b(2).FaceColor = 'k'; b(2).EdgeColor = 'none';
    b(2).DisplayName = sprintf('Sequential (N=%d)', length(seq_data)); 
    
    % 4. Add Median Lines
    med_sim = median(sim_data);
    med_seq = median(seq_data);
    
    if ~isempty(sim_data)
        xline(med_sim, '--k', 'LineWidth', 1.5, 'DisplayName', sprintf('Simultaneous Median: %.1fms', med_sim));
    end
    if ~isempty(seq_data)
        xline(med_seq, '-k', 'LineWidth', 1.5, 'DisplayName', sprintf('Sequential Median: %.1fms', med_seq));
    end
    
    % --- 5. ADD STATS BRACKET & STAR ---
    % Determine bracket position (Near Top of Plot)
    y_bracket = 23; % Fixed height since YLim is 25
    y_tick    = 1;  % Length of the bracket feet
    
    % Draw Bracket Line (Connecting the two medians)
    % Format: [x1 x1 x2 x2], [y_low y_high y_high y_low]
    plot([med_sim, med_sim, med_seq, med_seq], ...
         [y_bracket-y_tick, y_bracket, y_bracket, y_bracket-y_tick], ...
         '-k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Determine Star Text
    if p_val < 0.001, txt = '***';
    elseif p_val < 0.01, txt = '**';
    elseif p_val < 0.05, txt = '*';
    else, txt = 'n.s.';
    end
    
    % Draw Star above Bracket
    text((med_sim + med_seq)/2, y_bracket + 1.0, txt, ...
        'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    % --- Formatting ---
    xlabel('Time (ms)', 'FontSize', 22, 'FontName', 'Times New Roman');
    ylabel('Percentage (%)', 'FontSize', 22, 'FontName', 'Times New Roman');
    legend('Location','northeast', 'Box','off', 'FontSize', 14, 'FontName', 'Times New Roman');
    
    box off; 
    set(gca, 'FontSize', 13, 'FontName', 'Times New Roman', 'TickDir', 'out', 'LineWidth', 2);
    xlim([0, 17+0.5]);
    set(gca, 'XTick', 0 : 1 : 18);
    ylim([0, 25]); 
    set(gca, 'YTick', 0 : 5 : 25);
    axis square;
    if save_figure
        exportgraphics(gcf, fullfile(save_dir, [fig_name '.tiff']), 'ContentType', 'vector');
    end
end

% [NEW] Perform Global Statistic Test & Mean/SEM
p_global = ranksum(Global_Sim_Ch, Global_Seq_Ch);
mean_g_sim = mean(Global_Sim_Ch); sem_g_sim = std(Global_Sim_Ch)/sqrt(length(Global_Sim_Ch));
mean_g_seq = mean(Global_Seq_Ch); sem_g_seq = std(Global_Seq_Ch)/sqrt(length(Global_Seq_Ch));

fprintf('\n=== GLOBAL RESULT (All Data Pooled) ===\n');
fprintf('Total datasets: %d (Sim) vs %d (Seq)\n', length(Global_Sim_Ch), length(Global_Seq_Ch));
fprintf('Sim Mean: %.2f +/- %.2f ms\n', mean_g_sim, sem_g_sim);
fprintf('Seq Mean: %.2f +/- %.2f ms\n', mean_g_seq, sem_g_seq);
fprintf('Median Difference: %.2f ms\n', median(Global_Seq_Ch) - median(Global_Sim_Ch));
fprintf('Significance: p = %.10f\n', p_global);