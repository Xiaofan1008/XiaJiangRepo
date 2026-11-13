%% ============================================================
%   GRAND AVERAGE: PEAK LATENCY DISTRIBUTION (Pooled)
%   - Logic: 
%       1. Loads all single-dataset PeakLatency result files.
%       2. Pools ALL peak times together for each amplitude.
%       3. Plots Grouped Histograms (Sim vs Seq) per Amplitude.
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
save_figure = true;
save_dir    = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Group_Analysis/Peak_FR_Time/';
hist_bin_width = 1; % ms (Width of visual bars)

%% ================= 2. POOL DATA ACROSS DATASETS =================
fprintf('Pooling data from %d datasets...\n', length(file_paths));

% Structure to hold pooled data: Pooled.Amp_X.Sim = [all times...]
Pooled = struct();
All_Amps = [];

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file')
        warning('File not found: %s', file_paths{i}); continue; 
    end
    
    D = load(file_paths{i});
    if ~isfield(D, 'PeakData'), continue; end
    
    P = D.PeakData;
    Sets = fieldnames(P.Set); % Get sets (Set1, Set2...)
    
    % Loop through Sets in this file
    for s = 1:length(P.Set)
        if isempty(P.Set(s).Amp), continue; end
        
        % Loop through Amplitudes in this Set
        for a = 1:length(P.Set(s).Amp)
            ThisAmpData = P.Set(s).Amp(a);
            val = ThisAmpData.Val;
            
            % Create Field Name (e.g., A_5p0)
            fName = sprintf('A_%.1f', val);
            fName = strrep(fName, '.', 'p');
            
            % Initialize if new amplitude
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim = [];
                Pooled.(fName).Seq = [];
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            % Append Data (Pooling)
            Pooled.(fName).Sim = [Pooled.(fName).Sim; ThisAmpData.Sim(:)];
            Pooled.(fName).Seq = [Pooled.(fName).Seq; ThisAmpData.Seq(:)];
        end
    end
end

All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);

%% ================= 3. PLOT HISTOGRAMS (One Fig Per Amp) =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

fprintf('Generating plots for %d Amplitudes...\n', length(All_Amps));

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    sim_data = Pooled.(fName).Sim;
    seq_data = Pooled.(fName).Seq;
    
    % Skip empty amplitudes
    if isempty(sim_data) && isempty(seq_data), continue; end
    
    % --- Create Figure ---
    fig_name = sprintf('PeakLatency_Dist_%.1fuA', curr_amp);
    figure('Color','w', 'Position', [100 100 600 400], 'Name', fig_name);
    hold on;
    
    % 1. Define Bins
    % Find range to cover all data
    all_vals = [sim_data; seq_data];
    min_t = floor(min(all_vals)); 
    max_t = ceil(max(all_vals));
    if min_t < 2, min_t = 2; end % Clip to typical analysis window start
    
    edges = min_t : hist_bin_width : max_t;
    centers = edges(1:end-1) + hist_bin_width/2;
    
    % 2. Calculate Counts (Probability)
    N_sim = histcounts(sim_data, edges, 'Normalization', 'probability');
    N_seq = histcounts(seq_data, edges, 'Normalization', 'probability');
    
    % 3. Plot Grouped Bars (Side-by-Side)
    % Input must be [Bins x Groups]
    b = bar(centers, [N_sim; N_seq]', 'grouped');
    
    % Style
    b(1).FaceColor = [0 0.4 0.8]; % Sim (Blue)
    b(1).EdgeColor = 'none';
    b(1).DisplayName = sprintf('Simultaneous (N=%d)', length(sim_data));
    
    b(2).FaceColor = [0.8 0.1 0.1]; % Seq (Red)
    b(2).EdgeColor = 'none';
    b(2).DisplayName = sprintf('Sequential (N=%d)', length(seq_data));
    
    % 4. Add Median Lines (Dashed)
    if ~isempty(sim_data)
        xline(median(sim_data), '--', 'Color', [0 0.3 0.8], 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Sim Median: %.1fms', median(sim_data)));
    end
    if ~isempty(seq_data)
        xline(median(seq_data), '--', 'Color', [0.7 0 0], 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Seq Median: %.1fms', median(seq_data)));
    end
    
    % Format
    xlabel('Time from Pulse (ms)', 'FontSize', 10, 'FontWeight','bold');
    ylabel('Fraction of Channels', 'FontSize', 10, 'FontWeight','bold');
    title(sprintf('Peak Latency Distribution @ %.1f \\muA', curr_amp), 'FontSize', 10);
    legend('Location','best', 'Box','off');
    box off; set(gca, 'FontSize', 10);
    xlim([min_t-0.5, max_t+0.5]);
    ylim([0 1]); % Prob is 0-1
    
    % --- Save ---
    if save_figure
        saveas(gcf, fullfile(save_dir, [fig_name '.fig']));
        % saveas(gcf, fullfile(save_dir, [fig_name '.png'])); % Optional PNG
        fprintf('  Saved: %s\n', fig_name);
    end
end

fprintf('\n>>> All plots saved to: %s\n', save_dir);