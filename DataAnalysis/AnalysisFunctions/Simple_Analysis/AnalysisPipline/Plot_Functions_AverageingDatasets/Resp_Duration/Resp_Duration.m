%% ============================================================
%   GROUP DURATION ANALYSIS: Pooled Population
%   - Input: Result_Duration_*.mat files (from Single Analysis)
%   - Logic: Aggregates responding channel durations across all files.
%   - Outputs (Per Amplitude):
%       1. Grouped Box Plot (Distribution of Durations)
%       2. Line Plot (Mean +/- SEM with Statistics)
%   - Statistics: Wilcoxon Rank Sum Test (Sim vs Seq)
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
% List your single-dataset result files here
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX015/Result_Duration_Separated_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX015/Result_Duration_Separated_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX015/Result_Duration_Separated_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX015/Result_Duration_Separated_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX015/Result_Duration_Separated_Xia_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX015/Result_Duration_Separated_Xia_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX015/Result_Duration_Separated_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX014/Result_Duration_Separated_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX014/Result_Duration_Separated_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX014/Result_Duration_Separated_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX014/Result_Duration_Separated_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX014/Result_Duration_Separated_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX014/Result_Duration_Separated_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX013/Result_Duration_Separated_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX012/Result_Duration_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX012/Result_Duration_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX012/Result_Duration_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX011/Result_Duration_Separated_Xia_Exp1_Sim9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX010/Result_Duration_Separated_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX009/Result_Duration_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX005/Result_Duration_Separated_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX006/Result_Duration_Separated_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX006/Result_Duration_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX006/Result_Duration_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX006/Result_Duration_Separated_Xia_Exp1_Sim4.mat';
};

% Plot Settings
save_figures = false;
save_dir     = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Group_Analysis/Duration/';
fig_pos      = [100 100 1000 600];
color_sim    = [0 0.4 0.8]; % Blue
color_seq    = [0.8 0.1 0.1]; % Red

%% ================= 1. POOL DATA =================
fprintf('Pooling data from %d files...\n', length(file_paths));
Pooled = struct();
All_Amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), warning('Missing: %s', file_paths{f}); continue; end
    
    D = load(file_paths{f});
    if ~isfield(D, 'DurationData'), warning('No DurationData in %s', file_paths{f}); continue; end
    
    Data = D.DurationData;
    
    % Iterate through Sets
    for s = 1:length(Data.Set)
        if isempty(Data.Set(s).Label) || ~isfield(Data.Set(s), 'Amp'), continue; end
        
        % Safety check for the structure expansion bug
        if ~isstruct(Data.Set(s).Amp), continue; end
        
        % Iterate through Amplitudes
        for a = 1:length(Data.Set(s).Amp)
            AmpData = Data.Set(s).Amp(a);
            val = AmpData.Val;
            if val == 0, continue; end
            
            % Create Field Name (e.g. A_5p0 for 5.0uA) - No Tolerance, Exact Match
            fName = sprintf('A_%.1f', val); 
            fName = strrep(fName, '.', 'p');
            
            % Initialize if new amplitude
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_All = [];
                Pooled.(fName).Seq_All = [];
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            % Append Data
            if ~isempty(AmpData.Sim)
                Pooled.(fName).Sim_All = [Pooled.(fName).Sim_All; AmpData.Sim]; 
            end
            if ~isempty(AmpData.Seq)
                Pooled.(fName).Seq_All = [Pooled.(fName).Seq_All; AmpData.Seq]; 
            end
        end
    end
end

All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);

if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Found %d unique amplitudes.\n', length(All_Amps));

%% ================= 2. CALCULATE STATS & PREPARE VECTORS =================
% We prepare vectors for the summary line plot
Stats_Amp = [];
Stats_Mean_Sim = []; Stats_SEM_Sim = [];
Stats_Mean_Seq = []; Stats_SEM_Seq = [];
Stats_PVal = [];

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    s1 = Pooled.(fName).Sim_All;
    s2 = Pooled.(fName).Seq_All;
    
    % Calc Means & SEM
    mu1 = mean(s1); sem1 = std(s1)/sqrt(length(s1));
    mu2 = mean(s2); sem2 = std(s2)/sqrt(length(s2));
    
    % Rank Sum Test
    if ~isempty(s1) && ~isempty(s2)
        p = ranksum(s1, s2);
    else
        p = NaN;
    end
    
    % Store for Line Plot
    Stats_Amp = [Stats_Amp, curr_amp];
    Stats_Mean_Sim = [Stats_Mean_Sim, mu1]; Stats_SEM_Sim = [Stats_SEM_Sim, sem1];
    Stats_Mean_Seq = [Stats_Mean_Seq, mu2]; Stats_SEM_Seq = [Stats_SEM_Seq, sem2];
    Stats_PVal = [Stats_PVal, p];
end

%% ================= 3. PLOT 1: INDIVIDUAL BOX PLOTS (One Per Amplitude) =================
fprintf('Generating Individual Box Plots per Amplitude...\n');

% Define figure size for individual plots
fig_pos_single = [200 200 400 500]; 

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    s1 = Pooled.(fName).Sim_All;
    s2 = Pooled.(fName).Seq_All;
    
    % Skip if no data for this amplitude
    if isempty(s1) && isempty(s2), continue; end
    
    % 1. Create a NEW figure for THIS amplitude
    figName = sprintf('Duration_BoxPlot_%.1fuA', curr_amp);
    figure('Color','w', 'Position', fig_pos_single, 'Name', figName); hold on;
    
    % 2. Prepare Data for Boxplot (Concatenate & Group)
    data_plot = [s1; s2];
    groups    = [ones(size(s1)); 2*ones(size(s2))];
    
    % Define Labels dynamically based on what data exists
    labels = {};
    if ~isempty(s1), labels{end+1} = 'Sim'; end
    if ~isempty(s2), labels{end+1} = 'Seq'; end
    
    % 3. Plot (Updated with 'Notch', 'on')
    if ~isempty(data_plot)
         % [MODIFIED LINE] Added 'Notch', 'on'
         boxplot(data_plot, groups, 'Labels', labels, 'Widths', 0.5, 'Symbol', 'o', 'Notch', 'on');
         
         % 4. Color the boxes (using 'patch' to find box handles)
         h = findobj(gca,'Tag','Box');
         % Note: findobj returns handles in reverse order (Seq is usually first)
         
         if length(h) >= 1 && ~isempty(s2)
             patch(get(h(1),'XData'), get(h(1),'YData'), color_seq, 'FaceAlpha', 0.5); 
         end
         if length(h) >= 2 && ~isempty(s1)
             patch(get(h(2),'XData'), get(h(2),'YData'), color_sim, 'FaceAlpha', 0.5); 
         elseif length(h) == 1 && isempty(s2) && ~isempty(s1)
             patch(get(h(1),'XData'), get(h(1),'YData'), color_sim, 'FaceAlpha', 0.5);
         end
    end
    
    % 5. Add Significance Stars (Check P-Value calculated in Section 2)
    p = Stats_PVal(i);
    if ~isnan(p) && p < 0.05
        % Calculate Y-position for star (slightly above max data point)
        yl = ylim;
        y_star = yl(2) * 0.9; 
        
        txt = '';
        if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        
        % Place text in the center
        text(1.5, y_star, txt, 'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    % 6. Styling
    ylabel('Duration (ms)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Response Duration @ %.1f \\muA', curr_amp), 'FontSize', 14);
    grid on; box off;
    
    if save_figures, saveas(gcf, fullfile(save_dir, [figName '.fig'])); end
end

%% ================= 4. PLOT 2: LINE PLOT (Summary) =================
fprintf('Generating Line Plot...\n');
figNameB = 'Group_Duration_LinePlot';
figure('Color','w', 'Position', [150 150 800 600], 'Name', figNameB); hold on;

% Plot Lines with Error Bars
errorbar(Stats_Amp, Stats_Mean_Sim, Stats_SEM_Sim, '-o', 'Color', color_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'CapSize', 8, 'DisplayName', 'Simultaneous');
errorbar(Stats_Amp, Stats_Mean_Seq, Stats_SEM_Seq, '-s', 'Color', color_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'CapSize', 8, 'DisplayName', 'Sequential');

% Add Stars for Significance
for i = 1:length(Stats_Amp)
    p = Stats_PVal(i);
    if isnan(p), continue; end
    
    txt = '';
    if p < 0.001, txt = '***';
    elseif p < 0.01, txt = '**';
    elseif p < 0.05, txt = '*';
    end
    
    if ~isempty(txt)
        y_top = max(Stats_Mean_Sim(i)+Stats_SEM_Sim(i), Stats_Mean_Seq(i)+Stats_SEM_Seq(i));
        text(Stats_Amp(i), y_top + (y_top*0.05), txt, 'FontSize', 16, ...
             'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

xlabel('Amplitude (\muA)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean Duration (ms) \pm SEM', 'FontSize', 12, 'FontWeight', 'bold');
title('Effect of Sequential Stimulation on Duration', 'FontSize', 14);
grid on; legend('Location', 'northwest');
set(gca, 'XTick', Stats_Amp);
xlim([min(Stats_Amp)-0.5, max(Stats_Amp)+0.5]);

if save_figures, saveas(gcf, fullfile(save_dir, [figNameB '.fig'])); end

fprintf('\n>>> Group Analysis Complete.\n');