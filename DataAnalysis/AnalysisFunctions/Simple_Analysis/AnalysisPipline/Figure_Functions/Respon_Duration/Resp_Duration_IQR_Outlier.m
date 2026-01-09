%% ============================================================
%   GROUP DURATION ANALYSIS: Pooled Population (Refined)
%   - Input: Result_Duration_*.mat files
%   - Features: IQR Filtering, Effect Size Ribbon, Binned Line Plots
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

% Hard upper limit to catch obvious artifacts (before statistical filtering)
filter_time_hard_limit = 15; 

%% ================= 1. POOL DATA (With IQR Filter) =================
fprintf('Pooling data from %d files...\n', length(file_paths));
Pooled = struct();
All_Amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), warning('Missing: %s', file_paths{f}); continue; end
    
    D = load(file_paths{f});
    if ~isfield(D, 'DurationData'), warning('No DurationData in %s', file_paths{f}); continue; end
    Data = D.DurationData;
    
    for s = 1:length(Data.Set)
        if isempty(Data.Set(s).Label) || ~isfield(Data.Set(s), 'Amp'), continue; end
        if ~isstruct(Data.Set(s).Amp), continue; end
        
        for a = 1:length(Data.Set(s).Amp)
            AmpData = Data.Set(s).Amp(a);
            val = AmpData.Val;
            if val == 0, continue; end
            
            fName = sprintf('A_%.1f', val); fName = strrep(fName, '.', 'p');
            
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_All = [];
                Pooled.(fName).Seq_All = [];
                Pooled.(fName).Experiment_Diffs = []; 
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            % --- 1. Apply Hard Filter First ---
            curr_sim = AmpData.Sim; curr_sim(curr_sim > filter_time_hard_limit) = [];
            curr_seq = AmpData.Seq; curr_seq(curr_seq > filter_time_hard_limit) = [];
            
            % --- 2. Apply IQR Outlier Filter ---
            if ~isempty(curr_sim)
                iqr_val = iqr(curr_sim); q1 = prctile(curr_sim, 25); q3 = prctile(curr_sim, 75);
                lower_b = q1 - 1.5*iqr_val; upper_b = q3 + 1.5*iqr_val;
                % lower_b = q1 - 1*iqr_val; upper_b = q3 + 1*iqr_val;
                curr_sim(curr_sim < lower_b | curr_sim > upper_b) = [];
            end
            if ~isempty(curr_seq)
                iqr_val = iqr(curr_seq); q1 = prctile(curr_seq, 25); q3 = prctile(curr_seq, 75);
                lower_b = q1 - 1.5*iqr_val; upper_b = q3 + 1.5*iqr_val;
                % lower_b = q1 - 1*iqr_val; upper_b = q3 + 1*iqr_val;
                curr_seq(curr_seq < lower_b | curr_seq > upper_b) = [];
            end
            
            % --- 3. Store Data ---
            Pooled.(fName).Sim_All = [Pooled.(fName).Sim_All; curr_sim];
            Pooled.(fName).Seq_All = [Pooled.(fName).Seq_All; curr_seq];
        end
    end
end
All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Found %d unique amplitudes.\n', length(All_Amps));

%% ================= 2. CALCULATE STATS & PREPARE VECTORS =================
Stats_Amp = [];
Stats_Mean_Sim = []; Stats_SEM_Sim = [];
Stats_Mean_Seq = []; Stats_SEM_Seq = [];
Stats_PVal = [];

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    s1 = Pooled.(fName).Sim_All;
    s2 = Pooled.(fName).Seq_All;
    
    mu1 = mean(s1); sem1 = std(s1)/sqrt(length(s1));
    mu2 = mean(s2); sem2 = std(s2)/sqrt(length(s2));
    
    if ~isempty(s1) && ~isempty(s2), p = ranksum(s1, s2); else, p = NaN; end
    
    Stats_Amp = [Stats_Amp, curr_amp];
    Stats_Mean_Sim = [Stats_Mean_Sim, mu1]; Stats_SEM_Sim = [Stats_SEM_Sim, sem1];
    Stats_Mean_Seq = [Stats_Mean_Seq, mu2]; Stats_SEM_Seq = [Stats_SEM_Seq, sem2];
    Stats_PVal = [Stats_PVal, p];
end

%% ================= 3. PLOT 1: INDIVIDUAL BOX PLOTS (One Per Amplitude) =================
fprintf('Generating Individual Box Plots per Amplitude...\n');
fig_pos_single = [200 200 400 500]; 
for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    s1 = Pooled.(fName).Sim_All; s2 = Pooled.(fName).Seq_All;
    if isempty(s1) && isempty(s2), continue; end
    
    figName = sprintf('Duration_BoxPlot_%.1fuA', curr_amp);
    figure('Color','w', 'Position', fig_pos_single, 'Name', figName); hold on;
    
    data_plot = [s1; s2]; groups = [ones(size(s1)); 2*ones(size(s2))];
    labels = {}; if ~isempty(s1), labels{end+1} = 'Sim'; end; if ~isempty(s2), labels{end+1} = 'Seq'; end
    
    if ~isempty(data_plot)
         boxplot(data_plot, groups, 'Labels', labels, 'Widths', 0.5, 'Symbol', 'o', 'Notch', 'on');
         h = findobj(gca,'Tag','Box');
         if length(h) >= 1 && ~isempty(s2), patch(get(h(1),'XData'), get(h(1),'YData'), color_seq, 'FaceAlpha', 0.5); end
         if length(h) >= 2 && ~isempty(s1), patch(get(h(2),'XData'), get(h(2),'YData'), color_sim, 'FaceAlpha', 0.5); 
         elseif length(h) == 1 && isempty(s2) && ~isempty(s1), patch(get(h(1),'XData'), get(h(1),'YData'), color_sim, 'FaceAlpha', 0.5); end
    end
    
    p = Stats_PVal(i);
    if ~isnan(p) && p < 0.05
        yl = ylim; y_star = yl(2) * 0.9; 
        txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        text(1.5, y_star, txt, 'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    ylabel('Duration (ms)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Response Duration @ %.1f \\muA', curr_amp), 'FontSize', 14);
    grid on; box off;
    if save_figures, saveas(gcf, fullfile(save_dir, [figName '.fig'])); end
end

%% ================= 4. PLOT 2: LINE PLOT (Summary) =================
fprintf('Generating Line Plot...\n');
figNameB = 'Group_Duration_LinePlot';
figure('Color','w', 'Position', [150 150 800 600], 'Name', figNameB); hold on;
errorbar(Stats_Amp, Stats_Mean_Sim, Stats_SEM_Sim, '-o', 'Color', color_sim, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 8, 'DisplayName', 'Simultaneous');
errorbar(Stats_Amp, Stats_Mean_Seq, Stats_SEM_Seq, '-s', 'Color', color_seq, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 8, 'DisplayName', 'Sequential');

for i = 1:length(Stats_Amp)
    p = Stats_PVal(i);
    if isnan(p), continue; end
    txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
    if ~isempty(txt)
        y_top = max(Stats_Mean_Sim(i)+Stats_SEM_Sim(i), Stats_Mean_Seq(i)+Stats_SEM_Seq(i));
        text(Stats_Amp(i), y_top + (y_top*0.05), txt, 'FontSize', 16, 'HorizontalAlignment', 'center', 'Color', 'k');
    end
end
xlabel('Amplitude (\muA)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean Duration (ms) \pm SEM', 'FontSize', 12, 'FontWeight', 'bold');
title('Effect of Sequential Stimulation on Duration', 'FontSize', 14);
grid on; legend('Location', 'northwest');
set(gca, 'XTick', Stats_Amp); xlim([min(Stats_Amp)-0.5, max(Stats_Amp)+0.5]);
if save_figures, saveas(gcf, fullfile(save_dir, [figNameB '.fig'])); end

%% ================= 5. PLOT 3: EFFECT SIZE (Shaded Ribbon, No Dots) =================
% Shows "Seq - Sim" difference with a clean shaded error band.
fprintf('Generating Effect Size Plot (Clean Ribbon)...\n');
figure('Color','w', 'Position', [200 200 600 500], 'Name', 'Effect_Size_Plot'); hold on;

% 1. Calculate Statistics
diff_mean = Stats_Mean_Seq - Stats_Mean_Sim;
diff_sem  = sqrt(Stats_SEM_Seq.^2 + Stats_SEM_Sim.^2);

% 2. Create Shaded Error Ribbon (Polygon)
x_conf = [Stats_Amp, fliplr(Stats_Amp)];
y_conf = [diff_mean + diff_sem, fliplr(diff_mean - diff_sem)];
fill(x_conf, y_conf, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 

% 3. Plot Main Trend Line (Black)
plot(Stats_Amp, diff_mean, '-o', 'Color', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% 4. Style
yline(0, '--k', 'Alpha', 0.5); 
ylabel('\Delta Duration (Seq - Sim) [ms]', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Amplitude (\muA)', 'FontSize', 12, 'FontWeight', 'bold');
title('Added Duration by Sequence (Mean \pm SEM)', 'FontSize', 14);
grid on; box off;
xlim([min(Stats_Amp)-0.5, max(Stats_Amp)+0.5]); set(gca, 'XTick', Stats_Amp);
if save_figures, saveas(gcf, fullfile(save_dir, 'Group_Effect_Size.fig')); end

%% ================= 6. PLOT 4: BINNED AMPLITUDES (Line Plot) =================
% [IMPROVED] Replaces Swarm Chart with a simple Line Plot (Mean +/- SEM)
fprintf('Generating Binned Line Plot...\n');

% 1. Data Aggregation
bin_defs = { [1, 3.5],   [3.9, 6.5],    [7.5, 12] }; 
bin_names = {'Low (1-3\muA)', 'Med (4-6\muA)', 'High (8-10\muA)'};
BinStats = [];

for b = 1:length(bin_defs)
    s1 = []; s2 = [];
    fn = fieldnames(Pooled);
    for k = 1:length(fn)
        this_amp = Pooled.(fn{k}).Val;
        if this_amp >= bin_defs{b}(1) && this_amp <= bin_defs{b}(2)
            s1 = [s1; Pooled.(fn{k}).Sim_All];
            s2 = [s2; Pooled.(fn{k}).Seq_All];
        end
    end
    
    BinStats(b).Mean_Sim = mean(s1); BinStats(b).SEM_Sim = std(s1)/sqrt(length(s1));
    BinStats(b).Mean_Seq = mean(s2); BinStats(b).SEM_Seq = std(s2)/sqrt(length(s2));
    if ~isempty(s1) && ~isempty(s2), BinStats(b).P = ranksum(s1, s2); else, BinStats(b).P = NaN; end
end

% 2. Plot
figure('Color','w', 'Position', [300 300 600 500], 'Name', 'Binned_LinePlot'); hold on;
b_sim = [BinStats.Mean_Sim]; sem_sim = [BinStats.SEM_Sim];
b_seq = [BinStats.Mean_Seq]; sem_seq = [BinStats.SEM_Seq];
x_pts = 1:3;

errorbar(x_pts, b_sim, sem_sim, '-o', 'Color', color_sim, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 10, 'DisplayName', 'Simultaneous');
errorbar(x_pts, b_seq, sem_seq, '-s', 'Color', color_seq, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 10, 'DisplayName', 'Sequential');

% 3. Statistics
for b = 1:3
    p = BinStats(b).P;
    if isnan(p), continue; end
    txt = ''; if p<0.001, txt='***'; elseif p<0.01, txt='**'; elseif p<0.05, txt='*'; end
    if ~isempty(txt)
        y_top = max(b_sim(b)+sem_sim(b), b_seq(b)+sem_seq(b));
        text(b, y_top * 1.05, txt, 'FontSize', 18, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

% 4. Style
ylabel('Mean Duration (ms)', 'FontSize', 12, 'FontWeight', 'bold');
title('Response Duration (Binned)', 'FontSize', 14);
set(gca, 'XTick', 1:3, 'XTickLabel', bin_names);
legend('Location', 'northwest');
grid on; box off; xlim([0.5 3.5]);

if save_figures, saveas(gcf, fullfile(save_dir, 'Group_Binned_Line.fig')); end
fprintf('\n>>> Additional Summary Plots Generated.\n');