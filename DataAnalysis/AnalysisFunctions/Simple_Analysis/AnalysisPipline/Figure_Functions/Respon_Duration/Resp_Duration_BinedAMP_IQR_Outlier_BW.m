%% ============================================================
%   GROUP DURATION ANALYSIS: Pooled Population (Refined)
%   - Input: Result_Duration_*.mat files
%   - Features: IQR Filtering, Effect Size Ribbon, Binned Line Plots
%   - Style: Black & White (Publication Ready - Dashed/Solid Lines)
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
save_figures = true;
save_dir     = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE_EMBC/Figures/5_Respond_Duration';
fig_pos      = [100 100 800 600];
% --- Black and White Settings ---
color_sim    = 'k'; 
color_seq    = 'k'; 
filter_time_hard_limit = 14; 
IQR_out_parm = 1.5;
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
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            % --- 1. Apply Hard Filter First ---
            curr_sim = AmpData.Sim; curr_sim(curr_sim > filter_time_hard_limit) = [];
            curr_seq = AmpData.Seq; curr_seq(curr_seq > filter_time_hard_limit) = [];
            
            % --- 2. Apply IQR Outlier Filter ---
            if ~isempty(curr_sim)
                iqr_val = iqr(curr_sim); q1 = prctile(curr_sim, 25); q3 = prctile(curr_sim, 75);
                lower_b = q1 - IQR_out_parm*iqr_val; upper_b = q3 + IQR_out_parm*iqr_val;
                curr_sim(curr_sim < lower_b | curr_sim > upper_b) = [];
            end
            if ~isempty(curr_seq)
                iqr_val = iqr(curr_seq); q1 = prctile(curr_seq, 25); q3 = prctile(curr_seq, 75);
                lower_b = q1 - IQR_out_parm*iqr_val; upper_b = q3 + IQR_out_parm*iqr_val;
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

%% ================= 2. CALCULATE STATS (For Line Plots) =================
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

%% ================= 3. PLOT 1: BINNED BOX PLOTS (Low/Mid/High) =================
fprintf('Generating Binned Box Plot (Low/Mid/High)...\n');

% 1. Aggregate Data into Bins
Bins.Low.Sim = []; Bins.Low.Seq = [];   % 1-3 uA
Bins.Mid.Sim = []; Bins.Mid.Seq = [];   % 4-6 uA
Bins.High.Sim = []; Bins.High.Seq = []; % 8-10 uA

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    s_sim = Pooled.(fName).Sim_All;
    s_seq = Pooled.(fName).Seq_All;
    
    if curr_amp <= 3.1
        Bins.Low.Sim = [Bins.Low.Sim; s_sim];
        Bins.Low.Seq = [Bins.Low.Seq; s_seq];
    elseif curr_amp >= 3.9 && curr_amp <= 6.1
        Bins.Mid.Sim = [Bins.Mid.Sim; s_sim];
        Bins.Mid.Seq = [Bins.Mid.Seq; s_seq];
    elseif curr_amp >= 7.9
        Bins.High.Sim = [Bins.High.Sim; s_sim];
        Bins.High.Seq = [Bins.High.Seq; s_seq];
    end
end

% 2. Prepare Data for Boxplot (Grouping)
% We create a grouped plot: [LowSim, LowSeq, MidSim, MidSeq, HighSim, HighSeq]
% But to use 'boxplot' easily with grouping, we construct a data vector.

PlotData = [];
PlotGroups = [];
BinLabels = {'Low (1-3µA)', 'Medium (4-6µA)', 'High (8-10µA)'};
BinPvals = [];

fields = {'Low', 'Mid', 'High'};
pos_ctr = 1; 

figure('Color','w', 'Position', [100 100 800 500], 'Name', 'Duration_Binned_BoxPlot'); 
hold on;

for k = 1:length(fields)
    f = fields{k};
    d1 = Bins.(f).Sim;
    d2 = Bins.(f).Seq;
    
    % Positions for this bin (Side by Side)
    pos1 = pos_ctr - 0.15;
    pos2 = pos_ctr + 0.15;
    
    % --- Plot Simultaneous Box (White) ---
    if ~isempty(d1)
        % boxplot(d1, 'Positions', pos1, 'Widths', 0.25, 'Symbol', '', 'Colors', 'k', 'Notch', 'on');
        boxplot(d1, 'Positions', pos1, 'Widths', 0.22, 'Colors', 'k', ...
            'Symbol', '', 'Whisker', 0.75, 'Notch', 'on');
        
        % Hack to fill/style: find the last created objects
        h = findobj(gca,'Tag','Box'); 
        patch(get(h(1),'XData'), get(h(1),'YData'), 'w', 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    % --- Plot Sequential Box (Grey) ---
    if ~isempty(d2)
        % boxplot(d2, 'Positions', pos2, 'Widths', 0.25, 'Symbol', '', 'Colors', 'k', 'Notch', 'on');
        boxplot(d2, 'Positions', pos2, 'Widths', 0.22, 'Colors', 'k', ...
            'Symbol', '', 'Whisker', 0.75, 'Notch', 'on');
        h = findobj(gca,'Tag','Box'); 
        patch(get(h(1),'XData'), get(h(1),'YData'), [0.8 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    % --- Statistics (Wilcoxon Rank Sum) ---
    if ~isempty(d1) && ~isempty(d2)
        p = ranksum(d1, d2);
        BinPvals(k) = p;
        if p < 0.05
            txt = ''; 
            if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
            
            % Draw bracket or just text
            y_max = max([max(d1), max(d2)]);
            y_star = y_max + 1.0; 
            text(pos_ctr, y_star, txt, 'FontSize', 22, 'HorizontalAlignment', 'center', ...
                'FontName', 'Times New Roman', 'FontWeight', 'bold');
            
            % Optional: Draw a line connecting them?
            plot([pos1, pos2], [y_star-0.5, y_star-0.5], '-k', 'LineWidth', 1);
        end
    end
    
    pos_ctr = pos_ctr + 1;
end

% 3. Formatting
ylabel('Duration (ms)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Amplitude Range', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% title('Response Duration by Amplitude Group', 'FontSize', 16, 'FontName', 'Times New Roman');

set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'out');
xlim([0.5, 3.5]);
ylim([0 18]); % Adjust based on data
box off;
axis square;
% Add Dummy Legend
h1 = plot(NaN,NaN, 's', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
h2 = plot(NaN,NaN, 's', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
legend([h1, h2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 14, 'FontName', 'Times New Roman');

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'Duration_Binned_BoxPlot.tiff'), 'ContentType', 'vector');
end


%% ================= 4. PLOT 2: LINE PLOT (Summary B&W) =================
fprintf('Generating Line Plot...\n');
figNameB = 'Group_Duration_LinePlot';
figure('Color','w', 'Position', [150 150 800 600], 'Name', figNameB); hold on;
% Sim: Dashed Line ('--o'), Open Circle
errorbar(Stats_Amp, Stats_Mean_Sim, Stats_SEM_Sim, '--o', 'Color', color_sim, ...
    'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 10, 'DisplayName', 'Simultaneous', 'MarkerSize', 8);
% Seq: Solid Line ('-s'), Filled Square
errorbar(Stats_Amp, Stats_Mean_Seq, Stats_SEM_Seq, '-s', 'Color', color_seq, ...
    'LineWidth', 2, 'MarkerFaceColor', 'k', 'CapSize', 10, 'DisplayName', 'Sequential', 'MarkerSize', 8);

for i = 1:length(Stats_Amp)
    p = Stats_PVal(i);
    if isnan(p), continue; end
    txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
    if ~isempty(txt)
        y_top = max(Stats_Mean_Sim(i)+Stats_SEM_Sim(i), Stats_Mean_Seq(i)+Stats_SEM_Seq(i));
        text(Stats_Amp(i), y_top + 0.5, txt, 'FontSize', 20, 'HorizontalAlignment', 'center', 'Color', 'k', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    end
end
xlabel('Amplitude (\muA)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Mean Duration (ms) \pm SEM', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% title('Effect of Sequential Stimulation on Duration', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 16, 'FontName', 'Times New Roman');
set(gca, 'XTick', Stats_Amp, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'out');
xlim([min(Stats_Amp)-0.5, max(Stats_Amp)+0.5]);
axis square;
ylim([0,12]);
if save_figures
    exportgraphics(gcf, fullfile(save_dir, [figNameB '.tiff']),'ContentType', 'vector');
end


%% ================= 5. PLOT 3: EFFECT SIZE (Shaded Ribbon B&W) =================
fprintf('Generating Effect Size Plot (Clean Ribbon)...\n');
figure('Color','w', 'Position', [200 200 600 500], 'Name', 'Effect_Size_Plot'); hold on;

% 1. Calculate Statistics
diff_mean = Stats_Mean_Seq - Stats_Mean_Sim;
diff_sem  = sqrt(Stats_SEM_Seq.^2 + Stats_SEM_Sim.^2);

% 2. Create Shaded Error Ribbon (Light Gray)
x_conf = [Stats_Amp, fliplr(Stats_Amp)];
y_conf = [diff_mean + diff_sem, fliplr(diff_mean - diff_sem)];
fill(x_conf, y_conf, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 

% 3. Plot Main Trend Line (Black, Filled Circle)
plot(Stats_Amp, diff_mean, '-o', 'Color', 'k', 'LineWidth', 4, 'MarkerFaceColor', 'k', 'MarkerSize', 6);

% 4. Style
yline(0, '--k', 'Alpha', 0.5); 
ylabel('\Delta Duration (Seq - Sim) [ms]', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Amplitude (\muA)', 'FontSize', 18, 'FontWeight', 'bold');
% title('Added Duration by Sequence (Mean \pm SEM)', 'FontSize', 16);
box off;
xlim([min(Stats_Amp)-0.5, max(Stats_Amp)+0.5]); set(gca, 'XTick', Stats_Amp);
ylim([-3,3]);
axis square;
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'out');

fig_name = 'Group_Effect_Size.tiff';
if save_figures
    % saveas(gcf, fullfile(save_dir, 'Group_Binned_Line.fig')); 
    exportgraphics(gcf, fullfile(save_dir, fig_name),'ContentType', 'vector');
end

%% ================= 6. PLOT 4: BINNED AMPLITUDES (Line Plot B&W) =================
fprintf('Generating Binned Line Plot...\n');

% 1. Data Aggregation
bin_defs = { [1, 3.5],   [3.9, 6.5],    [7.5, 12] }; 
bin_names = {'Low (1-3\muA)', 'Medium (4-6\muA)', 'High (8-10\muA)'};
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

% [MODIFIED] Sim: Dashed Line ('--o'), Open Circle
errorbar(x_pts, b_sim, sem_sim, '--o', 'Color', color_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'CapSize', 10, 'DisplayName', 'Simultaneous');

% [MODIFIED] Seq: Solid Line ('-s'), Filled Square
errorbar(x_pts, b_seq, sem_seq, '-s', 'Color', color_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'k', 'CapSize', 10, 'DisplayName', 'Sequential');

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
ylabel('Mean Duration (ms)', 'FontSize', 18, 'FontWeight', 'bold');
xlabel('Amplitude Range', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% title('Response Duration (Binned)', 'FontSize', 16);
set(gca, 'XTick', 1:3, 'XTickLabel', bin_names);
set(gca, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'out');
legend('Location', 'northwest');
box off; xlim([0.5 3.5]);
axis square;

fig_name = 'Group_Binned_Line.tiff';
if save_figures
    % saveas(gcf, fullfile(save_dir, 'Group_Binned_Line.fig')); 
    exportgraphics(gcf, fullfile(save_dir, fig_name),'ContentType', 'vector');
    fprintf('\n>>> Plots saved.\n');
end
