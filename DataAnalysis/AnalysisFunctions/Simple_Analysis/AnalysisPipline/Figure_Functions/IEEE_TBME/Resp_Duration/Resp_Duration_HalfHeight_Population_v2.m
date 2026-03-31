%% ============================================================
%   GROUP DURATION ANALYSIS: Pooled Population (Stim-Pair Average)
%   - Input: Result_Duration_*.mat files (Using 50% FWHM metric)
%   - Pooling Logic: Calculates the median duration per Stimulation Pair (Set).
%   - Outlier Rejection: 1.5 * IQR filter applied to channels before pooling.
%   - Bins: "Operating Regimes" (Low: 1-3µA, Mid: 4-6µA, High: 8-10µA)
%   - Style: IEEE TBME Standard (Arial 9pt, clean B&W)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% List your single-dataset result files here
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX010_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX010_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX010_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX010_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX010_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX010_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX010_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX011_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX011_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX011_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX011_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX011_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX011_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX012_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX012_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX012_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_HalfHeight/Result_Duration_HalfHeight_DX016_Xia_Exp1_Seq_Full_4.mat';
};

% Plot Settings
save_figures = false;
save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Duration_Grouped';

% --- IEEE Format & Stats Settings ---
font_name = 'Arial';
font_size = 9;
color_sim = 'k'; 
color_seq = 'k'; 
IQR_multiplier = 1.5; % [MODIFIED] Defines how strictly to reject outliers
box_whisker    = 0.75; % [MODIFIED] Shrinks the whiskers on the box plot

%% ================= 2. POOL DATA (STIM-PAIR AVERAGE & IQR) =================
fprintf('Pooling data from %d files with %.1f IQR Filter...\n', length(file_paths), IQR_multiplier);
Pooled = struct();
All_Amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), warning('Missing: %s', file_paths{f}); continue; end
    D = load(file_paths{f});
    if ~isfield(D, 'DurationData'), warning('No DurationData in %s', file_paths{f}); continue; end
    Data = D.DurationData;
    
    for s = 1:length(Data.Set)
        if isempty(Data.Set(s).Label) || ~isfield(Data.Set(s), 'Amp'), continue; end
        
        for a = 1:length(Data.Set(s).Amp)
            AmpData = Data.Set(s).Amp(a); val = AmpData.Val;
            if val == 0, continue; end
            
            fName = sprintf('A_%.1f', val); fName = strrep(fName, '.', 'p');
            
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_SetMedians = []; 
                Pooled.(fName).Seq_SetMedians = [];
                All_Amps = [All_Amps, val]; 
            end
            
            curr_sim_channels = AmpData.Sim;
            curr_seq_channels = AmpData.Seq;
            
            % [MODIFIED] Apply IQR Filter to Simultaneous
            if ~isempty(curr_sim_channels)
                iqr_v = iqr(curr_sim_channels); q1 = prctile(curr_sim_channels, 25); q3 = prctile(curr_sim_channels, 75);
                curr_sim_channels(curr_sim_channels < (q1 - IQR_multiplier*iqr_v) | curr_sim_channels > (q3 + IQR_multiplier*iqr_v)) = [];
                set_median_sim = median(curr_sim_channels, 'omitnan');
                if ~isnan(set_median_sim), Pooled.(fName).Sim_SetMedians = [Pooled.(fName).Sim_SetMedians; set_median_sim]; end
            end
            
            % [MODIFIED] Apply IQR Filter to Sequential
            if ~isempty(curr_seq_channels)
                iqr_v = iqr(curr_seq_channels); q1 = prctile(curr_seq_channels, 25); q3 = prctile(curr_seq_channels, 75);
                curr_seq_channels(curr_seq_channels < (q1 - IQR_multiplier*iqr_v) | curr_seq_channels > (q3 + IQR_multiplier*iqr_v)) = [];
                set_median_seq = median(curr_seq_channels, 'omitnan');
                if ~isnan(set_median_seq), Pooled.(fName).Seq_SetMedians = [Pooled.(fName).Seq_SetMedians; set_median_seq]; end
            end
        end
    end
end
All_Amps = sort(unique(All_Amps));
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ================= 3. AGGREGATE INTO BINS & CALCULATE STATS =================
% [MODIFIED] Moved binning logic up so all plots use the exact same binned data
Bins.Low.Sim = []; Bins.Low.Seq = [];   % Peri-Threshold (1-3 uA)
Bins.Mid.Sim = []; Bins.Mid.Seq = [];   % Linear Range (4-6 uA)
Bins.High.Sim = []; Bins.High.Seq = []; % Saturation (8-10 uA)

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    s_sim = Pooled.(fName).Sim_SetMedians;
    s_seq = Pooled.(fName).Seq_SetMedians;
    
    if curr_amp <= 3.5
        Bins.Low.Sim = [Bins.Low.Sim; s_sim]; Bins.Low.Seq = [Bins.Low.Seq; s_seq];
    elseif curr_amp >= 3.9 && curr_amp <= 6.5
        Bins.Mid.Sim = [Bins.Mid.Sim; s_sim]; Bins.Mid.Seq = [Bins.Mid.Seq; s_seq];
    elseif curr_amp >= 7.5
        Bins.High.Sim = [Bins.High.Sim; s_sim]; Bins.High.Seq = [Bins.High.Seq; s_seq];
    end
end

% Calculate Binned Stats (for Line and Ribbon plots)
BinStats = struct();
fields = {'Low', 'Mid', 'High'};
for k = 1:3
    f = fields{k};
    s_sim = Bins.(f).Sim; s_seq = Bins.(f).Seq;
    BinStats(k).Mean_Sim = mean(s_sim); BinStats(k).SEM_Sim = std(s_sim)/sqrt(length(s_sim));
    BinStats(k).Mean_Seq = mean(s_seq); BinStats(k).SEM_Seq = std(s_seq)/sqrt(length(s_seq));
    if ~isempty(s_sim) && ~isempty(s_seq), BinStats(k).P = ranksum(s_sim, s_seq); else, BinStats(k).P = NaN; end
end

%% ================= 4. PLOT 1: BINNED BOX PLOTS =================
fprintf('Generating Binned Box Plot...\n');
BinLabels = {'1-3 \muA', '4-6 \muA', '8-10 \muA'};
pos_ctr = 1; 

figure('Units', 'centimeters', 'Position', [5, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Duration_Binned_BoxPlot'); hold on;
warning('off', 'all'); 

for k = 1:length(fields)
    f = fields{k};
    d1 = Bins.(f).Sim; d2 = Bins.(f).Seq;
    pos1 = pos_ctr - 0.15; pos2 = pos_ctr + 0.15;
    
    % Simultaneous Box (White)
    if ~isempty(d1)
        boxplot(d1, 'Positions', pos1, 'Widths', 0.22, 'Colors', 'k', 'Symbol', 'k.', 'Whisker', box_whisker);
        h = findobj(gca,'Tag','Box'); patch(get(h(1),'XData'), get(h(1),'YData'), 'w', 'FaceAlpha', 1, 'EdgeColor', 'k');
    end
    
    % Sequential Box (Grey)
    if ~isempty(d2)
        boxplot(d2, 'Positions', pos2, 'Widths', 0.22, 'Colors', 'k', 'Symbol', 'k.', 'Whisker', box_whisker);
        h = findobj(gca,'Tag','Box'); patch(get(h(1),'XData'), get(h(1),'YData'), [0.8 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'k');
    end
    
    % [MODIFIED] Force the Median lines to be visible on top of the patch
    hMed = findobj(gca, 'Tag', 'Median');
    set(hMed, 'Color', 'k', 'LineWidth', 1.5);
    % Loop through each median line individually to avoid the parent error
    for m = 1:length(hMed)
        uistack(hMed(m), 'top');
    end
    
    % Stats Star
    p = BinStats(k).P;
    if p < 0.05
        txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        y_max = max([max(d1), max(d2)]); y_star = y_max * 1.05; 
        plot([pos1, pos2], [y_star-0.5, y_star-0.5], '-k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        text(pos_ctr, y_star + 0.5, txt, 'FontSize', font_size, 'HorizontalAlignment', 'center', 'FontName', font_name);
    end
    pos_ctr = pos_ctr + 1;
end
warning('on', 'all');

ylabel('Duration (ms)', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels, 'FontSize', font_size, 'FontName', font_name, 'TickDir', 'out', 'LineWidth', 1.0);
xlim([0.5, 3.5]); ylim([0 max([Bins.High.Seq; Bins.High.Sim]) * 1.25]); box off; axis square;
h1 = plot(NaN,NaN, 's', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
h2 = plot(NaN,NaN, 's', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
legend([h1, h2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', font_size, 'FontName', font_name);

%% ================= 5. PLOT 2: LINE PLOT (Binned Tuning Curve) =================
fprintf('Generating Binned Line Plot...\n');
figure('Units', 'centimeters', 'Position', [10, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Duration_LinePlot'); hold on;

% [MODIFIED] Plotting exactly 3 points (x_pts = 1, 2, 3) instead of raw amplitudes
x_pts = 1:3;
b_sim = [BinStats.Mean_Sim]; sem_sim = [BinStats.SEM_Sim];
b_seq = [BinStats.Mean_Seq]; sem_seq = [BinStats.SEM_Seq];

e1 = errorbar(x_pts, b_sim, sem_sim, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'CapSize', 3, 'MarkerSize', 5);
e2 = errorbar(x_pts, b_seq, sem_seq, '--ks', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'CapSize', 3, 'MarkerSize', 5);

for b = 1:3
    p = BinStats(b).P;
    if isnan(p) || p >= 0.05, continue; end
    txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
    y_top = max(b_sim(b)+sem_sim(b), b_seq(b)+sem_seq(b));
    text(b, y_top + 0.5, txt, 'FontSize', font_size, 'HorizontalAlignment', 'center', 'FontName', font_name);
end

ylabel('Mean Duration (ms)', 'FontSize', font_size, 'FontName', font_name);
legend([e1, e2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels, 'FontSize', font_size, 'FontName', font_name, 'LineWidth', 1.0, 'TickDir', 'out');
xlim([0.5, 3.5]); box off; axis square;

%% ================= 6. PLOT 3: EFFECT SIZE (Binned Shaded Ribbon) =================
fprintf('Generating Binned Effect Size Plot...\n');
figure('Units', 'centimeters', 'Position', [15, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Effect_Size_Plot'); hold on;

% [MODIFIED] Using the binned differences
diff_mean = b_seq - b_sim;
diff_sem  = sqrt(sem_seq.^2 + sem_sim.^2);

x_conf = [x_pts, fliplr(x_pts)];
y_conf = [diff_mean + diff_sem, fliplr(diff_mean - diff_sem)];

fill(x_conf, y_conf, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
plot(x_pts, diff_mean, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
yline(0, '--k', 'Alpha', 0.5, 'HandleVisibility', 'off'); 

ylabel('\Delta Duration (Seq - Sim) [ms]', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels, 'FontSize', font_size, 'FontName', font_name, 'TickDir', 'out', 'LineWidth', 1.0);
xlim([0.5, 3.5]); box off; axis square;

%% ================= 7. SAVE & PRINT SUMMARY =================
if save_figures
    exportgraphics(figure(1), fullfile(save_dir, 'Duration_Binned_BoxPlot.tiff'), 'ContentType', 'vector');
    exportgraphics(figure(2), fullfile(save_dir, 'Duration_Binned_LinePlot.tiff'), 'ContentType', 'vector');
    exportgraphics(figure(3), fullfile(save_dir, 'Duration_Binned_EffectSize.tiff'), 'ContentType', 'vector');
end

fprintf('\n================================================================\n');
fprintf('             STATISTICAL SUMMARY (BINNED AMPLITUDES)            \n');
fprintf('================================================================\n');
fprintf('%-15s | %-20s | %-20s | %-10s\n', 'Range', 'Sim (Mean+/-SEM)', 'Seq (Mean+/-SEM)', 'P-Value');
fprintf('----------------------------------------------------------------\n');
for b = 1:length(fields)
    fprintf('%-15s | %6.2f +/- %.2f      | %6.2f +/- %.2f      | %.5f\n', ...
        BinLabels{b}, BinStats(b).Mean_Sim, BinStats(b).SEM_Sim, BinStats(b).Mean_Seq, BinStats(b).SEM_Seq, BinStats(b).P);
end
fprintf('================================================================\n');