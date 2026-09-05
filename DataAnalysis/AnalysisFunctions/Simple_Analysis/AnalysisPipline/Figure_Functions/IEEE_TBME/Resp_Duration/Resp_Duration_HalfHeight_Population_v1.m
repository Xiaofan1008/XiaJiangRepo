%% ============================================================
%   GROUP DURATION ANALYSIS: Pooled Population (Stim-Pair Average)
%   - Input: Result_Duration_*.mat files (Using 50% FWHM metric)
%   - Pooling Logic: Calculates the median duration per Stimulation Pair (Set).
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

% --- IEEE Format Settings ---
font_name = 'Arial';
font_size = 9;
color_sim = 'k'; 
color_seq = 'k'; 

%% ================= 2. POOL DATA (STIM-PAIR AVERAGE) =================
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
        
        for a = 1:length(Data.Set(s).Amp)
            AmpData = Data.Set(s).Amp(a);
            val = AmpData.Val;
            if val == 0, continue; end
            
            fName = sprintf('A_%.1f', val); fName = strrep(fName, '.', 'p');
            
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_SetMedians = []; % Store 1 value per Set
                Pooled.(fName).Seq_SetMedians = [];
                All_Amps = [All_Amps, val]; 
            end
            
            % --- STIM-PAIR AVERAGE LOGIC ---
            % Instead of pooling all channels, we take the median duration 
            % of all responsive channels for this specific stimulation set.
            curr_sim_channels = AmpData.Sim;
            curr_seq_channels = AmpData.Seq;
            
            % Only add to the pool if there was a valid response
            if ~isempty(curr_sim_channels)
                set_median_sim = median(curr_sim_channels, 'omitnan');
                if ~isnan(set_median_sim)
                    Pooled.(fName).Sim_SetMedians = [Pooled.(fName).Sim_SetMedians; set_median_sim];
                end
            end
            
            if ~isempty(curr_seq_channels)
                set_median_seq = median(curr_seq_channels, 'omitnan');
                if ~isnan(set_median_seq)
                    Pooled.(fName).Seq_SetMedians = [Pooled.(fName).Seq_SetMedians; set_median_seq];
                end
            end
        end
    end
end
All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
fprintf('Found %d unique amplitudes.\n', length(All_Amps));

%% ================= 3. CALCULATE STATS (For Line Plots) =================
Stats_Amp = [];
Stats_Mean_Sim = []; Stats_SEM_Sim = [];
Stats_Mean_Seq = []; Stats_SEM_Seq = [];
Stats_PVal = [];

for i = 1:length(All_Amps)
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    s1 = Pooled.(fName).Sim_SetMedians;
    s2 = Pooled.(fName).Seq_SetMedians;
    
    mu1 = mean(s1); sem1 = std(s1)/sqrt(length(s1));
    mu2 = mean(s2); sem2 = std(s2)/sqrt(length(s2));
    
    if ~isempty(s1) && ~isempty(s2), p = ranksum(s1, s2); else, p = NaN; end
    
    Stats_Amp = [Stats_Amp, curr_amp];
    Stats_Mean_Sim = [Stats_Mean_Sim, mu1]; Stats_SEM_Sim = [Stats_SEM_Sim, sem1];
    Stats_Mean_Seq = [Stats_Mean_Seq, mu2]; Stats_SEM_Seq = [Stats_SEM_Seq, sem2];
    Stats_PVal = [Stats_PVal, p];
end

%% ================= 4. PLOT 1: BINNED BOX PLOTS (Operating Regimes) =================
fprintf('Generating Binned Box Plot (Operating Regimes)...\n');
% 1. Aggregate Data into Bins
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

BinLabels = {'1-3 \muA', '4-6 \muA', '8-10 \muA'};
fields = {'Low', 'Mid', 'High'};
pos_ctr = 1; 

figure('Units', 'centimeters', 'Position', [5, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Duration_Binned_BoxPlot'); 
hold on;
warning('off', 'all'); % Mute boxplot warnings

for k = 1:length(fields)
    f = fields{k};
    d1 = Bins.(f).Sim; d2 = Bins.(f).Seq;
    pos1 = pos_ctr - 0.15; pos2 = pos_ctr + 0.15;
    
    % Simultaneous Box (White)
    if ~isempty(d1)
        boxplot(d1, 'Positions', pos1, 'Widths', 0.22, 'Colors', 'k', 'Symbol', 'k.', 'Whisker', 1.5);
        h = findobj(gca,'Tag','Box'); 
        patch(get(h(1),'XData'), get(h(1),'YData'), 'w', 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1.0);
    end
    
    % Sequential Box (Grey)
    if ~isempty(d2)
        boxplot(d2, 'Positions', pos2, 'Widths', 0.22, 'Colors', 'k', 'Symbol', 'k.', 'Whisker', 1.5);
        h = findobj(gca,'Tag','Box'); 
        patch(get(h(1),'XData'), get(h(1),'YData'), [0.8 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1.0);
    end
    
    % Statistics (Wilcoxon Rank Sum)
    if ~isempty(d1) && ~isempty(d2)
        p = ranksum(d1, d2);
        if p < 0.05
            txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
            y_max = max([max(d1), max(d2)]); y_star = y_max * 1.05; 
            plot([pos1, pos2], [y_star-0.5, y_star-0.5], '-k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
            text(pos_ctr, y_star + 0.5, txt, 'FontSize', font_size, 'HorizontalAlignment', 'center', 'FontName', font_name);
        end
    end
    pos_ctr = pos_ctr + 1;
end
warning('on', 'all');

% IEEE Formatting
ylabel('Duration (ms)', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels, 'FontSize', font_size, 'FontName', font_name, 'TickDir', 'out', 'LineWidth', 1.0);
xlim([0.5, 3.5]); ylim([0 max([Bins.High.Seq; Bins.High.Sim]) * 1.25]); box off; axis square;

% Dummy Legend
h1 = plot(NaN,NaN, 's', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
h2 = plot(NaN,NaN, 's', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
legend([h1, h2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', font_size, 'FontName', font_name);

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'Duration_Binned_BoxPlot.pdf'), 'ContentType', 'vector');
end

%% ================= 5. PLOT 2: LINE PLOT (Tuning Curve) =================
fprintf('Generating Line Plot...\n');
figure('Units', 'centimeters', 'Position', [10, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Duration_LinePlot'); hold on;

% Sim: Solid Line ('-ko'), Filled Circle
e1 = errorbar(Stats_Amp, Stats_Mean_Sim, Stats_SEM_Sim, '-ko', 'Color', color_sim, ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'CapSize', 3, 'MarkerSize', 5);
% Seq: Dashed Line ('--ks'), Open Square
e2 = errorbar(Stats_Amp, Stats_Mean_Seq, Stats_SEM_Seq, '--ks', 'Color', color_seq, ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'CapSize', 3, 'MarkerSize', 5);

for i = 1:length(Stats_Amp)
    p = Stats_PVal(i);
    if isnan(p) || p >= 0.05, continue; end
    txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
    y_top = max(Stats_Mean_Sim(i)+Stats_SEM_Sim(i), Stats_Mean_Seq(i)+Stats_SEM_Seq(i));
    text(Stats_Amp(i), y_top + 0.5, txt, 'FontSize', font_size, 'HorizontalAlignment', 'center', 'FontName', font_name);
end

xlabel('Stimulation Amplitude (\muA)', 'FontSize', font_size, 'FontName', font_name);
ylabel('Mean Duration (ms)', 'FontSize', font_size, 'FontName', font_name);
legend([e1, e2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', font_size, 'FontName', font_name);
set(gca, 'XTick', Stats_Amp, 'FontSize', font_size, 'FontName', font_name, 'LineWidth', 1.0, 'TickDir', 'out');
xlim([min(Stats_Amp)-1, max(Stats_Amp)+1]); box off; axis square;

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'Duration_LinePlot.pdf'),'ContentType', 'vector');
end

%% ================= 6. PLOT 3: EFFECT SIZE (Shaded Ribbon) =================
fprintf('Generating Effect Size Plot...\n');
figure('Units', 'centimeters', 'Position', [15, 5, 8.8, 8.8], 'Color', 'w', 'Name', 'Effect_Size_Plot'); hold on;

diff_mean = Stats_Mean_Seq - Stats_Mean_Sim;
diff_sem  = sqrt(Stats_SEM_Seq.^2 + Stats_SEM_Sim.^2);

x_conf = [Stats_Amp, fliplr(Stats_Amp)];
y_conf = [diff_mean + diff_sem, fliplr(diff_mean - diff_sem)];
fill(x_conf, y_conf, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
plot(Stats_Amp, diff_mean, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'MarkerSize', 5);

yline(0, '--k', 'Alpha', 0.5, 'HandleVisibility', 'off'); 
ylabel('\Delta Duration (Seq - Sim) [ms]', 'FontSize', font_size, 'FontName', font_name);
xlabel('Stimulation Amplitude (\muA)', 'FontSize', font_size, 'FontName', font_name);

box off; axis square;
xlim([min(Stats_Amp)-1, max(Stats_Amp)+1]); set(gca, 'XTick', Stats_Amp, 'FontSize', font_size, 'FontName', font_name, 'TickDir', 'out', 'LineWidth', 1.0);

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'Duration_EffectSize.pdf'),'ContentType', 'vector');
end

%% ================= 7. PRINT STATISTICS TO COMMAND WINDOW =================
fprintf('\n================================================================\n');
fprintf('             STATISTICAL SUMMARY (BINNED AMPLITUDES)            \n');
fprintf('================================================================\n');
fprintf('%-15s | %-20s | %-20s | %-10s\n', 'Range', 'Sim (Mean+/-SEM)', 'Seq (Mean+/-SEM)', 'P-Value');
fprintf('----------------------------------------------------------------\n');
bin_names = {'1-3uA', '4-6uA', '8-10uA'};
for b = 1:length(fields)
    f = fields{b};
    s_sim = Bins.(f).Sim; s_seq = Bins.(f).Seq;
    
    m_sim = mean(s_sim); sem_sim = std(s_sim)/sqrt(length(s_sim));
    m_seq = mean(s_seq); sem_seq = std(s_seq)/sqrt(length(s_seq));
    if ~isempty(s_sim) && ~isempty(s_seq), pval = ranksum(s_sim, s_seq); else, pval = NaN; end
    
    fprintf('%-15s | %6.2f +/- %.2f      | %6.2f +/- %.2f      | %.5f\n', ...
        bin_names{b}, m_sim, sem_sim, m_seq, sem_seq, pval);
end
fprintf('================================================================\n');