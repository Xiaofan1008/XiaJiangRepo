%% ============================================================
%   GROUP DURATION ANALYSIS: Pooled Population (Binned Medians)
%   - Metric: 3*SD Total Duration
%   - Logic: 1.5x IQR filter on raw trials -> Median per Set -> Binned
%   - Output: 3 Standardized Binned Plots (Box, Line, Effect Size)
%   - Style: Black & White (Publication Ready, Times New Roman)
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
% List your single-dataset result files here
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq2_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD/Result_Duration_3SD_DX016_Xia_Exp1_Seq_Full_4.mat';

    
    
};

% Plot Settings
save_figures = false;
save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Duration_3SD';

% --- Black and White Settings ---
color_sim    = 'k'; 
color_seq    = 'k'; 
IQR_out_parm = 1.5; % Strict IQR filter to drop intra-set noise

%% ================= 1. POOL DATA (Median per Set with IQR Filter) =================
fprintf('Pooling data from %d files (Calculating Medians)...\n', length(file_paths));
Pooled = struct();
All_Amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), continue; end
    D = load(file_paths{f});
    if ~isfield(D, 'DurationData'), continue; end
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
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            curr_sim = AmpData.Sim; 
            curr_seq = AmpData.Seq; 
            
            % [MODIFIED] Apply IQR Filter and then take MEDIAN
            if ~isempty(curr_sim)
                iqr_val = iqr(curr_sim); q1 = prctile(curr_sim, 25); q3 = prctile(curr_sim, 75);
                curr_sim(curr_sim < (q1 - IQR_out_parm*iqr_val) | curr_sim > (q3 + IQR_out_parm*iqr_val)) = [];
                med_sim = median(curr_sim, 'omitnan');
                if ~isnan(med_sim), Pooled.(fName).Sim_SetMedians = [Pooled.(fName).Sim_SetMedians; med_sim]; end
            end
            
            if ~isempty(curr_seq)
                iqr_val = iqr(curr_seq); q1 = prctile(curr_seq, 25); q3 = prctile(curr_seq, 75);
                curr_seq(curr_seq < (q1 - IQR_out_parm*iqr_val) | curr_seq > (q3 + IQR_out_parm*iqr_val)) = [];
                med_seq = median(curr_seq, 'omitnan');
                if ~isnan(med_seq), Pooled.(fName).Seq_SetMedians = [Pooled.(fName).Seq_SetMedians; med_seq]; end
            end
        end
    end
end
All_Amps = sort(unique(All_Amps));
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ================= 2. AGGREGATE INTO BINS & CALCULATE STATS =================
% [MODIFIED] Grouping the medians into Operating Regimes immediately
Bins.Low.Sim = []; Bins.Low.Seq = [];   % 1-3 uA
Bins.Mid.Sim = []; Bins.Mid.Seq = [];   % 4-6 uA
Bins.High.Sim = []; Bins.High.Seq = []; % 8-10 uA

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
BinStats = struct();

for k = 1:3
    f = fields{k};
    d1 = Bins.(f).Sim; d2 = Bins.(f).Seq;
    BinStats(k).Mean_Sim = mean(d1); BinStats(k).SEM_Sim = std(d1)/sqrt(length(d1));
    BinStats(k).Mean_Seq = mean(d2); BinStats(k).SEM_Seq = std(d2)/sqrt(length(d2));
    BinStats(k).Med_Sim  = median(d1); BinStats(k).IQR_Sim = iqr(d1);
    BinStats(k).Med_Seq  = median(d2); BinStats(k).IQR_Seq = iqr(d2);
    if ~isempty(d1) && ~isempty(d2), BinStats(k).P = ranksum(d1, d2); else, BinStats(k).P = NaN; end
end

%% ================= 3. PLOT 1: BINNED BOX PLOTS =================
fprintf('Generating Binned Box Plot...\n');
figure('Color','w', 'Position', [100 100 800 500], 'Name', 'Duration_Binned_BoxPlot'); hold on;
pos_ctr = 1; warning('off', 'all');

for k = 1:length(fields)
    f = fields{k};
    d1 = Bins.(f).Sim; d2 = Bins.(f).Seq;
    pos1 = pos_ctr - 0.15; pos2 = pos_ctr + 0.15;
    
    if ~isempty(d1)
        boxplot(d1, 'Positions', pos1, 'Widths', 0.22, 'Colors', 'k', 'Symbol', '', 'Whisker', 0.75, 'Notch', 'on');
        h = findobj(gca,'Tag','Box'); patch(get(h(1),'XData'), get(h(1),'YData'), 'w', 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    if ~isempty(d2)
        boxplot(d2, 'Positions', pos2, 'Widths', 0.22, 'Colors', 'k', 'Symbol', '', 'Whisker', 0.75, 'Notch', 'on');
        h = findobj(gca,'Tag','Box'); patch(get(h(1),'XData'), get(h(1),'YData'), [0.8 0.8 0.8], 'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    p = BinStats(k).P;
    if p < 0.05
        txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        y_max = max([max(d1), max(d2)]); y_star = y_max + 1.0; 
        text(pos_ctr, y_star, txt, 'FontSize', 22, 'HorizontalAlignment', 'center', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
        plot([pos1, pos2], [y_star-0.5, y_star-0.5], '-k', 'LineWidth', 1);
    end
    pos_ctr = pos_ctr + 1;
end
warning('on', 'all');

ylabel('Duration (ms)', 'FontSize', 20, 'FontName', 'Times New Roman');
xlabel('Amplitude Range', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels, 'FontSize', 18, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'out');
xlim([0.5, 3.5]); ylim([0 30]); box off; axis square;
h1 = plot(NaN,NaN, 's', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
h2 = plot(NaN,NaN, 's', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', 'k', 'MarkerSize', 12);
legend([h1, h2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 18, 'FontName', 'Times New Roman');
if save_figures, exportgraphics(gcf, fullfile(save_dir, 'Duration_Binned_BoxPlot.tiff'), 'ContentType', 'vector'); end

%% ================= 4. PLOT 2: LINE PLOT (Binned) =================
fprintf('Generating Binned Line Plot...\n');
figure('Color','w', 'Position', [150 150 800 600], 'Name', 'Group_Duration_LinePlot'); hold on;
x_pts = 1:3;
b_sim = [BinStats.Mean_Sim]; sem_sim = [BinStats.SEM_Sim];
b_seq = [BinStats.Mean_Seq]; sem_seq = [BinStats.SEM_Seq];

errorbar(x_pts, b_sim, sem_sim, '--o', 'Color', color_sim, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 10, 'DisplayName', 'Simultaneous', 'MarkerSize', 8);
errorbar(x_pts, b_seq, sem_seq, '-s', 'Color', color_seq, 'LineWidth', 2, 'MarkerFaceColor', 'k', 'CapSize', 10, 'DisplayName', 'Sequential', 'MarkerSize', 8);

for b = 1:3
    p = BinStats(b).P;
    if isnan(p) || p >= 0.05, continue; end
    txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
    y_top = max(b_sim(b)+sem_sim(b), b_seq(b)+sem_seq(b));
    text(b, y_top + 1.0, txt, 'FontSize', 20, 'HorizontalAlignment', 'center', 'Color', 'k', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
end

xlabel('Amplitude Range', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Mean Duration (ms) \pm SEM', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 16, 'FontName', 'Times New Roman');
set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'out');
xlim([0.5, 3.5]); ylim([0, 25]); axis square; box off;
if save_figures, exportgraphics(gcf, fullfile(save_dir, 'Group_Binned_LinePlot.tiff'),'ContentType', 'vector'); end

%% ================= 5. PLOT 3: EFFECT SIZE (Binned Ribbon) =================
fprintf('Generating Binned Effect Size Plot...\n');
figure('Color','w', 'Position', [200 200 600 500], 'Name', 'Effect_Size_Plot'); hold on;
diff_mean = b_seq - b_sim;
diff_sem  = sqrt(sem_seq.^2 + sem_sim.^2);

x_conf = [x_pts, fliplr(x_pts)];
y_conf = [diff_mean + diff_sem, fliplr(diff_mean - diff_sem)];
fill(x_conf, y_conf, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
plot(x_pts, diff_mean, '-o', 'Color', 'k', 'LineWidth', 4, 'MarkerFaceColor', 'k', 'MarkerSize', 6);
yline(0, '--k', 'Alpha', 0.5); 

ylabel('\Delta Duration (Seq - Sim) [ms]', 'FontSize', 22, 'FontName', 'Times New Roman');
xlabel('Amplitude Range', 'FontSize', 22, 'FontName', 'Times New Roman');
set(gca, 'XTick', 1:3, 'XTickLabel', BinLabels, 'FontSize', 16, 'FontName', 'Times New Roman', 'LineWidth', 2, 'TickDir', 'out');
xlim([0.5, 3.5]); ylim([-2, 6]); axis square; box off;
if save_figures, exportgraphics(gcf, fullfile(save_dir, 'Group_Effect_Size.tiff'),'ContentType', 'vector'); end

%% ================= 6. STATISTICAL TESTS (Between Amplitude Ranges) =================
% [NEW] Running tests to prove Dose-Response relationship
AmpStats = struct();

% Simultaneous tests
p_kwallis_sim = kruskalwallis([Bins.Low.Sim; Bins.Mid.Sim; Bins.High.Sim], ...
    [ones(size(Bins.Low.Sim)); 2*ones(size(Bins.Mid.Sim)); 3*ones(size(Bins.High.Sim))], 'off');
p_sim_L_vs_M = ranksum(Bins.Low.Sim, Bins.Mid.Sim);
p_sim_M_vs_H = ranksum(Bins.Mid.Sim, Bins.High.Sim);

% Sequential tests
p_kwallis_seq = kruskalwallis([Bins.Low.Seq; Bins.Mid.Seq; Bins.High.Seq], ...
    [ones(size(Bins.Low.Seq)); 2*ones(size(Bins.Mid.Seq)); 3*ones(size(Bins.High.Seq))], 'off');
p_seq_L_vs_M = ranksum(Bins.Low.Seq, Bins.Mid.Seq);
p_seq_M_vs_H = ranksum(Bins.Mid.Seq, Bins.High.Seq);

%% ================= 7. COMMAND WINDOW PRINTOUTS =================
fprintf('\n====================================================================\n');
fprintf('                TABLE 1: SIM VS SEQ (MEANS +/- SEM)                 \n');
fprintf('====================================================================\n');
fprintf('%-15s | %-18s | %-18s | %-10s\n', 'Range', 'Sim', 'Seq', 'P-Value');
fprintf('--------------------------------------------------------------------\n');
for b = 1:3
    fprintf('%-15s | %5.2f +/- %4.2f ms  | %5.2f +/- %4.2f ms  | %.5f\n', ...
        BinLabels{b}, BinStats(b).Mean_Sim, BinStats(b).SEM_Sim, BinStats(b).Mean_Seq, BinStats(b).SEM_Seq, BinStats(b).P);
end

fprintf('\n====================================================================\n');
fprintf('                TABLE 2: SIM VS SEQ (MEDIANS & IQR)                 \n');
fprintf('====================================================================\n');
fprintf('%-15s | %-18s | %-18s | %-10s\n', 'Range', 'Sim', 'Seq', 'P-Value');
fprintf('--------------------------------------------------------------------\n');
for b = 1:3
    fprintf('%-15s | %5.2f (IQR %4.2f) | %5.2f (IQR %4.2f) | %.5f\n', ...
        BinLabels{b}, BinStats(b).Med_Sim, BinStats(b).IQR_Sim, BinStats(b).Med_Seq, BinStats(b).IQR_Seq, BinStats(b).P);
end

fprintf('\n====================================================================\n');
fprintf('          TABLE 3: DOSE-RESPONSE (BETWEEN AMPLITUDE RANGES)         \n');
fprintf('====================================================================\n');
fprintf('SIMULTANEOUS (Kruskal-Wallis Overall P = %.5f)\n', p_kwallis_sim);
fprintf('  * Low (1-3) vs Mid (4-6): P = %.5f\n', p_sim_L_vs_M);
fprintf('  * Mid (4-6) vs High(8-10): P = %.5f\n', p_sim_M_vs_H);
fprintf('--------------------------------------------------------------------\n');
fprintf('SEQUENTIAL (Kruskal-Wallis Overall P = %.5f)\n', p_kwallis_seq);
fprintf('  * Low (1-3) vs Mid (4-6): P = %.5f\n', p_seq_L_vs_M);
fprintf('  * Mid (4-6) vs High(8-10): P = %.5f\n', p_seq_M_vs_H);
fprintf('====================================================================\n');