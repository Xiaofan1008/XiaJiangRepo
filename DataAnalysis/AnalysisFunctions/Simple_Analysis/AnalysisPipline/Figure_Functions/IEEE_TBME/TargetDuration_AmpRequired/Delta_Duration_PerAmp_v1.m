%% ============================================================
%   GROUP DURATION ANALYSIS: EXACT AMPLITUDE DELTA DURATION
%   - Metric: 3*SD Total Duration
%   - Logic: 1.5x IQR filter -> Strict Pairing -> Median per Set
%   - NEW: Use exact amplitudes directly (No inverse interpolation)
%   - Output 1: Mean duration vs exact amplitude
%   - Output 2: Delta duration (Seq - Sim) vs exact amplitude
%   - Style: IEEE Publication Ready (Arial, 8.8cm Width, 9pt Font)
% ============================================================
clear; 
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
stats_min_n_threshold =2;

if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ================= 1. POOL DATA (Strictly Paired Medians at Exact Amplitudes) =================
fprintf('Pooling data from %d files (Calculating paired medians at exact amplitudes)...\n', length(file_paths));
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
            AmpData = Data.Set(s).Amp(a); 
            val = AmpData.Val;
            if val == 0, continue; end
            
            fName = sprintf('A_%.1f', val); 
            fName = strrep(fName, '.', 'p');
            
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_SetMedians = [];
                Pooled.(fName).Seq_SetMedians = [];
                All_Amps = [All_Amps, val]; %#ok<AGROW>
            end
            
            curr_sim = AmpData.Sim; 
            curr_seq = AmpData.Seq; 
            
            med_sim = NaN; 
            med_seq = NaN;
            
            if ~isempty(curr_sim)
                iqr_val = iqr(curr_sim); 
                q1 = prctile(curr_sim, 25); 
                q3 = prctile(curr_sim, 75);
                curr_sim(curr_sim < (q1 - IQR_out_parm*iqr_val) | curr_sim > (q3 + IQR_out_parm*iqr_val)) = [];
                med_sim = median(curr_sim, 'omitnan');
            end
            
            if ~isempty(curr_seq)
                iqr_val = iqr(curr_seq); 
                q1 = prctile(curr_seq, 25); 
                q3 = prctile(curr_seq, 75);
                curr_seq(curr_seq < (q1 - IQR_out_parm*iqr_val) | curr_seq > (q3 + IQR_out_parm*iqr_val)) = [];
                med_seq = median(curr_seq, 'omitnan');
            end
            
            % Strict Pairing: Only keep if BOTH exist
            if ~isnan(med_sim) && ~isnan(med_seq)
                Pooled.(fName).Sim_SetMedians = [Pooled.(fName).Sim_SetMedians; med_sim];
                Pooled.(fName).Seq_SetMedians = [Pooled.(fName).Seq_SetMedians; med_seq];
            end
        end
    end
end

All_Amps = sort(unique(All_Amps));

%% ================= 2. AGGREGATE EXACT AMPLITUDES & CALCULATE STATS =================
AmpStats = struct();

for k = 1:length(All_Amps)
    curr_amp = All_Amps(k);
    fName = sprintf('A_%.1f', curr_amp); 
    fName = strrep(fName, '.', 'p');
    
    d1 = Pooled.(fName).Sim_SetMedians;
    d2 = Pooled.(fName).Seq_SetMedians;
    
    AmpStats(k).Amp = curr_amp;
    AmpStats(k).Sim = d1;
    AmpStats(k).Seq = d2;
    
    AmpStats(k).Mean_Sim = mean(d1, 'omitnan');
    AmpStats(k).SEM_Sim  = std(d1, 'omitnan') / sqrt(length(d1));
    AmpStats(k).Mean_Seq = mean(d2, 'omitnan');
    AmpStats(k).SEM_Seq  = std(d2, 'omitnan') / sqrt(length(d2));
    
    AmpStats(k).Med_Sim  = median(d1, 'omitnan');
    AmpStats(k).IQR_Sim  = iqr(d1);
    AmpStats(k).Med_Seq  = median(d2, 'omitnan');
    AmpStats(k).IQR_Seq  = iqr(d2);
    
    if ~isempty(d1) && ~isempty(d2)
        deltas = d2 - d1;   % Seq - Sim
        AmpStats(k).Delta = deltas;
        AmpStats(k).Mean_Delta = mean(deltas, 'omitnan');
        AmpStats(k).SEM_Delta  = std(deltas, 'omitnan') / sqrt(length(deltas));
        AmpStats(k).P = signrank(d1, d2);
        AmpStats(k).N = length(deltas);
    else
        AmpStats(k).Delta = [];
        AmpStats(k).Mean_Delta = NaN;
        AmpStats(k).SEM_Delta  = NaN;
        AmpStats(k).P = NaN;
        AmpStats(k).N = 0;
    end
end

% Vectors for plotting
AmpVals   = [AmpStats.Amp];
Mean_Sim  = [AmpStats.Mean_Sim];
SEM_Sim   = [AmpStats.SEM_Sim];
Mean_Seq  = [AmpStats.Mean_Seq];
SEM_Seq   = [AmpStats.SEM_Seq];
Mean_Del  = [AmpStats.Mean_Delta];
SEM_Del   = [AmpStats.SEM_Delta];

% Add origin point for plotting only
Plot_AmpVals  = [0, AmpVals];
Plot_Sim_Mean = [0, Mean_Sim];
Plot_Sim_SEM  = [0, SEM_Sim];
Plot_Seq_Mean = [0, Mean_Seq];
Plot_Seq_SEM  = [0, SEM_Seq];
Plot_Del_Mean = [0, Mean_Del];
Plot_Del_SEM  = [0, SEM_Del];

%% ================= 3. PLOT 1: MEAN DURATION VS EXACT AMPLITUDE =================
fprintf('Generating exact amplitude duration plot...\n');
figure('Color','w', 'Units', 'centimeters', 'Position', [5, 5, 8.8, 8.8], 'Name', 'ExactAmplitude_Duration'); hold on;

errorbar(Plot_AmpVals, Plot_Sim_Mean, Plot_Sim_SEM, '.', 'Color', color_sim, ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p1 = plot(Plot_AmpVals, Plot_Sim_Mean, '--o', 'Color', color_sim, 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5, 'DisplayName', 'Simultaneous');

errorbar(Plot_AmpVals, Plot_Seq_Mean, Plot_Seq_SEM, '.', 'Color', color_seq, ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p2 = plot(Plot_AmpVals, Plot_Seq_Mean, '-s', 'Color', color_seq, 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sequential');

% Stats stars
fprintf('\n=== EXACT AMPLITUDE DURATION STATS (Paired Signed Rank) ===\n');
fprintf('Min N = %d\n', stats_min_n_threshold);

for k = 1:length(AmpVals)
    curr_amp = AmpVals(k);
    n_pairs = AmpStats(k).N;
    
    if n_pairs < stats_min_n_threshold
        fprintf('Amp %.1f uA: Skipped (Low N=%d)\n', curr_amp, n_pairs);
        continue;
    end
    
    p = AmpStats(k).P;
    txt = '';
    if p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    end
    
    if ~isempty(txt)
        y_top = max(Mean_Sim(k) + SEM_Sim(k), Mean_Seq(k) + SEM_Seq(k));
        text(curr_amp, y_top + 0.35, txt, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontName', 'Arial');
        fprintf('Amp %.1f uA: MARKED with %s (p=%.5f)\n', curr_amp, txt, p);
    else
        fprintf('Amp %.1f uA: Not Significant (p=%.5f)\n', curr_amp, p);
    end
end

ylabel('Mean Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');
xlabel('Amplitude (\muA)', 'FontSize', 9, 'FontName', 'Arial');
legend([p1, p2], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

set(gca, 'FontSize', 9, 'FontName', 'Arial', 'LineWidth', 1.0, 'TickDir', 'out');
xlim([0, max(AmpVals)]);
set(gca, 'XTick', 0:2:max(AmpVals));
ylim([0, ceil(max([Plot_Sim_Mean + Plot_Sim_SEM, Plot_Seq_Mean + Plot_Seq_SEM]) + 1.2)]);
box off; axis square;

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'ExactAmplitude_Duration.tiff'), 'ContentType', 'vector', 'Resolution', 300);
end

%% ================= 4. PLOT 2: DELTA DURATION VS EXACT AMPLITUDE =================
fprintf('Generating exact amplitude delta duration plot...\n');
figure('Color','w', 'Units', 'centimeters', 'Position', [10, 5, 8.8, 8.8], 'Name', 'ExactAmplitude_DeltaDuration'); hold on;

% Scatter background
jitter_w = 0.20;
for k = 1:length(AmpVals)
    curr_amp = AmpVals(k);
    delta_vals = AmpStats(k).Delta;
    
    for i = 1:length(delta_vals)
        scatter(curr_amp + (rand-0.5)*jitter_w, delta_vals(i), 12, [0.7 0.7 0.7], 'o', ...
            'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
    end
end

errorbar(Plot_AmpVals, Plot_Del_Mean, Plot_Del_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p3 = plot(Plot_AmpVals, Plot_Del_Mean, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Seq - Sim');

yline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

% Stats stars
fprintf('\n=== EXACT AMPLITUDE DELTA DURATION STATS (Paired Signed Rank) ===\n');

for k = 1:length(AmpVals)
    curr_amp = AmpVals(k);
    n_pairs = AmpStats(k).N;
    
    if n_pairs < stats_min_n_threshold
        fprintf('Amp %.1f uA: Skipped (Low N=%d)\n', curr_amp, n_pairs);
        continue;
    end
    
    p = AmpStats(k).P;
    txt = '';
    if p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    end
    
    if ~isempty(txt)
        y_top = Mean_Del(k) + SEM_Del(k);
        text(curr_amp, y_top + 0.20, txt, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontName', 'Arial');
    end
end

ylabel('\Delta Duration (Sequential - Simultaneous) (ms)', 'FontSize', 9, 'FontName', 'Arial');
xlabel('Amplitude (\muA)', 'FontSize', 9, 'FontName', 'Arial');
legend([p3], {'Seq - Sim'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

set(gca, 'FontSize', 9, 'FontName', 'Arial', 'LineWidth', 1.0, 'TickDir', 'out');
xlim([0, max(AmpVals)]);
set(gca, 'XTick', 0:2:max(AmpVals));
ylim([floor(min([Plot_Del_Mean - Plot_Del_SEM, 0]) - 0.5), ceil(max([Plot_Del_Mean + Plot_Del_SEM]) + 0.8)]);
box off; axis square;

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'ExactAmplitude_DeltaDuration.tiff'), 'ContentType', 'vector', 'Resolution', 300);
end

%% ================= 5. COMMAND WINDOW PRINTOUTS =================
fprintf('\n====================================================================\n');
fprintf('        TABLE 1: EXACT AMPLITUDE DURATION (MEAN +/- SEM)            \n');
fprintf('====================================================================\n');
fprintf('%-12s | %-16s | %-16s | %-16s | %-10s\n', 'Amp', 'Sim', 'Seq', 'Delta (Seq-Sim)', 'P-Value');
fprintf('--------------------------------------------------------------------\n');

for k = 1:length(AmpVals)
    fprintf('%-12.1f | %5.2f +/- %4.2f  | %5.2f +/- %4.2f  | %5.2f +/- %4.2f  | %s\n', ...
        AmpVals(k), Mean_Sim(k), SEM_Sim(k), Mean_Seq(k), SEM_Seq(k), Mean_Del(k), SEM_Del(k), format_p(AmpStats(k).P));
end

fprintf('\n====================================================================\n');
fprintf('           TABLE 2: EXACT AMPLITUDE DURATION (MEDIAN & IQR)         \n');
fprintf('====================================================================\n');
fprintf('%-12s | %-18s | %-18s | %-10s\n', 'Amp', 'Sim', 'Seq', 'P-Value');
fprintf('--------------------------------------------------------------------\n');

for k = 1:length(AmpVals)
    fprintf('%-12.1f | %5.2f (IQR %4.2f) | %5.2f (IQR %4.2f) | %s\n', ...
        AmpVals(k), AmpStats(k).Med_Sim, AmpStats(k).IQR_Sim, AmpStats(k).Med_Seq, AmpStats(k).IQR_Seq, format_p(AmpStats(k).P));
end

%% ================= HELPER FUNCTIONS =================
% Formats P-values cleanly to prevent 0.00000 printouts
function p_str = format_p(p)
    if isnan(p)
        p_str = 'NaN';
    else
        p_str = sprintf('%.20f', p);
    end
end