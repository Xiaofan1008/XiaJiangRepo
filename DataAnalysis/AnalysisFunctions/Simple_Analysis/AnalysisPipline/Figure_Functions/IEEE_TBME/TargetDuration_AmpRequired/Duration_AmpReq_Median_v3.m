%% ============================================================
%   GROUP DURATION ANALYSIS: EXACT AMPLITUDE + SELECTED MATCHED DURATIONS
%   - Metric: 3*SD Total Duration
%   - Logic: 1.5x IQR filter -> Strict Pairing -> Median per Set
%   - Output 1: Mean duration vs exact amplitude
%   - Output 2: Delta duration (Seq - Sim) vs exact amplitude
%   - NEW Output 3: Required amplitude at selected matched durations
%   - NEW Output 4: Delta required amplitude at selected matched durations
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
stats_min_n_threshold = 2;

% --- Small swapped-axis test: selected matched durations ---
target_durations = [4:0.5:10];
enforce_monotonic = true;

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

%% ================= 5. SMALL SWAPPED-AXIS TEST: SELECTED MATCHED DURATIONS =================
fprintf('Generating small swapped-axis test at selected matched durations...\n');

Pool_Sim_Target = [];     % [target_duration, required_amp]
Pool_Seq_Target = [];     % [target_duration, required_amp]
Pool_Paired_Target = [];  % [target_duration, required_amp_sim, required_amp_seq]

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), continue; end
    D = load(file_paths{f});
    if ~isfield(D, 'DurationData'), continue; end
    Data = D.DurationData;
    
    FileCurve = struct();
    FileAmps = [];
    
    for s = 1:length(Data.Set)
        if isempty(Data.Set(s).Label) || ~isfield(Data.Set(s), 'Amp'), continue; end
        
        for a = 1:length(Data.Set(s).Amp)
            AmpData = Data.Set(s).Amp(a);
            val = AmpData.Val;
            if val == 0, continue; end
            
            fName = sprintf('A_%.1f', val);
            fName = strrep(fName, '.', 'p');
            
            if ~isfield(FileCurve, fName)
                FileCurve.(fName).Val = val;
                FileCurve.(fName).Sim_SetMedians = [];
                FileCurve.(fName).Seq_SetMedians = [];
                FileAmps = [FileAmps, val]; %#ok<AGROW>
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
            
            if ~isnan(med_sim) && ~isnan(med_seq)
                FileCurve.(fName).Sim_SetMedians = [FileCurve.(fName).Sim_SetMedians; med_sim];
                FileCurve.(fName).Seq_SetMedians = [FileCurve.(fName).Seq_SetMedians; med_seq];
            end
        end
    end
    
    FileAmps = sort(unique(FileAmps));
    amps_file = [];
    dur_sim_file = [];
    dur_seq_file = [];
    
    for i = 1:length(FileAmps)
        curr_amp = FileAmps(i);
        fName = sprintf('A_%.1f', curr_amp);
        fName = strrep(fName, '.', 'p');
        
        if ~isfield(FileCurve, fName), continue; end
        
        s_sim = FileCurve.(fName).Sim_SetMedians;
        s_seq = FileCurve.(fName).Seq_SetMedians;
        
        if isempty(s_sim) || isempty(s_seq), continue; end
        
        amps_file(end+1,1) = curr_amp; %#ok<SAGROW>
        dur_sim_file(end+1,1) = median(s_sim, 'omitnan'); %#ok<SAGROW>
        dur_seq_file(end+1,1) = median(s_seq, 'omitnan'); %#ok<SAGROW>
    end
    
    if length(amps_file) < 2
        continue;
    end
    
    [amps_file, idx_sort] = sort(amps_file);
    dur_sim_file = dur_sim_file(idx_sort);
    dur_seq_file = dur_seq_file(idx_sort);
    
    if enforce_monotonic
        dur_sim_file = cummax(dur_sim_file);
        dur_seq_file = cummax(dur_seq_file);
    end
    
    [dur_sim_u, ia_sim] = unique(dur_sim_file, 'stable');
    amps_sim_u = amps_file(ia_sim);
    
    [dur_seq_u, ia_seq] = unique(dur_seq_file, 'stable');
    amps_seq_u = amps_file(ia_seq);
    
    if length(dur_sim_u) < 2 || length(dur_seq_u) < 2
        continue;
    end
    
    for t = 1:length(target_durations)
        target_val = target_durations(t);
        
        req_amp_sim = NaN;
        req_amp_seq = NaN;
        
        if target_val >= min(dur_sim_u) && target_val <= max(dur_sim_u)
            req_amp_sim = interp1(dur_sim_u, amps_sim_u, target_val, 'linear');
        end
        
        if target_val >= min(dur_seq_u) && target_val <= max(dur_seq_u)
            req_amp_seq = interp1(dur_seq_u, amps_seq_u, target_val, 'linear');
        end
        
        if ~isnan(req_amp_sim)
            Pool_Sim_Target = [Pool_Sim_Target; target_val, req_amp_sim];
        end
        if ~isnan(req_amp_seq)
            Pool_Seq_Target = [Pool_Seq_Target; target_val, req_amp_seq];
        end
        if ~isnan(req_amp_sim) && ~isnan(req_amp_seq)
            Pool_Paired_Target = [Pool_Paired_Target; target_val, req_amp_sim, req_amp_seq];
        end
    end
end

Unique_Targets = unique([Pool_Sim_Target(:,1); Pool_Seq_Target(:,1)]);
Unique_Targets = sort(Unique_Targets);

Grand_Target_Sim_Mean = []; Grand_Target_Sim_SEM = [];
Grand_Target_Seq_Mean = []; Grand_Target_Seq_SEM = [];
Grand_Target_Del_Mean = []; Grand_Target_Del_SEM = [];
Grand_Target_N = [];

for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    
    vs = Pool_Sim_Target(abs(Pool_Sim_Target(:,1)-target_val)<0.001, 2);
    Grand_Target_Sim_Mean(k) = mean(vs, 'omitnan');
    Grand_Target_Sim_SEM(k)  = std(vs, 'omitnan') / sqrt(length(vs));
    
    vq = Pool_Seq_Target(abs(Pool_Seq_Target(:,1)-target_val)<0.001, 2);
    Grand_Target_Seq_Mean(k) = mean(vq, 'omitnan');
    Grand_Target_Seq_SEM(k)  = std(vq, 'omitnan') / sqrt(length(vq));
    
    current_pairs = Pool_Paired_Target(abs(Pool_Paired_Target(:,1)-target_val)<0.001, :);
    delta_vals = current_pairs(:,2) - current_pairs(:,3);   % Sim - Seq
    
    Grand_Target_Del_Mean(k) = mean(delta_vals, 'omitnan');
    Grand_Target_Del_SEM(k)  = std(delta_vals, 'omitnan') / sqrt(length(delta_vals));
    Grand_Target_N(k)        = length(delta_vals);
end

Plot_Targets2    = [0; Unique_Targets];
Plot_Target_Sim  = [0, Grand_Target_Sim_Mean];
Plot_Target_SimE = [0, Grand_Target_Sim_SEM];
Plot_Target_Seq  = [0, Grand_Target_Seq_Mean];
Plot_Target_SeqE = [0, Grand_Target_Seq_SEM];
Plot_Target_Del  = [0, Grand_Target_Del_Mean];
Plot_Target_DelE = [0, Grand_Target_Del_SEM];

%% ================= 6. PLOT 3: REQUIRED AMPLITUDE AT SELECTED MATCHED DURATIONS =================
figure('Color','w', 'Units', 'centimeters', 'Position', [15, 5, 8.8, 8.8], 'Name', 'SelectedMatchedDuration_RequiredAmplitude'); hold on;

jitter_w = 0.20;
for i = 1:size(Pool_Seq_Target, 1)
    scatter(Pool_Seq_Target(i,1)+(rand-0.5)*jitter_w, Pool_Seq_Target(i,2), 12, [0.7 0.7 0.7], 's', ...
        'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
for i = 1:size(Pool_Sim_Target, 1)
    scatter(Pool_Sim_Target(i,1)+(rand-0.5)*jitter_w, Pool_Sim_Target(i,2), 12, [0.8 0.8 0.8], 'o', ...
        'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end

errorbar(Plot_Targets2, Plot_Target_Sim, Plot_Target_SimE, '.', 'Color', color_sim, ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p4 = plot(Plot_Targets2, Plot_Target_Sim, '--o', 'Color', color_sim, 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'w', 'MarkerSize', 5, 'DisplayName', 'Simultaneous');

errorbar(Plot_Targets2, Plot_Target_Seq, Plot_Target_SeqE, '.', 'Color', color_seq, ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p5 = plot(Plot_Targets2, Plot_Target_Seq, '-s', 'Color', color_seq, 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sequential');

ylabel('Required Amplitude (\muA)', 'FontSize', 9, 'FontName', 'Arial');
xlabel('Matched Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');
legend([p4, p5], {'Simultaneous', 'Sequential'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

set(gca, 'FontSize', 9, 'FontName', 'Arial', 'LineWidth', 1.0, 'TickDir', 'out');
xlim([0, max(target_durations)]);
set(gca, 'XTick', 0:2:max(target_durations));
ylim([0, ceil(max([Plot_Target_Sim + Plot_Target_SimE, Plot_Target_Seq + Plot_Target_SeqE]) + 1.0)]);
box off; axis square;

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'SelectedMatchedDuration_RequiredAmplitude.tiff'), 'ContentType', 'vector', 'Resolution', 300);
end

%% ================= 7. PLOT 4: DELTA REQUIRED AMPLITUDE AT SELECTED MATCHED DURATIONS =================
figure('Color','w', 'Units', 'centimeters', 'Position', [20, 5, 8.8, 8.8], 'Name', 'SelectedMatchedDuration_DeltaAmplitude'); hold on;

jitter_w = 0.20;
for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired_Target(abs(Pool_Paired_Target(:,1)-target_val)<0.001, :);
    delta_vals = current_pairs(:,2) - current_pairs(:,3);   % Sim - Seq
    
    for i = 1:length(delta_vals)
        scatter(target_val+(rand-0.5)*jitter_w, delta_vals(i), 12, [0.7 0.7 0.7], 'o', ...
            'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
    end
end

errorbar(Plot_Targets2, Plot_Target_Del, Plot_Target_DelE, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p6 = plot(Plot_Targets2, Plot_Target_Del, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Sim - Seq');

yline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

ylabel('\Delta Required Amplitude (Sim - Seq, \muA)', 'FontSize', 9, 'FontName', 'Arial');
xlabel('Matched Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');
legend([p6], {'Sim - Seq'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

set(gca, 'FontSize', 9, 'FontName', 'Arial', 'LineWidth', 1.0, 'TickDir', 'out');
xlim([0, max(target_durations)]);
set(gca, 'XTick', 0:2:max(target_durations));
ylim([floor(min([Plot_Target_Del - Plot_Target_DelE, 0]) - 0.5), ceil(max([Plot_Target_Del + Plot_Target_DelE]) + 0.8)]);
box off; axis square;

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'SelectedMatchedDuration_DeltaAmplitude.tiff'), 'ContentType', 'vector', 'Resolution', 300);
end

%% ================= 8. COMMAND WINDOW PRINTOUTS =================
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
fprintf('        TABLE 2: SELECTED MATCHED DURATION REQUIRED AMPLITUDE       \n');
fprintf('====================================================================\n');
fprintf('%-16s | %-16s | %-16s | %-16s\n', 'Matched Duration', 'Sim', 'Seq', 'Delta (Sim-Seq)');
fprintf('--------------------------------------------------------------------\n');

for k = 1:length(Unique_Targets)
    fprintf('%-16.1f | %5.2f +/- %4.2f  | %5.2f +/- %4.2f  | %5.2f +/- %4.2f\n', ...
        Unique_Targets(k), ...
        Grand_Target_Sim_Mean(k), Grand_Target_Sim_SEM(k), ...
        Grand_Target_Seq_Mean(k), Grand_Target_Seq_SEM(k), ...
        Grand_Target_Del_Mean(k), Grand_Target_Del_SEM(k));
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