%% ============================================================
%   GROUP DURATION ANALYSIS: REQUIRED AMPLITUDE FOR MATCHED DURATION
%   - Metric: 3*SD Total Duration
%   - Logic: 1.5x IQR filter -> Strict Pairing -> Median per Set
%   - NEW: Build per-file exact-amplitude duration curves first
%   - NEW: Invert each file's duration curve (target duration -> required amp)
%   - Output 1: Required amplitude vs matched duration
%   - Output 2: Delta required amplitude (Sim - Seq) vs matched duration
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

% --- INVERSE ANALYSIS SETTINGS ---
target_durations = 5 : 0.5 : 10;   % Stage 1 target matched duration range (ms)
stats_min_n_threshold = 2;         % Minimum paired dataset number for stats
enforce_monotonic = true;          % Optional monotonic enforcement before inversion

if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% ================= 1. BUILD PER-FILE CURVES & INVERT =================
fprintf('Processing %d files (Building per-file exact-amplitude duration curves)...\n', length(file_paths));

Pool_Sim = [];         % [target_duration, required_amp]
Pool_Seq = [];         % [target_duration, required_amp]
Pool_Paired = [];      % [target_duration, required_amp_sim, required_amp_seq]

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), continue; end
    D = load(file_paths{f});
    if ~isfield(D, 'DurationData'), continue; end
    Data = D.DurationData;

    % --- Per-file temporary storage at exact amplitude level ---
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

            % Strict pairing: only keep if BOTH exist
            if ~isnan(med_sim) && ~isnan(med_seq)
                FileCurve.(fName).Sim_SetMedians = [FileCurve.(fName).Sim_SetMedians; med_sim];
                FileCurve.(fName).Seq_SetMedians = [FileCurve.(fName).Seq_SetMedians; med_seq];
            end
        end
    end

    % --- Collapse each file to one Sim/Seq median duration per exact amplitude ---
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

    % Need at least 2 amplitudes for interpolation
    if length(amps_file) < 2
        continue;
    end

    % Sort by amplitude
    [amps_file, idx_sort] = sort(amps_file);
    dur_sim_file = dur_sim_file(idx_sort);
    dur_seq_file = dur_seq_file(idx_sort);

    % Optional monotonic enforcement for inverse interpolation
    if enforce_monotonic
        dur_sim_file = cummax(dur_sim_file);
        dur_seq_file = cummax(dur_seq_file);
    end

    % Remove duplicate durations before interp1
    [dur_sim_u, ia_sim] = unique(dur_sim_file, 'stable');
    amps_sim_u = amps_file(ia_sim);

    [dur_seq_u, ia_seq] = unique(dur_seq_file, 'stable');
    amps_seq_u = amps_file(ia_seq);

    if length(dur_sim_u) < 2 || length(dur_seq_u) < 2
        continue;
    end

    % Invert each file: target duration -> required amplitude
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

        % Unpaired pools
        if ~isnan(req_amp_sim)
            Pool_Sim = [Pool_Sim; target_val, req_amp_sim];
        end
        if ~isnan(req_amp_seq)
            Pool_Seq = [Pool_Seq; target_val, req_amp_seq];
        end

        % Paired pool
        if ~isnan(req_amp_sim) && ~isnan(req_amp_seq)
            Pool_Paired = [Pool_Paired; target_val, req_amp_sim, req_amp_seq];
        end
    end
end

%% ================= 2. AGGREGATE ACROSS FILES =================
Unique_Targets = unique([Pool_Sim(:,1); Pool_Seq(:,1)]);
Unique_Targets = sort(Unique_Targets);

Grand_Sim_Mean = []; Grand_Sim_SEM = []; Grand_Sim_N = [];
Grand_Seq_Mean = []; Grand_Seq_SEM = []; Grand_Seq_N = [];

Grand_Delta_Mean = []; Grand_Delta_SEM = []; Grand_Delta_N = [];

for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);

    vs = Pool_Sim(abs(Pool_Sim(:,1)-target_val)<0.001, 2);
    Grand_Sim_Mean(k) = mean(vs, 'omitnan');
    Grand_Sim_SEM(k)  = std(vs, 'omitnan') / sqrt(length(vs));
    Grand_Sim_N(k)    = length(vs);

    vq = Pool_Seq(abs(Pool_Seq(:,1)-target_val)<0.001, 2);
    Grand_Seq_Mean(k) = mean(vq, 'omitnan');
    Grand_Seq_SEM(k)  = std(vq, 'omitnan') / sqrt(length(vq));
    Grand_Seq_N(k)    = length(vq);

    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);
    delta_vals = current_pairs(:,2) - current_pairs(:,3);   % Sim - Seq

    Grand_Delta_Mean(k) = mean(delta_vals, 'omitnan');
    Grand_Delta_SEM(k)  = std(delta_vals, 'omitnan') / sqrt(length(delta_vals));
    Grand_Delta_N(k)    = length(delta_vals);
end

% --- Add origin point for plotting only ---
Plot_Targets    = [0; Unique_Targets];

Plot_Sim_Mean   = [0, Grand_Sim_Mean];
Plot_Sim_SEM    = [0, Grand_Sim_SEM];

Plot_Seq_Mean   = [0, Grand_Seq_Mean];
Plot_Seq_SEM    = [0, Grand_Seq_SEM];

Plot_Delta_Mean = [0, Grand_Delta_Mean];
Plot_Delta_SEM  = [0, Grand_Delta_SEM];

%% ================= 3. PLOT 1: INVERSE CURVES =================
fprintf('Generating inverse duration curve plot...\n');
figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 8.8, 8.8], 'Name', 'Duration_InverseCurve'); hold on;

% A. Scatter (Background)
jitter_w = 0.08;
for i = 1:size(Pool_Seq, 1)
    scatter(Pool_Seq(i,1)+(rand-0.5)*jitter_w, Pool_Seq(i,2), 12, [0.7 0.7 0.7], 's', ...
        'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end
for i = 1:size(Pool_Sim, 1)
    scatter(Pool_Sim(i,1)+(rand-0.5)*jitter_w, Pool_Sim(i,2), 12, [0.8 0.8 0.8], 'o', ...
        'LineWidth', 0.7, 'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end

% B. Main Curves
errorbar(Plot_Targets, Plot_Sim_Mean, Plot_Sim_SEM, '.', 'Color', color_sim, ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p1 = plot(Plot_Targets, Plot_Sim_Mean, '--o', 'Color', color_sim, 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');

errorbar(Plot_Targets, Plot_Seq_Mean, Plot_Seq_SEM, '.', 'Color', color_seq, ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p2 = plot(Plot_Targets, Plot_Seq_Mean, '-s', 'Color', color_seq, 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

% C. Stats
fprintf('\n=== INVERSE DURATION STATS (Paired Signed Rank) ===\n');
fprintf('Min N = %d\n', stats_min_n_threshold);

for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);

    data_s = current_pairs(:,2);
    data_q = current_pairs(:,3);
    n_pairs = length(data_s);

    if n_pairs < stats_min_n_threshold
        fprintf('Target %.2f ms: Skipped (Low N=%d)\n', target_val, n_pairs);
        continue;
    end

    p = signrank(data_s, data_q);

    txt = '';
    if p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    end

    if ~isempty(txt)
        y_top = max(Grand_Sim_Mean(k) + Grand_Sim_SEM(k), Grand_Seq_Mean(k) + Grand_Seq_SEM(k));
        text(target_val, y_top + 0.35, txt, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontName', 'Arial');
        fprintf('Target %.2f ms: MARKED with %s (p=%.5f)\n', target_val, txt, p);
    else
        fprintf('Target %.2f ms: Not Significant (p=%.5f)\n', target_val, p);
    end
end

% --- Formatting ---
box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.0);
axis square;

xlabel('Matched Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Required Amplitude (\muA)', 'FontSize', 9, 'FontName', 'Arial');

legend([p1, p2], 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

xlim([0, max(Unique_Targets)]);
set(gca, 'XTick', 0:2:max(Unique_Targets));

ylim([0, ceil(max([Plot_Sim_Mean + Plot_Sim_SEM, Plot_Seq_Mean + Plot_Seq_SEM]) + 1.4)]);

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'Duration_InverseCurve.tiff'), 'ContentType', 'vector', 'Resolution', 300);
end

%% ================= 4. PLOT 2: DELTA AMPLITUDE =================
fprintf('Generating inverse duration delta plot...\n');
figure('Color', 'w', 'Units', 'centimeters', 'Position', [10, 5, 8.8, 8.8], 'Name', 'Duration_DeltaAmplitude'); hold on;

% A. Scatter (Background)
jitter_w = 0.08;
for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);
    delta_vals = current_pairs(:,2) - current_pairs(:,3);

    for i = 1:length(delta_vals)
        scatter(target_val+(rand-0.5)*jitter_w, delta_vals(i), 12, [0.7 0.7 0.7], 'o', ...
            'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
    end
end

% B. Main Curve
errorbar(Plot_Targets, Plot_Delta_Mean, Plot_Delta_SEM, '.', 'Color', 'k', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p3 = plot(Plot_Targets, Plot_Delta_Mean, '-o', 'Color', 'k', 'LineWidth', 1.5, ...
    'MarkerSize', 5, 'MarkerFaceColor', 'k', 'DisplayName', 'Sim - Seq');

yline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');

% C. Stats
fprintf('\n=== DELTA AMPLITUDE STATS (Paired Signed Rank) ===\n');
for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);

    data_s = current_pairs(:,2);
    data_q = current_pairs(:,3);
    n_pairs = length(data_s);

    if n_pairs < stats_min_n_threshold
        fprintf('Target %.2f ms: Skipped (Low N=%d)\n', target_val, n_pairs);
        continue;
    end

    p = signrank(data_s, data_q);

    txt = '';
    if p < 0.001
        txt = '***';
    elseif p < 0.01
        txt = '**';
    elseif p < 0.05
        txt = '*';
    end

    if ~isempty(txt)
        y_top = Grand_Delta_Mean(k) + Grand_Delta_SEM(k);
        text(target_val, y_top + 0.20, txt, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontName', 'Arial');
    end
end

% --- Formatting ---
box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.0);
axis square;

xlabel('Matched Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('\Delta Required Amplitude (Sim - Seq, \muA)', 'FontSize', 9, 'FontName', 'Arial');

legend([p3], 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

xlim([0, max(Unique_Targets)]);
set(gca, 'XTick', 0:2:max(Unique_Targets));

ylim([floor(min([Plot_Delta_Mean - Plot_Delta_SEM, 0]) - 0.5), ceil(max([Plot_Delta_Mean + Plot_Delta_SEM]) + 0.8)]);

if save_figures
    exportgraphics(gcf, fullfile(save_dir, 'Duration_DeltaAmplitude.tiff'), 'ContentType', 'vector', 'Resolution', 300);
end

%% ================= 5. COMMAND WINDOW PRINTOUTS =================
fprintf('\n====================================================================\n');
fprintf('     TABLE 1: REQUIRED AMPLITUDE FOR MATCHED DURATION (MEAN +/- SEM)\n');
fprintf('====================================================================\n');
fprintf('%-16s | %-16s | %-16s | %-16s\n', 'Matched Duration', 'Sim', 'Seq', 'Delta (Sim-Seq)');
fprintf('--------------------------------------------------------------------\n');

for k = 1:length(Unique_Targets)
    fprintf('%-16.2f | %5.2f +/- %4.2f  | %5.2f +/- %4.2f  | %5.2f +/- %4.2f\n', ...
        Unique_Targets(k), ...
        Grand_Sim_Mean(k), Grand_Sim_SEM(k), ...
        Grand_Seq_Mean(k), Grand_Seq_SEM(k), ...
        Grand_Delta_Mean(k), Grand_Delta_SEM(k));
end

fprintf('\n====================================================================\n');
fprintf('     TABLE 2: PAIRED SIGNRANK TEST FOR MATCHED DURATION LEVELS       \n');
fprintf('====================================================================\n');

for k = 1:length(Unique_Targets)
    target_val = Unique_Targets(k);
    current_pairs = Pool_Paired(abs(Pool_Paired(:,1)-target_val)<0.001, :);
    n_pairs = size(current_pairs,1);

    if n_pairs < stats_min_n_threshold
        fprintf('Target %.2f ms: Skipped (Low N=%d)\n', target_val, n_pairs);
        continue;
    end

    p = signrank(current_pairs(:,2), current_pairs(:,3));
    fprintf('Target %.2f ms: N=%d, P=%s\n', target_val, n_pairs, format_p(p));
end

%% ================= HELPER FUNCTIONS =================
function p_str = format_p(p)
    if isnan(p)
        p_str = 'NaN';
    else
        p_str = sprintf('%.20f', p);
    end
end