%% ============================================================
%   GRAND AVERAGE: IMPLANTATION COMPARISON (Vertical vs Parallel)
%   - Logic: Aggregates pre-calculated normalized recruitment curves.
%   - Grouping: Bins data by exact amplitude (No interpolation).
%   - Statistics: Mean +/- SEM across DATASETS.
%   - Comparison:
%         1) Simultaneous: Vertical vs Parallel
%         2) Sequential:   Vertical vs Parallel
%   - Statistics:
%         rank-sum test at each amplitude (unpaired)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% --- Statistical settings ---
stats_phys_threshold   = 1;   % ignore very low amplitudes if needed
stats_min_n_threshold  = 3;   % minimum dataset number per group
show_sig_star          = true;

% ---------------- VERTICAL implantation ----------------
file_paths_vertical = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim8.mat';    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX012/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX010/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX009/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX009/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX005/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';
};

% ---------------- PARALLEL implantation ----------------
file_paths_parallel = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';

    };

% --- Plot Settings ---
save_figure = false;
save_dir    = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/AnalysisPipline/Figure_Functions/IEEE_TBME/Implantation_Comparison';
fig_name    = 'Implantation_Comparison_Vertical_vs_Parallel.tiff';

%% ================= 2. AGGREGATE DATA =================
fprintf('Processing implantation comparison...\n');
fprintf('  Vertical datasets : %d\n', length(file_paths_vertical));
fprintf('  Parallel datasets : %d\n', length(file_paths_parallel));

% Pools for plotting
Pool_Sim_Vert = [];
Pool_Sim_Para = [];
Pool_Seq_Vert = [];
Pool_Seq_Para = [];

% Keep track of dataset counts
n_vert_valid = 0;
n_para_valid = 0;

% ---------------- VERTICAL ----------------
for i = 1:length(file_paths_vertical)
    if ~exist(file_paths_vertical{i}, 'file')
        fprintf('Missing vertical file: %s\n', file_paths_vertical{i});
        continue;
    end

    D = load(file_paths_vertical{i});
    if ~isfield(D, 'ResultNorm')
        fprintf('Skipped vertical file (no ResultNorm): %s\n', file_paths_vertical{i});
        continue;
    end

    R = D.ResultNorm;
    Amps = R.Amps;

    % dataset-level mean curve
    sim_ds_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));

    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_ds_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_ds_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end

    n_vert_valid = n_vert_valid + 1;

    for a = 1:length(Amps)
        if abs(Amps(a)) < 0.001
            v_sim = 0;
            v_seq = 0;
        else
            v_sim = sim_ds_mean(a);
            v_seq = seq_ds_mean(a);
        end

        if ~isnan(v_sim), Pool_Sim_Vert = [Pool_Sim_Vert; Amps(a), v_sim]; end
        if ~isnan(v_seq), Pool_Seq_Vert = [Pool_Seq_Vert; Amps(a), v_seq]; end
    end
end

% ---------------- PARALLEL ----------------
for i = 1:length(file_paths_parallel)
    if ~exist(file_paths_parallel{i}, 'file')
        fprintf('Missing parallel file: %s\n', file_paths_parallel{i});
        continue;
    end

    D = load(file_paths_parallel{i});
    if ~isfield(D, 'ResultNorm')
        fprintf('Skipped parallel file (no ResultNorm): %s\n', file_paths_parallel{i});
        continue;
    end

    R = D.ResultNorm;
    Amps = R.Amps;

    % dataset-level mean curve
    sim_ds_mean = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));

    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4
        seq_ds_mean = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else
        seq_ds_mean = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan'));
    end

    n_para_valid = n_para_valid + 1;

    for a = 1:length(Amps)
        if abs(Amps(a)) < 0.001
            v_sim = 0;
            v_seq = 0;
        else
            v_sim = sim_ds_mean(a);
            v_seq = seq_ds_mean(a);
        end

        if ~isnan(v_sim), Pool_Sim_Para = [Pool_Sim_Para; Amps(a), v_sim]; end
        if ~isnan(v_seq), Pool_Seq_Para = [Pool_Seq_Para; Amps(a), v_seq]; end
    end
end

fprintf('\nValid datasets loaded:\n');
fprintf('  Vertical : %d\n', n_vert_valid);
fprintf('  Parallel : %d\n', n_para_valid);

%% ================= 3. BUILD GRAND MEAN VECTORS =================
all_amps = [
    Pool_Sim_Vert(:,1); Pool_Sim_Para(:,1); ...
    Pool_Seq_Vert(:,1); Pool_Seq_Para(:,1)
];
Unique_Amps = unique(all_amps);
Unique_Amps = sort(Unique_Amps);
if Unique_Amps(1) ~= 0
    Unique_Amps = [0; Unique_Amps];
end

Grand_Sim_Vert_Mean = nan(size(Unique_Amps));
Grand_Sim_Vert_SEM  = nan(size(Unique_Amps));
Grand_Sim_Para_Mean = nan(size(Unique_Amps));
Grand_Sim_Para_SEM  = nan(size(Unique_Amps));

Grand_Seq_Vert_Mean = nan(size(Unique_Amps));
Grand_Seq_Vert_SEM  = nan(size(Unique_Amps));
Grand_Seq_Para_Mean = nan(size(Unique_Amps));
Grand_Seq_Para_SEM  = nan(size(Unique_Amps));

for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);

    % ---- Sim Vertical ----
    vals = Pool_Sim_Vert(abs(Pool_Sim_Vert(:,1)-amp)<0.001, 2);
    if isempty(vals) && amp == 0, vals = 0; end
    if ~isempty(vals)
        Grand_Sim_Vert_Mean(k) = mean(vals, 'omitnan');
        Grand_Sim_Vert_SEM(k)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    end

    % ---- Sim Parallel ----
    vals = Pool_Sim_Para(abs(Pool_Sim_Para(:,1)-amp)<0.001, 2);
    if isempty(vals) && amp == 0, vals = 0; end
    if ~isempty(vals)
        Grand_Sim_Para_Mean(k) = mean(vals, 'omitnan');
        Grand_Sim_Para_SEM(k)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    end

    % ---- Seq Vertical ----
    vals = Pool_Seq_Vert(abs(Pool_Seq_Vert(:,1)-amp)<0.001, 2);
    if isempty(vals) && amp == 0, vals = 0; end
    if ~isempty(vals)
        Grand_Seq_Vert_Mean(k) = mean(vals, 'omitnan');
        Grand_Seq_Vert_SEM(k)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    end

    % ---- Seq Parallel ----
    vals = Pool_Seq_Para(abs(Pool_Seq_Para(:,1)-amp)<0.001, 2);
    if isempty(vals) && amp == 0, vals = 0; end
    if ~isempty(vals)
        Grand_Seq_Para_Mean(k) = mean(vals, 'omitnan');
        Grand_Seq_Para_SEM(k)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    end
end

%% ================= 4. COMMAND WINDOW STATS =================
fprintf('\n============================================================\n');
fprintf('IMPLANTATION COMPARISON: VERTICAL vs PARALLEL\n');
fprintf('============================================================\n');
fprintf('Statistical test: rank-sum (unpaired)\n');
fprintf('Minimum N per group: %d\n', stats_min_n_threshold);
fprintf('Phys threshold: amp > %.1f uA\n', stats_phys_threshold);
fprintf('============================================================\n');

fprintf('\n---------------- SIMULTANEOUS ----------------\n');
for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    if amp == 0, continue; end

    data_vert = Pool_Sim_Vert(abs(Pool_Sim_Vert(:,1)-amp)<0.001, 2);
    data_para = Pool_Sim_Para(abs(Pool_Sim_Para(:,1)-amp)<0.001, 2);

    n_vert = length(data_vert);
    n_para = length(data_para);

    if amp <= stats_phys_threshold
        fprintf('Amp %.1f uA: Skipped (Below Phys Threshold)\n', amp);
        continue;
    end

    if n_vert < stats_min_n_threshold || n_para < stats_min_n_threshold
        fprintf('Amp %.1f uA: Skipped (Low N, Vert=%d, Para=%d)\n', amp, n_vert, n_para);
        continue;
    end

    p = ranksum(data_vert, data_para);

    fprintf(['Amp %.1f uA | Vert mean = %.4f ± %.4f (N=%d) | ', ...
             'Para mean = %.4f ± %.4f (N=%d) | p = %.6f\n'], ...
             amp, ...
             mean(data_vert,'omitnan'), std(data_vert,'omitnan')/sqrt(n_vert), n_vert, ...
             mean(data_para,'omitnan'), std(data_para,'omitnan')/sqrt(n_para), n_para, ...
             p);
end

fprintf('\n---------------- SEQUENTIAL ----------------\n');
for k = 1:length(Unique_Amps)
    amp = Unique_Amps(k);
    if amp == 0, continue; end

    data_vert = Pool_Seq_Vert(abs(Pool_Seq_Vert(:,1)-amp)<0.001, 2);
    data_para = Pool_Seq_Para(abs(Pool_Seq_Para(:,1)-amp)<0.001, 2);

    n_vert = length(data_vert);
    n_para = length(data_para);

    if amp <= stats_phys_threshold
        fprintf('Amp %.1f uA: Skipped (Below Phys Threshold)\n', amp);
        continue;
    end

    if n_vert < stats_min_n_threshold || n_para < stats_min_n_threshold
        fprintf('Amp %.1f uA: Skipped (Low N, Vert=%d, Para=%d)\n', amp, n_vert, n_para);
        continue;
    end

    p = ranksum(data_vert, data_para);

    fprintf(['Amp %.1f uA | Vert mean = %.4f ± %.4f (N=%d) | ', ...
             'Para mean = %.4f ± %.4f (N=%d) | p = %.6f\n'], ...
             amp, ...
             mean(data_vert,'omitnan'), std(data_vert,'omitnan')/sqrt(n_vert), n_vert, ...
             mean(data_para,'omitnan'), std(data_para,'omitnan')/sqrt(n_para), n_para, ...
             p);
end

%% ================= 5. PLOT =================
figure('Units', 'centimeters', 'Position', [2, 2, 17.8, 8.89], ...
    'Color', 'w', 'PaperPositionMode', 'auto');

% ============================================================
% PANEL 1: SIMULTANEOUS
% ============================================================
subplot(1,2,1); hold on;

% background scatter
jitter_w = 0.4;
for i = 1:size(Pool_Sim_Vert,1)
    scatter(Pool_Sim_Vert(i,1)+(rand-0.5)*jitter_w, Pool_Sim_Vert(i,2), ...
        12, [0.8 0.8 0.8], 'o', 'LineWidth', 0.7, ...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end
for i = 1:size(Pool_Sim_Para,1)
    scatter(Pool_Sim_Para(i,1)+(rand-0.5)*jitter_w, Pool_Sim_Para(i,2), ...
        12, [0.6 0.6 0.6], 's', 'filled', ...
        'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end

% mean curves
errorbar(Unique_Amps, Grand_Sim_Vert_Mean, Grand_Sim_Vert_SEM, '.', ...
    'Color', 'k', 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p1 = plot(Unique_Amps, Grand_Sim_Vert_Mean, '--o', ...
    'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Vertical');

errorbar(Unique_Amps, Grand_Sim_Para_Mean, Grand_Sim_Para_SEM, '.', ...
    'Color', 'k', 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p2 = plot(Unique_Amps, Grand_Sim_Para_Mean, '-s', ...
    'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Parallel');

% significance stars
if show_sig_star
    for k = 1:length(Unique_Amps)
        amp = Unique_Amps(k);
        if amp == 0 || amp <= stats_phys_threshold, continue; end

        data_vert = Pool_Sim_Vert(abs(Pool_Sim_Vert(:,1)-amp)<0.001, 2);
        data_para = Pool_Sim_Para(abs(Pool_Sim_Para(:,1)-amp)<0.001, 2);

        n_vert = length(data_vert);
        n_para = length(data_para);

        if n_vert < stats_min_n_threshold || n_para < stats_min_n_threshold
            continue;
        end

        p = ranksum(data_vert, data_para);

        txt = '';
        if p < 0.001
            txt = '***';
        elseif p < 0.01
            txt = '**';
        elseif p < 0.05
            txt = '*';
        end

        if ~isempty(txt)
            y_top = max(Grand_Sim_Vert_Mean(k)+Grand_Sim_Vert_SEM(k), ...
                        Grand_Sim_Para_Mean(k)+Grand_Sim_Para_SEM(k));
            text(amp, y_top + 0.08, txt, 'FontSize', 10, ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
                'FontName', 'Arial');
        end
    end
end

box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1);
axis square;
xlabel('Amplitude (µA)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Normalized Spike Count (a.u.)', 'FontSize', 9, 'FontName', 'Arial');
title('Simultaneous', 'FontSize', 9, 'FontName', 'Arial', 'FontWeight', 'normal');
legend([p1, p2], 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
xlim([0, max(Unique_Amps)]);
ylim([0 2]);
set(gca, 'YTick', 0:0.4:2);

% ============================================================
% PANEL 2: SEQUENTIAL
% ============================================================
subplot(1,2,2); hold on;

% background scatter
for i = 1:size(Pool_Seq_Vert,1)
    scatter(Pool_Seq_Vert(i,1)+(rand-0.5)*jitter_w, Pool_Seq_Vert(i,2), ...
        12, [0.8 0.8 0.8], 'o', 'LineWidth', 0.7, ...
        'MarkerEdgeAlpha', 0.5, 'HandleVisibility', 'off');
end
for i = 1:size(Pool_Seq_Para,1)
    scatter(Pool_Seq_Para(i,1)+(rand-0.5)*jitter_w, Pool_Seq_Para(i,2), ...
        12, [0.6 0.6 0.6], 's', 'filled', ...
        'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
end

% mean curves
errorbar(Unique_Amps, Grand_Seq_Vert_Mean, Grand_Seq_Vert_SEM, '.', ...
    'Color', 'k', 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p3 = plot(Unique_Amps, Grand_Seq_Vert_Mean, '--o', ...
    'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Vertical');

errorbar(Unique_Amps, Grand_Seq_Para_Mean, Grand_Seq_Para_SEM, '.', ...
    'Color', 'k', 'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
p4 = plot(Unique_Amps, Grand_Seq_Para_Mean, '-s', ...
    'Color', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Parallel');

% significance stars
if show_sig_star
    for k = 1:length(Unique_Amps)
        amp = Unique_Amps(k);
        if amp == 0 || amp <= stats_phys_threshold, continue; end

        data_vert = Pool_Seq_Vert(abs(Pool_Seq_Vert(:,1)-amp)<0.001, 2);
        data_para = Pool_Seq_Para(abs(Pool_Seq_Para(:,1)-amp)<0.001, 2);

        n_vert = length(data_vert);
        n_para = length(data_para);

        if n_vert < stats_min_n_threshold || n_para < stats_min_n_threshold
            continue;
        end

        p = ranksum(data_vert, data_para);

        txt = '';
        if p < 0.001
            txt = '***';
        elseif p < 0.01
            txt = '**';
        elseif p < 0.05
            txt = '*';
        end

        if ~isempty(txt)
            y_top = max(Grand_Seq_Vert_Mean(k)+Grand_Seq_Vert_SEM(k), ...
                        Grand_Seq_Para_Mean(k)+Grand_Seq_Para_SEM(k));
            text(amp, y_top + 0.08, txt, 'FontSize', 10, ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
                'FontName', 'Arial');
        end
    end
end

box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1);
axis square;
xlabel('Amplitude (µA)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Normalized Spike Count (a.u.)', 'FontSize', 9, 'FontName', 'Arial');
title('Sequential', 'FontSize', 9, 'FontName', 'Arial', 'FontWeight', 'normal');
legend([p3, p4], 'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
xlim([0, max(Unique_Amps)]);
ylim([0 2]);
set(gca, 'YTick', 0:0.4:2);

%% ================= 6. SAVE FIGURE =================
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    save_path = fullfile(save_dir, fig_name);
    print(gcf, save_path, '-dtiff', '-r300');
    fprintf('\nFigure saved to: %s\n', save_path);
end