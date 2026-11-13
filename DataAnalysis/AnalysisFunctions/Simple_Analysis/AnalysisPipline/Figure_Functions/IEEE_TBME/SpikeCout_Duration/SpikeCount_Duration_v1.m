%% ============================================================
%   POPULATION ANALYSIS: SPIKE COUNT vs DURATION
%   - Input: Result_Combined_SpikeDuration_*.mat files
%   - Logic:
%       1. Load all single-dataset combined result files.
%       2. Pool channel-level data across datasets.
%       3. Ignore amp = 0 points.
%       4. Use only valid points from the saved responsive pool.
%       5. Plot Spike Count vs Duration for Sim and Seq.
%       6. Run feasibility statistics:
%            a) Spearman correlation (Sim)
%            b) Spearman correlation (Seq)
%            c) Linear fit (Sim)
%            d) Linear fit (Seq)
%            e) Combined linear model:
%                 Duration ~ SpikeCount + Condition
%
%   - Purpose:
%       Check whether the Spike Count - Duration relationship is strong
%       enough to justify later matched-response analysis:
%           * same spike count -> compare duration
%           * same duration -> compare spike count
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
file_paths = {

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX005/Result_Combined_SpikeDuration_Separate_DX005_Xia_Exp1_Sim.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX006/Result_Combined_SpikeDuration_Separate_DX006_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX006/Result_Combined_SpikeDuration_Separate_DX006_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX006/Result_Combined_SpikeDuration_Separate_DX006_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX006/Result_Combined_SpikeDuration_Separate_DX006_Xia_Exp1_Sim4.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX009/Result_Combined_SpikeDuration_Separate_DX009_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX009/Result_Combined_SpikeDuration_Separate_DX009_Xia_Exp1_Sim5.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX010/Result_Combined_SpikeDuration_Separate_DX010_Xia_Exp1_Sim8.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX011/Result_Combined_SpikeDuration_Separate_DX011_Xia_Exp1_Sim9.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX012/Result_Combined_SpikeDuration_Separate_DX012_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX012/Result_Combined_SpikeDuration_Separate_DX012_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX012/Result_Combined_SpikeDuration_Separate_DX012_Xia_Exp1_Sim6.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX013/Result_Combined_SpikeDuration_DX013_Xia_Exp1_Seq_Sim8.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX014/Result_Combined_SpikeDuration_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX014/Result_Combined_SpikeDuration_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX014/Result_Combined_SpikeDuration_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX014/Result_Combined_SpikeDuration_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX014/Result_Combined_SpikeDuration_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX014/Result_Combined_SpikeDuration_DX014_Xia_Seq_Sim6.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX015/Result_Combined_SpikeDuration_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX015/Result_Combined_SpikeDuration_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX015/Result_Combined_SpikeDuration_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX015/Result_Combined_SpikeDuration_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX015/Result_Combined_SpikeDuration_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX015/Result_Combined_SpikeDuration_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX015/Result_Combined_SpikeDuration_DX015_Xia_Seq_Sim7.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX016/Result_Combined_SpikeDuration_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX016/Result_Combined_SpikeDuration_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX016/Result_Combined_SpikeDuration_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_Duration/DX016/Result_Combined_SpikeDuration_DX016_Xia_Exp1_Seq_Full_4.mat';

};

save_figures = false;
save_dir = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/SpikeCount_Duration';
dot_size = 14;
jitter_w = 0.12;

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% ================= 1. POOL CHANNEL-LEVEL DATA =================
fprintf('Pooling combined spike count + duration data from %d files...\n', length(file_paths));

Pool = struct();
Pool.Sim.Spike = [];
Pool.Sim.Duration = [];
Pool.Sim.Amp = [];
Pool.Sim.FileID = [];
Pool.Sim.SetID = [];
Pool.Sim.ChannelID = [];

Pool.Seq.Spike = [];
Pool.Seq.Duration = [];
Pool.Seq.Amp = [];
Pool.Seq.FileID = [];
Pool.Seq.SetID = [];
Pool.Seq.ChannelID = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file')
        fprintf('  Missing file: %s\n', file_paths{f});
        continue;
    end

    D = load(file_paths{f});
    if ~isfield(D, 'CombinedData')
        fprintf('  Skipped file (no CombinedData): %s\n', file_paths{f});
        continue;
    end

    Data = D.CombinedData;

    if ~isfield(Data, 'Set') || isempty(Data.Set)
        continue;
    end

    for ss = 1:length(Data.Set)
        if ~isfield(Data.Set(ss), 'Amp') || isempty(Data.Set(ss).Amp)
            continue;
        end

        for ai = 1:length(Data.Set(ss).Amp)
            if ~isfield(Data.Set(ss).Amp(ai), 'Val')
                continue;
            end

            curr_amp = Data.Set(ss).Amp(ai).Val;

            % --------------------------------------------------------
            % IMPORTANT FILTER:
            % Ignore amp = 0 baseline points for this population analysis
            % --------------------------------------------------------
            if curr_amp == 0
                continue;
            end

            % -------------------------
            % Sim condition
            % -------------------------
            if isfield(Data.Set(ss).Amp(ai), 'Sim') && ...
               isfield(Data.Set(ss).Amp(ai).Sim, 'SpikeCount_All') && ...
               isfield(Data.Set(ss).Amp(ai).Sim, 'Duration_All')

                spike_sim = Data.Set(ss).Amp(ai).Sim.SpikeCount_All(:);
                dur_sim   = Data.Set(ss).Amp(ai).Sim.Duration_All(:);

                valid_sim = ~isnan(spike_sim) & ~isnan(dur_sim);

                Pool.Sim.Spike    = [Pool.Sim.Spike; spike_sim(valid_sim)];
                Pool.Sim.Duration = [Pool.Sim.Duration; dur_sim(valid_sim)];
                Pool.Sim.Amp      = [Pool.Sim.Amp; curr_amp * ones(sum(valid_sim),1)];
                Pool.Sim.FileID   = [Pool.Sim.FileID; f * ones(sum(valid_sim),1)];
                Pool.Sim.SetID    = [Pool.Sim.SetID; ss * ones(sum(valid_sim),1)];
                Pool.Sim.ChannelID = [Pool.Sim.ChannelID; find(valid_sim)];
            end

            % -------------------------
            % Seq condition
            % -------------------------
            if isfield(Data.Set(ss).Amp(ai), 'Seq') && ...
               isfield(Data.Set(ss).Amp(ai).Seq, 'SpikeCount_All') && ...
               isfield(Data.Set(ss).Amp(ai).Seq, 'Duration_All')

                spike_seq = Data.Set(ss).Amp(ai).Seq.SpikeCount_All(:);
                dur_seq   = Data.Set(ss).Amp(ai).Seq.Duration_All(:);

                valid_seq = ~isnan(spike_seq) & ~isnan(dur_seq);

                Pool.Seq.Spike    = [Pool.Seq.Spike; spike_seq(valid_seq)];
                Pool.Seq.Duration = [Pool.Seq.Duration; dur_seq(valid_seq)];
                Pool.Seq.Amp      = [Pool.Seq.Amp; curr_amp * ones(sum(valid_seq),1)];
                Pool.Seq.FileID   = [Pool.Seq.FileID; f * ones(sum(valid_seq),1)];
                Pool.Seq.SetID    = [Pool.Seq.SetID; ss * ones(sum(valid_seq),1)];
                Pool.Seq.ChannelID = [Pool.Seq.ChannelID; find(valid_seq)];
            end
        end
    end
end

fprintf('Pooling complete.\n');
fprintf('  Sim valid points: %d\n', length(Pool.Sim.Spike));
fprintf('  Seq valid points: %d\n', length(Pool.Seq.Spike));

%% ================= 2. BASIC SANITY CHECK =================
if isempty(Pool.Sim.Spike) || isempty(Pool.Seq.Spike)
    error('No valid Sim or Seq data points found after filtering.');
end

%% ================= 3. SPEARMAN CORRELATION =================
fprintf('\n====================================================================\n');
fprintf('           POPULATION FEASIBILITY: SPIKE COUNT vs DURATION          \n');
fprintf('====================================================================\n');

[rho_sim, p_sim] = corr(Pool.Sim.Spike, Pool.Sim.Duration, 'Type', 'Spearman', 'Rows', 'complete');
[rho_seq, p_seq] = corr(Pool.Seq.Spike, Pool.Seq.Duration, 'Type', 'Spearman', 'Rows', 'complete');

fprintf('Spearman Correlation:\n');
fprintf('  Simultaneous: rho = %.4f, p = %s, N = %d\n', rho_sim, format_p(p_sim), length(Pool.Sim.Spike));
fprintf('  Sequential:   rho = %.4f, p = %s, N = %d\n', rho_seq, format_p(p_seq), length(Pool.Seq.Spike));

%% ================= 4. SIMPLE LINEAR FITS (VISUAL CHECK) =================
% ------------------------------------------------------------
% This is only a descriptive fit for quick visual inspection.
% It is NOT the final inferential test.
% ------------------------------------------------------------
fit_sim = polyfit(Pool.Sim.Spike, Pool.Sim.Duration, 1);
fit_seq = polyfit(Pool.Seq.Spike, Pool.Seq.Duration, 1);

xq_sim = linspace(min(Pool.Sim.Spike), max(Pool.Sim.Spike), 200);
yq_sim = polyval(fit_sim, xq_sim);

xq_seq = linspace(min(Pool.Seq.Spike), max(Pool.Seq.Spike), 200);
yq_seq = polyval(fit_seq, xq_seq);

mdl_sim = fitlm(Pool.Sim.Spike, Pool.Sim.Duration);
mdl_seq = fitlm(Pool.Seq.Spike, Pool.Seq.Duration);

fprintf('\nSimple Linear Fit:\n');
fprintf('  Simultaneous: Duration = %.4f * Spike + %.4f | R^2 = %.4f\n', ...
    fit_sim(1), fit_sim(2), mdl_sim.Rsquared.Ordinary);
fprintf('  Sequential:   Duration = %.4f * Spike + %.4f | R^2 = %.4f\n', ...
    fit_seq(1), fit_seq(2), mdl_seq.Rsquared.Ordinary);

%% ================= 5. COMBINED LINEAR MODEL =================
% ------------------------------------------------------------
% Main feasibility test:
% Does condition still explain duration after accounting for spike count?
%
% Model:
%   Duration ~ SpikeCount + Condition
%
% Condition coding:
%   Sim = 0
%   Seq = 1
% ------------------------------------------------------------
AllSpike = [Pool.Sim.Spike; Pool.Seq.Spike];
AllDur   = [Pool.Sim.Duration; Pool.Seq.Duration];
AllCond  = [zeros(length(Pool.Sim.Spike),1); ones(length(Pool.Seq.Spike),1)];

Tbl = table(AllSpike, AllDur, categorical(AllCond), ...
    'VariableNames', {'SpikeCount', 'Duration', 'Condition'});

mdl_all = fitlm(Tbl, 'Duration ~ SpikeCount + Condition');

fprintf('\nCombined Linear Model:\n');
disp(mdl_all);

fprintf('Model R^2 = %.4f\n', mdl_all.Rsquared.Ordinary);

coef_table = mdl_all.Coefficients;
fprintf('\nKey Coefficients:\n');
for i = 1:height(coef_table)
    fprintf('  %-20s | Beta = %8.4f | SE = %8.4f | t = %8.4f | p = %s\n', ...
        coef_table.Properties.RowNames{i}, ...
        coef_table.Estimate(i), ...
        coef_table.SE(i), ...
        coef_table.tStat(i), ...
        format_p(coef_table.pValue(i)));
end

%% ================= 6. FIGURE 1: POOLED SCATTER =================
fprintf('\nGenerating pooled scatter plot...\n');

fig1 = figure('Color','w', 'Units', 'centimeters', ...
    'Position', [5, 5, 8.8, 8.8], 'Name', 'Population_Spike_vs_Duration'); 
hold on;

% --- Raw pooled points ---
scatter(Pool.Sim.Spike + (rand(size(Pool.Sim.Spike))-0.5)*jitter_w, ...
        Pool.Sim.Duration, ...
        dot_size, [0.75 0.75 0.75], 'o', ...
        'LineWidth', 0.6, 'MarkerEdgeAlpha', 0.35, ...
        'HandleVisibility', 'off');

scatter(Pool.Seq.Spike + (rand(size(Pool.Seq.Spike))-0.5)*jitter_w, ...
        Pool.Seq.Duration, ...
        dot_size, [0.55 0.55 0.55], 's', ...
        'filled', 'MarkerFaceAlpha', 0.35, ...
        'HandleVisibility', 'off');

% --- Fitted lines ---
p1 = plot(xq_sim, yq_sim, '--k', 'LineWidth', 1.5, 'DisplayName', 'Simultaneous fit');
p2 = plot(xq_seq, yq_seq, '-k', 'LineWidth', 1.5, 'DisplayName', 'Sequential fit');

xlabel('Spike Count', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');

legend([p1, p2], {'Simultaneous fit', 'Sequential fit'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

set(gca, 'FontSize', 9, 'FontName', 'Arial', ...
    'LineWidth', 1.0, 'TickDir', 'out');
box off; axis square;

if save_figures
    exportgraphics(fig1, fullfile(save_dir, 'Population_Spike_vs_Duration.tiff'), ...
        'Resolution', 600);
end

%% ================= 7. FIGURE 2: CONDITION-COLORED MEAN BIN PLOT =================
% ------------------------------------------------------------
% Optional helper plot:
% Bin spike count and show mean duration ± SEM in each bin
% This is often easier to read than raw scatter alone.
% ------------------------------------------------------------
fprintf('Generating binned mean plot...\n');

bin_edges = 0:1:max(ceil(AllSpike));
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

sim_bin_mean = nan(size(bin_centers));
sim_bin_sem  = nan(size(bin_centers));
seq_bin_mean = nan(size(bin_centers));
seq_bin_sem  = nan(size(bin_centers));

for b = 1:length(bin_centers)
    % --- Sim ---
    idx_sim = Pool.Sim.Spike >= bin_edges(b) & Pool.Sim.Spike < bin_edges(b+1);
    y_sim = Pool.Sim.Duration(idx_sim);

    if ~isempty(y_sim)
        sim_bin_mean(b) = mean(y_sim, 'omitnan');
        % sim_bin_sem(b)  = std(y_sim, 'omitnan') / sqrt(sum(~isnan(y_sim)));
        sim_bin_sem(b)  = std(y_sim, 'omitnan') / sqrt(10);
    end

    % --- Seq ---
    idx_seq = Pool.Seq.Spike >= bin_edges(b) & Pool.Seq.Spike < bin_edges(b+1);
    y_seq = Pool.Seq.Duration(idx_seq);

    if ~isempty(y_seq)
        seq_bin_mean(b) = mean(y_seq, 'omitnan');
        % seq_bin_sem(b)  = std(y_seq, 'omitnan') / sqrt(sum(~isnan(y_seq)));
        seq_bin_sem(b)  = std(y_seq, 'omitnan') / sqrt(10);
    end
end

fig2 = figure('Color','w', 'Units', 'centimeters', ...
    'Position', [15, 5, 8.8, 8.8], 'Name', 'Population_Spike_vs_Duration_Binned'); 
hold on;

valid_sim = ~isnan(sim_bin_mean);
valid_seq = ~isnan(seq_bin_mean);

errorbar(bin_centers(valid_sim), sim_bin_mean(valid_sim), sim_bin_sem(valid_sim), ...
    '--ok', 'LineWidth', 1, 'MarkerFaceColor', 'w', 'CapSize', 8, ...
    'DisplayName', 'Simultaneous');

errorbar(bin_centers(valid_seq), seq_bin_mean(valid_seq), seq_bin_sem(valid_seq), ...
    '-sk', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'CapSize', 8, ...
    'DisplayName', 'Sequential');

xlabel('Spike Count Bin', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Mean Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

set(gca, 'FontSize', 9, 'FontName', 'Arial', ...
    'LineWidth', 1.0, 'TickDir', 'out');
box off; axis square;

if save_figures
    exportgraphics(fig2, fullfile(save_dir, 'Population_Spike_vs_Duration_Binned.tiff'), ...
        'Resolution', 600);
end

%% ================= 8. COMMAND WINDOW INTERPRETATION GUIDE =================
fprintf('\n====================================================================\n');
fprintf('INTERPRETATION GUIDE\n');
fprintf('====================================================================\n');
fprintf('1. Spearman rho:\n');
fprintf('   Shows whether spike count and duration are monotonically related.\n');
fprintf('   Positive rho means larger spike count tends to come with longer duration.\n');
fprintf('\n');
fprintf('2. Separate linear fits:\n');
fprintf('   These are descriptive only.\n');
fprintf('   They help visually compare the trend for Sim and Seq.\n');
fprintf('\n');
fprintf('3. Combined linear model:\n');
fprintf('   Duration ~ SpikeCount + Condition\n');
fprintf('   This is the key feasibility test.\n');
fprintf('   If the Condition term is significant, then Sim vs Seq differ in duration\n');
fprintf('   even after accounting for spike count.\n');
fprintf('====================================================================\n');

%% ================= HELPER FUNCTION =================
function p_str = format_p(p)
    if isnan(p)
        p_str = 'NaN';
    else
        p_str = sprintf('%.20f', p);
    end
end