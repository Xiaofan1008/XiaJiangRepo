%% ============================================================
%   POPULATION ANALYSIS: SPIKE COUNT vs DURATION (SET LEVEL)
%   - Input: Result_Combined_SpikeDuration_*.mat files
%   - Logic:
%       1. Load all single-dataset combined result files.
%       2. Pool SET-level data across datasets.
%       3. Ignore amp = 0 points.
%       4. Use set-level mean Spike Count and mean Duration
%          from the saved responsive pool.
%       5. Normalize Spike Count within each file/set
%          using Ref_Amp and max(Sim, Seq) reference response.
%       6. Plot Spike Count vs Duration for Sim and Seq.
%       7. Run feasibility statistics.
%
%   - NOTE:
%       This version uses SET-LEVEL means rather than channel-level points.
%       Spike Count is normalized within each file/set.
%       Duration remains unnormalized (ms).
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

save_figures = true;
save_dir = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/SpikeCount_Duration';
dot_size = 18;
jitter_w = 0.10;

% --- NEW: spike-count normalization settings ---
Ref_Amp = 5;             % reference amplitude for within-set normalization
min_ref_response = 0.01; % avoid dividing by tiny numbers

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% ================= 1. POOL SET-LEVEL DATA =================
fprintf('Pooling combined spike count + duration data from %d files (SET LEVEL)...\n', length(file_paths));

Pool = struct();
Pool.Sim.Spike = [];
Pool.Sim.Duration = [];
Pool.Sim.Amp = [];
Pool.Sim.FileID = [];
Pool.Sim.SetID = [];

Pool.Seq.Spike = [];
Pool.Seq.Duration = [];
Pool.Seq.Amp = [];
Pool.Seq.FileID = [];
Pool.Seq.SetID = [];

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

        % =========================================================
        % NEW: collect all set-level spike means first, so we can
        % calculate one normalization denominator per file/set
        % =========================================================
        amps_this_set = [];
        sim_spike_this_set = [];
        seq_spike_this_set = [];

        for ai = 1:length(Data.Set(ss).Amp)
            if ~isfield(Data.Set(ss).Amp(ai), 'Val')
                continue;
            end

            curr_amp = Data.Set(ss).Amp(ai).Val;

            if curr_amp == 0
                continue;
            end

            amps_this_set(end+1,1) = curr_amp; %#ok<SAGROW>

            if isfield(Data.Set(ss).Amp(ai), 'Sim') && ...
               isfield(Data.Set(ss).Amp(ai).Sim, 'SpikeCount_Mean')
                sim_spike_this_set(end+1,1) = Data.Set(ss).Amp(ai).Sim.SpikeCount_Mean; %#ok<SAGROW>
            else
                sim_spike_this_set(end+1,1) = NaN; %#ok<SAGROW>
            end

            if isfield(Data.Set(ss).Amp(ai), 'Seq') && ...
               isfield(Data.Set(ss).Amp(ai).Seq, 'SpikeCount_Mean')
                seq_spike_this_set(end+1,1) = Data.Set(ss).Amp(ai).Seq.SpikeCount_Mean; %#ok<SAGROW>
            else
                seq_spike_this_set(end+1,1) = NaN; %#ok<SAGROW>
            end
        end

        % --- determine within-set normalization denominator ---
        denom = NaN;

        if ~isempty(amps_this_set)
            [amps_this_set, sort_idx] = sort(amps_this_set);
            sim_spike_this_set = sim_spike_this_set(sort_idx);
            seq_spike_this_set = seq_spike_this_set(sort_idx);

            val_sim_ref = NaN;
            val_seq_ref = NaN;

            ref_idx = find(abs(amps_this_set - Ref_Amp) < 0.001);

            if ~isempty(ref_idx)
                val_sim_ref = sim_spike_this_set(ref_idx(1));
                val_seq_ref = seq_spike_this_set(ref_idx(1));
            else
                valid_sim_ref = ~isnan(sim_spike_this_set);
                if sum(valid_sim_ref) >= 2
                    val_sim_ref = interp1(amps_this_set(valid_sim_ref), sim_spike_this_set(valid_sim_ref), Ref_Amp, 'linear', 'extrap');
                elseif sum(valid_sim_ref) == 1
                    val_sim_ref = sim_spike_this_set(valid_sim_ref);
                end

                valid_seq_ref = ~isnan(seq_spike_this_set);
                if sum(valid_seq_ref) >= 2
                    val_seq_ref = interp1(amps_this_set(valid_seq_ref), seq_spike_this_set(valid_seq_ref), Ref_Amp, 'linear', 'extrap');
                elseif sum(valid_seq_ref) == 1
                    val_seq_ref = seq_spike_this_set(valid_seq_ref);
                end
            end

            ref_vals = [max(0, val_sim_ref), max(0, val_seq_ref)];
            if any(~isnan(ref_vals))
                denom = max(ref_vals, [], 'omitnan');
            end
        end

        % =========================================================
        % original pooling loop, now using normalized spike count
        % =========================================================
        for ai = 1:length(Data.Set(ss).Amp)
            if ~isfield(Data.Set(ss).Amp(ai), 'Val')
                continue;
            end

            curr_amp = Data.Set(ss).Amp(ai).Val;

            if curr_amp == 0
                continue;
            end

            % -------------------------
            % Sim condition (SET LEVEL)
            % -------------------------
            if isfield(Data.Set(ss).Amp(ai), 'Sim') && ...
               isfield(Data.Set(ss).Amp(ai).Sim, 'SpikeCount_Mean') && ...
               isfield(Data.Set(ss).Amp(ai).Sim, 'Duration_Mean')

                spike_sim = Data.Set(ss).Amp(ai).Sim.SpikeCount_Mean;
                dur_sim   = Data.Set(ss).Amp(ai).Sim.Duration_Mean;

                % NEW: normalize spike count within file/set
                if ~isnan(spike_sim) && ~isnan(denom) && denom > min_ref_response
                    spike_sim = spike_sim / denom;
                else
                    spike_sim = NaN;
                end

                if ~isnan(spike_sim) && ~isnan(dur_sim)
                    Pool.Sim.Spike    = [Pool.Sim.Spike; spike_sim];
                    Pool.Sim.Duration = [Pool.Sim.Duration; dur_sim];
                    Pool.Sim.Amp      = [Pool.Sim.Amp; curr_amp];
                    Pool.Sim.FileID   = [Pool.Sim.FileID; f];
                    Pool.Sim.SetID    = [Pool.Sim.SetID; ss];
                end
            end

            % -------------------------
            % Seq condition (SET LEVEL)
            % -------------------------
            if isfield(Data.Set(ss).Amp(ai), 'Seq') && ...
               isfield(Data.Set(ss).Amp(ai).Seq, 'SpikeCount_Mean') && ...
               isfield(Data.Set(ss).Amp(ai).Seq, 'Duration_Mean')

                spike_seq = Data.Set(ss).Amp(ai).Seq.SpikeCount_Mean;
                dur_seq   = Data.Set(ss).Amp(ai).Seq.Duration_Mean;

                % NEW: normalize spike count within file/set
                if ~isnan(spike_seq) && ~isnan(denom) && denom > min_ref_response
                    spike_seq = spike_seq / denom;
                else
                    spike_seq = NaN;
                end

                if ~isnan(spike_seq) && ~isnan(dur_seq)
                    Pool.Seq.Spike    = [Pool.Seq.Spike; spike_seq];
                    Pool.Seq.Duration = [Pool.Seq.Duration; dur_seq];
                    Pool.Seq.Amp      = [Pool.Seq.Amp; curr_amp];
                    Pool.Seq.FileID   = [Pool.Seq.FileID; f];
                    Pool.Seq.SetID    = [Pool.Seq.SetID; ss];
                end
            end
        end
    end
end

fprintf('Pooling complete.\n');
fprintf('  Sim valid SET-level points: %d\n', length(Pool.Sim.Spike));
fprintf('  Seq valid SET-level points: %d\n', length(Pool.Seq.Spike));

%% ================= 2. BASIC SANITY CHECK =================
if isempty(Pool.Sim.Spike) || isempty(Pool.Seq.Spike)
    error('No valid Sim or Seq set-level data points found after filtering.');
end

%% ================= 3. SPEARMAN CORRELATION =================
fprintf('\n====================================================================\n');
fprintf('       POPULATION FEASIBILITY: SPIKE COUNT vs DURATION (SET LEVEL)  \n');
fprintf('====================================================================\n');

[rho_sim, p_sim] = corr(Pool.Sim.Spike, Pool.Sim.Duration, 'Type', 'Spearman', 'Rows', 'complete');
[rho_seq, p_seq] = corr(Pool.Seq.Spike, Pool.Seq.Duration, 'Type', 'Spearman', 'Rows', 'complete');

fprintf('Spearman Correlation:\n');
fprintf('  Simultaneous: rho = %.4f, p = %s, N = %d\n', rho_sim, format_p(p_sim), length(Pool.Sim.Spike));
fprintf('  Sequential:   rho = %.4f, p = %s, N = %d\n', rho_seq, format_p(p_seq), length(Pool.Seq.Spike));

%% ================= 4. SIMPLE LINEAR FITS (VISUAL CHECK) =================
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
    'Position', [5, 5, 8.8, 8.8], 'Name', 'Population_Spike_vs_Duration_SetLevel'); 
hold on;

scatter(Pool.Sim.Spike + (rand(size(Pool.Sim.Spike))-0.5)*jitter_w, ...
        Pool.Sim.Duration, ...
        dot_size, [0.75 0.75 0.75], 'o', ...
        'LineWidth', 0.6, 'MarkerEdgeAlpha', 0.45, ...
        'HandleVisibility', 'off');

scatter(Pool.Seq.Spike + (rand(size(Pool.Seq.Spike))-0.5)*jitter_w, ...
        Pool.Seq.Duration, ...
        dot_size, [0.55 0.55 0.55], 's', ...
        'filled', 'MarkerFaceAlpha', 0.40, ...
        'HandleVisibility', 'off');

p1 = plot(xq_sim, yq_sim, '--k', 'LineWidth', 1.5, 'DisplayName', 'Simultaneous fit');
p2 = plot(xq_seq, yq_seq, '-k', 'LineWidth', 1.5, 'DisplayName', 'Sequential fit');

xlabel('Normalized Spike Count', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Duration (ms)', 'FontSize', 9, 'FontName', 'Arial');

legend([p1, p2], {'Simultaneous fit', 'Sequential fit'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

set(gca, 'FontSize', 9, 'FontName', 'Arial', ...
    'LineWidth', 1.0, 'TickDir', 'out');
box off; axis square;

if save_figures
    exportgraphics(fig1, fullfile(save_dir, 'Population_Spike_vs_Duration_Scatter.tiff'), ...
        'Resolution', 600);
end

%% ================= 7. FIGURE 2: CONDITION-COLORED MEAN BIN PLOT =================
fprintf('Generating binned mean plot...\n');

% bin_edges = 0:0.5:max(ceil(AllSpike*4))/4;
bin_edges = 0:0.2:max(ceil(1.2*4))/4;

bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

sim_bin_mean = nan(size(bin_centers));
sim_bin_sem  = nan(size(bin_centers));
seq_bin_mean = nan(size(bin_centers));
seq_bin_sem  = nan(size(bin_centers));

for b = 1:length(bin_centers)
    idx_sim = Pool.Sim.Spike >= bin_edges(b) & Pool.Sim.Spike < bin_edges(b+1);
    y_sim = Pool.Sim.Duration(idx_sim);

    if ~isempty(y_sim)
        sim_bin_mean(b) = mean(y_sim, 'omitnan');
        sim_bin_sem(b)  = std(y_sim, 'omitnan') / sqrt(sum(~isnan(y_sim)));
    end

    idx_seq = Pool.Seq.Spike >= bin_edges(b) & Pool.Seq.Spike < bin_edges(b+1);
    y_seq = Pool.Seq.Duration(idx_seq);

    if ~isempty(y_seq)
        seq_bin_mean(b) = mean(y_seq, 'omitnan');
        seq_bin_sem(b)  = std(y_seq, 'omitnan') / sqrt(sum(~isnan(y_seq)));
    end
end

fig2 = figure('Color','w', 'Units', 'centimeters', ...
    'Position', [15, 5, 8.8, 8.8], 'Name', 'Population_Spike_vs_Duration_Binned_SetLevel'); 
hold on;

valid_sim = ~isnan(sim_bin_mean);
valid_seq = ~isnan(seq_bin_mean);

errorbar(bin_centers(valid_sim), sim_bin_mean(valid_sim), sim_bin_sem(valid_sim), ...
    '--ok', 'LineWidth', 1, 'MarkerFaceColor', 'w', 'CapSize', 8, ...
    'DisplayName', 'Simultaneous');

errorbar(bin_centers(valid_seq), seq_bin_mean(valid_seq), seq_bin_sem(valid_seq), ...
    '-sk', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'CapSize', 8, ...
    'DisplayName', 'Sequential');

xlabel('Normalized Response Magnitude', 'FontSize', 9, 'FontName', 'Arial');
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
fprintf('\n');
fprintf('4. This version uses SET-LEVEL means, not channel-level points.\n');
fprintf('   Spike count is normalized within each file/set.\n');
fprintf('   Duration remains in ms.\n');
fprintf('====================================================================\n');

%% ================= HELPER FUNCTION =================
function p_str = format_p(p)
    if isnan(p)
        p_str = 'NaN';
    else
        p_str = sprintf('%.20f', p);
    end
end