%% ============================================================
%   ISO-RESPONSE PLOT: REQUIRED AMPLITUDE FOR TARGET OUTPUT
%   - Logic: Sigmoidal (Logistic) Regression per dataset.
%   - Style: IEEE TBME Publication Style (Percentage Axis)
%   - Statistics: Paired T-Test with Command Window Printout
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% --- ISO-RESPONSE SETTINGS ---
% Targets from 0.1 to 2.0 (10% to 200%)
Target_Spikes = 0.2 : 0.2 : 2.0; 

% --- STATISTICAL FILTERS ---
stats_min_n_threshold = 1; 
plot_min_n = 4; % Truncate plot if N drops below this

% List all result files
file_paths = {
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

%% ================= 2. AGGREGATE DATA =================
Dataset_Amps = {}; Dataset_SimVals = {}; Dataset_SeqVals = {};
Pool_Flipped_Raw = [];

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    D = load(file_paths{i}); R = D.ResultNorm; Amps = R.Amps;
    sim_ds = squeeze(mean(mean(R.Norm_Sim, 3, 'omitnan'), 1, 'omitnan'));
    seq_raw = R.Norm_Seq;
    if ndims(seq_raw) == 4, seq_ds = squeeze(mean(mean(mean(seq_raw, 4, 'omitnan'), 3, 'omitnan'), 1, 'omitnan'));
    else, seq_ds = squeeze(mean(mean(seq_raw, 3, 'omitnan'), 1, 'omitnan')); end
    
    Dataset_Amps{end+1} = Amps; Dataset_SimVals{end+1} = sim_ds; Dataset_SeqVals{end+1} = seq_ds;
    for a = 1:length(Amps)
        if ~isnan(sim_ds(a)), Pool_Flipped_Raw = [Pool_Flipped_Raw; sim_ds(a), Amps(a), 1]; end
        if ~isnan(seq_ds(a)), Pool_Flipped_Raw = [Pool_Flipped_Raw; seq_ds(a), Amps(a), 2]; end
    end
end

%% ================= 3. SIGMOIDAL FITTING ENGINE =================
fprintf('Fitting 4-Parameter Sigmoids to each dataset...\n');

% Sigmoid Function: y = Base + (Max-Base) / (1 + exp(-k * (x - V50)))
% beta = [Base, Max, k, V50]
sigmoid = @(beta, x) beta(1) + (beta(2) - beta(1)) ./ (1 + exp(-beta(3) * (x - beta(4))));

Cost_Sim = nan(length(Target_Spikes), length(Dataset_Amps));
Cost_Seq = nan(length(Target_Spikes), length(Dataset_Amps));

opts = statset('Display','off','MaxIter',500);

for d = 1:length(Dataset_Amps)
    amps = Dataset_Amps{d}(:);
    sim_s = Dataset_SimVals{d}(:);
    seq_s = Dataset_SeqVals{d}(:);
    
    % Initial guesses for [Base, Max, Slope, V50]
    beta0 = [0, max(sim_s), 1, mean(amps)];
    
    try
        % Fit Simultaneous
        beta_sim = nlinfit(amps, sim_s, sigmoid, beta0, opts);
        % Solve Inverse: Amp = V50 - (1/k) * log((Max-Base)/(S-Base) - 1)
        for t = 1:length(Target_Spikes)
            S = Target_Spikes(t);
            if S < beta_sim(2) && S > beta_sim(1)
                Cost_Sim(t, d) = beta_sim(4) - (1/beta_sim(3)) * log((beta_sim(2)-beta_sim(1))/(S-beta_sim(1)) - 1);
            end
        end
    catch, end
    
    try
        % Fit Sequential
        beta0_seq = [0, max(seq_s), 1, mean(amps)];
        beta_seq = nlinfit(amps, seq_s, sigmoid, beta0_seq, opts);
        for t = 1:length(Target_Spikes)
            S = Target_Spikes(t);
            if S < beta_seq(2) && S > beta_seq(1)
                Cost_Seq(t, d) = beta_seq(4) - (1/beta_seq(3)) * log((beta_seq(2)-beta_seq(1))/(S-beta_seq(1)) - 1);
            end
        end
    catch, end
end

% Summary Stats
Grand_Sim_Mean = mean(Cost_Sim, 2, 'omitnan');
Grand_Sim_SEM  = std(Cost_Sim, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(Cost_Sim), 2));
Grand_Seq_Mean = mean(Cost_Seq, 2, 'omitnan');
Grand_Seq_SEM  = std(Cost_Seq, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(Cost_Seq), 2));

% Apply Truncation Filter
N_Sim = sum(~isnan(Cost_Sim), 2);
N_Seq = sum(~isnan(Cost_Seq), 2);
Grand_Sim_Mean(N_Sim < plot_min_n) = NaN; Grand_Sim_SEM(N_Sim < plot_min_n) = NaN;
Grand_Seq_Mean(N_Seq < plot_min_n) = NaN; Grand_Seq_SEM(N_Seq < plot_min_n) = NaN;

%% ================= 4. PLOT & STATS =================
figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w'); hold on;
Target_P = Target_Spikes * 100;

% Background Scatter
scatter(Pool_Flipped_Raw(Pool_Flipped_Raw(:,3)==2,1)*100, Pool_Flipped_Raw(Pool_Flipped_Raw(:,3)==2,2), 8, [0.7 0.7 0.7], 's', 'filled', 'MarkerFaceAlpha', 0.2);
scatter(Pool_Flipped_Raw(Pool_Flipped_Raw(:,3)==1,1)*100, Pool_Flipped_Raw(Pool_Flipped_Raw(:,3)==1,2), 8, [0.8 0.8 0.8], 'o', 'MarkerEdgeAlpha', 0.2);

% Main Curves
errorbar(Target_P, Grand_Sim_Mean, Grand_Sim_SEM, '.', 'Color', 'k', 'CapSize', 4);
p1 = plot([0, Target_P], [0; Grand_Sim_Mean], '--o', 'Color', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous');
errorbar(Target_P, Grand_Seq_Mean, Grand_Seq_SEM, '.', 'Color', 'k', 'CapSize', 4);
p2 = plot([0, Target_P], [0; Grand_Seq_Mean], '-s', 'Color', 'k', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Sequential');

% Stats loop
fprintf('\n=== STATS ===\n');
for t = 1:length(Target_Spikes)
    [~, p] = ttest(Cost_Sim(t,:), Cost_Seq(t,:));
    if ~isnan(p)
        fprintf('Target %.0f%% | p = %.5f\n', Target_P(t), p);
        txt = ''; if p<0.001, txt='***'; elseif p<0.01, txt='**'; elseif p<0.05, txt='*'; end
        if ~isempty(txt)
            text(Target_P(t), (Grand_Sim_Mean(t)+Grand_Seq_Mean(t))/2, txt, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
end

% Formatting
box off; set(gca, 'FontSize', 9, 'FontName', 'Arial', 'LineWidth', 1.2);
xlabel('Relative Neural Recruitment (%)'); ylabel('Required Amplitude (\muA)');
legend([p1, p2], 'Location', 'northwest', 'Box', 'off');
xlim([0 200]); ylim([0 12]);