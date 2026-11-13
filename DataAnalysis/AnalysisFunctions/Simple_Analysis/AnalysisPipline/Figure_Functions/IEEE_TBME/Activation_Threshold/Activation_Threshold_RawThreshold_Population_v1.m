%% ============================================================
%   GRAND AVERAGE ACTIVATION THRESHOLD (RAW SPIKE METRIC)
%   - Metric: Amplitude (uA) required to reach X spikes per trial
%   - Logic: 
%       1. Uses RAW spike counts (not normalized) to avoid grid artifacts.
%       2. Anchor curves at 0 uA = 0 spikes.
%       3. Strict Exclusion: BOTH conditions must reach Target_Thresh.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS ============================
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

% [MODIFIED 1] Target is now in absolute units: Spikes per Trial
% A value of 1.0 is a very strong, standard activation threshold.
Target_Thresh = 0.35;   

% Plot Settings
dot_size = 25;
dot_color = [0 0 0]; % Deep professional blue
dot_alpha = 0.25;           

%% =================== 2. MASTER LOOP INITIALIZATION ====================
Thresh_Sim_Pool = [];
Thresh_Seq_Pool = [];
total_sets_processed = 0; 

fprintf('Starting Raw Spike Threshold Extraction (Target = %.1f spikes/trial)...\n', Target_Thresh);

for f_idx = 1:length(file_paths)
    data_file = file_paths{f_idx};
    if ~exist(data_file, 'file'), continue; end

    D = load(data_file);
    if ~isfield(D, 'ResultNorm'), continue; end
    
    R = D.ResultNorm;
    Amps_Raw = R.Amps(:)'; % [FIXED] Force horizontal
    
    % [MODIFIED 2] Switched from .Norm to .Raw matrices
    Raw_Sim_Data = R.Raw_Sim; 
    Raw_Seq_Data = R.Raw_Seq;

    [nCh_Total, nAmps_Raw, nSets] = size(Raw_Sim_Data);
    total_sets_processed = total_sets_processed + nSets;

    %% =================== 3. ANCHOR TO ZERO ====================
    if Amps_Raw(1) ~= 0
        Amps = [0, Amps_Raw];
        pad_zeros = zeros(nCh_Total, 1, nSets);
        % Append zeros to the RAW data
        Data_Sim = cat(2, pad_zeros, Raw_Sim_Data);
        Data_Seq = cat(2, pad_zeros, Raw_Seq_Data);
    else
        Amps = Amps_Raw;
        Data_Sim = Raw_Sim_Data;
        Data_Seq = Raw_Seq_Data;
    end

    %% =================== 4. CALCULATE THRESHOLDS ====================
    for ss = 1:nSets
        for ch = 1:nCh_Total
            curve_sim = Data_Sim(ch, :, ss);
            curve_seq = Data_Seq(ch, :, ss);
            
            if all(isnan(curve_sim)) || all(isnan(curve_seq)), continue; end
            
            t_sim = get_interpolated_threshold(Amps, curve_sim, Target_Thresh);
            t_seq = get_interpolated_threshold(Amps, curve_seq, Target_Thresh);
            
            % Strict Exclusion
            if isnan(t_sim) || isnan(t_seq), continue; end
            
            Thresh_Sim_Pool = [Thresh_Sim_Pool; t_sim];
            Thresh_Seq_Pool = [Thresh_Seq_Pool; t_seq];
        end
    end
end

%% =================== 5. SUMMARY STATS ====================
nPairs = length(Thresh_Sim_Pool);
if nPairs == 0, fprintf('No channels reached %.1f spikes.\n', Target_Thresh); return; end
p_val = signrank(Thresh_Sim_Pool, Thresh_Seq_Pool);

fprintf('\n================== FINAL SUMMARY ==================\n');
fprintf('Threshold Definition   : %.1f Mean Spikes/Trial\n', Target_Thresh);
fprintf('Total Stim Pairs (Sets): %d\n', total_sets_processed);
fprintf('Valid Paired Ch (N)    : %d\n', nPairs);
fprintf('---------------------------------------------------\n');
fprintf('Median Sim Threshold   : %.3f uA\n', median(Thresh_Sim_Pool));
fprintf('Median Seq Threshold   : %.3f uA\n', median(Thresh_Seq_Pool));
fprintf('Wilcoxon Paired p      : %.2e\n', p_val);
fprintf('===================================================\n\n');

%% =================== 6. UNITY SCATTER PLOT ====================
figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w', 'PaperPositionMode', 'auto'); 
hold on;

max_val = ceil(max([Thresh_Sim_Pool; Thresh_Seq_Pool]));
if max_val < 5, max_val = 5; end 

plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
scatter(Thresh_Sim_Pool, Thresh_Seq_Pool, dot_size, dot_color, 'filled', ...
    'MarkerFaceAlpha', dot_alpha,'MarkerEdgeColor', 'none');

% Significance Formatting
if p_val < 0.001, star = '***'; elseif p_val < 0.01, star = '**'; elseif p_val < 0.05, star = '*'; else, star = 'n.s.'; end
title(['Activation Threshold (p < 0.001', star, ')'], 'FontSize', 9, 'FontName', 'Arial', 'FontWeight', 'normal');

axis square; box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.0);
xlabel('Simultaneous Threshold (\muA)'); ylabel('Sequential Threshold (\muA)');
xlim([0 max_val]); ylim([0 max_val]);
hold off;

%% ==================== HELPER FUNCTIONS =========================
function thresh = get_interpolated_threshold(amps, curve, target)
    amps = squeeze(amps(:))'; curve = squeeze(curve(:))';
    valid_idx = ~isnan(curve);
    valid_amps = amps(valid_idx); valid_curve = curve(valid_idx);
    
    if isempty(valid_curve) || max(valid_curve) < target, thresh = NaN; return; end
    
    idx_above = find(valid_curve >= target, 1, 'first');
    if valid_curve(idx_above) == target, thresh = valid_amps(idx_above); return; end
    
    if idx_above == 1
        x0 = 0; y0 = 0; x1 = valid_amps(1); y1 = valid_curve(1);
    else
        x0 = valid_amps(idx_above - 1); x1 = valid_amps(idx_above);
        y0 = valid_curve(idx_above - 1); y1 = valid_curve(idx_above);
    end
    thresh = x0 + (target - y0) * (x1 - x0) / (y1 - y0);
end