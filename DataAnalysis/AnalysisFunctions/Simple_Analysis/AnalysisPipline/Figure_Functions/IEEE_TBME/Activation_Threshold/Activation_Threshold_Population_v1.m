%% ============================================================
%   GRAND AVERAGE ACTIVATION THRESHOLD EXTRACTOR
%   - Metric: Amplitude (uA) required to reach Target Normalized Response
%   - Logic: 
%       1. Loops through ALL datasets.
%       2. Anchor curves at 0 uA = 0 response.
%       3. Interpolate linearly to find first crossing of Target_Thresh.
%       4. Require BOTH Sim and Seq to cross target (Strict Exclusion).
%       5. Plot paired Unity Scatter & Print Summary Stats.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS ============================
% [MODIFIED 1] Replaced the single data_file with your full list of file paths
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX016/Result_SpikeNormGlobalRef_5uA_5ms_Xia_Exp1_Seq_Full_4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim7.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX014/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim5.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim7.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Seq_Sim8.mat';
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
    % NOTE: Add the rest of your files from your Grand Average script here!
};

% Threshold Definition
Target_Thresh = 0.5;   % Target normalized amplitude (e.g., 50% of the 5uA response)

% Plot Settings
dot_size = 30;
dot_color = [0.2 0.4 0.8]; 
dot_alpha = 0.4;           

%% =================== 2. MASTER LOOP INITIALIZATION ====================
% [MODIFIED 2] Master storage arrays and counters moved outside the loop
Thresh_Sim_Pool = [];
Thresh_Seq_Pool = [];
total_sets_processed = 0; % Tracks how many stimulation pairs we looked at

fprintf('Starting Grand Average Threshold Extraction (Target = %.2f)...\n', Target_Thresh);

% [MODIFIED 3] The Master Loop through all datasets
for f_idx = 1:length(file_paths)
    
    data_file = file_paths{f_idx};
    if ~exist(data_file, 'file')
        fprintf('Skipping missing file: %s\n', data_file);
        continue; 
    end

    D = load(data_file);
    if ~isfield(D, 'ResultNorm'), continue; end
    
    R = D.ResultNorm;
    Amps_Raw = R.Amps;
    Norm_Sim_Raw = R.Norm_Sim; % Dimensions: [nCh, nAmps, nSets]
    Norm_Seq_Raw = R.Norm_Seq;

    [nCh_Total, nAmps_Raw, nSets] = size(Norm_Sim_Raw);
    
    % Add to our running total of stimulation pairs
    total_sets_processed = total_sets_processed + nSets;

%% =================== 3. ANCHOR TO ZERO ====================
    % [MODIFIED] Force Amps_Raw to be a horizontal row vector (1D)
    % This prevents the 'horzcat' error if a dataset saved Amps as a column.
    Amps_Raw = Amps_Raw(:)'; 
    
    if Amps_Raw(1) ~= 0
        % Add 0 to the beginning of the Amps vector
        Amps = [0, Amps_Raw];
        
        % Pad the data matrices with zeros at the first amplitude index
        pad_zeros = zeros(nCh_Total, 1, nSets);
        Norm_Sim = cat(2, pad_zeros, Norm_Sim_Raw);
        Norm_Seq = cat(2, pad_zeros, Norm_Seq_Raw);
    else
        Amps = Amps_Raw;
        Norm_Sim = Norm_Sim_Raw;
        Norm_Seq = Norm_Seq_Raw;
    end
    %% =================== 4. CALCULATE THRESHOLDS ====================
    for ss = 1:nSets
        for ch = 1:nCh_Total
            
            curve_sim = Norm_Sim(ch, :, ss);
            curve_seq = Norm_Seq(ch, :, ss);
            
            if all(isnan(curve_sim)) || all(isnan(curve_seq)), continue; end
            
            t_sim = get_interpolated_threshold(Amps, curve_sim, Target_Thresh);
            t_seq = get_interpolated_threshold(Amps, curve_seq, Target_Thresh);
            
            % [MODIFIED 4] Reverted to STRICT EXCLUSION
            % If EITHER condition fails to reach 0.5, we throw it out.
            if isnan(t_sim) || isnan(t_seq)
                continue;
            end
            
            Thresh_Sim_Pool = [Thresh_Sim_Pool; t_sim];
            Thresh_Seq_Pool = [Thresh_Seq_Pool; t_seq];
        end
    end
end

%% =================== 5. SUMMARY STATS PRINTOUT ====================
% [MODIFIED 5] Calculates all the numbers you need for your paper's text
nPairs = length(Thresh_Sim_Pool);
if nPairs == 0
    fprintf('No channels reached the threshold. Try lowering the Target_Thresh.\n');
    return;
end

p_val = signrank(Thresh_Sim_Pool, Thresh_Seq_Pool);

fprintf('\n================== FINAL SUMMARY ==================\n');
fprintf('Target Threshold Level : %.2f a.u.\n', Target_Thresh);
fprintf('Total Datasets Loaded  : %d\n', length(file_paths));
fprintf('Total Stim Pairs (Sets): %d\n', total_sets_processed);
fprintf('Channels Plotted (N)   : %d (Strict Exclusion)\n', nPairs);
fprintf('---------------------------------------------------\n');
fprintf('Median Sim Threshold   : %.3f uA\n', median(Thresh_Sim_Pool));
fprintf('Median Seq Threshold   : %.3f uA\n', median(Thresh_Seq_Pool));
fprintf('Mean Sim Threshold     : %.3f uA\n', mean(Thresh_Sim_Pool));
fprintf('Mean Seq Threshold     : %.3f uA\n', mean(Thresh_Seq_Pool));
fprintf('Wilcoxon Paired p-value: %.2e\n', p_val);
fprintf('===================================================\n\n');

%% =================== 6. UNITY SCATTER PLOT ====================
figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w', 'PaperPositionMode', 'auto'); 
hold on;

max_val = ceil(max([Thresh_Sim_Pool; Thresh_Seq_Pool]));
if max_val < 5, max_val = 5; end 

plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

scatter(Thresh_Sim_Pool, Thresh_Seq_Pool, dot_size, dot_color, 'filled', ...
    'MarkerFaceAlpha', dot_alpha, 'MarkerEdgeColor', 'none');

if p_val < 0.001
    title_str = sprintf('Activation Threshold (p < 0.001***)');
elseif p_val < 0.01
    title_str = sprintf('Activation Threshold (p < 0.01**)');
elseif p_val < 0.05
    title_str = sprintf('Activation Threshold (p < 0.05*)');
else
    title_str = sprintf('Activation Threshold (n.s.)');
end
title(title_str, 'FontSize', 9, 'FontName', 'Arial', 'FontWeight', 'normal');

axis square;
box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.0);
xlabel('Simultaneous Threshold (\muA)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Sequential Threshold (\muA)', 'FontSize', 9,  'FontName', 'Arial');

xlim([0 max_val]);
ylim([0 max_val]);
xticks(0:2:max_val);
yticks(0:2:max_val);

hold off;

%% ==================== HELPER FUNCTIONS =========================
function thresh = get_interpolated_threshold(amps, curve, target)
    amps = squeeze(amps(:))';
    curve = squeeze(curve(:))';
    
    valid_idx = ~isnan(curve);
    valid_amps = amps(valid_idx);
    valid_curve = curve(valid_idx);
    
    if isempty(valid_curve) || max(valid_curve) < target
        thresh = NaN;
        return;
    end
    
    idx_above = find(valid_curve >= target, 1, 'first');
    
    if valid_curve(idx_above) == target
        thresh = valid_amps(idx_above);
        return;
    end
    
    if idx_above == 1
        if valid_amps(1) == 0
            thresh = 0; 
        else
            x0 = 0;
            y0 = 0;
            x1 = valid_amps(1);
            y1 = valid_curve(1);
            thresh = x0 + (target - y0) * (x1 - x0) / (y1 - y0);
        end
        return;
    end
    
    x0 = valid_amps(idx_above - 1);
    x1 = valid_amps(idx_above);
    y0 = valid_curve(idx_above - 1);
    y1 = valid_curve(idx_above);
    
    thresh = x0 + (target - y0) * (x1 - x0) / (y1 - y0);
end