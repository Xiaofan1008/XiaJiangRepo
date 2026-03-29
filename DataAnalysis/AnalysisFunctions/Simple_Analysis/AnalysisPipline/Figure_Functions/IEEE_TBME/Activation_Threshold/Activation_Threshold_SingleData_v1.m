%% ============================================================
%   ACTIVATION THRESHOLD EXTRACTOR (Single Dataset)
%   - Metric: Amplitude (uA) required to reach Target Normalized Response
%   - Logic: 
%       1. Anchor curves at 0 uA = 0 response.
%       2. Interpolate linearly to find first crossing of Target_Thresh.
%       3. Require BOTH Sim and Seq to cross target (Strict Exclusion).
%       4. Plot paired Unity Scatter (Sim vs. Seq).
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS ============================
% File path for the single dataset you want to test
data_file = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX006/Result_SpikeNormGlobalRef_5uA_Zeroed_5ms_Xia_Exp1_Sim.mat';

% Threshold Definition
Target_Thresh = 0.5;   % Target normalized amplitude (e.g., 50% of the 5uA response)

% Plot Settings
dot_size = 30;
dot_color = [0.2 0.4 0.8]; % A nice professional blue
dot_alpha = 0.4;           % Transparency to handle overlapping points

%% =================== 2. LOAD DATA ====================
fprintf('Loading data...\n');
if ~exist(data_file, 'file')
    error('File not found. Please check the path.');
end

D = load(data_file);
if ~isfield(D, 'ResultNorm')
    error('ResultNorm structure not found in the file.');
end

R = D.ResultNorm;
Amps_Raw = R.Amps;
Norm_Sim_Raw = R.Norm_Sim; % Dimensions: [nCh, nAmps, nSets]
Norm_Seq_Raw = R.Norm_Seq;

[nCh_Total, nAmps_Raw, nSets] = size(Norm_Sim_Raw);

%% =================== 3. ANCHOR TO ZERO ====================
% To handle "Hyper Responders" (channels that jump past 0.5 at the very first tested amp),
% we must ensure the curve starts at 0 uA = 0 a.u.
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
fprintf('Calculating thresholds (Target = %.2f a.u.)...\n', Target_Thresh);

% Storage arrays for the paired scatter plot
Thresh_Sim_Pool = [];
Thresh_Seq_Pool = [];

% Double Loop: Process Set by Set, then Channel by Channel
for ss = 1:nSets
    for ch = 1:nCh_Total
        
        % Extract the specific curves for this channel and set
        curve_sim = Norm_Sim(ch, :, ss);
        curve_seq = Norm_Seq(ch, :, ss);
        
        % Ignore channels that have no data (all NaNs) for this set
        if all(isnan(curve_sim)) || all(isnan(curve_seq)), continue; end
        
        % Run the interpolation helper function (defined at bottom of script)
        t_sim = get_interpolated_threshold(Amps, curve_sim, Target_Thresh);
        t_seq = get_interpolated_threshold(Amps, curve_seq, Target_Thresh);
        
        % STRICT EXCLUSION: Only keep if BOTH conditions reached the threshold
        if ~isnan(t_sim) && ~isnan(t_seq)
            Thresh_Sim_Pool = [Thresh_Sim_Pool; t_sim];
            Thresh_Seq_Pool = [Thresh_Seq_Pool; t_seq];
        end
    end
end

nPairs = length(Thresh_Sim_Pool);
fprintf('Found %d valid paired thresholds across %d Sets.\n', nPairs, nSets);

%% =================== 5. UNITY SCATTER PLOT ====================
if nPairs == 0
    fprintf('No channels reached the threshold. Try lowering the Target_Thresh.\n');
    return;
end

% Create IEEE formatted figure
figure('Units', 'centimeters', 'Position', [2, 2, 8.89, 8.89], 'Color', 'w', 'PaperPositionMode', 'auto'); 
hold on;

% Determine axis limits (give it a little buffer past the max threshold)
max_val = ceil(max([Thresh_Sim_Pool; Thresh_Seq_Pool]));
if max_val < 5, max_val = 5; end % Ensure at least 5uA axis

% 1. Plot the Unity Line (y = x)
% Anything BELOW this line means Sequential threshold < Simultaneous threshold
plot([0 max_val], [0 max_val], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% 2. Scatter the paired thresholds
scatter(Thresh_Sim_Pool, Thresh_Seq_Pool, dot_size, dot_color, 'filled', ...
    'MarkerFaceAlpha', dot_alpha, 'MarkerEdgeColor', 'none');

% 3. Statistics (Paired Wilcoxon Signed Rank)
p_val = signrank(Thresh_Sim_Pool, Thresh_Seq_Pool);
fprintf('\n=== STATS ===\n');
fprintf('Median Sim Threshold: %.2f uA\n', median(Thresh_Sim_Pool));
fprintf('Median Seq Threshold: %.2f uA\n', median(Thresh_Seq_Pool));
fprintf('Wilcoxon Signed Rank p-value: %.5f\n', p_val);

% Add significance star to plot if valid
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

% 4. Formatting to match your IEEE style
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
    % [MODIFIED 1] Ensure arrays are flat 1D vectors to prevent dimension mismatch
    amps = squeeze(amps(:))';
    curve = squeeze(curve(:))';

    % Remove NaNs from the curve
    valid_idx = ~isnan(curve);
    valid_amps = amps(valid_idx);
    valid_curve = curve(valid_idx);
    
    % If empty or the max response never reaches the target, return NaN (Exclude)
    if isempty(valid_curve) || max(valid_curve) < target
        thresh = NaN;
        return;
    end
    
    % Find the FIRST index where the curve crosses or equals the target
    idx_above = find(valid_curve >= target, 1, 'first');
    
    % Exact match case
    if valid_curve(idx_above) == target
        thresh = valid_amps(idx_above);
        return;
    end
    
    % [MODIFIED 2] BUG FIX: If the VERY FIRST valid point is already above the threshold
    if idx_above == 1
        if valid_amps(1) == 0
            thresh = 0; % It was somehow above the target exactly at 0 uA
        else
            % Interpolate assuming 0 response at 0 uA
            x0 = 0;
            y0 = 0;
            x1 = valid_amps(1);
            y1 = valid_curve(1);
            thresh = x0 + (target - y0) * (x1 - x0) / (y1 - y0);
        end
        return;
    end
    
    % Normal Interpolation case
    % Get the points immediately before and after the crossing
    x0 = valid_amps(idx_above - 1);
    x1 = valid_amps(idx_above);
    y0 = valid_curve(idx_above - 1);
    y1 = valid_curve(idx_above);
    
    % Linear interpolation formula: x = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
    thresh = x0 + (target - y0) * (x1 - x0) / (y1 - y0);
end