%% ============================================================
%   POPULATION POTENCY ANALYSIS: 32-Channel Dual-Shank Cluster
%   - Logic: 
%       1. Standard Normalization: Scaled to each Rat's Max.
%       2. Ratio Normalization: Sequential / Simultaneous.
%       3. Delta Change: Sequential - Simultaneous.
%   - Style: IEEE TBME Publication Style (B&W)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX005_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX009_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX009_Xia_Exp1_Seq5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq2_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX012_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim4.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim6.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_figs = true; 
save_dir  = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Active_Percentage';

%% =================== 1. AGGREGATE POOLED DATA =================
fprintf('Pooling Potency data from %d datasets...\n', length(file_paths));
Pooled = struct(); unique_amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), continue; end
    D = load(file_paths{f});
    
    % --- NEW: Find Max Response for this specific file/rat for Standard Normalization ---
    max_resp_rat = 0;
    for ss = 1:length(D.PotencyResults.Set)
        for ai = 1:length(D.PotencyResults.Set(ss).Amp)
            max_resp_rat = max([max_resp_rat, D.PotencyResults.Set(ss).Amp(ai).TotalPerc_Sim, D.PotencyResults.Set(ss).Amp(ai).TotalPerc_Seq]);
        end
    end
    if max_resp_rat == 0, max_resp_rat = 1; end % Avoid divide by zero

    for ss = 1:length(D.PotencyResults.Set)
        for ai = 1:length(D.PotencyResults.Set(ss).Amp)
            Data = D.PotencyResults.Set(ss).Amp(ai);
            if isempty(Data.Val) || Data.Val > 10, continue; end 
            
            fName = sprintf('A_%.1f', Data.Val); fName = strrep(fName, '.', 'p');
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = Data.Val; unique_amps = [unique_amps, Data.Val];
                Pooled.(fName).Sim_Raw = []; Pooled.(fName).Seq_Raw = [];
                Pooled.(fName).Sim_Std = []; Pooled.(fName).Seq_Std = [];
                Pooled.(fName).Delta   = []; Pooled.(fName).Ratio   = [];
            end
            
            % 1. Raw Data
            s_raw = Data.TotalPerc_Sim; q_raw = Data.TotalPerc_Seq;
            Pooled.(fName).Sim_Raw = [Pooled.(fName).Sim_Raw; s_raw];
            Pooled.(fName).Seq_Raw = [Pooled.(fName).Seq_Raw; q_raw];
            
            % 2. Standard Normalization (Scaled to Rat Max)
            Pooled.(fName).Sim_Std = [Pooled.(fName).Sim_Std; s_raw / max_resp_rat];
            Pooled.(fName).Seq_Std = [Pooled.(fName).Seq_Std; q_raw / max_resp_rat];
            
            % 3. Delta Change
            Pooled.(fName).Delta = [Pooled.(fName).Delta; (q_raw - s_raw)];
            
            % 4. Ratio (Seq / Sim) - Safety check for 0uA
            if s_raw > 0, r_val = q_raw / s_raw; else, r_val = 1; end
            Pooled.(fName).Ratio = [Pooled.(fName).Ratio; r_val];
        end
    end
end
unique_amps = sort(unique(unique_amps));

%% ================= 4. SUMMARY DIAGNOSTIC PLOTS (IEEE Single-Column) =================
sum_amps = []; n_counts = [];
mSimRaw = []; mSeqRaw = []; mSimStd = []; mSeqStd = []; mDelta = []; mRatio = [];
p_raw = []; p_std = []; p_delta = [];

% Setup arrays to collect raw data for ANOVA
anova_deltas = [];
anova_groups = [];

for i = 1:length(unique_amps)
    fName = sprintf('A_%.1f', unique_amps(i)); fName = strrep(fName, '.', 'p');
    P = Pooled.(fName);
    if length(P.Sim_Raw) < 5, continue; end
    sum_amps = [sum_amps; P.Val];
    n_counts = [n_counts; length(P.Sim_Raw)];
    
    % Collect raw data for ANOVA
    anova_deltas = [anova_deltas; P.Delta];
    anova_groups = [anova_groups; repmat(P.Val, length(P.Delta), 1)];
    
    % Statistics Calculation (Mean & SEM) - Within-Subject
    mSimRaw(end+1,:) = [mean(P.Sim_Raw), std(P.Sim_Raw)/sqrt(length(P.Sim_Raw))];
    mSeqRaw(end+1,:) = [mean(P.Seq_Raw), std(P.Seq_Raw)/sqrt(length(P.Seq_Raw))];
    mSimStd(end+1,:) = [mean(P.Sim_Std), std(P.Sim_Std)/sqrt(length(P.Sim_Std))];
    mSeqStd(end+1,:) = [mean(P.Seq_Std), std(P.Seq_Std)/sqrt(length(P.Seq_Std))];
    mDelta(end+1,:)  = [mean(P.Delta),   std(P.Delta)/sqrt(length(P.Delta))];
    mRatio(end+1,:)  = [mean(P.Ratio),   std(P.Ratio)/sqrt(length(P.Ratio))];
    
    [~, p_raw(end+1)]   = ttest2(P.Sim_Raw, P.Seq_Raw);
    [~, p_std(end+1)]   = ttest2(P.Sim_Std, P.Seq_Std);
    [~, p_delta(end+1)] = ttest(P.Delta, 0); 
end

% Perform Global One-Way ANOVA on Delta Change
[p_anova, tbl_anova, stats_anova] = anova1(anova_deltas, anova_groups, 'off');
F_stat = tbl_anova{2, 5};
df1 = tbl_anova{2, 3};
df2 = tbl_anova{3, 3};

% Determine ANOVA Star
if p_anova < 0.001, anova_star = '***';
elseif p_anova < 0.01, anova_star = '**';
elseif p_anova < 0.05, anova_star = '*';
else, anova_star = 'n.s.'; end

% --- Common IEEE Formatting Settings ---
fig_Width = 8.8;  % IEEE single column width in cm
fig_Height = 8.8; % Square aspect ratio
fSize = 9;        % IEEE standard font size
fName_Font = 'Arial';

% --- PANEL 1: Standard Normalization ---
figure('Color','w', 'Units', 'centimeters', 'Position', [2, 5, fig_Width, fig_Height], 'Name', 'Panel_A_Standard');
hold on;
errorbar(sum_amps, mSimStd(:,1), mSimStd(:,2), '--ok', 'LineWidth', 1, 'MarkerFaceColor','w', 'DisplayName', 'Simultaneous');
errorbar(sum_amps, mSeqStd(:,1), mSeqStd(:,2), '-sk', 'LineWidth', 1, 'MarkerFaceColor','k', 'DisplayName', 'Sequential');
ylabel('Normalized Active Percentage', 'FontName', fName_Font, 'FontSize', fSize); 
xlabel('Amplitude (µA)', 'FontName', fName_Font, 'FontSize', fSize);
% title('Standard (Max-Scaled)', 'FontName', fName_Font, 'FontSize', fSize); 
set(gca, 'FontSize', fSize, 'FontName', fName_Font, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1);
grid off; axis square;
add_sig_stars(sum_amps, mSimStd(:,1), mSeqStd(:,1), mSeqStd(:,2), p_std);
legend('Location', 'northwest', 'Box', 'off', 'FontSize', fSize);

% --- PANEL 2: Delta Change (Naka-Rushton Fit) ---
figure('Color','w', 'Units', 'centimeters', 'Position', [12, 5, fig_Width, fig_Height], 'Name', 'Panel_B_Delta');
hold on;
% Calculate Naka-Rushton Fit (Bio-Inspired Sigmoid)
% Response = (alpha * x^beta) / (x^beta + gamma^beta)
ft = fittype('a * x^b / (x^b + c^b)', 'independent', 'x');
opts = fitoptions(ft);
opts.StartPoint = [max(mDelta(:,1)), 2, 5]; % [alpha, beta, gamma]
opts.Lower = [0, 0.1, 0.1];
opts.Upper = [100, 15, 100];
[fitobj, gof] = fit(sum_amps, mDelta(:,1), ft, opts);
xq = linspace(0, 10, 100);
yq = feval(fitobj, xq);
plot(xq, yq, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5); % The Trend Line
% Plot Mean Markers (Black circles on top)
errorbar(sum_amps, mDelta(:,1), mDelta(:,2), 'ok', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'LineStyle', 'none');

% % Add ANOVA Spanning Bracket instead of point-by-point stars
% y_bracket = max(mDelta(:,1) + mDelta(:,2)) + 2.5; % Set bracket above highest error bar
% plot([min(sum_amps), max(sum_amps)], [y_bracket, y_bracket], '-k', 'LineWidth', 1); % Horizontal span
% plot([min(sum_amps), min(sum_amps)], [y_bracket-0.7, y_bracket], '-k', 'LineWidth', 1); % Left tick
% plot([max(sum_amps), max(sum_amps)], [y_bracket-0.7, y_bracket], '-k', 'LineWidth', 1); % Right tick
% text(mean([min(sum_amps), max(sum_amps)]), y_bracket + 1, anova_star, 'FontSize', 12, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', fName_Font);

% Add ANOVA Spanning Bracket instead of point-by-point stars
y_bracket = max(mDelta(:,1) + mDelta(:,2)) + 2.5; % Set bracket above highest error bar
% Define bracket edges with a slight indent (0.5 uA) so it doesn't touch the axes
x_left = min(sum_amps) + 0.5;  
x_right = max(sum_amps) - 0.5; 

plot([x_left, x_right], [y_bracket, y_bracket], '-k', 'LineWidth', 1); % Horizontal span
plot([x_left, x_left], [y_bracket-0.7, y_bracket], '-k', 'LineWidth', 1); % Left tick
plot([x_right, x_right], [y_bracket-0.7, y_bracket], '-k', 'LineWidth', 1); % Right tick
text(mean([x_left, x_right]), y_bracket + 1, anova_star, 'FontSize', 12, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', fName_Font);


ylim([0, y_bracket + 3.5]); % Dynamically adjust limits so bracket fits nicely

ylabel('Increase in Active Percentage', 'FontName', fName_Font, 'FontSize', fSize); 
xlabel('Amplitude (µA)', 'FontName', fName_Font, 'FontSize', fSize);
% title('Delta Change (Naka-Rushton)', 'FontName', fName_Font, 'FontSize', fSize); 
set(gca, 'FontSize', fSize, 'FontName', fName_Font, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1);
grid off; axis square;

% --- PANEL 3: Ratio (Shaded SEM) ---
figure('Color','w', 'Units', 'centimeters', 'Position', [22, 5, fig_Width, fig_Height], 'Name', 'Panel_C_Ratio');
hold on;
plot([0 10], [1 1], 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); % Baseline at 1.0
x_vec = sum_amps(:)'; m_vec = mRatio(:,1)'; s_vec = mRatio(:,2)';
fill([x_vec, fliplr(x_vec)], [m_vec+s_vec, fliplr(m_vec-s_vec)], 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x_vec, m_vec, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); 
ylabel('Sequential/Simultaneous', 'FontName', fName_Font, 'FontSize', fSize); 
xlabel('Amplitude (µA)', 'FontName', fName_Font, 'FontSize', fSize);
% title('Ratio (Efficiency)', 'FontName', fName_Font, 'FontSize', fSize); 
ylim([-2 5.0]); 
set(gca, 'FontSize', fSize, 'FontName', fName_Font, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1);
grid off; axis square;

%% ================= 5. MASTER STATISTICAL TABLE =================
fprintf('\n============================================================================\n');
fprintf('       GLOBAL STATISTICS & MODEL FIT RESULTS                                \n');
fprintf('============================================================================\n');
fprintf('1. ANOVA (Effect of Amplitude on Delta Change):\n');
fprintf('   F(%d, %d) = %.2f, p = %.4e\n', df1, df2, F_stat, p_anova);
if p_anova < 0.05
    fprintf('   -> SIGNIFICANT main effect of amplitude on sequential advantage.\n');
else
    fprintf('   -> NOT significant.\n');
end

fprintf('\n2. NAKA-RUSHTON MODEL FIT (Delta Change):\n');
fprintf('   Equation: (alpha * x^beta) / (x^beta + gamma^beta)\n');
fprintf('   Goodness of Fit: R^2 = %.3f, RMSE = %.3f\n', gof.rsquare, gof.rmse);
fprintf('   Coefficients:\n');
fprintf('      alpha (max gain) = %.3f\n', fitobj.a);
fprintf('      beta  (slope)    = %.3f\n', fitobj.b);
fprintf('      gamma (C50)      = %.3f\n', fitobj.c);
fprintf('\n----------------------------------------------------------------------------\n');
fprintf('       DIAGNOSTIC SUMMARY: COMPARING NORMALIZATION METHODS                  \n');
fprintf('----------------------------------------------------------------------------\n');
fprintf('%-6s | %-10s | %-10s | %-10s | %-10s\n', 'Amp', 'Raw Delta', 'Ratio', 'Std Sim', 'Std Seq');
fprintf('----------------------------------------------------------------------------\n');
for i = 1:length(sum_amps)
    fprintf('%3.1fuA  | +%4.1f%%     | %4.2fx      | %4.2f       | %4.2f\n', ...
        sum_amps(i), mDelta(i,1), mRatio(i,1), mSimStd(i,1), mSeqStd(i,1));
end
fprintf('============================================================================\n');

%% ================= 6. SAVE FIGURES (TIFF) =================
if save_figs
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    figHandles = findall(0, 'Type', 'figure');
    for i = 1:length(figHandles)
        f = figHandles(i);
        fName = get(f, 'Name'); if isempty(fName), fName = ['Fig_' num2str(i)]; end
        fName = strrep(fName, ' ', '_'); 
        exportgraphics(f, fullfile(save_dir, [fName '.tiff']), 'Resolution', 600);
        fprintf('Saved: %s.tiff\n', fName);
    end
end

function add_sig_stars(x, y1, y2, y2_err, pvals)
    for i = 1:length(x)
        p = pvals(i); if isnan(p), continue; end
        txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        if ~isempty(txt)
            y_max = max(y1(i), y2(i) + y2_err(i));
            text(x(i), y_max * 1.05, txt, 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
end