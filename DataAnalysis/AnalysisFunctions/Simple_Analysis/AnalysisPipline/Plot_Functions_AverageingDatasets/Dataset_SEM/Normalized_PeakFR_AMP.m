%% ============================================================
%   GRAND AVERAGE: PEAK FIRING RATE (Cross-Dataset)
%   - Loads multiple 'Result_PeakFiringRate_... .mat' files
%   - Aligns data by Amplitude
%   - Computes Grand Mean +/- SEM (N = Number of Datasets)
%   - PLOT: Shaded Error Bars + Mean Line
% ============================================================
clear;

%% ================= 1. DEFINE FILES TO LOAD =================
file_list = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set1_PeakFiringRate_FixedPop_5ms_Xia_Exp1_Sim1_251125_112055.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set4_PeakFiringRate_FixedPop_5ms_Xia_Exp1_Sim4_251125_152849.mat';
};

%% ================= 2. ALIGN DATA BY AMPLITUDE =================
All_Amps_Found = [];
Data_Holder = {}; 

for f = 1:length(file_list)
    if ~exist(file_list{f}, 'file'), error('File not found: %s', file_list{f}); end
    loaded = load(file_list{f});
    
    if isfield(loaded, 'ResultFR'), D = loaded.ResultFR;
    else, error('File %d does not contain ResultFR structure', f); end
    
    Data_Holder{f} = D;
    if isfield(D.Metadata, 'Amps'), amps = D.Metadata.Amps; else, amps = D.Amps; end
    All_Amps_Found = [All_Amps_Found; amps];
end

Master_Amps = unique(All_Amps_Found);
nAmps = length(Master_Amps);
nFiles = length(file_list);

Grand_Sim = nan(nFiles, nAmps);
Grand_Seq = nan(nFiles, nAmps);

%% ================= 3. EXTRACT AND ALIGN =================
for f = 1:length(file_list)
    D = Data_Holder{f};
    these_amps = D.Metadata.Amps;
    
    % Sim: Mean across Channels
    val_sim = mean(D.PeakFR.Sim_Z, 1, 'omitnan'); 
    
    % Seq: Mean across Channels (Dim 1) THEN Mean across Sets (Dim 3)
    step1 = mean(D.PeakFR.Seq_Z, 1, 'omitnan');
    step2 = mean(step1, 3, 'omitnan');
    val_seq = squeeze(step2)'; 
    
    % Align
    [~, loc_in_master] = ismember(these_amps, Master_Amps);
    Grand_Sim(f, loc_in_master) = val_sim;
    Grand_Seq(f, loc_in_master) = val_seq;
end

%% ================= 4. COMPUTE GRAND STATS =================
Grand_Mean_Sim = mean(Grand_Sim, 1, 'omitnan');
Grand_SEM_Sim  = std(Grand_Sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Grand_Sim), 1));

Grand_Mean_Seq = mean(Grand_Seq, 1, 'omitnan');
Grand_SEM_Seq  = std(Grand_Seq, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Grand_Seq), 1));

%% ================= 5. PLOT GRAND AVERAGE (SHADED) =================
figure('Color','w', 'Position',[500 500 700 500]); hold on;

% Define Colors
col_sim = [0 0.3 0.8];       % Blue
col_seq = [0.85 0.33 0.10];  % Orange

% --- Plot Shaded Areas using Helper Function ---
plot_shaded_error(Master_Amps, Grand_Mean_Sim, Grand_SEM_Sim, col_sim);
plot_shaded_error(Master_Amps, Grand_Mean_Seq, Grand_SEM_Seq, col_seq);

% --- Plot Mean Lines ---
plot(Master_Amps, Grand_Mean_Sim, '-o', 'Color', col_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');

plot(Master_Amps, Grand_Mean_Seq, '-s', 'Color', col_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

xlabel('Amplitude (µA)', 'FontWeight', 'bold');
ylabel('Normalized Peak Firing Rate (sps)', 'FontWeight', 'bold');
title(['Average Peak Firing Rate (N=' num2str(nFiles) ' Datasets)'], 'FontWeight', 'bold');
legend('Location', 'best', 'Box', 'off'); 
box off;

%% ================= 6. PRINT TABLE FOR EXCEL =================
fprintf('\n\n==========================================================\n');
fprintf('           GRAND AVERAGE RESULTS          \n');
fprintf('Amp(uA)\tSim_Mean\tSim_SEM\tSeq_Mean\tSeq_SEM\n');
for i = 1:nAmps
    fprintf('%.0f\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        Master_Amps(i), ...
        Grand_Mean_Sim(i), Grand_SEM_Sim(i), ...
        Grand_Mean_Seq(i), Grand_SEM_Seq(i));
end
fprintf('==========================================================\n');

%% ==================== HELPER FUNCTION ====================
function plot_shaded_error(x, y, se, col)
    % Ensures inputs are row vectors to avoid "fill" errors
    if numel(x) < 2, return; end
    
    x=x(:)'; y=y(:)'; se=se(:)'; 
    
    % Filter NaNs
    valid = ~isnan(y) & ~isnan(se); 
    x=x(valid); y=y(valid); se=se(valid);
    
    if isempty(x), return; end
    
    upper = y + se; 
    lower = y - se;
    
    % Plot shaded area
    fill([x fliplr(x)], [upper fliplr(lower)], col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end