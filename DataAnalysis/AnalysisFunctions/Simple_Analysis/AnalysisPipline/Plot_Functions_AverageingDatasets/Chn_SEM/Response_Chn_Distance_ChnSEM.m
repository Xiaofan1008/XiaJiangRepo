%% ============================================================
%   GRAND AVERAGE: SPATIAL DIVERSITY (Final Version)
%   - Fig 1: Response Probability vs Distance (Bar Plot - CLEAN)
%   - Fig 2: Activation Radius (d50) vs Amplitude (NaN -> 0)
% ============================================================
clear;

%% ================= 1. SETTINGS & FILES =================
file_list = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set1_Spatial_d50_Stats_5ms_Xia_Exp1_Sim1_251125_112055.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set4_Spatial_d50_Stats_5ms_Xia_Exp1_Sim4_251125_152849.mat';
};

Master_Dists = 0:50:600; 
nMasterDists = length(Master_Dists);
col_sim = [0.27 0.45 0.70]; % Slate Blue
col_seq = [0.85 0.50 0.30]; % Burnt Orange

%% ================= 2. DATA PREPARATION =================
All_Amps_Found = [];
Data_Holder = {}; 

for f = 1:length(file_list)
    if ~exist(file_list{f}, 'file'), error('File not found: %s', file_list{f}); end
    loaded = load(file_list{f});
    Data_Holder{f} = loaded.ResultSpatial;
    All_Amps_Found = [All_Amps_Found; loaded.ResultSpatial.d50_Table(:, 1)];
end
Master_Amps = unique(All_Amps_Found);
nAmps = length(Master_Amps);
nFiles = length(file_list);

%% ================= 3. EXTRACT & ALIGN DATA =================
Grand_d50_Sim = nan(nFiles, nAmps);
Grand_d50_Seq = nan(nFiles, nAmps);
Grand_Prob_Sim = nan(nFiles, nMasterDists);
Grand_Prob_Seq = nan(nFiles, nMasterDists);

for f = 1:length(file_list)
    D = Data_Holder{f};
    
    % --- PART A: d50 Curves ---
    file_amps = D.d50_Table(:, 1);
    file_sim  = D.d50_Table(:, 2);
    file_seq  = D.d50_Table(:, 3); % Seq Mean
    
    [~, loc_amp] = ismember(file_amps, Master_Amps);
    Grand_d50_Sim(f, loc_amp) = file_sim;
    Grand_d50_Seq(f, loc_amp) = file_seq;
    
    % --- PART B: Bar Plots ---
    file_dists = D.BarData.Distances;
    raw_probs  = D.BarData.Probs;
    sim_trace = raw_probs(:, 1);
    seq_trace = mean(raw_probs(:, 2:end), 2, 'omitnan');
    
    for i = 1:nMasterDists
        d_target = Master_Dists(i);
        idx = find(abs(file_dists - d_target) < 10, 1);
        if ~isempty(idx)
            Grand_Prob_Sim(f, i) = sim_trace(idx);
            Grand_Prob_Seq(f, i) = seq_trace(idx);
        else
            Grand_Prob_Sim(f, i) = 0; 
            Grand_Prob_Seq(f, i) = 0; 
        end
    end
end

%% ================= 4. REPLACE NaN WITH 0 (For d50) =================
% Treat failed fits as zero spread
nans_sim = isnan(Grand_d50_Sim);
if any(nans_sim(:)), Grand_d50_Sim(nans_sim) = 0; end

nans_seq = isnan(Grand_d50_Seq);
if any(nans_seq(:)), Grand_d50_Seq(nans_seq) = 0; end

%% ================= 5. COMPUTE STATS =================
% d50 Stats
Mean_d50_Sim = mean(Grand_d50_Sim, 1, 'omitnan');
SEM_d50_Sim  = std(Grand_d50_Sim, 0, 1, 'omitnan') ./ sqrt(nFiles);
Mean_d50_Seq = mean(Grand_d50_Seq, 1, 'omitnan');
SEM_d50_Seq  = std(Grand_d50_Seq, 0, 1, 'omitnan') ./ sqrt(nFiles);

% Bar Plot Stats
Mean_Prob_Sim = mean(Grand_Prob_Sim, 1, 'omitnan');
SEM_Prob_Sim  = std(Grand_Prob_Sim, 0, 1, 'omitnan') ./ sqrt(nFiles);
Mean_Prob_Seq = mean(Grand_Prob_Seq, 1, 'omitnan');
SEM_Prob_Seq  = std(Grand_Prob_Seq, 0, 1, 'omitnan') ./ sqrt(nFiles);

%% ================= 6. FIGURE 1: BAR PLOT (NO ERROR BARS) =================
figure('Color','w', 'Position',[100 500 700 500]); hold on;

% Plot Grouped Bars
b = bar(Master_Dists, [Mean_Prob_Sim; Mean_Prob_Seq]', 'grouped');

% Styling
b(1).FaceColor = col_sim; b(1).EdgeColor = 'none'; b(1).FaceAlpha = 0.8; b(1).DisplayName = 'Simultaneous';
b(2).FaceColor = col_seq; b(2).EdgeColor = 'none'; b(2).FaceAlpha = 0.8; b(2).DisplayName = 'Sequential';

xlabel('Distance from Stimulation (\mum)', 'FontWeight','bold');
ylabel('Probability of Response (%)', 'FontWeight','bold');
title(['Average Spatial Diversity'], 'FontWeight','bold');
legend('Location','northeast','Box','off'); box off;

xticks(Master_Dists);       % Force every 50um tick to show
% xtickangle(45);             % Rotate labels so they don't overlap
xlim([min(Master_Dists)-25, max(Master_Dists)+25]); 

%% ================= 7. FIGURE 2: d50 CURVE (Radius vs Amp) =================
figure('Color','w', 'Position',[850 500 600 500]); hold on;

% Shaded Errors
plot_shaded_error(Master_Amps, Mean_d50_Sim, SEM_d50_Sim, col_sim);
plot_shaded_error(Master_Amps, Mean_d50_Seq, SEM_d50_Seq, col_seq);

% Mean Lines
plot(Master_Amps, Mean_d50_Sim, '-o', 'Color', col_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');
plot(Master_Amps, Mean_d50_Seq, '-s', 'Color', col_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

xlabel('Amplitude (\muA)', 'FontWeight','bold');
ylabel('Activation Radius (d_{50}) [\mum]', 'FontWeight','bold');
title('Spatial Spread vs Amplitude', 'FontWeight','bold');
legend('Location','best','Box','off'); box off;

%% ================= 8. PRINT SUMMARY TABLES =================
fprintf('\n\n==========================================================\n');
fprintf('           SPATIAL PROBABILITY (Bar Plot Data)          \n');
fprintf('Dist(um)\tSim_Mean\tSim_SEM\tSeq_Mean\tSeq_SEM\n');
for i = 1:nMasterDists
    fprintf('%.0f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
        Master_Dists(i), ...
        Mean_Prob_Sim(i), SEM_Prob_Sim(i), ...
        Mean_Prob_Seq(i), SEM_Prob_Seq(i));
end
fprintf('\n==========================================================\n');
fprintf('           ACTIVATION RADIUS (d50 Curve Data)          \n');
fprintf('Amp(uA)\tSim_Mean\tSim_SEM\tSeq_Mean\tSeq_SEM\n');
for i = 1:nAmps
    fprintf('%.0f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
        Master_Amps(i), ...
        Mean_d50_Sim(i), SEM_d50_Sim(i), ...
        Mean_d50_Seq(i), SEM_d50_Seq(i));
end
fprintf('==========================================================\n');

%% ================= HELPER FUNCTION =================
function plot_shaded_error(x, y, se, col)
    if numel(x)<2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid=~isnan(y)&~isnan(se); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end