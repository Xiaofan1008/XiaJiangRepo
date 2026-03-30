%% ============================================================
%   GRAND AVERAGE: PEAK LATENCY TUNING CURVE (Line Graph)
%   - Logic: 
%       1. Loads all single-dataset PeakLatency result files.
%       2. Pools median channel peak times for each amplitude.
%       3. Calculates Mean and SEM across the population.
%       4. Plots an I/O Line Graph (Sim vs Seq) with Error Bars.
%       5. Adds statistical significance stars at each amplitude step.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =================
% List Result_PeakLatency_*.mat files (Same as your histogram script)
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/Result_PeakLatency_Separated_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX014/Result_PeakLatency_Separated_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX013/Result_PeakLatency_Separated_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX012/Result_PeakLatency_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX012/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX012/Result_PeakLatency_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX011/Result_PeakLatency_Separated_Xia_Exp1_Sim9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX010/Result_PeakLatency_Separated_Xia_Exp1_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX009/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX009/Result_PeakLatency_Separated_Xia_Exp1_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX006/Result_PeakLatency_Separated_Xia_Exp1_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX005/Result_PeakLatency_Separated_Xia_Exp1_Sim.mat';
};

% Save Settings
save_figure = false;
save_dir    = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure3/Peak_FR_Time_Amp';
fig_name    = 'PeakLatency_Amplitude_TuningCurve';

% Cutoffs
max_cutoff_sim = 14; 
max_cutoff_seq = 16;

%% ================= 2. POOL DATA ACROSS DATASETS =================
fprintf('Pooling data from %d datasets...\n', length(file_paths));
Pooled = struct();
All_Amps = [];

for i = 1:length(file_paths)
    if ~exist(file_paths{i}, 'file'), continue; end
    D = load(file_paths{i});
    if ~isfield(D, 'PeakData'), continue; end
    
    P = D.PeakData;
    for s = 1:length(P.Set)
        if isempty(P.Set(s).Amp), continue; end
        for a = 1:length(P.Set(s).Amp)
            ThisAmpData = P.Set(s).Amp(a);
            val = ThisAmpData.Val;
            if val == 0, continue; end % Skip 0 uA
            
            fName = sprintf('A_%.1f', val);
            fName = strrep(fName, '.', 'p');
            
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = val;
                Pooled.(fName).Sim_Ch = []; % Store Channel Medians
                Pooled.(fName).Seq_Ch = [];
                All_Amps = [All_Amps, val]; 
            end
            
            raw_sim = ThisAmpData.Sim(:);
            raw_seq = ThisAmpData.Seq(:);
            
            raw_sim = raw_sim(raw_sim <= max_cutoff_sim);
            raw_seq = raw_seq(raw_seq <= max_cutoff_seq);
            
            % Store median latency for each channel to prevent skewing
            if ~isempty(raw_sim)
                Pooled.(fName).Sim_Ch = [Pooled.(fName).Sim_Ch; median(raw_sim)];
            end
            if ~isempty(raw_seq)
                Pooled.(fName).Seq_Ch = [Pooled.(fName).Seq_Ch; median(raw_seq)];
            end
        end
    end
end
All_Amps = unique(All_Amps);
All_Amps = sort(All_Amps);

%% ================= 3. CALCULATE MEAN, SEM, AND STATS =================
num_amps = length(All_Amps);
Mean_Sim = zeros(1, num_amps);
Sem_Sim  = zeros(1, num_amps);
Mean_Seq = zeros(1, num_amps);
Sem_Seq  = zeros(1, num_amps);
P_Vals   = zeros(1, num_amps);

fprintf('\nCalculating Mean, SEM, and RankSum p-values per Amplitude...\n');
for i = 1:num_amps
    curr_amp = All_Amps(i);
    fName = sprintf('A_%.1f', curr_amp); fName = strrep(fName, '.', 'p');
    
    sim_ch = Pooled.(fName).Sim_Ch;
    seq_ch = Pooled.(fName).Seq_Ch;
    
    % Compute Mean and SEM
    Mean_Sim(i) = mean(sim_ch); 
    Sem_Sim(i)  = std(sim_ch) / sqrt(length(sim_ch));
    
    Mean_Seq(i) = mean(seq_ch); 
    Sem_Seq(i)  = std(seq_ch) / sqrt(length(seq_ch));
    
    % Compute Stats
    if ~isempty(sim_ch) && ~isempty(seq_ch)
        P_Vals(i) = ranksum(sim_ch, seq_ch);
    else
        P_Vals(i) = NaN;
    end
    
    fprintf('Amp %2.1f uA: Sim=%.2fms (n=%d), Seq=%.2fms (n=%d) | p=%.5f\n', ...
        curr_amp, Mean_Sim(i), length(sim_ch), Mean_Seq(i), length(seq_ch), P_Vals(i));
end

%% ================= 4. PLOT LINE GRAPH WITH ERROR BARS =================
figure('Units', 'centimeters', 'Position', [5, 5, 8.8, 8.8], 'Color', 'w');
hold on;

% % 1. Plot Error Bars & Lines
% % Simultaneous: Solid Black Line, Filled Black Circles
% e1 = errorbar(All_Amps, Mean_Sim, Sem_Sim, '-ko', ...
%     'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'CapSize', 3);
% 
% % Sequential: Dashed Black Line, Open Squares (White face, black edge)
% e2 = errorbar(All_Amps, Mean_Seq, Sem_Seq, '--ks', ...
%     'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'CapSize', 3);

% 1. Create Shaded Error Areas (Plotted first so they sit in the background)
x_fill = [All_Amps, fliplr(All_Amps)];

% Simultaneous Shading (Darker Gray, transparent)
y_sim_fill = [Mean_Sim + Sem_Sim, fliplr(Mean_Sim - Sem_Sim)];
fill(x_fill, y_sim_fill, [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');

% Sequential Shading (Lighter Gray, transparent)
y_seq_fill = [Mean_Seq + Sem_Seq, fliplr(Mean_Seq - Sem_Seq)];
fill(x_fill, y_seq_fill, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');

% 2. Plot Main Lines on Top
% Simultaneous: Solid Black Line, Filled Black Circles
e1 = plot(All_Amps, Mean_Sim, '-ko', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'k');

% Sequential: Dashed Black Line, Open Squares
e2 = plot(All_Amps, Mean_Seq, '--ks', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'w');

% 2. Add Statistical Stars above the points
for i = 1:num_amps
    p = P_Vals(i);
    if isnan(p) || p >= 0.05, continue; end
    
    if p < 0.001, txt = '***';
    elseif p < 0.01, txt = '**';
    elseif p < 0.05, txt = '*';
    end
    
    % Find the highest point at this amplitude to float the star above the error bar
    highest_y = max(Mean_Sim(i) + Sem_Sim(i), Mean_Seq(i) + Sem_Seq(i));
    
    text(All_Amps(i), highest_y + 0.4, txt, ...
        'FontSize', 12, 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontName', 'Arial');
end

% 3. IEEE Formatting
xlabel('Stimulation Amplitude (\muA)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Peak Latency (ms)', 'FontSize', 9, 'FontName', 'Arial');

lgd = legend([e1, e2], {'Simultaneous', 'Sequential'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');

box off;
set(gca, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', 'LineWidth', 1.0);

% Clean Axes
xlim([0,10]);
xticks(All_Amps); % Force ticks exactly at your tested amplitudes

% Dynamically adjust Y limits to fit the data and stars
max_y_limit = max(max(Mean_Sim + Sem_Sim), max(Mean_Seq + Sem_Seq)) + 1.5;
min_y_limit = min(min(Mean_Sim - Sem_Sim), min(Mean_Seq - Sem_Seq)) - 0.5;
if min_y_limit < 0, min_y_limit = 0; end
ylim([min_y_limit, max_y_limit]);
ylim([0,16]);

axis square;
hold off;

%% ================= 5. SAVE FIGURE =================
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    full_path = fullfile(save_dir, fig_name);
    exportgraphics(gcf, [full_path '.tiff'], 'ContentType', 'vector');
    fprintf('>>> Figure saved to: %s\n', full_path);
end