%% ============================================================
%   GRAND AVERAGE: RESPONSE DURATION
%   - Fig 1: Duration vs Amplitude (Line Plot with Shaded SEM)
%   - Fig 2: Box Plot of Population Distribution at Target Amp
%   - Input: Loads 'Result_Set*_Duration_... .mat' files
% ============================================================
clear;

%% ================= 1. SETTINGS & FILES =================
target_amp_for_box = 10; % Amplitude to use for Figure 2 (Box Plot)

file_list = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set1_Duration_SpecificPopulation_5ms_Xia_Exp1_Sim1_251125_112055.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set4_Duration_SpecificPopulation_5ms_Xia_Exp1_Sim4_251125_152849.mat';
};

%% ================= 2. DATA PREPARATION =================
% A. Scan for Amplitudes (for Line Plot Alignment)
All_Amps_Found = [];
Data_Holder = {}; 

for f = 1:length(file_list)
    if ~exist(file_list{f}, 'file'), error('File not found: %s', file_list{f}); end
    loaded = load(file_list{f});
    if isfield(loaded, 'ResultDur'), D = loaded.ResultDur;
    else, error('File %d does not contain ResultDur', f); end
    
    Data_Holder{f} = D;
    All_Amps_Found = [All_Amps_Found; D.Data.Amps];
end
Master_Amps = unique(All_Amps_Found);
nAmps = length(Master_Amps);
nFiles = length(file_list);

% Storage for Line Plot (Means per File)
Grand_Sim = nan(nFiles, nAmps);
Grand_Seq = nan(nFiles, nAmps);

% Storage for Box Plot (Raw Channels from all files)
Box_Data_Sim = [];
Box_Data_Seq = [];

%% ================= 3. EXTRACT DATA =================
for f = 1:length(file_list)
    D = Data_Holder{f};
    these_amps = D.Data.Amps;
    
    % --- PART A: Line Plot Data (Subject Means) ---
    % Sim: Mean across Channels (Dim 1) -> [1 x Amps]
    mean_sim_ch = mean(D.Data.Dur_Sim, 1, 'omitnan');
    
    % Seq: Mean across Channels (Dim 1) AND Sets (Dim 3) -> [1 x Amps]
    % Dur_Seq is [Channels x Amps x Sets x PTD]
    mean_seq_sets = mean(D.Data.Dur_Seq, 3, 'omitnan'); 
    mean_seq_ch   = mean(mean_seq_sets, 1, 'omitnan');
    val_seq       = squeeze(mean_seq_ch)';
    
    % Align and Store
    [~, loc] = ismember(these_amps, Master_Amps);
    Grand_Sim(f, loc) = mean_sim_ch;
    Grand_Seq(f, loc) = val_seq;
    
    % --- PART B: Box Plot Data (Raw Channels at Target Amp) ---
    idx_amp = find(these_amps == target_amp_for_box);
    if ~isempty(idx_amp)
        % Sim Raw Data
        raw_sim = D.Data.Dur_Sim(:, idx_amp);
        raw_sim = raw_sim(~isnan(raw_sim)); % Remove non-responders
        Box_Data_Sim = [Box_Data_Sim; raw_sim];
        
        % Seq Raw Data (Pool all sets)
        raw_seq = D.Data.Dur_Seq(:, idx_amp, :, :);
        raw_seq = raw_seq(:); % Flatten to vector
        raw_seq = raw_seq(~isnan(raw_seq));
        Box_Data_Seq = [Box_Data_Seq; raw_seq];
    end
end

%% ================= 4. COMPUTE STATS (Line Plot) =================
Grand_Mean_Sim = mean(Grand_Sim, 1, 'omitnan');
Grand_SEM_Sim  = std(Grand_Sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Grand_Sim), 1));

Grand_Mean_Seq = mean(Grand_Seq, 1, 'omitnan');
Grand_SEM_Seq  = std(Grand_Seq, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Grand_Seq), 1));

%% ================= 5. FIGURE 1: LINE PLOT (Duration vs Amp) =================
figure('Color','w', 'Position',[100 500 600 500]); hold on;
col_sim = [0 0.3 0.8];
col_seq = [0.85 0.33 0.10];

% Shaded Errors
plot_shaded_error(Master_Amps, Grand_Mean_Sim, Grand_SEM_Sim, col_sim);
plot_shaded_error(Master_Amps, Grand_Mean_Seq, Grand_SEM_Seq, col_seq);

% Mean Lines
plot(Master_Amps, Grand_Mean_Sim, '-o', 'Color', col_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');
plot(Master_Amps, Grand_Mean_Seq, '-s', 'Color', col_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

xlabel('Amplitude (\muA)', 'FontWeight', 'bold');
ylabel('Net Active Duration (ms)', 'FontWeight', 'bold');
title(['Grand Average Duration (N=' num2str(nFiles) ')'], 'FontWeight', 'bold');
legend('Location', 'best', 'Box', 'off'); grid on; box off;

%% ================= 6. FIGURE 2: BOX PLOT (Distribution at Target Amp) =================
figure('Color','w', 'Position',[800 500 400 500]); hold on;

% Combine data for boxplot
% Create Group Labels
g1 = repmat({'Simultaneous'}, length(Box_Data_Sim), 1);
g2 = repmat({'Sequential'}, length(Box_Data_Seq), 1);

if ~isempty(Box_Data_Sim) || ~isempty(Box_Data_Seq)
    boxplot([Box_Data_Sim; Box_Data_Seq], [g1; g2], 'Width', 0.5, 'Symbol', '.');
    ylabel('Duration (ms)', 'FontWeight', 'bold');
    title(sprintf('Population Distribution at %.0f µA', target_amp_for_box), 'FontWeight', 'bold');
    box off;
    
    % Add individual points (jittered)
    sw = 0.5; % scatter width
    x1 = (rand(size(Box_Data_Sim))-0.5)*0.2 + 1;
    scatter(x1, Box_Data_Sim, 15, col_sim, 'filled', 'MarkerFaceAlpha', 0.3);
    
    x2 = (rand(size(Box_Data_Seq))-0.5)*0.2 + 2;
    scatter(x2, Box_Data_Seq, 15, col_seq, 'filled', 'MarkerFaceAlpha', 0.3);
else
    text(0.5, 0.5, 'No Data at Target Amp', 'HorizontalAlignment', 'center');
end

%% ================= 7. PRINT TABLE =================
fprintf('\n\n==========================================================\n');
fprintf('           GRAND AVERAGE DURATION (ms)          \n');
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
    if numel(x) < 2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; 
    valid = ~isnan(y) & ~isnan(se); 
    x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    
    upper = y + se; 
    lower = y - se;
    fill([x fliplr(x)], [upper fliplr(lower)], col, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end