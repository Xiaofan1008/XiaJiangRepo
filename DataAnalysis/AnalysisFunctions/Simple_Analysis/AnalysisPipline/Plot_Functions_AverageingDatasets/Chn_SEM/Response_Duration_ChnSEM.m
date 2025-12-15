%% ============================================================
%   POOLED GRAND AVERAGE: RESPONSE DURATION
%   - Metric: Mean +/- SEM across ALL responding channels (Pooled)
%   - Fig 1: Duration vs Amp (Shaded Error Bars)
%   - Fig 2: Box Plot (Distribution at Target Amp)
% ============================================================
clear; clc;

%% ================= 1. SETTINGS & FILES =================
target_amp_for_box = 10; 

file_list = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set1_Duration_SpecificPopulation_5ms_Xia_Exp1_Sim1_251125_112055.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set4_Duration_SpecificPopulation_5ms_Xia_Exp1_Sim4_251125_152849.mat';
};

%% ================= 2. DATA PREPARATION =================
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

% Storage for POOLED Data (Cells to hold raw vectors)
Pooled_Sim = cell(1, nAmps);
Pooled_Seq = cell(1, nAmps);

%% ================= 3. EXTRACT AND POOL RAW DATA =================
fprintf('Pooling responding channels from %d files...\n', length(file_list));

for f = 1:length(file_list)
    D = Data_Holder{f};
    these_amps = D.Data.Amps;
    
    % Map this file's amps to Master list
    [~, loc_in_master] = ismember(these_amps, Master_Amps);
    
    for i = 1:length(these_amps)
        target_idx = loc_in_master(i);
        
        % --- 1. Simultaneous (Raw Channels) ---
        raw_sim = D.Data.Dur_Sim(:, i);
        raw_sim = raw_sim(~isnan(raw_sim)); % Filter Non-Responders
        Pooled_Sim{target_idx} = [Pooled_Sim{target_idx}; raw_sim];
        
        % --- 2. Sequential (Raw Channels, Pooled across Sets) ---
        % Data is [Channels x Amps x Sets x PTD]
        raw_seq = D.Data.Dur_Seq(:, i, :, :);
        raw_seq = raw_seq(:); % Flatten dimensions
        raw_seq = raw_seq(~isnan(raw_seq)); % Filter Non-Responders
        Pooled_Seq{target_idx} = [Pooled_Seq{target_idx}; raw_seq];
    end
end

%% ================= 4. COMPUTE POOLED STATISTICS =================
Grand_Mean_Sim = nan(1, nAmps); Grand_SEM_Sim = nan(1, nAmps);
Grand_Mean_Seq = nan(1, nAmps); Grand_SEM_Seq = nan(1, nAmps);
N_Sim_Ch = nan(1, nAmps); % To track sample size

for i = 1:nAmps
    % Sim Stats
    d = Pooled_Sim{i};
    if ~isempty(d)
        Grand_Mean_Sim(i) = mean(d);
        Grand_SEM_Sim(i)  = std(d) / sqrt(length(d)); % Channel SEM
        N_Sim_Ch(i) = length(d);
    end
    
    % Seq Stats
    d = Pooled_Seq{i};
    if ~isempty(d)
        Grand_Mean_Seq(i) = mean(d);
        Grand_SEM_Seq(i)  = std(d) / sqrt(length(d)); % Channel SEM
    end
end

%% ================= 5. FIGURE 1: LINE PLOT (Shaded Error) =================
figure('Color','w', 'Position',[100 500 600 500]); hold on;
col_sim = [0 0.3 0.8];
col_seq = [0.85 0.33 0.10];

% Plot Shaded Errors
plot_shaded_error(Master_Amps, Grand_Mean_Sim, Grand_SEM_Sim, col_sim);
plot_shaded_error(Master_Amps, Grand_Mean_Seq, Grand_SEM_Seq, col_seq);

% Plot Mean Lines
plot(Master_Amps, Grand_Mean_Sim, '-o', 'Color', col_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');
plot(Master_Amps, Grand_Mean_Seq, '-s', 'Color', col_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

xlabel('Amplitude (µA)', 'FontWeight', 'bold');
ylabel('Net Active Duration (ms)', 'FontWeight', 'bold');
title(['Duration Response'], 'FontWeight', 'bold');
legend('Location', 'best', 'Box', 'off'); box off;

%% ================= 6. FIGURE 2: BOX PLOT (Distribution) =================
figure('Color','w', 'Position',[800 500 400 500]); hold on;

% Find Index for Target Amp
idx_box = find(Master_Amps == target_amp_for_box);

if ~isempty(idx_box)
    box_sim = Pooled_Sim{idx_box};
    box_seq = Pooled_Seq{idx_box};
    
    if ~isempty(box_sim) || ~isempty(box_seq)
        g1 = repmat({'Simultaneous'}, length(box_sim), 1);
        g2 = repmat({'Sequential'}, length(box_seq), 1);
        
        boxplot([box_sim; box_seq], [g1; g2], 'Width', 0.5, 'Symbol', '.');
        ylabel('Duration (ms)', 'FontWeight', 'bold');
        title(sprintf('Population Duration at %.0f µA', target_amp_for_box), 'FontWeight', 'bold');
        box off;
        
        % Add Jitter Points
        x1 = (rand(size(box_sim))-0.5)*0.2 + 1;
        scatter(x1, box_sim, 15, col_sim, 'filled', 'MarkerFaceAlpha', 0.3);
        
        x2 = (rand(size(box_seq))-0.5)*0.2 + 2;
        scatter(x2, box_seq, 15, col_seq, 'filled', 'MarkerFaceAlpha', 0.3);
    else
        text(0.5, 0.5, 'No Data at Target Amp', 'HorizontalAlignment', 'center');
    end
end

%% ================= 7. PRINT TABLE =================
fprintf('\n\n==========================================================\n');
fprintf('           POOLED DURATION RESULTS (ms)          \n');
fprintf('Amp(uA)\tN_Ch\tSim_Mean\tSim_SEM\tSeq_Mean\tSeq_SEM\n');
for i = 1:nAmps
    fprintf('%.0f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        Master_Amps(i), N_Sim_Ch(i), ...
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