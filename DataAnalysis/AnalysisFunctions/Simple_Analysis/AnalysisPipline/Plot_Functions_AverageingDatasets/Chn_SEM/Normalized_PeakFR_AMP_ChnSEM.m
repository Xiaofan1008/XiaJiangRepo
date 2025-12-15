%% ============================================================
%   POOLED GRAND AVERAGE (Small N Strategy)
%   - Logic: Concatenates RAW channels from all datasets
%   - Metric: Mean +/- SEM across ALL recorded neurons
%   - Benefit: Reduces error bars when N_subjects is small
% ============================================================
clear; clc;

%% ================= 1. DEFINE FILES =================
file_list = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set1_PeakFiringRate_FixedPop_5ms_Xia_Exp1_Sim1_251125_112055.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set4_PeakFiringRate_FixedPop_5ms_Xia_Exp1_Sim4_251125_152849.mat';
};

%% ================= 2. SCAN FOR AMPLITUDES =================
% We need to know every amplitude tested to build our storage containers
All_Amps_Found = [];
for f = 1:length(file_list)
    loaded = load(file_list{f});
    if isfield(loaded, 'ResultFR'), D = loaded.ResultFR; 
    elseif isfield(loaded, 'ResultDur'), D = loaded.ResultDur; % Works for Duration too
    else, error('Structure not found'); end
    
    if isfield(D, 'Metadata'), amps = D.Metadata.Amps; else, amps = D.Data.Amps; end
    All_Amps_Found = [All_Amps_Found; amps];
end
Master_Amps = unique(All_Amps_Found);
nAmps = length(Master_Amps);

% Storage: Cell array to hold raw vectors for each amplitude
% {Amp1_Vector, Amp2_Vector, ...}
Pooled_Sim = cell(1, nAmps);
Pooled_Seq = cell(1, nAmps);

%% ================= 3. EXTRACT AND POOL RAW DATA =================
fprintf('Pooling data from %d files...\n', length(file_list));

for f = 1:length(file_list)
    loaded = load(file_list{f});
    
    % Handle different structure names (FR vs Duration)
    if isfield(loaded, 'ResultFR')
        D = loaded.ResultFR;
        Raw_Sim = D.PeakFR.Sim_Z;       % [Ch x Amp]
        Raw_Seq = D.PeakFR.Seq_Z;       % [Ch x Amp x Set]
        these_amps = D.Metadata.Amps;
    elseif isfield(loaded, 'ResultDur')
        D = loaded.ResultDur;
        Raw_Sim = D.Data.Dur_Sim;       % [Ch x Amp]
        Raw_Seq = D.Data.Dur_Seq;       % [Ch x Amp x Set]
        these_amps = D.Data.Amps;
    end
    
    % Align this file's amplitudes to Master list
    [~, loc_in_master] = ismember(these_amps, Master_Amps);
    
    % --- LOOP THROUGH AMPLITUDES IN THIS FILE ---
    for i = 1:length(these_amps)
        target_idx = loc_in_master(i);
        
        % 1. Get Sim Data (Vector of channels)
        vals_sim = Raw_Sim(:, i);
        vals_sim = vals_sim(~isnan(vals_sim)); % Remove NaNs
        Pooled_Sim{target_idx} = [Pooled_Sim{target_idx}; vals_sim];
        
        % 2. Get Seq Data (Vector of channels, pooled across sets)
        % Flatten the sets dimension so we just have a list of all seq responses
        vals_seq = Raw_Seq(:, i, :, :); 
        vals_seq = vals_seq(:); 
        vals_seq = vals_seq(~isnan(vals_seq)); % Remove NaNs
        Pooled_Seq{target_idx} = [Pooled_Seq{target_idx}; vals_seq];
    end
end

%% ================= 4. CALCULATE POOLED STATS =================
Grand_Mean_Sim = nan(1, nAmps); Grand_SEM_Sim = nan(1, nAmps);
Grand_Mean_Seq = nan(1, nAmps); Grand_SEM_Seq = nan(1, nAmps);
N_Sim = nan(1, nAmps); % Track how many neurons contributed

for i = 1:nAmps
    % Sim Stats
    d = Pooled_Sim{i};
    if ~isempty(d)
        Grand_Mean_Sim(i) = mean(d);
        Grand_SEM_Sim(i)  = std(d) / sqrt(length(d));
        N_Sim(i) = length(d);
    end
    
    % Seq Stats
    d = Pooled_Seq{i};
    if ~isempty(d)
        Grand_Mean_Seq(i) = mean(d);
        Grand_SEM_Seq(i)  = std(d) / sqrt(length(d));
    end
end

%% ================= 5. PLOT =================
figure('Color','w', 'Position',[500 500 700 500]); hold on;
col_sim = [0 0.3 0.8]; col_seq = [0.85 0.33 0.10];

% Helper function for shading
plot_shaded_error(Master_Amps, Grand_Mean_Sim, Grand_SEM_Sim, col_sim);
plot_shaded_error(Master_Amps, Grand_Mean_Seq, Grand_SEM_Seq, col_seq);

plot(Master_Amps, Grand_Mean_Sim, '-o', 'Color', col_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');
plot(Master_Amps, Grand_Mean_Seq, '-s', 'Color', col_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

xlabel('Amplitude (\muA)', 'FontWeight', 'bold');
ylabel('Normalized Peak Firing Rate (sps)', 'FontWeight', 'bold'); % Change label manually
title(['Normalized Peak Firing Rate'], 'FontWeight', 'bold');
legend('Location', 'best', 'Box', 'off'); box off;

%% ================= 6. PRINT TABLE =================
fprintf('\n\n=== POOLED POPULATION RESULTS ===\n');
fprintf('Amp\tN_Ch\tSim_Mean\tSim_SEM\tSeq_Mean\tSeq_SEM\n');
for i = 1:nAmps
    fprintf('%.0f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
        Master_Amps(i), N_Sim(i), ...
        Grand_Mean_Sim(i), Grand_SEM_Sim(i), ...
        Grand_Mean_Seq(i), Grand_SEM_Seq(i));
end

%% ================= HELPER =================
function plot_shaded_error(x, y, se, col)
    if numel(x) < 2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; 
    valid = ~isnan(y) & ~isnan(se); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end