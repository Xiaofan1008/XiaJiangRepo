%% ============================================================
%   POOLED GRAND AVERAGE (Small N Strategy)
%   - Logic: Concatenates RAW channels from all datasets
%   - Metric: Mean +/- SEM across ALL recorded neurons
%   - Added: Scatter points to visualize distribution
% ============================================================
clear; clc;

%% ================= 1. DEFINE FILES =================
file_list = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set1_SpikeNormSimMax_Zeroed_5ms_Xia_Exp1_Sim1_251125_112055.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set6_SpikeNormSimMax_Zeroed_5ms_Xia_Exp1_Sim6_251125_181554.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/Result_Set4_SpikeNormSimMax_Zeroed_5ms_Xia_Exp1_Sim4_251125_152849.mat'
};

%% ================= 2. SCAN FOR AMPLITUDES =================
All_Amps_Found = [];
for f = 1:length(file_list)
    loaded = load(file_list{f});
    
    if isfield(loaded, 'ResultNorm'), D = loaded.ResultNorm;
    elseif isfield(loaded, 'ResultFR'), D = loaded.ResultFR; 
    elseif isfield(loaded, 'ResultDur'), D = loaded.ResultDur;
    else, error('Structure not found'); end
    
    if isfield(D, 'Amps'), amps = D.Amps;
    elseif isfield(D, 'Metadata'), amps = D.Metadata.Amps; 
    else, amps = D.Data.Amps; end
    
    All_Amps_Found = [All_Amps_Found; amps];
end
Master_Amps = unique(All_Amps_Found);
nAmps = length(Master_Amps);
Pooled_Sim = cell(1, nAmps);
Pooled_Seq = cell(1, nAmps);

%% ================= 3. EXTRACT AND POOL DATA =================
fprintf('Pooling data from %d files...\n', length(file_list));
for f = 1:length(file_list)
    loaded = load(file_list{f});
    
    if isfield(loaded, 'ResultNorm')
        D = loaded.ResultNorm;
        Raw_Sim = D.Norm_Sim;
        Raw_Seq = D.Norm_Seq;
        these_amps = D.Amps;
    elseif isfield(loaded, 'ResultFR')
        D = loaded.ResultFR;
        Raw_Sim = D.PeakFR.Sim_Z;       
        Raw_Seq = D.PeakFR.Seq_Z;       
        these_amps = D.Metadata.Amps;
    elseif isfield(loaded, 'ResultDur')
        D = loaded.ResultDur;
        Raw_Sim = D.Data.Dur_Sim;       
        Raw_Seq = D.Data.Dur_Seq;       
        these_amps = D.Data.Amps;
    end
    
    [~, loc_in_master] = ismember(these_amps, Master_Amps);
    
    for i = 1:length(these_amps)
        target_idx = loc_in_master(i);
        
        vals_sim = Raw_Sim(:, i);
        vals_sim = vals_sim(~isnan(vals_sim)); 
        Pooled_Sim{target_idx} = [Pooled_Sim{target_idx}; vals_sim];
        
        vals_seq = Raw_Seq(:, i, :, :); 
        vals_seq = vals_seq(:); 
        vals_seq = vals_seq(~isnan(vals_seq)); 
        Pooled_Seq{target_idx} = [Pooled_Seq{target_idx}; vals_seq];
    end
end

%% ================= 4. CALCULATE STATS =================
Grand_Mean_Sim = nan(1, nAmps); Grand_SEM_Sim = nan(1, nAmps);
Grand_Mean_Seq = nan(1, nAmps); Grand_SEM_Seq = nan(1, nAmps);
N_Sim = nan(1, nAmps); 

for i = 1:nAmps
    d = Pooled_Sim{i};
    if ~isempty(d)
        Grand_Mean_Sim(i) = mean(d);
        Grand_SEM_Sim(i)  = std(d) / sqrt(length(d));
        N_Sim(i) = length(d);
    end
    d = Pooled_Seq{i};
    if ~isempty(d)
        Grand_Mean_Seq(i) = mean(d);
        Grand_SEM_Seq(i)  = std(d) / sqrt(length(d));
    end
end

%% ================= 5. PLOT WITH SCATTERS =================
figure('Color','w', 'Position',[500 500 700 500]); hold on;
col_sim = [0 0.3 0.8]; col_seq = [0.85 0.33 0.10];

% --- Plot Scatters (Added Block) ---
jitter = 0.2; % Horizontal spread
sz = 15;      % Dot size
alpha = 0.25; % Transparency

for i = 1:nAmps
    % Sim Scatters
    y = Pooled_Sim{i};
    if ~isempty(y)
        x = Master_Amps(i) + (rand(size(y)) - 0.5) * jitter;
        scatter(x, y, sz, col_sim, 'filled', ...
            'MarkerFaceAlpha', alpha, 'HandleVisibility', 'off');
    end
    
    % Seq Scatters
    y = Pooled_Seq{i};
    if ~isempty(y)
        x = Master_Amps(i) + (rand(size(y)) - 0.5) * jitter;
        scatter(x, y, sz, col_seq, 'filled', ...
            'MarkerFaceAlpha', alpha, 'HandleVisibility', 'off');
    end
end
% -----------------------------------

% Plot Curves
plot_shaded_error(Master_Amps, Grand_Mean_Sim, Grand_SEM_Sim, col_sim);
plot_shaded_error(Master_Amps, Grand_Mean_Seq, Grand_SEM_Seq, col_seq);

plot(Master_Amps, Grand_Mean_Sim, '-o', 'Color', col_sim, 'LineWidth', 2, ...
    'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');
plot(Master_Amps, Grand_Mean_Seq, '-s', 'Color', col_seq, 'LineWidth', 2, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Sequential');

xlabel('Amplitude (\muA)', 'FontWeight', 'bold');
ylabel('Normalized Spike Count per Trial', 'FontWeight', 'bold');
title(['Normalized Spike Count (Pooled)'], 'FontWeight', 'bold');
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