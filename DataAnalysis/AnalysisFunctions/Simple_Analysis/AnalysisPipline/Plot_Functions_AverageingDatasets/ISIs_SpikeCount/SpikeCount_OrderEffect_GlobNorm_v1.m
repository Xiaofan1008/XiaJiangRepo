%% ============================================================
%   Population Spike Count Analysis (Order 1 vs Order 2)
%   - Normalization: Each animal is divided by its own GLOBAL MAX
%   - Safety: Verifies Set 1 and Set 2 are reversed orders of the same electrodes.
%   - Metric: Matched Paired Sets -> Mean ± SEM across animals (N)
%   - FILTER: Wilcoxon Signed Rank Test across ISIs (Console Output)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% 1. Define your datasets here. 
dataset_files = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX016/Result_SpikeCount_FixWin_DX016_10uA_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX016/Result_SpikeCount_FixWin_DX016_10uA_Xia_Exp1_Seq_Full_4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX018/Result_SpikeCount_FixWin_DX018_5_10uA_Xia_ISI_SimSeq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX018/Result_SpikeCount_FixWin_DX018_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq2.mat';
    % Add more file paths here as needed...
};

% 2. Plotting Aesthetics
line_width = 2;
marker_size = 6;
cap_size = 3; 

% 3. Statistical Settings
stats_target_Amp = 10; % Focus the Plot and Stats on this specific Amplitude
stats_min_n_threshold = 6; % Minimum pairs required to run Wilcoxon test

%% =================== 1. SCOUT LOOP (Build Master Union) ====================
fprintf('Scanning datasets to build Union Map...\n');
Union_Amps = [];
Union_ISIs = [];
max_Sets   = 0;
nFiles     = length(dataset_files);

for k = 1:nFiles
    load(dataset_files{k}, 'ResultFR');
    Union_Amps = union(Union_Amps, ResultFR.Metadata.TargetAmps);
    Union_ISIs = union(Union_ISIs, ResultFR.Metadata.TargetISIs);
    
    if isfield(ResultFR.Metadata, 'Stimulation_Sets')
        max_Sets = max(max_Sets, size(ResultFR.Metadata.Stimulation_Sets, 1));
    else
        max_Sets = max(max_Sets, size(ResultFR.SpikeCounts.Sim, 3));
    end
end

fprintf('Union Amplitudes found: %s uA\n', num2str(Union_Amps));
fprintf('Union ISIs found: %s ms\n', num2str(Union_ISIs));

% Pre-allocate the Master Population Array: [Amplitudes × ISIs × Sets × Animals]
Pop_Data = nan(length(Union_Amps), length(Union_ISIs), max_Sets, nFiles);

%% ================= 2. GATHER & NORMALIZE =================
fprintf('\nExtracting and Normalizing Data...\n');
for k = 1:nFiles
    load(dataset_files{k}, 'ResultFR');
    sim_data = ResultFR.SpikeCounts.Sim;
    seq_data = ResultFR.SpikeCounts.Seq;
    nSets_this = size(sim_data, 3);
    
    % --- THE BOUNCER: Check Dimension 3 (Sets) ---
    if nSets_this < 2
        continue; % Silently skip files without paired sets
    end
    
    % --- THE BOUNCER: Electrode Verification ---
    try
        sets_meta = ResultFR.Metadata.Stimulation_Sets;
        set1_ch = sets_meta(1, :); set1_ch = sort(set1_ch(set1_ch > 0));
        set2_ch = sets_meta(2, :); set2_ch = sort(set2_ch(set2_ch > 0));
        
        if ~isequal(set1_ch, set2_ch)
            fprintf('  Skipping File %d: Sets are NOT reverse orders.\n', k);
            continue;
        end
    catch
        fprintf('  Warning File %d: Metadata missing. Assuming reversed pairs.\n', k);
    end
    
    this_Amps = ResultFR.Metadata.TargetAmps;
    this_ISIs = ResultFR.Metadata.TargetISIs;
    
    % Extract Animal Means (Averaging across channels)
    animal_means = nan(length(Union_Amps), length(Union_ISIs), nSets_this);
    
    for ss = 1:nSets_this
        for a = 1:length(Union_Amps)
            amp_val = Union_Amps(a);
            ai_local = find(abs(this_Amps - amp_val) < 0.001);
            if isempty(ai_local), continue; end 
            
            for p = 1:length(Union_ISIs)
                isi_val = Union_ISIs(p);
                if isi_val == 0
                    ch_data = sim_data(:, ai_local, ss);
                else
                    pi_local = find(abs(this_ISIs - isi_val) < 0.001);
                    if isempty(pi_local), continue; end 
                    ch_data = seq_data(:, ai_local, ss, pi_local);
                end
                animal_means(a, p, ss) = nanmean(ch_data); 
            end
        end
    end
    
    % Global Max Normalization (across ALL valid Sets to maintain relative scale)
    valid_means = animal_means(~isnan(animal_means));
    if ~isempty(valid_means)
        Global_Max = max(valid_means);
        if Global_Max > 0
            Pop_Data(:, :, 1:nSets_this, k) = animal_means / Global_Max;
        end
    end
end

%% =================== 3. POPULATION STATISTICS =================
fprintf('\nCalculating Paired Population Mean and SEM...\n');

% Create a logical mask ensuring a rat has valid data for BOTH Order 1 and Order 2
% This guarantees perfect pairing for our statistical and graphical comparison
valid_pair_mask = ~isnan(Pop_Data(:,:,1,:)) & ~isnan(Pop_Data(:,:,2,:));

% Apply mask so unpaired data points are ignored
Pop_Data_Paired = Pop_Data;
Pop_Data_Paired(~repmat(valid_pair_mask, [1,1,max_Sets,1])) = NaN;

% Calculate Means for Order 1 and Order 2 across Dimension 4 (Animals)
Pop_Ord1_Mean = nanmean(Pop_Data_Paired(:,:,1,:), 4);
Pop_Ord2_Mean = nanmean(Pop_Data_Paired(:,:,2,:), 4);

% Dynamic N (same for both orders due to strict pairing mask)
Pop_N = sum(valid_pair_mask, 4);

% Calculate SEM
Pop_Ord1_SEM = nanstd(Pop_Data_Paired(:,:,1,:), 0, 4) ./ sqrt(Pop_N);
Pop_Ord2_SEM = nanstd(Pop_Data_Paired(:,:,2,:), 0, 4) ./ sqrt(Pop_N);

%% ===================== 4. PLOT (Target Amplitude Only) ======================
target_a_idx = find(abs(Union_Amps - stats_target_Amp) < 0.001);

if isempty(target_a_idx)
    error('Target amplitude (%.1f uA) not found in dataset union.', stats_target_Amp);
end

figure('Color','w', 'Position',[100 100 800 600]); 
hold on;

% Extract data strictly for the target amplitude
y1_mean = squeeze(Pop_Ord1_Mean(target_a_idx, :, 1, :));
y1_sem  = squeeze(Pop_Ord1_SEM(target_a_idx, :, 1, :));
y2_mean = squeeze(Pop_Ord2_Mean(target_a_idx, :, 1, :));
y2_sem  = squeeze(Pop_Ord2_SEM(target_a_idx, :, 1, :));
n_vals  = squeeze(Pop_N(target_a_idx, :, 1, :));

% The Bridge: Only plot points with valid pairs
valid_idx = (n_vals > 0);

if any(valid_idx)
    plot_x = Union_ISIs(valid_idx)';
    
    % Plot Order 1
    errorbar(plot_x, y1_mean(valid_idx), y1_sem(valid_idx), '-o', 'Color', 'k', 'LineWidth', line_width, ...
        'MarkerFaceColor', 'w', 'MarkerSize', marker_size, 'CapSize', cap_size, 'DisplayName', 'Order 1 (A \rightarrow B)');
    
    % Plot Order 2
    errorbar(plot_x, y2_mean(valid_idx), y2_sem(valid_idx), '--s', 'Color', 'k', 'LineWidth', line_width, ...
        'MarkerFaceColor', 'k', 'MarkerSize', marker_size, 'CapSize', cap_size, 'DisplayName', 'Order 2 (B \rightarrow A)');
end

% Figure Formatting
xlabel('Inter-Stimulus Interval (ms)', 'FontSize', 10); 
ylabel('Normalized Spike Count', 'FontSize', 10);
title_str = sprintf('Stimulation Order Effect across ISIs (%.1f µA)', stats_target_Amp);
title(title_str, 'FontSize', 12);
xticks(sort(Union_ISIs));
box off;
legend('Location','best','Box','off', 'FontSize', 10); 

%% ================= 5. CONSOLE STATISTICS TABLE =================
fprintf('\n=================================================================\n');
fprintf('   WILCOXON SIGNED-RANK TEST: ISI DIRECTIONALITY (%.1f uA)\n', stats_target_Amp);
fprintf('=================================================================\n');
fprintf(' ISI (ms) |   N pairs  |  Ord1 Mean |  Ord2 Mean |  p-value \n');
fprintf('-----------------------------------------------------------------\n');

for p = 1:length(Union_ISIs)
    isi = Union_ISIs(p);
    
    % Extract paired data for this specific ISI and Target Amp across all animals
    data_1 = squeeze(Pop_Data_Paired(target_a_idx, p, 1, :));
    data_2 = squeeze(Pop_Data_Paired(target_a_idx, p, 2, :));
    
    % Remove NaNs to get clean paired vectors
    valid_animals = ~isnan(data_1) & ~isnan(data_2);
    data_1 = data_1(valid_animals);
    data_2 = data_2(valid_animals);
    n_pairs = length(data_1);
    
    % --- FILTER: Check Min N ---
    if n_pairs < stats_min_n_threshold
        fprintf(' %8.1f | %10d |       Skipped (Low N < %d)\n', isi, n_pairs, stats_min_n_threshold);
        continue;
    end
    
    % --- STATS: Signed Rank Test ---
    pval = signrank(data_1, data_2);
    m1 = mean(data_1);
    m2 = mean(data_2);
    
    sig_flag = '';
    if pval < 0.05
        sig_flag = '<- SIGNIFICANT';
    else
        sig_flag = '(n.s.)';
    end
    
    fprintf(' %8.1f | %10d | %10.3f | %10.3f |  %.4f %s\n', isi, n_pairs, m1, m2, pval, sig_flag);
end
fprintf('=================================================================\n');
fprintf('Done! Directionality ISI figure generated.\n\n');