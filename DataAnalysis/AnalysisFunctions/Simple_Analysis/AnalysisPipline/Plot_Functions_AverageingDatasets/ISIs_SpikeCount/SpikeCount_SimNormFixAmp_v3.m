%% ============================================================
%   Population Spike Count Analysis (Universal Normalizer)
%   - Normalization: Dynamically subtracts or divides by a target baseline
%   - Logic: Union of Amplitudes and ISIs from CONDENSED .mat files
%   - Metric: Pooled Sets -> Mean ± SEM across animals (N)
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
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX019/Result_SpikeCount_FixWin_DX019_5_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX022/Result_SpikeCount_FixWin_DX022_10uA_Xia_ISI_10uA_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX023/Result_SpikeCount_FixWin_DX023_10uA_Xia_ISI_10uA_SinSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX024/Result_SpikeCount_FixWin_DX024_10uA_Xia_ISI_10uA_SimSeq1.mat';
};

% 2. Plotting Aesthetics
line_width = 2.5;
marker_size = 5;
cap_size = 3.5; % Size of the error bar caps

% 3. Statistical Settings
stats_target_Amp = 10; % Select the specific Amplitude to run the LMM on

% [MODIFIED] 4. Universal Normalization Settings
norm_target_Amp = 10;           % Amplitude to use for baseline normalization (e.g., 10)
norm_target_ISI = 0;            % ISI time (ms) to use for baseline normalization (e.g., 0)
norm_method     = 'subtract';   % CHOOSE ONE: 'subtract' or 'divide'

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

% Pre-allocate the Master Population Array
Pop_Data = nan(length(Union_Amps), length(Union_ISIs), max_Sets, nFiles);

%% ================= 2. GATHER & NORMALIZE (The Core Engine) =============
fprintf('\nExtracting and Normalizing Data using method: [%s]\n', upper(norm_method));
for k = 1:nFiles
    fprintf('  Processing File %d/%d...\n', k, nFiles);
    
    load(dataset_files{k}, 'ResultFR');
    sim_data = ResultFR.SpikeCounts.Sim;
    seq_data = ResultFR.SpikeCounts.Seq;
    nSets_this = size(sim_data, 3);
    
    this_Amps = ResultFR.Metadata.TargetAmps;
    this_ISIs = ResultFR.Metadata.TargetISIs;
    
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
    
    % [MODIFIED] 3. Dynamic Normalization Engine
    idx_normAmp = find(abs(Union_Amps - norm_target_Amp) < 0.001);
    idx_normISI = find(abs(Union_ISIs - norm_target_ISI) < 0.001);
    
    if ~isempty(idx_normAmp) && ~isempty(idx_normISI)
        for ss = 1:nSets_this
            % Get the target baseline spike count
            baseline_val = animal_means(idx_normAmp, idx_normISI, ss);
            
            % If we have a valid baseline, apply the chosen mathematical operation
            if ~isnan(baseline_val)
                if strcmpi(norm_method, 'subtract')
                    Pop_Data(:, :, ss, k) = animal_means(:, :, ss) - baseline_val;
                    
                elseif strcmpi(norm_method, 'divide')
                    if baseline_val > 0
                        Pop_Data(:, :, ss, k) = animal_means(:, :, ss) / baseline_val;
                    else
                        fprintf('    -> Warning: Baseline is 0. Cannot perform division.\n');
                    end
                    
                else
                    error('Unknown norm_method: %s. Please use ''subtract'' or ''divide''.', norm_method);
                end
            else
                fprintf('    -> Warning: No valid %.1fuA at %.1fms baseline found. Excluded from average.\n', norm_target_Amp, norm_target_ISI);
            end
        end
    else
        fprintf('    -> Warning: %.1fuA or %.1fms missing from dataset union. Cannot normalize.\n', norm_target_Amp, norm_target_ISI);
    end
end

%% =================== 3. POPULATION STATISTICS =================
fprintf('\nCalculating Population Mean and SEM...\n');
Pop_Data_Pooled = nanmean(Pop_Data, 3);
Pop_Mean = nanmean(Pop_Data_Pooled, 4);
Pop_N = sum(~isnan(Pop_Data_Pooled), 4);
Pop_Std = nanstd(Pop_Data_Pooled, 0, 4);
Pop_SEM = Pop_Std ./ sqrt(Pop_N);

%% =================== 3.5 STATISTICAL TEST (LMM) =================
target_a_idx = find(abs(Union_Amps - stats_target_Amp) < 0.001);
if ~isempty(target_a_idx)
    tbl_Spikes = []; tbl_ISI = []; tbl_Animal = [];
    
    for k = 1:nFiles
        for p = 1:length(Union_ISIs)
            val = Pop_Data_Pooled(target_a_idx, p, 1, k);
            if ~isnan(val)
                tbl_Spikes = [tbl_Spikes; val];
                tbl_ISI    = [tbl_ISI; Union_ISIs(p)];
                tbl_Animal = [tbl_Animal; k];
            end
        end
    end
    
    T = table(tbl_Spikes, categorical(tbl_ISI), categorical(tbl_Animal), ...
        'VariableNames', {'SpikeCount', 'ISI', 'AnimalID'});
    
    fprintf('\n--------------------------------------------------\n');
    fprintf('Running Linear Mixed-Effects Model for %.1f uA\n', stats_target_Amp);
    fprintf('Formula: SpikeCount ~ 1 + ISI + (1 | AnimalID)\n');
    fprintf('--------------------------------------------------\n');
    
    lme = fitlme(T, 'SpikeCount ~ 1 + ISI + (1|AnimalID)');
    disp(anova(lme));
else
    fprintf('\nWarning: Target amplitude for stats (%.1f uA) not found in datasets.\n', stats_target_Amp);
end

%% ===================== 4. PLOT (Pooled Population Dose-Response) ======================
num_amps = length(Union_Amps);
gray_shades = linspace(0.6, 0, num_amps)'; 
colors = [gray_shades, gray_shades, gray_shades]; 

line_styles = {'-', '--', ':', '-.'}; 
markers = {'s', 'o', '^', 'd'};

figure('Color','w', 'Position',[100 100 800 600]); 
hold on;

for a = 1:length(Union_Amps)
    current_amp = Union_Amps(a);
    y_mean = squeeze(Pop_Mean(a, :));
    y_sem  = squeeze(Pop_SEM(a, :));
    n_vals = squeeze(Pop_N(a, :));
    
    valid_idx = (n_vals > 0) & ~isnan(y_mean);
    
    if any(valid_idx)
        col = colors(a, :);
        ls  = line_styles{mod(a-1, length(line_styles))+1};
        mk  = markers{mod(a-1, length(markers))+1};
        lbl = sprintf('%.1f uA', current_amp); 
        
        plot_x   = Union_ISIs(valid_idx);
        plot_y   = y_mean(valid_idx);
        plot_err = y_sem(valid_idx);
        
        errorbar(plot_x, plot_y, plot_err, 'LineStyle', ls, 'Marker', mk, 'Color', col, ...
            'LineWidth', line_width, 'MarkerFaceColor', 'w', 'MarkerSize', marker_size, ...
            'CapSize', cap_size, 'DisplayName', lbl);
    end
end

% Figure Formatting
xlabel('Inter-Stimulus Interval (ms)', 'FontSize', 12); 

% [MODIFIED] Dynamic Y-axis label based on user settings
if strcmpi(norm_method, 'subtract')
    ylabel_str = sprintf('\\Delta Spike Count (relative to %.1f \\muA at %.1f ms)', norm_target_Amp, norm_target_ISI);
elseif strcmpi(norm_method, 'divide')
    ylabel_str = sprintf('Normalized Spike Count (Fraction of %.1f \\muA at %.1f ms)', norm_target_Amp, norm_target_ISI);
end
ylabel(ylabel_str, 'FontSize', 12);

title_str = sprintf('ISI Population Average');
title(title_str, 'FontSize', 14);

xticks(sort(Union_ISIs));
box off;
lgd = legend('Location','best','Box','off'); 
title(lgd, 'Amplitudes'); 

fprintf('Figure generated.\n');