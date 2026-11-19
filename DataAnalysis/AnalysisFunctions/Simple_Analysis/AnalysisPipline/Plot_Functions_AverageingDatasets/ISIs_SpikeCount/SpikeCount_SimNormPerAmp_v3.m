%% ============================================================
%   Population Spike Count Analysis (Within-Amplitude Baseline)
%   - Normalization: Each Amplitude is divided by its own ISI=0 baseline
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
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX022/Result_SpikeCount_FixWin_DX022_10uA_Xia_ISI_10uA_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX023/Result_SpikeCount_FixWin_DX023_10uA_Xia_ISI_10uA_SinSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX024/Result_SpikeCount_FixWin_DX024_10uA_Xia_ISI_10uA_SimSeq1.mat';
    % Add more file paths here as needed...
};

% 2. Plotting Aesthetics
line_width = 2.5;
marker_size = 5;
cap_size = 3.5; % Size of the error bar caps

% 3. Statistical Settings
stats_target_Amp = 10; % Select the specific Amplitude to run the LMM on

%% =================== 1. SCOUT LOOP (Build Master Union) ====================
fprintf('Scanning datasets to build Union Map...\n');
Union_Amps = [];
Union_ISIs = [];
max_Sets   = 0;
nFiles     = length(dataset_files);

for k = 1:nFiles
    load(dataset_files{k}, 'ResultFR');
    
    % Collect all unique Amps and ISIs found across all files
    Union_Amps = union(Union_Amps, ResultFR.Metadata.TargetAmps);
    Union_ISIs = union(Union_ISIs, ResultFR.Metadata.TargetISIs);
    
    % Track the maximum number of Stimulation Sets
    if isfield(ResultFR.Metadata, 'Stimulation_Sets')
        max_Sets = max(max_Sets, size(ResultFR.Metadata.Stimulation_Sets, 1));
    else
        max_Sets = max(max_Sets, size(ResultFR.SpikeCounts.Sim, 3));
    end
end
fprintf('Union Amplitudes found: %s uA\n', num2str(Union_Amps));
fprintf('Union ISIs found: %s ms\n', num2str(Union_ISIs));

% Pre-allocate the Master Population Array
% Dimensions: [Amplitudes × ISIs × Sets × Animals]
Pop_Data = nan(length(Union_Amps), length(Union_ISIs), max_Sets, nFiles);

%% ================= 2. GATHER & NORMALIZE (The Core Engine) =============
fprintf('\nExtracting and Normalizing Data...\n');
for k = 1:nFiles
    fprintf('  Processing File %d/%d...\n', k, nFiles);
    
    % 1. Load Result Data
    load(dataset_files{k}, 'ResultFR');
    sim_data = ResultFR.SpikeCounts.Sim;
    seq_data = ResultFR.SpikeCounts.Seq;
    nSets_this = size(sim_data, 3);
    
    % Extract local metadata instead of loading raw hardware maps
    this_Amps = ResultFR.Metadata.TargetAmps;
    this_ISIs = ResultFR.Metadata.TargetISIs;
    
    % 2. Extract Animal Means (Averaging across channels)
    animal_means = nan(length(Union_Amps), length(Union_ISIs), nSets_this);
    
    for ss = 1:nSets_this
        for a = 1:length(Union_Amps)
            amp_val = Union_Amps(a);
            
            % Match against the local condensed metadata
            ai_local = find(abs(this_Amps - amp_val) < 0.001);
            if isempty(ai_local), continue; end % Rat didn't test this Amp
            
            for p = 1:length(Union_ISIs)
                isi_val = Union_ISIs(p);
                
                if isi_val == 0
                    ch_data = sim_data(:, ai_local, ss);
                else
                    % Match against the local condensed metadata
                    pi_local = find(abs(this_ISIs - isi_val) < 0.001);
                    if isempty(pi_local), continue; end % Rat didn't test this ISI
                    ch_data = seq_data(:, ai_local, ss, pi_local);
                end
                
                % nanmean safely ignores scrubbed artifacts!
                animal_means(a, p, ss) = nanmean(ch_data); 
            end
        end
    end
    
    % [MODIFIED] 3. Within-Amplitude Baseline Normalization (ISI = 0)
    % Find the index corresponding to ISI = 0 ms in the Union Map
    idx_0ms = find(abs(Union_ISIs - 0) < 0.001);
    
    if ~isempty(idx_0ms)
        for ss = 1:nSets_this
            for a = 1:length(Union_Amps)
                % Get the baseline spike count at 0 ms for this specific Amp and Set
                baseline_val = animal_means(a, idx_0ms, ss);
                
                % If we have a valid baseline > 0, divide all ISIs for this Amp by it
                if ~isnan(baseline_val) && baseline_val > 0
                    Pop_Data(a, :, ss, k) = animal_means(a, :, ss) / baseline_val;
                end
            end
        end
    end
end

%% =================== 3. POPULATION STATISTICS =================
fprintf('\nCalculating Population Mean and SEM...\n');

% POOLING: Average across the Sets (Dimension 3) FIRST
% This smooths directional variance and yields [Amplitudes x ISIs x 1 x Animals]
Pop_Data_Pooled = nanmean(Pop_Data, 3);

% Calculate Mean across the 4th dimension (Animals) using POOLED data
Pop_Mean = nanmean(Pop_Data_Pooled, 4);

% Calculate dynamic N (number of valid animals per specific point) using POOLED data
Pop_N = sum(~isnan(Pop_Data_Pooled), 4);

% Calculate SEM = std / sqrt(N) using POOLED data
Pop_Std = nanstd(Pop_Data_Pooled, 0, 4);
Pop_SEM = Pop_Std ./ sqrt(Pop_N);

%% =================== 3.5 STATISTICAL TEST (LMM) =================
% Run Linear Mixed-Effects Model for the target amplitude
target_a_idx = find(abs(Union_Amps - stats_target_Amp) < 0.001);

if ~isempty(target_a_idx)
    % Initialize empty arrays for Long-Format table
    tbl_Spikes = [];
    tbl_ISI    = [];
    tbl_Animal = [];
    
    % Extract valid data points into the lists
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
    
    % Create MATLAB Table. ISI is categorical to test its specific steps.
    % Animal is categorical to act as the random effect grouping variable.
    T = table(tbl_Spikes, categorical(tbl_ISI), categorical(tbl_Animal), ...
        'VariableNames', {'SpikeCount', 'ISI', 'AnimalID'});
    
    % Fit the LMM: SpikeCount predicted by ISI (Fixed), accounting for AnimalID (Random)
    fprintf('\n--------------------------------------------------\n');
    fprintf('Running Linear Mixed-Effects Model for %.1f uA\n', stats_target_Amp);
    fprintf('Formula: SpikeCount ~ 1 + ISI + (1 | AnimalID)\n');
    fprintf('--------------------------------------------------\n');
    
    lme = fitlme(T, 'SpikeCount ~ 1 + ISI + (1|AnimalID)');
    
    % Display the ANOVA table which shows the main effect of ISI
    disp(anova(lme));
else
    fprintf('\nWarning: Target amplitude for stats (%.1f uA) not found in datasets.\n', stats_target_Amp);
end

%% ===================== 4. PLOT (Pooled Population Dose-Response) ======================
amp_colors = lines(length(Union_Amps)); 

% Removed the 'for ss = 1:max_Sets' loop. We only generate ONE master figure.
figure('Color','w', 'Position',[100 100 800 600]); 
hold on;

% Loop through each Union Amplitude
for a = 1:length(Union_Amps)
    current_amp = Union_Amps(a);
    
    % Extract data for this Amp (Sets dimension is already pooled out)
    y_mean = squeeze(Pop_Mean(a, :));
    y_sem  = squeeze(Pop_SEM(a, :));
    n_vals = squeeze(Pop_N(a, :));
    
    % [THE BRIDGE]: Only plot points where N >= 1
    valid_idx = (n_vals > 0) & ~isnan(y_mean);
    
    if any(valid_idx)
        col = amp_colors(a, :);
        lbl = sprintf('%.1f uA', current_amp); 
        
        plot_x   = Union_ISIs(valid_idx);
        plot_y   = y_mean(valid_idx);
        plot_err = y_sem(valid_idx);
        
        % Use errorbar to plot Mean ± SEM
        errorbar(plot_x, plot_y, plot_err, '-o', 'Color', col, 'LineWidth', line_width, ...
            'MarkerFaceColor', 'w', 'MarkerSize', marker_size, 'CapSize', cap_size, 'DisplayName', lbl);
    end
end

% Figure Formatting
xlabel('Inter-Stimulus Interval (ms)', 'FontSize', 12); 
% [MODIFIED] Updated Y-axis label to reflect within-amplitude normalization
ylabel('Normalizaed Spike Count (Fraction of Simultaneous Response)', 'FontSize', 12);

% Updated Title
title_str = sprintf('ISI Population Average');
title(title_str, 'FontSize', 12);

xticks(sort(Union_ISIs));
box off;
lgd = legend('Location','best','Box','off'); 
title(lgd, 'Amplitudes'); 
% axis square;

fprintf('Figure generated.\n');


%% ===================== 4. PLOT (Pooled Population Dose-Response) ======================
% 1. Style Dictionary for Grayscale Distinctiveness
num_amps = length(Union_Amps);
% Generate shades from 0.6 (Medium Gray) to 0 (Solid Black)
gray_shades = linspace(0.6, 0, num_amps)'; 
colors = [gray_shades, gray_shades, gray_shades]; 

% Rotate through distinct line styles and markers
line_styles = {'-', '--', ':', '-.'}; 
markers = {'s', 'o', '^', 'd'};

% Removed the 'for ss = 1:max_Sets' loop. We only generate ONE master figure.
figure('Color','w', 'Position',[100 100 800 600]); 
hold on;

% Loop through each Union Amplitude
for a = 1:length(Union_Amps)
    current_amp = Union_Amps(a);
    
    % Extract data for this Amp (Sets dimension is already pooled out)
    y_mean = squeeze(Pop_Mean(a, :));
    y_sem  = squeeze(Pop_SEM(a, :));
    n_vals = squeeze(Pop_N(a, :));
    
    % [THE BRIDGE]: Only plot points where N >= 1
    valid_idx = (n_vals > 0) & ~isnan(y_mean);
    
    if any(valid_idx)
        % Select aesthetics for this amplitude
        col = colors(a, :);
        ls  = line_styles{mod(a-1, length(line_styles))+1};
        mk  = markers{mod(a-1, length(markers))+1};
        lbl = sprintf('%.1f uA', current_amp); 
        
        plot_x   = Union_ISIs(valid_idx);
        plot_y   = y_mean(valid_idx);
        plot_err = y_sem(valid_idx);
        
        % Use classic errorbar with distinct grayscale styling
        errorbar(plot_x, plot_y, plot_err, 'LineStyle', ls, 'Marker', mk, 'Color', col, ...
            'LineWidth', line_width, 'MarkerFaceColor', 'w', 'MarkerSize', marker_size, ...
            'CapSize', cap_size, 'DisplayName', lbl);
    end
end

% Figure Formatting
xlabel('Inter-Stimulus Interval (ms)', 'FontSize', 12); 
% [MODIFIED] Updated Y-axis label to reflect within-amplitude normalization
ylabel('Normalized Spike Count (Fraction of Simultaneous Response)', 'FontSize', 12);

% Updated Title
title_str = sprintf('ISI Population Average');
title(title_str, 'FontSize', 12);

xticks(sort(Union_ISIs));
box off;
lgd = legend('Location','best','Box','off'); 
title(lgd, 'Amplitudes'); 

fprintf('Figure generated.\n');