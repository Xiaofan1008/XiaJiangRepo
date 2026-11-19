%% ============================================================
%   Population Spike Count Analysis (Global Max Normalized)
%   - Normalization: Each animal is divided by its own GLOBAL MAX
%   - Logic: Union of Amplitudes and ISIs from CONDENSED .mat files
%   - Metric: Mean ± SEM across animals (N)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% 1. Define your datasets here. 
% [MODIFIED] Since the .mat files are now fully self-contained,
% you ONLY need to provide the paths to the result files.

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
marker_size = 8;
cap_size = 5; % Size of the error bar caps

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
    
    % [MODIFIED] Extract local metadata instead of loading raw hardware maps
    this_Amps = ResultFR.Metadata.TargetAmps;
    this_ISIs = ResultFR.Metadata.TargetISIs;
    
    % 2. Extract Animal Means (Averaging across channels)
    animal_means = nan(length(Union_Amps), length(Union_ISIs), nSets_this);
    
    for ss = 1:nSets_this
        for a = 1:length(Union_Amps)
            amp_val = Union_Amps(a);
            
            % [MODIFIED] Match against the local condensed metadata
            ai_local = find(abs(this_Amps - amp_val) < 0.001);
            if isempty(ai_local), continue; end % Rat didn't test this Amp
            
            for p = 1:length(Union_ISIs)
                isi_val = Union_ISIs(p);
                
                if isi_val == 0
                    ch_data = sim_data(:, ai_local, ss);
                else
                    % [MODIFIED] Match against the local condensed metadata
                    pi_local = find(abs(this_ISIs - isi_val) < 0.001);
                    if isempty(pi_local), continue; end % Rat didn't test this ISI
                    ch_data = seq_data(:, ai_local, ss, pi_local);
                end
                
                % nanmean safely ignores scrubbed artifacts!
                animal_means(a, p, ss) = nanmean(ch_data); 
            end
        end
    end
    
    % 3. Global Max Normalization
    % Find absolute max value across all valid Amps/ISIs/Sets for THIS animal
    valid_means = animal_means(~isnan(animal_means));
    if ~isempty(valid_means)
        Global_Max = max(valid_means);
        
        % Divide all data by max and slot into Master Array
        if Global_Max > 0
            Pop_Data(:, :, 1:nSets_this, k) = animal_means / Global_Max;
        end
    end
end

%% =================== 3. POPULATION STATISTICS =================
fprintf('\nCalculating Population Mean and SEM...\n');

% Calculate Mean across the 4th dimension (Animals)
Pop_Mean = nanmean(Pop_Data, 4);

% Calculate dynamic N (number of valid animals per specific point)
Pop_N = sum(~isnan(Pop_Data), 4);

% Calculate SEM = std / sqrt(N)
Pop_Std = nanstd(Pop_Data, 0, 4);
Pop_SEM = Pop_Std ./ sqrt(Pop_N);

%% ===================== 4. PLOT (Population Dose-Response) ======================
amp_colors = lines(length(Union_Amps)); 

for ss = 1:max_Sets
    figure('Color','w', 'Position',[100 100 800 600]); 
    hold on;
    
    % Loop through each Union Amplitude
    for a = 1:length(Union_Amps)
        current_amp = Union_Amps(a);
        
        % Extract data for this Set and Amp
        y_mean = squeeze(Pop_Mean(a, :, ss));
        y_sem  = squeeze(Pop_SEM(a, :, ss));
        n_vals = squeeze(Pop_N(a, :, ss));
        
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
    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold', 'FontSize', 12); 
    ylabel('Normalized Net Mean Spikes (Max = 1.0)', 'FontWeight','bold', 'FontSize', 12);
    
    title_str = sprintf('Population Average: Set %d (N = %d Datasets)', ss, nFiles);
    title(title_str, 'FontWeight','bold', 'FontSize', 14);
    
    xticks(sort(Union_ISIs));
    box off;
    lgd = legend('Location','best','Box','off'); 
    title(lgd, 'Amplitudes'); 
end

fprintf('Done! Population figures generated.\n');