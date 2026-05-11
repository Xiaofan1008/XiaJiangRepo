%% ============================================================
%   Population Spike Count Analysis (Within-Set Max Normalized)
%   - Normalization: Each SET is divided independently by its own MAX response
%   - Logic: Union of Amplitudes and ISIs from CONDENSED .mat files
%   - Metric: Pooled Sets -> Mean ± SEM across animals (N)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% 1. Define your datasets here. 
dataset_files = {
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX014/Result_SpikeCount_FixWin_DX014_5_10uA_Xia_Seq_Sim1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX014/Result_SpikeCount_FixWin_DX014_5_10uA_Xia_Seq_Sim3.mat';
    % 
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX016/Result_SpikeCount_FixWin_DX016_10uA_Xia_Exp1_Seq_Full_3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX016/Result_SpikeCount_FixWin_DX016_10uA_Xia_Exp1_Seq_Full_4.mat';
    % 
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX018/Result_SpikeCount_FixWin_DX018_5_10uA_Xia_ISI_SimSeq2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX018/Result_SpikeCount_FixWin_DX018_10uA_Xia_ISI_SimSeq1.mat';
    % 
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX019/Result_SpikeCount_FixWin_DX019_5_10uA_Xia_ISI_SimSeq1.mat';
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_10uA_Xia_ISI_10uA_SimSeq3.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX020/Result_SpikeCount_FixWin_DX020_5_10uA_Xia_ISI_SimSeq2.mat';
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX021/Result_SpikeCount_FixWin_DX021_5_10uA_Xia_ISI_SimSeq2.mat';
    
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX022/Result_SpikeCount_FixWin_DX022_10uA_Xia_ISI_10uA_SimSeq1.mat';
    % 
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX023/Result_SpikeCount_FixWin_DX023_10uA_Xia_ISI_10uA_SinSeq1.mat';
    % 
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX024/Result_SpikeCount_FixWin_DX024_10uA_Xia_ISI_10uA_SimSeq1.mat';
};


% 2. Plotting Aesthetics
line_width = 1.5;
marker_size = 4;
cap_size = 6;

% --- Paper Figure Settings ---
fig_width_cm  = 8.8;     % IEEE single-column width
fig_height_cm = 8.8;     % slightly shorter than square
axis_fontsize = 9;
label_fontsize = 9;
title_fontsize = 10;
legend_fontsize = 8;
axis_linewidth = 1.0;
font_name = 'Arial';

% 3. Statistical Settings
stats_target_Amp = 5; 
friedman_target_ISIs = [0 5 9 10 12 15];

% 4. Normalization Settings
% Each set will be divided by its own maximum response across all
% available amplitudes and ISIs (including ISI = 0 sim condition).
min_valid_denominator = 0;

% 5. Figure Saving
save_figures = false;
save_dir = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure5/ISI/ISI_SingeAnimal';
save_name = 'ISI_DX021_Average_MaxNorm';
save_dpi = 600;

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

Pop_Data = nan(length(Union_Amps), length(Union_ISIs), max_Sets, nFiles);

%% ================= 2. GATHER & NORMALIZE (The Core Engine) =============
fprintf('\nExtracting and Normalizing Data...\n');
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
    
    % ========================================================
    % MODIFIED: Within-set max normalization
    % Each set is divided by its own maximum response across
    % all amplitudes and all ISIs available in that set.
    % ========================================================
    for ss = 1:nSets_this
        this_set_matrix = animal_means(:, :, ss);
        this_set_max = max(this_set_matrix(:), [], 'omitnan');
        
        if ~isnan(this_set_max) && this_set_max > min_valid_denominator
            Pop_Data(:, :, ss, k) = this_set_matrix / this_set_max;
        else
            fprintf('    -> Warning: Set %d has invalid or 0 maximum response. Skipping normalization for this set.\n', ss);
        end
    end
end

%% =================== 3. POPULATION STATISTICS =================
fprintf('\nCalculating Population Mean and SEM...\n');

% First average across sets within each file
Pop_Data_Pooled = nanmean(Pop_Data, 3);

% Then average across files
Pop_Mean = nanmean(Pop_Data_Pooled, 4);
Pop_N    = sum(~isnan(Pop_Data_Pooled), 4);
Pop_Std  = nanstd(Pop_Data_Pooled, 0, 4);
% Pop_SEM  = Pop_Std ./ sqrt(Pop_N);
% Pop_SEM  = Pop_Std ./ sqrt(5); 
Pop_SEM  = Pop_Std./sqrt(1);

%% =================== 3.5 STATISTICAL TEST (REDUCED FRIEDMAN) =================
target_a_idx = find(abs(Union_Amps - stats_target_Amp) < 0.001);

if ~isempty(target_a_idx)
    fprintf('\n--------------------------------------------------\n');
    fprintf('Running Reduced Friedman Test for %.1f uA\n', stats_target_Amp);
    fprintf('Matched unit: File x Set\n');
    fprintf('Requested ISIs: %s ms\n', num2str(friedman_target_ISIs));
    fprintf('--------------------------------------------------\n');
    
    % --------------------------------------------------
    % Map requested ISIs onto the Union_ISIs index
    % --------------------------------------------------
    friedman_idx = [];
    friedman_ISIs_found = [];
    
    for ii = 1:length(friedman_target_ISIs)
        idx = find(abs(Union_ISIs - friedman_target_ISIs(ii)) < 0.001);
        if ~isempty(idx)
            friedman_idx(end+1) = idx; %#ok<SAGROW>
            friedman_ISIs_found(end+1) = Union_ISIs(idx); %#ok<SAGROW>
        else
            fprintf('  Warning: ISI %.1f ms not found in Union_ISIs. Skipped.\n', friedman_target_ISIs(ii));
        end
    end
    
    if length(friedman_idx) < 2
        fprintf('Not enough ISIs found for Friedman test.\n');
    else
        fprintf('Using reduced ISIs: %s ms\n', num2str(friedman_ISIs_found));
        
        % ========================================================
        % Build matched matrix:
        %   rows    = File x Set
        %   columns = reduced ISIs only
        % ========================================================
        Friedman_Matrix = [];
        Friedman_FileID = [];
        Friedman_SetID  = [];
        
        for k = 1:nFiles
            this_file_data = Pop_Data(target_a_idx, :, :, k);   % [1 x ISI x Set]
            this_file_data = squeeze(this_file_data);           % [ISI x Set] or [ISI x 1]
            
            if isempty(this_file_data)
                continue;
            end
            
            if isvector(this_file_data)
                this_file_data = this_file_data(:);
            end
            
            nSets_this = size(this_file_data, 2);
            
            for ss = 1:nSets_this
                this_row_full = this_file_data(:, ss)';                  % [1 x all ISIs]
                this_row_red  = this_row_full(friedman_idx);             % [1 x reduced ISIs]
                
                % Keep only complete matched rows for reduced ISI set
                if all(~isnan(this_row_red))
                    Friedman_Matrix = [Friedman_Matrix; this_row_red]; %#ok<SAGROW>
                    Friedman_FileID = [Friedman_FileID; k]; %#ok<SAGROW>
                    Friedman_SetID  = [Friedman_SetID; ss]; %#ok<SAGROW>
                end
            end
        end
        
        fprintf('Total matched File x Set units included: %d\n', size(Friedman_Matrix, 1));
        fprintf('Number of reduced ISIs compared: %d\n', size(Friedman_Matrix, 2));
        
        if size(Friedman_Matrix, 1) >= 2 && size(Friedman_Matrix, 2) >= 2
            % --------------------------------------------------
            % Overall Friedman test
            % --------------------------------------------------
            [p_friedman, tbl_friedman, stats_friedman] = friedman(Friedman_Matrix, 1, 'off');
            
            chi2_stat = tbl_friedman{2, 5};
            df_friedman = tbl_friedman{2, 3};
            
            fprintf('\nFriedman Test Result:\n');
            fprintf('  chi-square(%d) = %.4f, p = %s\n', ...
                df_friedman, chi2_stat, format_p(p_friedman));
            
            if p_friedman < 0.05
                fprintf('  -> Significant overall ISI effect at %.1f uA.\n', stats_target_Amp);
            else
                fprintf('  -> No significant overall ISI effect at %.1f uA.\n', stats_target_Amp);
            end
            
            % --------------------------------------------------
            % Post hoc Wilcoxon signed-rank tests vs ISI = 0 ms
            % --------------------------------------------------
            idx_isi0_local = find(abs(friedman_ISIs_found - 0) < 0.001);
            
            if ~isempty(idx_isi0_local)
                fprintf('\nPost hoc Wilcoxon signed-rank tests (each reduced ISI vs 0 ms):\n');
                fprintf('%-10s | %-10s | %-10s | %-10s\n', 'ISI (ms)', 'N', 'p-value', 'Sig');
                fprintf('----------------------------------------------------------\n');
                
                base_vals = Friedman_Matrix(:, idx_isi0_local);
                
                for p = 1:length(friedman_ISIs_found)
                    if p == idx_isi0_local
                        continue;
                    end
                    
                    compare_vals = Friedman_Matrix(:, p);
                    valid_pair = ~isnan(base_vals) & ~isnan(compare_vals);
                    
                    if sum(valid_pair) < 2
                        fprintf('%-10.1f | %-10d | %-10s | %-10s\n', ...
                            friedman_ISIs_found(p), sum(valid_pair), 'NaN', 'Skipped');
                        continue;
                    end
                    
                    try
                        p_pair = signrank(base_vals(valid_pair), compare_vals(valid_pair));
                    catch
                        p_pair = NaN;
                    end
                    
                    sig_str = '';
                    if ~isnan(p_pair)
                        if p_pair < 0.001
                            sig_str = '***';
                        elseif p_pair < 0.01
                            sig_str = '**';
                        elseif p_pair < 0.05
                            sig_str = '*';
                        else
                            sig_str = 'n.s.';
                        end
                    else
                        sig_str = 'NaN';
                    end
                    
                    fprintf('%-10.1f | %-10d | %-10s | %-10s\n', ...
                        friedman_ISIs_found(p), sum(valid_pair), format_p(p_pair), sig_str);
                end
            else
                fprintf('\nWarning: ISI = 0 ms not included in reduced ISI set, so post hoc comparisons vs 0 ms were skipped.\n');
            end
        else
            fprintf('Not enough matched data for reduced Friedman test.\n');
        end
    end
else
    fprintf('\nWarning: Target amplitude for stats (%.1f uA) not found in datasets.\n', stats_target_Amp);
end

% --------------------------------------------------
% Additional planned pairwise Wilcoxon signed-rank tests
% Focus on peak-region comparisons
% --------------------------------------------------
fprintf('\nPlanned pairwise Wilcoxon signed-rank tests:\n');
fprintf('%-12s | %-12s | %-10s | %-10s | %-10s\n', 'ISI 1 (ms)', 'ISI 2 (ms)', 'N', 'p-value', 'Sig');
fprintf('--------------------------------------------------------------------------\n');

planned_pairs = [
    9   5;
    10  5;
    9   12;
    10  12;
    9   15;
    10  15;
    9   10;
];

for ii = 1:size(planned_pairs, 1)
    isi_1 = planned_pairs(ii, 1);
    isi_2 = planned_pairs(ii, 2);
    
    idx1 = find(abs(friedman_ISIs_found - isi_1) < 0.001);
    idx2 = find(abs(friedman_ISIs_found - isi_2) < 0.001);
    
    if isempty(idx1) || isempty(idx2)
        fprintf('%-12.1f | %-12.1f | %-10s | %-10s | %-10s\n', ...
            isi_1, isi_2, 'NaN', 'NaN', 'Skipped');
        continue;
    end
    
    vals1 = Friedman_Matrix(:, idx1);
    vals2 = Friedman_Matrix(:, idx2);
    valid_pair = ~isnan(vals1) & ~isnan(vals2);
    
    if sum(valid_pair) < 2
        fprintf('%-12.1f | %-12.1f | %-10d | %-10s | %-10s\n', ...
            isi_1, isi_2, sum(valid_pair), 'NaN', 'Skipped');
        continue;
    end
    
    try
        p_pair = signrank(vals1(valid_pair), vals2(valid_pair));
    catch
        p_pair = NaN;
    end
    
    sig_str = '';
    if ~isnan(p_pair)
        if p_pair < 0.001
            sig_str = '***';
        elseif p_pair < 0.01
            sig_str = '**';
        elseif p_pair < 0.05
            sig_str = '*';
        else
            sig_str = 'n.s.';
        end
    else
        sig_str = 'NaN';
    end
    
    fprintf('%-12.1f | %-12.1f | %-10d | %-10s | %-10s\n', ...
        isi_1, isi_2, sum(valid_pair), format_p(p_pair), sig_str);
end

%% ===================== 4. PLOT (Pooled Population Dose-Response) ======================
num_amps = length(Union_Amps);

% --- MODIFIED: use black only, distinguish by line/marker style ---
colors = repmat([0 0 0], num_amps, 1);
line_styles = {'--','-'};     % 5 uA dashed, 10 uA solid
markers = {'o','o'};          % both circles
marker_faces = {'w','k'};     % 5 uA open, 10 uA filled

figure('Color','w', 'Units','centimeters', ...
    'Position',[5 5 fig_width_cm fig_height_cm]);
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
        mfc = marker_faces{mod(a-1, length(marker_faces))+1};
        lbl = sprintf('%.1f µA', current_amp);
        
        plot_x   = Union_ISIs(valid_idx);
        plot_y   = y_mean(valid_idx);
        plot_err = y_sem(valid_idx);
        
        % --- MODIFIED: standard error bars instead of shaded band ---
        errorbar(plot_x, plot_y, plot_err, ...
            'LineStyle', ls, 'Marker', mk, 'Color', col, ...
            'LineWidth', 1.2, 'MarkerFaceColor', mfc, ...
            'MarkerEdgeColor', col, ...
            'MarkerSize', marker_size, ...
            'CapSize', 5, ...
            'DisplayName', lbl);
    end
end

xlabel('Inter-Stimulus Interval (ms)', 'FontSize', label_fontsize, 'FontName', font_name);
ylabel('Normalized Spike Count Response', 'FontSize', label_fontsize, 'FontName', font_name);
% title('ISI Population Average', 'FontSize', title_fontsize, 'FontName', font_name, 'FontWeight', 'bold');

xticks([0 3 5 7 9 10 12 15 20]);

xlim([0 20]);

% --- MODIFIED: tighten upper limit so top does not look open ---
ylim([0.1 1.1]);
yticks(0.1:0.2:1.1);

axis square;

set(gca, 'FontSize', axis_fontsize, ...
         'FontName', font_name, ...
         'LineWidth', axis_linewidth, ...
         'TickDir', 'out', ...
         'XTickLabelRotation', 0);

box off;

lgd = legend('Location','southeast', 'Box','off', ...
    'FontSize', legend_fontsize, 'FontName', font_name);
% title(lgd, 'Amplitudes');
lgd.Title.FontSize = legend_fontsize;
lgd.Title.FontName = font_name;

fprintf('Figure generated.\n');


if save_figures
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    exportgraphics(gcf, fullfile(save_dir, [save_name '.tiff']), ...
        'Resolution', save_dpi);
    
    fprintf('Figure saved to: %s\n', fullfile(save_dir, [save_name '.tiff']));
end



%% =================== Functions =================
function p_str = format_p(p)
    if isnan(p)
        p_str = 'NaN';
    else
        p_str = sprintf('%.20f', p);
    end
end

function plot_shaded_error(x, y, se, col)
    if numel(x) < 2
        return;
    end
    
    x  = x(:)';
    y  = y(:)';
    se = se(:)';
    
    valid = ~isnan(x) & ~isnan(y) & ~isnan(se);
    x = x(valid);
    y = y(valid);
    se = se(valid);
    
    if isempty(x)
        return;
    end
    
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end