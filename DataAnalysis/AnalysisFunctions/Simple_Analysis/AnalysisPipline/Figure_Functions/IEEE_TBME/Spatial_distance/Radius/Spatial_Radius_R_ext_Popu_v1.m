%% ============================================================
%   POPULATION SPATIAL ANALYSIS: R_ext Tuning Curve
%   - Input: Result_Spatial_RawRadius_*.mat files
%   - Logic: 
%       1. Dynamic Radius: Calculates R_ext (Maximal Responding Distance) 
%          from raw saved arrays for every individual stimulation set.
%       2. Statistics: Paired Wilcoxon signed-rank test (Sim vs. Seq).
%       3. Plotting: Amplitude vs. R_ext (Mean ± SEM) with significance markers.
%       4. Style: IEEE Black & White formatting.
% ============================================================
clear;

%% ================= USER SETTINGS =================
% Paste your specific list of saved .mat files here
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX005_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX006_Xia_Exp1_Seq4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX009_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX009_Xia_Exp1_Seq5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq2_5ms.mat';
    % % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq5_5ms.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq6_5ms.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq7.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq8.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX012_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_2.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius/Result_Spatial_RawRadius_DX016_Xia_Exp1_Seq_Full_4.mat';


    };

save_dir     = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Spatial_Radius';
save_figures = true; % Set to true to export .tiff files
tiff_dpi     = 600;  % High resolution for publication

%% ================= 1. DYNAMIC AMPLITUDE HARVESTING =================
fprintf('Scanning explicitly provided .mat files...\n');
num_files = length(file_paths);

if num_files == 0
    error('The file_paths cell array is empty.');
end

% Find all unique amplitudes across the explicitly listed files
all_amps_collected = [];
for f = 1:num_files
    tmp = load(file_paths{f}, 'Amps');
    all_amps_collected = [all_amps_collected; tmp.Amps(:)];
end
MasterAmps = unique(all_amps_collected);
MasterAmps(MasterAmps == 0) = []; % Remove baseline if it exists

fprintf('Loaded %d files. Identified Amplitudes: ', num_files);
fprintf('%.1f ', MasterAmps); fprintf('uA\n');

%% ================= 2. MASTER BUCKET POOLING (N = Sets) =================
% Preallocate cell arrays to hold the R_ext values for every set at each amplitude
Master_Rext_Sim = cell(length(MasterAmps), 1);
Master_Rext_Seq = cell(length(MasterAmps), 1);

total_sets_pooled = 0;

for f = 1:num_files
    data = load(file_paths{f});
    nSets = length(data.SpatialResults.Set);
    
    for ss = 1:nSets
        total_sets_pooled = total_sets_pooled + 1;
        
        for ai = 1:length(data.SpatialResults.Set(ss).Amp)
            D = data.SpatialResults.Set(ss).Amp(ai);
            if isempty(D.Val) || D.Val == 0, continue; end
            
            % Find where this amplitude belongs in the Master list
            m_idx = find(abs(MasterAmps - D.Val) < 0.001);
            if isempty(m_idx), continue; end
            
            % --- CALCULATION LOGIC: R_ext (Maximal Responding Distance) ---
            % We calculate this on the fly from the raw cleaned arrays
            dist = D.Dist_to_Boundary;
            
            % Find max distance where channel was active (Sim)
            act_sim = find(D.Clean_Sim > 0);
            if isempty(act_sim), r_sim = 0; else, r_sim = max(dist(act_sim)); end
            
            % Find max distance where channel was active (Seq)
            act_seq = find(D.Clean_Seq > 0);
            if isempty(act_seq), r_seq = 0; else, r_seq = max(dist(act_seq)); end
            
            % SAFETY: If no channels fired at all (potency = 0%), we treat as NaN 
            % to avoid dragging the average to zero if the electrode was non-functional.
            if D.TotalPerc_Sim == 0 && D.TotalPerc_Seq == 0
                r_sim = NaN; r_seq = NaN;
            end
            
            % Append to Master Buckets
            Master_Rext_Sim{m_idx} = [Master_Rext_Sim{m_idx}; r_sim];
            Master_Rext_Seq{m_idx} = [Master_Rext_Seq{m_idx}; r_seq];
        end
    end
end
fprintf('Total unique stimulation sets pooled (N): %d\n', total_sets_pooled);

%% ================= 3. POPULATION MATH & STATISTICS =================
PopResults = struct();

for a = 1:length(MasterAmps)
    PopResults.Amp(a).Val = MasterAmps(a);
    
    r_sim = Master_Rext_Sim{a};
    r_seq = Master_Rext_Seq{a};
    
    if isempty(r_sim), continue; end
    
    % Valid N-count (ignoring NaNs)
    N_count = sum(~isnan(r_sim) & ~isnan(r_seq));
    PopResults.Amp(a).N = N_count;
    
    % --- MEANS & SEM (OmitNaN) ---
    PopResults.Amp(a).Sim_Mean = mean(r_sim, 'omitnan');
    PopResults.Amp(a).Seq_Mean = mean(r_seq, 'omitnan');
    
    PopResults.Amp(a).Sim_SEM  = std(r_sim, 'omitnan') / sqrt(sum(~isnan(r_sim)));
    PopResults.Amp(a).Seq_SEM  = std(r_seq, 'omitnan') / sqrt(sum(~isnan(r_seq)));
    
    % --- PAIRED STATISTICS (Wilcoxon Signed-Rank) ---
    try 
        [pval, ~] = signrank(r_sim, r_seq); 
    catch 
        pval = NaN; 
    end
    PopResults.Amp(a).pval = pval;
end

%% ================= 4. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     POPULATION STATISTICS: R_ext (Maximal Responding Distance)     \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-4s | %-15s | %-15s | %-10s | %-6s\n', 'Amp', 'N', 'R_ext Sim (um)', 'R_ext Seq (um)', 'p-value', 'Sig');
fprintf('--------------------------------------------------------------------\n');
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N == 0, continue; end
    
    % sig_str = ''; 
    % if P.pval < 0.05, sig_str = '*'; end
    % if P.pval < 0.01, sig_str = '**'; end
    % if P.pval < 0.001, sig_str = '***'; end

    % Multi-level star logic for Command Window
    sig_str = ''; 
    if P.pval < 0.001, sig_str = '***';
    elseif P.pval < 0.01, sig_str = '**';
    elseif P.pval < 0.05, sig_str = '*'; end
    
    fprintf('%4.1fuA | %-4d | %5.1f ± %-5.1f | %5.1f ± %-5.1f | %-10.10f | %s\n', ...
        P.Val, P.N, P.Sim_Mean, P.Sim_SEM, P.Seq_Mean, P.Seq_SEM, P.pval, sig_str);
end
fprintf('\n');

%% ================= 5. PLOTTING: R_ext TUNING CURVE =================
% fig = figure('Units', 'centimeters', 'Position', [5, 5, 12, 11], 'Color', 'w', 'Name', 'Population R_ext');
fig = figure('Units', 'centimeters', 'Position', [2, 2, 8.8, 8.8], 'Color', 'w', 'Name', 'Population R_ext');
hold on; 

% Extract valid data for plotting
v_amps = []; sim_m = []; sim_s = []; seq_m = []; seq_s = []; p_list = [];
for a = 1:length(MasterAmps)
    P = PopResults.Amp(a);
    if isempty(P.N) || P.N < 2, continue; end
    v_amps = [v_amps; P.Val];
    sim_m = [sim_m; P.Sim_Mean]; sim_s = [sim_s; P.Sim_SEM];
    seq_m = [seq_m; P.Seq_Mean]; seq_s = [seq_s; P.Seq_SEM];
    p_list = [p_list; P.pval];
end

% --- PLOT LINES WITH ERROR BARS ---
% Sim: Dashed, White Marker
errorbar(v_amps, sim_m, sim_s, '--ok', 'LineWidth', 1, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Simultaneous', 'CapSize', 8);

% Seq: Solid, Black Marker
errorbar(v_amps, seq_m, seq_s, '-sk', 'LineWidth', 1, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'k', 'DisplayName', 'Sequential', 'CapSize', 8);

% --- ADD SIGNIFICANCE ASTERISKS ---
% for i = 1:length(v_amps)
%     if p_list(i) < 0.05
%         % Place asterisk above the highest error bar
%         y_pos = max([sim_m(i)+sim_s(i), seq_m(i)+seq_s(i)]) + 30; 
%         text(v_amps(i), y_pos, '*', 'FontSize', 20, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
%     end
% end

for i = 1:length(v_amps)
    sig_str = '';
    if p_list(i) < 0.001, sig_str = '***';
    elseif p_list(i) < 0.01, sig_str = '**';
    elseif p_list(i) < 0.05, sig_str = '*';
    end
    
    if ~isempty(sig_str)
        % Place stars 10% higher than the top of the highest error bar
        y_max_local = max([sim_m(i)+sim_s(i), seq_m(i)+seq_s(i)]);
        % y_pos = y_max_local + (max(sim_m)*0.05); 
        y_pos = max([sim_m(i)+sim_s(i), seq_m(i)+seq_s(i)]) + (max(sim_m)*0.08); 
        text(v_amps(i), y_pos, sig_str, 'FontSize', 10, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontName', 'Arial');
    end
end

% --- FORMATTING ---
% xlabel('Current Amplitude (\muA)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Maximal Responding Distance R_{ext} (\mum)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Effective Spatial Containment', 'FontSize', 13);
% set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 11);
% legend('Location', 'northwest', 'Box', 'off');
% xlim([0, max(v_amps)+1]); ylim([0, max([sim_m; seq_m])+100]);
% axis square;

% --- IEEE FORMATTING ---
xlabel('Amplitude (µA)', 'FontSize', 9, 'FontName', 'Arial');
ylabel('Maximal Activation Spread (µm)', 'FontSize', 9, 'FontName', 'Arial'); 

set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 1, 'FontSize', 9, 'FontName', 'Arial');
lgd = legend('Location', 'northwest', 'Box', 'off', 'FontSize', 9, 'FontName', 'Arial');
lgd.Position(2) = 0.82;
lgd.Position(1) = 0.15;
% xlim([0, max(v_amps)+1]); 
% ylim([0, 550]);
xlim([0, 10]); 
ylim([0, 450]);
axis square;

% --- SAVE TIFF ---
if save_figures
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    filename = fullfile(save_dir, 'Population_R_ext_TuningCurve.tiff');
    exportgraphics(fig, filename, 'Resolution', tiff_dpi);
    fprintf('>>> Population Figure saved to: %s\n', filename);
end