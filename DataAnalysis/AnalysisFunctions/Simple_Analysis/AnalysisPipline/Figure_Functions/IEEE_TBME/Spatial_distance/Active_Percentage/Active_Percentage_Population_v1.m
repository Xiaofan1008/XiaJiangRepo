%% ============================================================
%   POPULATION POTENCY ANALYSIS: 32-Channel Dual-Shank Cluster
%   - Logic: 
%       1. Pooling: Aggregates Active % from 32-ch dual-shank files.
%       2. N-Filter: Only plots amplitudes with N >= 5 datasets.
%       3. Save Choice: High-res TIFF export for Potency Figure.
%   - Style: IEEE TBME Publication Style (B&W)
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX005_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX006_Xia_Exp1_Seq4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX009_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX009_Xia_Exp1_Seq5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq2_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX012_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage/Result_ActivePercentage_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_figs = true; 
save_dir  = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Active_Percentage';

%% =================== 1. AGGREGATE POOLED DATA =================
fprintf('Pooling Potency data from %d datasets...\n', length(file_paths));
Pooled = struct(); unique_amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), continue; end
    D = load(file_paths{f});
    
    for ss = 1:length(D.PotencyResults.Set)
        for ai = 1:length(D.PotencyResults.Set(ss).Amp)
            Data = D.PotencyResults.Set(ss).Amp(ai);
            if isempty(Data.Val) || Data.Val > 10, continue; end 
            
            fName = sprintf('A_%.1f', Data.Val); fName = strrep(fName, '.', 'p');
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = Data.Val; unique_amps = [unique_amps, Data.Val];
                Pooled.(fName).Sim_TotalPerc = []; 
                Pooled.(fName).Seq_TotalPerc = [];
            end
            Pooled.(fName).Sim_TotalPerc = [Pooled.(fName).Sim_TotalPerc; Data.TotalPerc_Sim];
            Pooled.(fName).Seq_TotalPerc = [Pooled.(fName).Seq_TotalPerc; Data.TotalPerc_Seq];
        end
    end
end
unique_amps = sort(unique(unique_amps));

%% ================= 4. SUMMARY POTENCY TREND (N >= 5) =================
sum_amps = []; sum_sim_perc = []; sum_seq_perc = []; 
p_perc = []; n_counts = [];

for i = 1:length(unique_amps)
    fName = sprintf('A_%.1f', unique_amps(i)); fName = strrep(fName, '.', 'p');
    P = Pooled.(fName);
    
    % --- N-FILTER ---
    if length(P.Sim_TotalPerc) < 5, continue; end
    
    sum_amps = [sum_amps; P.Val];
    n_counts = [n_counts; length(P.Sim_TotalPerc)];
    
    % Calculate Mean and SEM for N=32 cluster
    sum_sim_perc = [sum_sim_perc; mean(P.Sim_TotalPerc), std(P.Sim_TotalPerc)/sqrt(length(P.Sim_TotalPerc))];
    sum_seq_perc = [sum_seq_perc; mean(P.Seq_TotalPerc), std(P.Seq_TotalPerc)/sqrt(length(P.Seq_TotalPerc))];
    
    [~, p_perc(end+1)] = ttest2(P.Sim_TotalPerc, P.Seq_TotalPerc);
end

figure('Color','w', 'Units', 'centimeters', 'Position', [15, 5, 8.8, 8.8], 'Name', 'Active Percentage'); hold on;
hold on;
errorbar(sum_amps, sum_sim_perc(:,1), sum_sim_perc(:,2), '--ok', 'MarkerFaceColor','w', 'LineWidth', 1, 'CapSize', 8, 'DisplayName', 'Sim');
errorbar(sum_amps, sum_seq_perc(:,1), sum_seq_perc(:,2), '-sk', 'MarkerFaceColor','k', 'LineWidth', 1, 'CapSize', 8, 'DisplayName', 'Seq');

ylabel('Active Channel Percentage %','FontSize', 9, 'FontName', 'Arial'); 
xlabel('Amplitude (µA)','FontSize', 9, 'FontName', 'Arial'); 
% title('Active Percentage'); 
axis square; 
% grid on;
add_sig_stars(sum_amps, sum_sim_perc(:,1), sum_seq_perc(:,1), sum_seq_perc(:,2), p_perc);
legend('Location', 'northwest', 'Box', 'off');

%% ================= 5. POTENCY STATISTICAL TABLE =================
fprintf('\n============================================================================\n');
fprintf('       MASTER POTENCY SUMMARY (Dual-Shank Cluster N=32, N >= 5 Rats)        \n');
fprintf('============================================================================\n');
fprintf('%-6s | %-3s | %-10s | %-10s | %-10s | %-6s\n', 'Amp', 'N', 'Mean Sim', 'Mean Seq', 'P-value', 'Delta');
fprintf('----------------------------------------------------------------------------\n');
for i = 1:length(sum_amps)
    delta = sum_seq_perc(i,1) - sum_sim_perc(i,1);
    fprintf('%3.1fuA  | %2d  | %5.1f%%    | %5.1f%%    | p=%.4f   | +%4.1f%%\n', ...
        sum_amps(i), n_counts(i), sum_sim_perc(i,1), sum_seq_perc(i,1), p_perc(i), delta);
end

%% ================= 6. SAVE FIGURES (TIFF) =================
if save_figs
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    figHandles = findall(0, 'Type', 'figure');
    for i = 1:length(figHandles)
        f = figHandles(i);
        fName = get(f, 'Name'); if isempty(fName), fName = ['Fig_' num2str(i)]; end
        fName = strrep(fName, ' ', '_'); 
        exportgraphics(f, fullfile(save_dir, [fName '.tiff']), 'Resolution', 300);
        fprintf('Saved: %s.tiff at 300 DPI\n', fName);
    end
end

function add_sig_stars(x, y1, y2, y2_err, pvals)
    for i = 1:length(x)
        p = pvals(i); if isnan(p), continue; end
        txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        if ~isempty(txt)
            y_max = max(y1(i), y2(i) + y2_err(i));
            text(x(i), y_max * 1.10, txt, 'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
end