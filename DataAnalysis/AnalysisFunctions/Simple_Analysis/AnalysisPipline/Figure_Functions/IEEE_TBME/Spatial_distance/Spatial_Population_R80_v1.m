%% ============================================================
%   POPULATION SPATIAL ANALYSIS: Pooled Recruitment & R80 Reach
%   - Logic: 
%       1. HD 50um Mapping: Handles high-res input matrices.
%       2. N-Filter: Only plots amplitudes with N >= 5 datasets.
%       3. Stats: T-test comparisons for Potency and Reach (R80).
%   - Style: IEEE TBME Publication Style (B&W)
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% Update these paths to point to your new v3_R80 result files
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX010_Xia_Exp1_Seq1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX010_Xia_Exp1_Seq2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX010_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX010_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX010_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX010_Xia_Exp1_Seq7_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX010_Xia_Exp1_Seq8_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq2_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq4_5ms_new.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq5_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX011_Xia_Exp1_Seq9.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX012_Xia_Exp1_Seq1_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX012_Xia_Exp1_Seq4_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX012_Xia_Exp1_Seq6_5ms.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2/Result_Spatial_R80_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_dir = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Spatial_analysis';

%% =================== 1. AGGREGATE POOLED DATA =================
fprintf('Pooling data from %d datasets...\n', length(file_paths));
Pooled = struct(); unique_amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), continue; end
    D = load(file_paths{f});
    % Supports both single-struct and Sim/Seq split formats
    for ss = 1:length(D.SpatialResults.Set)
        for ai = 1:length(D.SpatialResults.Set(ss).Amp)
            Data = D.SpatialResults.Set(ss).Amp(ai);
            if isempty(Data.Val) || Data.Val > 10, continue; end % Cap at 10uA
            
            fName = sprintf('A_%.1f', Data.Val); fName = strrep(fName, '.', 'p');
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = Data.Val; unique_amps = [unique_amps, Data.Val];
                Pooled.(fName).Sim_Prob_Global = []; Pooled.(fName).Seq_Prob_Global = [];
                Pooled.(fName).Sim_TotalPerc   = []; Pooled.(fName).Seq_TotalPerc   = [];
                Pooled.(fName).Sim_R80         = []; Pooled.(fName).Seq_R80         = [];
                Pooled.(fName).Sim_Norm_Global = []; Pooled.(fName).Seq_Norm_Global = [];
            end
            
            % Append HD Profiles and R80 Reach
            Pooled.(fName).Sim_Prob_Global = [Pooled.(fName).Sim_Prob_Global; Data.Prob_Global_Sim];
            Pooled.(fName).Seq_Prob_Global = [Pooled.(fName).Seq_Prob_Global; Data.Prob_Global_Seq];
            Pooled.(fName).Sim_TotalPerc   = [Pooled.(fName).Sim_TotalPerc;   Data.TotalPerc_Sim];
            Pooled.(fName).Seq_TotalPerc   = [Pooled.(fName).Seq_TotalPerc;   Data.TotalPerc_Seq];
            Pooled.(fName).Sim_R80         = [Pooled.(fName).Sim_R80;         Data.R80_Sim];
            Pooled.(fName).Seq_R80         = [Pooled.(fName).Seq_R80;         Data.R80_Seq];
            Pooled.(fName).Sim_Norm_Global = [Pooled.(fName).Sim_Norm_Global; Data.Norm_Global_Sim];
            Pooled.(fName).Seq_Norm_Global = [Pooled.(fName).Seq_Norm_Global; Data.Norm_Global_Seq];
        end
    end
end
unique_amps = sort(unique(unique_amps));

%% ================= 2. GENERATE POPULATION PROFILES =================
% (Section 2 remains functionally similar, but now utilizes 50um resolution)
for i = 1:length(unique_amps)
    fName = sprintf('A_%.1f', unique_amps(i)); fName = strrep(fName, '.', 'p');
    P = Pooled.(fName);
    if size(P.Sim_Prob_Global, 1) < 5, continue; end % Only plot if N >= 5
    
    figure('Units', 'centimeters', 'Position', [2, 2, 10, 8], 'Color', 'w', 'Name', ['Profile_' fName]);
    hold on;
    mSim = mean(P.Sim_Prob_Global, 1, 'omitnan'); semSim = std(P.Sim_Prob_Global, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(P.Sim_Prob_Global),1));
    mSeq = mean(P.Seq_Prob_Global, 1, 'omitnan'); semSeq = std(P.Seq_Prob_Global, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(P.Seq_Prob_Global),1));
    x_plot = linspace(0, 800, length(mSim));
    
    fill([x_plot(:).', fliplr(x_plot(:).')], [mSim(:).'-semSim(:).', fliplr(mSim(:).'+semSim(:).')], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill([x_plot(:).', fliplr(x_plot(:).')], [mSeq(:).'-semSeq(:).', fliplr(mSeq(:).'+semSeq(:).')], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x_plot, mSim, '--ok', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
    plot(x_plot, mSeq, '-sk', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
    title(['Global Profile @ ' num2str(unique_amps(i)) ' \muA (N=' num2str(size(P.Sim_Prob_Global,1)) ')']);
    xlabel('Distance (\mum)'); ylabel('Probability'); grid on; xlim([0 800]); ylim([0 1.1]);
end

%% ================= 4. SUMMARY TREND PLOTS (N >= 5) =================
sum_amps = []; sum_sim_perc = []; sum_seq_perc = []; sum_sim_r80 = []; sum_seq_r80 = [];
p_perc = []; p_r80 = []; n_counts = [];

for i = 1:length(unique_amps)
    fName = sprintf('A_%.1f', unique_amps(i)); fName = strrep(fName, '.', 'p');
    P = Pooled.(fName);
    
    % --- CRITICAL N-FILTER ---
    if length(P.Sim_TotalPerc) < 5, continue; end
    
    sum_amps = [sum_amps; P.Val];
    n_counts = [n_counts; length(P.Sim_TotalPerc)];
    
    % Potency Stats
    sum_sim_perc = [sum_sim_perc; mean(P.Sim_TotalPerc), std(P.Sim_TotalPerc)/sqrt(length(P.Sim_TotalPerc))];
    sum_seq_perc = [sum_seq_perc; mean(P.Seq_TotalPerc), std(P.Seq_TotalPerc)/sqrt(length(P.Seq_TotalPerc))];
    [~, p_perc(end+1)] = ttest2(P.Sim_TotalPerc, P.Seq_TotalPerc);
    
    % Reach Stats (R80)
    sum_sim_r80 = [sum_sim_r80; mean(P.Sim_R80, 'omitnan'), std(P.Sim_R80, 'omitnan')/sqrt(sum(~isnan(P.Sim_R80)))];
    sum_seq_r80 = [sum_seq_r80; mean(P.Seq_R80, 'omitnan'), std(P.Seq_R80, 'omitnan')/sqrt(sum(~isnan(P.Seq_R80)))];
    [~, p_r80(end+1)] = ttest2(P.Sim_R80, P.Seq_R80);
end

figure('Units', 'centimeters', 'Position', [5, 5, 20, 10], 'Color', 'w'); t = tiledlayout(1,2);
% Panel A: Potency
nexttile; hold on;
errorbar(sum_amps, sum_sim_perc(:,1), sum_sim_perc(:,2), '--ok', 'MarkerFaceColor','w', 'LineWidth', 1.5, 'CapSize', 8, 'DisplayName', 'Sim');
errorbar(sum_amps, sum_seq_perc(:,1), sum_seq_perc(:,2), '-sk', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 8, 'DisplayName', 'Seq');
ylabel('Active Channel %'); xlabel('Amplitude (\muA)'); title('Network Potency'); axis square; grid on;
add_sig_stars(sum_amps, sum_sim_perc(:,1), sum_seq_perc(:,1), sum_seq_perc(:,2), p_perc);
legend('Location', 'northwest', 'Box', 'off');

% Panel B: Reach (R80)
nexttile; hold on;
errorbar(sum_amps, sum_sim_r80(:,1), sum_sim_r80(:,2), '--ok', 'MarkerFaceColor','w', 'LineWidth', 1.5, 'CapSize', 8);
errorbar(sum_amps, sum_seq_r80(:,1), sum_seq_r80(:,2), '-sk', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 8);
ylabel('Effective Reach (R_{80}, \mum)'); xlabel('Amplitude (\muA)'); title('Spatial Reach'); axis square; grid on;
add_sig_stars(sum_amps, sum_sim_r80(:,1), sum_seq_r80(:,1), sum_seq_r80(:,2), p_r80);

%% ================= 5. MASTER STATISTICAL TABLE =================
fprintf('\n=========================================================================================\n');
fprintf('       MASTER POPULATION SUMMARY (N >= 5 Filter applied, Capped at 10uA)                 \n');
fprintf('=========================================================================================\n');
fprintf('%-6s | %-4s | %-15s | %-15s | %-12s | %-10s\n', 'Amp', 'N', 'Active% (Sim/Seq)', 'R80 (Sim/Seq)', 'P (Potency)', 'Delta');
fprintf('-----------------------------------------------------------------------------------------\n');
for i = 1:length(sum_amps)
    delta = sum_seq_perc(i,1) - sum_sim_perc(i,1);
    fprintf('%3.1fuA  | %2d   | %4.1f%% vs %4.1f%% | %4.0fum vs %4.0fum | p=%.4f    | +%4.1f%%\n', ...
        sum_amps(i), n_counts(i), sum_sim_perc(i,1), sum_seq_perc(i,1), sum_sim_r80(i,1), sum_seq_r80(i,1), p_perc(i), delta);
end

function add_sig_stars(x, y1, y2, y2_err, pvals)
    for i = 1:length(x)
        p = pvals(i); if isnan(p), continue; end
        txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        if ~isempty(txt)
            y_max = max(y1(i), y2(i) + y2_err(i));
            text(x(i), y_max * 1.12, txt, 'FontSize', 12, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
end