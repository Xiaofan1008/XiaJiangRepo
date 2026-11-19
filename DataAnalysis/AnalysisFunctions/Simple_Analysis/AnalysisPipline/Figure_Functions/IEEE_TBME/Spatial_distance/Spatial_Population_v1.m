%% ============================================================
%   POPULATION SPATIAL ANALYSIS: Pooled Recruitment & Reach (100um)
%   - Logic: 
%       1. Dynamic X-mapping for 100um input matrices.
%       2. Profile Plots: Global, Inner, and Outer with Shaded SEM.
%       3. Stats: T-test comparisons per Amplitude.
%   - Style: IEEE TBME Publication Style (B&W)
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
file_paths = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX013_Xia_Exp1_Seq_Sim8.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX014_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX014_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX014_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX014_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX014_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX014_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX015_Xia_Seq_Sim1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX015_Xia_Seq_Sim2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX015_Xia_Seq_Sim3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX015_Xia_Seq_Sim4.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX015_Xia_Seq_Sim5.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX015_Xia_Seq_Sim6.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX015_Xia_Seq_Sim7.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX016_Xia_Exp1_Seq_Full_1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX016_Xia_Exp1_Seq_Full_2.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX016_Xia_Exp1_Seq_Full_3.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank/Result_Spatial_SameShank_DX016_Xia_Exp1_Seq_Full_4.mat';
};

save_dir = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure4/Spatial_analysis';

%% =================== 1. AGGREGATE POOLED DATA =================
fprintf('Pooling data from %d datasets...\n', length(file_paths));
Pooled = struct(); unique_amps = [];

for f = 1:length(file_paths)
    if ~exist(file_paths{f}, 'file'), continue; end
    D = load(file_paths{f});
    for ss = 1:length(D.SpatialResults.Set)
        for ai = 1:length(D.SpatialResults.Set(ss).Amp)
            Data = D.SpatialResults.Set(ss).Amp(ai);
            if isempty(Data.Val) || Data.Val == 0, continue; end
            fName = sprintf('A_%.1f', Data.Val); fName = strrep(fName, '.', 'p');
            if ~isfield(Pooled, fName)
                Pooled.(fName).Val = Data.Val; unique_amps = [unique_amps, Data.Val];
                Pooled.(fName).Sim_Prob_Global = []; Pooled.(fName).Seq_Prob_Global = [];
                Pooled.(fName).Sim_Prob_Inner  = []; Pooled.(fName).Seq_Prob_Inner  = [];
                Pooled.(fName).Sim_Prob_Outer  = []; Pooled.(fName).Seq_Prob_Outer  = [];
                Pooled.(fName).Sim_Norm_Global = []; Pooled.(fName).Seq_Norm_Global = [];
                Pooled.(fName).Sim_TotalPerc   = []; Pooled.(fName).Seq_TotalPerc   = [];
                Pooled.(fName).Sim_R50         = []; Pooled.(fName).Seq_R50         = [];
            end
            Pooled.(fName).Sim_Prob_Global = [Pooled.(fName).Sim_Prob_Global; Data.Prob_Global_Sim];
            Pooled.(fName).Seq_Prob_Global = [Pooled.(fName).Seq_Prob_Global; Data.Prob_Global_Seq];
            Pooled.(fName).Sim_Prob_Inner  = [Pooled.(fName).Sim_Prob_Inner;  Data.Prob_Inner_Sim];
            Pooled.(fName).Seq_Prob_Inner  = [Pooled.(fName).Seq_Prob_Inner;  Data.Prob_Inner_Seq];
            Pooled.(fName).Sim_Prob_Outer  = [Pooled.(fName).Sim_Prob_Outer;  Data.Prob_Outer_Sim];
            Pooled.(fName).Seq_Prob_Outer  = [Pooled.(fName).Seq_Prob_Outer;  Data.Prob_Outer_Seq];
            Pooled.(fName).Sim_Norm_Global = [Pooled.(fName).Sim_Norm_Global; Data.Norm_Global_Sim];
            Pooled.(fName).Seq_Norm_Global = [Pooled.(fName).Seq_Norm_Global; Data.Norm_Global_Seq];
            Pooled.(fName).Sim_TotalPerc   = [Pooled.(fName).Sim_TotalPerc;   Data.TotalPerc_Sim];
            Pooled.(fName).Seq_TotalPerc   = [Pooled.(fName).Seq_TotalPerc;   Data.TotalPerc_Seq];
            Pooled.(fName).Sim_R50         = [Pooled.(fName).Sim_R50;         Data.R50_Sim];
            Pooled.(fName).Seq_R50         = [Pooled.(fName).Seq_R50;         Data.R50_Seq];
        end
    end
end
unique_amps = sort(unique(unique_amps));

%% ================= 2. GENERATE POPULATION PROFILES (Shaded SEM) =================
for i = 1:length(unique_amps)
    fName = sprintf('A_%.1f', unique_amps(i)); fName = strrep(fName, '.', 'p');
    P = Pooled.(fName);
    figure('Units', 'centimeters', 'Position', [2, 2, 26, 9], 'Color', 'w', 'Name', ['Profiles_' fName]);
    t = tiledlayout(1,3, 'TileSpacing', 'compact');
    
    zone_names = {'Global', 'Inner (200µm Gap)', 'Outer'};
    sim_data_sets = {P.Sim_Prob_Global, P.Sim_Prob_Inner, P.Sim_Prob_Outer};
    seq_data_sets = {P.Seq_Prob_Global, P.Seq_Prob_Inner, P.Seq_Prob_Outer};
    
    for z = 1:3
        nexttile; hold on;
        currSim = sim_data_sets{z}; currSeq = seq_data_sets{z};
        
        mSim = mean(currSim, 1, 'omitnan'); semSim = std(currSim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(currSim),1));
        mSeq = mean(currSeq, 1, 'omitnan'); semSeq = std(currSeq, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(currSeq),1));
        
        % --- MODIFICATION: Dynamic X-Axis to match 100um (8 bins) or 50um (16 bins) ---
        nBins = length(mSim);
        if z == 2, x_plot = linspace(0, 200, nBins); else, x_plot = linspace(0, 800, nBins); end
        
        % --- SHADED SEM ---
        fill([x_plot(:).', fliplr(x_plot(:).')], [mSim(:).'-semSim(:).', fliplr(mSim(:).'+semSim(:).')], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        fill([x_plot(:).', fliplr(x_plot(:).')], [mSeq(:).'-semSeq(:).', fliplr(mSeq(:).'+semSeq(:).')], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % --- MEAN LINES ---
        plot(x_plot, mSim, '--ok', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
        plot(x_plot, mSeq, '-sk', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
        
        title(zone_names{z}); ylim([0 1.1]);
        if z == 1, ylabel('Recruitment Probability'); end
        if z == 2, xlim([0 200]); else, xlim([0 800]); end
    end
    xlabel(t, 'Distance (µm)', 'FontName', 'Arial');
    legend('Location', 'northeast', 'Box', 'off');
end

%% ================= 3. PLOT POPULATION DENSITY (Normalized Shank) =================
figure('Units', 'centimeters', 'Position', [5, 5, 12, 10], 'Color', 'w', 'Name', 'Population_Density'); hold on;
% Plotting the density for the highest amplitude available
fHigh = sprintf('A_%.1f', max(unique_amps)); fHigh = strrep(fHigh, '.', 'p');
mSimNorm = mean(Pooled.(fHigh).Sim_Norm_Global, 1, 'omitnan');
mSeqNorm = mean(Pooled.(fHigh).Seq_Norm_Global, 1, 'omitnan');
x_norm = linspace(0, 800, length(mSimNorm));

plot(x_norm, mSimNorm, '--ok', 'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
plot(x_norm, mSeqNorm, '-sk', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
xlabel('Distance (µm)'); ylabel('Fraction of Total Shank Activity');
title(['Total Density @ ' num2str(max(unique_amps)) ' µA']); xlim([0 800]); ylim([0 0.5]); axis square; 

%% ================= 4. SUMMARY PLOTS WITH STATS =================
sum_sim_perc = []; sum_seq_perc = []; sum_sim_r50 = []; sum_seq_r50 = []; p_vals_perc = []; p_vals_r50 = [];

for i = 1:length(unique_amps)
    fName = sprintf('A_%.1f', unique_amps(i)); fName = strrep(fName, '.', 'p');
    P = Pooled.(fName);
    sum_sim_perc(i,:) = [mean(P.Sim_TotalPerc), std(P.Sim_TotalPerc)/sqrt(length(P.Sim_TotalPerc))];
    sum_seq_perc(i,:) = [mean(P.Seq_TotalPerc), std(P.Seq_TotalPerc)/sqrt(length(P.Seq_TotalPerc))];
    sum_sim_r50(i,:)  = [mean(P.Sim_R50, 'omitnan'), std(P.Sim_R50, 'omitnan')/sqrt(sum(~isnan(P.Sim_R50)))];
    sum_seq_r50(i,:)  = [mean(P.Seq_R50, 'omitnan'), std(P.Seq_R50, 'omitnan')/sqrt(sum(~isnan(P.Seq_R50)))];
    [~, p_vals_perc(i)] = ttest2(P.Sim_TotalPerc, P.Seq_TotalPerc);
    [~, p_vals_r50(i)]  = ttest2(P.Sim_R50, P.Seq_R50);
end

figure('Units', 'centimeters', 'Position', [5, 5, 20, 10], 'Color', 'w'); t = tiledlayout(1,2);
nexttile; hold on; % Active %
errorbar(unique_amps, sum_sim_perc(:,1), sum_sim_perc(:,2), '--ok', 'MarkerFaceColor','w', 'LineWidth', 1.5, 'CapSize', 8, 'DisplayName', 'Sim');
errorbar(unique_amps, sum_seq_perc(:,1), sum_seq_perc(:,2), '-sk', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 8, 'DisplayName', 'Seq');
ylabel('Active Channel %'); xlabel('Amplitude (µA)'); title('Active Percentage'); axis square; 
add_sig_stars(unique_amps, sum_sim_perc(:,1), sum_seq_perc(:,1), sum_seq_perc(:,2), p_vals_perc);

nexttile; hold on; % R50
errorbar(unique_amps, sum_sim_r50(:,1), sum_sim_r50(:,2), '--ok', 'MarkerFaceColor','w', 'LineWidth', 1.5, 'CapSize', 8);
errorbar(unique_amps, sum_seq_r50(:,1), sum_seq_r50(:,2), '-sk', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 8);
ylabel('R50 Distance (µm)'); xlabel('Amplitude (µA)'); title('R_{50}'); axis square; 
add_sig_stars(unique_amps, sum_sim_r50(:,1), sum_seq_r50(:,1), sum_seq_r50(:,2), p_vals_r50);

%% ================= 5. COMMAND WINDOW STATS =================
fprintf('\n====================================================================\n');
fprintf('         POPULATION STATISTICAL SUMMARY (N = %d Datasets)            \n', length(file_paths));
fprintf('====================================================================\n');
for i = 1:length(unique_amps)
    fprintf('%3.1fuA | Active%% Sim: %5.1f%%, Seq: %5.1f%% | R50 Seq: %5.1fum | p=%.3f\n', ...
        unique_amps(i), sum_sim_perc(i,1), sum_seq_perc(i,1), sum_seq_r50(i,1), p_vals_perc(i));
end

function add_sig_stars(x, y1, y2, y2_err, pvals)
    for i = 1:length(x)
        p = pvals(i); if isnan(p), continue; end
        txt = ''; if p < 0.001, txt = '***'; elseif p < 0.01, txt = '**'; elseif p < 0.05, txt = '*'; end
        if ~isempty(txt)
            y_max = max(y1(i), y2(i) + y2_err(i));
            text(x(i), y_max * 1.15, txt, 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
end