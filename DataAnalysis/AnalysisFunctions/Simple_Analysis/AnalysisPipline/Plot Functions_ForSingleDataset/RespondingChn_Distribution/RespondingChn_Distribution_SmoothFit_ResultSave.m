%% ============================================================
%   Spatial Analysis: Distribution (Bar), Focality (R80), & CPF
%   - Fig 1: Spatial Distribution Bar Plot (Target Amp)
%   - Fig 2: Focality Radius (R80) vs Amplitude (Mean +/- SEM)
%   - Fig 3: Cumulative Probability Function (CPF) Curves
%   - Stats: Pooled ANOVA on R80
%   - Output: Saves Results & Prints Excel-Ready Tables
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';
Electrode_Type = 1; 

target_amp = 10; % For Bar Plot and CPF
target_ptd_seq = 5; 
containment_threshold = 0.80; % R80

%% =================== LOAD DATA ====================
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% Extract Params
Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim;
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP; 
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq); for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v  = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

Coords = get_electrode_coordinates(Electrode_Type);
if Electrode_Type == 1, valid_channels = 1:32; else, valid_channels = [1:16, 33:48, 49:64, 17:32]; end
idx_ptd = find(abs(PTDs_ms - target_ptd_seq) < 0.1, 1); if isempty(idx_ptd), idx_ptd=1; end

%% ================= STEP 1: CALCULATE BAR PLOT DATA (Target Amp) =================
% A. SIMULTANEOUS
idx_sim_amp = find(Amps_sim == target_amp, 1);
stim_locs_sim  = Coords(Rsim.set(1).stimChannels, :); 
Sim_Results_Plot = []; 
chan_structs = Rsim.set(1).amp(idx_sim_amp).ptd(1).channel;
for ch = valid_channels
    rec_loc = Coords(ch, :);
    dists = sqrt(sum((stim_locs_sim - rec_loc).^2, 2));
    min_dist = round(min(dists));
    if min_dist < 1, continue; end 
    is_resp = 0; if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive, is_resp=1; end
    Sim_Results_Plot = [Sim_Results_Plot; min_dist, is_resp];
end

% B. SEQUENTIAL
idx_seq_amp = find(Amps_seq == target_amp, 1);
Seq_Results_Plot = cell(nSets_seq, 1);
for ss = 1:nSets_seq
    stim_chans_seq = uniqueComb_seq(ss, :); stim_chans_seq = stim_chans_seq(stim_chans_seq>0);
    stim_locs_seq  = Coords(stim_chans_seq, :);
    chan_structs = Rseq.set(ss).amp(idx_seq_amp).ptd(idx_ptd).channel;
    
    for ch = valid_channels
        rec_loc = Coords(ch, :);
        dists = sqrt(sum((stim_locs_seq - rec_loc).^2, 2));
        min_dist = round(min(dists));
        if min_dist < 1, continue; end
        is_resp = 0; if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive, is_resp=1; end
        Seq_Results_Plot{ss} = [Seq_Results_Plot{ss}; min_dist, is_resp];
    end
end

%% ================= STEP 2: CALCULATE R80 FOR ALL AMPLITUDES =================
R80_Results = []; % [Amp, Sim_R80, Seq_Mean, Seq_SEM, Set1, Set2...]
ANOVA_Data_Y = []; ANOVA_Data_X = []; ANOVA_Group = {}; 

for ai = 1:length(Amps_sim)
    amp_val = Amps_sim(ai);
    
    % --- Sim ---
    chan_structs = Rsim.set(1).amp(ai).ptd(1).channel;
    stim_locs = Coords(Rsim.set(1).stimChannels, :);
    resp_dists = [];
    for ch = valid_channels
        if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive
            rec_loc = Coords(ch, :);
            d = round(min(sqrt(sum((stim_locs - rec_loc).^2, 2))));
            resp_dists = [resp_dists; d];
        end
    end
    sim_r80 = NaN;
    if ~isempty(resp_dists)
        sorted_dists = sort(resp_dists);
        cutoff_idx = ceil(containment_threshold * length(sorted_dists));
        sim_r80 = sorted_dists(cutoff_idx);
        
        ANOVA_Data_Y = [ANOVA_Data_Y; sim_r80];
        ANOVA_Data_X = [ANOVA_Data_X; amp_val];
        ANOVA_Group = [ANOVA_Group; {'Simultaneous'}];
    end
    
    % --- Seq (Pooled) ---
    seq_r80_list = nan(1, nSets_seq);
    ai_seq = find(Amps_seq == amp_val);
    if ~isempty(ai_seq)
        for ss = 1:nSets_seq
            chan_structs = Rseq.set(ss).amp(ai_seq).ptd(idx_ptd).channel;
            stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
            stim_locs = Coords(stimCh, :);
            resp_dists = [];
            for ch = valid_channels
                if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive
                    rec_loc = Coords(ch, :);
                    d = round(min(sqrt(sum((stim_locs - rec_loc).^2, 2))));
                    resp_dists = [resp_dists; d];
                end
            end
            if ~isempty(resp_dists)
                sorted_dists = sort(resp_dists);
                cutoff_idx = ceil(containment_threshold * length(sorted_dists));
                val = sorted_dists(cutoff_idx);
                seq_r80_list(ss) = val;
                
                ANOVA_Data_Y = [ANOVA_Data_Y; val];
                ANOVA_Data_X = [ANOVA_Data_X; amp_val];
                ANOVA_Group = [ANOVA_Group; {'Sequential'}];
            end
        end
    end
    
    % Stats for Seq
    valid_seq = seq_r80_list(~isnan(seq_r80_list));
    if ~isempty(valid_seq)
        seq_avg = mean(valid_seq);
        seq_sem = std(valid_seq) / sqrt(length(valid_seq));
    else
        seq_avg = NaN; seq_sem = NaN;
    end
    
    % Store: [Amp, Sim, Seq_Mean, Seq_SEM, Set1, Set2...]
    row_data = [amp_val, sim_r80, seq_avg, seq_sem, seq_r80_list];
    R80_Results = [R80_Results; row_data];
end

%% ================= STEP 3: FIGURE 1 - BAR PLOT (DISTRIBUTION) =================
all_dists = Sim_Results_Plot(:, 1);
for ss=1:nSets_seq, all_dists = [all_dists; Seq_Results_Plot{ss}(:,1)]; end
Unique_Dists = unique(all_dists); Unique_Dists = Unique_Dists(Unique_Dists <= 600);
Plot_Data = zeros(length(Unique_Dists), 1 + nSets_seq);

Total_Sim = sum(Sim_Results_Plot(:,2));
for i = 1:length(Unique_Dists)
    d = Unique_Dists(i); mask = (Sim_Results_Plot(:,1) == d);
    if Total_Sim > 0, Plot_Data(i, 1) = sum(Sim_Results_Plot(mask, 2)) / Total_Sim * 100; end
end
for ss = 1:nSets_seq
    dat = Seq_Results_Plot{ss};
    Total_Seq = sum(dat(:,2));
    for i = 1:length(Unique_Dists)
        d = Unique_Dists(i); mask = (dat(:,1) == d);
        if Total_Seq > 0, Plot_Data(i, 1+ss) = sum(dat(mask, 2)) / Total_Seq * 100; end
    end
end

figure('Color','w', 'Position',[100 600 800 500]); hold on;
b = bar(Unique_Dists, Plot_Data, 'grouped');
prof_colors = [0.25 0.45 0.65; 0.80 0.45 0.30; 0.55 0.40 0.65];
for i = 1:length(b)
    if i <= size(prof_colors,1), b(i).FaceColor = prof_colors(i,:); end
    b(i).EdgeColor = 'none'; b(i).FaceAlpha = 0.85;
end
ylabel('% of Total Responding Population', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Distance (µm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Spatial Distribution (at %.0f µA)', target_amp), 'FontSize', 14);
legend('Simultaneous', 'Seq Set 1', 'Seq Set 2', 'Location','northeast','Box','off'); box off;

%% ================= STEP 4: FIGURE 2 - R80 vs AMPLITUDE (AVG) =================
figure('Color','w', 'Position',[100 100 600 500]); hold on;
col_sim = [0 0.3 0.8];
col_seq = [0.85 0.33 0.10]; 

% 1. Plot Sim Curve
plot(R80_Results(:,1), R80_Results(:,2), '-o', 'Color', col_sim, ...
    'LineWidth', 2, 'MarkerFaceColor', col_sim, 'DisplayName', 'Simultaneous');

% 2. Plot Seq Average Curve with Error Bars
errorbar(R80_Results(:,1), R80_Results(:,3), R80_Results(:,4), '-s', ...
    'Color', col_seq, 'LineWidth', 2, 'MarkerFaceColor', 'w', ...
    'CapSize', 0, 'DisplayName', 'Sequential (Avg)');

xlabel('Amplitude (µA)', 'FontWeight','bold'); 
ylabel(sprintf('Focality Radius (R_{%.0f}) [µm]', containment_threshold*100), 'FontWeight','bold');
title('Spatial Spread (R80) vs Amplitude', 'FontWeight','bold'); 
grid on; box off; legend('Location','best');

%% ================= STEP 5: FIGURE 3 - CPF CURVES (Target Amp) =================
% Cumulative Probability Function
figure('Color','w', 'Position',[750 100 600 500]); hold on;
dist_fine = 0:10:800;

% A. Sim CPF
Sim_CDF = zeros(size(dist_fine));
if Total_Sim > 0
    % Get raw distances of responders
    sim_resp_dists = Sim_Results_Plot(Sim_Results_Plot(:,2)==1, 1);
    for i = 1:length(dist_fine)
        Sim_CDF(i) = sum(sim_resp_dists <= dist_fine(i)) / Total_Sim;
    end
    plot(dist_fine, Sim_CDF, 'Color', col_sim, 'LineWidth', 3, 'DisplayName', 'Simultaneous');
end

% B. Seq CPF (Average across sets)
Seq_CDF_Matrix = [];
for ss = 1:nSets_seq
    dat = Seq_Results_Plot{ss};
    tot = sum(dat(:,2));
    if tot > 0
        seq_resp_dists = dat(dat(:,2)==1, 1);
        cdf_curve = zeros(size(dist_fine));
        for i = 1:length(dist_fine)
            cdf_curve(i) = sum(seq_resp_dists <= dist_fine(i)) / tot;
        end
        Seq_CDF_Matrix = [Seq_CDF_Matrix; cdf_curve];
    end
end

if ~isempty(Seq_CDF_Matrix)
    seq_cdf_mean = mean(Seq_CDF_Matrix, 1);
    seq_cdf_sem  = std(Seq_CDF_Matrix, 0, 1) / sqrt(size(Seq_CDF_Matrix, 1));
    
    % Plot Shaded Error
    fill([dist_fine, fliplr(dist_fine)], [seq_cdf_mean + seq_cdf_sem, fliplr(seq_cdf_mean - seq_cdf_sem)], ...
        col_seq, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(dist_fine, seq_cdf_mean, 'Color', col_seq, 'LineWidth', 3, 'DisplayName', 'Sequential (Avg)');
end

yline(containment_threshold, '--k', sprintf('R_{%.0f}', containment_threshold*100), 'HandleVisibility', 'off');
xlabel('Distance (µm)', 'FontWeight','bold'); 
ylabel('Cumulative Fraction of Activity', 'FontWeight','bold');
title(sprintf('Cumulative Spatial Distribution (at %.0f µA)', target_amp), 'FontWeight','bold'); 
grid on; box off; legend('Location','southeast'); xlim([0 600]);

%% ================= STEP 6: POOLED ANOVA (R80) =================
fprintf('\n--- POOLED ANOVA (R80 Comparison) ---\n');
if ~isempty(ANOVA_Data_Y)
    [p, tbl, stats] = anovan(ANOVA_Data_Y, {ANOVA_Group, ANOVA_Data_X}, ...
        'model', 'interaction', 'varnames', {'StimType', 'Amplitude'}, 'display', 'off');
    fprintf('StimType P-Value: %.5f\n', p(1));
    fprintf('Amplitude P-Value: %.5f\n', p(2));
    fprintf('Interaction P-Value: %.5f\n', p(3));
else
    fprintf('Insufficient data for ANOVA.\n');
end

%% ============================================================
%   7. OUTPUT RESULTS (COMMAND WINDOW)
% ============================================================
fprintf('\n\n=================================================================================\n');
fprintf('             FOCALITY RADIUS (R%d) DATA\n', containment_threshold*100);

fprintf('Amp(uA)\tSim_R80\tSeq_Mean\tSeq_SEM');
for ss=1:nSets_seq, fprintf('\tSet%d', ss); end
fprintf('\n');

for i = 1:size(R80_Results, 1)
    fprintf('%.0f\t%.2f\t%.2f\t%.2f', R80_Results(i,1), R80_Results(i,2), R80_Results(i,3), R80_Results(i,4));
    for k = 5:size(R80_Results, 2), fprintf('\t%.2f', R80_Results(i,k)); end
    fprintf('\n');
end
%% ============================================================
%   8. SAVE RESULTS
% ============================================================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Set1_Spatial_Full_5ms_' exp_id '.mat']);

ResultSpatial = struct();
ResultSpatial.R80.Table = R80_Results; 
ResultSpatial.R80.Stats.P = p;
ResultSpatial.BarData = Plot_Data;
ResultSpatial.CPF.Sim = Sim_CDF;
ResultSpatial.CPF.Seq = Seq_CDF_Matrix;
save(out_filename, 'ResultSpatial');
fprintf('Success! Results saved to:\n  %s\n', out_filename);

%% ==================== HELPERS ====================
function Coords = get_electrode_coordinates(type)
    Coords = zeros(64, 2); 
    if type == 1 
        for ch = 1:32, Coords(ch, 1) = 0; Coords(ch, 2) = (ch-1)*50; end
    elseif type == 2
        for i = 1:16, ch=i; Coords(ch,1)=0; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=32+i; Coords(ch,1)=200; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=48+i; Coords(ch,1)=400; Coords(ch,2)=(i-1)*50; end
        for i = 1:16, ch=16+i; Coords(ch,1)=600; Coords(ch,2)=(i-1)*50; end
    end
end
function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    S = load(dir('*_exp_datafile_*.mat').name);
    sp=[]; trig=[];
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
end