%% ============================================================
%   Spatial Diversity Analysis (Refined: Gaussian Fit + Peak Plot)
%   - Metric: Spatial Probability of Response (N_resp / N_total)
%   - Figures: 
%       1. Grouped Bar Chart (Spatial Profile at Target Amp)
%       2. Gaussian Fit Curves (Spatial Profile at Target Amp)
%       3. Peak Probability vs Amplitude (Intensity Summary)
%   - Analysis: Sigma (Spread) & Peak_Prob (Intensity)
%   - Output: Saves Results & Prints Detailed Excel Tables
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));
%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim6_251125_181554';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq6_5ms_251125_182437';
Electrode_Type = 1; 
target_amp = 6; % For the spatial profile figures (Fig 1 & 2)
target_ptd_seq = 5; 
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
% Geometry
Coords = get_electrode_coordinates(Electrode_Type);
if Electrode_Type == 1, valid_channels = 1:32; else, valid_channels = [1:16, 33:48, 49:64, 17:32]; end
idx_ptd = find(abs(PTDs_ms - target_ptd_seq) < 0.1, 1); if isempty(idx_ptd), idx_ptd=1; end

%% ================= STEP 1: PREPARE DATA FOR PLOTTING (Target Amp) =================
% A. SIMULTANEOUS
idx_sim_amp = find(Amps_sim == target_amp, 1);
stim_chans_sim = Rsim.set(1).stimChannels; stim_locs_sim = Coords(stim_chans_sim, :); 
Sim_Results_Plot = []; 
chan_structs = Rsim.set(1).amp(idx_sim_amp).ptd(1).channel;
for ch = valid_channels
    rec_loc = Coords(ch, :); d = round(min(sqrt(sum((stim_locs_sim - rec_loc).^2, 2))));
    if d<1, continue; end; is_resp=0; if isfield(chan_structs(ch),'is_responsive')&&chan_structs(ch).is_responsive, is_resp=1; end
    Sim_Results_Plot = [Sim_Results_Plot; d, is_resp];
end
% B. SEQUENTIAL
idx_seq_amp = find(Amps_seq == target_amp, 1);
Seq_Results_Plot = cell(nSets_seq, 1);
for ss = 1:nSets_seq
    stim_chans_seq = uniqueComb_seq(ss,:); stim_chans_seq = stim_chans_seq(stim_chans_seq>0); stim_locs_seq = Coords(stim_chans_seq, :);
    chan_structs = Rseq.set(ss).amp(idx_seq_amp).ptd(idx_ptd).channel;
    for ch = valid_channels
        rec_loc = Coords(ch, :); d = round(min(sqrt(sum((stim_locs_seq - rec_loc).^2, 2))));
        if d<1, continue; end; is_resp=0; if isfield(chan_structs(ch),'is_responsive')&&chan_structs(ch).is_responsive, is_resp=1; end
        Seq_Results_Plot{ss} = [Seq_Results_Plot{ss}; d, is_resp];
    end
end

%% ================= STEP 2: CALCULATE GAUSSIAN FIT FOR ALL AMPS =================
% Model: P(x) = Peak * exp(-(x/Sigma)^2)
gauss_eq = @(b,x) b(1) .* exp(-(x./b(2)).^2);
Sigma_Results = []; 
Peak_Results  = []; 
warning('off', 'stats:nlinfit:ModelConstantWRTParam'); warning('off', 'stats:nlinfit:IllConditionedJacobian');

for ai = 1:length(Amps_sim)
    amp_val = Amps_sim(ai);
    beta0 = [0.5, 100]; 
    
    % --- Fit Sim ---
    [X, Y] = get_XY_data(Rsim, Coords, 1, ai, 1, valid_channels, stim_locs_sim);
    try
        beta_sim = nlinfit(X, Y, gauss_eq, beta0);
        P_peak_sim = min(max(beta_sim(1), 0), 1); 
        Sigma_sim  = min(max(abs(beta_sim(2)), 0), 2000); 
    catch
        P_peak_sim = NaN; Sigma_sim = NaN;
    end
    
    % --- Fit Seq (Individual Sets) ---
    seq_sigma_list = nan(1, nSets_seq); seq_peak_list = nan(1, nSets_seq);
    for ss = 1:nSets_seq
        ai_seq = find(Amps_seq == amp_val); if isempty(ai_seq), continue; end
        stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0); stimLoc = Coords(stimCh, :);
        [X, Y] = get_XY_data(Rseq, Coords, ss, ai_seq, idx_ptd, valid_channels, stimLoc);
        try
            beta_seq = nlinfit(X, Y, gauss_eq, beta0);
            p = min(max(beta_seq(1), 0), 1); s = min(max(abs(beta_seq(2)), 0), 2000);
            seq_sigma_list(ss) = s; seq_peak_list(ss) = p;
        catch; end
    end
    
    % Statistics
    valid_seq_s = seq_sigma_list(~isnan(seq_sigma_list));
    if ~isempty(valid_seq_s), seq_sigma_avg = mean(valid_seq_s); seq_sigma_sem = std(valid_seq_s)/sqrt(length(valid_seq_s)); else, seq_sigma_avg = NaN; seq_sigma_sem = NaN; end
    
    valid_seq_p = seq_peak_list(~isnan(seq_peak_list));
    if ~isempty(valid_seq_p), seq_peak_avg = mean(valid_seq_p); seq_peak_sem = std(valid_seq_p)/sqrt(length(valid_seq_p)); else, seq_peak_avg = NaN; seq_peak_sem = NaN; end
    
    % Store Rows
    Sigma_Results = [Sigma_Results; amp_val, Sigma_sim, seq_sigma_avg, seq_sigma_sem, seq_sigma_list];
    Peak_Results  = [Peak_Results;  amp_val, P_peak_sim, seq_peak_avg, seq_peak_sem, seq_peak_list];
end
warning('on', 'stats:nlinfit:ModelConstantWRTParam'); warning('on', 'stats:nlinfit:IllConditionedJacobian');

%% ================= STEP 4: PLOTTING =================
prof_colors = [0.25 0.45 0.65; 0.80 0.45 0.30; 0.55 0.40 0.65]; 

% --- PREPARE BAR DATA ---
all_dists = Sim_Results_Plot(:, 1); for ss=1:nSets_seq, all_dists = [all_dists; Seq_Results_Plot{ss}(:,1)]; end
Unique_Dists = unique(all_dists); Unique_Dists = Unique_Dists(Unique_Dists <= 600);
Plot_Data = zeros(length(Unique_Dists), 1 + nSets_seq);
for i = 1:length(Unique_Dists)
    d = Unique_Dists(i); mask = (Sim_Results_Plot(:,1) == d);
    if sum(mask)>0, Plot_Data(i, 1) = mean(Sim_Results_Plot(mask, 2)) * 100; end
end
for ss = 1:nSets_seq
    dat = Seq_Results_Plot{ss};
    for i = 1:length(Unique_Dists)
        d = Unique_Dists(i); mask = (dat(:,1) == d);
        if sum(mask)>0, Plot_Data(i, 1+ss) = mean(dat(mask, 2)) * 100; end
    end
end

% --- FIG 1: BAR CHART ---
figure('Color','w', 'Position',[200 600 800 500]); hold on;
leg_names = {'Simultaneous'}; for ss = 1:nSets_seq, leg_names{end+1}=sprintf('Seq Set %d', ss); end
b = bar(Unique_Dists, Plot_Data, 'grouped');
for i = 1:length(b), if i<=size(prof_colors,1), b(i).FaceColor=prof_colors(i,:); end; b(i).EdgeColor='none'; b(i).FaceAlpha=0.85; end
ylabel('% Responding', 'FontSize', 12, 'FontWeight', 'bold'); xlabel('Distance (µm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Spatial Spread at %.0f µA', target_amp), 'FontSize', 14); legend(leg_names, 'Location','northeast','Box','off'); box off;

% --- FIG 2: GAUSSIAN CURVES (Target Amp) ---
figure('Color','w', 'Position',[200 100 800 500]); hold on;
x_fit = linspace(0, 600, 200)';
try
    beta_sim = nlinfit(Sim_Results_Plot(:,1), Sim_Results_Plot(:,2), gauss_eq, [0.5, 100]);
    y_sim = gauss_eq(beta_sim, x_fit);
    plot(x_fit, y_sim, '-', 'Color', prof_colors(1,:), 'LineWidth', 3, 'DisplayName', leg_names{1});
catch; end
for ss = 1:nSets_seq
    dat = Seq_Results_Plot{ss}; if isempty(dat), continue; end
    col_idx = 1 + ss; if col_idx > size(prof_colors,1), col_idx = size(prof_colors,1); end
    try
        beta_seq = nlinfit(dat(:,1), dat(:,2), gauss_eq, [0.5, 100]);
        y_seq = gauss_eq(beta_seq, x_fit);
        plot(x_fit, y_seq, '-', 'Color', prof_colors(col_idx,:), 'LineWidth', 3, 'DisplayName', leg_names{col_idx});
    catch; end
end
ylabel('Probability', 'FontSize', 10, 'FontWeight', 'bold'); xlabel('Distance (µm)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Gaussian Spatial Fit at %.0f µA', target_amp), 'FontSize', 14); legend('Location','northeast','Box','off'); box off;

% --- FIG 3: PEAK PROBABILITY vs AMPLITUDE (New Figure) ---
figure('Color','w', 'Position',[800 100 600 500]); hold on;
amps_x = Peak_Results(:, 1);
% 1. Simultaneous
plot(amps_x, Peak_Results(:, 2), '-o', 'Color', prof_colors(1,:), 'LineWidth', 2, 'MarkerFaceColor', prof_colors(1,:), 'DisplayName', 'Simultaneous');
% 2. Sequential (Individual Sets)
for ss = 1:nSets_seq
    col_idx = 1 + ss; if col_idx > size(prof_colors,1), col_idx = size(prof_colors,1); end
    col = prof_colors(col_idx, :);
    set_data = Peak_Results(:, 4 + ss); % Col 5 is Set 1, Col 6 is Set 2, etc.
    plot(amps_x, set_data, '-s', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', sprintf('Seq Set %d', ss));
end
xlabel('Amplitude (µA)', 'FontWeight','bold'); ylabel('Peak Probability (P_{peak})', 'FontWeight','bold');
title('Response Intensity vs Amplitude', 'FontWeight','bold');
legend('Location','best','Box','off'); grid on; box off; ylim([0 1.1]);

%% ================= OUTPUT & SAVE =================
fprintf('\n\n====\n   FINAL RESULT\n========\n');

% --- TABLE 1: SIGMA (SPREAD) ---
fprintf('\n--- ACTIVATION SPREAD (Sigma) ---\n');
fprintf('Amp(uA)\tSim\tMean\tSEM');
for ss=1:nSets_seq, fprintf('\tSet%d', ss); end; fprintf('\n');
for i = 1:size(Sigma_Results, 1)
    fprintf('%.0f\t%.2f\t%.2f\t%.2f', Sigma_Results(i,1), Sigma_Results(i,2), Sigma_Results(i,3), Sigma_Results(i,4));
    for k = 5:size(Sigma_Results, 2), fprintf('\t%.2f', Sigma_Results(i,k)); end; fprintf('\n');
end

% --- TABLE 2: PEAK PROBABILITY (INTENSITY) ---
fprintf('\n--- PEAK PROBABILITY (Intensity) ---\n');
fprintf('Amp(uA)\tSim\tMean\tSEM');
for ss=1:nSets_seq, fprintf('\tSet%d', ss); end; fprintf('\n');
for i = 1:size(Peak_Results, 1)
    fprintf('%.0f\t%.2f\t%.2f\t%.2f', Peak_Results(i,1), Peak_Results(i,2), Peak_Results(i,3), Peak_Results(i,4));
    for k = 5:size(Peak_Results, 2), fprintf('\t%.2f', Peak_Results(i,k)); end; fprintf('\n');
end

save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Set6_Spatial_Sigma_Stats_5ms_' exp_id '.mat']); 

ResultSpatial = struct();
ResultSpatial.Sigma_Table = Sigma_Results; 
ResultSpatial.Peak_Table  = Peak_Results; 
ResultSpatial.BarData.Distances = Unique_Dists;
ResultSpatial.BarData.Probs = Plot_Data; 
ResultSpatial.BarData.TargetAmp = target_amp;
save(out_filename, 'ResultSpatial');
fprintf('\nResults saved to %s\n', out_filename);

%% ================= HELPER FUNCTIONS =================
function [X, Y] = get_XY_data(R, Coords, set_idx, amp_idx, ptd_idx, valid_ch, stim_locs)
    X = []; Y = []; chan_structs = R.set(set_idx).amp(amp_idx).ptd(ptd_idx).channel;
    for ch = valid_ch
        rec_loc = Coords(ch, :); d = round(min(sqrt(sum((stim_locs - rec_loc).^2, 2))));
        if d < 1, continue; end
        is_resp = 0; if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive, is_resp = 1; end
        X = [X; d]; Y = [Y; is_resp];
    end
end
function Coords = get_electrode_coordinates(type)
    Coords = zeros(64, 2); 
    if type == 1, for ch = 1:32, Coords(ch, 1) = 0; Coords(ch, 2) = (ch-1)*50; end
    elseif type == 2, for i = 1:16, ch=i; Coords(ch,1)=0; Coords(ch,2)=(i-1)*50; end, for i = 1:16, ch=32+i; Coords(ch,1)=200; Coords(ch,2)=(i-1)*50; end, for i = 1:16, ch=48+i; Coords(ch,1)=400; Coords(ch,2)=(i-1)*50; end, for i = 1:16, ch=16+i; Coords(ch,1)=600; Coords(ch,2)=(i-1)*50; end
    end
end
function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder); f = dir('*RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file'); end
    S_sp = load(f(1).name); if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);
    QC.BadCh = []; QC.BadTrials = []; f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end; f_bt = dir('*.BadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end