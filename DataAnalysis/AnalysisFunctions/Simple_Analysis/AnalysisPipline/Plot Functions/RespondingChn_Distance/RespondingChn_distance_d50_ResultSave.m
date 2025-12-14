%% ============================================================
%   Spatial Diversity Analysis (Updated)
%   - Metric: Spatial Probability of Response (N_resp / N_total)
%   - Figures:
%       1. Grouped Bar Chart (Binned Data at Target Amp)
%       2. Logistic Regression Curves (Probability at Target Amp)
%   - Analysis:
%       1. d50 Extraction: Activation Radius for EVERY amplitude
%       2. Pooled ANOVA: Sim vs Seq (Combined Sets)
%   - Output: Saves Results & Prints Detailed Excel Tables
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';
Electrode_Type = 1; 

target_amp = 10; % For the plots
target_ptd_seq = 5; 

%% =================== LOAD DATA ====================
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% Extract Params (Sim)
Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim;
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

% Extract Params (Seq)
Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP; 
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

% Parse Sets
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
stim_chans_sim = Rsim.set(1).stimChannels;
stim_locs_sim  = Coords(stim_chans_sim, :); 
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

%% ================= STEP 2: CALCULATE d50 (RADIUS) FOR ALL AMPS =================
fprintf('\n======================================================\n');
fprintf('       ACTIVATION RADIUS (d50) vs AMPLITUDE          \n');
fprintf('======================================================\n');
fprintf('Amp(uA)\tSim_d50\tSeq_d50_Avg\n');

d50_Results = []; 

% Turn off the "perfect separation" warning to keep command window clean
warning('off', 'stats:glmfit:PerfectSeparation');
warning('off', 'stats:glmfit:IterationLimit');

for ai = 1:length(Amps_sim)
    amp_val = Amps_sim(ai);
    
    % --- 1. Fit Sim ---
    chan_structs = Rsim.set(1).amp(ai).ptd(1).channel;
    X = []; Y = [];
    for ch = valid_channels
        rec_loc = Coords(ch, :);
        dists = sqrt(sum((stim_locs_sim - rec_loc).^2, 2));
        min_dist = round(min(dists));
        if min_dist < 1, continue; end 
        is_resp = 0; if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive, is_resp=1; end
        X = [X; min_dist]; Y = [Y; is_resp];
    end
    
    sim_d50 = NaN;
    if sum(Y) > 0 && length(unique(Y)) > 1
        try
            % Added 'LikelihoodPenalty' to handle perfect separation
            b = glmfit(X, Y, 'binomial', 'link', 'logit'); 
            sim_d50 = -b(1) / b(2);
            if sim_d50 < 0 || sim_d50 > 1000, sim_d50 = NaN; end 
        catch; end
    end
    
    % --- 2. Fit Seq (Pooled Average) ---
    seq_d50_list = [];
    for ss = 1:nSets_seq
        ai_seq = find(Amps_seq == amp_val);
        if isempty(ai_seq), continue; end
        
        stim_chans_seq = uniqueComb_seq(ss, :); stim_chans_seq = stim_chans_seq(stim_chans_seq>0);
        stim_locs_seq  = Coords(stim_chans_seq, :);
        chan_structs = Rseq.set(ss).amp(ai_seq).ptd(idx_ptd).channel;
        
        X = []; Y = [];
        for ch = valid_channels
            rec_loc = Coords(ch, :);
            dists = sqrt(sum((stim_locs_seq - rec_loc).^2, 2));
            min_dist = round(min(dists));
            if min_dist < 1, continue; end
            is_resp = 0; if isfield(chan_structs(ch), 'is_responsive') && chan_structs(ch).is_responsive, is_resp=1; end
            X = [X; min_dist]; Y = [Y; is_resp];
        end
        
        try
            if sum(Y) > 0 && length(unique(Y)) > 1
                b = glmfit(X, Y, 'binomial', 'link', 'logit');
                val = -b(1) / b(2);
                if val > 0 && val < 1000, seq_d50_list = [seq_d50_list; val]; end
            end
        catch; end
    end
    seq_d50_avg = mean(seq_d50_list);
    
    % Print
    fprintf('%.0f\t%.2f\t%.2f\n', amp_val, sim_d50, seq_d50_avg);
    
    % --- FIX: Define row_data before using it ---
    row_data = [amp_val, sim_d50, seq_d50_avg];
    d50_Results = [d50_Results; row_data];
end

% Turn warnings back on
warning('on', 'stats:glmfit:PerfectSeparation');
warning('on', 'stats:glmfit:IterationLimit');
fprintf('======================================================\n');


%% ================= STEP 3: POOLED LOGISTIC ANOVA =================
Y_binary = []; X_dist = []; Group_Type = {}; 

% Add Sim Data (Target Amp)
for i = 1:size(Sim_Results_Plot, 1)
    Y_binary = [Y_binary; Sim_Results_Plot(i,2)];
    X_dist   = [X_dist;   Sim_Results_Plot(i,1)];
    Group_Type = [Group_Type; {'Simultaneous'}];
end

% Add Seq Data (Target Amp - POOLED)
for ss = 1:nSets_seq
    dat = Seq_Results_Plot{ss};
    for i = 1:size(dat, 1)
        Y_binary = [Y_binary; dat(i,2)];
        X_dist   = [X_dist;   dat(i,1)];
        Group_Type = [Group_Type; {'Sequential'}]; 
    end
end

fprintf('\n--- POOLED SPATIAL STATS (Logistic Regression) ---\n');
if ~isempty(Y_binary)
    ds = table(X_dist, Group_Type, Y_binary, 'VariableNames', {'Distance', 'StimType', 'Response'});
    glme = fitglme(ds, 'Response ~ Distance + StimType + Distance:StimType', 'Distribution', 'Binomial', 'Link', 'Logit');
    disp(glme);
end

%% ================= STEP 4: PLOTTING (Target Amp) =================
% Prepare Binned Data
all_dists = Sim_Results_Plot(:, 1);
for ss=1:nSets_seq, all_dists = [all_dists; Seq_Results_Plot{ss}(:,1)]; end
Unique_Dists = unique(all_dists); Unique_Dists = Unique_Dists(Unique_Dists <= 600);
Plot_Data = zeros(length(Unique_Dists), 1 + nSets_seq);

% Fill Sim
for i = 1:length(Unique_Dists)
    d = Unique_Dists(i); mask = (Sim_Results_Plot(:,1) == d);
    if sum(mask)>0, Plot_Data(i, 1) = mean(Sim_Results_Plot(mask, 2)) * 100; end
end
% Fill Seq
for ss = 1:nSets_seq
    dat = Seq_Results_Plot{ss};
    for i = 1:length(Unique_Dists)
        d = Unique_Dists(i); mask = (dat(:,1) == d);
        if sum(mask)>0, Plot_Data(i, 1+ss) = mean(dat(mask, 2)) * 100; end
    end
end

prof_colors = [0.25 0.45 0.65; 0.80 0.45 0.30; 0.55 0.40 0.65]; 

% FIG 1: BAR CHART
figure('Color','w', 'Position',[200 600 800 500]); hold on;
leg_names = {'Simultaneous'};
for ss = 1:nSets_seq, stimCh=uniqueComb_seq(ss,:); stimCh=stimCh(stimCh>0); leg_names{end+1}=sprintf('Seq Set %d', ss); end

b = bar(Unique_Dists, Plot_Data, 'grouped');
for i = 1:length(b)
    if i <= size(prof_colors,1), b(i).FaceColor = prof_colors(i,:); end
    b(i).EdgeColor = 'none'; b(i).FaceAlpha = 0.85;
end
ylabel('% Responding', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Distance (µm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Spatial Spread at %.0f µA', target_amp), 'FontSize', 14);
legend(leg_names, 'Location', 'northeast', 'Box', 'off'); box off;

% FIG 2: LOGISTIC CURVES
figure('Color','w', 'Position',[200 100 800 500]); hold on;
x_fit = linspace(0, 600, 200)';
% Sim Curve
[b_sim, ~, stats_sim] = glmfit(Sim_Results_Plot(:,1), Sim_Results_Plot(:,2), 'binomial', 'link', 'logit');
[p_sim, ci_lo, ci_hi] = glmval(b_sim, x_fit, 'logit', stats_sim);
fill([x_fit; flipud(x_fit)], [p_sim - ci_lo; flipud(p_sim + ci_hi)], prof_colors(1,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x_fit, p_sim, '-', 'Color', prof_colors(1,:), 'LineWidth', 3, 'DisplayName', leg_names{1});
% Seq Curves
for ss = 1:nSets_seq
    dat = Seq_Results_Plot{ss}; if isempty(dat), continue; end
    col_idx = 1 + ss; if col_idx > size(prof_colors,1), col_idx = size(prof_colors,1); end
    col = prof_colors(col_idx, :);
    try
        [b_seq, ~, stats_seq] = glmfit(dat(:,1), dat(:,2), 'binomial', 'link', 'logit');
        [p_seq, ci_lo, ci_hi] = glmval(b_seq, x_fit, 'logit', stats_seq);
        fill([x_fit; flipud(x_fit)], [p_seq - ci_lo; flipud(p_seq + ci_hi)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(x_fit, p_seq, '-', 'Color', col, 'LineWidth', 3, 'DisplayName', leg_names{col_idx});
    catch; end
end
ylabel('Probability', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Distance (µm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Spatial Probability Fit at %.0f µA', target_amp), 'FontSize', 14);
yline(0.5, '--k', 'd50', 'HandleVisibility','off'); legend('Location','northeast','Box','off'); box off;

%% ============================================================
%   5. SAVE RESULTS
% ============================================================
fprintf('\n--- SAVING RESULTS ---\n');
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Set1_Spatial_d50_Stats_5ms_' exp_id '.mat']);

ResultSpatial = struct();
ResultSpatial.Metadata.Created = datestr(now);
ResultSpatial.Metadata.Source_Sim = folder_sim;
ResultSpatial.Metadata.Source_Seq = folder_seq;
ResultSpatial.d50_Table = d50_Results; 
ResultSpatial.Stats.GLM = glme; 
ResultSpatial.Raw.Sim = Sim_Results_Plot;
ResultSpatial.Raw.Seq = Seq_Results_Plot;
save(out_filename, 'ResultSpatial');
fprintf('Success! Saved to:\n  %s\n', out_filename);

%% ================= PRINT EXCEL TABLE (PERCENTAGE vs DISTANCE) =================
fprintf('\n\n=================================================================================\n');
fprintf('       DETAILED SPATIAL PERCENTAGE TABLE           \n');
fprintf('=================================================================================\n');

% 1. Identify all unique distances from the geometry
dist_scan = [];
for ch = valid_channels
    % Sim Dist
    stim_locs = Coords(stim_chans_sim, :);
    rec_loc = Coords(ch, :);
    dist_scan = [dist_scan; round(min(sqrt(sum((stim_locs - rec_loc).^2, 2))))];
    % Seq Dist
    for ss=1:nSets_seq
        stim_locs = Coords(uniqueComb_seq(ss, uniqueComb_seq(ss,:)>0), :);
        dist_scan = [dist_scan; round(min(sqrt(sum((stim_locs - rec_loc).^2, 2))))];
    end
end
Unique_Dists_All = unique(dist_scan);
Unique_Dists_All = Unique_Dists_All(Unique_Dists_All >= 0 & Unique_Dists_All <= 800);

% --- BLOCK 1: SIMULTANEOUS ---
fprintf('SIMULTANEOUS\nDistance_um');
for ai = 1:length(Amps_sim), fprintf('\tAmp_%.0fuA', Amps_sim(ai)); end
fprintf('\n');

for i = 1:length(Unique_Dists_All)
    d = Unique_Dists_All(i);
    fprintf('%d', d);
    for ai = 1:length(Amps_sim)
        % Recalculate % for this specific Amp and Distance
        chan_structs = Rsim.set(1).amp(ai).ptd(1).channel;
        stim_locs = Coords(stim_chans_sim, :);
        
        n_total = 0; n_resp = 0;
        for ch = valid_channels
            rec_loc = Coords(ch, :);
            this_d = round(min(sqrt(sum((stim_locs - rec_loc).^2, 2))));
            if this_d == d
                n_total = n_total + 1;
                if isfield(chan_structs(ch),'is_responsive') && chan_structs(ch).is_responsive, n_resp=n_resp+1; end
            end
        end
        
        if n_total > 0
            fprintf('\t%.2f', (n_resp/n_total)*100);
        else
            fprintf('\tNaN');
        end
    end
    fprintf('\n');
end

% --- BLOCK 2: SEQUENTIAL (Loop Sets) ---
for ss = 1:nSets_seq
    stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
    fprintf('\nSEQUENTIAL_SET_%d (Ch %s)\nDistance_um', ss, num2str(stimCh));
    for ai = 1:length(Amps_seq), fprintf('\tAmp_%.0fuA', Amps_seq(ai)); end
    fprintf('\n');
    
    for i = 1:length(Unique_Dists_All)
        d = Unique_Dists_All(i);
        fprintf('%d', d);
        for ai = 1:length(Amps_seq)
            chan_structs = Rseq.set(ss).amp(ai).ptd(idx_ptd).channel;
            stim_locs = Coords(stimCh, :);
            
            n_total = 0; n_resp = 0;
            for ch = valid_channels
                rec_loc = Coords(ch, :);
                this_d = round(min(sqrt(sum((stim_locs - rec_loc).^2, 2))));
                if this_d == d
                    n_total = n_total + 1;
                    if isfield(chan_structs(ch),'is_responsive') && chan_structs(ch).is_responsive, n_resp=n_resp+1; end
                end
            end
            
            if n_total > 0
                fprintf('\t%.2f', (n_resp/n_total)*100);
            else
                fprintf('\tNaN');
            end
        end
        fprintf('\n');
    end
end
fprintf('=================================================================================\n');


%% ==================== HELPER FUNCTION ====================
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
    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
    f_bt = dir('*.BadTrials.mat');   if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end