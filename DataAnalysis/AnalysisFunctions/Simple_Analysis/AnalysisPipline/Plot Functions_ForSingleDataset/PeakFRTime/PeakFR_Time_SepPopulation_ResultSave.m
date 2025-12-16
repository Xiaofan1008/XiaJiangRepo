%% ============================================================
%   Peak Latency Analysis (Refined: Only Responsive Channels)
%   - Metric: Time of Peak Firing Rate (ms)
%   - Logic: Latency is ONLY calculated if the channel is 
%            responsive to that SPECIFIC condition.
%   - Figures:
%       1. Grouped Histograms (Distribution at Target Amp)
%       2. Average Latency vs Amplitude (Curves)
%   - Output: Saves Result & Prints Summary for Excel
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim6_251125_181554';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq6_5ms_251125_182437';
Electrode_Type = 1;

% Analysis
search_win_ms = [2 20]; % Start >0 to avoid stim artifact 
FS = 30000;            

% Smoothing (SYMMETRIC Gaussian)
bin_ms = 1; sigma_ms = 3; sigma_bins = sigma_ms / bin_ms;

% PLOTTING SETTINGS
target_amp = 10; % Amplitude for Histogram
target_ptds_for_hist = [5]; 
hist_edges = 0:1:20; % Bins

% Kernel
edges = -100:bin_ms:100; bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g = gausswin(kernel_size); g = g / sum(g); 

%% =================== 1. LOAD DATA & QC INFO ====================
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
stim_chans_sim = Rsim.set(1).stimChannels; 
d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

%% ================= 2. IDENTIFY ALL POTENTIAL CHANNELS =============
% We keep the union list just for indexing, but we will filter data later.
resp_sim = false(nCh_Total, 1); for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive')&&this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1); for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive')&&this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end
resp_channels = find(resp_sim | resp_seq); 

%% =================== 3. COMPUTE LATENCY (SEPARATE POPULATIONS) =================
LatAmp_sim = nan(length(resp_channels), length(Amps_sim));
LatAmp_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq); 

% --- SIMULTANEOUS ---
for ai = 1:length(Amps_sim)
    tr_ids = find(ampIdx_sim == ai); 
    if isempty(tr_ids), continue; end
    
    for ci = 1:length(resp_channels)
        ch_idx = resp_channels(ci); recCh = d(ch_idx);
        
        % 1. CHECK RESPONSIVENESS (Condition Specific)
        try 
            is_resp = Rsim.set(1).amp(ai).ptd(1).channel(ch_idx).is_responsive;
        catch
            is_resp = false;
        end
        
        if ~is_resp
            LatAmp_sim(ci, ai) = NaN; % Explicitly exclude non-responders
            continue; 
        end
        
        % 2. QC Filter
        bad_trs = []; if ~isempty(QC_Sim.BadTrials) && ch_idx<=numel(QC_Sim.BadTrials), bad_trs = QC_Sim.BadTrials{ch_idx}; end
        tr_ids_clean = setdiff(tr_ids, bad_trs);
        
        S_ch = sp_sim{recCh};
        LatAmp_sim(ci,ai) = get_peak_latency(tr_ids_clean, trig_sim, S_ch, search_win_ms, bin_ms, g, FS);
    end
end

% --- SEQUENTIAL ---
p_idx_trend = find(abs(PTDs_ms - target_ptds_for_hist(1)) < 0.1, 1);
for ss = 1:nSets_seq
    for ai = 1:length(Amps_seq)
        tr_ids = find(combClass_seq==ss & ptdIdx_seq==p_idx_trend & ampIdx_seq==ai);
        if isempty(tr_ids), continue; end
        
        for ci = 1:length(resp_channels)
            ch_idx = resp_channels(ci); recCh = d(ch_idx);
            
            % 1. CHECK RESPONSIVENESS (Condition Specific)
            try
                is_resp = Rseq.set(ss).amp(ai).ptd(p_idx_trend).channel(ch_idx).is_responsive;
            catch
                is_resp = false;
            end
            
            if ~is_resp
                LatAmp_seq(ci, ai, ss) = NaN; % Explicitly exclude non-responders
                continue;
            end
            
            % 2. QC Filter
            bad_trs = []; if ~isempty(QC_Seq.BadTrials) && ch_idx<=numel(QC_Seq.BadTrials), bad_trs = QC_Seq.BadTrials{ch_idx}; end
            tr_ids_clean = setdiff(tr_ids, bad_trs);
            
            S_ch = sp_seq{recCh};
            LatAmp_seq(ci,ai,ss) = get_peak_latency(tr_ids_clean, trig_seq, S_ch, search_win_ms, bin_ms, g, FS);
        end
    end
end

%% ===================== FIGURE 1: GROUPED HISTOGRAMS ======================
idx_sim_t = find(Amps_sim == target_amp, 1);
idx_seq_t = find(Amps_seq == target_amp, 1);

% Prepare Data (NaNs are already filtered out above)
Lat_sim_dist = LatAmp_sim(:, idx_sim_t); Lat_sim_dist = Lat_sim_dist(~isnan(Lat_sim_dist));
Lat_seq_dist = cell(nSets_seq, 1);
Lat_seq_all = []; % Storage for combined Seq data

for ss=1:nSets_seq
    temp = squeeze(LatAmp_seq(:, idx_seq_t, ss));
    temp = temp(~isnan(temp));
    Lat_seq_dist{ss} = temp;
    Lat_seq_all = [Lat_seq_all; temp];
end

figure('Color','w', 'Position',[100 100 900 600]); hold on;
bar_data = []; group_names = {}; group_colors = {};

% Sim Data
if isempty(Lat_sim_dist)
    bar_data(:, 1) = zeros(length(hist_edges)-1, 1);
else
    [counts_sim, ~] = histcounts(Lat_sim_dist, hist_edges);
    bar_data(:, 1) = counts_sim(:) / length(Lat_sim_dist) * 100;
end
group_names{1} = 'Simultaneous'; group_colors{1} = [0 0.3 0.8];

% Seq Data
cols_seq = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 
col_idx = 2;
for ss = 1:nSets_seq
    data = Lat_seq_dist{ss}; 
    if isempty(data)
        bar_data(:, col_idx) = zeros(length(hist_edges)-1, 1);
    else
        [counts_seq, ~] = histcounts(data, hist_edges);
        bar_data(:, col_idx) = counts_seq(:) / length(data) * 100;
    end
    stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
    group_names{col_idx} = sprintf('Seq Set %d', ss);
    group_colors{col_idx} = cols_seq(mod(ss-1,3)+1, :);
    col_idx = col_idx + 1;
end

% Plot Bars
b = bar(hist_edges(1:end-1)+diff(hist_edges)/2, bar_data, 'grouped');
for i = 1:length(b), b(i).FaceColor = group_colors{i}; b(i).EdgeColor = 'none'; b(i).FaceAlpha = 0.3; end

% Define scale globally
scale = 100 * (hist_edges(2)-hist_edges(1)); 

% Plot Smooth Curves & Medians (Sim)
if length(Lat_sim_dist) > 1 
    [f_sim, xi_sim] = ksdensity(Lat_sim_dist, 'Bandwidth', 1.5);
    plot(xi_sim, f_sim*scale, 'Color', group_colors{1}, 'LineWidth', 2.5);
    xline(median(Lat_sim_dist), '--', 'Color', group_colors{1}, 'LineWidth', 2);
end

% Plot Smooth Curves & Medians (Seq)
for i = 2:length(group_names)
    ss_idx = i - 1; data = Lat_seq_dist{ss_idx};
    if length(data) > 1
        [f_seq, xi_seq] = ksdensity(data, 'Bandwidth', 1.5);
        plot(xi_seq, f_seq*scale, 'Color', group_colors{i}, 'LineWidth', 2.5);
        xline(median(data), '--', 'Color', group_colors{i}, 'LineWidth', 2);
    end
end
ylabel('% of Channels', 'FontSize', 10, 'FontWeight', 'bold'); 
xlabel('Peak Latency (ms)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Latency Distribution at %.0f µA', target_amp), 'FontSize', 14);
legend(b, group_names, 'Location','northeast'); box off; xlim([0 20]);

%% ===================== FIGURE 2: AVERAGE LATENCY VS AMPLITUDE ======================
figure('Color','w', 'Position',[900 100 600 500]); hold on;

% 1. Simultaneous Curve
mu_sim = mean(LatAmp_sim, 1, 'omitnan'); 
err_sim = std(LatAmp_sim,0,1,'omitnan') ./ sqrt(sum(~isnan(LatAmp_sim),1));

errorbar(Amps_sim, mu_sim, err_sim, '-o', 'Color', group_colors{1}, ...
    'LineWidth', 2, 'MarkerFaceColor', group_colors{1}, 'CapSize', 0, ...
    'DisplayName', 'Simultaneous');

% 2. Sequential Curves (Per Set)
for ss = 1:nSets_seq
    data = squeeze(LatAmp_seq(:, :, ss));
    mu_seq = mean(data, 1, 'omitnan'); 
    err_seq = std(data,0,1,'omitnan') ./ sqrt(sum(~isnan(data),1));
    
    col = cols_seq(mod(ss-1,3)+1, :);
    stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
    lbl = sprintf('Seq Set %d (Ch:%s)', ss, num2str(stimCh));
    
    errorbar(Amps_seq, mu_seq, err_seq, '-s', 'Color', col, ...
        'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 0, ...
        'DisplayName', lbl);
end

xlabel('Amplitude (µA)', 'FontWeight','bold'); 
ylabel('Mean Peak Latency (ms)', 'FontWeight','bold');
title('Average Peak Firing Rate Time vs. Amplitude', 'FontWeight','bold'); 
legend('Location','best'); box off;

%% ===================== STATS: TWO-WAY ANOVA (POOLED) ======================
y_values = []; g_stim = []; g_amp = [];

% 1. Sim
for ai = 1:length(Amps_sim)
    lats = LatAmp_sim(:, ai); lats = lats(~isnan(lats));
    if isempty(lats), continue; end
    y_values = [y_values; lats];
    g_stim   = [g_stim; repmat({'Simultaneous'}, length(lats), 1)];
    g_amp    = [g_amp;  repmat(Amps_sim(ai), length(lats), 1)];
end

% 2. Seq (MERGED)
for ss = 1:nSets_seq
    % --- POOLED GROUP NAME ---
    group_name = 'Sequential'; 
    % -------------------------
    
    for ai = 1:length(Amps_seq)
        lats = squeeze(LatAmp_seq(:, ai, ss)); lats = lats(~isnan(lats));
        if isempty(lats), continue; end
        y_values = [y_values; lats];
        g_stim   = [g_stim; repmat({group_name}, length(lats), 1)];
        g_amp    = [g_amp;  repmat(Amps_seq(ai), length(lats), 1)];
    end
end

fprintf('\n--- ANOVA RESULTS (POOLED) ---\n');
if length(unique(g_stim)) > 1
    [p_vals, ~, stats] = anovan(y_values, {g_stim, g_amp}, 'model', 'interaction', 'varnames', {'StimType', 'Amplitude'}, 'display', 'on');
    fprintf('StimType P: %.5f\nAmp P: %.5f\nInteraction P: %.5f\n', p_vals(1), p_vals(2), p_vals(3));
    figure('Color','w', 'Name', 'Post-Hoc Comparison'); multcompare(stats, 'Dimension', 1);
else
    fprintf('Not enough groups for ANOVA.\n');
end

%% ============================================================
%   5. SAVE RESULTS
% ============================================================
fprintf('\n--- SAVING RESULTS ---\n');
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Set6_Latency_SepPop_5ms_' exp_id '.mat']);

ResultLat = struct();
ResultLat.Metadata.Created = datestr(now);
ResultLat.Metadata.Source_Sim = folder_sim;
ResultLat.Metadata.Source_Seq = folder_seq;
ResultLat.Metadata.StimGroups.Sim = stim_chans_sim;
ResultLat.Metadata.StimGroups.Seq = uniqueComb_seq(:, uniqueComb_seq(1,:) > 0);
ResultLat.Metadata.Amps = Amps_sim;
ResultLat.Metadata.RespChannels = resp_channels;
% Store Latency Data
ResultLat.Latency.Sim = LatAmp_sim; 
ResultLat.Latency.Seq = LatAmp_seq; 
% Store QC Data
ResultLat.QC.Sim = QC_Sim;
ResultLat.QC.Seq = QC_Seq;
save(out_filename, 'ResultLat');
fprintf('Latency Results saved to:\n  %s\n', out_filename);

%% ================= PRINT SUMMARY STATISTICS (EXCEL READY) =================
fprintf('\n\n========================================================================\n');
fprintf('       SUMMARY STATISTICS: PEAK LATENCY          \n');
fprintf('========================================================================\n');

% 1. Print Header Row
fprintf('Condition      \t');
for ai = 1:length(Amps_sim)
    fprintf('Amp_%.0fuA_Mean   SEM\t', Amps_sim(ai), Amps_sim(ai));
end
fprintf('\n');

% 2. Simultaneous Row
fprintf('Simultaneous\t');
for ai = 1:length(Amps_sim)
    d = LatAmp_sim(:, ai);
    mu = mean(d, 'omitnan');
    sem = std(d, 0, 'omitnan') / sqrt(sum(~isnan(d)));
    if isnan(mu), fprintf('NaN\tNaN\t'); else, fprintf('%.4f\t%.4f\t', mu, sem); end
end
fprintf('\n');

% 3. Sequential Rows (Per Set)
for ss = 1:nSets_seq
    stimCh_List = uniqueComb_seq(ss, :); stimCh_List = stimCh_List(stimCh_List > 0);
    fprintf('Seq_Set%d(Ch%s)\t', ss, num2str(stimCh_List));
    
    for ai = 1:length(Amps_seq)
        d = squeeze(LatAmp_seq(:, ai, ss)); 
        mu = mean(d, 'omitnan');
        sem = std(d, 0, 'omitnan') / sqrt(sum(~isnan(d)));
        if isnan(mu), fprintf('NaN\tNaN\t'); else, fprintf('%.4f\t%.4f\t', mu, sem); end
    end
    fprintf('\n');
end
fprintf('========================================================================\n');

%% ==================== HELPER FUNCTIONS ====================
function peak_time_ms = get_peak_latency(tr_ids, trig, sp_data, search_win, bin_ms, kern, FS)
    nTr = numel(tr_ids); if nTr == 0, peak_time_ms = NaN; return; end
    all_spikes = [];
    for k = 1:nTr, tr=tr_ids(k); t0=trig(tr)/FS*1000; tt=sp_data(:,1)-t0; all_spikes=[all_spikes; tt(tt>=-20 & tt<=80)]; end
    
    bin_s = bin_ms/1000; edges = -20:bin_ms:80; t_centers = edges(1:end-1)+bin_ms/2;
    h = histcounts(all_spikes, edges); rate_smooth = conv(h/(nTr*bin_s), kern, 'same');
    
    valid_idx = find(t_centers >= search_win(1) & t_centers <= search_win(2));
    valid_rate = rate_smooth(valid_idx); valid_time = t_centers(valid_idx);
    
    if isempty(valid_rate) || max(valid_rate)==0, peak_time_ms = NaN; else
        [~, max_idx] = max(valid_rate); peak_time_ms = valid_time(max_idx);
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