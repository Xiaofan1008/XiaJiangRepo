%% ============================================================
%   Response Duration Analysis (FINAL)
%   - Metric: Response Span (End Time - Start Time)
%   - Method: Dual Threshold (Detect > 3.0 SD, Measure > 1.5 SD)
%   - Filter: Zero-Phase Symmetric Gaussian
%   - Plot: Line+Scatter (Trend) AND Box Plot (Stats at 10uA)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type = 1;

% ANALYSIS SETTINGS
% We look for the response start/end within this window
search_win_ms = [2 20];  
FS = 30000;            

% Smoothing (SYMMETRIC / ZERO-PHASE)
bin_ms     = 1;
sigma_ms   = 3;        
sigma_bins = sigma_ms / bin_ms;

% ---- DUAL THRESHOLDS ----
baseline_win_ms = [-90 -10];
baseDur_s       = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;

z_detect   = 3.0; % Peak must reach this (Reliability)
z_boundary = 1.5; % Width measured at this level (Full Span)

% Plotting
jitter_width  = 0.2;     
scatter_alpha = 0.4;    
dot_size      = 25;          

% Create Symmetric Kernel
edges = -100:bin_ms:100; 
bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g = gausswin(kernel_size); g = g / sum(g); 

%% =================== LOAD DATA ====================
% (Standard Loading Block)
cd(folder_sim);
Rsim_file = dir('*RespondingChannels.mat'); Rsim = load(Rsim_file(1).name).Responding;
sp_sim = load_ssd_spikes(folder_sim);   
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig_sim = loadTrig(0);
Ssim = load(dir('*_exp_datafile_*.mat').name); Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim;
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

cd(folder_seq);
Rseq_file = dir('*RespondingChannels.mat'); Rseq = load(Rseq_file(1).name).Responding;
sp_seq = load_ssd_spikes(folder_seq);
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig_seq = loadTrig(0);
Sseq = load(dir('*_exp_datafile_*.mat').name); Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP; 
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);
stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq); for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v  = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

%% ================= DETERMINE RESPONDING CHANNELS =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_sim = false(nCh_Total, 1);
for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1);
for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end
resp_channels = find(resp_sim | resp_seq); 

%% =================== COMPUTE DURATION =================
Dur_sim = nan(length(resp_channels), length(Amps_sim));
Dur_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch = resp_channels(ci); recCh = d(ch);
    
    % --- 1. SIMULTANEOUS ---
    S_ch = sp_sim{recCh};
    for ai = 1:length(Amps_sim)
        tr_ids = find(ampIdx_sim == ai);
        if isempty(tr_ids), continue; end
        
        Dur_sim(ci,ai) = get_duration_dual(tr_ids, trig_sim, S_ch, ...
            baseline_win_ms, baseDur_s, search_win_ms, bin_ms, g, z_detect, z_boundary, FS);
    end
    
    % --- 2. SEQUENTIAL ---
    S_ch = sp_seq{recCh};
    for p = 1:numel(PTDs_ms)
        for ss = 1:nSets_seq
            for ai = 1:length(Amps_seq)
                tr_ids = find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai);
                if isempty(tr_ids), continue; end
                
                Dur_seq(ci,ai,ss,p) = get_duration_dual(tr_ids, trig_seq, S_ch, ...
                    baseline_win_ms, baseDur_s, search_win_ms, bin_ms, g, z_detect, z_boundary, FS);
            end
        end
    end
end 

%% ===================== PLOT: 2-PANEL FIGURE ======================
figure('Color','w', 'Position',[100 100 1100 550]); 
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- PANEL 1: LINE PLOT (Trend vs Amplitude) ---
ax1 = nexttile([1 2]); hold(ax1, 'on');

% 1. Simultaneous
AvgSim = mean(Dur_sim, 1, 'omitnan');
SEMSim = std(Dur_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Dur_sim), 1));
sim_col = [0 0.3 0.8]; 

if any(~isnan(AvgSim))
    for i = 1:length(Amps_sim)
        x_jit = (rand(size(Dur_sim(:,i))) - 0.5) * jitter_width + Amps_sim(i);
        scatter(ax1, x_jit, Dur_sim(:,i), dot_size, sim_col, 'filled', ...
            'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
    end
    plot_shaded_error(Amps_sim, AvgSim, SEMSim, sim_col);
    plot(ax1, Amps_sim, AvgSim, '-o', 'Color', sim_col, 'LineWidth', 2, ...
        'MarkerFaceColor', sim_col, 'DisplayName', 'Simultaneous');
end

% 2. Sequential
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 

for p = 1:numel(PTDs_ms)
    for ss = 1:nSets_seq
        data_set = squeeze(Dur_seq(:, :, ss, p));
        AvgSeq = mean(data_set, 1, 'omitnan');
        SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
        
        col = set_colors(mod(ss-1,3)+1, :);
        stimCh_List = uniqueComb_seq(ss, :); stimCh_List = stimCh_List(stimCh_List > 0);
        lbl = sprintf('Seq Set %d (Ch: %s)', ss, num2str(stimCh_List));
        
        if any(~isnan(AvgSeq))
            for i = 1:length(Amps_seq)
                x_jit = (rand(size(data_set(:,i))) - 0.5) * jitter_width + Amps_seq(i);
                scatter(ax1, x_jit, data_set(:,i), dot_size, col, 'filled', ...
                    'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
            end
            plot_shaded_error(Amps_seq, AvgSeq, SEMSeq, col);
            plot(ax1, Amps_seq, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, ...
                'MarkerFaceColor', 'w', 'DisplayName', lbl);
        end
    end
end

xlabel(ax1, 'Amplitude (µA)', 'FontWeight', 'bold');
ylabel(ax1, 'Response Duration (ms)', 'FontWeight', 'bold');
title(ax1, 'Duration vs. Amplitude', 'FontWeight', 'bold');
legend(ax1, 'Location','northwest');

% --- PANEL 2: BOX PLOT (At 10 uA) ---
ax2 = nexttile; hold(ax2, 'on');
target_amp = 10;
idx_sim = find(Amps_sim == target_amp, 1);
idx_seq = find(Amps_seq == target_amp, 1);

% Prepare Data for Boxplot
bp_data = [];
bp_groups = {};
colors_bp = {};

% Sim Data
if ~isempty(idx_sim)
    d_sim = Dur_sim(:, idx_sim);
    bp_data = [bp_data; d_sim];
    bp_groups = [bp_groups; repmat({'Sim'}, length(d_sim), 1)];
end

% Seq Data
for ss = 1:nSets_seq
    d_seq = squeeze(Dur_seq(:, idx_seq, ss, 1));
    bp_data = [bp_data; d_seq];
    bp_groups = [bp_groups; repmat({sprintf('Seq Set %d', ss)}, length(d_seq), 1)];
end

% Plot Boxplot
boxplot(ax2, bp_data, bp_groups, 'Width', 0.5, 'Symbol', 'k.');
ylabel(ax2, 'Duration (ms)');
title(ax2, sprintf('Distribution at %.0f µA', target_amp), 'FontWeight', 'bold');

%% ==================== HELPER FUNCTION (DUAL THRESHOLD) ====================
function duration_ms = get_duration_dual(tr_ids, trig, sp_data, ...
    base_win, base_dur, search_win, bin_ms, kern, z_detect, z_boundary, FS)
    
    nTr = numel(tr_ids);
    all_spikes_base = [];
    all_spikes_post = [];
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        all_spikes_base = [all_spikes_base; tt(tt >= base_win(1) & tt < base_win(2))];
        all_spikes_post = [all_spikes_post; tt(tt >= -20 & tt <= 80)];
    end
    
    if nTr == 0, duration_ms = NaN; return; end
    
    % Compute Smoothed PSTH
    bin_s = bin_ms/1000;
    base_rate_avg = numel(all_spikes_base) / (nTr * base_dur);
    edges_psth = -20 : bin_ms : 80;
    t_centers  = edges_psth(1:end-1) + bin_ms/2;
    h = histcounts(all_spikes_post, edges_psth);
    rate_raw = h / (nTr * bin_s);
    rate_smooth = conv(rate_raw, kern, 'same');
    rate_smooth(t_centers < 0) = base_rate_avg; % Clamp Pre-Stim
    
    % Statistics
    baseline_part = rate_smooth(t_centers >= base_win(1) & t_centers < base_win(2));
    if isempty(baseline_part), sigma_base = 5; else, sigma_base = std(baseline_part); end
    sigma_base = max(sigma_base, 2.0); % Noise floor
    
    % --- DUAL THRESHOLD LOGIC ---
    thresh_high = base_rate_avg + z_detect * sigma_base;
    thresh_low  = base_rate_avg + z_boundary * sigma_base;
    
    % 1. VALIDITY CHECK
    valid_idx = find(t_centers >= search_win(1) & t_centers <= search_win(2));
    if isempty(valid_idx), duration_ms = 0; return; end
    
    peak_in_window = max(rate_smooth(valid_idx));
    
    if peak_in_window < thresh_high
        duration_ms = 0; % Peak didn't reach 3 SD
        return;
    end
    
    % 2. DURATION MEASUREMENT
    sig_bins = valid_idx(rate_smooth(valid_idx) > thresh_low);
    
    if isempty(sig_bins)
        duration_ms = 0;
    else
        t_start = t_centers(sig_bins(1));
        t_end   = t_centers(sig_bins(end));
        duration_ms = t_end - t_start;
        if duration_ms < 2, duration_ms = 0; end
    end
end

function sp = load_ssd_spikes(folder)
    cd(folder); f=dir('*sp_xia_SSD.mat'); if isempty(f), error('No SSD file'); end
    S=load(f(1).name); if isfield(S,'sp_corr'), sp=S.sp_corr; elseif isfield(S,'sp_SSD'), sp=S.sp_SSD; else, sp=S.sp_in; end
end

function plot_shaded_error(x, y, se, col)
    if numel(x)<2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid=~isnan(y); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end