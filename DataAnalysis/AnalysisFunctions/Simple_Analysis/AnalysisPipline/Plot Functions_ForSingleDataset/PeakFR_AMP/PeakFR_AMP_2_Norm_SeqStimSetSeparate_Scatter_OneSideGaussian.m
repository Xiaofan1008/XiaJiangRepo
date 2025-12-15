%% ============================================================
%   Peak Firing Rate (Normalized) vs Amplitude
%   Union of Responding Channels
%   Simultaneous vs Sequential (SEPARATE SETS)
%   + SCATTER PLOTS
%   + ONE-SIDED (CAUSAL) SMOOTHING KERNEL
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type = 1;

% Post-stim window for PSTH / peak FR (ms)
post_win_ms = [2 20];
FS = 30000;            % Hz

% PSTH parameters
bin_ms     = 1;
sigma_bins = 3;        % Gaussian smoothing kernel (in bins)

% Baseline window (relative to pulse 1 trigger, ms)
baseline_win_ms = [-90 -10];
baselineDur_s   = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;

% Modified Z-score parameters
k_MAD   = 1.4826;                % scaling constant
epsilon = 1 / baselineDur_s;     % soft noise floor
MAX_Z   = 50;                    % Artifact ceiling

% Plotting Aesthetics
jitter_width = 0.2;     
scatter_alpha = 0.4;    
dot_size = 20;          

% ---- PSTH bins and smoothing kernel (MODIFIED: ONE-SIDED) ----
edges = post_win_ms(1):bin_ms:post_win_ms(2);
bin_s = bin_ms/1000;

% Causal (One-Sided) Kernel
% Start at 0, decay over 3 sigmas
kernel_len = ceil(3 * sigma_bins); 
g = exp(-0.5 * ((0:kernel_len-1) / sigma_bins).^2);
g = g / sum(g);

%% =================== LOAD RESPONDING CHANNELS ================
cd(folder_sim);
Rsim_file = dir('*RespondingChannels.mat');
if isempty(Rsim_file), error('No RespondingChannels.mat in simultaneous folder'); end
Rsim = load(Rsim_file(1).name).Responding;

cd(folder_seq);
Rseq_file = dir('*RespondingChannels.mat');
if isempty(Rseq_file), error('No RespondingChannels.mat in sequential folder'); end
Rseq = load(Rseq_file(1).name).Responding;

%% =================== LOAD SIM SPIKES & STIMPARAMS ============
cd(folder_sim);
sp_sim = load_ssd_spikes(folder_sim);   
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_sim = loadTrig(0);
Ssim = load(dir('*_exp_datafile_*.mat').name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');

Stim_sim = Ssim.StimParams;
simN_sim = Ssim.simultaneous_stim;
E_MAP_sim = Ssim.E_MAP;
nTr_sim   = Ssim.n_Trials;

amps_all_sim  = cell2mat(Stim_sim(2:end,16));
trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim);
Amps_sim(Amps_sim==-1) = 0;

%% =================== LOAD SEQ SPIKES & STIMPARAMS ============
cd(folder_seq);
sp_seq = load_ssd_spikes(folder_seq);
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_seq = loadTrig(0);
Sseq = load(dir('*_exp_datafile_*.mat').name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');

Stim_seq = Sseq.StimParams;
simN_seq = Sseq.simultaneous_stim;
E_MAP_seq = Sseq.E_MAP;

% ROBUST n_Trials Check
if isfield(Sseq, 'n_Trials')
    nTr_seq = Sseq.n_Trials;
else
    nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq;
end

amps_all_seq  = cell2mat(Stim_seq(2:end,16));
trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq);
Amps_seq(Amps_seq==-1) = 0;

PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6));
PTD_all_ms = PTD_all_us / 1000;
[PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

% Parse Sets (Orders)
stimNames = Stim_seq(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq
    rr = (t-1)*simN_seq + (1:simN_seq);
    v  = idx_all(rr); v = v(v>0);
    comb_seq(t,1:numel(v)) = v(:).';
end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable');
nSets_seq = size(uniqueComb_seq,1);

%% ================= DETERMINE RESPONDING CHANNELS (UNION) =====
d = Depth_s(Electrode_Type); % moved here to get channel count
nCh_Total = length(d);

resp_sim = false(nCh_Total, 1);
for si = 1:numel(Rsim.set)
    for ai = 1:numel(Rsim.set(si).amp)
        for pi = 1:numel(Rsim.set(si).amp(ai).ptd)
            this = Rsim.set(si).amp(ai).ptd(pi).channel;
            for ch = 1:min(length(this), nCh_Total)
                if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                    resp_sim(ch) = true;
                end
            end
        end
    end
end

resp_seq = false(nCh_Total, 1);
for si = 1:numel(Rseq.set)
    for ai = 1:numel(Rseq.set(si).amp)
        for pi = 1:numel(Rseq.set(si).amp(ai).ptd)
            this = Rseq.set(si).amp(ai).ptd(pi).channel;
            for ch = 1:min(length(this), nCh_Total)
                if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                    resp_seq(ch) = true;
                end
            end
        end
    end
end

resp_channels = find(resp_sim | resp_seq); 
fprintf('Found %d Responding Channels (Union of Sim/Seq)\n', length(resp_channels));

%% =================== COMPUTE PEAK + NORMALIZED PEAK FR =======
% Storage Arrays
NormPeakFR_sim = nan(length(resp_channels), length(Amps_sim));
NormPeakFR_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch = resp_channels(ci);
    recCh = d(ch); % The actual electrode ID    
    
    % --- 1. SIMULTANEOUS STIMULATION ---
    S_ch = sp_sim{recCh};
    for ai = 1:length(Amps_sim)
        tr_ids = find(ampIdx_sim == ai);
        if isempty(tr_ids), continue; end
        
        [peak_vals, FR_baseline] = get_trial_data(tr_ids, trig_sim, S_ch, ...
            baseline_win_ms, baselineDur_s, post_win_ms, edges, g, bin_s, FS);
            
        NormPeakFR_sim(ci,ai) = normalize_data(peak_vals, FR_baseline, k_MAD, epsilon, MAX_Z);
    end
    
    % --- 2. SEQUENTIAL STIMULATION ---
    S_ch = sp_seq{recCh};
    for p = 1:numel(PTDs_ms)
        for ss = 1:nSets_seq
            for ai = 1:length(Amps_seq)
                tr_ids = find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai);
                if isempty(tr_ids), continue; end
                
                [peak_vals, FR_baseline] = get_trial_data(tr_ids, trig_seq, S_ch, ...
                    baseline_win_ms, baselineDur_s, post_win_ms, edges, g, bin_s, FS);
                
                NormPeakFR_seq(ci,ai,ss,p) = normalize_data(peak_vals, FR_baseline, k_MAD, epsilon, MAX_Z);
            end
        end
    end
end 

%% ===================== AVERAGE ACROSS CHANNELS & PLOT ======================
figure('Color','w', 'Position',[200 200 800 600]); hold on;

% --- 1. Plot Simultaneous ---
AvgSim_Pop     = mean(NormPeakFR_sim, 1, 'omitnan');
AvgSim_Pop_SEM = std(NormPeakFR_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(NormPeakFR_sim), 1));
sim_col = [0 0.3 0.8];

x = Amps_sim; y = AvgSim_Pop; se = AvgSim_Pop_SEM;

if any(~isnan(y))
    % A. Scatter Plot (Individual Channels)
    for i = 1:length(x)
        amp_val = x(i);
        chan_vals = NormPeakFR_sim(:, i);
        x_jittered = (rand(size(chan_vals)) - 0.5) * jitter_width + amp_val;
        scatter(x_jittered, chan_vals, dot_size, sim_col, 'filled', ...
            'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility', 'off');
    end
    % B. Line Plot
    plot_shaded_error(x, y, se, sim_col);
    plot(x, y, '-o', 'Color', sim_col, 'LineWidth', 2, ...
        'MarkerFaceColor', sim_col, 'DisplayName', 'Simultaneous');
end

% --- 2. Plot Sequential (SEPARATE SETS) ---
% Color Palette for different sets
set_colors = [0.85 0.33 0.10;  % Orange
              0.60 0.20 0.60;  % Purple
              0.20 0.60 0.20]; % Green

for p = 1:numel(PTDs_ms)
    for ss = 1:nSets_seq
        
        % Extract data for this specific Set
        % Data Shape: [Channels x Amps]
        data_this_set = squeeze(NormPeakFR_seq(:, :, ss, p));
        
        AvgSeq_Set     = mean(data_this_set, 1, 'omitnan');
        AvgSeq_Set_SEM = std(data_this_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_this_set), 1));
        
        x = Amps_seq; y = AvgSeq_Set; se = AvgSeq_Set_SEM;
        col = set_colors(mod(ss-1, size(set_colors,1))+1, :);
        
        if any(~isnan(y))
            % --- Create Legend Label with Channels ---
            stimCh_List = uniqueComb_seq(ss, :);
            stimCh_List = stimCh_List(stimCh_List > 0);
            stimCh_Str  = num2str(stimCh_List);
            
            lbl = sprintf('Seq Set %d (Ch: %s)', ss, stimCh_Str);
            if numel(PTDs_ms) > 1
                lbl = sprintf('Seq Set %d (PTD %.1f) (Ch: %s)', ss, PTDs_ms(p), stimCh_Str);
            end
            
            % A. Scatter Plot (Individual Channels for this Set)
            for i = 1:length(x)
                amp_val = x(i);
                chan_vals = data_this_set(:, i);
                x_jittered = (rand(size(chan_vals)) - 0.5) * jitter_width + amp_val;
                scatter(x_jittered, chan_vals, dot_size, col, 'filled', ...
                    'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility', 'off');
            end
            
            % B. Line Plot
            plot_shaded_error(x, y, se, col);
            plot(x, y, '-s', 'Color', col, 'LineWidth', 2, ...
                'MarkerFaceColor', 'w', 'DisplayName', lbl);
        end
    end
end

xlabel('Amplitude (µA)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Normalized Peak Firing Rate (sps)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Normalized Peak Firing Rate (%d Channels)', length(resp_channels)), 'FontWeight', 'bold');
legend('Location','best', 'Box', 'off');
box off;

%% ==================== HELPER FUNCTIONS =========================
function [peak_val, baseline_rates] = get_trial_data(tr_ids, trig, sp_data, ...
    base_win, base_dur, post_win, edges, kern, bin_s, FS)
    
    nTr = numel(tr_ids);
    baseline_rates = zeros(nTr, 1);
    all_post_spikes = []; 
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        mask_b = tt >= base_win(1) & tt < base_win(2);
        baseline_rates(k) = sum(mask_b) / base_dur;
        
        tt_post = tt(tt >= post_win(1) & tt <= post_win(2));
        all_post_spikes = [all_post_spikes; tt_post]; 
    end
    
    if isempty(all_post_spikes)
        peak_val = 0;
    else
        h = histcounts(all_post_spikes, edges);
        % Divide by (N_Trials * Bin_Size_sec) to get Hz
        rate_curve = h / (nTr * bin_s);
        
        % Filter with the One-Sided Kernel
        rate_smooth = conv(rate_curve, kern, 'same');
        
        % For one-sided kernels, 'same' output might shift the peak.
        % If accurate latency is needed, sometimes 'full' and trimming is better,
        % but 'same' is standard if the kernel is short.
        
        peak_val = max(rate_smooth);
    end
end

function Z_val = normalize_data(peak_val, baseline_rates, k_MAD, epsilon, max_z)
    medB = median(baseline_rates);
    MADB = median(abs(baseline_rates - medB));
    denom = k_MAD * MADB + epsilon;
    Z_val = (peak_val - medB) / denom;
    if abs(Z_val) > max_z, Z_val = NaN; end
end

function sp = load_ssd_spikes(folder)
    cd(folder); f=dir('*sp_xia_SSD.mat'); if isempty(f), error('No SSD file'); end
    S=load(f(1).name); if isfield(S,'sp_corr'), sp=S.sp_corr; elseif isfield(S,'sp_SSD'), sp=S.sp_SSD; else, sp=S.sp_in; end
end

function plot_shaded_error(x, y, se, col)
    if numel(x) < 2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid=~isnan(y); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    upper = y + se; lower = y - se;
    fill([x fliplr(x)], [upper fliplr(lower)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end