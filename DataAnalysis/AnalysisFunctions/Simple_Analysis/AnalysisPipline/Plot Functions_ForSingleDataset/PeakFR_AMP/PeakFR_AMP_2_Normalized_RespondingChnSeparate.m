%% ============================================================
%   Peak Firing Rate (Normalized) vs Amplitude
%   CHECK: Condition-Specific Responding Channels (No Union)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim4_251125_152849';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq4_5ms_251125_154235';
Electrode_Type = 1;

% Post-stim window for PSTH / peak FR (ms)
post_win_ms = [2 20];
FS = 30000;            % Hz

% PSTH parameters
bin_ms     = 1;
sigma_bins = 2;        % Gaussian smoothing kernel (in bins)

% Baseline window (relative to pulse 1 trigger, ms)
baseline_win_ms = [-90 -10];
baselineDur_s   = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;

% Modified Z-score parameters
k_MAD   = 1.4826;                % scaling constant
epsilon = 1 / baselineDur_s;     % soft noise floor
MAX_Z   = 50;                    % Artifact ceiling

% ---- PSTH bins and smoothing kernel ----
edges = post_win_ms(1):bin_ms:post_win_ms(2);
bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g = gausswin(kernel_size);
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
nTr_seq   = Sseq.n_Trials;

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

%% ================= DETERMINE RESPONDING CHANNELS (SEPARATE) =====
d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

% 1. Find Simultaneous Responders
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
resp_channels_SIM = find(resp_sim); % ONLY Sim responders

% 2. Find Sequential Responders
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
resp_channels_SEQ = find(resp_seq); % ONLY Seq responders

fprintf('--- SEPARATE LISTS CHECK ---\n');
fprintf('Simultaneous Responders: %d channels\n', length(resp_channels_SIM));
fprintf('Sequential Responders:   %d channels\n', length(resp_channels_SEQ));


%% =================== COMPUTE PEAK FR (SEPARATE LOOPS) =======

% --- 1. SIMULTANEOUS CALCULATION ---
NormPeakFR_sim = nan(length(resp_channels_SIM), length(Amps_sim));

for ci = 1:length(resp_channels_SIM)
    ch = resp_channels_SIM(ci);
    recCh = d(ch); 
    
    S_ch = sp_sim{recCh};
    for ai = 1:length(Amps_sim)
        tr_ids = find(ampIdx_sim == ai);
        if isempty(tr_ids), continue; end
        
        [peak_vals, FR_baseline] = get_trial_data(tr_ids, trig_sim, S_ch, ...
            baseline_win_ms, baselineDur_s, post_win_ms, edges, g, bin_s, FS);
            
        NormPeakFR_sim(ci,ai) = normalize_data(peak_vals, FR_baseline, k_MAD, epsilon, MAX_Z);
    end
end

% --- 2. SEQUENTIAL CALCULATION ---
% Note: Using resp_channels_SEQ list here
NormPeakFR_seq = nan(length(resp_channels_SEQ), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels_SEQ)
    ch = resp_channels_SEQ(ci);
    recCh = d(ch);
    
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


%% ===================== AVERAGE ACROSS CHANNELS ======================

% 1. SIMULTANEOUS POPULATION
AvgSim_Pop     = mean(NormPeakFR_sim, 1, 'omitnan');
AvgSim_Pop_SEM = std(NormPeakFR_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(NormPeakFR_sim), 1));


% 2. SEQUENTIAL POPULATION (Combined Sets)
% First: Average the Sets (Orders) within each channel
NormPeakFR_seq_Combined = mean(NormPeakFR_seq, 3, 'omitnan');

% Second: Average across Channels
AvgSeq_Pop     = squeeze(mean(NormPeakFR_seq_Combined, 1, 'omitnan'));
AvgSeq_Pop_SEM = squeeze(std(NormPeakFR_seq_Combined, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(NormPeakFR_seq_Combined), 1)));


%% ===================== PLOT AVERAGE POPULATION CURVES ======================
figure('Color','w', 'Position',[200 200 750 520]); hold on;

% --- Plot Simultaneous ---
sim_col = [0 0.3 0.8];
x = Amps_sim;
y = AvgSim_Pop;
se = AvgSim_Pop_SEM;

if any(~isnan(y))
    plot_shaded_error(x, y, se, sim_col);
    plot(x, y, '-o', 'Color', sim_col, 'LineWidth', 2, ...
        'MarkerFaceColor', sim_col, 'DisplayName', sprintf('Simultaneous (N=%d)', length(resp_channels_SIM)));
end

% --- Plot Sequential (Combined Order) ---
seq_cols = cool(numel(PTDs_ms)+2); 

for p = 1:numel(PTDs_ms)
    x = Amps_seq;
    
    if numel(PTDs_ms) == 1
        y = AvgSeq_Pop(:)'; 
        se = AvgSeq_Pop_SEM(:)';
        col = [0.85 0.33 0.1]; 
    else
        y = AvgSeq_Pop(:,p)';
        se = AvgSeq_Pop_SEM(:,p)';
        col = seq_cols(p+1,:);
    end
    
    if any(~isnan(y))
        lbl = sprintf('Sequential (PTD %.1f ms | N=%d)', PTDs_ms(p), length(resp_channels_SEQ));
        plot_shaded_error(x, y, se, col);
        plot(x, y, '-s', 'Color', col, 'LineWidth', 2, ...
            'MarkerFaceColor', 'w', 'DisplayName', lbl);
    end
end

xlabel('Amplitude (µA)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Peak FR (Z_{robust})', 'FontSize', 12, 'FontWeight', 'bold');
title('CHECK: Condition-Specific Responding Channels', 'FontWeight', 'bold');
legend('Location','best', 'Box', 'off');
grid on; box off;


%% ==================== HELPER FUNCTIONS =========================
function [peak_vals, FR_baseline] = get_trial_data(tr_ids, trig, sp_data, ...
    base_win, base_dur, post_win, edges, kern, bin_s, FS)
    % Helper to extract raw peak and baseline values for list of trials
    nTr = numel(tr_ids);
    peak_vals = zeros(nTr,1);
    FR_baseline = zeros(nTr,1);
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        % Baseline
        mask_b = tt >= base_win(1) & tt < base_win(2);
        FR_baseline(k) = sum(mask_b)/base_dur;
        
        % Peak
        tt_post = tt(tt >= post_win(1) & tt <= post_win(2));
        if isempty(tt_post)
            peak_vals(k) = 0;
        else
            h = histcounts(tt_post, edges);
            rate_curve = conv(h/bin_s, kern, 'same');
            peak_vals(k) = max(rate_curve);
        end
    end
end

function Z_mean = normalize_data(peak_vals, baseline_FR, k_MAD, epsilon, max_z)
    % Calculates trial-averaged Robust Z-score
    medB = median(baseline_FR);
    MADB = median(abs(baseline_FR - medB));
    denom = k_MAD * MADB + epsilon;
    
    Z_trials = (peak_vals - medB) ./ denom;
    
    % Artifact Check
    Z_trials(abs(Z_trials) > max_z) = NaN;
    
    Z_mean = mean(Z_trials, 'omitnan');
end

function sp = load_ssd_spikes(folder)
    cd(folder);
    f = dir('*sp_xia_SSD.mat');
    if isempty(f), error('No sp_xia_SSD.mat in %s', folder); end
    S = load(f(1).name);
    if isfield(S,'sp_corr'), sp = S.sp_corr;
    elseif isfield(S,'sp_SSD'), sp = S.sp_SSD;
    elseif isfield(S,'sp_in'),  sp = S.sp_in;
    else, error('SSD file missing variables.'); end
end

function plot_shaded_error(x, y, se, col)
    if numel(x) < 2, return; end
    % Ensure rows
    x=x(:)'; y=y(:)'; se=se(:)';
    valid = ~isnan(y);
    x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    
    upper = y + se;
    lower = y - se;
    fill([x fliplr(x)], [upper fliplr(lower)], ...
        col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end