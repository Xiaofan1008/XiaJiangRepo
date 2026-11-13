%% ============================================================
%   Input-Output Curve: MEAN FIRING RATE (Total Spikes)
%   - Separate Curves for Each Sequential Set
%   - Legend includes Stimulation Channel Numbers
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type = 1;

% ANALYSIS WINDOW (Integration Time for AUC)
post_win_ms = [2 20];  
window_dur_s = (post_win_ms(2) - post_win_ms(1)) / 1000;

FS = 30000;            % Hz

% Baseline window
baseline_win_ms = [-90 -10];
baselineDur_s   = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;

% Modified Z-score parameters
k_MAD   = 1.4826;                
epsilon = 1 / baselineDur_s;     
MAX_Z   = 50;                    

% Plotting Aesthetics
jitter_width = 0.2;     
scatter_alpha = 0.35;    
dot_size = 20;          

%% =================== LOAD DATA ====================
% --- Load Simultaneous ---
cd(folder_sim);
Rsim_file = dir('*RespondingChannels.mat'); Rsim = load(Rsim_file(1).name).Responding;
sp_sim = load_ssd_spikes(folder_sim);   
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig_sim = loadTrig(0);
Ssim = load(dir('*_exp_datafile_*.mat').name); Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim;
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

% --- Load Sequential ---
cd(folder_seq);
Rseq_file = dir('*RespondingChannels.mat'); Rseq = load(Rseq_file(1).name).Responding;
sp_seq = load_ssd_spikes(folder_seq);
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig_seq = loadTrig(0);
Sseq = load(dir('*_exp_datafile_*.mat').name); Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP; 

% ROBUST n_Trials Check
if isfield(Sseq, 'n_Trials')
    nTr_seq = Sseq.n_Trials;
else
    nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq;
end

amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

% --- Parse Sets (Identify Channel Orders) ---
stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v  = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable');
nSets_seq = size(uniqueComb_seq,1);

%% ================= DETERMINE RESPONDING CHANNELS =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_sim = false(nCh_Total, 1);
for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1);
for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end
resp_channels = find(resp_sim | resp_seq); 

%% =================== COMPUTE MEAN FR (Total Spikes) ==========
NormFR_sim = nan(length(resp_channels), length(Amps_sim));
NormFR_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch = resp_channels(ci); recCh = d(ch);
    
    % --- 1. SIMULTANEOUS ---
    S_ch = sp_sim{recCh};
    for ai = 1:length(Amps_sim)
        tr_ids = find(ampIdx_sim == ai);
        if isempty(tr_ids), continue; end
        
        [mean_rate, FR_baseline] = get_mean_rate(tr_ids, trig_sim, S_ch, ...
            baseline_win_ms, baselineDur_s, post_win_ms, window_dur_s, FS);
        NormFR_sim(ci,ai) = normalize_data(mean_rate, FR_baseline, k_MAD, epsilon, MAX_Z);
    end
    
    % --- 2. SEQUENTIAL ---
    S_ch = sp_seq{recCh};
    for p = 1:numel(PTDs_ms)
        for ss = 1:nSets_seq
            for ai = 1:length(Amps_seq)
                tr_ids = find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai);
                if isempty(tr_ids), continue; end
                
                [mean_rate, FR_baseline] = get_mean_rate(tr_ids, trig_seq, S_ch, ...
                    baseline_win_ms, baselineDur_s, post_win_ms, window_dur_s, FS);
                NormFR_seq(ci,ai,ss,p) = normalize_data(mean_rate, FR_baseline, k_MAD, epsilon, MAX_Z);
            end
        end
    end
end 

%% ===================== AVERAGE & PLOT ======================
figure('Color','w', 'Position',[200 200 800 600]); hold on;

% --- 1. Plot Simultaneous ---
AvgSim_Pop     = mean(NormFR_sim, 1, 'omitnan');
AvgSim_Pop_SEM = std(NormFR_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(NormFR_sim), 1));
sim_col = [0 0.3 0.8]; 

x = Amps_sim; y = AvgSim_Pop; se = AvgSim_Pop_SEM;
if any(~isnan(y))
    % Jittered Scatter
    for i = 1:length(x)
        x_jit = (rand(size(NormFR_sim(:,i))) - 0.5) * jitter_width + x(i);
        scatter(x_jit, NormFR_sim(:,i), dot_size, sim_col, 'filled', 'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
    end
    plot_shaded_error(x, y, se, sim_col);
    plot(x, y, '-o', 'Color', sim_col, 'LineWidth', 2, 'MarkerFaceColor', sim_col, 'DisplayName', 'Simultaneous');
end

% --- 2. Plot Sequential (SEPARATE SETS) ---
% Color Scheme for sets (Warm colors)
set_colors = [0.85 0.33 0.10;  % Orange
              0.60 0.20 0.60;  % Purple
              0.20 0.60 0.20]; % Green (if needed)

for p = 1:numel(PTDs_ms)
    for ss = 1:nSets_seq
        
        % Extract Data for this specific Set & PTD
        % Data Shape: [Channels x Amps]
        data_this_set = squeeze(NormFR_seq(:, :, ss, p));
        
        % Average across channels
        AvgSeq_Set     = mean(data_this_set, 1, 'omitnan');
        AvgSeq_Set_SEM = std(data_this_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_this_set), 1));
        
        x = Amps_seq;
        y = AvgSeq_Set;
        se = AvgSeq_Set_SEM;
        
        % Select Color (Cycle if > 3 sets)
        col = set_colors(mod(ss-1, size(set_colors,1))+1, :);
        
        if any(~isnan(y))
            % --- Generate Legend Label with Channels ---
            stimCh_List = uniqueComb_seq(ss, :);
            stimCh_List = stimCh_List(stimCh_List > 0);
            stimCh_Str  = num2str(stimCh_List);
            
            lbl = sprintf('Seq Set %d (Ch: %s)', ss, stimCh_Str);
            if numel(PTDs_ms) > 1
                lbl = sprintf('Seq Set %d (PTD %.1f) (Ch: %s)', ss, PTDs_ms(p), stimCh_Str);
            end
            
            % Jittered Scatter
            for i = 1:length(x)
                x_jit = (rand(size(data_this_set(:,i))) - 0.5) * jitter_width + x(i);
                scatter(x_jit, data_this_set(:,i), dot_size, col, 'filled', 'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
            end
            
            plot_shaded_error(x, y, se, col);
            plot(x, y, '-s', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', lbl);
        end
    end
end

xlabel('Amplitude (µA)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Normalized Mean Firing Rate (sps)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Total Activation (Mean Rate 2-20ms)\n %d Channels', length(resp_channels)), 'FontWeight', 'bold');
legend('Location','best', 'Box', 'off'); box off;

%% ==================== HELPER FUNCTIONS =========================
function [mean_rate, baseline_rates] = get_mean_rate(tr_ids, trig, sp_data, ...
    base_win, base_dur, post_win, win_dur, FS)
    
    nTr = numel(tr_ids);
    baseline_rates = zeros(nTr, 1);
    total_post_spikes = 0; 
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        % Baseline Rate (Hz) per trial
        mask_b = tt >= base_win(1) & tt < base_win(2);
        baseline_rates(k) = sum(mask_b) / base_dur;
        
        % Post-stim Count
        total_post_spikes = total_post_spikes + sum(tt >= post_win(1) & tt <= post_win(2));
    end
    
    % Mean Firing Rate (Hz) = Total Spikes / Total Time Monitored
    if nTr > 0
        mean_rate = total_post_spikes / (nTr * win_dur);
    else
        mean_rate = 0;
    end
end

function Z_val = normalize_data(val, baseline_rates, k_MAD, epsilon, max_z)
    medB = median(baseline_rates);
    MADB = median(abs(baseline_rates - medB));
    denom = k_MAD * MADB + epsilon;
    Z_val = (val - medB) / denom;
    if abs(Z_val) > max_z, Z_val = NaN; end
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