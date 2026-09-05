%% ============================================================
%   Peak Firing Rate (Normalized) + PSTH Storage
%   - Calculates: Normalized Peak FR (Z-score) & Full PSTH
%   - Output: Saves strictly to 'Result_PeakFiringRate_... .mat'
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim4_251125_152849';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq4_5ms_251125_154235';
Electrode_Type = 1;

% 1. Analysis Window (For calculating the SINGLE Peak Value)
post_win_ms = [2 20]; 

% 2. PSTH Storage Window (For saving the full time course)
psth_win_ms = [-50 100]; 

FS = 30000;            

% Smoothing
bin_ms     = 1;
sigma_bins = 3;        % Gaussian smoothing (bins)
baseline_win_ms = [-90 -10];
baselineDur_s   = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;

% Normalization
k_MAD   = 1.4826;              
epsilon = 1 / baselineDur_s;    
MAX_Z   = 50;                  

% Plotting Aesthetics
jitter_width = 0.2; scatter_alpha = 0.4; dot_size = 20;          

% ---- Kernels ----
% A. For Peak Detection (Narrow window, Causal/One-sided is often better, but keeping Gauss)
edges_peak = post_win_ms(1):bin_ms:post_win_ms(2);
bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g_peak = gausswin(kernel_size); g_peak = g_peak / sum(g_peak);

% B. For PSTH Storage (Wide window)
edges_psth = psth_win_ms(1):bin_ms:psth_win_ms(2);
time_vector_psth = edges_psth(1:end-1) + bin_ms/2;
g_sym = gausswin(kernel_size); g_sym = g_sym / sum(g_sym);

%% =================== 1. LOAD DATA & QC INFO ====================
% Helper to load data and QC files (defined at bottom)
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% Extract Stim Params (Sim)
Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim;
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

% Extract Stim Params (Seq)
Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP;
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

% Parse Seq Sets
stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

% Extract Sim Channels for Metadata
stim_chans_sim = Rsim.set(1).stimChannels; 

%% ================= 2. IDENTIFY RESPONDING CHANNELS =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_sim = false(nCh_Total, 1);
for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1);
for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end
resp_channels = find(resp_sim | resp_seq); 
fprintf('Analyzing %d Responding Channels\n', length(resp_channels));

%% =================== 3. COMPUTE PEAK FR & PSTH =================
% Init Output Arrays
NormPeakFR_sim = nan(length(resp_channels), length(Amps_sim));
NormPeakFR_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

% PSTH Matrices: [Channels x TimeBins x Amps]
nBins = length(edges_psth)-1;
PSTH_Sim = zeros(length(resp_channels), nBins, length(Amps_sim));
PSTH_Seq = zeros(length(resp_channels), nBins, length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch_idx = resp_channels(ci); recCh = d(ch_idx);    
    
    % --- SIMULTANEOUS ---
    S_ch = sp_sim{recCh};
    
    % Check if BadTrials exists for this channel index
    bad_trs = [];
    if ~isempty(QC_Sim.BadTrials) && ch_idx <= length(QC_Sim.BadTrials)
        bad_trs = QC_Sim.BadTrials{ch_idx}; 
    end
    
    for ai = 1:length(Amps_sim)
        tr_ids = find(ampIdx_sim == ai);
        tr_ids = setdiff(tr_ids, bad_trs); % QC Filter
        
        if isempty(tr_ids), continue; end
        
        [peak_val, FR_baseline, psth_curve] = get_trial_metrics(tr_ids, trig_sim, S_ch, ...
            baseline_win_ms, baselineDur_s, post_win_ms, edges_peak, ...
            edges_psth, g_peak, g_sym, bin_s, FS);
            
        NormPeakFR_sim(ci,ai) = normalize_data(peak_val, FR_baseline, k_MAD, epsilon, MAX_Z);
        PSTH_Sim(ci, :, ai)   = psth_curve;
    end
    
    % --- SEQUENTIAL ---
    S_ch = sp_seq{recCh};
    bad_trs = [];
    if ~isempty(QC_Seq.BadTrials) && ch_idx <= length(QC_Seq.BadTrials)
        bad_trs = QC_Seq.BadTrials{ch_idx};
    end
    
    for p = 1:numel(PTDs_ms)
        for ss = 1:nSets_seq
            for ai = 1:length(Amps_seq)
                tr_ids = find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai);
                tr_ids = setdiff(tr_ids, bad_trs); % QC Filter
                
                if isempty(tr_ids), continue; end
                
                [peak_val, FR_baseline, psth_curve] = get_trial_metrics(tr_ids, trig_seq, S_ch, ...
                    baseline_win_ms, baselineDur_s, post_win_ms, edges_peak, ...
                    edges_psth, g_peak, g_sym, bin_s, FS);
                
                NormPeakFR_seq(ci,ai,ss,p) = normalize_data(peak_val, FR_baseline, k_MAD, epsilon, MAX_Z);
                PSTH_Seq(ci, :, ai, ss, p) = psth_curve;
            end
        end
    end
end 

%% ===================== 4. PLOT (Full Figure) ======================
figure('Color','w', 'Position',[100 100 800 600]); hold on;

% --- Plot Simultaneous ---
AvgSim = mean(NormPeakFR_sim, 1, 'omitnan');
SEMSim = std(NormPeakFR_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(NormPeakFR_sim), 1));
sim_col = [0 0.3 0.8];

if any(~isnan(AvgSim))
    for i = 1:length(Amps_sim)
        x_jit = (rand(size(NormPeakFR_sim(:,i))) - 0.5) * jitter_width + Amps_sim(i);
        scatter(x_jit, NormPeakFR_sim(:,i), dot_size, sim_col, 'filled', ...
            'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
    end
    plot_shaded_error(Amps_sim, AvgSim, SEMSim, sim_col);
    plot(Amps_sim, AvgSim, '-o', 'Color', sim_col, 'LineWidth', 2, ...
        'MarkerFaceColor', sim_col, 'DisplayName', 'Simultaneous');
end

% --- Plot Sequential Sets ---
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 

for p = 1:numel(PTDs_ms)
    for ss = 1:nSets_seq
        data_set = squeeze(NormPeakFR_seq(:, :, ss, p));
        AvgSeq = mean(data_set, 1, 'omitnan');
        SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
        col = set_colors(mod(ss-1,3)+1, :);
        
        stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
        lbl = sprintf('Seq Set %d (Ch:%s)', ss, num2str(stimCh));
        
        if any(~isnan(AvgSeq))
            for i = 1:length(Amps_seq)
                x_jit = (rand(size(data_set(:,i))) - 0.5) * jitter_width + Amps_seq(i);
                scatter(x_jit, data_set(:,i), dot_size, col, 'filled', ...
                    'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
            end
            plot_shaded_error(Amps_seq, AvgSeq, SEMSeq, col);
            plot(Amps_seq, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, ...
                'MarkerFaceColor', 'w', 'DisplayName', lbl);
        end
    end
end

xlabel('Amplitude (uA)', 'FontWeight','bold'); ylabel('Norm Peak FR (sps)', 'FontWeight','bold');
title('Peak Firing Rate vs Amplitude', 'FontWeight','bold');
legend('Location','best','Box','off'); box off;

%% ============================================================
%   5. SAVE RESULTS (Specific File for Firing Rate)
% ============================================================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

parts = split(folder_sim, filesep); 
exp_id = parts{end};
% File Name: Result_PeakFiringRate_[ExpID].mat
out_filename = fullfile(save_dir, ['Result_Set2_PeakFiringRate5_5ms_' exp_id '.mat']);

% Construct Data Structure
ResultFR = struct();
ResultFR.Metadata.Created = datestr(now);
ResultFR.Metadata.Source_Sim = folder_sim;
ResultFR.Metadata.Source_Seq = folder_seq;
ResultFR.Metadata.StimGroups.Sim = stim_chans_sim;
ResultFR.Metadata.StimGroups.Seq = uniqueComb_seq(:, uniqueComb_seq(1,:) > 0);
ResultFR.Metadata.RespChannels = resp_channels;
ResultFR.Metadata.Amps = Amps_sim;

% 1. Peak Firing Rate (Normalized)
ResultFR.PeakFR.Sim_Z = NormPeakFR_sim;
ResultFR.PeakFR.Seq_Z = NormPeakFR_seq;

% 2. Full PSTH Data
ResultFR.PSTH.Sim = PSTH_Sim;
ResultFR.PSTH.Seq = PSTH_Seq;
ResultFR.PSTH.TimeVector = time_vector_psth;
ResultFR.PSTH.Parameters.Window = psth_win_ms;
ResultFR.PSTH.Parameters.BinSize = bin_ms;

% 3. QC Data Used
ResultFR.QC.Sim = QC_Sim;
ResultFR.QC.Seq = QC_Seq;

% Save
save(out_filename, 'ResultFR');
fprintf('\n>>> Results Saved to: %s\n', out_filename);


%% ==================== HELPER FUNCTIONS =========================
function [peak_val, baseline_rates, psth_trace] = get_trial_metrics(tr_ids, trig, sp_data, ...
    base_win, base_dur, peak_win, peak_edges, psth_edges, g_peak, g_sym, bin_s, FS)
    
    nTr = numel(tr_ids);
    baseline_rates = zeros(nTr, 1);
    
    % Collect spikes for both windows
    all_peak_spikes = [];
    all_psth_spikes = [];
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        % Baseline Rate (per trial)
        mask_b = tt >= base_win(1) & tt < base_win(2);
        baseline_rates(k) = sum(mask_b) / base_dur;
        
        % Collect for Aggregate Peak (Narrow Window)
        all_peak_spikes = [all_peak_spikes; tt(tt >= peak_win(1) & tt <= peak_win(2))]; 
        
        % Collect for Full PSTH (Wide Window)
        all_psth_spikes = [all_psth_spikes; tt(tt >= psth_edges(1) & tt <= psth_edges(end))];
    end
    
    % A. Calculate Peak (Narrow Window)
    if isempty(all_peak_spikes)
        peak_val = 0;
    else
        h = histcounts(all_peak_spikes, peak_edges);
        rate = h / (nTr * bin_s);
        smooth = conv(rate, g_peak, 'same');
        peak_val = max(smooth);
    end
    
    % B. Calculate Full PSTH (Wide Window)
    h_psth = histcounts(all_psth_spikes, psth_edges);
    rate_psth = h_psth / (nTr * bin_s);
    psth_trace = conv(rate_psth, g_sym, 'same');
end

function Z_val = normalize_data(peak_val, baseline_rates, k_MAD, epsilon, max_z)
    medB = median(baseline_rates);
    MADB = median(abs(baseline_rates - medB));
    denom = k_MAD * MADB + epsilon;
    Z_val = (peak_val - medB) / denom;
    if abs(Z_val) > max_z, Z_val = NaN; end
end

function plot_shaded_error(x, y, se, col)
    if numel(x) < 2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid = ~isnan(y); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    upper = y + se; lower = y - se;
    fill([x fliplr(x)], [upper fliplr(lower)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    % Responding
    f = dir('*RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    
    % Spikes (Check SSD first)
    f = dir('*sp_xia_SSD.mat'); 
    if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr;
    elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD;
    else, sp = S_sp.sp_in; end
    
    % Triggers
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    
    % Params
    S = load(dir('*_exp_datafile_*.mat').name);
    
    % QC: Load BadChannels/Trials if they exist
    QC.BadCh = []; QC.BadTrials = [];
    
    f_bc = dir('*.BadChannels.mat');
    if ~isempty(f_bc)
        tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; 
    end
    
    f_bt = dir('*.BadTrials.mat');
    if ~isempty(f_bt)
        tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials;
    end
end