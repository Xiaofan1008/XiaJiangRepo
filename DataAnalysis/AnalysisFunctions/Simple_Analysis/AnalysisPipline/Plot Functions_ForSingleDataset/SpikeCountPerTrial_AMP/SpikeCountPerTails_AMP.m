%% ============================================================
%   Input-Output Curve: SPIKE COUNT PER TRIAL (Total Activation)
%   - AVERAGED Sequential Sets (One Curve for Seq, One for Sim)
%   - Metric: Raw Spike Count (No Normalization)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type = 1;

% ANALYSIS WINDOW (Integration Time)
post_win_ms = [2 20];  
FS = 30000;            % Hz

% Baseline window (kept for ref, not used for normalization here)
baseline_win_ms = [-90 -10];
baselineDur_s   = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;

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

%% =================== COMPUTE SPIKE COUNT (Per Trial) ==========
% NO Z-SCORING. JUST RAW COUNTS.
Count_sim = nan(length(resp_channels), length(Amps_sim));
Count_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch = resp_channels(ci); recCh = d(ch);
    
    % --- 1. SIMULTANEOUS ---
    S_ch = sp_sim{recCh};
    for ai = 1:length(Amps_sim)
        tr_ids = find(ampIdx_sim == ai);
        if isempty(tr_ids), continue; end
        
        % Get Average Spike Count per Trial
        Count_sim(ci,ai) = get_spike_count(tr_ids, trig_sim, S_ch, post_win_ms, FS);
    end
    
    % --- 2. SEQUENTIAL ---
    S_ch = sp_seq{recCh};
    for p = 1:numel(PTDs_ms)
        for ss = 1:nSets_seq
            for ai = 1:length(Amps_seq)
                tr_ids = find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai);
                if isempty(tr_ids), continue; end
                
                Count_seq(ci,ai,ss,p) = get_spike_count(tr_ids, trig_seq, S_ch, post_win_ms, FS);
            end
        end
    end
end 

%% ===================== AVERAGE & PLOT ======================
figure('Color','w', 'Position',[200 200 800 600]); hold on;

% --- 1. Plot Simultaneous ---
AvgSim_Pop     = mean(Count_sim, 1, 'omitnan');
AvgSim_Pop_SEM = std(Count_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Count_sim), 1));
sim_col = [0 0.3 0.8]; 

x = Amps_sim; y = AvgSim_Pop; se = AvgSim_Pop_SEM;
if any(~isnan(y))
    % Jittered Scatter
    for i = 1:length(x)
        x_jit = (rand(size(Count_sim(:,i))) - 0.5) * jitter_width + x(i);
        scatter(x_jit, Count_sim(:,i), dot_size, sim_col, 'filled', 'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
    end
    plot_shaded_error(x, y, se, sim_col);
    plot(x, y, '-o', 'Color', sim_col, 'LineWidth', 2, 'MarkerFaceColor', sim_col, 'DisplayName', 'Simultaneous');
end

% --- 2. Plot Sequential (AVERAGED ACROSS SETS) ---
% First, average across Sets (Dim 3) to combine Set 1, Set 2, etc.
% Result: [Channels x Amps x 1 x PTDs]
Count_seq_Combined = mean(Count_seq, 3, 'omitnan');

% Color: Orange
seq_col = [0.85 0.33 0.10]; 

for p = 1:numel(PTDs_ms)
    
    % Extract data for this PTD (now averaged across sets)
    % Shape: [Channels x Amps]
    data_this_ptd = squeeze(Count_seq_Combined(:, :, 1, p));
    
    AvgSeq_Pop     = mean(data_this_ptd, 1, 'omitnan');
    AvgSeq_Pop_SEM = std(data_this_ptd, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_this_ptd), 1));
    
    x = Amps_seq; y = AvgSeq_Pop; se = AvgSeq_Pop_SEM;
    
    if any(~isnan(y))
        lbl = 'Sequential (Avg Sets)';
        if numel(PTDs_ms) > 1
            lbl = sprintf('Sequential (PTD %.1f ms)', PTDs_ms(p));
        end
        
        % Jittered Scatter (Plotting the set-averaged channel responses)
        for i = 1:length(x)
            x_jit = (rand(size(data_this_ptd(:,i))) - 0.5) * jitter_width + x(i);
            scatter(x_jit, data_this_ptd(:,i), dot_size, seq_col, 'filled', 'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
        end
        
        plot_shaded_error(x, y, se, seq_col);
        plot(x, y, '-s', 'Color', seq_col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', lbl);
    end
end

xlabel('Amplitude (µA)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Mean Spikes per Trial (Total Count)', 'FontSize', 10, 'FontWeight', 'bold');
title(sprintf('Total Activation (Window: %d-%d ms)\n %d Channels', ...
    post_win_ms(1), post_win_ms(2), length(resp_channels)), 'FontWeight', 'bold');
legend('Location','best', 'Box', 'off'); box off;

%% ==================== HELPER FUNCTIONS =========================
function avg_count = get_spike_count(tr_ids, trig, sp_data, post_win, FS)
    % Calculates the AVERAGE number of spikes per trial in the window
    nTr = numel(tr_ids);
    total_spikes = 0; 
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        % Strictly Count
        count_post = sum(tt >= post_win(1) & tt <= post_win(2));
        total_spikes = total_spikes + count_post;
    end
    
    if nTr > 0
        avg_count = total_spikes / nTr; 
    else
        avg_count = 0;
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