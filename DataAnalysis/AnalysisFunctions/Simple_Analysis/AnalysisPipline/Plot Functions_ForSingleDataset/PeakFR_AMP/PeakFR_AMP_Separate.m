%% ============================================================
%   Peak Firing Rate vs Amplitude — Only Responding Channels
%   Loads both Sim & Seq responding-channel files
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type = 1;

% The FR window in ms (should match your responding-channel detection)
post_win_ms = [2 20];    
FS = 30000;

% PSTH parameter
bin_ms = 1;
sigma_bins = 2; % Gaussian smoothing kernel
       
edges = post_win_ms(1):bin_ms:post_win_ms(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;        
% Gaussian smoothing kernel
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

% SSD spikes
sp_sim = load_ssd_spikes(folder_sim);

% Triggers
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_sim = loadTrig(0);

% StimParams
Ssim = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim_sim = Ssim.StimParams;
simN_sim  = Ssim.simultaneous_stim;
E_MAP_sim = Ssim.E_MAP;
nTr_sim   = Ssim.n_Trials;

amps_all_sim  = cell2mat(Stim_sim(2:end,16));
trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim);
Amps_sim(Amps_sim==-1) = 0;

stimNames = Stim_sim(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP_sim(2:end));
stimChPerTrial = cell(nTr_sim,1);
for t = 1:nTr_sim
    rr = (t-1)*simN_sim + (1:simN_sim);
    v  = idx_all(rr); v = v(v>0);
    stimChPerTrial{t} = v(:).';
end

comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim
    v = stimChPerTrial{t};
    comb_sim(t,1:numel(v)) = v;
end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable');
nSets_sim = size(uniqueComb_sim,1);

%% =================== LOAD SEQ SPIKES & STIMPARAMS ============
cd(folder_seq);

sp_seq = load_ssd_spikes(folder_seq);

if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_seq = loadTrig(0);

Sseq = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim_seq = Sseq.StimParams;
simN_seq  = Sseq.simultaneous_stim;
E_MAP_seq = Sseq.E_MAP;
nTr_seq   = Sseq.n_Trials;

amps_all_seq  = cell2mat(Stim_seq(2:end,16));
trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq);
Amps_seq(Amps_seq==-1) = 0;

PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6));
PTD_all_ms = PTD_all_us / 1000;
[PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

stimNames = Stim_seq(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));

stimChPerTrial = cell(nTr_seq,1);
for t = 1:nTr_seq
    rr = (t-1)*simN_seq + (1:simN_seq);
    v = idx_all(rr); v = v(v>0);
    stimChPerTrial{t} = v(:).';
end

comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq
    v = stimChPerTrial{t};
    comb_seq(t,1:numel(v)) = v;
end

[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable');
nSets_seq = size(uniqueComb_seq,1);

%% ================= DETERMINE RESPONDING CHANNELS =============
% Find all responding channels in SIM
resp_sim = false(32,1);
for si = 1:numel(Rsim.set)
    for ai = 1:numel(Rsim.set(si).amp)
        for pi = 1:numel(Rsim.set(si).amp(ai).ptd)
            this = Rsim.set(si).amp(ai).ptd(pi).channel;
            for ch = 1:length(this)
                if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                    resp_sim(ch) = true;
                end
            end
        end
    end
end

% Find all responding channels in SEQ
resp_seq = false(32,1);
for si = 1:numel(Rseq.set)
    for ai = 1:numel(Rseq.set(si).amp)
        for pi = 1:numel(Rseq.set(si).amp(ai).ptd)
            this = Rseq.set(si).amp(ai).ptd(pi).channel;
            for ch = 1:length(this)
                if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                    resp_seq(ch) = true;
                end
            end
        end
    end
end

% UNION of responding channels
resp_channels = find(resp_sim | resp_seq);
%% =================== COMPUTE PEAK FR ==========================

d = Depth_s(Electrode_Type);

% Output matrices
PeakFR_sim = nan(length(resp_channels), length(Amps_sim));
PeakFR_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)

    ch = resp_channels(ci);
    recCh_sim = d(ch);
    recCh_seq = d(ch);

    %% -------------------- SIMULTANEOUS ------------------------
    S_ch = sp_sim{recCh_sim};

    for ai = 1:length(Amps_sim)
        tr_ids = find(ampIdx_sim==ai);
        if isempty(tr_ids), continue; end
        % Build PSTH across all trials
        counts = zeros(1,length(edges)-1);    
        for k = 1:numel(tr_ids)
            t0 = trig_sim(tr_ids(k))/FS*1000;
            tt = S_ch(:,1) - t0;       % align spikes
            tt = tt(tt>=post_win_ms(1) & tt<=post_win_ms(2));        
            counts = counts + histcounts(tt, edges);
        end       
        rate_curve = counts / (numel(tr_ids)*bin_s);
        rate_curve = conv(rate_curve, g, "same");   % smooth PSTH
       
        PeakFR_sim(ci, ai) = max(rate_curve);
    end

    %% -------------------- SEQUENTIAL --------------------------
    S_ch = sp_seq{recCh_seq};
    for ss = 1:nSets_seq
        for p = 1:numel(PTDs_ms)
            for ai = 1:length(Amps_seq)
                % trials for THIS set, THIS PTD, THIS amplitude
                tr_ids = find(combClass_seq==ss & ...
                              ptdIdx_seq==p       & ...
                              ampIdx_seq==ai);
                if isempty(tr_ids), continue; end
                % Build PSTH across these trials
                counts = zeros(1, length(edges)-1);
                for k = 1:numel(tr_ids)
                    t0 = trig_seq(tr_ids(k))/FS*1000;   % ms
                    tt = S_ch(:,1) - t0;                % align to trigger
                    tt = tt(tt >= post_win_ms(1) & tt <= post_win_ms(2));
                    counts = counts + histcounts(tt, edges);
                end
                rate_curve = counts / (numel(tr_ids) * bin_s);  % sp/s
                rate_curve = conv(rate_curve, g, 'same');       % smooth
                PeakFR_seq(ci, ai, ss, p) = max(rate_curve);
            end
        end
    end
    
end

%% ===================== PLOT RESULTS ======================

sim_col = [0 0.3 0.8];
seq_color_map = lines(nSets_seq);

for ci = 1:length(resp_channels)

    ch = resp_channels(ci);

    figure('Color','w','Position',[200 200 600 450]);
    hold on;

    % ---- SIMULTANEOUS ----
    stimCh_sim = uniqueComb_sim(1,uniqueComb_sim(1,:)>0);
    sim_label = sprintf('Sim Ch %s', mat2str(stimCh_sim));

    plot(Amps_sim, PeakFR_sim(ci,:), '-o', ...
        'Color', sim_col, 'LineWidth', 2, 'MarkerFaceColor', sim_col, ...
        'DisplayName', sim_label);

    % ---- SEQUENTIAL ----
    for ss = 1:nSets_seq
        for p = 1:numel(PTDs_ms)

            stimCh = uniqueComb_seq(ss, uniqueComb_seq(ss,:)>0);
            lbl = sprintf('Seq Ch %s | PTD %g ms', mat2str(stimCh), PTDs_ms(p));

            plot(Amps_seq, PeakFR_seq(ci,:,ss,p), '-o', ...
                'Color', seq_color_map(ss,:), 'LineWidth', 1.6, ...
                'DisplayName', lbl);
        end
    end

    is_resp_sim = resp_sim(ch);
    is_resp_seq = resp_seq(ch);   
    if is_resp_sim && is_resp_seq
       resp_label = 'RESP: Sim + Seq';
    elseif is_resp_sim
        resp_label = 'RESP: Sim only';
    elseif is_resp_seq
        resp_label = 'RESP: Seq only';
    else
        resp_label = 'RESP: (unexpected none)'; % Should never happen
    end    
    stimCh_sim = uniqueComb_sim(1, uniqueComb_sim(1,:)>0);
    title(sprintf('Channel %d — Peak FR vs Amplitude\n%s', ch, resp_label), 'FontWeight','bold');
    xlabel('Amplitude (µA)');
    ylabel('Peak firing rate (sp/s)');
    legend('Location','eastoutside','Box','off');
    box off;
end
%% ==================== HELPER FUNCTION =========================
function sp = load_ssd_spikes(folder)
    cd(folder);
    f = dir('*sp_xia_SSD.mat');
    if isempty(f), error('No sp_xia_SSD.mat in %s', folder); end
    S = load(f(1).name);
    if isfield(S,'sp_corr'), sp=S.sp_corr;
    elseif isfield(S,'sp_SSD'), sp=S.sp_SSD;
    elseif isfield(S,'sp_in'), sp=S.sp_in;
    else, error('SSD file missing variables.'); end
end