%% ============================================================
%   Response Magnitude: Normalized to Simultaneous Max (with Noise Floor)
%   - Metric: Mean Spike Count (2-20ms)
%   - RULE: If Channel is NOT responsive -> Count = 0 (Ignore noise)
%   - Normalization: Seq / max(Sim_Max, Noise_Floor)
%   - Output: Prints tables for both Raw and Normalized data.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim4_251125_152849';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq4_5ms_251125_154235';
Electrode_Type = 1;

% 1. Analysis Window (Total Spike Count captures BOTH peaks)
post_win_ms = [2 20]; 

% 2. NOISE FLOOR (The Fix)
noise_floor = 0.1; 

% 3. PSTH Settings
FS = 30000; bin_ms = 1; sigma_bins = 3;        
psth_win_ms = [-50 100]; 
edges_psth = psth_win_ms(1):bin_ms:psth_win_ms(2);
bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g_sym = gausswin(kernel_size); g_sym = g_sym / sum(g_sym);

% Plotting
jitter_width = 0.2; dot_size = 20;          

%% =================== 1. LOAD DATA ====================
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% Extract Stim Params
Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim;
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP;
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq); for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);
stim_chans_sim = Rsim.set(1).stimChannels; 

%% ================= 2. IDENTIFY UNION POPULATION =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_sim = false(nCh_Total, 1); for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1); for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end
resp_channels = find(resp_sim | resp_seq); 

%% =================== 3. COMPUTE RAW COUNTS & NORMALIZE =================
% Initialize arrays for Normalized Data
Norm_Sim = nan(length(resp_channels), length(Amps_sim));
Norm_Seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

% Initialize arrays for Raw Data (Unnormalized)
Raw_Sim_All = nan(length(resp_channels), length(Amps_sim));
Raw_Seq_All = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch_idx = resp_channels(ci); recCh = d(ch_idx);    
    
    % --- A. Calculate Raw Sim ---
    raw_sim_vec = nan(1, length(Amps_sim));
    S_ch = sp_sim{recCh};
    bad_trs_sim = []; if ~isempty(QC_Sim.BadTrials) && ch_idx <= length(QC_Sim.BadTrials), bad_trs_sim = QC_Sim.BadTrials{ch_idx}; end
    
    for ai = 1:length(Amps_sim)
        % --- RULE: Check Responsiveness ---
        is_resp = 0; 
        try is_resp = Rsim.set(1).amp(ai).ptd(1).channel(ch_idx).is_responsive; catch, end
        
        if ~is_resp
            raw_sim_vec(ai) = 0; % Force 0 if not responsive
            continue;
        end
        % ----------------------------------
        
        tr_ids = setdiff(find(ampIdx_sim == ai), bad_trs_sim); 
        if ~isempty(tr_ids), raw_sim_vec(ai) = get_spike_count(tr_ids, trig_sim, S_ch, post_win_ms, FS); end
    end
    Raw_Sim_All(ci, :) = raw_sim_vec; % Store Raw
    
    % --- B. Calculate Raw Seq ---
    raw_seq_mat = nan(length(Amps_seq), nSets_seq, numel(PTDs_ms));
    S_ch = sp_seq{recCh};
    bad_trs_seq = []; if ~isempty(QC_Seq.BadTrials) && ch_idx <= length(QC_Seq.BadTrials), bad_trs_seq = QC_Seq.BadTrials{ch_idx}; end
    
    for p = 1:numel(PTDs_ms)
        for ss = 1:nSets_seq
            for ai = 1:length(Amps_seq)
                % --- RULE: Check Responsiveness ---
                is_resp = 0;
                try is_resp = Rseq.set(ss).amp(ai).ptd(p).channel(ch_idx).is_responsive; catch, end
                
                if ~is_resp
                    raw_seq_mat(ai,ss,p) = 0; % Force 0 if not responsive
                    continue;
                end
                % ----------------------------------

                tr_ids = setdiff(find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai), bad_trs_seq);
                if ~isempty(tr_ids), raw_seq_mat(ai,ss,p) = get_spike_count(tr_ids, trig_seq, S_ch, post_win_ms, FS); end
            end
        end
    end
    Raw_Seq_All(ci, :, :, :) = raw_seq_mat; % Store Raw
    
    % --- C. NORMALIZE TO SIMULTANEOUS MAX (WITH NOISE FLOOR) ---
    sim_max = max(raw_sim_vec, [], 'omitnan');
    
    % Apply Noise Floor Logic: Use max(Actual Sim Max, Noise Floor)
    denominator = max(sim_max, noise_floor);
    
    if ~isnan(denominator) && denominator > 0
        Norm_Sim(ci, :) = raw_sim_vec / denominator;
        Norm_Seq(ci, :, :, :) = raw_seq_mat / denominator;
    end
end 

%% ===================== 4. PLOT NORMALIZED RESULTS ======================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% Plot Sim
AvgSim = mean(Norm_Sim, 1, 'omitnan');
SEMSim = std(Norm_Sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Norm_Sim), 1));
sim_col = [0 0.3 0.8];
if any(~isnan(AvgSim))
    for i = 1:length(Amps_sim)
        valid = ~isnan(Norm_Sim(:,i));
        if any(valid), x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps_sim(i); scatter(x_jit, Norm_Sim(valid,i), dot_size, sim_col, 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off'); end
    end
    plot_shaded_error(Amps_sim, AvgSim, SEMSim, sim_col);
    plot(Amps_sim, AvgSim, '-o', 'Color', sim_col, 'LineWidth', 2, 'MarkerFaceColor', sim_col, 'DisplayName', 'Simultaneous');
end

% Plot Seq
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 
for p = 1:numel(PTDs_ms)
    for ss = 1:nSets_seq
        data_set = squeeze(Norm_Seq(:, :, ss, p));
        AvgSeq = mean(data_set, 1, 'omitnan');
        SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
        col = set_colors(mod(ss-1,3)+1, :);
        
        stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
        lbl = sprintf('Seq Set %d', ss);
        
        if any(~isnan(AvgSeq))
            for i = 1:length(Amps_seq)
                valid = ~isnan(data_set(:,i));
                if any(valid), x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps_seq(i); scatter(x_jit, data_set(valid,i), dot_size, col, 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off'); end
            end
            plot_shaded_error(Amps_seq, AvgSeq, SEMSeq, col);
            plot(Amps_seq, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', lbl);
        end
    end
end
yline(1.0, '--k', 'Sim Max', 'HandleVisibility','off');
xlabel('Amplitude (uA)', 'FontWeight','bold'); ylabel('Normalized Spike/Trial', 'FontWeight','bold');
title(sprintf('Spike per Trial'), 'FontWeight','bold');
legend('Location','best','Box','off'); box off; ylim([0 3.0]);


%% ================= PRINT SUMMARY TABLES =================
fprintf('\n\n========================================================================\n');
fprintf('                 TABLE 1: UNNORMALIZED (RAW SPIKE COUNTS)               \n');
fprintf('========================================================================\n');
fprintf('Condition      \t');
for ai = 1:length(Amps_sim), fprintf('Amp_%.0fuA_Mean\tSEM\t', Amps_sim(ai)); end; fprintf('\n');

% 1. Simultaneous (Raw)
fprintf('Simultaneous   \t');
for ai = 1:length(Amps_sim)
    d = Raw_Sim_All(:, ai); mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan')/sqrt(sum(~isnan(d)));
    fprintf('%.4f\t%.4f\t', mu, sem);
end
fprintf('\n');

% 2. Sequential Sets (Raw)
for ss = 1:nSets_seq
    stimCh = uniqueComb_seq(ss,:); stimCh=stimCh(stimCh>0);
    fprintf('Seq_Set%d(Ch%s)\t', ss, num2str(stimCh));
    for ai = 1:length(Amps_seq)
        d = squeeze(Raw_Seq_All(:, ai, ss, 1)); mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan')/sqrt(sum(~isnan(d)));
        fprintf('%.4f\t%.4f\t', mu, sem);
    end
    fprintf('\n');
end

fprintf('\n\n========================================================================\n');
fprintf('                 TABLE 2: NORMALIZED (TO SIM MAX w/ FLOOR)              \n');
fprintf('========================================================================\n');
fprintf('Condition      \t');
for ai = 1:length(Amps_sim), fprintf('Amp_%.0fuA_Mean\tSEM\t', Amps_sim(ai)); end; fprintf('\n');

% 1. Simultaneous (Norm)
fprintf('Simultaneous   \t');
for ai = 1:length(Amps_sim)
    d = Norm_Sim(:, ai); mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan')/sqrt(sum(~isnan(d)));
    fprintf('%.4f\t(%.4f)\t', mu, sem);
end
fprintf('\n');

% 2. Sequential Sets (Norm)
for ss = 1:nSets_seq
    stimCh = uniqueComb_seq(ss,:); stimCh=stimCh(stimCh>0);
    fprintf('Seq_Set%d(Ch%s)\t', ss, num2str(stimCh));
    for ai = 1:length(Amps_seq)
        d = squeeze(Norm_Seq(:, ai, ss, 1)); mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan')/sqrt(sum(~isnan(d)));
        fprintf('%.4f\t(%.4f)\t', mu, sem);
    end
    fprintf('\n');
end
fprintf('========================================================================\n');

%% ================= SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/DX012';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Set4_SpikeNormSimMax_Zeroed_5ms_' exp_id '.mat']);

ResultNorm = struct();
ResultNorm.Raw_Sim = Raw_Sim_All; % Save Raw
ResultNorm.Raw_Seq = Raw_Seq_All; % Save Raw
ResultNorm.Norm_Sim = Norm_Sim;
ResultNorm.Norm_Seq = Norm_Seq;
ResultNorm.Amps = Amps_seq;
save(out_filename, 'ResultNorm');
fprintf('\n>>> Results Saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function count_val = get_spike_count(tr_ids, trig, sp_data, count_win, FS)
    nTr = numel(tr_ids);
    total_spikes = 0;
    for k = 1:nTr
        tr = tr_ids(k); t0 = trig(tr)/FS*1000; tt = sp_data(:,1) - t0;
        mask = tt >= count_win(1) & tt <= count_win(2);
        total_spikes = total_spikes + sum(mask);
    end
    if nTr > 0, count_val = total_spikes / nTr; else, count_val = NaN; end
end
function plot_shaded_error(x, y, se, col)
    if numel(x) < 2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid = ~isnan(y); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
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
    f_bt = dir('*.BadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end