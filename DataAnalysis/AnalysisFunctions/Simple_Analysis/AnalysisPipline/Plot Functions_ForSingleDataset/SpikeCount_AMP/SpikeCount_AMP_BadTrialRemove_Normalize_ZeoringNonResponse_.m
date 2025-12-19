%% ============================================================
%   Response Magnitude: GLOBAL MAX AT REF AMP NORMALIZATION
%   - Metric: Mean Spike Count (2-20ms)
%   - Filter 1: If Channel NOT responsive -> Count = 0
%   - Filter 2: If Channel BAD -> Count = NaN 
%   - Filter 3: Exclude BAD TRIALS
%   - Normalization: Value / MAX(Sim@5uA, Seq@5uA)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim5';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq5_5ms';
Electrode_Type = 1;

% 1. Analysis Window
post_win_ms = [2 20]; 

% 2. NORMALIZATION SETTINGS
Ref_Amp = 10;            % Look at this amplitude (uA)
min_ref_response = 0.1; % Floor to avoid dividing by noise

% 3. Plotting
FS = 30000;         
jitter_width = 0.2; dot_size = 20;          

%% =================== 1. LOAD DATA ====================
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

stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq); for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

%% ================= 2. IDENTIFY UNION POPULATION =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_sim = false(nCh_Total, 1); for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1); for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end
resp_channels = find(resp_sim | resp_seq); 

%% =================== 3. COMPUTE RAW COUNTS & NORMALIZE =================
% Init Output Arrays
Norm_Sim = nan(length(resp_channels), length(Amps_sim));
Norm_Seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));
Raw_Sim_All = nan(length(resp_channels), length(Amps_sim));
Raw_Seq_All = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch_idx = resp_channels(ci); recCh = d(ch_idx);    
    
    % --- A. Calculate Raw Sim ---
    raw_sim_vec = nan(1, length(Amps_sim));
    S_ch = sp_sim{recCh};
    bad_trs_sim = []; if ~isempty(QC_Sim.BadTrials) && ch_idx <= length(QC_Sim.BadTrials), bad_trs_sim = QC_Sim.BadTrials{ch_idx}; end
    is_bad_ch_sim = false; if ~isempty(QC_Sim.BadCh) && ~isempty(QC_Sim.BadCh{1}), if ismember(ch_idx, QC_Sim.BadCh{1}), is_bad_ch_sim = true; end; end
    
    for ai = 1:length(Amps_sim)
        if is_bad_ch_sim, raw_sim_vec(ai) = NaN; continue; end
        is_resp = 0; try is_resp = Rsim.set(1).amp(ai).ptd(1).channel(ch_idx).is_responsive; catch, end
        if ~is_resp, raw_sim_vec(ai) = 0; continue; end
        tr_ids = setdiff(find(ampIdx_sim == ai), bad_trs_sim); 
        if ~isempty(tr_ids), raw_sim_vec(ai) = get_spike_count(tr_ids, trig_sim, S_ch, post_win_ms, FS); end
    end
    Raw_Sim_All(ci, :) = raw_sim_vec;
    
    % --- B. Calculate Raw Seq ---
    raw_seq_mat = nan(length(Amps_seq), nSets_seq, numel(PTDs_ms));
    S_ch = sp_seq{recCh};
    bad_trs_seq = []; if ~isempty(QC_Seq.BadTrials) && ch_idx <= length(QC_Seq.BadTrials), bad_trs_seq = QC_Seq.BadTrials{ch_idx}; end
    
    for p = 1:numel(PTDs_ms)
        for ss = 1:nSets_seq
            is_bad_ch_seq = false; if ~isempty(QC_Seq.BadCh) && ss <= length(QC_Seq.BadCh), if ismember(ch_idx, QC_Seq.BadCh{ss}), is_bad_ch_seq = true; end; end
            for ai = 1:length(Amps_seq)
                if is_bad_ch_seq, raw_seq_mat(ai,ss,p) = NaN; continue; end
                is_resp = 0; try is_resp = Rseq.set(ss).amp(ai).ptd(p).channel(ch_idx).is_responsive; catch, end
                if ~is_resp, raw_seq_mat(ai,ss,p) = 0; continue; end
                tr_ids = setdiff(find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai), bad_trs_seq);
                if ~isempty(tr_ids), raw_seq_mat(ai,ss,p) = get_spike_count(tr_ids, trig_seq, S_ch, post_win_ms, FS); end
            end
        end
    end
    Raw_Seq_All(ci, :, :, :) = raw_seq_mat;
    
    % --- C. NORMALIZE TO GLOBAL MAX AT REFERENCE AMPLITUDE (5uA) ---
    norm_denom = NaN;
    
    % 1. Get Sim value at Ref Amp (5uA)
    ref_idx_sim = find(abs(Amps_sim - Ref_Amp) < 0.1);
    val_sim_ref = [];
    if ~isempty(ref_idx_sim), val_sim_ref = raw_sim_vec(ref_idx_sim); end
    
    % 2. Get Seq values at Ref Amp (5uA) across all sets
    ref_idx_seq = find(abs(Amps_seq - Ref_Amp) < 0.1);
    val_seq_ref = [];
    if ~isempty(ref_idx_seq)
        slice = raw_seq_mat(ref_idx_seq, :, :); % [Amp x Set x PTD] -> Extract Set dim
        val_seq_ref = slice(:); % Flatten all Seq sets at 5uA
    end
    
    % 3. Find GLOBAL MAX at 5uA (Highest of Sim or Seq)
    all_ref_vals = [val_sim_ref; val_seq_ref];
    if ~isempty(all_ref_vals)
        norm_denom = max(all_ref_vals, [], 'omitnan');
    end
    
    % 4. Apply Normalization
    if ~isnan(norm_denom) && norm_denom > min_ref_response
        Norm_Sim(ci, :) = raw_sim_vec / norm_denom;
        Norm_Seq(ci, :, :, :) = raw_seq_mat / norm_denom;
    else
        norm_denom = min_ref_response;
        Norm_Sim(ci, :) = raw_sim_vec / norm_denom;
        Norm_Seq(ci, :, :, :) = raw_seq_mat / norm_denom;

        % Norm_Sim(ci, :) = NaN;
        % Norm_Seq(ci, :, :, :) = NaN;
    end
end 

%% ===================== 4. PLOT RESULTS ======================
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
        if any(~isnan(AvgSeq))
            for i = 1:length(Amps_seq)
                valid = ~isnan(data_set(:,i));
                if any(valid), x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps_seq(i); scatter(x_jit, data_set(valid,i), dot_size, col, 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off'); end
            end
            plot_shaded_error(Amps_seq, AvgSeq, SEMSeq, col);
            plot(Amps_seq, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', sprintf('Seq Set %d', ss));
        end
    end
end
yline(1.0, '--k', sprintf('Ref Max @ %.0fuA', Ref_Amp), 'HandleVisibility','off');
xlabel('Amplitude (uA)', 'FontWeight','bold'); 
ylabel(sprintf('Normalized (%.0fuA)', Ref_Amp), 'FontWeight','bold');
title(sprintf('Response Magnitude (Global Max Norm @ %.0fuA)', Ref_Amp), 'FontWeight','bold');
legend('Location','best','Box','off'); box off; 
ylim([0 3.0]);

%% ================= SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_SpikeNormGlobalRef_' num2str(Ref_Amp) 'uA_Zeroed_5ms_' exp_id '.mat']);
ResultNorm = struct();
ResultNorm.Raw_Sim = Raw_Sim_All; ResultNorm.Raw_Seq = Raw_Seq_All; 
ResultNorm.Norm_Sim = Norm_Sim; ResultNorm.Norm_Seq = Norm_Seq; ResultNorm.Amps = Amps_seq;
save(out_filename, 'ResultNorm');
fprintf('\n>>> Results Saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function count_val = get_spike_count(tr_ids, trig, sp_data, count_win, FS)
    nTr = numel(tr_ids); total_spikes = 0;
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