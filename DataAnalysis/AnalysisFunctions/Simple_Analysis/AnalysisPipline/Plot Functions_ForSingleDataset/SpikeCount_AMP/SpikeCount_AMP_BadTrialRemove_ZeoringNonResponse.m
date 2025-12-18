%% ============================================================
%   Spike Count Analysis (FIXED POPULATION / UNION)
%   - Metric: Mean Spike Count per Trial (Raw)
%   - Logic: Counts spikes in window [2 20]ms
%   - REFINEMENT 1: Non-responsive channels -> Count = 0
%   - REFINEMENT 2: Bad Channels -> Count = NaN (Excluded from mean)
%   - REFINEMENT 3: Bad Trials -> Removed from trial average
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim3';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq3';
Electrode_Type = 1;

% 1. Analysis Window (Spike Counting)
post_win_ms = [2 20]; 

% 2. Plotting & PSTH
psth_win_ms = [-50 100]; 
FS = 30000;            
bin_ms     = 1;
sigma_bins = 3;        
jitter_width = 0.2; scatter_alpha = 0.4; dot_size = 20;          

% Kernels
edges_peak = post_win_ms(1):bin_ms:post_win_ms(2);
bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g_sym = gausswin(kernel_size); g_sym = g_sym / sum(g_sym);
edges_psth = psth_win_ms(1):bin_ms:psth_win_ms(2);

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

% Parse Seq Sets
stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

%% ================= 2. IDENTIFY UNION POPULATION =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_sim = false(nCh_Total, 1);
for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1);
for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive') && this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end

resp_channels = find(resp_sim | resp_seq); 
fprintf('Analyzing Fixed Population (Union): %d Channels\n', length(resp_channels));

%% =================== 3. COMPUTE SPIKE COUNTS (FIXED POPULATION) =================
% Init Output Arrays
SpikeCount_sim = nan(length(resp_channels), length(Amps_sim));
SpikeCount_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq, numel(PTDs_ms));
nBins = length(edges_psth)-1;
PSTH_Sim = nan(length(resp_channels), nBins, length(Amps_sim));
PSTH_Seq = nan(length(resp_channels), nBins, length(Amps_seq), nSets_seq, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch_idx = resp_channels(ci); recCh = d(ch_idx);    
    
    % --- SIMULTANEOUS ---
    S_ch = sp_sim{recCh};
    
    % Get Bad Trials
    bad_trs = []; 
    if ~isempty(QC_Sim.BadTrials) && ch_idx <= length(QC_Sim.BadTrials)
        bad_trs = QC_Sim.BadTrials{ch_idx}; 
    end
    
    % Get Bad Channels (Simultaneous usually Set 1 or flat list)
    is_bad_ch = false;
    if ~isempty(QC_Sim.BadCh)
        if iscell(QC_Sim.BadCh)
             % Check if bad in ANY set (conservative)
             all_bad = [QC_Sim.BadCh{:}];
             if ismember(ch_idx, all_bad), is_bad_ch = true; end
        else
             if ismember(ch_idx, QC_Sim.BadCh), is_bad_ch = true; end
        end
    end

    for ai = 1:length(Amps_sim)
        % RULE 1: If Bad Channel -> NaN (Ignore from Mean)
        if is_bad_ch
            SpikeCount_sim(ci,ai) = NaN; 
            PSTH_Sim(ci, :, ai)   = nan(1, nBins);
            continue;
        end
        
        % RULE 2: Check Responsiveness
        is_resp = 0; 
        try is_resp = Rsim.set(1).amp(ai).ptd(1).channel(ch_idx).is_responsive; catch, end
        if ~is_resp
            SpikeCount_sim(ci,ai) = 0; 
            PSTH_Sim(ci, :, ai)   = zeros(1, nBins); 
            continue; 
        end
        
        % RULE 3: Filter Bad Trials
        tr_ids = find(ampIdx_sim == ai);
        tr_ids = setdiff(tr_ids, bad_trs); 
        
        if isempty(tr_ids), continue; end
        
        [count_val, psth_curve] = get_spike_count(tr_ids, trig_sim, S_ch, ...
            post_win_ms, edges_psth, g_sym, bin_s, FS);
            
        SpikeCount_sim(ci,ai) = count_val; 
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
            
            % Check Bad Channel for THIS Set
            is_bad_ch_seq = false;
            if ~isempty(QC_Seq.BadCh) && ss <= length(QC_Seq.BadCh)
                if ismember(ch_idx, QC_Seq.BadCh{ss}), is_bad_ch_seq = true; end
            end

            for ai = 1:length(Amps_seq)
                
                % RULE 1: Bad Channel -> NaN
                if is_bad_ch_seq
                    SpikeCount_seq(ci,ai,ss,p) = NaN;
                    PSTH_Seq(ci, :, ai, ss, p) = nan(1, nBins);
                    continue;
                end
                
                % RULE 2: Check Responsiveness
                is_resp = 0;
                try is_resp = Rseq.set(ss).amp(ai).ptd(p).channel(ch_idx).is_responsive; catch, end
                if ~is_resp
                    SpikeCount_seq(ci,ai,ss,p) = 0; 
                    PSTH_Seq(ci, :, ai, ss, p) = zeros(1, nBins); 
                    continue; 
                end
                
                % RULE 3: Filter Bad Trials
                tr_ids = find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai);
                tr_ids = setdiff(tr_ids, bad_trs); 
                
                if isempty(tr_ids), continue; end
                
                [count_val, psth_curve] = get_spike_count(tr_ids, trig_seq, S_ch, ...
                    post_win_ms, edges_psth, g_sym, bin_s, FS);
                
                SpikeCount_seq(ci,ai,ss,p) = count_val; 
                PSTH_Seq(ci, :, ai, ss, p) = psth_curve;
            end
        end
    end
end 

%% ===================== 4. PLOT (Full Figure) ======================
figure('Color','w', 'Position',[100 100 800 600]); hold on;

% --- Plot Simultaneous ---
AvgSim = mean(SpikeCount_sim, 1, 'omitnan');
SEMSim = std(SpikeCount_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(SpikeCount_sim), 1));
sim_col = [0 0.3 0.8];
if any(~isnan(AvgSim))
    for i = 1:length(Amps_sim)
        valid = ~isnan(SpikeCount_sim(:,i));
        if any(valid)
            x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps_sim(i);
            scatter(x_jit, SpikeCount_sim(valid,i), dot_size, sim_col, 'filled', ...
                'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
        end
    end
    plot_shaded_error(Amps_sim, AvgSim, SEMSim, sim_col);
    plot(Amps_sim, AvgSim, '-o', 'Color', sim_col, 'LineWidth', 2, ...
        'MarkerFaceColor', sim_col, 'DisplayName', 'Simultaneous');
end

% --- Plot Sequential Sets ---
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 
for p = 1:numel(PTDs_ms)
    for ss = 1:nSets_seq
        data_set = squeeze(SpikeCount_seq(:, :, ss, p));
        AvgSeq = mean(data_set, 1, 'omitnan');
        SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
        col = set_colors(mod(ss-1,3)+1, :);
        
        stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
        lbl = sprintf('Seq Set %d (Ch:%s)', ss, num2str(stimCh));
        
        if any(~isnan(AvgSeq))
            for i = 1:length(Amps_seq)
                valid = ~isnan(data_set(:,i));
                if any(valid)
                    x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps_seq(i);
                    scatter(x_jit, data_set(valid,i), dot_size, col, 'filled', ...
                        'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
                end
            end
            plot_shaded_error(Amps_seq, AvgSeq, SEMSeq, col);
            plot(Amps_seq, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, ...
                'MarkerFaceColor', 'w', 'DisplayName', lbl);
        end
    end
end

xlabel('Amplitude (uA)', 'FontWeight','bold'); 
ylabel('Mean Spike Count (2-20 ms)', 'FontWeight','bold');
title('Response Magnitude (Total Spike Count - Cleaned)', 'FontWeight','bold');
legend('Location','best','Box','off'); box off;

%% ============================================================
%   5. STATISTICS: POOLED ANOVA (Sim vs All Seq)
% ============================================================
y_val  = []; g_stim = []; g_amp  = [];
% Sim
for ai = 1:length(Amps_sim)
    d = SpikeCount_sim(:, ai); d = d(~isnan(d));
    if isempty(d), continue; end
    y_val  = [y_val; d];
    g_stim = [g_stim; repmat({'Simultaneous'}, length(d), 1)];
    g_amp  = [g_amp;  repmat(Amps_sim(ai), length(d), 1)];
end
% Seq (Pooled)
for ss = 1:nSets_seq
    for ai = 1:length(Amps_seq)
        d = squeeze(SpikeCount_seq(:, ai, ss, 1)); d = d(~isnan(d));
        if isempty(d), continue; end
        y_val  = [y_val; d];
        g_stim = [g_stim; repmat({'Sequential'}, length(d), 1)];
        g_amp  = [g_amp;  repmat(Amps_seq(ai), length(d), 1)];
    end
end
if ~isempty(y_val)
    fprintf('\n=== ANOVA RESULTS (Pooled) ===\n');
    [p, tbl, stats] = anovan(y_val, {g_stim, g_amp}, ...
        'model', 'interaction', 'varnames', {'StimType', 'Amplitude'}, 'display', 'on');
    fprintf('StimType P: %.5f\nAmplitude P: %.5f\nInteract P: %.5f\n', p(1), p(2), p(3));
    
    figure('Color','w','Name','SpikeCount Comparison');
    multcompare(stats, 'Dimension', 1);
end

%% =============== PRINT SUMMARY STATISTICS ===============
fprintf('\n\n========================================================================\n');
fprintf('       SUMMARY STATISTICS: SPIKE COUNT         \n');
fprintf('========================================================================\n');
fprintf('Condition       \t');
for ai = 1:length(Amps_sim), fprintf('%.0fuA_Mean SEM\t', Amps_sim(ai)); end
fprintf('\n');
% 2. Simultaneous Row
fprintf('Simultaneous       \t');
for ai = 1:length(Amps_sim)
    d = SpikeCount_sim(:, ai);
    mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan') / sqrt(sum(~isnan(d)));
    if isnan(mu), fprintf('NaN\tNaN\t'); else, fprintf('%.4f %.4f\t', mu, sem); end
end
fprintf('\n');
% 3. Sequential Rows
for ss = 1:nSets_seq
    stimCh_List = uniqueComb_seq(ss, :); stimCh_List = stimCh_List(stimCh_List > 0);
    fprintf('Seq_Set%d(Ch%s)\t', ss, num2str(stimCh_List));
    for ai = 1:length(Amps_seq)
        d = squeeze(SpikeCount_seq(:, ai, ss, 1)); 
        mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan') / sqrt(sum(~isnan(d)));
        if isnan(mu), fprintf('NaN\tNaN\t'); else, fprintf('%.4f %.4f\t', mu, sem); end
    end
    fprintf('\n');
end
fprintf('========================================================================\n');

%% ============================================================
%   6. SAVE RESULTS
% ============================================================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX011/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_SpikeCount_FixedPop_Zeroed_5ms_' exp_id '.mat']);

ResultFR = struct();
ResultFR.Metadata.Created = datestr(now);
ResultFR.Metadata.Metric = 'Mean Spike Count per Trial (Cleaned)';
ResultFR.Metadata.Window = post_win_ms;
ResultFR.Metadata.Amps = Amps_sim;
ResultFR.PeakFR.Sim_Z = SpikeCount_sim; 
ResultFR.PeakFR.Seq_Z = SpikeCount_seq;
ResultFR.PSTH.Sim = PSTH_Sim;
ResultFR.PSTH.Seq = PSTH_Seq;
save(out_filename, 'ResultFR');
fprintf('\n>>> Results Saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function [count_val, psth_trace] = get_spike_count(tr_ids, trig, sp_data, ...
    count_win, psth_edges, g_sym, bin_s, FS)
    
    nTr = numel(tr_ids);
    all_psth_spikes = [];
    total_spikes_in_window = 0;
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        % 1. Count spikes for Metric
        mask_count = tt >= count_win(1) & tt <= count_win(2);
        total_spikes_in_window = total_spikes_in_window + sum(mask_count);
        
        % 2. Collect spikes for PSTH Visual
        all_psth_spikes = [all_psth_spikes; tt(tt >= psth_edges(1) & tt <= psth_edges(end))];
    end
    
    if nTr > 0
        count_val = total_spikes_in_window / nTr; % Mean spikes per trial
    else
        count_val = NaN;
    end
    
    % Compute PSTH for plotting
    h_psth = histcounts(all_psth_spikes, psth_edges);
    rate_psth = h_psth / (nTr * bin_s);
    psth_trace = conv(rate_psth, g_sym, 'same');
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