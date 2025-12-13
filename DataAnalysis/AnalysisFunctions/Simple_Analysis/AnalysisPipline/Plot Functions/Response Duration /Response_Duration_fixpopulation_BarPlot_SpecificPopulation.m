%% ============================================================
%   Response Duration Analysis (FINAL)
%   - Metric: Response Span (End - Start)
%   - Method: Dual Threshold (Detect > 3.0 SD, Measure > 1.5 SD)
%   - Stats: N-Way ANOVA with Post-Hoc Comparison
%   - Output: Saves to 'Result_ResponseDuration_... .mat'
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';
Electrode_Type = 1;

% Analysis window
search_win_ms = [2 20];  
FS = 30000;  

% Thresholds (Dual)
baseline_win_ms = [-90 -10];
baseDur_s       = (baseline_win_ms(2) - baseline_win_ms(1)) / 1000;
z_detect   = 3.0; % Peak must reach this to be "Significant"
z_boundary = 1.5; % Onset/Offset measured at this lower level

% Smoothing
bin_ms     = 1;
sigma_ms   = 3;        
sigma_bins = sigma_ms / bin_ms;

% Plotting Aesthetics
target_amp = 10; 
jitter_width = 0.2;     
scatter_alpha = 0.4;    
dot_size = 25;          

% Kernel
edges = -100:bin_ms:100; 
bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g = gausswin(kernel_size); g = g / sum(g); 

%% =================== 1. LOAD DATA & QC INFO ====================
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% Extract Params (Sim)
Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim;
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

% Extract Params (Seq)
Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP; 
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

% Parse Sets
stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq); for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v  = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

stim_chans_sim = Rsim.set(1).stimChannels; 
d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

%% =================== 2. DEFINE POPULATIONS =================
% A. FIXED POPULATION (Union of all responses across all amps)
is_fixed_sim = false(nCh_Total, 1);
for ai=1:length(Amps_sim), for ch=1:nCh_Total, try if Rsim.set(1).amp(ai).ptd(1).channel(ch).is_responsive, is_fixed_sim(ch)=true; end; catch; end; end; end

is_fixed_seq = false(nCh_Total, nSets_seq);
for ss=1:nSets_seq, for ai=1:length(Amps_seq), for ch=1:nCh_Total, try if Rseq.set(ss).amp(ai).ptd(1).channel(ch).is_responsive, is_fixed_seq(ch,ss)=true; end; catch; end; end; end; end

is_fixed_Resp_Chn = find(is_fixed_seq(:,1) | is_fixed_sim | is_fixed_seq(:,2)); % Union Indices

%% =================== 3. COMPUTE DURATION =================
% Arrays for Data
Dur_sim = nan(nCh_Total, length(Amps_sim));
Dur_seq = nan(nCh_Total, length(Amps_seq), nSets_seq, numel(PTDs_ms));

% Init Timing Arrays
Onset_sim = nan(nCh_Total, length(Amps_sim)); Offset_sim = nan(nCh_Total, length(Amps_sim));
Onset_seq = nan(nCh_Total, length(Amps_seq), nSets_seq, numel(PTDs_ms)); Offset_seq = nan(nCh_Total, length(Amps_seq), nSets_seq, numel(PTDs_ms));

% --- SIMULTANEOUS ---
for ai = 1:length(Amps_sim)
    tr_ids = find(ampIdx_sim == ai);
    if isempty(tr_ids), continue; end
    
    for ch = 1:nCh_Total
        % QC Filter
        bad_trs = []; if ~isempty(QC_Sim.BadTrials) && ch<=numel(QC_Sim.BadTrials), bad_trs = QC_Sim.BadTrials{ch}; end
        tr_ids_clean = setdiff(tr_ids, bad_trs);
        
        recCh = d(ch); S_ch = sp_sim{recCh};
        
        [dur, t_on, t_off] = get_duration_dual(tr_ids_clean, trig_sim, S_ch, ...
            baseline_win_ms, baseDur_s, search_win_ms, bin_ms, g, z_detect, z_boundary, FS);
        
        % Fixed Population Metric
        if ismember(ch, is_fixed_Resp_Chn)
            Dur_sim(ch, ai) = dur; 
        end
        
        Onset_sim(ch,ai) = t_on; Offset_sim(ch,ai) = t_off;
    end
end

% --- SEQUENTIAL ---
for ss = 1:nSets_seq
    for p = 1:numel(PTDs_ms)
        for ai = 1:length(Amps_seq)
            tr_ids = find(combClass_seq==ss & ptdIdx_seq==p & ampIdx_seq==ai);
            if isempty(tr_ids), continue; end
            
            for ch = 1:nCh_Total
                bad_trs = []; if ~isempty(QC_Seq.BadTrials) && ch<=numel(QC_Seq.BadTrials), bad_trs = QC_Seq.BadTrials{ch}; end
                tr_ids_clean = setdiff(tr_ids, bad_trs);
                
                recCh = d(ch); S_ch = sp_seq{recCh};
                
                [dur, t_on, t_off] = get_duration_dual(tr_ids_clean, trig_seq, S_ch, ...
                    baseline_win_ms, baseDur_s, search_win_ms, bin_ms, g, z_detect, z_boundary, FS);
                
                % Fixed Population Metric
                if ismember(ch, is_fixed_Resp_Chn)
                    Dur_seq(ch, ai, ss, p) = dur; 
                end
                
                Onset_seq(ch,ai,ss,p) = t_on; Offset_seq(ch,ai,ss,p) = t_off;
            end
        end
    end
end

%% ===================== 4. PLOT FIGURES ======================
figure('Color','w', 'Position',[100 100 1100 550]); 
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- PANEL 1: TREND ---
ax1 = nexttile([1 2]); hold(ax1, 'on');
% Sim
AvgSim = mean(Dur_sim, 1, 'omitnan');
SEMSim = std(Dur_sim, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Dur_sim), 1));
if any(~isnan(AvgSim))
    for i = 1:length(Amps_sim)
        vals = Dur_sim(:,i); valid = ~isnan(vals);
        if any(valid)
            x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps_sim(i);
            scatter(ax1, x_jit, vals(valid), dot_size, [0 0.3 0.8], 'filled', ...
                'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
        end
    end
    plot_shaded_error(Amps_sim, AvgSim, SEMSim, [0 0.3 0.8]);
    plot(ax1, Amps_sim, AvgSim, '-o', 'Color', [0 0.3 0.8], 'LineWidth', 2, 'DisplayName', 'Simultaneous');
end
% Seq
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 
for ss = 1:nSets_seq
    data_set = squeeze(Dur_seq(:, :, ss, 1));
    AvgSeq = mean(data_set, 1, 'omitnan');
    SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
    col = set_colors(mod(ss-1,3)+1, :);
    if any(~isnan(AvgSeq))
        for i = 1:length(Amps_seq)
            vals = data_set(:,i); valid = ~isnan(vals);
            if any(valid)
                x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps_seq(i);
                scatter(ax1, x_jit, vals(valid), dot_size, col, 'filled', ...
                    'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
            end
        end
        plot_shaded_error(Amps_seq, AvgSeq, SEMSeq, col);
        plot(ax1, Amps_seq, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, 'DisplayName', sprintf('Seq Set %d', ss));
    end
end
xlabel(ax1,'Amplitude (µA)', 'FontWeight','bold'); ylabel(ax1,'Duration (ms)', 'FontWeight','bold'); 
title(ax1, 'Response Duration', 'FontWeight','bold'); legend(ax1,'Location','best');

% --- PANEL 2: BOX PLOT ---
ax2 = nexttile; hold(ax2, 'on');
idx_sim = find(Amps_sim == target_amp, 1); idx_seq = find(Amps_seq == target_amp, 1);
bp_data = []; bp_groups = {};
if ~isempty(idx_sim)
    d = Dur_sim(:, idx_sim); d = d(~isnan(d)); 
    bp_data = [bp_data; d]; bp_groups = [bp_groups; repmat({'Sim'}, length(d), 1)];
end
for ss = 1:nSets_seq
    d = squeeze(Dur_seq(:, idx_seq, ss, 1)); d = d(~isnan(d));
    bp_data = [bp_data; d]; bp_groups = [bp_groups; repmat({sprintf('Seq%d',ss)}, length(d), 1)];
end
if ~isempty(bp_data)
    boxplot(ax2, bp_data, bp_groups, 'Width', 0.5, 'Symbol', 'k.');
end
ylabel(ax2, 'Duration (ms)'); title(ax2, sprintf('Distribution at %.0f µA', target_amp), 'FontWeight', 'bold');

%% ============================================================
%   STATISTICS: N-Way ANOVA (Duration vs StimType & Amplitude)
% ============================================================
y_dur = []; g_stim = []; g_amp = [];

% 1. Collect Simultaneous Data
for ai = 1:length(Amps_sim)
    d = Dur_sim(:, ai); d = d(~isnan(d)); % Remove NaNs for ANOVA
    if isempty(d), continue; end
    y_dur = [y_dur; d];
    g_stim = [g_stim; repmat({'Simultaneous'}, length(d), 1)];
    g_amp  = [g_amp;  repmat(Amps_sim(ai), length(d), 1)];
end

% 2. Collect Sequential Data (Loop all sets)
for ss = 1:nSets_seq
    group_label = sprintf('Seq Set %d', ss);
    for ai = 1:length(Amps_seq)
        d = squeeze(Dur_seq(:, ai, ss, 1)); d = d(~isnan(d));
        if isempty(d), continue; end
        y_dur = [y_dur; d];
        g_stim = [g_stim; repmat({group_label}, length(d), 1)];
        g_amp  = [g_amp;  repmat(Amps_seq(ai), length(d), 1)];
    end
end

% 3. Run ANOVA
fprintf('\n=== ANOVA RESULTS ===\n');
[p, tbl, stats] = anovan(y_dur, {g_stim, g_amp}, 'model', 'interaction', ...
    'varnames', {'StimType', 'Amplitude'}, 'display', 'off');
fprintf('P-value (StimType):  %.5f\n', p(1));
fprintf('P-value (Amplitude): %.5f\n', p(2));
fprintf('P-value (Interaction): %.5f\n', p(3));

% 4. Post-Hoc Comparison (Visual)
figure('Color','w', 'Name', 'Post-Hoc Comparison');
multcompare(stats, 'Dimension', 1);
title('Pairwise Comparison: Duration by Stim Type');


%% ==================== HELPER FUNCTIONS ====================
function [duration_ms, onset_ms, offset_ms] = get_duration_dual(tr_ids, trig, sp_data, ...
    base_win, base_dur, search_win, bin_ms, kern, z_detect, z_boundary, FS)
    
    nTr = numel(tr_ids);
    all_spikes_base = []; all_spikes_post = [];
    for k = 1:nTr
        tr = tr_ids(k); t0 = trig(tr)/FS*1000; tt = sp_data(:,1) - t0;
        all_spikes_base = [all_spikes_base; tt(tt >= base_win(1) & tt < base_win(2))];
        all_spikes_post = [all_spikes_post; tt(tt >= -20 & tt <= 80)];
    end
    if nTr == 0, duration_ms=0; onset_ms=NaN; offset_ms=NaN; return; end
    
    bin_s = bin_ms/1000;
    base_rate_avg = numel(all_spikes_base) / (nTr * base_dur);
    edges_psth = -20 : bin_ms : 80; t_centers = edges_psth(1:end-1) + bin_ms/2;
    h = histcounts(all_spikes_post, edges_psth);
    rate_raw = h / (nTr * bin_s);
    rate_smooth = conv(rate_raw, kern, 'same');
    rate_smooth(t_centers < 0) = base_rate_avg; 
    
    baseline_part = rate_smooth(t_centers >= base_win(1) & t_centers < base_win(2));
    if isempty(baseline_part), sigma_base = 5; else, sigma_base = std(baseline_part); end
    sigma_base = max(sigma_base, 2.0); 
    
    thresh_high = base_rate_avg + z_detect * sigma_base;   
    thresh_low  = base_rate_avg + z_boundary * sigma_base; 
    
    valid_idx = find(t_centers >= search_win(1) & t_centers <= search_win(2));
    if isempty(valid_idx) || max(rate_smooth(valid_idx)) < thresh_high
        duration_ms=0; onset_ms=NaN; offset_ms=NaN; return;
    end
    
    sig_bins = valid_idx(rate_smooth(valid_idx) > thresh_low);
    if isempty(sig_bins)
        duration_ms=0; onset_ms=NaN; offset_ms=NaN;
    else
        t_start = t_centers(sig_bins(1)); t_end = t_centers(sig_bins(end));
        duration_ms = t_end - t_start; onset_ms = t_start; offset_ms = t_end;
        if duration_ms < 2, duration_ms=0; onset_ms=NaN; offset_ms=NaN; end
    end
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
    f_bt = dir('*.BadTrials.mat');   if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end
function plot_shaded_error(x, y, se, col)
    if numel(x)<2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid=~isnan(y); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end