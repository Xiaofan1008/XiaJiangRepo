%% ============================================================
%   Peak Latency Analysis (Grouped Bars + Smooth Curves)
%   Fig 1: Grouped Histograms + KDE Curves + Medians
%   Fig 2: Latency vs Amplitude Trend
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type = 1;
search_win_ms = [0 20];  
FS = 30000;            

% Smoothing (SYMMETRIC Gaussian)
bin_ms = 1; sigma_ms = 3; sigma_bins = sigma_ms / bin_ms;

% PLOTTING SETTINGS
target_amp = 10; % Amplitude for Histogram
target_ptds_for_hist = [5]; 
hist_edges = 0:1:20; % Bins

% Kernel
edges = -100:bin_ms:100; bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g = gausswin(kernel_size); g = g / sum(g); 

%% =================== LOAD DATA (Collapsed) ====================
cd(folder_sim); Rsim=load(dir('*RespondingChannels.mat').name).Responding; sp_sim=load_ssd_spikes(folder_sim); if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig_sim=loadTrig(0); Ssim=load(dir('*_exp_datafile_*.mat').name); Stim_sim=Ssim.StimParams; simN_sim=Ssim.simultaneous_stim; amps_all_sim=cell2mat(Stim_sim(2:end,16)); trialAmps_sim=amps_all_sim(1:simN_sim:end); [Amps_sim,~,ampIdx_sim]=unique(trialAmps_sim); Amps_sim(Amps_sim==-1)=0;
cd(folder_seq); Rseq=load(dir('*RespondingChannels.mat').name).Responding; sp_seq=load_ssd_spikes(folder_seq); if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig_seq=loadTrig(0); Sseq=load(dir('*_exp_datafile_*.mat').name); Stim_seq=Sseq.StimParams; simN_seq=Sseq.simultaneous_stim; E_MAP_seq=Sseq.E_MAP; if isfield(Sseq,'n_Trials'), nTr_seq=Sseq.n_Trials; else, nTr_seq=(size(Stim_seq,1)-1)/simN_seq; end; amps_all_seq=cell2mat(Stim_seq(2:end,16)); trialAmps_seq=amps_all_seq(1:simN_seq:end); [Amps_seq,~,ampIdx_seq]=unique(trialAmps_seq); Amps_seq(Amps_seq==-1)=0;
PTD_all_us=cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms=PTD_all_us/1000; [PTDs_ms,~,ptdIdx_seq]=unique(PTD_all_ms); stimNames=Stim_seq(2:end,1); [~,idx_all]=ismember(stimNames,E_MAP_seq(2:end)); comb_seq=zeros(nTr_seq,simN_seq); for t=1:nTr_seq, rr=(t-1)*simN_seq+(1:simN_seq); v=idx_all(rr); v=v(v>0); comb_seq(t,1:numel(v))=v(:).'; end; [uniqueComb_seq,~,combClass_seq]=unique(comb_seq,'rows','stable'); nSets_seq=size(uniqueComb_seq,1);

d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_sim = false(nCh_Total, 1); for si=1:numel(Rsim.set), for ai=1:numel(Rsim.set(si).amp), for pi=1:numel(Rsim.set(si).amp(ai).ptd), this=Rsim.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive')&&this(ch).is_responsive, resp_sim(ch)=true; end; end; end; end; end
resp_seq = false(nCh_Total, 1); for si=1:numel(Rseq.set), for ai=1:numel(Rseq.set(si).amp), for pi=1:numel(Rseq.set(si).amp(ai).ptd), this=Rseq.set(si).amp(ai).ptd(pi).channel; for ch=1:min(length(this),nCh_Total), if isfield(this(ch),'is_responsive')&&this(ch).is_responsive, resp_seq(ch)=true; end; end; end; end; end
resp_channels = find(resp_sim | resp_seq); 

%% =================== 1. CALC LATENCY (Distributions at 10uA) =================
idx_sim_t = find(Amps_sim == target_amp, 1);
idx_seq_t = find(Amps_seq == target_amp, 1);

Lat_sim = []; Lat_seq = cell(nSets_seq, 1); 

% Sim
if ~isempty(idx_sim_t)
    for ci=1:length(resp_channels), ch=resp_channels(ci); recCh=d(ch);
        tr_ids=find(ampIdx_sim==idx_sim_t); if isempty(tr_ids), continue; end
        lat=get_peak_latency(tr_ids, trig_sim, sp_sim{recCh}, search_win_ms, bin_ms, g, FS);
        Lat_sim=[Lat_sim; lat];
    end
end
% Seq (Target PTD)
ptd_idx = find(abs(PTDs_ms - target_ptds_for_hist(1)) < 0.1, 1);
for ss=1:nSets_seq
    temp=[];
    for ci=1:length(resp_channels), ch=resp_channels(ci); recCh=d(ch);
        tr_ids=find(combClass_seq==ss & ptdIdx_seq==ptd_idx & ampIdx_seq==idx_seq_t); if isempty(tr_ids), continue; end
        lat=get_peak_latency(tr_ids, trig_seq, sp_seq{recCh}, search_win_ms, bin_ms, g, FS);
        temp=[temp; lat];
    end
    Lat_seq{ss} = temp;
end

%% ===================== FIGURE 1: GROUPED BARS + SMOOTH + MEDIAN ======================
figure('Color','w', 'Position',[100 100 900 600]); hold on;

% --- 1. Prepare Data for Grouped Bar ---
num_bins = length(hist_edges) - 1;
bar_data = [];
group_names = {};
group_colors = {};

% Sim Data
[counts_sim, ~] = histcounts(Lat_sim, hist_edges);
bar_data(:, 1) = counts_sim(:) / length(Lat_sim) * 100;
group_names{1} = 'Simultaneous';
group_colors{1} = [0 0.3 0.8]; % Blue

% Seq Data
cols_seq = [0.85 0.33 0.10; 0.60 0.20 0.60]; % Orange, Purple
col_idx = 2;
for ss = 1:nSets_seq
    data = Lat_seq{ss};
    if isempty(data), continue; end
    
    [counts_seq, ~] = histcounts(data, hist_edges);
    bar_data(:, col_idx) = counts_seq(:) / length(data) * 100;
    
    stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
    group_names{col_idx} = sprintf('Seq Set %d (Ch:%s)', ss, num2str(stimCh));
    group_colors{col_idx} = cols_seq(mod(ss-1,2)+1, :);
    col_idx = col_idx + 1;
end

% --- 2. Plot Grouped Bars ---
b = bar(hist_edges(1:end-1)+diff(hist_edges)/2, bar_data, 'grouped');

for i = 1:length(b)
    b(i).FaceColor = group_colors{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 0.3; 
    b(i).DisplayName = [group_names{i} ' (Bars)'];
end

% --- 3. Overlay Smooth Curves ---
% Sim Curve
[f_sim, xi_sim] = ksdensity(Lat_sim, 'Bandwidth', 1.5);
scale = 100 * (hist_edges(2)-hist_edges(1));
plot(xi_sim, f_sim*scale, 'Color', group_colors{1}, 'LineWidth', 2.5, 'DisplayName', 'Sim (Smooth)');

% Seq Curves
for i = 2:length(group_names) % Skip Sim (index 1)
    % Find correct data for this group (matching the loop order above)
    % Note: The order in 'group_names' matches 'Lat_seq' iteration
    ss_idx = i - 1; 
    data = Lat_seq{ss_idx}; 
    
    [f_seq, xi_seq] = ksdensity(data, 'Bandwidth', 1.5);
    plot(xi_seq, f_seq*scale, 'Color', group_colors{i}, 'LineWidth', 2.5, 'DisplayName', 'Seq (Smooth)');
end

% --- 4. Plot Median Lines (ON TOP) ---
% Sim Median
med_sim = median(Lat_sim, 'omitnan');
xline(med_sim, '--', 'Color', group_colors{1}, 'LineWidth', 2, 'DisplayName', 'Sim Median');
% Add text label at the top
text(med_sim, max(bar_data(:))*1.05, sprintf('Sim: %.1f', med_sim), ...
    'Color', group_colors{1}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Seq Medians
for i = 2:length(group_names)
    ss_idx = i - 1;
    data = Lat_seq{ss_idx};
    med_seq = median(data, 'omitnan');
    
    xline(med_seq, '--', 'Color', group_colors{i}, 'LineWidth', 2, 'HandleVisibility','off');
    text(med_seq, max(bar_data(:))*1.05, sprintf('%.1f', med_seq), ...
        'Color', group_colors{i}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

ylabel('% of Channels', 'FontSize', 12, 'FontWeight', 'bold'); 
xlabel('Peak Latency (ms)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Latency Distribution at %.0f µA', target_amp), 'FontSize', 14);
legend(b, group_names, 'Location','northeast'); % Legend only shows bars to keep it clean
box off; xlim([0 20]);

%% =================== 2. LATENCY vs AMPLITUDE TREND =================
LatAmp_sim = nan(length(resp_channels), length(Amps_sim));
LatAmp_seq = nan(length(resp_channels), length(Amps_seq), nSets_seq); 

% Calc Sim
for ai=1:length(Amps_sim), for ci=1:length(resp_channels), ch=resp_channels(ci); recCh=d(ch);
    tr_ids=find(ampIdx_sim==ai); if isempty(tr_ids), continue; end
    LatAmp_sim(ci,ai) = get_peak_latency(tr_ids, trig_sim, sp_sim{recCh}, search_win_ms, bin_ms, g, FS);
end; end

% Calc Seq
p_idx_trend = find(abs(PTDs_ms - target_ptds_for_hist(1)) < 0.1, 1);
for ss=1:nSets_seq, for ai=1:length(Amps_seq), for ci=1:length(resp_channels), ch=resp_channels(ci); recCh=d(ch);
    tr_ids=find(combClass_seq==ss & ptdIdx_seq==p_idx_trend & ampIdx_seq==ai); if isempty(tr_ids), continue; end
    LatAmp_seq(ci,ai,ss) = get_peak_latency(tr_ids, trig_seq, sp_seq{recCh}, search_win_ms, bin_ms, g, FS);
end; end; end

%% ===================== FIGURE 2: TREND PLOT ======================
figure('Color','w', 'Position',[900 100 600 500]); hold on;

% 1. Sim
mu = mean(LatAmp_sim, 1, 'omitnan'); err = std(LatAmp_sim,0,1,'omitnan')./sqrt(sum(~isnan(LatAmp_sim),1));
errorbar(Amps_sim, mu, err, '-o', 'Color', group_colors{1}, 'LineWidth', 2, 'MarkerFaceColor', group_colors{1}, 'CapSize', 0, 'DisplayName', 'Simultaneous');

% 2. Seq
for ss = 1:nSets_seq
    data = squeeze(LatAmp_seq(:, :, ss));
    mu = mean(data, 1, 'omitnan'); err = std(data,0,1,'omitnan')./sqrt(sum(~isnan(data),1));
    
    col = cols_seq(mod(ss-1,2)+1, :);
    stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
    lbl = sprintf('Seq Set %d (Ch:%s)', ss, num2str(stimCh));
    
    errorbar(Amps_seq, mu, err, '-s', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 0, 'DisplayName', lbl);
end

xlabel('Amplitude (µA)', 'FontWeight','bold'); ylabel('Mean Peak Latency (ms)', 'FontWeight','bold');
title('Latency vs. Amplitude', 'FontWeight','bold');
legend('Location','best'); box off;


%% ============================================================
%   STATISTICS: Two-Way ANOVA for Peak Latency
%   Factors: Stimulation Type (Sim vs Seq) & Amplitude
% ============================================================

% 1. Prepare Variables for ANOVA Table
% We need to flatten the data into long columns:
% Y (Latency) | Group (Sim/Seq1/Seq2) | Amp (3,5,6,10)

y_values = [];
g_stim   = []; % Grouping variable for Stim Type
g_amp    = []; % Grouping variable for Amplitude

% --- A. Collect Simultaneous Data ---
for ai = 1:length(Amps_sim)
    % Extract latencies for this amplitude (remove NaNs)
    lats = LatAmp_sim(:, ai);
    lats = lats(~isnan(lats));
    
    % Append to master lists
    y_values = [y_values; lats];
    g_stim   = [g_stim; repmat({'Simultaneous'}, length(lats), 1)];
    g_amp    = [g_amp;  repmat(Amps_sim(ai), length(lats), 1)];
end

% --- B. Collect Sequential Data (Loop Sets) ---
% Note: Uses the same 'p_idx_trend' (PTD) defined in the Trend Plot section
for ss = 1:nSets_seq
    group_name = sprintf('Seq Set %d', ss);
    
    for ai = 1:length(Amps_seq)
        % Extract latencies
        lats = squeeze(LatAmp_seq(:, ai, ss));
        lats = lats(~isnan(lats));
        
        % Append
        y_values = [y_values; lats];
        g_stim   = [g_stim; repmat({group_name}, length(lats), 1)];
        g_amp    = [g_amp;  repmat(Amps_seq(ai), length(lats), 1)];
    end
end

% 2. Run N-Way ANOVA
% 'varnames' labels the axes in the output figure
[p_vals, tbl, stats] = anovan(y_values, {g_stim, g_amp}, ...
    'model', 'interaction', ...
    'varnames', {'StimType', 'Amplitude'}, ...
    'display', 'on');

% 3. Display P-Values in Command Window
fprintf('\n--- ANOVA P-VALUES ---\n');
fprintf('Stimulation Type:  %.5f  (Is Seq different from Sim?)\n', p_vals(1));
fprintf('Amplitude:         %.5f  (Does Amp affect Latency?)\n', p_vals(2));
fprintf('Interaction:       %.5f  (Does the gap change with Amp?)\n', p_vals(3));

% 4. Post-Hoc Pairwise Comparison (The "Proof" Plot)
% This figure shows exactly which groups are significantly different.
% Non-overlapping lines = Significant difference.
figure('Color','w', 'Name', 'Pairwise Comparison');
results = multcompare(stats, 'Dimension', 1);
title('Pairwise Comparison: Mean Latency by Group');
ylabel('Group Name'); xlabel('Mean Latency (ms)');

%% ==================== HELPER ====================
function peak_time_ms = get_peak_latency(tr_ids, trig, sp_data, search_win, bin_ms, kern, FS)
    nTr = numel(tr_ids); if nTr == 0, peak_time_ms = NaN; return; end
    all_spikes = [];
    for k = 1:nTr, tr=tr_ids(k); t0=trig(tr)/FS*1000; tt=sp_data(:,1)-t0; all_spikes=[all_spikes; tt(tt>=-20 & tt<=80)]; end
    
    bin_s = bin_ms/1000; edges = -20:bin_ms:80; t_centers = edges(1:end-1)+bin_ms/2;
    h = histcounts(all_spikes, edges); rate_smooth = conv(h/(nTr*bin_s), kern, 'same');
    
    valid_idx = find(t_centers >= search_win(1) & t_centers <= search_win(2));
    valid_rate = rate_smooth(valid_idx); valid_time = t_centers(valid_idx);
    
    if isempty(valid_rate) || max(valid_rate)==0, peak_time_ms = NaN; else
        [~, max_idx] = max(valid_rate); peak_time_ms = valid_time(max_idx);
    end
end
function sp = load_ssd_spikes(folder)
    cd(folder); f=dir('*sp_xia_SSD.mat'); if isempty(f), error('No SSD file'); end
    S=load(f(1).name); if isfield(S,'sp_corr'), sp=S.sp_corr; elseif isfield(S,'sp_SSD'), sp=S.sp_SSD; else, sp=S.sp_in; end
end