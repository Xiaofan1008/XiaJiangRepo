%% ============================================================
%   Spike Count Analysis (Combined Dataset: PTD 0=Sim, PTD X=Seq)
%   - Metric: Mean Spike Count per Trial (Raw)
%   - Logic: Counts spikes in window [2 20]ms
%   - FILTER: Extracts Sim (PTD=0) and Selected Seq ISIs from single dataset
%   - MODIFICATION: Uses RAW counts for all Union channels (No forced 0)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% Single dataset folder containing both Sim (PTD=0) and Seq
data_folder = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_1'; 
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% [NEW] Select which ISIs (PTDs) you want to analyze (e.g., [1, 5, 10])
target_ISIs = [0,5]; 

% 1. Analysis Window (Spike Counting)
post_win_ms = [2 25]; 
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

%% =================== 1. LOAD DATA (Single Source) ====================
[R, sp, trig, S, QC] = load_experiment_data(data_folder);

% --- Extract Stim Params ---
Stim = S.StimParams; 
simN = S.simultaneous_stim; 
E_MAP = S.E_MAP;
if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end

% Amplitudes
amps_all  = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
[Amps,~,ampIdx] = unique(trialAmps); Amps(Amps==-1) = 0;

% PTDs (Pulse Train Delays)
if simN > 1
    PTD_all_us = cell2mat(Stim(3:simN:end,6)); 
else
    PTD_all_us = zeros(nTr,1);
end
PTD_all_ms = PTD_all_us / 1000; 
[PTDs_ms,~,ptdIdx] = unique(PTD_all_ms);

% Parse Sets (Stimulation Combinations)
stimNames = Stim(2:end,1); 
[~, idx_all] = ismember(stimNames, E_MAP(2:end));
comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[uniqueComb,~,combClass] = unique(comb,'rows','stable'); 
nSets = size(uniqueComb,1);

%% ================= 2. IDENTIFY UNION POPULATION =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_channels_mask = false(nCh_Total, 1);

% Consolidate responsive channels from the single R structure
for si=1:numel(R.set)
    for ai=1:numel(R.set(si).amp)
        for pi=1:numel(R.set(si).amp(ai).ptd)
            % Check if this PTD is relevant (0 or target_ISIs)
            curr_ptd = R.set(si).amp(ai).ptd(pi).PTD_ms;
            
            % [MODIFIED] Check against target_ISIs list
            is_target_seq = any(abs(target_ISIs - curr_ptd) < 0.001);
            
            if abs(curr_ptd - 0) > 0.001 && ~is_target_seq
                continue; 
            end
            
            this = R.set(si).amp(ai).ptd(pi).channel; 
            for ch=1:min(length(this),nCh_Total)
                if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                    resp_channels_mask(ch)=true; 
                end
            end
        end
    end
end
resp_channels = find(resp_channels_mask); 
fprintf('Analyzing Fixed Population (Union of Sim & Seq %s ms): %d Channels\n', num2str(target_ISIs), length(resp_channels));

%% =================== 3. COMPUTE SPIKE COUNTS =================
% Init Output Arrays
SpikeCount_sim = nan(length(resp_channels), length(Amps), nSets);
nBins = length(edges_psth)-1;
PSTH_Sim = nan(length(resp_channels), nBins, length(Amps), nSets);

SpikeCount_seq = nan(length(resp_channels), length(Amps), nSets, numel(PTDs_ms));
PSTH_Seq = nan(length(resp_channels), nBins, length(Amps), nSets, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch_idx = resp_channels(ci); recCh = d(ch_idx);    
    S_ch = sp{recCh};
    
    % Get Bad Trials for this channel
    bad_trs = []; 
    if ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
        bad_trs = QC.BadTrials{ch_idx}; 
    end
    
    % --- SIMULTANEOUS (PTD = 0) ---
    ptd_sim_idx = find(abs(PTDs_ms - 0) < 0.001);
                  
    if ~isempty(ptd_sim_idx)
        for ss = 1:nSets
            % Bad Channel Check
            is_bad_ch = false;
            if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
                is_bad_ch = true;
            end
            
            for ai = 1:length(Amps)
                if is_bad_ch
                    SpikeCount_sim(ci,ai,ss) = NaN; PSTH_Sim(ci, :, ai, ss) = nan(1, nBins); continue;
                end
                
                % [MODIFIED] Removed "if ~is_resp" check.
                % We calculate RAW counts for all channels in the Union Population.
                
                % Filter trials: Set + Amp + PTD=0
                tr_ids = find(combClass == ss & ampIdx == ai & ptdIdx == ptd_sim_idx);
                tr_ids = setdiff(tr_ids, bad_trs); 
                
                if isempty(tr_ids), continue; end
                
                [count_val, psth_curve] = get_spike_count(tr_ids, trig, S_ch, ...
                    post_win_ms, edges_psth, g_sym, bin_s, FS);
                    
                SpikeCount_sim(ci,ai,ss) = count_val; 
                PSTH_Sim(ci, :, ai, ss)  = psth_curve;
            end
        end
    end
    
    % --- SEQUENTIAL (Target ISIs) ---
    for p = 1:numel(PTDs_ms)
        
        % [MODIFIED] FILTER: Process target ISIs only
        if ~any(abs(target_ISIs - PTDs_ms(p)) < 0.001)
            continue; 
        end
        
        for ss = 1:nSets
            is_bad_ch = false;
            if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
                is_bad_ch = true;
            end
            
            for ai = 1:length(Amps)
                if is_bad_ch
                    SpikeCount_seq(ci,ai,ss,p) = NaN; PSTH_Seq(ci, :, ai, ss, p) = nan(1, nBins); continue;
                end
                
                % [MODIFIED] Removed "if ~is_resp" check.
                
                % Filter trials: Set + Amp + PTD=p
                tr_ids = find(combClass==ss & ptdIdx==p & ampIdx==ai);
                tr_ids = setdiff(tr_ids, bad_trs); 
                
                if isempty(tr_ids), continue; end
                
                [count_val, psth_curve] = get_spike_count(tr_ids, trig, S_ch, ...
                    post_win_ms, edges_psth, g_sym, bin_s, FS);
                
                SpikeCount_seq(ci,ai,ss,p) = count_val; 
                PSTH_Seq(ci, :, ai, ss, p) = psth_curve;
            end
        end
    end
end 

%% ===================== 4. PLOT (Full Figure) ======================
figure('Color','w', 'Position',[100 100 800 600]); hold on;

% --- Plot Simultaneous (PTD=0) ---
sim_base_col = [0 0.3 0.8]; 
for ss = 1:nSets
    data_set = squeeze(SpikeCount_sim(:, :, ss));
    if all(all(isnan(data_set) | data_set==0)), continue; end 
    
    AvgSim = mean(data_set, 1, 'omitnan');
    SEMSim = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
    
    if nSets > 1, col = sim_base_col * (0.5 + 0.5*(ss/nSets)); else, col = sim_base_col; end
    
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    lbl = sprintf('Sim Set %d (Ch:%s)', ss, num2str(stimCh));
    
    if any(~isnan(AvgSim))
        for i = 1:length(Amps)
            valid = ~isnan(data_set(:,i));
            if any(valid)
                x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(i);
                scatter(x_jit, data_set(valid,i), dot_size, col, 'filled', 'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
            end
        end
        plot_shaded_error(Amps, AvgSim, SEMSim, col);
        plot(Amps, AvgSim, '-o', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, 'DisplayName', lbl);
    end
end

% --- Plot Sequential (Target ISIs) ---
% [MODIFIED] Define colors for different ISIs
% If you analyze many ISIs, this creates a distinct color for each.
isi_colors = lines(length(target_ISIs)); 

for p_idx = 1:length(target_ISIs)
    target_isi = target_ISIs(p_idx);
    
    % Find corresponding index in PTDs_ms array
    p = find(abs(PTDs_ms - target_isi) < 0.001);
    if isempty(p), continue; end
    
    for ss = 1:nSets
        data_set = squeeze(SpikeCount_seq(:, :, ss, p));
        AvgSeq = mean(data_set, 1, 'omitnan');
        SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
        
        % Assign color based on ISI (p_idx)
        col = isi_colors(p_idx, :);
        
        stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
        lbl = sprintf('Seq %dms Set %d (Ch:%s)', round(target_isi), ss, num2str(stimCh));
        
        if any(~isnan(AvgSeq))
            for i = 1:length(Amps)
                valid = ~isnan(data_set(:,i));
                if any(valid)
                    x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(i);
                    scatter(x_jit, data_set(valid,i), dot_size, col, 'filled', 'MarkerFaceAlpha', scatter_alpha, 'HandleVisibility','off');
                end
            end
            plot_shaded_error(Amps, AvgSeq, SEMSeq, col);
            plot(Amps, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', lbl);
        end
    end
end

xlabel('Amplitude (uA)', 'FontWeight','bold'); 
ylabel('Mean Spike Count (2-20 ms)', 'FontWeight','bold');
title('Response Magnitude (Total Spike Count)', 'FontWeight','bold');
legend('Location','best','Box','off'); box off;

%% ============================================================
%   5. STATISTICS: POOLED ANOVA (Sim vs Seq Selected ISIs)
% ============================================================
y_val  = []; g_stim = []; g_amp  = [];

% Sim Data
for ss = 1:nSets
    for ai = 1:length(Amps)
        d = squeeze(SpikeCount_sim(:, ai, ss)); d = d(~isnan(d));
        if isempty(d), continue; end
        y_val  = [y_val; d];
        g_stim = [g_stim; repmat({'Simultaneous'}, length(d), 1)];
        g_amp  = [g_amp;  repmat(Amps(ai), length(d), 1)];
    end
end

% Seq Data (Target ISIs)
for ss = 1:nSets
    for p = 1:numel(PTDs_ms)
        % [MODIFIED] Check against target ISIs
        if ~any(abs(target_ISIs - PTDs_ms(p)) < 0.001), continue; end
        
        label_str = sprintf('Seq_%dms', round(PTDs_ms(p)));
        
        for ai = 1:length(Amps)
            d = squeeze(SpikeCount_seq(:, ai, ss, p)); d = d(~isnan(d));
            if isempty(d), continue; end
            y_val  = [y_val; d];
            g_stim = [g_stim; repmat({label_str}, length(d), 1)];
            g_amp  = [g_amp;  repmat(Amps(ai), length(d), 1)];
        end
    end
end

if ~isempty(y_val)
    fprintf('\n=== ANOVA RESULTS (Pooled) ===\n');
    [p_anova, tbl, stats] = anovan(y_val, {g_stim, g_amp}, ...
        'model', 'interaction', 'varnames', {'StimType', 'Amplitude'}, 'display', 'on');
    fprintf('StimType P: %.5f\nAmplitude P: %.5f\nInteract P: %.5f\n', p_anova(1), p_anova(2), p_anova(3));
    
    figure('Color','w','Name','SpikeCount Comparison');
    multcompare(stats, 'Dimension', 1);
end

%% =============== PRINT SUMMARY STATISTICS ===============
fprintf('\n\n========================================================================\n');
fprintf('       SUMMARY STATISTICS: SPIKE COUNT         \n');
fprintf('========================================================================\n');
fprintf('Condition       \t');
for ai = 1:length(Amps), fprintf('%.0fuA_Mean SEM\t', Amps(ai)); end
fprintf('\n');

% Sim Rows
for ss = 1:nSets
    stimCh_List = uniqueComb(ss, :); stimCh_List = stimCh_List(stimCh_List > 0);
    fprintf('Sim_Set%d(Ch%s)\t', ss, num2str(stimCh_List));
    for ai = 1:length(Amps)
        d = squeeze(SpikeCount_sim(:, ai, ss));
        mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan') / sqrt(sum(~isnan(d)));
        if isnan(mu), fprintf('NaN\tNaN\t'); else, fprintf('%.4f %.4f\t', mu, sem); end
    end
    fprintf('\n');
end

% Seq Rows (Target ISIs)
for ss = 1:nSets
    for p = 1:numel(PTDs_ms)
        % [MODIFIED] Check against target ISIs
        if ~any(abs(target_ISIs - PTDs_ms(p)) < 0.001), continue; end
        
        stimCh_List = uniqueComb(ss, :); stimCh_List = stimCh_List(stimCh_List > 0);
        fprintf('Seq_Set%d(Ch%s_%dms)\t', ss, num2str(stimCh_List), round(PTDs_ms(p)));
        for ai = 1:length(Amps)
            d = squeeze(SpikeCount_seq(:, ai, ss, p)); 
            mu = mean(d, 'omitnan'); sem = std(d, 0, 'omitnan') / sqrt(sum(~isnan(d)));
            if isnan(mu), fprintf('NaN\tNaN\t'); else, fprintf('%.4f %.4f\t', mu, sem); end
        end
        fprintf('\n');
    end
end
fprintf('========================================================================\n');

%% ============================================================
%   6. SAVE RESULTS
% ============================================================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_PTD_SpikeCount/DX016/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};

% [MODIFIED] Generate filename based on selected ISIs
isi_str = strjoin(string(target_ISIs), '_');
out_filename = fullfile(save_dir, ['Result_SpikeCount_ISIs_' char(isi_str) 'ms_' exp_id '.mat']);

ResultFR = struct();
ResultFR.Metadata.Created = datestr(now);
ResultFR.Metadata.Metric = 'Mean Spike Count per Trial (Cleaned)';
ResultFR.Metadata.Window = post_win_ms;
ResultFR.Metadata.Amps = Amps;
ResultFR.Metadata.TargetISIs = target_ISIs;
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
        count_val = total_spikes_in_window / nTr; 
    else
        count_val = NaN;
    end
    
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