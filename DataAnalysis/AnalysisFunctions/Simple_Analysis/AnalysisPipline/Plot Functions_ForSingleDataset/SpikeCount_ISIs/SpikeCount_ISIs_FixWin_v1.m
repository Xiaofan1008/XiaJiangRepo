%% ============================================================
%   Spike Count Analysis (Fixed Macro-Window Method)
%   - Metric: Mean Spike Count per Trial (Raw)
%   - Logic: Counts spikes in a single, wide window [start, end] for ALL ISIs.
%   - FILTER: Evaluates responding channels PER SET for TARGET AMP only.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% Single dataset folder containing both Sim (PTD=0) and Seq
data_folder = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_1'; 
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% Select which ISIs (PTDs) you want to analyze (e.g., [0, 5, 10, 15])
target_ISIs = [0, 5, 8, 10, 12, 15, 17, 20, 25]; 

% [MODIFIED] Choose the FIXED macro window [start, end]
% This wide window is used for EVERY condition to ensure 100% fairness
fixed_macro_win_ms = [2, 50]; 

% Choose which amplitude to run the Statistical Comparison on
target_Amp_Stats = 10; % uA

% Choose whether to involve bad trials
exclude_bad_trials = false; % true = remove them, false = keep all trials

% 2. Plotting & PSTH
psth_win_ms = [-50 100]; 
FS = 30000;             
bin_ms     = 1;
sigma_bins = 3;        
jitter_width = 0.2; scatter_alpha = 0.4; dot_size = 20;           

% Kernels
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

%% ================= 2. IDENTIFY POPULATION PER SET =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
% Store responding channels individually for each set
resp_channels_per_set = cell(nSets, 1);
% Find the index for the target amplitude to filter channels strictly
ai_target = find(abs(Amps - target_Amp_Stats) < 0.001);
if isempty(ai_target)
    error('Target Amplitude %.1f uA not found in dataset.', target_Amp_Stats);
end
fprintf('Analyzing Responding Channels (PER SET, %.1fuA ONLY) for ISIs %s ms:\n', target_Amp_Stats, num2str(target_ISIs));
for ss = 1:nSets
    local_resp_mask = false(nCh_Total, 1);
    
    if ss <= numel(R.set)
        % Only loop through the specific target amplitude (ai_target)
        if ai_target <= numel(R.set(ss).amp)
            for pi=1:numel(R.set(ss).amp(ai_target).ptd)
                curr_ptd = R.set(ss).amp(ai_target).ptd(pi).PTD_ms;
                
                % Check against target_ISIs list
                is_target_seq = any(abs(target_ISIs - curr_ptd) < 0.001);
                if ~is_target_seq
                    continue; 
                end
                
                this = R.set(ss).amp(ai_target).ptd(pi).channel; 
                for ch=1:min(length(this),nCh_Total)
                    if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch)=true; 
                    end
                end
            end
        end
    end
    
    resp_channels_per_set{ss} = find(local_resp_mask); 
    
    % Print responding channels specifically for this set
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    fprintf('  Set %d (Ch:%s) -> %d Channels\n', ss, num2str(stimCh), length(resp_channels_per_set{ss}));
end

%% =================== 3. COMPUTE SPIKE COUNTS =================
% Init Output Arrays based on total channels, NaN will handle non-responding
SpikeCount_sim = nan(nCh_Total, length(Amps), nSets);
nBins = length(edges_psth)-1;
PSTH_Sim = nan(nCh_Total, nBins, length(Amps), nSets);
SpikeCount_seq = nan(nCh_Total, length(Amps), nSets, numel(PTDs_ms));
PSTH_Seq = nan(nCh_Total, nBins, length(Amps), nSets, numel(PTDs_ms));
% Loop through sets first, then only process that set's responding channels
for ss = 1:nSets
    current_set_channels = resp_channels_per_set{ss};
    
    for ci = 1:length(current_set_channels)
        ch_idx = current_set_channels(ci); 
        recCh = d(ch_idx);    
        S_ch = sp{recCh};
        
        % Logic to toggle Bad Trials
        bad_trs = []; 
        if exclude_bad_trials && ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
            bad_trs = QC.BadTrials{ch_idx}; 
        end
        
        is_bad_ch = false;
        if exclude_bad_trials && ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
            is_bad_ch = true;
        end
        if is_bad_ch, continue; end
        
        % --- SIMULTANEOUS (PTD = 0) ---
        ptd_sim_idx = find(abs(PTDs_ms - 0) < 0.001);
                      
        if ~isempty(ptd_sim_idx) && any(abs(target_ISIs - 0) < 0.001)
            for ai = 1:length(Amps)
                tr_ids = find(combClass == ss & ampIdx == ai & ptdIdx == ptd_sim_idx);
                tr_ids = setdiff(tr_ids, bad_trs); 
                
                if isempty(tr_ids), continue; end
                
                % [MODIFIED] Pass the exact same fixed window for 0ms
                [count_val, psth_curve] = get_spike_count_macro(tr_ids, trig, S_ch, ...
                    fixed_macro_win_ms, edges_psth, g_sym, bin_s, FS);
                    
                SpikeCount_sim(ch_idx,ai,ss) = count_val; 
                PSTH_Sim(ch_idx, :, ai, ss)  = psth_curve;
            end
        end
        
        % --- SEQUENTIAL (Target ISIs > 0) ---
        for p = 1:numel(PTDs_ms)
            
            if PTDs_ms(p) == 0, continue; end % Handled above
            if ~any(abs(target_ISIs - PTDs_ms(p)) < 0.001), continue; end
            
            for ai = 1:length(Amps)
                tr_ids = find(combClass==ss & ptdIdx==p & ampIdx==ai);
                tr_ids = setdiff(tr_ids, bad_trs); 
                
                if isempty(tr_ids), continue; end
                
                % [MODIFIED] Pass the exact same fixed window for Seq ISIs
                [count_val, psth_curve] = get_spike_count_macro(tr_ids, trig, S_ch, ...
                    fixed_macro_win_ms, edges_psth, g_sym, bin_s, FS);
                
                SpikeCount_seq(ch_idx,ai,ss,p) = count_val; 
                PSTH_Seq(ch_idx, :, ai, ss, p) = psth_curve;
            end
        end
    end
end 

%% ===================== 4. PLOT (ISI Tuning Curve) ======================
figure('Color','w', 'Position',[100 100 800 600]); hold on;
% (ai_target is already found in Section 2)
% Define distinct colors for each SET
set_colors = lines(nSets); 
for ss = 1:nSets
    y_mean = nan(1, length(target_ISIs));
    y_sem  = nan(1, length(target_ISIs));
    
    % Gather the data across all chosen ISIs for this Set
    for p_idx = 1:length(target_ISIs)
        target_isi = target_ISIs(p_idx);
        
        if target_isi == 0
            data_set = squeeze(SpikeCount_sim(:, ai_target, ss));
        else
            p = find(abs(PTDs_ms - target_isi) < 0.001);
            if isempty(p)
                data_set = NaN;
            else
                data_set = squeeze(SpikeCount_seq(:, ai_target, ss, p));
            end
        end
        
        if ~all(isnan(data_set))
            y_mean(p_idx) = mean(data_set, 'omitnan');
            y_sem(p_idx)  = std(data_set, 0, 'omitnan') / sqrt(sum(~isnan(data_set)));
        end
    end
    
    % Plot the line for this Set
    if any(~isnan(y_mean))
        col = set_colors(ss, :);
        stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
        lbl = sprintf('Set %d (Ch:%s)', ss, num2str(stimCh));
        
        % Using standard errorbar for discrete intervals
        errorbar(target_ISIs, y_mean, y_sem, '-o', 'Color', col, 'LineWidth', 2, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', lbl);
    end
end
% Formatting
xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold', 'FontSize', 12); 
% [MODIFIED] Updated Y label to reflect the Macro Window
ylabel(sprintf('Spike Count (%.1f uA, Fixed Win: %.1f to %.1f ms)', target_Amp_Stats, fixed_macro_win_ms(1), fixed_macro_win_ms(2)), 'FontWeight','bold', 'FontSize', 12);
title('ISI Tuning Curve (Fixed Macro-Window Method)', 'FontWeight','bold', 'FontSize', 14);
% Force X-ticks to exactly match the tested ISIs
xticks(sort(target_ISIs));
legend('Location','best','Box','off'); box off;

%% ============================================================
%   5. STATISTICS: ANOVA (Compare ISIs at TARGET AMPLITUDE)
% ============================================================
y_val  = []; g_isi = []; 
fprintf('\n=== ANOVA RESULTS (Amplitude: %.1f uA) ===\n', target_Amp_Stats);
if ~isempty(ai_target)
    % Collect Sim Data
    if any(abs(target_ISIs - 0) < 0.001)
        for ss = 1:nSets
            d = squeeze(SpikeCount_sim(:, ai_target, ss)); d = d(~isnan(d));
            if ~isempty(d)
                y_val  = [y_val; d];
                g_isi = [g_isi; repmat({'ISI_0ms'}, length(d), 1)];
            end
        end
    end
    % Collect Seq Data
    for ss = 1:nSets
        for p = 1:numel(PTDs_ms)
            if PTDs_ms(p) == 0, continue; end
            if ~any(abs(target_ISIs - PTDs_ms(p)) < 0.001), continue; end
            
            label_str = sprintf('ISI_%dms', round(PTDs_ms(p)));
            d = squeeze(SpikeCount_seq(:, ai_target, ss, p)); d = d(~isnan(d));
            
            if ~isempty(d)
                y_val  = [y_val; d];
                g_isi = [g_isi; repmat({label_str}, length(d), 1)];
            end
        end
    end
    if ~isempty(y_val)
        [p_anova, tbl, stats] = anovan(y_val, {g_isi}, ...
            'varnames', {'ISI'}, 'display', 'on');
        fprintf('ISI Effect P-value: %.5f\n', p_anova(1));
        
        figure('Color','w','Name',sprintf('ISI Comparison @ %.1fuA', target_Amp_Stats));
        multcompare(stats, 'Dimension', 1);
    end
end

%% ============================================================
%   6. SAVE RESULTS
% ============================================================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_PTD_SpikeCount/DX016/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
isi_str = strjoin(string(target_ISIs), '_');
% [MODIFIED] Rename save file to reflect Fixed Macro-Window
out_filename = fullfile(save_dir, ['Result_SpikeCount_MacroWin_' char(isi_str) 'ms_' exp_id '.mat']);
ResultFR = struct();
ResultFR.Metadata.Created = datestr(now);
ResultFR.Metadata.Metric = 'Fixed Macro-Window Spike Count';
ResultFR.Metadata.FixedMacroWin = fixed_macro_win_ms;
ResultFR.Metadata.Amps = Amps;
ResultFR.Metadata.TargetISIs = target_ISIs;
ResultFR.Metadata.BadTrialsExcluded = exclude_bad_trials;
ResultFR.PeakFR.Sim = SpikeCount_sim; 
ResultFR.PeakFR.Seq = SpikeCount_seq;
ResultFR.PSTH.Sim = PSTH_Sim;
ResultFR.PSTH.Seq = PSTH_Seq;
save(out_filename, 'ResultFR');
fprintf('\n>>> Results Saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
% [MODIFIED] Simplified counting logic for the Fixed Macro-Window
function [count_val, psth_trace] = get_spike_count_macro(tr_ids, trig, sp_data, ...
    count_win, psth_edges, g_sym, bin_s, FS)
    
    nTr = numel(tr_ids);
    all_psth_spikes = [];
    total_spikes_in_window = 0;
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        % Count spikes strictly inside the provided static window
        mask_count = (tt >= count_win(1) & tt <= count_win(2));
        total_spikes_in_window = total_spikes_in_window + sum(mask_count);
        
        % Collect spikes for PSTH Visual
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