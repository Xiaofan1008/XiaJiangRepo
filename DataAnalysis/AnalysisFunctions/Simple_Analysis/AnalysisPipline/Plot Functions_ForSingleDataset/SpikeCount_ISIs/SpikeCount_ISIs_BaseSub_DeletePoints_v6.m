%% ============================================================
%   Spike Count Analysis (Fixed Macro-Window Method)
%   - Metric: Mean Spike Count per Trial per Channel (Classic)
%   - Logic: Counts spikes in a single, wide window [start, end] for ALL ISIs.
%   - FILTER: Evaluates responding channels PER SET using a UNION mask 
%             across all TARGET AMPS.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% Single dataset folder containing both Sim (PTD=0) and Seq
data_folder = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_4'; 
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% Select which ISIs (PTDs) you want to analyze. 
% [MODIFIED] Leave as [] to automatically process ALL tested ISIs.
target_ISIs = []; 

% Choose the FIXED macro window [start, end]
fixed_macro_win_ms = [0, 40]; 

% Added baseline window for spontaneous noise subtraction
baseline_win_ms = [-50, -10]; 

% Array of amplitudes to run the Comparison on. 
% [MODIFIED] Leave as [] to automatically process ALL tested Amplitudes.
target_Amps = [10]; % uA

% ARTIFACT SCRUBBER: List known noisy points to delete and bridge.
% Format: [Set, Amplitude (uA), ISI (ms)]; Add as many rows as needed.
manual_artifacts = [ 
  
]; 

% Choose whether to involve bad trials
exclude_bad_trials = true; % true = remove them, false = keep all trials

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

% =====================================================================
% [NEW MODIFICATION] AUTO-DETECT EMPTY TARGET ARRAYS
% =====================================================================
if isempty(target_Amps)
    target_Amps = Amps(:)';
    fprintf('Auto-detected target_Amps: %s uA\n', num2str(target_Amps));
end

if isempty(target_ISIs)
    % Ensure 0 is always included for Sim trials alongside sequential ISIs
    target_ISIs = sort(unique([0; PTDs_ms(:)]))'; 
    fprintf('Auto-detected target_ISIs: %s ms\n', num2str(target_ISIs));
end

%% ================= 2. IDENTIFY POPULATION PER SET =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_channels_per_set = cell(nSets, 1);

% Find indices for ALL target amplitudes
[~, ai_targets] = ismember(target_Amps, Amps);
if any(ai_targets == 0)
    error('One or more Target Amplitudes not found in dataset.');
end

fprintf('Analyzing Responding Channels (UNION MASK across %s uA) for ISIs %s ms:\n', num2str(target_Amps), num2str(target_ISIs));

for ss = 1:nSets
    local_resp_mask = false(nCh_Total, 1);
    
    if ss <= numel(R.set)
        for a_idx = 1:length(ai_targets)
            ai_targ = ai_targets(a_idx);
            
            if ai_targ <= numel(R.set(ss).amp)
                for pi=1:numel(R.set(ss).amp(ai_targ).ptd)
                    curr_ptd = R.set(ss).amp(ai_targ).ptd(pi).PTD_ms;
                    
                    if ~any(abs(target_ISIs - curr_ptd) < 0.001), continue; end
                    
                    this = R.set(ss).amp(ai_targ).ptd(pi).channel; 
                    for ch=1:min(length(this),nCh_Total)
                        if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                            local_resp_mask(ch)=true; 
                        end
                    end
                end
            end
        end
    end
    
    resp_channels_per_set{ss} = find(local_resp_mask); 
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    fprintf('  Set %d (Ch:%s) -> %d Channels\n', ss, num2str(stimCh), length(resp_channels_per_set{ss}));
end

%% =================== 3. COMPUTE SPIKE COUNTS =================
% Matrices are tightly condensed to match your Targets exactly
SpikeCount_sim = nan(nCh_Total, length(target_Amps), nSets);
SpikeCount_seq = nan(nCh_Total, length(target_Amps), nSets, length(target_ISIs));

nBins = length(edges_psth)-1;
PSTH_Sim = nan(nCh_Total, nBins, length(target_Amps), nSets);
PSTH_Seq = nan(nCh_Total, nBins, length(target_Amps), nSets, length(target_ISIs));

for ss = 1:nSets
    current_set_channels = resp_channels_per_set{ss};
    
    for ci = 1:length(current_set_channels)
        ch_idx = current_set_channels(ci); 
        recCh = d(ch_idx);    
        S_ch = sp{recCh};
        
        bad_trs = []; 
        if exclude_bad_trials && ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
            bad_trs = QC.BadTrials{ch_idx}; 
        end
        if exclude_bad_trials && ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
            continue; 
        end
        
        % --- SIMULTANEOUS (PTD = 0) ---
        p_sim_idx = find(target_ISIs == 0); 
        if ~isempty(p_sim_idx)
            ptd_sim_raw = find(abs(PTDs_ms - 0) < 0.001);
            if ~isempty(ptd_sim_raw)
                for a_idx = 1:length(target_Amps)
                    ai_raw = ai_targets(a_idx);
                    
                    tr_ids = find(combClass == ss & ampIdx == ai_raw & ptdIdx == ptd_sim_raw);
                    tr_ids = setdiff(tr_ids, bad_trs); 
                    
                    if isempty(tr_ids), continue; end
                    
                    [count_val, psth_curve] = get_spike_count_macro(tr_ids, trig, S_ch, ...
                        fixed_macro_win_ms, baseline_win_ms, edges_psth, g_sym, bin_s, FS);
                        
                    SpikeCount_sim(ch_idx, a_idx, ss) = count_val;
                    PSTH_Sim(ch_idx, :, a_idx, ss)    = psth_curve;
                end
            end
        end
        
        % --- SEQUENTIAL (Target ISIs > 0) ---
        for p_idx = 1:length(target_ISIs)
            curr_isi = target_ISIs(p_idx);
            if curr_isi == 0, continue; end 
            
            p_raw = find(abs(PTDs_ms - curr_isi) < 0.001);
            if isempty(p_raw), continue; end
            
            for a_idx = 1:length(target_Amps)
                ai_raw = ai_targets(a_idx);
                
                tr_ids = find(combClass==ss & ptdIdx==p_raw & ampIdx==ai_raw);
                tr_ids = setdiff(tr_ids, bad_trs); 
                
                if isempty(tr_ids), continue; end
                
                [count_val, psth_curve] = get_spike_count_macro(tr_ids, trig, S_ch, ...
                    fixed_macro_win_ms, baseline_win_ms, edges_psth, g_sym, bin_s, FS);
                
                SpikeCount_seq(ch_idx, a_idx, ss, p_idx) = count_val;
                PSTH_Seq(ch_idx, :, a_idx, ss, p_idx)    = psth_curve;
            end
        end
    end
end 

%% =================== 3.5 ARTIFACT SCRUBBER ===================
if ~isempty(manual_artifacts)
    fprintf('\nApplying Artifact Scrubber...\n');
    for r = 1:size(manual_artifacts, 1)
        scrub_set = manual_artifacts(r, 1);
        scrub_amp = manual_artifacts(r, 2);
        scrub_isi = manual_artifacts(r, 3);
        
        a_scrub_idx = find(abs(target_Amps - scrub_amp) < 0.001);
        p_scrub_idx = find(abs(target_ISIs - scrub_isi) < 0.001);
        
        if isempty(a_scrub_idx) || isempty(p_scrub_idx), continue; end
        
        if scrub_isi == 0
            if scrub_set <= nSets
                SpikeCount_sim(:, a_scrub_idx, scrub_set) = NaN;
                fprintf('  Scrubbed: Set %d | %.1f uA | %d ms\n', scrub_set, scrub_amp, scrub_isi);
            end
        else
            if scrub_set <= nSets
                SpikeCount_seq(:, a_scrub_idx, scrub_set, p_scrub_idx) = NaN;
                fprintf('  Scrubbed: Set %d | %.1f uA | %d ms\n', scrub_set, scrub_amp, scrub_isi);
            end
        end
    end
    fprintf('\n');
end

%% ===================== 4. PLOT (ISI Tuning Curve) ======================
amp_colors = lines(length(target_Amps)); 

for ss = 1:nSets
    figure('Color','w', 'Position',[100 100 800 600]); 
    hold on;
    
    stimCh = uniqueComb(ss,:); 
    stimCh = stimCh(stimCh>0);
    
    for a_idx = 1:length(target_Amps)
        current_amp = target_Amps(a_idx);
        
        y_mean = nan(1, length(target_ISIs));
        y_sem  = nan(1, length(target_ISIs));
        
        for p_idx = 1:length(target_ISIs)
            target_isi = target_ISIs(p_idx);
            
            if target_isi == 0
                data_set = squeeze(SpikeCount_sim(:, a_idx, ss));
            else
                data_set = squeeze(SpikeCount_seq(:, a_idx, ss, p_idx));
            end
            
            data_set = data_set(~isnan(data_set)); 
            
            if ~isempty(data_set)
                y_mean(p_idx) = mean(data_set);
                y_sem(p_idx)  = std(data_set) / sqrt(length(data_set));
            end
        end
        
        if any(~isnan(y_mean))
            col = amp_colors(a_idx, :);
            lbl = sprintf('%.1f uA', current_amp); 
            
            valid_idx = ~isnan(y_mean);
            plot_x = target_ISIs(valid_idx);
            plot_y = y_mean(valid_idx);
            
            plot(plot_x, plot_y, '-o', 'Color', col, 'LineWidth', 2, ...
                'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', lbl);
        end
    end
    
    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold', 'FontSize', 12); 
    ylabel('Net Mean Spikes', 'FontWeight','bold', 'FontSize', 12);
    title(sprintf('Set %d (Stim Channels: %s)', ss, num2str(stimCh)), 'FontWeight','bold', 'FontSize', 14);
    xticks(sort(target_ISIs));
    box off;
    lgd = legend('Location','best','Box','off'); 
    title(lgd, 'Amplitudes'); 
end

%% ============================================================
%   6. SAVE RESULTS
% ============================================================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX016/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};

% [MODIFIED] Smart Filename Generator to prevent OS Path Limit Crashes
if length(target_Amps) > 4
    amp_str = 'AllAmps';
else
    amp_str = strjoin(string(target_Amps), '_');
end

if length(target_ISIs) > 4
    isi_str = 'AllISIs';
else
    isi_str = strjoin(string(target_ISIs), '_');
end

out_filename = fullfile(save_dir, ['Result_SpikeCount_FixWin_DX016_' char(amp_str) 'uA_' exp_id '.mat']);

ResultFR = struct();
ResultFR.Metadata.Created = datestr(now);
ResultFR.Metadata.Metric = 'Classic Fixed Macro-Window Spike Count (Net Mean per Channel)';
ResultFR.Metadata.FixedMacroWin = fixed_macro_win_ms;
ResultFR.Metadata.BaselineWin = baseline_win_ms; 
ResultFR.Metadata.TargetAmps = target_Amps; 
ResultFR.Metadata.TargetISIs = target_ISIs;
ResultFR.Metadata.BadTrialsExcluded = exclude_bad_trials;
ResultFR.Metadata.Stimulation_Sets = uniqueComb;
ResultFR.Metadata.Responding_Channels = resp_channels_per_set;
ResultFR.Metadata.Manual_Artifacts_Scrubbed = manual_artifacts; 
ResultFR.Metadata.HardwareAmps = Amps;
ResultFR.Metadata.HardwarePTDs = PTDs_ms;

ResultFR.SpikeCounts.Sim = SpikeCount_sim; 
ResultFR.SpikeCounts.Seq = SpikeCount_seq;
ResultFR.PSTH.Sim = PSTH_Sim;
ResultFR.PSTH.Seq = PSTH_Seq;

save(out_filename, 'ResultFR');
fprintf('\n>>> Results Saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function [count_val, psth_trace] = get_spike_count_macro(tr_ids, trig, sp_data, ...
    count_win, base_win, psth_edges, g_sym, bin_s, FS)
    
    nTr = numel(tr_ids);
    all_psth_spikes = [];
    total_net_spikes_in_window = 0;
    
    dur_evoked = count_win(2) - count_win(1);
    dur_base   = base_win(2) - base_win(1);
    
    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        
        mask_count = (tt >= count_win(1) & tt <= count_win(2));
        evoked_count = sum(mask_count);
        
        mask_base = (tt >= base_win(1) & tt <= base_win(2));
        base_count = sum(mask_base);
        
        net_count = max(0, evoked_count - (base_count * (dur_evoked / dur_base))); 
        total_net_spikes_in_window = total_net_spikes_in_window + net_count;
        
        all_psth_spikes = [all_psth_spikes; tt(tt >= psth_edges(1) & tt <= psth_edges(end))];
    end
    
    if nTr > 0
        count_val = total_net_spikes_in_window / nTr; 
    else
        count_val = NaN;
    end
    
    h_psth = histcounts(all_psth_spikes, psth_edges);
    rate_psth = h_psth / (nTr * bin_s);
    psth_trace = conv(rate_psth, g_sym, 'same');
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*MultiISIRespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    
    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end
    
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);
    
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.MultiISIsBadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
    f_bt = dir('*.MultiISIsBadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end