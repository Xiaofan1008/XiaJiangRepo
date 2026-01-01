%% Spike injection: Data-Driven Thinning (Simultaneous Proxy)
%  - Uses 1-5uA Simultaneous data to learn recruitment ratios.
%  - FIX: Enforces MONOTONICITY (cummax) so curves always rise.
%  - Extrapolates to est. 7uA max using slope logic.
%  - Thins 7uA Single spikes to match target density.
%  - FIX: Skips injection for 0uA/Baseline trials.
clear all;
%% User parameters
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Electrode_test'; % Source (7uA)
sequential_folder = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Exp1_Seq4';    % Target (1-5uA)
simultaneous_folder = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Exp1_Sim4'; % Reference (1-5uA)
Source_Amp_uA     = 7;          % The amplitude of the SINGLE pulse data (Source)
pulse_offset_ms   = 0;          % shift for injection in sequential data
use_fallback      = true;       % if no first spike, use PTD + 2 ms
FS                = 30000;
Electrode_type    = 1;
win_start_ms      = 0;          % injection window always starts at 2 ms
calc_win_ms       = [2 10];     % Window to calculate ratios
%% ============================================================
%   STEP 1: Load Simultaneous Data (To Learn Ratios)
% ============================================================
fprintf('Loading Simultaneous data to learn recruitment curves...\n');
cd(simultaneous_folder);
sim_sp_file = dir('*sp_xia_SSD.mat');
if isempty(sim_sp_file), error('No Simultaneous sp_xia file found'); end
load(sim_sp_file(1).name, 'sp_corr'); sp_sim = sp_corr;
S_sim = load(dir('*_exp_datafile_*.mat').name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_sim = S_sim.StimParams; simN_sim = S_sim.simultaneous_stim; E_MAP_sim = S_sim.E_MAP;
if isfield(S_sim, 'n_Trials'), nTr_sim = S_sim.n_Trials; else, nTr_sim = (size(StimParams_sim, 1)-1)/simN_sim; end
% Extract Amps
trialAmps_sim = cell2mat(StimParams_sim(2:end,16)); trialAmps_sim = trialAmps_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;
% Extract Sets
stimNames_sim = StimParams_sim(2:end,1); [~, idx_all_sim] = ismember(stimNames_sim, E_MAP_sim(2:end));
comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim, rr = (t-1)*simN_sim + (1:simN_sim); v = idx_all_sim(rr); v = v(v>0); comb_sim(t,1:numel(v)) = v(:).'; end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable'); nSets_sim = size(uniqueComb_sim,1);
% Load Triggers
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_sim = loadTrig(0);
% --- Calculate Ratios per Channel/Set/Amp ---
% Matrix: [nCh x nAmps x nSets]
Ratio_Map = zeros(numel(sp_sim), length(Amps_sim), nSets_sim);
fprintf('Computing ratios with slope extrapolation (Monotonicity Enforced)...\n');
for ss = 1:nSets_sim
    for ch = 1:numel(sp_sim)
        % 1. Get Mean Count for each Amp (1-5uA)
        mean_counts = zeros(1, length(Amps_sim));
        for ai = 1:length(Amps_sim)
            tr_ids = find(combClass_sim == ss & ampIdx_sim == ai);
            if isempty(tr_ids), continue; end
            
            % Count spikes
            tot = 0;
            for k = 1:numel(tr_ids)
                t0 = trig_sim(tr_ids(k))/FS*1000;
                tt = sp_sim{ch}(:,1) - t0;
                tot = tot + sum(tt >= calc_win_ms(1) & tt <= calc_win_ms(2));
            end
            mean_counts(ai) = tot / numel(tr_ids);
        end
        
        % --- [FIX] FORCE MONOTONICITY ---
        % Ensure curve never drops (e.g. [10 9 12] -> [10 10 12])
        mean_counts = cummax(mean_counts);
        
        % 2. Estimate 7uA Count (Slope Extrapolation)
        % Assume Amps_sim is [1, 2, 3, 4, 5]
        idx_4uA = find(abs(Amps_sim - 4) < 0.1);
        idx_5uA = find(abs(Amps_sim - 5) < 0.1);
        
        est_7uA = max(mean_counts); % Default fallback
        
        if ~isempty(idx_4uA) && ~isempty(idx_5uA)
            val4 = mean_counts(idx_4uA);
            val5 = mean_counts(idx_5uA);
            
            % Slope is guaranteed >= 0 because of cummax
            slope = (val5 - val4) / (5 - 4);
            est_7uA = val5 + slope * (7 - 5);
        end
        if est_7uA < 0.1, est_7uA = 0.1; end % Avoid div by zero
        
        % 3. Store Ratios
        Ratio_Map(ch, :, ss) = mean_counts / est_7uA;
    end
end
% Clip Ratios > 1
Ratio_Map(Ratio_Map > 1) = 1;
%% ============================================================
%   STEP 2: Load Single-Pulse Data (Source)
% ============================================================
cd(single_folder);
single_sp_file = dir('*sp_xia_SSD.mat');
assert(~isempty(single_sp_file),'No single sp_xia file found');
load(single_sp_file(1).name, 'sp_corr'); sp_single = sp_corr;
nChn = numel(sp_single); trig_single = loadTrig(0);
S1 = load(dir('*_exp_datafile_*.mat').name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_single = S1.StimParams; E_MAP = S1.E_MAP; n_Trials_single = S1.n_Trials; simultaneous_stim1 = S1.simultaneous_stim;
trialAmps_single = cell2mat(StimParams_single(2:end,16)); trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);
stimNames_single = StimParams_single(2:end,1); [~, idx_all_single] = ismember(stimNames_single, E_MAP(2:end));
stimChPerTrial_single = cell(n_Trials_single,1);
for t = 1:n_Trials_single, rr = (t-1)*simultaneous_stim1 + (1:simultaneous_stim1); v  = idx_all_single(rr); v  = v(v>0); stimChPerTrial_single{t} = v(:).'; end
%% ============================================================
%   STEP 3: Load Sequential Data (Target)
% ============================================================
cd(sequential_folder);
seq_sp_file = dir('*sp_xia.mat');
assert(~isempty(seq_sp_file),'No sequential sp_xia file found');
% load(seq_sp_file(1).name, 'sp'); sp_seq = sp;
load(seq_sp_file(1).name, 'sp_clipped');sp_seq = sp_clipped;
trig_seq = loadTrig(0);
S2 = load(dir('*_exp_datafile_*.mat').name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_seq = S2.StimParams; n_Trials_seq = S2.n_Trials; simultaneous_stim2 = S2.simultaneous_stim;
trialAmps_seq = cell2mat(StimParams_seq(2:end,16)); trialAmps_seq = trialAmps_seq(1:simultaneous_stim2:end);
stimNames_seq = StimParams_seq(2:end,1); [~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));
stimChPerTrial_seq = cell(n_Trials_seq,1);
for t = 1:n_Trials_seq, rr = (t-1)*simultaneous_stim2 + (1:simultaneous_stim2); v  = idx_all_seq(rr); v  = v(v>0); stimChPerTrial_seq{t} = v(:).'; end
comb_seq = zeros(n_Trials_seq, simultaneous_stim2); for t = 1:n_Trials_seq, v = stimChPerTrial_seq{t}; comb_seq(t, 1:numel(v)) = v; end
[uniqueComb_seq, ~, combClass_seq] = unique(comb_seq, 'rows', 'stable'); nSeqSets = size(uniqueComb_seq,1);
PTD_us_all = cell2mat(StimParams_seq(2:end, 6)); PTD_us = PTD_us_all(2:simultaneous_stim2:end); PTD_ms = PTD_us / 1000;
isSimultaneous = (PTD_ms == 0);
fst_file = dir('*FirstSpikeTimes.mat');
assert(~isempty(fst_file), 'Missing FirstSpikeTimes.mat');
load(fst_file(1).name, 'firstSpikeTimes');
%% ============================================================
%   STEP 4: Injection Loop
% ============================================================
unique_amps = unique(trialAmps_seq);
unique_PTDs = unique(PTD_ms);
total_spikes_added_all = 0;
for a = 1:numel(unique_amps)
    amp_val = unique_amps(a);
    
    % [FIX] Skip injection for 0uA (Baseline) or -1 (No Stim)
    if amp_val <= 0
        fprintf('Skipping Amplitude %g µA (Baseline/Control) - No injection needed.\n', amp_val);
        continue;
    end
    
    % Find Ratio Index for this amplitude
    ratio_amp_idx = find(abs(Amps_sim - amp_val) < 0.1);
    if isempty(ratio_amp_idx)
        fprintf('Warning: Seq Amp %g uA not found in Sim data. Using P_keep=1.0\n', amp_val);
    end
    
    total_spikes_added_amp = 0;
    
    % Source is always 7uA
    mask_single_amp = (trialAmps_single == Source_Amp_uA);
    
    for set_id = 1:nSeqSets
        stimVec = uniqueComb_seq(set_id,:);
        stimVec = stimVec(stimVec > 0);
        if isempty(stimVec), continue; end        
        ch_first_stim = stimVec(1);
        
        mask_single_chan = cellfun(@(x) ismember(ch_first_stim, x), stimChPerTrial_single);
        single_trials_for_group = find(mask_single_amp & mask_single_chan);
        
        if isempty(single_trials_for_group), continue; end
        
        for p = 1:numel(unique_PTDs)
            ptd_val = unique_PTDs(p);
            if ptd_val == 0, continue; end
            
            mask_seq_group = (trialAmps_seq == amp_val) & (combClass_seq == set_id) & (PTD_ms == ptd_val) & ~isSimultaneous;
            seq_trials_group = find(mask_seq_group);
            if isempty(seq_trials_group), continue; end          
            
            for gi = 1:numel(seq_trials_group)
                tr_seq = seq_trials_group(gi);
                idx_single = mod(gi-1, numel(single_trials_for_group)) + 1;
                tr_single  = single_trials_for_group(idx_single);
                
                t0_seq_ms    = trig_seq(tr_seq)/FS*1000;
                t0_single_ms = trig_single(tr_single)/FS*1000;
                
                for rec_ch = 1:nChn
                    spikes_src = sp_single{rec_ch};
                    if isempty(spikes_src), continue; end
                    
                    % Determine Window
                    first_ms = NaN;
                    if rec_ch <= numel(firstSpikeTimes) && tr_seq <= numel(firstSpikeTimes{rec_ch})
                        first_ms = firstSpikeTimes{rec_ch}(tr_seq);
                    end
                    win_start = win_start_ms;
                    if isfinite(first_ms) && first_ms > win_start
                        win_end = first_ms;
                    else
                        if ~use_fallback, continue; end
                        win_end = ptd_val + 2;
                    end
                    if win_end <= win_start, continue; end
                    
                    % 1. Delete Existing
                    seq_spikes_ch = sp_seq{rec_ch};
                    if ~isempty(seq_spikes_ch)
                        rel_seq = seq_spikes_ch(:,1) - t0_seq_ms;
                        del_mask = (rel_seq >= win_start) & (rel_seq < win_end);
                        if any(del_mask)
                            seq_spikes_ch(del_mask,:) = [];
                            sp_seq{rec_ch} = seq_spikes_ch;
                        end
                    end
                    
                    % 2. Select Source
                    t_start = t0_single_ms + win_start;
                    t_end   = t0_single_ms + win_end;
                    in_win = (spikes_src(:,1) >= t_start) & (spikes_src(:,1) < t_end);
                    spikes_potential = spikes_src(in_win,:);
                    if isempty(spikes_potential), continue; end
                    
                    % 3. Apply Thinning
                    p_keep = 1.0;
                    if ~isempty(ratio_amp_idx)
                        p_keep = Ratio_Map(rec_ch, ratio_amp_idx, set_id);
                    end
                    
                    keep_mask = rand(size(spikes_potential,1), 1) < p_keep;
                    spikes_add = spikes_potential(keep_mask, :);
                    
                    if isempty(spikes_add), continue; end
                    
                    % 4. Inject
                    t_inject = t0_seq_ms + pulse_offset_ms;
                    spikes_add(:,1) = spikes_add(:,1) - t0_single_ms + t_inject;
                    sp_seq{rec_ch} = sortrows([sp_seq{rec_ch}; spikes_add], 1);
                    total_spikes_added_amp = total_spikes_added_amp + size(spikes_add,1);
                end
            end
        end
    end
    fprintf('Amplitude %g µA: %d spikes added.\n', amp_val, total_spikes_added_amp);
    total_spikes_added_all = total_spikes_added_all + total_spikes_added_amp;
end
%% Save
base_sp_name = seq_sp_file(1).name;
new_name = strrep(base_sp_name, '.sp_xia.mat', '.sp_xia_FirstPulse.mat');
save(fullfile(sequential_folder, new_name), 'sp_seq', '-v7.3');
fprintf('Saved injected sequential file: %s\n', new_name);