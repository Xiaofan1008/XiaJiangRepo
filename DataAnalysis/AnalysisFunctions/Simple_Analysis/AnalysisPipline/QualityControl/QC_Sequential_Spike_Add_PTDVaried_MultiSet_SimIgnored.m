%% Spike injection: grouped by amplitude, sequence set, PTD
%  FIXED: Matches Stimulation Channels by NAME, not ID (Prevents cross-mapping errors)
clear all;
%% User parameters
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Single2';
sequential_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim2';
pulse_offset_ms   = 0;          % shift for injection in sequential data
use_fallback      = true;       % if no first spike, use PTD + 2 ms
FS                = 30000;
Electrode_type    = 2;
win_start_ms      = 0;          % injection window always starts at 2 ms

%% Load single-pulse data (source spikes)
cd(single_folder);
single_sp_file = dir('*sp_xia_SSD.mat');
assert(~isempty(single_sp_file),'No single sp_xia file found');
% load(single_sp_file(1).name, 'sp_clipped');
% sp_single   = sp_clipped;
load(single_sp_file(1).name, 'sp_corr');
sp_single   = sp_corr;
nChn        = numel(sp_single);
trig_single = loadTrig(0);
S1 = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_single  = S1.StimParams;
E_MAP_single       = S1.E_MAP; % [MODIFIED] Store as E_MAP_single
n_Trials_single    = S1.n_Trials;
simultaneous_stim1 = S1.simultaneous_stim;
trialAmps_single = cell2mat(StimParams_single(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);
stimNames_single  = StimParams_single(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP_single(2:end));
stimChPerTrial_single = cell(n_Trials_single,1);
for t = 1:n_Trials_single
    rr = (t-1)*simultaneous_stim1 + (1:simultaneous_stim1);
    v  = idx_all_single(rr);
    v  = v(v>0);
    stimChPerTrial_single{t} = v(:).';
end

%% Load sequential data (target spikes)
cd(sequential_folder);
seq_sp_file = dir('*sp_xia.mat');
assert(~isempty(seq_sp_file),'No sequential sp_xia file found');
load(seq_sp_file(1).name, 'sp_clipped');
sp_seq   = sp_clipped;
trig_seq = loadTrig(0);
S2 = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_seq    = S2.StimParams;
E_MAP_seq         = S2.E_MAP; % [MODIFIED] Store as E_MAP_seq
n_Trials_seq      = S2.n_Trials;
simultaneous_stim2 = S2.simultaneous_stim;
trialAmps_seq = cell2mat(StimParams_seq(2:end,16));
trialAmps_seq = trialAmps_seq(1:simultaneous_stim2:end);
stimNames_seq = StimParams_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP_seq(2:end));
stimChPerTrial_seq = cell(n_Trials_seq,1);
for t = 1:n_Trials_seq
    rr = (t-1)*simultaneous_stim2 + (1:simultaneous_stim2);
    v  = idx_all_seq(rr);
    v  = v(v>0);
    stimChPerTrial_seq{t} = v(:).';
end
comb_seq = zeros(n_Trials_seq, simultaneous_stim2);
for t = 1:n_Trials_seq
    v = stimChPerTrial_seq{t};
    comb_seq(t, 1:numel(v)) = v;
end
[uniqueComb_seq, ~, combClass_seq] = unique(comb_seq, 'rows', 'stable');
nSeqSets = size(uniqueComb_seq,1);
PTD_us_all = cell2mat(StimParams_seq(2:end, 6));
PTD_us     = PTD_us_all(2:simultaneous_stim2:end);
PTD_ms     = PTD_us / 1000;
isSimultaneous = (PTD_ms == 0);

%% Load first-spike times (from sequential data)
fst_file = dir('*FirstSpikeTimes.mat');
assert(~isempty(fst_file), 'Missing FirstSpikeTimes.mat');
load(fst_file(1).name, 'firstSpikeTimes');

%% Channel mapping (if you need depth mapping)
d = Depth_s(Electrode_type);

%% Grouped injection: amp -> set -> PTD -> trial -> channel
unique_amps = unique(trialAmps_seq);
unique_PTDs = unique(PTD_ms);
total_spikes_added_all = 0;

for a = 1:numel(unique_amps)
    amp_val = unique_amps(a);
    total_spikes_added_amp = 0;
    
    mask_single_amp = (trialAmps_single == amp_val);
    
    for set_id = 1:nSeqSets
        stimVec = uniqueComb_seq(set_id,:);
        stimVec = stimVec(stimVec > 0);
        if isempty(stimVec)
            continue;
        end        
        
        % 1. Get the Index in Seq E_MAP
        ch_first_idx_seq = stimVec(1);
        
        % 2. Get the Actual String Name (e.g., 'E-006')
        % Note: E_MAP usually has an empty first cell, so idx maps to E_MAP{idx+1}
        ch_name = E_MAP_seq{ch_first_idx_seq + 1};
        
        % 3. Find the Correct Index in Single E_MAP
        % Returns the index into E_MAP_single(2:end)
        ch_idx_in_single_map = find(strcmp(E_MAP_single(2:end), ch_name));
        
        if isempty(ch_idx_in_single_map)
            fprintf('Warning: Seq Channel %s not found in Single Data Map.\n', ch_name);
            continue;
        end
        
        % 4. Select Single Trials using the CORRECTED index
        % The idx_all_single logic above mapped names to 1..N indices matching E_MAP(2:end)
        mask_single_chan = cellfun(@(x) ismember(ch_idx_in_single_map, x), stimChPerTrial_single);
        
        single_trials_for_group = find(mask_single_amp & mask_single_chan);
        
        if isempty(single_trials_for_group)
            fprintf('No single trials for amp %g, Stim %s (SeqID %d)\n', amp_val, ch_name, set_id);
            continue;
        end
        
        for p = 1:numel(unique_PTDs)
            ptd_val = unique_PTDs(p);
            
            % skip PTD = 0 (simultaneous stimulation)
            if ptd_val == 0
                continue;
            end
            
            % sequential trials in this (amp, set, PTD) group
            mask_seq_group = (trialAmps_seq == amp_val) & ...
                             (combClass_seq == set_id) & ...
                             (PTD_ms == ptd_val) & ...
                             ~isSimultaneous;
            seq_trials_group = find(mask_seq_group);
            
            if isempty(seq_trials_group)
                continue;
            end          
            
            % loop trials in this group
            for gi = 1:numel(seq_trials_group)
                tr_seq = seq_trials_group(gi);
                
                % pick a matching single trial using deterministic cycling
                idx_single = mod(gi-1, numel(single_trials_for_group)) + 1;
                tr_single  = single_trials_for_group(idx_single);
                
                t0_seq_ms    = trig_seq(tr_seq)/FS*1000;
                t0_single_ms = trig_single(tr_single)/FS*1000;
                
                % per-channel spike injection for this trial
                for rec_ch = 1:nChn
                    spikes_src = sp_single{rec_ch};
                    if isempty(spikes_src)
                        continue;
                    end
                    
                    % trial-specific first spike latency in sequential data
                    first_ms = NaN;
                    if rec_ch <= numel(firstSpikeTimes) && ...
                       tr_seq <= numel(firstSpikeTimes{rec_ch})
                        first_ms = firstSpikeTimes{rec_ch}(tr_seq);
                    end
                    
                    win_start = win_start_ms;
                    if isfinite(first_ms) && first_ms > win_start
                        win_end = first_ms;
                    else
                        if ~use_fallback
                            continue;
                        end
                        win_end = ptd_val + 2;
                    end
                    
                    if win_end <= win_start
                        continue;
                    end
                    
                    % delete any existing sequential spikes in [2, win_end] relative to first pulse
                    seq_spikes_ch = sp_seq{rec_ch};
                    if ~isempty(seq_spikes_ch)
                        rel_seq = seq_spikes_ch(:,1) - t0_seq_ms;
                        del_mask = (rel_seq >= win_start) & (rel_seq < win_end);
                        if any(del_mask)
                            seq_spikes_ch(del_mask,:) = [];
                            sp_seq{rec_ch} = seq_spikes_ch;
                        end
                    end
                    
                    % injection window in single dataset
                    t_start = t0_single_ms + win_start;
                    t_end   = t0_single_ms + win_end;
                    
                    in_win = (spikes_src(:,1) >= t_start) & (spikes_src(:,1) < t_end);
                    spikes_add = spikes_src(in_win,:);
                    if isempty(spikes_add)
                        continue;
                    end
                    
                    % shift single spikes into sequential timeline
                    t_inject = t0_seq_ms + pulse_offset_ms;
                    spikes_add(:,1) = spikes_add(:,1) - t0_single_ms + t_inject;
                    
                    % append and sort
                    sp_seq{rec_ch} = sortrows([sp_seq{rec_ch}; spikes_add], 1);
                    total_spikes_added_amp = total_spikes_added_amp + size(spikes_add,1);
                end
            end
        end
    end
    
    fprintf('Amplitude %g µA: %d spikes added.\n', amp_val, total_spikes_added_amp);
    total_spikes_added_all = total_spikes_added_all + total_spikes_added_amp;
end
% fprintf('\nTotal spikes added across all amplitudes: %d\n', total_spikes_added_all);
%% Save new spike file
base_sp_name = seq_sp_file(1).name;
new_name = strrep(base_sp_name, '.sp_xia.mat', '.sp_xia_FirstPulse.mat');
save(fullfile(sequential_folder, new_name), 'sp_seq', '-v7.3');
fprintf('Saved injected sequential file: %s\n', new_name);