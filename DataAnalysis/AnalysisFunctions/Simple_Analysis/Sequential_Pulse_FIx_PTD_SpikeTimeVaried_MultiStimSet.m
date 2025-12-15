%% === SPIKE INJECTION: Correct ORDER-SENSITIVE version ===
clear all;

%% === USER PARAMETERS === %%
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Single6_251125_180744';
sequential_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq6_5ms_251125_182437';
pulse_offset_ms   = 0;     
fallback_win_end  = 7;     
use_fallback      = true;
FS = 30000;
Electrode_type = 1;

%% ================================
%  LOAD SINGLE-PULSE DATA (SOURCE)
% ================================
cd(single_folder);
load(dir('*sp_xia.mat').name, 'sp_clipped'); 
sp_single = sp_clipped;
nChn = numel(sp_single);
trig_single = loadTrig(0);

S1 = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials','n_REP');
StimParams_single  = S1.StimParams;
E_MAP              = S1.E_MAP;
n_Trials_single    = S1.n_Trials;
simultaneous_stim1 = S1.simultaneous_stim;

trialAmps_single = cell2mat(StimParams_single(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);

stimNames_single  = StimParams_single(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP(2:end));
stimChPerTrial_single = cell(n_Trials_single,1);
for t = 1:n_Trials_single
    rr = (t-1)*simultaneous_stim1 + (1:simultaneous_stim1);
    stimChPerTrial_single{t} = unique(idx_all_single(rr)); 
end

%% ================================
%  LOAD SEQUENTIAL DATA (TARGET)
% ================================
cd(sequential_folder);
load(dir('*sp_xia.mat').name, 'sp_clipped'); 
sp_seq  = sp_clipped;
trig_seq = loadTrig(0);

S2 = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials','n_REP');
StimParams_seq  = S2.StimParams;
n_Trials_seq    = S2.n_Trials;
simultaneous_stim2 = S2.simultaneous_stim;

trialAmps_seq   = cell2mat(StimParams_seq(2:end,16));
trialAmps_seq   = trialAmps_seq(1:simultaneous_stim2:end);

stimNames_seq   = StimParams_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));
stimChPerTrial_seq = cell(n_Trials_seq,1);
for t = 1:n_Trials_seq
    rr = (t-1)*simultaneous_stim2 + (1:simultaneous_stim2);
    v  = idx_all_seq(rr);     % keep order
    v  = v(v > 0);            % drop any zeros
    stimChPerTrial_seq{t} = v(:).';
end

% Determine stimulation set ID for each sequential trial (order-sensitive)
comb_seq = zeros(n_Trials_seq, simultaneous_stim2);
for t = 1:n_Trials_seq
    v = stimChPerTrial_seq{t};
    comb_seq(t,1:numel(v)) = v;  % preserve ORDER
end
[uniqueComb_seq, ~, combClass_seq] = unique(comb_seq, 'rows', 'stable');
nSeqSets = size(uniqueComb_seq,1);

%% === Extract POST-TRIGGER DELAY (PTD) for sequential pulses ===
PTD_us_all = cell2mat(StimParams_seq(2:end, 6));   % µs for each pulse row
PTD_us     = PTD_us_all(2:simultaneous_stim2:end); % one PTD per TRIAL (first pulse row)
PTD_ms     = PTD_us / 1000;                        % convert to ms

%% Load first-spike time file
fst_file = dir('*FirstSpikeTimes.mat');
assert(~isempty(fst_file), 'Missing FirstSpikeTimes.mat');
load(fst_file(1).name, 'firstSpikeTimes');

%% Channel mapping
d = Depth_s(Electrode_type);

%% ================================
%  MAIN SPIKE INJECTION
% ================================
unique_amps = unique(trialAmps_seq);
total_spikes_added_all = 0;
singleIdxCounter = ones(nSeqSets,1);   % 1 counter per seq stimulus set
for a = 1:numel(unique_amps)
    amp = unique_amps(a);
    trials_seq_this_amp = find(trialAmps_seq == amp);
    total_spikes_added_amp = 0;

    for ti = 1:numel(trials_seq_this_amp)
        tr_seq = trials_seq_this_amp(ti);

        % Determine FIRST electrode of sequential stimulation (ORDER)
        stimSet_seq = stimChPerTrial_seq{tr_seq};
        if isempty(stimSet_seq), continue; end
        ch_first_stim = stimSet_seq(1);   % ALWAYS first in order
        % Find matching SINGLE trial for same amplitude + same channel
    
        match_single = find( cellfun(@(x) ismember(ch_first_stim, x), stimChPerTrial_single) & ...
                             (trialAmps_single == amp) );
        if isempty(match_single)
            warning('No single trial for stim %d at %duA (seq trial %d).', ...
                     ch_first_stim, amp, tr_seq);
            continue;
        end
        % Safe cyclic reuse of single trials
        % tr_single = match_single( mod(ti-1, numel(match_single)) + 1 );
        % Determine which seq set this trial belongs to
        seqSetID = combClass_seq(tr_seq);
        % Cycle INSIDE this stim set only
        k = singleIdxCounter(seqSetID);
        tr_single = match_single(k);        
        % update counter for this set only
        k = k + 1;
        if k > numel(match_single)
            k = 1;
        end
        singleIdxCounter(seqSetID) = k;
        
        % Trigger times
        t0_seq_ms    = trig_seq(tr_seq)/FS*1000;
        t0_single_ms = trig_single(tr_single)/FS*1000;
    
        for rec_ch = 1:nChn
            % get channel-specific first-spike time in seq
            first_ms = NaN;
            if rec_ch <= numel(firstSpikeTimes) && tr_seq <= numel(firstSpikeTimes{rec_ch})
                first_ms = firstSpikeTimes{rec_ch}(tr_seq);
            end
            if ~(isfinite(first_ms) && first_ms > 1)
                if use_fallback
                    % win_end = fallback_win_end;
                    win_end = PTD_ms(tr_seq) + 1;
                else
                    continue;
                end
            else
                win_end = first_ms;
            end

            % extract from SINGLE
            t_start = t0_single_ms;
            t_end   = t0_single_ms + win_end;

            spikes_src = sp_single{rec_ch};
            if isempty(spikes_src), continue; end

            in_win = spikes_src(:,1) >= t_start & spikes_src(:,1) < t_end;
            spikes_add = spikes_src(in_win,:);
            if isempty(spikes_add), continue; end

            % shift into seq timeline
            t_inject = t0_seq_ms + pulse_offset_ms;
            spikes_add(:,1) = spikes_add(:,1) - t0_single_ms + t_inject;

            % append
            sp_seq{rec_ch} = sortrows([sp_seq{rec_ch}; spikes_add], 1);
            total_spikes_added_amp = total_spikes_added_amp + size(spikes_add,1);
        end
    end

    fprintf('Amplitude %d µA: %d spikes added.\n', amp, total_spikes_added_amp);
    total_spikes_added_all = total_spikes_added_all + total_spikes_added_amp;
end

fprintf('\nTOTAL SPIKES ADDED: %d\n', total_spikes_added_all);

%% save
new_name = strrep(dir('*sp_xia.mat').name, '.sp_xia.mat', '.sp_xia_FirstPulse.mat');
save(fullfile(sequential_folder, new_name), 'sp_seq', '-v7.3');
fprintf('Saved injected seq file: %s\n', new_name);