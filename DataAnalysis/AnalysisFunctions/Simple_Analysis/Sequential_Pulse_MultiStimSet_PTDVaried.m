%% === SPIKE INJECTION: ORDER-SENSITIVE + PTD-DEPENDENT WINDOW (REFINED) ===
clear all;

%% === USER PARAMETERS === %%
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Single4_251128_150137';
sequential_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim4_251128_150648';

pulse_offset_ms    = 0;     % extra shift if you want (usually 0)
blank_margin_ms    = 1;     % start injecting at least 1 ms AFTER PTD
fallback_len_ms    = 10;    % window length (after PTD) if no first spike
use_fallback       = true;  % still inject even if no first spike
FS                 = 30000;
Electrode_type     = 2;

%% ===================== LOAD SINGLE-PULSE DATA ===================== %%
cd(single_folder);
load(dir('*sp_xia.mat').name, 'sp_clipped'); 
sp_single = sp_clipped;
nChn = numel(sp_single);
trig_single = loadTrig(0);

S1 = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_single  = S1.StimParams;
E_MAP              = S1.E_MAP;
n_Trials_single    = S1.n_Trials;
simultaneous_stim1 = S1.simultaneous_stim;

% amplitudes in SINGLE
trialAmps_single = cell2mat(StimParams_single(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);

% which channel is stimulated in SINGLE
stimNames_single  = StimParams_single(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP(2:end));

stimChPerTrial_single = cell(n_Trials_single,1);
for t = 1:n_Trials_single
    rr = (t-1)*simultaneous_stim1 + (1:simultaneous_stim1);
    v  = idx_all_single(rr);
    v  = v(v > 0);
    % single stimulation: usually one channel in v
    stimChPerTrial_single{t} = unique(v); 
end


%% ===================== LOAD SEQUENTIAL DATA ===================== %%
cd(sequential_folder);
load(dir('*sp_xia.mat').name, 'sp_clipped'); 
sp_seq  = sp_clipped;
trig_seq = loadTrig(0);

S2 = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_seq  = S2.StimParams;
n_Trials_seq    = S2.n_Trials;
simultaneous_stim2 = S2.simultaneous_stim;

% amplitudes in SEQUENTIAL
trialAmps_seq = cell2mat(StimParams_seq(2:end,16));
trialAmps_seq = trialAmps_seq(1:simultaneous_stim2:end);

% stimulation set (ordered) in SEQUENTIAL
stimNames_seq = StimParams_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));

stimChPerTrial_seq = cell(n_Trials_seq,1);
for t = 1:n_Trials_seq
    rr = (t-1)*simultaneous_stim2 + (1:simultaneous_stim2);
    v  = idx_all_seq(rr);
    v  = v(v > 0);
    % ORDER preserved: e.g. [20 24] vs [24 20]
    stimChPerTrial_seq{t} = v(:).';
end

% order-sensitive stimulation sets for SEQ
comb_seq = zeros(n_Trials_seq, simultaneous_stim2);
for t = 1:n_Trials_seq
    v = stimChPerTrial_seq{t};
    comb_seq(t,1:numel(v)) = v;
end
[uniqueComb_seq, ~, combClass_seq] = unique(comb_seq, 'rows', 'stable');
nSeqSets = size(uniqueComb_seq,1);

% PTD PER TRIAL (µs → ms)
PTD_us_all = cell2mat(StimParams_seq(2:end, 6));
PTD_us     = PTD_us_all(2:simultaneous_stim2:end);
PTD_ms     = PTD_us / 1000;


%% ===================== LOAD FIRST-SPIKE TIMES ===================== %%
% firstSpikeTimes{ch}(trial) is the FIRST spike time relative to trigger (ms),
% found in [PTD, PTD+something] for that trial.
fst_file = dir('*FirstSpikeTimes.mat');
assert(~isempty(fst_file), 'Missing FirstSpikeTimes.mat');
load(fst_file(1).name, 'firstSpikeTimes');

% ===================== DEPTH MAPPING (if needed later) ===================== %%
d = Depth_s(Electrode_type); 


%% ===================== MAIN SPIKE INJECTION ===================== %%
unique_amps = unique(trialAmps_seq);
total_spikes_added_all = 0;

singleIdxCounter = ones(nSeqSets,1);

fprintf('\n=== SPIKE INSERTION (first-spike OR PTD+1 ms fallback) ===\n');

for a = 1:numel(unique_amps)
    amp = unique_amps(a);
    seq_trials_this_amp = find(trialAmps_seq == amp);
    spikes_added_amp = 0;

    for ti = 1:numel(seq_trials_this_amp)
        tr_seq = seq_trials_this_amp(ti);

        % FIRST stimulated channel (ordered SEQ)
        stimSet_seq = stimChPerTrial_seq{tr_seq};
        if isempty(stimSet_seq), continue; end
        ch_first_stim = stimSet_seq(1);

        % Find matching SINGLE trials (same amp + same first stim channel)
        match_single = find(trialAmps_single == amp & ...
                             cellfun(@(x) ismember(ch_first_stim, x), stimChPerTrial_single) );

        if isempty(match_single)
            warning('No matching single trial for seq trial %d (stim %d)', ...
                     tr_seq, ch_first_stim);
            continue;
        end

        % Cycle through matching single trials inside each SEQ stim set
        seqSetID = combClass_seq(tr_seq);
        k = singleIdxCounter(seqSetID);
        tr_single = match_single(k);

        k = k + 1; 
        if k > numel(match_single), k = 1; end
        singleIdxCounter(seqSetID) = k;

        % Trigger times
        t0_seq_ms    = trig_seq(tr_seq)/FS*1000;
        t0_single_ms = trig_single(tr_single)/FS*1000;
        % PTD of THIS seq trial
        ptd_ms = PTD_ms(tr_seq);

        % LOOP OVER RECORDING CHANNELS (correct position)
        for rec_ch = 1:nChn
            % Determine window_end_ms
            first_ms = NaN;

            if rec_ch <= numel(firstSpikeTimes)
                if tr_seq <= numel(firstSpikeTimes{rec_ch})
                    first_ms = firstSpikeTimes{rec_ch}(tr_seq);
                end
            end

            if isfinite(first_ms)
                window_end_ms = first_ms;
            else
                window_end_ms = ptd_ms + 1;   % fallback
            end

            % Extract SINGLE spikes from [0 → window_end_ms]
            spikes_src = sp_single{rec_ch};
            if isempty(spikes_src), continue; end

            t_start_single = t0_single_ms;
            t_end_single   = t0_single_ms + window_end_ms;

            in_win = spikes_src(:,1) >= t_start_single & spikes_src(:,1) < t_end_single;
            spikes_add = spikes_src(in_win,:);
            if isempty(spikes_add), continue; end
            t_inject_start = t0_seq_ms + pulse_offset_ms;

            spikes_add(:,1) = spikes_add(:,1) - t0_single_ms + t_inject_start;

            % Store
            sp_seq{rec_ch} = sortrows([sp_seq{rec_ch}; spikes_add], 1);
            spikes_added_amp = spikes_added_amp + size(spikes_add,1);
        end
    end

    fprintf('Amp %d µA: %d spikes injected.\n', amp, spikes_added_amp);
    total_spikes_added_all = total_spikes_added_all + spikes_added_amp;
end

fprintf('\nTOTAL SPIKES INJECTED: %d\n', total_spikes_added_all);
%% ===================== SAVE RESULT ===================== %%
new_name = strrep(dir('*sp_xia.mat').name, '.sp_xia.mat', '.sp_xia_FirstPulse.mat');
save(fullfile(sequential_folder, new_name), 'sp_seq', '-v7.3');
fprintf('Saved injected sequential file: %s\n', new_name);