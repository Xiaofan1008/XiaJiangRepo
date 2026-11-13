%% === SPIKE INJECTION: use sequential first-spike per (channel, trial);
%                       inject ONLY single-pulse spikes from the first stim electrode ===
clear all;

%% === USER PARAMETERS === %%
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Single8';
sequential_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq8';
pulse_offset_ms   = 0;        % shift injected spikes relative to the 1st sequential pulse (ms)
fallback_win_end  = 5;        % used if first-spike time is missing/invalid (ms)
use_fallback      = true;    % set true to use fallback; false to skip those trials/channels
FS = 30000;

%% === LOAD SINGLE STIM DATA (SOURCE) === %%
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

% single: per-trial amplitude & stim-channel mapping
trialAmps_single = cell2mat(StimParams_single(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);
stimNames_single  = StimParams_single(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP(2:end));
stimChPerTrial_single = arrayfun(@(t) unique(idx_all_single((t-1)*simultaneous_stim1 + (1:simultaneous_stim1))), ...
    (1:n_Trials_single)', 'UniformOutput', false);

%% === LOAD SEQUENTIAL STIM DATA (TARGET) === %%
cd(sequential_folder);
load(dir('*sp_xia.mat').name, 'sp_clipped'); 
sp_seq  = sp_clipped;
trig_seq = loadTrig(0);

S2 = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials','n_REP');
StimParams_seq  = S2.StimParams;
n_Trials_seq    = S2.n_Trials;
simultaneous_stim2 = S2.simultaneous_stim;

% sequential: per-trial amplitude & stim-channel mapping
trialAmps_seq   = cell2mat(StimParams_seq(2:end,16));
trialAmps_seq   = trialAmps_seq(1:simultaneous_stim2:end);
stimNames_seq   = StimParams_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));
stimChPerTrial_seq = arrayfun(@(t) unique(idx_all_seq((t-1)*simultaneous_stim2 + (1:simultaneous_stim2))), ...
    (1:n_Trials_seq)', 'UniformOutput', false);

% optional: POST-TRIGGER DELAY of 2nd pulse present in column 6 (not used for 1st pulse injection)
% post_trigger_delay_us = cell2mat(StimParams_seq(3:simultaneous_stim2:end,6));

% channel index mapping
d = Depth_s(1);

%% === LOAD FIRST SPIKE TIMES FROM THE SEQUENTIAL DATA === %%
% Expect firstSpikeTimes to be a 1×nChn cell: firstSpikeTimes{ch}(trial) in ms (relative to seq trigger)
fst_file = dir('*FirstSpikeTimes.mat');
assert(~isempty(fst_file), 'No *FirstSpikeTimes.mat found in the sequential folder.');
load(fst_file(1).name, 'firstSpikeTimes');
fprintf('Loaded sequential first-spike times: %s\n', fst_file(1).name);

%% === MAIN INJECTION LOOP === %%
unique_amps = unique(trialAmps_seq);
total_spikes_added_all = 0;

fprintf('\n=== Injecting single-pulse spikes using seq first-spike per (ch,trial); single source = FIRST stim electrode ===\n');

for a = 1:numel(unique_amps)
    amp = unique_amps(a);
    trials_seq_this_amp = find(trialAmps_seq == amp);
    total_spikes_added_amp = 0;

    for ti = 1:numel(trials_seq_this_amp)
        tr_seq = trials_seq_this_amp(ti);

        % --- Determine first stim electrode for this sequential trial (e.g., [24,28] → 24) ---
        stimSet_seq = stimChPerTrial_seq{tr_seq};
        if isempty(stimSet_seq), continue; end
        ch_first_stim = stimSet_seq(1);  % ALWAYS use the first electrode of the sequential pair

        % --- Find matching single trials where ch_first_stim was stimulated at the same amplitude ---
        match_single = find( cellfun(@(x) ismember(ch_first_stim, x), stimChPerTrial_single) & ...
                             (trialAmps_single == amp) );
        if isempty(match_single)
            warning('No single trial found for stim ch %d at %duA (seq trial %d).', ch_first_stim, amp, tr_seq);
            continue;
        end

        % safe cyclic reuse
        tr_single = match_single( mod(ti-1, numel(match_single)) + 1 );

        % anchors (ms)
        t0_seq_ms    = trig_seq(tr_seq)   / FS * 1000;
        t0_single_ms = trig_single(tr_single) / FS * 1000;

        % For each RECORDING channel, define its OWN window end from sequential first spike time
        for rec_ch = 1:nChn
            % get channel-specific first-spike time in SEQ data for THIS trial
            first_ms = NaN;
            if rec_ch <= numel(firstSpikeTimes) && tr_seq <= numel(firstSpikeTimes{rec_ch})
                first_ms = firstSpikeTimes{rec_ch}(tr_seq);
            end

            if ~(isfinite(first_ms) && first_ms > 1)
                if use_fallback
                    win_end = fallback_win_end;   % use small default window
                else
                    continue;                      % skip this channel/trial
                end
            else
                win_end = first_ms;                % dynamic end
            end

            % extract window in SINGLE data (relative to t0_single)
            % t_start = t0_single_ms + 1;        % 1 ms
            t_start = t0_single_ms;        % 1 ms
            t_end   = t0_single_ms + win_end;  % channel-specific
            % spikes_src = sp_single{d(rec_ch)};
            spikes_src = sp_single{rec_ch};
            if isempty(spikes_src), continue; end

            in_win = spikes_src(:,1) >= t_start & spikes_src(:,1) < t_end;
            spikes_add = spikes_src(in_win,:);
            if isempty(spikes_add), continue; end

            % shift to SEQ timeline; inject at first pulse time (no PTD added)
            t_inject = t0_seq_ms + pulse_offset_ms;
            spikes_add(:,1) = spikes_add(:,1) - t0_single_ms + t_inject;

            % append and keep time-ordered
            sp_seq{d(rec_ch)} = sortrows([sp_seq{d(rec_ch)}; spikes_add], 1);
            total_spikes_added_amp = total_spikes_added_amp + size(spikes_add,1);
        end
    end

    fprintf('Amplitude %d µA: %d spikes added.\n', amp, total_spikes_added_amp);
    total_spikes_added_all = total_spikes_added_all + total_spikes_added_amp;
end

fprintf('\nTOTAL SPIKES ADDED TO SEQ DATA: %d\n', total_spikes_added_all);

%% === SAVE UPDATED SEQ FILE === %%
new_name = strrep(dir('*sp_xia.mat').name, '.sp_xia.mat', '.sp_xia_FirstPulse.mat');
save(fullfile(sequential_folder, new_name), 'sp_seq', '-v7.3');
fprintf('Saved updated sequential spikes as: %s\n', new_name);