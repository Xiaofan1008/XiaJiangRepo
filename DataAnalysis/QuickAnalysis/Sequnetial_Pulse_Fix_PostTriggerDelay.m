%% === SPIKE INJECTION PROCESS WITH POST-TRIGGER DELAY SUPPORT === %%
clear all

% === USER PARAMETERS === %
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Single6_251106_162313';
sequential_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq6_5ms_251106_164125';
single_win_ms     = [1 6.5];     % spike window after trigger (ms)
pulse_offset_ms   = 0;         % when to inject relative to 1st pulse (ms)
FS = 30000;                    % Sampling rate

%% === LOAD SINGLE STIMULATION DATA === %%
cd(single_folder);
load(dir('*sp_xia.mat').name, 'sp_clipped'); 
sp_single = sp_clipped; 
nChn = numel(sp_single);
trig_single = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams', 'simultaneous_stim', 'E_MAP', 'n_Trials', 'n_REP');
StimParams = S.StimParams; 
E_MAP = S.E_MAP;
n_Trials_single = S.n_Trials; 
n_REP_single = S.n_REP;
simultaneous_stim = S.simultaneous_stim;

trialAmps_single = cell2mat(StimParams(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim:end);
stimNames_single = StimParams(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP(2:end));
stimChPerTrial_single = arrayfun(@(t) unique(idx_all_single((t-1)*simultaneous_stim + (1:simultaneous_stim))), ...
    (1:n_Trials_single)', 'UniformOutput', false);

%% === LOAD SEQUENTIAL STIMULATION DATA === %%
cd(sequential_folder);
load(dir('*sp_xia.mat').name, 'sp_clipped'); 
sp_seq = sp_clipped;
trig_seq = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams', 'simultaneous_stim', 'E_MAP', 'n_Trials', 'n_REP');
StimParams = S.StimParams; 
E_MAP = S.E_MAP;
n_Trials_seq = S.n_Trials; 
n_REP_seq = S.n_REP;
simultaneous_stim = S.simultaneous_stim;

trialAmps_seq = cell2mat(StimParams(2:end,16));
trialAmps_seq = trialAmps_seq(1:simultaneous_stim:end);
stimNames_seq = StimParams(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));
stimChPerTrial_seq = arrayfun(@(t) unique(idx_all_seq((t-1)*simultaneous_stim + (1:simultaneous_stim))), ...
    (1:n_Trials_seq)', 'UniformOutput', false);

% Load POST-TRIGGER DELAY (column 6)
post_trigger_delay_us = cell2mat(StimParams(3:simultaneous_stim:end,6));

d = Depth_s(1);
total_spikes_added_all = 0;
unique_amps = unique(trialAmps_seq);

fprintf('\n=== Counting Added Spikes to Sequential Data ===\n');

for a = 1:numel(unique_amps)
    amp = unique_amps(a);
    total_spikes_added = 0;

    trials_seq_this_amp = find(trialAmps_seq == amp);
    t0_seq_all = trig_seq(trials_seq_this_amp) / FS * 1000;

    for ti = 1:numel(trials_seq_this_amp)
        tr_seq = trials_seq_this_amp(ti);
        stimSet_seq = stimChPerTrial_seq{tr_seq};
        ch_first = stimSet_seq(1);  % Always use first channel in seq set

        % --- Find matching single-stim trials for this channel --- %
        matching_single_trials = find(cellfun(@(x) ismember(ch_first, x), stimChPerTrial_single) ...
                                      & trialAmps_single == amp);

        if isempty(matching_single_trials)
            warning('No single-trial found for channel %d at %duA', ch_first, amp);
            continue;
        end

        % --- Safe cyclic reuse of single trials --- %
        tr_single_idx = mod(ti-1, length(matching_single_trials)) + 1;
        tr_single = matching_single_trials(tr_single_idx);

        % Trigger times
        t0_seq = t0_seq_all(ti);
        t0_single = trig_single(tr_single) / FS * 1000;

        % Compute injection time (account for post-trigger delay)
        PTD_ms = post_trigger_delay_us(tr_seq) / 1000;
        % t_inject = t0_seq + PTD_ms + pulse_offset_ms;
        t_inject = t0_seq + pulse_offset_ms;
        t_extract = t0_single + single_win_ms;

        % Inject spikes from single data
        for ch = 1:nChn
            spikes = sp_single{d(ch)};
            if isempty(spikes), continue; end

            in_win = spikes(:,1) >= t_extract(1) & spikes(:,1) < t_extract(2);
            spikes_add = spikes(in_win,:);
            if isempty(spikes_add), continue; end

            spikes_add(:,1) = spikes_add(:,1) - t0_single + t_inject;
            sp_seq{d(ch)} = sortrows([sp_seq{d(ch)}; spikes_add], 1);
            total_spikes_added = total_spikes_added + size(spikes_add,1);
        end
    end

    fprintf('Amplitude %d uA: %d spikes added\n', amp, total_spikes_added);
    total_spikes_added_all = total_spikes_added_all + total_spikes_added;
end

fprintf('=== TOTAL SPIKES ADDED TO SEQ DATA: %d ===\n', total_spikes_added_all);

% === SAVE NEW FILE === %
new_name = strrep(dir('*sp_xia.mat').name, '.sp_xia.mat', '.sp_xia_FirstPulse.mat');
save(fullfile(sequential_folder, new_name), 'sp_seq', '-v7.3');
fprintf('Injection complete. Saved as: %s\n', new_name);