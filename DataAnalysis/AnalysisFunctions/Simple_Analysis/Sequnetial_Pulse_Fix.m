clear all


%% === USER-DEFINED PARAMETERS === %%
single_folder     = '/Volumes/MACData/Data/Data_Xia/CJ268/XIa_Exp1_Single8_001_251023_030306';        % folder with single stim data
sequential_folder = '/Volumes/MACData/Data/Data_Xia/CJ268/Xia_Exp1_Seq8_5ms_001_251023_032314';    % folder with sequential stim data
single_win_ms     = [2 5];     % spike window after trigger (in ms)
pulse_offset_ms   = 2;       % when to inject relative to 1st pulse (in ms)
FS = 30000;                    % Sampling rate

%% === LOAD SINGLE STIMULATION DATA === %%
cd(single_folder);

% Load sp_clipped
sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No *_sp_xia.mat found in single folder.');
load(sp_file(1).name, 'sp_clipped');  % variable: sp_clipped
sp_single = sp_clipped;
nChn = numel(sp_single);

% Load trig
trig_single = loadTrig(0); 

% Load StimParams and decode amp, set
S = load(dir('*_exp_datafile_*.mat').name, ...
         'StimParams', 'simultaneous_stim', 'CHN', 'E_MAP', 'n_Trials','n_REP');
StimParams = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP = S.E_MAP;
n_Trials_single = S.n_Trials;
n_REP_single = S.n_REP;

% Trial amplitudes and stim sets
trialAmps_single = cell2mat(StimParams(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim:end);
stimNames_single = StimParams(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP(2:end));
stimChPerTrial_single = arrayfun(@(t) unique(idx_all_single((t-1)*simultaneous_stim+(1:simultaneous_stim))), ...
                            (1:n_Trials_single)', 'UniformOutput', false);

%% === LOAD SEQUENTIAL STIMULATION DATA === %%
cd(sequential_folder);

% Load sp_clipped
sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No *_sp_xia.mat found in sequential folder.');
load(sp_file(1).name, 'sp_clipped');  % variable reused: sp_clipped
sp_seq = sp_clipped;
% Load trig
trig_seq =loadTrig(0); 

% Load StimParams and decode
S = load(dir('*_exp_datafile_*.mat').name, ...
         'StimParams', 'simultaneous_stim', 'CHN', 'E_MAP', 'n_Trials','n_REP');
StimParams = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP = S.E_MAP;
n_Trials_seq = S.n_Trials;
CHN = S.CHN;
n_REP_seq = S.n_REP;


% Trial amps and stim sets
trialAmps_seq = cell2mat(StimParams(2:end,16));
trialAmps_seq = trialAmps_seq(1:simultaneous_stim:end);
stimNames_seq = StimParams(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));
stimChPerTrial_seq = arrayfun(@(t) unique(idx_all_seq((t-1)*simultaneous_stim+(1:simultaneous_stim))), ...
                          (1:n_Trials_seq)', 'UniformOutput', false);


d = Depth_s(1); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex
%% === SPIKE INJECTION PROCESS === %%
fprintf('\n=== Counting Added Spikes to Sequential Data ===\n');
total_spikes_added_all = 0;

unique_amps = unique(trialAmps_seq);
for a = 1:numel(unique_amps)
    amp = unique_amps(a);
    total_spikes_added = 0;

    % --- Find all trials with this amplitude --- %
    trials_seq_this_amp    = find(trialAmps_seq    == amp);
    trials_single_this_amp = find(trialAmps_single == amp);

    % === Remove overlapping stim channels === %
    keep_idx = false(size(trials_single_this_amp));
    for k = 1:numel(trials_single_this_amp)
        stimSet_single = stimChPerTrial_single{trials_single_this_amp(k)};
        overlap = false;
        for j = 1:numel(trials_seq_this_amp)
            stimSet_seq = stimChPerTrial_seq{trials_seq_this_amp(j)};
            if isempty(intersect(stimSet_single, stimSet_seq))
                overlap = true;
                break;
            end
        end
        keep_idx(k) = overlap;
    end
    trials_single_this_amp = trials_single_this_amp(keep_idx);

    if numel(trials_seq_this_amp) ~= n_REP_seq
        warning('Unexpected number of seq trials for amp %d', amp);
    end
    if numel(trials_single_this_amp) ~= n_REP_single
        warning('Unexpected number of single trials for amp %d', amp);
    end

    % --- Extract trigger times --- %
    t0_seq_all    = trig_seq(trials_seq_this_amp)    / FS * 1000;
    t0_single_all = trig_single(trials_single_this_amp) / FS * 1000;

    % === Loop over trials === %
    for ti = 1:n_REP_seq
        tr_seq    = trials_seq_this_amp(ti);
        tr_single = trials_single_this_amp(ti);

        stimSet_seq    = stimChPerTrial_seq{tr_seq};
        stimSet_single = stimChPerTrial_single{tr_single};

        % Get injection channel (non-overlapping)
        ch_to_inject = setdiff(stimSet_single, stimSet_seq);
        if isempty(ch_to_inject), continue; end
        if numel(ch_to_inject) > 1
            warning('Multiple extra stim channels in single trial %d. Using first.', tr_single);
            ch_to_inject = ch_to_inject(1);
        end

        % === Time windows === %
        t0_seq    = t0_seq_all(ti);
        t0_single = t0_single_all(ti);
        t_extract = t0_single + single_win_ms;
        t_inject  = t0_seq + single_win_ms + pulse_offset_ms;

        % === Loop over channel(s) === %
        for ch = 1:nChn
            spikes = sp_single{d(ch)};
            if isempty(spikes), continue; end

            in_win = spikes(:,1) >= t_extract(1) & spikes(:,1) < t_extract(2);
            spikes_add = spikes(in_win,:);
            if isempty(spikes_add), continue; end

            % Shift spike timing
            spikes_add(:,1) = spikes_add(:,1) - t0_single + t0_seq + pulse_offset_ms;

            % Inject and count
            sp_seq{d(ch)} = sortrows([sp_seq{d(ch)}; spikes_add], 1);
            total_spikes_added = total_spikes_added + size(spikes_add,1);
        end
    end

    fprintf('Amplitude %d µA: %d spikes added\n', amp, total_spikes_added);
    total_spikes_added_all = total_spikes_added_all + total_spikes_added;
end

fprintf('=== TOTAL SPIKES ADDED TO SEQ DATA: %d ===\n', total_spikes_added_all);

fprintf(' Injection complete! Saving to new file...\n');

%% === SAVE UPDATED sp_clipped FILE === %%
base_seq = sp_file(1).name;
new_name = strrep(base_seq, '.sp_xia.mat', '.sp_xia_FirstPulse.mat');
save(fullfile(sequential_folder, new_name), 'sp_seq', '-v7.3');

