% ============================================================
%  Extract First Post-Stimulation Spike Time per Trial
% 
%  This script loads filtered spike data (.sp_xia.mat) and the experiment
%  parameters to find, for each channel and each trial, the *first spike*
%  that occurs after stimulation onset.
%
%  The result (firstSpikeTimes) is saved into a .mat file.
% ============================================================

clear;

%% === User Input ===
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq8_251106_190122';
response_window_ms = [5, 15];  % Time window after stim to search for first spike (ms)
FS = 30000;  % Sampling rate (Hz)

%% === Check folder and load files ===
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Working in: %s\n', data_folder);

% --- Load spike data (filtered result) ---
sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No .sp_xia.mat file found.');
load(sp_file(1).name, 'sp_clipped');
fprintf('Loaded filtered spike file: %s\n', sp_file(1).name);

% --- Load trigger file ---
if isempty(dir('*.trig.dat'))
    error('No trigger file found in this folder.');
end
trig = loadTrig(0);  % function from MASSIVE toolbox
fprintf('Loaded trigger data.\n');

% --- Load StimParams (for trial count etc.) ---
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat file found.');
S = load(param_file(1).name, 'StimParams', 'simultaneous_stim', 'n_Trials');
StimParams = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
n_Trials = S.n_Trials;
fprintf('Loaded StimParams: %d total trials.\n', n_Trials);

%% === Initialize result storage ===
nChn = numel(sp_clipped);
firstSpikeTimes = cell(nChn, 1);  % each cell: 1 × n_Trials vector (ms)
hasSpike = zeros(nChn, n_Trials); % 1 if at least one spike in window

fprintf('\n===== Extracting First Spike Times =====\n');

%% === Loop through channels and trials ===
for ch = 1:nChn
    S_ch = sp_clipped{ch};
    if isempty(S_ch)
        firstSpikeTimes{ch} = nan(n_Trials, 1);
        continue;
    end

    spike_timestamps = S_ch(:,1); % spike times in samples
    firstSpikeTimes{ch} = nan(n_Trials, 1);

    for tr = 1:n_Trials
        t0 = trig(tr)/FS*1000;  % convert to ms
        rel_spk = spike_timestamps - t0; % spike times relative to trigger

        % find spikes within post-stimulation window
        post_spk = rel_spk(rel_spk >= response_window_ms(1) & rel_spk <= response_window_ms(2));

        if ~isempty(post_spk)
            firstSpikeTimes{ch}(tr) = post_spk(1);
            hasSpike(ch, tr) = 1;
        end
    end

    fprintf('Ch %2d: found first spikes in %d / %d trials.\n', ...
        ch, sum(hasSpike(ch,:)), n_Trials);
end

fprintf('=======================================\n\n');

%% === Save results ===
save_name = sprintf('%s_FirstSpikeTimes.mat', erase(sp_file(1).name, '.mat'));
save(save_name, 'firstSpikeTimes', 'response_window_ms', 'hasSpike', 'FS', 'n_Trials');
fprintf('Saved results to: %s\n', fullfile(data_folder, save_name));

%% === Optional summary ===
total_detected = sum(hasSpike(:));
fprintf('Total detected first spikes: %d across all channels and trials.\n', total_detected);