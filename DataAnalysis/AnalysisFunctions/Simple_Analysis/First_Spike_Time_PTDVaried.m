% ============================================================
%  Extract First Post-Stimulation Spike Time per Trial
%
%  NOW SUPPORTS MULTIPLE PTD VALUES IN THE SAME FOLDER
%
%  For each trial:
%      PTD_ms = StimParams column 6 (µs → ms)
%      First-spike window = [PTD_ms, PTD_ms + win_ms]
%
%  Saves: firstSpikeTimes, PTD_ms, hasSpike
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% === User Input ===
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim4_251128_150648';

win_ms = 10;       % window length after PTD (ms)
FS = 30000;        % sampling rate

%% === Load folder ===
if ~isfolder(data_folder)
    error('The specified folder does not exist.');
end
cd(data_folder);
fprintf('Working in: %s\n', data_folder);

%% --- Load filtered spikes ---
sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No .sp_xia.mat file found.');
load(sp_file(1).name, 'sp_clipped');
fprintf('Loaded filtered spikes: %s\n', sp_file(1).name);

%% --- Load triggers ---
if isempty(dir('*.trig.dat'))
    error('No trigger file found.');
end
trig = loadTrig(0);
fprintf('Loaded trigger file.\n');

%% --- Load StimParams ---
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat found.');
S = load(param_file(1).name, ...
         'StimParams','simultaneous_stim','n_Trials');

StimParams = S.StimParams;
simN       = S.simultaneous_stim;
n_Trials   = S.n_Trials;

fprintf('Loaded StimParams: %d trials.\n', n_Trials);

%%       EXTRACT PTD PER TRIAL

if simN > 1
    % PTD for each trial is row 2 of every block of simN
    PTD_us = cell2mat(StimParams(3:simN:end, 6));   % <-- IMPORTANT
else
    PTD_us = zeros(n_Trials,1);
end

PTD_ms = PTD_us / 1000;
fprintf('Detected PTD values (ms): '); disp(unique(PTD_ms)');

%%        FIND FIRST SPIKE PER TRIAL BASED ON THAT PTD

nChn = numel(sp_clipped);
firstSpikeTimes = cell(nChn,1);
hasSpike = zeros(nChn, n_Trials);

fprintf('\n===== Extracting First Spike Times (PTD-based windows) =====\n');

for ch = 1:nChn
    S_ch = sp_clipped{ch};
    if isempty(S_ch)
        firstSpikeTimes{ch} = nan(n_Trials,1);
        continue;
    end

    spike_times_ms = S_ch(:,1);   % already in ms
    firstSpikeTimes{ch} = nan(n_Trials,1);

    for tr = 1:n_Trials
        t0 = trig(tr)/FS * 1000;   % trigger time in ms
        
        PTD_this = PTD_ms(tr);     % <-- trial-specific PTD
        win_start = PTD_this;
        win_end   = PTD_this + win_ms;

        rel_spk = spike_times_ms - t0;
        % window: [PTD, PTD+win_ms]
        spk_in_win = rel_spk(rel_spk >= win_start & rel_spk <= win_end);

        if ~isempty(spk_in_win)
            firstSpikeTimes{ch}(tr) = spk_in_win(1);
            hasSpike(ch,tr) = 1;
        end
    end

    fprintf('Ch %2d: first spikes found in %d/%d trials.\n', ...
        ch, sum(hasSpike(ch,:)), n_Trials);
end

fprintf('==============================================================\n');

%% === Save ===
save_name = sprintf('%s_FirstSpikeTimes.mat', erase(sp_file(1).name, '.mat'));
save(save_name, 'firstSpikeTimes', 'PTD_ms', 'win_ms', 'hasSpike', ...
                'FS','n_Trials');

fprintf('Saved to: %s\n', fullfile(data_folder, save_name));