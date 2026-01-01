% ============================================================
%  Extract First Post-Stimulation Spike Time per Trial
%  (Supports multiple PTDs, and now detects simultaneous stim)
%
%  If PTD = 0  → simultaneous stimulation → skip spike search
%
%  Output:
%       firstSpikeTimes{ch}(trial)
%       hasSpike(ch,trial)
%       PTD_ms(trial)
%       isSimultaneous(trial)
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% === User Input ===
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim1';

win_ms = 10;       % window length after PTD (ms)
FS = 30000;        % sampling rate

%% === Load folder ===
cd(data_folder);

%% --- Load filtered spikes ---
sp_file = dir('*sp_xia.mat');
load(sp_file(1).name, 'sp_clipped');

%% --- Load triggers ---
trig = loadTrig(0);

%% --- Load StimParams ---
param_file = dir('*_exp_datafile_*.mat');
S = load(param_file(1).name, 'StimParams','simultaneous_stim','n_Trials');

StimParams = S.StimParams;
simN       = S.simultaneous_stim;
n_Trials   = S.n_Trials;

%% --- Extract PTD for each trial ---
if simN > 1
    % PTD is the 3rd row of each SIM block
    PTD_us = cell2mat(StimParams(3:simN:end, 6));
else
    PTD_us = zeros(n_Trials,1);
end

PTD_ms = PTD_us / 1000;

%% --- Identify Simultaneous-Stimulation Trials ---
% PTD = 0 → simultaneous stimulation
isSimultaneous = (PTD_ms == 0);

fprintf('\nDetected PTD values (ms): '); disp(unique(PTD_ms)');

%% --- Initialize outputs ---
nChn = numel(sp_clipped);
firstSpikeTimes = cell(nChn,1);
hasSpike = zeros(nChn, n_Trials);

fprintf('\nExtracting First Spike Times\n');

%% ===== MAIN LOOP =====
for ch = 1:nChn

    S_ch = sp_clipped{ch};
    firstSpikeTimes{ch} = nan(n_Trials,1);

    if isempty(S_ch)
        continue;
    end

    spike_times_ms = S_ch(:,1);

    for tr = 1:n_Trials

        %% --- Skip simultaneous trials ---
        if isSimultaneous(tr)
            % Do NOT detect first spike for simultaneous stimulation
            continue;
        end

        %% --- Time window based on PTD ---
        t0 = trig(tr)/FS * 1000;       % trigger time ms
        PTD_this = PTD_ms(tr);
        win_start = PTD_this;
        win_end   = PTD_this + win_ms;

        rel_t = spike_times_ms - t0;

        spk_in_win = rel_t(rel_t >= win_start & rel_t <= win_end);

        if ~isempty(spk_in_win)
            firstSpikeTimes{ch}(tr) = spk_in_win(1);
            hasSpike(ch,tr) = 1;
        end
    end

    fprintf('Ch %2d: first spikes found in %d/%d sequential trials.\n', ...
        ch, sum(hasSpike(ch,~isSimultaneous)), sum(~isSimultaneous));
end

fprintf('==============================================================\n');

%% === Save ===
save_name = sprintf('%s_FirstSpikeTimes.mat', erase(sp_file(1).name, '.mat'));
save(save_name, ...
     'firstSpikeTimes', 'PTD_ms', 'isSimultaneous', ...
     'win_ms', 'hasSpike', 'FS','n_Trials');

fprintf('Saved to: %s\n', fullfile(data_folder, save_name));