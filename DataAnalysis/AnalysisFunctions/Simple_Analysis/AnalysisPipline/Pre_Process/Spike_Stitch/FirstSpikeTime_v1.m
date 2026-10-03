%% ============================================================
% FirstSpikeTime_Sequential_v2.m
% Step 1 of the sequential-stimulation artifact recovery pipeline.
%
% For each SEQUENTIAL trial only (PTD > 0), find the time of the first
% spike occurring after the second pulse -- i.e. where real,
% uncontaminated recording resumes after the second pulse's artifact.
% This is used by the NEXT script to decide how far into each sequential
% trial the "blank and substitute with single-pulse data" window should
% extend (falling back to a fixed duration when no spike is detected).
%
% This script is READ-ONLY with respect to the spike data -- it never
% modifies sp_clipped, just inspects it. SIMULTANEOUS trials (PTD == 0)
% are explicitly skipped in the per-trial loop below and left as NaN:
% nothing is computed, estimated, or written for them, since they never
% need this recovery step at all.
%
% Reads the FILTERED spike-times file (<base>.sp_xia_filtered.mat,
% produced by SpikeFiltering_Evoked_Cleanup.m) directly -- entered below,
% not guessed from the folder name. The trigger file and the stim-
% protocol file are still located via the data folder, since their
% naming is fixed/unambiguous (one *.trig.dat set, one
% *_exp_datafile_*.mat) -- this isn't the kind of guessing the pipeline
% moved away from earlier (that was about multiple possible spike
% variable names/suffixes, not these).
% ============================================================
clear all; clc;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1';
                    % folder holding this dataset's trigger
                    % (*.trig.dat) and stim-protocol
                    % (*_exp_datafile_*.mat) files.

times_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_filtered.mat';
                    % the FILTERED spike-times file for this SAME
                    % dataset (sp_clipped) -- entered directly.

win_ms = 10;        % how far past PTD to search for the first real spike (ms)
FS     = 30000;     % sampling rate (overridden below if saved in the file)

%% ================= LOAD DATA =================
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

assert(isfile(times_file), 'Cannot find %s.', times_file);
S_t = load(times_file);
assert(isfield(S_t, 'sp_clipped'), 'Variable "sp_clipped" not found in %s.', times_file);
sp_clipped = S_t.sp_clipped;
if isfield(S_t, 'fs'), FS = S_t.fs; end

if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);

param_file = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(param_file), 'No *_exp_datafile_*.mat found in %s.', data_folder);
S = load(fullfile(data_folder, param_file(1).name), 'StimParams', 'simultaneous_stim', 'n_Trials');
StimParams = S.StimParams;
simN       = S.simultaneous_stim;
n_Trials   = S.n_Trials;

%% ================= PTD PER TRIAL =================
if simN > 1
    PTD_us = cell2mat(StimParams(3:simN:end, 6));   % PTD stored on the 2nd-pulse row of each block
else
    PTD_us = zeros(n_Trials,1);
end
PTD_ms = PTD_us / 1000;

% Simultaneous trials (PTD == 0): never processed below -- left as NaN.
isSimultaneous = (PTD_ms == 0);
fprintf('\nDetected PTD values (ms): '); disp(unique(PTD_ms)');
fprintf('%d/%d trials are simultaneous (PTD=0) and will be left untouched.\n', ...
    sum(isSimultaneous), n_Trials);

%% ================= FIRST SPIKE TIME PER SEQUENTIAL TRIAL =================
nChn = numel(sp_clipped);
firstSpikeTimes = cell(nChn,1);
hasSpike = zeros(nChn, n_Trials);

fprintf('\nExtracting first post-2nd-pulse spike time (sequential trials only)...\n');
for ch = 1:nChn
    firstSpikeTimes{ch} = nan(n_Trials,1);
    S_ch = sp_clipped{ch};
    if isempty(S_ch), continue; end
    spike_times_ms = S_ch(:,1);

    for tr = 1:n_Trials
        if isSimultaneous(tr)
            continue;   % simultaneous trial -- do not touch or estimate anything
        end

        t0 = trig(tr) / FS * 1000;
        win_start = PTD_ms(tr);
        win_end   = PTD_ms(tr) + win_ms;

        rel_t = spike_times_ms - t0;
        spk_in_win = rel_t(rel_t >= win_start & rel_t <= win_end);

        if ~isempty(spk_in_win)
            firstSpikeTimes{ch}(tr) = spk_in_win(1);
            hasSpike(ch,tr) = 1;
        end
    end

    fprintf('Ch %2d: first spike found in %d/%d sequential trials.\n', ...
        ch, sum(hasSpike(ch, ~isSimultaneous)), sum(~isSimultaneous));
end
fprintf('==============================================================\n');

%% ================= SAVE =================
% Output name derived directly from the INPUT times_file, not the folder
% name, matching the rest of the pipeline's convention.
out_name = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_FirstSpikeTimes.mat');
if strcmp(out_name, times_file)
    [p,n,e] = fileparts(times_file); out_name = fullfile(p, [n '_FirstSpikeTimes' e]);
end
save(out_name, 'firstSpikeTimes', 'hasSpike', 'PTD_ms', 'isSimultaneous', 'win_ms', 'FS', 'n_Trials');
fprintf('\nSaved to: %s\n', out_name);