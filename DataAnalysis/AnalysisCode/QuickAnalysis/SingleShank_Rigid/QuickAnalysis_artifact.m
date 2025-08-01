%% Quick Analysis for Artifact
clear all
close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Sabrina/Experimental_Design'));

%% Load Intan parameters
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

dName='amplifier';
% info = dir([filepath filesep dName '.dat']);
% info = info.bytes/2;
% nL = (ceil(info / (nChn*FS*double(T)))+1);
vFID = fopen([filepath filesep dName '.dat'],'r'); % read data file

%% Quick Analysis data parameters
nTrials = 500; % number of trials used for analysis

artifact_window_ms = [-1, 5];  % artifact window in ms
artifact_window_samp = round(artifact_window_ms / 1000 * FS);
artifact_samples = diff(artifact_window_samp) + 1;
artifact_time = (artifact_window_samp(1):artifact_window_samp(2)) / FS * 1000;

Artifact_all = zeros(nChn, artifact_samples, nTrials);  % unfiltered


%% Load trigger and parameters
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);
trig = trig(1:nTrials);
nTrig = length(trig); % length of trigger
d = Depth;
% load trial parameters
TrialParams = loadTrialParams;
trialIDs = cell2mat(TrialParams(:,2));
trialIDs = trialIDs(1:nTrials);

%% Load data
for tr = 1:nTrials
    fprintf('Artifact trial %d/%d\n', tr, nTrials);
    stim_index = trig(tr);
    start_index = stim_index + artifact_window_samp(1);
    start_byte = start_index * nChn * 2;
    fseek(vFID, start_byte, 'bof');
    raw_data = fread(vFID, [nChn, artifact_samples], 'int16') * 0.195;

    for ch = 1:nChn
        Artifact_all(ch, :, tr) = raw_data(ch, :);
    end
end

%% Plot Artifact traces (averaged across all trials)
Artifact_avg = mean(Artifact_all, 3);
abs_max = max(abs(Artifact_avg(:)));
y_lim = ceil(abs_max / 10) * 10;

figure('Name', 'Stimulation Artifact (Raw)', 'NumberTitle', 'off', 'Color', 'w');
tl = tiledlayout(8, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Average Raw Artifact around Stim (−1 ms to +5 ms)', 'FontWeight', 'bold');

for row = 1:8
    for col = 1:4
        ch = (col - 1) * 8 + row;
        if ch > nChn
            break;
        end
        ax = nexttile((row - 1) * 4 + col);
        plot(artifact_time, Artifact_avg(ch, :), 'r');
        title(sprintf('Ch %d', ch), 'FontSize', 8);
        xlim([artifact_time(1), artifact_time(end)]);
        ylim([-y_lim, y_lim]);
        xticks(linspace(artifact_time(1), artifact_time(end), 3));
        yticks([-y_lim, 0, y_lim]);
        xlabel('Time (ms)', 'FontSize', 7);
        ylabel('μV', 'FontSize', 7);
    end
end