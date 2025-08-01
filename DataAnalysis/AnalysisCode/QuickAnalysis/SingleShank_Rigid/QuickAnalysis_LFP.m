%% Quick Analysis Code
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
%% Quick Analysis parameters
nTrials = 500; % number of trials used for analysis
window_ms = [-100,100]; % time window for data extract
window_samps = round(window_ms/1000 * FS);
total_sample = diff(window_samps) + 1;
trialLength = diff(window_samps)+1;
blank_ms = 0;%[-1,1]; % artifact blanking window ±1 ms
blank_samps = round(blank_ms/1000 * FS);

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
%% Bandpass filter design for MUA (300-3000 Hz)
lfpFilt = designfilt('lowpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency', 300, ...
    'SampleRate', FS);

% Time vector for plotting
t_axis = (window_samps(1):window_samps(2)) / FS * 1000;

% Preallocate
LFP_all = zeros(nChn, trialLength, nTrials);


for tr = 1:nTrials
    fprintf('Reading trial %d/%d\n', tr, nTrials);
    start_index = trig(tr) + window_samps(1);
    start_byte = start_index * nChn * 2; % 2 bytes per int16
    fseek(vFID, start_byte, 'bof');
    data = fread(vFID, [nChn, total_sample], 'int16') * 0.195; % in uV

    stim_idx = -window_samps(1) + 1;
    blank_range = max(1, stim_idx - blank_samps) : min(total_sample, stim_idx + blank_samps);
    keep_range = setdiff(1:total_sample, blank_range);

    for ch = 1:nChn
        
        thisChn = d(ch);
        raw = data(thisChn, :);

        % Interpolate blank
        raw_interp = raw;
        raw_interp(blank_range) = interp1(keep_range, raw(keep_range), blank_range, 'linear', 'extrap');
        % raw(blank_range(1):blank_range(end)) = interpolate(raw(blank_range(1):blank_range(end)), 1);
        % Filter & rectify
        filtered = filtfilt(lfpFilt, raw_interp);
        % Store
        LFP_all(ch, :, tr) = filtered;

    end
end

fclose(vFID);

 
%% Average LFP across trials
LFP_avg = mean(LFP_all, 3);
y_max = max(LFP_avg(:)) * 1.1;
y_min = min(LFP_avg(:)) * 1.1;

% Set plot layout
nRows = 8;
nCols = 4;
figure('Name', 'Average LFP across all trials', 'NumberTitle', 'off', 'Color', 'w');
tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Average LFP across all trials', 'FontWeight', 'bold');

for row = 1:nRows
    for col = 1:nCols
        ch = (col - 1) * nRows + row;
        if ch > nChn
            break;
        end

        ax = nexttile((row - 1) * nCols + col);
        plot(t_axis, LFP_avg(ch, :), 'b');
        title(sprintf('Ch %d', ch), 'FontSize', 8);
        xlim([t_axis(1), t_axis(end)]);
        ylim([floor(y_min/100)*100, ceil(y_max/100)*100]);

        xticks(linspace(t_axis(1), t_axis(end), 5));
        yticks([round(y_min), round((y_min + y_max)/2), round(y_max)]);
        xlabel('Time (ms)', 'FontSize', 7);
        ylabel('LFP (μV)', 'FontSize', 7);
    end
end


%% Average LFP for each trial
uniqueIDs = unique(trialIDs);
nCond = numel(uniqueIDs);  % number of unique trial types
nRows = 8;
nCols = 4;

for i = 1:nCond
    trialID = uniqueIDs(i);
    trial_idx = find(trialIDs == trialID);  % find trials with this ID

    % Average across trials with this ID
    LFP_cond_avg = mean(LFP_all(:, :, trial_idx), 3);
    y_max = max(LFP_cond_avg(:)) * 1.1;
    y_min = min(LFP_cond_avg(:)) * 1.1;

    % Create figure
    figure('Name', sprintf('LFP - Trial Type %d', trialID), 'NumberTitle', 'off', 'Color', 'w');
    tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Average LFP | Trial ID = %d', trialID), 'FontWeight', 'bold');

    % Plot channels in column-major order
    for row = 1:nRows
        for col = 1:nCols
            ch = (col - 1) * nRows + row;  % column-major indexing
            if ch > nChn
                break;
            end
            nexttile((row - 1) * nCols + col);
            plot(t_axis, LFP_cond_avg(ch, :), 'k');
            title(sprintf('Ch %d', ch));
            xlim([t_axis(1), t_axis(end)]);
            ylim([floor(y_min/100)*100, ceil(y_max/100)*100]);
            % Only show x-axis on bottom row
            if row == nRows; xlabel('Time (ms)');end
            % Only show y-axis on left column
            if col == 1; ylabel('µV'); end
        end
    end
end


