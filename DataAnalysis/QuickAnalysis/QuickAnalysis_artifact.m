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
vFID = fopen([filepath filesep dName '.dat'],'r'); % read data file

%% ------ Quick Analysis data parameters ------ %%
nTrials = 220; % number of trials used for analysis

artifact_window_ms = [-1, 10];  % artifact window in ms
artifact_window_samp = round(artifact_window_ms / 1000 * FS);
artifact_samples = diff(artifact_window_samp) + 1;
artifact_time = (artifact_window_samp(1):artifact_window_samp(2)) / FS * 1000;

Artifact_all = zeros(nChn, artifact_samples, nTrials);  % unfiltered


%% ------ Load trigger and parameters ------ %%
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

% load StimParameters
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
fileDIR = fileDIR(1).name;
S = load(fileDIR,'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
CHN = S.CHN;
E_MAP = S.E_MAP;
n_Trials = S.n_Trials;

%% ------- Load data ------ %
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

%% ------ Plot Artifact traces (averaged across all trials) ------ %
% Artifact_avg = mean(Artifact_all, 3);
% abs_max = max(abs(Artifact_avg(:)));
% y_lim = ceil(abs_max / 10) * 10;
% 
% figure('Name', 'Stimulation Artifact (Raw)', 'NumberTitle', 'off', 'Color', 'w');
% tl = tiledlayout(8, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% title(tl, 'Average Raw Artifact around Stim', 'FontWeight', 'bold');
% 
% for row = 1:8
%     for col = 1:4
%         ch = (col - 1) * 8 + row;
%         if ch > nChn
%             break;
%         end
%         ax = nexttile((row - 1) * 4 + col);
%         plot(artifact_time, Artifact_avg(ch, :), 'r');
%         title(sprintf('Ch %d', ch), 'FontSize', 8);
%         xlim([artifact_time(1), artifact_time(end)]);
%         ylim([-y_lim, y_lim]);
%         xticks(linspace(artifact_time(1), artifact_time(end), 3));
%         yticks([-y_lim, 0, y_lim]);
%         xlabel('Time (ms)', 'FontSize', 7);
%         ylabel('μV', 'FontSize', 7);
%     end
% end

%% ------ Plot Artifact (across Amplitudes) ------ %%
% Load trial AMP
fileDIR = dir('*_exp_datafile_*.mat');
if isempty(fileDIR)
    return
end
fileDIR = fileDIR.name;
load(fileDIR,'StimParams');
load(fileDIR,'simultaneous_stim');
load(fileDIR,'CHN');

trialAmps = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps(1:simultaneous_stim:(nTrials*simultaneous_stim));

% Group trials by amplitude 
[uniqAmps, ~, ampIdx] = unique(trialAmps(:));
n_AMP = numel(uniqAmps);
Amps = uniqAmps;
if any(Amps == -1)
    Amps(Amps == -1) = 0;
end
% ------ Average per amplitude ------ %
Artifact_avg_amp = zeros(nChn, artifact_samples, n_AMP); 
nPerAmp = zeros(n_AMP,1);
for k = 1: n_AMP
    idx = (ampIdx == k); 
    nPerAmp(k) = nnz(idx);
    if nPerAmp(k) > 0
        Artifact_avg_amp(:,:,k) = mean(Artifact_all(:,:,idx), 3);
    end
end

% ------ Stimulation Channels ------ %
E_NAME = E_MAP(2:end); % Channel name e.g. "A-001"
if exist('StimParams','var') && ~isempty(StimParams)
    % Extract first column for channel 
    stimNames = StimParams(2:end,1);
    [tf, idx_all] = ismember(stimNames, E_NAME);
end
% group channels per trial
stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v  = unique(idx_all(rr));                                  % dedupe within trial
    v  = v(v>=0).';                                             % drop zeros, row vector
    stimChPerTrial_all{t} = v;
end
% find unique channels 
uniqueCh = unique([stimChPerTrial_all{:}]); 
% find unique channel set for each trial
comb = zeros(n_Trials, simultaneous_stim);
for i = 1:n_Trials
    v = stimChPerTrial_all{i};  
    comb(i,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb, 'rows');
combClass_win = combClass(1 : nTrials);

sets_here = unique(combClass_win(:))';
nSets     = numel(sets_here);
% Label for each set and its stim channels
setLabel   = strings(1,nSets);
stimIdxSet = cell(1,nSets);
for s = 1:nSets
    cc = sets_here(s);
    stimIdx = uniqueComb(cc,:); 
    stimIdx = stimIdx(stimIdx>0);
    stimIdxSet{s} = stimIdx;
    % setLabel(s) = strjoin(E_NAME(stimIdx), ' + ');
    setLabel(s) = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx,'UniformOutput', false), ' + ');
end

% ------ Average artifact per (Amplitude × Set) ------ %
Artifact_avg_amp_set = zeros(nChn, artifact_samples, n_AMP, nSets);
nPerAmpSet           = zeros(n_AMP, nSets);

for s = 1:nSets
    cc = sets_here(s);
    tr_in_set = find(combClass_win == cc);    % trials belonging to this set

    for k = 1:n_AMP
        tr_sel = tr_in_set(ampIdx(tr_in_set) == k);  % trials in this set at this amp
        nPerAmpSet(k,s) = numel(tr_sel);
        if nPerAmpSet(k,s) > 0
            Artifact_avg_amp_set(:,:,k,s) = mean(Artifact_all(:,:,tr_sel), 3);
        end
    end
end


%% -------- plot -------- %%

% ------ (1): Average across amp ------ %
% y-axis 
abs_max = max(abs(Artifact_avg_amp), [], 'all');
y_lim   = ceil(abs_max / 1000) * 1000;

cmap = turbo(n_AMP);

% 1) Subplot function 
figure('Name', 'Stimulation Artifact (Raw) by Amplitude', ...
       'NumberTitle', 'off', 'Color', 'w');
tl = tiledlayout(8, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Artifact Waveforms across Stim Amplitudes'), ...
      'FontWeight', 'bold');
        
for ch = 1:nChn
        % subplot(4,8,ch)
        subplot(8,8,ch)
        hold on
        % Soft red background if this is a stimulation channel
        if ismember(ch, uniqueCh)
            set(gca, 'Color', [1 0.95 0.95]);    % light pink
            baseLW = 1.8;                       % slightly thicker for stim chans
        else
            baseLW = 1.2;
        end

        % overlay one curve per amplitude
        for k = 1:n_AMP
            if nPerAmp(k) > 0
                h(k) = plot(artifact_time, squeeze(Artifact_avg_amp(ch,:,k)), ...
                     'Color', cmap(k,:), 'LineWidth', 2);
            end
        end

        xline(0, 'r-');  % stim time
        title(sprintf('Ch %d', ch), 'FontSize', 8);
        xlim([artifact_time(1), artifact_time(end)]);
        ylim([-y_lim, y_lim]);
        xticks(linspace(artifact_time(1), artifact_time(end), 5));
        yticks([-y_lim, 0, y_lim]);
        xlabel('Time (ms)', 'FontSize', 8);
        ylabel('μV', 'FontSize', 8);
        box('off');

        hold off
end

% Single legend for amplitudes (placed below the grid)
legLabels = arrayfun(@(a) sprintf('%g \\muA', a), Amps, 'UniformOutput', false);
lgd = legend(h, legLabels, 'Orientation','horizontal', 'Box','off','FontSize',10);
lgd.Position = [0.1, 0.01, 0.8, 0.05];   % centered below all subplots


% 2) Handset function 
figure('Name', 'Stimulation Artifact (Raw) by Amplitude', ...
       'NumberTitle', 'off', 'Color', 'w');
tl = tiledlayout(8, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Artifact Waveforms across Stim Amplitudes'), ...
      'FontWeight', 'bold');

for row = 1:8
    for col = 1:4
        ch = (col - 1) * 8 + row;
        if ch > nChn, break; end
        ax = nexttile((row - 1) * 4 + col); hold(ax,'on');
        
        ax = nexttile(ch); hold(ax,'on');
        % Soft red background if this is a stimulation channel
        if ismember(ch, uniqueCh)
            set(ax, 'Color', [1 0.95 0.95]);    % light pink
            baseLW = 1.8;                       % slightly thicker for stim chans
        else
            baseLW = 1.2;
        end

        % overlay one curve per amplitude
        for k = 1:n_AMP
            if nPerAmp(k) > 0
                plot(ax, artifact_time, squeeze(Artifact_avg_amp(ch,:,k)), ...
                     'Color', cmap(k,:), 'LineWidth', 2);
            end
        end
        xline(ax, 0, 'r-');  % stim time
        title(ax, sprintf('Ch %d', ch), 'FontSize', 8);
        xlim(ax, [artifact_time(1), artifact_time(end)]);
        ylim(ax, [-y_lim, y_lim]);
        xticks(ax, linspace(artifact_time(1), artifact_time(end), 5));
        yticks(ax, [-y_lim, 0, y_lim]);
        xlabel(ax, 'Time (ms)', 'FontSize', 8);
        ylabel(ax, 'μV', 'FontSize', 8);
        box(ax,'off');

     end
end

% Single legend for amplitudes (placed below the grid)
lg = legend(compose('%.0f µA', Amps), ...
            'NumColumns', min(n_AMP, 6),'Orientation', 'horizontal','Location','southoutside','FontSize',10);
lg.Box = 'off';



% ------ (2) Average across Amp, different Stim Set ------ %
% Global y-limit across all sets & amps
abs_max = max(abs(Artifact_avg_amp_set), [], 'all');
y_lim   = max(100, 100 * ceil(abs_max/100));      % nice rounded

cmap = turbo(n_AMP);
legLabels = arrayfun(@(a) sprintf('%g \\muA', a), Amps, 'UniformOutput', false);

for s = 1:nSets
    figName = sprintf('Artifact | Set %d: %s', sets_here(s), setLabel(s));
    figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w');
    tl = tiledlayout(8, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Artifact (raw) | %s', setLabel(s)), 'FontWeight','bold');

    hL = gobjects(1,n_AMP);  % for legend handles once

    for ch = 1:nChn
        ax = nexttile(tl); hold(ax,'on');

        % Light background if this channel is a stim channel for this set
        if ismember(ch, stimIdxSet{s})
            set(ax, 'Color', [1 0.95 0.95]);
        end

        % Overlays by amplitude for THIS set
        for k = 1:n_AMP
            if nPerAmpSet(k,s) > 0
                h = plot(ax, artifact_time, squeeze(Artifact_avg_amp_set(ch,:,k,s)), ...
                         'Color', cmap(k,:), 'LineWidth', 1.8);
                if ~isgraphics(hL(k)), hL(k) = h; end  % keep one handle per amp
            end
        end

        xline(ax, 0, 'r-');
        title(ax, sprintf('Ch %d', ch), 'FontSize', 8);
        xlim(ax, [artifact_time(1), artifact_time(end)]);
        ylim(ax, [-y_lim, y_lim]);
        xticks(ax, linspace(artifact_time(1), artifact_time(end), 5));
        yticks(ax, [-y_lim, 0, y_lim]);
        xlabel(ax, 'Time (ms)', 'FontSize', 8);
        ylabel(ax, '\muV', 'FontSize', 8);
        box(ax,'off');
    end

    % One legend per figure (set)
    lgd = legend(hL(isgraphics(hL)), legLabels(isgraphics(hL)), ...
        'Orientation','horizontal','Box','off','FontSize',10, ...
        'NumColumns', min(n_AMP, 6), 'Location','southoutside');
    lgd.Layout.Tile = 'south';
end