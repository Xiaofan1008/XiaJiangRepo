%% ============================================================
% NN_PSTH_Random_v1.m
% Raster + PSTH for the RANDOM-subset recovery (NN_SpikeRecovery_Random_v1).
% READ-ONLY: only plots, saves nothing. Based on NN_RasterPSTH_v1.
%
% In the Random recovery only a random share of the trials of each level is recovered; the other
% trials keep the original data with the blank gap. This script shows, for every chosen condition
% and channel, the RECOVERED trials and the NOT-recovered trials separately, so you can compare them:
%   raster : recovered trials on top (substituted spikes 0 .. blank_end in red, green dot = window end),
%            then the not-recovered trials (all black), one row per trial
%   PSTH   : recovered trials (solid line) and not-recovered trials (dashed line), mean over the trials
%            of each group, number of trials in the legend
%
% trialSet = 'both'         recovered and not recovered, separately (default)
%            'recovered'    recovered trials only
%            'notRecovered' not-recovered trials only
%            'all'          all trials of the condition together (one PSTH)
% mode     = 'overlay'  one figure per channel, all chosen conditions together (rasters stacked)
%            'separate' one figure per condition and channel (raster on top, PSTH below)
% condList = [] prints the table of conditions with their recovered / total trials and stops.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
spike_file  = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_TrainRecovery_Random.mat';

condList = [15 45];        % condition numbers to plot ([] = print the table and stop)
chList   = 1:10;           % channels to plot (index in sp_clipped)
mode     = 'overlay';      % 'overlay' or 'separate'
trialSet = 'both';         % 'both' | 'recovered' | 'notRecovered' | 'all'

xRange     = [-20 100];        % [from to] ms after the trigger; [] = automatic (see next two)
xPre       = 10;        % automatic range: ms before the trigger
xPostLast  = 50;        % automatic range: ms after the LAST pulse of the longest chosen condition
binMs      = 1;         % PSTH bin width (ms)
smoothMs   = 3;         % Gaussian smoothing of the PSTH (ms, std); 0 = none
yMaxPSTH   = [];        % fixed top of the PSTH axis (spikes/s); [] = automatic
showPulses = true;      % dashed lines at the pulse times
maxTrials  = Inf;       % show at most this many trials per group (the first ones)
colorRecovered = true;  % substituted spikes (0 to blank_end) in red
FS = 30000;

%% ================= LOAD =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(isfile(spike_file), 'Cannot find %s.', spike_file);
assert(ismember(trialSet, {'both','recovered','notRecovered','all'}), 'trialSet: both, recovered, notRecovered or all.');
assert(ismember(mode, {'overlay','separate'}), 'mode: overlay or separate.');
S = load(spike_file, 'sp_clipped', 'blank_end_ms', 'isRecovered', 'Recovery_params', 'FS');
assert(isfield(S,'sp_clipped') && isfield(S,'Recovery_params'), '%s is not a TrainRecovery file.', spike_file);
sp = S.sp_clipped;  if isfield(S,'FS') && ~isempty(S.FS), FS = double(S.FS); end
nChn = numel(sp);
trialCond = S.Recovery_params.trialCond(:);
condTable = S.Recovery_params.condTable;
if isfield(S,'isRecovered'), isRec = logical(S.isRecovered(:));
elseif isfield(S.Recovery_params,'isRecovered'), isRec = logical(S.Recovery_params.isRecovered(:));
else, error('%s has no isRecovered: use the file written by NN_SpikeRecovery_Random_v1.', spike_file); end
be = [];
if colorRecovered && isfield(S,'blank_end_ms'), be = S.blank_end_ms; end
clear S
n_Trials = numel(trialCond);
nCond = numel(condTable);
assert(numel(isRec) == n_Trials, 'isRecovered has %d entries but there are %d trials.', numel(isRec), n_Trials);

f = dir(fullfile(data_folder, '*.trig.dat'));  f = f(~startsWith({f.name}, '.'));
assert(~isempty(f), 'No *.trig.dat in %s.', data_folder);
fid = fopen(fullfile(data_folder, f(1).name), 'r');  trig = fread(fid, [1 inf], 'double');  fclose(fid);
nUse = min(numel(trig), n_Trials);
t0 = nan(n_Trials,1);  t0(1:nUse) = trig(1:nUse)/FS*1000;

fe = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));  fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fullfile(data_folder, fe(1).name), 'StimParams', 'simultaneous_stim');
nEl = double(E.simultaneous_stim);

labels = arrayfun(@(k) sprintf('cond %d  %s', k, char(condTable(k).label)), 1:nCond, 'UniformOutput', false);

% pulse times per condition (ms), from the exp file
pulseT = cell(nCond,1);  lastP = zeros(nCond,1);
for k = 1:nCond
    t1 = find(trialCond == k, 1);
    if isempty(t1), continue; end
    for e = 1:nEl
        row = 1 + (t1-1)*nEl + e;
        if double(E.StimParams{row,16}) > 0
            n = double(E.StimParams{row,8});  per = double(E.StimParams{row,9});  p0 = double(E.StimParams{row,6});
            pulseT{k} = [pulseT{k}, (p0 + (0:n-1)*per)/1000];
        end
    end
    pulseT{k} = unique(pulseT{k});
    if ~isempty(pulseT{k}), lastP(k) = max(pulseT{k}); end
end

%% ================= TABLE =================
if isempty(condList)
    fprintf(' cond  recovered / valid trials  lastPulse(ms)  label\n');
    for k = 1:nCond
        ok = trialCond == k & isfinite(t0);
        fprintf(' %3d       %3d / %3d          %8.1f     %s\n', k, sum(ok & isRec), sum(ok), lastP(k), labels{k});
    end
    fprintf('\nSet condList (and chList) at the top, then run again.\n');
    return;
end
assert(all(condList >= 1 & condList <= nCond), 'condList must be between 1 and %d.', nCond);
assert(all(chList >= 1 & chList <= nChn), 'chList must be between 1 and %d.', nChn);

fprintf('Chosen conditions (recovered trials = trial number within the condition):\n');
for k = condList(:)'
    allk = find(trialCond == k);
    fprintf(' %s : %d recovered, %d not recovered | recovered trials %s\n', labels{k}, sum(isRec(allk)), sum(~isRec(allk)), fmtList_(find(isRec(allk))'));
end

if isempty(xRange), xRange = [-xPre, max(lastP(condList)) + xPostLast]; end
bins = xRange(1):binMs:xRange(2);
ctr  = bins(1:end-1) + binMs/2;
switch trialSet
    case 'both', setList = {'recovered','notRecovered'};
    otherwise,   setList = {trialSet};
end

%% ================= PLOT =================
cols = lines(max(numel(condList), 2));
for ch = chList(:)'
    s = sp{ch}(:,1);
    if strcmp(mode, 'overlay')
        figure('Color','w','Position',[100 100 950 800],'Name',sprintf('ch %d overlay', ch));
        ax1 = subplot(3,1,[1 2]);  hold(ax1,'on');  ax2 = subplot(3,1,3);  hold(ax2,'on');
        yOff = 0;  yt = [];  ytl = {};  hh = [];  hn = {};
        for ci = 1:numel(condList)
            k = condList(ci);
            for si = 1:numel(setList)
                tr = trialsOf_(trialCond, t0, isRec, k, setList{si}, maxTrials);
                if isempty(tr), continue; end
                beK = nan(numel(tr),1);  if ~isempty(be), beK = be{ch}(tr); end
                raster_(ax1, s, t0(tr), beK, xRange, yOff, cols(ci,:), [0.85 0.15 0.15]);
                yt(end+1) = yOff + numel(tr)/2;  ytl{end+1} = sprintf('%d %s', k, tag_(setList{si})); %#ok<SAGROW>
                yOff = yOff + numel(tr) + 1;
                ls = '-';  if strcmp(setList{si}, 'notRecovered'), ls = '--'; end
                hh(end+1) = plot(ax2, ctr, rate_(s, t0(tr), bins, binMs, smoothMs), ls, 'Color', cols(ci,:), 'LineWidth', 1.3); %#ok<SAGROW>
                hn{end+1} = sprintf('cond %d %s (n=%d)', k, tag_(setList{si}), numel(tr)); %#ok<SAGROW>
            end
        end
        if yOff == 0, close(gcf); fprintf('ch %d: no trials to show.\n', ch); continue; end
        set(ax1, 'YTick', yt, 'YTickLabel', ytl); ylabel(ax1, 'condition and group (trials stacked)'); ylim(ax1,[0 yOff]); xlim(ax1, xRange);
        title(ax1, sprintf('channel %d  |  %s  |  red = substituted spikes', ch, shortName_(spike_file)), 'Interpreter','none');
        xlabel(ax2,'time from trigger (ms)'); ylabel(ax2,'spikes/s'); xlim(ax2, xRange);
        if ~isempty(yMaxPSTH), ylim(ax2,[0 yMaxPSTH]); end
        legend(ax2, hh, hn, 'Location','northeast','Box','off','FontSize',8,'Interpreter','none');
        if showPulses
            allP = unique([pulseT{condList}]);
            for ax = [ax1 ax2], drawPulses_(ax, allP); end
        end
        linkaxes([ax1 ax2], 'x');
    else
        for k = condList(:)'
            figure('Color','w','Position',[100 100 900 650],'Name',sprintf('cond %d ch %d', k, ch));
            ax1 = subplot(2,1,1);  hold(ax1,'on');  ax2 = subplot(2,1,2);  hold(ax2,'on');
            yOff = 0;  yt = [];  ytl = {};  hh = [];  hn = {};
            for si = 1:numel(setList)
                tr = trialsOf_(trialCond, t0, isRec, k, setList{si}, maxTrials);
                if isempty(tr), continue; end
                beK = nan(numel(tr),1);  if ~isempty(be), beK = be{ch}(tr); end
                raster_(ax1, s, t0(tr), beK, xRange, yOff, [0 0 0], [0.85 0.15 0.15]);
                yt(end+1) = yOff + numel(tr)/2;  ytl{end+1} = tag_(setList{si}); %#ok<SAGROW>
                yOff = yOff + numel(tr) + 1;
                if strcmp(setList{si}, 'notRecovered'), ls = '--'; c = [0 0.45 0.74]; else, ls = '-'; c = [0 0 0]; end
                hh(end+1) = plot(ax2, ctr, rate_(s, t0(tr), bins, binMs, smoothMs), ls, 'Color', c, 'LineWidth', 1.5); %#ok<SAGROW>
                hn{end+1} = sprintf('%s (n=%d)', tag_(setList{si}), numel(tr)); %#ok<SAGROW>
            end
            if yOff == 0, close(gcf); fprintf('cond %d ch %d: no trials to show.\n', k, ch); continue; end
            set(ax1, 'YTick', yt, 'YTickLabel', ytl); ylim(ax1,[0 yOff]); xlim(ax1, xRange); ylabel(ax1,'trials');
            title(ax1, sprintf('%s  |  ch %d  |  %s', labels{k}, ch, shortName_(spike_file)), 'Interpreter','none');
            xlabel(ax2,'time from trigger (ms)'); ylabel(ax2,'spikes/s'); xlim(ax2, xRange);
            if ~isempty(yMaxPSTH), ylim(ax2,[0 yMaxPSTH]); end
            legend(ax2, hh, hn, 'Location','northeast','Box','off');
            if showPulses
                for ax = [ax1 ax2], drawPulses_(ax, pulseT{k}); end
            end
            linkaxes([ax1 ax2], 'x');
        end
    end
end

%% ================= HELPERS =================
function tr = trialsOf_(trialCond, t0, isRec, k, setName, maxTrials)
ok = trialCond == k & isfinite(t0);
switch setName
    case 'recovered',    ok = ok & isRec;
    case 'notRecovered', ok = ok & ~isRec;
end
tr = find(ok);
if numel(tr) > maxTrials, tr = tr(1:maxTrials); end
end

function t = tag_(setName)
switch setName
    case 'recovered',    t = 'recovered';
    case 'notRecovered', t = 'not recovered';
    otherwise,           t = 'all trials';
end
end

function raster_(ax, s, t0s, bes, xl, yOff, colOrig, colRec)
for i = 1:numel(t0s)
    rel = s(s >= t0s(i)+xl(1) & s <= t0s(i)+xl(2)) - t0s(i);
    if isempty(rel), continue; end
    isRec = isfinite(bes(i)) & rel >= 0 & rel < bes(i);
    for g = 1:2
        if g == 1, m = ~isRec; c = colOrig; else, m = isRec; c = colRec; end
        x = rel(m);
        if isempty(x), continue; end
        n = numel(x);
        X = [x(:)'; x(:)'; nan(1,n)];
        Y = [yOff+i-0.4+zeros(1,n); yOff+i+0.4+zeros(1,n); nan(1,n)];
        plot(ax, X(:), Y(:), 'Color', c, 'LineWidth', 0.8);
    end
    if isfinite(bes(i)), plot(ax, bes(i), yOff+i, '.', 'Color', [0 0.7 0], 'MarkerSize', 8); end
end
end

function r = rate_(s, t0s, bins, binMs, smoothMs)
cnt = zeros(1, numel(bins)-1);
for i = 1:numel(t0s)
    rel = s(s >= t0s(i)+bins(1) & s < t0s(i)+bins(end)) - t0s(i);
    cnt = cnt + histcounts(rel, bins);
end
r = cnt / numel(t0s) / binMs * 1000;
if smoothMs > 0
    sd = smoothMs / binMs;
    x = -ceil(4*sd):ceil(4*sd);
    g = exp(-0.5*(x/sd).^2);  g = g / sum(g);
    r = conv(r, g, 'same');
end
end

function drawPulses_(ax, p)
yl = ylim(ax);
for tp = p(:)'
    h = plot(ax, [tp tp], yl, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
end
ylim(ax, yl);
end

function n = shortName_(p)
[~, n, e] = fileparts(p);  n = [n e];
end

function s = fmtList_(v)
% compact MATLAB-style list, e.g. [1 2 4:7]
v = sort(v(:)');
if isempty(v), s = '[]'; return; end
parts = {};  i = 1;
while i <= numel(v)
    j = i;
    while j < numel(v) && v(j+1) == v(j) + 1, j = j + 1; end
    if j == i,        parts{end+1} = sprintf('%d', v(i));          %#ok<AGROW>
    elseif j == i+1,  parts{end+1} = sprintf('%d %d', v(i), v(j)); %#ok<AGROW>
    else,             parts{end+1} = sprintf('%d:%d', v(i), v(j)); %#ok<AGROW>
    end
    i = j + 1;
end
s = ['[' strjoin(parts, ' ') ']'];
end
