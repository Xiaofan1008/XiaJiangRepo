%% ============================================================
% NN_RasterPSTH_v1.m
% Raster + PSTH for any NN spike file, you choose the file, conditions and channels.
% READ-ONLY: only plots, saves nothing.
%
% spike_file can be any file that holds sp_clipped, for example
%   *.sp_xia.mat                 (unfiltered)
%   *.sp_xia_filtered.mat        (after waveform filtering)
%   *.sp_xia_TrainRecovery.mat   (after recovery; substituted spikes can be drawn in red)
%
% mode = 'separate' : one figure per condition and channel (raster on top, PSTH below)
% mode = 'overlay'  : one figure per channel, all chosen conditions together
%                     (rasters stacked by condition, PSTHs on the same axes)
% condList = [] prints the condition table and stops.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
spike_file  = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_TrainRecovery.mat';

condList = [15 45 49];          % condition numbers to plot ([] = print the table and stop)
chList   = [1:10];       % channels to plot (index in sp_clipped)
mode     = 'overlay';  % 'separate' or 'overlay'

xRange     = [];        % [from to] ms after the trigger; [] = automatic (see next two)
xPre       = 10;        % automatic range: ms before the trigger
xPostLast  = 30;        % automatic range: ms after the LAST pulse of the longest chosen condition
binMs      = 1;         % PSTH bin width (ms)
smoothMs   = 0;         % Gaussian smoothing of the PSTH (ms, std); 0 = none
yMaxPSTH   = [];        % fixed top of the PSTH axis (spikes/s); [] = automatic
showPulses = true;      % dashed lines at the pulse times
maxTrials  = Inf;       % show at most this many trials per condition (the first ones)
colorRecovered = true;  % TrainRecovery file only: substituted spikes (0 to blank_end) in red
FS = 30000;

%% ================= LOAD =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
cd(data_folder);
assert(isfile(spike_file), 'Cannot find %s.', spike_file);
S = load(spike_file);
assert(isfield(S,'sp_clipped'), 'sp_clipped not found in %s.', spike_file);
sp = S.sp_clipped;  if isfield(S,'fs'), FS = S.fs; elseif isfield(S,'FS'), FS = S.FS; end
nChn = numel(sp);
be = [];
if colorRecovered && isfield(S,'blank_end_ms'), be = S.blank_end_ms; end
clear S

f = dir('*.trig.dat'); f = f(~startsWith({f.name}, '.'));
assert(~isempty(f), 'No *.trig.dat in %s.', data_folder);
fid = fopen(f(1).name, 'r'); trig = fread(fid, [1 inf], 'double'); fclose(fid);

fe = dir('*_exp_datafile_*.mat'); fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fe(1).name, 'TrialParams', 'StimParams', 'simultaneous_stim', 'n_Trials');
nEl = E.simultaneous_stim;
condAll = cell2mat(E.TrialParams(2:end,2));
trialCond = condAll(1:nEl:end);  trialCond = trialCond(:);
n_Trials = numel(trialCond);
nUse = min(numel(trig), n_Trials);
t0 = nan(n_Trials,1);  t0(1:nUse) = trig(1:nUse)/FS*1000;
nCond = max(trialCond);

labels = arrayfun(@(k) sprintf('cond %d', k), 1:nCond, 'UniformOutput', false);
fc = dir('*.conds.mat'); fc = fc(~startsWith({fc.name}, '.'));
if ~isempty(fc)
    K = load(fc(1).name);
    if isfield(K,'cond') && numel(K.cond) >= nCond && isfield(K.cond,'label')
        for k = 1:nCond, labels{k} = sprintf('cond %d  %s', k, K.cond(k).label); end
    end
end

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

if isempty(condList)
    fprintf(' cond  trials  lastPulse(ms)  label\n');
    for k = 1:nCond
        fprintf(' %3d   %4d     %8.1f     %s\n', k, sum(trialCond == k), lastP(k), labels{k});
    end
    fprintf('\nSet condList (and chList) at the top, then run again.\n');
    return;
end
assert(all(condList >= 1 & condList <= nCond), 'condList must be between 1 and %d.', nCond);
assert(all(chList >= 1 & chList <= nChn), 'chList must be between 1 and %d.', nChn);

if isempty(xRange), xRange = [-xPre, max(lastP(condList)) + xPostLast]; end
bins = xRange(1):binMs:xRange(2);
ctr  = bins(1:end-1) + binMs/2;

%% ================= PLOT =================
cols = lines(max(numel(condList), 2));
for ch = chList(:)'
    s = sp{ch}(:,1);
    if strcmp(mode, 'overlay')
        figure('Color','w','Position',[100 100 950 800],'Name',sprintf('ch %d overlay', ch));
        ax1 = subplot(3,1,[1 2]);  hold(ax1,'on');  ax2 = subplot(3,1,3);  hold(ax2,'on');
        yOff = 0;  yt = [];  ytl = {};  hh = [];
        for ci = 1:numel(condList)
            k = condList(ci);
            tr = trialsOf_(trialCond, t0, k, maxTrials);
            if isempty(tr), continue; end
            beK = nan(numel(tr),1);  if ~isempty(be), beK = be{ch}(tr); end
            raster_(ax1, s, t0(tr), beK, xRange, yOff, cols(ci,:), [0.85 0.15 0.15]);
            yt(end+1) = yOff + numel(tr)/2;  ytl{end+1} = sprintf('%d', k); %#ok<SAGROW>
            yOff = yOff + numel(tr) + 1;
            hh(end+1) = plot(ax2, ctr, rate_(s, t0(tr), bins, binMs, smoothMs), 'Color', cols(ci,:), 'LineWidth', 1.3); %#ok<SAGROW>
        end
        set(ax1, 'YTick', yt, 'YTickLabel', ytl); ylabel(ax1, 'condition (trials stacked)'); ylim(ax1,[0 yOff]); xlim(ax1, xRange);
        title(ax1, sprintf('channel %d  |  %s', ch, shortName_(spike_file)), 'Interpreter','none');
        xlabel(ax2,'time from trigger (ms)'); ylabel(ax2,'spikes/s'); xlim(ax2, xRange);
        if ~isempty(yMaxPSTH), ylim(ax2,[0 yMaxPSTH]); end
        legend(ax2, hh, arrayfun(@(k) sprintf('cond %d', k), condList(1:numel(hh)), 'UniformOutput', false), 'Location','northeast','Box','off');
        if showPulses
            allP = unique([pulseT{condList}]);
            for ax = [ax1 ax2], drawPulses_(ax, allP); end
        end
        linkaxes([ax1 ax2], 'x');
    else
        for k = condList(:)'
            tr = trialsOf_(trialCond, t0, k, maxTrials);
            if isempty(tr), fprintf('cond %d has no trials.\n', k); continue; end
            beK = nan(numel(tr),1);  if ~isempty(be), beK = be{ch}(tr); end
            figure('Color','w','Position',[100 100 900 650],'Name',sprintf('cond %d ch %d', k, ch));
            ax1 = subplot(2,1,1);  hold(ax1,'on');
            raster_(ax1, s, t0(tr), beK, xRange, 0, [0 0 0], [0.85 0.15 0.15]);
            ylim(ax1,[0 numel(tr)+1]); xlim(ax1, xRange); ylabel(ax1,'trial');
            title(ax1, sprintf('%s  |  ch %d  |  %s', labels{k}, ch, shortName_(spike_file)), 'Interpreter','none');
            ax2 = subplot(2,1,2);  hold(ax2,'on');
            plot(ax2, ctr, rate_(s, t0(tr), bins, binMs, smoothMs), 'k', 'LineWidth', 1.5);
            xlabel(ax2,'time from trigger (ms)'); ylabel(ax2,'spikes/s'); xlim(ax2, xRange);
            if ~isempty(yMaxPSTH), ylim(ax2,[0 yMaxPSTH]); end
            if showPulses
                for ax = [ax1 ax2], drawPulses_(ax, pulseT{k}); end
            end
            linkaxes([ax1 ax2], 'x');
        end
    end
end

%% ================= HELPERS =================
function tr = trialsOf_(trialCond, t0, k, maxTrials)
tr = find(trialCond == k & isfinite(t0));
if numel(tr) > maxTrials, tr = tr(1:maxTrials); end
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
