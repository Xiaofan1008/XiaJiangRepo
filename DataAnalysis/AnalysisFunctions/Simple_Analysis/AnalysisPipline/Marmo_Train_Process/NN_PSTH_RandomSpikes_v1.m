%% ============================================================
% NN_PSTH_RandomSpikes_v1.m
% Raster + PSTH for several RECOVERY FILES of the same dataset, drawn together, so the result of the
% partial (random-spikes) recovery can be compared with the full recovery (or with any other file).
% READ-ONLY: only plots and prints, saves nothing.
%
% spike_files : cell array of TrainRecovery files of the SAME dataset (same trials), e.g.
%                 *.sp_xia_TrainRecovery.mat               (v1, full recovery)
%                 *.sp_xia_TrainRecovery_RandomSpikes.mat  (NN_SpikeRecovery_RandomSpikes_v1)
%               file_labels gives the name of each file in the legends.
%
% mode = 'perChannel' : one figure per condition and channel
%                       top: raster, one block of trials per file (substituted spikes in red)
%                       middle: PSTH of every file on the same axes
%                       bottom (showDiff): PSTH of each file minus the PSTH of the FIRST file
%        'average'    : one figure per condition, PSTH averaged over the channels in chList
%                       (mean of the channels, shaded = SEM across channels), no raster
%
% COUNT TABLE (printed): for every condition and file, the mean spike count per trial and channel
%   in the window [last pulse + countWin_ms(1), last pulse + countWin_ms(2)] (the region where the
%   recovery differs between files), and the ratio to the first file.
% condList = [] prints the condition table and stops.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
spike_files = { ...
    '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_TrainRecovery.mat', ...
    '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_TrainRecovery_RandSpikes.mat'};
file_labels = {'v1 full recovery', 'RandomSpikes'};

condList = [15 45];        % condition numbers to plot ([] = print the table and stop)
chList   = 1:32;           % channels (index in sp_clipped); used for the figures and for the count table
mode     = 'average';   % 'perChannel' or 'average'
showDiff = true;           % extra panel: each file minus the first file

xRange     = [-20 180];        % [from to] ms after the trigger; [] = automatic (see next two)
xPre       = 10;        % automatic range: ms before the trigger
xPostLast  = 30;        % automatic range: ms after the LAST pulse of the longest chosen condition
binMs      = 1;         % PSTH bin width (ms)
smoothMs   = 3;         % Gaussian smoothing of the PSTH (ms, std); 0 = none
yMaxPSTH   = [];        % fixed top of the PSTH axis (spikes/s); [] = automatic
showPulses = true;      % dashed lines at the pulse times
maxTrials  = Inf;       % raster: at most this many trials per file (the first ones)
colorRecovered = true;  % raster: substituted spikes (0 to blank_end) in red
countWin_ms = [0 10];   % count table window relative to the last pulse of the condition (ms)
FS = 30000;

%% ================= LOAD =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(ismember(mode, {'perChannel','average'}), 'mode: perChannel or average.');
if ischar(spike_files), spike_files = {spike_files}; end
nF = numel(spike_files);
assert(numel(file_labels) == nF, 'file_labels must have one entry per file (%d).', nF);
SP = cell(nF,1);  BE = cell(nF,1);
for fi = 1:nF
    assert(isfile(spike_files{fi}), 'Cannot find %s.', spike_files{fi});
    S = load(spike_files{fi}, 'sp_clipped', 'blank_end_ms', 'Recovery_params', 'FS');
    assert(isfield(S,'sp_clipped') && isfield(S,'Recovery_params'), '%s is not a TrainRecovery file.', spike_files{fi});
    SP{fi} = S.sp_clipped;
    if colorRecovered && isfield(S,'blank_end_ms'), BE{fi} = S.blank_end_ms; end
    if fi == 1
        trialCond = S.Recovery_params.trialCond(:);
        condTable = S.Recovery_params.condTable;
        if isfield(S,'FS') && ~isempty(S.FS), FS = double(S.FS); end
    else
        assert(isequal(S.Recovery_params.trialCond(:), trialCond), 'File %d has a different trial list than file 1.', fi);
    end
    clear S
end
nChn = numel(SP{1});
assert(all(cellfun(@numel, SP) == nChn), 'The files have different numbers of channels.');
n_Trials = numel(trialCond);
nCond = numel(condTable);

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

%% ================= TABLE OF CONDITIONS =================
if isempty(condList)
    fprintf(' cond  trials  lastPulse(ms)  label\n');
    for k = 1:nCond
        fprintf(' %3d   %4d     %8.1f     %s\n', k, sum(trialCond == k & isfinite(t0)), lastP(k), labels{k});
    end
    fprintf('\nSet condList (and chList) at the top, then run again.\n');
    return;
end
assert(all(condList >= 1 & condList <= nCond), 'condList must be between 1 and %d.', nCond);
assert(all(chList >= 1 & chList <= nChn), 'chList must be between 1 and %d.', nChn);

%% ================= COUNT TABLE =================
fprintf('\n--- MEAN SPIKE COUNT PER TRIAL AND CHANNEL in [last pulse %+g, last pulse %+g] ms (channels %s) ---\n', ...
    countWin_ms(1), countWin_ms(2), mat2str(chList(:)'));
fprintf(' %-34s', 'condition');
for fi = 1:nF, fprintf(' | %-18s', file_labels{fi}); end
fprintf('\n');
for k = condList(:)'
    tr = find(trialCond == k & isfinite(t0));
    if isempty(tr), continue; end
    m = zeros(1, nF);
    for fi = 1:nF
        c = 0;
        for ch = chList(:)'
            s = SP{fi}{ch}(:,1);
            for i = 1:numel(tr)
                a = t0(tr(i)) + lastP(k) + countWin_ms(1);  b = t0(tr(i)) + lastP(k) + countWin_ms(2);
                c = c + (lb_(s, b) - lb_(s, a));
            end
        end
        m(fi) = c / (numel(tr) * numel(chList));
    end
    fprintf(' %-34s', labels{k}(1:min(end,34)));
    for fi = 1:nF
        if fi == 1 || m(1) == 0, fprintf(' | %8.3f          ', m(fi));
        else, fprintf(' | %8.3f (x%5.2f)  ', m(fi), m(fi)/m(1)); end
    end
    fprintf('\n');
end

if isempty(xRange), xRange = [-xPre, max(lastP(condList)) + xPostLast]; end
bins = xRange(1):binMs:xRange(2);
ctr  = bins(1:end-1) + binMs/2;
colF = lines(max(nF, 2));

%% ================= PLOT =================
if strcmp(mode, 'perChannel')
    for k = condList(:)'
        tr = find(trialCond == k & isfinite(t0));
        if isempty(tr), fprintf('cond %d has no trials.\n', k); continue; end
        trR = tr(1:min(numel(tr), maxTrials));
        for ch = chList(:)'
            figure('Color','w','Position',[100 80 950 880],'Name',sprintf('cond %d ch %d', k, ch));
            if showDiff && nF > 1
                ax1 = subplot(4,1,[1 2]);  ax2 = subplot(4,1,3);  ax3 = subplot(4,1,4);
            else
                ax1 = subplot(3,1,[1 2]);  ax2 = subplot(3,1,3);  ax3 = [];
            end
            hold(ax1,'on');  hold(ax2,'on');  if ~isempty(ax3), hold(ax3,'on'); end
            yOff = 0;  yt = [];  ytl = {};  hh = [];  R = cell(nF,1);
            for fi = 1:nF
                s = SP{fi}{ch}(:,1);
                beK = nan(numel(trR),1);  if ~isempty(BE{fi}), beK = BE{fi}{ch}(trR); end
                raster_(ax1, s, t0(trR), beK, xRange, yOff, [0 0 0], [0.85 0.15 0.15]);
                yt(end+1) = yOff + numel(trR)/2;  ytl{end+1} = file_labels{fi}; %#ok<SAGROW>
                yOff = yOff + numel(trR) + 1;
                R{fi} = rate_(s, t0(tr), bins, binMs, smoothMs);
                hh(end+1) = plot(ax2, ctr, R{fi}, 'Color', colF(fi,:), 'LineWidth', 1.4); %#ok<SAGROW>
            end
            set(ax1, 'YTick', yt, 'YTickLabel', ytl, 'TickLabelInterpreter', 'none'); ylim(ax1, [0 yOff]); xlim(ax1, xRange);
            ylabel(ax1, 'trials per file');
            title(ax1, sprintf('%s  |  ch %d  |  red = substituted spikes', labels{k}, ch), 'Interpreter','none');
            ylabel(ax2, 'spikes/s'); xlim(ax2, xRange);
            if ~isempty(yMaxPSTH), ylim(ax2, [0 yMaxPSTH]); end
            legend(ax2, hh, file_labels, 'Location','northeast','Box','off','Interpreter','none');
            if ~isempty(ax3)
                for fi = 2:nF, plot(ax3, ctr, R{fi} - R{1}, 'Color', colF(fi,:), 'LineWidth', 1.4); end
                yline(ax3, 0, ':', 'Color', [0.5 0.5 0.5]);
                ylabel(ax3, 'minus first file (spikes/s)'); xlim(ax3, xRange);
                xlabel(ax3, 'time from trigger (ms)');
            else
                xlabel(ax2, 'time from trigger (ms)');
            end
            if showPulses
                axs = [ax1 ax2 ax3];
                for ax = axs, drawPulses_(ax, pulseT{k}); end
            end
            linkaxes([ax1 ax2 ax3], 'x');
        end
    end
else
    for k = condList(:)'
        tr = find(trialCond == k & isfinite(t0));
        if isempty(tr), fprintf('cond %d has no trials.\n', k); continue; end
        nCh = numel(chList);
        RR = cell(nF,1);
        for fi = 1:nF
            RR{fi} = zeros(nCh, numel(ctr));
            for ci = 1:nCh
                RR{fi}(ci,:) = rate_(SP{fi}{chList(ci)}(:,1), t0(tr), bins, binMs, smoothMs);
            end
        end
        figure('Color','w','Position',[100 100 950 620],'Name',sprintf('cond %d average', k));
        if showDiff && nF > 1, ax2 = subplot(2,1,1);  ax3 = subplot(2,1,2); else, ax2 = subplot(1,1,1);  ax3 = []; end
        hold(ax2,'on');  if ~isempty(ax3), hold(ax3,'on'); end
        hh = [];
        for fi = 1:nF
            hh(end+1) = band_(ax2, ctr, RR{fi}, colF(fi,:)); %#ok<SAGROW>
        end
        ylabel(ax2, 'spikes/s (mean over channels)'); xlim(ax2, xRange);
        if ~isempty(yMaxPSTH), ylim(ax2, [0 yMaxPSTH]); end
        legend(ax2, hh, file_labels, 'Location','northeast','Box','off','Interpreter','none');
        title(ax2, sprintf('%s  |  mean of %d channels (shaded = SEM across channels)', labels{k}, nCh), 'Interpreter','none');
        if ~isempty(ax3)
            for fi = 2:nF, band_(ax3, ctr, RR{fi} - RR{1}, colF(fi,:)); end
            yline(ax3, 0, ':', 'Color', [0.5 0.5 0.5]);
            ylabel(ax3, 'minus first file (spikes/s)'); xlim(ax3, xRange);
            xlabel(ax3, 'time from trigger (ms)');
        else
            xlabel(ax2, 'time from trigger (ms)');
        end
        if showPulses
            for ax = [ax2 ax3], drawPulses_(ax, pulseT{k}); end
        end
        linkaxes([ax2 ax3], 'x');
    end
end

%% ================= HELPERS =================
function i = lb_(s, x)
% first index with s(i) >= x (numel(s)+1 if none); s sorted ascending
lo = 1; hi = numel(s) + 1;
while lo < hi
    m = floor((lo + hi)/2);
    if s(m) < x, lo = m + 1; else, hi = m; end
end
i = lo;
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

function h = band_(ax, x, M, col)
% mean over rows of M with a shaded SEM band; returns the handle of the mean line
mu = mean(M, 1);
se = std(M, 0, 1) / sqrt(max(size(M,1), 1));
patch(ax, [x fliplr(x)], [mu+se fliplr(mu-se)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h = plot(ax, x, mu, 'Color', col, 'LineWidth', 1.5);
end

function drawPulses_(ax, p)
if isempty(ax) || ~isgraphics(ax), return; end
yl = ylim(ax);
for tp = p(:)'
    h = plot(ax, [tp tp], yl, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    set(get(get(h,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
end
ylim(ax, yl);
end
