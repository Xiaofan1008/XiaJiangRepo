%% ============================================================
% NN_PSTH_RandomEvent_v1.m
% Check of the NN_SpikeRecovery_RandomEvent_v1 simulation: raster / PSTH of several recovery files drawn
% together, PLUS a per-pulse check that the data were thinned by the percentage you set for each level.
% READ-ONLY: only plots and prints, saves nothing.
%
% spike_files : TrainRecovery files of the SAME dataset, e.g.
%                 *.sp_xia_TrainRecovery.mat                  (v1, full recovery = reference, FIRST in the list)
%                 *.sp_xia_TrainRecovery_RandomEvent.mat      (NN_SpikeRecovery_RandomEvent_v1)
%               file_labels = names in the legends; paramFile = index of the file written by
%               NN_SpikeRecovery_RandomEvent_v1 (it holds the percent and event_ms of every level).
%
% FIGURE GROUP 1 (Plot_PSTH)  PSTH of every file on the same axes + (showDiff) each file minus the first
%   mode = 'average'    : PSTH averaged over the channels in chList (shaded = SEM across channels)
%          'perChannel' : one figure per condition and channel with a raster block per file
%
% FIGURE GROUP 2 (Plot_PerPulse)  per-pulse check, one figure per condition in condList
%   For every pulse j of the (final) train: the mean spike count per trial and channel in the window
%   after that pulse, for each file, and the RATIO file / first file (in %). The ratio of the
%   RandomEvent file should follow the percent you set for the level of that pulse (black line).
%   RAW counts are used (no baseline subtraction), because the simulation thins ALL spikes in the window.
%   Window: pulseWin_ms = 'event' -> [0, event_ms of that level] taken from the RandomEvent file
%           pulseWin_ms = [a b]   -> the same window [a, b] ms after every pulse
%   Pulse j belongs to level j of the chain (level 1 = first pulse, never thinned = 100 %).
%
% condList = [] prints the table of conditions and stops.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
spike_files = { ...
    '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_TrainRecovery.mat', ...
    '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_TrainRecovery_RandomEvent.mat'};
file_labels = {'v1 full recovery', 'RandomEvent'};
paramFile   = 2;           % which file holds percent_k / event_ms_k (the RandomEvent file)

condList = [15 45];        % FINAL condition numbers to plot ([] = print the table and stop)
chList   = 1:32;           % channels used for everything
mode     = 'perChannel';      % PSTH figures: 'average' or 'perChannel'
Plot_PSTH     = true;      % figure group 1
Plot_PerPulse = true;      % figure group 2
showDiff = true;           % PSTH figures: extra panel, each file minus the first file

% ---- figure group 1 ----
xRange     = [-20 180];        % [from to] ms after the trigger; [] = automatic (see next two)
xPre       = 10;        % automatic range: ms before the trigger
xPostLast  = 30;        % automatic range: ms after the LAST pulse of the longest chosen condition
binMs      = 1;         % PSTH bin width (ms)
smoothMs   = 3;         % Gaussian smoothing of the PSTH (ms, std); 0 = none
yMaxPSTH   = [];        % fixed top of the PSTH axis (spikes/s); [] = automatic
showPulses = true;      % dashed lines at the pulse times
maxTrials  = Inf;       % perChannel raster: at most this many trials per file
colorRecovered = true;  % perChannel raster: substituted spikes (0 to blank_end) in red

% ---- figure group 2 ----
pulseWin_ms = 'event';  % 'event' or [a b] (ms after each pulse)
xAxisPulse  = 'pulse';  % x axis of the per-pulse figure: 'pulse' (pulse number) or 'time' (ms)
FS = 30000;

%% ================= LOAD =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(ismember(mode, {'perChannel','average'}), 'mode: perChannel or average.');
assert(ismember(xAxisPulse, {'pulse','time'}), 'xAxisPulse: pulse or time.');
if ischar(spike_files), spike_files = {spike_files}; end
nF = numel(spike_files);
assert(numel(file_labels) == nF, 'file_labels must have one entry per file (%d).', nF);
assert(paramFile >= 1 && paramFile <= nF, 'paramFile must be between 1 and %d.', nF);
SP = cell(nF,1);  BE = cell(nF,1);  RP = cell(nF,1);
for fi = 1:nF
    assert(isfile(spike_files{fi}), 'Cannot find %s.', spike_files{fi});
    S = load(spike_files{fi}, 'sp_clipped', 'blank_end_ms', 'Recovery_params', 'FS');
    assert(isfield(S,'sp_clipped') && isfield(S,'Recovery_params'), '%s is not a TrainRecovery file.', spike_files{fi});
    SP{fi} = S.sp_clipped;  RP{fi} = S.Recovery_params;
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
nCh = numel(chList);

if isempty(xRange), xRange = [-xPre, max(lastP(condList)) + xPostLast]; end
bins = xRange(1):binMs:xRange(2);
ctr  = bins(1:end-1) + binMs/2;
colF = lines(max(nF, 2));

%% ================= FIGURE GROUP 1: PSTH =================
if Plot_PSTH
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
                yOff = 0;  yt = [];  ytl = {};  hh = [];  Rt = cell(nF,1);
                for fi = 1:nF
                    s = SP{fi}{ch}(:,1);
                    beK = nan(numel(trR),1);  if ~isempty(BE{fi}), beK = BE{fi}{ch}(trR); end
                    raster_(ax1, s, t0(trR), beK, xRange, yOff, [0 0 0], [0.85 0.15 0.15]);
                    yt(end+1) = yOff + numel(trR)/2;  ytl{end+1} = file_labels{fi}; %#ok<SAGROW>
                    yOff = yOff + numel(trR) + 1;
                    Rt{fi} = rate_(s, t0(tr), bins, binMs, smoothMs);
                    hh(end+1) = plot(ax2, ctr, Rt{fi}, 'Color', colF(fi,:), 'LineWidth', 1.4); %#ok<SAGROW>
                end
                set(ax1, 'YTick', yt, 'YTickLabel', ytl, 'TickLabelInterpreter', 'none'); ylim(ax1, [0 yOff]); xlim(ax1, xRange);
                ylabel(ax1, 'trials per file');
                title(ax1, sprintf('%s  |  ch %d  |  red = substituted spikes', labels{k}, ch), 'Interpreter','none');
                ylabel(ax2, 'spikes/s'); xlim(ax2, xRange);
                if ~isempty(yMaxPSTH), ylim(ax2, [0 yMaxPSTH]); end
                legend(ax2, hh, file_labels, 'Location','northeast','Box','off','Interpreter','none');
                if ~isempty(ax3)
                    for fi = 2:nF, plot(ax3, ctr, Rt{fi} - Rt{1}, 'Color', colF(fi,:), 'LineWidth', 1.4); end
                    yline(ax3, 0, ':', 'Color', [0.5 0.5 0.5]);
                    ylabel(ax3, 'minus first file (spikes/s)'); xlim(ax3, xRange);
                    xlabel(ax3, 'time from trigger (ms)');
                else
                    xlabel(ax2, 'time from trigger (ms)');
                end
                if showPulses
                    for ax = [ax1 ax2 ax3], drawPulses_(ax, pulseT{k}); end
                end
                linkaxes([ax1 ax2 ax3], 'x');
            end
        end
    else
        for k = condList(:)'
            tr = find(trialCond == k & isfinite(t0));
            if isempty(tr), fprintf('cond %d has no trials.\n', k); continue; end
            RR = cell(nF,1);
            for fi = 1:nF
                RR{fi} = zeros(nCh, numel(ctr));
                for ci = 1:nCh
                    RR{fi}(ci,:) = rate_(SP{fi}{chList(ci)}(:,1), t0(tr), bins, binMs, smoothMs);
                end
            end
            figure('Color','w','Position',[100 100 950 620],'Name',sprintf('cond %d average PSTH', k));
            if showDiff && nF > 1, ax2 = subplot(2,1,1);  ax3 = subplot(2,1,2); else, ax2 = subplot(1,1,1);  ax3 = []; end
            hold(ax2,'on');  if ~isempty(ax3), hold(ax3,'on'); end
            hh = [];
            for fi = 1:nF, hh(end+1) = band_(ax2, ctr, RR{fi}, colF(fi,:)); end %#ok<SAGROW>
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
end

%% ================= FIGURE GROUP 2: PER-PULSE CHECK =================
if Plot_PerPulse
    rp = RP{paramFile};
    havePct = isfield(rp,'percent_k') && isfield(rp,'event_ms_k');
    if ~havePct
        warning('RE:params', 'File %d has no percent_k / event_ms_k (not a RandomEvent file): no expected line, and the window must be numeric.', paramFile);
        assert(isnumeric(pulseWin_ms), 'pulseWin_ms = ''event'' needs the RandomEvent file as paramFile; set pulseWin_ms = [a b].');
    end
    for k = condList(:)'
        tr = find(trialCond == k & isfinite(t0));
        if isempty(tr), fprintf('cond %d has no trials.\n', k); continue; end
        % chain of levels ending at k: pulse j <-> level j
        ch_ = k;
        while isfinite(condTable(ch_(1)).sourceCond), ch_ = [condTable(ch_(1)).sourceCond, ch_]; end %#ok<AGROW>
        L = numel(ch_);
        tj = arrayfun(@(c) double(condTable(c).lastPulse_ms), ch_);
        pct = 100*ones(1, L);  ev = nan(1, L);
        if havePct
            for j = 2:L
                pct(j) = double(rp.percent_k(ch_(j)));
                ev(j)  = double(rp.event_ms_k(ch_(j)));
            end
            ev(1) = ev(2);                                  % level 1 has no setting: use the one of level 2
        end
        if isnumeric(pulseWin_ms)
            w1 = pulseWin_ms(1)*ones(1,L);  w2 = pulseWin_ms(2)*ones(1,L);
        else
            assert(strcmp(pulseWin_ms, 'event'), 'pulseWin_ms must be ''event'' or [a b].');
            w1 = zeros(1,L);  w2 = ev;
        end
        % counts per file and pulse (mean per trial and channel; SEM across channels)
        M = nan(nF, L);  SE = nan(nF, L);  T = nan(nF, L);
        for fi = 1:nF
            for j = 1:L
                cc = zeros(nCh, 1);
                for ci = 1:nCh
                    s = SP{fi}{chList(ci)}(:,1);
                    for i = 1:numel(tr)
                        a = t0(tr(i)) + tj(j) + w1(j);  b = t0(tr(i)) + tj(j) + w2(j);
                        cc(ci) = cc(ci) + (lb_(s, b) - lb_(s, a));
                    end
                end
                T(fi,j) = sum(cc);
                cc = cc / numel(tr);
                M(fi,j) = mean(cc);  SE(fi,j) = std(cc) / sqrt(nCh);
            end
        end
        ratio = 100 * T ./ T(1,:);                          % total counts, in % of the first file
        if strcmp(xAxisPulse, 'pulse'), xv = 1:L;  xl = 'pulse number (= level)'; else, xv = tj;  xl = 'time of the pulse (ms)'; end

        fprintf('\n--- PER-PULSE CHECK: %s (channels %s, %d trials) ---\n', labels{k}, mat2str(chList(:)'), numel(tr));
        fprintf('  pulse  t(ms)  window(ms)        expected%%   ');
        for fi = 1:nF, fprintf(' | %-16s', file_labels{fi}); end
        fprintf('\n');
        for j = 1:L
            fprintf('  %4d  %6.1f  [%4.1f, %5.1f]    %6.1f   ', j, tj(j), w1(j), w2(j), pct(j));
            for fi = 1:nF
                if fi == 1, fprintf(' | %7.3f          ', M(fi,j));
                else, fprintf(' | %7.3f (%5.1f%%) ', M(fi,j), ratio(fi,j)); end
            end
            fprintf('\n');
        end

        figure('Color','w','Position',[100 80 1000 700],'Name',sprintf('cond %d per pulse', k));
        axA = subplot(2,1,1);  hold(axA,'on');  hA = [];
        for fi = 1:nF
            hA(end+1) = errorbar(axA, xv, M(fi,:), SE(fi,:), 'o-', 'Color', colF(fi,:), 'MarkerFaceColor', colF(fi,:), 'LineWidth', 1.3); %#ok<SAGROW>
        end
        ylabel(axA, 'spikes per trial and channel'); xlabel(axA, xl);
        title(axA, sprintf('%s  |  count in the window after each pulse (mean +/- SEM across %d channels)', labels{k}, nCh), 'Interpreter','none');
        legend(axA, hA, file_labels, 'Location','best','Box','off','Interpreter','none'); box(axA,'off');
        axB = subplot(2,1,2);  hold(axB,'on');  hB = [];  nB = {};
        for fi = 2:nF
            hB(end+1) = plot(axB, xv, ratio(fi,:), 'o-', 'Color', colF(fi,:), 'MarkerFaceColor', colF(fi,:), 'LineWidth', 1.3); %#ok<SAGROW>
            nB{end+1} = [file_labels{fi} ' / ' file_labels{1}]; %#ok<SAGROW>
        end
        if havePct
            hB(end+1) = plot(axB, xv, pct, 's--', 'Color', [0 0 0], 'LineWidth', 1.2); %#ok<SAGROW>
            nB{end+1} = 'percent set for the level'; %#ok<SAGROW>
        end
        yline(axB, 100, ':', 'Color', [0.5 0.5 0.5]);
        ylabel(axB, 'count, % of the first file'); xlabel(axB, xl);
        title(axB, 'ratio of the counts (summed over trials and channels) vs the percentage set');
        if ~isempty(hB), legend(axB, hB, nB, 'Location','best','Box','off','Interpreter','none'); end
        box(axB,'off');
        if strcmp(xAxisPulse, 'pulse'), set(axA, 'XTick', xv); set(axB, 'XTick', xv); end
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
