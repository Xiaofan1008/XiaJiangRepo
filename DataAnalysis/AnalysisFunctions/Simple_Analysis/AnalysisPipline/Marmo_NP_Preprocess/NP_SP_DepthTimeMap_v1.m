%% ============================================================
% NP_SP_DepthTimeMap_v1.m
% NP probe | any NP spike file (NP_Spikes / NP_Spikes_TrainRecovery...) | DEPTH x TIME FIRING-RATE MAP | v1
% READ-ONLY: figures only, nothing is saved.
%
% PURPOSE
%   Where along the probe, and when, does a condition drive activity? One figure per condition:
%     x = time (ms) after the trigger (visual conditions: after the photodiode onset, alignVisual)
%     y = depth from the probe TIP (um), every channel at its real position (settings.xml), shallow at the top
%     colour = firing rate per channel (blank-aware PSTH, as NP_SP_RasterPSTH), measure set by colorMode
%   + a side panel: each channel's mean value in profileWinMs (depth profile).
%   Pulse times are marked above the map; for visual conditions a bar shows when the spot was on.
%   conds = [] -> prints the condition list and stops.
%
% colorMode
%   'baseSub' (default) : rate - the channel's baseline rate (baseWinMs, same trials, relative to the TRIGGER)
%   'rate'              : spikes/s
%   'z'                 : (rate - baseline) / SD of the channel's PSTH in the baseline window
% Grey = no data (channel skipped in step 2, or a bin blanked in more than (1 - minDataFrac) of the trials).
%
% INPUTS (recording1 folder): <short>.NP_TrialTable.mat, the spike file (default: the recovered file if it exists,
%         else <short>.NP_Spikes.mat)
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile  = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery_randomEvent.mat';            % '' = <short>.NP_Spikes_TrainRecovery.mat if it exists, else <short>.NP_Spikes.mat
conds      = [15,45,52];            % conditions, one figure each (TT.cond numbers); [] = print the list and stop

winMs       = [-50 250];    % time axis (ms) around the alignment point
alignVisual = 'photodiode'; % visual trials: 'photodiode' | 'trigger'  (electrical / blank: always the trigger)
binMs       = 1;            % PSTH bin (ms)
smoothSigmaMs = 3;          % Gaussian smoothing SD along time (ms); 0 = none
minDataFrac = 0.5;          % a bin is shown only if at least this fraction of the trials had data in it

colorMode   = 'rate';    % 'baseSub' | 'rate' | 'z'   (see header)
baseWinMs   = [-300 -5];    % baseline window, ms relative to the TRIGGER (inside the detection window of step 2)
depthBinCh  = 1;            % channels averaged per row (1 = every channel; 4 = 80 um rows)
yRange      = [];           % [min max] depth from the tip (um) to show; [] = the whole probe
sharedClim  = true;         % true = one colour scale for all figures of this run; false = one per figure
clim        = [];           % fixed colour limits, e.g. [-100 100]; [] = automatic (99th percentile)
groupLines  = [];           % depths (um) drawn as horizontal lines, e.g. [4800] for a group / area boundary

showProfile  = true;        % side panel: mean value per channel in profileWinMs
profileWinMs = [5 150];     % ms (same alignment as the map)
%% =================================================

assert(any(strcmp(colorMode, {'baseSub', 'rate', 'z'})), 'colorMode must be ''baseSub'', ''rate'' or ''z''.');

%% ---- files -----------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L

if isempty(conds)
    fprintf('\nConditions of %s  (set conds = [...] to plot):\n  cond  kind     level  trials  label\n', short);
    for c = 1:numel(TT.cond)
        fprintf('  %4d  %-7s  %5d  %6d  %s\n', c, TT.cond(c).kind, TT.cond(c).level, TT.cond(c).nTrials, TT.cond(c).label);
    end
    return
end

if isempty(spikeFile)
    spikeFile = fullfile(recDir, [short '.NP_Spikes_TrainRecovery.mat']);
    if ~isfile(spikeFile), spikeFile = fullfile(recDir, [short '.NP_Spikes.mat']); end
end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
assert(isfile(spikeFile), 'Spike file not found: %s', spikeFile);
SP = load(spikeFile, 'sp_clipped', 'sp_trial', 'blank_ms', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed');
fs = SP.FS;  nCh = numel(SP.sp_clipped);  nTr = SP.nTrialsProcessed;
[~, spName] = fileparts(spikeFile);
assert(all(conds >= 1 & conds <= numel(TT.cond)), 'conds must be between 1 and %d.', numel(TT.cond));
if baseWinMs(1) < -SP.preMs, warning('NPDT:base', 'baseWinMs starts before the detection window (-%g ms).', SP.preMs); end
fprintf('\nspike file: %s\n', spikeFile);

tr = TT.trial;  trialCond = tr.cond(1:nTr);  kinds = {TT.cond.kind};
trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;
off = zeros(nTr, 1);                                         % alignment point, ms after each trial's trigger
if strcmp(alignVisual, 'photodiode')
    for k = 1:nTr
        if strcmp(kinds{trialCond(k)}, 'visual') && ~isnan(tr.diodeOnMs(k)), off(k) = tr.diodeOnMs(k); end
    end
end

%% ---- rows: channels sorted by depth, optionally binned -------------------------------------------------
y = TT.chan.y_um(:);
o = TT.chan.orderTipToTop(:);
if ~isempty(yRange), o = o(y(o) >= yRange(1) & y(o) <= yRange(2)); end
nRow = ceil(numel(o)/depthBinCh);
rowCh = cell(nRow, 1);  rowY = nan(nRow, 1);
for r = 1:nRow
    rowCh{r} = o((r-1)*depthBinCh + 1 : min(numel(o), r*depthBinCh));
    rowY(r) = mean(y(rowCh{r}));
end

%% ---- spikes per trial (relative to the trigger) for every channel ----------------------------------------
fprintf('sorting spikes by trial ...\n');
Bt = cell(nCh, 1);
for ch = 1:nCh, Bt{ch} = byTrial_(SP.sp_clipped{ch}, SP.sp_trial{ch}, trigMs, nTr); end

edges = winMs(1):binMs:winMs(2);  ctr = edges(1:end-1) + binMs/2;
bEdges = baseWinMs(1):binMs:baseWinMs(2);
if smoothSigmaMs > 0
    sb = smoothSigmaMs/binMs;  kx = -ceil(3*sb):ceil(3*sb);  ker = exp(-0.5*(kx/sb).^2);  ker = ker/sum(ker);
else
    ker = 1;
end

%% ---- compute every condition first (for the shared colour scale) -------------------------------------------
fprintf('computing %d condition(s) ...\n', numel(conds));
Mall = cell(numel(conds), 1);
for i = 1:numel(conds)
    c = conds(i);  ks = find(trialCond == c);  ks = ks(:)';
    Mch = nan(nCh, numel(ctr));
    for ch = find(SP.goodCh(:))'
        rt = psth_(Bt{ch}, SP.blank_ms(ch, :), ks, off, edges, ker, minDataFrac);       % aligned
        if ~strcmp(colorMode, 'rate')
            bl = psth_(Bt{ch}, SP.blank_ms(ch, :), ks, zeros(nTr, 1), bEdges, ker, minDataFrac);   % baseline, trigger-aligned
            mB = mean(bl, 'omitnan');
            if strcmp(colorMode, 'baseSub')
                rt = rt - mB;
            else
                sB = std(bl, 'omitnan');  if ~(sB > 0), sB = NaN; end
                rt = (rt - mB)/sB;
            end
        end
        Mch(ch, :) = rt;
    end
    M = nan(nRow, numel(ctr));
    for r = 1:nRow, M(r, :) = mean(Mch(rowCh{r}, :), 1, 'omitnan'); end
    Mall{i} = M;
end
if isempty(clim)
    v = cell2mat(cellfun(@(m) m(:), Mall, 'UniformOutput', false));  v = v(~isnan(v));
    if strcmp(colorMode, 'rate')
        cl = [0, pct_(v, 99)];
    else
        a = pct_(abs(v), 99);  cl = [-a a];
    end
    if any(~isfinite(cl)) || diff(cl) <= 0, cl = [-1 1]; end
else
    cl = clim;
end
switch colorMode
    case 'rate',    cmap = parula(256);  cbl = 'spikes/s';
    case 'baseSub', cmap = redblue_();   cbl = 'rate - baseline (spikes/s)';
    case 'z',       cmap = redblue_();   cbl = 'z (baseline SD)';
end

%% ---- figures -----------------------------------------------------------------------------------------------------
for i = 1:numel(conds)
    c = conds(i);  ks = find(trialCond == c);  M = Mall{i};
    if ~sharedClim && isempty(clim)
        v = M(~isnan(M));
        if strcmp(colorMode, 'rate'), cl_ = [0, pct_(v, 99)]; else, a = pct_(abs(v), 99); cl_ = [-a a]; end
        if any(~isfinite(cl_)) || diff(cl_) <= 0, cl_ = [-1 1]; end
    else
        cl_ = cl;
    end
    fh = figure('Name', sprintf('%s depth x time cond %d', short, c), 'Color', 'w', 'Position', [60 60 1100 750]);
    if showProfile, ax = axes(fh, 'Position', [0.08 0.10 0.62 0.78]); else, ax = axes(fh, 'Position', [0.08 0.10 0.80 0.78]); end
    hI = imagesc(ax, ctr, rowY, M);  set(hI, 'AlphaData', ~isnan(M));
    set(ax, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'CLim', cl_, 'TickDir', 'out');  colormap(ax, cmap);  hold(ax, 'on');
    cb = colorbar(ax, 'Location', 'westoutside');  cb.Label.String = cbl;
    xlabel(ax, sprintf('ms after the %s', ternary_(strcmp(kinds{c}, 'visual') && strcmp(alignVisual, 'photodiode'), 'photodiode onset', 'trigger')));
    ylabel(ax, 'depth from the tip (um)');
    yl = [min(rowY) - 10, max(rowY) + 10];
    xline(ax, 0, 'k:');
    for g = groupLines(:)', yline(ax, g, 'k-', 'LineWidth', 1); end
    if ~isempty(ks)
        pt = tr.pulseMs{ks(1)} - off(ks(1));
        top = yl(2) + 0.02*diff(yl);
        if ~isempty(pt), plot(ax, [pt(:)'; pt(:)'], repmat([yl(2); top + 0.03*diff(yl)], 1, numel(pt)), 'r-', 'LineWidth', 0.8, 'Clipping', 'off'); end
        if strcmp(kinds{c}, 'visual')
            on_ = median(tr.diodeOnMs(ks) - off(ks), 'omitnan');  of_ = median(tr.diodeOffMs(ks) - off(ks), 'omitnan');
            if ~isnan(on_), plot(ax, [on_ of_], [top top] + 0.015*diff(yl), 'b-', 'LineWidth', 4, 'Clipping', 'off'); end
        end
    end
    ylim(ax, yl);  xlim(ax, winMs);
    if showProfile
        axP = axes(fh, 'Position', [0.74 0.10 0.22 0.78]);  hold(axP, 'on');
        sel = ctr >= profileWinMs(1) & ctr < profileWinMs(2);
        prof = mean(M(:, sel), 2, 'omitnan');
        plot(axP, prof, rowY, 'k-', 'LineWidth', 1);
        if ~strcmp(colorMode, 'rate'), xline(axP, 0, ':'); end
        for g = groupLines(:)', yline(axP, g, 'k-', 'LineWidth', 1); end
        ylim(axP, yl);  grid(axP, 'on');  set(axP, 'YTickLabel', []);
        xlabel(axP, sprintf('mean %g..%g ms', profileWinMs));
    end
    sgtitle(fh, sprintf('%s   cond %d: %s   (n = %d trials)\n%s; colour = %s; %d channel(s) per row; file %s', ...
            strrep(short, '_', '\_'), c, strrep(TT.cond(c).label, '_', '\_'), numel(ks), ...
            ternary_(strcmp(kinds{c}, 'visual') && strcmp(alignVisual, 'photodiode'), 'aligned to the photodiode onset', 'aligned to the trigger'), ...
            colorMode, depthBinCh, strrep(spName, '_', '\_')), 'FontSize', 10);
end
fprintf('%d figure(s) drawn; colour limits %s\n', numel(conds), mat2str(round(cl, 1)));

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function B = byTrial_(tAbs, k, trigMs, nTr)
B = cell(nTr, 1);  for t = 1:nTr, B{t} = zeros(0, 1); end
if isempty(k), return; end
k = k(:);  rel = tAbs(:) - trigMs(k);
[ks, o] = sort(k);  rel = rel(o);
e = [0; find(diff(ks)); numel(ks)];
for i = 1:numel(e) - 1
    rows = e(i)+1:e(i+1);  B{ks(rows(1))} = sort(rel(rows));
end
end

function rt = psth_(B, blankRow, ks, off, edges, ker, minFrac)
% blank-aware PSTH (spikes/s) of trials ks; times shifted by off(trial); bins with data in < minFrac of the trials = NaN
binMs = edges(2) - edges(1);  nB = numel(edges) - 1;
cnt = zeros(1, nB);  expo = zeros(1, nB);
for t = ks
    cnt = cnt + histcounts(B{t} - off(t), edges);
    covr = ones(1, nB);  iv = blankRow{t};
    for b = 1:size(iv, 1)
        lo = max(edges(1:end-1), iv(b,1) - off(t));  hi = min(edges(2:end), iv(b,2) - off(t));
        covr = covr - max(0, hi - lo)/binMs;
    end
    expo = expo + max(covr, 0);
end
cS = conv(cnt, ker, 'same');  eS = conv(expo, ker, 'same');
rt = cS ./ (eS*binMs/1e3);
rt(expo < minFrac*max(numel(ks), 1) | eS <= 0) = NaN;
end

function v = pct_(x, p)
x = sort(x(~isnan(x(:))));
if isempty(x), v = NaN; return; end
r = 1 + (numel(x) - 1)*p/100;  lo = floor(r);  hi = ceil(r);
v = x(lo) + (r - lo)*(x(hi) - x(lo));
end

function cm = redblue_()
n = 128;  a = linspace(0, 1, n)';
cm = [[a, a, ones(n,1)]; [ones(n,1), flipud(a), flipud(a)]];
cm(n, :) = [1 1 1];
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
