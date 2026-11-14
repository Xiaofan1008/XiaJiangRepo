%% ============================================================
% NP_Train_FR_RespChnDepthPlot_v1.m
% NP probe | pulse TRAIN sessions | FIRING-RATE analysis | RESPONSIVE-CHANNEL DEPTH PLOT | v1
%
% PURPOSE
%   Plots the SAVED responsive channels (<short>.NP_RespChn.mat, from NP_Train_FR_RespChnEntry_v1) along the probe:
%   one panel per pattern, depth from the tip vs the check's z of the decision condition. Only plots and prints:
%   NOTHING IS SAVED, no file is changed. Needs only the NP_RespChn file (saved by the entry version that stores the
%   numbers: RespChn.pattern(p).score / tested; re-save with the entry if they are missing).
%
% SAME X AXIS FOR ALL PANELS (to compare patterns)
%   xMax = 'auto': the upper limit is set so that at most nOutMax points (all panels together) lie beyond it, i.e. a few
%          extreme channels do not squeeze the rest; rounded up to a multiple of xStep. If there are no more than nOutMax
%          points in total, the largest value is used.
%   xMax = a number: fixed upper limit.
%   Points beyond the limit are drawn ON the limit as OPEN markers, with their value written next to them.
%
% MARKERS
%   green dot    = final responsive set         grey dot       = tested for this pattern, not chosen
%   grey bar     = not tested for this pattern (no z; drawn at the left edge)   dashed = the check's threshold
%   (channels added by hand are simply part of the final set; channels removed by hand count as 'not chosen')
%
% Y RANGE: yRange = [from to] (um from the tip) shows only that part of the probe; the automatic x limit then uses only
%   the channels inside it. yTickStep sets the spacing of the depth ticks.
%
% PRINTED: per pattern the final count, counts per depth band (final of tested), and the channels beyond the x limit.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
respFile   = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_RespChn.mat';            % '' = <short>.NP_RespChn.mat in the recording folder (or a full path)
patterns   = 'all';         % 'all' | names, e.g. {'sim', 'seq B->A', 'visual bright', 'visual dark'}

% ---- x axis (shared by all panels) ----
xMax       = 'auto';        % 'auto' | a number (fixed upper limit of z)
nOutMax    = 2;             % auto: at most this many points (all panels together) may lie beyond the limit
xStep      = 5;             % auto: the limit is rounded up to a multiple of this

% ---- what to draw ----
showTested    = true;       % grey dots: tested for the pattern, not in the final set
showUntested  = true;       % grey bar at the left edge: processed but not tested for the pattern
showThreshold = true;       % dashed line at the check's threshold (z = 3)
showOutValues = true;       % write the value next to the points beyond the x limit
groupLines    = [2200 3000 6000];   % depth lines (um) and bands for the printed counts; [] = none
yRange        = [];         % [] = the whole probe | [from to] in um from the tip, e.g. [0 6000]
yTickStep     = 1000;       % spacing of the depth ticks (um); [] = MATLAB's choice

% ---- style ----
colFinal   = [0 0.6 0];     % final set
colTested  = [0.6 0.6 0.6]; % tested / not tested
markerSize = 4;             % size of the dots
fontSize   = 9;
figPos     = [40 40 1500 800];
%% =================================================

%% ---- load ---------------------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
if isempty(respFile), respFile = [short '.NP_RespChn.mat']; end
if ~isfile(respFile), respFile = fullfile(recDir, respFile); end
assert(isfile(respFile), 'Responsive-channel file not found: %s (save it with NP_Train_FR_RespChnEntry_v1).', respFile);
L = load(respFile, 'RespChn');  RC = L.RespChn;  clear L
assert(isfield(RC.pattern, 'score') && isfield(RC.pattern, 'tested'), ...
       ['%s has no z values / tested lists (saved by an older entry version).\n' ...
        'Run NP_Train_FR_RespChnEntry_v1 again with doSave = true.'], respFile);
y = double(RC.y_um(:));  nCh = numel(y);  chProc = reshape(RC.chTest, 1, []);   % processed channels
thr = 3;  if isfield(RC, 'scoreThr'), thr = RC.scoreThr; end
sName = 'z (response vs baseline)';  if isfield(RC, 'scoreName'), sName = RC.scoreName; end

% patterns to plot
P = RC.pattern;
if ischar(patterns) && strcmpi(patterns, 'all')
    sel = 1:numel(P);
else
    if ischar(patterns), patterns = {patterns}; end
    sel = zeros(1, numel(patterns));
    for i = 1:numel(patterns)
        k = find(strcmpi(strtrim({P.name}), strtrim(patterns{i})), 1);
        assert(~isempty(k), 'Pattern "%s" not in the file. Patterns: %s', patterns{i}, strjoin({P.name}, ', '));
        sel(i) = k;
    end
end
nP = numel(sel);

%% ---- y range and shared x limits -------------------------------------------------------------------------------------
if isempty(yRange), yl = [min(y) max(y)]; else, yl = sort(yRange); end
inY = @(v) v(y(v) >= yl(1) & y(v) <= yl(2));                      % channels inside the y range
S = cell(1, nP);  vals = [];
for i = 1:nP
    p = P(sel(i));  s = nan(nCh, 1);  sc = double(p.score(:));  s(1:min(nCh, numel(sc))) = sc(1:min(nCh, numel(sc)));
    S{i} = s;
    shown = p.final;  if showTested, shown = union(shown, p.tested); end
    v = s(inY(shown));  vals = [vals; v(isfinite(v))]; %#ok<AGROW>
end
if isempty(vals), vals = [0; thr]; end
if ischar(xMax) && strcmpi(xMax, 'auto')
    v = sort(vals, 'descend');
    if numel(v) > nOutMax, xHi = v(nOutMax + 1); else, xHi = v(1); end
    xHi = ceil(max(xHi, thr)/xStep)*xStep;
else
    xHi = xMax;
end
xLo = min([vals; 0]);
xl  = [xLo - 0.08*(xHi - xLo), xHi + 0.03*(xHi - xLo)];
xEdge = xl(1) + 0.012*diff(xl);                                   % x of the 'not tested' bar

%% ---- figure ----------------------------------------------------------------------------------------------------------
fh = figure('Name', [short ' responsive channels along the probe'], 'Color', 'w', 'Position', figPos);
gL = 0.06;  gR = 0.015;  gap = 0.025;  bot = 0.13;  hgt = 0.77;
w = (1 - gL - gR - (nP - 1)*gap)/nP;
outTxt = cell(1, nP);
for i = 1:nP
    p = P(sel(i));  s = S{i};
    ax = axes(fh, 'Position', [gL + (i - 1)*(w + gap), bot, w, hgt]);  hold(ax, 'on');
    for g = groupLines(:)', plot(ax, xl, [g g], '-', 'Color', [0.75 0.75 0.75]); end
    if showThreshold, plot(ax, [thr thr], yl, '--', 'Color', [0.2 0.8 0.2], 'LineWidth', 1); end
    % groups of channels
    fin  = p.final;
    tst  = setdiff(p.tested, fin);
    unt  = setdiff(chProc, union(p.tested, fin));
    if showUntested, plot(ax, repmat(xEdge, 1, numel(unt)), y(unt), '.', 'Color', colTested, 'MarkerSize', markerSize*2); end
    if showTested, dots_(ax, s, y, tst, xHi, xEdge, '.', colTested, markerSize*2.5, false); end
    dots_(ax, s, y, fin, xHi, xEdge, 'o', colFinal, markerSize, true);
    % values beyond the limit
    shown = fin;  if showTested, shown = union(shown, tst); end
    shown = inY(shown);  out = shown(s(shown) > xHi);  outTxt{i} = out;
    if showOutValues
        for ch = out(:)'
            text(ax, xHi - 0.01*diff(xl), y(ch), sprintf('%.0f', s(ch)), 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'middle', 'FontSize', fontSize - 2, 'Color', [0.3 0.3 0.3]);
        end
    end
    xlim(ax, xl);  ylim(ax, yl);  grid(ax, 'on');  set(ax, 'FontSize', fontSize, 'Layer', 'top');
    if ~isempty(yTickStep), set(ax, 'YTick', ceil(yl(1)/yTickStep)*yTickStep : yTickStep : yl(2)); end
    xlabel(ax, sName, 'FontSize', fontSize);
    if i == 1, ylabel(ax, 'depth from the tip (\mum)'); else, set(ax, 'YTickLabel', []); end
    title(ax, sprintf('%s\n%d of %d tested', p.name, numel(fin), numel(p.tested)), 'FontSize', fontSize + 1, 'Interpreter', 'none');
    if i == 1, ax1 = ax; end
end
% legend from dummy lines (always complete), centred below the panels
hL = plot(ax1, NaN, NaN, 'o', 'Color', colFinal, 'MarkerFaceColor', colFinal, 'MarkerSize', markerSize);  lab = {'final set'};
if showTested,   hL(end+1) = plot(ax1, NaN, NaN, '.', 'Color', colTested, 'MarkerSize', markerSize*2.5);  lab{end+1} = 'tested, not chosen'; end
if showUntested, hL(end+1) = plot(ax1, NaN, NaN, '.', 'Color', colTested, 'MarkerSize', markerSize*2);    lab{end+1} = 'not tested (left edge)'; end
hL(end+1) = plot(ax1, NaN, NaN, 'o', 'Color', colFinal, 'MarkerSize', markerSize + 1);  lab{end+1} = sprintf('beyond the x limit (%g)', xHi);
if showThreshold, hL(end+1) = plot(ax1, NaN, NaN, '--', 'Color', [0.2 0.8 0.2], 'LineWidth', 1);  lab{end+1} = sprintf('threshold (%g)', thr); end
lg = legend(ax1, hL, lab, 'Orientation', 'horizontal', 'FontSize', fontSize, 'Box', 'off');
lg.Units = 'normalized';  lp = lg.Position;  lg.Position = [0.5 - lp(3)/2, 0.01, lp(3), lp(4)];
sgtitle(fh, sprintf('%s: responsive channels along the probe (saved %s)', strrep(short, '_', '\_'), RC.date), 'FontSize', fontSize + 3);

%% ---- printout ----------------------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP RESPONSIVE-CHANNEL DEPTH PLOT (v1): %s\n', short);
fprintf('============================================================\n');
fprintf('file    : %s (saved %s)\ny range : %g-%g um\nx limit : %s (%s)\n', respFile, RC.date, yl(1), yl(2), num2str(xHi), ...
        ternary_(ischar(xMax), sprintf('auto, at most %d points beyond', nOutMax), 'fixed'));
bandE = [-Inf, groupLines(:)', Inf];
for i = 1:nP
    p = P(sel(i));
    fprintf('\n--- %s (conds %s, decision cond %d): %d of %d tested ---\n', p.name, mat2str(p.conds), p.decCond, numel(p.final), numel(p.tested));
    for b = numel(bandE) - 1:-1:1
        inB = @(v) v(y(v) >= bandE(b) & y(v) < bandE(b + 1));
        fprintf('    depth %5s-%5s um: %3d of %3d tested\n', num2s_(bandE(b)), num2s_(bandE(b + 1)), numel(inB(p.final)), numel(inB(p.tested)));
    end
    if ~isempty(outTxt{i})
        fprintf('  beyond the x limit: %s\n', strjoin(arrayfun(@(ch) sprintf('ch %d (z %.1f)', ch, S{i}(ch)), outTxt{i}, 'UniformOutput', false), ', '));
    end
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function dots_(ax, s, y, chs, xHi, xEdge, mk, col, ms, filled)
% channels chs at x = s (NaN -> left edge; beyond xHi -> on xHi as an open marker)
if isempty(chs), return; end
x = s(chs);  x = x(:);  yy = y(chs);  yy = yy(:);
x(isnan(x) | x == -Inf) = xEdge;
out = x > xHi;  x(out) = xHi;
if filled, fc = col; else, fc = 'none'; end
if any(~out), plot(ax, x(~out), yy(~out), mk, 'Color', col, 'MarkerFaceColor', fc, 'MarkerSize', ms, 'LineWidth', 1); end
if any(out),  plot(ax, x(out), yy(out), mk, 'Color', col, 'MarkerFaceColor', 'none', 'MarkerSize', ms + 1, 'LineWidth', 1); end
end

function s = num2s_(v)
if isnan(v), s = '-'; elseif isinf(v), s = ternary_(v > 0, 'top', 'tip'); elseif v == round(v), s = sprintf('%d', v); else, s = sprintf('%.1f', v); end
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
