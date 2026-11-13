%% ============================================================
% NP_SP_RasterPSTH_v1.m
% NP probe | after NP_2_SpikeDetect (or a later recovered spike file) | RASTER + PSTH | v1   (read-only: nothing is saved)
%
% PURPOSE
%   Raster and PSTH of chosen channels and conditions, to look at the detected spikes.
%   conds = []  -> only PRINTS the condition list of the session (number, label, kind, level, trials) and stops.
%
% FIGURES (displayed only, not saved)
%   One figure per CHANNEL and PATTERN:  'sim'  |  'seq' (A->B or B->A)  |  'blank + visual'
%   Inside a figure: one column per chosen condition (raster on top, PSTH below); after maxCols conditions the
%   panels continue on a new row. All PSTHs of a figure share the same y axis.
%   Raster : one row per trial (repetition 1 at the bottom); grey = blanked stretch of that trial and channel;
%            red ticks above = stimulation pulses; blue bar above = visual spot on (photodiode).
%   PSTH   : spikes/s; each bin is divided by the number of trials that had data in it (blanked parts do not count),
%            bins where less than minDataFrac of the trials had data are left empty (gap, not zero).
%            Optional Gaussian smoothing that does not bridge the gaps.
%
% ALIGNMENT
%   electrical and blank trials: the trigger (alignElec = 'firstPulse' -> the first pulse, ~0.25 ms later)
%   visual trials             : the photodiode onset (alignVisual = 'trigger' -> the trigger)
%
% INPUTS (recording1 folder of the session): <short>.NP_TrialTable.mat and the spike file (default <short>.NP_Spikes.mat)
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile  = '';          % '' = <short>.NP_Spikes.mat; or a file name in recording1 / a full path (e.g. the recovered file later)
conds      = [10:15 40:45 49 52];          % conditions to plot (TT.cond numbers); [] = print the condition list and stop
channels   = [];          % AP channels (1-based); [] = automatic choice (see below)

% automatic channel choice (only used when channels = [])
autoN       = 8;          % channels spread evenly along the ACTIVE part of the probe
autoRateMin = 20;         % spikes/s: a channel counts as active above this rate

% time axis and alignment
winMs       = [-20 200];      % ms around the alignment point
alignElec   = 'trigger';      % 'trigger' | 'firstPulse'   (electrical and blank trials)
alignVisual = 'photodiode';   % 'photodiode' | 'trigger'   (visual trials)

% PSTH
binMs         = 1;        % bin width (ms)
smoothSigmaMs = 2;        % Gaussian smoothing SD (ms); 0 = no smoothing
minDataFrac   = 0.5;      % a bin is drawn only if at least this fraction of the trials had (unblanked) data in it

% layout
maxCols    = 5;           % conditions per row of panels; more conditions continue on the next row(s)
dotSize    = 3;           % raster dot size (points)
%% =================================================

%% ---- files -----------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
nTrAll = numel(TT.trial.trigAP);

%% ---- condition list only ------------------------------------------------------------------------------
if isempty(conds)
    fprintf('\nConditions of %s  (set conds = [...] to plot):\n', short);
    fprintf('  cond  kind     level  final  trials  label\n');
    for c = 1:numel(TT.cond)
        C = TT.cond(c);
        fprintf('  %4d  %-7s  %5d  %5d  %6d  %s\n', c, C.kind, C.level, C.isFinal, C.nTrials, C.label);
    end
    return
end

%% ---- spike file -------------------------------------------------------------------------------------------
if isempty(spikeFile), spikeFile = [short '.NP_Spikes.mat']; end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
assert(isfile(spikeFile), 'Spike file not found: %s', spikeFile);
SP = load(spikeFile, 'sp_clipped', 'sp_trial', 'blank_ms', 'nSpk', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed');
fs = SP.FS;  nCh = numel(SP.sp_clipped);  nTr = SP.nTrialsProcessed;
assert(all(conds >= 1 & conds <= numel(TT.cond)), 'conds must be between 1 and %d.', numel(TT.cond));
fprintf('\nspike file: %s\n', spikeFile);

tr     = TT.trial;
trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;               % trigger of each trial, ms from the AP file start
kinds  = {TT.cond.kind};
y      = TT.chan.y_um(:);

% alignment point of every trial, in ms AFTER ITS TRIGGER
off = zeros(nTr, 1);
for k = 1:nTr
    if strcmp(kinds{tr.cond(k)}, 'visual')
        if strcmp(alignVisual, 'photodiode') && ~isnan(tr.diodeOnMs(k)), off(k) = tr.diodeOnMs(k); end
    elseif strcmp(alignElec, 'firstPulse') && ~isempty(tr.pulseMs{k})
        off(k) = tr.pulseMs{k}(1);
    end
end
if winMs(1) + min(off) < -SP.preMs || winMs(2) + max(off) > SP.postMs
    warning('NPRP:win', 'Part of the window lies outside the detection window (-%g..+%g ms): those parts have no spikes.', SP.preMs, SP.postMs);
end

%% ---- channels ------------------------------------------------------------------------------------------
rate = sum(SP.nSpk(:, 1:nTr), 2)/nTr/((SP.preMs + SP.postMs)/1e3);
if isempty(channels)
    o   = TT.chan.orderTipToTop(:);  o = o(SP.goodCh(o));
    act = o(rate(o) >= autoRateMin);
    assert(~isempty(act), 'No channel above autoRateMin = %g spikes/s.', autoRateMin);
    channels = act(unique(round(linspace(1, numel(act), min(autoN, numel(act))))));
end
[~, q] = sort(y(channels), 'descend');  channels = channels(q);
fprintf('channels (top -> tip): %s\n', mat2str(channels(:)'));

%% ---- patterns of the chosen conditions ---------------------------------------------------------------
grpOf = cell(size(conds));
for i = 1:numel(conds)
    switch kinds{conds(i)}
        case 'sim',              grpOf{i} = 'sim';
        case {'AB', 'BA'},       grpOf{i} = 'seq';
        case {'blank', 'visual'}, grpOf{i} = 'blank + visual';
        otherwise,               grpOf{i} = kinds{conds(i)};
    end
end
grpList = unique(grpOf, 'stable');

%% ---- PSTH bins and smoothing kernel ---------------------------------------------------------------------
edges = winMs(1):binMs:winMs(2);  ctr = edges(1:end-1) + binMs/2;  nB = numel(ctr);
if smoothSigmaMs > 0
    sb = smoothSigmaMs/binMs;  kx = -ceil(3*sb):ceil(3*sb);  ker = exp(-0.5*(kx/sb).^2);  ker = ker/sum(ker);
else
    ker = 1;
end

%% ---- figures ------------------------------------------------------------------------------------------------
for ch = channels(:)'
    % spikes of this channel, ms after each trial's trigger, grouped by trial
    if isempty(SP.sp_clipped{ch})
        relBy = cell(nTr, 1);
    else
        kk  = SP.sp_trial{ch}(:);  rel = SP.sp_clipped{ch}(:) - trigMs(kk);
        relBy = accumarray(kk, rel, [nTr 1], @(v) {v});
    end
    for gI = 1:numel(grpList)
        cs = conds(strcmp(grpOf, grpList{gI}));
        nC = numel(cs);  nCol = min(maxCols, nC);  nRB = ceil(nC/nCol);
        fh = figure('Name', sprintf('%s ch %d %s', short, ch, grpList{gI}), 'Color', 'w', ...
                    'Position', [40 40 min(1800, 300*nCol) min(1000, 420*nRB)]);
        axP = gobjects(nC, 1);  yMax = 0;
        for i = 1:nC
            c = cs(i);  ks = find(tr.cond(1:nTr) == c);  ks = ks(:)';  nK = numel(ks);
            row = ceil(i/nCol);  col = i - (row-1)*nCol;
            pR = [(3*row-3)*nCol + col, (3*row-2)*nCol + col];  pP = (3*row-1)*nCol + col;

            % ---- raster ----
            axR = subplot(3*nRB, nCol, pR);  hold(axR, 'on');
            bx = [];  by = [];  sx = [];  sy = [];
            for r = 1:nK
                k = ks(r);
                iv = SP.blank_ms{ch, k};
                for b = 1:size(iv, 1)
                    bx = [bx, [iv(b,1); iv(b,2); iv(b,2); iv(b,1)] - off(k)]; %#ok<AGROW>
                    by = [by, [r-0.5; r-0.5; r+0.5; r+0.5]];                  %#ok<AGROW>
                end
                t = relBy{k} - off(k);  t = t(t >= winMs(1) & t < winMs(2));
                sx = [sx; t(:)];  sy = [sy; r*ones(numel(t), 1)]; %#ok<AGROW>
            end
            if ~isempty(bx), patch(axR, 'XData', bx, 'YData', by, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); end
            plot(axR, sx, sy, 'k.', 'MarkerSize', dotSize);
            if nK > 0
                k1 = ks(1);
                if ~isempty(tr.pulseMs{k1})
                    pt = tr.pulseMs{k1} - off(k1);
                    plot(axR, pt, (nK + 1.2)*ones(size(pt)), 'r|', 'MarkerSize', 4);
                end
                if strcmp(kinds{c}, 'visual')
                    on_ = median(tr.diodeOnMs(ks) - off(ks), 'omitnan');  of_ = median(tr.diodeOffMs(ks) - off(ks), 'omitnan');
                    if ~isnan(on_), plot(axR, [on_ of_], (nK + 1.2)*[1 1], 'b-', 'LineWidth', 3); end
                end
            end
            xline(axR, 0, ':');  xlim(axR, winMs);  ylim(axR, [0.5 max(nK, 1) + 1.8]);
            set(axR, 'XTickLabel', []);  box(axR, 'on');
            if col == 1, ylabel(axR, 'trial'); end
            title(axR, sprintf('%d: %s (n=%d)', c, TT.cond(c).label, nK), 'FontSize', 8, 'Interpreter', 'none');

            % ---- PSTH (blank-aware) ----
            cnt = zeros(1, nB);  expo = zeros(1, nB);
            for r = 1:nK
                k = ks(r);
                t = relBy{k} - off(k);
                cnt = cnt + histcounts(t, edges);
                cov = ones(1, nB);                                % fraction of each bin with data
                iv = SP.blank_ms{ch, k};
                for b = 1:size(iv, 1)
                    lo = max(edges(1:end-1), iv(b,1) - off(k));  hi = min(edges(2:end), iv(b,2) - off(k));
                    cov = cov - max(0, hi - lo)/binMs;
                end
                expo = expo + max(cov, 0);
            end
            cntS = conv(cnt, ker, 'same');  expS = conv(expo, ker, 'same');
            rt = cntS ./ (expS*binMs/1e3);
            rt(expo < minDataFrac*max(nK, 1) | expS <= 0) = NaN;
            axP(i) = subplot(3*nRB, nCol, pP);  hold(axP(i), 'on');
            plot(axP(i), ctr, rt, 'k', 'LineWidth', 1.2);
            xline(axP(i), 0, ':');  xlim(axP(i), winMs);  box(axP(i), 'on');
            yMax = max([yMax, rt(~isnan(rt))]);
            if col == 1, ylabel(axP(i), 'spikes/s'); end
            if row == nRB || i + nCol > nC, xlabel(axP(i), 'ms'); end
        end
        if yMax > 0, set(axP, 'YLim', [0 1.1*yMax]); end
        alignTxt = sprintf('aligned to: electrical/blank = %s, visual = %s', alignElec, alignVisual);
        sgtitle(fh, sprintf('%s  ch %d (y %d um from the tip)  %s\n%s; bin %g ms, smoothing %g ms; grey = blanked', ...
                strrep(short, '_', '\_'), ch, y(ch), grpList{gI}, alignTxt, binMs, smoothSigmaMs), 'FontSize', 10);
    end
end
fprintf('%d figure(s) drawn.\n', numel(channels)*numel(grpList));
