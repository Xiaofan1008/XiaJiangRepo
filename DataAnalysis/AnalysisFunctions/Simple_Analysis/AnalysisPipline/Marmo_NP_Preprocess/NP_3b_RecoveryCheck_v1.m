%% ============================================================
% NP_3b_RecoveryCheck_v1.m
% NP probe | after NP_3_SpikeRecovery | RECOVERY CHECK: automatic test + count table + figures | v1
% READ-ONLY: nothing is saved or modified.   (NP versions of NN_RecoveryCheck_v1 + NN_RecoveryCountCheck_v1)
%
% PART 1  Automatic check (all channels, all recovered trials; runs whenever runPart1 = true):
%   OUTSIDE the replaced stretches (rec_ms) the recovered spikes must be IDENTICAL to the original ones.
%   INSIDE  the replaced stretch the recovered spikes must equal (count AND times) those of the paired SOURCE trial
%           in the same stretch (recovered data).
%   With simMode 'randomEvent' the thinning window [last pulse, last pulse + event_ms) of each thinned level is left out
%   of both tests (spikes there were removed on purpose).
% PART 2  Count table (conditions in conds, channels in channels): spikes per trial and channel
%   right after each pulse  [pulse, pulse + gap_ms)   and in the rest of the interval  [pulse + gap_ms, next pulse),
%   in the ORIGINAL data (blanked -> ~0 right after the pulse), the RECOVERED data and the SOURCE level (recovered).
% PART 3  Figures (conds x channels), one figure per channel and pattern (sim / seq):
%   raster : black = original spikes, orange = copied spikes (inside rec_ms), grey = still without data,
%            red ticks = pulses
%   PSTH   : black = recovered, grey = original (gaps where blanked), dashed blue = SOURCE level (recovered).
%            Inside the replaced stretch the black line should follow the dashed one; after it, black = grey.
%   showLevelOverlay: one more figure per channel and pattern with the recovered PSTHs of the chosen levels overlaid
%            (light = low level, dark = high level): they should lie on top of each other up to each level's last pulse.
% conds = []  -> prints the condition list and runs Part 1 only.
%
% INPUTS (recording1 folder): <short>.NP_TrialTable.mat, <short>.NP_Spikes.mat, the recovered file (recFile)
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
recFile    = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery.mat';          % '' = <short>.NP_Spikes_TrainRecovery.mat; or a file name in recording1 / a full path
conds      = [10 14 15 30 40 45];          % conditions for Part 2 and 3 (TT.cond numbers); [] = print the list, run Part 1 only
channels   = [];          % AP channels for Part 2 and 3; [] = automatic (autoN active channels along the probe)
autoN       = 8;          % automatic choice: channels spread along the ACTIVE part of the probe
autoRateMin = 20;         % spikes/s: a channel counts as active above this rate
runPart1   = true;        % false = skip the automatic check (it takes a little while)

gap_ms     = 1;           % Part 2: 'right after the pulse' = [pulse, pulse + gap_ms)

winMs         = [-20 200];% Part 3: ms around the trigger
binMs         = 1;        % PSTH bin (ms)
smoothSigmaMs = 3;        % PSTH Gaussian smoothing SD (ms); 0 = none
minDataFrac   = 0.5;      % a PSTH bin is drawn only if at least this fraction of the trials had data in it
maxCols       = 5;        % conditions per row of panels
dotSize       = 3;        % raster dot size
showLevelOverlay = true;  % extra figure per channel and pattern: recovered PSTHs of the chosen levels overlaid
%% =================================================

%% ---- load ----------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
if isempty(recFile), recFile = [short '.NP_Spikes_TrainRecovery.mat']; end
if ~isfile(recFile), recFile = fullfile(recDir, recFile); end
assert(isfile(recFile), 'Recovered file not found: %s', recFile);
RC = load(recFile);
oFile = RC.spikeFile;  if ~isfile(oFile), oFile = fullfile(recDir, [short '.NP_Spikes.mat']); end
assert(isfile(oFile), 'Original spike file not found: %s', oFile);
SO = load(oFile, 'sp_clipped', 'sp_trial', 'blank_ms', 'nTrialsProcessed', 'FS');
fs = RC.FS;  nCh = numel(RC.sp_clipped);  nTr = RC.nTrialsProcessed;
hS = 0.5/fs*1e3;                                        % half a sample: replaced stretches include their first and last sample
assert(SO.nTrialsProcessed == nTr && numel(SO.sp_clipped) == nCh, 'Original and recovered files do not match.');
tr = TT.trial;  trialCond = tr.cond(1:nTr);
trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;
CT = RC.condTable;  depth = RC.depth;  RP = RC.Recovery_params;
isSim = isfield(RP, 'simMode') && strcmp(RP.simMode, 'randomEvent');
thinWin = nan(nTr, 2);                                  % [start end) of the thinning window of each trial (randomEvent)
if isSim
    for t = 1:nTr
        c = trialCond(t);
        if depth(c) > 0 && RP.percent(c) < 100, thinWin(t, :) = tr.lastPulseMs(t) + [0 RP.event_ms(c)]; end
    end
end
fprintf('\n============================================================\nNP RECOVERY CHECK: %s\n', short);
fprintf('original : %s\nrecovered: %s  (boundary %s, simulation %s)\n', oFile, recFile, RP.boundaryMode, RP.simMode);

%% ---- condition list only ------------------------------------------------------------------------------------
if isempty(conds)
    fprintf('\nConditions (depth 0 = not recovered):\n  cond  kind   level  source  depth  trials  label\n');
    for c = 1:numel(CT)
        fprintf('  %4d  %-5s  %5d  %6s  %5d  %6d  %s\n', c, CT(c).kind, CT(c).level, n2s_(CT(c).sourceCond), depth(c), ...
                sum(trialCond == c), CT(c).label);
    end
end

%% ---- PART 1: automatic check ---------------------------------------------------------------------------------
if runPart1
    fprintf('\nPART 1: automatic check over %d channels ...\n', nCh);
    tol = 1e-6;  nOut = 0;  badOut = 0;  nIn = 0;  badIn = 0;  badChOut = zeros(nCh, 1);  badChIn = zeros(nCh, 1);
    for ch = 1:nCh
        O = byTrial_(SO.sp_clipped{ch}, SO.sp_trial{ch}, trigMs, nTr);
        R = byTrial_(RC.sp_clipped{ch}, RC.sp_trial{ch}, trigMs, nTr);
        for t = 1:nTr
            iv = RC.rec_ms{ch, t};  ex = iv;
            if ~isnan(thinWin(t, 1)), ex = [ex; thinWin(t, :)]; end %#ok<AGROW>
            % outside: identical to the original
            a = O{t}(~inIv_(O{t}, ex, hS));  b = R{t}(~inIv_(R{t}, ex, hS));
            nOut = nOut + 1;
            if ~same_(a, b, tol), badOut = badOut + 1; badChOut(ch) = badChOut(ch) + 1; end
            % inside: identical to the source trial
            if RC.isRecovered(t) && ~isempty(iv)
                s = RC.srcTrial(t);
                x = R{t}(inIv_(R{t}, iv, hS));  y = R{s}(inIv_(R{s}, iv, hS));
                if ~isnan(thinWin(t, 1))                     % leave out this trial's thinning window
                    x = x(~inIv_(x, thinWin(t, :), 0));  y = y(~inIv_(y, thinWin(t, :), 0));
                end
                nIn = nIn + 1;
                if ~same_(x, y, tol), badIn = badIn + 1; badChIn(ch) = badChIn(ch) + 1; end
            end
        end
    end
    fprintf('  outside the replaced stretches = original : %d of %d trial x channel pairs %s\n', nOut - badOut, nOut, passTxt_(badOut));
    if badOut > 0, fprintf('     channels with mismatches: %s\n', mat2str(find(badChOut)')); end
    fprintf('  inside  the replaced stretch   = source   : %d of %d trial x channel pairs %s\n', nIn - badIn, nIn, passTxt_(badIn));
    if badIn > 0, fprintf('     channels with mismatches: %s\n', mat2str(find(badChIn)')); end
    if isSim, fprintf('  (randomEvent: the thinning windows were left out of both tests)\n'); end
end
if isempty(conds), return; end
assert(all(conds >= 1 & conds <= numel(CT)), 'conds must be between 1 and %d.', numel(CT));
skip = conds(depth(conds) == 0);
if ~isempty(skip), fprintf('\nconds %s are not recovered (level 1, blank, visual): left out of Part 2 and 3.\n', mat2str(skip)); end
conds = conds(depth(conds) > 0);
if isempty(conds), return; end

%% ---- channels ----------------------------------------------------------------------------------------------------
y = TT.chan.y_um(:);
if isempty(channels)
    rate = sum(RC.nSpk(:, 1:nTr), 2)/nTr/((RC.preMs + RC.postMs)/1e3);
    o = TT.chan.orderTipToTop(:);  o = o(RC.goodCh(o));
    act = o(rate(o) >= autoRateMin);
    assert(~isempty(act), 'No channel above autoRateMin = %g spikes/s.', autoRateMin);
    channels = act(unique(round(linspace(1, numel(act), min(autoN, numel(act))))));
end
[~, q] = sort(y(channels), 'descend');  channels = channels(q);
fprintf('\nchannels (top -> tip): %s\n', mat2str(channels(:)'));
Oc = cell(nCh, 1);  Rc = cell(nCh, 1);                    % per-trial spikes of the chosen channels
for ch = channels(:)'
    Oc{ch} = byTrial_(SO.sp_clipped{ch}, SO.sp_trial{ch}, trigMs, nTr);
    Rc{ch} = byTrial_(RC.sp_clipped{ch}, RC.sp_trial{ch}, trigMs, nTr);
end

%% ---- PART 2: count table -----------------------------------------------------------------------------------------
fprintf('\nPART 2: spikes per trial and channel (mean over the %d channels): right after each pulse [p, p+%g ms) | rest of the interval\n', ...
        numel(channels), gap_ms);
for c = conds(:)'
    ks = find(trialCond == c);  if isempty(ks), continue; end
    pl = tr.pulseMs{ks(1)};  sp = RP.spacing_ms(c);
    nxt = [pl(2:end), pl(end) + sp];
    fprintf('\n cond %d %s (source %d), %d trials\n', c, CT(c).label, CT(c).sourceCond, numel(ks));
    fprintf('   pulse  t(ms)   | after: original  recovered   source | rest: original  recovered   source\n');
    for p = 1:numel(pl)
        w1 = [pl(p), pl(p) + gap_ms];  w2 = [pl(p) + gap_ms, nxt(p)];
        v = zeros(1, 6);
        for ch = channels(:)'
            for t = ks(:)'
                s = RC.srcTrial(t);
                v = v + [cnt_(Oc{ch}{t}, w1), cnt_(Rc{ch}{t}, w1), cnt_(Rc{ch}{s}, w1), ...
                         cnt_(Oc{ch}{t}, w2), cnt_(Rc{ch}{t}, w2), cnt_(Rc{ch}{s}, w2)];
            end
        end
        v = v/(numel(ks)*numel(channels));
        fprintf('   %4d  %6.1f   |        %6.3f    %6.3f   %6.3f |       %6.3f    %6.3f   %6.3f\n', p, pl(p), v);
    end
end

%% ---- PART 3: figures ---------------------------------------------------------------------------------------------
edges = winMs(1):binMs:winMs(2);  ctr = edges(1:end-1) + binMs/2;
if smoothSigmaMs > 0
    sb = smoothSigmaMs/binMs;  kx = -ceil(3*sb):ceil(3*sb);  ker = exp(-0.5*(kx/sb).^2);  ker = ker/sum(ker);
else
    ker = 1;
end
grpOf = cell(size(conds));
for i = 1:numel(conds)
    if strcmp(CT(conds(i)).kind, 'sim'), grpOf{i} = 'sim'; else, grpOf{i} = 'seq'; end
end
grpList = unique(grpOf, 'stable');
orange = [0.95 0.5 0.05];
for ch = channels(:)'
    for gI = 1:numel(grpList)
        cs = conds(strcmp(grpOf, grpList{gI}));
        nC = numel(cs);  nCol = min(maxCols, nC);  nRB = ceil(nC/nCol);
        fh = figure('Name', sprintf('%s ch %d %s recovery', short, ch, grpList{gI}), 'Color', 'w', ...
                    'Position', [40 40 min(1800, 300*nCol) min(1000, 420*nRB)]);
        axP = gobjects(nC, 1);  yMax = 0;
        for i = 1:nC
            c = cs(i);  ks = find(trialCond == c);  ks = ks(:)';  nK = numel(ks);
            row = ceil(i/nCol);  col = i - (row-1)*nCol;
            pR = [(3*row-3)*nCol + col, (3*row-2)*nCol + col];  pP = (3*row-1)*nCol + col;
            % raster
            axR = subplot(3*nRB, nCol, pR);  hold(axR, 'on');
            bx = [];  by = [];  xo = [];  yo = [];  xc = [];  yc = [];
            for r = 1:nK
                t = ks(r);  iv = RC.blank_ms{ch, t};
                for b = 1:size(iv, 1)
                    bx = [bx, [iv(b,1); iv(b,2); iv(b,2); iv(b,1)]]; by = [by, [r-0.5; r-0.5; r+0.5; r+0.5]]; %#ok<AGROW>
                end
                s_ = Rc{ch}{t};  s_ = s_(s_ >= winMs(1) & s_ < winMs(2));
                cp = inIv_(s_, RC.rec_ms{ch, t}, hS);
                xo = [xo; s_(~cp)];  yo = [yo; r*ones(sum(~cp), 1)]; %#ok<AGROW>
                xc = [xc; s_(cp)];   yc = [yc; r*ones(sum(cp), 1)];  %#ok<AGROW>
            end
            if ~isempty(bx), patch(axR, 'XData', bx, 'YData', by, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); end
            plot(axR, xo, yo, 'k.', 'MarkerSize', dotSize);
            plot(axR, xc, yc, '.', 'Color', orange, 'MarkerSize', dotSize);
            if nK > 0 && ~isempty(tr.pulseMs{ks(1)})
                pt = tr.pulseMs{ks(1)};
                plot(axR, [pt(:)'; pt(:)'], repmat([nK + 0.8; nK + 1.6], 1, numel(pt)), 'r-', 'LineWidth', 0.8);
            end
            xline(axR, 0, ':');  xlim(axR, winMs);  ylim(axR, [0.5 max(nK, 1) + 1.8]);  set(axR, 'XTickLabel', []);  box(axR, 'on');
            if col == 1, ylabel(axR, 'trial'); end
            title(axR, sprintf('%d: %s <- %d (n=%d)', c, CT(c).label, CT(c).sourceCond, nK), 'FontSize', 8, 'Interpreter', 'none');
            % PSTHs
            rR = psth_(Rc{ch}, RC.blank_ms(ch, :), ks, edges, ker, minDataFrac);
            rO = psth_(Oc{ch}, SO.blank_ms(ch, :), ks, edges, ker, minDataFrac);
            sk = find(trialCond == CT(c).sourceCond);  sk = sk(:)';
            rS = psth_(Rc{ch}, RC.blank_ms(ch, :), sk, edges, ker, minDataFrac);
            axP(i) = subplot(3*nRB, nCol, pP);  hold(axP(i), 'on');
            plot(axP(i), ctr, rO, 'Color', [0.6 0.6 0.6], 'LineWidth', 2.5);
            plot(axP(i), ctr, rS, '--', 'Color', [0.1 0.3 0.9], 'LineWidth', 1.2);
            plot(axP(i), ctr, rR, 'k', 'LineWidth', 1.0);
            xline(axP(i), 0, ':');  xlim(axP(i), winMs);  box(axP(i), 'on');
            yMax = max([yMax, rR(~isnan(rR)), rO(~isnan(rO)), rS(~isnan(rS))]);
            if col == 1, ylabel(axP(i), 'spikes/s'); end
            if i == 1, legend(axP(i), {'original', 'source', 'recovered'}, 'FontSize', 7, 'Location', 'northwest'); end
            if row == nRB || i + nCol > nC, xlabel(axP(i), 'ms after the trigger'); end
        end
        if yMax > 0, set(axP, 'YLim', [0 1.1*yMax]); end
        sgtitle(fh, sprintf('%s  ch %d (y %d um)  %s: recovery check\nraster: black = original, orange = copied, grey = no data | PSTH: grey = original, blue dashed = source level, black = recovered', ...
                strrep(short, '_', '\_'), ch, y(ch), grpList{gI}), 'FontSize', 9);

        if showLevelOverlay && nC > 1
            [~, o] = sort([CT(cs).level]);  cs2 = cs(o);
            fo = figure('Name', sprintf('%s ch %d %s levels', short, ch, grpList{gI}), 'Color', 'w', 'Position', [80 80 900 450]);
            ax = axes(fo);  hold(ax, 'on');
            cm = parula(numel(cs2) + 2);  cm = flipud(cm(1:numel(cs2), :));     % light (low) -> dark (high)
            lab = cell(1, numel(cs2));
            for i = 1:numel(cs2)
                c = cs2(i);  ks = find(trialCond == c);  ks = ks(:)';
                plot(ax, ctr, psth_(Rc{ch}, RC.blank_ms(ch, :), ks, edges, ker, minDataFrac), 'Color', cm(i,:), 'LineWidth', 1.2);
                lab{i} = sprintf('%d: %s', c, CT(c).label);
            end
            xline(ax, 0, ':');  xlim(ax, winMs);  grid(ax, 'on');  box(ax, 'on');
            xlabel(ax, 'ms after the trigger');  ylabel(ax, 'spikes/s (recovered)');
            legend(ax, lab, 'FontSize', 7, 'Location', 'eastoutside', 'Interpreter', 'none');
            title(ax, sprintf('%s  ch %d (y %d um)  %s: recovered PSTH of each level (should overlap up to each level''s last pulse)', ...
                  strrep(short, '_', '\_'), ch, y(ch), grpList{gI}), 'FontSize', 9);
        end
    end
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function B = byTrial_(tAbs, k, trigMs, nTr)
% spikes in ms after each trial's trigger, one sorted column per trial
B = cell(nTr, 1);  for t = 1:nTr, B{t} = zeros(0, 1); end
if isempty(k), return; end
k = k(:);  rel = tAbs(:) - trigMs(k);
[ks, o] = sort(k);  rel = rel(o);
e = [0; find(diff(ks)); numel(ks)];
for i = 1:numel(e) - 1
    rows = e(i)+1:e(i+1);  B{ks(rows(1))} = sort(rel(rows));
end
end

function m = inIv_(t, R, h)
% times inside any [start end] row of R, widened by h on both sides
m = false(size(t));
for b = 1:size(R, 1), m = m | (t >= R(b,1) - h & t <= R(b,2) + h); end
end

function tf = same_(a, b, tol)
tf = numel(a) == numel(b) && (isempty(a) || max(abs(sort(a) - sort(b))) < tol);
end

function n = cnt_(t, w)
n = sum(t >= w(1) & t < w(2));
end

function rt = psth_(B, blankRow, ks, edges, ker, minFrac)
% blank-aware PSTH (spikes/s) of trials ks; bins with data in fewer than minFrac of the trials are NaN
binMs = edges(2) - edges(1);  nB = numel(edges) - 1;
cnt = zeros(1, nB);  expo = zeros(1, nB);
for t = ks
    cnt = cnt + histcounts(B{t}, edges);
    covr = ones(1, nB);  iv = blankRow{t};
    for b = 1:size(iv, 1)
        lo = max(edges(1:end-1), iv(b,1));  hi = min(edges(2:end), iv(b,2));
        covr = covr - max(0, hi - lo)/binMs;
    end
    expo = expo + max(covr, 0);
end
cS = conv(cnt, ker, 'same');  eS = conv(expo, ker, 'same');
rt = cS ./ (eS*binMs/1e3);
rt(expo < minFrac*max(numel(ks), 1) | eS <= 0) = NaN;
end

function s = n2s_(v)
if isnan(v), s = '-'; else, s = sprintf('%d', v); end
end

function s = passTxt_(nBad)
if nBad == 0, s = '-> PASS'; else, s = sprintf('-> %d MISMATCHES', nBad); end
end
