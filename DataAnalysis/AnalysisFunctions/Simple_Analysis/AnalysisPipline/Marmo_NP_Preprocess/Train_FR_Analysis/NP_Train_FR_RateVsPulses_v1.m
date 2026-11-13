%% ============================================================
% NP_Train_FR_RateVsPulses_v1.m
% NP probe | pulse TRAIN data | FIRING-RATE analysis | RATE (and COUNT) vs NUMBER OF PULSES per depth group | v1
% (NP version of NN_Train_FR_CountVsPulses_v1 / _Norm_v1, without normalisation)
% READ-ONLY: one session (one spike file) per run; prints and plots, NOTHING IS SAVED.
%
% PER LEVEL (every recovered level is its own condition), PER CHANNEL
%   Window   : first pulse (+trainStart_ms) ... last pulse of the level + tail_ms
%   Count    : spikes in the window - baseline spikes (baseWinMs, same trial) scaled to the window length,
%              mean over the (good) trials of the level
%   RATE     : count / window length (spikes/s above baseline)                                 <- MAIN MEASURE
%   ADDED    : count(level) - count(previous plotted level)              (added spikes)
%   ADDED RATE = added count / (window(level) - window(previous level))  = rate of the newly added stretch of the
%              train (spikes/s); the first plotted level: its own rate. Negative = suppression.
%   Group value = mean over the channels of the group, error bar = SEM over channels (channels are the samples).
%   x axis   : pulses per electrode (sim n; seq [n n] = n, half-steps [n n-1] = n - 0.5), or the time of the last pulse.
%
% CHANNEL GROUPS (chGroups)
%   'all'                              : one group, every included channel
%   {'name', [ymin ymax]; ...}         : depth ranges from the probe TIP (um), one curve per group
%   Included channels inside a group:
%     - processed in step 2 (goodCh) and active: rate over the whole window >= autoRateMin (spikes/s);
%     - respFilter = true: and RESPONSIVE: in the final level of at least one pattern the train-window rate is above
%       baseline by > respZ standard errors (over trials). The same channels are used for sim and seq (union).
%     - if <short>.NP_RespChn.mat exists and useRespChn = true, its channels (RespChn.respChn{cond}, union over the
%       final conditions) replace this automatic test (file written by a later script).
%   Bad trials: if <short>.NP_BadTrials.mat exists and useBadTrials = true, BadTrials.badRel{cond} (trial numbers
%   within the condition) are left out; otherwise every trial is kept.
%
% FIGURES
%   1. FIRING RATE : rate vs x  and  added rate vs x      (colour = depth group; sim solid / seq dashed)
%   2. SPIKE COUNT : count vs x and  added spikes vs x     (Plot_CountFigure)
%   3. CHANNELS    : final-level rate of every channel vs depth, the groups and the included channels marked (Plot_ChannelFigure)
%   Printed: the groups and their channels, and a table of every level.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile  = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery_randomEvent.mat';            % '' = <short>.NP_Spikes_TrainRecovery.mat; or another recovered file (e.g. a _RandomEvent one)

% ---- channel groups (depth from the probe TIP, um). Session-specific: read them off NP_SP_DepthTimeMap. ----
chGroups = { ...
    'deep (0-2200)',     [ 0 2200]; ...
    'middle (3000-4500)',  [3000 4500]; ...
    % 'upper (4500-5800)',   [4500 5800]; ...
    };                      % or chGroups = 'all';
autoRateMin = 20;           % spikes/s over the whole window: a channel counts as active above this
respFilter  = true;         % true = keep only channels responsive in the final level of sim or seq (see header)
respZ       = 2;            % responsive if (train-window rate - baseline) / SE over trials > respZ
useRespChn   = false;        % use <short>.NP_RespChn.mat when it exists 
useBadTrials = false;        % use <short>.NP_BadTrials.mat when it exists

% ---- counting window and baseline ----
tail_ms       = 15;         % window = first pulse + trainStart_ms ... last pulse + tail_ms
trainStart_ms = 0;
baseWinMs     = [-300 -5];  % baseline window relative to the trigger (ms)

% ---- levels and axes ----
seqLevels = 'pair';         % 'pair' ([n n] only, comparable with sim n) | 'all' (also the half-steps)
xAxis     = 'n';            % 'n' = pulses per electrode | 'time' = time of the last pulse (ms)

% ---- figures ----
Plot_CountFigure   = true;  % spike count and added spikes
Plot_ChannelFigure = true;  % final-level rate of every channel vs depth
Show_Channels      = false; % thin line per channel behind the group means (figure 1)
Marker_Size = 6;  Line_Width = 1.5;
%% =================================================

assert(any(strcmp(seqLevels, {'pair', 'all'})), 'seqLevels: ''pair'' or ''all''.');
assert(any(strcmp(xAxis, {'n', 'time'})), 'xAxis: ''n'' or ''time''.');

%% ---- files -------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
if isempty(spikeFile), spikeFile = [short '.NP_Spikes_TrainRecovery.mat']; end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
assert(isfile(spikeFile), 'Spike file not found: %s (run NP_3_SpikeRecovery first).', spikeFile);
SP = load(spikeFile, 'sp_clipped', 'sp_trial', 'nSpk', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed');
fs = SP.FS;  nCh = numel(SP.sp_clipped);  nTr = SP.nTrialsProcessed;
tr = TT.trial;  trialCond = tr.cond(1:nTr);  trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;
y = TT.chan.y_um(:);
assert(baseWinMs(1) >= -SP.preMs, 'baseWinMs starts before the detection window (-%g ms).', SP.preMs);

% optional files
badRel = cell(numel(TT.cond), 1);  bf = fullfile(recDir, [short '.NP_BadTrials.mat']);
if useBadTrials && isfile(bf)
    Bf = load(bf, 'BadTrials');  badRel = Bf.BadTrials.badRel;  fprintf('bad trials from %s\n', bf);
else
    fprintf('no bad-trial file used: all trials kept\n');
end
rf = fullfile(recDir, [short '.NP_RespChn.mat']);  RespC = [];
if useRespChn && isfile(rf), Rf = load(rf, 'RespChn'); RespC = Rf.RespChn; fprintf('responsive channels from %s\n', rf); end

%% ---- spikes per trial -----------------------------------------------------------------------------------
B = cell(nCh, 1);
for ch = 1:nCh, B{ch} = byTrial_(SP.sp_clipped{ch}, SP.sp_trial{ch}, trigMs, nTr); end

%% ---- patterns and levels to plot ---------------------------------------------------------------------------
P = TT.pattern;  nP = numel(P);
assert(nP > 0, 'No train pattern in the trial table.');
Lv = struct('name', {}, 'conds', {}, 'x', {}, 'win', {}, 'ks', {});
for p = 1:nP
    cs = P(p).conds;  keep = true(size(cs));  xv = nan(size(cs));  wv = nan(size(cs));
    for i = 1:numel(cs)
        C = TT.cond(cs(i));  nsp = C.nsp;  amp = C.amp;  nsp(amp <= 0) = 0;
        if strcmp(P(p).kind, 'sim')
            xv(i) = max(nsp);
        else
            if nsp(1) == nsp(2), xv(i) = nsp(1); else, xv(i) = max(nsp) - 0.5; keep(i) = strcmp(seqLevels, 'all'); end
        end
        ks = find(trialCond == cs(i));
        if isempty(ks), keep(i) = false; continue; end
        if strcmp(xAxis, 'time'), xv(i) = tr.lastPulseMs(ks(1)); end
        wv(i) = (tr.lastPulseMs(ks(1)) + tail_ms) - (tr.firstPulseMs(ks(1)) + trainStart_ms);
    end
    cs = cs(keep);  xv = xv(keep);  wv = wv(keep);
    kk = cell(1, numel(cs));
    for i = 1:numel(cs)
        ks = find(trialCond == cs(i))';  bad = badRel{cs(i)};
        if ~isempty(bad), ks(bad(bad >= 1 & bad <= numel(ks))) = []; end
        kk{i} = ks;
    end
    switch P(p).kind, case 'sim', nm = 'sim'; case 'AB', nm = 'seq A->B'; case 'BA', nm = 'seq B->A'; otherwise, nm = P(p).kind; end
    Lv(end+1) = struct('name', nm, 'conds', cs, 'x', xv, 'win', wv, 'ks', {kk}); %#ok<SAGROW>
end

%% ---- per channel and level: count and rate ------------------------------------------------------------------
cntC = cell(nP, 1);  zC = cell(nP, 1);
baseDur = diff(baseWinMs)/1e3;
for p = 1:nP
    nL = numel(Lv(p).conds);  cntC{p} = nan(nCh, nL);  zC{p} = nan(nCh, nL);
    for i = 1:nL
        ks = Lv(p).ks{i};  k1 = ks(1);
        w = [tr.firstPulseMs(k1) + trainStart_ms, tr.lastPulseMs(k1) + tail_ms];  wd = diff(w)/1e3;
        for ch = find(SP.goodCh(:))'
            v = zeros(numel(ks), 1);
            for q = 1:numel(ks)
                s = B{ch}{ks(q)};
                nb = sum(s >= baseWinMs(1) & s < baseWinMs(2));
                v(q) = sum(s >= w(1) & s < w(2)) - nb/baseDur*wd;              % baseline-subtracted count
            end
            cntC{p}(ch, i) = mean(v);
            se = std(v)/sqrt(numel(v));
            if se > 0, zC{p}(ch, i) = mean(v)/se; elseif mean(v) > 0, zC{p}(ch, i) = Inf; else, zC{p}(ch, i) = 0; end
        end
    end
end

%% ---- channels per group ---------------------------------------------------------------------------------------------
rateAll = sum(SP.nSpk(:, 1:nTr), 2)/nTr/((SP.preMs + SP.postMs)/1e3);
incl = SP.goodCh(:) & rateAll >= autoRateMin;
if ~isempty(RespC)
    rc = [];
    for p = 1:nP, rc = [rc, RespC.respChn{Lv(p).conds(end)}]; end %#ok<AGROW>
    incl = incl & ismember((1:nCh)', unique(rc));
elseif respFilter
    resp = false(nCh, 1);
    for p = 1:nP, resp = resp | zC{p}(:, end) > respZ; end
    incl = incl & resp;
end
if ischar(chGroups) && strcmp(chGroups, 'all')
    G = struct('name', 'all channels', 'ch', find(incl)', 'range', [min(y) max(y)]);
else
    G = struct('name', {}, 'ch', {}, 'range', {});
    for g = 1:size(chGroups, 1)
        r = chGroups{g, 2};
        G(end+1) = struct('name', chGroups{g, 1}, 'ch', find(incl & y >= r(1) & y <= r(2))', 'range', r); %#ok<SAGROW>
    end
end
nG = numel(G);

%% ---- group statistics -------------------------------------------------------------------------------------------------
S = struct();
for p = 1:nP
    win = Lv(p).win(:)';  dW = [win(1), diff(win)];
    cm  = cntC{p};                                         % channels x levels
    add = [cm(:, 1), diff(cm, 1, 2)];
    S(p).rate    = bsxfun(@rdivide, cm, win/1e3);
    S(p).addRate = bsxfun(@rdivide, add, dW/1e3);
    S(p).cnt = cm;  S(p).add = add;
end
fld = {'rate', 'addRate', 'cnt', 'add'};

%% ---- printout ---------------------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP TRAIN FR: RATE vs PULSES  %s\n', short);
fprintf('============================================================\n');
fprintf('spike file : %s\nwindow     : first pulse + %g ms ... last pulse + %g ms; baseline %s ms\n', spikeFile, trainStart_ms, tail_ms, mat2str(baseWinMs));
fprintf('channels   : goodCh, rate >= %g /s%s\n', autoRateMin, ternary_(~isempty(RespC), ', NP_RespChn file', ternary_(respFilter, sprintf(', responsive (z > %g in a final level)', respZ), '')));
for g = 1:nG
    fprintf('  group %-22s y %4d-%4d um: %3d channels\n', G(g).name, G(g).range(1), G(g).range(2), numel(G(g).ch));
end
for g = 1:nG
    fprintf('\n--- %s (%d channels) ---\n  pattern    cond      x   trials   rate(Hz)  +/-SEM   addedRate +/-SEM  |   count  +/-SEM   added  +/-SEM\n', G(g).name, numel(G(g).ch));
    for p = 1:nP
        for i = 1:numel(Lv(p).conds)
            v = zeros(1, 8);
            for f = 1:4, [v(2*f-1), v(2*f)] = ms_(S(p).(fld{f})(G(g).ch, i)); end
            fprintf('  %-9s  %4d  %5.1f   %6d   %7.1f  %6.1f   %8.1f  %6.1f  | %7.2f  %6.2f  %6.2f  %6.2f\n', Lv(p).name, Lv(p).conds(i), ...
                    Lv(p).x(i), numel(Lv(p).ks{i}), v);
        end
    end
end

%% ---- figures -----------------------------------------------------------------------------------------------------------
colG = lines(max(nG, 1));
styP = {'-o', '--s', ':^'};  mkP = {'o', 's', '^'};   % line style / marker per pattern
xlab = ternary_(strcmp(xAxis, 'n'), 'pulses per electrode', 'time of the last pulse (ms)');
figs = {'Firing rate', 'rate', 'addRate', 'rate above baseline (spikes/s)', 'added rate (spikes/s of the new stretch)'};
if Plot_CountFigure, figs(end+1, :) = {'Spike count', 'cnt', 'add', 'count above baseline (spikes/trial)', 'added spikes (count n - count n-1)'}; end
for fi = 1:size(figs, 1)
    fh = figure('Name', sprintf('%s: %s', short, figs{fi, 1}), 'Color', 'w', 'Position', [60 60 1100 480]);
    for sp = 1:2
        ax = subplot(1, 2, sp);  hold(ax, 'on');  hl = [];  lab = {};
        for g = 1:nG
            if isempty(G(g).ch), continue; end
            for p = 1:nP
                M = S(p).(figs{fi, 1 + sp})(G(g).ch, :);
                if Show_Channels && fi == 1
                    plot(ax, Lv(p).x, M', '-', 'Color', [colG(g,:) 0.15], 'LineWidth', 0.5, 'HandleVisibility', 'off');
                end
                mu = mean(M, 1, 'omitnan');  se = std(M, 0, 1, 'omitnan')./sqrt(sum(~isnan(M), 1));
                h = errorbar(ax, Lv(p).x, mu, se, styP{min(p, 3)}, 'Color', colG(g,:), 'MarkerFaceColor', ternary_(p == 1, colG(g,:), 'w'), ...
                             'MarkerSize', Marker_Size, 'LineWidth', Line_Width, 'CapSize', 3);
                hl(end+1) = h; lab{end+1} = sprintf('%s | %s (%d ch)', G(g).name, Lv(p).name, numel(G(g).ch)); %#ok<SAGROW>
            end
        end
        yline(ax, 0, ':');  grid(ax, 'on');  xlabel(ax, xlab);  ylabel(ax, figs{fi, 3 + sp});
        if strcmp(xAxis, 'n'), xticks(ax, 0:1:ceil(max([Lv.x]))); end
        if sp == 1 && ~isempty(hl), legend(ax, hl, lab, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'none'); end
    end
    sgtitle(fh, sprintf('%s: %s vs number of pulses (mean +/- SEM over channels); window first pulse ... last + %g ms', ...
            strrep(short, '_', '\_'), lower(figs{fi, 1}), tail_ms), 'FontSize', 10);
end
if Plot_ChannelFigure
    fh = figure('Name', sprintf('%s: channels', short), 'Color', 'w', 'Position', [100 100 600 750]);
    ax = axes(fh);  hold(ax, 'on');
    for g = 1:nG
        r = G(g).range;  patch(ax, [-1e4 1e4 1e4 -1e4], [r(1) r(1) r(2) r(2)], colG(g,:), 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hl = [];  lab = {};
    for p = 1:nP
        v = S(p).rate(:, end);
        hl(end+1) = plot(ax, v, y, mkP{min(p, 3)}, 'Color', [0.5 0.5 0.5], 'MarkerSize', 3); %#ok<SAGROW>
        lab{end+1} = sprintf('%s final level, all channels', Lv(p).name); %#ok<SAGROW>
        inc = find(incl);
        plot(ax, v(inc), y(inc), mkP{min(p, 3)}, 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'HandleVisibility', 'off');
    end
    xline(ax, 0, ':');  grid(ax, 'on');  ylim(ax, [min(y) max(y)]);
    if strcmp(xAxis, 'n'), xticks(ax, 0:1:ceil(max([Lv.x]))); end
    vv = cell2mat(arrayfun(@(s) s.rate(:, end), S(:), 'UniformOutput', false));  vv = vv(isfinite(vv));
    if ~isempty(vv), xlim(ax, [min(0, min(vv)) max(vv)]*1.05 + [-1 1]); end
    xlabel(ax, 'rate above baseline, final level (spikes/s)');  ylabel(ax, 'depth from the tip (um)');
    legend(ax, hl, lab, 'Location', 'best', 'FontSize', 8);
    title(ax, sprintf('%s: black = included channels; shaded = groups', strrep(short, '_', '\_')), 'FontSize', 10);
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function Bt = byTrial_(tAbs, k, trigMs, nTr)
Bt = cell(nTr, 1);  for t = 1:nTr, Bt{t} = zeros(0, 1); end
if isempty(k), return; end
k = k(:);  rel = tAbs(:) - trigMs(k);
[ks, o] = sort(k);  rel = rel(o);
e = [0; find(diff(ks)); numel(ks)];
for i = 1:numel(e) - 1
    rows = e(i)+1:e(i+1);  Bt{ks(rows(1))} = sort(rel(rows));
end
end

function [m, s] = ms_(v)
v = v(~isnan(v));
if isempty(v), m = NaN; s = NaN; return; end
m = mean(v);  s = std(v)/sqrt(numel(v));
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
