%% ============================================================
% NP_Train_FR_RateVsPulses_AreaBad_v2.m
% NP probe | pulse TRAIN data | FIRING-RATE analysis | RATE (and COUNT) vs NUMBER OF PULSES per depth group | v2 (PER-AREA BAD TRIALS)
% (NP version of NN_Train_FR_CountVsPulses_v1 / _Norm_v1, without normalisation)
% READ-ONLY: one session (one spike file) per run; prints and plots, NOTHING IS SAVED.
% v1 = NP_Train_FR_RateVsPulses_v1 is NOT changed. The ONLY difference of v2: the bad trials are loaded PER AREA.
%
% BAD TRIALS (v2)
%   File <short>.NP_BadTrialsArea.mat (struct BadTrialsArea, written by NP_Train_FR_BadTrialEntry_v2). Each channel
%   uses the list of ITS OWN area: BadTrialsArea.badRel.(area){cond} = the 'all' trials + that area's own trials.
%     - channel in V1 / V2 (label from NP_GetAreas)        -> the list of that area
%     - border channel with borderAs = 'V1' or 'V2'          -> the list of that area (pooled into it)
%     - border ('own'), unassigned, or chGroups = 'all' / depth ranges (no area split) -> the UNION list
%       (BadTrialsArea.badRelUnion: a trial bad in any area is left out; conservative)
%   A condition not reviewed for an area keeps all its trials (empty list). Trial numbers are WITHIN the condition.
%   If the v2 file does not exist and fallbackV1 = true, the old <short>.NP_BadTrials.mat (one list for all channels)
%   is used for every channel; otherwise all trials are kept. Because the trials differ between areas, the number of
%   trials per level can differ between V1 and V2 (printed in the table).
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
%   'area'                             : one curve per AREA (V1, V2; border as set by borderAs), same axes, area colours.
%                                        Needs areaSource = 'file' or 'manual' (labels from NP_GetAreas, RF folder).
%                                        'outside' and 'unassigned' channels are never in an area group.
%   excludeBad = true (and an area source): channels labelled 'bad' are left out of EVERY grouping.
%   borderAs ('area' only): 'exclude' (default) | 'V1' / 'V2' (border channels pooled into that area: sensitivity
%                           check) | 'own' (border as its own curve)
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
%   1. FIRING RATE : rate vs x  and  added rate vs x      (colour = group; sim solid / seq dashed)
%      groupLayout = 'separate': one row per group (V1 / V2 ... each in its own panels, same limits if sameYLim)
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
% chGroups = { ...
    % 'deep (0-3000)',     [ 0 3000]; ...
    % 'middle (3020-4500)',  [3020 4500]; ...
    % 'upper (4500-5800)',   [4500 5800]; ...
    % };                      % or chGroups = 'all';
chGroups = 'area';           % 'all' | 'area' | depth-range table (above)

% ---- areas (labels from NP_GetAreas.m) ----
areaFuncDir = '/Volumes/LaCie/Marmoset_Data/AnalysisFunctions/NP_Analysis/RF_Detection/NP_RF';   % where NP_GetAreas.m is
areaSource  = 'file';       % 'file' (area file of NP_AreaEntry_v1) | 'manual' (ranges below) | 'none'
areaFile    = '/Volumes/LaCie/Marmoset_Data/CJ287/CJ287_NP_Areas_20261110_0506.mat';
areasManual = { ...         % used when areaSource = 'manual': rows {from_um to_um label}
    0       2600   'V1';
    2600    2800   'border';
    2800    5400   'V2';
    5400    6000   'V2';
    6000    7660   'outside';
    };
badManual   = [83];         % used when areaSource = 'manual'
excludeBad  = true;         % true = channels labelled 'bad' are never included (any chGroups)
borderAs    = 'exclude';    % chGroups = 'area': 'exclude' | 'V1' | 'V2' | 'own'
autoRateMin = 20;           % spikes/s over the whole window: a channel counts as active above this
respFilter  = true;         % true = keep only channels responsive in the final level of sim or seq (see header)
respZ       = 2;            % responsive if (train-window rate - baseline) / SE over trials > respZ
useRespChn   = true;        % use <short>.NP_RespChn.mat when it exists 
useBadTrials = true;        % use the bad-trial file(s) when they exist
fallbackV1   = true;        % if <short>.NP_BadTrialsArea.mat is missing: use the old <short>.NP_BadTrials.mat for all channels

% ---- counting window and baseline ----
tail_ms       = 15;         % window = first pulse + trainStart_ms ... last pulse + tail_ms
trainStart_ms = 0;
baseWinMs     = [-100 -5];  % baseline window relative to the trigger (ms)

% ---- levels and axes ----
seqLevels = 'pair';         % 'pair' ([n n] only, comparable with sim n) | 'all' (also the half-steps)
xAxis     = 'n';            % 'n' = pulses per electrode | 'time' = time of the last pulse (ms)

% ---- figures ----
Plot_CountFigure   = true;  % spike count and added spikes
Plot_ChannelFigure = true;  % final-level rate of every channel vs depth
Show_Channels      = false; % thin line per channel behind the group means (figure 1)
groupLayout        = 'separate';  % 'separate' = one row of panels per group (e.g. V1 above, V2 below) | 'together' = one axes
sameYLim           = true;  % 'separate': the same axis limits in every group, so they can be compared
Marker_Size = 6;  Line_Width = 1.5;
%% =================================================

assert(any(strcmp(seqLevels, {'pair', 'all'})), 'seqLevels: ''pair'' or ''all''.');
assert(any(strcmp(xAxis, {'n', 'time'})), 'xAxis: ''n'' or ''time''.');
assert(any(strcmp(borderAs, {'exclude', 'V1', 'V2', 'own'})), 'borderAs: ''exclude'', ''V1'', ''V2'' or ''own''.');
byArea = ischar(chGroups) && strcmp(chGroups, 'area');

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
useArea = ~strcmpi(areaSource, 'none');
assert(~byArea || useArea, 'chGroups = ''area'' needs areaSource = ''file'' or ''manual''.');
if useArea, addpath(areaFuncDir); end
[aLab, aInfo] = NP_GetAreas(areaSource, areaFile, areasManual, badManual, y);
assert(baseWinMs(1) >= -SP.preMs, 'baseWinMs starts before the detection window (-%g ms).', SP.preMs);

% optional files: bad trials PER AREA (lists: one per area in BTA.areas, then 'union')
nCondT = numel(TT.cond);  emptyC_ = cell(nCondT, 1);
listNames = {};  listBad = {};  haveAreaBad = false;
bfA = fullfile(recDir, [short '.NP_BadTrialsArea.mat']);  bf1 = fullfile(recDir, [short '.NP_BadTrials.mat']);
if useBadTrials && isfile(bfA)
    Bf = load(bfA, 'BadTrialsArea');  BTA = Bf.BadTrialsArea;  haveAreaBad = true;
    for a = 1:numel(BTA.areas), listNames{end+1} = BTA.areas{a};  listBad{end+1} = BTA.badRel.(BTA.areas{a}); end %#ok<SAGROW>
    listNames{end+1} = 'union';  listBad{end+1} = BTA.badRelUnion;
    fprintf('bad trials PER AREA from %s (areas: %s)\n', bfA, strjoin(BTA.areas, ', '));
elseif useBadTrials && fallbackV1 && isfile(bf1)
    Bf = load(bf1, 'BadTrials');
    listNames = {'union'};  listBad = {Bf.BadTrials.badRel};
    fprintf('NO per-area file; the old bad-trial file %s is used for ALL channels\n', bf1);
else
    listNames = {'union'};  listBad = {emptyC_};
    fprintf('no bad-trial file used: all trials kept\n');
end
for l = 1:numel(listBad), if numel(listBad{l}) < nCondT, listBad{l}(nCondT) = {[]}; end, end
iUnion = find(strcmp(listNames, 'union'));
% list used by every channel
chList = iUnion*ones(nCh, 1);
for ch = 1:nCh
    if ~byArea, break; end                       % no area split: every channel uses the union list
    lb = aLab{ch};
    if strcmp(lb, 'border') && any(strcmp(borderAs, {'V1', 'V2'})), lb = borderAs; end
    j = find(strcmp(listNames, lb), 1);
    if byArea && ~isempty(j) && ~strcmp(lb, 'union'), chList(ch) = j; end
end
rf = fullfile(recDir, [short '.NP_RespChn.mat']);  RespC = [];
if useRespChn && isfile(rf), Rf = load(rf, 'RespChn'); RespC = Rf.RespChn; fprintf('responsive channels from %s\n', rf); end

%% ---- spikes per trial -----------------------------------------------------------------------------------
B = cell(nCh, 1);
for ch = 1:nCh, B{ch} = byTrial_(SP.sp_clipped{ch}, SP.sp_trial{ch}, trigMs, nTr); end

%% ---- patterns and levels to plot ---------------------------------------------------------------------------
P = TT.pattern;  nP = numel(P);
assert(nP > 0, 'No train pattern in the trial table.');
Lv = struct('name', {}, 'conds', {}, 'x', {}, 'win', {}, 'ks', {}, 'ks0', {});
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
    kk = cell(numel(listNames), numel(cs));  k0 = cell(1, numel(cs));       % kk{list, level}: trials kept for that list
    for i = 1:numel(cs)
        ks = find(trialCond == cs(i))';  k0{i} = ks;
        for l = 1:numel(listNames)
            bad = listBad{l}{cs(i)};  kl = ks;
            if ~isempty(bad), kl(bad(bad >= 1 & bad <= numel(ks))) = []; end
            kk{l, i} = kl;
        end
    end
    switch P(p).kind, case 'sim', nm = 'sim'; case 'AB', nm = 'seq A->B'; case 'BA', nm = 'seq B->A'; otherwise, nm = P(p).kind; end
    Lv(end+1) = struct('name', nm, 'conds', cs, 'x', xv, 'win', wv, 'ks', {kk}, 'ks0', {k0}); %#ok<SAGROW>
end

%% ---- per channel and level: count and rate ------------------------------------------------------------------
cntC = cell(nP, 1);  zC = cell(nP, 1);
baseDur = diff(baseWinMs)/1e3;
for p = 1:nP
    nL = numel(Lv(p).conds);  cntC{p} = nan(nCh, nL);  zC{p} = nan(nCh, nL);
    for i = 1:nL
        k1 = Lv(p).ks0{i}(1);                              % timing is the same in every trial of the level
        w = [tr.firstPulseMs(k1) + trainStart_ms, tr.lastPulseMs(k1) + tail_ms];  wd = diff(w)/1e3;
        for ch = find(SP.goodCh(:))'
            ks = Lv(p).ks{chList(ch), i};                  % this channel's own (area) trial list
            if isempty(ks), continue; end
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
isBad = false(nCh, 1);
if useArea && excludeBad, isBad = strcmp(aLab, 'bad');  incl = incl & ~isBad; end
areaCol = struct('V1', [0.2 0.4 0.9], 'V2', [0.9 0.45 0.1], 'border', [0.45 0.45 0.45]);
G = struct('name', {}, 'ch', {}, 'range', {}, 'bands', {}, 'col', {});
if ischar(chGroups) && strcmp(chGroups, 'all')
    G(1).name = 'all channels';  G(1).ch = find(incl)';  G(1).range = [min(y) max(y)];  G(1).bands = G(1).range;  G(1).col = [0 0 0];
elseif byArea
    labG = aLab;  if ~excludeBad, labG = aInfo.rangeLabel; end
    rg = aInfo.ranges;  isBorderR = strcmp(rg(:, 3), 'border');
    names = {'V1', 'V2'};  if strcmp(borderAs, 'own'), names{end+1} = 'border'; end
    for a = 1:numel(names)
        lb = names{a};  inA = strcmp(labG, lb);  rIdx = strcmp(rg(:, 3), lb);  gName = lb;
        if strcmp(borderAs, lb)                                    % border pooled into this area
            inA = inA | strcmp(labG, 'border');  rIdx = rIdx | isBorderR;  gName = [lb ' (+border)'];
        end
        ch = find(incl & inA)';
        g = numel(G) + 1;
        G(g).name = gName;  G(g).ch = ch;  G(g).bands = cell2mat(rg(rIdx, 1:2));  G(g).col = areaCol.(lb);
        if isempty(ch), G(g).range = [NaN NaN]; else, G(g).range = [min(y(ch)) max(y(ch))]; end
    end
else
    for g = 1:size(chGroups, 1)
        r = chGroups{g, 2};
        G(g).name = chGroups{g, 1};  G(g).ch = find(incl & y >= r(1) & y <= r(2))';  G(g).range = r;  G(g).bands = r;  G(g).col = [];
    end
end
nG = numel(G);
% list (bad-trial set) of each group, for the trial counts in the table
gListIdx = iUnion*ones(1, nG);
for g = 1:nG
    ch = G(g).ch;
    if ~isempty(ch) && all(chList(ch) == chList(ch(1))), gListIdx(g) = chList(ch(1)); end
end
gList_ = @(g) gListIdx(g);
cl = lines(max(nG, 1));
for g = 1:nG, if isempty(G(g).col), G(g).col = cl(g, :); end, end

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
fprintf('%s\n', aInfo.text);
if any(isBad), fprintf('bad        : channels %s excluded\n', mat2str(find(isBad)')); end
if byArea, fprintf('grouping   : by area, border %s\n', ternary_(strcmp(borderAs, 'exclude'), 'excluded', ternary_(strcmp(borderAs, 'own'), 'as its own group', ['pooled into ' borderAs]))); end
for g = 1:nG
    fprintf('  group %-22s y %4g-%4g um: %3d channels, bad-trial list: %s\n', G(g).name, G(g).range(1), G(g).range(2), numel(G(g).ch), listNames{gListIdx(g)});
end
for g = 1:nG
    fprintf('\n--- %s (%d channels) ---\n  pattern    cond      x   trials(kept)   rate(Hz)  +/-SEM   addedRate +/-SEM  |   count  +/-SEM   added  +/-SEM\n', G(g).name, numel(G(g).ch));
    for p = 1:nP
        for i = 1:numel(Lv(p).conds)
            v = zeros(1, 8);
            for f = 1:4, [v(2*f-1), v(2*f)] = ms_(S(p).(fld{f})(G(g).ch, i)); end
            fprintf('  %-9s  %4d  %5.1f   %6d   %7.1f  %6.1f   %8.1f  %6.1f  | %7.2f  %6.2f  %6.2f  %6.2f\n', Lv(p).name, Lv(p).conds(i), ...
                    Lv(p).x(i), numel(Lv(p).ks{gList_(g), i}), v);
        end
    end
end

%% ---- figures -----------------------------------------------------------------------------------------------------------
colG = cat(1, G.col);  if isempty(colG), colG = lines(1); end
styP = {'-o', '--s', ':^'};  mkP = {'o', 's', '^'};   % line style / marker per pattern
xlab = ternary_(strcmp(xAxis, 'n'), 'pulses per electrode', 'time of the last pulse (ms)');
figs = {'Firing rate', 'rate', 'addRate', 'rate above baseline (spikes/s)', 'added rate (spikes/s of the new stretch)'};
if Plot_CountFigure, figs(end+1, :) = {'Spike count', 'cnt', 'add', 'count above baseline (spikes/trial)', 'added spikes (count n - count n-1)'}; end
for fi = 1:size(figs, 1)
    sepG  = strcmp(groupLayout, 'separate') && nG > 1;          % one row of panels per group
    nRowG = ternary_(sepG, nG, 1);
    fh = figure('Name', sprintf('%s: %s', short, figs{fi, 1}), 'Color', 'w', 'Position', [60 60 1100 min(980, 150 + 330*nRowG)]);
    axAll = gobjects(nRowG, 2);
    for sp = 1:2
      for gr = 1:nRowG
        ax = subplot(nRowG, 2, (gr - 1)*2 + sp);  hold(ax, 'on');  hl = [];  lab = {};  axAll(gr, sp) = ax;
        gList = ternary_(sepG, gr, 1:nG);
        for g = gList
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
        yline(ax, 0, ':');   xlabel(ax, xlab);  ylabel(ax, figs{fi, 3 + sp});
        if strcmp(xAxis, 'n'), xticks(ax, 0:1:ceil(max([Lv.x]))); end
        if sp == 1 && ~isempty(hl), legend(ax, hl, lab, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'none'); end
        if sepG, title(ax, sprintf('%s (%d ch)', G(gr).name, numel(G(gr).ch)), 'Color', G(gr).col, 'FontSize', 10, 'Interpreter', 'none'); end
      end
      if sepG && sameYLim, linkaxes(axAll(:, sp), 'xy'); end         % same scale in every group (per column)
    end
    % sgtitle(fh, sprintf('%s: %s vs number of pulses (mean +/- SEM over channels); window first pulse ... last + %g ms', ...
            % strrep(short, '_', '\_'), lower(figs{fi, 1}), tail_ms), 'FontSize', 10);
end
if Plot_ChannelFigure
    fh = figure('Name', sprintf('%s: channels', short), 'Color', 'w', 'Position', [100 100 600 750]);
    ax = axes(fh);  hold(ax, 'on');
    for g = 1:nG
        for bb = 1:size(G(g).bands, 1)
            r = G(g).bands(bb, :);  patch(ax, [-1e4 1e4 1e4 -1e4], [r(1) r(1) r(2) r(2)], colG(g,:), 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    hl = [];  lab = {};
    for p = 1:nP
        v = S(p).rate(:, end);
        hl(end+1) = plot(ax, v, y, mkP{min(p, 3)}, 'Color', [0.5 0.5 0.5], 'MarkerSize', 3); %#ok<SAGROW>
        lab{end+1} = sprintf('%s final level, all channels', Lv(p).name); %#ok<SAGROW>
        inc = find(incl);
        plot(ax, v(inc), y(inc), mkP{min(p, 3)}, 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'HandleVisibility', 'off');
    end
    xline(ax, 0, ':'); ylim(ax, [min(y) max(y)]);
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
