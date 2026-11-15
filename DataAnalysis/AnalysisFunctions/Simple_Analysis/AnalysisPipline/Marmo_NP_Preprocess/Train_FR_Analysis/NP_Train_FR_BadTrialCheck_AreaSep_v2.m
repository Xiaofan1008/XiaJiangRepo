%% ============================================================
% NP_Train_FR_BadTrialCheck_v2.m
% NP probe | pulse TRAIN sessions | FIRING-RATE analysis | BAD-TRIAL CHECK PER AREA (V1 / V2) | v2
% (v1 = NP_Train_FR_BadTrialCheck_v1: the same rules on all channels together; v1 is unchanged)
% READ-ONLY: prints and plots, NOTHING IS SAVED.
%
% WHAT IS NEW IN v2
%   The channels are split by AREA (labels from NP_GetAreas.m in the RF folder: an area file of NP_AreaEntry_v1 or
%   ranges typed here). Every rule, every display value and every flag is computed SEPARATELY for each area in areaList
%   (default V1 and V2), using only that area's channels. 'bad' channels are left out (excludeBad = true); border /
%   outside channels are not in V1 or V2 (add 'border' to areaList to check them as their own group).
%   Figures: every v1 figure is drawn once PER AREA (area name in the title), plus an AREA SCATTER: each trial's
%   ratio in the first area vs the second area. A trial far out along the diagonal is a GLOBAL outlier (both areas);
%   far out on one axis only = an AREA-SPECIFIC outlier.
%   Printed: the v1 table per area, a comparison per condition (which areas flagged each trial), and two paste-ready
%   lists: (a) the UNION of the areas in the v1 entry format (NP_Train_FR_BadTrialEntry_v1), (b) per area, with an
%   area column (for an area-aware entry).
%
% CONDITIONS CHECKED
%   every train level (sim, seq; seqLevels = 'pair' -> [n n] only, 'all' -> also the half-steps) and every VISUAL
%   condition. Blank trials are not checked (they still count in the session-wide baseline reference).
%
% WINDOWS (per trial)
%   train level : first pulse (+trainStart_ms) ... last pulse + tail_ms   (relative to the trigger)
%   visual      : visWinMs after the photodiode onset (trigger if no photodiode onset)
%   baseline    : baseWinMs relative to the trigger
%
% WHAT THE RULES TEST (per area; as v1, but the channels are the ACTIVE channels OF THAT AREA)
%   For every channel, a trial's spike count is divided by that channel's MEAN count (denominator >= minMeanCount):
%     RESP ratio  = mean over the area's channels of  count(trial) / mean count of the channel over the condition's trials
%     MAXCH ratio = largest single-channel ratio of the trial (in the area)
%     BASE ratio  = mean over the area's channels of baseline count(trial) / the channel's mean baseline count (session)
%   R+ / R-  RESP  : |ratio - median of the condition| > respMAD_k  x max(1.4826 x MAD, respMadFloor)
%   M        MAXCH : ratio - median of the condition   > maxChMAD_k x max(1.4826 x MAD, maxChMadFloor)   (high side)
%   B        BASE  : ratio - median of the SESSION     > baseMAD_k  x max(1.4826 x MAD, baseMadFloor)    (high side)
%
% DISPLAY VALUES (per area): RESPONSE / MAXCH over the area's RESPONSIVE channels (NP_RespChn file or z > respZ),
%   BASELINE over the area's active channels.
% Trial numbers = number WITHIN the condition (TT.trial.rep); absolute trial in brackets.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile  = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery_randomEvent.mat';            % '' = <short>.NP_Spikes_TrainRecovery.mat

% ---- channels ----
autoRateMin  = 20;          % spikes/s over the whole window: ACTIVE channels (used by the rules)
respZ        = 2;           % responsive channel (figures only): (window rate - baseline)/SE over trials > respZ
useRespChn   = true;        % use <short>.NP_RespChn.mat for the figures when it exists
minMeanCount = 1;           % spikes: smallest denominator of a channel's ratio (stops near-silent channels from exploding)

% ---- windows ----
tail_ms       = 15;         % train window = first pulse + trainStart_ms ... last pulse + tail_ms
trainStart_ms = 0;
visWinMs      = [0 200];    % visual window, ms after the photodiode onset
baseWinMs     = [-100 -10];  % baseline window, ms relative to the trigger
seqLevels     = 'pair';     % 'pair' ([n n]) | 'all' (also the half-steps)
checkVisual   = true;       % also check the visual conditions

% ---- rules (PLACEHOLDERS) ----
Use_RespRule  = true;   respMAD_k  = 5;   respMadFloor  = 0.05;   Resp_FlagLow = true;
Use_MaxChRule = true;   maxChMAD_k = 5;   maxChMadFloor = 0.5;
Use_BaseRule  = true;   baseMAD_k  = 5;   baseMadFloor  = 0.05;
minTrialsForMAD         = 5;      % R and M need at least this many trials in the condition
Warn_Exclusion_Fraction = 0.30;   % warn when a condition loses more than this share of its trials

% ---- areas (labels from NP_GetAreas.m) ----
areaFuncDir = '/Volumes/LaCie/Marmoset_Data/AnalysisFunctions/NP_Analysis/RF_Detection/NP_RF';   % where NP_GetAreas.m is
areaSource  = 'file';       % 'file' (area file of NP_AreaEntry_v1) | 'manual' (ranges below)
areaFile    = '/Volumes/LaCie/Marmoset_Data/CJ287/CJ287_NP_Areas_20261110_0506.mat';
areasManual = { ...         % used when areaSource = 'manual': rows {from_um to_um label}
    0       2600   'V1';
    2600    2800   'border';
    2800    5400   'V2';
    5400    6000   'V2';
    6000    7660   'outside';
    };
badManual   = [83];         % used when areaSource = 'manual'
excludeBad  = true;         % true = channels labelled 'bad' are never used
areaList    = {'V1', 'V2'}; % areas checked separately (add 'border' to check it as its own group)
Plot_AreaScatter = true;    % trial ratio in areaList{1} vs areaList{2}

% ---- figures ----
Trial_Measure       = [1 2 3];    % panel figures: 1 RESPONSE | 2 MAXCH | 3 BASELINE
Plot_PanelFigures   = true;
Plot_SummaryFigure  = true;
Plot_SessionBaseline = true;
Plot_RatioFigure    = false;      % the ratios the rules tested (one panel figure)
nPanelCols          = 5;
Label_All_Trials    = true;       % true = number on EVERY trial (flagged ones red + rule); false = only the flagged trials
Label_FontSize      = 9;          % font size of the flagged-trial labels
Series_Shift        = 0.30;       % sim / seq shifted left / right (trial-number units)
Marker_Size_Kept    = 6;
Marker_Size_Flagged = 10;
%% =================================================

assert(any(strcmp(seqLevels, {'pair', 'all'})), 'seqLevels: ''pair'' or ''all''.');

%% ---- load -----------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
if isempty(spikeFile), spikeFile = [short '.NP_Spikes_TrainRecovery.mat']; end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
assert(isfile(spikeFile), 'Spike file not found: %s', spikeFile);
SP = load(spikeFile, 'sp_clipped', 'sp_trial', 'nSpk', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed');
fs = SP.FS;  nCh = numel(SP.sp_clipped);  nTr = SP.nTrialsProcessed;
tr = TT.trial;  trialCond = tr.cond(1:nTr);  rep = tr.rep(1:nTr);
trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;
kinds = {TT.cond.kind};
RespC = [];  rf = fullfile(recDir, [short '.NP_RespChn.mat']);
if useRespChn && isfile(rf), Rf = load(rf, 'RespChn'); RespC = Rf.RespChn; end

%% ---- checked conditions --------------------------------------------------------------------------------------
Q = struct('cond', {}, 'pat', {}, 'x', {}, 'isVis', {});          % one entry per checked condition
patName = @(k) ternary_(strcmp(k, 'sim'), 'sim', ternary_(strcmp(k, 'AB'), 'seq A->B', ternary_(strcmp(k, 'BA'), 'seq B->A', k)));
for p = 1:numel(TT.pattern)
    for c = TT.pattern(p).conds(:)'
        C = TT.cond(c);  nsp = C.nsp;  nsp(C.amp <= 0) = 0;
        if strcmp(C.kind, 'sim'), x = max(nsp);
        elseif nsp(1) == nsp(2), x = nsp(1);
        else, if strcmp(seqLevels, 'pair'), continue; end, x = max(nsp) - 0.5;
        end
        if ~any(trialCond == c), continue; end
        Q(end+1) = struct('cond', c, 'pat', patName(C.kind), 'x', x, 'isVis', false); %#ok<SAGROW>
    end
end
if checkVisual
    for c = find(strcmp(kinds, 'visual'))
        if any(trialCond == c), Q(end+1) = struct('cond', c, 'pat', 'visual', 'x', NaN, 'isVis', true); end %#ok<SAGROW>
    end
end
nQ = numel(Q);
assert(nQ > 0, 'No condition to check.');

%% ---- counts per channel and trial ---------------------------------------------------------------------------------
winOf = nan(nTr, 2);                                        % analysis window of every trial (ms after the trigger)
for q = 1:nQ
    ks = find(trialCond == Q(q).cond)';
    for t = ks
        if Q(q).isVis
            on = tr.diodeOnMs(t);  if isnan(on), on = 0; end
            winOf(t, :) = on + visWinMs;
        else
            winOf(t, :) = [tr.firstPulseMs(t) + trainStart_ms, tr.lastPulseMs(t) + tail_ms];
        end
    end
end
Cw = nan(nCh, nTr);  Cb = nan(nCh, nTr);                    % spikes in the window / in the baseline
for ch = find(SP.goodCh(:))'
    k = SP.sp_trial{ch}(:);  rel = SP.sp_clipped{ch}(:) - trigMs(k);
    Cb(ch, :) = accumarray(k(rel >= baseWinMs(1) & rel < baseWinMs(2)), 1, [nTr 1])';
    inW = rel >= winOf(k, 1) & rel < winOf(k, 2);
    Cw(ch, :) = accumarray(k(inW), 1, [nTr 1])';
end
Cw(:, isnan(winOf(:, 1))) = NaN;
baseDur = diff(baseWinMs)/1e3;  winDur = diff(winOf, 1, 2)'/1e3;            % s
Csub = Cw - bsxfun(@times, Cb/baseDur, winDur);                               % baseline-subtracted count

%% ---- channel sets -----------------------------------------------------------------------------------------------------
rateAll = sum(SP.nSpk(:, 1:nTr), 2)/nTr/((SP.preMs + SP.postMs)/1e3);
act = SP.goodCh(:) & rateAll >= autoRateMin;
zOf = @(ks) mean(Csub(:, ks), 2) ./ (std(Csub(:, ks), 0, 2)/sqrt(numel(ks)));
respTrain = false(nCh, 1);  respVis = false(nCh, 1);
if ~isempty(RespC)
    for q = 1:nQ
        if Q(q).isVis, respVis(RespC.respChn{Q(q).cond}) = true; else, respTrain(RespC.respChn{Q(q).cond}) = true; end
    end
    respSrc = 'NP_RespChn file';
else
    for p = 1:numel(TT.pattern)
        fc = TT.pattern(p).finalCond;  ks = find(trialCond == fc);
        if numel(ks) > 1, respTrain = respTrain | zOf(ks) > respZ; end
    end
    for q = find([Q.isVis])
        ks = find(trialCond == Q(q).cond);
        if numel(ks) > 1, respVis = respVis | zOf(ks) > respZ; end
    end
    respSrc = sprintf('automatic, z > %g', respZ);
end
respTrain = respTrain & act;  respVis = respVis & act;

%% ---- areas -----------------------------------------------------------------------------------------------------------
assert(any(strcmpi(areaSource, {'file', 'manual'})), 'areaSource must be ''file'' or ''manual'' in v2 (v1 checks all channels together).');
addpath(areaFuncDir);
[aLab, aInfo] = NP_GetAreas(areaSource, areaFile, areasManual, badManual, TT.chan.y_um(:));
labUse = aLab;  if ~excludeBad, labUse = aInfo.rangeLabel; end
nA = numel(areaList);

%% ---- per area: per-trial values (ratios for the rules, display values) and the rules -----------------------------------
isChk = false(nTr, 1);  for q = 1:nQ, isChk(trialCond == Q(q).cond) = true; end
A = struct('name', {}, 'act', {}, 'respTrain', {}, 'respVis', {}, 'Rr', {}, 'Mr', {}, 'Br', {}, 'Rd', {}, 'Md', {}, 'Bd', {}, ...
           'flagR', {}, 'flagRlow', {}, 'flagM', {}, 'flagB', {}, 'flagAny', {}, 'mdB', {}, 'scB', {}, 'codeOf', {});
for ia = 1:nA
    inA = strcmp(labUse, areaList{ia});
    actA = act & inA;  rTA = respTrain & inA;  rVA = respVis & inA;
    if ~any(actA), warning('NPBT:area', 'Area %s has no active channel: nothing is flagged there.', areaList{ia}); end
    Rr = nan(nTr, 1);  Mr = nan(nTr, 1);  Br = nan(nTr, 1);       % ratios
    Rd = nan(nTr, 1);  Md = nan(nTr, 1);  Bd = nan(nTr, 1);       % display
    if any(actA)
        bMean = mean(Cb(actA, :), 2);                              % each channel's mean baseline count, whole session
        Br(:) = mean(bsxfun(@rdivide, Cb(actA, :), max(bMean, minMeanCount)), 1)';
        Bd(:) = mean(Cb(actA, :), 1)'/baseDur;
        for q = 1:nQ
            ks = find(trialCond == Q(q).cond)';
            m = mean(Cw(actA, ks), 2);
            rat = bsxfun(@rdivide, Cw(actA, ks), max(m, minMeanCount));
            Rr(ks) = mean(rat, 1)';  Mr(ks) = max(rat, [], 1)';
            rs = ternary_(Q(q).isVis, rVA, rTA);
            if any(rs), Rd(ks) = mean(Csub(rs, ks), 1)'; Md(ks) = max(Csub(rs, ks), [], 1)'; end
        end
    end
    fR = false(nTr, 1);  fRl = false(nTr, 1);  fM = false(nTr, 1);  fB = false(nTr, 1);  mdB = NaN;  scB = NaN;
    for q = 1:nQ
        ks = find(trialCond == Q(q).cond);
        if numel(ks) < minTrialsForMAD || ~any(actA), continue; end
        if Use_RespRule
            [md, sc] = madScale_(Rr(ks), respMadFloor);
            fR(ks) = Rr(ks) - md > respMAD_k*sc;
            if Resp_FlagLow, fRl(ks) = md - Rr(ks) > respMAD_k*sc; end
        end
        if Use_MaxChRule
            [md, sc] = madScale_(Mr(ks), maxChMadFloor);
            fM(ks) = Mr(ks) - md > maxChMAD_k*sc;
        end
    end
    if Use_BaseRule && any(actA)
        [mdB, scB] = madScale_(Br, baseMadFloor);                  % whole session (every trial, blank included)
        fB = Br - mdB > baseMAD_k*scB;
    end
    fB = fB & isChk;                                               % blank trials are not flagged
    fAny = fR | fRl | fM | fB;
    A(ia) = struct('name', areaList{ia}, 'act', actA, 'respTrain', rTA, 'respVis', rVA, 'Rr', Rr, 'Mr', Mr, 'Br', Br, ...
                   'Rd', Rd, 'Md', Md, 'Bd', Bd, 'flagR', fR, 'flagRlow', fRl, 'flagM', fM, 'flagB', fB, 'flagAny', fAny, ...
                   'mdB', mdB, 'scB', scB, 'codeOf', @(t) codeStr_(t, fR, fRl, fM, fB));
end
flagUnion = false(nTr, 1);  for ia = 1:nA, flagUnion = flagUnion | A(ia).flagAny; end

%% ---- printout -----------------------------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP TRAIN FR BAD-TRIAL CHECK PER AREA (v2): %s\n', short);
fprintf('============================================================\n');
fprintf('spike file : %s\n%s%s\n', spikeFile, aInfo.text, ternary_(excludeBad && ~isempty(aInfo.badCh), sprintf(' | bad channels excluded: %s', mat2str(aInfo.badCh)), ''));
for ia = 1:nA
    fprintf('  %-7s: %3d ACTIVE channels for the rules (rate >= %g /s) | figures: %d responsive (trains), %d (visual) [%s]\n', ...
            A(ia).name, sum(A(ia).act), autoRateMin, sum(A(ia).respTrain), sum(A(ia).respVis), respSrc);
end
fprintf('windows    : train first pulse + %g ... last pulse + %g ms | visual %s ms after the photodiode | baseline %s ms\n', ...
        trainStart_ms, tail_ms, mat2str(visWinMs), mat2str(baseWinMs));
fprintf('rules      : R %s (k %g, floor %g, low side %d) | M %s (k %g, floor %g) | B %s (k %g, floor %g, session-wide) | MAD needs >= %d trials\n', ...
        onoff_(Use_RespRule), respMAD_k, respMadFloor, Resp_FlagLow, onoff_(Use_MaxChRule), maxChMAD_k, maxChMadFloor, ...
        onoff_(Use_BaseRule), baseMAD_k, baseMadFloor, minTrialsForMAD);
patList = unique({Q.pat}, 'stable');
for ia = 1:nA
    a = A(ia);
    fprintf('\n=============== %s ===============\n', a.name);
    if Use_BaseRule, fprintf('session baseline ratio: median %.3f, limit %.3f\n', a.mdB, a.mdB + baseMAD_k*a.scB); end
    fprintf('  pattern     cond     x  trials  median resp  flagged trials (rel [abs]: rules)\n');
    nFlagPat = zeros(size(patList));  nTrPat = zeros(size(patList));
    for q = 1:nQ
        ks = find(trialCond == Q(q).cond)';  fl = ks(a.flagAny(ks));
        s = strjoin(arrayfun(@(t) sprintf('%d [%d]: %s', rep(t), t, a.codeOf(t)), fl, 'UniformOutput', false), ', ');
        if isempty(s), s = '-'; end
        fprintf('  %-10s  %4d  %4s  %6d  %11.2f  %s\n', Q(q).pat, Q(q).cond, num2strNaN_(Q(q).x), numel(ks), median(a.Rd(ks), 'omitnan'), s);
        if numel(fl)/numel(ks) > Warn_Exclusion_Fraction
            fprintf('     WARNING: %d of %d trials flagged in this condition\n', numel(fl), numel(ks));
        end
        ip = find(strcmp(patList, Q(q).pat));  nFlagPat(ip) = nFlagPat(ip) + numel(fl);  nTrPat(ip) = nTrPat(ip) + numel(ks);
    end
    fprintf('exclusion rate per pattern (%s):\n', a.name);
    for ip = 1:numel(patList), fprintf('  %-10s %4d of %4d trials (%.1f %%)\n', patList{ip}, nFlagPat(ip), nTrPat(ip), 100*nFlagPat(ip)/nTrPat(ip)); end
end
% which areas flagged each trial
fprintf('\n=============== AREAS COMPARED (trial within condition: areas that flagged it) ===============\n');
nBoth = 0;  nOne = zeros(1, nA);
for q = 1:nQ
    ks = find(trialCond == Q(q).cond)';  fl = ks(flagUnion(ks));
    if isempty(fl), continue; end
    parts = cell(1, numel(fl));
    for j = 1:numel(fl)
        t = fl(j);  who = areaList(arrayfun(@(ia) A(ia).flagAny(t), 1:nA));
        parts{j} = sprintf('%d [%s]', rep(t), strjoin(who, '+'));
        if numel(who) == nA, nBoth = nBoth + 1; elseif numel(who) == 1, nOne(strcmp(areaList, who{1})) = nOne(strcmp(areaList, who{1})) + 1; end
    end
    fprintf('  %-10s  %4d  %4s  %s\n', Q(q).pat, Q(q).cond, num2strNaN_(Q(q).x), strjoin(parts, ', '));
end
fprintf('  flagged in all %d areas: %d trials', nA, nBoth);
for ia = 1:nA, fprintf(' | only %s: %d', areaList{ia}, nOne(ia)); end
fprintf('\n');
% paste-ready lists
fprintf('\nPASTE-READY (a) UNION of the areas, format of NP_Train_FR_BadTrialEntry_v1 ([] = reviewed, none bad):\nentries = { ...\n');
for q = 1:nQ
    ks = find(trialCond == Q(q).cond)';  fl = rep(ks(flagUnion(ks)))';
    if isempty(fl), fls = '[]'; else, fls = mat2str(fl); end
    fprintf('    %-12s %5s, %3d, %s; ...\n', ['''' Q(q).pat ''','], num2strNaN_(ternary_(Q(q).isVis, 0, Q(q).x)), Q(q).cond, fls);
end
fprintf('    };\n');
fprintf('\nPASTE-READY (b) PER AREA, with an area column (area, pattern, n, cond, trials):\nentries = { ...\n');
for ia = 1:nA
    for q = 1:nQ
        ks = find(trialCond == Q(q).cond)';  fl = rep(ks(A(ia).flagAny(ks)))';
        if isempty(fl), fls = '[]'; else, fls = mat2str(fl); end
        fprintf('    %-7s %-12s %5s, %3d, %s; ...\n', ['''' A(ia).name ''','], ['''' Q(q).pat ''','], num2strNaN_(ternary_(Q(q).isVis, 0, Q(q).x)), Q(q).cond, fls);
    end
end
fprintf('    };\n');

%% ---- figures (one set per area) ------------------------------------------------------------------------------------------------
mNames = {'RESPONSE: baseline-subtracted count per channel (responsive channels)', ...
          'MAX CHANNEL: largest single-channel baseline-subtracted count (responsive channels)', ...
          'BASELINE: mean baseline rate (spikes/s, active channels)'};
mYLab = {'count per channel', 'count (one channel)', 'spikes/s'};
xs = unique([Q(~[Q.isVis]).x]);
panelsE = struct('title', {}, 'qs', {});  panelsV = panelsE;
for x = xs, panelsE(end+1) = struct('title', sprintf('%g pulse(s) per electrode', x), 'qs', find(~[Q.isVis] & [Q.x] == x)); end %#ok<SAGROW>
for q = find([Q.isVis]), panelsV(end+1) = struct('title', TT.cond(Q(q).cond).label, 'qs', q); end %#ok<SAGROW>
patCol = containers.Map({'sim', 'seq A->B', 'seq B->A', 'visual'}, {[0 0.45 0.74], [0.47 0.67 0.19], [0.47 0.67 0.19], [0.49 0.18 0.56]});
patMk  = containers.Map({'sim', 'seq A->B', 'seq B->A', 'visual'}, {'o', 's', 's', 'd'});
for ia = 1:nA
    a = A(ia);  dispV = {a.Rd, a.Md, a.Bd};  tag = sprintf('%s  [%s]', short, a.name);
    pArgs = {Q, trialCond, rep, a.flagAny, a.codeOf, patCol, patMk, nPanelCols, Series_Shift, Label_All_Trials, Label_FontSize, Marker_Size_Kept, Marker_Size_Flagged};
    if Plot_PanelFigures
        for mI = Trial_Measure(:)'
            if ~isempty(panelsE), drawPanels_(panelsE, dispV{mI}, pArgs{:}, sprintf('%s  ELECTRICAL  %s', tag, mNames{mI}), mYLab{mI}); end
            if ~isempty(panelsV), drawPanels_(panelsV, dispV{mI}, pArgs{:}, sprintf('%s  VISUAL  %s', tag, mNames{mI}), mYLab{mI}); end
        end
    end
    if Plot_RatioFigure
        if ~isempty(panelsE), drawPanels_(panelsE, a.Rr, pArgs{:}, sprintf('%s  ELECTRICAL  RESP ratio tested by the rules', tag), 'ratio'); end
        if ~isempty(panelsV), drawPanels_(panelsV, a.Rr, pArgs{:}, sprintf('%s  VISUAL  RESP ratio tested by the rules', tag), 'ratio'); end
    end
    if Plot_SummaryFigure
        fh = figure('Name', [tag ' summary'], 'Color', 'w', 'Position', [80 80 1100 300*numel(Trial_Measure)]);
        for r = 1:numel(Trial_Measure)
            mI = Trial_Measure(r);  ax = subplot(numel(Trial_Measure), 1, r);  hold(ax, 'on');  hl = [];  lab = {};
            for ip = 1:numel(patList)
                qq = find(strcmp({Q.pat}, patList{ip}) & ~[Q.isVis]);  if isempty(qq), continue; end
                xv = [];  yv = [];  fv = [];
                for q = qq
                    ks = find(trialCond == Q(q).cond);
                    xv = [xv; Q(q).x*ones(numel(ks), 1) + (ip - 1.5)*0.15]; yv = [yv; dispV{mI}(ks)]; fv = [fv; a.flagAny(ks)]; %#ok<AGROW>
                end
                hl(end+1) = plot(ax, xv, yv, patMk(patList{ip}), 'Color', patCol(patList{ip}), 'MarkerSize', 4); %#ok<SAGROW>
                plot(ax, xv(fv > 0), yv(fv > 0), 'o', 'Color', 'r', 'MarkerSize', 9, 'LineWidth', 1.2, 'HandleVisibility', 'off');
                lab{end+1} = patList{ip}; %#ok<SAGROW>
            end
            grid(ax, 'on');  ylabel(ax, sprintf('measure %d', mI));  if ~isempty(hl), legend(ax, hl, lab, 'Location', 'best'); end
            title(ax, mNames{mI}, 'FontSize', 9);
            if r == numel(Trial_Measure), xlabel(ax, 'pulses per electrode (red ring = flagged)'); end
        end
        sgtitle(fh, strrep(tag, '_', '\_'));
    end
end
if Plot_SessionBaseline
    fh = figure('Name', [short ' session baseline per area'], 'Color', 'w', 'Position', [100 100 1200 260*nA + 60]);
    for ia = 1:nA
        a = A(ia);
        ax = subplot(nA, 2, 2*ia - 1);  hold(ax, 'on');
        plot(ax, 1:nTr, a.Bd, '.', 'Color', [0.5 0.5 0.5]);  plot(ax, find(a.flagAny), a.Bd(a.flagAny), 'ro', 'MarkerSize', 6);
        grid(ax, 'on');  ylabel(ax, 'baseline (spikes/s)');  title(ax, sprintf('%s: baseline of every trial (red = flagged in %s)', a.name, a.name), 'FontSize', 9);
        ax = subplot(nA, 2, 2*ia);  hold(ax, 'on');
        plot(ax, 1:nTr, a.Br, '.', 'Color', [0.5 0.5 0.5]);  plot(ax, find(a.flagB), a.Br(a.flagB), 'ro', 'MarkerSize', 6);
        if Use_BaseRule && ~isnan(a.mdB), yline(ax, a.mdB, 'k:'); yline(ax, a.mdB + baseMAD_k*a.scB, 'r--'); end
        grid(ax, 'on');  ylabel(ax, 'BASE ratio');  title(ax, sprintf('%s: BASE ratio (red dashed = limit; red = flagged by B)', a.name), 'FontSize', 9);
        if ia == nA, xlabel(subplot(nA, 2, 2*ia - 1), 'trial (recording order)'); xlabel(ax, 'trial (recording order)'); end
    end
    sgtitle(fh, strrep(short, '_', '\_'));
end
if Plot_AreaScatter && nA >= 2
    a1 = A(1);  a2 = A(2);
    fh = figure('Name', sprintf('%s %s vs %s', short, a1.name, a2.name), 'Color', 'w', 'Position', [120 120 1100 520]);
    sc = {'Rr', 'RESP ratio (rules R+ / R-)'; 'Br', 'BASE ratio (rule B)'};
    for k = 1:2
        ax = subplot(1, 2, k);  hold(ax, 'on');  hl = [];  lab = {};
        v1a = a1.(sc{k, 1});  v2a = a2.(sc{k, 1});  chk = isChk & isfinite(v1a) & isfinite(v2a);
        for ip = 1:numel(patList)
            tt = find(chk & ismember(trialCond(:), [Q(strcmp({Q.pat}, patList{ip})).cond]));
            if isempty(tt), continue; end
            hl(end+1) = plot(ax, v1a(tt), v2a(tt), patMk(patList{ip}), 'Color', patCol(patList{ip}), 'MarkerSize', 4); %#ok<SAGROW>
            lab{end+1} = patList{ip}; %#ok<SAGROW>
        end
        both = chk & a1.flagAny & a2.flagAny;  only1 = chk & a1.flagAny & ~a2.flagAny;  only2 = chk & ~a1.flagAny & a2.flagAny;
        if any(both),  hl(end+1) = plot(ax, v1a(both),  v2a(both),  'o', 'Color', 'r', 'MarkerSize', 10, 'LineWidth', 1.5); lab{end+1} = 'flagged in both'; end %#ok<SAGROW>
        if any(only1), hl(end+1) = plot(ax, v1a(only1), v2a(only1), 'o', 'Color', [0.2 0.4 0.9], 'MarkerSize', 10, 'LineWidth', 1.5); lab{end+1} = ['only ' a1.name]; end %#ok<SAGROW>
        if any(only2), hl(end+1) = plot(ax, v1a(only2), v2a(only2), 'o', 'Color', [0.9 0.45 0.1], 'MarkerSize', 10, 'LineWidth', 1.5); lab{end+1} = ['only ' a2.name]; end %#ok<SAGROW>
        lim = [min([v1a(chk); v2a(chk)]) max([v1a(chk); v2a(chk)])];
        if ~isempty(lim) && all(isfinite(lim)) && diff(lim) > 0
            plot(ax, lim, lim, 'k:', 'HandleVisibility', 'off');  xlim(ax, lim + [-0.05 0.05]*diff(lim));  ylim(ax, lim + [-0.05 0.05]*diff(lim));
        end
        axis(ax, 'square');  grid(ax, 'on');  xlabel(ax, sprintf('%s  %s', a1.name, sc{k, 2}));  ylabel(ax, sprintf('%s  %s', a2.name, sc{k, 2}));
        if ~isempty(hl), legend(ax, hl, lab, 'Location', 'best', 'FontSize', 8); end
    end
    sgtitle(fh, sprintf('%s: every checked trial, %s vs %s (on the diagonal = both areas alike; off it = one area only)', ...
            strrep(short, '_', '\_'), a1.name, a2.name), 'FontSize', 10);
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function drawPanels_(panels, V, Q, trialCond, rep, flagAny, codeOf, patCol, patMk, nCols, shift, labAll, fsz, msK, msF, ttl, yLab)
nPan = numel(panels);  nC = min(nCols, nPan);  nR = ceil(nPan/nC);
fh = figure('Name', ttl, 'Color', 'w', 'Position', [40 40 min(1800, 340*nC) min(1000, 280*nR + 60)]);
hLeg = gobjects(0);  legNames = {};  ax1 = [];
for i = 1:nPan
    ax = subplot(nR, nC, i);  hold(ax, 'on');  if i == 1, ax1 = ax; end
    qs = panels(i).qs;  nRepMax = 1;  LB = struct('x', {}, 'y', {}, 'txt', {}, 'flag', {}, 'up', {});
    for j = 1:numel(qs)
        q = qs(j);  ks = find(trialCond == Q(q).cond)';  nRepMax = max(nRepMax, max(rep(ks)));
        off = (j - (numel(qs) + 1)/2)*shift;
        x = rep(ks)' + off;  v = V(ks)';  col = patCol(Q(q).pat);
        h = plot(ax, x, v, patMk(Q(q).pat), 'LineStyle', 'none', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', msK);
        if ~any(strcmp(legNames, Q(q).pat)), hLeg(end+1) = h; legNames{end+1} = Q(q).pat; end %#ok<AGROW>
        yline(ax, median(v, 'omitnan'), ':', 'Color', col, 'LineWidth', 1, 'HandleVisibility', 'off');
        fl = flagAny(ks)';
        if any(fl)
            hf = plot(ax, x(fl), v(fl), 'o', 'LineStyle', 'none', 'Color', 'r', 'MarkerSize', msF, 'LineWidth', 1.5);
            if ~any(strcmp(legNames, 'flagged')), hLeg(end+1) = hf; legNames{end+1} = 'flagged'; end %#ok<AGROW>
        end
        up = mod(j, 2) == 1;                                  % 1st series: labels ABOVE the points, 2nd series: BELOW
        for r = 1:numel(ks)
            if fl(r)
                LB(end+1) = struct('x', x(r), 'y', v(r), 'txt', sprintf('%d %s', rep(ks(r)), codeOf(ks(r))), 'flag', true, 'up', up); %#ok<AGROW>
            elseif labAll
                LB(end+1) = struct('x', x(r), 'y', v(r), 'txt', sprintf('%d', rep(ks(r))), 'flag', false, 'up', up); %#ok<AGROW>
            end
        end
    end
    % labels: above / below the point, after the axis range is known
    yl = ylim(ax);  dy = 0.05*diff(yl);  ylim(ax, yl + [-2.5 2.5]*dy);
    for b = 1:numel(LB)
        if LB(b).up, yy = LB(b).y + dy; va = 'bottom'; else, yy = LB(b).y - dy; va = 'top'; end
        if LB(b).flag
            text(ax, LB(b).x, yy, LB(b).txt, 'FontSize', fsz, 'FontWeight', 'bold', 'Color', 'r', 'BackgroundColor', 'w', ...
                 'EdgeColor', [1 0.7 0.7], 'Margin', 0.5, 'HorizontalAlignment', 'center', 'VerticalAlignment', va);
        else
            text(ax, LB(b).x, yy, LB(b).txt, 'FontSize', max(fsz - 2, 6), 'Color', [0.35 0.35 0.35], ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', va);
        end
    end
    xlim(ax, [0.5 nRepMax + 0.5]);  xticks(ax, 1:nRepMax);  grid(ax, 'on');  box(ax, 'on');
    title(ax, panels(i).title, 'FontSize', 10, 'Interpreter', 'none');
    if mod(i - 1, nC) == 0, ylabel(ax, yLab); end
    if i > nPan - nC, xlabel(ax, 'trial within condition'); end
end
if ~isempty(hLeg)
    legend(ax1, hLeg, legNames, 'Location', 'best', 'FontSize', 8, 'Interpreter', 'none');   % inside the first panel
end
sgtitle(fh, sprintf('%s   (dotted = median; numbers = trial within condition: sim / 1st series above, seq below; red = flagged + rule)', strrep(ttl, '_', '\_')), 'FontSize', 11);
end

function [md, sc] = madScale_(v, floorV)
v = v(~isnan(v));
if isempty(v), md = NaN; sc = NaN; return; end
md = median(v);  sc = max(1.4826*median(abs(v - md)), floorV);
end

function s = num2strNaN_(v)
if isnan(v), s = '-'; elseif v == round(v), s = sprintf('%d', v); else, s = sprintf('%.1f', v); end
end

function s = onoff_(b)
if b, s = 'on'; else, s = 'off'; end
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end


function s = codeStr_(t, fR, fRl, fM, fB)
% rule codes of trial t, e.g. 'R+B'
s = strjoin([repmat({'R+'}, 1, double(fR(t))), repmat({'R-'}, 1, double(fRl(t))), repmat({'M'}, 1, double(fM(t))), repmat({'B'}, 1, double(fB(t)))], '');
end