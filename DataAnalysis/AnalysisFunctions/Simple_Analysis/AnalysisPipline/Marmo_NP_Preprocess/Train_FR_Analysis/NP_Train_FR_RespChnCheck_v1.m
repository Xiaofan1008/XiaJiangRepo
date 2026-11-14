%% ============================================================
% NP_Train_FR_RespChnCheck_v1.m
% NP probe | pulse TRAIN sessions | FIRING-RATE analysis | RESPONSIVE-CHANNEL CHECK | v1   (NP version of NN_Train_FR_RespChnCheck_v1)
% Prints and plots. Saves NOTHING unless saveProposal = true: then only <short>.NP_RespChn_proposal.mat (the automatic
% lists and the per-channel numbers), which NP_Train_FR_RespChnEntry_v1 reads. The proposal is simply rewritten each time.
%
% RULE, ruleType = 'rate' (DEFAULT): FIRING RATE over the whole response window, paired over trials
%   Per trial: response rate = spikes in the window / its length (blanked parts of the window are not counted in the length)
%              window: electrical = first pulse ... last pulse + tail_ms; visual = visRespWinMs after the photodiode onset
%              baseline rate = spikes in baseWinMs (relative to the trigger) / its length
%              d = response rate - baseline rate (same trial: the brain state of the trial cancels out)
%   Per channel (one value): rateSDmode = 'paired' (default): z = mean(d) / (SD(d)/sqrt(n trials))
%                            rateSDmode = 'baselineSD'      : z = (mean response - mean baseline) / SD of the baseline rates
%   RESPONSIVE if z > rateZ.  z < -rateZ = SUPPRESSED (listed for information, not responsive).
%   The number of trials with d > 0 is printed next to it.
%
% RULE, ruleType = 'psth' (as NN: baseline mean + k x SD, applied to the PSTH of each channel)
%   PSTH bins of binMs, blank-aware (a bin counts only the trials that had data in it; bins blanked in more than
%   (1 - minDataFrac) of the trials are ignored). Baseline = the bins in baseWinMs (relative to the TRIGGER): mean mu,
%   SD sd of the bin rates. Threshold = mu + kSD x max(sd, sdFloor_Hz).
%   A channel RESPONDS if at least minConsec CONSECUTIVE bins in the response window are above the threshold.
%   Response window: electrical = respWinMs after the FIRST pulse; visual = visRespWinMs after the photodiode onset.
%   Latency = start of the first such run; peak = highest bin in the window.
%   Visual conditions use their own (wider) bin, visBinMs, because their window is long (many bins = many chances).
%
% FALSE-POSITIVE CHECK (fpCheck = true)
%   rate rule: the same test with a STIMULUS-FREE window of the same length as each trial's response window, placed just
%   before the baseline window (ending fpGapMs before it). psth rule: the same rule, same bins and same threshold,
%   applied to a STIMULUS-FREE window of the same length before the trigger
%   (electrical: starting at fpStartMs; visual: starting at fpStartVisMs; not overlapping the baseline window).
%   Channels that 'respond' there are false positives: the printed percentage is the chance rate of the current settings.
%
% WHICH CONDITION DECIDES (one channel set per pattern; sim and seq separate; every visual condition separate)
%   decideBy = 'first' (default, as NN): the first level with ONE pulse per electrode (sim n = 1, seq [1 1]):
%              real data, almost nothing copied, no saturation of long trains.
%   decideBy = 'final': the final level of each pattern (window after the FIRST pulse as well).
%   Information (does not change the recommendation): z of the whole-train count of the FINAL level vs baseline
%   (mean difference / SE over trials; window first pulse ... last pulse + tail_ms).
%
% CHANNELS   all channels processed in step 2 (goodCh). TESTED is decided PER PATTERN: with onlyActive = true a channel is
%            tested only if its mean rate in the RESPONSE window of the decision level is >= autoRateMin
%            (electrical: trigger (0) ... last pulse + tail_ms; visual: visRespWinMs after the photodiode; blanked time
%            not counted). A channel can therefore be tested for one pattern and 'not tested' for another.
% BAD TRIALS <short>.NP_BadTrials.mat is used when it exists (useBadTrials = true).
%
% FIGURES
%   1. ALONG THE PROBE: per decision condition, (peak - threshold) of every TESTED channel vs depth (green = responsive);
%      visual columns only when plotVisual = true;
%      last column: final-level train z vs depth. Depth groups (groupLines) can be drawn.
%   2. TILES (plotMode = 'psth' or 'rasterPsth'), for plotConds (default: the FINAL condition of each pattern, plus the
%      visual conditions when plotVisual = true). tileChannels = 'all': EVERY recorded (processed) channel, responsive or
%      not, tested or not ('tested' = channels tested in at least one pattern). Channels not tested for the page's pattern
%      have a grey title marked '(not tested)'.
%      tileChannels can also be a list of channel numbers, e.g. [2 68 70]. The page title gives the depth range of the page,
%      and a channel -> page index is printed.
%      One page per tilesPerPage channels (sorted by depth, top of the probe first); each tile = PSTH (black, smoothed
%      for display only) [+ raster dots behind it], dotted = baseline mean, dashed = threshold of THAT condition,
%      yellow = response window, red lines = pulses. Green title = recommended for that pattern.
%      plotMode = 'none' -> no tiles.
% PRINTED: per pattern the rule numbers of every channel (or only the responsive ones), counts per depth band, and a
%   paste-ready list for NP_Train_FR_RespChnEntry_v1.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile  = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery.mat';            % '' = <short>.NP_Spikes_TrainRecovery.mat (the unthinned recovery)
useBadTrials = true;        % leave out the trials in <short>.NP_BadTrials.mat when it exists
onlyActive   = true;        % true = per pattern, test only channels with a response-window rate >= autoRateMin; false = every processed channel
autoRateMin  = 20;          % spikes/s in the response window (0 ... last pulse + tail_ms; visual: visRespWinMs) of the decision level

% ---- rule ----
decideBy     = 'final';     % 'first' (sim n = 1, seq [1 1]; as NN) | 'final'
ruleType     = 'rate';      % 'rate' (default): whole-window firing rate vs baseline, paired over trials | 'psth': NN peak rule
rateZ        = 3;           % rate rule: responsive if z > rateZ (z < -rateZ = suppressed)
rateSDmode   = 'paired';    % rate rule: 'paired' (mean difference / SE over trials) | 'baselineSD' (difference / SD of baseline)
fpGapMs      = 10;          % rate rule: the stimulus-free window ends this many ms before the baseline window starts
checkVisual  = true;        % also test every visual condition
baseWinMs    = [-100 -10];   % baseline, ms relative to the trigger
respWinMs    = [1 20];      % psth rule only: electrical response window, ms after the FIRST pulse (rate rule: whole train)
visRespWinMs = [20 200];    % visual: response window, ms after the photodiode onset
binMs        = 1;           % PSTH bin (ms), electrical conditions
visBinMs     = 5;           % PSTH bin (ms), visual conditions (long window -> wider bins)
kSD          = 3;           % threshold = mu + kSD x max(SD, sdFloor_Hz)
minConsec    = 2;           % consecutive bins above the threshold
sdFloor_Hz   = 10;          % floor on the SD of the baseline bins (spikes/s)
minDataFrac  = 0.5;         % a bin is used only if at least this fraction of the trials had data in it
tail_ms      = 15;          % information: final-level train window = first pulse ... last pulse + tail_ms
fpCheck      = true;        % false-positive check: the rule on a stimulus-free window before the trigger
fpStartMs    = -200;        % electrical: start of the fake window (ms, trigger); its length = the response window
fpStartVisMs = -290;        % visual: start of the fake window (ms, trigger); keep start + length < baseline start

% ---- printout ----
printAllChannels = false;   % true = one line per channel; false = only the responsive channels + summary

% ---- proposal for NP_Train_FR_RespChnEntry_v1 ----
saveProposal = true;       % true = save <short>.NP_RespChn_proposal.mat in the recording folder (overwritten each run)

% ---- figures ----
Plot_ProbeFigure = true;
plotVisual   = false;       % true = the visual conditions also get columns in the probe figure and tile pages (many figures)
tileChannels = 'all';       % 'all' = every recorded (processed) channel | 'tested' = only the tested channels | a list, e.g. [2 68 70]
groupLines   = [2200 3000 4500 5800];   % depths (um) drawn as lines in the probe figure; [] = none
plotMode     = 'rasterPsth';      % 'none' | 'psth' | 'rasterPsth'   (tiles)
plotConds    = [15 45];          % conditions to tile; [] = the FINAL condition of each pattern (add visual ones by number)
tilesPerPage = 48;          % channels per page (8 x 6)
nTileCols    = 8;
xRange       = [-100 250];    % ms around the alignment point (electrical: trigger; visual: photodiode onset)
xRangeVis    = [-50 250];
curveSmooth_ms = 3;         % display only: Gaussian smoothing of the drawn PSTH (ms); the rule always uses raw bins
Tile_FontSize  = 7;
%% =================================================

assert(any(strcmp(decideBy, {'first', 'final'})), 'decideBy must be ''first'' or ''final''.');
assert(any(strcmp(plotMode, {'none', 'psth', 'rasterPsth'})), 'plotMode must be ''none'', ''psth'' or ''rasterPsth''.');

%% ---- load ---------------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
if isempty(spikeFile), spikeFile = [short '.NP_Spikes_TrainRecovery.mat']; end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
assert(isfile(spikeFile), 'Spike file not found: %s', spikeFile);
SP = load(spikeFile, 'sp_clipped', 'sp_trial', 'blank_ms', 'nSpk', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed');
fs = SP.FS;  nCh = numel(SP.sp_clipped);  nTr = SP.nTrialsProcessed;  tr = TT.trial;
trialCond = tr.cond(1:nTr);  rep = tr.rep(1:nTr);  trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;
y = TT.chan.y_um(:);  kindC = {TT.cond.kind};
useTr = true(nTr, 1);  bf = fullfile(recDir, [short '.NP_BadTrials.mat']);  badTxt = 'none (all trials kept)';
if useBadTrials && isfile(bf)
    Bf = load(bf, 'BadTrials');
    for c = 1:numel(TT.cond)
        b = Bf.BadTrials.badRel{c};
        if ~isempty(b), ks = find(trialCond == c); useTr(ks(b(b <= numel(ks)))) = false; end
    end
    badTxt = sprintf('%s (%d trials left out)', bf, sum(~useTr));
end
chTest = find(SP.goodCh(:))';              % candidates = every processed channel; the TESTED set is decided per pattern below

%% ---- patterns: decision and final conditions ------------------------------------------------------------------------
Pt = struct('name', {}, 'decCond', {}, 'finalCond', {}, 'isVis', {}, 'conds', {});   % conds = every level of the pattern
for p = 1:numel(TT.pattern)
    P = TT.pattern(p);  nm = ternary_(strcmp(P.kind, 'sim'), 'sim', ternary_(strcmp(P.kind, 'AB'), 'seq A->B', 'seq B->A'));
    dc = P.finalCond;
    if strcmp(decideBy, 'first')
        dc = NaN;
        for c = P.conds(:)'
            nsp = TT.cond(c).nsp;  nsp(TT.cond(c).amp <= 0) = 0;
            if all(nsp == 1), dc = c; break; end
        end
        assert(~isnan(dc), 'Pattern %s has no level with one pulse per electrode.', nm);
    end
    Pt(end+1) = struct('name', nm, 'decCond', dc, 'finalCond', P.finalCond, 'isVis', false, 'conds', P.conds(:)'); %#ok<SAGROW>
end
if checkVisual
    for c = find(strcmp(kindC, 'visual'))
        if any(trialCond == c), Pt(end+1) = struct('name', sprintf('visual %d', c), 'decCond', c, 'finalCond', c, 'isVis', true, 'conds', c); end %#ok<SAGROW>
    end
end
nP = numel(Pt);

%% ---- spikes per trial (relative to the trigger) for the tested channels ------------------------------------------------
Bt = cell(nCh, 1);
for ch = chTest, Bt{ch} = byTrial_(SP.sp_clipped{ch}, SP.sp_trial{ch}, trigMs, nTr); end
alignOf = @(c, ks) ternary_(strcmp(kindC{c}, 'visual'), nanTo0_(tr.diodeOnMs(ks)), cellfun(@(v) firstOr0_(v), tr.pulseMs(ks)));

%% ---- rule for every pattern and channel -------------------------------------------------------------------------------------
assert(any(strcmp(ruleType, {'rate', 'psth'})), 'ruleType must be ''rate'' or ''psth''.');
assert(any(strcmp(rateSDmode, {'paired', 'baselineSD'})), 'rateSDmode must be ''paired'' or ''baselineSD''.');
isRate = strcmp(ruleType, 'rate');  fpClipped = false;
R = struct('resp', {}, 'mu', {}, 'sd', {}, 'thr', {}, 'peak', {}, 'lat', {}, 'zFinal', {}, 'fp', {}, 'bin', {}, 'fpWin', {}, ...
           'rBase', {}, 'rResp', {}, 'z', {}, 'nPos', {}, 'nT', {}, 'supp', {}, 'win', {}, 'tested', {}, 'actRate', {});
for p = 1:nP
    c = Pt(p).decCond;  ks = find(trialCond == c & useTr);  off = alignOf(c, ks);
    w = ternary_(Pt(p).isVis, visRespWinMs, respWinMs);
    bm = ternary_(Pt(p).isVis, visBinMs, binMs);
    fpW = ternary_(Pt(p).isVis, fpStartVisMs, fpStartMs) + [0 diff(w)];
    res = struct('resp', false(nCh, 1), 'mu', nan(nCh, 1), 'sd', nan(nCh, 1), 'thr', nan(nCh, 1), 'peak', nan(nCh, 1), ...
                 'lat', nan(nCh, 1), 'zFinal', nan(nCh, 1), 'fp', false(nCh, 1), 'bin', bm, 'fpWin', fpW, ...
                 'rBase', nan(nCh, 1), 'rResp', nan(nCh, 1), 'z', nan(nCh, 1), 'nPos', nan(nCh, 1), 'nT', nan(nCh, 1), ...
                 'supp', false(nCh, 1), 'win', [], 'tested', false(nCh, 1), 'actRate', nan(nCh, 1));
    kf = find(trialCond == Pt(p).finalCond & useTr);
    if isRate                                                       % per-trial windows (relative to the trigger)
        Wr = rateWin_(tr, Pt(p).isVis, ks, tail_ms, visRespWinMs);
        L = diff(Wr, 1, 2);  s0 = baseWinMs(1) - fpGapMs - L;       % stimulus-free window right before the baseline
        if any(s0 < -SP.preMs), fpClipped = true; s0 = max(s0, -SP.preMs); end
        Wfp = [s0, min(s0 + L, baseWinMs(1) - fpGapMs)];
        res.win = Wr;  res.fpWin = [min(Wfp(:,1)) max(Wfp(:,2))];
    end
    % activity: mean rate in the response window of the decision level (electrical from the trigger, 0 ms)
    Wa = rateWin_(tr, Pt(p).isVis, ks, tail_ms, visRespWinMs);  if ~Pt(p).isVis, Wa(:, 1) = 0; end
    for ch = chTest
        [~, res.actRate(ch)] = rateCh_(Bt{ch}, SP.blank_ms(ch, :), ks, Wa, baseWinMs, rateSDmode);
        res.tested(ch) = ~onlyActive || res.actRate(ch) >= autoRateMin;
    end
    for ch = find(res.tested)'
        if isRate
            [mB, mR, z, nPos, nT] = rateCh_(Bt{ch}, SP.blank_ms(ch, :), ks, Wr, baseWinMs, rateSDmode);
            res.rBase(ch) = mB;  res.rResp(ch) = mR;  res.z(ch) = z;  res.nPos(ch) = nPos;  res.nT(ch) = nT;
            res.resp(ch) = z > rateZ;  res.supp(ch) = z < -rateZ;
            if fpCheck
                [~, ~, zf] = rateCh_(Bt{ch}, SP.blank_ms(ch, :), ks, Wfp, baseWinMs, rateSDmode);
                res.fp(ch) = zf > rateZ;
            end
        else
            [mu, sd, thr, pk, lat, ok] = ruleCh_(Bt{ch}, SP.blank_ms(ch, :), ks, off, baseWinMs, w, bm, kSD, sdFloor_Hz, minConsec, minDataFrac);
            if fpCheck                                              % same rule on a stimulus-free window (relative to the trigger)
                [~, ~, ~, ~, ~, res.fp(ch)] = ruleCh_(Bt{ch}, SP.blank_ms(ch, :), ks, zeros(size(off)), baseWinMs, fpW, bm, kSD, sdFloor_Hz, minConsec, minDataFrac);
            end
            res.resp(ch) = ok;  res.mu(ch) = mu;  res.sd(ch) = sd;  res.thr(ch) = thr;  res.peak(ch) = pk;  res.lat(ch) = lat;
        end
        res.zFinal(ch) = trainZ_(Bt{ch}, kf, tr, baseWinMs, tail_ms, Pt(p).isVis, visRespWinMs);
    end
    R(p) = res;
end

%% ---- printout ------------------------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP TRAIN FR RESPONSIVE-CHANNEL CHECK (v1): %s\n', short);
fprintf('============================================================\n');
fprintf('spike file : %s\nbad trials : %s\nchannels   : %d processed; tested per pattern: %s\n', spikeFile, badTxt, numel(chTest), ...
        ternary_(onlyActive, sprintf('rate >= %g spikes/s in the response window of the decision level (0 ... last pulse + %g ms; visual %s ms)', ...
                 autoRateMin, tail_ms, mat2str(visRespWinMs)), 'all processed'));
if isRate
    fprintf('rule       : FIRING RATE, window first pulse ... last pulse + %g ms (visual %s ms after the photodiode) vs baseline %s ms;\n', ...
            tail_ms, mat2str(visRespWinMs), mat2str(baseWinMs));
    fprintf('             responsive if z > %g (%s); z < -%g = suppressed\n', rateZ, ...
            ternary_(strcmp(rateSDmode, 'paired'), 'z = mean(response - baseline) / SE over trials', 'z = (mean response - mean baseline) / SD of the baseline rates'), rateZ);
    if fpCheck
        fprintf('false-pos. : the same test on a stimulus-free window of the same length, ending %g ms before the baseline window\n', fpGapMs);
        if fpClipped, warning('NPRC:fpclip', 'Some stimulus-free windows had to be shortened to fit after -%g ms (start of the detection window).', SP.preMs); end
    end
else
fprintf('rule       : %g ms bins (visual %g ms), baseline %s ms, threshold = mu + %g x max(SD, %g Hz), >= %d consecutive bins in %s ms after the first pulse (visual: %s ms after the photodiode)\n', ...
        binMs, visBinMs, mat2str(baseWinMs), kSD, sdFloor_Hz, minConsec, mat2str(respWinMs), mat2str(visRespWinMs));
end
if fpCheck && ~isRate
    fprintf('false-pos. : the same rule on a stimulus-free window (electrical %s ms, visual %s ms, relative to the trigger)\n', ...
            mat2str(fpStartMs + [0 diff(respWinMs)]), mat2str(fpStartVisMs + [0 diff(visRespWinMs)]));
    if fpStartVisMs + diff(visRespWinMs) > baseWinMs(1) || fpStartMs + diff(respWinMs) > baseWinMs(1)
        warning('NPRC:fp', 'The false-positive window overlaps the baseline window: move fpStartMs / fpStartVisMs earlier.');
    end
end
fprintf('decides    : %s level\n', decideBy);
bandsE = [-Inf, groupLines(:)', Inf];
for p = 1:nP
    r = R(p);  rc = find(r.resp)';  chT = find(r.tested)';          % tested channels of THIS pattern
    fprintf('\n--- %s: decision cond %d (%s, %d trials) | final cond %d ---\n', Pt(p).name, Pt(p).decCond, TT.cond(Pt(p).decCond).label, ...
            sum(trialCond == Pt(p).decCond & useTr), Pt(p).finalCond);
    fprintf('  responsive: %d of %d tested channels (%d of %d processed below %g spikes/s: not tested)', numel(rc), numel(chT), ...
            numel(chTest) - numel(chT), numel(chTest), autoRateMin);
    if fpCheck
        fprintf('   |  FALSE-POSITIVE CHECK: %d of %d (%.1f %%) "respond" in the stimulus-free window', sum(r.fp(chT)), numel(chT), ...
                100*mean(r.fp(chT)));
    end
    fprintf('\n');
    for b = 1:numel(bandsE) - 1
        inB = chT(y(chT) >= bandsE(b) & y(chT) < bandsE(b + 1));
        fprintf('    depth %5s-%5s um: %3d of %3d\n', num2s_(bandsE(b)), num2s_(bandsE(b + 1)), sum(r.resp(inB)), numel(inB));
    end
    if printAllChannels, chP = chT; else, chP = rc; end
    if isRate && ~isempty(chP)
        fprintf('   ch   y(um)  base(Hz)  resp(Hz)  diff(Hz)       z   trials d>0   RESP\n');
        [~, o] = sort(y(chP), 'descend');
        for ch = chP(o)
            fprintf('  %3d  %6d  %8.1f  %8.1f  %8.1f  %6.1f      %2d/%2d     %s\n', ch, y(ch), r.rBase(ch), r.rResp(ch), ...
                    r.rResp(ch) - r.rBase(ch), r.z(ch), r.nPos(ch), r.nT(ch), ternary_(r.resp(ch), 'YES', '-'));
        end
    end
    if isRate
        sc = find(r.supp)';
        fprintf('  suppressed (z < -%g, information only): %d channels %s\n', rateZ, numel(sc), chList_(sc));
    end
    if ~isRate && ~isempty(chP)
        fprintf('   ch   y(um)  base(Hz)  SD(Hz)  thr(Hz)  peak(Hz)  latency(ms)  final z  RESP\n');
        [~, o] = sort(y(chP), 'descend');
        for ch = chP(o)
            fprintf('  %3d  %6d  %8.1f  %6.1f  %7.1f  %8.1f  %10s  %7.1f   %s\n', ch, y(ch), r.mu(ch), r.sd(ch), r.thr(ch), r.peak(ch), ...
                    num2s_(r.lat(ch)), r.zFinal(ch), ternary_(r.resp(ch), 'YES', '-'));
        end
    end
end
fprintf('\nPASTE-READY for NP_Train_FR_RespChnEntry_v1 (one set per pattern; edit before saving):\nentries = { ...\n');
for p = 1:nP
    rc = find(R(p).resp)';
    if Pt(p).isVis
        fprintf('    %-12s %3d, %s; ...\n', '''visual'',', Pt(p).decCond, chList_(rc));
    else
        fprintf('    %-12s %5s, %s; ...\n', ['''' Pt(p).name ''','], '''all''', chList_(rc));
    end
end
fprintf('    };\n');

%% ---- proposal file for NP_Train_FR_RespChnEntry_v1 (only when saveProposal = true) ------------------------------------------
if saveProposal
    Prop = struct();
    Prop.short = short;  Prop.recDir = recDir;  Prop.spikeFile = spikeFile;  Prop.badTrials = badTxt;
    Prop.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');  Prop.code = 'NP_Train_FR_RespChnCheck_v1';
    Prop.settings = struct('ruleType', ruleType, 'rateZ', rateZ, 'rateSDmode', rateSDmode, 'decideBy', decideBy, ...
        'baseWinMs', baseWinMs, 'visRespWinMs', visRespWinMs, 'tail_ms', tail_ms, 'onlyActive', onlyActive, ...
        'autoRateMin', autoRateMin, 'respWinMs', respWinMs, 'binMs', binMs, 'visBinMs', visBinMs, 'kSD', kSD, 'minConsec', minConsec);
    Prop.nCh = nCh;  Prop.y_um = y;  Prop.chTest = chTest;               % depth of every channel; processed channels (tested: per pattern)
    Prop.scoreName = ternary_(isRate, 'z (response vs baseline)', 'peak - threshold (Hz)');   % the number on the probe figures
    Prop.scoreThr  = ternary_(isRate, rateZ, 0);
    pr = struct('name', {}, 'decCond', {}, 'finalCond', {}, 'conds', {}, 'isVis', {}, 'auto', {}, 'score', {}, ...
                'rBase', {}, 'rResp', {}, 'nPos', {}, 'nT', {}, 'tested', {}, 'actRate', {});
    for p = 1:nP
        r = R(p);
        if isRate, sc = r.z; else, sc = r.peak - r.thr; end
        pr(p).name = Pt(p).name;  pr(p).decCond = Pt(p).decCond;  pr(p).finalCond = Pt(p).finalCond;  pr(p).conds = Pt(p).conds;
        pr(p).isVis = Pt(p).isVis;  pr(p).auto = find(r.resp)';  pr(p).score = sc;
        pr(p).rBase = r.rBase;  pr(p).rResp = r.rResp;  pr(p).nPos = r.nPos;  pr(p).nT = r.nT;
        pr(p).tested = find(r.tested)';  pr(p).actRate = r.actRate;      % tested channels of this pattern; response-window rate
    end
    Prop.pattern = pr;
    propFile = fullfile(recDir, [short '.NP_RespChn_proposal.mat']);
    save(propFile, 'Prop');
    fprintf('\nPROPOSAL saved (rewritten each run): %s\n', propFile);
end

%% ---- figure 1: along the probe -----------------------------------------------------------------------------------------------
if Plot_ProbeFigure
    pShow = find(~[Pt.isVis] | plotVisual);  nS = numel(pShow);
    fh = figure('Name', [short ' responsive channels along the probe'], 'Color', 'w', 'Position', [40 40 min(1800, 260*(nS + 1)) 750]);
    for ip = 1:nS
        p = pShow(ip);  ax = subplot(1, nS + 1, ip);  hold(ax, 'on');  r = R(p);
        if isRate, d = r.z; else, d = r.peak - r.thr; end
        chT = find(r.tested)';  nr = chT(~r.resp(chT));  rr = chT(r.resp(chT));
        plot(ax, d(nr), y(nr), '.', 'Color', [0.6 0.6 0.6]);  plot(ax, d(rr), y(rr), 'o', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 3);
        xline(ax, 0, 'k:');  for g = groupLines(:)', yline(ax, g, 'k-', 'Color', [0.7 0.7 0.7]); end
        if isRate, xline(ax, rateZ, 'g--'); xline(ax, -rateZ, 'b--'); end
        grid(ax, 'on');  ylim(ax, [min(y) max(y)]);  xlabel(ax, ternary_(isRate, 'z (response vs baseline)', 'peak - threshold (Hz)'));
        title(ax, sprintf('%s\ncond %d: %d responsive', Pt(p).name, Pt(p).decCond, sum(r.resp)), 'FontSize', 8, 'Interpreter', 'none');
        if ip == 1, ylabel(ax, 'depth from the tip (um)'); else, set(ax, 'YTickLabel', []); end
    end
    ax = subplot(1, nS + 1, nS + 1);  hold(ax, 'on');  cl = lines(nS);
    for ip = 1:nS, chT = find(R(pShow(ip)).tested)'; plot(ax, R(pShow(ip)).zFinal(chT), y(chT), '.', 'Color', cl(ip,:)); end
    xline(ax, 0, 'k:');  for g = groupLines(:)', yline(ax, g, 'k-', 'Color', [0.7 0.7 0.7]); end
    grid(ax, 'on');  ylim(ax, [min(y) max(y)]);  set(ax, 'YTickLabel', []);  xlabel(ax, 'final-level train z');
    legend(ax, {Pt(pShow).name}, 'FontSize', 7, 'Location', 'best', 'Interpreter', 'none');  title(ax, 'information: final level', 'FontSize', 8);
    sgtitle(fh, sprintf('%s: responsive channels (green) along the probe; decision = %s level', strrep(short, '_', '\_'), decideBy));
end

%% ---- figure 2: tiles ----------------------------------------------------------------------------------------------------------------
if ~strcmp(plotMode, 'none')
    if isempty(plotConds)
        plotConds = unique([Pt(~[Pt.isVis]).finalCond], 'stable');
        if plotVisual, plotConds = [plotConds, [Pt([Pt.isVis]).decCond]]; end
    end
    if smoothSigmaOK_(curveSmooth_ms, binMs)
        sb = curveSmooth_ms/binMs;  kx = -ceil(3*sb):ceil(3*sb);  ker = exp(-0.5*(kx/sb).^2);  ker = ker/sum(ker);
    else
        ker = 1;
    end
    if isnumeric(tileChannels)
        chTile = intersect(tileChannels(:)', find(SP.goodCh(:))');
        if numel(chTile) < numel(tileChannels), warning('NPRC:tilech', 'Some channels in tileChannels were not processed in step 2: left out.'); end
    elseif strcmp(tileChannels, 'all')
        chTile = find(SP.goodCh(:))';
    else
        chTile = find(any([R.tested], 2))';                    % 'tested' = tested in at least one pattern
    end
    for ch = setdiff(chTile, chTest)                           % spikes of channels that were not tested
        Bt{ch} = byTrial_(SP.sp_clipped{ch}, SP.sp_trial{ch}, trigMs, nTr);
    end
    [~, o] = sort(y(chTile), 'descend');  chP = chTile(o);
    nPage = ceil(numel(chP)/tilesPerPage);  nRows = ceil(tilesPerPage/nTileCols);
    fprintf('\nTILE PAGES (channels sorted by depth, top of the probe first):\n');
    for pg = 1:nPage
        chs = chP((pg - 1)*tilesPerPage + 1 : min(numel(chP), pg*tilesPerPage));
        fprintf('  page %d: %4d-%4d um | ch %s\n', pg, max(y(chs)), min(y(chs)), strjoin(arrayfun(@num2str, chs, 'UniformOutput', false), ' '));
    end
    % tile positions (normalized figure units): faster than subplot
    gL = 0.03;  gR = 0.01;  gB = 0.035;  gT = 0.075;
    tw = (1 - gL - gR)/nTileCols;  th = (1 - gB - gT)/nRows;
    for c = plotConds(:)'
        ks = find(trialCond == c & useTr);  off = alignOf(c, ks);  isV = strcmp(kindC{c}, 'visual');
        p = find([Pt.decCond] == c | [Pt.finalCond] == c, 1);            % pattern of this condition (recommendation, threshold)
        if isempty(p), pp = find(arrayfun(@(q) ~Pt(q).isVis && any(TT.pattern(q).conds == c), 1:nP), 1); else, pp = p; end
        if isempty(pp), pName = '-'; else, pName = Pt(pp).name; end
        if isV, pTag = 'VISUAL'; elseif isempty(pp), pTag = upper(kindC{c}); else, pTag = upper(pName); end   % big label on the page
        xr = ternary_(isV, xRangeVis, xRange);  w = ternary_(isV, visRespWinMs, respWinMs);
        ed = xr(1):binMs:xr(2);  ct = ed(1:end-1) + binMs/2;
        bmT = ternary_(isV, visBinMs, binMs);                       % bin of the RULE (threshold line)
        if isRate                                                   % the rate window of this condition, in plot (aligned) time
            Wt = rateWin_(tr, isV, ks, tail_ms, visRespWinMs);
            if isempty(ks), wT = w; else, wT = Wt(1, :) - off(1); end
        else
            wT = w;
        end
        pl = [];  if ~isempty(ks) && ~isempty(tr.pulseMs{ks(1)}), pl = tr.pulseMs{ks(1)} - off(1); end
        for pg = 1:nPage
            chs = chP((pg - 1)*tilesPerPage + 1 : min(numel(chP), pg*tilesPerPage));
            fh = figure('Name', sprintf('%s | cond %d | %s | page %d/%d', pTag, c, short, pg, nPage), 'Color', 'w', 'Position', [20 20 1700 980], 'Visible', 'off');
            for i = 1:numel(chs)
                ch = chs(i);  rI = ceil(i/nTileCols);  cI = i - (rI - 1)*nTileCols;
                ax = axes(fh, 'Position', [gL + (cI - 1)*tw + 0.012, gB + (nRows - rI)*th + 0.018, tw - 0.02, th - 0.045]);  hold(ax, 'on');
                if isRate                                                   % dotted = baseline rate, dashed = mean rate in the window
                    [mu, thr] = rateCh_(Bt{ch}, SP.blank_ms(ch, :), ks, Wt, baseWinMs, rateSDmode);
                else
                    [mu, ~, thr] = ruleCh_(Bt{ch}, SP.blank_ms(ch, :), ks, off, baseWinMs, w, bmT, kSD, sdFloor_Hz, minConsec, minDataFrac);
                end
                rt = psthAligned_(Bt{ch}, SP.blank_ms(ch, :), ks, off, ed, minDataFrac);
                rs = rt;  if numel(ker) > 1, rs = nanConv_(rt, ker); end
                yTop = max([rs(:); thr; 1], [], 'omitnan')*1.15;
                patch(ax, [wT(1) wT(2) wT(2) wT(1)], [0 0 yTop yTop], [1 1 0.6], 'EdgeColor', 'none');
                if strcmp(plotMode, 'rasterPsth')
                    xs = [];  ys = [];
                    for q = 1:numel(ks)
                        s_ = Bt{ch}{ks(q)} - off(q);  s_ = s_(s_ >= xr(1) & s_ < xr(2));
                        xs = [xs; s_(:)];  ys = [ys; (q/numel(ks))*yTop*0.95*ones(numel(s_), 1)]; %#ok<AGROW>
                    end
                    plot(ax, xs, ys, '.', 'Color', [0.65 0.65 0.65], 'MarkerSize', 2);
                end
                if ~isempty(pl)                                             % all pulses as ONE line object (fast)
                    px = [pl(:)'; pl(:)'; nan(1, numel(pl))];  py = repmat([0; yTop; NaN], 1, numel(pl));
                    plot(ax, px(:), py(:), 'r-', 'LineWidth', 0.3);
                end
                plot(ax, ct, rs, 'k-', 'LineWidth', 0.8);
                plot(ax, xr, [mu mu], 'b:');
                if isRate, plot(ax, wT, [thr thr], 'g--', 'LineWidth', 1); else, plot(ax, xr, [thr thr], 'b--'); end
                xlim(ax, xr);  ylim(ax, [0 yTop]);  set(ax, 'FontSize', Tile_FontSize, 'TickLength', [0.02 0.02]);
                isR = ~isempty(pp) && R(pp).resp(ch);  if isempty(pp), isT = any(cellfun(@(v) v(ch), {R.tested})); else, isT = R(pp).tested(ch); end
                title(ax, sprintf('ch %d  (%d um)%s', ch, y(ch), ternary_(isT, '', ' (not tested)')), 'FontSize', Tile_FontSize + 1, ...
                      'Color', ternary_(isR, [0 0.55 0], [0.3 0.3 0.3]), 'FontWeight', ternary_(isR, 'bold', 'normal'));
            end
            sgtitle(fh, sprintf(['{\\bf\\fontsize{14}%s  |  cond %d: %s}\n%s  (%d trials)  page %d/%d (%d-%d um)   |  green title = recommended for %s; dotted = baseline mean, ' ...
                    '%s, yellow = response window, red = pulses; %s'], pTag, c, strrep(TT.cond(c).label, '_', '\_'), strrep(short, '_', '\_'), ...
                    numel(ks), pg, nPage, max(y(chs)), min(y(chs)), pName, ...
                    ternary_(isRate, 'green dashed = mean rate in the window', 'dashed = threshold'), ...
                    ternary_(isV, 'aligned to the photodiode', 'aligned to the first pulse')), 'FontSize', 9);
            set(fh, 'Visible', 'on');                                % draw the finished page once
        end
    end
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function B = byTrial_(tAbs, k, trigMs, nTr)
B = cell(nTr, 1);  for t = 1:nTr, B{t} = zeros(0, 1); end
if isempty(k), return; end
k = k(:);  rel = tAbs(:) - trigMs(k);
[ks, o] = sort(k);  rel = rel(o);
e = [0; find(diff(ks)); numel(ks)];
for i = 1:numel(e) - 1, rows = e(i)+1:e(i+1);  B{ks(rows(1))} = sort(rel(rows)); end
end

function rt = psthAligned_(B, blankRow, ks, off, ed, minFrac)
% blank-aware PSTH (spikes/s); times shifted by off(q) for trial ks(q); NaN where < minFrac of the trials had data
binMs = ed(2) - ed(1);  nB = numel(ed) - 1;  cnt = zeros(1, nB);  expo = zeros(1, nB);
for q = 1:numel(ks)
    t = ks(q);  cnt = cnt + histcounts(B{t} - off(q), ed);
    covr = ones(1, nB);  iv = blankRow{t};
    for b = 1:size(iv, 1)
        lo = max(ed(1:end-1), iv(b,1) - off(q));  hi = min(ed(2:end), iv(b,2) - off(q));
        covr = covr - max(0, hi - lo)/binMs;
    end
    expo = expo + max(covr, 0);
end
rt = cnt./(expo*binMs/1e3);
rt(expo < minFrac*max(numel(ks), 1)) = NaN;
end

function [mu, sd, thr, pk, lat, ok] = ruleCh_(B, blankRow, ks, off, bw, w, binMs, kSD, sdFloor, minC, minFrac)
% baseline bins: relative to the TRIGGER; response bins: relative to the alignment point (off)
mu = NaN; sd = NaN; thr = NaN; pk = NaN; lat = NaN; ok = false;
if isempty(ks), return; end
bEd = bw(1):binMs:bw(2);  rb = psthAligned_(B, blankRow, ks, zeros(size(ks)), bEd, minFrac);
mu = mean(rb, 'omitnan');  sd = std(rb, 'omitnan');  thr = mu + kSD*max(sd, sdFloor);
rEd = w(1):binMs:w(2);  rr = psthAligned_(B, blankRow, ks, off, rEd, minFrac);
pk = max(rr, [], 'omitnan');
above = rr > thr;  above(isnan(rr)) = false;
run = 0;
for i = 1:numel(above)
    if above(i), run = run + 1; else, run = 0; end
    if run >= minC, ok = true; lat = rEd(i - minC + 1); return; end
end
end

function W = rateWin_(tr, isVis, ks, tail, visW)
% response window of every trial (ms after the trigger): electrical = first pulse ... last pulse + tail; visual = visW after the photodiode
W = nan(numel(ks), 2);
for q = 1:numel(ks)
    t = ks(q);
    if isVis
        on = tr.diodeOnMs(t);  if isnan(on), on = 0; end
        W(q, :) = on + visW;
    else
        W(q, :) = [tr.firstPulseMs(t), tr.lastPulseMs(t) + tail];
    end
end
end

function [mB, mR, z, nPos, nT] = rateCh_(B, blankRow, ks, W, bw, mode)
% per trial: response rate in W(q,:) (blanked time not counted) and baseline rate in bw; one z per channel
nT = numel(ks);  rB = nan(nT, 1);  rR = nan(nT, 1);
for q = 1:nT
    t = ks(q);  s = B{t};
    rB(q) = sum(s >= bw(1) & s < bw(2))/(diff(bw)/1e3);
    len = diff(W(q, :));  iv = blankRow{t};
    for b = 1:size(iv, 1), len = len - max(0, min(W(q,2), iv(b,2)) - max(W(q,1), iv(b,1))); end
    if len > 0, rR(q) = sum(s >= W(q,1) & s < W(q,2) & ~inBlank_(s, iv))/(len/1e3); end
end
ok = ~isnan(rR);  rB = rB(ok);  rR = rR(ok);  nT = sum(ok);
mB = mean(rB);  mR = mean(rR);  d = rR - rB;  nPos = sum(d > 0);
if nT < 2, z = NaN; return; end
if strcmp(mode, 'paired'), sc = std(d)/sqrt(nT); else, sc = std(rB); end
if sc > 0, z = (mR - mB)/sc; elseif mR > mB, z = Inf; elseif mR < mB, z = -Inf; else, z = 0; end
end

function m = inBlank_(s, iv)
m = false(size(s));
for b = 1:size(iv, 1), m = m | (s >= iv(b,1) & s <= iv(b,2)); end
end

function z = trainZ_(B, ks, tr, bw, tail, isVis, visW)
% information: whole-train (or visual window) count vs baseline over trials, mean difference / SE
z = NaN;  if numel(ks) < 2, return; end
d = zeros(numel(ks), 1);
for q = 1:numel(ks)
    t = ks(q);  s = B{t};
    if isVis
        on = tr.diodeOnMs(t);  if isnan(on), on = 0; end
        w = on + visW;
    else
        w = [tr.firstPulseMs(t), tr.lastPulseMs(t) + tail];
    end
    nb = sum(s >= bw(1) & s < bw(2));
    d(q) = sum(s >= w(1) & s < w(2)) - nb*diff(w)/diff(bw);
end
se = std(d)/sqrt(numel(d));
if se > 0, z = mean(d)/se; end
end

function y = nanConv_(x, k)
m = ~isnan(x);  x0 = x;  x0(~m) = 0;
y = conv(x0, k, 'same')./conv(double(m), k, 'same');  y(~m) = NaN;
end

function v = nanTo0_(v), v(isnan(v)) = 0; end
function v = firstOr0_(p), if isempty(p), v = 0; else, v = p(1); end, end
function ok = smoothSigmaOK_(s, b), ok = s > 0 && s >= b/2; end

function s = chList_(v)
if isempty(v), s = '[]'; return; end
v = sort(v(:)');  parts = {};  i = 1;
while i <= numel(v)
    j = i;  while j < numel(v) && v(j+1) == v(j) + 1, j = j + 1; end
    if j == i, parts{end+1} = sprintf('%d', v(i)); else, parts{end+1} = sprintf('%d:%d', v(i), v(j)); end %#ok<AGROW>
    i = j + 1;
end
s = ['[' strjoin(parts, ' ') ']'];
end

function s = num2s_(v)
if isnan(v), s = '-'; elseif isinf(v), s = ternary_(v > 0, 'top', 'tip'); elseif v == round(v), s = sprintf('%d', v); else, s = sprintf('%.1f', v); end
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
