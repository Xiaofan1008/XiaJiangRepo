%% ============================================================
% NP_Train_FR_RespChnCheck_v1.m
% NP probe | pulse TRAIN sessions | FIRING-RATE analysis | RESPONSIVE-CHANNEL CHECK | v1   (NP version of NN_Train_FR_RespChnCheck_v1)
% READ-ONLY: prints and plots, NOTHING IS SAVED. The final channels are typed into NP_Train_FR_RespChnEntry_v1.
%
% RULE (as NN: baseline mean + k x SD, applied to the PSTH of each channel)
%   PSTH bins of binMs, blank-aware (a bin counts only the trials that had data in it; bins blanked in more than
%   (1 - minDataFrac) of the trials are ignored). Baseline = the bins in baseWinMs (relative to the TRIGGER): mean mu,
%   SD sd of the bin rates. Threshold = mu + kSD x max(sd, sdFloor_Hz).
%   A channel RESPONDS if at least minConsec CONSECUTIVE bins in the response window are above the threshold.
%   Response window: electrical = respWinMs after the FIRST pulse; visual = visRespWinMs after the photodiode onset.
%   Latency = start of the first such run; peak = highest bin in the window.
%   Visual conditions use their own (wider) bin, visBinMs, because their window is long (many bins = many chances).
%
% FALSE-POSITIVE CHECK (fpCheck = true)
%   The same rule, same bins and same threshold, applied to a STIMULUS-FREE window of the same length before the trigger
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
% CHANNELS   all channels processed in step 2 (goodCh); onlyActive = true -> only channels with a rate >= autoRateMin
% BAD TRIALS <short>.NP_BadTrials.mat is used when it exists (useBadTrials = true).
%
% FIGURES
%   1. ALONG THE PROBE: per decision condition, (peak - threshold) of every TESTED channel vs depth (green = responsive);
%      visual columns only when plotVisual = true;
%      last column: final-level train z vs depth. Depth groups (groupLines) can be drawn.
%   2. TILES (plotMode = 'psth' or 'rasterPsth'), for plotConds (default: the FINAL condition of each pattern, plus the
%      visual conditions when plotVisual = true). tileChannels = 'all': EVERY recorded (processed) channel, responsive or
%      not, tested or not ('tested' = only the tested channels). Untested channels have a grey title marked '(not tested)'.
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
onlyActive   = true;        % true = test only active channels (rate >= autoRateMin); false = every processed channel
autoRateMin  = 20;          % spikes/s over the whole window

% ---- rule ----
decideBy     = 'final';     % 'first' (sim n = 1, seq [1 1]; as NN) | 'final'
checkVisual  = true;        % also test every visual condition
baseWinMs    = [-100 -10];   % baseline, ms relative to the trigger
respWinMs    = [1 20];      % electrical: response window, ms after the FIRST pulse
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

% ---- figures ----
Plot_ProbeFigure = true;
plotVisual   = false;       % true = the visual conditions also get columns in the probe figure and tile pages (many figures)
tileChannels = 'all';       % 'all' = every recorded (processed) channel | 'tested' = only the tested channels | a list, e.g. [2 68 70]
groupLines   = [2200 3000 4500 5800];   % depths (um) drawn as lines in the probe figure; [] = none
plotMode     = 'rasterPsth';      % 'none' | 'psth' | 'rasterPsth'   (tiles)
plotConds    = [15 45];          % conditions to tile; [] = the FINAL condition of each pattern (add visual ones by number)
tilesPerPage = 48;          % channels per page (8 x 6)
nTileCols    = 8;
xRange       = [-50 200];    % ms around the alignment point (electrical: trigger; visual: photodiode onset)
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
rateAll = sum(SP.nSpk(:, 1:nTr), 2)/nTr/((SP.preMs + SP.postMs)/1e3);
chTest = find(SP.goodCh(:) & (~onlyActive | rateAll >= autoRateMin))';

%% ---- patterns: decision and final conditions ------------------------------------------------------------------------
Pt = struct('name', {}, 'decCond', {}, 'finalCond', {}, 'isVis', {});
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
    Pt(end+1) = struct('name', nm, 'decCond', dc, 'finalCond', P.finalCond, 'isVis', false); %#ok<SAGROW>
end
if checkVisual
    for c = find(strcmp(kindC, 'visual'))
        if any(trialCond == c), Pt(end+1) = struct('name', sprintf('visual %d', c), 'decCond', c, 'finalCond', c, 'isVis', true); end %#ok<SAGROW>
    end
end
nP = numel(Pt);

%% ---- spikes per trial (relative to the trigger) for the tested channels ------------------------------------------------
Bt = cell(nCh, 1);
for ch = chTest, Bt{ch} = byTrial_(SP.sp_clipped{ch}, SP.sp_trial{ch}, trigMs, nTr); end
alignOf = @(c, ks) ternary_(strcmp(kindC{c}, 'visual'), nanTo0_(tr.diodeOnMs(ks)), cellfun(@(v) firstOr0_(v), tr.pulseMs(ks)));

%% ---- rule for every pattern and channel -------------------------------------------------------------------------------------
R = struct('resp', {}, 'mu', {}, 'sd', {}, 'thr', {}, 'peak', {}, 'lat', {}, 'zFinal', {}, 'fp', {}, 'bin', {}, 'fpWin', {});
for p = 1:nP
    c = Pt(p).decCond;  ks = find(trialCond == c & useTr);  off = alignOf(c, ks);
    w = ternary_(Pt(p).isVis, visRespWinMs, respWinMs);
    bm = ternary_(Pt(p).isVis, visBinMs, binMs);
    fpW = ternary_(Pt(p).isVis, fpStartVisMs, fpStartMs) + [0 diff(w)];
    res = struct('resp', false(nCh, 1), 'mu', nan(nCh, 1), 'sd', nan(nCh, 1), 'thr', nan(nCh, 1), 'peak', nan(nCh, 1), ...
                 'lat', nan(nCh, 1), 'zFinal', nan(nCh, 1), 'fp', false(nCh, 1), 'bin', bm, 'fpWin', fpW);
    kf = find(trialCond == Pt(p).finalCond & useTr);
    for ch = chTest
        [mu, sd, thr, pk, lat, ok] = ruleCh_(Bt{ch}, SP.blank_ms(ch, :), ks, off, baseWinMs, w, bm, kSD, sdFloor_Hz, minConsec, minDataFrac);
        if fpCheck                                                  % same rule on a stimulus-free window (relative to the trigger)
            [~, ~, ~, ~, ~, res.fp(ch)] = ruleCh_(Bt{ch}, SP.blank_ms(ch, :), ks, zeros(size(off)), baseWinMs, fpW, bm, kSD, sdFloor_Hz, minConsec, minDataFrac);
        end
        res.resp(ch) = ok;  res.mu(ch) = mu;  res.sd(ch) = sd;  res.thr(ch) = thr;  res.peak(ch) = pk;  res.lat(ch) = lat;
        res.zFinal(ch) = trainZ_(Bt{ch}, kf, tr, baseWinMs, tail_ms, Pt(p).isVis, visRespWinMs);
    end
    R(p) = res;
end

%% ---- printout ------------------------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP TRAIN FR RESPONSIVE-CHANNEL CHECK (v1): %s\n', short);
fprintf('============================================================\n');
fprintf('spike file : %s\nbad trials : %s\nchannels   : %d tested (%s)\n', spikeFile, badTxt, numel(chTest), ...
        ternary_(onlyActive, sprintf('active, >= %g spikes/s', autoRateMin), 'all processed'));
fprintf('rule       : %g ms bins (visual %g ms), baseline %s ms, threshold = mu + %g x max(SD, %g Hz), >= %d consecutive bins in %s ms after the first pulse (visual: %s ms after the photodiode)\n', ...
        binMs, visBinMs, mat2str(baseWinMs), kSD, sdFloor_Hz, minConsec, mat2str(respWinMs), mat2str(visRespWinMs));
if fpCheck
    fprintf('false-pos. : the same rule on a stimulus-free window (electrical %s ms, visual %s ms, relative to the trigger)\n', ...
            mat2str(fpStartMs + [0 diff(respWinMs)]), mat2str(fpStartVisMs + [0 diff(visRespWinMs)]));
    if fpStartVisMs + diff(visRespWinMs) > baseWinMs(1) || fpStartMs + diff(respWinMs) > baseWinMs(1)
        warning('NPRC:fp', 'The false-positive window overlaps the baseline window: move fpStartMs / fpStartVisMs earlier.');
    end
end
fprintf('decides    : %s level\n', decideBy);
bandsE = [-Inf, groupLines(:)', Inf];
for p = 1:nP
    r = R(p);  rc = find(r.resp)';
    fprintf('\n--- %s: decision cond %d (%s, %d trials) | final cond %d ---\n', Pt(p).name, Pt(p).decCond, TT.cond(Pt(p).decCond).label, ...
            sum(trialCond == Pt(p).decCond & useTr), Pt(p).finalCond);
    fprintf('  responsive: %d of %d tested channels', numel(rc), numel(chTest));
    if fpCheck
        fprintf('   |  FALSE-POSITIVE CHECK: %d of %d (%.1f %%) "respond" in the stimulus-free window', sum(r.fp(chTest)), numel(chTest), ...
                100*mean(r.fp(chTest)));
    end
    fprintf('\n');
    for b = 1:numel(bandsE) - 1
        inB = chTest(y(chTest) >= bandsE(b) & y(chTest) < bandsE(b + 1));
        fprintf('    depth %5s-%5s um: %3d of %3d\n', num2s_(bandsE(b)), num2s_(bandsE(b + 1)), sum(r.resp(inB)), numel(inB));
    end
    if printAllChannels, chP = chTest; else, chP = rc; end
    if ~isempty(chP)
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

%% ---- figure 1: along the probe -----------------------------------------------------------------------------------------------
if Plot_ProbeFigure
    pShow = find(~[Pt.isVis] | plotVisual);  nS = numel(pShow);
    fh = figure('Name', [short ' responsive channels along the probe'], 'Color', 'w', 'Position', [40 40 min(1800, 260*(nS + 1)) 750]);
    for ip = 1:nS
        p = pShow(ip);  ax = subplot(1, nS + 1, ip);  hold(ax, 'on');  r = R(p);
        d = r.peak - r.thr;  nr = chTest(~r.resp(chTest));  rr = chTest(r.resp(chTest));
        plot(ax, d(nr), y(nr), '.', 'Color', [0.6 0.6 0.6]);  plot(ax, d(rr), y(rr), 'o', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 3);
        xline(ax, 0, 'k:');  for g = groupLines(:)', yline(ax, g, 'k-', 'Color', [0.7 0.7 0.7]); end
        grid(ax, 'on');  ylim(ax, [min(y) max(y)]);  xlabel(ax, 'peak - threshold (Hz)');
        title(ax, sprintf('%s\ncond %d: %d responsive', Pt(p).name, Pt(p).decCond, sum(r.resp)), 'FontSize', 8, 'Interpreter', 'none');
        if ip == 1, ylabel(ax, 'depth from the tip (um)'); else, set(ax, 'YTickLabel', []); end
    end
    ax = subplot(1, nS + 1, nS + 1);  hold(ax, 'on');  cl = lines(nS);
    for ip = 1:nS, plot(ax, R(pShow(ip)).zFinal(chTest), y(chTest), '.', 'Color', cl(ip,:)); end
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
        chTile = chTest;
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
        xr = ternary_(isV, xRangeVis, xRange);  w = ternary_(isV, visRespWinMs, respWinMs);
        ed = xr(1):binMs:xr(2);  ct = ed(1:end-1) + binMs/2;
        bmT = ternary_(isV, visBinMs, binMs);                       % bin of the RULE (threshold line)
        pl = [];  if ~isempty(ks) && ~isempty(tr.pulseMs{ks(1)}), pl = tr.pulseMs{ks(1)} - off(1); end
        for pg = 1:nPage
            chs = chP((pg - 1)*tilesPerPage + 1 : min(numel(chP), pg*tilesPerPage));
            fh = figure('Name', sprintf('%s cond %d page %d/%d', short, c, pg, nPage), 'Color', 'w', 'Position', [20 20 1700 980], 'Visible', 'off');
            for i = 1:numel(chs)
                ch = chs(i);  rI = ceil(i/nTileCols);  cI = i - (rI - 1)*nTileCols;
                ax = axes(fh, 'Position', [gL + (cI - 1)*tw + 0.012, gB + (nRows - rI)*th + 0.018, tw - 0.02, th - 0.045]);  hold(ax, 'on');
                [mu, ~, thr] = ruleCh_(Bt{ch}, SP.blank_ms(ch, :), ks, off, baseWinMs, w, bmT, kSD, sdFloor_Hz, minConsec, minDataFrac);
                rt = psthAligned_(Bt{ch}, SP.blank_ms(ch, :), ks, off, ed, minDataFrac);
                rs = rt;  if numel(ker) > 1, rs = nanConv_(rt, ker); end
                yTop = max([rs(:); thr; 1], [], 'omitnan')*1.15;
                patch(ax, [w(1) w(2) w(2) w(1)], [0 0 yTop yTop], [1 1 0.6], 'EdgeColor', 'none');
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
                plot(ax, xr, [mu mu], 'b:');  plot(ax, xr, [thr thr], 'b--');
                xlim(ax, xr);  ylim(ax, [0 yTop]);  set(ax, 'FontSize', Tile_FontSize, 'TickLength', [0.02 0.02]);
                isR = ~isempty(pp) && R(pp).resp(ch);  isT = ismember(ch, chTest);
                title(ax, sprintf('ch %d  (%d um)%s', ch, y(ch), ternary_(isT, '', ' (not tested)')), 'FontSize', Tile_FontSize + 1, ...
                      'Color', ternary_(isR, [0 0.55 0], [0.3 0.3 0.3]), 'FontWeight', ternary_(isR, 'bold', 'normal'));
            end
            sgtitle(fh, sprintf(['%s  cond %d: %s  (%d trials)  page %d/%d (%d-%d um)   |  green title = recommended for %s; dotted = baseline mean, ' ...
                    'dashed = threshold, yellow = response window, red = pulses; %s'], strrep(short, '_', '\_'), c, ...
                    strrep(TT.cond(c).label, '_', '\_'), numel(ks), pg, nPage, max(y(chs)), min(y(chs)), pName, ...
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
