%% ============================================================
% NN_Train_FR_RespChnCheck_v1.m
% NN probe | pulse TRAIN data | FIRING-RATE analysis | RESPONSIVE-CHANNEL CHECK | v1
%
% PURPOSE
%   Show every channel (raster + PSTH, 32 channels in ONE figure per condition) and print the response
%   numbers of every channel, plus a RECOMMENDED list of responsive channels, so that you can check it.
%   NOTHING IS SAVED and no existing file is changed. The final channels are typed into the entry
%   script that comes next.
%
% RULE (baseline mean + k * SD, applied to the PSTH)
%   PSTH bins of binMs. Baseline = the PSTH bins inside baseWin_ms (before the trigger): mean mu and SD sd
%   of the bin rates (spikes/s). Threshold = mu + kSD * max(sd, sdFloor_Hz).
%   A channel RESPONDS if at least minConsec CONSECUTIVE bins inside respWin_ms (after the trigger) are above
%   the threshold. Latency = start of the first such run; peak = highest bin in the window.
%   sdFloor_Hz stops a tiny SD (almost no baseline spikes) from making single spikes "significant" (placeholder).
%
% WHICH CONDITIONS DECIDE THE RECOMMENDATION
%   The first level with ONE pulse per electrode of each pattern: sim n=1, and for seq the first pair level [1 1]
%   (the very first seq level has only one electrode, so it is skipped). Real data, almost no copied parts. Choosing channels by the first pulse keeps the selection independent of what is analysed later
%   (late pulses / saturation). Recommended set = one set per pattern.
%   (This dataset has sim and seq B->A only.) Higher levels are only printed as information (they do NOT change the recommendation):
%     "+" = the last pulse of that level evokes a response (same rule, window respWin_ms after the last pulse)
%     "." = it does not.   Last column = mean baseline-subtracted whole-train count of the FINAL level.
%
% DATA
%   rec_tag = ''           : *.sp_xia_TrainRecovery.mat  (unthinned, v1 full recovery)  <- default
%   rec_tag = 'RandomEvent': *.sp_xia_TrainRecovery_RandomEvent.mat. NOTE: the decision condition of seq is level 2
%                            ([1 1]), which that file may have thinned, so make the DECISION on the unthinned file;
%                            use the RandomEvent file only to look at the information columns.
%   Trials with isRecovered = false in levels >= 2 are left out. excludeRel lists trials to leave out.
%
% FIGURES (Plot_Figures): default = the decision condition of each pattern; per condition one figure, 8 x 4 tiles (one per channel):
%   PSTH curve (black) with the raster overlaid behind it (grey ticks), dotted = mu, dashed = threshold,
%   yellow = response window, red lines = pulse times.
%   Green title = recommended.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
rec_tag     = '';          % '' = unthinned v1 file; 'RandomEvent' etc. = *.sp_xia_TrainRecovery_<rec_tag>.mat
rec_file    = '';          % full path (overrides rec_tag)
FS          = 30000;       % used only if the recovery file has no FS

patternConds = [];         % FINAL condition number(s) of the patterns; [] = all sim / AB / BA patterns found
chRange      = 1:32;       % channels used (also the channels in the figures)

% ---- rule ----
baseWin_ms  = [-50 -5];    % baseline window relative to the trigger (ms)
respWin_ms  = [1 20];      % response window after the pulse (ms)
binMs       = 1;           % PSTH bin width (ms)
kSD         = 3;           % threshold = mu + kSD * SD
minConsec   = 2;           % consecutive bins above the threshold that make a response
sdFloor_Hz  = 10;          % floor on the SD of the baseline bins (spikes/s)  (placeholder)
tail_ms     = 15;          % whole-train count window = first pulse ... last pulse + tail_ms (information column)

% ---- trials to leave out (manual; the bad-trial entry script will replace this later) ----
excludeRel = [];           % rows [cond  relTrialNumber], e.g. [15 3; 15 7]; [] = none

% ---- figures ----
Plot_Figures = true;
plotConds    = [];         % [] = the decision condition of every pattern; or a list of condition numbers
extraConds   = [];         % conditions added to the figures (e.g. a higher level of a pattern)
xRange       = [-20 80];   % ms around the trigger
nTileCols    = 8;          % tiles per row (32 channels -> 8 x 4)
yMaxPSTH     = [];         % fixed top of the PSTH axes (spikes/s); [] = automatic per tile (the y scale differs per channel)
curveSmooth_ms = 3;        % display only: Gaussian smoothing of the drawn PSTH curve (ms, SD = the smoothing factor); 0 = raw bins joined by a line. The rule always uses the raw bins.
Tile_FontSize = 7;

%% ================= CHECK SETTINGS =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(numel(baseWin_ms) == 2 && baseWin_ms(2) > baseWin_ms(1), 'baseWin_ms must be [start end].');
assert(numel(respWin_ms) == 2 && respWin_ms(2) > respWin_ms(1), 'respWin_ms must be [start end].');
assert(binMs > 0 && minConsec >= 1, 'binMs > 0 and minConsec >= 1.');
nBaseBins = floor(diff(baseWin_ms)/binMs + 1e-9);
nRespBins = floor(diff(respWin_ms)/binMs + 1e-9);
assert(nBaseBins >= 5, 'The baseline window holds only %d bins: widen baseWin_ms or lower binMs.', nBaseBins);
assert(nRespBins >= minConsec, 'The response window holds only %d bins: widen respWin_ms or lower binMs / minConsec.', nRespBins);

%% ================= LOAD =================
if isempty(rec_file)
    if isempty(rec_tag), pat_ = '*.sp_xia_TrainRecovery.mat'; else, pat_ = ['*.sp_xia_TrainRecovery_' rec_tag '.mat']; end
    d = dir(fullfile(data_folder, pat_));  d = d(~startsWith({d.name}, '.'));
    assert(~isempty(d), 'No %s in %s.', pat_, data_folder);
    assert(numel(d) == 1, 'More than one %s in %s - set rec_file.', pat_, data_folder);
    rec_file = fullfile(data_folder, d(1).name);
end
assert(isfile(rec_file), 'Cannot find %s', rec_file);
vn_ = who('-file', rec_file);
R = load(rec_file, vn_{ismember(vn_, {'sp_clipped','Recovery_params','FS','isRecovered'})});
assert(isfield(R, 'sp_clipped') && isfield(R, 'Recovery_params'), '%s is not a TrainRecovery file.', rec_file);
sp = R.sp_clipped;
if isfield(R, 'FS') && ~isempty(R.FS), FS = double(R.FS); end
trialCond = R.Recovery_params.trialCond(:);
condTable = R.Recovery_params.condTable;
nChn = numel(sp);  nTrialsAll = numel(trialCond);  nCond = numel(condTable);
if isfield(R, 'isRecovered') && numel(R.isRecovered) == nTrialsAll
    isRec = logical(R.isRecovered(:));  haveRec = true;
else
    isRec = true(nTrialsAll, 1);  haveRec = false;
end
clear R

ft = dir(fullfile(data_folder, '*.trig.dat'));  ft = ft(~startsWith({ft.name}, '.'));
assert(~isempty(ft), 'No *.trig.dat in %s.', data_folder);
fid = fopen(fullfile(data_folder, ft(1).name), 'r');  trig = fread(fid, [1 inf], 'double');  fclose(fid);
t0 = nan(nTrialsAll, 1);  nUse = min(numel(trig), nTrialsAll);  t0(1:nUse) = trig(1:nUse)/FS*1000;
if numel(trig) < nTrialsAll
    warning('RC:trig', '%d triggers for %d trials: the last trials have no trigger and are ignored.', numel(trig), nTrialsAll);
end

fe = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));  fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fullfile(data_folder, fe(1).name), 'StimParams', 'simultaneous_stim');
nEl = double(E.simultaneous_stim);
peK = zeros(nCond, nEl);  pulseT = cell(nCond, 1);
for k = 1:nCond
    t1 = find(trialCond == k, 1);
    if isempty(t1), continue; end
    for e = 1:nEl
        row = 1 + (t1-1)*nEl + e;
        if double(E.StimParams{row,16}) > 0
            n_ = double(E.StimParams{row,8});  per = double(E.StimParams{row,9});  p0 = double(E.StimParams{row,6});
            peK(k,e) = n_;
            pulseT{k} = [pulseT{k}, (p0 + (0:n_-1)*per)/1000];
        end
    end
    pulseT{k} = unique(pulseT{k});
end
clear E

chUse = chRange(:)';
assert(~isempty(chUse) && all(chUse >= 1 & chUse <= nChn), 'chRange must lie between 1 and %d.', nChn);
nChUse = numel(chUse);

%% ================= PATTERNS AND CHAINS =================
kinds = arrayfun(@(c) char(c.kind), condTable, 'UniformOutput', false);
isFinal = arrayfun(@(c) isfield(c,'isFinal') && ~isempty(c.isFinal) && logical(c.isFinal), condTable);
if isempty(patternConds)
    patternConds = find(isFinal(:)' & ismember(kinds(:)', {'sim','AB','BA'}));
end
assert(~isempty(patternConds), 'No sim / AB / BA pattern found in the recovery file.');
nPat = numel(patternConds);
chain = cell(nPat,1);  patName = cell(nPat,1);
for p = 1:nPat
    c = patternConds(p);  ch_ = c;
    while isfinite(condTable(ch_(1)).sourceCond), ch_ = [condTable(ch_(1)).sourceCond, ch_]; end %#ok<AGROW>
    chain{p} = ch_;
    switch kinds{c}
        case 'sim', patName{p} = 'sim';
        case 'AB',  patName{p} = 'seq A->B';
        case 'BA',  patName{p} = 'seq B->A';
        otherwise,  patName{p} = kinds{c};
    end
end

%% ================= TRIAL SETS =================
validOf = @(cond) validTrials_(cond, trialCond, t0, isRec, condTable, excludeRel);

%% ================= DECISION CONDITION AND INFORMATION PER LEVEL =================
decision = cell(nPat, 1);          % struct array over channels (decision condition)
decCond  = zeros(nPat, 1);         % decision condition of each pattern (first level with one pulse per electrode)
levelN   = cell(nPat, 1);          % pulses per electrode of each pair level
lastResp = cell(nPat, 1);          % [nChUse x nLevels] logical: last pulse evokes a response
wholeCnt = cell(nPat, 1);          % [nChUse x nLevels] mean baseline-subtracted whole-train count
nTrLvl   = cell(nPat, 1);
for p = 1:nPat
    cc = chain{p};  L = numel(cc);
    % pair levels only ([n n] for seq)
    isPair = false(1, L);
    for n = 1:L, isPair(n) = all(peK(cc(n),:) == peK(cc(n),1)); end
    lv = find(isPair);
    levelN{p} = arrayfun(@(n) max(peK(cc(n),:)), lv);
    lastResp{p} = false(nChUse, numel(lv));  wholeCnt{p} = nan(nChUse, numel(lv));  nTrLvl{p} = zeros(1, numel(lv));
    for q = 1:numel(lv)
        n = lv(q);  cond = cc(n);
        if q == 1, decCond(p) = cond; end
        tr = validOf(cond);  tt = t0(tr);  nTrLvl{p}(q) = numel(tr);
        if numel(tr) < 3, continue; end
        tL = double(condTable(cond).lastPulse_ms);
        b = tL + tail_ms;  scaleF = b / diff(baseWin_ms);
        for ci = 1:nChUse
            s = sp{chUse(ci)}(:,1);
            if q == 1
                decision{p}(ci) = analyzeCh_(s, tt, baseWin_ms, binMs, kSD, minConsec, sdFloor_Hz, respWin_ms); %#ok<SAGROW>
            end
            rl = analyzeCh_(s, tt, baseWin_ms, binMs, kSD, minConsec, sdFloor_Hz, tL + respWin_ms);
            lastResp{p}(ci, q) = rl.resp;
            cs = zeros(numel(tt), 1);
            for i = 1:numel(tt)
                cs(i) = (lb_(s, tt(i)+b) - lb_(s, tt(i))) - (lb_(s, tt(i)+baseWin_ms(2)) - lb_(s, tt(i)+baseWin_ms(1))) * scaleF;
            end
            wholeCnt{p}(ci, q) = mean(cs);
        end
    end
    assert(~isempty(decision{p}), 'Pattern %s has no usable trials in its decision condition.', patName{p});
end

%% ================= PRINT =================
fprintf('\n============================================================\n');
fprintf('NN TRAIN FR RESPONSIVE-CHANNEL CHECK (v1)\n');
fprintf('============================================================\n');
fprintf('Folder   : %s\nRecovery : %s\n', data_folder, rec_file);
fprintf('Channels : %s (%d)\n', fmtList_(chUse), nChUse);
fprintf('Rule     : PSTH bins %g ms | baseline [%g, %g] ms: threshold = mu + %g * max(SD, %g Hz) | response = >= %d consecutive bins above it in [%g, %g] ms\n', ...
    binMs, baseWin_ms, kSD, sdFloor_Hz, minConsec, respWin_ms);
fprintf('Decision : first level with ONE pulse per electrode of each pattern (sim n=1; seq [1 1], i.e. the first pair level). Higher levels are information only.\n');
if ~haveRec, fprintf('Note     : no isRecovered in the file; all trials used.\n'); end
if ~isempty(excludeRel), fprintf('Excluded : %d manual trial(s) (excludeRel).\n', size(excludeRel,1)); end

recSet = cell(nPat, 1);
for p = 1:nPat
    D = decision{p};  cond1 = decCond(p);  tr1 = validOf(cond1);
    recSet{p} = chUse([D.resp]);
    fprintf('\n--- %s | decision condition %d | %d trials ---\n', upper(patName{p}), cond1, numel(tr1));
    fprintf('  ch   base(Hz)  SD(Hz)  thr(Hz)  peak(Hz)  latency(ms)  added/trial   RESPONSIVE\n');
    for ci = 1:nChUse
        x = D(ci);
        fprintf(' %3d   %7.1f  %6.1f  %7.1f  %8.1f   %9s   %9.2f      %s\n', chUse(ci), x.mu, x.sdUsed, x.thr, x.peak, ...
            ternary_(isfinite(x.lat), sprintf('%.1f', x.lat), '-'), x.addCount, ternary_(x.resp, 'YES', '-'));
    end
    fprintf('  recommended (%d of %d): %s\n', numel(recSet{p}), nChUse, fmtList_(recSet{p}));
end

fprintf('\n--- LAST-PULSE RESPONSE OF HIGHER LEVELS (information only; + response, . none) ---\n');
for p = 1:nPat
    fprintf('\n%s  (columns = pulses per electrode: %s)\n', upper(patName{p}), strtrim(sprintf('%d ', levelN{p})));
    fprintf('  ch  rec  levels%s  count(final)\n', repmat(' ', 1, max(0, 2*numel(levelN{p}) - 6 + 1)));
    for ci = 1:nChUse
        c_ = ternary_cell_(lastResp{p}(ci,:));
        str = sprintf('%s ', c_{:});
        fprintf(' %3d   %s    %s  %8.2f\n', chUse(ci), ternary_(decision{p}(ci).resp, 'Y', '-'), str, wholeCnt{p}(ci, end));
    end
    fprintf('  (trials per level: %s)\n', strtrim(sprintf('%d ', nTrLvl{p})));
end

fprintf('\n--- SUMMARY ---\n');
allSet = chUse;  anySet = [];
for p = 1:nPat, allSet = intersect(allSet, recSet{p}); anySet = union(anySet, recSet{p}); end
fprintf(' recommended in ALL patterns : %s (%d)\n', fmtList_(allSet), numel(allSet));
fprintf(' recommended in ANY pattern  : %s (%d)\n', fmtList_(anySet), numel(anySet));
fprintf('\n--- PASTE-READY LIST for the responsive-channel entry script ---\n');
fprintf('%%   pattern      channels\n');
for p = 1:nPat
    fprintf('    ''%s'',%s %s;\n', patName{p}, repmat(' ', 1, max(1, 11-numel(patName{p}))), fmtList_(recSet{p}));
end
fprintf('\n(Nothing was saved. Check the figures, then edit the lists before using them.)\n');

%% ================= FIGURES =================
if Plot_Figures
    if isempty(plotConds), plotConds = decCond'; end
    plotConds = unique([plotConds(:)', extraConds(:)'], 'stable');
    assert(all(plotConds >= 1 & plotConds <= nCond), 'plotConds / extraConds must lie between 1 and %d.', nCond);
    nR = ceil(nChUse / nTileCols);
    for k = plotConds
        tr = validOf(k);
        if numel(tr) < 3, fprintf('cond %d has fewer than 3 usable trials: no figure.\n', k); continue; end
        tt = t0(tr);
        fig = figure('Color','w','Units','normalized','Position',[0.01 0.03 0.98 0.90], ...
            'Name', sprintf('cond %d %s', k, char(condTable(k).label)));
        mL = 0.035; mR = 0.01; mB = 0.05; mT = 0.06; gx = 0.014; gy = 0.03;
        tw = (1 - mL - mR - (nTileCols-1)*gx) / nTileCols;
        th = (1 - mB - mT - (nR-1)*gy) / nR;
        edgesP = xRange(1):binMs:xRange(2);  ctr = edgesP(1:end-1) + binMs/2;
        for ci = 1:nChUse
            cc_ = mod(ci-1, nTileCols);  rr = floor((ci-1)/nTileCols);
            x0 = mL + cc_*(tw+gx);  y0 = 1 - mT - (rr+1)*th - rr*gy;
            ax = axes('Position', [x0, y0, tw, 0.86*th], 'Parent', fig);  hold(ax, 'on');
            s = sp{chUse(ci)}(:,1);
            % PSTH (raw bins decide the rule; the drawn curve may be smoothed for display)
            cnt = cnt_(s, tt, edgesP);
            rate = cnt / numel(tt) / binMs * 1000;
            r = analyzeCh_(s, tt, baseWin_ms, binMs, kSD, minConsec, sdFloor_Hz, respWin_ms);
            rateD = rate;
            if curveSmooth_ms > 0
                sdb = curveSmooth_ms / binMs;  xg = -ceil(4*sdb):ceil(4*sdb);
                g = exp(-0.5*(xg/sdb).^2);  g = g / sum(g);
                rateD = conv(rate, g, 'same');
            end
            ymax = max([rateD, r.thr, 1]) * 1.1;  if ~isempty(yMaxPSTH), ymax = yMaxPSTH; end
            patch(ax, [respWin_ms(1) respWin_ms(2) respWin_ms(2) respWin_ms(1)], [0 0 ymax ymax], [1 0.85 0.2], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
            % raster, drawn behind the curve over the full height (trial i occupies the band i of nTrials)
            X = [];  Y = [];  nT_ = numel(tt);
            for i = 1:nT_
                a_ = lb_(s, tt(i)+xRange(1));  b_ = lb_(s, tt(i)+xRange(2));
                rel = s(a_:b_-1) - tt(i);
                y1 = ymax*(i-1+0.1)/nT_;  y2 = ymax*(i-0.1)/nT_;
                X = [X, [rel(:)'; rel(:)'; nan(1, numel(rel))]]; %#ok<AGROW>
                Y = [Y, [y1+zeros(1,numel(rel)); y2+zeros(1,numel(rel)); nan(1, numel(rel))]]; %#ok<AGROW>
            end
            if ~isempty(X), plot(ax, X(:), Y(:), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.5); end
            for tp = pulseT{k}(:)'
                if tp >= xRange(1) && tp <= xRange(2), plot(ax, [tp tp], [0 ymax], '-', 'Color', [0.8 0.5 0.5], 'LineWidth', 0.5); end
            end
            plot(ax, xRange, r.mu*[1 1], ':', 'Color', [0.1 0.4 0.9], 'LineWidth', 1);
            plot(ax, xRange, r.thr*[1 1], '--', 'Color', [0.85 0.1 0.1], 'LineWidth', 1);
            plot(ax, ctr, rateD, '-', 'Color', [0 0 0], 'LineWidth', 1.3);
            xlim(ax, xRange);  ylim(ax, [0 ymax]);
            set(ax, 'FontSize', Tile_FontSize, 'Box', 'off', 'TickLength', [0.02 0.02]);
            if rr ~= nR-1 && ci + nTileCols <= nChUse, set(ax, 'XTickLabel', []); end
            col = ternary_(r.resp, [0 0.55 0.1], [0.3 0.3 0.3]);
            title(ax, sprintf('ch %d%s', chUse(ci), ternary_(r.resp, '  RESP', '')), 'FontSize', Tile_FontSize+2, 'Color', col, 'FontWeight', 'bold');
        end
        annotation(fig, 'textbox', [0 0.965 1 0.035], 'String', sprintf(['cond %d  %s  |  %d trials  |  black = PSTH (%g ms bins%s), grey ticks = raster (trials stacked over the full height)  |  ' ...
            'dotted = baseline mean, dashed = mean + %g SD, yellow = response window [%g, %g] ms, red lines = pulses  |  green = responsive'], ...
            k, char(condTable(k).label), numel(tt), binMs, ternary_(curveSmooth_ms > 0, sprintf(', smoothed %g ms', curveSmooth_ms), ''), kSD, respWin_ms), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Interpreter', 'none');
    end
end

%% ================= HELPERS =================
function tr = validTrials_(cond, trialCond, t0, isRec, condTable, excludeRel)
all_ = find(trialCond == cond);
isL1 = ~isfinite(condTable(cond).sourceCond);
ok = isfinite(t0(all_)) & (isRec(all_) | isL1);
if ~isempty(excludeRel)
    rows = excludeRel(excludeRel(:,1) == cond, 2);
    ok(rows(rows >= 1 & rows <= numel(all_))) = false;
end
tr = all_(ok);
end

function r = analyzeCh_(s, tt, bw, binMs, kSD, minC, sdFloor, w)
% baseline mean / SD of the PSTH bins in bw; response = >= minC consecutive bins above mu + kSD*SD inside w
bE = edgesOf_(bw, binMs);  rE = edgesOf_(w, binMs);
n = numel(tt);
cb = cnt_(s, tt, bE);  rb = cnt_(s, tt, rE);
baseRate = cb / n / binMs * 1000;  respRate = rb / n / binMs * 1000;
mu = mean(baseRate);  sd = std(baseRate);  sdU = max(sd, sdFloor);
thr = mu + kSD*sdU;
above = respRate > thr;
run = 0;  best = 0;  startIdx = NaN;  cur = 0;
for i = 1:numel(above)
    if above(i), cur = cur + 1; else, cur = 0; end
    if cur > best, best = cur; end
    if cur == minC && isnan(startIdx), startIdx = i - minC + 1; end
end
r.mu = mu;  r.sd = sd;  r.sdUsed = sdU;  r.thr = thr;
r.peak = max(respRate);  r.run = best;  r.resp = best >= minC;
r.lat = NaN;  if isfinite(startIdx), r.lat = rE(startIdx); end
r.addCount = sum(rb)/n - mu*numel(rb)*binMs/1000;
end

function e = edgesOf_(w, binMs)
nb = floor(diff(w)/binMs + 1e-9);
e = w(1) + (0:nb)*binMs;
end

function c = cnt_(s, tt, edges)
c = zeros(1, numel(edges)-1);
for i = 1:numel(tt)
    a = lb_(s, tt(i)+edges(1));  b = lb_(s, tt(i)+edges(end));
    if b > a, c = c + histcounts(s(a:b-1) - tt(i), edges); end
end
end

function i = lb_(s, x)
lo = 1;  hi = numel(s) + 1;
while lo < hi
    m = floor((lo + hi)/2);
    if s(m) < x, lo = m + 1; else, hi = m; end
end
i = lo;
end

function s = fmtList_(v)
v = sort(v(:)');
if isempty(v), s = '[]'; return; end
parts = {};  i = 1;
while i <= numel(v)
    j = i;
    while j < numel(v) && v(j+1) == v(j) + 1, j = j + 1; end
    if j == i,        parts{end+1} = sprintf('%d', v(i));          %#ok<AGROW>
    elseif j == i+1,  parts{end+1} = sprintf('%d %d', v(i), v(j)); %#ok<AGROW>
    else,             parts{end+1} = sprintf('%d:%d', v(i), v(j)); %#ok<AGROW>
    end
    i = j + 1;
end
s = ['[' strjoin(parts, ' ') ']'];
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end

function c = ternary_cell_(tf)
c = repmat({'.'}, 1, numel(tf));
c(tf) = {'+'};
end
