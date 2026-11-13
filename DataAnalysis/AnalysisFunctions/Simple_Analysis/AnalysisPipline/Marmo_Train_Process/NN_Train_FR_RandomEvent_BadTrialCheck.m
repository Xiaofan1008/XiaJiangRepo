%% ============================================================
% NN_Train_FR_RandomEvent_BadTrialCheck.m
% NN probe | pulse TRAIN data | FIRING-RATE (spike count) analysis | BAD-TRIAL CHECK | RandomEvent version of NN_Train_FR_BadTrialCheck_v2  (v1, v2 untouched)
%
% This script = v2 for the file made by NN_SpikeRecovery_RandomEvent_v1 (*.sp_xia_TrainRecovery_RandomEvent.mat).
%   Differences to v2:
%   - the recovery file is found by rec_tag (default 'RandomEvent'), or set rec_file.
%   - the recovery settings of every level (percent, event_ms) are printed in the table.
%   - trials that were NOT recovered (isRecovered = false in a level >= 2, e.g. no usable source trial)
%     keep a blank gap: they are listed as 'unrec', are NOT used for the MAD statistics, and are not drawn.
%     Level 1 is never recovered by design and is always used.
%   - the RandomEvent thinning is random (and pooled over trials), so some trials lose more spikes than others
%     by chance. A flag on a thinned level may therefore come from the simulation and not from the trial.
%     The same rules are used as in v2; compare with the flags of the unthinned run (v2) before trusting them.
%
% PURPOSE
%   (v2 text) Look at every trial of every recovered level (1, 2, 3 ... N pulses) of the simultaneous
%   and the sequential pattern, flag suspicious trials, and show them so you can decide which
%   trials to remove. Each recovered level is treated as its OWN, independent condition.
%   NOTHING IS SAVED and no existing file is changed (final bad trials are typed into the
%   entry script that comes next).
%
% PER-TRIAL VALUES (over the channels in chRange)
%   Window : whole train, from the first pulse (+trainStart_ms) to the last pulse of the level
%            + tail_ms. The window therefore grows with the number of pulses.
%   Metric : mean baseline-subtracted spike count per channel         (measure 1 = RESPONSE)
%            count in the window - count in baseWin_ms * (window length / baseline length)
%   MaxCh  : largest single-channel baseline-subtracted count          (measure 2 = MAX single chn)
%   Base   : mean baseline firing rate (Hz) of the channels            (measure 3 = BASELINE)
%
% LEVELS SHOWN
%   simultaneous : every level n (n pulses on both electrodes)
%   sequential   : the [n n] levels (n pulses per electrode), comparable with simultaneous n.
%                  seqLevels = 'all' also adds the half-steps [n n-1] (one more pulse on the
%                  electrode that goes first), drawn as triangles.
%
% RULES (all relative to the other trials of the SAME level; thresholds are placeholders)
%   R+  RESP-MAD high : metric - median  > respMAD_k * max(1.4826*MAD, madFloor)
%   R-  RESP-MAD low  : median - metric  > respMAD_k * max(1.4826*MAD, madFloor)   (Resp_FlagLow)
%   M   MAXCH-MAD     : max-channel count - median > maxChMAD_k * max(1.4826*MAD, maxChFloor)  (high side only)
%   B   BASE-MAD      : baseline rate - median > baseMAD_k * max(1.4826*MAD, baseMadFloorHz)    (high side only)
%   MAD groups = one condition (one level of one pattern). Each rule can be switched off.
%   No absolute limits: you judge extreme values by eye from the figures.
%   The label of a flagged trial shows the codes, e.g. "11 B" or "4 R+M".
%
% TRIAL NUMBERS
%   "rel" = trial number WITHIN the condition (1 ... number of repetitions), as printed on the figures.
%   The absolute trial number in the recording is printed next to each flagged trial.
%
% OUTPUT (nothing is saved)
%   1. table: every level, trials, median metric, flagged trials and the rule(s) that flagged them
%   2. exclusion rate per pattern and level (warning above Warn_Exclusion_Fraction)
%   3. responsive channels per level (z of the baseline-subtracted count > zThr)
%   4. a paste-ready list (pattern, level, trials) for the bad-trial entry script
%   5. figure per measure : one panel per pulse number, trial number vs value, every point labelled
%   6. summary figure     : value of every trial vs number of pulses (scatter, all levels together)
%   7. channel figure     : z of every channel for every level
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
rec_tag     = 'RandomEvent';   % file = *.sp_xia_TrainRecovery_<rec_tag>.mat in data_folder  ('' = the untagged v1 file)
rec_file    = '';          % '' = find by rec_tag; otherwise the full path (overrides rec_tag)
FS          = 30000;       % used only if the recovery file has no FS

patternConds = [];         % FINAL condition number(s) of the patterns to check; [] = all sim / AB / BA patterns found
chRange      = 1:32;       % channels used for the per-trial values, e.g. 1:32, 5:12, [3 7 9]

% ---- counting window and baseline ----
tail_ms       = 15;        % window = first pulse + trainStart_ms  ...  last pulse + tail_ms
trainStart_ms = 0;         % start of the window relative to the first pulse (ms)
baseWin_ms    = [-50 -5];  % baseline window relative to the trigger (ms); CHANGE HERE if needed

% ---- which levels ----
seqLevels = 'pair';        % sequential patterns: 'pair' ([n n] only) | 'all' (also the half-steps [n n-1])

% ---- spike-count rules (MAD only; PLACEHOLDERS) ----
Use_RespMADRule  = true;   respMAD_k       = 4;     madFloor = 1.0;       % floor on the MAD scale (counts)
Resp_FlagLow     = true;   % also flag trials whose response is unusually LOW (R-)
Use_MaxChMADRule = true;   maxChMAD_k      = 4;     maxChFloor = 2.0;     % floor on the MAD scale (counts), high side only
Use_BaseMADRule  = true;   baseMAD_k       = 4;     baseMadFloorHz = 10;  % floor on the MAD scale (Hz), high side only
minTrialsForMAD  = 5;      % MAD rules need at least this many valid trials in the condition
Warn_Exclusion_Fraction = 0.30;   % warn when a level loses more than this share of its trials

% ---- responsive channels (guide for the responsive-channel entry) ----
zThr  = 3;                 % channel counts as responsive if mean(baseline-subtracted count)/SEM over trials > zThr
zClip = 30;                % z values are clipped to +-zClip in the figure

% ---- figures ----
Trial_Measure        = [1 2 3]; % panel figures to draw: 1 RESPONSE | 2 MAX single chn | 3 BASELINE
Plot_PanelFigures    = true;
Plot_SummaryFigure   = true;    % value of every trial vs number of pulses (one row per measure in Trial_Measure)
Plot_ChannelResponse = true;
nPanelCols           = 5;       % panels per row in the panel figures
Label_All_Trials     = true;    % trial number next to every point
Label_FontSize       = 11;      % trial number labels (black; flagged ones bold)
Series_Shift         = 0.30;    % sim and seq points are shifted left/right by this much (in trial-number units) so they do not overlap
Marker_Size_Kept     = 6;
Marker_Size_Flagged  = 10;     % flagged = series colour and shape + thick red ring
Show_Median_Lines    = true;

%% ================= CHECK SETTINGS =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(ismember(seqLevels, {'pair','all'}), 'seqLevels: pair or all.');
assert(numel(baseWin_ms) == 2 && baseWin_ms(2) > baseWin_ms(1), 'baseWin_ms must be [start end] with end > start.');
assert(tail_ms >= 0, 'tail_ms must be >= 0.');
Trial_Measure = unique(Trial_Measure(:)');
assert(~isempty(Trial_Measure) && all(ismember(Trial_Measure, 1:3)), 'Trial_Measure may only contain 1, 2, 3.');

%% ================= LOAD =================
if isempty(rec_file)
    if isempty(rec_tag), pat_ = '*.sp_xia_TrainRecovery.mat'; else, pat_ = ['*.sp_xia_TrainRecovery_' rec_tag '.mat']; end
    d = dir(fullfile(data_folder, pat_));
    d = d(~startsWith({d.name}, '.'));
    assert(~isempty(d), 'No %s in %s (run the recovery script first).', pat_, data_folder);
    assert(numel(d) == 1, 'More than one %s in %s - set rec_file.', pat_, data_folder);
    rec_file = fullfile(data_folder, d(1).name);
end
assert(isfile(rec_file), 'Cannot find %s', rec_file);
R = load(rec_file, 'sp_clipped', 'Recovery_params', 'FS', 'isRecovered');
assert(isfield(R, 'sp_clipped') && isfield(R, 'Recovery_params'), '%s is not a TrainRecovery file.', rec_file);
sp = R.sp_clipped;
if isfield(R, 'FS') && ~isempty(R.FS), FS = double(R.FS); end
trialCond = R.Recovery_params.trialCond(:);
condTable = R.Recovery_params.condTable;
if isfield(R, 'isRecovered') && numel(R.isRecovered) == numel(trialCond)
    isRec = logical(R.isRecovered(:));  haveRec = true;
else
    isRec = true(numel(trialCond), 1);  haveRec = false;
    warning('BT:rec', 'No isRecovered in this file: all trials are treated as recovered.');
end
rp_ = R.Recovery_params;
pctK = 100*ones(numel(condTable), 1);  evK = nan(numel(condTable), 1);
if isfield(rp_, 'percent_k'), pctK = double(rp_.percent_k(:)); end
if isfield(rp_, 'event_ms_k'), evK = double(rp_.event_ms_k(:)); end
nChn = numel(sp);  nTrialsAll = numel(trialCond);  nCond = numel(condTable);
clear R

ft = dir(fullfile(data_folder, '*.trig.dat'));  ft = ft(~startsWith({ft.name}, '.'));
assert(~isempty(ft), 'No *.trig.dat in %s.', data_folder);
fid = fopen(fullfile(data_folder, ft(1).name), 'r');  trig = fread(fid, [1 inf], 'double');  fclose(fid);
t0 = nan(nTrialsAll, 1);  nUse = min(numel(trig), nTrialsAll);  t0(1:nUse) = trig(1:nUse)/FS*1000;   % trigger time (ms)
if numel(trig) < nTrialsAll
    warning('BT:trig', '%d triggers for %d trials: the last %d trials have no trigger and are ignored.', numel(trig), nTrialsAll, nTrialsAll-numel(trig));
end

fe = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));  fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fullfile(data_folder, fe(1).name), 'StimParams', 'simultaneous_stim');
nEl = double(E.simultaneous_stim);
% pulses per electrode of every condition (rows of the exp file with amplitude > 0)
peK = zeros(nCond, nEl);
for k = 1:nCond
    t1 = find(trialCond == k, 1);
    if isempty(t1), continue; end
    for e = 1:nEl
        row = 1 + (t1-1)*nEl + e;
        if double(E.StimParams{row,16}) > 0, peK(k,e) = double(E.StimParams{row,8}); end
    end
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
chain = cell(nPat,1);  tn = cell(nPat,1);  patName = cell(nPat,1);
for p = 1:nPat
    c = patternConds(p);
    ch_ = c;
    while isfinite(condTable(ch_(1)).sourceCond)
        ch_ = [condTable(ch_(1)).sourceCond, ch_]; %#ok<AGROW>
    end
    chain{p} = ch_;
    tn{p} = arrayfun(@(k) condTable(k).lastPulse_ms, ch_);
    assert(all(diff(tn{p}) > 0), 'cond %d: pulse times of the levels are not increasing.', c);
    switch kinds{c}
        case 'sim', patName{p} = 'sim';
        case 'AB',  patName{p} = 'seq A->B';
        case 'BA',  patName{p} = 'seq B->A';
        otherwise,  patName{p} = kinds{c};
    end
end

%% ================= PER-LEVEL, PER-TRIAL VALUES =================
ruleNames = {'RESP+','RESP-','MAXCH','BASE'};     % rows of F (long names, for the table)
codeNames = {'R+','R-','M','B'};                  % short codes for the figure labels
baseLen = diff(baseWin_ms);
Lv = struct([]);
iL = 0;
for p = 1:nPat
    cc = chain{p};  L = numel(cc);  t = tn{p};  pe = peK(cc,:);
    for n = 1:L
        isPair = all(pe(n,:) == pe(n,1));
        if ~isPair && strcmp(seqLevels, 'pair'), continue; end
        nn   = max(pe(n,:));
        cond = cc(n);
        trAll = find(trialCond == cond);
        nT = numel(trAll);
        tt = t0(trAll);
        isL1 = ~isfinite(condTable(cond).sourceCond);        % level 1: never recovered by design
        recOK = isRec(trAll) | isL1;
        valid = isfinite(tt) & recOK(:);
        nUnrec = sum(isfinite(tt) & ~recOK(:));
        a = t(1) + trainStart_ms;  b = t(n) + tail_ms;
        scaleF = (b - a)/baseLen;

        CNT = nan(nT, nChUse);  BAS = nan(nT, nChUse);
        for ci = 1:nChUse
            s = sp{chUse(ci)}(:,1);
            for i = find(valid(:))'
                CNT(i,ci) = lb_(s, tt(i)+b) - lb_(s, tt(i)+a);
                BAS(i,ci) = lb_(s, tt(i)+baseWin_ms(2)) - lb_(s, tt(i)+baseWin_ms(1));
            end
        end
        CORR   = CNT - BAS*scaleF;
        metric = mean(CORR, 2);
        maxch  = max(CORR, [], 2);
        baseHz = mean(BAS, 2) / (baseLen/1000);
        vi = find(valid(:));

        % ---- rules (rows of F: R+, R-, M, B) ----
        F = false(4, nT);
        if numel(vi) >= minTrialsForMAD
            med = median(metric(vi));  sc = max(1.4826*median(abs(metric(vi)-med)), madFloor);
            F(1,:) = Use_RespMADRule & ((metric' - med) > respMAD_k*sc);
            F(2,:) = Use_RespMADRule & Resp_FlagLow & ((med - metric') > respMAD_k*sc);
            medM = median(maxch(vi));  scM = max(1.4826*median(abs(maxch(vi)-medM)), maxChFloor);
            F(3,:) = Use_MaxChMADRule & ((maxch' - medM) > maxChMAD_k*scM);
            medB = median(baseHz(vi));  scB = max(1.4826*median(abs(baseHz(vi)-medB)), baseMadFloorHz);
            F(4,:) = Use_BaseMADRule & ((baseHz' - medB) > baseMAD_k*scB);
        end

        % ---- responsive channels: z of the baseline-subtracted count over trials ----
        z = nan(1, nChUse);
        for ci = 1:nChUse
            dd = CORR(vi, ci);
            if numel(dd) >= 3
                sd = std(dd);  mu = mean(dd);
                if sd > 0, z(ci) = mu / (sd/sqrt(numel(dd)));
                elseif mu ~= 0, z(ci) = sign(mu)*Inf;
                else, z(ci) = 0; end
            end
        end

        iL = iL + 1;
        Lv(iL).pat = p;  Lv(iL).lvl = n;  Lv(iL).cond = cond;  Lv(iL).nn = nn;
        Lv(iL).isHalf = ~isPair;  Lv(iL).pe = pe(n,:);  Lv(iL).t = t(n);
        Lv(iL).trAll = trAll(:)';  Lv(iL).rel = 1:nT;  Lv(iL).valid = valid(:)';
        Lv(iL).metric = metric(:)';  Lv(iL).maxch = maxch(:)';  Lv(iL).baseHz = baseHz(:)';
        Lv(iL).F = F;  Lv(iL).flag = any(F, 1) & valid(:)';  Lv(iL).z = z;
        Lv(iL).win = [a b];
        Lv(iL).nUnrec = nUnrec;  Lv(iL).pct = ternary_(isL1, 100, pctK(cond));  Lv(iL).ev = evK(cond);
        Lv(iL).name = sprintf('%s n=%d%s', patName{p}, nn, ternary_(~isPair, ' (half)', ''));
    end
end
assert(~isempty(Lv), 'No level to check.');
nLv = numel(Lv);

% series (pattern x pair/half) for colours, markers, legends
keys = [Lv.pat] + 100*[Lv.isHalf];
uk = unique(keys);  nS = numel(uk);
sid = arrayfun(@(k) find(uk == k), keys);
for j = 1:nLv, Lv(j).sid = sid(j); end
serName = cell(1,nS);  serHalf = false(1,nS);  serPat = zeros(1,nS);
for s_ = 1:nS
    j = find(sid == s_, 1);
    serPat(s_) = Lv(j).pat;  serHalf(s_) = Lv(j).isHalf;
    if strcmp(patName{serPat(s_)}, 'sim'), serName{s_} = 'sim';
    else, serName{s_} = sprintf('%s %s', patName{serPat(s_)}, ternary_(serHalf(s_), '[n n-1]', '[n n]')); end
end
mk = {'o','d','s','v','p','h'};
baseCols = lines(max(nPat,2));
serCol = zeros(nS,3);  serMk = cell(1,nS);
for s_ = 1:nS
    c = baseCols(serPat(s_),:);
    if serHalf(s_), c = 0.5*c + 0.5*[1 1 1]; serMk{s_} = '^'; else, serMk{s_} = mk{mod(serPat(s_)-1, numel(mk))+1}; end
    serCol(s_,:) = c;
end

%% ================= PRINT: SETTINGS AND TABLE =================
fprintf('\n============================================================\n');
fprintf('NN TRAIN FR BAD-TRIAL CHECK (RandomEvent)\n');
fprintf('============================================================\n');
fprintf('Folder   : %s\nRecovery : %s\n', data_folder, rec_file);
fprintf('Channels : %s (%d) | window: first pulse%+g ms ... last pulse +%g ms | baseline [%g, %g] ms\n', ...
    fmtList_(chUse), nChUse, trainStart_ms, tail_ms, baseWin_ms);
fprintf('Rules    : RESP-MAD(k=%g, floor %g, low side %s) %s | MAXCH-MAD(k=%g, floor %g) %s | BASE-MAD(k=%g, floor %g Hz) %s\n', ...
    respMAD_k, madFloor, onoff_(Resp_FlagLow), onoff_(Use_RespMADRule), maxChMAD_k, maxChFloor, onoff_(Use_MaxChMADRule), ...
    baseMAD_k, baseMadFloorHz, onoff_(Use_BaseMADRule));
fprintf('Patterns : '); for p = 1:nPat, fprintf('%s (final cond %d)   ', patName{p}, patternConds(p)); end; fprintf('\n');
fprintf('Recovered: %s\n', ternary_(haveRec, sprintf('%d of %d trials flagged isRecovered (level 1 is never recovered)', sum(isRec), numel(isRec)), 'no isRecovered in the file'));
fprintf('Levels   : %s (seqLevels = ''%s''). Each level is checked on its own.\n', ternary_(nPat>0,'sim: all n; seq: [n n]',''), seqLevels);

fprintf('\n--- TRIALS PER LEVEL (rel = trial number within the condition) ---\n');
fprintf(' %-9s  n  cond  pulses[E1 E2]  trials unrec  pct  event(ms)  median   flagged   rel trial (abs): reason\n', 'pattern');
for j = 1:nLv
    Lj = Lv(j);  vi = find(Lj.valid);
    fl = find(Lj.flag);
    if isempty(fl), txt = '-';
    else
        pieces = cell(1, numel(fl));
        for q = 1:numel(fl)
            pieces{q} = sprintf('%d (%d): %s', Lj.rel(fl(q)), Lj.trAll(fl(q)), strjoin(ruleNames(Lj.F(:,fl(q))'), '+'));
        end
        txt = strjoin(pieces, ' | ');
    end
    fprintf(' %-9s %2d  %3d    [%2d %2d]       %4d  %4d  %4g   %6.2f   %6.2f   %4d      %s\n', patName{Lj.pat}, Lj.nn, Lj.cond, ...
        Lj.pe(1), Lj.pe(min(2,end)), numel(vi), Lj.nUnrec, Lj.pct, ternary_(isfinite(Lj.ev), Lj.ev, NaN), median(Lj.metric(vi)), numel(fl), ...
        [ternary_(Lj.isHalf,'(half) ','') txt]);
end

if any([Lv.nUnrec] > 0)
    fprintf('\n  NOTE: "unrec" = trials with isRecovered = false (blank gap kept, not used above and not flagged).\n');
end
if any([Lv.pct] < 100)
    fprintf('  NOTE: levels with pct < 100 are randomly thinned: flags there may be caused by the simulation (compare with the unthinned run).\n');
end

fprintf('\n--- EXCLUSION RATE (flagged / valid trials) ---\n');
for p = 1:nPat
    jj = find([Lv.pat] == p);
    nFl = sum(arrayfun(@(j) sum(Lv(j).flag), jj));  nTot = sum(arrayfun(@(j) sum(Lv(j).valid), jj));
    fprintf(' %-9s overall %d / %d (%.1f%%)\n', patName{p}, nFl, nTot, 100*nFl/max(nTot,1));
    for j = jj
        fr = sum(Lv(j).flag) / max(sum(Lv(j).valid), 1);
        if fr > Warn_Exclusion_Fraction
            fprintf('   WARNING: %s loses %.0f%% of its trials (%d / %d)\n', Lv(j).name, 100*fr, sum(Lv(j).flag), sum(Lv(j).valid));
        end
    end
end

fprintf('\n--- RESPONSIVE CHANNELS (z > %g; baseline-subtracted whole-train count, over trials) ---\n', zThr);
for j = 1:nLv
    rc = chUse(Lv(j).z > zThr);
    fprintf(' %-9s n=%2d%s cond %3d  %2d / %d  %s\n', patName{Lv(j).pat}, Lv(j).nn, ternary_(Lv(j).isHalf,'h',' '), ...
        Lv(j).cond, numel(rc), nChUse, fmtList_(rc));
end

fprintf('\n--- PASTE-READY LIST for the bad-trial entry script (rel trial numbers; all trials flagged above) ---\n');
fprintf('%%   pattern      n   cond   trials\n');
for j = 1:nLv
    fprintf('    ''%s'',%s %2d,  %3d,   %s;\n', patName{Lv(j).pat}, repmat(' ', 1, max(1, 11-numel(patName{Lv(j).pat}))), ...
        Lv(j).nn, Lv(j).cond, fmtList_(Lv(j).rel(Lv(j).flag)));
end
fprintf('\n(Nothing was saved. Edit the lists above before using them: they are only what the placeholder rules flagged.)\n');

%% ================= FIGURES: ONE PANEL PER PULSE NUMBER =================
yField = {'metric','maxch','baseHz'};
yLabel = {'mean baseline-subtracted count per channel', 'max single-channel count', 'baseline rate (Hz)'};
mName  = {'RESPONSE', 'MAX single channel', 'BASELINE'};
nnAll = unique([Lv.nn]);  nP = numel(nnAll);
nTmax = max(arrayfun(@(l) numel(l.rel), Lv));
if Plot_PanelFigures
    for m = Trial_Measure
        nC = min(nPanelCols, nP);  nR = ceil(nP/nC);
        fig = figure('Color','w','Units','normalized','Position',[0.01 0.04 0.98 0.88], 'Name', sprintf('Trials per level: %s', mName{m}));
        tl = tiledlayout(fig, nR, nC, 'TileSpacing','compact', 'Padding','compact');
        for ip = 1:nP
            ax = nexttile(tl);  hold(ax, 'on');
            jj = find([Lv.nn] == nnAll(ip));
            hl = [];  ln = {};
            LB = struct('x',{},'y',{},'str',{},'flag',{},'va',{});     % labels, drawn after the axis limits are known
            for j = jj
                Lj = Lv(j);  s_ = Lj.sid;  c = serCol(s_,:);
                off = ((s_-1) - (nS-1)/2) * Series_Shift;               % sim left, seq right
                x = Lj.rel;  y = Lj.(yField{m});
                kept = Lj.valid & ~Lj.flag;  flg = Lj.flag;
                h = plot(ax, x(kept)+off, y(kept), serMk{s_}, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', Marker_Size_Kept, 'LineStyle','none');
                plot(ax, x(flg)+off, y(flg), serMk{s_}, 'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', [1 0 0], ...
                    'MarkerSize', Marker_Size_Flagged, 'LineWidth', 2.5, 'LineStyle','none');
                if Show_Median_Lines && any(Lj.valid)
                    plot(ax, [0.5 nTmax+0.5], median(y(Lj.valid))*[1 1], ':', 'Color', c, 'LineWidth', 1);
                end
                va = ternary_(mod(s_,2) == 1, 'bottom', 'top');           % first series above its points, second below
                for i = find(Lj.valid)
                    if ~Label_All_Trials && ~Lj.flag(i), continue; end
                    str = sprintf('%d', x(i));
                    if Lj.flag(i), str = sprintf('%d %s', x(i), strjoin(codeNames(Lj.F(:,i)'), '')); end
                    LB(end+1) = struct('x', x(i)+off, 'y', y(i), 'str', str, 'flag', Lj.flag(i), 'va', va); %#ok<AGROW>
                end
                hl(end+1) = h(1); ln{end+1} = serName{s_}; %#ok<SAGROW>
            end
            hF = plot(ax, nan, nan, 'o', 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 1 1], 'MarkerSize', Marker_Size_Flagged, 'LineWidth', 2.5, 'LineStyle','none');
            yl = ylim(ax);  r = diff(yl);  if r <= 0, r = max(abs(yl(1)), 1); end
            dy = 0.05*r;
            ylim(ax, [yl(1)-0.10*r, yl(2)+0.10*r]);
            for q = 1:numel(LB)
                sg = ternary_(strcmp(LB(q).va,'bottom'), 1, -1);
                text(ax, LB(q).x, LB(q).y + sg*dy, LB(q).str, 'FontSize', Label_FontSize, 'FontWeight', ternary_(LB(q).flag,'bold','normal'), ...
                    'Color', [0 0 0], 'HorizontalAlignment','center', 'VerticalAlignment', LB(q).va, 'Clipping','on');
            end
            xlim(ax, [0.3 nTmax+0.7]);  box(ax, 'off');
            title(ax, sprintf('%d pulse(s) per electrode', nnAll(ip)), 'FontSize', 10);
            if ip > (nR-1)*nC, xlabel(ax, 'trial number'); end
            if mod(ip-1, nC) == 0, ylabel(ax, yLabel{m}, 'FontSize', 9); end
            if ip == 1, legend(ax, [hl hF], [ln {'flagged'}], 'Box','off', 'FontSize', 8, 'Location','best'); end
        end
        title(tl, sprintf('%s per trial | red ring = flagged (R+ R- response, M max channel, B baseline) | dotted = median | window = whole train + %g ms tail', mName{m}, tail_ms), 'Interpreter','none');
    end
end

%% ================= FIGURE: SUMMARY SCATTER (value of every trial vs number of pulses) =================
if Plot_SummaryFigure
    nM = numel(Trial_Measure);
    figure('Color','w','Position',[60 60 1100 min(900, 330*nM+80)], 'Name', 'Summary: every trial vs number of pulses');
    for mi = 1:nM
        m = Trial_Measure(mi);
        ax = subplot(nM, 1, mi);  hold(ax, 'on');
        hl = gobjects(1,nS);
        for s_ = 1:nS
            off = ((s_-1) - (nS-1)/2) * 0.22;
            c = serCol(s_,:);
            hl(s_) = plot(ax, nan, nan, serMk{s_}, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', Marker_Size_Kept, 'LineStyle','none');
            for j = find([Lv.sid] == s_)
                Lj = Lv(j);  y = Lj.(yField{m});
                x = Lj.nn + off + ((Lj.rel - mean(Lj.rel)) / max(nTmax,1)) * 0.16;
                kept = Lj.valid & ~Lj.flag;
                plot(ax, x(kept), y(kept), serMk{s_}, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', Marker_Size_Kept-1, 'LineStyle','none');
                plot(ax, x(Lj.flag), y(Lj.flag), serMk{s_}, 'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', [1 0 0], 'MarkerSize', Marker_Size_Flagged-1, 'LineWidth', 2, 'LineStyle','none');
                if any(Lj.valid)
                    plot(ax, Lj.nn + off + [-0.1 0.1], median(y(Lj.valid))*[1 1], '-', 'Color', c*0.7, 'LineWidth', 1.8);
                end
            end
        end
        set(ax, 'XTick', nnAll);  xlim(ax, [min(nnAll)-0.6, max(nnAll)+0.6]);
        xlabel(ax, 'number of pulses (per electrode)');  ylabel(ax, yLabel{m});
        title(ax, sprintf('%s of every trial (red ring = flagged, bar = median)', mName{m}));
        box(ax, 'off');
        if mi == 1, legend(ax, hl, serName, 'Box','off', 'Location','northwest'); end
    end
end

%% ================= FIGURE: CHANNEL RESPONSE (z per channel and level) =================
if Plot_ChannelResponse
    figure('Color','w','Position',[80 80 min(1700, 600*nPat+100) 480], 'Name', 'Channel response per level');
    for p = 1:nPat
        ax = subplot(1, nPat, p);  hold(ax, 'on');
        jj = find([Lv.pat] == p);
        for j = jj
            zz = max(min(Lv(j).z, zClip), -zClip);
            scatter(ax, chUse, zz, 22, repmat(Lv(j).nn, 1, nChUse), 'filled', 'MarkerFaceAlpha', 0.7);
        end
        yline(ax, zThr, 'k--', sprintf('z = %g', zThr));  yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
        colormap(ax, parula);  cb = colorbar(ax);  ylabel(cb, 'pulses per electrode');
        caxis(ax, [min(nnAll) max(nnAll)]);
        xlabel(ax, 'channel');  ylabel(ax, 'z of baseline-subtracted whole-train count');
        title(ax, sprintf('%s: one dot per channel and level', patName{p}));
        xlim(ax, [min(chUse)-1 max(chUse)+1]);  box(ax, 'off');
    end
end

%% ================= HELPERS =================
function i = lb_(s, x)
% index of the first element of the sorted vector s that is >= x
lo = 1;  hi = numel(s) + 1;
while lo < hi
    m = floor((lo + hi)/2);
    if s(m) < x, lo = m + 1; else, hi = m; end
end
i = lo;
end

function s = fmtList_(v)
% compact MATLAB-style list, e.g. [1 2 4:7]
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

function s = onoff_(tf)
if tf, s = 'on'; else, s = 'off'; end
end
