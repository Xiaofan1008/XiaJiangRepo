%% ============================================================
% NP_Train_FR_RespChnEntry_v1.m
% NP probe | pulse TRAIN sessions | FIRING-RATE analysis | RESPONSIVE-CHANNEL ENTRY | v1
%
% PURPOSE
%   Reads the automatic responsive lists of NP_Train_FR_RespChnCheck_v1 (<short>.NP_RespChn_proposal.mat, written by the
%   check when saveProposal = true), applies YOUR corrections per pattern, prints and plots the result and, when
%   doSave = true, saves ONE file that the FR analysis scripts read (NP_Train_FR_RateVsPulses_v1, useRespChn = true).
%   Nothing is recomputed here: edit "resp", re-run (instant), look at the printout / figure, repeat, then save.
%   The recording, the spike files and the proposal are never changed.
%
% PATTERNS (one row each in "resp")
%   'sim', 'seq B->A', 'seq A->B'  : the names the check prints; EVERY level of the pattern gets the same set
%   'visual bright', 'visual dark' : one shared set for all contrasts of that polarity (bright / dark taken from the
%                                    condition label)
%        visStart = 'final' : 'auto' = the automatic list of the HIGHEST contrast of that polarity (default)
%        visStart = 'union' : 'auto' = the union of the automatic lists of all contrasts of that polarity
%   Blank condition(s) get an empty set. A pattern that is not typed gets an empty set (isReviewed = false, warning).
%
% STEPS, per pattern, in this order
%   1. start  : 'auto' (the proposal) | 'none' (empty) | 'all' (every channel TESTED for this pattern) | a list, e.g. [10:44 60]
%   'Tested' is per pattern (the check: response-window rate >= autoRateMin on the decision level); visual bright / dark use
%   the tested set of the decision condition (visStart = 'final') or the union over the contrasts (visStart = 'union').
%   2. keep   : keep only channels whose depth (um from the tip) is inside one of the [from to] rows (inclusive);
%               [] = no depth limit
%   3. add    : channels added; they are KEPT even outside the keep ranges
%   4. remove : channels removed; remove always wins
%   Channel lists accept MATLAB ranges (100:120). Channels that were not tested by the check can be added (a note is printed).
%
% PRINTED per pattern: start count, channels dropped by depth, added, removed, final count per depth band, the final list,
%   and the difference (+ / -) from the saved NP_RespChn file if it exists.
% FIGURE (plotFigure = true): one panel per pattern, depth vs the check's number (z) of the decision condition:
%   grey = tested, not chosen | green = final set | orange ring = dropped by depth | blue triangle = added |
%   red x = removed | blue shading = keep ranges | dashed = the check's threshold.
%
% OUTPUT (doSave = true)  <recording1>/<short>.NP_RespChn.mat  with the struct RespChn:
%   RespChn.respChn{cond}     final channels of every condition (sorted row)   <- read by the FR scripts
%   RespChn.isReviewed(cond)  true = the condition's pattern was typed in resp
%   RespChn.pattern(p) also keeps the check's numbers of the decision condition, per channel (index = channel):
%                             score (z), rBase, rResp (Hz), nPos (trials with d > 0), tested (channel list);
%                             RespChn.scoreName / scoreThr describe the score  <- read by NP_Train_FR_RespChnDepthPlot_v1
%   RespChn.pattern(p)        name, conds, decCond, how, auto, start, keep, add, remove, startList, droppedDepth,
%                             added, removed, final, isReviewed
%   RespChn.resp              the table as typed;  visStart, proposalFile, proposalDate, ruleSettings, chTest, y_um,
%                             note, ttFile, date, version
%   An existing file is first copied to <name>_backup_<yyyymmdd_HHMMSS>.mat, then replaced.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;
codeVersion = 'NP_Train_FR_RespChnEntry_v1 (2026-11-14)';

%% ================= USER SETTINGS =================
sessionDir   = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
proposalFile = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_RespChn_proposal.mat';          % '' = <short>.NP_RespChn_proposal.mat in the recording folder (or a full path)
visStart     = 'final';     % 'final' = 'auto' is the highest contrast of each polarity | 'union' = all contrasts of the polarity

% pattern | start | keep depth ranges (um, rows of [from to]) | add | remove
resp = { ...
    'sim',           'auto',  [],   [314,282,303],   [233,251];
    'seq B->A',      'auto',  [],   [264,262,260,190,182,94,90,78],   [];
    'visual bright', 'auto',  [],   [],   [];
    'visual dark',   'auto',  [],   [],   [];
    };

doSave     = true;         % false = DRY RUN (printout + figure only) | true = save <short>.NP_RespChn.mat
note       = '';            % free text stored in the file
plotFigure = true;          % true = the probe figure (one panel per pattern)
groupLines = [2200 3000 6000];   % depth bands for the printed counts and the lines on the figure (um); [] = none
%% =================================================

assert(any(strcmp(visStart, {'final', 'union'})), 'visStart must be ''final'' or ''union''.');
assert(iscell(resp) && size(resp, 2) == 5, 'resp must have 5 columns: pattern, start, keep, add, remove.');

%% ---- load ---------------------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');  ttFile = fullfile(recDir, tf(1).name);
L = load(ttFile, 'TT');  TT = L.TT;  clear L
if isempty(proposalFile), proposalFile = [short '.NP_RespChn_proposal.mat']; end
if ~isfile(proposalFile), proposalFile = fullfile(recDir, proposalFile); end
assert(isfile(proposalFile), ['Proposal not found: %s\nRun NP_Train_FR_RespChnCheck_v1 with saveProposal = true first ' ...
       '(with the settings you want, e.g. rateSDmode = ''paired'').'], proposalFile);
L = load(proposalFile, 'Prop');  Prop = L.Prop;  clear L
assert(strcmp(Prop.short, short), 'The proposal belongs to %s, not to %s.', Prop.short, short);
nCond = numel(TT.cond);  nCh = Prop.nCh;  y = double(Prop.y_um(:));  chTest = reshape(Prop.chTest, 1, []);
PP = Prop.pattern;
outFile = fullfile(recDir, [short '.NP_RespChn.mat']);

%% ---- patterns of this session (electrical as in the check; visual grouped by polarity) -------------------------------
Q = struct('name', {}, 'conds', {}, 'decCond', {}, 'auto', {}, 'score', {}, 'how', {}, 'tested', {}, 'rBase', {}, 'rResp', {}, 'nPos', {});
hasT = isfield(PP, 'tested');                                      % older proposals: one tested set for the whole session
if ~hasT, warning('NPRE:oldProp', 'Old proposal without per-pattern tested sets: every pattern uses Prop.chTest. Re-run the check.'); end
for p = find(~[PP.isVis])
    q = struct();  q.name = PP(p).name;  q.conds = reshape(PP(p).conds, 1, []);  q.decCond = PP(p).decCond;
    q.auto = rowv_(PP(p).auto);  q.score = PP(p).score(:);  q.how = sprintf('automatic list of cond %d', PP(p).decCond);
    if hasT, q.tested = rowv_(PP(p).tested); else, q.tested = chTest; end
    q.rBase = PP(p).rBase(:);  q.rResp = PP(p).rResp(:);  q.nPos = PP(p).nPos(:);
    Q(end+1) = q; %#ok<SAGROW>
end
pv = find([PP.isVis]);
if ~isempty(pv)
    pol = cell(1, numel(pv));
    for i = 1:numel(pv)
        c = PP(pv(i)).decCond;  lb = lower(TT.cond(c).label);
        if contains(lb, 'bright'), pol{i} = 'bright';
        elseif contains(lb, 'dark'), pol{i} = 'dark';
        else, error('Cannot tell bright or dark from the label of cond %d: "%s".', c, TT.cond(c).label);
        end
    end
    for pn = {'bright', 'dark'}
        ii = pv(strcmp(pol, pn{1}));
        if isempty(ii), continue; end
        cs = [PP(ii).decCond];  ctr = arrayfun(@(c) double(TT.cond(c).contrast), cs);
        [~, k] = max(ctr);  dc = cs(k);                                   % highest contrast of this polarity
        q = struct();  q.name = ['visual ' pn{1}];  q.conds = sort(cs);  q.decCond = dc;  q.auto = [];  q.score = PP(ii(k)).score(:);
        if strcmp(visStart, 'final')
            q.auto = rowv_(PP(ii(k)).auto);  q.how = sprintf('automatic list of cond %d (highest contrast)', dc);
            if hasT, q.tested = rowv_(PP(ii(k)).tested); else, q.tested = chTest; end
        else
            q.auto = rowv_(unique([PP(ii).auto]));  q.how = sprintf('union of the automatic lists of conds %s', mat2str(sort(cs)));
            if hasT, q.tested = rowv_(unique([PP(ii).tested])); else, q.tested = chTest; end
        end
        q.rBase = PP(ii(k)).rBase(:);  q.rResp = PP(ii(k)).rResp(:);  q.nPos = PP(ii(k)).nPos(:);   % numbers of the decision cond
        Q(end+1) = q; %#ok<SAGROW>
    end
end
nQ = numel(Q);

% match the typed rows to the patterns
names = lower(strtrim({Q.name}));  rowOf = zeros(1, nQ);
for i = 1:size(resp, 1)
    k = find(strcmp(names, lower(strtrim(resp{i, 1}))));
    if isempty(k)
        error('resp row %d: pattern "%s" not found. Patterns of this session: %s', i, resp{i, 1}, ...
              strjoin(cellfun(@(s) ['''' s ''''], {Q.name}, 'UniformOutput', false), ', '));
    end
    assert(rowOf(k) == 0, 'Pattern "%s" is typed twice in resp.', Q(k).name);
    rowOf(k) = i;
end

%% ---- apply the steps -------------------------------------------------------------------------------------------------
res = struct('name', {}, 'conds', {}, 'decCond', {}, 'how', {}, 'auto', {}, 'start', {}, 'keep', {}, 'add', {}, ...
             'remove', {}, 'startList', {}, 'droppedDepth', {}, 'added', {}, 'removed', {}, 'final', {}, 'isReviewed', {}, 'tested', {}, 'score', {}, 'rBase', {}, 'rResp', {}, 'nPos', {});
for k = 1:nQ
    q = Q(k);  i = rowOf(k);
    if i == 0
        warning('NPRE:notTyped', 'Pattern "%s" is not typed in resp: its conditions get an EMPTY set (not reviewed).', q.name);
        st = 'none';  kp = [];  ad = [];  rm = [];  rev = false;
    else
        st = resp{i, 2};  kp = resp{i, 3};  ad = resp{i, 4};  rm = resp{i, 5};  rev = true;
    end
    % 1. start
    if ischar(st) || isstring(st)
        switch lower(char(st))
            case 'auto', s0 = q.auto;
            case 'none', s0 = zeros(1, 0);
            case 'all',  s0 = q.tested;
            otherwise, error('Pattern "%s": start must be ''auto'', ''none'', ''all'' or a channel list.', q.name);
        end
    else
        s0 = chkCh_(st, nCh, q.name, 'start');
    end
    s0 = rowv_(unique(s0));
    % 2. keep depth ranges
    if isempty(kp)
        inK = true(1, numel(s0));
    else
        assert(isnumeric(kp) && size(kp, 2) == 2, 'Pattern "%s": keep must be [] or rows of [from to] in um.', q.name);
        inK = reshape(inRanges_(y(s0), kp), 1, []);
    end
    dropD = s0(~inK);  s1 = s0(inK);
    % 3. add (kept even outside the keep ranges)
    ad = chkCh_(ad, nCh, q.name, 'add');  added = rowv_(setdiff(ad, s1));  s2 = rowv_(union(s1, ad));
    % 4. remove (always wins)
    rm = chkCh_(rm, nCh, q.name, 'remove');  removed = rowv_(intersect(s2, rm));  fin = rowv_(setdiff(s2, rm));

    r = struct();  r.name = q.name;  r.conds = q.conds;  r.decCond = q.decCond;  r.how = q.how;  r.auto = q.auto;
    r.start = st;  r.keep = kp;  r.add = ad;  r.remove = rm;  r.startList = s0;  r.droppedDepth = rowv_(dropD);
    r.added = added;  r.removed = removed;  r.final = fin;  r.isReviewed = rev;
    r.tested = q.tested;  r.score = q.score;  r.rBase = q.rBase;  r.rResp = q.rResp;  r.nPos = q.nPos;   % for the depth plot
    res(k) = r;
end

% one set per condition
respChn = repmat({zeros(1, 0)}, 1, nCond);  isReviewed = false(1, nCond);  owner = zeros(1, nCond);
for k = 1:nQ
    for c = res(k).conds
        assert(owner(c) == 0, 'Condition %d belongs to two patterns (%s, %s).', c, res(max(owner(c), 1)).name, res(k).name);
        owner(c) = k;  respChn{c} = res(k).final;  isReviewed(c) = res(k).isReviewed;
    end
end

%% ---- printout --------------------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP TRAIN FR RESPONSIVE-CHANNEL ENTRY (v1): %s\n', short);
fprintf('============================================================\n');
S = Prop.settings;
fprintf('proposal   : %s (made %s)\n', proposalFile, Prop.date);
if strcmp(S.ruleType, 'rate')
    fprintf('rule       : rate, %s, z > %g, %s level, baseline %s ms\n', S.rateSDmode, S.rateZ, S.decideBy, mat2str(S.baseWinMs));
else
    fprintf('rule       : psth peak rule, %s level\n', S.decideBy);
end
fprintf('processed  : %d channels (tested: per pattern, see below) | visual start: %s\n', numel(chTest), visStart);
bandE = [-Inf, groupLines(:)', Inf];
for k = 1:nQ
    r = res(k);
    fprintf('\n--- %s (conds %s) ---\n', r.name, mat2str(r.conds));
    if ~r.isReviewed, fprintf('  NOT TYPED in resp -> EMPTY set\n'); continue; end
    if ischar(r.start) || isstring(r.start), stTxt = char(r.start); else, stTxt = 'list'; end
    fprintf('  start  : %-5s -> %3d channels   (%s)\n', stTxt, numel(r.startList), ternary_(strcmp(stTxt, 'auto'), r.how, 'typed'));
    if isempty(r.keep), fprintf('  keep   : no depth limit\n');
    else, fprintf('  keep   : %s um -> %d dropped by depth: %s\n', keepTxt_(r.keep), numel(r.droppedDepth), chList_(r.droppedDepth)); end
    fprintf('  add    : %3d new %s\n', numel(r.added), chList_(r.added));
    nt = setdiff(r.added, r.tested);
    if ~isempty(nt), fprintf('           note: not tested for this pattern (response rate too low / not processed): %s\n', chList_(nt)); end
    fprintf('  remove : %3d     %s\n', numel(r.removed), chList_(r.removed));
    fprintf('  FINAL  : %3d channels\n', numel(r.final));
    for b = numel(bandE) - 1:-1:1
        inB = @(v) v(y(v) >= bandE(b) & y(v) < bandE(b + 1));
        fprintf('    depth %5s-%5s um: %3d of %3d tested\n', num2s_(bandE(b)), num2s_(bandE(b + 1)), numel(inB(r.final)), numel(inB(r.tested)));
    end
    fprintf('  list   : %s\n', chList_(r.final));
end
cNone = find(owner == 0);
if ~isempty(cNone)
    fprintf('\nconditions without a pattern (empty set): %s\n', strjoin(arrayfun(@(c) sprintf('%d (%s)', c, TT.cond(c).label), ...
            cNone, 'UniformOutput', false), ', '));
end

% difference from the saved file
if isfile(outFile)
    O = load(outFile, 'RespChn');  old = O.RespChn.respChn;
    fprintf('\nCOMPARED WITH THE SAVED FILE (%s):\n', outFile);
    for k = 1:nQ
        c = res(k).decCond;
        if c > numel(old), fprintf('  %-14s not in the saved file\n', res(k).name); continue; end
        pl = setdiff(res(k).final, old{c});  mi = setdiff(old{c}, res(k).final);
        if isempty(pl) && isempty(mi), fprintf('  %-14s unchanged\n', res(k).name);
        else, fprintf('  %-14s + %s   - %s\n', res(k).name, chList_(pl), chList_(mi)); end
    end
end

%% ---- figure ----------------------------------------------------------------------------------------------------------
if plotFigure
    fh = figure('Name', [short ' responsive-channel entry'], 'Color', 'w', 'Position', [40 40 min(1800, 330*nQ + 80) 780]);
    hL = gobjects(0);  lab = {};
    for k = 1:nQ
        r = res(k);  s = double(Q(k).score(:));  ax = subplot(1, nQ, k);  hold(ax, 'on');
        if numel(s) < nCh, s(end+1:nCh) = NaN; end
        fs = s(r.tested);  fs = fs(isfinite(fs));
        if isempty(fs), xl = [-1 1]; else, xl = [min([fs; 0]) max([fs; Prop.scoreThr])]; end
        if diff(xl) == 0, xl = xl + [-1 1]; end
        xl = xl + [-1 1]*0.05*diff(xl);
        s(s == Inf) = xl(2);  s(s == -Inf) = xl(1);  s(isnan(s)) = xl(1);   % off-scale / untested at the edges
        yl = [min(y) max(y)];
        for rr = 1:size(r.keep, 1)                                        % keep ranges
            a = max(min(r.keep(rr, :)), yl(1));  b = min(max(r.keep(rr, :)), yl(2));
            if b > a, patch(ax, xl([1 2 2 1]), [a a b b], [0.85 0.92 1], 'EdgeColor', 'none'); end
        end
        for g = groupLines(:)', plot(ax, xl, [g g], '-', 'Color', [0.75 0.75 0.75]); end
        plot(ax, [Prop.scoreThr Prop.scoreThr], yl, 'g--');
        other = setdiff(chTest, [r.final r.droppedDepth r.removed]);
        plot(ax, s(other), y(other), '.', 'Color', [0.6 0.6 0.6]);
        plot(ax, s(r.final), y(r.final), 'o', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 3);
        plot(ax, s(r.droppedDepth), y(r.droppedDepth), 'o', 'Color', [1 0.5 0], 'MarkerSize', 5);
        plot(ax, s(r.added), y(r.added), '^', 'Color', [0 0.3 1], 'MarkerFaceColor', [0 0.3 1], 'MarkerSize', 5);
        plot(ax, s(r.removed), y(r.removed), 'x', 'Color', [0.9 0 0], 'MarkerSize', 7, 'LineWidth', 1.2);
        if k == 1                                                   % legend from dummy lines: always all 5 entries,
            hL = [plot(ax, NaN, NaN, '.', 'Color', [0.6 0.6 0.6]), ...  % even when a group is empty in this panel
                  plot(ax, NaN, NaN, 'o', 'Color', [0 0.6 0], 'MarkerFaceColor', [0 0.6 0], 'MarkerSize', 3), ...
                  plot(ax, NaN, NaN, 'o', 'Color', [1 0.5 0], 'MarkerSize', 5), ...
                  plot(ax, NaN, NaN, '^', 'Color', [0 0.3 1], 'MarkerFaceColor', [0 0.3 1], 'MarkerSize', 5), ...
                  plot(ax, NaN, NaN, 'x', 'Color', [0.9 0 0], 'MarkerSize', 7, 'LineWidth', 1.2)];
            lab = {'tested, not chosen', 'final set', 'dropped by depth', 'added', 'removed'};
        end
        xlim(ax, xl);  ylim(ax, yl);  grid(ax, 'on');
        xlabel(ax, sprintf('%s, cond %d', Prop.scoreName, r.decCond), 'FontSize', 8);
        if k == 1, ylabel(ax, 'depth from the tip (um)'); else, set(ax, 'YTickLabel', []); end
        title(ax, sprintf('%s\n%d channels%s', r.name, numel(r.final), ternary_(r.isReviewed, '', ' (NOT TYPED)')), ...
              'FontSize', 9, 'Interpreter', 'none');
    end
    if ~isempty(hL), legend(subplot(1, nQ, 1), hL, lab, 'FontSize', 7, 'Location', 'southeast'); end
    sgtitle(fh, sprintf('%s: responsive channels per pattern (%s)', strrep(short, '_', '\_'), ternary_(doSave, 'SAVED', 'DRY RUN')));
end

%% ---- save ------------------------------------------------------------------------------------------------------------
if doSave
    if isfile(outFile)
        bak = strrep(outFile, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
        copyfile(outFile, bak);  fprintf('\nOld file copied to %s\n', bak);
    end
    RespChn = struct();
    RespChn.respChn = respChn;  RespChn.isReviewed = isReviewed;  RespChn.pattern = res;  RespChn.resp = resp;
    RespChn.visStart = visStart;  RespChn.proposalFile = proposalFile;  RespChn.proposalDate = Prop.date;
    RespChn.scoreName = Prop.scoreName;  RespChn.scoreThr = Prop.scoreThr;
    RespChn.ruleSettings = Prop.settings;  RespChn.chTest = chTest;  RespChn.y_um = y;  RespChn.note = note;
    RespChn.ttFile = ttFile;  RespChn.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');  RespChn.version = codeVersion;
    save(outFile, 'RespChn');
    fprintf('\nSAVED: %s\n', outFile);
else
    fprintf('\nDRY RUN: nothing saved. Set doSave = true to save %s\n', outFile);
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function v = rowv_(v)
v = reshape(v, 1, []);
end

function v = chkCh_(v, nCh, nm, what)
% channel list check: integers 1..nCh, sorted, unique, row
if isempty(v), v = zeros(1, 0); return; end
assert(isnumeric(v) && all(v(:) == round(v(:))) && all(v(:) >= 1 & v(:) <= nCh), ...
       'Pattern "%s": %s must be channel numbers between 1 and %d.', nm, what, nCh);
v = reshape(unique(v(:)), 1, []);
end

function m = inRanges_(v, R)
% true where v is inside any [from to] row of R (inclusive)
m = false(size(v));
for r = 1:size(R, 1), m = m | (v >= min(R(r, :)) & v <= max(R(r, :))); end
end

function s = keepTxt_(R)
s = strjoin(arrayfun(@(r) sprintf('%g-%g', R(r, 1), R(r, 2)), 1:size(R, 1), 'UniformOutput', false), ', ');
end

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
