%% ============================================================
% NN_Train_FR_CountVsPulses_UnionRes_v1.m
% NN probe | pulse TRAIN data | FIRING-RATE analysis A | spike COUNT and ADDED SPIKES vs NUMBER OF PULSES | UNION of responsive channels | v1
%
% = NN_Train_FR_CountVsPulses_v1 (untouched) with the channel set changed: by default (chSetMode = 'union') sim and seq are
%   BOTH analysed on the UNION of the responsive channels of all patterns and levels, so that sim and seq are compared
%   on the same channels (a channel responsive in one pattern only is kept for both).
%
% PURPOSE
%   How does the response grow with the number of pulses, and where does it saturate (or drop)?
%   ONE recovery file is analysed per run (change rec_tag / rec_file by hand to look at another file).
%   NOTHING IS SAVED; no existing file is changed.
%
% PER LEVEL (each recovered level is its own condition)
%   Window  : first pulse (+trainStart_ms) ... last pulse of the level + tail_ms        (tail_ms changeable)
%   Count   : spikes in the window minus the baseline count (baseWin_ms) scaled to the window length,
%             mean over the (good) trials of the level, per channel
%   Group   : mean over channels, error bar = SEM over channels (channels are the samples)
%   ADDED   : count(level) - count(previous level), taken per channel, then mean / SEM over channels.
%             First level: its own count. With seqLevels = 'pair' a seq step adds one pulse on EACH electrode;
%             with 'all' every step adds one pulse. Negative added spikes = suppression.
%   RATE    : count / window length (spikes/s)         SPIKES PER PULSE : count / pulses per electrode
%   x axis  : pulses per electrode (sim n; seq [n n]; half-steps n-0.5), or the time of the last pulse (xAxis = 'time')
%
% OPTIONAL INPUT FILES (searched in data_folder; used when present, otherwise everything is kept)
%   <folder name>.FR_RespChn.mat    from NN_Train_FR_RespChnEntry_v1   -> channels per condition
%                                   (not found / pattern not set -> chRange)
%   <folder name>.FR_BadTrials.mat  struct BadTrials with BadTrials.badRel{cond} = trial numbers WITHIN the
%                                   condition (1 ... reps) to leave out; a missing file or an empty / missing entry
%                                   keeps all trials of that condition. (Written by the bad-trial entry script, not yet made.)
%   chSetMode = 'union'    : all patterns use every channel that is responsive in ANY pattern / level   <- default of this script
%               'intersect': all patterns use only the channels responsive in ALL patterns / levels
%               'pattern'  : every pattern uses its own channel set (sim and seq may differ)
%
% FIGURES
%   1. count vs x  and  added spikes vs x   (mean +/- SEM over channels; sim and seq in the same axes)
%   2. rate vs x   and  spikes per pulse vs x      (Plot_SecondFigure)
%   Printed: table of every level (trials used, channels used, count, added spikes, rate).
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
rec_tag     = '';          % '' = *.sp_xia_TrainRecovery.mat (unthinned); 'RandomEvent' etc. = the tagged file
rec_file    = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA/Xia_Train_15Pulses_3_100Hz_30uA.sp_xia_TrainRecovery_RandomEvent.mat';          % full path (overrides rec_tag)
FS          = 30000;       % used only if the recovery file has no FS

patternConds = [];         % FINAL condition number(s) of the patterns; [] = all sim / AB / BA patterns found
chRange      = 1:32;       % manual channel list: used when no responsive-channel file / pattern entry exists
useRespChn   = true;       % true = use <folder>.FR_RespChn.mat when it exists
chSetMode    = 'union';    % 'union' | 'intersect' | 'pattern'
useBadTrials = true;       % true = use <folder>.FR_BadTrials.mat when it exists
resp_file    = '';         % '' = <data_folder>/<folder name>.FR_RespChn.mat
bad_file     = '';         % '' = <data_folder>/<folder name>.FR_BadTrials.mat

% ---- counting window and baseline ----
tail_ms       = 15;        % window = first pulse + trainStart_ms ... last pulse + tail_ms
trainStart_ms = 0;
baseWin_ms    = [-50 -5];  % baseline window relative to the trigger (ms)

% ---- levels and axes ----
seqLevels = 'pair';        % sequential patterns: 'pair' ([n n] only) | 'all' (also the half-steps [n n-1])
xAxis     = 'n';           % 'n' = pulses per electrode | 'time' = time of the last pulse (ms)

% ---- figures ----
Plot_SecondFigure = true;  % rate and spikes per pulse
Show_Channels     = false; % thin line per channel behind the group mean
Marker_Size       = 7;
Line_Width        = 1.6;
Label_FontSize    = 11;

%% ================= CHECK SETTINGS =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(ismember(seqLevels, {'pair','all'}), 'seqLevels: pair or all.');
assert(ismember(xAxis, {'n','time'}), 'xAxis: n or time.');
assert(ismember(chSetMode, {'union','intersect','pattern'}), 'chSetMode: union, intersect or pattern.');
assert(numel(baseWin_ms) == 2 && baseWin_ms(2) > baseWin_ms(1), 'baseWin_ms must be [start end].');
assert(tail_ms >= 0, 'tail_ms must be >= 0.');

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
    isRec = logical(R.isRecovered(:));
else
    isRec = true(nTrialsAll, 1);
end
clear R

ft = dir(fullfile(data_folder, '*.trig.dat'));  ft = ft(~startsWith({ft.name}, '.'));
assert(~isempty(ft), 'No *.trig.dat in %s.', data_folder);
fid = fopen(fullfile(data_folder, ft(1).name), 'r');  trig = fread(fid, [1 inf], 'double');  fclose(fid);
t0 = nan(nTrialsAll, 1);  nUse = min(numel(trig), nTrialsAll);  t0(1:nUse) = trig(1:nUse)/FS*1000;

fe = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));  fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fullfile(data_folder, fe(1).name), 'StimParams', 'simultaneous_stim');
nEl = double(E.simultaneous_stim);
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

[~, folderName] = fileparts(regexprep(data_folder, '[\\/]+$', ''));
% ---- optional: responsive channels ----
haveResp = false;  RespChn = [];
if useRespChn
    if isempty(resp_file), resp_file = fullfile(data_folder, [folderName '.FR_RespChn.mat']); end
    if isfile(resp_file)
        S_ = load(resp_file, 'RespChn');  RespChn = S_.RespChn;  haveResp = true;  clear S_
        assert(numel(RespChn.respChn) == nCond, 'The responsive-channel file has %d conditions, the recovery file %d.', numel(RespChn.respChn), nCond);
    end
end
% ---- optional: bad trials ----
haveBad = false;  BadTrials = [];
if useBadTrials
    if isempty(bad_file), bad_file = fullfile(data_folder, [folderName '.FR_BadTrials.mat']); end
    if isfile(bad_file)
        S_ = load(bad_file, 'BadTrials');  BadTrials = S_.BadTrials;  haveBad = true;  clear S_
        assert(isfield(BadTrials, 'badRel'), 'The bad-trial file must hold BadTrials.badRel{cond}.');
    end
end

%% ================= PATTERNS AND LEVELS =================
kinds = arrayfun(@(c) char(c.kind), condTable, 'UniformOutput', false);
isFinal = arrayfun(@(c) isfield(c,'isFinal') && ~isempty(c.isFinal) && logical(c.isFinal), condTable);
if isempty(patternConds)
    patternConds = find(isFinal(:)' & ismember(kinds(:)', {'sim','AB','BA'}));
end
assert(~isempty(patternConds), 'No sim / AB / BA pattern found in the recovery file.');
nPat = numel(patternConds);
chain = cell(nPat,1);  patName = cell(nPat,1);  lvPos = cell(nPat,1);
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
    isPair = arrayfun(@(k) all(peK(k,:) == peK(k,1)) && peK(k,1) > 0, ch_);
    if ~strcmp(kinds{c}, 'sim') && strcmp(seqLevels, 'all')
        lvPos{p} = find(sum(peK(ch_,:), 2)' > 0);       % every step (half-steps included)
    else
        lvPos{p} = find(isPair);                          % sim: every level; seq: the [n n] levels
    end
end

%% ================= CHANNEL SETS PER LEVEL =================
chSet = cell(nPat, 1);          % chSet{p}{q} = channels of level q
usedFile = false(nPat, 1);
for p = 1:nPat
    chSet{p} = cell(1, numel(lvPos{p}));
    for q = 1:numel(lvPos{p})
        k = chain{p}(lvPos{p}(q));
        if haveResp && RespChn.isSet(k)
            chSet{p}{q} = intersect(double(RespChn.respChn{k}(:)'), 1:nChn);  usedFile(p) = true;
        else
            chSet{p}{q} = intersect(chRange(:)', 1:nChn);
        end
    end
end
if strcmp(chSetMode, 'intersect')
    common = 1:nChn;
    for p = 1:nPat, for q = 1:numel(chSet{p}), common = intersect(common, chSet{p}{q}); end, end
    common = common(:)';
    for p = 1:nPat, for q = 1:numel(chSet{p}), chSet{p}{q} = common; end, end
elseif strcmp(chSetMode, 'union')
    common = zeros(1, 0);
    for p = 1:nPat, for q = 1:numel(chSet{p}), common = union(common, chSet{p}{q}); end, end
    common = common(:)';                                   % row vector (union of an empty start returns a column)
    for p = 1:nPat, for q = 1:numel(chSet{p}), chSet{p}{q} = common; end, end
end
if haveResp && ismember(chSetMode, {'union','intersect'}) && ~all(usedFile)
    warning('UR:notset', 'Pattern(s) %s have no entry in the responsive-channel file: chRange was used for them, which widens the %s set.', ...
        strjoin(patName(~usedFile), ', '), chSetMode);
end

%% ================= PER-LEVEL VALUES =================
baseLen = diff(baseWin_ms);
Lv = struct([]);
for p = 1:nPat
    cc = chain{p};  t = arrayfun(@(k) double(condTable(k).lastPulse_ms), cc);
    for q = 1:numel(lvPos{p})
        n = lvPos{p}(q);  cond = cc(n);
        all_ = find(trialCond == cond);
        isL1 = ~isfinite(condTable(cond).sourceCond);
        ok = isfinite(t0(all_)) & (isRec(all_) | isL1);
        nBad = 0;
        if haveBad && numel(BadTrials.badRel) >= cond && ~isempty(BadTrials.badRel{cond})
            br = BadTrials.badRel{cond}(:)';  br = br(br >= 1 & br <= numel(all_));
            nBad = sum(ok(br));  ok(br) = false;
        end
        tr = all_(ok);  tt = t0(tr);
        a = t(1) + trainStart_ms;  b = t(n) + tail_ms;  scaleF = (b - a)/baseLen;
        chs = chSet{p}{q};
        cm = nan(1, nChn);
        if numel(tr) >= 1
            for ch = chs
                s = sp{ch}(:,1);  v = zeros(numel(tt), 1);
                for i = 1:numel(tt)
                    v(i) = (lb_(s, tt(i)+b) - lb_(s, tt(i)+a)) - (lb_(s, tt(i)+baseWin_ms(2)) - lb_(s, tt(i)+baseWin_ms(1))) * scaleF;
                end
                cm(ch) = mean(v);
            end
        end
        Lv(p).cond(q) = cond;                     %#ok<SAGROW>
        Lv(p).pe(q)   = mean(peK(cond,:));        %#ok<SAGROW>
        Lv(p).t(q)    = t(n);                     %#ok<SAGROW>
        Lv(p).win(q)  = b - a;                    %#ok<SAGROW>
        Lv(p).nTr(q)  = numel(tr);                %#ok<SAGROW>
        Lv(p).nBad(q) = nBad;                     %#ok<SAGROW>
        Lv(p).cm(q,:) = cm;                       %#ok<SAGROW>
    end
    % added spikes: difference of consecutive levels, per channel
    cmP = Lv(p).cm;
    Lv(p).add = [cmP(1,:); cmP(2:end,:) - cmP(1:end-1,:)];
    Lv(p).rate = Lv(p).cm ./ (Lv(p).win(:) / 1000);
    Lv(p).spp  = Lv(p).cm ./ Lv(p).pe(:);
    if strcmp(xAxis, 'n'), Lv(p).x = Lv(p).pe; else, Lv(p).x = Lv(p).t; end
end

%% ================= PRINT =================
fprintf('\n============================================================\n');
fprintf('NN TRAIN FR | COUNT AND ADDED SPIKES VS PULSES | UNION RES (v1)\n');
fprintf('============================================================\n');
fprintf('Folder   : %s\nRecovery : %s\n', data_folder, rec_file);
fprintf('Window   : first pulse%+g ms ... last pulse +%g ms | baseline [%g, %g] ms | levels: seq ''%s'' | x axis: %s\n', ...
    trainStart_ms, tail_ms, baseWin_ms, seqLevels, xAxis);
if haveResp, fprintf('Channels : responsive-channel file %s (chSetMode ''%s'')\n', resp_file, chSetMode);
    if ~strcmp(chSetMode, 'pattern'), fprintf('           common set used for every pattern and level (%d channels): %s\n', numel(chSet{1}{1}), fmtList_(chSet{1}{1})); end
else,        fprintf('Channels : no responsive-channel file used -> chRange %s (chSetMode ''%s'')\n', fmtList_(chRange), chSetMode); end
if haveBad,  fprintf('Bad trials: file %s\n', bad_file);
else,        fprintf('Bad trials: no bad-trial file used -> all trials kept\n'); end
for p = 1:nPat
    fprintf('\n--- %s (final cond %d) ---\n', upper(patName{p}), patternConds(p));
    fprintf('  cond   x     trials(bad)  chn   count    +/-SEM    added    +/-SEM    rate(Hz)  spikes/pulse\n');
    for q = 1:numel(Lv(p).cond)
        [mC, sC, nC] = ms_(Lv(p).cm(q,:));  [mA, sA] = ms_(Lv(p).add(q,:));  mR = ms_(Lv(p).rate(q,:));  mS = ms_(Lv(p).spp(q,:));
        fprintf('  %3d  %5.1f   %4d (%d)    %3d  %7.2f  %7.2f  %7.2f  %7.2f   %8.1f   %8.3f\n', Lv(p).cond(q), Lv(p).x(q), Lv(p).nTr(q), Lv(p).nBad(q), nC, mC, sC, mA, sA, mR, mS);
    end
end
fprintf('\n(Nothing was saved.)\n');

%% ================= FIGURES =================
cols = [0.00 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556];
mks = {'o','s','d','^'};
panels = {'cm', 'add'; 'rate', 'spp'};
ylab = struct('cm', 'baseline-subtracted count per trial (whole train + tail)', 'add', 'added spikes per level (count n - count n-1)', ...
    'rate', 'rate (spikes/s)', 'spp', 'spikes per pulse (count / pulses per electrode)');
xlab = ternary_(strcmp(xAxis, 'n'), 'pulses per electrode', 'time of the last pulse (ms)');
nFig = 1 + double(Plot_SecondFigure);
for fi = 1:nFig
    figure('Color','w','Position',[60+40*fi 60+30*fi 1150 460], 'Name', ternary_(fi == 1, 'Count and added spikes', 'Rate and spikes per pulse'));
    for pi_ = 1:2
        fld = panels{fi, pi_};
        ax = subplot(1, 2, pi_);  hold(ax, 'on');  hl = gobjects(1, nPat);
        for p = 1:nPat
            c = cols(mod(p-1, size(cols,1))+1, :);  mk = mks{mod(p-1, numel(mks))+1};
            Y = Lv(p).(fld);  x = Lv(p).x(:)';
            if Show_Channels
                plot(ax, x, Y', '-', 'Color', [c 0.25], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            end
            m = nan(1, numel(x));  se = nan(1, numel(x));
            for q = 1:numel(x), [m(q), se(q)] = ms_(Y(q,:)); end
            hl(p) = errorbar(ax, x, m, se, ['-' mk], 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', Marker_Size, 'LineWidth', Line_Width, 'CapSize', 4);
        end
        if strcmp(fld, 'add') || strcmp(fld, 'cm'), yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]); end
        xlabel(ax, xlab, 'FontSize', Label_FontSize);  ylabel(ax, ylab.(fld), 'FontSize', Label_FontSize);
        if strcmp(xAxis, 'n'), xlim(ax, [min(cellfun(@min, {Lv.x}))-0.5, max(cellfun(@max, {Lv.x}))+0.5]); end
        legend(ax, hl, cellfun(@(nm, p) sprintf('%s (%d chn)', nm, round(mean(sum(isfinite(Lv(p).cm), 2)))), patName(:)', num2cell(1:nPat), 'UniformOutput', false), ...
            'Box', 'off', 'Location', 'best');
        box(ax, 'off');  set(ax, 'FontSize', Label_FontSize);
    end
    sgtitle(sprintf('%s | mean +/- SEM over channels | channels: %s | bad trials: %s', strrep(rec_file(find(rec_file == filesep, 1, 'last')+1:end), '_', '\_'), ...
        ternary_(haveResp, 'responsive-channel file', ['chRange ' fmtList_(chRange)]), ternary_(haveBad, 'file', 'none')), 'FontSize', 10);
end

%% ================= HELPERS =================
function [m, se, n] = ms_(v)
v = v(isfinite(v));  n = numel(v);
if n == 0, m = NaN;  se = NaN;  return; end
m = mean(v);
if n > 1, se = std(v)/sqrt(n); else, se = NaN; end
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
