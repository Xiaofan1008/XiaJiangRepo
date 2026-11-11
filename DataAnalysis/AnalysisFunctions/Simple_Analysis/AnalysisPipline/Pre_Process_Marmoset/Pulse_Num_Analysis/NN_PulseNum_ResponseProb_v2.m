%% ============================================================
% NN_PulseNum_ResponseProb_v2.m      (v2 = v1 + seqMode for sequential patterns)
% Response PROBABILITY against pulse number.
%
%   P(level n) = fraction of trials with at least one spike in [t_n + winStart_ms, t_n + winStart_ms + W]
%   (t_n = time of the last pulse of level n), for one channel.
% A baseline probability p0 is measured on the same trials in a window of the same length just before the
% trigger ([-5 - W, -5] ms). With correctBaseline = true the plotted value is (P - p0) / (1 - p0): the
% probability that the pulses caused a spike, with chance spikes removed (0 = no response, 1 = every trial).
%
% SEQUENTIAL patterns (A pulse, then B pulse 5 ms later) - seqMode:
%   'AandB'      (default) the A response and the B response are measured SEPARATELY and then averaged:
%                  A response: level where A is the last pulse (level 2n-1), window after that A pulse
%                  B response: level where B is the last pulse (level 2n),   window after that B pulse
%                  plotted value = mean of the two probabilities (probabilities are averaged, not added,
%                  because a sum could exceed 1); error = sqrt(seA^2 + seB^2)/2.
%                  x = pulses per electrode, like simultaneous n.
%   'pairWindow' one window on level 2n anchored at the START of the pair (the A pulse):
%                  [t_A + winStart_ms, t_A + winStart_ms + W]; for simultaneous the same with t_n.
%                  It contains both responses, like the simultaneous window.
%   'lastPulse'  window after the last pulse only (the B pulse for [n n]); seqLevels picks 'matched' or 'all'.
% Simultaneous patterns are not affected (one window after each simultaneous pair).
%
% Only the genuine part of each level is used (after the blank), so levels are independent of each other.
% It is the response after the n-th pulse in the context of the n-1 before it.
% A spike COUNT can stay constant while the probability of responding falls (bursts), so this can show
% fatigue that counts hide.
%
% ERROR BARS  per channel: binomial standard error sqrt(P(1-P)/N) (combined with that of p0 when corrected);
%             average    : standard error ACROSS CHANNELS.
% With 15 trials a single channel's probability moves in steps of 1/15, so rely on the population average.
%
% Levels come from the TrainRecovery file (chain of sourceCond back from a FINAL condition).
% READ-ONLY: plots and prints, saves nothing.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses';
rec_file    = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses/Xia_Train_Freq4_100Hz_20Pulses.sp_xia_TrainRecovery.mat';

finalConds = [60];       % FINAL condition number(s) of the pattern(s) ([] = print the list of train patterns and stop)
levelRange = [];    % Sim 1 12 17 24 28 32 Seq 5 12 21    % levels to analyse, e.g. 1:10 ([] = all levels of the pattern)
chRange    = [5 12 21];     % channels to analyse: a range like 1:32 or 5:12, or a list like [3 7 9]
plotModes  = {'average'};  % any of 'perChannel' | 'grid' | 'average'
xAxis      = 'perElectrode';  % 'perElectrode' (pulses per electrode), 'time' (time of the last pulse, ms) or 'level'
seqMode    = 'AandB';         % sequential patterns: 'AandB' | 'pairWindow' | 'lastPulse' (see header)
seqLevels  = 'matched';       % only for seqMode = 'lastPulse': 'matched' ([n n] levels) or 'all'

winStart_ms = 2;       % window start after the last pulse (ms); should be longer than the blank
winLen_ms   = 20;      % window length W (ms); a vector, e.g. [5 10 20], draws one curve each
correctBaseline = true;% plot (P - p0)/(1 - p0) instead of the raw P
baseEnd_ms  = -5;      % the baseline window is [baseEnd_ms - W, baseEnd_ms] before the trigger

normalize   = 'none';  % 'none' or 'toFirst': extra figure divided by each channel's value at the first plotted level
responsiveOnly = false; % average only over channels with a clear response at the first plotted level
zThr        = 3;       % z-score needed to count as responsive
FS = 30000;

%% ================= LOAD =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
cd(data_folder);
assert(isfile(rec_file), 'Cannot find %s (run NN_SpikeRecovery_v1 first).', rec_file);
R = load(rec_file);
assert(all(isfield(R, {'sp_clipped','Recovery_params'})), '%s must be a TrainRecovery file.', rec_file);
sp = R.sp_clipped;
if isfield(R,'FS'), FS = R.FS; end
trialCond = R.Recovery_params.trialCond(:);
condTable = R.Recovery_params.condTable;
nChn = numel(sp);  n_Trials = numel(trialCond);  nCond = numel(condTable);
clear R

f = dir('*.trig.dat'); f = f(~startsWith({f.name}, '.'));
assert(~isempty(f), 'No *.trig.dat in %s.', data_folder);
fid = fopen(f(1).name, 'r'); trig = fread(fid, [1 inf], 'double'); fclose(fid);
t0 = nan(n_Trials,1);  nUse = min(numel(trig), n_Trials);  t0(1:nUse) = trig(1:nUse)/FS*1000;

fe = dir('*_exp_datafile_*.mat'); fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fe(1).name, 'StimParams', 'simultaneous_stim');
nEl = E.simultaneous_stim;
peK = zeros(nCond, nEl);
for k = 1:nCond
    t1 = find(trialCond == k, 1);
    if isempty(t1), continue; end
    for e = 1:nEl
        row = 1 + (t1-1)*nEl + e;
        if double(E.StimParams{row,16}) > 0, peK(k,e) = double(E.StimParams{row,8}); end
    end
end

chUse = chRange(:)';
assert(~isempty(chUse) && all(chUse >= 1 & chUse <= nChn), 'chRange must lie between 1 and %d.', nChn);
if ischar(plotModes), plotModes = {plotModes}; end
assert(all(ismember(plotModes, {'perChannel','grid','average'})), 'plotModes: use perChannel, grid and/or average.');
assert(ismember(xAxis, {'perElectrode','time','level'}), 'xAxis: perElectrode, time or level.');
assert(ismember(seqLevels, {'matched','all'}), 'seqLevels: matched or all.');
assert(ismember(seqMode, {'AandB','pairWindow','lastPulse'}), 'seqMode: AandB, pairWindow or lastPulse.');
isFinal = arrayfun(@(c) isfield(c,'isFinal') && ~isempty(c.isFinal) && logical(c.isFinal), condTable);
nLev = arrayfun(@(c) c.nLevels, condTable);

if isempty(finalConds)
    fprintf(' Pulse-train patterns (use the cond number in finalConds):\n cond  kind     levels  trials  label\n');
    for k = find(isFinal(:)' & nLev(:)' > 1)
        fprintf(' %3d   %-7s  %3d     %4d    %s\n', k, condTable(k).kind, condTable(k).nLevels, sum(trialCond==k), condTable(k).label);
    end
    fprintf('\nSet finalConds at the top, then run again.\n');
    return;
end

%% ================= CHAINS =================
nPat = numel(finalConds);
chain = cell(nPat,1);  tn = cell(nPat,1);
for p = 1:nPat
    c = finalConds(p);
    ch_ = c;
    while isfinite(condTable(ch_(1)).sourceCond)
        ch_ = [condTable(ch_(1)).sourceCond, ch_]; %#ok<AGROW>
    end
    if numel(ch_) ~= condTable(c).nLevels
        warning('RP:chain', 'cond %d: chain has %d levels but the table says %d.', c, numel(ch_), condTable(c).nLevels);
    end
    chain{p} = ch_;
    tn{p} = arrayfun(@(k) condTable(k).lastPulse_ms, ch_);
end

%% ================= COMPUTE =================
cols = lines(max(nPat,2));
styles = {'o-','s--','^:','d-.','v-'};
nW = numel(winLen_ms);
res = cell(nPat,1);
for p = 1:nPat
    cc = chain{p};  t = tn{p};  L = numel(cc);
    pe = peK(cc,:);
    prev = [zeros(1,nEl); pe(1:L-1,:)];
    inc = pe > prev;
    isSingleStep = arrayfun(@(n) sum(inc(n,:)) == 1, 1:L);
    if all(~isSingleStep), ptype = 'simultaneous';
    elseif all(isSingleStep), ptype = 'sequential';
    else, ptype = 'mixed'; end
    if nEl == 2, matched = (pe(:,1) == pe(:,2))'; else, matched = true(1,L); end
    isSeq2 = strcmp(ptype,'sequential') && nEl == 2;
    if isSeq2 && (strcmp(seqMode,'AandB') || strcmp(seqMode,'pairWindow') || strcmp(seqLevels,'matched'))
        levels = find(matched);
    else
        levels = 1:L;
    end
    if isSeq2 && ~strcmp(seqMode,'lastPulse') && any(levels < 2), levels = levels(levels >= 2); end
    if ~isempty(levelRange), levels = intersect(levels, levelRange(:)'); end
    assert(~isempty(levels), 'Pattern %d: no level left after levelRange / seqLevels.', p);
    nPulsesPerEl = sum(pe,2)'/nEl;
    xLevel = nan(1,L);
    for n = 1:L
        switch xAxis
            case 'perElectrode', xLevel(n) = nPulsesPerEl(n);
            case 'time',         xLevel(n) = t(n);
            case 'level',        xLevel(n) = n;
        end
    end

    ser = struct('name',{},'x',{},'Y',{},'SE',{},'style',{},'W',{});
    zFirst = nan(nChn,1);
    winStr = cell(1, numel(levels));
    for iw = 1:nW
        W = winLen_ms(iw);
        Y = nan(nChn, numel(levels));  SE = nan(nChn, numel(levels));
        for ch = chUse
            s = sp{ch}(:,1);
            for j = 1:numel(levels)
                n = levels(j);
                if isSeq2 && strcmp(seqMode,'AandB')
                    % A response (level n-1, window after its last pulse) and B response (level n, window after its last pulse)
                    trA = t0(trialCond == cc(n-1) & isfinite(t0));
                    trB = t0(trialCond == cc(n)   & isfinite(t0));
                    if isempty(trA) || isempty(trB), continue; end
                    [vA, seA, zA] = probStats_(s, trA, t(n-1)+winStart_ms, W, baseEnd_ms, correctBaseline);
                    [vB, seB, zB] = probStats_(s, trB, t(n)  +winStart_ms, W, baseEnd_ms, correctBaseline);
                    Y(ch,j)  = (vA + vB)/2;
                    SE(ch,j) = sqrt(seA^2 + seB^2)/2;
                    zj = max(zA, zB);
                    if iw == 1, winStr{j} = sprintf('A[%.1f,%.1f] B[%.1f,%.1f]', t(n-1)+winStart_ms, t(n-1)+winStart_ms+W, t(n)+winStart_ms, t(n)+winStart_ms+W); end
                else
                    if isSeq2 && strcmp(seqMode,'pairWindow'), tA = t(n-1); else, tA = t(n); end
                    tr = t0(trialCond == cc(n) & isfinite(t0));
                    if isempty(tr), continue; end
                    [Y(ch,j), SE(ch,j), zj] = probStats_(s, tr, tA+winStart_ms, W, baseEnd_ms, correctBaseline);
                    if iw == 1, winStr{j} = sprintf('[%.1f,%.1f]', tA+winStart_ms, tA+winStart_ms+W); end
                end
                if iw == 1 && j == 1, zFirst(ch) = zj; end
            end
        end
        ser(end+1) = struct('name', sprintf('W=%g ms', W), 'x', xLevel(levels), 'Y', Y, 'SE', SE, ...
            'style', styles{min(iw,numel(styles))}, 'W', W); %#ok<AGROW>
    end
    isResp = zFirst > zThr;
    useCh = intersect(chUse(:)', find(isResp(:))');
    if ~responsiveOnly, useCh = chUse(:)'; end
    res{p} = struct('chain',cc,'t',t,'pe',pe,'ptype',ptype,'levels',levels,'ser',ser,'useCh',useCh,'winStr',{winStr});
    modeStr = ptype; if isSeq2, modeStr = [ptype ', seqMode ' seqMode]; end
    fprintf('\nPattern %d: final cond %d (%s), %d levels, %s, %d levels plotted. Channels used: %d of %d (%s)\n', p, finalConds(p), ...
        condTable(finalConds(p)).label, L, modeStr, numel(levels), numel(useCh), numel(chUse), ternary_(responsiveOnly,'responsive only','all chosen'));
end

%% ================= POPULATION TABLE =================
yname = ternary_(correctBaseline, 'baseline-corrected probability', 'probability');
for p = 1:nPat
    r = res{p};  u = r.useCh;
    for q = 1:numel(r.ser)
        S = r.ser(q);
        fprintf('\nPattern %d (%s), W = %g ms, mean over %d channels: %s\n', p, r.ptype, S.W, numel(u), yname);
        fprintf('  n  cond  pulses[E1 E2]  t_n(ms)  window(ms)                        trials  probability  sem(ch)\n');
        for j = 1:numel(r.levels)
            n = r.levels(j);
            v = S.Y(u,j);
            fprintf(' %2d  %3d    [%2d %2d]      %7.1f  %-32s  %4d   %9.3f   %7.3f\n', n, r.chain(n), r.pe(n,1), r.pe(n,min(2,end)), ...
                r.t(n), r.winStr{j}, sum(trialCond == r.chain(n) & isfinite(t0)), mean(v,'omitnan'), ...
                std(v,'omitnan')/sqrt(max(sum(~isnan(v)),1)));
        end
    end
end

%% ================= FIGURE: GRID =================
if any(strcmp(plotModes,'grid'))
    for p = 1:nPat
        r = res{p};
        nr = ceil(sqrt(numel(chUse)));  nc = ceil(numel(chUse)/nr);
        figure('Color','w','Position',[50 50 1300 800],'Name',sprintf('pattern %d grid', p));
        for i = 1:numel(chUse)
            ch = chUse(i);
            subplot(nr, nc, i);  hold on
            for q = 1:numel(r.ser)
                S = r.ser(q);
                errorbar(S.x, S.Y(ch,:), S.SE(ch,:), S.style, 'MarkerSize', 3, 'LineWidth', 0.8);
            end
            yline(0, ':', 'Color', [0.5 0.5 0.5]);  ylim([-0.2 1.1]);
            ttl = sprintf('ch %d', ch); if ~ismember(ch, r.useCh), ttl = [ttl ' (not resp.)']; end %#ok<AGROW>
            title(ttl, 'FontSize', 8); box off
            if i > (nr-1)*nc, xlabel(xlabel_(xAxis)); end
            if mod(i-1, nc) == 0, ylabel('probability'); end
            if i == 1 && numel(r.ser) > 1, legend({r.ser.name}, 'Box','off', 'FontSize', 6); end
        end
        sgtitle(sprintf('%s  |  response probability (mean +/- binomial SE)', condTable(finalConds(p)).label), 'Interpreter','none');
    end
end

%% ================= FIGURE: ONE FIGURE PER CHANNEL =================
if any(strcmp(plotModes,'perChannel'))
    for ch = chUse(:)'
        figure('Color','w','Position',[100 100 1200 420],'Name',sprintf('channel %d', ch));
        ax1 = axes; hold(ax1,'on');
        lg = {};  hl = [];
        for p = 1:nPat
            r = res{p};
            for q = 1:numel(r.ser)
                S = r.ser(q);
                hl(end+1) = errorbar(ax1, S.x, S.Y(ch,:), S.SE(ch,:), S.style, 'Color', cols(p,:), 'LineWidth', 1.2); %#ok<SAGROW>
                lg{end+1} = sprintf('cond %d: %s', finalConds(p), S.name); %#ok<SAGROW>
            end
        end
        yline(ax1, 0, ':', 'Color', [0.5 0.5 0.5]);
        ylabel(ax1, yname);  xlabel(ax1, xlabel_(xAxis));
        title(ax1, sprintf('channel %d: response probability (mean +/- binomial SE)', ch));
        legend(ax1, hl, lg, 'Location', 'northeast', 'Box', 'off', 'Interpreter','none');
        box(ax1,'off');
    end
end

%% ================= FIGURE: POPULATION AVERAGE =================
if any(strcmp(plotModes,'average'))
    figure('Color','w','Position',[100 100 1200 420],'Name','population response probability');
    ax1 = axes; hold(ax1,'on');
    doNorm = strcmp(normalize, 'toFirst');
    ax2 = [];
    if doNorm
        figure('Color','w','Position',[100 560 1200 420],'Name','population response probability, normalised');
        ax2 = axes; hold(ax2,'on');
    end
    lg = {};  hl = [];  hl2 = [];
    for p = 1:nPat
        r = res{p};  u = r.useCh;
        if isempty(u), warning('RP:nochan','Pattern %d: no channel left to average.', p); continue; end
        for q = 1:numel(r.ser)
            S = r.ser(q);
            A = S.Y(u,:);
            mu = mean(A,1,'omitnan');  se = std(A,0,1,'omitnan')./sqrt(max(sum(~isnan(A),1),1));
            hl(end+1) = errorbar(ax1, S.x, mu, se, S.style, 'Color', cols(p,:), 'LineWidth', 1.3); %#ok<SAGROW>
            lg{end+1} = sprintf('cond %d: %s (%d ch)', finalConds(p), S.name, numel(u)); %#ok<SAGROW>
            if doNorm
                ok = A(:,1) > 0;
                An = A(ok,:) ./ A(ok,1);
                mn = mean(An,1,'omitnan');  sn = std(An,0,1,'omitnan')./sqrt(max(sum(~isnan(An),1),1));
                hl2(end+1) = errorbar(ax2, S.x, mn, sn, S.style, 'Color', cols(p,:), 'LineWidth', 1.3); %#ok<SAGROW>
            end
        end
    end
    yline(ax1, 0, ':', 'Color', [0.5 0.5 0.5]);
    ylabel(ax1, yname);  xlabel(ax1, xlabel_(xAxis));
    title(ax1, 'population: response probability (mean +/- SEM across channels)');
    if ~isempty(hl), legend(ax1, hl, lg, 'Location', 'northeast', 'Box', 'off', 'Interpreter','none'); end
    box(ax1, 'off');
    if doNorm
        yline(ax2, 1, ':', 'Color', [0.5 0.5 0.5]);
        ylabel(ax2, 'probability / value at first level');  xlabel(ax2, xlabel_(xAxis));
        title(ax2, 'population: normalised to each channel''s first plotted level (mean +/- SEM across channels)');
        if ~isempty(hl2), legend(ax2, hl2, lg(1:numel(hl2)), 'Location', 'northeast', 'Box', 'off', 'Interpreter','none'); end
        box(ax2, 'off');
    end
end

%% ================= HELPERS =================
function s = xlabel_(xAxis)
switch xAxis
    case 'perElectrode', s = 'pulses per electrode';
    case 'time',         s = 'time of the last pulse (ms)';
    otherwise,           s = 'level (pulse number)';
end
end

function [val, se, z] = probStats_(s, tr, a, W, baseEnd, correct)
% probability of >=1 spike in [a, a+W] (ms after trigger) over the trials tr, with baseline correction
N  = numel(tr);
P  = mean(counts_(s, tr, a, a + W) > 0);
P0 = mean(counts_(s, tr, baseEnd - W, baseEnd) > 0);
seP = sqrt(P*(1-P)/N);  seP0 = sqrt(P0*(1-P0)/N);
if correct && P0 < 1
    val = (P - P0)/(1 - P0);  se = sqrt(seP^2 + seP0^2)/(1 - P0);
else
    val = P;  se = seP;
end
den = sqrt(seP^2 + seP0^2);
if den > 0, z = (P - P0)/den; else, z = sign(P - P0)*Inf; end
end

function c = counts_(s, t0s, a, b)
% spike count in [t0+a, t0+b) for each trigger (s sorted ascending column)
c = zeros(numel(t0s),1);
for i = 1:numel(t0s)
    c(i) = lb_(s, t0s(i)+b) - lb_(s, t0s(i)+a);
end
end

function i = lb_(s, x)
lo = 1; hi = numel(s) + 1;
while lo < hi
    m = floor((lo + hi)/2);
    if s(m) < x, lo = m + 1; else, hi = m; end
end
i = lo;
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
