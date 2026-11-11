%% ============================================================
% NN_PulseNum_SpikesPerPulse_v1.m
% Spikes per pulse against pulse number, using the whole RECOVERED signal of every level.
%
%   spikes per pulse (level n) = [ spikes in [0, t_n + W] - baseline rate * duration ] / number of pulses
%
%   t_n = time of the last pulse of level n, W = trailing window after it (trailLen_ms),
%   baseline = pre-trigger rate of the same trials (baseWin_ms), number of pulses = pulses per
%   electrode (divideBy = 'perElectrode') or all pulses (divideBy = 'totalPulses').
% If spikes per pulse falls with n, later pulses add less than the average pulse (adaptation /
% saturation). Dividing by the number of pulses removes the trivial growth with train length.
%
% NOTES
%  - In-train parts of a level are copies of lower levels (trial by trial), so the levels are NESTED and
%    the points at different n are not independent. Error bars are per level and do not account for that.
%  - The ~2 ms right after each pulse is filled from the level below (approximation, same for all levels).
%  - SEQUENTIAL patterns: seqLevels = 'matched' uses only the levels where both electrodes have the same
%    number of pulses ([n n]), so sequential [n n] is comparable with simultaneous n.
%  - ERROR BARS  per channel: standard error across trials of the per-trial spikes per pulse.
%                average    : standard error ACROSS CHANNELS.
%
% Levels come from the TrainRecovery file (chain of sourceCond back from a FINAL condition).
% READ-ONLY: plots and prints, saves nothing.
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses';
rec_file    = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses/Xia_Train_Freq4_100Hz_20Pulses.sp_xia_TrainRecovery.mat';

finalConds = [20];       % FINAL condition number(s) of the pattern(s) to analyse ([] = print the list and stop)
chRange    = 1:10;     % channels to analyse: a range like 1:32 or 5:12, or a list like [3 7 9]
plotModes  = {'perChannel'};  % any of:
                       %   'perChannel' : one figure for each channel in chRange
                       %   'grid'       : all channels of chRange in ONE figure (one small panel each)
                       %   'average'    : average over the channels of chRange (raw and normalised, separate figures)
xAxis      = 'perElectrode';  % 'perElectrode' (pulses per electrode), 'time' (time of the last pulse, ms) or 'level'
seqLevels  = 'matched';       % sequential patterns: 'matched' ([n n] levels only) or 'all' (every level)

trailLen_ms = 10;      % trailing window W after the last pulse (ms); a vector, e.g. [20 50 100], draws one curve each
divideBy    = 'perElectrode';  % 'perElectrode' or 'totalPulses'
subtractBaseline = true;       % subtract the pre-trigger rate (times the window duration)
baseWin_ms  = [-50 -5];        % pre-trigger window for the baseline rate

normalize   = 'none';   % 'none' or 'toFirst': second figure divides by each channel's value at the first plotted level
responsiveOnly = false;     % average only over channels with a clear response at the first plotted level
zThr        = 3;           % z-score (mean / SEM of net spikes over trials) needed to count as responsive
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
peK = zeros(nCond, nEl);          % pulses per electrode of every condition
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
assert(ismember(divideBy, {'perElectrode','totalPulses'}), 'divideBy: perElectrode or totalPulses.');
isFinal = arrayfun(@(c) isfield(c,'isFinal') && ~isempty(c.isFinal) && logical(c.isFinal), condTable);

if isempty(finalConds)
    fprintf(' FINAL conditions (use these in finalConds):\n cond  kind     levels  trials  label\n');
    for k = find(isFinal(:))'
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
        warning('SPP:chain', 'cond %d: chain has %d levels but the table says %d.', c, numel(ch_), condTable(c).nLevels);
    end
    chain{p} = ch_;
    tn{p} = arrayfun(@(k) condTable(k).lastPulse_ms, ch_);
end

%% ================= COMPUTE =================
cols = lines(max(nPat,2));
styles = {'o-','s--','^:','d-.','v-'};
nW = numel(trailLen_ms);
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
    if strcmp(ptype,'sequential') && strcmp(seqLevels,'matched'), levels = find(matched); else, levels = 1:L; end
    nPulsesPerEl = sum(pe,2)'/nEl;
    switch divideBy
        case 'perElectrode', divisor = nPulsesPerEl;
        otherwise,           divisor = sum(pe,2)';
    end
    xLevel = nan(1,L);
    for n = 1:L
        switch xAxis
            case 'perElectrode', xLevel(n) = nPulsesPerEl(n);
            case 'time',         xLevel(n) = t(n);
            case 'level',        xLevel(n) = n;
        end
    end

    ser = struct('name',{},'x',{},'Y',{},'SE',{},'style',{},'W',{});
    zFirst = nan(nChn,1);  baseFirst = nan(nChn,1);
    for iw = 1:nW
        W = trailLen_ms(iw);
        Y = nan(nChn, numel(levels));  SE = nan(nChn, numel(levels));
        for ch = chUse
            s = sp{ch}(:,1);
            for j = 1:numel(levels)
                n = levels(j);
                tr = t0(trialCond == cc(n) & isfinite(t0));
                if isempty(tr), continue; end
                dur = t(n) + W;
                c = counts_(s, tr, 0, dur);
                if subtractBaseline
                    b = counts_(s, tr, baseWin_ms(1), baseWin_ms(2)) * (dur / diff(baseWin_ms));
                else
                    b = zeros(size(c));
                end
                net = (c - b) / divisor(n);
                Y(ch,j)  = mean(net);
                SE(ch,j) = std(net) / sqrt(numel(net));
                if iw == 1 && j == 1
                    netc = c - b;
                    se = std(netc)/sqrt(numel(netc));
                    if se > 0, zFirst(ch) = mean(netc)/se; else, zFirst(ch) = sign(mean(netc))*Inf; end
                end
            end
        end
        ser(end+1) = struct('name', sprintf('W=%g ms', W), 'x', xLevel(levels), 'Y', Y, 'SE', SE, ...
            'style', styles{min(iw,numel(styles))}, 'W', W); %#ok<AGROW>
    end
    isResp = zFirst > zThr;
    useCh = intersect(chUse(:)', find(isResp(:))');
    if ~responsiveOnly, useCh = chUse(:)'; end
    res{p} = struct('chain',cc,'t',t,'pe',pe,'ptype',ptype,'levels',levels,'ser',ser,'useCh',useCh,'z',zFirst,'divisor',divisor);
    fprintf('\nPattern %d: final cond %d (%s), %d levels, %s, %d levels plotted. Channels used: %d of %d (%s)\n', p, finalConds(p), ...
        condTable(finalConds(p)).label, L, ptype, numel(levels), numel(useCh), numel(chUse), ternary_(responsiveOnly,'responsive only','all chosen'));
end

%% ================= POPULATION TABLE =================
for p = 1:nPat
    r = res{p};  u = r.useCh;
    for q = 1:numel(r.ser)
        S = r.ser(q);
        fprintf('\nPattern %d (%s), W = %g ms, mean over %d channels: net spikes per pulse (divided by %s)\n', p, r.ptype, S.W, numel(u), divideBy);
        fprintf('  n  cond  pulses[E1 E2]  t_n(ms)  window(ms)         trials  spikes/pulse  sem(ch)\n');
        for j = 1:numel(r.levels)
            n = r.levels(j);
            v = S.Y(u,j);
            fprintf(' %2d  %3d    [%2d %2d]      %7.1f  [0, %7.1f]    %4d   %9.3f   %7.3f\n', n, r.chain(n), r.pe(n,1), r.pe(n,min(2,end)), ...
                r.t(n), r.t(n)+S.W, sum(trialCond == r.chain(n) & isfinite(t0)), mean(v,'omitnan'), ...
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
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            ttl = sprintf('ch %d', ch); if ~ismember(ch, r.useCh), ttl = [ttl ' (not resp.)']; end %#ok<AGROW>
            title(ttl, 'FontSize', 8); box off
            if i > (nr-1)*nc, xlabel(xlabel_(xAxis)); end
            if mod(i-1, nc) == 0, ylabel('spikes per pulse'); end
            if i == 1 && numel(r.ser) > 1, legend({r.ser.name}, 'Box','off', 'FontSize', 6); end
        end
        sgtitle(sprintf('%s  |  spikes per pulse (mean +/- SEM over trials)', condTable(finalConds(p)).label), 'Interpreter','none');
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
                h = errorbar(ax1, S.x, S.Y(ch,:), S.SE(ch,:), S.style, 'Color', cols(p,:), 'LineWidth', 1.2);
                hl(end+1) = h; lg{end+1} = sprintf('cond %d: %s', finalConds(p), S.name); %#ok<SAGROW>
            end
        end
        yline(ax1, 0, ':', 'Color', [0.5 0.5 0.5]);
        ylabel(ax1, 'spikes per pulse');  xlabel(ax1, xlabel_(xAxis));
        title(ax1, sprintf('channel %d: spikes per pulse (mean +/- SEM over trials)', ch));
        legend(ax1, hl, lg, 'Location', 'northeast', 'Box', 'off', 'Interpreter','none');
        box(ax1,'off');
    end
end

%% ================= FIGURE: POPULATION AVERAGE =================
if any(strcmp(plotModes,'average'))
    figure('Color','w','Position',[100 100 1200 420],'Name','population spikes per pulse');
    ax1 = axes; hold(ax1,'on');
    doNorm = strcmp(normalize, 'toFirst');
    ax2 = [];
    if doNorm
        figure('Color','w','Position',[100 560 1200 420],'Name','population spikes per pulse, normalised');
        ax2 = axes; hold(ax2,'on');
    end
    lg = {};  hl = [];  hl2 = [];
    for p = 1:nPat
        r = res{p};  u = r.useCh;
        if isempty(u), warning('SPP:nochan','Pattern %d: no channel left to average.', p); continue; end
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
    ylabel(ax1, 'spikes per pulse');  xlabel(ax1, xlabel_(xAxis));
    title(ax1, 'population: net spikes per pulse (mean +/- SEM across channels)');
    if ~isempty(hl), legend(ax1, hl, lg, 'Location', 'northeast', 'Box', 'off', 'Interpreter','none'); end
    box(ax1, 'off');
    if doNorm
        yline(ax2, 1, ':', 'Color', [0.5 0.5 0.5]);
        ylabel(ax2, 'spikes per pulse / value at first level');  xlabel(ax2, xlabel_(xAxis));
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
