%% ============================================================
% NN_PulseNum_AddedSpikes_v1.m      (measure "B")
% How many spikes does the n-th pulse ADD, and until which n does it stop adding?
%
%   added(n) = [spikes of level n in a window after pulse n] - [spikes of level n-1 in the SAME
%              window, same time after the trigger]
% Level n-1 has no n-th pulse, so the difference is what pulse n added, given the n-1 before it.
% The window starts at  t_n + winStart_ms  (t_n = time of the n-th pulse) and lasts winLen_ms.
% Both levels are genuine original data there (after the blank), so this does not depend on
% how the recovery filled the train. Spontaneous activity cancels in the difference.
% Level 1 has no lower level: its "added" value is its response in the window minus the
% pre-trigger baseline (switch: includePulse1).
%
% SIMULTANEOUS vs SEQUENTIAL are handled differently, detected from the pulses of each level:
%   simultaneous : every step adds a pulse to BOTH electrodes at once  -> one curve "steps"
%   sequential   : every step adds ONE pulse, alternating electrodes   -> curve "E1 steps" (odd
%                  levels, e.g. A), curve "E2 steps" (even levels, e.g. B), and curve "pair [n n]":
%                  the two steps that complete [n n] added together (comparable with simultaneous n)
% Cumulative added spikes are given at the levels where both electrodes have the same number of
% pulses ([n n]; every level for simultaneous), so simultaneous n and sequential [n n] match.
%
% CHANNELS / PLOTS: chRange picks the channels (e.g. 1:32, 5:12, [3 7 9]); plotModes picks
% 'perChannel' (one figure each), 'grid' (all in one figure) and/or 'average' (mean over chRange).
% Levels are taken from the TrainRecovery file (chain of sourceCond back from a FINAL condition).
% READ-ONLY: plots and prints, saves nothing.
%
% WINDOW OPTIONS
%   winStart_ms : start after each pulse (ms), should be longer than the blank (see diagnostics)
%   winLen_ms   : length of the window for each pulse number n, ONE of
%                   scalar            same length for every pulse, e.g. 10
%                   vector            entry n = length for pulse n (entry 1 only used if includePulse1;
%                                     if shorter than the number of levels the last entry repeats)
%                   'nextPulse'       until the next pulse of the pattern (t_{n+1} - t_n - winStart_ms)
%                                     for the last pulse winLastPulse_ms is used
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses';
rec_file    = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses/Xia_Train_Freq4_100Hz_20Pulses.sp_xia_TrainRecovery.mat';

finalConds = [20 60];       % FINAL condition number(s) of the pattern(s) to analyse ([] = print the list and stop)
chRange    = 1:32;     % channels to analyse: a range like 1:32 or 5:12, or a list like [3 7 9]
plotModes  = {'perChannel'};  % any of:
                       %   'perChannel' : one figure for each channel in chRange
                       %   'grid'       : all channels of chRange in ONE figure (one small panel each)
                       %   'average'    : one figure with the average over the channels of chRange
                       % e.g. {'perChannel','average'}
xAxis      = 'perElectrode';  % x axis: 'perElectrode' (pulses per electrode), 'time' (time of the pulse, ms) or 'level'

winStart_ms     = 2;   % window start after each pulse (ms)
winLen_ms       = 10;  % scalar | vector (per pulse number) | 'nextPulse'
winLastPulse_ms = 10;  % window length for the last pulse when winLen_ms = 'nextPulse'
includePulse1   = true;% also count the response to the first pulse (level 1 minus baseline); needed so that
                       % cumulative curves of simultaneous n and sequential [n n] match

normalize   = 'toPulse1';  % 'none' or 'toPulse1': also show added / (channel's response to the pattern's first step)
normLen_ms  = 10;          % window after the first pulse used for that response
baseWin_ms  = [-50 -5];    % pre-trigger window for the baseline rate
responsiveOnly = true;     % average only over channels responding to the first pulse
respLen_ms  = 10;          % window after pulse 1 (starting at winStart_ms) used for the test
zThr        = 3;           % z-score above baseline needed to count as responsive

showDiagnostics = true;    % figure with: blank_end distribution, level-1 PSTH, difference curves
nExplore   = [2 3 5 10 20];% pulse numbers shown in the difference curves
binMs      = 1;            % PSTH bin for the diagnostics
FS = 30000;

%% ================= LOAD =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
cd(data_folder);
assert(isfile(rec_file), 'Cannot find %s (run NN_SpikeRecovery_v1 first).', rec_file);
R = load(rec_file);
need = {'sp_clipped','blank_end_ms','Recovery_params'};
assert(all(isfield(R, need)), '%s must be a TrainRecovery file.', rec_file);
sp = R.sp_clipped;  be = R.blank_end_ms;
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

chUse = chRange(:)';
assert(~isempty(chUse) && all(chUse >= 1 & chUse <= nChn), 'chRange must lie between 1 and %d.', nChn);
if ischar(plotModes), plotModes = {plotModes}; end
assert(all(ismember(plotModes, {'perChannel','grid','average'})), 'plotModes: use perChannel, grid and/or average.');
assert(ismember(xAxis, {'perElectrode','time','level'}), 'xAxis: perElectrode, time or level.');
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
        warning('PN:chain', 'cond %d: chain has %d levels but the table says %d.', c, numel(ch_), condTable(c).nLevels);
    end
    chain{p} = ch_;
    tn{p} = arrayfun(@(k) condTable(k).lastPulse_ms, ch_);     % time of the n-th pulse (ms)
    assert(all(diff(tn{p}) > 0), 'cond %d: pulse times of the levels are not increasing.', c);
end

%% ================= COMPUTE =================
cols = lines(max(nPat,2));
res = cell(nPat,1);
for p = 1:nPat
    cc = chain{p};  t = tn{p};  L = numel(cc);
    W = nan(1,L);
    for n = 1:L, W(n) = winLen_(n, L, t, winStart_ms, winLen_ms, winLastPulse_ms); end
    nt = zeros(1, L);
    for n = 1:L, nt(n) = sum(trialCond == cc(n) & isfinite(t0)); end

    % which electrode(s) does each step add a pulse to?
    pe   = peK(cc,:);                          % L x nEl pulses per electrode
    prev = [zeros(1,nEl); pe(1:L-1,:)];
    inc  = pe > prev;
    stepElec = zeros(1,L);                     % 0 = both electrodes, e = electrode e only
    for n = 1:L
        if sum(inc(n,:)) == 1, stepElec(n) = find(inc(n,:)); end
        if sum(inc(n,:)) == 0, warning('PN:step', 'pattern %d level %d adds no pulse.', p, n); end
    end
    if all(stepElec == 0), ptype = 'simultaneous';
    elseif all(stepElec > 0), ptype = 'sequential';
    else, ptype = 'mixed'; end
    firstLevel = 1 + ~includePulse1;

    add = nan(nChn, L);  sem = nan(nChn, L);
    for ch = chUse
        s = sp{ch}(:,1);
        for n = firstLevel:L
            if ~(W(n) > 0), continue; end
            a = t(n) + winStart_ms;  b = a + W(n);
            cn = counts_(s, t0(trialCond == cc(n) & isfinite(t0)), a, b);
            if isempty(cn), continue; end
            if n == 1
                bl = counts_(s, t0(trialCond == cc(1) & isfinite(t0)), baseWin_ms(1), baseWin_ms(2)) * (W(1)/diff(baseWin_ms));
                add(ch,1) = mean(cn) - mean(bl);
                sem(ch,1) = sqrt(var(cn)/numel(cn) + var(bl)/numel(bl));
            else
                cn1 = counts_(s, t0(trialCond == cc(n-1) & isfinite(t0)), a, b);
                if isempty(cn1), continue; end
                add(ch,n) = mean(cn) - mean(cn1);
                sem(ch,n) = sqrt(var(cn)/numel(cn) + var(cn1)/numel(cn1));
            end
        end
    end

    % response to the first step (normalisation / responsive channels)
    tr1 = t0(trialCond == cc(1) & isfinite(t0));
    R1 = nan(nChn,1);  zResp = nan(nChn,1);
    for ch = chUse
        s = sp{ch}(:,1);
        a = t(1) + winStart_ms;
        rN = counts_(s, tr1, a, a + normLen_ms);
        rZ = counts_(s, tr1, a, a + respLen_ms);
        bc = counts_(s, tr1, baseWin_ms(1), baseWin_ms(2));
        bL = diff(baseWin_ms);
        if isempty(rN), continue; end
        R1(ch) = mean(rN) - mean(bc * (normLen_ms / bL));
        bz = bc * (respLen_ms / bL);
        dz = mean(rZ) - mean(bz);
        se = sqrt(var(rZ)/numel(rZ) + var(bz)/numel(bz));
        if se > 0, zResp(ch) = dz/se; else, zResp(ch) = sign(dz)*Inf; end
    end
    isResp = zResp > zThr;
    useCh = intersect(chUse(:)', find(isResp(:))');
    if ~responsiveOnly, useCh = chUse(:)'; end

    % ---- x values for a level, series, cumulative ----
    xLevel = nan(1,L);
    for n = 1:L
        switch xAxis
            case 'perElectrode', if stepElec(n) > 0, xLevel(n) = pe(n, stepElec(n)); else, xLevel(n) = pe(n,1); end
            case 'time',         xLevel(n) = t(n);
            case 'level',        xLevel(n) = n;
        end
    end
    ser = struct('name',{},'x',{},'Y',{},'SE',{},'style',{});
    switch ptype
        case 'simultaneous'
            idx = firstLevel:L;
            ser(end+1) = struct('name','steps', 'x',xLevel(idx), 'Y',add(:,idx), 'SE',sem(:,idx), 'style','o-');
        otherwise
            sty = {'o-','s--','^:'};
            for e = unique(stepElec(stepElec>0))
                idx = find(stepElec == e & (1:L) >= firstLevel);
                ser(end+1) = struct('name',sprintf('E%d steps', e), 'x',xLevel(idx), 'Y',add(:,idx), 'SE',sem(:,idx), 'style',sty{min(e,3)}); %#ok<AGROW>
            end
            if strcmp(ptype,'sequential') && nEl == 2
                Lc = find(pe(:,1) == pe(:,2) & (1:L)' >= firstLevel+1)';
                Lc = Lc(stepElec(Lc-1) > 0);
                if ~isempty(Lc)
                    xp = nan(1,numel(Lc));
                    for j = 1:numel(Lc)
                        switch xAxis
                            case 'perElectrode', xp(j) = pe(Lc(j),1);
                            case 'time',         xp(j) = t(Lc(j));
                            case 'level',        xp(j) = Lc(j);
                        end
                    end
                    ser(end+1) = struct('name','pair [n n]', 'x',xp, ...
                        'Y',add(:,Lc-1)+add(:,Lc), 'SE',sqrt(sem(:,Lc-1).^2+sem(:,Lc).^2), 'style','d-'); %#ok<AGROW>
                end
            end
    end
    % cumulative added spikes, at levels with equal pulses on both electrodes ([n n]; every level if simultaneous)
    if nEl == 2, complete = find(pe(:,1) == pe(:,2))'; else, complete = 1:L; end
    complete = complete(complete >= firstLevel);
    cumY = nan(nChn, numel(complete));
    for j = 1:numel(complete)
        cumY(:,j) = sum(add(:, firstLevel:complete(j)), 2, 'omitnan');
    end
    cumX = xLevel(complete);
    if ~strcmp(ptype,'simultaneous') && strcmp(xAxis,'perElectrode')
        cumX = pe(complete,1)';
    end
    res{p} = struct('chain',cc,'t',t,'W',W,'nTrials',nt,'add',add,'sem',sem,'R1',R1,'z',zResp,'useCh',useCh, ...
        'pe',pe,'stepElec',stepElec,'ptype',ptype,'ser',ser,'cumX',cumX,'cumY',cumY,'complete',complete,'firstLevel',firstLevel);
    fprintf('\nPattern %d: final cond %d (%s), %d levels, %s. Channels used: %d of %d (%s)\n', p, finalConds(p), ...
        condTable(finalConds(p)).label, L, ptype, numel(useCh), numel(chUse), ternary_(responsiveOnly,'responsive only','all chosen'));
end

%% ================= POPULATION TABLE =================
for p = 1:nPat
    r = res{p};  L = numel(r.chain);  u = r.useCh;
    fprintf('\nPattern %d (%s), population mean over %d channels (spikes added per trial and channel)\n', p, r.ptype, numel(u));
    fprintf('  n  cond  pulses[E1 E2]  step  t_n(ms)  window(ms)          trials  added   sem(ch)   cumulative\n');
    cum = 0;
    for n = r.firstLevel:L
        v = r.add(u,n);
        mu = mean(v, 'omitnan');  se = std(v, 'omitnan')/sqrt(max(sum(~isnan(v)),1));
        cum = cum + mu;
        stp = 'both'; if r.stepElec(n) > 0, stp = sprintf('E%d', r.stepElec(n)); end
        fprintf(' %2d  %3d    [%2d %2d]       %-4s  %7.1f  [%6.1f, %6.1f]   %4d   %6.3f  %6.3f    %7.3f\n', n, r.chain(n), ...
            r.pe(n,1), r.pe(n,min(2,end)), stp, r.t(n), r.t(n)+winStart_ms, r.t(n)+winStart_ms+r.W(n), r.nTrials(n), mu, se, cum);
    end
end

%% ================= FIGURE: GRID (all channels of chRange in one figure) =================
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
            if mod(i-1, nc) == 0, ylabel('added spikes'); end
            if i == 1, legend({r.ser.name}, 'Box','off', 'FontSize', 6); end
        end
        sgtitle(sprintf('%s  |  spikes added by pulse n (mean +/- SEM over trials)', condTable(finalConds(p)).label), 'Interpreter','none');
    end
end

%% ================= FIGURE: ONE FIGURE PER CHANNEL =================
if any(strcmp(plotModes,'perChannel'))
    for ch = chUse(:)'
        figure('Color','w','Position',[100 100 800 650],'Name',sprintf('channel %d', ch));
        ax1 = subplot(2,1,1); hold(ax1,'on');  ax2 = subplot(2,1,2); hold(ax2,'on');
        lg = {};  hl = [];
        for p = 1:nPat
            r = res{p};
            for q = 1:numel(r.ser)
                S = r.ser(q);
                h = errorbar(ax1, S.x, S.Y(ch,:), S.SE(ch,:), S.style, 'Color', cols(p,:), 'LineWidth', 1.2);
                hl(end+1) = h; lg{end+1} = sprintf('cond %d: %s', finalConds(p), S.name); %#ok<SAGROW>
            end
            plot(ax2, r.cumX, r.cumY(ch,:), 'o-', 'Color', cols(p,:), 'LineWidth', 1.2);
        end
        yline(ax1, 0, ':', 'Color', [0.5 0.5 0.5]);
        ylabel(ax1, 'added spikes of pulse n');  ylabel(ax2, 'cumulative added spikes ([n n] levels)');  xlabel(ax2, xlabel_(xAxis));
        title(ax1, sprintf('channel %d: spikes added by pulse n (mean +/- SEM over trials)', ch));
        legend(ax1, hl, lg, 'Location', 'northeast', 'Box', 'off', 'Interpreter','none');
        box(ax1,'off'); box(ax2,'off');
    end
end

%% ================= FIGURE: POPULATION AVERAGE =================
if any(strcmp(plotModes,'average'))
    figure('Color','w','Position',[100 100 1100 850],'Name','population added spikes');
    ax1 = subplot(3,1,1); hold(ax1,'on');  ax2 = subplot(3,1,2); hold(ax2,'on');  ax3 = subplot(3,1,3); hold(ax3,'on');
    lg = {};  hl = [];
    for p = 1:nPat
        r = res{p};  u = r.useCh;
        if isempty(u), warning('PN:nochan','Pattern %d: no channel left to average.', p); continue; end
        for q = 1:numel(r.ser)
            S = r.ser(q);
            A = S.Y(u,:);
            mu = mean(A,1,'omitnan');  se = std(A,0,1,'omitnan')./sqrt(max(sum(~isnan(A),1),1));
            h = errorbar(ax1, S.x, mu, se, S.style, 'Color', cols(p,:), 'LineWidth', 1.3);
            hl(end+1) = h; lg{end+1} = sprintf('cond %d: %s (%d ch)', finalConds(p), S.name, numel(u)); %#ok<SAGROW>
            if strcmp(normalize, 'toPulse1')
                ok = r.R1(u) > 0;
                An = A(ok,:) ./ r.R1(u(ok));
                mn = mean(An,1,'omitnan');  sn = std(An,0,1,'omitnan')./sqrt(max(sum(~isnan(An),1),1));
                errorbar(ax2, S.x, mn, sn, S.style, 'Color', cols(p,:), 'LineWidth', 1.3);
            end
        end
        C = r.cumY(u,:);
        mc = mean(C,1,'omitnan');  sc = std(C,0,1,'omitnan')./sqrt(max(sum(~isnan(C),1),1));
        errorbar(ax3, r.cumX, mc, sc, 'o-', 'Color', cols(p,:), 'LineWidth', 1.3);
    end
    yline(ax1, 0, ':', 'Color', [0.5 0.5 0.5]);  yline(ax2, 0, ':', 'Color', [0.5 0.5 0.5]);
    ylabel(ax1, 'added spikes of pulse n');  title(ax1, 'population: spikes added by pulse n (mean +/- SEM over channels)');
    ylabel(ax2, 'added / response to first step');  title(ax2, sprintf('normalised to each channel''s response to the pattern''s first step (%s)', normalize));
    ylabel(ax3, 'cumulative added spikes');  xlabel(ax3, xlabel_(xAxis));
    title(ax3, 'cumulative, at levels with equal pulses on both electrodes ([n n]; every level if simultaneous)');
    if ~isempty(hl), legend(ax1, hl, lg, 'Location', 'northeast', 'Box', 'off', 'Interpreter','none'); end
    for ax = [ax1 ax2 ax3], box(ax, 'off'); end
end

%% ================= FIGURE: DIAGNOSTICS =================
if showDiagnostics
    figure('Color','w','Position',[150 100 1100 850],'Name','diagnostics');
    d = [];
    for ch = 1:nChn
        ok = find(isfinite(be{ch}));
        d = [d; be{ch}(ok) - arrayfun(@(k) condTable(k).lastPulse_ms, trialCond(ok))]; %#ok<AGROW>
    end
    subplot(3,1,1); histogram(d, 0:0.25:12, 'FaceColor', [0.5 0.5 0.5]); hold on
    xline(winStart_ms, 'r-', 'winStart');  xline(prctile(d,95), 'b--', '95th pct');
    xlabel('recovery window end minus start of the last pulse (ms)'); ylabel('windows');
    title(sprintf('blank end after the last pulse: median %.2f ms, 95th percentile %.2f ms (set winStart_ms above this)', median(d), prctile(d,95)));
    box off
    p = 1;  r = res{p};  cc = r.chain;  t = r.t;  L = numel(cc);  u = r.useCh;
    if isempty(u), u = chUse(:)'; end
    tau = -10 : binMs : 40;
    subplot(3,1,2); hold on
    pre1 = psthRel_(sp, u, t0(trialCond == cc(1) & isfinite(t0)), t(1), tau, binMs);
    plot(tau(1:end-1)+binMs/2, pre1, 'k', 'LineWidth', 1.5);
    xregion_(gca, winStart_ms, winStart_ms + normLen_ms);
    xlabel('ms after pulse 1'); ylabel('spikes/s'); title('level 1: response to the first pulse (shaded = normalisation window)'); box off
    subplot(3,1,3); hold on
    ne = nExplore(nExplore >= 2 & nExplore <= L);
    cm = parula(max(numel(ne),2));
    for i = 1:numel(ne)
        n = ne(i);
        rn  = psthRel_(sp, u, t0(trialCond == cc(n)   & isfinite(t0)), t(n), tau, binMs);
        rn1 = psthRel_(sp, u, t0(trialCond == cc(n-1) & isfinite(t0)), t(n), tau, binMs);
        plot(tau(1:end-1)+binMs/2, rn - rn1, 'Color', cm(i,:), 'LineWidth', 1.2);
    end
    xline(winStart_ms, 'r-');  yline(0, ':');
    xlabel('ms after pulse n'); ylabel('level n - level n-1 (spikes/s)');
    title('what the n-th pulse adds over time (red line = winStart_ms); pattern 1');
    legend(arrayfun(@(n) sprintf('n = %d', n), ne, 'UniformOutput', false), 'Box', 'off'); box off
end

%% ================= HELPERS =================
function s = xlabel_(xAxis)
switch xAxis
    case 'perElectrode', s = 'pulses per electrode';
    case 'time',         s = 'time of the pulse (ms)';
    otherwise,           s = 'level (pulse number)';
end
end

function W = winLen_(n, L, t, g, winLen, lastLen)
if ischar(winLen) || isstring(winLen)
    assert(strcmpi(winLen, 'nextPulse'), 'winLen_ms must be a number, a vector or ''nextPulse''.');
    if n < L, W = t(n+1) - t(n) - g; else, W = lastLen; end
elseif isscalar(winLen)
    W = winLen;
else
    W = winLen(min(n, numel(winLen)));
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

function r = psthRel_(sp, chs, t0s, tPulse, tau, binMs)
% rate (spikes/s) vs time after a pulse at tPulse, averaged over channels chs
cnt = zeros(1, numel(tau)-1);
for ch = chs(:)'
    s = sp{ch}(:,1);
    for i = 1:numel(t0s)
        rel = s(s >= t0s(i)+tPulse+tau(1) & s < t0s(i)+tPulse+tau(end)) - t0s(i) - tPulse;
        cnt = cnt + histcounts(rel, tau);
    end
end
r = cnt / max(numel(t0s)*numel(chs),1) / binMs * 1000;
end

function xregion_(ax, a, b)
yl = ylim(ax);
h = patch(ax, [a b b a], [yl(1) yl(1) yl(2) yl(2)], [1 0.9 0.5], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
uistack(h, 'bottom');
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
