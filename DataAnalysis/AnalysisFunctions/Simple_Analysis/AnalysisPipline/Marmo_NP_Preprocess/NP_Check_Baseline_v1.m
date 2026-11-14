%% ============================================================
% NP_Check_Baseline_v1.m
% NP probe | after NP_2 / NP_3 | CHECK OF THE PRE-TRIGGER BASELINE by trial type, depth band and time | v1
% READ-ONLY: figures and printout only, nothing is saved.
%
% Trial types: sim, seq (A->B / B->A), blank, visual (from TT.cond.kind).
% Channels: processed (goodCh) and ACTIVE (rate over the whole window >= autoRateMin), split into depth bands.
% Baseline of a trial = spikes in baseWinMs (relative to the TRIGGER) / its length, mean over the band's channels.
%
% PART A  pre-trigger rate (preWinMs, binMs bins) per trial type, one panel per depth band (+ all active channels)
%         -> visual higher EVENLY over the window = chance / state;  rising towards 0 = something happens before the trial
% PART B  setup before the trial: photodiode (NI AI0) averaged over diodeWinMs per trial type (a screen change before
%         the trigger would show here), and the time since the previous trigger per trial type
% PART C  distribution of the per-trial baseline (log10 scale) per trial type and depth band  (two states -> two humps)
% PART D  baseline of every trial in recording order, coloured by trial type, per depth band (runs of active trials?)
% PART E  response vs baseline, trial by trial, for respConds: is the response different in active-state trials?
%         response window: train = first pulse ... last pulse + tail_ms; visual = visWinMs after the photodiode onset
% PRINTED: median baseline per type and band, share of 'active-state' trials (baseline > activeFactor x the band's
%          median over all trials), permutation test of visual vs electrical (difference of medians), Spearman r of E.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir  = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFile   = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery_randomEvent.mat';            % '' = <short>.NP_Spikes_TrainRecovery.mat if it exists (needed for Part E), else <short>.NP_Spikes.mat
bands = { ...                % depth bands from the probe TIP (um)
    '0-2200',     [   0 2200]; ...
    '2200-3600',  [2200 3600]; ...
    '3600-4500',  [3600 4500]; ...
    '4500-5800',  [4500 5800]; ...
    };
autoRateMin  = 20;           % spikes/s: active channels
baseWinMs    = [-100 -5];    % baseline window (ms, relative to the trigger)
activeFactor = 2;            % 'active-state' trial: baseline > activeFactor x median baseline of the band (all trials)

doA = true;  preWinMs = [-300 0];  binMs = 10;
doB = true;  diodeWinMs = [-1000 300];
doC = true;
doD = true;
doE = true;  respConds = [];  % conditions for Part E; [] = final sim + final seq + the visual conditions with contrast 1
tail_ms = 15;  visWinMs = [0 200];
nPerm = 5000;                 % permutations for the visual-vs-electrical test
%% =================================================

%% ---- load -------------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
if isempty(spikeFile)
    spikeFile = fullfile(recDir, [short '.NP_Spikes_TrainRecovery.mat']);
    if ~isfile(spikeFile), spikeFile = fullfile(recDir, [short '.NP_Spikes.mat']); end
end
if ~isfile(spikeFile), spikeFile = fullfile(recDir, spikeFile); end
SP = load(spikeFile, 'sp_clipped', 'sp_trial', 'nSpk', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed');
fs = SP.FS;  nTr = SP.nTrialsProcessed;  tr = TT.trial;
trialCond = tr.cond(1:nTr);  trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;
assert(min([preWinMs baseWinMs]) >= -SP.preMs, 'A window starts before the detection window (-%g ms).', SP.preMs);
kindC = {TT.cond.kind};
typ = cell(nTr, 1);
for t = 1:nTr
    k = kindC{trialCond(t)};
    if any(strcmp(k, {'AB', 'BA'})), typ{t} = 'seq'; else, typ{t} = k; end
end
types = intersect({'sim', 'seq', 'blank', 'visual'}, unique(typ), 'stable');  types = types(:)';   % always a row
tcol = containers.Map({'sim', 'seq', 'blank', 'visual'}, {[0 0.45 0.74], [0.47 0.67 0.19], [0.4 0.4 0.4], [0.49 0.18 0.56]});
y = TT.chan.y_um(:);
rateAll = sum(SP.nSpk(:, 1:nTr), 2)/nTr/((SP.preMs + SP.postMs)/1e3);
act = SP.goodCh(:) & rateAll >= autoRateMin;
nB = size(bands, 1);  bandCh = cell(nB + 1, 1);  bandName = [bands(:, 1)', {'all active'}];
for b = 1:nB, r = bands{b, 2}; bandCh{b} = find(act & y >= r(1) & y < r(2)); end
bandCh{nB + 1} = find(act);
fprintf('\n============================================================\nNP BASELINE CHECK: %s\n', short);
fprintf('spike file: %s\nbaseline %s ms (trigger) | active channels: %d\n', spikeFile, mat2str(baseWinMs), sum(act));
for b = 1:nB + 1, fprintf('  band %-11s %3d channels\n', bandName{b}, numel(bandCh{b})); end

%% ---- spikes relative to the trigger, per trial (active channels only) ---------------------------------------
Rel = cell(numel(SP.sp_clipped), 1);
for ch = find(act)'
    k = SP.sp_trial{ch}(:);  Rel{ch} = [k, SP.sp_clipped{ch}(:) - trigMs(k)];
end
cntIn = @(ch, w) accumarray(Rel{ch}(Rel{ch}(:,2) >= w(1) & Rel{ch}(:,2) < w(2), 1), 1, [nTr 1]);
Bt = nan(nTr, nB + 1);                                   % baseline rate per trial and band
for b = 1:nB + 1
    if isempty(bandCh{b}), continue; end
    M = zeros(nTr, numel(bandCh{b}));
    for q = 1:numel(bandCh{b}), M(:, q) = cntIn(bandCh{b}(q), baseWinMs); end
    Bt(:, b) = mean(M, 2)/(diff(baseWinMs)/1e3);
end

%% ---- printout: medians, active share, permutation test --------------------------------------------------------
fprintf('\nmedian baseline (spikes/s) | share of active-state trials (> %g x band median)\n  %-11s', activeFactor, 'band');
for i = 1:numel(types), fprintf(' %16s', types{i}); end
fprintf('   visual vs electrical: diff of medians, p (perm.)\n');
isV = strcmp(typ, 'visual');  isE = strcmp(typ, 'sim') | strcmp(typ, 'seq');
rng(1);
for b = 1:nB + 1
    v = Bt(:, b);  if all(isnan(v)), continue; end
    thrA = activeFactor*median(v, 'omitnan');
    fprintf('  %-11s', bandName{b});
    for i = 1:numel(types), m = strcmp(typ, types{i}); fprintf('   %6.1f (%3.0f%%)', median(v(m)), 100*mean(v(m) > thrA)); end
    if any(isV) && any(isE)
        d0 = median(v(isV)) - median(v(isE));  pool = v(isV | isE);  nVv = sum(isV);  cnt = 0;
        for p = 1:nPerm
            pp = pool(randperm(numel(pool)));  cnt = cnt + (abs(median(pp(1:nVv)) - median(pp(nVv+1:end))) >= abs(d0));
        end
        fprintf('   %+7.1f /s, p = %.4f', d0, (cnt + 1)/(nPerm + 1));
    end
    fprintf('\n');
end

%% ---- PART A: pre-trigger rate per trial type --------------------------------------------------------------------
if doA
    ed = preWinMs(1):binMs:preWinMs(2);  ct = ed(1:end-1) + binMs/2;
    fh = figure('Name', [short ' A pre-trigger rate'], 'Color', 'w', 'Position', [40 40 1500 380]);
    for b = 1:nB + 1
        ax = subplot(1, nB + 1, b);  hold(ax, 'on');
        for i = 1:numel(types)
            m = find(strcmp(typ, types{i}));  R = zeros(1, numel(ct));
            for ch = bandCh{b}(:)'
                s = Rel{ch};  s = s(ismember(s(:,1), m), 2);
                R = R + histcounts(s, ed);
            end
            R = R/(numel(m)*max(numel(bandCh{b}), 1)*binMs/1e3);
            plot(ax, ct, R, '-', 'Color', tcol(types{i}), 'LineWidth', 1.5);
        end
        xline(ax, baseWinMs(1), ':');  xline(ax, baseWinMs(2), ':');  grid(ax, 'on');  xlim(ax, preWinMs);
        title(ax, sprintf('%s (%d ch)', bandName{b}, numel(bandCh{b})), 'FontSize', 9);  xlabel(ax, 'ms before the trigger');
        if b == 1, ylabel(ax, 'spikes/s per channel'); legend(ax, types, 'Location', 'best'); end
    end
    sgtitle(fh, sprintf('%s  A: pre-trigger rate per trial type (dotted = baseline window)', strrep(short, '_', '\_')));
end

%% ---- PART B: photodiode and time since the previous trigger -------------------------------------------------------
if doB
    I = TT.info;  niFile = I.niFile;  if ~isfile(niFile), niFile = [recDir strrep(I.niFile, I.recDir, '')]; end
    oe = jsondecode(fileread(fullfile(recDir, 'structure.oebin')));
    cont = oe.continuous;  if isstruct(cont), cont = num2cell(cont); end
    iNI = find(cellfun(@(c) contains(c.folder_name, 'NI-DAQ'), cont), 1);
    nNI = double(cont{iNI}.num_channels);  chI = cont{iNI}.channels;  if isstruct(chI), chI = num2cell(chI); end
    dAI = I.settings.diodeAI;  niV = double(chI{dAI + 1}.bit_volts);  fsNI = I.fsNI;
    mmNI = memmapfile(niFile, 'Format', {'int16', [nNI I.nSampNI], 'x'});
    o1 = round(diodeWinMs(1)*fsNI/1e3);  o2 = round(diodeWinMs(2)*fsNI/1e3);  tD = (o1:o2)/fsNI*1e3;
    Dm = zeros(numel(types), numel(tD));  Dpre = nan(nTr, 1);
    for i = 1:numel(types)
        m = find(strcmp(typ, types{i}))';  acc = zeros(1, numel(tD));  n = 0;
        for t = m
            a = tr.trigNI(t) + o1;  bb = tr.trigNI(t) + o2;
            if a < 1 || bb > I.nSampNI, continue; end
            v = double(mmNI.Data.x(dAI + 1, a:bb))*niV;  acc = acc + v;  n = n + 1;
            Dpre(t) = mean(v(tD < 0));
        end
        Dm(i, :) = acc/max(n, 1);
    end
    gap = [NaN; diff(double(tr.trigAP(1:nTr)))/fs];        % s since the previous trigger
    fh = figure('Name', [short ' B setup'], 'Color', 'w', 'Position', [60 60 1300 420]);
    ax = subplot(1, 3, 1);  hold(ax, 'on');
    for i = 1:numel(types), plot(ax, tD, Dm(i, :), 'Color', tcol(types{i}), 'LineWidth', 1.2); end
    xline(ax, 0, ':');  grid(ax, 'on');  xlabel(ax, 'ms after the trigger');  ylabel(ax, 'photodiode AI0 (V), mean');
    legend(ax, types, 'Location', 'best');  title(ax, 'photodiode, mean per trial type', 'FontSize', 9);
    ax = subplot(1, 3, 2);  hold(ax, 'on');
    for i = 1:numel(types), plot(ax, tD(tD < 0), Dm(i, tD < 0), 'Color', tcol(types{i}), 'LineWidth', 1.2); end
    grid(ax, 'on');  xlabel(ax, 'ms after the trigger');  title(ax, 'same, before the trigger only (zoomed)', 'FontSize', 9);
    ax = subplot(1, 3, 3);  hold(ax, 'on');
    for i = 1:numel(types)
        m = strcmp(typ, types{i});  g = gap(m);  g = g(~isnan(g));
        plot(ax, i + 0.15*randn(numel(g), 1), g, '.', 'Color', tcol(types{i}));
        plot(ax, i + [-0.3 0.3], median(g)*[1 1], 'k-', 'LineWidth', 2);
    end
    set(ax, 'XTick', 1:numel(types), 'XTickLabel', types);  grid(ax, 'on');  ylabel(ax, 's since the previous trigger');
    title(ax, 'time since the previous trigger (black = median)', 'FontSize', 9);
    sgtitle(fh, sprintf('%s  B: what the setup does before the trial', strrep(short, '_', '\_')));
    fprintf('\nPART B: mean photodiode before the trigger (V) and median time since the previous trigger (s):\n');
    for i = 1:numel(types)
        m = strcmp(typ, types{i});
        fprintf('  %-7s photodiode %.4f V   gap %.2f s\n', types{i}, mean(Dpre(m), 'omitnan'), median(gap(m), 'omitnan'));
    end
end

%% ---- PART C: baseline distribution --------------------------------------------------------------------------------
if doC
    fh = figure('Name', [short ' C baseline distribution'], 'Color', 'w', 'Position', [80 80 1500 380]);
    ed = linspace(-0.5, 3, 36);                                   % log10(spikes/s + 1)
    for b = 1:nB + 1
        ax = subplot(1, nB + 1, b);  hold(ax, 'on');
        for i = 1:numel(types)
            v = log10(Bt(strcmp(typ, types{i}), b) + 1);
            stairs(ax, ed(1:end-1), histcounts(v, ed, 'Normalization', 'probability'), 'Color', tcol(types{i}), 'LineWidth', 1.3);
        end
        xline(ax, log10(activeFactor*median(Bt(:, b), 'omitnan') + 1), 'k--');
        grid(ax, 'on');  xlabel(ax, 'log_{10}(baseline + 1)');  title(ax, bandName{b}, 'FontSize', 9);
        if b == 1, ylabel(ax, 'fraction of trials'); legend(ax, types, 'Location', 'best'); end
    end
    sgtitle(fh, sprintf('%s  C: per-trial baseline (dashed = active-state limit, %g x band median)', strrep(short, '_', '\_'), activeFactor));
end

%% ---- PART D: baseline in recording order ---------------------------------------------------------------------------
if doD
    fh = figure('Name', [short ' D recording order'], 'Color', 'w', 'Position', [100 100 1300 180*(nB + 1)]);
    for b = 1:nB + 1
        ax = subplot(nB + 1, 1, b);  hold(ax, 'on');
        for i = 1:numel(types)
            m = find(strcmp(typ, types{i}));
            plot(ax, m, Bt(m, b), '.', 'Color', tcol(types{i}), 'MarkerSize', 8);
        end
        plot(ax, 1:nTr, movmedian(Bt(:, b), 15, 'omitnan'), 'k-', 'LineWidth', 1);
        yline(ax, activeFactor*median(Bt(:, b), 'omitnan'), 'k--');
        grid(ax, 'on');  ylabel(ax, 'spikes/s');  title(ax, bandName{b}, 'FontSize', 9);  xlim(ax, [0 nTr + 1]);
        if b == 1, legend(ax, [types, {'running median (15 trials)', 'active-state limit'}], 'Location', 'eastoutside'); end
    end
    xlabel(ax, 'trial (recording order)');
    sgtitle(fh, sprintf('%s  D: baseline of every trial in recording order', strrep(short, '_', '\_')));
end

%% ---- PART E: response vs baseline -------------------------------------------------------------------------------------
if doE
    if isempty(respConds)
        respConds = [TT.pattern.finalCond];
        for c = find(strcmp(kindC, 'visual')), if abs(TT.cond(c).contrast - 1) < 1e-9, respConds(end+1) = c; end, end %#ok<AGROW>
    end
    nC = numel(respConds);
    fh = figure('Name', [short ' E response vs baseline'], 'Color', 'w', 'Position', [120 40 340*nC 260*(nB + 1)]);
    fprintf('\nPART E: Spearman r between baseline and response (per trial)\n  %-11s', 'band');
    for c = respConds, fprintf('  %22s', TT.cond(c).label); end
    fprintf('\n');
    rs = nan(nB + 1, nC);
    for j = 1:nC
        c = respConds(j);  ks = find(trialCond == c);
        W = nan(numel(ks), 2);
        for q = 1:numel(ks)
            t = ks(q);
            if strcmp(kindC{c}, 'visual')
                on = tr.diodeOnMs(t);  if isnan(on), on = 0; end
                W(q, :) = on + visWinMs;
            else
                W(q, :) = [tr.firstPulseMs(t), tr.lastPulseMs(t) + tail_ms];
            end
        end
        for b = 1:nB + 1
            if isempty(bandCh{b}), continue; end
            Rr = zeros(numel(ks), 1);
            for ch = bandCh{b}(:)'
                s = Rel{ch};
                for q = 1:numel(ks), w = s(s(:,1) == ks(q), 2);  Rr(q) = Rr(q) + sum(w >= W(q,1) & w < W(q,2)); end
            end
            Rr = Rr./(numel(bandCh{b})*diff(W, 1, 2)/1e3);          % spikes/s per channel in the response window
            x = Bt(ks, b);
            rs(b, j) = spearman_(x, Rr);
            ax = subplot(nB + 1, nC, (b - 1)*nC + j);  hold(ax, 'on');
            plot(ax, x, Rr, 'o', 'Color', tcol(typ{ks(1)}), 'MarkerFaceColor', tcol(typ{ks(1)}), 'MarkerSize', 4);
            xline(ax, activeFactor*median(Bt(:, b), 'omitnan'), 'k--');  grid(ax, 'on');
            title(ax, sprintf('%s | %s  (r = %.2f)', TT.cond(c).label, bandName{b}, rs(b, j)), 'FontSize', 8, 'Interpreter', 'none');
            if j == 1, ylabel(ax, 'response (spikes/s)'); end
            if b == nB + 1, xlabel(ax, 'baseline (spikes/s)'); end
        end
    end
    for b = 1:nB + 1, fprintf('  %-11s', bandName{b}); fprintf('  %22.2f', rs(b, :)); fprintf('\n'); end
    sgtitle(fh, sprintf('%s  E: response vs baseline per trial (dashed = active-state limit)', strrep(short, '_', '\_')));
end

%% ---- helpers ------------------------------------------------------------------------------------------------------------
function r = spearman_(x, y)
ok = ~isnan(x) & ~isnan(y);  x = x(ok);  y = y(ok);
if numel(x) < 3, r = NaN; return; end
cc = corrcoef(rank_(x), rank_(y));  r = cc(1, 2);
end

function rk = rank_(v)
% ranks with ties averaged (no Statistics Toolbox needed)
[s, o] = sort(v(:));  rk = zeros(numel(v), 1);  i = 1;
while i <= numel(s)
    j = i;  while j < numel(s) && s(j+1) == s(i), j = j + 1; end
    rk(o(i:j)) = (i + j)/2;  i = j + 1;
end
end
