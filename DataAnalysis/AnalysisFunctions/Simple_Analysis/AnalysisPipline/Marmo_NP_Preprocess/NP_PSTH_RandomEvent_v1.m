%% ============================================================
% NP_PSTH_RandomEvent_v1.m
% NP probe | check of the RandomEvent simulation of NP_3_SpikeRecovery (simMode = 'randomEvent')  | v1
% READ-ONLY: only prints and plots, saves nothing.   (NP version of NN_PSTH_RandomEvent_v1)
%
% FILES  spikeFiles = recovered files of the SAME session:
%          1st  = the normal recovery (reference = 100 %), e.g. <short>.NP_Spikes_TrainRecovery.mat
%          then = one or more RandomEvent files, e.g. <short>.NP_Spikes_TrainRecovery_RandomEvent.mat
%        The percent and event_ms of every level are read from each file's own Recovery_params.
%        spikeFiles = {} -> the reference + every *.NP_Spikes_TrainRecovery*RandomEvent*.mat in recording1.
%
% PART 1  Exact check (all channels; always runs)
%   For every thinned level L (percent < 100): spikes in [t_L, t_L + event_ms) after its last pulse t_L, over all
%   trials of the level, in the reference (N) and in the RandomEvent file (kept). The simulation keeps exactly
%   round(percent/100 x N) per channel -> prints the % kept pooled over channels, how many channels match exactly,
%   and the range over channels.
% PART 2  Per-pulse check (conds; channels)   one figure per condition
%   For every pulse j of the condition's train: mean spikes per trial and channel in the window after pulse j
%   (pulseWin_ms = 'event': [0, event_ms of level j]; or a fixed [a b] ms after every pulse), for every file, and the
%   RATIO file / reference in %. Pulse j belongs to level j of the chain (its data were thinned at level j and copied
%   upwards), so the ratio should follow the percent set for level j (black line). Level 1 is never thinned (100 %).
%   Raw counts (no baseline subtraction): the simulation thins ALL spikes in the window.
% PART 3  PSTH (conds)
%   psthMode = 'average'    : one figure per condition, PSTH of every file averaged over the channels (shaded = SEM
%                             across channels) + (showDiff) each file minus the reference
%   psthMode = 'perChannel' : one figure per condition and plotChannels channel: a raster block per file + PSTHs
% conds = [] -> prints the condition list (with the percent / event_ms of each file) and runs Part 1 only.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir  = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
spikeFiles  = {'/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery.mat',...
    '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29/Record Node 101/experiment1/recording1/CJ287_15Pulses3_100Hz.NP_Spikes_TrainRecovery_randomEvent.mat'};         % {} = automatic (see header); or names in recording1 / full paths, REFERENCE FIRST
fileLabels  = {'v1 full recovery', 'RandomEvent'};         % legend names, one per file; {} = 'reference' + the file suffixes
conds       = [15,45];         % conditions for Part 2 and 3 (e.g. the final levels: [15 45]); [] = list + Part 1 only

channels    = [];         % channels averaged in Part 2 and 3; [] = all ACTIVE channels (rate >= autoRateMin)
autoRateMin = 20;         % spikes/s (reference file, whole window): a channel counts as active above this

pulseWin_ms = 'event';    % Part 2 window after each pulse: 'event' = [0, event_ms of that level]; or [a b] in ms

psthMode     = 'perChannel'; % Part 3: 'average' | 'perChannel'
plotChannels = [];        % 'perChannel' only: channels to draw; [] = 4 active channels spread along the probe
showDiff     = true;      % extra panel = file minus reference ('average': below the PSTH; 'perChannel': below raster + PSTH)
winMs        = [-20 200]; % ms around the trigger
binMs        = 1;         % PSTH bin (ms) 
smoothSigmaMs = 3;        % PSTH Gaussian smoothing SD (ms); 0 = none
minDataFrac  = 0.5;       % a PSTH bin is drawn only if at least this fraction of the trials had data in it
dotSize      = 2;         % raster dot size ('perChannel')
%% =================================================

%% ---- files -------------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
if isempty(spikeFiles)
    ref = fullfile(recDir, [short '.NP_Spikes_TrainRecovery.mat']);
    d = dir(fullfile(recDir, [short '.NP_Spikes_TrainRecovery_*RandomEvent*.mat']));
    d = d(~startsWith({d.name}, '._') & ~contains({d.name}, '_backup_') & ~contains({d.name}, '_waveforms'));
    spikeFiles = [{ref}, fullfile({d.folder}, {d.name})];
end
for i = 1:numel(spikeFiles)
    if ~isfile(spikeFiles{i}), spikeFiles{i} = fullfile(recDir, spikeFiles{i}); end
    assert(isfile(spikeFiles{i}), 'File not found: %s', spikeFiles{i});
end
nF = numel(spikeFiles);
assert(nF >= 2, 'Need the reference and at least one RandomEvent file (found %d file).', nF);
if isempty(fileLabels)
    fileLabels = cell(1, nF);  fileLabels{1} = 'reference';
    for i = 2:nF
        [~, n_] = fileparts(spikeFiles{i});  t_ = regexp(n_, 'TrainRecovery_(.*)$', 'tokens', 'once');
        if isempty(t_), fileLabels{i} = n_; else, fileLabels{i} = t_{1}; end
    end
end
assert(numel(fileLabels) == nF, 'fileLabels must have one name per file.');
F = cell(nF, 1);
for i = 1:nF
    F{i} = load(spikeFiles{i}, 'sp_clipped', 'sp_trial', 'blank_ms', 'nSpk', 'goodCh', 'FS', 'preMs', 'postMs', ...
                'nTrialsProcessed', 'Recovery_params', 'condTable', 'depth');
end
fs = F{1}.FS;  nCh = numel(F{1}.sp_clipped);  nTr = F{1}.nTrialsProcessed;
for i = 2:nF
    assert(F{i}.nTrialsProcessed == nTr && numel(F{i}.sp_clipped) == nCh, 'File %d does not match the reference.', i);
end
if ~strcmp(F{1}.Recovery_params.simMode, 'none'), warning('NPRE:ref', 'The FIRST file is not a normal recovery (simMode = %s).', F{1}.Recovery_params.simMode); end
tr = TT.trial;  trialCond = tr.cond(1:nTr);  trigMs = (double(tr.trigAP(1:nTr)) - 1)/fs*1e3;
CT = F{1}.condTable;  depth = F{1}.depth;  nCond = numel(CT);
pct = nan(nCond, nF);  evm = nan(nCond, nF);
for i = 1:nF
    RP = F{i}.Recovery_params;  pct(:, i) = RP.percent(:);  evm(:, i) = RP.event_ms(:);
end
pct(depth == 0, :) = 100;

fprintf('\n============================================================\nNP RANDOM-EVENT CHECK: %s\n', short);
for i = 1:nF, fprintf('  file %d (%s): %s\n', i, fileLabels{i}, spikeFiles{i}); end

%% ---- condition list ---------------------------------------------------------------------------------------
if isempty(conds)
    fprintf('\n  cond  kind   level  depth');
    for i = 2:nF, fprintf('  | %s: percent event(ms)', fileLabels{i}); end
    fprintf('  label\n');
    for c = 1:nCond
        if ~CT(c).isElec, continue; end
        fprintf('  %4d  %-5s  %5d  %5d', c, CT(c).kind, CT(c).level, depth(c));
        for i = 2:nF, fprintf('  | %7.1f %8.2f', pct(c, i), evm(c, i)); end
        fprintf('  %s\n', CT(c).label);
    end
end

%% ---- spikes per trial, every file, every channel ---------------------------------------------------------------
fprintf('\nsorting spikes by trial ...\n');
B = cell(nF, nCh);
for i = 1:nF
    for ch = 1:nCh, B{i, ch} = byTrial_(F{i}.sp_clipped{ch}, F{i}.sp_trial{ch}, trigMs, nTr); end
end

%% ---- PART 1: exact check ------------------------------------------------------------------------------------------
fprintf('\nPART 1: thinning windows [last pulse, last pulse + event_ms), all %d channels\n', nCh);
for i = 2:nF
    fprintf('\n file %d (%s)\n  cond  level  percent  event(ms) | spikes ref    kept   kept %%  | channels exact   min %%   max %%\n', i, fileLabels{i});
    th = find(depth > 0 & pct(:, i) < 100)';
    if isempty(th), fprintf('  no thinned level in this file.\n'); continue; end
    nOK_all = 0;  nCh_all = 0;
    for c = th
        ks = find(trialCond == c)';
        nRef = zeros(nCh, 1);  nKeep = zeros(nCh, 1);
        for ch = 1:nCh
            for t = ks
                w = tr.lastPulseMs(t) + [0 evm(c, i)];
                nRef(ch)  = nRef(ch)  + cnt_(B{1, ch}{t}, w);
                nKeep(ch) = nKeep(ch) + cnt_(B{i, ch}{t}, w);
            end
        end
        exact = nKeep == round(pct(c, i)/100*nRef);
        r = 100*nKeep(nRef > 0)./nRef(nRef > 0);
        fprintf('  %4d  %5d  %7.1f  %9.2f | %10d  %6d  %6.1f  | %6d of %3d   %6.1f  %6.1f\n', c, CT(c).level, pct(c, i), evm(c, i), ...
                sum(nRef), sum(nKeep), 100*sum(nKeep)/max(sum(nRef), 1), sum(exact), nCh, min([r; NaN]), max([r; NaN]));
        nOK_all = nOK_all + sum(exact);  nCh_all = nCh_all + nCh;
    end
    if nOK_all == nCh_all, fprintf('  -> PASS: every channel of every thinned level kept exactly round(percent x N)\n');
    else, fprintf('  -> %d of %d channel x level pairs differ from round(percent x N)\n', nCh_all - nOK_all, nCh_all); end
end
if isempty(conds), return; end
assert(all(conds >= 1 & conds <= nCond), 'conds must be between 1 and %d.', nCond);

%% ---- channels ------------------------------------------------------------------------------------------------------
y = TT.chan.y_um(:);
rate = sum(F{1}.nSpk(:, 1:nTr), 2)/nTr/((F{1}.preMs + F{1}.postMs)/1e3);
o = TT.chan.orderTipToTop(:);  o = o(F{1}.goodCh(o));
act = o(rate(o) >= autoRateMin);
if isempty(channels), channels = act; end
channels = channels(:)';
fprintf('\n%d channels used for the averages (%s)\n', numel(channels), ternary_(isequal(channels(:), act(:)), ...
        sprintf('all active, rate >= %g /s', autoRateMin), 'typed list'));
col = [0 0 0; lines(nF - 1)];

%% ---- PART 2: per-pulse check ---------------------------------------------------------------------------------------
for c = conds(:)'
    ks = find(trialCond == c)';
    if isempty(ks) || isempty(tr.pulseMs{ks(1)}), fprintf('cond %d has no pulses: skipped in Part 2.\n', c); continue; end
    pl = tr.pulseMs{ks(1)};  nP = numel(pl);
    lvC = nan(1, nP);                                          % condition at level j of the same pattern
    for j = 1:nP
        q = find(strcmp({CT.kind}, CT(c).kind) & [CT.level] == j, 1);
        if ~isempty(q), lvC(j) = q; end
    end
    M = nan(nP, nF);  expP = nan(nP, nF);
    for i = 1:nF
        for j = 1:nP
            if ischar(pulseWin_ms)
                if isnan(lvC(j)), continue; end
                w = pl(j) + [0 evm(lvC(j), max(i, 2))];       % event window of level j (from the RandomEvent file)
            else
                w = pl(j) + pulseWin_ms;
            end
            s = 0;
            for ch = channels, for t = ks, s = s + cnt_(B{i, ch}{t}, w); end, end
            M(j, i) = s/(numel(ks)*numel(channels));
            if ~isnan(lvC(j)), expP(j, i) = pct(lvC(j), i); end
        end
    end
    R = 100*M./M(:, 1);
    fprintf('\nPART 2: cond %d %s, %d trials, %d channels  (window %s)\n   pulse  level-cond |', c, CT(c).label, numel(ks), numel(channels), ...
            ternary_(ischar(pulseWin_ms), '[0, event_ms of the level]', mat2str(pulseWin_ms)));
    for i = 1:nF, fprintf(' %12s', fileLabels{i}); end
    for i = 2:nF, fprintf(' | %s: ratio%% (set%%)', fileLabels{i}); end
    fprintf('\n');
    for j = 1:nP
        fprintf('   %5d  %10s |', j, n2s_(lvC(j)));
        fprintf(' %12.3f', M(j, :));
        for i = 2:nF, fprintf(' | %8.1f (%5.1f)', R(j, i), expP(j, i)); end
        fprintf('\n');
    end
    fh = figure('Name', sprintf('%s per-pulse cond %d', short, c), 'Color', 'w', 'Position', [60 60 900 650]);
    ax1 = subplot(2, 1, 1);  hold(ax1, 'on');
    for i = 1:nF, plot(ax1, 1:nP, M(:, i), '-o', 'Color', col(i,:), 'MarkerFaceColor', col(i,:), 'MarkerSize', 4); end
    ylabel(ax1, 'spikes per trial and channel');  grid(ax1, 'on');  xlim(ax1, [0.5 nP + 0.5]);
    legend(ax1, fileLabels, 'Location', 'best', 'Interpreter', 'none');
    title(ax1, sprintf('%s: counts in the window after each pulse', CT(c).label), 'Interpreter', 'none');
    ax2 = subplot(2, 1, 2);  hold(ax2, 'on');
    for i = 2:nF
        stairs(ax2, (1:nP+1) - 0.5, [expP(:, i); expP(end, i)], 'k-', 'LineWidth', 1.5);
        plot(ax2, 1:nP, R(:, i), '-o', 'Color', col(i,:), 'MarkerFaceColor', col(i,:), 'MarkerSize', 4);
    end
    yline(ax2, 100, ':');  grid(ax2, 'on');  xlim(ax2, [0.5 nP + 0.5]);
    xlabel(ax2, 'pulse number (= level of the chain)');  ylabel(ax2, 'file / reference (%)');
    title(ax2, 'ratio to the reference; black = percent set for that level', 'FontSize', 9);
    sgtitle(fh, sprintf('%s  cond %d: per-pulse RandomEvent check (%d channels)', strrep(short, '_', '\_'), c, numel(channels)));
end

%% ---- PART 3: PSTH ----------------------------------------------------------------------------------------------------------
edges = winMs(1):binMs:winMs(2);  ctr = edges(1:end-1) + binMs/2;
if smoothSigmaMs > 0
    sb = smoothSigmaMs/binMs;  kx = -ceil(3*sb):ceil(3*sb);  ker = exp(-0.5*(kx/sb).^2);  ker = ker/sum(ker);
else
    ker = 1;
end
if strcmp(psthMode, 'average')
    for c = conds(:)'
        ks = find(trialCond == c)';  if isempty(ks), continue; end
        mu = nan(nF, numel(ctr));  se = mu;  P1 = [];
        for i = 1:nF
            P = nan(numel(channels), numel(ctr));
            for q = 1:numel(channels)
                ch = channels(q);  P(q, :) = psth_(B{i, ch}, F{i}.blank_ms(ch, :), ks, edges, ker, minDataFrac);
            end
            mu(i, :) = mean(P, 1, 'omitnan');  se(i, :) = std(P, 0, 1, 'omitnan')/sqrt(size(P, 1));
            if i == 1, P1 = P; end
            DF{i} = P - P1; %#ok<SAGROW>
        end
        fh = figure('Name', sprintf('%s PSTH cond %d', short, c), 'Color', 'w', 'Position', [80 80 1000 300 + 250*showDiff]);
        ax = subplot(1 + showDiff, 1, 1);  hold(ax, 'on');  h = gobjects(nF, 1);
        for i = 1:nF
            fill(ax, [ctr fliplr(ctr)], [mu(i,:) - se(i,:), fliplr(mu(i,:) + se(i,:))], col(i,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            h(i) = plot(ax, ctr, mu(i, :), 'Color', col(i,:), 'LineWidth', 1.3);
        end
        pulseTicks_(ax, tr.pulseMs{ks(1)});
        xline(ax, 0, ':');  xlim(ax, winMs);  grid(ax, 'on');  ylabel(ax, 'spikes/s (mean over channels)');
        legend(ax, h, fileLabels, 'Interpreter', 'none', 'Location', 'northeast');
        title(ax, sprintf('%s (n = %d trials, %d channels; shaded = SEM over channels)', CT(c).label, numel(ks), numel(channels)), 'Interpreter', 'none');
        if showDiff
            ax = subplot(2, 1, 2);  hold(ax, 'on');
            for i = 2:nF, plot(ax, ctr, mean(DF{i}, 1, 'omitnan'), 'Color', col(i,:), 'LineWidth', 1.3); end
            pulseTicks_(ax, tr.pulseMs{ks(1)});
            yline(ax, 0, ':');  xline(ax, 0, ':');  xlim(ax, winMs);  grid(ax, 'on');
            xlabel(ax, 'ms after the trigger');  ylabel(ax, 'file - reference (spikes/s)');
            legend(ax, fileLabels(2:end), 'Interpreter', 'none', 'Location', 'northeast');
        end
        sgtitle(fh, sprintf('%s  cond %d: PSTH of every file', strrep(short, '_', '\_'), c));
    end
else
    if isempty(plotChannels), plotChannels = act(unique(round(linspace(1, numel(act), min(4, numel(act)))))); end
    for c = conds(:)'
        ks = find(trialCond == c)';  nK = numel(ks);  if nK == 0, continue; end
        for ch = plotChannels(:)'
            nRow = 3 + showDiff;                                        % raster (2 rows) + PSTH (+ difference)
            fh = figure('Name', sprintf('%s cond %d ch %d', short, c, ch), 'Color', 'w', 'Position', [100 100 800 650 + 180*showDiff]);
            axR = subplot(nRow, 1, [1 2]);  hold(axR, 'on');
            for i = 1:nF
                off = (nF - i)*(nK + 2);  xs = [];  ys = [];
                for r = 1:nK
                    s_ = B{i, ch}{ks(r)};  s_ = s_(s_ >= winMs(1) & s_ < winMs(2));
                    xs = [xs; s_(:)];  ys = [ys; (off + r)*ones(numel(s_), 1)]; %#ok<AGROW>
                end
                plot(axR, xs, ys, '.', 'Color', col(i,:), 'MarkerSize', dotSize + 4);
                text(axR, winMs(2), off + nK/2, [' ' fileLabels{i}], 'Interpreter', 'none', 'FontSize', 8, 'Color', col(i,:));
                if i < nF, yline(axR, off - 0.5, '-', 'Color', [0.7 0.7 0.7]); end
            end
            pulseTicks_(axR, tr.pulseMs{ks(1)});
            xline(axR, 0, ':');  xlim(axR, winMs);  ylim(axR, [0 nF*(nK + 2)]);  set(axR, 'XTickLabel', [], 'YTick', []);  box(axR, 'on');
            title(axR, sprintf('%s  ch %d (y %d um): one raster block per file', CT(c).label, ch, y(ch)), 'Interpreter', 'none');
            axP = subplot(nRow, 1, 3);  hold(axP, 'on');
            Pc = nan(nF, numel(ctr));
            for i = 1:nF
                Pc(i, :) = psth_(B{i, ch}, F{i}.blank_ms(ch, :), ks, edges, ker, minDataFrac);
                plot(axP, ctr, Pc(i, :), 'Color', col(i,:), 'LineWidth', 1.2);
            end
            pulseTicks_(axP, tr.pulseMs{ks(1)});
            xline(axP, 0, ':');  xlim(axP, winMs);  grid(axP, 'on');  ylabel(axP, 'spikes/s');
            legend(axP, fileLabels, 'Interpreter', 'none', 'Location', 'northeast');
            if showDiff                                                 % each file minus the reference
                set(axP, 'XTickLabel', []);
                axD = subplot(nRow, 1, 4);  hold(axD, 'on');
                for i = 2:nF, plot(axD, ctr, Pc(i, :) - Pc(1, :), 'Color', col(i,:), 'LineWidth', 1.2); end
                pulseTicks_(axD, tr.pulseMs{ks(1)});
                yline(axD, 0, ':');  xline(axD, 0, ':');  xlim(axD, winMs);  grid(axD, 'on');
                ylabel(axD, 'file - reference (spikes/s)');  xlabel(axD, 'ms after the trigger');
                legend(axD, fileLabels(2:end), 'Interpreter', 'none', 'Location', 'northeast');
            else
                xlabel(axP, 'ms after the trigger');
            end
            sgtitle(fh, strrep(short, '_', '\_'));
        end
    end
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function Bt = byTrial_(tAbs, k, trigMs, nTr)
Bt = cell(nTr, 1);  for t = 1:nTr, Bt{t} = zeros(0, 1); end
if isempty(k), return; end
k = k(:);  rel = tAbs(:) - trigMs(k);
[ks, o] = sort(k);  rel = rel(o);
e = [0; find(diff(ks)); numel(ks)];
for i = 1:numel(e) - 1
    rows = e(i)+1:e(i+1);  Bt{ks(rows(1))} = sort(rel(rows));
end
end

function n = cnt_(t, w)
% spikes in [w1, w2); the edges are widened by 1e-6 ms so that rounding of the times does not move a spike across
n = sum(t >= w(1) - 1e-6 & t < w(2) - 1e-6);
end

function rt = psth_(Bt, blankRow, ks, edges, ker, minFrac)
binMs = edges(2) - edges(1);  nB = numel(edges) - 1;
cnt = zeros(1, nB);  expo = zeros(1, nB);
for t = ks
    cnt = cnt + histcounts(Bt{t}, edges);
    covr = ones(1, nB);  iv = blankRow{t};
    for b = 1:size(iv, 1)
        lo = max(edges(1:end-1), iv(b,1));  hi = min(edges(2:end), iv(b,2));
        covr = covr - max(0, hi - lo)/binMs;
    end
    expo = expo + max(covr, 0);
end
cS = conv(cnt, ker, 'same');  eS = conv(expo, ker, 'same');
rt = cS ./ (eS*binMs/1e3);
rt(expo < minFrac*max(numel(ks), 1) | eS <= 0) = NaN;
end

function pulseTicks_(ax, pt)
if isempty(pt), return; end
yl = ylim(ax);
for p = pt(:)', xline(ax, p, '-', 'Color', [1 0.6 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off'); end
ylim(ax, yl);
end

function s = n2s_(v)
if isnan(v), s = '-'; else, s = sprintf('%d', v); end
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
