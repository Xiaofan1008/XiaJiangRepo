%% ============================================================
% NP_SP_Waveforms_v1.m
% NP probe | after NP_2_SpikeDetect | SPIKE WAVEFORM CHECK | v1   (read-only: nothing is saved)
%
% PURPOSE
%   Are the detected events real spikes? For chosen channels it overlays the spike waveforms, ALIGNED ON THEIR
%   NEGATIVE PEAK, split by WHEN the spike happened in the trial:
%       pre        : before the trigger (spontaneous activity = the reference shape)
%       afterBlank : within afterBlankMs after the end of the artifact blank (where an artifact remnant would show)
%       post       : every other spike after the trigger (blank and visual trials: everything after the trigger)
%   plus a histogram of the peak-to-peak amplitude of each group, and a summary of all channels along the probe.
%
% ALIGNMENT
%   The detection already puts the spike time on the most negative point, and the snippet is centred there, so the
%   trough should be the middle sample. The script reports how many waveforms have their minimum elsewhere.
%   realign = true shifts every waveform so that its minimum (searched within +/- realignMs of the centre) sits on
%   the centre sample; the samples shifted in from outside the snippet are left empty (NaN).
%
% INPUTS (in the recording1 folder of the session, found from the short name of the trial table)
%   <short>.NP_TrialTable.mat, <short>.NP_Spikes.mat, <short>.NP_Spikes_waveforms.mat (saveWaveforms = true in step 2)
%
% FIGURES
%   1..n  waveform figures, up to chPerFig channels each: [pre | afterBlank | post | peak-to-peak histogram]
%         grey = nShow random waveforms, black = mean, dashed = mean +/- SD
%   last  all channels along the probe: spike rate, mean waveform peak-to-peak, noise; the plotted channels are marked
%
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
channels   = [];          % AP channels (1-based) to plot; [] = automatic choice (see below)

% automatic channel choice (only used when channels = [])
autoN       = 20;          % channels spread evenly along the ACTIVE part of the probe
autoRateMin = 20;         % spikes/s: a channel counts as active above this rate
autoQuietN  = 2;          % plus this many channels from the quietest part at the TOP of the probe (noise reference)

% waveform groups
afterBlankMs = 10;        % ms after the end of the blank that count as 'afterBlank'
conds        = [];        % only spikes of these conditions (TT.cond numbers); [] = all conditions

% plotting
nShow     = 200;          % random waveforms drawn per panel
chPerFig  = 4;            % channels per figure (one row each)
realign   = true;         % true = put each waveform's minimum on the centre sample (see ALIGNMENT)
realignMs = 0.2;          % ms: the minimum is searched within +/- this of the centre
ptpEdges  = 0:5:300;      % uV: bins of the peak-to-peak histogram
rngSeed   = 1;            % random choice of the drawn waveforms (same seed = same picture)
%% =================================================

rng(rngSeed);

%% ---- files -----------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;
short  = erase(tf(1).name, '.NP_TrialTable.mat');
spFile = fullfile(recDir, [short '.NP_Spikes.mat']);
wfFile = fullfile(recDir, [short '.NP_Spikes_waveforms.mat']);
assert(isfile(spFile), 'Missing %s: run NP_2_SpikeDetect first.', spFile);
assert(isfile(wfFile), 'Missing %s: run NP_2_SpikeDetect with saveWaveforms = true.', wfFile);

L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;
SP = load(spFile, 'sp_clipped', 'sp_trial', 'nSpk', 'blankEnd_ms', 'noise_uV', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed');
fprintf('loading waveforms (%s) ...\n', wfFile);
WF = load(wfFile, 'sp_waveforms', 'wfScale_uV', 'wfTime_ms');
fs = SP.FS;  nCh = numel(SP.sp_clipped);  nTr = SP.nTrialsProcessed;
tWf = WF.wfTime_ms(:)';  nSmp = numel(tWf);  cS = find(abs(tWf) == min(abs(tWf)), 1);   % centre sample

trigMs  = (double(TT.trial.trigAP(1:nTr)) - 1)/fs*1e3;      % trigger of every trial, ms from AP file start
condTr  = TT.trial.cond(1:nTr);
useTr   = true(nTr, 1);  if ~isempty(conds), useTr = ismember(condTr, conds); end
durS    = (SP.preMs + SP.postMs)/1e3;
rate    = sum(SP.nSpk(:, 1:nTr), 2)/nTr/durS;                % spikes/s per channel
y       = TT.chan.y_um(:);

%% ---- channel choice ----------------------------------------------------------------------------------
if isempty(channels)
    o  = TT.chan.orderTipToTop(:);  o = o(SP.goodCh(o));
    act = o(rate(o) >= autoRateMin);
    assert(~isempty(act), 'No channel above autoRateMin = %g spikes/s.', autoRateMin);
    pick = act(unique(round(linspace(1, numel(act), min(autoN, numel(act))))));
    top  = o(y(o) > max(y(act)));                            % above the active part
    if ~isempty(top)
        [~, q] = sort(rate(top));  quiet = top(q(1:min(autoQuietN, numel(top))));
    else
        quiet = [];
    end
    channels = [pick(:); quiet(:)];
end
[~, q] = sort(y(channels), 'descend');  channels = channels(q);   % top of the probe first
fprintf('channels (top -> tip): %s\n', mat2str(channels(:)'));

%% ---- per channel: groups, alignment ------------------------------------------------------------------
grpName = {'pre', 'afterBlank', 'post'};
res = struct('ch', {}, 'W', {}, 'g', {}, 'n', {}, 'offCentre', {}, 'ptp', {});
for c = channels(:)'
    tAbs = SP.sp_clipped{c};  k = SP.sp_trial{c};
    W = double(WF.sp_waveforms{c}) * WF.wfScale_uV;          % spikes x samples, uV
    assert(size(W, 1) == numel(tAbs), 'ch %d: %d waveforms but %d spike times.', c, size(W, 1), numel(tAbs));
    sel = useTr(k);  tAbs = tAbs(sel);  k = k(sel);  W = W(sel, :);
    rel = tAbs - trigMs(k);                                   % ms after the trigger
    bE  = SP.blankEnd_ms(c, k)';  bE = bE(:);
    g = 3*ones(numel(rel), 1);                                % post
    g(rel < 0) = 1;                                           % pre
    g(~isnan(bE) & rel >= bE & rel < bE + afterBlankMs) = 2;  % afterBlank
    [~, mi] = min(W, [], 2);
    offC = mean(mi ~= cS)*100;
    if realign
        r = round(realignMs*fs/1e3);
        win = max(1, cS - r):min(nSmp, cS + r);
        [~, m2] = min(W(:, win), [], 2);  sh = win(m2)' - cS;   % shift needed for each waveform
        W2 = nan(size(W));
        for s = unique(sh)'
            rows = sh == s;
            if s >= 0, W2(rows, 1:nSmp-s) = W(rows, 1+s:nSmp); else, W2(rows, 1-s:nSmp) = W(rows, 1:nSmp+s); end
        end
        W = W2;
    end
    ptp = max(W, [], 2) - min(W, [], 2);
    nG = [sum(g == 1), sum(g == 2), sum(g == 3)];
    res(end+1) = struct('ch', c, 'W', {W}, 'g', {g}, 'n', nG, 'offCentre', offC, 'ptp', {ptp}); %#ok<SAGROW>
end

%% ---- printout --------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP SPIKE WAVEFORMS: %s\n', short);
fprintf('waveform %.2f ms (%d samples), centre sample %d; realign = %d\n', tWf(end) - tWf(1), nSmp, cS, realign);
fprintf('  ch    y(um)  rate(/s)  noise(uV)   n pre   n afterBlank   n post   min not at centre   median p-p (uV)\n');
for q = 1:numel(res)
    c = res(q).ch;
    fprintf(' %3d  %6d   %7.1f    %6.1f   %7d   %9d   %8d       %5.1f %%          %6.1f\n', c, y(c), rate(c), ...
            SP.noise_uV(c), res(q).n(1), res(q).n(2), res(q).n(3), res(q).offCentre, median(res(q).ptp, 'omitnan'));
end

%% ---- waveform figures -------------------------------------------------------------------------------------
col = [0.65 0.65 0.65];
nFig = ceil(numel(res)/chPerFig);
for f = 1:nFig
    rr = (f-1)*chPerFig + 1 : min(numel(res), f*chPerFig);
    figure('Name', sprintf('Waveforms %s (%d/%d)', short, f, nFig), 'Color', 'w', 'Position', [60 60 1300 230*numel(rr)]);
    for i = 1:numel(rr)
        R = res(rr(i));  c = R.ch;
        lim = prctile_(R.W(:), [0.5 99.5]);  lim = lim + [-1 1]*0.1*diff(lim);
        if any(~isfinite(lim)) || diff(lim) <= 0, lim = [-100 50]; end
        for gI = 1:3
            ax = subplot(numel(rr), 4, (i-1)*4 + gI);  hold(ax, 'on');
            Wg = R.W(R.g == gI, :);
            if ~isempty(Wg)
                d = randperm(size(Wg, 1), min(nShow, size(Wg, 1)));
                plot(ax, tWf, Wg(d, :)', 'Color', [col 0.25], 'LineWidth', 0.5);
                mu = mean(Wg, 1, 'omitnan');  sd = std(Wg, 0, 1, 'omitnan');
                plot(ax, tWf, mu, 'k', 'LineWidth', 2);
                plot(ax, tWf, mu + sd, 'k--', tWf, mu - sd, 'k--', 'LineWidth', 0.8);
            end
            xline(ax, 0, ':');  ylim(ax, lim);  xlim(ax, tWf([1 end]));  grid(ax, 'on');
            title(ax, sprintf('ch %d (y %d um) %s: n = %d', c, y(c), grpName{gI}, size(Wg, 1)), 'FontSize', 8);
            if gI == 1, ylabel(ax, 'uV'); end
            if i == numel(rr), xlabel(ax, 'ms from the trough'); end
        end
        ax = subplot(numel(rr), 4, (i-1)*4 + 4);  hold(ax, 'on');
        cl = lines(3);  shown = {};
        for gI = 1:3
            p = R.ptp(R.g == gI);  p = p(~isnan(p));
            if isempty(p), continue; end
            h = histcounts(p, ptpEdges, 'Normalization', 'probability');
            stairs(ax, ptpEdges(1:end-1), h, 'Color', cl(gI,:), 'LineWidth', 1.2);
            shown{end+1} = grpName{gI}; %#ok<SAGROW>
        end
        if ~isempty(shown), legend(ax, shown, 'FontSize', 7, 'Location', 'northeast'); end
        grid(ax, 'on');
        title(ax, sprintf('peak-to-peak; min not at centre %.1f %%', R.offCentre), 'FontSize', 8);
        if i == numel(rr), xlabel(ax, 'peak-to-peak (uV)'); end
    end
    sgtitle(sprintf('%s: spike waveforms aligned on the trough%s', strrep(short, '_', '\_'), ternary_(realign, ' (realigned)', '')));
end

%% ---- summary along the probe -----------------------------------------------------------------------------
mPtp = nan(nCh, 1);
for c = find(SP.goodCh(:))'
    if isempty(WF.sp_waveforms{c}), continue; end
    Wc = double(WF.sp_waveforms{c}) * WF.wfScale_uV;
    if size(Wc, 1) > 2000, Wc = Wc(randperm(size(Wc, 1), 2000), :); end
    mw = mean(Wc, 1);  mPtp(c) = max(mw) - min(mw);
end
figure('Name', ['Probe summary ' short], 'Color', 'w', 'Position', [80 80 1000 700]);
vals = {rate, mPtp, SP.noise_uV(:)};  xl = {'spike rate (/s, whole window)', 'peak-to-peak of the mean waveform (uV)', 'noise (uV)'};
for i = 1:3
    ax = subplot(1, 3, i);  hold(ax, 'on');
    plot(ax, vals{i}, y, '.', 'Color', [0.4 0.4 0.4]);
    plot(ax, vals{i}(channels), y(channels), 'ro', 'MarkerFaceColor', 'r');
    xlabel(ax, xl{i});  grid(ax, 'on');  ylim(ax, [min(y) max(y)]);
    if i == 1, ylabel(ax, 'y from the tip (um)'); end
end
sgtitle(sprintf('%s: all channels along the probe (red = plotted channels)', strrep(short, '_', '\_')));

%% ---- helpers ------------------------------------------------------------------------------------------------
function v = prctile_(x, p)
x = sort(x(~isnan(x(:))));
v = nan(size(p));
if isempty(x), return; end
for i = 1:numel(p)
    r = 1 + (numel(x) - 1)*p(i)/100;  lo = floor(r);  hi = ceil(r);
    v(i) = x(lo) + (r - lo)*(x(hi) - x(lo));
end
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
