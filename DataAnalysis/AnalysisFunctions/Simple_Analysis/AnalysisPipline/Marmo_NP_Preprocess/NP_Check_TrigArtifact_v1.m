%% ============================================================
% NP_Check_TrigArtifact_v1.m
% NP probe | after NP_1_TrialTable | CHECK: trigger (NI AI3) vs stimulation artifact (AP) on one time axis | v1
% READ-ONLY: figures and printout only, nothing is saved.
%
% For every chosen trial one figure, all rows on the SAME time axis (ms after the detected trigger, probe clock):
%   row 1  TRIGGER   : NI AI3 voltage (each NI sample converted to probe time with the sync fit of step 1),
%                      the 2 V detection threshold (trigThrV of step 1) and the detected trigger (t = 0)
%   row 2  RAW AP    : raw signal of the chosen channels (minus their median before the trigger)
%   row 3  FILTERED  : the same channels band-passed (bpHz, zero phase), as used for spike detection - no blanking here
%   row 4  COMMON STEP: median over ALL channels of |sample-to-sample step| of the raw signal (the artifact signal),
%                      its threshold (artK x pre-trigger median + artAbsUV) and the detected artifact onset
%   red dashed lines: the pulse times stored in the trial table (design + measured delay)
%   two time ranges: zoomWinMs (around the first pulse) and, if showTrain = true, fullWinMs (the whole train)
% conds = [] -> prints the condition list and stops.
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir   = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
conds        = [1 2 15 45];          % conditions to look at; [] = print the list and stop
trialsPerCond = 1;          % how many trials per condition (the first ones), one figure each
channels     = [];          % AP channels for rows 2 and 3; [] = 4 channels spread over 3600-5800 um from the tip
zoomWinMs    = [-5 10];      % ms around the trigger (first figure column)
showTrain    = true;        % also show the whole train (second figure column)
fullWinMs    = [-10 160];   % ms around the trigger for the whole-train column
bpHz         = [300 6000];  % band-pass for row 3 (as step 2)
artK         = 15;  artAbsUV = 20;   % artifact threshold, as the artifact check of step 1
%% =================================================

tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
L = load(fullfile(recDir, tf(1).name), 'TT');  TT = L.TT;  clear L
if isempty(conds)
    fprintf('\nConditions of %s:\n  cond  kind     level  trials  label\n', short);
    for c = 1:numel(TT.cond), fprintf('  %4d  %-7s  %5d  %6d  %s\n', c, TT.cond(c).kind, TT.cond(c).level, TT.cond(c).nTrials, TT.cond(c).label); end
    return
end

%% ---- files and streams -------------------------------------------------------------------------------
I = TT.info;  fsAP = I.fsAP;  fsNI = I.fsNI;  nCh = I.nChAP;  uv = I.uVperBit;
apFile = I.apFile;  niFile = I.niFile;
if ~isfile(apFile), apFile = [recDir strrep(I.apFile, I.recDir, '')]; end
if ~isfile(niFile), niFile = [recDir strrep(I.niFile, I.recDir, '')]; end
oe = jsondecode(fileread(fullfile(recDir, 'structure.oebin')));
cont = oe.continuous;  if isstruct(cont), cont = num2cell(cont); end
iNI = find(cellfun(@(c) contains(c.folder_name, 'NI-DAQ'), cont), 1);
nNI = double(cont{iNI}.num_channels);  chI = cont{iNI}.channels;  if isstruct(chI), chI = num2cell(chI); end
trigAI = I.settings.trigAI;  niV = double(chI{trigAI + 1}.bit_volts);  trigThrV = I.settings.trigThrV;
mmAP = memmapfile(apFile, 'Format', {'int16', [nCh I.nSampAP], 'x'});
mmNI = memmapfile(niFile, 'Format', {'int16', [nNI I.nSampNI], 'x'});
[bF, aF] = butter(3, bpHz/(fsAP/2), 'bandpass');

y = TT.chan.y_um(:);
if isempty(channels)
    o = TT.chan.orderTipToTop(:);  o = o(y(o) >= 3600 & y(o) <= 5800);
    channels = o(round(linspace(1, numel(o), 4)));
end
channels = channels(:)';
col = lines(numel(channels));
tr = TT.trial;

%% ---- figures -----------------------------------------------------------------------------------------------
cols = {zoomWinMs};  if showTrain, cols{end+1} = fullWinMs; end
fprintf('\n  trial  cond  label                          trigger width   artifact onset (ms after trigger)\n');
for c = conds(:)'
    ks = find(tr.cond == c, trialsPerCond);
    for t = ks(:)'
        fh = figure('Name', sprintf('%s trial %d cond %d', short, tr.trialNum(t), c), 'Color', 'w', 'Position', [40 40 650*numel(cols) 850]);
        onsetMs = NaN;
        for ci = 1:numel(cols)
            W = cols{ci};  padMs = 50;
            % ---- AP window (with padding for the filter) ----
            a1 = tr.trigAP(t) + round((W(1) - padMs)*fsAP/1e3);  a2 = tr.trigAP(t) + round((W(2) + padMs)*fsAP/1e3);
            a1 = max(a1, 1);  a2 = min(a2, I.nSampAP);
            X = double(mmAP.Data.x(:, a1:a2))*uv;
            tA = ((a1:a2) - tr.trigAP(t))/fsAP*1e3;                      % ms after the trigger
            pre = tA < -0.5 & tA >= W(1) - padMs;
            X = bsxfun(@minus, X, median(X(:, pre), 2));
            Xf = filtfilt(bF, aF, X(channels, :)')';
            st = [0, median(abs(diff(X, 1, 2)), 1)];                       % common step signal
            thr = artK*median(st(pre)) + artAbsUV;
            if ci == 1
                k = find(tA >= -0.5 & st > thr, 1);  if ~isempty(k), onsetMs = tA(k); end
            end
            % ---- NI window, converted to probe time ----
            n1 = tr.trigNI(t) + round(W(1)*fsNI/1e3) - 5;  n2 = tr.trigNI(t) + round(W(2)*fsNI/1e3) + 5;
            n1 = max(n1, 1);  n2 = min(n2, I.nSampNI);
            v = double(mmNI.Data.x(trigAI + 1, n1:n2))*niV;
            tN = (polyval(TT.sync.P, ((n1:n2) - 1)/fsNI) - (tr.trigAP(t) - 1)/fsAP)*1e3;   % ms after the trigger (probe clock)
            inW = tA >= W(1) & tA <= W(2);
            pt = tr.pulseMs{t};
            % ---- rows ----
            ax = subplot(4, numel(cols), ci);  hold(ax, 'on');
            stairs(ax, tN, v, 'k', 'LineWidth', 1);  yline(ax, trigThrV, 'b--');  xline(ax, 0, 'b-', 'LineWidth', 1);
            ylabel(ax, 'AI3 (V)');  title(ax, sprintf('trigger (NI %g kHz): blue = 2 V threshold and detected trigger', fsNI/1e3), 'FontSize', 9);
            ax(2) = subplot(4, numel(cols), numel(cols) + ci);  hold(ax(2), 'on');
            for q = 1:numel(channels), plot(ax(2), tA(inW), X(channels(q), inW), 'Color', col(q,:)); end
            ylabel(ax(2), 'raw (uV)');  title(ax(2), 'raw AP (median before the trigger removed)', 'FontSize', 9);
            legend(ax(2), arrayfun(@(ch) sprintf('ch %d (%d um)', ch, y(ch)), channels, 'UniformOutput', false), 'FontSize', 7, 'Location', 'best');
            ax(3) = subplot(4, numel(cols), 2*numel(cols) + ci);  hold(ax(3), 'on');
            for q = 1:numel(channels), plot(ax(3), tA(inW), Xf(q, inW), 'Color', col(q,:)); end
            ylabel(ax(3), 'filtered (uV)');  title(ax(3), sprintf('band-passed %g-%g Hz (no blanking)', bpHz), 'FontSize', 9);
            ax(4) = subplot(4, numel(cols), 3*numel(cols) + ci);  hold(ax(4), 'on');
            plot(ax(4), tA(inW), st(inW), 'k');  yline(ax(4), thr, 'r--');
            if ~isnan(onsetMs), xline(ax(4), onsetMs, 'r-', 'LineWidth', 1); end
            ylabel(ax(4), 'step (uV)');  xlabel(ax(4), 'ms after the trigger (probe clock)');
            title(ax(4), 'common step (median over all channels); red = threshold and detected artifact onset', 'FontSize', 9);
            for a = ax
                for p = pt(:)', xline(a, p, 'r:', 'HandleVisibility', 'off'); end
                xline(a, 0, 'b-', 'HandleVisibility', 'off');  xlim(a, W);  grid(a, 'on');
            end
            linkaxes(ax, 'x');
        end
        w = NaN;  if isfield(tr, 'trigWidthMs'), w = tr.trigWidthMs(t); end
        fprintf('  %5d  %4d  %-28s  %8.1f ms      %6.3f\n', tr.trialNum(t), c, TT.cond(c).label, w, onsetMs);
        sgtitle(fh, sprintf('%s  trial %d, cond %d: %s   (artifact onset %.3f ms after the trigger)', strrep(short, '_', '\_'), ...
                tr.trialNum(t), c, strrep(TT.cond(c).label, '_', '\_'), onsetMs), 'FontSize', 10);
    end
end
