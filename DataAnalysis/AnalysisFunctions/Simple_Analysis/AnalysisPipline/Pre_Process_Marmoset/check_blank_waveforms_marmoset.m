%% CHECK_BLANK_WAVEFORMS_MARMOSET  (script; needs R2018b or newer)
% Plots, for chosen channels and conditions, the SAME trials from the three
% files written by the NN pipeline, so you can see what blanking and filtering did:
%   1) RAW        amplifier.dat            (before any processing)
%   2) BLANKED    amplifier_dn_sab.dat     (after artifact blanking + step removal)
%   3) MU FILTER  <name>.mu_sab.dat        (after the MU filter, stored x10 in the file)
% One figure per (condition, channel). Left column = wide window around the
% trigger, right column = zoom on the artifact. Several trials are overlaid
% in grey, the first one in colour. Black line = trigger (0 ms), magenta
% dashed = each pulse (from StimParams). Red dashed (bottom row) = approximate
% spike threshold.
%
% Usage: edit the settings, then press Run. Nothing is written or changed.
% Run it AFTER run_preprocess_marmoset has finished on that folder.

%% ============================ SETTINGS ============================
folderPath = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses';  % <-- EDIT
channels   = [9];       % amplifier channel NUMBERS to plot (1..nChn = position in amplifier.dat)
plotCond   = [];          % condition IDs to plot; [] = all conditions
nTrialsPerCond = 5;       % trials overlaid per condition
trialPick  = 'random';     % 'first' = first N trials of the condition, 'random' = N random trials
winMs      = [-20 100];   % wide window around the trigger, ms  [before after]
zoomMs     = [-2 15];     % zoom window around the trigger, ms
yLimRaw    = [];          % y limits in uV for the raw row, e.g. [-1500 1500]; [] = automatic
yLimBlank  = [];          % y limits for the blanked row, e.g. [-300 300]; [] = automatic
yLimMU     = [];          % y limits for the MU row, e.g. [-200 200]; [] = automatic
showThr    = true;        % draw an approximate spike threshold on the MU row
thrFac     = -4.5;        % same factor as allExtract_sab_1_stopAtMu (threshfac)
maxFigures = 40;          % safety stop: never open more than this many figures
muScale    = 10;          % SCALEFACTOR used when the MU file was written (allExtract: 10)

%% ============================ LOAD ============================
origDir = pwd;
cleanupObj = onCleanup(@() cd(origDir)); %#ok<NASGU>
cd(folderPath);

[amplifier_channels, fp] = read_Intan_RHS2000_file;
nChn = size(amplifier_channels,2);
FS   = fp.amplifier_sample_rate;
assert(all(channels >= 1 & channels <= nChn), 'channels must be between 1 and %d', nChn);

rawName = 'amplifier.dat';
dnName  = 'amplifier_dn_sab.dat';
muList  = dir('*.mu_sab.dat');
assert(~isempty(dir(rawName)), 'No %s in %s', rawName, folderPath);
assert(~isempty(dir(dnName)),  'No %s - run run_preprocess_marmoset first.', dnName);
assert(~isempty(muList),       'No *.mu_sab.dat - run run_preprocess_marmoset first.');
muName = muList(1).name;

fRaw = fopen(rawName,'r');
fDn  = fopen(dnName,'r');
fMu  = fopen(muName,'r');
closeFiles = onCleanup(@() fclose('all')); %#ok<NASGU>

% triggers (absolute sample numbers; -500 = no trigger on that trial)
tf = dir('*.trig.dat');
assert(~isempty(tf), 'No *.trig.dat in the folder.');
ft = fopen(tf(1).name,'r');
trig = fread(ft,[1 inf],'double');
fclose(ft);

% trial -> condition (same rule as the blanking function)
ef = dir('*_exp_datafile_*.mat');
assert(~isempty(ef), 'No *_exp_datafile_*.mat in the folder.');
S = load(ef(1).name,'TrialParams','StimParams');
condAll = cell2mat(S.TrialParams(2:end,2))';
nElect  = min(diff(find(diff(condAll)~=0)));
trialCond = condAll(1:nElect:end);
nTr = min(numel(trialCond), numel(trig));
if numel(trialCond) ~= numel(trig)
    warning('check:ntrials','%d trials in the exp file but %d triggers; using the first %d.', ...
        numel(trialCond), numel(trig), nTr);
end
trialCond = trialCond(1:nTr);
trig      = trig(1:nTr);

condIDs = unique(trialCond);
if ~isempty(plotCond), condIDs = intersect(condIDs, plotCond); end

% per-condition label and pulse times (ms after the trigger), from StimParams
label = containers.Map('KeyType','double','ValueType','char');
evMs  = containers.Map('KeyType','double','ValueType','any');
for c = condIDs
    t0   = find(trialCond == c, 1);
    rows = S.StimParams(1 + (t0-1)*nElect + (1:nElect), :);
    parts = {}; ev = [];
    for r = 1:size(rows,1)
        amp = rows{r,16};
        if isnumeric(amp) && ~isempty(amp) && amp > 0
            ptd = rows{r,6}/1000;  nP = rows{r,8};  per = rows{r,9}/1000;
            ev  = [ev, ptd + (0:nP-1)*per]; %#ok<AGROW>
            parts{end+1} = sprintf('%s %gms x%d %guA', toStr_(rows{r,1}), ptd, nP, amp); %#ok<AGROW>
        end
    end
    if isempty(parts), label(c) = sprintf('cond %d: no stim', c);
    else,              label(c) = sprintf('cond %d: %s', c, strjoin(parts,' + ')); end
    evMs(c) = ev;
end

% approximate spike threshold per channel, from the first <=20 s of the MU file
thr = nan(1,nChn);
if showThr
    dMu = dir(muName);
    nB = min(20*FS, floor(dMu.bytes/(2*nChn)));
    fseek(fMu,0,'bof');
    m = double(fread(fMu,[nChn nB],'int16'))/muScale;
    thr = thrFac * median(abs(m),2)' / 0.6745;
    clear m
end

%% ============================ PLOT ============================
nWide = round(diff(winMs)*FS/1000)+1;
rng(1);
nFig = 0;
for c = condIDs
    idx = find(trialCond == c & trig > 0);
    if isempty(idx), continue; end
    if strcmp(trialPick,'random')
        idx = idx(randperm(numel(idx), min(nTrialsPerCond,numel(idx))));
    else
        idx = idx(1:min(nTrialsPerCond,numel(idx)));
    end

    % read the wide window of every picked trial from the three files (all channels)
    W = nan(3, nChn, nWide, numel(idx));
    for k = 1:numel(idx)
        s0 = round(trig(idx(k)) + winMs(1)*FS/1000);
        W(1,:,:,k) = readSeg_(fRaw, nChn, s0, nWide, 0.195);
        W(2,:,:,k) = readSeg_(fDn,  nChn, s0, nWide, 0.195);
        W(3,:,:,k) = readSeg_(fMu,  nChn, s0, nWide, 1/muScale);
    end
    tms = winMs(1) + (0:nWide-1)*1000/FS;

    for ch = channels
        if nFig >= maxFigures
            warning('check:max','Stopped at %d figures (maxFigures).', maxFigures);
            return
        end
        nFig = nFig + 1;
        figure('Name', sprintf('cond %d | ch %d', c, ch), 'Color','w', 'Position',[80 60 1200 760]);
        rowName = {'RAW (uV)','BLANKED (uV)','MU FILTER (uV)'};
        yl      = {yLimRaw, yLimBlank, yLimMU};
        for row = 1:3
            for col = 1:2
                subplot(3,2,(row-1)*2+col); hold on;
                for k = numel(idx):-1:1
                    y = squeeze(W(row,ch,:,k));
                    if k == 1, plot(tms, y, 'Color',[0.85 0.2 0.1], 'LineWidth',1.1);
                    else,      plot(tms, y, 'Color',[0.65 0.65 0.65]); end
                end
                if col == 1, xlim(winMs); else, xlim(zoomMs); end
                if ~isempty(yl{row}), ylim(yl{row}); end
                xline(0,'k-');
                for e = evMs(c), xline(e,'m--'); end
                if row == 3 && showThr && ~isnan(thr(ch)), yline(thr(ch),'r--'); end
                if col == 1, ylabel(rowName{row}); end
                if row == 3, xlabel('ms from trigger'); end
                if row == 1
                    if col == 1, title('wide window'); else, title('zoom'); end
                end
                box off
            end
        end
        sgtitle(sprintf('%s   |   channel %d (%s)   |   %d trials: %s', label(c), ch, ...
            toStr_(amplifier_channels(ch).native_channel_name), numel(idx), mat2str(idx(:)')), ...
            'Interpreter','none','FontSize',10);
    end
end
fprintf('Done: %d figure(s).\n', nFig);

%% ============================ HELPERS ============================
function X = readSeg_(fid, nChn, s0, n, scale)
% Read n samples starting at sample s0 (1-based) for all channels of an
% interleaved int16 file; samples outside the file come back as NaN.
X = nan(nChn, n);
first = max(s0,1);
skip  = first - s0;                       % samples cut off at the start
if first - 1 < 0, return; end
fseek(fid, (first-1)*nChn*2, 'bof');
d = fread(fid, [nChn, n-skip], 'int16');
if isempty(d), return; end
X(:, skip+(1:size(d,2))) = double(d)*scale;
end

function s = toStr_(x)
if isnumeric(x), s = num2str(x); else, s = char(string(x)); end
end
