%% QUICK_COMPARE_ELEC_VIS_NN  (script; needs R2018b or newer)
% Quick look: raster + PSTH for ANY set of conditions on the same axes, so you
% can put an electrical condition (e.g. the one-pulse case) next to visual
% condition(s), or plot visual-only conditions on their own. Nothing is saved.
%
% Per channel: one figure, raster on top (conditions stacked in blocks, one
% colour per condition), PSTH below. An extra figure shows the population
% PSTH (mean rate per channel over all chosen channels).
%
% HOW TO USE
%   1) Set folderPath, then run once with condList = [] : it prints every
%      condition (ID, kind, number of trials, description) and stops.
%   2) Put the IDs you want in condList, e.g. one electrical + visual ones,
%      and run again.
%
% NEEDS (in folderPath): *.sp_xia.mat (variable sp_clipped), *.trig.dat,
% *_exp_datafile_*.mat. To know which "no stimulation" conditions are VISUAL
% (and their contrast / polarity) it also needs the *.conds.mat that the
% stimulus PC saved next to the Neurostim file: copy it into folderPath or set
% condFile. Without it, conditions are listed by their electrical content only.
%
% NOTES
%  * Time zero = the trigger for every condition. For visual trials the screen
%    adds a delay after the trigger (not corrected here; this is a quick look).
%  * Grey band = period blanked by the artifact removal, estimated from the
%    pulse times (first pulse - 1 ms to last pulse + blankAfterLastMs). The
%    PSTH is left empty (gap) where more than 80% of a bin is blanked.

%% ============================ SETTINGS ============================
folderPath   = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq1_50Hz';   % <-- EDIT
analysisPath = '';          % optional: folder added with addpath(genpath(..)) so Depth_s and
                            % cleanTrig_sabquick are found, e.g. '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'
condFile     = '';          % '' = look for *.conds.mat inside folderPath
Electrode_Type = 1;         % same as your raster code: 0 single shank rigid, 1 single shank flex, 2 four shank flex

condList     = [7 29];          % condition IDs to plot together, e.g. [3 52 53]. [] = print the list and stop
chanList     = 20:32;       % channel indices into d = Depth_s(Electrode_Type) (same idea as raster_chn_start:raster_chn_end)

ras_win      = [-40 200];   % ms around the trigger
binMs        = 5;           % PSTH bin, ms
smoothMs     = 10;          % Gaussian smoothing sigma for the PSTH, ms (0 = none)
maxTrialsPerCond = Inf;     % use at most this many trials of each condition (first ones)
blankAfterLastMs = 6;       % blank is assumed to end this long after the last pulse (artifact recovery cap)
showPopulation = true;      % extra figure: mean rate over all chosen channels
maxFigures   = 20;          % never open more than this many per-channel figures
FS           = 30000;       % sampling rate used for the trigger file

%% ============================ LOAD ============================
if ~isempty(analysisPath), addpath(genpath(analysisPath)); end
assert(isfolder(folderPath), 'Folder not found: %s', folderPath);
origDir = pwd;
cleanupObj = onCleanup(@() cd(origDir)); %#ok<NASGU>
cd(folderPath);

% spikes
f = listClean_('*.sp_xia.mat');
assert(~isempty(f), 'No *.sp_xia.mat in %s', folderPath);
SP = load(f(1).name, 'sp_clipped');
assert(isfield(SP,'sp_clipped'), 'Variable sp_clipped not found in %s', f(1).name);
sp_clipped = SP.sp_clipped;

% triggers (absolute samples; -500 = no trigger on that trial)
f = listClean_('*.trig.dat');
if isempty(f)
    cleanTrig_sabquick;
    f = listClean_('*.trig.dat');
end
assert(~isempty(f), 'No *.trig.dat in %s', folderPath);
ft = fopen(f(1).name,'r');
trig = fread(ft,[1 inf],'double');
fclose(ft);

% trial -> condition
f = listClean_('*_exp_datafile_*.mat');
assert(~isempty(f), 'No *_exp_datafile_*.mat in %s', folderPath);
E = load(f(1).name, 'TrialParams', 'StimParams', 'simultaneous_stim');
nEl = E.simultaneous_stim;
condAll   = cell2mat(E.TrialParams(2:end,2));
trialCond = condAll(1:nEl:end)';
nTr = min(numel(trialCond), numel(trig));
if numel(trialCond) ~= numel(trig)
    warning('quick:ntr','%d trials in the exp file but %d triggers; using the first %d.', ...
        numel(trialCond), numel(trig), nTr);
end
trialCond = trialCond(1:nTr);
trig      = trig(1:nTr);

% conds file (kind / label of every condition)
if isempty(condFile)
    cf = listClean_('*.conds.mat');
    if ~isempty(cf), condFile = fullfile(folderPath, cf(1).name); end
end
condKind  = {};   condLabel = {};
if ~isempty(condFile) && exist(condFile,'file') == 2
    K = load(condFile);
    nC = numel(K.cond);
    condKind = repmat({'electrical'}, 1, nC);
    if isfield(K,'info')
        condKind = {K.info.kind};
    elseif isfield(K,'isVisual')
        condKind(K.isVisual(:)' ~= 0) = {'visual'};
    end
    condLabel = repmat({''}, 1, nC);
    if isfield(K.cond,'label')
        condLabel = {K.cond.label};
    end
else
    warning('quick:nocond','No conds file found: visual and blank trials cannot be told apart (both show as "no stim").');
end

% what each condition delivers (from StimParams) and where its blank is
ids = unique(trialCond);
descr = containers.Map('KeyType','double','ValueType','char');
blank = containers.Map('KeyType','double','ValueType','any');
nTrOf = containers.Map('KeyType','double','ValueType','double');
for k = ids
    t0   = find(trialCond == k, 1);
    rows = E.StimParams(1 + (t0-1)*nEl + (1:nEl), :);
    parts = {}; ev = [];
    for r = 1:size(rows,1)
        amp = rows{r,16};
        if isnumeric(amp) && ~isempty(amp) && amp > 0
            ptd = rows{r,6}/1000;  nP = rows{r,8};  per = rows{r,9}/1000;
            ev  = [ev, ptd + (0:nP-1)*per]; %#ok<AGROW>
            parts{end+1} = sprintf('%s %gms x%d %guA', toStr_(rows{r,1}), ptd, nP, amp); %#ok<AGROW>
        end
    end
    if isempty(parts), descr(k) = 'no electrical stim'; blank(k) = [];
    else,              descr(k) = strjoin(parts,' + ');  blank(k) = [min(ev)-1, max(ev)+blankAfterLastMs]; end
    nTrOf(k) = sum(trialCond == k & trig > 0);
end
kindOf = @(k) ternary_(~isempty(condKind) && k <= numel(condKind), @() condKind{k}, @() '?');
labOf  = @(k) ternary_(~isempty(condLabel) && k <= numel(condLabel) && ~isempty(condLabel{k}), @() condLabel{k}, @() '');

fprintf('\nConditions in %s\n', folderPath);
for k = ids
    fprintf('  [%3d] %-10s n=%3d  %s  %s\n', k, kindOf(k), nTrOf(k), labOf(k), descr(k));
end
if isempty(condList)
    fprintf('\ncondList is empty: pick condition IDs from the list above and run again.\n');
    return
end
assert(all(ismember(condList, ids)), 'condList contains an ID that is not in the list above.');

d = Depth_s(Electrode_Type);
chans = d(chanList);

%% ============================ PREPARE ============================
nSel  = numel(condList);
cmap  = lines(max(nSel,2));
edges = ras_win(1):binMs:ras_win(2);
ctr   = edges(1:end-1) + binMs/2;
nb    = numel(edges) - 1;
bin_s = binMs/1000;
if smoothMs > 0
    sb = smoothMs/binMs;  half = ceil(3*sb);
    g  = exp(-0.5*((-half:half)/sb).^2);  g = g/sum(g);
else
    g = 1;
end

trialsOf = cell(1,nSel);  validOf = cell(1,nSel);  shortLab = cell(1,nSel);
for j = 1:nSel
    k  = condList(j);
    tr = find(trialCond == k & trig > 0);
    tr = tr(1:min(numel(tr), maxTrialsPerCond));
    trialsOf{j} = tr;
    v = ones(1,nb);
    bw = blank(k);
    if ~isempty(bw)
        ov = max(0, min(edges(2:end), bw(2)) - max(edges(1:end-1), bw(1))) / binMs;
        v = 1 - ov;
    end
    validOf{j} = v;
    lab = labOf(k);
    if isempty(lab), lab = descr(k); end
    shortLab{j} = sprintf('c%d %s (n=%d)', k, lab, numel(tr));
end

countsCh = cell(1,nSel);                 % per condition: nChan x nb spike counts summed over trials
for j = 1:nSel, countsCh{j} = zeros(numel(chans), nb); end

%% ============================ PER-CHANNEL FIGURES ============================
nFig = 0;
for ci = 1:numel(chans)
    ch = chans(ci);
    spk = [];
    if ch >= 1 && ch <= numel(sp_clipped) && ~isempty(sp_clipped{ch})
        spk = sp_clipped{ch}(:,1);
    end
    doFig = nFig < maxFigures;
    if doFig
        figure('Color','w','Name',sprintf('chan idx %d (amp ch %d)', chanList(ci), ch), 'Position',[100 80 900 700]);
        ax1 = subplot(4,1,1:3); hold(ax1,'on'); box(ax1,'off');
        ax2 = subplot(4,1,4);   hold(ax2,'on'); box(ax2,'off');
        nFig = nFig + 1;
    end
    yCur = 0;  yTick = zeros(1,nSel);
    for j = 1:nSel
        k = condList(j);  tr = trialsOf{j};  col = cmap(j,:);
        XS = cell(1,numel(tr));  YS = cell(1,numel(tr));
        for t = 1:numel(tr)
            t0 = trig(tr(t))/FS*1000;
            tt = [];
            if ~isempty(spk)
                tt = spk(spk >= t0+ras_win(1) & spk <= t0+ras_win(2)) - t0;
            end
            countsCh{j}(ci,:) = countsCh{j}(ci,:) + histcounts(tt, edges);
            y = yCur + t;
            tt = tt(:)';
            XS{t} = [tt; tt; nan(1,numel(tt))];
            YS{t} = [repmat(y-0.4,1,numel(tt)); repmat(y+0.4,1,numel(tt)); nan(1,numel(tt))];
        end
        if doFig
            if ~isempty(tr)
                X = [XS{:}];  Y = [YS{:}];
                if ~isempty(X), plot(ax1, X(:), Y(:), 'Color', col, 'LineWidth', 1); end
                bw = blank(k);
                if ~isempty(bw)
                    patch(ax1, [bw(1) bw(2) bw(2) bw(1)], [yCur yCur yCur+numel(tr) yCur+numel(tr)], [0.5 0.5 0.5], ...
                        'FaceAlpha',0.15, 'EdgeColor','none');
                end
            end
            yTick(j) = yCur + max(numel(tr),1)/2;
            if j < nSel, yline(ax1, yCur + numel(tr), 'Color',[0.8 0.8 0.8]); end
        end
        yCur = yCur + numel(tr);
        if doFig
            rate = smoothRate_(countsCh{j}(ci,:), validOf{j}, numel(tr), g, bin_s, 1);
            plot(ax2, ctr, rate, 'Color', col, 'LineWidth', 1.6);
        end
    end
    if doFig
        xline(ax1, 0, 'k--');  xlim(ax1, ras_win);  ylim(ax1, [0 max(yCur,1)]);
        yticks(ax1, yTick);  yticklabels(ax1, shortLab);  set(ax1,'TickLabelInterpreter','none');
        title(ax1, sprintf('Channel index %d (amplifier channel %d)', chanList(ci), ch), 'Interpreter','none');
        xline(ax2, 0, 'k--');  xlim(ax2, ras_win);
        xlabel(ax2,'ms from trigger');  ylabel(ax2,'Rate (sp/s)');
    end
end
if numel(chans) > maxFigures
    fprintf('Only the first %d channel figures were drawn (maxFigures); the population plot uses all %d channels.\n', ...
        maxFigures, numel(chans));
end

%% ============================ POPULATION ============================
if showPopulation
    figure('Color','w','Name','population PSTH','Position',[120 100 900 380]);
    ax = axes; hold(ax,'on'); box(ax,'off');
    for j = 1:nSel
        rate = smoothRate_(sum(countsCh{j},1), validOf{j}, numel(trialsOf{j}), g, bin_s, numel(chans));
        plot(ax, ctr, rate, 'Color', cmap(j,:), 'LineWidth', 1.8);
    end
    xline(ax, 0, 'k--');  xlim(ax, ras_win);
    xlabel(ax,'ms from trigger');  ylabel(ax,'Rate per channel (sp/s)');
    title(ax, sprintf('Population: mean over %d channels', numel(chans)));
    legend(ax, shortLab, 'Box','off', 'Location','northeast', 'Interpreter','none');
end
fprintf('Done.\n');

%% ============================ HELPERS ============================
function f = listClean_(pat)
% dir() without macOS "._" junk files
f = dir(pat);
if ~isempty(f), f = f(~startsWith({f.name}, '.')); end
end

function s = toStr_(x)
if isnumeric(x), s = num2str(x); else, s = char(string(x)); end
end

function out = ternary_(cond, a, b)
if cond, out = a(); else, out = b(); end
end

function rate = smoothRate_(counts, valid, nTr, g, bin_s, nCh)
% Smoothed rate in sp/s per channel. Counts and valid (non-blanked) time are
% smoothed separately, so blanked bins do not pull the rate towards zero.
num = conv(counts, g, 'same');
den = conv(valid*nTr*bin_s*nCh, g, 'same');
rate = num ./ den;
rate(conv(valid, g, 'same') < 0.2) = NaN;
end
