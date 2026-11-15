%% ============================================================
% NP_AreaEntry_v1.m
% NP probe | AREA ENTRY (V1 / V2 / border / outside / bad) for ONE PENETRATION | v1
%
% PURPOSE
%   You type the depth ranges of the areas the probe crossed (decided from the noise-grid RFs, later from histology).
%   The script checks them, gives every channel a label, prints / plots the result and, when doSave = true, saves
%   ONE file per penetration in the animal folder, plus a readable .txt copy:
%       <animalDir>/<animal>_NP_Areas_<yyyymmdd>_<HHMM>.mat   (+ .txt)
%   The date and time come from the evidence recording's folder name (e.g. ..._2026-11-10_05-06-09 -> 20261110_0506).
%   Every session recorded at this probe position uses this file: paste its path into the analysis script.
%   Nothing else is changed (recordings, spike files, RF files are only read).
%
% RANGES  rows of  { from_um  to_um  label  confidence }, depth from the probe tip
%   - labels allowed: 'V1', 'V2', 'border', 'outside' ('bad' is given by badCh only)
%   - a range is [from, to): the shared edge belongs to the range ABOVE it; the topmost range also includes its 'to'
%   - ranges must not overlap; gaps and depths not covered are reported (those channels are 'unassigned')
%   - badCh: channels labelled 'bad' whatever their depth (dead / noisy channels)
%
% USING THE FILE IN ANOTHER SCRIPT
%   A = load(areaFile);  lab = A.Areas.label;          % lab{ch} = 'V1', 'V2', 'border', 'outside', 'bad', 'unassigned'
%   v1 = find(strcmp(lab, 'V1'));                       % channel numbers in V1
%
% SAVED  Areas (current version) and History (all earlier versions, oldest first):
%   Areas.label{ch}, Areas.y_um(ch)        label and depth of every channel (1 ... nCh)
%   Areas.ranges, badCh, surface_um, status, note     as typed
%   Areas.animal, penetration, evidenceDir, xmlFile, version, date, code
%   An existing file is first copied to <name>_backup_<yyyymmdd_HHMMSS>.mat; the .txt is rewritten.
%
% FIGURE (plotFigure = true): depth on the y axis, area bands shaded in every panel:
%   areas | MUA rate (noise grid) | RF strength (peak/null ratio) | RF azimuth | RF elevation | RF peak lag
%   The RF panels use noisegrid_rf.mat / noisegrid_mua.mat in the evidence recording when they exist
%   (channels with RF strength >= rfSnrMin; centre = weighted centroid of the smoothed map at the peak lag).
% Written without access to MATLAB: please report any error.
% ============================================================
clear; clc;
codeVersion = 'NP_AreaEntry_v1 (2026-11-15)';

%% ================= USER SETTINGS =================
animalDir   = '/Volumes/LaCie/Marmoset_Data/CJ287';
evidenceDir = '/Volumes/LaCie/Marmoset_Data/CJ287/CJ287_noisegrid_050605_2026-11-10_05-06-09';  % noise grid of this penetration

% one row per range:   from (um)   to (um)   label   confidence
areas = { ...
    0       2500   'V1'       'high';
    2500    3400   'border'   '-';
    3400    5400   'V2'       'high';
    5400    5900   'V2'       'medium';
    5900    7660   'outside'  '-';
    };
badCh      = [83 250];            % always 'bad', whatever their depth
surface_um = 5900;                % brain surface (um from the tip); NaN = unknown
status     = 'provisional (RF)';  % e.g. 'provisional (RF)' | 'histology'
note       = 'Border from the reversal of the RF azimuth progression (noise grid 05:06).';

doSave     = false;               % false = DRY RUN (printout + figure) | true = save the .mat and .txt
fileName   = '';                  % '' = <animal>_NP_Areas_<yyyymmdd>_<HHMM> from the evidence folder | or your own name

% ---- figure ----
plotFigure = true;
rfSnrMin   = 1.05;                % RF panels: only channels with peak/null ratio >= this
rfPol      = 'on';                % 'on' | 'off' map for the RF centres
flipY      = true;                % as FLIP_Y in NP_5 (applied on top of NP_4's coordinates) so the elevation matches NP_5
smoothSD   = 1;                   % Gaussian SD (randels) of the map before the centroid; 0 = none
fontSize   = 9;
%% =================================================

allowed = {'V1', 'V2', 'border', 'outside'};
cols = struct('V1', [0.35 0.55 0.95], 'V2', [0.95 0.6 0.25], 'border', [0.75 0.75 0.75], ...
              'outside', [0.9 0.9 0.9], 'bad', [0.9 0.1 0.1], 'unassigned', [1 1 1]);

%% ---- channel depths from the evidence recording's settings.xml -------------------------------------------------------
assert(isfolder(evidenceDir), 'Evidence folder not found: %s', evidenceDir);
[y, xmlFile] = readDepths_(evidenceDir);
nCh = numel(y);
fprintf('\n============================================================\nNP AREA ENTRY (v1)\n============================================================\n');
fprintf('evidence : %s\nprobe    : %d channels, %g to %g um from the tip (%s)\n', evidenceDir, nCh, min(y), max(y), xmlFile);

%% ---- names ---------------------------------------------------------------------------------------------------------------
[~, animal] = fileparts(animalDir);
if isempty(fileName)
    [~, ev] = fileparts(evidenceDir);
    t = regexp(ev, '(\d{4})-(\d{2})-(\d{2})_(\d{2})-(\d{2})-\d{2}$', 'tokens', 'once');
    assert(~isempty(t), 'Cannot read the date/time from the evidence folder name "%s": set fileName yourself.', ev);
    fileName = sprintf('%s_NP_Areas_%s%s%s_%s%s', animal, t{1}, t{2}, t{3}, t{4}, t{5});
end
outMat = fullfile(animalDir, [fileName '.mat']);  outTxt = fullfile(animalDir, [fileName '.txt']);
fprintf('file     : %s\n', outMat);

%% ---- check the ranges ------------------------------------------------------------------------------------------------
assert(iscell(areas) && size(areas, 2) == 4, 'areas must have 4 columns: from, to, label, confidence.');
nR = size(areas, 1);  fr = zeros(nR, 1);  to = zeros(nR, 1);
for r = 1:nR
    fr(r) = areas{r, 1};  to(r) = areas{r, 2};
    assert(to(r) > fr(r), 'Range %d: "to" (%g) must be larger than "from" (%g).', r, to(r), fr(r));
    assert(any(strcmp(areas{r, 3}, allowed)), 'Range %d: label "%s" not allowed (use %s).', r, areas{r, 3}, strjoin(allowed, ', '));
end
[~, o] = sort(fr);  areas = areas(o, :);  fr = fr(o);  to = to(o);
for r = 1:nR - 1
    assert(fr(r + 1) >= to(r), 'Ranges overlap: %g-%g and %g-%g.', fr(r), to(r), fr(r + 1), to(r + 1));
    if fr(r + 1) > to(r), warning('NPAE:gap', 'Gap between %g and %g um: channels there are UNASSIGNED.', to(r), fr(r + 1)); end
end
if min(y) < fr(1), warning('NPAE:cover', 'Depths %g-%g um are below the first range: UNASSIGNED.', min(y), fr(1)); end
if max(y) > to(end), warning('NPAE:cover', 'Depths %g-%g um are above the last range: UNASSIGNED.', to(end), max(y)); end
badCh = unique(badCh(:)');
assert(all(badCh >= 1 & badCh <= nCh & badCh == round(badCh)), 'badCh must be channel numbers 1..%d.', nCh);

%% ---- label every channel -------------------------------------------------------------------------------------------
label = repmat({'unassigned'}, nCh, 1);  rangeOf = zeros(nCh, 1);
for r = 1:nR
    if r == nR, in = y >= fr(r) & y <= to(r); else, in = y >= fr(r) & y < to(r); end
    label(in) = areas(r, 3);  rangeOf(in) = r;
end
label(badCh) = {'bad'};

%% ---- printout -------------------------------------------------------------------------------------------------------
fprintf('status   : %s\nsurface  : %s um\n\nRANGES (um from the tip)\n', status, num2str(surface_um));
for r = 1:nR
    ch = find(rangeOf == r & ~strcmp(label, 'bad'))';
    fprintf('  %5g - %5g  %-8s %-7s %3d ch  %s\n', fr(r), to(r), areas{r, 3}, areas{r, 4}, numel(ch), compact_(ch));
end
fprintf('\nCHANNELS PER LABEL\n');
for lb = [allowed, {'bad', 'unassigned'}]
    ch = find(strcmp(label, lb{1}))';
    if isempty(ch) && strcmp(lb{1}, 'unassigned'), continue; end
    fprintf('  %-10s %3d  %s\n', lb{1}, numel(ch), compact_(ch));
end

%% ---- figure ---------------------------------------------------------------------------------------------------------
if plotFigure
    % evidence from the noise grid (optional)
    rate = nan(nCh, 1);  snr = nan(nCh, 1);  lag = nan(nCh, 1);  cx = nan(nCh, 1);  cy = nan(nCh, 1);
    mf = dir(fullfile(evidenceDir, '**', 'noisegrid_mua.mat'));  mf = mf(~startsWith({mf.name}, '._'));
    if ~isempty(mf)
        M = load(fullfile(mf(1).folder, mf(1).name), 'rate', 'channels');
        rate(M.channels) = M.rate;
    end
    rfF = dir(fullfile(evidenceDir, '**', 'noisegrid_rf.mat'));  rfF = rfF(~startsWith({rfF.name}, '._'));
    if ~isempty(rfF)
        fn = fullfile(rfF(1).folder, rfF(1).name);  zName = ternary_(strcmpi(rfPol, 'off'), 'zOff', 'zOn');
        fprintf('\nreading RF maps (%s) from %s ...\n', zName, fn);
        Rf = load(fn, 'channels', zName, 'peak_bin', 'peak_lag_ms', 'snr', 'size_h', 'size_w', 'x_deg', 'y_deg');
        H = double(Rf.size_h);  W = double(Rf.size_w);  xC = Rf.x_deg(:)';  yC = Rf.y_deg(:)';
        if flipY, yC = fliplr(yC); end
        Z = Rf.(zName);  nB = size(Z, 2);  chs = Rf.channels(:)';
        snr(chs) = Rf.snr;  lag(chs) = Rf.peak_lag_ms;
        ker = 1;  if smoothSD > 0, kx = -ceil(3*smoothSD):ceil(3*smoothSD); g = exp(-0.5*(kx/smoothSD).^2); ker = (g'*g)/sum(g)^2; end
        for k = find(Rf.snr(:)' >= rfSnrMin)
            b = Rf.peak_bin(k);  w = max(1, b - 1):min(nB, b + 1);
            m = reshape(mean(double(Z(:, w, k)), 2), H, W);
            if numel(ker) > 1, m = conv2(m, ker, 'same'); end
            [pk, ip] = max(m(:));  [r0, c0] = ind2sub([H W], ip);
            rr = max(1, r0 - 4):min(H, r0 + 4);  cc = max(1, c0 - 4):min(W, c0 + 4);
            sub = m(rr, cc);  wgt = sub.*(sub >= pk/2);  [CC, RR] = meshgrid(cc, rr);
            cx(chs(k)) = sum(wgt(:).*reshape(xC(CC), [], 1))/sum(wgt(:));
            cy(chs(k)) = sum(wgt(:).*reshape(yC(RR), [], 1))/sum(wgt(:));
        end
        clear Z Rf
        lag(snr < rfSnrMin) = NaN;
    end

    fh = figure('Name', [fileName ' areas'], 'Color', 'w', 'Position', [40 40 1500 800]);
    nPan = 6;  gL = 0.05;  gR = 0.015;  gap = 0.02;  w = (1 - gL - gR - (nPan - 1)*gap)/nPan;
    yl = [min(y) max(y)];  ok = ~strcmp(label, 'bad');
    pan = {'areas', 'MUA rate (spikes/s)', 'RF strength (peak/null)', 'RF azimuth (deg)', 'RF elevation (deg)', 'RF peak lag (ms)'};
    vals = {[], rate, snr, cx, cy, lag};
    for i = 1:nPan
        ax = axes(fh, 'Position', [gL + (i - 1)*(w + gap), 0.08, w, 0.82]);  hold(ax, 'on');
        if i == 1, xl = [0 1]; else
            v = vals{i}(ok);  v = v(isfinite(v));
            if isempty(v), xl = [0 1]; else, xl = [min(v) max(v)] + [-0.05 0.05]*max(max(v) - min(v), eps); end
            if i == 2 && ~isempty(v), set(ax, 'XScale', 'log'); xl = [max(min(v), 0.1) max(v)].*[0.8 1.25]; end
        end
        for r = 1:nR                                               % area bands
            c = cols.(areas{r, 3});  a = max(fr(r), yl(1));  b = min(to(r), yl(2));
            patch(ax, xl([1 2 2 1]), [a a b b], c, 'EdgeColor', 'none', 'FaceAlpha', ternary_(i == 1, 1, 0.25));
            if i == 1, text(ax, 0.5, (a + b)/2, sprintf('%s\n%g-%g', areas{r, 3}, fr(r), to(r)), 'HorizontalAlignment', 'center', ...
                            'FontSize', fontSize, 'FontWeight', 'bold'); end
        end
        if ~isnan(surface_um), plot(ax, xl, [surface_um surface_um], 'k--', 'LineWidth', 1); end
        if i == 1
            plot(ax, repmat(0.92, 1, numel(badCh)), y(badCh), 'x', 'Color', cols.bad, 'MarkerSize', 7, 'LineWidth', 1.5);
            set(ax, 'XTick', []);
        else
            plot(ax, vals{i}(ok), y(ok), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 8);
            if i == 3, plot(ax, [rfSnrMin rfSnrMin], yl, ':', 'Color', [0.3 0.3 0.3]); end
            if i >= 4 && all(isnan(vals{i})), text(ax, mean(xl), mean(yl), 'no RF file', 'HorizontalAlignment', 'center'); end
        end
        xlim(ax, xl);  ylim(ax, yl);  grid(ax, 'on');  set(ax, 'FontSize', fontSize, 'Layer', 'top');
        title(ax, pan{i}, 'FontSize', fontSize + 1);
        if i == 1, ylabel(ax, 'depth from the tip (\mum)'); else, set(ax, 'YTickLabel', []); end
    end
    sgtitle(fh, sprintf('%s: areas along the probe (%s, %s)', strrep(fileName, '_', '\_'), status, ternary_(doSave, 'SAVED', 'DRY RUN')), ...
            'FontSize', fontSize + 3);
end

%% ---- save -------------------------------------------------------------------------------------------------------------
if doSave
    History = struct([]);  ver = 1;
    if isfile(outMat)
        O = load(outMat);
        if isfield(O, 'History'), History = O.History; end
        if isempty(History), History = O.Areas; else, History(end + 1) = O.Areas; end
        ver = O.Areas.version + 1;
        bak = strrep(outMat, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
        copyfile(outMat, bak);  fprintf('\nOld file copied to %s\n', bak);
    end
    Areas = struct();
    Areas.label = label;  Areas.y_um = y;  Areas.ranges = areas;  Areas.badCh = badCh;  Areas.surface_um = surface_um;
    Areas.status = status;  Areas.note = note;  Areas.animal = animal;  Areas.penetration = fileName;
    Areas.evidenceDir = evidenceDir;  Areas.xmlFile = xmlFile;  Areas.version = ver;
    Areas.date = datestr(now, 'yyyy-mm-dd HH:MM:SS');  Areas.code = codeVersion;
    if ~isempty(History), History = orderfields_(History, Areas); end
    save(outMat, 'Areas', 'History');
    writeTxt_(outTxt, Areas, History, rangeOf, nR);
    fprintf('\nSAVED: %s (version %d)\n       %s\n', outMat, ver, outTxt);
else
    fprintf('\nDRY RUN: nothing saved. Set doSave = true to save %s (+ .txt)\n', outMat);
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function [y, xmlFile] = readDepths_(root)
% electrode y positions (um from the tip) from the Open Ephys settings.xml (as in NP_5)
fx = dir(fullfile(root, '**', 'settings.xml'));  fx = fx(~startsWith({fx.name}, '._'));
y = [];  xmlFile = '';
for q = 1:numel(fx)
    cand = fullfile(fx(q).folder, fx(q).name);
    try, txt = fileread(cand); catch, continue; end
    tok = regexp(txt, '<ELECTRODE_YPOS\s([^/>]*)/?>', 'tokens', 'once');
    if isempty(tok), continue; end
    t = regexp(tok{1}, 'CH(\d+)\s*=\s*"(-?[\d.]+)"', 'tokens');
    if isempty(t), continue; end
    iy = cellfun(@(c) str2double(c{1}), t) + 1;  vy = cellfun(@(c) str2double(c{2}), t);
    y = nan(max(iy), 1);  y(iy) = vy;  xmlFile = cand;  break;
end
assert(~isempty(y) && ~any(isnan(y)), 'Could not read the electrode depths from settings.xml below %s.', root);
end

function s = compact_(v)
% channel list in MATLAB form, runs with a constant step of 1 or 2 written as a:b or a:2:b
if isempty(v), s = '[]'; return; end
v = sort(v(:)');  parts = {};  i = 1;
while i <= numel(v)
    j = i;
    if i < numel(v) && any(v(i + 1) - v(i) == [1 2])
        d = v(i + 1) - v(i);  while j < numel(v) && v(j + 1) - v(j) == d, j = j + 1; end
        if j - i < 2 && d == 2, j = i; end                         % a step-2 run needs at least 3 values
    end
    if j == i, parts{end+1} = sprintf('%d', v(i)); %#ok<AGROW>
    elseif v(j) - v(j - 1) == 1, parts{end+1} = sprintf('%d:%d', v(i), v(j)); %#ok<AGROW>
    else, parts{end+1} = sprintf('%d:2:%d', v(i), v(j)); %#ok<AGROW>
    end
    i = j + 1;
end
s = ['[' strjoin(parts, ' ') ']'];
end

function H = orderfields_(H, A)
% make the history entries have the same fields as the current version (missing fields = [])
f = fieldnames(A);
for q = 1:numel(f), if ~isfield(H, f{q}), [H.(f{q})] = deal([]); end, end
extra = setdiff(fieldnames(H), f);  if ~isempty(extra), H = rmfield(H, extra); end
H = orderfields(H, A);
end

function writeTxt_(fn, A, Hst, rangeOf, nR)
fid = fopen(fn, 'w');  assert(fid > 0, 'Cannot write %s', fn);
c = onCleanup(@() fclose(fid));
fprintf(fid, 'NP AREAS: %s   (version %d, %s)\n', A.penetration, A.version, A.date);
fprintf(fid, 'status   : %s\nevidence : %s\nsurface  : %s um from the tip\nnote     : %s\n\n', A.status, A.evidenceDir, num2str(A.surface_um), A.note);
fprintf(fid, 'RANGES (um from the tip; [from, to), the top range includes its end)\n');
for r = 1:nR
    ch = find(rangeOf == r & ~strcmp(A.label, 'bad'))';
    fprintf(fid, '  %5g - %5g  %-8s %-7s %3d channels  %s\n', A.ranges{r, 1}, A.ranges{r, 2}, A.ranges{r, 3}, A.ranges{r, 4}, numel(ch), compact_(ch));
end
fprintf(fid, '\nCHANNELS PER LABEL\n');
for lb = {'V1', 'V2', 'border', 'outside', 'bad', 'unassigned'}
    ch = find(strcmp(A.label, lb{1}))';
    if ~isempty(ch) || ~strcmp(lb{1}, 'unassigned'), fprintf(fid, '  %-10s %3d  %s\n', lb{1}, numel(ch), compact_(ch)); end
end
fprintf(fid, '\nUSE:  A = load(''<this file>.mat'');  lab = A.Areas.label;  v1 = find(strcmp(lab, ''V1''));\n');
if ~isempty(Hst)
    fprintf(fid, '\nEARLIER VERSIONS\n');
    for q = 1:numel(Hst)
        rg = Hst(q).ranges;  txt = strjoin(arrayfun(@(r) sprintf('%g-%g %s', rg{r, 1}, rg{r, 2}, rg{r, 3}), 1:size(rg, 1), ...
                                  'UniformOutput', false), ' | ');
        fprintf(fid, '  v%d  %s  %s: %s\n', Hst(q).version, Hst(q).date, Hst(q).status, txt);
    end
end
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
