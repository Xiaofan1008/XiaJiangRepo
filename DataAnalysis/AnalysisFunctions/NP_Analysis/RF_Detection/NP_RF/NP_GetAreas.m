function [lab, info] = NP_GetAreas(areaSource, areaFile, areasManual, badManual, y_um)
%% NP_GetAreas  area label of every NP channel (V1 / V2 / border / outside / bad / unassigned)
%
%   [lab, info] = NP_GetAreas(areaSource, areaFile, areasManual, badManual, y_um)
%
%   areaSource  'file'   : the ranges and bad channels of an area file saved by NP_AreaEntry_v1
%                          (<animal>_NP_Areas_<yyyymmdd>_<HHMM>.mat)
%               'manual' : areasManual / badManual typed in the calling script
%               'none'   : no areas: lab = [] and info.source = 'none'
%   areaFile    full path of the area file ('file' only)
%   areasManual rows {from_um  to_um  label} (a 4th column, e.g. confidence, is allowed and ignored); 'manual' only
%   badManual   channel numbers labelled 'bad' ('manual' only); [] = none
%   y_um        depth of every channel from the probe tip (um), e.g. TT.chan.y_um; its length = number of channels
%
%   lab{ch}     'V1' | 'V2' | 'border' | 'outside' | 'bad' | 'unassigned'   (nCh x 1 cell)
%   info        .source ('file' | 'manual' | 'none'), .file, .status, .text (one line for printouts / titles),
%               .ranges (rows {from to label}, sorted), .badCh, .rangeLabel{ch} (the label WITHOUT the bad override,
%               for scripts that do not exclude bad channels)
%
%   RULES (the same as NP_AreaEntry_v1): labels 'V1', 'V2', 'border', 'outside'; a range is [from, to) and the topmost
%   range also includes its 'to'; overlapping ranges are an error; depths in no range are 'unassigned' (warning).
%   With an area file, the depths saved in the file are compared with y_um: a mismatch (wrong penetration / probe map)
%   gives a warning, and the labels follow y_um (the session's own depths).
%   Written without access to MATLAB: please report any error.

y = double(y_um(:));  nCh = numel(y);
info = struct('source', lower(areaSource), 'file', '', 'status', '', 'text', 'areas: none', 'ranges', {{}}, ...
              'badCh', [], 'rangeLabel', {{}});
switch lower(areaSource)
    case 'none'
        lab = [];  return
    case 'file'
        assert(isfile(areaFile), 'Area file not found: %s', areaFile);
        L = load(areaFile, 'Areas');  A = L.Areas;
        rg = A.ranges(:, 1:3);  bad = A.badCh;
        if numel(A.y_um) ~= nCh
            warning('NPGA:nch', 'The area file has %d channels, the session %d: check that it is the right penetration.', numel(A.y_um), nCh);
        elseif any(abs(double(A.y_um(:)) - y) > 1)
            warning('NPGA:depth', 'Channel depths in the area file differ from the session''s: check the penetration / probe map.');
        end
        [~, nm] = fileparts(areaFile);
        info.file = areaFile;  info.status = A.status;
        info.text = sprintf('areas: %s (v%d, %s)', nm, A.version, A.status);
    case 'manual'
        assert(iscell(areasManual) && size(areasManual, 2) >= 3, 'areasManual must be rows {from to label}.');
        rg = areasManual(:, 1:3);  bad = badManual;
        info.status = 'manual';  info.text = 'areas: typed in the script (manual)';
    otherwise
        error('areaSource must be ''file'', ''manual'' or ''none''.');
end

% ---- ranges -> labels ----
allowed = {'V1', 'V2', 'border', 'outside'};
nR = size(rg, 1);  fr = cell2mat(rg(:, 1));  to = cell2mat(rg(:, 2));
for r = 1:nR
    assert(to(r) > fr(r), 'Area range %d: "to" (%g) must be larger than "from" (%g).', r, to(r), fr(r));
    assert(any(strcmp(rg{r, 3}, allowed)), 'Area range %d: label "%s" not allowed (%s).', r, rg{r, 3}, strjoin(allowed, ', '));
end
[~, o] = sort(fr);  rg = rg(o, :);  fr = fr(o);  to = to(o);
for r = 1:nR - 1
    assert(fr(r + 1) >= to(r), 'Area ranges overlap: %g-%g and %g-%g.', fr(r), to(r), fr(r + 1), to(r + 1));
end
lab = repmat({'unassigned'}, nCh, 1);
for r = 1:nR
    if r == nR, in = y >= fr(r) & y <= to(r); else, in = y >= fr(r) & y < to(r); end
    lab(in) = rg(r, 3);
end
nU = sum(strcmp(lab, 'unassigned'));
if nU > 0, warning('NPGA:unassigned', '%d channel(s) are in no area range: labelled ''unassigned''.', nU); end
info.rangeLabel = lab;  info.ranges = rg;

% ---- bad channels ----
bad = unique(bad(:)');
if ~isempty(bad)
    assert(all(bad >= 1 & bad <= nCh & bad == round(bad)), 'Bad channels must be channel numbers 1..%d.', nCh);
    lab(bad) = {'bad'};
end
info.badCh = bad;
end
