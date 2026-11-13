%% ============================================================
% NN_Train_FR_RespChnEntry_v1.m
% NN probe | pulse TRAIN data | FIRING-RATE analysis | RESPONSIVE-CHANNEL ENTRY | v1
%
% PURPOSE
%   You type the FINAL responsive channels (from NN_Train_FR_RespChnCheck_v1 and your own judgement) for each
%   pattern; the script checks the entries, turns them into a channel list for EVERY condition of the dataset,
%   prints a summary and saves ONE new file. The analysis scripts read that file later.
%   The recording data and all existing files are never changed.
%
% ENTRY FORMAT  (rows of  pattern, levels, channels)
%   entries = { ...
%     'sim',      'all',   [2 4 7 9 16:19];      % one set for every level of the pattern
%     'sim',      10:15,   [2 4 7 9 10 21];      % OVERRIDE for these levels (pulses per electrode)
%     'seq B->A', 'all',   [];                    % [] = reviewed, NO responsive channels
%   };
%   - pattern names: 'sim', 'seq A->B', 'seq B->A'  (only the ones present in the recovery file)
%   - levels: 'all', or a list of pulses per electrode (the pair levels: sim n, seq [n n])
%   - 'all' row: at most one per pattern; override rows win over the 'all' row; a level may be overridden once.
%   - a pattern without any row is NOT set (the analysis then falls back to its own manual channel list).
%   - half-steps of seq ([n n-1]) always get the 'all' set of their pattern.
%   - the channel set is one per dataset and pattern: it is used for the unthinned file and for every thinned
%     (Random / RandomSpikes / RandomEvent) variant, so the channels do not change between percent settings.
%
% OUTPUT FILE  <data_folder>/<folder name>.FR_RespChn.mat   with the struct  RespChn:
%   RespChn.respChn{cond}   channels of condition cond ([] = none or not set; see isSet)
%   RespChn.isSet(cond)     true if an entry covers this condition
%   RespChn.patName / patCond / patChain   patterns found, their final conditions and condition chains
%   RespChn.entries, chRange, rec_file, date   what was typed and where it came from
%   Use in an analysis script:   S = load(file, 'RespChn');  chans = S.RespChn.respChn{cond};
%   An existing file is NOT overwritten unless overwrite = true (then the old one is first copied to a dated backup).
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
rec_tag     = '';          % '' = *.sp_xia_TrainRecovery.mat (only used to read the condition table and trial list)
rec_file    = '';          % full path (overrides rec_tag)
chRange     = 1:32;        % valid channels

% PLACEHOLDER = what NN_Train_FR_RespChnCheck_v1 recommended for the CJ287 100 Hz 30 uA data. EDIT before saving.
entries = { ...
    'sim',      'all',  [2 4 7:10 14 16:23 26 27 29:32]; ...
    'seq B->A', 'all',  [2 4 6:20 22:27 29:32]; ...
    };

out_file  = '';            % '' = <data_folder>/<folder name>.FR_RespChn.mat
overwrite = true;         % true = replace an existing file (a dated backup of the old one is made first)
do_save   = true;          % false = only check and print

%% ================= LOAD CONDITION TABLE =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
if isempty(rec_file)
    if isempty(rec_tag), pat_ = '*.sp_xia_TrainRecovery.mat'; else, pat_ = ['*.sp_xia_TrainRecovery_' rec_tag '.mat']; end
    d = dir(fullfile(data_folder, pat_));  d = d(~startsWith({d.name}, '.'));
    assert(~isempty(d), 'No %s in %s.', pat_, data_folder);
    assert(numel(d) == 1, 'More than one %s in %s - set rec_file.', pat_, data_folder);
    rec_file = fullfile(data_folder, d(1).name);
end
assert(isfile(rec_file), 'Cannot find %s', rec_file);
R = load(rec_file, 'Recovery_params');
trialCond = R.Recovery_params.trialCond(:);
condTable = R.Recovery_params.condTable;
nCond = numel(condTable);
clear R

fe = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));  fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fullfile(data_folder, fe(1).name), 'StimParams', 'simultaneous_stim');
nEl = double(E.simultaneous_stim);
peK = zeros(nCond, nEl);                 % pulses per electrode of every condition
for k = 1:nCond
    t1 = find(trialCond == k, 1);
    if isempty(t1), continue; end
    for e = 1:nEl
        row = 1 + (t1-1)*nEl + e;
        if double(E.StimParams{row,16}) > 0, peK(k,e) = double(E.StimParams{row,8}); end
    end
end
clear E

%% ================= PATTERNS, CHAINS, PAIR LEVELS =================
kinds = arrayfun(@(c) char(c.kind), condTable, 'UniformOutput', false);
isFinal = arrayfun(@(c) isfield(c,'isFinal') && ~isempty(c.isFinal) && logical(c.isFinal), condTable);
patCond = find(isFinal(:)' & ismember(kinds(:)', {'sim','AB','BA'}));
assert(~isempty(patCond), 'No sim / AB / BA pattern found in the recovery file.');
nPat = numel(patCond);
patName = cell(1, nPat);  patChain = cell(1, nPat);  pairLvl = cell(1, nPat);  pairN = cell(1, nPat);
for p = 1:nPat
    c = patCond(p);  ch_ = c;
    while isfinite(condTable(ch_(1)).sourceCond), ch_ = [condTable(ch_(1)).sourceCond, ch_]; end %#ok<AGROW>
    patChain{p} = ch_;
    switch kinds{c}
        case 'sim', patName{p} = 'sim';
        case 'AB',  patName{p} = 'seq A->B';
        case 'BA',  patName{p} = 'seq B->A';
        otherwise,  patName{p} = kinds{c};
    end
    isPair = arrayfun(@(k) all(peK(k,:) == peK(k,1)) && peK(k,1) > 0, ch_);
    pairLvl{p} = find(isPair);                                    % positions in the chain
    pairN{p}   = arrayfun(@(k) max(peK(k,:)), ch_(isPair));       % pulses per electrode of those levels
end

%% ================= CHECK AND RESOLVE ENTRIES =================
assert(iscell(entries) && size(entries,2) == 3, 'entries must be a cell array with 3 columns: pattern, levels, channels.');
nE = size(entries, 1);
assert(isnumeric(chRange) && ~isempty(chRange), 'chRange must be a numeric list.');
respChn = cell(nCond, 1);  isSet = false(nCond, 1);
for k = 1:nCond, respChn{k} = zeros(1, 0); end
allRow = zeros(1, nPat);                      % entry row of the 'all' set per pattern
ovRow  = cell(1, nPat);                       % entry row of each override, by pulses per electrode
for p = 1:nPat, ovRow{p} = zeros(1, numel(pairN{p})); end
for r = 1:nE
    pn = entries{r,1};  lv = entries{r,2};  chn = entries{r,3};
    assert(ischar(pn) || isstring(pn), 'entries row %d: the pattern must be text.', r);
    pn = char(pn);
    p = find(strcmp(patName, pn), 1);
    assert(~isempty(p), 'entries row %d: pattern ''%s'' not found. Available: %s', r, pn, strjoin(patName, ', '));
    assert(isnumeric(chn), 'entries row %d: the channels must be a numeric list (use [] for none).', r);
    chn = unique(double(chn(:)'));
    assert(all(ismember(chn, chRange)), 'entries row %d (%s): channels outside chRange: %s', r, pn, mat2str(setdiff(chn, chRange)));
    if (ischar(lv) || isstring(lv))
        assert(strcmp(char(lv), 'all'), 'entries row %d: levels must be ''all'' or a list of numbers.', r);
        assert(allRow(p) == 0, 'entries row %d: pattern ''%s'' already has an ''all'' row (row %d).', r, pn, allRow(p));
        allRow(p) = r;
    else
        assert(isnumeric(lv) && ~isempty(lv), 'entries row %d: levels must be ''all'' or a list of numbers.', r);
        lv = unique(double(lv(:)'));
        assert(all(ismember(lv, pairN{p})), 'entries row %d (%s): levels %s do not exist. Pair levels are: %s', ...
            r, pn, mat2str(setdiff(lv, pairN{p})), mat2str(pairN{p}));
        for n = lv
            j = find(pairN{p} == n, 1);
            assert(ovRow{p}(j) == 0, 'entries row %d: level %d of ''%s'' is already overridden by row %d.', r, n, pn, ovRow{p}(j));
            ovRow{p}(j) = r;
        end
    end
end
for p = 1:nPat
    if allRow(p) == 0 && any(ovRow{p} > 0)
        warning('RCE:noall', 'Pattern ''%s'' has override rows but no ''all'' row: the other levels stay NOT set.', patName{p});
    end
    cc = patChain{p};
    for n = 1:numel(cc)
        k = cc(n);  r = 0;
        j = find(pairLvl{p} == n, 1);
        if ~isempty(j) && ovRow{p}(j) > 0, r = ovRow{p}(j); elseif allRow(p) > 0, r = allRow(p); end
        if r > 0
            respChn{k} = unique(double(entries{r,3}(:)'));  isSet(k) = true;
        end
    end
end

%% ================= PRINT =================
fprintf('\n============================================================\n');
fprintf('NN TRAIN FR RESPONSIVE-CHANNEL ENTRY (v1)\n');
fprintf('============================================================\n');
fprintf('Folder   : %s\nCond. table from : %s\n', data_folder, rec_file);
fprintf('Patterns found : %s\n', strjoin(cellfun(@(n, c) sprintf('%s (final cond %d)', n, c), patName, num2cell(patCond), 'UniformOutput', false), ', '));
for p = 1:nPat
    fprintf('\n--- %s ---\n', upper(patName{p}));
    if allRow(p) == 0 && ~any(ovRow{p} > 0)
        fprintf('  NOT SET (no entry): the analysis will use its own channel list.\n');
        continue;
    end
    if allRow(p) > 0
        c_ = unique(double(entries{allRow(p),3}(:)'));
        fprintf('  all levels   : %2d channels %s\n', numel(c_), fmtList_(c_));
    end
    for j = find(ovRow{p} > 0)
        c_ = unique(double(entries{ovRow{p}(j),3}(:)'));
        fprintf('  level n = %2d : %2d channels %s   (override)\n', pairN{p}(j), numel(c_), fmtList_(c_));
    end
    fprintf('  conditions   : ');
    for n = 1:numel(patChain{p})
        k = patChain{p}(n);
        fprintf('%d:%s ', k, ternary_(isSet(k), sprintf('%d', numel(respChn{k})), 'x'));
    end
    fprintf('\n  (cond:number of channels; x = not set)\n');
end

%% ================= SAVE =================
[~, folderName] = fileparts(regexprep(data_folder, '[\\/]+$', ''));
if isempty(out_file), out_file = fullfile(data_folder, [folderName '.FR_RespChn.mat']); end
RespChn = struct('respChn', {respChn}, 'isSet', isSet, 'patName', {patName}, 'patCond', patCond, 'patChain', {patChain}, ...
    'pairLevels_n', {pairN}, 'entries', {entries}, 'chRange', chRange, 'rec_file', rec_file, 'date', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    'version', 'NN_Train_FR_RespChnEntry_v1');
fprintf('\nOutput file : %s\n', out_file);
if ~do_save
    fprintf('do_save = false: nothing saved.\n');
elseif isfile(out_file) && ~overwrite
    error('The file already exists and overwrite = false. Nothing was changed. Set overwrite = true to replace it (a dated backup is made first).');
else
    if isfile(out_file)
        bak = strrep(out_file, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
        copyfile(out_file, bak);
        fprintf('Old file copied to %s\n', bak);
    end
    save(out_file, 'RespChn');
    fprintf('Saved.\n');
end

%% ================= HELPERS =================
function s = fmtList_(v)
v = sort(v(:)');
if isempty(v), s = '[]'; return; end
parts = {};  i = 1;
while i <= numel(v)
    j = i;
    while j < numel(v) && v(j+1) == v(j) + 1, j = j + 1; end
    if j == i,        parts{end+1} = sprintf('%d', v(i));          %#ok<AGROW>
    elseif j == i+1,  parts{end+1} = sprintf('%d %d', v(i), v(j)); %#ok<AGROW>
    else,             parts{end+1} = sprintf('%d:%d', v(i), v(j)); %#ok<AGROW>
    end
    i = j + 1;
end
s = ['[' strjoin(parts, ' ') ']'];
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
