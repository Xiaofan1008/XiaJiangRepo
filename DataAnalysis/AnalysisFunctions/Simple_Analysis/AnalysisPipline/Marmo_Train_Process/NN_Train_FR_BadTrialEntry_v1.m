%% ============================================================
% NN_Train_FR_BadTrialEntry_v1.m
% NN probe | pulse TRAIN data | FIRING-RATE analysis | BAD-TRIAL ENTRY | v1
%
% PURPOSE
%   You type the FINAL bad trials of every condition (level) you reviewed with NN_Train_FR_BadTrialCheck_v2
%   (or the RandomEvent version); the script checks the entries, prints a summary and saves ONE new file that the
%   analysis scripts read. The recording data and all existing files are never changed.
%
% ENTRY FORMAT  (the same rows the bad-trial check prints as its paste-ready list)
%   entries = { ...
%     'sim',      3,  3,   [4 9];        % pattern, n (pulses per electrode), condition number, bad trials
%     'sim',      4,  4,   [];           % [] = reviewed, NO bad trials
%     'seq B->A', 2,  17,  [1 6 12];
%   };
%   - the CONDITION NUMBER is the key; pattern and n are checked against it (error if they do not agree)
%   - trial numbers = trial number WITHIN the condition (1 ... repetitions), as printed on the check figures
%   - every condition is independent: a bad trial removes the whole trial from THAT condition only
%   - a condition without a row = NOT reviewed (nothing is removed there)
%   - one row per condition; a list REPLACES whatever was saved for that condition before
%
% MODE
%   'replace' : the file is written from the entries typed here. An existing file is first copied to a dated backup,
%               then replaced (conditions that are not typed here become "not reviewed").            <- default
%   'new'     : refuses to run when the file already exists.
%
% OUTPUT FILE  <data_folder>/<folder name>.FR_BadTrials.mat   with the struct  BadTrials:
%   BadTrials.badRel{cond}      bad trial numbers within the condition
%   BadTrials.badAbs{cond}      the same trials as absolute trial numbers in the recording (for reference)
%   BadTrials.isReviewed(cond)  true if the condition has a row (an empty list means reviewed, none bad)
%   BadTrials.nTrials(cond)     repetitions of the condition
%   BadTrials.entries, rec_file, date, version
%   The list is the same for the unthinned and every thinned recovery file (same trial numbers).
%   (The older <dataset>.BadTrials.mat in the data folder is a different file and is not touched.)
%   Use in an analysis script:   S = load(file, 'BadTrials');  bad = S.BadTrials.badRel{cond};
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/CJ287/Xia_Train_15Pulses_3_100Hz_30uA';
rec_tag     = '';          % '' = *.sp_xia_TrainRecovery.mat (only used to read the condition table and trial list)
rec_file    = '';          % full path (overrides rec_tag)

% paste the rows printed by the bad-trial check here, then edit them
entries = {  
    'sim',          1,    1,   [11 12 13];
    'sim',          2,    2,   [7 8 11];
    'sim',          3,    3,   [5 8 10 11];
    'sim',          4,    4,   [6 13];
    'sim',          5,    5,   [2 3 9 13];
    'sim',          6,    6,   [6];
    'sim',          7,    7,   [12 15];
    'sim',          8,    8,   [7 8 10];
    'sim',          9,    9,   [2 7 9];
    'sim',         10,   10,   [];
    'sim',         11,   11,   [5 9 14];
    'sim',         12,   12,   [2 6 9 15];
    'sim',         13,   13,   [2 15];
    'sim',         14,   14,   [7];
    'sim',         15,   15,   [4 9];
    'seq B->A',     1,   17,   [2 14];
    'seq B->A',     2,   19,   [11 13 14 15];
    'seq B->A',     3,   21,   [];
    'seq B->A',     4,   23,   [];
    'seq B->A',     5,   25,   [5 7];
    'seq B->A',     6,   27,   [2 10 14];
    'seq B->A',     7,   29,   [11];
    'seq B->A',     8,   31,   [];
    'seq B->A',     9,   33,   [10 11];
    'seq B->A',    10,   35,   [9 11 14];
    'seq B->A',    11,   37,   [1 2 4 5 6];
    'seq B->A',    12,   39,   [9 14 15];
    'seq B->A',    13,   41,   [7];
    'seq B->A',    14,   43,   [4 11];
    'seq B->A',    15,   45,   [15];
    };

mode      = 'replace';     % 'replace' | 'new'
out_file  = '';            % '' = <data_folder>/<folder name>.FR_BadTrials.mat
do_save   = true;          % false = only check and print
Warn_Exclusion_Fraction = 0.30;   % warn when a condition loses more than this share of its trials
Warn_MinTrialsLeft      = 5;      % warn when fewer trials than this remain

%% ================= LOAD CONDITION TABLE =================
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
assert(ismember(mode, {'replace','new'}), 'mode: replace or new.');
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
nTr = arrayfun(@(k) sum(trialCond == k), 1:nCond);

fe = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));  fe = fe(~startsWith({fe.name}, '.'));
assert(~isempty(fe), 'No *_exp_datafile_*.mat in %s.', data_folder);
E = load(fullfile(data_folder, fe(1).name), 'StimParams', 'simultaneous_stim');
nEl = double(E.simultaneous_stim);
peK = zeros(nCond, nEl);
for k = 1:nCond
    t1 = find(trialCond == k, 1);
    if isempty(t1), continue; end
    for e = 1:nEl
        row = 1 + (t1-1)*nEl + e;
        if double(E.StimParams{row,16}) > 0, peK(k,e) = double(E.StimParams{row,8}); end
    end
end
clear E

%% ================= PATTERNS =================
kinds = arrayfun(@(c) char(c.kind), condTable, 'UniformOutput', false);
isFinal = arrayfun(@(c) isfield(c,'isFinal') && ~isempty(c.isFinal) && logical(c.isFinal), condTable);
patCond = find(isFinal(:)' & ismember(kinds(:)', {'sim','AB','BA'}));
assert(~isempty(patCond), 'No sim / AB / BA pattern found in the recovery file.');
nPat = numel(patCond);
patName = cell(1, nPat);  patChain = cell(1, nPat);
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
end

%% ================= CHECK ENTRIES =================
if isempty(entries), entries = cell(0, 4); end
assert(iscell(entries) && size(entries,2) == 4, 'entries must be a cell array with 4 columns: pattern, n, condition, trials.');
nE = size(entries, 1);
badRel = cell(nCond, 1);  badAbs = cell(nCond, 1);  isReviewed = false(nCond, 1);
for k = 1:nCond, badRel{k} = zeros(1, 0);  badAbs{k} = zeros(1, 0); end
rowOf = zeros(nCond, 1);
for r = 1:nE
    pn = entries{r,1};  nn = entries{r,2};  cond = entries{r,3};  tl = entries{r,4};
    assert(ischar(pn) || isstring(pn), 'entries row %d: the pattern must be text.', r);
    pn = char(pn);
    p = find(strcmp(patName, pn), 1);
    assert(~isempty(p), 'entries row %d: pattern ''%s'' not found. Available: %s', r, pn, strjoin(patName, ', '));
    assert(isnumeric(cond) && isscalar(cond) && cond == round(cond) && cond >= 1 && cond <= nCond, ...
        'entries row %d: the condition must be a whole number between 1 and %d.', r, nCond);
    assert(ismember(cond, patChain{p}), 'entries row %d: condition %d is not a level of ''%s'' (its conditions: %s).', r, cond, pn, fmtList_(patChain{p}));
    assert(isnumeric(nn) && isscalar(nn) && nn == max(peK(cond,:)), ...
        'entries row %d: n = %s does not agree with condition %d (%d pulse(s) per electrode at most).', r, num2str(nn), cond, max(peK(cond,:)));
    assert(rowOf(cond) == 0, 'entries row %d: condition %d already has a row (row %d). One row per condition.', r, cond, rowOf(cond));
    assert(isnumeric(tl), 'entries row %d: the trials must be a numeric list (use [] for none).', r);
    tl = unique(double(tl(:)'));
    assert(all(tl == round(tl)) && all(tl >= 1 & tl <= nTr(cond)), ...
        'entries row %d (condition %d): trial numbers must be whole numbers between 1 and %d (found %s).', r, cond, nTr(cond), mat2str(tl(~(tl >= 1 & tl <= nTr(cond)))));
    rowOf(cond) = r;  isReviewed(cond) = true;
    badRel{cond} = tl;
    all_ = find(trialCond == cond);
    badAbs{cond} = all_(tl)';
end

%% ================= PRINT =================
fprintf('\n============================================================\n');
fprintf('NN TRAIN FR BAD-TRIAL ENTRY (v1)\n');
fprintf('============================================================\n');
fprintf('Folder   : %s\nCond. table from : %s\n', data_folder, rec_file);
fprintf('Mode     : %s | entries typed: %d | conditions reviewed: %d of %d\n', mode, nE, sum(isReviewed), nCond);
nWarn = 0;
for p = 1:nPat
    fprintf('\n--- %s ---\n', upper(patName{p}));
    fprintf('  cond   n   reps   state           bad  left   rel trials (abs trials)\n');
    for k = patChain{p}
        nb = numel(badRel{k});  left = nTr(k) - nb;
        if ~isReviewed(k)
            st = 'not reviewed';  txt = '';
        elseif nb == 0
            st = 'reviewed';  txt = '[]';
        else
            st = 'reviewed';  txt = sprintf('%s (%s)', fmtList_(badRel{k}), fmtList_(badAbs{k}));
        end
        fl = '';
        if isReviewed(k) && nb/max(nTr(k),1) > Warn_Exclusion_Fraction, fl = sprintf('  WARNING: %.0f%% removed', 100*nb/nTr(k));  nWarn = nWarn + 1; end
        if isReviewed(k) && left < Warn_MinTrialsLeft, fl = [fl sprintf('  WARNING: only %d trials left', left)];  nWarn = nWarn + 1; end %#ok<AGROW>
        fprintf('  %3d  %2d   %3d    %-13s   %3d   %3d   %s%s\n', k, max(peK(k,:)), nTr(k), st, nb, left, txt, fl);
    end
end
fprintf('\nTotal bad trials: %d in %d reviewed conditions. Warnings: %d.\n', sum(cellfun(@numel, badRel)), sum(isReviewed), nWarn);

%% ================= SAVE =================
[~, folderName] = fileparts(regexprep(data_folder, '[\\/]+$', ''));
if isempty(out_file), out_file = fullfile(data_folder, [folderName '.FR_BadTrials.mat']); end
BadTrials = struct('badRel', {badRel}, 'badAbs', {badAbs}, 'isReviewed', isReviewed, 'nTrials', nTr(:), 'entries', {entries}, ...
    'rec_file', rec_file, 'date', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'version', 'NN_Train_FR_BadTrialEntry_v1');
fprintf('\nOutput file : %s\n', out_file);
if ~do_save
    fprintf('do_save = false: nothing saved.\n');
elseif isfile(out_file) && strcmp(mode, 'new')
    error('The file already exists and mode = ''new''. Nothing was changed. Use mode = ''replace'' to replace it (a dated backup is made first).');
else
    if isfile(out_file)
        bak = strrep(out_file, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
        copyfile(out_file, bak);
        fprintf('Old file copied to %s and replaced.\n', bak);
    end
    save(out_file, 'BadTrials');
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
