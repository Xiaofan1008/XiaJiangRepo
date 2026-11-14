%% ============================================================
% NP_Train_FR_BadTrialEntry_v1.m
% NP probe | pulse TRAIN sessions | FIRING-RATE analysis | BAD-TRIAL ENTRY | v1   (NP version of NN_Train_FR_BadTrialEntry_v1)
%
% PURPOSE
%   You type the FINAL bad trials of every condition you reviewed with NP_Train_FR_BadTrialCheck_v1; the script checks
%   the entries, prints a summary and saves ONE file that the FR analysis scripts read (NP_Train_FR_RateVsPulses ...).
%   The recording, the spike files and all other files are never changed.
%
% ENTRY FORMAT  (the same rows the bad-trial check prints as its PASTE-READY list)
%   entries = { ...
%     'sim',        3,  3,  [4 9];        % pattern, n (pulses per electrode), condition number, bad trials
%     'seq B->A',   2, 17,  [];           % [] = reviewed, NO bad trials
%     'visual',     0, 49,  [1];          % visual conditions: n = 0
%   };
%   - the CONDITION NUMBER is the key; pattern and n are checked against it (error if they do not agree)
%   - pattern names: 'sim', 'seq A->B', 'seq B->A', 'visual'; half-steps of seq have n = x.5 (e.g. 2.5 for [3 2])
%   - trial numbers = trial number WITHIN the condition (1 ... repetitions), as on the check figures
%   - a bad trial removes the whole trial (all channels) from THAT condition only
%   - a condition without a row = NOT reviewed (all its trials are kept by the analyses)
%   - one row per condition
%
% MODE
%   'replace' (default): the file is written from the entries typed here. An existing file is first copied to a dated
%             backup, then replaced (conditions that are not typed here become "not reviewed").
%   'new'   : refuses to run when the file already exists.
%
% OUTPUT  <recording1>/<short>.NP_BadTrials.mat   with the struct BadTrials:
%   BadTrials.badRel{cond}      bad trial numbers within the condition   <- read by the FR scripts
%   BadTrials.badAbs{cond}      the same trials as rows of TT.trial (recording order), for reference
%   BadTrials.isReviewed(cond)  true if the condition has a row ([] = reviewed, none bad)
%   BadTrials.nTrials(cond)     repetitions of the condition
%   BadTrials.entries, note, ttFile, date, version
%   The list is valid for every spike file of the session (original, recovered, RandomEvent): the trials are the same.
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';

% PLACEHOLDER rows: paste the list printed by NP_Train_FR_BadTrialCheck_v1 here and edit it.
entries = { ...
    'sim',           1,   1, [7,12 14]; ...
    'sim',           2,   2, [8 15]; ...
    'sim',           3,   3, [1 4 7 10 14]; ...
    'sim',           4,   4, 6; ...
    'sim',           5,   5, [2 3 5 9]; ...
    'sim',           6,   6, [2 9 10 11]; ...
    'sim',           7,   7, [5 10 15]; ...
    'sim',           8,   8, [2 10 5 8]; ...
    'sim',           9,   9, [2 15]; ...
    'sim',          10,  10, [1 11 12 15]; ...
    'sim',          11,  11, [4 9 10 ]; ...
    'sim',          12,  12, [2 8 10 13 15]; ...
    'sim',          13,  13, []; ...
    'sim',          14,  14, 1; ...
    'sim',          15,  15, [4 15]; ...
    'seq B->A',      1,  17, [4,6 12,13 14]; ...
    'seq B->A',      2,  19, [4 11 13 14 15]; ...
    'seq B->A',      3,  21, [5 7 8 13]; ...
    'seq B->A',      4,  23, [7 15]; ...
    'seq B->A',      5,  25, [1 7 11]; ...
    'seq B->A',      6,  27, [9 14]; ...
    'seq B->A',      7,  29, [4 14]; ...
    'seq B->A',      8,  31, [8 11]; ...
    'seq B->A',      9,  33, [2 7]; ...
    'seq B->A',     10,  35, [6 9 11 13 14]; ...
    'seq B->A',     11,  37, [1 2 4 5 6]; ...
    'seq B->A',     12,  39, [3 7 8 9 10 15]; ...
    'seq B->A',     13,  41, [5]; ...
    'seq B->A',     14,  43, [10]; ...
    'seq B->A',     15,  45, [];
    };

note    = '';         % free text saved in the file, e.g. 'BadTrialCheck v1, k = 5 for R/M/B, baseline -100..-5 ms'
mode    = 'replace';  % 'replace' (default) | 'new'
do_save = true;       % false = only check and print
%% =================================================

codeVersion = 'NP_Train_FR_BadTrialEntry_v1 (2026-11-14)';
assert(any(strcmp(mode, {'replace', 'new'})), 'mode must be ''replace'' or ''new''.');

%% ---- trial table ----------------------------------------------------------------------------------------
tf = dir(fullfile(sessionDir, '**', '*.NP_TrialTable.mat'));
tf = tf(~startsWith({tf.name}, '._') & ~contains({tf.name}, '_backup_'));
assert(numel(tf) == 1, 'Found %d trial tables below %s (need exactly 1).', numel(tf), sessionDir);
recDir = tf(1).folder;  short = erase(tf(1).name, '.NP_TrialTable.mat');
ttFile = fullfile(recDir, tf(1).name);
L = load(ttFile, 'TT');  TT = L.TT;  clear L
nCond = numel(TT.cond);  trialCond = TT.trial.cond(:);

% pattern name and n of every condition (as the bad-trial check prints them)
patOf = cell(nCond, 1);  nOf = nan(nCond, 1);
for c = 1:nCond
    C = TT.cond(c);
    switch C.kind
        case 'sim',    patOf{c} = 'sim';
        case 'AB',     patOf{c} = 'seq A->B';
        case 'BA',     patOf{c} = 'seq B->A';
        case 'visual', patOf{c} = 'visual';  nOf(c) = 0;
        otherwise,     patOf{c} = C.kind;
    end
    if any(strcmp(C.kind, {'sim', 'AB', 'BA'}))
        nsp = C.nsp;  nsp(C.amp <= 0) = 0;
        if strcmp(C.kind, 'sim') || nsp(1) == nsp(2), nOf(c) = max(nsp); else, nOf(c) = max(nsp) - 0.5; end
    end
end

%% ---- check the entries ------------------------------------------------------------------------------------
assert(iscell(entries) && size(entries, 2) == 4, 'entries must be a cell array with 4 columns: pattern, n, cond, trials.');
nE = size(entries, 1);
badRel = cell(nCond, 1);  badAbs = cell(nCond, 1);  isReviewed = false(nCond, 1);  nTrials = zeros(nCond, 1);
for c = 1:nCond, nTrials(c) = sum(trialCond == c); badRel{c} = zeros(1, 0); badAbs{c} = zeros(1, 0); end
rowOf = zeros(nCond, 1);
for r = 1:nE
    pn = entries{r, 1};  n = entries{r, 2};  c = entries{r, 3};  bt = entries{r, 4};
    assert(ischar(pn) || isstring(pn), 'entries row %d: the pattern must be text.', r);  pn = char(pn);
    assert(isnumeric(c) && isscalar(c) && c == round(c) && c >= 1 && c <= nCond, 'entries row %d: condition must be 1..%d.', r, nCond);
    assert(strcmp(pn, patOf{c}), 'entries row %d: condition %d is ''%s'' (%s), not ''%s''.', r, c, patOf{c}, TT.cond(c).label, pn);
    assert(isnumeric(n) && isscalar(n) && abs(n - nOf(c)) < 1e-9, 'entries row %d: condition %d has n = %g, not %g.', r, c, nOf(c), n);
    assert(rowOf(c) == 0, 'entries row %d: condition %d is already entered in row %d.', r, c, rowOf(c));
    assert(isnumeric(bt), 'entries row %d: the bad trials must be a numeric list ([] = none).', r);
    bt = unique(double(bt(:)'));
    assert(all(bt == round(bt) & bt >= 1 & bt <= nTrials(c)), 'entries row %d (cond %d): trials must be 1..%d, got %s.', ...
           r, c, nTrials(c), mat2str(bt));
    rowOf(c) = r;  isReviewed(c) = true;
    ks = find(trialCond == c);
    badRel{c} = bt;  badAbs{c} = ks(bt)';
end

%% ---- summary --------------------------------------------------------------------------------------------------
fprintf('\n============================================================\nNP TRAIN FR BAD-TRIAL ENTRY (v1): %s\n', short);
fprintf('============================================================\n');
fprintf('  cond  pattern      n   trials  bad  label\n');
for c = find(isReviewed)'
    fprintf('  %4d  %-10s %4s  %6d  %3d  %-28s %s\n', c, patOf{c}, num2s_(nOf(c)), nTrials(c), numel(badRel{c}), ...
            TT.cond(c).label, ternary_(isempty(badRel{c}), '(none)', mat2str(badRel{c})));
end
pats = unique(patOf(isReviewed), 'stable');
fprintf('\nper pattern:\n');
for p = 1:numel(pats)
    cs = find(isReviewed & strcmp(patOf, pats{p}));
    nb = sum(cellfun(@numel, badRel(cs)));  nt = sum(nTrials(cs));
    fprintf('  %-10s %3d conditions reviewed, %4d of %4d trials bad (%.1f %%)\n', pats{p}, numel(cs), nb, nt, 100*nb/max(nt, 1));
end
notRev = find(~isReviewed & ismember(patOf, {'sim', 'seq A->B', 'seq B->A', 'visual'}) & nTrials > 0);
if ~isempty(notRev)
    fprintf('\nNOT reviewed (all trials kept): conds %s\n', mat2str(notRev'));
end

%% ---- save ----------------------------------------------------------------------------------------------------------
outFile = fullfile(recDir, [short '.NP_BadTrials.mat']);
fprintf('\nOutput: %s\n', outFile);
if ~do_save, fprintf('do_save = false: nothing saved.\n'); return; end
if isfile(outFile)
    if strcmp(mode, 'new'), error('%s already exists and mode = ''new''. Nothing was changed.', outFile); end
    bak = strrep(outFile, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
    copyfile(outFile, bak);  fprintf('Old file copied to %s\n', bak);
end
BadTrials = struct('badRel', {badRel}, 'badAbs', {badAbs}, 'isReviewed', isReviewed, 'nTrials', nTrials, ...
    'entries', {entries}, 'note', note, 'ttFile', ttFile, 'date', datestr(now, 'yyyy-mm-dd HH:MM:SS'), 'version', codeVersion);
save(outFile, 'BadTrials');
fprintf('Saved.\n');

%% ---- helpers --------------------------------------------------------------------------------------------------------
function s = num2s_(v)
if isnan(v), s = '-'; elseif v == round(v), s = sprintf('%d', v); else, s = sprintf('%.1f', v); end
end

function o = ternary_(c, a, b)
if c, o = a; else, o = b; end
end
