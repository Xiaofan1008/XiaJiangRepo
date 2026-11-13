%% MISSING_TRIG_3_APPLY
%
% Step 3 of 3 in the missing-trigger repair workflow:
%   1) missing_trig_1_diagnose.m   (recount pulses + locate gaps)
%   2) missing_trig_2_recover.m    (work out which trials are missing)
%   3) missing_trig_3_apply.m      <-- this script (back up and rewrite .trig.dat)
%
% THIS IS THE ONLY STEP THAT CHANGES A FILE.
%
% What it does:
%   1) Refuses to run unless step 2 says PASS and was run after the latest
%      step 1, and unless the current .trig.dat is still exactly the file
%      step 1 recorded (so a file can never be patched twice by accident).
%   2) Copies the original .trig.dat into missing_trig_backup/ under a name
%      that does NOT end in .trig.dat (loadTrig would otherwise find two).
%   3) Overwrites .trig.dat in place with step 2's proposed trig (double,
%      no header -- same format cleanTrig_sabquick writes), missing trials
%      as -500, then reads it back with loadTrig to confirm.
%   4) Saves trig_info.mat (per-trial table incl. continuous recording
%      after each stimulus) and missing_trig_summary.txt in the folder,
%      and records everything in report.apply.
%
% To undo: copy the backup file back over the .trig.dat (instructions are
% printed at the end and written to the summary).

clear; clc;

%% ======================= SETTINGS =======================
data_folder    = '/Volumes/MACData/Data/Data_Xia/DX039/Xia_Linearity_S3_100_400um_SimSeq1_261009_190544';   % <-- same folder as steps 1 and 2
save_trig_info = true;    % write trig_info.mat with per-trial info for analysis scripts

%% ======================= CHECKS =======================
assert(isfolder(data_folder), 'data_folder not found:\n%s', data_folder);
cd(data_folder);
report_file = fullfile(data_folder, 'missing_trig_report.mat');
assert(isfile(report_file), 'No missing_trig_report.mat here -- run steps 1 and 2 first.');
load(report_file, 'report');
assert(isfield(report, 'diagnose'), 'No step 1 results in the report.');
assert(isfield(report, 'recover'),  'No step 2 results in the report -- run missing_trig_2_recover.m first.');
info = report.info;  dg = report.diagnose;  rc = report.recover;

assert(isfield(rc, 'version') && rc.version >= 2, ...
    'Step 2 results are from the old version -- rerun the new missing_trig_2_recover.m first.');
assert(isequal(rc.diagnose_run_time, dg.run_time), ...
    'Step 1 was rerun after step 2 -- rerun missing_trig_2_recover.m first.');
assert(strcmp(rc.status, 'pass'), ...
    'Step 2 status is "%s", not PASS -- nothing written.', rc.status);
if isfield(report, 'apply')
    error(['This folder was already patched on %s (backup: %s). Nothing written.\n' ...
           'To redo it, restore the backup and run steps 1-3 again.'], ...
          char(report.apply.run_time), report.apply.backup_file);
end

trigDIR = dir(fullfile(data_folder, '*.trig.dat'));
assert(numel(trigDIR) == 1, 'Expected exactly one *.trig.dat in the folder, found %d.', numel(trigDIR));
trig_name = trigDIR(1).name;
trig_path = fullfile(data_folder, trig_name);
assert(strcmp(trig_name, info.trig_file), ...
    'The .trig.dat is now "%s" but step 1 recorded "%s".', trig_name, info.trig_file);

fid = fopen(trig_path, 'r');
assert(fid > 0, 'Could not open %s', trig_name);
current = fread(fid, inf, 'double');
fclose(fid);
assert(isequal(current(:), double(info.old_trig(:))), ...
    '%s has changed since step 1 ran -- rerun steps 1 and 2 before applying.', trig_name);

new_trig = double(rc.proposed_trig(:))';
n_exp    = info.n_expected;
assert(numel(new_trig) == n_exp, 'Proposed trig has %d entries, expected %d.', numel(new_trig), n_exp);
assert(~any(isnan(new_trig)), 'Proposed trig still has unresolved trials.');
assert(all(new_trig == -500 | (new_trig > 0 & new_trig <= info.nSam)), 'Proposed trig has invalid values.');
pos = new_trig(new_trig > 0);
assert(all(diff(pos) > 0), 'Proposed trig times are not in increasing order.');

missing_trials = find(new_trig == -500);
n_changed = sum(current(:)' ~= new_trig);
fprintf('Folder: %s\n', data_folder);
fprintf('%s: %d entries, %d will change (%d -> %d placeholders)\n', trig_name, n_exp, n_changed, ...
    sum(current == -500), numel(missing_trials));
fprintf('Missing trials (-500): %s\n\n', mat2str(missing_trials));

%% ======================= BACKUP =======================
backup_dir = fullfile(data_folder, 'missing_trig_backup');
if ~isfolder(backup_dir), mkdir(backup_dir); end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
backup_file = fullfile(backup_dir, sprintf('%s.original_%s.bak', trig_name, stamp));
[ok, msg] = copyfile(trig_path, backup_file);
assert(ok, 'Backup failed: %s', msg);
fid = fopen(backup_file, 'r');  bk = fread(fid, inf, 'double');  fclose(fid);
assert(isequal(bk(:), current(:)), 'Backup does not match the original -- nothing written.');
fprintf('Backup:  %s\n', backup_file);

%% ======================= WRITE + VERIFY =======================
fid = fopen(trig_path, 'w');
assert(fid > 0, 'Could not open %s for writing (backup is safe).', trig_name);
cnt = fwrite(fid, new_trig, 'double');
fclose(fid);
assert(cnt == n_exp, 'Only %d of %d values written -- restore from the backup.', cnt, n_exp);

check = loadTrig(0);
assert(isequal(double(check(:)), new_trig(:)), ...
    'Read-back with loadTrig does not match -- restore from the backup:\n  %s', backup_file);
assert(numel(dir(fullfile(data_folder, '*.trig.dat'))) == 1, 'More than one *.trig.dat now in the folder.');
fprintf('Written: %s (verified with loadTrig)\n', trig_path);

%% ======================= trig_info.mat =======================
trial    = (1:n_exp)';
trig_col = new_trig(:);
is_missing = trig_col == -500;
trig_index = nan(n_exp, 1);
known = ~isnan(rc.trial_of);
trig_index(rc.trial_of(known)) = find(known);
amp = info.trialAmps(:);  if numel(amp) < n_exp, amp(end+1:n_exp) = NaN; end
cont_ms = rc.cont_ms(:);
trial_table = table(trial, trig_col, is_missing, trig_index, amp(1:n_exp), cont_ms, ...
    'VariableNames', {'trial', 'trig', 'missing', 'detected_trigger', 'amp', 'cont_ms'});

c = cont_ms(~isnan(cont_ms));
if save_trig_info
    trig_info = struct( ...
        'created', datetime('now'), 'data_folder', data_folder, 'trig_file', trig_name, ...
        'n_trials', n_exp, 'missing_trials', missing_trials, ...
        'trials', trial_table, ...
        'note', ['cont_ms = continuous recording after each stimulus before the Intan restarts ' ...
                 'for the next trial. Analysis windows longer than this cross a restart.']);
    save(fullfile(data_folder, 'trig_info.mat'), 'trig_info');
    fprintf('Saved:   trig_info.mat (per-trial table; cont_ms median %.0f ms, min %.0f ms)\n', median(c), min(c));
end

%% ======================= SUMMARY TEXT =======================
sum_file = fullfile(data_folder, 'missing_trig_summary.txt');
fid = fopen(sum_file, 'w');
fprintf(fid, 'Missing-trigger repair\n');
fprintf(fid, 'Date:        %s\n', char(datetime('now')));
fprintf(fid, 'Folder:      %s\n', data_folder);
fprintf(fid, 'Trig file:   %s (overwritten)\n', trig_name);
fprintf(fid, 'Backup:      %s\n\n', backup_file);
fprintf(fid, 'Expected trials: %d, detected triggers: %d, missing: %d\n', n_exp, dg.n_found, numel(missing_trials));
fprintf(fid, 'Method: one time.dat counter reset per trial, solved between anchors (pauses after trial 499/999/..., short final pulse)\n');
fprintf(fid, 'Trig convention: raw %s edge + %d samples (same as the original file)\n\n', dg.edge_used, rc.trig_offset);
fprintf(fid, 'Missing trials (written as -500; TTL never fired, so not stimulated):\n');
for i = 1:height(rc.missing)
    fprintf(fid, '  trial %d: between trig %d and %d, restart at %.3f s, amp %g, %s\n', ...
        rc.missing.trial(i), rc.missing.after_trig(i), rc.missing.before_trig(i), ...
        rc.missing.reset_time_s(i), rc.missing.amp(i), rc.missing.electrodes{i});
end
if ~isempty(c)
    fprintf(fid, '\nContinuous recording after each stimulus: median %.0f ms, 5th percentile %.0f ms, min %.0f ms\n', ...
        median(c), prctile(c, 5), min(c));
end
fprintf(fid, '\nTo undo: copy the backup over %s, then delete missing_trig_report.mat.\n', trig_name);
fclose(fid);
fprintf('Saved:   missing_trig_summary.txt\n');

%% ======================= RECORD IN REPORT =======================
report.apply = struct( ...
    'run_time', datetime('now'), 'recover_run_time', rc.run_time, ...
    'trig_file', trig_path, 'backup_file', backup_file, ...
    'old_trig', current(:)', 'new_trig', new_trig, ...
    'missing_trials', missing_trials, 'n_changed', n_changed, ...
    'trig_info_saved', save_trig_info);
save(report_file, 'report');

fprintf('\nDone. Downstream scripts using loadTrig(0) now get the corrected trig.\n');
fprintf('To undo:\n  copyfile(''%s'', ''%s'');\n', backup_file, trig_path);