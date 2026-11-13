%% PATCH_TRIG_PREPEND_TRIAL1
%
% One-time fix for this specific dataset's .trig.dat:
%
%   Trial 1's real TTL pulse never registered on digitalin.dat (confirmed:
%   trigOffset=-1 works as a single global constant across the whole file,
%   with no doubled gap anywhere in trig(1:809), and no recoverable raw
%   artifact for trial 1 in amplifier.dat channels 5/15 -- see
%   recover_trial1_artifact.m). Everything from trial 2 onward is a real,
%   correctly-spaced pulse, just shifted one slot earlier than its true
%   trial number.
%
%   Fix: prepend a single -500 sentinel (the same "count restored, no
%   real value recoverable" convention cleanTrig_sabquick.m already uses
%   for its own catch-trial corrections) to the FRONT of trig. After this,
%   trig(tr) correctly indexes trial tr for every tr, and trial 1 can be
%   detected/skipped downstream via trig(tr)==-500 instead of a
%   trigOffset/trigIdx<1 hack.
%
% Run this ONCE from inside the dataset folder (same folder loadTrig/
% cleanTrig_sabquick already operate in). It backs up the original file
% before writing anything.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_300_500um_SimSeq1';
cd(data_folder);

%% ====================== LOCATE + BACK UP THE FILE ======================
tmp = dir('*.trig.dat');
assert(~isempty(tmp), 'No *.trig.dat file found in %s', data_folder);
trig_path = fullfile(data_folder, tmp.name);

backup_path = [trig_path '.bak_pre_trial1fix'];
if isfile(backup_path)
    fprintf('Backup already exists at %s -- not overwriting it. Delete it manually if you want a fresh backup.\n', backup_path);
else
    copyfile(trig_path, backup_path);
    fprintf('Backed up original file to: %s\n', backup_path);
end

%% ====================== READ CURRENT TRIG (TRUE LENGTH, NOT THE BUGGY sz) ======================
true_count = tmp.bytes / 8;   % double = 8 bytes; loadTrig.m's own sz (bytes/2) is wrong but harmless
fid = fopen(trig_path, 'r');
trig_old = fread(fid, [1, true_count], 'double');
fclose(fid);

fprintf('Read %d values from %s.\n', numel(trig_old), tmp.name);
fprintf('First 3: %s | Last 3: %s\n', mat2str(trig_old(1:3)), mat2str(trig_old(end-2:end)));

assert(numel(trig_old) == 810, ...
    'Expected 810 values (809 real + 1 trailing -500 sentinel) but got %d -- STOPPING, check before proceeding.', numel(trig_old));
assert(trig_old(end) == -500, ...
    'Expected the last value to be the known trailing -500 sentinel but got %g -- STOPPING, check before proceeding.', trig_old(end));

%% ====================== PREPEND THE SENTINEL ======================
trig_new = [-500, trig_old];

fprintf('\nNew array: %d values (was %d).\n', numel(trig_new), numel(trig_old));
fprintf('trig_new(1:5)  = %s   (trial 1 = sentinel, trial 2 = old trial 1''s pulse, ...)\n', mat2str(trig_new(1:5)));
fprintf('trig_new(end)  = %g   (unchanged trailing sentinel, now trial 811)\n', trig_new(end));

%% ====================== WRITE BACK OUT ======================
fid = fopen(trig_path, 'w');
fwrite(fid, trig_new, 'double');
fclose(fid);

fprintf('\nWrote %d values back to %s.\n', numel(trig_new), trig_path);
fprintf('Verifying by re-reading with loadTrig(0)...\n');

trig_check = loadTrig(0);
fprintf('loadTrig(0) now returns %d values.\n', numel(trig_check));
fprintf('trig_check(1:5) = %s\n', mat2str(trig_check(1:5)));
assert(isequal(trig_check, trig_new), 'Mismatch after re-read -- something is wrong, restore from the .bak file.');
fprintf('\nOK: trig(tr) now correctly indexes trial tr for tr = 1..811 (tr=1 and tr=811 are -500 sentinels).\n');