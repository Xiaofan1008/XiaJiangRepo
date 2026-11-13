%% PATCH_TRIG_THREE_GAPS
%
% This dataset expects 270 trials (n_Trials=270, simultaneous_stim=1) but
% only 201 real triggers were ever detected on the digital line -- no
% detection-method bug (bit-collision ruled out) and no periodic artifact
% evidence in the quiet stretches (ruled out recoverable-but-unlogged
% delivery), so the 69 missing triggers are treated as genuinely
% unrecoverable and filled with -500 sentinels, at BEST-ESTIMATE
% positions (gap duration / typical ITI, rounded to the nearest trial):
%
%   ~24 sentinels BEFORE the first real pulse   (9.69s quiet / 397ms ITI)
%   ~21 sentinels IN THE MIDDLE, between real pulse 62 and 63 (8895.8ms gap)
%   ~24 sentinels AFTER the last real pulse     (9.42s quiet / 397ms ITI)
%   24 + 21 + 24 = 69, matching the full shortfall (270 - 201) exactly.
%
% These counts are ESTIMATES, not exact recoveries -- real inter-trial
% timing jitters by ~30ms (MAD) around the 397ms median, so each count
% could plausibly be off by a trial or two. There is no independent log
% to pin these down further, per discussion.
%
% Run this ONCE from inside the dataset folder. Backs up the original
% .trig.dat before writing anything.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_100_150um_Single1_260929_144315';
cd(data_folder);

n_Trials_expected = 270;
n_before = 24;
n_mid    = 21;
n_after  = 24;
gap_after_pulse = 62;   % the mid-recording gap sits between real pulse 62 and 63

%% ====================== LOCATE + BACK UP THE FILE ======================
tmp = dir('*.trig.dat');
assert(~isempty(tmp), 'No *.trig.dat file found in %s', data_folder);
trig_path = fullfile(data_folder, tmp.name);

backup_path = [trig_path '.bak_pre_threegapfix'];
if isfile(backup_path)
    fprintf('Backup already exists at %s -- not overwriting it.\n', backup_path);
else
    copyfile(trig_path, backup_path);
    fprintf('Backed up original file to: %s\n', backup_path);
end

%% ====================== READ CURRENT TRIG, KEEP ONLY REAL VALUES ======================
true_count = tmp.bytes / 8;   % double = 8 bytes
fid = fopen(trig_path, 'r');
trig_old = fread(fid, [1, true_count], 'double');
fclose(fid);

fprintf('Read %d values from %s.\n', numel(trig_old), tmp.name);

real_pulses = trig_old(trig_old ~= -500);
fprintf('Of those, %d are real (non-sentinel) pulse times.\n', numel(real_pulses));

assert(numel(real_pulses) == 201, ...
    'Expected 201 real pulses but found %d -- STOPPING, check before proceeding.', numel(real_pulses));
assert(gap_after_pulse < numel(real_pulses), 'gap_after_pulse out of range.');

%% ====================== BUILD THE PATCHED ARRAY ======================
before_block = real_pulses(1:gap_after_pulse);
after_block  = real_pulses(gap_after_pulse+1:end);

trig_new = [ ...
    -500*ones(1,n_before), ...
    before_block, ...
    -500*ones(1,n_mid), ...
    after_block, ...
    -500*ones(1,n_after) ...
    ];

fprintf('\nPatched array length: %d (expected %d)\n', numel(trig_new), n_Trials_expected);
assert(numel(trig_new) == n_Trials_expected, ...
    'Patched length %d does not match expected n_Trials %d -- STOPPING, check the sentinel counts.', ...
    numel(trig_new), n_Trials_expected);

fprintf('trig_new(1:5)              = %s   (all -500, before-block sentinels)\n', mat2str(trig_new(1:5)));
fprintf('trig_new(%d:%d)          = %s ... (start of real before-block)\n', n_before+1, n_before+3, mat2str(trig_new(n_before+1:n_before+3)));
fprintf('trig_new(%d)               = %g   (last real pulse before the mid-gap)\n', n_before+gap_after_pulse, trig_new(n_before+gap_after_pulse));
fprintf('trig_new(%d:%d)          = %s   (mid-gap sentinels)\n', n_before+gap_after_pulse+1, n_before+gap_after_pulse+3, mat2str(trig_new(n_before+gap_after_pulse+1:n_before+gap_after_pulse+3)));
fprintf('trig_new(end-%d:end)      = %s   (all -500, after-block sentinels)\n', n_after-1, mat2str(trig_new(end-n_after+1:end)));

%% ====================== WRITE BACK OUT ======================
fid = fopen(trig_path, 'w');
fwrite(fid, trig_new, 'double');
fclose(fid);

fprintf('\nWrote %d values back to %s.\n', numel(trig_new), trig_path);

trig_check = loadTrig(0);
fprintf('Verifying: loadTrig(0) now returns %d values.\n', numel(trig_check));
assert(isequal(trig_check, trig_new), 'Mismatch after re-read -- restore from the .bak file.');
fprintf('OK: trig(tr) now approximates trial tr for tr = 1..%d (sentinel blocks at 1:%d, %d:%d, and the last %d).\n', ...
    n_Trials_expected, n_before, n_before+gap_after_pulse+1, n_before+gap_after_pulse+n_mid, n_after);