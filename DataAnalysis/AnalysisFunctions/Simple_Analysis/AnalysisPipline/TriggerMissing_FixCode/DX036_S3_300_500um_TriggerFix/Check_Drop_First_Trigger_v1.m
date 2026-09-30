%% CHECK_DROPPED_FIRST_TRIGGER
%
% cleanTrig_sabquick.m unconditionally deletes the very first detected
% pulse if trig(1) < 30000 samples (i.e. within the first ~1 second),
% assuming it's spurious pre-experiment noise:
%
%   if trig(1)<30000
%       trig(1)=[];
%   end
%
% If that first pulse was actually trial 1's real trigger, deleting it
% shifts every later trial's index down by one for the rest of the file --
% which would explain the "pause lands one index early" observation
% without needing any anomaly earlier in the recording.
%
% This redoes cleanTrig_sabquick's OWN detection method (exact
% digital_in==1 match, same de-dup trick) WITHOUT the trig(1)<30000
% deletion, so we can see exactly what that discarded first pulse was,
% and how many pulses this gives in total.
%
% Diagnostic only -- does not touch your .trig.dat.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_300_500um_SimSeq1';
cd(data_folder);

[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;

%% ====================== REPLICATE cleanTrig_sabquick's OWN DETECTION ======================
fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);

stimDig = flip(find(digital_in == 1));
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
stimDig = flip(stimDig);   % chronological order, BEFORE the trig(1)<30000 deletion

fprintf('Total pulses detected (before trig(1)<30000 deletion): %d\n', numel(stimDig));
fprintf('First pulse: sample %d (t = %.4f s)\n', stimDig(1), stimDig(1)/FS);
fprintf('Second pulse: sample %d (t = %.4f s)\n', stimDig(2), stimDig(2)/FS);
fprintf('Gap between 1st and 2nd pulse: %.1f ms\n', (stimDig(2)-stimDig(1))/FS*1000);

fprintf(['\nWould cleanTrig_sabquick delete this first pulse? %s (threshold: sample < 30000, i.e. t < %.3f s)\n'], ...
    string(stimDig(1) < 30000), 30000/FS);

fprintf(['\nFor reference, RunMultiElect_exp.m does pause(20) before the loop starts, then\n' ...
    'begins trial 1 fairly promptly. If stimDig(1) sits somewhere around t ~ 20-21 s\n' ...
    '(sample ~600000-630000), it is almost certainly the genuine trial-1 trigger, NOT\n' ...
    'spurious noise near t=0 -- meaning trig(1)<30000 (t<1s) would only ever fire on truly\n' ...
    'spurious pre-recording noise, not this. If instead stimDig(1) really is within the\n' ...
    'first second (t < 1s), it likely IS spurious (e.g. a DAQ reset blip) and the deletion\n' ...
    'is correct, and this theory does not explain the observed shift.\n']);

%% ====================== TOTAL COUNT WITHOUT THE DELETION ======================
fprintf('\nIf this pulse is genuinely trial 1''s trigger and should NOT be deleted,\n');
fprintf('total pulse count (undeleted) = %d, vs. n_Trials = 810.\n', numel(stimDig));
if numel(stimDig) == 810
    fprintf(['This matches exactly -- strong support that keeping this pulse\n' ...
        '(i.e. not applying the trig(1)<30000 deletion) fully accounts for all\n' ...
        '810 trials with no missing trigger anywhere, including at the pause.\n']);
else
    fprintf('This does NOT match 810 on its own -- there may still be another separate issue.\n');
end