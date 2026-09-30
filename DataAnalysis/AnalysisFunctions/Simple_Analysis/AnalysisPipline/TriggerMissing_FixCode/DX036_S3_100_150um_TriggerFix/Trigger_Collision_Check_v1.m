%% CHECK_VIS_COLLISION
%
% Tests whether the 69 "missing" pulses are actually real stim events that
% the exact-match detector (digital_in == 1) missed because a separate
% "vis" digital line (bit 1) happened to be high at the same sample,
% making digital_in read 3 instead of 1.
%
% Robust detection uses bitget(digital_in,1) -- true whenever the stim bit
% is set, regardless of any other bit -- and compares the count/timing
% against the fragile exact-match method.
%
% Diagnostic only -- does not modify any files.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_100_150um_Single1_260929_144315';
cd(data_folder);

[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;

fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);

fprintf('Unique digital_in values present: %s\n', mat2str(unique(digital_in)'));

%% ====================== EXACT-MATCH METHOD (current pipeline) ======================
stimDig_exact = flip(find(digital_in == 1));
dt = diff(stimDig_exact);
kill = dt == -1;
stimDig_exact(kill) = [];
stimDig_exact = flip(stimDig_exact);
fprintf('Exact-match (digital_in==1) pulse count: %d\n', numel(stimDig_exact));

%% ====================== ROBUST METHOD (bit 0, regardless of other bits) ======================
stim_bit = bitget(digital_in, 1);
rising = find(diff([0; stim_bit]) == 1);   % rising edges of bit 0
fprintf('Robust (bitget bit0) pulse count: %d  (expected 270)\n', numel(rising));

%% ====================== HOW MANY OF THE ROBUST PULSES DID EXACT-MATCH MISS? ======================
missed_by_exact = setdiff(rising, stimDig_exact);
fprintf('\nPulses found by the robust method but NOT by exact-match: %d\n', numel(missed_by_exact));

if ~isempty(missed_by_exact)
    fprintf('digital_in value at each missed pulse (should be 3 if the vis-collision theory is right):\n');
    vals_at_missed = digital_in(missed_by_exact);
    tabulate_vals = unique(vals_at_missed);
    for v = tabulate_vals'
        fprintf('  digital_in == %d at %d of the missed pulses\n', v, sum(vals_at_missed==v));
    end
end

%% ====================== SANITY: DOES ROBUST COUNT MATCH THE JITTER PATTERN TOO? ======================
iti_ms = diff(rising)/FS*1000;
fprintf('\nRobust-method ITI: median %.2f ms, MAD %.2f ms\n', median(iti_ms), median(abs(iti_ms-median(iti_ms))));