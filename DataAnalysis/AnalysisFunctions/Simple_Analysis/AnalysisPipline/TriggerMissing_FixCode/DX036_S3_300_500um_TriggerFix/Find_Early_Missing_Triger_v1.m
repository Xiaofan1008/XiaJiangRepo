%% FIND_EARLIER_MISSING_TRIGGER
%
% We've confirmed the pause-every-500-trials gap sits at array positions
% 498/499 instead of the expected 499/500 -- meaning the array is already
% short by one trial by that point, even though the trig(1)<30000
% deletion isn't the cause. This looks specifically at intervals BEFORE
% index 498, with a lower/more sensitive threshold, to find that earlier
% dropped pulse (a doubled gap that may not have been an extreme outlier
% given how jittered this recording's trial timing is).
%
% Diagnostic only.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_300_500um_SimSeq1';
cd(data_folder);

[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;

fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);

stimDig = flip(find(digital_in == 1));
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
stimDig = flip(stimDig);   % chronological, 809 entries, matches current trig(1:809)

fprintf('Total pulses: %d\n', numel(stimDig));

%% ====================== INTERVALS, RESTRICTED TO BEFORE INDEX 498 ======================
iti_samples = diff(stimDig);
iti_ms = iti_samples / FS * 1000;

region = 1:497;   % intervals between pulse(1)-pulse(2) up to pulse(497)-pulse(498)
iti_region = iti_ms(region);

med_iti = median(iti_region);
mad_iti = median(abs(iti_region - med_iti));

fprintf('Median ITI in trials 1-498: %.2f ms (MAD %.2f ms)\n', med_iti, mad_iti);

robust_z = (iti_region - med_iti) / (1.4826 * max(mad_iti, eps));

% Show the TOP 10 largest gaps in this region regardless of a hard
% threshold, so we don't miss a borderline one.
[sorted_z, order] = sort(robust_z, 'descend');
fprintf('\nTop 10 largest inter-trigger intervals in trials 1-498 (by robust z):\n');
for k = 1:10
    i = region(order(k));
    fprintf('  Between pulse(%d) and pulse(%d): %.1f ms (robust z = %.2f) [samples %d -> %d]\n', ...
        i, i+1, iti_ms(i), sorted_z(k), stimDig(i), stimDig(i+1));
end

fprintf(['\nA genuinely doubled gap (one missing pulse) should print here at roughly\n' ...
    '2x the median (~%.0f ms), even if it is not dramatically larger than everything\n' ...
    'else. Compare the top few candidates above to that -- whichever one is closest\n' ...
    'to double the median, and isolated (not near the pause region), is the most\n' ...
    'likely spot for the second missing trigger.\n'], 2*med_iti);