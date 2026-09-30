%% DIAGNOSE_LARGE_TRIGGER_GAP
%
% This dataset expects 270 triggers (n_Trials=270, simultaneous_stim=1)
% but cleanTrig_sabquick only found 201 -- a shortfall of 69, much larger
% than the single-missing-pulse case fixed earlier. This script looks at
% the raw digitalin.dat pulses directly (independent of trig.dat/
% cleanTrig_sabquick's own corrections) to see:
%   1) how many raw pulses are actually there before any correction,
%   2) whether the missing ones show up as real gaps (likely, given how
%      many are missing) and where those gaps are concentrated,
%   3) what the overall digitalin.dat trace looks like, to spot a dead/
%      disconnected stretch.
%
% Diagnostic only -- does not modify any files.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_100_150um_Single1_260929_144315';
cd(data_folder);

[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;

fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
fprintf('Total recording length: %.1f s\n', nSam/FS);

digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);

% Same exact-match detection cleanTrig_sabquick itself uses (bit 0 == 1
% exactly, i.e. digital_in == 1)
stimDig = flip(find(digital_in == 1));
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
stimDig = flip(stimDig);   % chronological

fprintf('Raw pulses detected (exact-match method): %d  (expected 270)\n', numel(stimDig));

%% ====================== GAP ANALYSIS OVER THE WHOLE RECORDING ======================
iti_samples = diff(stimDig);
iti_ms = iti_samples / FS * 1000;

med_iti = median(iti_ms);
mad_iti = median(abs(iti_ms - med_iti));

fprintf('Median ITI: %.2f ms (MAD %.2f ms)\n', med_iti, mad_iti);

robust_z = (iti_ms - med_iti) / (1.4826 * max(mad_iti, eps));

% Print every gap that looks like it could hide one or more missing
% pulses (>= ~1.5x median), not just the top 10 -- there should be many,
% given 69 pulses are missing.
suspect = find(iti_ms > 1.5 * med_iti);
[~, order] = sort(iti_ms(suspect), 'descend');
suspect = suspect(order);

fprintf('\n%d intervals exceed 1.5x the median ITI (sorted, largest first):\n', numel(suspect));
total_implied_missing = 0;
for k = 1:numel(suspect)
    i = suspect(k);
    n_missing_here = round(iti_ms(i)/med_iti) - 1;
    total_implied_missing = total_implied_missing + max(n_missing_here,0);
    fprintf('  pulse(%d)->pulse(%d): %.1f ms (%.2fx median, z=%.2f) [samples %d -> %d] -- implies ~%d missing here\n', ...
        i, i+1, iti_ms(i), iti_ms(i)/med_iti, robust_z(i), stimDig(i), stimDig(i+1), n_missing_here);
end
fprintf('\nSum of implied-missing across all flagged gaps: %d  (need 69 to fully explain the shortfall)\n', total_implied_missing);

%% ====================== FULL-RECORDING RAW TRACE (DOWNSAMPLED) ======================
% Downsample for plotting so this doesn't choke on a long recording
step = max(1, floor(nSam / 2e6));
figure('Color','w','Position',[100 100 1400 400]);
plot((0:step:nSam-1)/FS, digital_in(1:step:end));
xlabel('Time (s)'); ylabel('digital\_in raw value');
title('Full digitalin.dat trace (downsampled)');