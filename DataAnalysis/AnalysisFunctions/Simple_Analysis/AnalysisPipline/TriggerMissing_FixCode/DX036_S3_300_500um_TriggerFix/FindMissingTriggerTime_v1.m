%% RECOVER_MISSING_TRIGGER_FROM_DIGITALIN
%
% cleanTrig_sabquick.m detects the stim trigger line via:
%     stimDig = find(digital_in == 1)
% which requires ALL other digital input bits to be exactly 0 at that
% sample. If any other digital line (e.g. the "vis" line, bit 1, or noise
% on an unused pin) was also high at the same sample as a real stim pulse,
% that pulse is invisible to this method -- it looks "missing" even though
% it really happened on time.
%
% This re-detects the stim line using the actual bit (bitget), which is
% robust to other bits being high at the same time, and compares it
% against the exact-match method cleanTrig_sabquick uses, to find the
% pulse(s) that were silently dropped.
%
% Diagnostic only -- prints/plots, does not touch your .trig.dat.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_300_500um_SimSeq1';
cd(data_folder);

[amplifier_channels, frequency_parameters] = read_Intan_RHS2000_file;
FS = frequency_parameters.amplifier_sample_rate;

%% ====================== LOAD RAW DIGITAL LINE ======================
fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);

%% ====================== OLD (FRAGILE) METHOD ======================
stimDig_old = flip(find(digital_in == 1));
dt = diff(stimDig_old);
stimDig_old(dt == -1) = [];
stimDig_old = flip(stimDig_old);   % back to chronological order

fprintf('Old method (digital_in == 1, exact match): %d pulses detected.\n', numel(stimDig_old));

%% ====================== NEW (BIT-ROBUST) METHOD ======================
stimBit = bitget(digital_in, 1);   % bit 0 = stim line, regardless of other bits
% Rising edges: 0 -> 1 transitions
rising = find(diff([0; stimBit]) == 1);
% Falling edges: 1 -> 0 transitions (to mirror the "falling edge" convention
% cleanTrig_sabquick actually uses)
falling = find(diff([stimBit; 0]) == -1) ;

fprintf('New method (bitget, robust to other bits): %d rising edges, %d falling edges.\n', ...
    numel(rising), numel(falling));

%% ====================== FIND WHAT'S NEW ======================
% Pulses present in the new (robust) detection but not near any pulse in
% the old detection -- these are the ones the old method silently missed.
tol_samples = round(1e-3 * FS);   % 1 ms tolerance for "same pulse"
isNew = true(size(falling));
for i = 1:numel(falling)
    if any(abs(stimDig_old - falling(i)) <= tol_samples)
        isNew(i) = false;
    end
end
recovered = falling(isNew);

fprintf('\n%d pulse(s) found by the robust method that the old method missed:\n', numel(recovered));
for i = 1:numel(recovered)
    samp = recovered(i);
    % where does it sit relative to its neighbours in the old detection?
    before = stimDig_old(find(stimDig_old < samp, 1, 'last'));
    after  = stimDig_old(find(stimDig_old > samp, 1, 'first'));
    fprintf('  Sample %d (t = %.3f s). Sits between old pulses at sample %d and %d ', ...
        samp, samp/FS, before, after);
    if ~isempty(before)
        fprintf('(%.1f ms after previous, ', (samp-before)/FS*1000);
    end
    if ~isempty(after)
        fprintf('%.1f ms before next)', (after-samp)/FS*1000);
    end
    fprintf('\n');
end

fprintf(['\nCompare the gap timings above to your recording''s usual inter-trial\n' ...
    'spacing (~420 ms, from the earlier diagnostic). If a recovered sample here\n' ...
    'splits the old ~2850 ms gap (between trig(498) and trig(499)) into two\n' ...
    'normal-looking ~420 ms gaps, that is very strong confirmation this is the\n' ...
    'true missing pulse, and its exact sample index is what should replace the\n' ...
    '-500 sentinel cleanTrig_sabquick inserted.\n']);