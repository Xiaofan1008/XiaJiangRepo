%% RECOVER_TRIAL1_ARTIFACT
%
% Trial 1's TTL pulse never registered on digitalin.dat (trig(1) is
% actually trial 2's pulse -- confirmed via trigOffset=-1 working
% uniformly across the whole file, with no doubled-gap anywhere in
% trig(1:20) because a missing FIRST pulse leaves no gap to detect).
%
% This checks whether trial 1's STIMULATION itself actually happened
% (just wasn't logged digitally) by looking for a real artifact on the
% blanking-reference channels (5, 15 -- same channels cleanTrig_sabquick
% uses for its own catch-trial checkpointing) in the window where
% trial 1's pulse should sit if delivery was on-schedule: roughly one
% median-ITI before trig(1).
%
% Diagnostic only -- does not modify any files.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_300_500um_SimSeq1';
cd(data_folder);

[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;
nChn = numel(amp_channels);

trig = loadTrig(0);

% Use the same "trials 1-498" ITI stats you already computed
med_iti_ms = 418.53;              % from find_earlier_missing_trigger.m
med_iti_samp = round(med_iti_ms/1000 * FS);

expected_center = trig(1) - med_iti_samp;   % where trial 1's pulse "should" be
win_samp = round(0.15 * FS);                % +/- 150 ms search window around it

search_start = max(1, expected_center - win_samp);
search_end   = trig(1) - 1;                 % up to (not including) trig(1) itself

fprintf('Searching for a trial-1 artifact between sample %d and %d (t = %.2f - %.2f s)\n', ...
    search_start, search_end, search_start/FS, search_end/FS);
fprintf('(trig(1) = %d, i.e. t = %.2f s; expected trial-1 center = %d, i.e. t = %.2f s)\n', ...
    trig(1), trig(1)/FS, expected_center, expected_center/FS);

%% ====================== READ RAW AMPLIFIER DATA IN THIS WINDOW ======================
fid = fopen('amplifier.dat','r');
byte_pos = search_start * nChn * 2;
fseek(fid, byte_pos, 'bof');
nSampWin = search_end - search_start + 1;
data_block = fread(fid, [nChn, nSampWin], 'int16') * 0.195;
fclose(fid);

%% ====================== CHECK BLANKING-REFERENCE CHANNELS (5, 15) ======================
% Same threshold cleanTrig_sabquick uses for catch-trial artifact detection
ref_channels = [5, 15];
thresh_uV = 2000;

figure('Color','w','Position',[150 150 1200 600]);
for i = 1:numel(ref_channels)
    ch = ref_channels(i);
    subplot(2,1,i);
    t_axis = (search_start:search_end)/FS;
    plot(t_axis, data_block(ch,:));
    hold on;
    yline(thresh_uV,'r--'); yline(-thresh_uV,'r--');
    xline(trig(1)/FS,'g--','trig(1)');
    xline(expected_center/FS,'m--','expected trial-1 center');
    title(sprintf('Channel %d (raw amplifier.dat), search window', ch));
    xlabel('Time (s)'); ylabel('\muV');

    exceeds = find(abs(data_block(ch,:)) > thresh_uV);
    if ~isempty(exceeds)
        fprintf('Channel %d: found %d samples exceeding +/-%d uV, first at sample %d (t=%.4f s)\n', ...
            ch, numel(exceeds), thresh_uV, search_start+exceeds(1)-1, (search_start+exceeds(1)-1)/FS);
    else
        fprintf('Channel %d: NO samples exceed +/-%d uV in this window.\n', ch, thresh_uV);
    end
end

fprintf(['\nIf either channel shows a clear artifact blip inside this window, trial 1''s\n' ...
    'stimulation was actually delivered -- we can recover its true sample index and give\n' ...
    'trial 1 real data instead of discarding it. If nothing shows up at all, the pulse\n' ...
    'was likely never delivered (not just unlogged), and trial 1 should stay excluded.\n']);