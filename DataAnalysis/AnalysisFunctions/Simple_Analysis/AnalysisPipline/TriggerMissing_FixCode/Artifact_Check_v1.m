%% CHECK_ARTIFACT_IN_QUIET_WINDOWS
%
% The digital line (digitalin.dat) is missing 69 of the expected 270
% pulses, in three blocks: ~24 trials' worth of quiet time BEFORE the
% first detected pulse, ~24 AFTER the last detected pulse, and ~21 in one
% large gap mid-recording (between pulse 62 and 63). Gap analysis and the
% bit-collision check both ruled out a detection-method bug -- so the
% question is whether stimulation was actually DELIVERED in these three
% windows (just not logged on the digital line), by looking for real,
% periodic stim artifacts directly in the raw amplifier.dat trace, the
% same way recover_trial1_artifact.m checked the single missing pulse in
% the earlier dataset.
%
% Diagnostic only -- does not modify any files.

clear; clc;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_100_150um_Single1_260929_144315';
cd(data_folder);

[amp_channels, freq_params] = read_Intan_RHS2000_file;
FS = freq_params.amplifier_sample_rate;
nChn = numel(amp_channels);

% Same reference channels cleanTrig_sabquick itself uses for artifact
% checkpointing
ref_channels = [48, 52];
thresh_uV = 2000;   % start with the same threshold used before; loosen if nothing shows up

windows = struct( ...
    'label',  {'BEFORE first pulse (~24 trials implied)', 'AFTER last pulse (~24 trials implied)', 'MID-RECORDING gap (~21 trials implied)'}, ...
    'sample_start', {1,        2912182,  1004109}, ...
    'sample_end',   {290613,   3194880,  1270983} );

fid = fopen('amplifier.dat', 'r');

for w = 1:numel(windows)
    s0 = windows(w).sample_start;
    s1 = windows(w).sample_end;
    nSampWin = s1 - s0 + 1;

    fseek(fid, (s0-1) * nChn * 2, 'bof');
    data_block = fread(fid, [nChn, nSampWin], 'int16') * 0.195;

    figure('Color','w','Position',[100 100 1300 500],'Name',windows(w).label);
    for i = 1:numel(ref_channels)
        ch = ref_channels(i);
        subplot(numel(ref_channels),1,i);
        t_axis = (s0:s1)/FS;
        plot(t_axis, data_block(ch,:));
        hold on; yline(thresh_uV,'r--'); yline(-thresh_uV,'r--');
        title(sprintf('%s -- Channel %d', windows(w).label, ch), 'Interpreter','none');
        xlabel('Time (s)'); ylabel('\muV');

        exceeds = find(abs(data_block(ch,:)) > thresh_uV);
        if ~isempty(exceeds)
            % Count distinct events (require > 100 samples/~3ms apart to
            % avoid counting one artifact's ringing as several)
            gaps = diff(exceeds);
            n_events = 1 + sum(gaps > 100);
            fprintf('%s | Ch %d: %d threshold-crossing samples, ~%d distinct artifact events\n', ...
                windows(w).label, ch, numel(exceeds), n_events);
        else
            fprintf('%s | Ch %d: no samples exceed +/-%d uV\n', windows(w).label, ch, thresh_uV);
        end
    end
end

fclose(fid);

fprintf(['\nIf periodic artifacts show up at roughly 397ms spacing in these windows,\n' ...
    'stimulation was genuinely delivered there -- we can recover real trigger times\n' ...
    'directly from these artifacts instead of relying on digitalin.dat for these blocks.\n' ...
    'If nothing shows up, these trials likely were never delivered/recorded at all.\n']);