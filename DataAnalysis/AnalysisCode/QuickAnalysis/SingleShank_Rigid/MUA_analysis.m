%cleanTrig_sabquick;

% === Load trigger timestamps and trial info ===
trig = loadTrig(0);
TrialParams = loadTrialParams;
trialIDs = cell2mat(TrialParams(:,2));

% === Restrict to first N trials ===
nKeep = 160;
trig = trig(1:nKeep);
trialIDs = trialIDs(1:nKeep);

% === Artifact parameters ===
zero_start_ms = -10;
zero_end_ms = 10;
zero_idx = (time_axis >= zero_start_ms) & (time_axis < zero_end_ms);

% === Settings ===
nTrials = length(trig);
nChn = 32;
FS = 30000;
window_ms = [-100, 100];  % in milliseconds
nSamps = round(diff(window_ms)/1000 * FS);  % 100 ms window

% === Open amplifier file ===
filepath = pwd; % Or set manually
amplifier_file = fullfile(filepath, 'amplifier.dat');
fid = fopen(amplifier_file, 'r');

% === Bandpass filter design (for MUA) ===
d = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', 300, ...
    'HalfPowerFrequency2', 3000, ...
    'SampleRate', FS);

% === Preallocate array ===
all_mua = zeros(nChn, nSamps, nTrials);

for i = 1:nTrials
    t0 = trig(i) + round(window_ms(1)/1000 * FS);
    offset = int64(nChn * 2 * t0); % byte offset
    fseek(fid, offset, 'bof');
    data = fread(fid, [nChn, nSamps], 'int16');
    if size(data,2) < nSamps
        warning("Trial %d too short. Skipping.", i);
        continue;
    end
    uVdata = double(data) * 0.195;

    % Extract MUA: bandpass → rectify → smooth
    mua_data = zeros(size(uVdata));
    for ch = 1:nChn
        filtered = filtfilt(d, uVdata(ch, :));
        rectified = abs(filtered);
        smoothed = movmean(rectified, 15);  % ~0.5 ms smoothing
        mua_data(ch, :) = smoothed;
    end

    all_mua(:,:,i) = mua_data;
end
fclose(fid);

% === Time axis in ms ===
time_axis = linspace(window_ms(1), window_ms(2), nSamps);

% === Get all unique stimulation IDs ===
uniqueIDs = unique(trialIDs);
shankMap = [24; 8; 7; 23; 26; 10; 5; 21; 0; 16; 1; 11; 31; 15; 25; 20; 30; 17; 6; 12; 2; 14; 3; 19; 29; 9; 4; 13; 28; 18; 27; 22];  % 32x1
nRows = 32; nCols = 1;

% === Loop through each stimulation site ID ===
for id = uniqueIDs(:)'  % force row vector
    % === Find matching trials ===
    idx = trialIDs == id;

    % === Average across trials ===
    mean_mua = mean(all_mua(:,:,idx), 3, 'omitnan');
    mean_mua(:, zero_idx) = 0;  % zero artifact period if needed

    % === Set y-axis limits ===
    ylims = [0, max(mean_mua(:)) * 1.1];

    % === Create figure for this ID ===
    figure('Name', sprintf('MUA - ID %d', id), ...
           'Color', 'w', 'Position', [100, 100, 400, 1600]);
    sgtitle(sprintf('Stimulation at ID %d electrode site', id), 'FontWeight', 'bold');
    
    for row = 1:nRows
        ch = shankMap(row) + 1;
        subplot(nRows, nCols, row);
        plot(time_axis, mean_mua(ch,:), 'k', 'LineWidth', 1);
        text(time_axis(1)-5, ylims(2)*0.5, sprintf('Ch %d', ch-1), ...
             'HorizontalAlignment', 'right', 'FontSize', 6);
        ylim(ylims);
        xlim([window_ms(1), window_ms(2)]);
        % axis off;
    end
end