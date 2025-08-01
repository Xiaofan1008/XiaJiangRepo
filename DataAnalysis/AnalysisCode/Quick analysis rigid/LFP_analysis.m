%cleanTrig_sabquick;

% === Load trigger timestamps and trial info ===
trig = loadTrig(0);
TrialParams = loadTrialParams;
trialIDs = cell2mat(TrialParams(:,2));

% === Restrict to first N trials ===
nKeep = 160;
trig = trig(1:nKeep);
trialIDs = trialIDs(1:nKeep);

% === Settings ===
nTrials = length(trig);
nChn = 32;
FS = 30000;
window_ms = [-100, 100];  % in milliseconds
nSamps = round(diff(window_ms)/1000 * FS);  % 100 ms window

% === Time axis ===
time_axis = linspace(window_ms(1), window_ms(2), nSamps);

% === Artifact zeroing window ===
zero_start_ms = 0;
zero_end_ms = 0;
zero_idx = (time_axis >= zero_start_ms) & (time_axis < zero_end_ms);

% === Open amplifier file ===
filepath = pwd; % Or set manually
amplifier_file = fullfile(filepath, 'amplifier.dat');
fid = fopen(amplifier_file, 'r');

% === Low-pass filter for LFP ===
lfpFilt = designfilt('lowpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency', 300, ...
    'SampleRate', FS);

% === Preallocate array ===
all_lfp = zeros(nChn, nSamps, nTrials);

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

    % Extract LFP: low-pass filter
    lfp_data = zeros(size(uVdata));
    for ch = 1:nChn
        lfp_data(ch,:) = filtfilt(lfpFilt, uVdata(ch,:));
    end

    all_lfp(:,:,i) = lfp_data;
end
fclose(fid);

% === Get all unique stimulation IDs ===
uniqueIDs = unique(trialIDs);
shankMap = [24; 8; 7; 23; 26; 10; 5; 21; 0; 16; 1; 11; 31; 15; 25; 20; 30; 17; 6; 12; 2; 14; 3; 19; 29; 9; 4; 13; 28; 18; 27; 22];  % 32x1
nRows = 32; nCols = 1;

% === Loop over each ID ===
for id = uniqueIDs(:)'  % Make sure it's a row vector
    % === Find trials with this ID ===
    idx = trialIDs == id;
    
    % === Average LFP across trials ===
    mean_lfp = mean(all_lfp(:,:,idx), 3, 'omitnan');
    mean_lfp(:, zero_idx) = 0;  % Artifact blanking
    
    % === Y-limits consistent for each figure ===
    ylims = [min(mean_lfp(:)) * 1.1, max(mean_lfp(:)) * 1.1];
    
    % === Plotting ===
    figure('Name', sprintf('Stimulation ID %d', id), ...
           'Color', 'w', 'Position', [100, 100, 400, 1600]);
    sgtitle(sprintf('Stimulation at ID %d electrode site', id), 'FontWeight', 'bold');
    
    for row = 1:nRows
        ch = shankMap(row) + 1;
        subplot(nRows, nCols, row);
        plot(time_axis, mean_lfp(ch,:), 'k', 'LineWidth', 1);
        text(time_axis(1)-5, ylims(2)*0.5, sprintf('Ch %d', ch-1), ...
             'HorizontalAlignment', 'right', 'FontSize', 6);
        ylim(ylims);
        xlim([window_ms(1), window_ms(2)]);
        axis off;
    end
end
