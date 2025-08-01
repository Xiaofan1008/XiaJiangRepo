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
nSamps = round(0.005 * FS); % 5 ms = 150 samples

% === Open amplifier file ===
filepath = pwd; % Or set manually
amplifier_file = fullfile(filepath, 'amplifier.dat');
fid = fopen(amplifier_file, 'r');

% === Preallocate array ===
all_uVdata = zeros(nChn, nSamps, nTrials); % Channels x Time x Trials

% === Loop to extract artifact traces ===
for i = 1:nTrials
    offset = int64(nChn * 2 * trig(i)); % 2 bytes per int16
    fseek(fid, offset, 'bof');
    data = fread(fid, [nChn, nSamps], 'int16');
    if size(data,2) < nSamps
        warning("Trial %d too short. Skipping.", i);
        continue;
    end
    uVdata = data .* 0.195; % Convert to µV
    all_uVdata(:,:,i) = uVdata; % Store entire uV trace
end
fclose(fid);

% === Get all unique stimulation IDs ===
uniqueIDs = unique(trialIDs);
shankMap = [24; 8; 7; 23; 26; 10; 5; 21; 0; 16; 1; 11; 31; 15; 25; 20; 30; 17; 6; 12; 2; 14; 3; 19; 29; 9; 4; 13; 28; 18; 27; 22];  % 32x1
ylims = [-6000, 6000];  % µV fixed scale
nRows = 32; nCols = 1;

for id = uniqueIDs(:)'  % loop through each unique trial ID
    % === Find trials matching this ID ===
    idx = trialIDs == id;

    % === Average uVdata across trials for this ID ===
    mean_uVdata = mean(all_uVdata(:,:,idx), 3, 'omitnan'); % 32 x nSamps

    % === Plot the averaged LFP traces ===
    figure('Name', sprintf('ID %d', id), 'Color', 'w', 'Position', [100, 100, 400, 1600]);
    sgtitle(sprintf('Stimulation at ID %d electrode site', id), 'FontWeight', 'bold');

    for row = 1:nRows
        ch = shankMap(row) + 1;  % convert to 1-based indexing
        subplot(nRows, nCols, row);
        plot(mean_uVdata(ch,:), 'k', 'LineWidth', 1);
        text(1, 0, sprintf('Ch %d', ch-1), 'HorizontalAlignment', 'right', 'FontSize', 6);
        ylim(ylims);
        axis off;
    end
end