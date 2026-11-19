function [t_ms,block] = Artifact_check(trialIdx, win_ms, filepath)
% Artifact_check  Plot raw artifact for ONE trial from amplifier.dat
%
%   Artifact_check(trialIdx)
%   Artifact_check(trialIdx, win_ms)
%   Artifact_check(trialIdx, win_ms, folder)
%
%   Inputs
%   trialIdx : which trigger index to plot (1-based)
%   win_ms   : [t0 t1] window in ms around trigger (default [-2 10])
%   filepath   : path containing info.rhs and amplifier.dat (default pwd)
%
%   Outputs
%   t_ms     : time vector for the window (ms)
%   block    : [nChn × time] raw data (µV) for that trial
%
%   Notes
%   - Uses ONLY amplifier.dat
%   - Row-major plotting (4×8 grid: Ch1..Ch8 left→right, next row Ch9..Ch16, etc.)
%   - Returns the snippet so you can reuse it if needed


if nargin < 2 || isempty(win_ms), win_ms = [-2, 8]; end
if nargin < 3 || isempty(filepath), filepath = pwd; end

% ---- Intan Files ----
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

% ---- Data File ----
dName='amplifier';
datFile = fullfile(filepath, [dName '.dat']);
vFID = fopen([filepath filesep dName '.dat'],'r'); % read data file
assert(vFID~=-1, 'Data file not found: %s', datFile);

% For bounds checking
info = dir(datFile);
bytes_per_sample = 2;            % int16
total_samples_per_ch = info.bytes / (nChn * bytes_per_sample);

% Triggers (sample indices)
trig = loadTrig(0); 

%% ---- Artifact extract for ONE specific trial ----
Trial_to_extract = trialIdx;                  % trial index
assert(Trial_to_extract >= 1 && Trial_to_extract <= numel(trig), 'Trial index out of range.');

% Window around trigger
win_samp = round(win_ms/1000 * FS);      % samples relative to trigger
win_len  = diff(win_samp) + 1;
t_ms     = (win_samp(1):win_samp(2)) / FS * 1000;

% byte offset and read block
t0   = trig(Trial_to_extract);           % trigger (in samples)
samp_start = t0 + win_samp(1);           % absolute start sample (per channel)
samp_end   = t0 + win_samp(2);

% Bounds check (skip/clip if near file edges)
if samp_start < 0 || samp_end >= total_samples_per_ch
    error('Selected window [%d,%d] samples is out of file bounds for trial %d.', samp_start, samp_end, Trial_to_extract);
end

start_byte = int64(samp_start) * nChn * bytes_per_sample;
fseek(vFID, start_byte, 'bof');
block = fread(vFID, [nChn, win_len], 'int16') * 0.195;   % µV

%% ---- Plot: all channels, one trial ----
nRows = 8; nCols = 4;
figure('Name', sprintf('Artifact | Trial %d', Trial_to_extract), 'Color', 'w');
tl = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% Compute a symmetric y-limit per figure
absmax = max(abs(block), [], 'all');
y_lim  = ceil(absmax/10)*10;   % round to nearest 10 µV

for row = 1:nRows
    for col = 1:nCols
        ch = (col-1)*nRows + row; 
        if ch > nChn, break; end
        ax = nexttile((row-1)*nCols + col);
        plot(t_ms, block(ch, :), 'k'); hold on;
        xline(0, '--r', '0 ms', 'LabelVerticalAlignment','bottom', 'LabelOrientation','horizontal');
        yline(0, ':');

        title(sprintf('Ch %d', ch), 'FontSize', 10);
        xlabel('Time (ms)', 'FontSize', 7);
        ylabel('µV', 'FontSize', 7);
        xlim([t_ms(1) t_ms(end)]);
        if y_lim > 0, ylim([-y_lim y_lim]); end
        box off;
    end
end



