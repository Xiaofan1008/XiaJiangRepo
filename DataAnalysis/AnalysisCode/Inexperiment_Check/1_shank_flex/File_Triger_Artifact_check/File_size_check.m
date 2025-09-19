% File size and duration check 

% ---- Intan Files ----
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

% ---- Data File ----
dName='amplifier';
datFile = fullfile(filepath, [dName '.dat']);
vFID = fopen([filepath filesep dName '.dat'],'r'); % read data file

if ~exist(datFile,'file')
    error('Data file not found: %s', datFile);
end

info = dir(datFile);
bytes = info.bytes;
bytes_per_sample = 2;  % int16
total_samples_per_channel = bytes / (nChn * bytes_per_sample);
duration_s = total_samples_per_channel / FS;
duration_ms = duration_s*1000;
duration_mins = duration_s/60;

fprintf('\n=== FILE / DURATION SUMMARY ===\n');
fprintf('File: %s\n', datFile);
fprintf('Size: %.2f MB (%.2f GB)\n', bytes/1e6, bytes/1e9);
fprintf('Channels: %d, FS: %.0f Hz, datatype: int16 (%d bytes)\n', nChn, FS, bytes_per_sample);
fprintf('Total samples per channel: %.0f\n', total_samples_per_channel);
fprintf('Recording duration: %.3f s (%.3f mins)\n', duration_s, duration_mins);

