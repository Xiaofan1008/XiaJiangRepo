%% PRE-DEFINED VARIABLES
DATALIMIT = 60*(1e3); % How much data to load at maximum, in seconds
if ~exist('filepath','var') && ~exist('file','var')
    [file, path, ~] = ...
        uigetfile('*.rhs', 'Select an RHS2000 Data File', 'MultiSelect', 'off');
else    
    if ~exist('filepath','var')
        filepath = file;
    end
    file = 'info.rhs';
    path = [filepath '\'];    
end
%% MAIN FUNCTION
info = read_Intan_RHS2000_file(path,file);
nChn = length(info.amplifier_channels);
FS = info.frequency_parameters.amplifier_sample_rate;
fileinfo = dir([path 'stim.dat']);
if ~isempty(fileinfo)
    num_samples = fileinfo.bytes/(nChn * 2); % uint16 = 2 bytes    
    if num_samples > ((FS/1e3)*DATALIMIT)
        num_samples = ((FS/1e3)*DATALIMIT);
    end
    fid = fopen([path 'stim.dat'], 'r');
    data = fread(fid, [nChn, num_samples], 'uint16');
    fclose(fid);
    i = bitand(data, 255) * info.stim_parameters.stim_step_size; % current magnitude
    sign = (128 - bitand(data, 256))/128; % convert sign bit to 1 or -1
    i = i .* sign; % signed current in Amps
end
fileinfo = dir([path 'amplifier.dat']);
num_samples = fileinfo.bytes/(nChn * 2); % int16 = 2 bytes
if num_samples > ((FS/1e3)*DATALIMIT)
    num_samples = ((FS/1e3)*DATALIMIT);
end
fid = fopen([path 'amplifier.dat'], 'r');
v = fread(fid, [nChn, num_samples], 'int16');
fclose(fid);
v = v * 0.195; % convert to microvolts
fileinfo = dir([path 'digitalin.dat']);
if ~isempty(fileinfo)
    num_samples = fileinfo.bytes/2;
    if num_samples > ((FS/1e3)*DATALIMIT)
        num_samples = ((FS/1e3)*DATALIMIT);
    end
    fid = fopen([path 'digitalin.dat'],'r');
    digital_in = fread(fid, num_samples, 'uint16');
    fclose(fid);
    stimDig = flip(find(digital_in == 1)); % Fix for finding the falling line instead of the rising line
    visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
    dt = diff(stimDig);
    kill = dt == -1;
    stimDig(kill) = [];
    dt = diff(visDig);
    kill = dt == -1;
    visDig(kill) = [];
    nStamps = max([length(stimDig); length(visDig)]);
    time_stamps = nan(2,nStamps);
    time_stamps(1,1:length(stimDig)) = flip(stimDig);
    time_stamps(2,1:length(visDig)) = flip(visDig);
    time_stamps = time_stamps ./ (FS/1e3); % Convert to milliseconds
    if ~isempty(time_stamps)
        if isnan(time_stamps(2,2))
            time_stamps = time_stamps(1,:);
        else
            time_stamps = time_stamps(2,:);
        end
    end
end
%% Post-processing
[Mufilt,Lfpfilt] = generate_Filters;
MuNf = length(Mufilt);
LfpNf = length(Lfpfilt);
% Threshold
v_info = dir([path, 'amplifier.dat']);
threshBYTES = min([(FS/1e3)*(100*1000),v_info.bytes/64]);
fid = fopen([path, 'amplifier.dat'], 'r');
v_thresh = fread(fid,[nChn, threshBYTES],'int16') * 0.195;
thresh = zeros(2,nChn);
if ~exist('threshfac','var')
    threshfac = -3.5;
end
for c = 1:nChn
    tmp = conv(v_thresh(c,:),Mufilt);
    filt_v = tmp(1,MuNf/2:threshBYTES+MuNf/2-1);
    thresh(1,c) = threshfac*median(abs(filt_v))./0.6795;
    thresh(2,c) = -threshfac*median(abs(filt_v))./0.6795;
end
MU = zeros(nChn,length(v(1,:)));
X_HF = 1/30:1/30:length(MU(1,:))/30;
nData = length(X_HF);
X_LF = 1:1:length((LfpNf/2:FS/1e3:nData+LfpNf/2-1));
RAW = zeros(nChn,length(X_HF));
LFP = zeros(nChn,length(X_HF));
for c = 1:nChn    
    RAW(c,:) = v(c,:);
    tmp = conv(v(c,:),Mufilt);
    MU(c,:) = tmp(1,MuNf/2:nData+MuNf/2-1);
    tmp = conv(v(c,:),Lfpfilt);
    LFP(c,:) = tmp(1,LfpNf/2:nData+LfpNf/2-1);
end