%% High-Frequency analysis function
function HF_D(filepath)
%% Deal with input error
if ~(nargin)
    disp('I need a filepath!');
    return;
elseif nargin < 2
    datafile = dir([filepath, '\*_exp_datafile*.mat']);
    if isempty(datafile)
        disp('I need a datafile!');
        return;
    end
    datafile = ['\', datafile.name];
end
if nargin < 3
    threshfac = -3.5;
end
%% Load in Intan data
info = read_Intan_RHS2000_file(filepath,'\info.rhs',0);
clearvars -except info filepath datafile threshfac
FS = info.frequency_parameters.amplifier_sample_rate;       % Sampling frequency
nChn = length(info.amplifier_channels);
%% Load in the digital lines
fileinfo = dir([filepath, '\digitalin.dat']);
if ~isempty(fileinfo)
    num_samples = fileinfo.bytes/2;
    digin_fid = fopen([filepath, '\digitalin.dat'],'r');
    % I should note that the digitalin.dat file in this case has two rows,
    % one for electrical stimulation, and one for visual stimulation
    digital_in = fread(digin_fid, num_samples, 'uint16');
    fclose(digin_fid);
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
    if isnan(time_stamps(2,2))
        time_stamps = time_stamps(1,:);
    else
        time_stamps = time_stamps(2,:);
    end
end
%% Load in the datafile and select the Trial ID for analysis
AMP = [];DUR = [];CHN = [];
load([filepath, datafile],'AMP','DUR','n_REP','CHN','StimParams','n_Trials','TrialParams');
if (n_Trials ~= nStamps)
    disp('Warning: Digital Lines do not match expected number of trials');
    return;
end
uniqueTrials = length(AMP)*length(DUR)*length(CHN);
%% Build out the filters
% Generate a MUA filter
[Mufilt,~] = generate_Filters;
MuNf = length(Mufilt);
LOOP = 1;
v_fid = fopen([filepath, '\amplifier.dat'],'r');
BIN = [-1000, 2000];    % Extraction window for reading in data
nData = diff(BIN)*30+1;  % Number of datapoints
% Load in a little data for thresholding
v_info = dir([filepath, '\amplifier.dat']);
threshBYTES = min([(FS/1e3)*(100*1000),v_info.bytes/2]);
v_thresh = fread(v_fid,[nChn, threshBYTES],'int16') * 0.195;
% Generate a threshold
thresh = zeros(1,nChn);
for c = 1:nChn
    thresh(c) = threshfac*median(abs(v_thresh(c,:)))./0.6745;
end
depth = Depth+1;
while (LOOP ~= 0)
    % Administration
    disp('Event loop. Set LOOP = 0 to exit. Set LOOP = n to target that trial ID');
    disp(['AMP: ' num2str(AMP(1)) ':' num2str(AMP(end))]);
    disp(['DUR: ' num2str(DUR(1)) ':' num2str(DUR(end))]);
    disp(['CHN: ' num2str(CHN(1)) '|' num2str(CHN(end))]);
    LOOP = input(['Please select a trial identifier between 1 and ' num2str(uniqueTrials) '.\n'],'s');
    tic;
    LOOP = str2double(LOOP);
    if (LOOP == 0)
        break;
    elseif (LOOP > n_REP)
        continue;
    end
    allTrials = cell2mat(TrialParams(2:(uniqueTrials*n_REP)+1,2));
    selectedTrials = zeros(2,n_REP);
    selectedTrials(1,:) = find(allTrials == LOOP);
    selectedTrials(2,:) = time_stamps(1,selectedTrials(1,:));
    disp(['Display PSTHs for Trial ID n = ' num2str(LOOP)]);
    % Loading, filtering, and blanking.
    X = 1/30:1/30:nData/30;
    for r = 1:n_REP        
        OFFSET = cast(nChn*2*(FS/1e3)*(selectedTrials(2,r)+BIN(1)),'int64');
        fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
        v = fread(v_fid,[nChn,(FS/1e3)*diff(BIN)],'int16') * 0.195;
        %% Extract spikes from v. OUTPUT will contain n_REP cell arrays
        % Blank the artefact.
        v(:,(abs(BIN(1))-1)*30:(abs(BIN(1))+3)*30) = 0;
        filt_v = zeros(nChn,nData);
        for c = 1:nChn
            chan = depth(c);
            % Filter data
            tmp = conv(v(chan,:),Mufilt);
            filt_v(c,:) = tmp(1,MuNf/2:nData+MuNf/2-1);            
        end
    end
end