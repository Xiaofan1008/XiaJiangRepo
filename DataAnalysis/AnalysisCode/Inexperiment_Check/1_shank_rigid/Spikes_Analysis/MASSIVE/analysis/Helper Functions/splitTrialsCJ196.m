function splitTrialsCJ196(WHICH)
%% Fix function for splitting trials around off-channel stimulation
filepath = '\\ad.monash.edu\home\User030\tall0003\Desktop\CJ196\CJ196_Intan_datafile001_180831_095228';
datafile = dir([filepath, '\*_exp_datafile*.mat']);
if isempty(datafile)
    disp('I need a datafile!');
    return;
end
datafile = ['\', datafile.name];
threshfac = -3.5;
info = read_Intan_RHS2000_file(filepath,'\info.rhs',0);
clearvars -except info filepath datafile threshfac WHICH
FS = info.frequency_parameters.amplifier_sample_rate;       % Sampling frequency
nChn = length(info.amplifier_channels);
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
StimParams = [];AMP = [];DUR = [];CHN = [];
load([filepath, datafile],'AMP','DUR','n_REP','CHN','StimParams','n_Trials','TrialParams');
if (n_Trials ~= nStamps)
    disp('Warning: Digital Lines do not match expected number of trials');
    return;
end
uniqueTrials = length(AMP)*length(DUR)*length(CHN);

BIN = [-100, 200];
OFFSET = cast(nChn*2*(FS/1e3)*(time_stamps(2)+BIN(1)),'int64');


% KEY
STIM_COND = zeros(3,n_Trials);
for i = 1:500
    if strcmp(StimParams{i+1,1},'A-001')
        STIM_COND(1,i) = StimParams{i+1,16};
        STIM_COND(2,i) = STIM_COND(2,i-1);
    else
        STIM_COND(2,i) = StimParams{i+1,16};
        if i == 1
            STIM_COND(1,i) = 0;
        else
            STIM_COND(1,i) = STIM_COND(1,i-1);
        end
    end
end
for i = 1:n_Trials
    A = find(AMP == STIM_COND(1,i),1);
    B = find(AMP == STIM_COND(2,i),1);
    STIM_COND(3,i) = A + (B-1)*length(AMP);
end

selectedTrials(1,:) = find(STIM_COND(3,:) == WHICH);
selectedTrials(2,:) = time_stamps(1,selectedTrials(1,:));

%% Build out the filters
% Generate a MUA filter
[Mufilt,Lfpfilt] = generate_Filters;
MuNf = length(Mufilt);
LfpNf = length(Lfpfilt);
v_fid = fopen([filepath, '\amplifier.dat'],'r');
BIN = [-100, 200];    % Extraction window for reading in data
X = BIN(1):BIN(2);
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
if (nChn == 32)
    ROW = 6;
    COL = 6;
    A = zeros(6,6);
    for i = 1:32
        A(i+3) = i;
    end
    for i = 1:28
        A(i+2) = i;
    end
    A(1,6) = 0;
    for i = 1:4
        A(i+1) = i;
    end
    A(6,1) = 0;
end
depth = Depth+1;
XLABEL = [25,30,32,33,34,35];
YLABEL = [2,7,13,19,25,32];
ARTIFACT_CUTOFF = 100;
PSTH = figure('units','normalized','outerposition',[0 0 1 1]);
if (STIM_COND(1,selectedTrials(1,1)) ~= 0)
    STIMCHN = 'A-001';
elseif (STIM_COND(2,selectedTrials(1,1)) ~= 0)
    STIMCHN = 'A-017';
else
    STIMCHN = 'A-094';
end
STIMCHN = STIMCHN(4:5);
lfp = cell(1,nChn);
MAX = 0;
for c = 1:nChn
    lfp{c} = zeros(1,301);
end
LFP = figure('units','normalized','outerposition',[0 0 1 1]);
hold on
for r = 1:length(selectedTrials(1,:))
    OFFSET = cast(nChn*2*(FS/1e3)*(selectedTrials(2,r)+BIN(1)),'int64');
    fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
    v = fread(v_fid,[nChn,(FS/1e3)*diff(BIN)],'int16') * 0.195;
    %% Extract spikes from v. OUTPUT will contain n_REP cell arrays
    % Blank the artefact.
    v(:,(abs(BIN(1))-1)*30:(abs(BIN(1))+3)*30) = 0;
    for c = 1:nChn
        chan = depth(c);
        % Filter data
        tmp = conv(v(chan,:),Mufilt);
        filt_v = tmp(1,MuNf/2:nData+MuNf/2-1);
        % Extract spikes
        [WAVES,Sp_tmp] = spikeextract(filt_v,thresh(c),FS);
        % Remove artifact outliers
        if sum(Sp_tmp) ~= 0
            for l = 1:length(Sp_tmp)
                if min(WAVES(l,:)) < -ARTIFACT_CUTOFF || max(WAVES(l,:)) > ARTIFACT_CUTOFF
                    Sp_tmp(l) = NaN;
                end
            end
            Sp_tmp(isnan(Sp_tmp)) = [];
        else
            Sp_tmp(Sp_tmp == 0) = [];
        end
        OUTPUT{c,r} = Sp_tmp;        
    end
    for c = 1:nChn
        chan = depth(c);
        % Filter data
        tmp = conv(v(chan,:),Lfpfilt);
        lfp{c} = lfp{c} + tmp(1,LfpNf/2:FS/1e3:nData+LfpNf/2-1);
        if max(lfp{c}) > MAX
            MAX = max(lfp{c});
        end
    end
end
SMOOTHING = 5;
H = figure;
MAXRATE = 0;
for c = 1:nChn
    INDEX = find(A' == c);
    subplot(ROW,COL,INDEX)
    [rate,~] = psth(OUTPUT(c,:),[-100,200]);
    if max(rate) > MAXRATE
        MAXRATE = max(rate);
    end
end
close(H);
for c = 1:nChn
    % Work out the subplot position
    chan = depth(c);
    INDEX = find(A' == c);
    figure(PSTH);
    subplot(ROW,COL,INDEX);
    psth(OUTPUT(c,:),[-100,200],SMOOTHING,MAXRATE);
    if (chan-1) < 10
        if ~strcmp(STIMCHN(2),num2str(chan-1))
            box off
        end
    elseif ~strcmp(STIMCHN,num2str(chan-1))
        box off
    elseif strcmp(STIMCHN,num2str(chan-1))
        box on
    end
    if find(XLABEL == INDEX,1)
        xlabel('Time (msec)');
    end
    if find(YLABEL == INDEX,1)
        ylabel('Firing rate (Spikes/sec)');
    end
    axis([-100 200 0 MAXRATE])
    title(['CHN ' num2str(chan-1) ' | Depth ' num2str(c)]);
end
figure(LFP);
hold on
SEP = MAX ./ n_REP;
if str2double(STIMCHN) < 10
    STIMCHN = str2double(STIMCHN(2));
else
    STIMCHN = str2double(STIMCHN);
end
STIMCHN = find(depth-1 == STIMCHN,1);
if isempty(STIMCHN)
    STIMCHN = 50;
end
for ch = 1:nChn
    % VERY IMPORTANT. FILTERING JUST ADDS EACH TRIAL TOGETHER
    lfp{ch} = lfp{ch} ./ n_REP;
    if ch == STIMCHN
        plot(X,lfp{ch}+(SEP*(STIMCHN-ch)),'Color','r');
    else
        plot(X,lfp{ch}+(SEP*(STIMCHN-ch)),'Color','k');
    end
end
xlabel('Time (msec)');
ylabel('Voltage (uV) - with Offset');
yl = ylim;
line([0 0],[yl(1) yl(2)],'Color','blue');
end