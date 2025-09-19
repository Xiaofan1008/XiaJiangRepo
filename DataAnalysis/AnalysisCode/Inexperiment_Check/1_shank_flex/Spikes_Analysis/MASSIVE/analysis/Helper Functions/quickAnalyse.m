%% A series of scripts designed to quickly generate example plots of response profiles from large Intan datafiles
function quickAnalyse(t,filepath,TRIAL)
% t determines the type of analysis produced.
% 1 for a PSTH based around a small sample of trials
% 2 for tuning curves built piecemeal from small samples of trials
% 3 for a few LFP traces
% 4 for a spectrogram of a sample set of trials
close all
%% Global variables
if (nargin < 1)
    fprintf('No TYPE of analysis provided.\n t = 1: RAW waveforms.\n t = 2: High-pass filtered response.\n t = 3: Low-pass filtered response.\n t = 4: PSTH responses.\n t = 5: Spectrogram responses.\n');
    return
end
if (nargin < 3)
    % No trial specified
    TRIAL = 20;
    if (nargin < 2)
        % No filepath provided
        filepath = uigetdir;
    end
end
filepath = [filepath '\'];
info = read_Intan_RHS2000_file(filepath,'info.rhs',0);
clearvars -except info t filepath TRIAL
FS = info.frequency_parameters.amplifier_sample_rate;       % Sampling frequency
nChn = length(info.amplifier_channels);
home = pwd;
cd(filepath) % n.b. This line is actually slow. If I really wanted to optimise, I'd remove this, and use the filepath to locate everything.
%% Clean and save the digital line
fileinfo = dir('digitalin.dat');
if ~isempty(fileinfo)
    num_samples = fileinfo.bytes/2;
    digin_fid = fopen('digitalin.dat','r');
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
else
    if isempty(dir('time_stamps.mat'))
        fileinfo = dir('stim.dat');
        num_samples = fileinfo.bytes/(nChn * 2); % uint16 = 2 bytes
        fid = fopen('stim.dat', 'r');
        data = fread(fid, [nChn, num_samples], 'uint16');
        fclose(fid);
        i = bitand(data, 255) * info.stim_parameters.stim_step_size; % current magnitude
        sign = (128 - bitand(data, 256))/128; % convert sign bit to 1 or -1
        i = i .* sign; % signed current in Amps
        % Isolate the channel with timestamps on it
        [m,~] = find(i ~= 0,1); % s is the stimulating channel
        i = i(m,:);
        % Loop through i and take only the timestamps at the beginning of each
        % stimulus
        time_stamps = flip(find(i < 0));
        dt = diff(time_stamps);
        kill = dt == -1;
        time_stamps(kill) = [];
        time_stamps = flip(time_stamps);
        clearvars i sign data % Save on working memory
        save('time_stamps.mat','time_stamps','m'); % Save time stamps to disk to save on computation time later
    else % time_stamps already calculated previously
        load('time_stamps.mat','time_stamps','m');
    end
    n_REP = length(time_stamps);
end
%% Associate time_stamps with trial types
% Preload variables to avoid annoying warnings
TrialParams = [];
datafile = dir('*_exp_datafile_*.mat');
if ~isempty(datafile)
    names = {datafile.name};
    if length(names) > 1
        disp('Please ensure you only have the one corresponding datafile in this subdirectory.');
        return;
    end
    n_REP = 0; % Sometimes I hate how MATLAB demands strings to load variables in from files, then can't parse them
    load(names{1,1},'AMP','DUR','n_REP','CHN','StimParams','n_Trials','TrialParams');
    if (n_Trials ~= nStamps)
        disp('Warning: Digital Lines do not match expected number of trials');
    end
    uniqueTrials = length(AMP)*length(DUR)*length(CHN);
    trial = input(['Please select a trial identifier between 1 and ' num2str(uniqueTrials) '.\n'],'s');
    allTrials = cell2mat(TrialParams(2:(uniqueTrials*n_REP)+1,2));
    selectedTrials = zeros(2,n_REP);
    selectedTrials(1,:) = find(allTrials == str2double(trial));
    % selectedTrials(:,28:end) = [];
    % n_REP = 27;
    selectedTrials(2,:) = time_stamps(1,selectedTrials(1,:));
end
%% PSTHs
if (t == 1)
    % Generate a MUA filter
    [Mufilt,~] = generate_Filters;
    MuNf = length(Mufilt);
    v_fid = fopen('amplifier.dat','r');
    BIN = [-1000, 1000];    % Extraction window for reading in data
    nData = diff([-100 200])*30+1;  % Number of datapoints
    threshfac = 4.2;       % Threshold factor for spike extraction. Increase to be harsher.
    SPIKES = cell(nChn,n_REP); % Will hold spike times for PSTH
    OUTPUT = cell(nChn,n_REP);
    X = 1/30:1/30:49/30;
    for r = 1:n_REPq
        OFFSET = cast(nChn*2*(FS/1e3)*(selectedTrials(2,r)+BIN(1)),'int64');
        fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
        v = fread(v_fid,[nChn,(FS/1e3)*6000],'int16') * 0.195;
        fseek(v_fid,0,'bof'); % Return to the start of the file
        %% Extract spikes from v. OUTPUT will contain n_REP cell arrays
        % Blank. Note that BlankArte always takes a window of [-500, 500].
        blank_v = v;
        blank_v(:,(501*30):(1500*30)) = BlankArte(v(:,(501*30):(1500*30)),nChn);
        % Downsample blank-v to a smaller window
        blank_v(:,999*30:1004*30) = 0;
        blank_v = blank_v(:,(901*30):(1200*30));
        % Filter
        for c = 16
            tmp = conv(blank_v(c,:),Mufilt);
            filt_v = tmp(1,MuNf/2:nData+MuNf/2-1);
            % Generate a threshold
            sd = median(abs(filt_v(200*30:299*30)))./0.6745;
            thresh = threshfac*sd;
            %thresh = 35;
            % Plot for debugging
            %X = 1/30:1/30:length(filt_v)/30;
            %figure;
            %plot(X,filt_v);
            %line([0 X(end)],[thresh thresh],'Color','red');
            % Extract spikes
            [WAVES,Sp_tmp] = spikeextract(filt_v,thresh,FS);
            SPIKES{c,r} = [Sp_tmp WAVES];
            OUTPUT{c,r} = Sp_tmp;
            % Plot the spike waveform
            PRE = figure(1);
            hold on
            POST = figure(2);
            hold on
            for sw = 1:length(Sp_tmp)
                if (Sp_tmp(sw) > 80) && (Sp_tmp(sw) < 101)
                    figure(PRE);
                    plot(X,WAVES(sw,:),'Color',[.7 .7 .7]);
                    xlabel('Time (msec)');
                    ylabel('Voltage (uV)');
                    title('20 msec before STIM');
                end
                if (Sp_tmp(sw) > 103) && (Sp_tmp(sw) < 121)
                    figure(POST);
                    plot(X,WAVES(sw,:),'Color',[.7 .7 .7]);
                    xlabel('Time (msec)');
                    ylabel('Voltage (uV)');
                    title('20 msec after STIM');
                    %axis([0 50 -60 40]);
                end
            end
        end
    end
    MEAN = zeros(1,49);
    count = 0;
    AMPLITUDE = cell(50,1);
    WIDTH = cell(50,1);
    for i = 1:50
        for n = 1:length(SPIKES{c,i}(:,1))
            if length(SPIKES{c,i}(n,2:end)) == 49
                if SPIKES{c,i}(n,1) > 103 && SPIKES{c,i}(n,1) < 121
                    MEAN = MEAN + SPIKES{c,i}(n,2:end);
                    count = count + 1;
%                     AMPLITUDE{i}(n) = max(SPIKES{c,i}(n,2:end));
%                     wdtmp(1) = find(SPIKES{c,i}(n,2:end) >= 0.5*AMPLITUDE{i}(n),1,'first')-1;
%                     wdtmp(2) = find(SPIKES{c,i}(n,2:end) >= 0.5*AMPLITUDE{i}(n),1,'last')+1;
%                     WIDTH{i}(n) = (wdtmp(2) - wdtmp(1))/30;
                end
            end
        end
    end
    MEAN = MEAN ./ count;
    figure(POST);
    hold on
    plot(X,MEAN,'Color','k','LineWidth',2)
    spikes = ['n = ' num2str(count)];
    line([0.0335 1.6334],[0 0],'Color','red');
    text(1.2, 250, spikes,'FontSize',16);
    axis([0 1.6 -40 60]);
    text(1.2, 50, spikes,'FontSize',16);
    MEAN = zeros(1,49);
    count = 0;
    for i = 1:50
        for n = 1:length(SPIKES{c,i}(:,1))
            if length(SPIKES{c,i}(n,2:end)) == 49
                if SPIKES{c,i}(n,1) > 80 && SPIKES{c,i}(n,1) < 101
                    MEAN = MEAN + SPIKES{c,i}(n,2:end);
                    count = count + 1;
                end
            end
        end
    end
    MEAN = MEAN ./ count;
    figure(PRE);
    hold on
    plot(X,MEAN,'Color','k','LineWidth',2)
    spikes = ['n = ' num2str(count)];
    line([0.0335 1.6334],[0 0],'Color','red');
    text(1.2, 30, spikes,'FontSize',16);
%     Plot the PSTH
%         stimCHN = find(cell2mat(TrialParams(2:(uniqueTrials*n_REP)+1,2)) == str2double(trial),1);
%         stimCHN = TrialParams{stimCHN+1,3};
%         channel = input(['Which channel would you like to plot?\nThe stimulating channel is: ' num2str(stimCHN) '\n'],'s');
%         channel = str2double(channel);
%     figure
%     hold on
%     for i = 1:50        
%         scatter(WIDTH{i},AMPLITUDE{i},'filled');
%     end
    for channel = 16
        SMOOTHING = 5;
        MAXRATE = 500; % Would require a spike every 2ms, which isn't feasible with real data. Therefore, will always overestimate.
        figure;
        psth(OUTPUT(channel,:),[-100,200],SMOOTHING,MAXRATE);
        if strcmp(trial,'1')
            title('Stimulation: 0 uA | Channel: 15');
        elseif strcmp(trial,'5')
            title('Stimulation: 20 uA | Channel: 15');
        elseif strcmp(trial,'6')
            title('Stimulation: 0 uA | Channel: 16');
        elseif strcmp(trial,'10')
            title('Stimulation: 20 uA | Channel: 16');
        end
        axis([-100 200 0 MAXRATE])
        xlabel('Time (msec)');
        ylabel('Trials (10*#)');
    end    
    cd(home)
end
%% Periodograms and Spectrograms
if (t == 5)
    %% Generate filters
    [~,Lfpfilt] = generate_Filters;
    LfpNf = length(Lfpfilt);
    %% Set up the channels we want to look at
    BIN = [-500, 1000];
    depth = Depth;
    depth = depth + 1;
    % Stimulating channel is m
    stim = m;
    % Adjacent channel is . . .
    stimchn = find(depth == stim);
    adja1 = depth(stimchn-1);
    adja2 = depth(stimchn+1);
    far = depth(32-stimchn);
    CHAN = [stim; adja1; adja2; far];
    %% Set up a loop through the actual data. There are n_time_stamps trials
    v_id = fopen('amplifier.dat', 'r');
    filt_lfp = zeros(length(CHAN),diff(BIN));
    d_v = zeros(length(CHAN),diff(BIN)*30);
    for n = TRIAL
        % Zoom to where the trial takes place
        bytes = (time_stamps(n)+BIN(1)*30)*2*nChn; % The number of bytes between bof and the timestamp. nChn*2 bytes per sample
        fseek(v_id,(bytes),'bof');
        v = fread(v_id, [nChn, diff(BIN)*30], 'int16') * 0.195;
        nData = size(v,2);
        if (nData) < diff(BIN)*30
            continue;
        end
        for c = 1:length(CHAN)
            % Detrend v
            d_v(c,:) = Detrend(v(CHAN(c),:),4);
            % Low pass filter
            tmp = conv(d_v(c,:),Lfpfilt);
            filt_lfp(c,:) = tmp(1,LfpNf/2:FS/1e3:nData+LfpNf/2-1);
        end
    end
    FS = 1e3;
    % Periodograms
    TITLE = {'Stimulating Channel','Adjacent Channel','Adjacent Channel','Distant Channel'};
    figure;
    hold on
    for c = 1:length(CHAN)
        xdft = fft(filt_lfp(c,:));
        xdft = xdft(1:floor((diff(BIN)+1)/2+1));
        psdx = (1/(FS*(diff(BIN)+1))) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        lfp_freq = 0:FS/(diff(BIN)+1):FS/2;
        subplot(2,2,c)
        plot(lfp_freq,10*log10(psdx));
        xlim([0 70]);
        xlabel('Frequency (Hz)');
        ylabel('Power (dB)');
        title(TITLE{c});
        grid on
    end
    figure;
    hold on
    for c = 1:length(CHAN)
        subplot(2,2,c)
        spectrogram(filt_lfp(c,:),'yaxis');
        xlabel('Time (msec)');
        title(TITLE{c});
    end
    cd(home)
end