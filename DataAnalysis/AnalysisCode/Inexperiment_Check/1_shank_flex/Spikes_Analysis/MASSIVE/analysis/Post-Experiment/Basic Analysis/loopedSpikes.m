function loopedSpikes(filepath,threshsign,ARTIFACT_CUTOFF,BIN,spWINDOW,savepath)
%% Function takes in a filepath and experimental datafile, prompts for Trial
%% ID and outputs two plots - one of spikes before stimulation, and one of spikes after
%% Initialisation
% Deal with input error
if ~(nargin)
    disp('I need a filepath!');
    return;
end
if nargin < 2
    threshsign = 1;
end
if nargin < 3
    ARTIFACT_CUTOFF = 300;
end
if nargin < 4
    BIN = [-100 200];
end
if nargin < 5
    spWINDOW = [1 11];
end
if nargin < 6
    savepath = 'C:\Users\tall0003\Google Drive\Monash PhD\Data\';
end
% Load in the Intan data
generateIntanData;
% Load in the digital lines
generateTimeStamps;
% Load in the datafile and select the Trial ID for analysis
StimParams = [];AMP = [];DUR = [];CHN = [];
load([filepath, datafile],'AMP','DUR','n_REP','CHN','StimParams','n_Trials','TrialParams');
if (n_Trials ~= nStamps)
    disp('Warning: Digital Lines do not match expected number of trials');
    %return;
end
uniqueTrials = length(AMP)*length(DUR)*length(CHN);
% Build out the filters
generateFiltersandThresholds;
% Setup basic graph logic
setupGraphsBasic;
X = 1/30:1/30:49/30;
%% Loop logic
for LOOP = 1:uniqueTrials
    % Administration
    disp('Spikes Event loop.');
    disp(['AMP: ' num2str(AMP(1)) ':' num2str(AMP(end))]);
    disp(['DUR: ' num2str(DUR(1)) ':' num2str(DUR(end))]);
    disp(['CHN: ' num2str(CHN(1)) '|' num2str(CHN(end))]);
    allTrials = cell2mat(TrialParams(2:(uniqueTrials*n_REP)+1,2));
    selectedTrials = zeros(2,n_REP);
    selectedTrials(1,:) = find(allTrials == LOOP);
    selectedTrials(2,:) = time_stamps(1,selectedTrials(1,:));
    STIMCHN = StimParams{selectedTrials(1,1)+1,1};
    if str2double(STIMCHN(4:5)) < 10
        STIMCHN = STIMCHN(5);
    else
        STIMCHN = STIMCHN(4:5);
    end
    disp(['Display Spikes for Trial ID n = ' num2str(LOOP)]);
    thisAmp = num2str(StimParams{selectedTrials(1,1)+1,16});
    if length(thisAmp) == 1
        thisAmp = ['00' thisAmp];
    elseif length(thisAmp) == 2
        thisAmp = ['0' thisAmp];
    end
    thisDur = num2str(StimParams{selectedTrials(1,1)+1,13});
    thisChn = StimParams{selectedTrials(1,1)+1,1};
    disp(['AMP: ' thisAmp]);
    disp(['DUR: ' thisDur]);
    disp(['CHN: ' thisChn]);
    % Administration
    disp('Display spikes surrounding all repeats');
    POST = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    % Loading, filtering, and blanking.
    for r = 1:n_REP        
        OFFSET = cast(nChn*2*(FS/1e3)*(selectedTrials(2,r)+BIN(1)),'int64');
        fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
        v = fread(v_fid,[nChn,(FS/1e3)*diff(BIN)],'int16') * 0.195;
        %% Extract and plot spikes from v. OUTPUT will contain n_REP cell arrays
        for c = 1:nChn
            % Work out the subplot position
            INDEX = find(A' == c);
            % Blank the data
            chan = depth(c);
            if ~(strcmp(thisAmp,'0-1')) && ~(strcmp(thisDur,'0-1'))
                b_v = BlankArte(v(chan,:),BIN);
            else
                b_v = v(chan,:);
            end
            % Filter data
            tmp = conv(b_v,Mufilt);
            filt_v = tmp(1,MuNf/2:nData+MuNf/2-1);
            % Extract spikes
            [WAVES,Sp] = spikeextract(filt_v,thresh(threshsign,c),FS);                                               
            % Sort spikes to remove artefact outliers
            if sum(Sp) ~= 0
                for n = 1:length(Sp)
                    if max(WAVES(n,:)) > ARTIFACT_CUTOFF || min(WAVES(n,:)) < -ARTIFACT_CUTOFF
                        Sp(n) = NaN;
                    end
                    if Sp(n) < abs(BIN(1))+spWINDOW(1) || Sp(n) > abs(BIN(1))+spWINDOW(2)
                        Sp(n) = NaN;
                    end
                end
                WAVES(isnan(Sp),:) = [];
                Sp(isnan(Sp)) = [];
            else
                Sp(Sp == 0) = [];
            end
            % Window the spikes with the denoise function
            [WAVES,~] = windowSpikes(WAVES,Sp,threshsign);
            % Plot the spike waveform
            figure(POST);
            ax = subplot(ROW,COL,INDEX);
            hold on
            if ~isempty(WAVES)
                plot(X,WAVES');
            end
            if strcmp(STIMCHN,num2str(chan-1))
                box on
                ax.XColor = 'r';
                ax.YColor = 'r';
            end
            if find(XLABEL == INDEX,1)
                xlabel('Time (msec)');
            end
            if find(YLABEL == INDEX,1)
                ylabel('Voltage (uV)');
            end            
            title(['CHN ' num2str(chan-1) ' | Depth ' num2str(c)]);            
        end
        disp(['The current loop is: ' num2str(r)]);
    end
    saveSPIKEPLOTS;
end
% Exit gracefully
fclose(v_fid);
end
