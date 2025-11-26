%% Function takes in a filepath and returns LFP plots
function loopedLFP(filepath,savepath)
%% Deal with input error
if ~(nargin)
    disp('I need a filepath!');
    return;
end
if nargin < 2
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
X = BIN(1):BIN(2);
%% Loop logic
for LOOP = 1:uniqueTrials
    % Administration    
    disp('LFP Event loop.');
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
    disp(['Display LFPs for Trial ID n = ' num2str(LOOP)]);    
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
    LFP = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on
    % Loading, filtering, and blanking.
    lfp = cell(1,nChn);
    MAX = 0;
    for c = 1:nChn
        lfp{c} = zeros(1,diff(BIN)+1);
    end   
    for r = 1:n_REP
        OFFSET = cast(nChn*2*(FS/1e3)*(selectedTrials(2,r)+BIN(1)),'int64');
        fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
        v = fread(v_fid,[nChn,(FS/1e3)*diff(BIN)],'int16') * 0.195;
        %% Extract spikes from v. OUTPUT will contain n_REP cell arrays        
        for c = 1:nChn            
            % Blank and detrend the artefact.
            chan = depth(c);
            if ~(strcmp(thisAmp,'0-1')) && ~(strcmp(thisDur,'0-1'))
                b_v = BlankArte(v(chan,:));
            else
                b_v = v(chan,:);
            end            
            % Filter data
            tmp = conv(b_v,Lfpfilt);
            lfp{c} = lfp{c} + tmp(1,LfpNf/2:FS/1e3:nData+LfpNf/2-1);
            if (max(lfp{c}) > MAX) && (r == 1)
                MAX = max(lfp{c});
            end
        end
    end
    % VERY IMPORTANT. FILTERING JUST ADDS EACH TRIAL TOGETHER
    for ch = 1:nChn
        lfp{ch} = lfp{ch} ./ n_REP;
    end
    figure(LFP);
    hold on
    SEP = MAX ./ (4); % Magic number to improve legibility
    for ch = 1:nChn        
        if ch == STIMCHN
            plot(X,lfp{ch}+(SEP*(STIMCHN-ch)),'Color','r','LineWidth',2);
        else
            plot(X,lfp{ch}+(SEP*(STIMCHN-ch)),'Color','k','LineWidth',2);
        end
    end
    set(gca,'FontSize',24);
    xlabel('Time (msec)','FontSize',36);
    ylabel('Voltage (uV) - with Offset','FontSize',36);
    yl = ylim;
    line([0 0],[yl(1) yl(2)],'Color','blue','LineWidth',2);
    title(['LFP - Trial ' num2str(LOOP) ' | Arranged by depth'],'FontSize',36);    
    ylim([SEP*(-34+STIMCHN) SEP*(STIMCHN+2)]);
    xlim([-5 120]);
    % Save the output
    saveLFPPLOTS;
end
% Exit gracefully
fclose(v_fid);
end