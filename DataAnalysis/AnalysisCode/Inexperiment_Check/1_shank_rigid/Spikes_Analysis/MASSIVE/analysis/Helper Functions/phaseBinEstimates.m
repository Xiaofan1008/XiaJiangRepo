%% Tool for assessing reliability of phase bin estimation
function phaseBinEstimates(phaseOFFSET,stimulation,proPhaseEstimate)
dbstop if error
%% Setup
% Base directory
base = 'X:/Tim';
% Recordings included in the global analysis
listOfDirectories = {
                     %'\RAT0017\estim_pen1_003_190401_134749',3,-1,-1;...
                     %'\RAT0017\estim_pen1_004_190401_153246',3,-1,-1;...
                     %'\RAT0017\estim_pen1_009_190401_193112',3,-1,-1;...
                     '\RAT0017\estim_pen1_010_190401_201325',3,-1,-1;...
                     %'\RAT0017\estim_pen1_011_190401_204831',3,-1,-1;...
                     };
%% Basic Variables
BAND = [4,6,10,20,80]; % The frequency bands to include in analysis
chanOFFSET = 0; % The channel to assess phase on, relative to the spike channel
phaseBINS = 8; % Number of phase-bins to use
coloured = 1; % Controls the colormap of the histograms
switch nargin        
    case 2
        proPhaseEstimate = 0; % Controls whether to output a series of plots showing the modulation index
    case 1
        stimulation = 1; % Controls whether to assess baseline or stimulation events
        proPhaseEstimate = 0; % Controls whether to output a series of plots showing the modulation index
    case 0
        phaseOFFSET = [-100 -100 -100 -100 -25]; % The temporal offset from stimulation to measure phase
        stimulation = 1; % Controls whether to assess baseline or stimulation events
        proPhaseEstimate = 0; % Controls whether to output a series of plots showing the modulation index
end
%% Conversion
REC = listOfDirectories(:,1)';
pMap = listOfDirectories(:,2)';
ID = listOfDirectories(:,3)';
nT = listOfDirectories(:,4)';
%% Processing Logic
nRec = size(REC,2); gMean = cell(1,nRec);
for n = 1:nRec    
    cd([base REC{n}]);
    disp(['Recording: ' num2str(n)]);    
    gMean{n} = singleRecording(ID{n},BAND,phaseOFFSET,chanOFFSET,phaseBINS,pMap{n},nT{n},stimulation,proPhaseEstimate);
end
%% Collapsing Statistics
thisChnMean = zeros(nRec,length(BAND),phaseBINS);
nChan = 0;
for n = 1:nRec
    chnMean = gMean{n};   
    cd([base REC{n}]);
    load('collapse.mat','collapse');
    collapse(collapse(:,1) ~= 12,:) = [];
    chk = cast(collapse(:,2),'logical');
    chnMap = collapse(chk,1);
    nChn = size(chnMap,1);
    thisMean = zeros(nChn,length(BAND),phaseBINS);   
    for c = 1:nChn
        thisMean(c,:,:) = chnMean{chnMap(c)};
        nChan = nChan + 1;
    end    
    thisChnMean(n,:,:) = reshape(nansum(thisMean,1),[5 8]);        
end
gMean = reshape(nansum(thisChnMean,1),[5 8]);
%% Final Statistics

%% Output Plots
%% Plotting
C = cell(1,8);C{1} = 'b';C{2} = 'r';C{3} = 'g';C{4} = 'y';
C{5} = 'c';C{6} = 'm';C{7} = [0.5 0.5 0.5];C{8} = 'k';
%% Plot MEAN
MAX = zeros(length(BAND),1);
figure('Position',[0 0 1920 1080]); hold on;
for b = 1:length(BAND)
    subplot(1,length(BAND),b); hold on;
    if (coloured)
        for i = 1:phaseBINS
            bar(i,100*gMean(b,i)/sum(gMean(b,:)),'FaceColor',C{i});
        end
    else
        bar(1:phaseBINS,100*gMean(b,i)/sum(gMean(b,:))); %#ok<*UNRCH>
    end       
    set(gca,'FontSize',10);
    title([num2str(BAND(b)) ' Hz'],'FontSize',24);
    xticks(0.5:1:phaseBINS+0.5)
    XTL = cell(1,phaseBINS+1);
    for i = 1:phaseBINS+1
        XTL{i} = [num2str((i-1)*45) char(176)];
    end
    xticklabels(XTL);
    if (b == 1)
        ylabel('Percentage Occurrence of Phasebins (%)','FontSize',24);
    end
    if (b == 3)
        xlabel('Centre Phase Angle (Deg)','FontSize',24);
    end
    xlim([0 phaseBINS+1]);
    MAX(b) = max(100*gMean(b,:)/sum(gMean(b,:)));
end
if (stimulation)
    sgtitle([num2str(nRec) ' Recordings | ' num2str(nChan) ' Channels | Stimulation Condition'],'FontSize',30);
else
    sgtitle([num2str(nRec) ' Recordings | ' num2str(nChan) ' Channels | Baseline Condition'],'FontSize',30);
end
for b = 1:length(BAND)
    subplot(1,length(BAND),b);
    ylim([0 max(MAX)*1.1]);
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
end
end

function chnMean = singleRecording(ID,BAND,phaseOFFSET,chanOFFSET,phaseBINS,pMap,nTrig,stimulation,pPE)
%% Set up
% Load in the collapse matrix
if isempty(dir('collapse.mat'))
    loopedResponseCurve;
end
load('collapse.mat','collapse');
chk = cast(collapse(:,2),'logical');
chnMap = collapse(chk,1);
%loadStimChn; chnMap = stimChn;
nChn = size(chnMap,1);
% Set up the depth matrix
if (pMap > 2)
    depth = Depth(pMap);
else
    depth = Depth(pMap)+1;
end
%% Loading in data
lfp = loadLFP; trig = loadTrig; TrialParams = loadTrialParams;
%% Processing Logic
if (stimulation)
    if (ID ~= -1)
        thisTrial = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == ID));
        trig = cast(trig(thisTrial),'int32');
        if nTrig ~= -1
            trig = trig(1:nTrig);
        end
    end    
end
chnMean = cell(nChn,1);
for n = 1:nChn
    if (chnMap(n) + chanOFFSET) > 0 && (chnMap(n) + chanOFFSET) <= size(lfp,1)
        thisChn = depth(chnMap(n) + chanOFFSET);
    else
        thisChn = depth(chnMap(n));
    end
    thisLfp = lfp(thisChn,:);
    disp(['Channel: ' num2str(chnMap(n))]);
    chnMean{n} = singleChannel(thisLfp,trig,BAND,stimulation,phaseOFFSET,phaseBINS,pPE);    
end
end

function mean = singleChannel(lfp,trig,BAND,stimulation,phaseOFFSET,phaseBINS,pPE)
%% Processing Logic
[phase] = generatePhaseVector(lfp,BAND);
if (stimulation)
    phasebin_vector = zeros(length(BAND),length(trig));
    for b = 1:length(BAND)
        for t = 1:length(trig)
            phasebin_vector(b,t) = ceil(phase{b}(trig(t)+phaseOFFSET(b))./(360/phaseBINS));            
        end
        phasebin_vector(b,phasebin_vector(b,:) <= 0) = phasebin_vector(b,phasebin_vector(b,:) <= 0) + phaseBINS;
    end
else
    phasebin_vector = zeros(length(BAND),length(lfp));
    for b = 1:length(BAND)
        phasebin_vector(b,:) = ceil(phase{b}./(360/phaseBINS));
        phasebin_vector(b,phasebin_vector(b,:) <= 0) = phasebin_vector(b,phasebin_vector(b,:) <= 0) + phaseBINS;
    end
end
%% Initial Statistics
mean = zeros(length(BAND),phaseBINS);
for b = 1:length(BAND)
    for c = 1:phaseBINS
        mean(b,c) = sum(phasebin_vector(b,:) == c);
    end
end

if (pPE)
    progressivePhaseEstimate(phase,trig,BAND,phaseBINS);
end
end

function progressivePhaseEstimate(phase,trig,BAND,phaseBINS)
% Produces a plot of the moduulation index for a moving window relative to
% stimulation
%% Variables
DUR = 2e3; % Amount of time to analyse
BN = ceil(2e3.*(1./BAND)); % Bin size in msec
%% Logic
WINDOW = cell(1,length(BAND));
for b = 1:length(BAND)
    WINDOW{b} = floor(-DUR/2 + BN(b)/2):ceil(BN(b)/10):ceil(DUR/2 - BN(b)/2);
end
phasebin_vector = zeros(length(BAND),length(phase{1}));
for b = 1:length(BAND)
    phasebin_vector(b,:) = ceil(phase{b}./(360/phaseBINS));
    phasebin_vector(b,phasebin_vector(b,:) <= 0) = phasebin_vector(b,phasebin_vector(b,:) <= 0) + phaseBINS;
end
progPhase = cell(1,length(BAND));
for b = 1:length(BAND)
    progPhase{b} = NaN(1,length(trig),BN(b)+1);
    for w = 1:length(WINDOW{b})
        for t = 1:length(trig)
            progPhase{b}(w,t,:) = phasebin_vector(b,trig(t)+WINDOW{b}(w) - BN(b)/2:trig(t)+WINDOW{b}(w) + BN(b)/2);
        end
    end
end
%% Collapse across the time dimension
mean = cell(1,length(BAND));
for b = 1:length(BAND)   
    mean{b} = zeros(phaseBINS,length(WINDOW{b}));
    for w = 1:length(WINDOW{b})
        for c = 1:phaseBINS
            for t = 1:BN(b)+1
                mean{b}(c,w) = mean{b}(c,w) + sum(progPhase{b}(w,:,t) == c);
            end
        end
    end
end
%% Collapse phasebins counts to a modulation index
MI = cell(1,length(BAND));
for b = 1:length(BAND)
    MI{b} = zeros(1,length(WINDOW{b}));
    for w = 1:length(WINDOW{b})
        mean{b}(:,w) = 100 .* mean{b}(:,w) ./ sum(mean{b}(:,w));
        MI{b}(w) = max(mean{b}(:,w)) - min(mean{b}(:,w));
    end
end
%% Plot the data
figure; hold on;
for b = 1:length(BAND)    
    subplot(5,1,b); hold on;
    plot(WINDOW{b},MI{b}(1,:),'Color','k')
    grid on    
    xticks(WINDOW{b}(1):50:WINDOW{b}(end))
    xlabel('Time (ms)');
    ylabel('MI (% Difference)');
    line([0 0],[-100 100],'Color','r'); ylim([0 5]);
    title([num2str(BAND(b)) ' Hz']);
end
sgtitle('Modulation Index over time in each frequency band');
drawnow;
end