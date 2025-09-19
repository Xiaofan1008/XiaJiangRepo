function genRaster(thisChn,trial,RATE)
dbstop if error
%% Variables
%% Loading
TrialParams = loadTrialParams; trig = loadTrig; sp = loadSpikes; AMP = loadAMP; loadNTrials; d = Depth; loadCHN;
%% Initialisation of data structures
if ~isempty(TrialParams)
    nTr = size(TrialParams,1)/n_REP;
else
    nTr = length(trig);
end
nT = zeros(n_REP,nTr);
TrialParams = cell2mat(TrialParams);
for i = 1:nTr
    if ~isempty(TrialParams)
        nT(:,i) = TrialParams(TrialParams(:,2) == i,1);        
    else
        nT(i,:) = i;
    end
end
chn = d(thisChn);
%sp = cast(sp{chn},'int32');
sp = sp{chn};
sp = denoiseSpikes(sp);
%% PSTH Methodology for window adjustment
SPIKES = cell(1,nTr);
figure; hold on;
MR = 250;
XLIM = [-20 150];
for r = trial
    for t = 1:n_REP
        thisSp = sp(sp(:,1) >= (trig(nT(t,r)) - 500) & sp(:,1) <= (trig(nT(t,r)) + 500),:);
        SPIKES{t} = thisSp(:,1) - (trig(nT(t,r)) - 500);
    end
    if r > 1
        psth(SPIKES(:),[-500 500],1,MR,[],RATE,MR*(r-2));
        line(XLIM,[MR*(r-2) MR*(r-2)],'Color','b');
    else
        psth(SPIKES(:),[-500 500],1,MR,[],RATE,MR*(r-1));        
        line(XLIM,[MR*(r-1) MR*(r-1)],'Color','b');
    end
end
xlim(XLIM);
set(gca,'YColor',[1 1 1]);
xlabel('Time from Stimulation (ms)')
beautifyPlot;
%title(['Raster Plot | Channel ' num2str(thisChn) ' | ' num2str(AMP(trial)) ' uA']);
line([3 3],[0 MR*(r-1)],'Color','r');
line([8 8],[0 MR*(r-1)],'Color','r');
patch([-2.99 2.99 2.99 -2.99],[0.4 0.4 MR*(r-1) MR*(r-1)],[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end