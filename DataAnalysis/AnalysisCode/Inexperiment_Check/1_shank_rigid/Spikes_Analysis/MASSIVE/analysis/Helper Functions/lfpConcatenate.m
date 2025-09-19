function lfpConcatenate
dbstop if error
%% This function concatenates twenty-one milliseconds of high-pass filtered data from all stimulation events
%% Variables
BIN = [-1000 2000];
chan = [1,3:8,10:32];
trialID = 4;
eC = 0.9;
%% Loading in data
d = Depth(3); trig = []; TrialParams = []; lfp = [];
loadTrig;
FREQ = mean(diff(trig));
loadLFP;
nChn = size(lfp,1);
loadTrialParams;
if ~isempty(TrialParams)
    TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == trialID));
else
    TrialParams = 1:length(trig);
end
trig = trig(TrialParams);
nTrig = length(trig);
SPACE = 1.2 * sqrt(nTrig);
%% Initialise data matrices
LFP = nan(nChn,length(chan),diff(BIN)+1);
mLFP = nan(length(chan),diff(BIN)+1);
sLFP = nan(length(chan),diff(BIN)+1);
figure('Position',[0,0,400,1200]); hold on;
X = BIN(1):1:BIN(2);
for c = chan
    for t = 1:nTrig
        try
%        if max(abs(lfp(d(c),trig(t)+BIN(1):trig(t)+BIN(2)))) <= 6*std(lfp(d(c),trig(t)+BIN(1):trig(t)+BIN(2)))
            LFP(c,t,:) = lfp(d(c),trig(t)+BIN(1):trig(t)+BIN(2));         
%        end
        catch
        end
    end
    sLFP(c,:) = nanstd(squeeze(LFP(c,:,:))) ./ sqrt(nTrig);
    % Remove DC offset
    mLFP(c,:) = nanmean(squeeze(LFP(c,:,:)));
%     m = nanmean(mLFP(c,1:100));
%     mLFP(c,:) = mLFP(c,:) - m;
end
for c = chan
    errorbar(X,mLFP(c,:)+(c-1)*SPACE,sLFP(c,:),'Color',[eC eC eC],'Linestyle','none');
end
for c = chan
    plot(X,mLFP(c,:)+(c-1)*SPACE,'Color','k');
    text(BIN(1)-45,(c-1)*SPACE,num2str(c));
end
line([0 0],[-100 (c-1)*SPACE+100],'Color','b');
line([FREQ FREQ],[-100 (c-1)*SPACE+100],'Color','b');
ylim([-(4*SPACE) (c)*SPACE+250]);
xlim([BIN(1)-50 BIN(2)+250]);
xlabel('Time (msec)');
set(gca,'FontSize',12);
yticks('');
xticks(BIN(1):100:BIN(2));
set(gca,'ycolor',[1 1 1 0]);
line([BIN(2)+50 BIN(2)+50],[((c)*SPACE)/2 (((c)*SPACE)/2)+200],'Color','k','LineWidth',3);
text(BIN(2)+80,(((c)*SPACE)/2)+100,'200 uV','FontSize',16);
line([100 100],[-100 (c-1)*SPACE+100],'Color','r','LineStyle','--');
%% Also plot a few single trials
if (length(chan) == 1) % Only do this for a single channel
    lfp = squeeze(LFP(c,:,:));
    nSam = 10; % Number of trials to plot
    r = randperm(nTrig,nSam);
    figure; hold on;
    for t = 1:nSam
        plot(lfp(r(t),:)+(t-1)*SPACE*sqrt(nTrig),'Color','k');
    end
    line([abs(BIN(1)) abs(BIN(1))],[-SPACE*sqrt(nTrig) t*SPACE*sqrt(nTrig)],'Color','b');
    line([abs(BIN(1))+FREQ abs(BIN(1))+FREQ],[-SPACE*sqrt(nTrig) t*SPACE*sqrt(nTrig)],'Color','b');
end
%% Save the results to a .mat file for CSD analysis
% Remove NaNs
mLFP(isnan(mLFP)) = 0;
save('CSD.mat', 'mLFP');
end