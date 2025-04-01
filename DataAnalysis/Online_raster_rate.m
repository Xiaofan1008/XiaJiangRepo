function [rate,xdata,ydata] = Online_raster_rate(ID,Spike_array)

sp = Spike_array;
trig = loadTrig(0);
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
TP = loadTrialParams;
tID = find(cell2mat(TP(numsim_elect:numsim_elect:end,2)) == ID);
theseTrig = trig(tID)./30;
nT=length(theseTrig);
%% Set up the raster data structure
BIN = [-200 200]; 
SMOOTHING = 2; MAX = 400;
xdata = {};
ydata = {};
xdatahold = [];
ydatahold = [];

for tr = 1:nT

    theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
    for i = 1:length(theseSp)
        xdatahold = [xdatahold, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
        ydatahold = [ydatahold, tr*(MAX/nT)];
    end
    xdata{end+1} = xdatahold;
    ydata{end+1} = ydatahold;
    xdatahold = [];
    ydatahold = [];
end

    % figure; hold on; axM = gca;
    % yScale = MAX./nT;
    % for i = 1:size(xdata,2)
    %     line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)+yScale/2.3],'Color','k','LineWidth',2);
    % end
    % line([200 200],[0 MAX],'Color','b');

    % Add the convolved spikerate
    Z = hist(xdatahold,0:400); %#ok<HIST>
    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
    rate = (1000/nT)*conv(Z,window);
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
    % plot(axM,rate,'b','LineWidth',2);
    % ylim(axM,[0 MAX]); xlim(axM,[150 300]);
    % % Tidy Up
    % xlabel(axM,'Time (ms)');
    % ylabel(axM,'Firing Rate (Sp/s)');
    % xticks(axM,150:50:300);
    % xticklabels(axM,-50:50:100);
    % beautifyPlot(30,axM);
end