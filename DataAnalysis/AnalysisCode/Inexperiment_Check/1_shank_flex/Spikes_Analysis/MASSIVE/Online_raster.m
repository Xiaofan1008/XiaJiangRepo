function resp=Online_raster(ID,Spike_array)

%sp = Spike_array;
trig = loadTrig(0);
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
TP = loadTrialParams;
tID = find(cell2mat(TP(numsim_elect:numsim_elect:end,2)) == ID);
theseTrig = trig(tID)./30;
nT=length(theseTrig);
%% Set up the raster data structure
BIN = [-99 99]; 
SMOOTHING = 2; MAX = 400;
xdata = [];
ydata = [];


for tr = 1:nT
    theseSp = (Spike_array(Spike_array(:,1) > theseTrig(tr)+BIN(1) & Spike_array(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
    for i = 1:length(theseSp)
        xdata = [xdata, (theseSp(i))]; %#ok<*AGROW>
        %ydata = [ydata, tr*(MAX/nT)];
    end
end


    if size(xdata,2)>5
    % Add the convolved spikerate
    Z = hist(xdata,BIN(1):BIN(2)); %#ok<HIST>
    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
    rate = (1000/nT)*conv((Z),window);
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
    meanbase=mean(rate(1:abs(ceil(BIN(1)/2))));
    sd=std(rate(1:abs(ceil(BIN(1)/2))));
    resp=find(rate>meanbase+5*sd & rate>5 & sum(rate>meanbase+5*sd)>1);%5sp/s is the smoothed rate for a single spike
    resp(resp>190)=[];
    %if ~isempty(resp)
        MAX=max(rate);
        for tr = 1:nT
            theseSp = (Spike_array(Spike_array(:,1) > theseTrig(tr)+BIN(1) & Spike_array(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
            for i = 1:length(theseSp)    
                ydata = [ydata, tr*(MAX/nT)];
            end
        end
   figure; hold on; axM = gca;
        yScale = MAX/nT;%MAX./nT;
        for i = 1:size(xdata,2)
            line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)+yScale/2.3],'Color','k','LineWidth',2);
        end
        line([0 0],[0 MAX],'Color','b');
        
        plot(axM,BIN(1):BIN(2),rate,'b','LineWidth',2);

    %ylim(axM,[0 MAX]);
    % Tidy Up
    xlabel(axM,'Time (ms)');
    ylabel(axM,'Firing Rate (Sp/s)');
    %xticks(axM,150:50:300);
    %xticklabels(axM,0:50:400);
    beautifyPlot(30,axM);
   % end
    else
        resp=[];
    end
end