function [resp,rate,xdata,ydata]=ANS_raster(ID,Spike_array)

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
BIN = [-90 90];
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
    temp=rate>meanbase+5*sd;
    resp=find(temp==1 & rate>5 & sum(temp)>1);%5sp/s is the smoothed rate for a single spike
    %resp(resp>190)=[];
    if ~isempty(resp)
        MAX=max(rate);
        for tr = 1:nT
            theseSp = (Spike_array(Spike_array(:,1) > theseTrig(tr)+BIN(1) & Spike_array(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
            for i = 1:length(theseSp)
                ydata = [ydata, tr*(MAX/nT)];
            end
        end
    end
else
    resp=[];
    rate=[];
end
end