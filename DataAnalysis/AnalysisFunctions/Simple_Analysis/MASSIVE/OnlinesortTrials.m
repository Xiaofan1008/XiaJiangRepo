function [Spike_trialstruct,baslinespike_trialstruct] = OnlinesortTrials(trig,sp,chn,startpointms,mstoanalyse)
%Sort into trial IDs and plot example spikes
%OUTPUT - the structure containing spiking information for each trial/repeat where
%each cell relates to one trial ID. The array in this cell is in the format
%of (electrodes,trials/repeats) where the total number of spikes per
%repeat is calculated for the given time

%INPUT - starting point after trigegr to begin counting spikes
%(startpointms) in ms, end point after trigger to stop counting spikes
%(mstoanalyse) in ms, trigger information (trig)
FS=30000;
TrialParams = loadTrialParams;
TrialParams=cell2mat(TrialParams);
maxtid=max(TrialParams(:,2));
nospI=[];
Spike_trialstruct=[];
nChn=length(sp);
portnunelect=32;
for i=1:nChn/portnunelect
    C1 = intersect(sp{1+(portnunelect*(i-1))}(:,1),sp{2+(portnunelect*(i-1))}(:,1));%checks later electrodes
    C2 = intersect(sp{24+(portnunelect*(i-1))}(:,1),sp{6+(portnunelect*(i-1))}(:,1)); %early electrodes
    C3 = intersect(sp{21+(portnunelect*(i-1))}(:,1),sp{5+(portnunelect*(i-1))}(:,1)); %middle electrodes
    test1= intersect(C1,C2);
    test2= intersect(C2,C3);
    check=intersect(test1,test2);
    for chncount=1+(portnunelect*(i-1)):portnunelect*i
        [C,r,c]=(intersect(sp{chncount}(:,1),check));
        sp{chncount}(r,:)=[];
    end
end
dispstat('','init');
dispstat(sprintf('Working through trial: '),'keepthis','n');
% Sort into trial IDs
for tID=1:maxtid
    dispstat(sprintf('%d',tID));
    if isempty(TrialParams(TrialParams(1:end,2) == tID))
        fprintf('No trial ID: %d\n',tID)
    else
        TrialParamstID = find(TrialParams(1:end,2) == tID); %identifies trial row matching trial ID
        num_elect=min(diff(find(diff(TrialParamstID)~=1)));
        if num_elect>1
        TrialParamstID=TrialParamstID(num_elect:num_elect:end);
        TrialParamstID=TrialParamstID./num_elect;
        end
        trigtID = trig(TrialParamstID)./(FS/1000);
        trigtID(trigtID==-500/(FS/1000))=[];
        nTrig = length(trigtID);
        for indT=1:nTrig
            for chsp=1:1:length(chn)
                v = sp{chn(chsp)};
                if ~isempty(v)
                    spikedetailstrig=v(((v(:,1)>(trigtID(indT)+startpointms))&(v(:,1)<(trigtID(indT)+mstoanalyse))),:); %v(((v(:,1)>0)&(v(:,1)<20000)),:);
                    timems=mstoanalyse-startpointms;
                    avgtimebs=1;
                    baslinespiketrig=v(((v(:,1)>(trigtID(indT)-5-avgtimebs*timems))&(v(:,1)<(trigtID(indT)-5))),:); %v(((v(:,1)>0)&(v(:,1)<20000)),:);
                    flag=checkTrials(trig(indT),(avgtimebs*timems+5)/1000,mstoanalyse/1000,nChn,FS,chsp,spikedetailstrig,baslinespiketrig);
                if flag==1
                    spikedetailstrig=[];
                    baslinespiketrig=[];
                end

                    
                    if ~isempty(spikedetailstrig)
                        spike=size(spikedetailstrig,1);
                    else
                        spike=0;
                    end
                    if ~isempty(baslinespiketrig)
                        spikebaseline=size(baslinespiketrig,1)/avgtimebs;
                    else
                        spikebaseline=0;
                    end
                else
                    spikedetailstrig=[];
                    baslinespiketrig=[];
                end
                mchsp=chn(chsp)-1;
                cnum=['Chn_' num2str(mchsp)];
                IDnum=['ID_' num2str(tID)];
                tnum=['Trial_' num2str(indT)];
                if isfield(Spike_trialstruct,cnum)&& isfield(Spike_trialstruct.(cnum),(IDnum)) && isfield(Spike_trialstruct.(cnum).(IDnum),(indT))
                    Spike_trialstruct.(cnum).(IDnum).(tnum)=[Spike_trialstruct.cnum.IDnum.tnum; spikedetailstrig];
                    baslinespike_trialstruct.(cnum).(IDnum).(tnum)=[ baslinespike_trialstruct.cnum.IDnum.tnum;  baslinespiketrig];
                else
                    Spike_trialstruct.(cnum).(IDnum).(tnum)=spikedetailstrig;
                    baslinespike_trialstruct.(cnum).(IDnum).(tnum)=baslinespiketrig;
                end

                nospI(chsp,indT)=spike;
                basespI(chsp,indT)=spikebaseline;
            end
        end
        basespI=[];
        nospI=[];
    end
end

end

