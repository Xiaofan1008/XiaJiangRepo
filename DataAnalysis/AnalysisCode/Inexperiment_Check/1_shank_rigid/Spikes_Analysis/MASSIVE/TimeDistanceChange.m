function TimeDistanceChange(TimeBegin,TimeEnd, TimeStep,starttrial, endtrial, depthdriven)
%for plotting distance and time of electrodes
trig = loadTrig(0);
theseTrig = trig;
trialinfo=loadTrialInfo;
trialinfo(1,:)=[];
trialjump= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition 
TrialParams = loadTrialParams;
maxtid=max(cell2mat(TrialParams(:,2)));
loadNREP;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=loadMapNum;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
try
    loadoriginalEND;
    if ((originalEND-1)/(n_REP_true*2))<maxtid
        maxtid=((originalEND-1)/(n_REP_true*2));
    end
catch
    %%code was made before loadoriginalEND
end
Loopnum=floor((TimeEnd-TimeBegin)/TimeStep);
loadStimChn;
x=zeros(Loopnum,nChn);
y=zeros(Loopnum,nChn);
z=zeros(Loopnum,nChn);
AMP=loadAMP;
overallcounter=0;
figure
for Ampchange=starttrial:trialjump:endtrial
    overallcounter=overallcounter+1;
for counter=1:Loopnum
        [IDstruct]=sortTrials_SM(TimeBegin+1+TimeStep*(counter-1),TimeStep*(counter)+TimeBegin,trig,0,starttrial,trialjump,endtrial);
        [avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct);
        z(counter,:)=(TimeStep*(counter)+TimeBegin).*ones(1,length(avgnospT(:,Ampchange)));
        x(counter,:)=depthdriven:-50:depthdriven-(nChn-1)*50;%1:length(avgnospT(:,Ampchange));
        y(counter,:)=(1000/(TimeStep-TimeBegin)).*avgnospT(:,Ampchange);
end
        subplot((length(AMP)-1)/2,2,overallcounter)
        hold on
        surf(z,x,y)
        yline(depthdriven-(cell2mat(trialinfo(Ampchange*2-1,2))-1)*50,'Color','r','LineWidth',2)
        if cell2mat(trialinfo(Ampchange*2,2))~=0
           yline(depthdriven-(cell2mat(trialinfo(Ampchange*2,2))-1)*50,'Color','r','LineWidth',2) 
        end
        zlabel('Sp/s')
        xlabel('Time (ms)')
        ylabel('Estimated Depth (um)')
        title([num2str(cell2mat(trialinfo(Ampchange*2,18))) 'uA'])
        caxis([0 200])
        ylim([depthdriven-(nChn-1)*50 depthdriven])
        xlim([TimeBegin+TimeStep TimeEnd])
        set(gca, 'YDir','reverse')
        
end
colorbar('Position',[0.93 0.1 0.03 0.85])
end