function ActivationDepth(AMPInterestSingleLinePlot,avgnospT,startpointseconds, secondstoanalyse,depthdriven,cutoffsp)


trialinfo=loadTrialInfo(0);
loadNORECORDELECT;
loadStimChn;
loadVarAmp;
loadNREP;
TrialParams=loadTrialParams;
loadAMP_all;
AMP=AMP_all;
AMP_original=loadAMP;
avgnostim=mean(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1),2);%average reponse without stimulation
filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
    originalEND=size(TrialParams,1);
    AMP=AMP_original;
else
    E_Mapnumber=loadMapNum;
    loadoriginalEND;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
 loopcounter=0;
 trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;%trials to jump over
 endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions 
 if isempty(endtrialelect)
     endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
 end
maxid=(originalEND-1)/(n_REP_true*2);
checkconsecutive=1;
singleLineplotval=zeros(nChn,trialjump);
depthactivated=zeros(6,7);
 for group_related=1:endtrialelect*2:maxid*2 %group_related is used to go through groups of related trials
     for TJ_related=1:2:trialjump*2 %goes through trials related by trial jump
         desiredchanneltrial_one=(find((cell2mat(trialinfo(:,2))==cell2mat(trialinfo(TJ_related+(group_related-1),2))))+1)/2; %finds trials with desired initial electrode
         desiredchanneltrial_two=find(cell2mat(trialinfo(:,2))==cell2mat(trialinfo(TJ_related+1+(group_related-1),2)))/2; %finds trials with desired second electrode
         chosen_trials=intersect(desiredchanneltrial_one,desiredchanneltrial_two); % finds the trials that intersect with two trials of interest
         Desired_trialequal=zeros(1,length(AMP));
         for Amploop=1:length(AMP) %used to identify existing 75/25 amplitudes for chn 2
             try
                 Desired_trialequal(Amploop)=(chosen_trials(cell2mat(trialinfo(chosen_trials*2,18))==AMP(Amploop)));%array of mathcning trial number
             catch
                 Desired_trialequal(Amploop)=0;
             end
         end
                  
         if all(chosen_trials<=maxid)&&(VarAmp==1)&&(~isempty(chosen_trials(diff(chosen_trials(chosen_trials<maxid))==1))) %deals with varying amplitude having the same two electrodes for all middle trials
             chosen_trials=chosen_trials(checkconsecutive:3:length(chosen_trials));
             checkconsecutive=checkconsecutive+1; %flag for middle trials
             %find closest AMP level to your chosen amp level to analyse (in case you
             %pick a level not tested)
             [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP_original(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
         else
             chosen_trials=Desired_trialequal;
             checkconsecutive=1;
             %find closest AMP level to your chosen amp level to analyse (in case you
             %pick a level not tested)
             [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
         end
         chosen_trials(chosen_trials>endtrialelect*(group_related))=[]; % removes any trials that are greater than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         chosen_trials(chosen_trials<group_related/2)=[]; % removes any trials that are smaller than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         normalisedAvgspikingT=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,chosen_trials)-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
         loopcounter=loopcounter+1;
         singleLineplotval(:,loopcounter)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUAL);
     end
         SMOOTHING=1;
         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
         if  (NORECORDELECT(1)==0) && (length(NORECORDELECT)==1)
             for p=1:1
                 rate = conv(singleLineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                 rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                 depthElect=depthdriven:-50:depthdriven-(nChn-1)*50;
                 depthactivated(1,3)=any(depthElect(rate>=cutoffsp)<150);% layer 1 %activated_Regions=find(rate>=cutoffsp)
                 depthactivated(2,3)=any(150<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=580);% layer 2/3
                 depthactivated(3,3)=any(580<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=800);% layer 4
                 depthactivated(4,3)=any(800<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=1130);% layer 5
                 depthactivated(5,3)=any(1130<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=1415);% layer 6a
                 depthactivated(6,3)=any(1415<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=1490);% layer 6b
             end
                
                stimdepth=(depthdriven-50*(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2))-1));
                Estim=cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2));
                Estimtwo=0;
                depthactivated(1,1)=any((stimdepth)<150);% layer 1 %activated_Regions=find(rate>=cutoffsp)
                 depthactivated(2,1)=any(150<(stimdepth) & (stimdepth)<=580);% layer 2/3
                 depthactivated(3,1)=any(580<(stimdepth) & (stimdepth)<=800);% layer 4
                 depthactivated(4,1)=any(800<(stimdepth) & (stimdepth)<=1130);% layer 5
                 depthactivated(5,1)=any(1130<(stimdepth) & (stimdepth)<=1415);% layer 6a
                 depthactivated(6,1)=any(1415<(stimdepth) & (stimdepth)<=1490);% layer 6b
                 depthactivated(1:6,2)=0;% no second elect
             else
             for p=1:size(singleLineplotval,2)
                 rate = conv(singleLineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                 rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                 depthElect=depthdriven:-50:depthdriven-(nChn-1)*50;
                 depthactivated(1,p+2)=any(depthElect(rate>=cutoffsp)<150);% layer 1 %activated_Regions=find(rate>=cutoffsp)
                 depthactivated(2,p+2)=any(150<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=580);% layer 2/3
                 depthactivated(3,p+2)=any(580<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=800);% layer 4
                 depthactivated(4,p+2)=any(800<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=1130);% layer 5
                 depthactivated(5,p+2)=any(1130<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=1415);% layer 6a
                 depthactivated(6,p+2)=any(1415<depthElect(rate>=cutoffsp) & depthElect(rate>=cutoffsp)<=1490);% layer 6b
             end
             

             if VarAmp==1
                 stimdepth=[(depthdriven-50*(cell2mat(trialinfo((chosen_trials(1)-2)*2,2))-1)) (depthdriven-50*(cell2mat(trialinfo((chosen_trials(1)-2)*2-1,2))-1))];
                 Estim=cell2mat(trialinfo((chosen_trials(1)-2)*2,2));
                 Estimtwo=cell2mat(trialinfo((chosen_trials(1)-2)*2-1,2));
             elseif (length(NORECORDELECT)==1) && (VarAmp==0)
                 stimdepth=[(depthdriven-50*(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2))-1)) (depthdriven-50*(cell2mat(trialinfo((chosen_trials(3)-2)*2,2))-1))];
                 Estim=cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2));
                 Estimtwo=cell2mat(trialinfo((chosen_trials(3)-2)*2,2));
             else
                 for looprecordelect=1:length(NORECORDELECT)+1
                     stimdepth(looprecordelect)=(depthdriven-50*((cell2mat(trialinfo((chosen_trials(1)-length(NORECORDELECT)+(looprecordelect-1))*2-1,2)))));
                 end
                 Estim=(cell2mat(trialinfo((chosen_trials(1)-length(NORECORDELECT))*2-1,2)));
                 Estimtwo=(cell2mat(trialinfo((chosen_trials(1)-length(NORECORDELECT)+1)*2-1,2)));
             end

             
             depthactivated(1,1)=any((stimdepth(1))<150);% layer 1
             depthactivated(2,1)=any(150<(stimdepth(1)) & (stimdepth(1))<=580);% layer 2/3
             depthactivated(3,1)=any(580<(stimdepth(1)) & (stimdepth(1))<=800);% layer 4
             depthactivated(4,1)=any(800<(stimdepth(1)) & (stimdepth(1))<=1130);% layer 5
             depthactivated(5,1)=any(1130<(stimdepth(1)) & (stimdepth(1))<=1415);% layer 6a
             depthactivated(6,1)=any(1415<(stimdepth(1)) & (stimdepth(1))<=1490);% layer 6b
             depthactivated(1,2)=any((stimdepth(2))<150);% layer 1
             depthactivated(2,2)=any(150<(stimdepth(2)) & (stimdepth(2))<=580);% layer 2/3
             depthactivated(3,2)=any(580<(stimdepth(2)) & (stimdepth(2))<=800);% layer 4
             depthactivated(4,2)=any(800<(stimdepth(2)) & (stimdepth(2))<=1130);% layer 5
             depthactivated(5,2)=any(1130<(stimdepth(2)) & (stimdepth(2))<=1415);% layer 6a
             depthactivated(6,2)=any(1415<(stimdepth(2)) & (stimdepth(2))<=1490);% layer 6b
             
         end
         loopcounter=0;
         depthtitles=[{'E1'} {'E2'} {'E1_100'} {'E1_75'} {'50_50'} {'E2_75'} {'E2_100'}];
         layertitles=[{'layer'}; {'L1'}; {'L2/3'}; {'L4'}; {'L5'}; {'L6A'}; {'L6B'}];

         depthactivated=[depthtitles; num2cell(depthactivated)];
         depthactivated=[layertitles depthactivated];
         
         file=['ActivatedDepth_' num2str(Estim) '_' num2str(Estimtwo) '.mat'];
         save(file,'depthactivated')
         depthactivated=zeros(6,7);
 end

end

