function truedatastruct=TrueData_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,depthdriven)
%Creates Heatmaps of dual electrode stimulation and a linecut at the
%amplitude of interest input into the function
plotdata=0;
trialinfo=loadTrialInfo(0);
lastwarn('', '');
loadNORECORDELECT;
[warnMsg, warnId] = lastwarn();
if(~isempty(warnId))
    NORECORDELECT=[];
end


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
truedatastruct=[];
 loopcounter=0;
 trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;%trials to jump over
 endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions 
 if isempty(endtrialelect)
     endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
 end
maxid=(originalEND-1)/(n_REP_true*2);
checkconsecutive=1;
singleLineplotval=zeros(nChn,trialjump);
singleLineplotvalstd=zeros(nChn,trialjump);

for group_related=1:endtrialelect*2:maxid*2 %group_related is used to go through groups of related trials
    if plotdata==1
    figure
    fignum=gcf;
    fignum=fignum.Number;
    end
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
         stdsp=stderrspktrial(:,chosen_trials); %finds the standard deviation of chosen trials for plotting
         loopcounter=loopcounter+1;
         singleLineplotval(:,loopcounter)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUAL);
         singleLineplotvalstd(:,loopcounter)=(1000/(secondstoanalyse-startpointseconds)).*stdsp(:,AMPInterestSingleLinePlotINDEXDUAL);

         
         if length(chosen_trials)==length(AMP)
             chosen_trials_amp=AMP';
         elseif length(chosen_trials)==length(AMP_original)
             chosen_trials_amp=AMP_original;
         end
         if cell2mat(trialinfo((chosen_trials(2)*2),2))==0 
             check=['T100_' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
             if TJ_related==1 && ~isempty(NORECORDELECT)
                 subplot(3,2,1)
             elseif ~isempty(NORECORDELECT)
                 subplot(3,2,2)
             end
         else
             check=['T' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))*100/(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2),18))*100/((cell2mat(trialinfo((chosen_trials(2)*2)-1,18)))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
             subplot(2,3,loopcounter+2)
         end
         if plotdata==1
             DepthChangeingSpiking_SM(normalisedAvgspikingT, chosen_trials,fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp); %plots heat maps
             
             yline(1450-(16-1)*50,'k')
             yline(1450-(32-1)*50,'k')
             yline(1450-(32+16-1)*50,'k')
         end
         truedatastruct.(check) = normalisedAvgspikingT;
         %          if VarAmp==1 && cell2mat(trialinfo(TJ_related+1+(TJ_related-1),2))==0
%              xlim([0 AMP(end)/2]) % makes axes comparable on single stim trials
%          end
     end
     if plotdata==1
         figure %plotting real data line plot
         hold on
         SMOOTHING=1;
         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
         if (length(NORECORDELECT)==1)&& (NORECORDELECT(1)==0)
             for p=1:1
                 rate = conv(singleLineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                 rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                 errorbar(rate,singleLineplotvalstd(:,p))
             end
                 xline(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)),'k')
                 legend('100/0','Stim E1')
                 title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)))])
         else
             for p=1:size(singleLineplotval,2)
                 rate = conv(singleLineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                 rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                 errorbar(rate,singleLineplotvalstd(:,p))
             end
             if VarAmp==1
                 xline(cell2mat(trialinfo((chosen_trials(1)-2)*2,2)),'k')
                 xline(cell2mat(trialinfo((chosen_trials(1)-2)*2-1,2)),'-.k')
                 legend('100/0','75/25','50/50','25/75','0/100', 'Stim E1','Stim E2')
                 title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(1)-2)*2,2))) ' & ' num2str(cell2mat(trialinfo((chosen_trials(1)-2)*2-1,2)))])
             elseif (length(NORECORDELECT)==1) && (VarAmp==0) 
                 legend('50/50','0/100','100/0','Stim E1','Stim E2')
                 xline(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)),'-.k')
                 xline(cell2mat(trialinfo((chosen_trials(3)-2)*2,2)),'k')
                 title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2,2))) ' & ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)))])
             else
                 for looprecordelect=1:length(NORECORDELECT)
                     lgnd_names{looprecordelect}=['50/50 w/ Stimchn ' num2str(cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(looprecordelect-1))*2,2))) ', ' num2str(cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(looprecordelect-1))*2-1,2)))];
                 end
                 for looprecordelect=1:(length(NORECORDELECT)+1)
                     lgnd_names{(looprecordelect+length(NORECORDELECT))}=['100/0 w/ Stimchn ' num2str(cell2mat(trialinfo((chosen_trials(1)-length(NORECORDELECT)+(looprecordelect-1))*2-1,2)))];
                 end
                 linS = {'-','--',':','-.'};
                 for looprecordelect=1:length(NORECORDELECT)+1
                     lgnd_names{(looprecordelect+length(NORECORDELECT)*2+1)}=['Stim E' num2str((looprecordelect))];
                     xline((cell2mat(trialinfo((chosen_trials(1)-length(NORECORDELECT)+(looprecordelect-1))*2-1,2))),'LineStyle', linS{looprecordelect}, 'Color', 'k')
                 end
                 title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA'])
                 legend(lgnd_names)
             end
         end
         ylabel('Sp/s')
         xlabel('Depth (um)')
         xlim([1 nChn])
         xticklabel_depth=(depthdriven-50*4):-50*5:depthdriven-(nChn-1)*50;
         xticklabels(cellstr(string((xticklabel_depth))));
         loopcounter=0;
     end
end
 if plotdata==1
 figure(fignum)
 ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
 ax.Label.String='Sp/s';
 ax.Label.Rotation=270;
 end
 save('truedatastruct.mat','truedatastruct')
end

