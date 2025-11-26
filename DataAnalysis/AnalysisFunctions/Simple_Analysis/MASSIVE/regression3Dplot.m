function [bf]=regression3Dplot(avgnospT,stderrspktrial,startpointseconds, secondstoanalyse, AMPInterestSingleLinePlotchn1, AMPInterestSingleLinePlotchn2)
%Creates Heatmaps of dual electrode stimulation and a linecut at the
%amplitude of interest input into the function

trialinfo=loadTrialInfo(0);
lastwarn('', '');
loadNORECORDELECT;
[~, warnId] = lastwarn();
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
singleLineplotvalchn1=zeros(nChn,trialjump);
singleLineplotvalchn2=zeros(nChn,trialjump);
singleLineplotvalstd=zeros(nChn,trialjump);

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
             [~, AMPInterestSingleLinePlotINDEXDUALchn1]=min(abs(AMP_original(2:end)-AMPInterestSingleLinePlotchn1)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             AMPInterestSingleLinePlotINDEXDUALchn1=AMPInterestSingleLinePlotINDEXDUALchn1+1; %since you did not check -1 condition add one to the index AMP(2:end)
            [~, AMPInterestSingleLinePlotINDEXDUALchn2]=min(abs(AMP_original(2:end)-AMPInterestSingleLinePlotchn1)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             AMPInterestSingleLinePlotINDEXDUALchn2=AMPInterestSingleLinePlotINDEXDUALchn2+1; %since you did not check -1 condition add one to the index AMP(2:end)
             
         else
             chosen_trials=Desired_trialequal;
             checkconsecutive=1;
             %find closest AMP level to your chosen amp level to analyse (in case you
             %pick a level not tested)
             [~, AMPInterestSingleLinePlotINDEXDUALchn1]=min(abs(AMP(2:end)-AMPInterestSingleLinePlotchn1)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             AMPInterestSingleLinePlotINDEXDUALchn1=AMPInterestSingleLinePlotINDEXDUALchn1+1; %since you did not check -1 condition add one to the index AMP(2:end)
            [~, AMPInterestSingleLinePlotINDEXDUALchn2]=min(abs(AMP(2:end)-AMPInterestSingleLinePlotchn2)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             AMPInterestSingleLinePlotINDEXDUALchn2=AMPInterestSingleLinePlotINDEXDUALchn2+1; %since you did not check -1 condition add one to the index AMP(2:end)
             
         end
         chosen_trials(chosen_trials>endtrialelect*(group_related))=[]; % removes any trials that are greater than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         chosen_trials(chosen_trials<group_related/2)=[]; % removes any trials that are smaller than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         normalisedAvgspikingT=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,chosen_trials));%-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
         stdsp=stderrspktrial(:,chosen_trials); %finds the standard deviation of chosen trials for plotting
         loopcounter=loopcounter+1;
         singleLineplotvalchn1(:,loopcounter)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUALchn1);
         singleLineplotvalstd(:,loopcounter)=(1000/(secondstoanalyse-startpointseconds)).*stdsp(:,AMPInterestSingleLinePlotINDEXDUALchn1);

         singleLineplotvalchn2(:,loopcounter)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUALchn2);
         
         if cell2mat(trialinfo((chosen_trials(2)*2),2))==0 
             check=['T100_' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
             if TJ_related==1 && ~isempty(NORECORDELECT)

                 E1=num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)));
             elseif ~isempty(NORECORDELECT)
                E2=num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)));

             end
         else
             check=['T' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))*100/(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2),18))*100/((cell2mat(trialinfo((chosen_trials(2)*2)-1,18)))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
         end
         %DepthChangeingSpiking_SM(normalisedAvgspikingT, chosen_trials,fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp); %plots heat maps
        
         
         truedatastruct.(check) = normalisedAvgspikingT;
         %          if VarAmp==1 && cell2mat(trialinfo(TJ_related+1+(TJ_related-1),2))==0
%              xlim([0 AMP(end)/2]) % makes axes comparable on single stim trials
%          end
     end
%          figure %plotting real data line plot
%          hold on
%          scatter3(singleLineplotvalchn1(:,1),singleLineplotvalchn1(:,5),singleLineplotvalchn1(:,3))
%          xlabel(['Channel ' E1 ' Sp/s'])
%          ylabel(['Channel ' E2 ' Sp/s'])
%          zlabel('E1+E2')
%          title('Chn1 AMP')
%          X_line=best_fit_3D_line([singleLineplotvalchn2(:,1),singleLineplotvalchn2(:,5),singleLineplotvalchn2(:,3)]);
%          plot3(X_line(:,1),X_line(:,2),X_line(:,3),'-r','LineWidth',3) % best fit line 

         
         figure 
         ax=scatter(singleLineplotvalchn2(:,1),singleLineplotvalchn1(:,5));
         bf = singleLineplotvalchn2(:,1)\singleLineplotvalchn1(:,5);
         ycalc=bf*singleLineplotvalchn2(:,1);
         hold on
         plot(singleLineplotvalchn2(:,1),ycalc)
         xlabel(['Channel ' E1 ' Sp/s'])
         ylabel(['Channel ' E2 ' Sp/s'])
         title(['chn ' E2 ' ' num2str(AMPInterestSingleLinePlotchn1) 'uA ' 'chn ' E1 ' ' num2str(AMPInterestSingleLinePlotchn2) 'uA'])
         xLimits = get(gca,'XLim');  %# Get the range of the x axis
         yLimits = get(gca,'YLim');  %# Get the range of the y axis
         if xLimits(2)>yLimits(2)
             upperlim=xLimits(2);
         else
             upperlim=yLimits(2);
         end
         if xLimits(1)>yLimits(1)
             lowerlim=xLimits(1);
         else
             lowerlim=yLimits(1);
         end
         xlim([lowerlim upperlim])
         ylim([lowerlim upperlim])
         loopcounter=0;
%          
%          if bf<1
%             points=singleLineplotvalchn2(:,1)-singleLineplotvalchn1(:,5);
%          else
%              points=singleLineplotvalchn1(:,5)-singleLineplotvalchn2(:,1);
%          end
%          figure
%          scatter(1:length(points),points)
%          x=1:24;
%          f = @(b,x) b(1).*exp(b(2).*(x.^2));                                     % Objective Function
%          B = fminsearch(@(b) (sum((points(1:24)' - f(b,x)).^2)), [0; 0]);                % Estimate Parameters
%          hold on
%          plot(x, f(B,x'), '-r')
%          [B, x, f]=equationregression;
         
         
         
end
end

