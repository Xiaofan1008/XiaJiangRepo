function [stimshankcentroid,truedatastruct,stackedacross]=TrueData_heatmapLinecutFOURSHANK(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,IDstruct,startpointseconds, secondstoanalyse,depthdriven)
%Creates Heatmaps of dual electrode stimulation and a linecut at the
%amplitude of interest input into the function
plottingnorm=0;
plotheat=0;
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
E_MAP = Depth(E_Mapnumber);

stimshankcentroid=zeros(20,1);
peakshank=zeros(20,1);
nostim=[];
trialnostim=find(cell2mat(trialinfo(1:2:end,18))==-1);
for Tnum=1:length(trialnostim)
    nostim=[nostim IDstruct.(['T', num2str(trialnostim(Tnum))])];
end
avgnostim=mean(nostim,2);%average reponse without stimulation
avgnostim=avgnostim(E_MAP);
varnostim=var(nostim,0,2)';%variance without stimulation
varnostim=varnostim(E_MAP);
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
         if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
             laminar=1;
         else
             laminar=0;
         end
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
         normalisedAvgspikingT=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,chosen_trials));%-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
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
             if plotheat==1
                 if TJ_related==1 && ~isempty(NORECORDELECT)
                     %subplot(3,2,1)
                     figure
                     fignum=gcf;
                     fignum=fignum.Number;
                 elseif ~isempty(NORECORDELECT)
                     %subplot(3,2,2)
                     figure
                     fignum=gcf;
                     fignum=fignum.Number;
                 end
             end
         else
             check=['T' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))*100/(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2),18))*100/((cell2mat(trialinfo((chosen_trials(2)*2)-1,18)))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
             %subplot(2,3,loopcounter+2)
             if plotheat==1
                 figure
                 fignum=gcf;
                 fignum=fignum.Number;
             end
         end
         if plotheat==1
             subplot(1,4,1)
             DepthChangeingSpiking_SM(normalisedAvgspikingT(1:16,:), chosen_trials,fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp); %plots heat maps
             if (cell2mat(trialinfo(chosen_trials(1)*2,2)))<17 && (cell2mat(trialinfo(chosen_trials(1)*2,2))>0)
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2,2))-1)*50,'r')
             end
             if (cell2mat(trialinfo(chosen_trials(1)*2-1,2)))<17 && (cell2mat(trialinfo(chosen_trials(1)*2-1,2))>0)
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2-1,2))-1)*50,'r')
             end
             subplot(1,4,4)
             DepthChangeingSpiking_SM(normalisedAvgspikingT(17:32,:), chosen_trials,fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp); %plots heat maps
             if ((cell2mat(trialinfo(chosen_trials(1)*2,2)))<33) && cell2mat(trialinfo(chosen_trials(1)*2,2))>16
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2,2))-16-1)*50,'r')
             end
             if ((cell2mat(trialinfo(chosen_trials(1)*2-1,2)))<33) && (cell2mat(trialinfo(chosen_trials(1)*2-1,2))>16)
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2-1,2))-16-1)*50,'r')
             end
             
             subplot(1,4,2)
             DepthChangeingSpiking_SM(normalisedAvgspikingT(33:48,:), chosen_trials,fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp); %plots heat maps
             if ((cell2mat(trialinfo(chosen_trials(1)*2,2)))<49) && cell2mat(trialinfo(chosen_trials(1)*2,2))>32
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2,2))-32-1)*50,'r')
             end
             if ((cell2mat(trialinfo(chosen_trials(1)*2-1,2)))<49) && (cell2mat(trialinfo(chosen_trials(1)*2-1,2))>32)
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2-1,2))-32-1)*50,'r')
             end
             subplot(1,4,3)
             DepthChangeingSpiking_SM(normalisedAvgspikingT(49:64,:), chosen_trials,fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp); %plots heat maps
             if ((cell2mat(trialinfo(chosen_trials(1)*2,2)))<65) && (cell2mat(trialinfo(chosen_trials(1)*2,2))>48)
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2,2))-48-1)*50,'r')
             end
             if ((cell2mat(trialinfo(chosen_trials(1)*2-1,2)))<65) && (cell2mat(trialinfo(chosen_trials(1)*2-1,2))>48)
                 yline(depthdriven-(cell2mat(trialinfo(chosen_trials(1)*2-1,2))-48-1)*50,'r')
             end
             ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
             ax.Label.String='Sp/s';
             ax.Label.Rotation=270;
         end
         truedatastruct.(check) = normalisedAvgspikingT;
         %          if VarAmp==1 && cell2mat(trialinfo(TJ_related+1+(TJ_related-1),2))==0
%              xlim([0 AMP(end)/2]) % makes axes comparable on single stim trials
%          end
     end
     electmean=zeros(1,nChn);
     electvar=zeros(1,nChn);
     chosenstimchn=1;
     for elect=1:nChn
         check=['T100_' num2str(CHN(chosenstimchn))];
         electrodeampint100=truedatastruct.(check)(elect,:);
         check=['T75_25_' num2str(CHN(chosenstimchn))];
         electrodeampint75=truedatastruct.(check)(elect,:);
         check=['T25_75_'  num2str(CHN(chosenstimchn))];
         electrodeampint25=truedatastruct.(check)(elect,:);
         check=['T50_50_'  num2str(CHN(chosenstimchn))];
         electrodeampint50=truedatastruct.(check)(elect,:);
         
         if laminar==1
             check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
         else
            check=['T100_' num2str(NORECORDELECT(1))];
         end
         electrodeampint0=truedatastruct.(check)(elect,:);
         electall=[electrodeampint100 electrodeampint75 electrodeampint25 electrodeampint50 electrodeampint0];
         electmean(elect)=mean(electall);
         electvar(elect)=var(electall);
     end
     
     if plotdata==1

         figure %plotting real data line plot
         hold on
         SMOOTHING=1;
         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
         stackedacross=[];
         acrossshankplot=zeros(5,4);
         if (length(NORECORDELECT)==1)&& (NORECORDELECT(1)==0)
             for p=1:1
                 subplot(1,4,1)
                 rate = conv(singleLineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                 rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                 errorbar(rate,singleLineplotvalstd(:,p))
             end
             xline(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)),'k')
             legend('100/0','Stim E1')
             title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)))])
         else
             for shankplot=1:4
                 if shankplot==1
                     subplot(1,4,1)
                     title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(1)])
                     if (cell2mat(trialinfo(2*2,2)))<17 && (cell2mat(trialinfo(2*2,2))>0)
                         yline((cell2mat(trialinfo(2*2,2))),'k')

                     end
                     if (cell2mat(trialinfo(2*2-1,2)))<17 && (cell2mat(trialinfo(2*2-1,2))>0)
                         yline((cell2mat(trialinfo(2*2-1,2))),'k')
                     end
                 elseif shankplot==2
                     subplot(1,4,4)
                     title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(4)])
                     if ((cell2mat(trialinfo(2*2,2)))<33) && cell2mat(trialinfo(2*2,2))>16
                         yline((cell2mat(trialinfo(2*2,2))-16),'k')

                     end
                     if ((cell2mat(trialinfo(2*2-1,2)))<33) && (cell2mat(trialinfo(2*2-1,2))>16)
                         yline((cell2mat(trialinfo(2*2-1,2))-16),'k')
                     end
                 elseif shankplot==3
                     subplot(1,4,shankplot-1)
                     title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(2)])
                     if ((cell2mat(trialinfo(2*2,2)))<49) && cell2mat(trialinfo(2*2,2))>32
                         yline((cell2mat(trialinfo(2*2,2))-32),'k')

                     end
                     if ((cell2mat(trialinfo(2*2-1,2)))<49) && (cell2mat(trialinfo(2*2-1,2))>32)
                         yline((cell2mat(trialinfo(2*2-1,2))-32),'k')
                     end
                     
                 else
                     subplot(1,4,shankplot-1)
                     title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(3)])
                     if ((cell2mat(trialinfo(2*2,2)))<65) && (cell2mat(trialinfo(2*2,2))>48)
                         yline((cell2mat(trialinfo(2*2,2))-48),'k')
                         stimplot=4;
                     end
                     if ((cell2mat(trialinfo(2*2-1,2)))<65) && (cell2mat(trialinfo(2*2-1,2))>48)
                         yline((cell2mat(trialinfo(2*2-1,2))-48),'k')
                     end
                 end
                 for p=1:size(singleLineplotval,2)
%                                           if p==2
%                                               ax = gca;
%                                               ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
%                                           end
                     if plottingnorm==1
                         maxplot=max(singleLineplotval,[],2);
                         normspk=singleLineplotval(1+((shankplot-1)*16):(shankplot*16),p)./maxplot(1+((shankplot-1)*16):(shankplot*16));%varnostim(1+((shankplot-1)*16):(shankplot*16))';
                         normspk(isnan(normspk))=0;
                         rate = conv(normspk,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                         rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                         acrossshankplot(p,shankplot)=mean(rate);
                         if p==5
                             stackedacross(:,shankplot)=singleLineplotval(1+((shankplot-1)*16):(shankplot*16),p);
                             itamount=0.5;
                         end
                         hold on
                         plot(rate,1:16)
                         rate(rate<0)=0;
                         A = trapz(1:16, rate);
                         B(1)=0;

                         for lims = 2:16
                             B(lims) =  trapz(1:lims, rate(1:lims));
                         end

                         [~,electrodecentroid]=min(abs(B-(A/2)));
                         stimshankcentroid(p+(shankplot-1)*5)=electrodecentroid;
                         %[~,electrodemax]=max(rate);
                         %stimshankcentroid(p+(shankplot-1)*5)=electrodemax;
                         xlabel('Normalised spike rate')
                     else
                         rate = conv(singleLineplotval(1+((shankplot-1)*16):(shankplot*16),p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                         xlabel('Sp/s')
                         rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                         acrossshankplot(p,shankplot)=rate(8);
                         if p==5
                             stackedacross(:,shankplot)=rate;
                             itamount=100;
                         end
                         hold on
                         errorbar(rate,1:16,singleLineplotvalstd(1+((shankplot-1)*16):(shankplot*16),p),'horizontal')
                     end
                 end
                 if shankplot==4
                     lgd = legend;
                     lgd.Position=[0.922,0.7036,0.06,0.149];
                     lgd.String={'Stim Elect','Stim Elect','0/100','25/75','50/50','75/25','100/0'};
                     %lgd.String={'Stim Elect','Stim Elect','0/100','50/50','100/0'};
                     %lgd.String={'Stim Elect','Stim Elect','25/75','75/25'};
                     %lgd.String={'Stim Elect','Stim Elect','25/75','50/50','75/25'};
                 end
                 %set(gca, 'YDir','reverse')
                 ylabel('Electrode')
%                  newcolors = {'#330000' '#990000', '#FF0000', '#FF6666', '#CBA8A8'};
%                  colororder(newcolors)
                 ylim([1 16])
             end
             
             
             figure
             hold on
             for i=1:16
                 plot(stackedacross(i,:)+(i-1)*itamount, 'k')
                 yline(i*itamount,'-.k')
             end
             if CHN(1)>48
                 scchn=CHN(1)-48;
                 shankk=3;
             elseif CHN(1)>32 && CHN(1)<49
                 scchn=CHN(1)-32;
                 shankk=2;
             elseif CHN(1)>16 && CHN(1)<33
                 scchn=CHN(1)-16;
                 shankk=4;
             else
                 scchn=CHN(1);
                 shankk=1;
             end
             scatter(shankk,scchn*itamount, 'r', 'filled')
             xlabel('shank')
             title('Average accross all electrodes on each shank. Dotted line is 50% max spkrate')
            
         end
         
         %              temp=singleLineplotval;
         %              singleLineplotvalN=singleLineplotval;
         %              singleLineplotvalN(49:64,:)=temp(17:32,:);
         %              singleLineplotvalN(17:32,:)=temp(33:48,:);
         %              singleLineplotvalN(33:48,:)=temp(49:64,:);
         %              shankspikes(1,:)=[mean(singleLineplotvalN(1:16,1)),mean(singleLineplotvalN(17:32,1)),mean(singleLineplotvalN(33:48,1)),mean(singleLineplotvalN(49:64,1))];
         %              shankspikes(2,:)=[mean(singleLineplotvalN(1:16,2)),mean(singleLineplotvalN(17:32,2)),mean(singleLineplotvalN(33:48,2)),mean(singleLineplotvalN(49:64,2))];
         %              shankspikes(3,:)=[mean(singleLineplotvalN(1:16,3)),mean(singleLineplotvalN(17:32,3)),mean(singleLineplotvalN(33:48,3)),mean(singleLineplotvalN(49:64,3))];
         %              shankspikes(1,:)=[mean(singleLineplotvalN(1:16,4)),mean(singleLineplotvalN(17:32,4)),mean(singleLineplotvalN(33:48,4)),mean(singleLineplotvalN(49:64,4))];
         %              shankspikes(5,:)=[mean(singleLineplotvalN(1:16,5)),mean(singleLineplotvalN(17:32,5)),mean(singleLineplotvalN(33:48,5)),mean(singleLineplotvalN(49:64,5))];
         %              hold on
         %              plot(1:4,shankspikes)
         %
         %              plot(acrossshankplot')
         %              temp=acrossshankplot;
         %              acrossshankplot(:,4)=temp(:,2);
         %              acrossshankplot(:,2)=temp(:,3);
         %              acrossshankplot(:,3)=temp(:,4);
         %              plot(acrossshankplot')
     end
end
if plotheat==1
    figure(fignum)
    ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
    ax.Label.String='Sp/s';
    ax.Label.Rotation=270;
end

%%no stim

%     check=['T100_' num2str(CHN(chosenstimchn))];
%     electrodeampint100ns=truedatastruct.(check)(:,1);
%     check=['T75_25_' num2str(CHN(chosenstimchn))];
%     electrodeampint75ns=truedatastruct.(check)(:,1);
%     check=['T25_75_'  num2str(CHN(chosenstimchn))];
%     electrodeampint25ns=truedatastruct.(check)(:,1);
%     check=['T50_50_'  num2str(CHN(chosenstimchn))];
%     electrodeampint50ns=truedatastruct.(check)(:,1);
%     if laminar==1
%         check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
%     else
%         check=['T100_' num2str(NORECORDELECT(1))];
%     end
%     electrodeampint0ns=truedatastruct.(check)(:,1);
% ns=[electrodeampint100ns electrodeampint75ns electrodeampint25ns electrodeampint50ns electrodeampint0ns];
% nsE=mean(ns,2);
%     check=['T100_' num2str(CHN(chosenstimchn))];
%     truedatastruct.(check)(:,1)=nsE;
%     check=['T75_25_' num2str(CHN(chosenstimchn))];
%     truedatastruct.(check)(:,1)=nsE;
%     check=['T25_75_'  num2str(CHN(chosenstimchn))];
%     truedatastruct.(check)(:,1)=nsE;
%     check=['T50_50_'  num2str(CHN(chosenstimchn))];
%     truedatastruct.(check)(:,1)=nsE;
%     if laminar==1
%         check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
%     else
%         check=['T100_' num2str(NORECORDELECT(1))];
%     end
%     truedatastruct.(check)(:,1)=nsE;
%%

save('truedatastruct.mat','truedatastruct')
end

