function electrodepreference = Depth_heatmap_and_linecut(AMPInterestSingleLinePlot,trialinfo,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,printFigures,SavetoPPT)
%%used to plot heatmaps outputs electrode preference deep or shallow as
%%well as the electrode of preference in 50/50 trial
%SavetoPPT=1;
%AMPInterestSingleLinePlot=8;%input in uA
AMP=loadAMP;
[PlatChnArray, AmpPlatVal]=AllElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,printFigures);
close all
for createindex=1:length(AMP)
    PlatChnArray(PlatChnArray==AMP(createindex))=createindex;
end
[x, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot));
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1;
[x, AMPInterestSingleLinePlotINDEXSINGLE]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot/2));%x is a throw away variable but it will tell you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXSINGLE=AMPInterestSingleLinePlotINDEXSINGLE+1;
avgnostim=mean(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1),2);%average reponse without stimulation
if SavetoPPT==1
    import mlreportgen.ppt.* %need this to import ppt save format
    TemplateFile = 'C:\Users\smei0006\Documents\myRasterTemplate.pptx'; %template where you can alter slide master and selection pane names layout etc.
    presentationPath = 'depthSpike.pptx'; %saving file
    presentationObj = Presentation(presentationPath,TemplateFile);%create presentation with the specified template
end
loadVarAmp;
loadNORECORDELECT;
 trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;%trials to jump over
 endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions 
 if isempty(endtrialelect)
     endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
 end
lgnd_names{1,(length(NORECORDELECT)*2+1)}=[];
checkconsecutive=1;
checkifsecondloop=0;
numrecordelectinterest=0;
checkNE=0;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
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

if VarAmp==1
    singleLineplotval=zeros(nChn,5);
    singleLineplotvalstd=zeros(nChn,5);
    ampratio=find(AMP<=AMP(end)/2);
    indexcomparabletrials=zeros(length(ampratio),1);
    failtrials=zeros(length(ampratio),1);
    electrodepreference=zeros((size(trialinfo,1)/2)/endtrialelect,4);
else
    singleLineplotval=zeros(nChn,3);
    singleLineplotvalstd=zeros(nChn,3);
    electrodepreference=zeros((size(trialinfo,1)/2)/(endtrialelect/length(NORECORDELECT)),4);
end
estimatelineplotval=zeros(nChn,3);
Largeloopcounter=0;

 for i=1:endtrialelect*2:maxid*2 %i is used to go through groups of related trials
     if SavetoPPT==1
         pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
         replace(pictureSlide, 'Title', 'Depth Heatmaps'); %replace title
     end
     innerLOOPcounter=0;
     Largeloopcounter=Largeloopcounter+1;
     flag=1;
     for j=1:2:trialjump*2 %goes through trials related by trial jump
         desiredchanneltrial_one=(find((cell2mat(trialinfo(:,2))==cell2mat(trialinfo(j+(i-1),2))))+1)/2; %finds trials with desired initial electrode
         desiredchanneltrial_two=find(cell2mat(trialinfo(:,2))==cell2mat(trialinfo(j+1+(i-1),2)))/2; %finds trials with desired second electrode
         chosen_trials=intersect(desiredchanneltrial_one,desiredchanneltrial_two); % finds the trials that intersect with two trials of interest
         if (VarAmp==1)&&(~isempty(chosen_trials(diff(chosen_trials)==1))) %deals with varying amplitude having the same two electrodes for all middle trials
             chosen_trials=chosen_trials(checkconsecutive:3:length(chosen_trials));
             checkconsecutive=checkconsecutive+1; %flag for middle trials
         else
             checkconsecutive=1;
         end
         chosen_trials(chosen_trials>endtrialelect*(i))=[]; % removes any trials that are greater than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         chosen_trials(chosen_trials<i/2)=[]; % removes any trials that are smaller than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         normalisedAvgspikingT=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,chosen_trials)-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
         stdsp=stderrspktrial(:,chosen_trials); %finds the standard deviation of chosen trials for plotting
         if printFigures==1
             fignum=j+(i-1);
             DepthChangeingSpiking_SM(normalisedAvgspikingT, chosen_trials,fignum); %plots heat maps
             if VarAmp==1 && cell2mat(trialinfo(j+1+(i-1),2))==0
                 xlim([0 AMP(end)/2]) % makes axes comparable on single stim trials
             end
             if SavetoPPT==1 %saves to powerpoint
                 if (length(NORECORDELECT)>1)&&(j/length(NORECORDELECT)>ceil(trialjump/length(NORECORDELECT)*(flag)))
                     pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
                     replace(pictureSlide,'Title', 'Depth Heatmaps'); %replace title
                     flag=flag+1;
                     innerLOOPcounter=0;
                 end
                 innerLOOPcounter=innerLOOPcounter+1;
                 num=num2str(innerLOOPcounter);
                 saveas(gcf,['Depth_plot' num2str(j+(i-1)) '.png'])
                 pichandle=Picture(['Depth_plot' num2str(j+(i-1)) '.png']);%save figure as picture
                 replace(pictureSlide,['Picture ' num],pichandle); %replace picture with name Picture X
             end
         end
         singleLineplotval(:,(j+1)/2)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUAL); %finds amplitude of interest trial i.e. 3rd AMP level for dual trials
         singleLineplotvalstd(:,(j+1)/2)=(1000/(secondstoanalyse-startpointseconds))*stdsp(:,AMPInterestSingleLinePlotINDEXDUAL); %standard devation for line plot converted to Sp/s
         
         if (cell2mat(trialinfo(j+1+(i-1),2)))==0
             singleLineplotval(:,(j+1)/2)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXSINGLE);% matches second AMP level of single trials
             singleLineplotvalstd(:,(j+1)/2)=(1000/(secondstoanalyse-startpointseconds))*stdsp(:,AMPInterestSingleLinePlotINDEXSINGLE); %single amp trial standard deviation converted to Sp/s
             checkifsecondloop=checkifsecondloop+1;%used to determine if values from two single electrode trials are available
             if checkifsecondloop>=2
                 if SavetoPPT==1
                     pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
                     replace(pictureSlide, 'Title', 'Depth Heatmaps'); %replace title
                 end
                 
                 %if length(NORECORDELECT)>1&&checkNE==1 %used for if NORECORDELECT>1 since there is only one set of single electrode trials for both second electrode trials
                    % checkNE=0;
                     %checkifsecondloop=0;
                 if length(NORECORDELECT)>1 && checkNE==0
                     if numrecordelectinterest==length(NORECORDELECT)-1
                        %checkNE=1;
                        checkifsecondloop=0;
                     else 
                         checkifsecondloop=1;
                     end
                     numrecordelectinterest=numrecordelectinterest+1;
                 else
                     checkifsecondloop=0;
                 end
                 
                 if VarAmp==1
                     for ratioplot=0.25:0.25:0.75 %ratio values to iterate through
                         additivespikingsingchn=normalisedAvgspikingTsingchn1*2*(1-ratioplot)+normalisedAvgspikingT*2*(ratioplot);%adds together ratio of spike rates for prediction based on single electrode trials
                         
                         %                      for chncount=1:nChn
                         %                          for replaceval=1:size(additivespikingsingchn,2)
                         %                              if additivespikingsingchn(chncount,replaceval)>AmpPlatVal(chncount)+0.1*AmpPlatVal(chncount)
                         %                                  additivespikingsingchn(chncount,replaceval)=AmpPlatVal(chncount)+0.1*AmpPlatVal(chncount);
                         %                              end
                         %                          end
                         %                      end
                         if printFigures==1
                             fignum=ratioplot/0.25+(i-1)+500;
                             DepthChangeingSpiking_SM(additivespikingsingchn, chosen_trials,fignum);%plots heatmap of estimated plots
                             if VarAmp==1 && cell2mat(trialinfo(j+1+(i-1),2))==0
                                 xlim([0 AMP(end)/2])%compensates single electrode axes
                                 chnT1=j+(i-1)-1;
                                 xlabel_halvamp=[0:2:AMP(end)];
                                 xticklabels({xlabel_halvamp});
                             end
                             title(['Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo(j+(i-1),2))) ' ' num2str(cell2mat(trialinfo(chnT1,2))) ' ESTIMATED '  num2str((ratioplot)*100) '/' num2str((1-ratioplot)*100)])
                             yline(cell2mat(trialinfo(chnT1,2)),'Color','r','Linewidth',3*(1-ratioplot),'Alpha',1)
                             yline(cell2mat(trialinfo(j+(i-1),2)),'Color','r','Linewidth',3*(ratioplot),'Alpha',1)
                             if SavetoPPT==1
                                 saveas(gcf,['Depth_plot_estimate' num2str(i*100+j*5+3*ratioplot/0.25) '.png'])
                                 pichandle=Picture(['Depth_plot_estimate' num2str(i*100+j*5+3*ratioplot/0.25) '.png']);%save figure as picture
                                 replace(pictureSlide,['Picture ' num2str((ratioplot/0.25)+1)],pichandle); %replace picture with name Picture X
                             end
                         end
                         estimatelineplotval(:,(ratioplot/0.25))=additivespikingsingchn(:,AMPInterestSingleLinePlotINDEXDUAL);
                     end
                     if SavetoPPT==1
                         pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
                         replace(pictureSlide, 'Title', 'Depth Heatmaps ERROR'); %replace title
                     end
                     for ratioplot=0.25:0.25:0.75 %ratio values to iterate through
                         %PLOTTING OF ERROR GRAPH
                         indexcomparabletrials(1)=1;
                         for l=2:length(ampratio)
                             try
                                indexcomparabletrials(l)=find(AMP==AMP(ampratio(l))*2);%%%%%%%%fix this
                             catch
                                 indexcomparabletrials(l)=0;
                                 failtrials(l)=1;
                             end
                         end
                         indexcomparabletrials(indexcomparabletrials==0)=[];
                         if printFigures==1
                             fignum=(ratioplot*2/0.25)+1+(i-1);
                             figure(fignum)
                             object_handles = findall(gcf,'Type','Surface');
                             ZdataResultplot=object_handles.ZData;
                             fignum=ratioplot/0.25+(i-1)+500;
                             figure(fignum)
                             object_handles = findall(gcf,'Type','Surface');
                             ZdataEstimateplot=object_handles.ZData;
                             fignum=ratioplot/0.25+(i-1)+1000;
                             DepthChangeingSpiking_SM(ZdataResultplot(:,indexcomparabletrials)-ZdataEstimateplot(:,ampratio(failtrials==0)), chosen_trials(ampratio(failtrials==0)),fignum);%plots heatmap of estimated plots
                             title(['ERROR in estimation. Stimchn: ' num2str(cell2mat(trialinfo(j+(i-1),2))) ' ' num2str(cell2mat(trialinfo(chnT1,2))) ' @ '  num2str((ratioplot)*100) '/' num2str((1-ratioplot)*100)])
                             yline(cell2mat(trialinfo(chnT1,2)),'Color','r','Linewidth',3*(1-ratioplot),'Alpha',1)
                             yline(cell2mat(trialinfo(j+(i-1),2)),'Color','r','Linewidth',3*(ratioplot),'Alpha',1)
                             x=get(gca,'xticklabels');%just temp variable
                             xlabel_halvamp=[0:AMP(indexcomparabletrials(end))/(length(x)-1):AMP(indexcomparabletrials(end))];
                             xticklabels({xlabel_halvamp});
                             caxis([-150 150])
                             if SavetoPPT==1
                                 saveas(gcf,['Depth_plot_ERROR' num2str(i*100+j*5+3*ratioplot/0.25) '.png'])
                                 pichandle=Picture(['Depth_plot_ERROR' num2str(i*100+j*5+3*ratioplot/0.25) '.png']);%save figure as picture
                                 replace(pictureSlide,['Picture ' num2str((ratioplot/0.25)+1)],pichandle); %replace picture with name Picture X
                             end
                         end
                     end
                     errorlineplot=singleLineplotval(:,2:4)-estimatelineplotval;
                     if printFigures==1
                         if SavetoPPT==1
                             pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
                             replace(pictureSlide, 'Title', 'Depth Heatmaps'); %replace title
                         end
                         figure %plotting estimation line plot
                         ax = gca;
                         ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
                         hold on
                         SMOOTHING=1;
                         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
                         colorOrder = get(gca, 'ColorOrder');
                         for p=1:size(estimatelineplotval,2)
                             rate = conv(estimatelineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                             rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                             plot(rate)
                         end
                         ylabel('Sp/s')
                         xlabel('Channel number')
                         xlim([1 32])
                         xline(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)),'-.k')
                         xline(cell2mat(trialinfo((chosen_trials(3)-2)*2,2)),'k')
                         legend('75/25','50/50','25/75','Stim E1','Stim E2')
                         title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2,2))) ' & ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2))) ' ESTIMATION'])
                         if SavetoPPT==1
                             saveas(gcf,['LineplotEstimation' num2str(i*100+j*5+3*ratioplot/0.25) '.png'])
                             pichandle=Picture(['LineplotEstimation' num2str(i*100+j*5+3*ratioplot/0.25) '.png']);%save figure as picture
                             replace(pictureSlide,['Picture ' num2str(3)],pichandle); %replace picture with name Picture X
                         end

                         figure %plotting error line plot
                         ax = gca;
                         ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
                         hold on
                         SMOOTHING=1;
                         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
                         colorOrder = get(gca, 'ColorOrder');
                         for p=1:size(errorlineplot,2)
                             rate = conv(errorlineplot(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                             rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                             plot(rate)
                         end
                         ylabel('Sp/s')
                         xlabel('Channel number')
                         xlim([1 nChn])
                         xline(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)),'-.k')
                         xline(cell2mat(trialinfo((chosen_trials(3)-2)*2,2)),'k')
                         legend('75/25','50/50','25/75','Stim E1','Stim E2')
                         title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2,2))) ' & ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2))) ' ERROR'])
                         if SavetoPPT==1
                             saveas(gcf,['LineplotERROR' num2str(i*100+j*5+3*ratioplot/0.25) '.png'])
                             pichandle=Picture(['LineplotERROR' num2str(i*100+j*5+3*ratioplot/0.25) '.png']);%save figure as picture
                             replace(pictureSlide,['Picture ' num2str(4)],pichandle); %replace picture with name Picture X
                         end
                     end
                 else
                     additivespikingsingchn=normalisedAvgspikingTsingchn1+normalisedAvgspikingT;%adds together ratio of spike rates for prediction based on single electrode trials
                     estimatelineplotval=additivespikingsingchn(:,AMPInterestSingleLinePlotINDEXDUAL);
                     errorlineplot=singleLineplotval(:,numrecordelectinterest)-estimatelineplotval;
                     if printFigures==1
                         figure
                         fignum=gcf;
                         fignum=fignum.Number;
                         DepthChangeingSpiking_SM(additivespikingsingchn, chosen_trials,fignum);%plots heatmap of estimated plots
                         title(['Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo(j+(i-1),2))) ' ' num2str(cell2mat(trialinfo(chnT1,2))) ' ESTIMATED '  '50/50'])
                         yline(cell2mat(trialinfo(chnT1,2)),'Color','r','Linewidth',1.5,'Alpha',1)
                         yline(cell2mat(trialinfo(j+(i-1),2)),'Color','r','Linewidth',1.5,'Alpha',1)
                         if SavetoPPT==1
                             saveas(gcf,['Depth_plot_estimate' num2str(i*100+j*5+3) '.png'])
                             pichandle=Picture(['Depth_plot_estimate' num2str(i*100+j*5+3) '.png']);%save figure as picture
                             replace(pictureSlide,['Picture ' num2str(2)],pichandle); %replace picture with name Picture X
                         end
                         ZdataEstimateplot=additivespikingsingchn;
                         %PLOTTING OF ERROR GRAPH
                         if length(NORECORDELECT)>1
                             fignum=j+(i-1)-4-length(NORECORDELECT);
                         else
                             fignum=j+(i-1)-4;
                         end
                         figure(fignum)
                         object_handles = findall(gcf,'Type','Surface');
                         ZdataResultplot=object_handles.ZData;
                         figure
                         fignum=gcf;
                         fignum=fignum.Number;
                         DepthChangeingSpiking_SM(ZdataResultplot-ZdataEstimateplot, chosen_trials,fignum);%plots heatmap of estimated plots
                         title(['ERROR in estimation. Stimchn: ' num2str(cell2mat(trialinfo(j+(i-1),2))) ' ' num2str(cell2mat(trialinfo(chnT1,2))) ' @ '  '50/50'])
                         yline(cell2mat(trialinfo(chnT1,2)),'Color','r','Linewidth',1.5,'Alpha',1)
                         yline(cell2mat(trialinfo(j+(i-1),2)),'Color','r','Linewidth',1.5,'Alpha',1)
                         caxis([-150 150])
                         if SavetoPPT==1
                             saveas(gcf,['Depth_plot_ERROR' num2str(i*100+j*5+3) '.png'])
                             pichandle=Picture(['Depth_plot_ERROR' num2str(i*100+j*5+3) '.png']);%save figure as picture
                             replace(pictureSlide,['Picture ' num2str(3)],pichandle); %replace picture with name Picture X
                         end
                         if SavetoPPT==1
                             pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
                             replace(pictureSlide, 'Title', 'Depth Heatmaps'); %replace title
                         end
                         figure %plotting estimation line plot
                         hold on
                         SMOOTHING=1;
                         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
                         for p=1:size(estimatelineplotval,2)
                             rate = conv(estimatelineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                             rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                             plot(rate)
                         end
                         ylabel('Sp/s')
                         xlabel('Channel number')
                         xlim([1 nChn])
%                          xline(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)),'-.k')
%                          xline(cell2mat(trialinfo((chosen_trials(3)-2)*2,2)),'k')
                         dualtrial=floor( chosen_trials(1) / trialjump ) * trialjump +1;
                         xline(cell2mat(trialinfo(dualtrial*2-1,2)),'-.k')
                         xline(cell2mat(trialinfo(dualtrial*2,2)),'k')
                         legend('50/50', 'Stim E1','Stim E2')
                         title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(1))*2-1,2))) ' & ' num2str(cell2mat(trialinfo((chnT1),2))) ' ESTIMATION'])
                         if SavetoPPT==1
                             saveas(gcf,['LineplotEstimation' num2str(i*100+j*5+3) '.png'])
                             pichandle=Picture(['LineplotEstimation' num2str(i*100+j*5+3) '.png']);%save figure as picture
                             replace(pictureSlide,['Picture ' num2str(3)],pichandle); %replace picture with name Picture X
                         end
                         
                         
                         figure %plotting error line plot
                         
                         hold on
                         SMOOTHING=1;
                         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
                         for p=1:size(errorlineplot,2)
                             rate = conv(errorlineplot(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                             rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                             plot(rate)
                         end
                         ylabel('Sp/s')
                         xlabel('Channel number')
                         xlim([1 nChn])
                         xline(cell2mat(trialinfo((chnT1),2)),'-.k')
                         xline(cell2mat(trialinfo((chosen_trials(1))*2-1,2)),'k')
                         legend('50/50', 'Stim E1','Stim E2')
                         title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(1))*2-1,2))) ' & ' num2str(cell2mat(trialinfo((chnT1),2)))  ' ERROR'])
                         if SavetoPPT==1
                             saveas(gcf,['LineplotERROR' num2str(i*100+j*5+3) '.png'])
                             pichandle=Picture(['LineplotERROR' num2str(i*100+j*5+3) '.png']);%save figure as picture
                             replace(pictureSlide,['Picture ' num2str(4)],pichandle); %replace picture with name Picture X
                         end
                     end
                 end
             else
                 numrecordelectinterest=0;
                 normalisedAvgspikingTsingchn1=normalisedAvgspikingT;%if this is the first single electrode data, then keep it for the next loop with the second electrode
                 chnT1=j+(i-1);%note which channel number this is so we can use it later
             end
         end
     end
     if printFigures==1
         figure %plotting real data line plot
         hold on
         SMOOTHING=1;
         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
         for p=1:size(singleLineplotval,2)
             rate = conv(singleLineplotval(:,p),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
             rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
             errorbar(rate,singleLineplotvalstd(:,p))
         end
         if VarAmp==1
             legend('100/0','75/25','50/50','25/75','0/100', 'Stim E1','Stim E2')
             xline(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)),'-.k')
             xline(cell2mat(trialinfo((chosen_trials(3)-2)*2,2)),'k')
             title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2,2))) ' & ' num2str(cell2mat(trialinfo((chosen_trials(3)-2)*2-1,2)))])
         elseif (length(NORECORDELECT)==1) && (VarAmp==0)
             legend('50/50','0/100','100/0','Stim E1','Stim E2')
             ratioplot=0.25;
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
             ratioplot=0.25;
         end
         ylabel('Sp/s')
         xlabel('Channel number')
         xlim([1 nChn])         
         if SavetoPPT==1
             saveas(gcf,['LineplotReal' num2str(i*100+j*5+3*ratioplot/0.25) '.png'])
             pichandle=Picture(['LineplotReal' num2str(i*100+j*5+3*ratioplot/0.25) '.png']);%save figure as picture
             replace(pictureSlide,['Picture ' num2str(2)],pichandle); %replace picture with name Picture X
         end
     end
     %find where the average peak occurs 
     if VarAmp==1
        averagespElect=find(mean(singleLineplotval(:,2:4),2)==max(mean(singleLineplotval(:,2:4),2)));
     else
         averagespElect=find(singleLineplotval(:,1)==max(singleLineplotval(:,1)));
     end
     if length(NORECORDELECT)>1 && VarAmp==0
         for norecrodElect=1:length(NORECORDELECT)
             averagespElect=find(singleLineplotval(:,norecrodElect)==max(singleLineplotval(:,norecrodElect)));
         if abs(cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(norecrodElect-1))*2-1,2))-averagespElect(1))<abs(cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(norecrodElect-1))*2,2))-averagespElect(1))
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,1)='D';%for deep, note this converts to ascii
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,2)=averagespElect(1);%for deep
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,3)=cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(norecrodElect-1))*2-1,2));
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,4)=cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(norecrodElect-1))*2,2));
         else
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,1)='S';%for shallow
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,2)=averagespElect(1);%for deep
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,3)=cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(norecrodElect-1))*2-1,2));
             electrodepreference((Largeloopcounter-1)*length(NORECORDELECT)+norecrodElect,4)=cell2mat(trialinfo((chosen_trials(1)-trialjump+1+(norecrodElect-1))*2,2));
         end
         end
     else
         if abs(cell2mat(trialinfo(chosen_trials(1)*2-1,2))-averagespElect)<abs(cell2mat(trialinfo(chosen_trials(1)*2-2,2))-averagespElect)
             electrodepreference(Largeloopcounter,1)='D';%for deep, note this converts to ascii
             electrodepreference(Largeloopcounter,2)=averagespElect;%for deep
             electrodepreference(Largeloopcounter,3)=cell2mat(trialinfo(chosen_trials(1)*2-1,2));
             electrodepreference(Largeloopcounter,4)=cell2mat(trialinfo(chosen_trials(1)*2-2,2));
         else
             electrodepreference(Largeloopcounter,1)='S';%for shallow
             electrodepreference(Largeloopcounter,2)=averagespElect;%for deep
             electrodepreference(Largeloopcounter,3)=cell2mat(trialinfo(chosen_trials(1)*2-1,2));
             electrodepreference(Largeloopcounter,4)=cell2mat(trialinfo(chosen_trials(1)*2-2,2));
         end
     end
 end
 if SavetoPPT==1
     replace(pictureSlide,'Title', 'Depth Heatmaps'); %replace title
     close(presentationObj); %close presentation to keep changes
     if ispc
         winopen('depthSpike.pptx'); %open presentation (WINDOWS ONLY FUNCTION)
     end
 end
 
end

