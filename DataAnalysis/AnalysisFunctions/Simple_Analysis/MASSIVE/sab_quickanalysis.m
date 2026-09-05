%% For quick analysis during experiment

% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-500; %removes spikes below this threshold
artefact_high=500; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=8; %How long after the trigger do you want to analyse spikes for(ms)? 
printspiking=0;
par=0;

%% 1. Blank stimulus

FS=30000;
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
dName='amplifier';
vFID = fopen([filepath filesep dName '.dat'],'r');
mem_check=dir('amplifier.dat');
T = mem_check.bytes ./ (2 * nChn * 30000);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * 30000);
if t_len < (Overall_time_to_analyse+Startpoint_analyse)
    error('Time to analyse exceeds the data recording period')
end
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
elseif Overall_time_to_analyse~=0
    T=Overall_time_to_analyse;
end
info = fileinfo.bytes/2;
nL = (ceil(info / (nChn*FS*double(T)))+1);
vblank=[];
BREAK = 1;
N=1;
denoiseIntan_sab(filepath, dName, T, par, Startpoint_analyse, Overall_time_to_analyse);
trig = loadTrig(0);
theseTrig = trig;

%% 2. Thresholds & Mu
allExtract_sab_1(dName,filepath,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');

%% Extract trial information
%trialinfo=TrialInfo_sab;

%% 4. Calculate Structure of sorted trials according to IDs
SavetoPPT=0;
starttrial=1;
trialjump=1;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
[IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
save('IDstruct.mat','IDstruct','baslinespikestruct')
%%
if printspiking==1 && SavetoPPT==1
    import mlreportgen.ppt.* %need this to import ppt save format
    TemplateFile = 'C:\Users\smei0006\Documents\myRasterTemplate.pptx'; %template where you can alter slide master and selection pane names layout etc.
    presentationPath = 'Spiking.pptx'; %saving file
    presentationObj = Presentation(presentationPath,TemplateFile);%create presentation with the specified template
    for j=1:floor(nChn/10)
        pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
        Loopcounter=0;
        for i=2:2:10
             if ishandle(i+(j-1)*10)
                 Loopcounter=Loopcounter+1;
                 figure(i+(j-1)*10)
                 fig=gcf;
                 saveas(gcf,['spike_plot' num2str(i+(j-1)*10) '.png'])
                 pichandle=Picture(['spike_plot' num2str(i+(j-1)*10) '.png']);%save figure as picture
                 replace(pictureSlide,['Picture ' num2str(Loopcounter)],pichandle); %replace picture with name Picture X
             end
        end
        replace(pictureSlide,'Title', 'Spiking'); %replace title
    end
    close(presentationObj); %close presentation to keep changes
    if ispc
        winopen('Spiking.pptx'); %open presentation (WINDOWS ONLY FUNCTION)
    end
end
%% 5. Calculates template of trials and spiking responses (Output in true electrode order)
[avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct, baslinespikestruct);
 save('Averagetrialresponse.mat','avgnospT','stderrspktrial')
 %% 8. plotting average electrode response for all electrodes classes as significant and not significant
 AMP=loadAMP;
 loadStimChn;
 avgnospT_sps=(1000/(secondstoanalyse-startpointseconds)).*avgnospT;
 cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
 avgnostim=mean(avgnospT_sps(:,cell2mat(trialinfo(1:2:end,18))==-1),2);
 numofsigperchn=zeros(nChn,1);
whichchn=zeros(nChn,1);
 %testifsignificant=[(avgnospT_sps(:,cell2mat(trialinfo(16:2:28,18))==AMP(end))) (avgnospT_sps(:,cell2mat(trialinfo(14:2:28,18))==AMP(end-1))) (avgnospT_sps(:,cell2mat(trialinfo(14:2:28,18))==AMP(end-2)))];
 %testsig_maxamp=testifsignificant(:,1:size(avgnostim,2));
 chn_sigboth=zeros(2,nChn);
 for Chosenstimchn=1:length(stimChn)
    desiredchanneltrial=find(cell2mat(trialinfo(:,2))==stimChn(Chosenstimchn)); %finds trials with desired initial electrode
    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds trials with matching recording electrode spacing
    Desired_trialfirst=cell2mat(trialinfo(desiredchannel__singleampmath,1));%array of mathcning trial number
     testifsignificant=(avgnospT_sps(:,Desired_trialfirst(end-2:end)));
     chn_sig=zeros(1,nChn);
     for i=1:nChn
         [h,p]=ttest(testifsignificant(i,:)',avgnostim(i)','tail','right','Alpha',0.05);%one-tailed ttest determineing if the electrode has a mean significantly larger than the no stimulation trials with a 95% confidence level
         chn_sig(i)=h;
     end
     chn_sig(isnan(chn_sig))=0;
     numofsigperchn(stimChn(Chosenstimchn))=sum(chn_sig);
     chn_sigboth(Chosenstimchn,:)=chn_sig;
     whichchn(stimChn(Chosenstimchn))=1;
 end
 save('sigchn.mat','numofsigperchn','whichchn','chn_sigboth')
 

 %AllElectrodeResponseCurve_SM(trialinfo,avgnospT_sps,stderrspktrial,1,h);
% AllElectrodeResponseCurve_SM(trialinfo,avgnospT_sps,(1000/(secondstoanalyse-startpointseconds)).*stderrspktrial,1,h);
%% Plots Heatmaps
depthdriven=1000-50;%in um -50 frm tip
SavetoPPT=0; %%if single trial, you cannot automatically save to ppt.
AMPInterestSingleLinePlot=4;%input in uA
ERROR_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,depthdriven)

TrueData_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,depthdriven);
AdditivePrediction_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT, stderrspktrial,startpointseconds, secondstoanalyse, depthdriven)

 figure 
 plot(1:nChn,chn_sigboth)
  

printFigures=1;
%note electrode preference depends on amplitude
electrodepreference = Depth_heatmap_and_linecut(AMPInterestSingleLinePlot,trialinfo,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse, printFigures,SavetoPPT); 




%% Raster or histogram of firing rate
 chn=16;
generate_StackedRaster_sab(chn);
 %% Plots unfiltered data
tID=[7]; %trial ID of interest
Chan=[15];%channels of interest - recommend only one or two for raw
Chosen_trig=15;
Binstart=-250;
Binend=400;
stimConcatenate(tID,Chosen_trig,Chan,'RAW',Binstart,Binend);
title(['Chn ' ,num2str(Chan), '. Stimchn ' ,num2str(cell2mat(trialinfo(tID*2-1,2))),' & ' num2str(cell2mat(trialinfo(tID*2,2))), ' at ' ,num2str(cell2mat(trialinfo(tID*2-1,18))), 'uA filtered'])
xlabel('Time (ms)')

%% Plots filtered data for selected channels and trial ID
Binstart=-250;
Binend=400;
tID=[7]; %trial ID of interest
Chan=[15];%channels of interest 
Chosen_trig=15;
stimConcatenate(tID,Chosen_trig,Chan,'MU',Binstart,Binend);
title(['Chn ' ,num2str(Chan), '. Stimchn ' ,num2str(cell2mat(trialinfo(tID*2-1,2))),' & ' num2str(cell2mat(trialinfo(tID*2,2))), ' at ' ,num2str(cell2mat(trialinfo(tID*2-1,18))), 'uA filtered'])
xlabel('Time (ms)')

%% 7. Plot stimulating electrode response curve
StimElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial);


 %% hypothesis
loadStimChn;
loadNORECORDELECT;
outputMaxHist=[];
 for i=1:length(CHN)
     for j=1:NORECORDELECT(end)+2
         chn=CHN(i)+j-1;
         outputMaxHist=hypothplot_sab(chn,outputMaxHist);
     end
 end
 save('outputMaxHist.mat','outputMaxHist');
 %%
 loadVarAmp;
 AMP=loadAMP;
 AMP(1)=0;
 chosenfield='Chn_5_StimChn_8_3';
 ddddaata=outputMaxHist.(chosenfield);
 ddddaata(:,1)=mean(ddddaata(:,1));
 figure
 hold on
 if VarAmp==1
     for i=1:5
         plot(AMP(1:4),ddddaata(i,1:4)')
     end
     legend('100/0','75/25','50/50','25/75','0/100')
 else
     for i=1:3
         plot(AMP(1:4),ddddaata(i,1:4)')
     end
     legend('100/0','50/50','0/100')
 end

 title(chosenfield,'Interpreter','none')
 xlabel('Current uA')
 ylabel('Sp/s')
 ylim([0 450])

 %% 8. plotting average electrode response for all electrodes classes as significant and not significant
 chnlist=zeros(nChn,1);
 trialinfo=loadTrialInfo;
 trialinfo(1,:)=[];
 for channel=1:nChn
 chn=Depth(1);
 chn=chn(channel);
 loadNREP;
  cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
  counter=0;
  currents=zeros(n_REP_true,n_REP_true/cond);
 for trial=1:cond:n_REP_true
     counter=counter+1;
    currents(:,counter)=IDstruct.(['T' num2str(trial)])(chn,:)';
 end
 [p,h,stats] = ranksum(sum(currents(:,1:4),2),sum(currents(:,5:8),2));
 %p = kruskalwallis(currents);
 %[p,tbl,stats] = anova1(currents);
 if p<0.05
     chnlist(channel)=1;
 end

 end
  chnlist=zeros(nChn,1);

 AMP_orig=loadAMP;
 AMP_orig=AMP_orig./2;
  avgnostim=(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1)).*1000/(secondstoanalyse-startpointseconds); 
  testifsignificantchn1=[(avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end)) + (cell2mat(trialinfo(2:2:end,18))==AMP_orig(end))+ (cell2mat(trialinfo(1:2:end,2))==stimChn(1)))==3)) ...
      (avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end-1)) + (cell2mat(trialinfo(2:2:end,18))==AMP_orig(end-1))+ (cell2mat(trialinfo(1:2:end,2))==stimChn(1)))==3))...
      (avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end-2)) + (cell2mat(trialinfo(2:2:end,18))==AMP_orig(end-2))+ (cell2mat(trialinfo(1:2:end,2))==stimChn(1)))==3))];
 testifsignificantchn2=[(avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end)) + (cell2mat(trialinfo(2:2:end,18))==AMP_orig(end))+ (cell2mat(trialinfo(1:2:end,2))==stimChn(2)))==3)) ...
      (avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end-1)) + (cell2mat(trialinfo(2:2:end,18))==AMP_orig(end-1))+ (cell2mat(trialinfo(1:2:end,2))==stimChn(2)))==3))...
      (avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end-2)) + (cell2mat(trialinfo(2:2:end,18))==AMP_orig(end-2))+ (cell2mat(trialinfo(1:2:end,2))==stimChn(2)))==3))];
% 
%    testifsignificantchnb=[(avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end)) + (cell2mat(trialinfo(1:2:end,2))==stimChn(1))+ (cell2mat(trialinfo(2:2:end,2))==stimChn(2)))==3)) ...
%       (avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end-1)) + (cell2mat(trialinfo(1:2:end,2))==stimChn(1))+ (cell2mat(trialinfo(2:2:end,2))==stimChn(2)))==3))...
%       (avgnospT(:,((cell2mat(trialinfo(1:2:end,18))==AMP_orig(end-2)) + (cell2mat(trialinfo(1:2:end,2))==stimChn(1))+ (cell2mat(trialinfo(2:2:end,2))==stimChn(2)))==3))];

  testsig_maxamp=testifsignificant(:,1:size(avgnostim,2)).*1000/(secondstoanalyse-startpointseconds);
 testsigch1=testifsignificantchn1.*1000/(secondstoanalyse-startpointseconds);
  testsigch2=testifsignificantchn2.*1000/(secondstoanalyse-startpointseconds);
  testsigchb=testifsignificantchnb.*1000/(secondstoanalyse-startpointseconds);
 for channel=1:nChn
     [p1,h,stats] = ranksum(testsigch1(channel,:),avgnostim(channel,:));
     [p2,h,stats] = ranksum(testsigch2(channel,:),avgnostim(channel,:));
     if (p1<0.05) && (p2<0.05)
         chnlist(channel)=1;
     end
 end
 [h,p]=ttest(testsig_maxamp',avgnostim','tail','right','Alpha',0.05);%one-tailed ttest determineing if the electrode has a mean significantly larger than the no stimulation trials with a 95% confidence level
 AllElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,1,h);

%% 6. Print all trial details
Legend=cell(nChn,1);
loadStimChn;
loadStimParams;
stimamplitudes=unique(cell2mat(StimParams(2:end,16)));
fprintf('The current amplitude in ascending order is: \n');
fprintf('%d\n',stimamplitudes);
fprintf('The stimulation channels were: \n');
fprintf('%d\n',stimChn);

%% 8. plotting average electrode response for all electrodes classes as significant and not significant
AllElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,1,h);

%% 9. plotting electrodes between stim and stim electrodes
InbetweenElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial);

%% Extracts responses from chosen single electrode stimulation for specific hannel
Chan=[14 15 16 17];%can input as many channels of interest as you want
stimChn=15;%only one stim chan
singleElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,Chan,stimChn);

%% Curve according to chosen trial with number of trials to skip
Chan=[15 16];%can input as many as you want
starttrial=3;%starting trial to analyse
trialjump=5;%number of trials to jump for next consecutive response
singleElectrodeResponseCurveAMPLITUDE_SM(trialinfo,avgnospT,stderrspktrial,Chan,starttrial,trialjump);

%% Plot response curve change over time (3D graph)
 TimeBegin=0; %How long after the trigger do you want skip spike analysis(ms)? 
 TimeEnd=14; %How long after the trigger do you want to stop spike analysis(ms)? 
 TimeStep=2; %Timestep of spike analysis between TimeBegin and TimeEnd(ms)? 
 starttrial=3;%starting trial to analyse
 endtrial=40;%trial to stop analysing
 electrode=14;
 trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %number of trials to jump for next consecutive response equal to cond as default
 TimeChangeinSpiking_SM(TimeBegin, TimeEnd, TimeStep);%, electrode, starttrial, trialjump, endtrial); %can input end trial but haven't tested
 
 %%

%% PCA spike sorting attempt
% CHNinterest=25;
% plotPCAspiking_sab(CHNinterest);

%%     %for i=1:length(chosen_trials)%trialjump
      %   starttrial=chosen_trials(1);%i+(j-1)*endtrialelect;
      %   endtrial=chosen_trials(end);%j*endtrialelect;

%% Deleting non-original files from filepath
%Delete_non_orig_sab;



%% MUA response

MUA_response_SM;
Binstart=-250;
Binend=400;
tID=[10]; %trial ID of interest
Chan=[20 21 22 23];%channels of interest
Chosen_trig=1;
stimConcatenate(tID,Chosen_trig,Chan,'MUA',Binstart,Binend);
title(['Chn ', num2str(Chan),' at ' string(trialinfo(tID*2,18)) 'uA MUA (filtered)'])
xlabel('Time (ms)')

%%
 SavetoPPT=1;
if SavetoPPT==1
    import mlreportgen.ppt.* %need this to import ppt save format
    TemplateFile = 'C:\Users\smei0006\Documents\myRasterTemplate.pptx'; %template where you can alter slide master and selection pane names layout etc.
    presentationPath = 'channelActivationSpread.pptx'; %saving file
    presentationObj = Presentation(presentationPath,TemplateFile);%create presentation with the specified template
end

loadVarAmp;
loadNORECORDELECT;

trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/(2);
endtrialelect=75;
checkconsecutive=1;
num_archive=NORECORDELECT(1);
SMOOTHING=0.6;
Z=[];
avgnostim=mean(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1),2);%avgnostim=mean(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1),2);
for i=1:endtrialelect*2:maxid*2
    innerLOOPcounter=0;
    flag=1;
    for K=1:length(NORECORDELECT)
        if SavetoPPT==1
            pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
            replace(pictureSlide, 'Title', 'Spread of activation along probe'); %replace title
        end
        if length(NORECORDELECT)==1
            temp=trialjump*2;
        else 
            temp=trialjump;
        end
        for j=1:2:temp%(maxid/endtrialelect)
            if (K>1)&&(j>=(trialjump-1))
            desiredchanneltrial_one=(find((cell2mat(trialinfo(:,2))==cell2mat(trialinfo(-1+(i-1)+(K-1)*(trialjump+1),2))))+1)/2; %finds trials with desired initial electrode
            desiredchanneltrial_two=find(cell2mat(trialinfo(:,2))==cell2mat(trialinfo((i-1)+(K-1)*(trialjump+1),2)))/2;
            chosen_trials=intersect(desiredchanneltrial_one,desiredchanneltrial_two);
            k=1;
            else
            desiredchanneltrial_one=(find((cell2mat(trialinfo(:,2))==cell2mat(trialinfo(j+(i-1)+(K-1)*(trialjump+1),2))))+1)/2; %finds trials with desired initial electrode
            desiredchanneltrial_two=find(cell2mat(trialinfo(:,2))==cell2mat(trialinfo(j+1+(i-1)+(K-1)*(trialjump+1),2)))/2;
            chosen_trials=intersect(desiredchanneltrial_one,desiredchanneltrial_two);
            k=0;
            end
            if (VarAmp==1)&&(~isempty(chosen_trials(diff(chosen_trials)==1)))
                chosen_trials=chosen_trials(checkconsecutive:3:length(chosen_trials));
                checkconsecutive=checkconsecutive+1;
            else
                checkconsecutive=1;
            end
            
            if (j==1)||(j==5)||(j==9)
                for l=2:length(chosen_trials)
                    Z(:, l-1)=avgnospT(:,chosen_trials(l))-avgnostim;
                    window = normpdf(-2*SMOOTHING:2*SMOOTHING,0,SMOOTHING);
                    rate = conv(Z(1:nChn,l-1),window);%1000ms over 12-2ms number of trials
                    rate = rate(round(2*SMOOTHING+1):end-round(2*SMOOTHING));

                    rate = conv(flipud(rate),window);%1000ms over 12-2ms number of trials
                    rate=flipud(rate);
                    rate = rate(round(2*SMOOTHING+1):end-round(2*SMOOTHING));
                    figure(l+600+(i-1)*10+(K-1)*100)
                    hold on
                    plot(1:nChn, rate,'LineWidth',1.5);
                    if (j>=(temp-1))
                        legend('100/0','50/50','0/100')
                        ylim([0 160])
                        xlabel('Channel')
                        ylabel('average number sp/trial 2-12m')
                        if K==1
                            title(['Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo(((chosen_trials(2)-3)*2),2))) ' & ' num2str(cell2mat(trialinfo(((chosen_trials(2)-3)*2)-1,2))) ' @ ' num2str(cell2mat(trialinfo((chosen_trials(l)*2)-1,18))) 'uA']);
                        else
                            title(['Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo(((chosen_trials(1)+2)*2),2))) ' & ' num2str(cell2mat(trialinfo(((chosen_trials(1)+3)*2)-1,2))) ' @ ' num2str(cell2mat(trialinfo((chosen_trials(l)*2)-1,18))) 'uA']);
                        end
                        if SavetoPPT==1
%                             if (length(NORECORDELECT)>1)&&(j/length(NORECORDELECT)>ceil(trialjump/length(NORECORDELECT)*(flag)))
%                                 pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
%                                 replace(pictureSlide,'Title', 'Depth Heatmaps'); %replace title
%                                 flag=flag+1;
%                                 OVERALLOOPcounter=0;
%                             end
                            innerLOOPcounter=innerLOOPcounter+1;
                            if innerLOOPcounter>5
                                pictureSlide = add(presentationObj,'PictureLayout2'); %Create picture slide - custom layout
                                replace(pictureSlide, 'Title', 'Spread of activation along probe'); %replace title
                                innerLOOPcounter=1;
                            end
                            num=num2str(innerLOOPcounter);
                            saveas(gcf,['Activation_spread' num2str(l+(j-1)*trialjump+(i-1)+K*trialjump) '.png'])
                            pichandle=Picture(['Activation_spread' num2str(l+(j-1)*trialjump+(i-1)+K*trialjump) '.png']);%save figure as picture
                            replace(pictureSlide,['Picture ' num],pichandle); %replace picture with name Picture X
                        end
                    end
                end
            end
        end
    end
end

 if SavetoPPT==1
     close(presentationObj); %close presentation to keep changes
     if ispc
         winopen('channelActivationSpread.pptx'); %open presentation (WINDOWS ONLY FUNCTION)
     end
 end
 