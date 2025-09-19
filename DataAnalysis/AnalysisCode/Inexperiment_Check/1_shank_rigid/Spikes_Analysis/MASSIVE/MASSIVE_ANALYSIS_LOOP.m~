function MASSIVE_ANALYSIS_LOOP(SubDir_Path)
%%for analysing data on massive and saving denoised files

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

folder = fileparts(which('MASSIVE_ANALYSIS_LOOP')); % Determines filepath to folder containing your .m file.
addpath(genpath(folder)); % Add that folder plus all subfolders to the path.
%fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd(SubDir_Path)%change working directory to where data is stored - currently manually input
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
dName='amplifier';
mem_check=dir('amplifier.dat');
T = mem_check.bytes ./ (2 * nChn * FS);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * FS);
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

denoiseIntan_sab(filepath, dName, T, par, Startpoint_analyse, Overall_time_to_analyse,FS);
cleanTrig_sabquick;
trig = loadTrig(0);
%% 2. Thresholds & Mu
allExtract_sab_1(dName,filepath,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');

%% 4. Calculate Structure of sorted trials according to IDs
starttrial=1;
trialjump=1;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
[IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
save('IDstruct.mat', 'IDstruct','baslinespikestruct')

%% 5. Calculates template of trials and spiking responses (Output in true electrode order)
[avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct, baslinespikestruct);
save('Averagetrialresponse.mat','avgnospT','stderrspktrial')


%% Heatmap of each condition
Spk_array=loadSpikes;
for ID=1:maxid
[array,arraycolmod,xvalues,yvalues]=OnlineHeatmap(nChn,ID,Spk_array);
saveidname=['HeatmapData_ID_',num2str(ID),'.mat'];
save(saveidname,'array','arraycolmod','xvalues','yvalues')
end
colormap(cmap)
maxsp=600;
for ID=144:256%maxid
    saveidname=['HeatmapData_ID_',num2str(ID),'.mat'];
    load(saveidname,'array','arraycolmod','xvalues','yvalues')

    figure (1)
    h=heatmap(xvalues,yvalues,arraycolmod);
    h.ColorLimits=[0 maxsp];
    pause(0.2)
end
%% 6. Work out where the stim electrode was and what layers were activated
% AMPInterestSingleLinePlot=4;
% depthdriven=1500-50;
% cutoffsp=50;
% ActivationDepth(AMPInterestSingleLinePlot,avgnospT,startpointseconds, secondstoanalyse,depthdriven,cutoffsp)
% %%
% 
%  AMP=loadAMP;
%  loadStimChn;
%  avgnospT_sps=(1000/(secondstoanalyse-startpointseconds)).*avgnospT;
%  cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
%  avgnostim=mean(avgnospT_sps(:,cell2mat(trialinfo(1:2:end,18))==-1),2);
%  numofsigperchn=zeros(nChn,1);
%  whichchn=zeros(nChn,1);
%  %testifsignificant=[(avgnospT_sps(:,cell2mat(trialinfo(16:2:28,18))==AMP(end))) (avgnospT_sps(:,cell2mat(trialinfo(14:2:28,18))==AMP(end-1))) (avgnospT_sps(:,cell2mat(trialinfo(14:2:28,18))==AMP(end-2)))];
%  %testsig_maxamp=testifsignificant(:,1:size(avgnostim,2));
%  for Chosenstimchn=1:length(stimChn)
%     desiredchanneltrial=find(cell2mat(trialinfo(:,2))==stimChn(Chosenstimchn)); %finds trials with desired initial electrode
%     desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds trials with matching recording electrode spacing
%     Desired_trialfirst=cell2mat(trialinfo(desiredchannel__singleampmath,1));%array of mathcning trial number
%      testifsignificant=(avgnospT_sps(:,Desired_trialfirst(end-2:end)));
%      chn_sig=zeros(1,nChn);
%      for i=1:nChn
%          [h,p]=ttest(testifsignificant(i,:)',avgnostim(i)','tail','right','Alpha',0.05);%one-tailed ttest determineing if the electrode has a mean significantly larger than the no stimulation trials with a 95% confidence level
%          chn_sig(i)=h;
%      end
%      chn_sig(isnan(chn_sig))=0;
%      numofsigperchn(stimChn(Chosenstimchn))=sum(chn_sig);
%      chn_sigboth(Chosenstimchn,:)=chn_sig;
%      whichchn(stimChn(Chosenstimchn))=1;
%  end
% save('sigchn.mat','numofsigperchn','whichchn','chn_sigboth')

% depthdriven=1000-50;%in um -50 frm tip
% 
% AMPInterestSingleLinePlot=4;%input in uA
% 
% truedatastruct=TrueData_heatmapLinecutFOURSHANK(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,IDstruct,startpointseconds, secondstoanalyse,depthdriven);
% 
% ratioALLnormalise;




cd(folder)
fprintf(['End of Analysis for: ' SubDir_Path newline])
end

