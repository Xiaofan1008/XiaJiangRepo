%% For quick analysis during experiment

% Parameters to alter
artefact=-100; %removes spikes below this threshold
artefact_high=100; %removes spikes above this threshold
startpointseconds=11; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=20; %How long after the trigger do you want to analyse spikes for(ms)? 
printspiking=0;
par=0;

%% 1. Blank stimulus
nChn=32;
FS=30000;
filepath = pwd;
dName='amplifier';
vFID = fopen([filepath '\' dName '.dat'],'r');
mem = memory;
T = mem.MaxPossibleArrayBytes ./ (2 * 32 * 30000);
fileinfo = dir([filepath '\' dName '.dat']);
t_len = fileinfo.bytes/(32 * 2 * 30000);
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256; 
end
info = fileinfo.bytes/2;
nL = (ceil(info / (nChn*FS*double(T)))+1);
vblank=[];
BREAK = 1;
N=1;
denoiseIntan_sab(filepath,dName,T,par)
trig = loadTrig(0);
theseTrig = trig;

%% 2. Thresholds & Mu
allExtract_sab_1(dName,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');
%% 4. Calculate Structure of sorted trials according to IDs
[IDstruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking);

%% 5. Calculates template of trials and spiking responses (Output in true electrode order)
[avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct);
cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition 
avgnostim=(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1));
testifsignificant=(avgnospT(:,cell2mat(trialinfo(1:2:end,18))~=-1));
checknumnovsstim=size(testifsignificant,2)/size(avgnostim,2);
[h,p]=ttest(avgnostim',testifsignificant(1:nChn,1:size(avgnostim,2))');
h(isnan(h))=0;

%% Plot response curve change over time (3D graph)
TimeBegin=2; %How long after the trigger do you want skip spike analysis(ms)? 
TimeEnd=30; %How long after the trigger do you want to stop spike analysis(ms)? 
TimeStep=2; %Timestep of spike analysis between TimeBegin and TimeEnd(ms)? 
TimeChangeinSpiking_SM(TimeBegin, TimeEnd, TimeStep);

%% 6. Print all trial details
Legend=cell(nChn,1);
loadStimChn;
loadStimParams;
stimamplitudes=unique(cell2mat(StimParams(2:end,16)));
fprintf('The current amplitude in ascending order is: \n');
fprintf('%d\n',stimamplitudes);
fprintf('The stimulation channels were: \n');
fprintf('%d\n',stimChn);

%% 7. Plot stimulating electrode response curve
StimElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial);

%% 8. plotting average electrode response for all electrodes classes as significant and not significant
AllElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,h)

%% 9. plotting electrodes between stim and stim electrodes
InbetweenElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial);









% lenstimamplitudes=length(stimamplitudes);
% chnstep=1;
% %%
% for nextamp=1:maxtid/cond
%     figure
%     for condchange=1:cond
%         subplot(1,cond,condchange)
%         hold on
%         plot(10:20,(avgnospT(10:20,condchange+(cond*(nextamp-1)))),'b')
%         ylim([0 yvalmax])
%         xlabel('Electrode number')
%         ylabel('Average number of threshold crossings')
%         if (cell2mat(trialinfo((nextamp-1)*(cond*2)+(condchange-1)*2+2,2))==0)
%             title(['Stimulating Electrode: ', (trialinfo((nextamp-1)*(cond*2)+(condchange-1)*2+1,2)), ' @ ', (trialinfo((nextamp-1)*(cond*2)+(condchange-1)*2+1,18)),'uA'])
%         else
%             title(['Stimulating Electrode: ', (trialinfo((nextamp-1)*(cond*2)+(condchange-1)*2+1,2)), ' @ ', (trialinfo((nextamp-1)*(cond*2)+(condchange-1)*2+1,18)),'uA & ', (trialinfo((nextamp-1)*(cond*2)+(condchange-1)*2+2,2)),' @ ', (trialinfo((nextamp-1)*(cond*2)+(condchange-1)*2+2,18)), 'uA'])
%         end
%     end
% 
% end



%%
%
% figure
% for trialid=1:length(reorder16)
%     %subplot(1,length(reorder16),trialid)
%     hold on
%     scatter(reorder16(trialid,1),avgnospT(electrode,reorder16(trialid,2)),'b')
% end
%%
% for i=1:length(stimChn)
%     for chn=1:chnstep:nChn
%         if chn<4
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'r')
%         elseif chn<7
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'k')
%         elseif chn<10
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'m')
%         elseif chn<13
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'g')
%         elseif chn<16
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'c')
%         elseif chn<19
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'b')
%         elseif chn<22
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'k')
%         elseif chn<25
%             subplot(1,length(stimChn),i)
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'r')
%         elseif (chn<28)
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'k')
%         else
%             subplot(1,length(stimChn),i)
%             hold on
%             plot((avgnospT(emap(chn),i*length(stimamplitudes)-length(stimamplitudes)+1:length(stimamplitudes)*i)),'b')
%         end
%     end
%     ylim([0 yvalmax])
%     xlim([1 lenstimamplitudes])
%     xticks(1:lenstimamplitudes)
%     xticklabels(num2str(stimamplitudes))
%     if (cell2mat(TrialParams(2,1))==cell2mat(TrialParams(3,1)))
%         if i<length(stimChn)/2+1
%             title(['Stimulating Electrode: ', (num2str(stimChn(i*2-1)))])
%         else
%             title(['Stimulating Electrode: ', (num2str(stimChn(i*2-1-length(stimChn)))),' & ', num2str(stimChn(i*2-length(stimChn)))])
%         end
%     else
%         title(['Stimulating Electrode: ', (num2str(stimChn(i)))])
%     end
%     xlabel('Stim current (uA)')
%     ylabel('Average number of threshold crossings')
% end
%  Legend=string(num2cell(1:chnstep:nChn));
% legend(Legend)



%%
%electrode=16;

%%plotting stimulating electrodes
% plotcounter=1;
% yvalmax=30;
% figure
% hold on
% for count=1:length(stimChn)
%     electrode=stimChn(count);
% 
%     equal16=find(cell2mat(trialinfo(:,2))==electrode);
%     %single stim
% %     if (equal16(end)==length(trialinfo))
% %         equal16(end)=[];
% %     end
% %     equal16l=equal16((cell2mat(trialinfo(equal16+1,2))==0),1);
%     %
%     equal16l=equal16;
%     amp16=cell2mat(trialinfo(equal16l,18));
%     trials16=cell2mat(trialinfo(equal16l,1));
%     reorder16=[amp16 trials16];
%     reorder16=sortrows(reorder16);
%    % trials in order of amplitude
% 
%     ampnumber=unique(reorder16(:,1));
%     for trialid=1:length(ampnumber)
%         for i=1:length(reorder16)
%             if (ampnumber(trialid,1)==reorder16(i,1))
%                 trialamp(i)=avgnospT(electrode,reorder16(i,2));
%                 errortrial(i)=stderrspktrial(electrode,reorder16(i,2));
%             end
%         end
%         trialamp(trialamp==0)=[];
%         errortrial(errortrial==0)=[];
%         avgtrialamp(trialid)=mean(trialamp);
%         avgerrortrial(trialid)=mean(errortrial);
%         clear trialamp
%         clear errortrial
%     end
%     
%     hold on
%     subplot(1,length(stimChn),plotcounter)
%     errorbar(ampnumber,(avgtrialamp),avgerrortrial) 
%     ylabel('Average number of spikes')
%     xlabel('Stimulating current(uA)')
%     title(['Electrode ', (num2str(electrode))])
%     ylim([0 yvalmax])
%     xlim([-2 17])
%     plotcounter=plotcounter+1;
%     clear trialamp
%     clear errortrial
%     clear avgtrialamp
%     clear avgerrortrial
%     
% end


% 
% %threshold crossings in individual trials add then average
% TrialParams = loadTrialParams;
% loadNREP;%number repeats
% spike =0;
% 
% %fileID=fopen('amplifier_dn.dat','r');
% 
% 
% 
% name = pwd;
% name = strsplit(name,'\');
% name = name{end};
% name = name(1:end-14);
% load([name '.sp.mat'])
% chk=1;
% maxtid=max(cell2mat(TrialParams(:,2)));
% spike_trig_chn=[];
% nospI=[];
% IDstruct=[];
% loopcount=0;
% flag=0;
% emap=[22,27,18,28,13,4,9,29,19,3,14,2,12,6,17,30,20,25,15,31,11,1,16,0,21,5,10,26,23,7,8,24];
% emap=emap+1;
% %% remove values of common spike times across the array
% C1 = intersect(sp{1}(:,1),sp{2}(:,1));
% C2 = intersect(sp{24}(:,1),sp{6}(:,1));
% C3 = intersect(sp{3}(:,1),sp{5}(:,1));
% test1= intersect(C1,C2);
% test2= intersect(C2,C3);
% check=intersect(test1,test2);
% for chncount=1:nChn
%     [C,r,c]=(intersect(sp{chncount}(:,1),check));
%     sp{chncount}(r,:)=[];
% end
% %% Sort into trial IDs and plot spikes for trial 1
% 
% for tID=1:maxtid
%     loopcount=1+loopcount;
%     fprintf('%d\n',loopcount)
%     if isempty(TrialParams(cell2mat(TrialParams(1:end,2)) == tID))
%         fprintf('No trial ID: %d\n',tID)
%     else
%         TrialParamstID = find(cell2mat(TrialParams(1:end,2)) == tID); %identifies trial row matching trial ID
%         TrialParamstID(1:2:end)=[];
%         TrialParamstID=TrialParamstID./2;
%         %                 if TrialParamstID(end)>=2697
%         %                     TrialParamstID(end)=[];
%         %                 end
%         trigtID = trig(TrialParamstID)./(FS/1000);
%         %                 trigtID=trigtID(((loopcount-1)*T*FS-FS*secondstoanalyse)<trigtID);
%         %                 trigtID=trigtID(trigtID<=(T*FS*loopcount-FS*secondstoanalyse));
%         nTrig = length(trigtID);
%         for indT=1:nTrig
%             %                     if (trigtID(indT)-((T*FS)*(loopcount-1))+FS*secondstoanalyse)<=length(vblankmu)
%             %                         if loopcount~=1
%             %                             v = tempblankmu(1:nChn, (trigtID(indT)-(T*FS*(loopcount-1)-FS*secondstoanalyse)+FS*startpointseconds):trigtID(indT)-(T*FS*(loopcount-1)-FS*secondstoanalyse)+FS*secondstoanalyse); % reading (secondstoanalyse)s following trial
%             %                         else
%             %                             v = tempblankmu(1:nChn, (trigtID(indT)+FS*startpointseconds):(trigtID(indT)+FS*secondstoanalyse)); % reading (secondstoanalyse)s following trial
%             %                         end
%             %                     else
%             %                         error('should not occur');
%             %                     end
%             for chsp=1:1:nChn
%                 v = sp{chsp};
%                 spikedetailstrig=v(((v(:,1)>trigtID(indT)+startpointseconds*(FS/1000))&(v(:,1)<trigtID(indT)+secondstoanalyse*(FS/1000))),:);
%                 if ~isempty(spikedetailstrig)
%                     if any(spikedetailstrig(:,2:end)>(150))
%                         
%                         fileID=fopen('rat_001_005_005.mu.dat','r');
%                         shortbytes=2;
%                         offset=trig(TrialParamstID(indT))*nChn*shortbytes;%trig(TrialParamstID(tID))*nChn*shortbytes;
%                         ftell(fileID)
%                         fseek(fileID,offset,'bof');
%                         ftell(fileID)
%                         vblankmu = fread(fileID,[nChn, 1*FS],'short')./10;
%                         figure (400)
%                         plot (vblankmu(chsp,:))
%                         title(['channel' num2str(chsp)])
%                         xlabel('Samples')
%                         ylabel('uV')
%                         fclose(fileID);
%                         %
%                         %                                     fileID=fopen('rat_001_005_005.mu.dat','r');
%                         %                                     shortbytes=2;
%                         %                                     offset=0*nChn*shortbytes;%trig(TrialParamstID(tID))*nChn*shortbytes;
%                         %                                     ftell(fileID)
%                         %                                     fseek(fileID,offset,'bof');
%                         %                                     ftell(fileID)
%                         %                                     vblankmu = fread(fileID,[nChn, 30*FS],'short')./10;
%                         %                                     figure (403)
%                         %                                     plot (vblankmu(32,:))
%                         %                                     title(['channel' num2str(chsp)])
%                         %                                     xlabel('Samples')
%                         %                                     ylabel('uV')
%                         %                                     fclose(fileID);
%                         
%                     end
%                     spike=size(spikedetailstrig,1);
%                     if (tID<3) && (indT<=10)
%                         figure(emap(chsp))
%                         hold on
%                         plot(1*1000/FS:1000/FS:49*1000/FS,spikedetailstrig(:,2:50))
%                         title(['Channel ' num2str(emap(chsp))])
%                         ylabel('Spike amplitude (uV)')
%                         xlabel('Time (ms)')
%                     end
%                 else
%                     spike=0;
%                 end
%                 nospI(chsp,indT)=spike;%(length(r1)+spike);
%                 spike=0;
%             end
%         end
%         IDstruct=StructureIDgeneration(nospI, tID, IDstruct);
%         nospI=[];
%     end
%     
% end



% avgnospT=[];
% stderrspktrial=[];
% trialinfo={};
% TrialParams=loadTrialParams;
% loadStimParams;
% for ID=1:maxtid
%     avgnospT=meanstruct(avgnospT, ID, IDstruct);
%     stderrspktrial=stderrorstruct(stderrspktrial, ID, IDstruct);
%     num=find(cell2mat(TrialParams(1:end,2)) == ID);
%     trialinfo(ID+(ID-1),:)=[{ID},TrialParams(num(1),3), StimParams(num(1)+1,:)];
%     trialinfo((ID+1)+(ID-1),:)=[{ID},TrialParams(num(2),3), StimParams(num(2)+1,:)];
% end
% cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;
% tempavg=avgnospT;
% avgnospT=avgnospT(emap,:);
% stderrspktrial=stderrspktrial(emap,:);
% avgnostim=mean(avgnospT(:,1:cond),2);

% 
% [Mufilt,Lfpfilt] = generate_Filters;
% MuNf = length(Mufilt);
% LfpNf = length(Lfpfilt);
% artchk = zeros(1,nChn);
% munoise = cell(1,nChn);
% nRun = 0;
% nSam = fileinfo.bytes/(nChn * 2);
% tRun = ceil(nSam / FS / 10) + 1; % Number of 10 second chunks of time in the data
% thresh = cell(1, nChn);
% ARTMAX = 1.5e3;
% fileID = fopen([filepath '\' dName '_dn_sab.dat'],'r');
% %fileID=fopen('rat_001_005_005.mu.dat','r');
% threshfac = -4.5;
% dispstat(sprintf('\n Thresholding . . .'),'keepthis','n')
% 
% while sum(artchk) < nChn
%     v = fread(fileID,[32, (20*FS)], 'int16');
%     %v = vblank(1:nChn, 1:FS*20); % reading 20s
%     nRun = nRun + 1;
%     dispstat(sprintf('Progress %03.2f%%',(100*(nRun/tRun))),'timestamp');
%     if ~isempty(v)
%         % Checks for artifact
%         for iChn = 1:nChn
%             if isempty(thresh{iChn})
%                 munoise{iChn} = [];
%                 mu = conv(v(iChn,:),Mufilt);
%                 if max(abs(mu(MuNf+1:end-MuNf))) < ARTMAX
%                     artchk(iChn) = 1;
%                     sd = median(abs(mu(MuNf+1:end-MuNf)))./0.6745;
%                     thresh{iChn} = threshfac*sd;
%                 else
%                     munoise{iChn} = [munoise{iChn} mu(MuNf+1:end-MuNf)];
%                 end
%             end
%         end
%     else
%         for iChn = 1:nChn
%             if isempty(thresh{iChn}) % If threshold is still 0 - noisy channel
%                 artchk(iChn) = 1;
%                 sd = median(abs(munoise{iChn}))./0.6745;
%                 thresh{iChn} = threshfac*sd;
%             end
%         end
%     end
%     
%     disp(['Total recording time: ' num2str(nSam / FS) ' seconds.']);
%     disp(['Time analysed per loop: ' num2str(T) ' seconds.']);
%     
% end
%% 3. Filtering MU
% N=1;
% Ninput=1e6; %why this value
% vblankmu=[];
% chk=1;
% ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);
% dispstat('','init');
% dispstat(sprintf('Processing MU . . .'),'keepthis','n');
%     while (chk && N < Ninput)
%         if (FS*N*T) >= length(vblank)
%           data = vblank(1:nChn, 1+FS*(N-1)*T:end);
%             chk=0;
%         else 
%             data = vblank(1:nChn, 1+FS*(N-1)*T:FS*N*T); 
%         end
%         Ndata = size(data,2);
%         mu=zeros(nChn, Ndata);
%         for iChn = 1:nChn
%             flip_data = fliplr(data(iChn,:));
%             tmp = conv(flip_data,Mufilt);
%             mu(iChn,:) = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
%         end
%         vblankmu = [vblankmu, mu];
%         
% 
%         dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
%         N = N + 1; 
%         if (size(data,2) < FS * T)
%             chk = 0;
%         end
%     end

% %% 4. Response curves
% 
% loadThreshold;
% 
% %single electrode stim first in trial order
% startpointseconds=2/1000;
% secondstoanalyse=11/1000;
% %threshold crossings in individual trials add then average
% TrialParams = loadTrialParams;
% loadNREP;%number repeats
% spike =0;
% %fileID=fopen('amplifier_dn.dat','r');
% %fileID=fopen('rat_001_005_005.mu.dat','r');
% name = pwd;
% name = strsplit(name,'\');
% name = name{end};
% name = name(1:end-14);
% fileID = fopen([name '.mu_sab3.dat'],'r');
% chk=1;
% maxtid=max(cell2mat(TrialParams(:,2)));
% spike_trig_chn=[];
% nospI=[];
% IDstruct=[];
% loopcount=0;
% flag=0;
% while (chk) 
%     loopcount=1+loopcount;
%     ftell(fileID)
%     if loopcount ~=1
%         tempvblankmu=vblankmu(1:32, end-FS*secondstoanalyse:end);
%     else
%         tempvblankmu=[];
%     end
%     vblankmu = fread(fileID,[nChn, FS*T],'short')./10;
%     tempblankmu=[tempvblankmu vblankmu];
%     %tempblankmu=tempblankmu-smoothdata(tempblankmu);
%     % C = textscan(fileID,'%s %s %f32 %d8 %u %f %f %s %f');
%     if ~isempty(vblankmu)
%         for tID=1:maxtid
%             if isempty(TrialParams(cell2mat(TrialParams(1:end,2)) == tID))
%                 fprintf('No trial ID: %d\n',tID)
%             else
%                 TrialParamstID = find(cell2mat(TrialParams(1:end,2)) == tID); %identifies trial row matching trial ID
%                 TrialParamstID(1:2:end)=[];
%                 TrialParamstID=TrialParamstID./2;
% %                 if TrialParamstID(end)>=2697
% %                     TrialParamstID(end)=[];
% %                 end
%                 trigtID = trig(TrialParamstID);
%                 trigtID=trigtID(((loopcount-1)*T*FS-FS*secondstoanalyse)<trigtID);
%                 trigtID=trigtID(trigtID<=(T*FS*loopcount-FS*secondstoanalyse));
%                 nTrig = length(trigtID);
%                 for indT=1:nTrig
%                     if (trigtID(indT)-((T*FS)*(loopcount-1))+FS*secondstoanalyse)<=length(vblankmu)
%                         if loopcount~=1
%                             v = tempblankmu(1:nChn, (trigtID(indT)-(T*FS*(loopcount-1)-FS*secondstoanalyse)+FS*startpointseconds):trigtID(indT)-(T*FS*(loopcount-1)-FS*secondstoanalyse)+FS*secondstoanalyse); % reading (secondstoanalyse)s following trial
%                         else
%                             v = tempblankmu(1:nChn, (trigtID(indT)+FS*startpointseconds):(trigtID(indT)+FS*secondstoanalyse)); % reading (secondstoanalyse)s following trial
%                         end
%                     else
%                         error('should not occur');
%                     end
%                     for chsp=1:1:nChn
%                             [Sp,Spktimes] = spikeextract(v(chsp,:),thresh{chsp},FS,-200);
%                             if (~isempty(Sp) && any(Sp(:)<(-200)))
%                                 flag=flag+1;
%                             end
%                             if Spktimes~=0
%                                 spike=length(Spktimes);
%                                 if tID==27
%                                     figure(chsp)
%                                     hold on
%                                     plot(1:49,Sp(:,1:49))
% 
%                                 end
%                             else
%                                 spike=0;
%                             end
%                             nospI(chsp,indT)=spike;%(length(r1)+spike);
%                             spike=0;
%                     end
%                 end
%                 IDstruct=StructureIDgeneration(nospI, tID, IDstruct);
%                  nospI=[];
%             end
%             
%         end
%     else
%         chk=0;
%     end
% end
% fclose(fileID);