
%% Make sigmoid from chosen channels
%ONLY if there is single electrode stimulation
%can also plot raster - need to uncomment
electrodes='A-030,C-009';%"A-009'"% for all electrodes = Depth(1)
amplifier_channels=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);
electrodes='all';
% % all electrodes
% electrodes='';
% temp = ProbeMAP;
% for i=2:nChn+1
%     electrodes=strcat(electrodes,temp{i,6},',');
% end
% electrodes(end)=[];

Spk_array=OnlineSigmoidGenerator(electrodes,nChn);


%% firing rate response all channels 
trialinfo=loadTrialInfo;%use this to investigate which ID you want
ID=9;
amplifier_channels=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);
Chns='A-030,B-004';
Spk_array=OnlineFRresponse(Chns,nChn,ID);

%% heatmap
trialinfo=loadTrialInfo;%use this to investigate which ID you want
ID=5;
amplifier_channels=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);
array=OnlineHeatmap(nChn,ID);
 %% Raster threshold
 VA=2;%V1=1, V2=2
        if VA==2
            chn_range=[1 64];
        elseif VA==1
            chn_range=[65 128];
        end
uniqueStim=zeros(64,1);
filepath=pwd;
[filepathm,name,ext] = fileparts(filepath);
  ti=loadTrialInfo;
    load([name(1:end-14) '.sp.mat'],'sp')
    respsave=[];
    order=Depth(1);
   for ID=10:10:size(ti,1)-1%5:5:size(ti,1)-1
        VAsigcount=0;
        for chn=chn_range(1):chn_range(end) %=1:64%65:128
            resp=Online_raster(ID,sp{chn});
            if ~isempty(resp)&& any(resp>100)
                title(['chn: ' num2str(find(chn==order)) ' ID: ' num2str(ID)])
                %saveas(gcf,['Plot' num2str(chn) '_' num2str(ID) '.fig']);
                respsave=[respsave resp];
                VAsigcount=VAsigcount+1;
            end
            close all
        end
   end
   %% 
   for ID=1:10%5:5:size(ti,1)-1 
       for chn=chn_range(1):chn_range(end) %=1:64%65:128
           try
           openfig(['Plot' num2str(chn) '_' num2str(ID) '.fig']);
           catch
           end
       end
   end
 
 
%% Raster threshold %loop folders
VA=2;%V1=1, V2=2
respsave=[];
uniqueStim=zeros(64,1);
folders=dir;
for loopfold=3:size(folders,1)
    try
        cd([folders(loopfold).folder filesep folders(loopfold).name])
    catch
        continue
    end
    
    try
        ti=loadTrialInfo;
    load([folders(loopfold).name(1:end-14) '.sp.mat'],'sp')
    catch
        continue
    end
    %respsave=[];
    nostimID=cell2mat(ti([false; (cell2mat(ti(2:end,18))==-1)],1));
    order=Depth(1);
    for ID=5:5:size(ti,1)-1
        if VA==2
            chn_range=[1 64];
        elseif VA==1
            chn_range=[65 128];
        end
        VAsigcount=0;
        for chn=chn_range(1):chn_range(end) %=1:64%65:128
            resp=Online_raster(ID,sp{chn});
            if ~isempty(resp) && any(resp>100)
                title(['chn: ' num2str(find(chn==order)) ' ID: ' num2str(ID)])
                respsave=[respsave resp];
                VAsigcount=VAsigcount+1;
            end
            close all
        end
        if ti{ID+1,2}>64
             uniqueStim(ti{ID+1,2}-64)=VAsigcount;
        else
            uniqueStim(ti{ID+1,2})=VAsigcount;
        end
    end
end
%% percentage error
VA=2;%V1=1, V2=2

uniqueStim=zeros(64,1);
folders=dir;
for loopfold=3:size(folders,1)
    try
        cd([folders(loopfold).folder filesep folders(loopfold).name])
    catch
        continue
    end
    ti=loadTrialInfo;
    nostimID=cell2mat(ti([false; (cell2mat(ti(2:end,18))==-1)],1));
    load([folders(loopfold).name(1:end-14) '.sp.mat'],'sp')
    respsave=[];
    order=Depth(1);
    for ID=nostimID
        if VA==2
            chn_range=[1 64];
        elseif VA==1
            chn_range=[65 128];
        end
        VAsigcount=0;
        parfor chn=chn_range(1):chn_range(end) %=1:64%65:128
            resp=Online_raster(ID,sp{chn});
            if ~isempty(resp) && any(resp>100)
                title(['chn: ' num2str(find(chn==order)) ' ID: ' num2str(ID)])
                respsave=[respsave resp];
                VAsigcount=VAsigcount+1;
            end
            close all
        end
        if ti{ID+1,2}>64
             uniqueStim(ti{ID+1,2}-64)=VAsigcount;
        else
            uniqueStim(ti{ID+1,2})=VAsigcount;
        end
    end
end

%%
figure
 plot(-99:100,hist(respsave,1:200))
xlabel('Time(ms)')
ylabel('# of electrodes with FR above 5xSD')
legend('V2','V1')
%%
uniquestimsorted=reshape(uniqueStim,16,4); uniquestimsorted=uniquestimsorted(:,[1 3 4 2]);
figure
hm=heatmap(uniquestimsorted);
hm.YDisplayData=flip(hm.YDisplayData);
xlabel('shank #')
ylabel('electrode #')
title('# V1 electrodes active to V2 stimulation at that position')
%%
ID=50;
trig = loadTrig(0);
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
TP = loadTrialParams;
tID = find(cell2mat(TP(numsim_elect:numsim_elect:end,2)) == ID);
theseTrig = trig(tID);

filepath=pwd;
[filepathm,name,ext] = fileparts(filepath);
name = name(1:end-14);
 %fileID=fopen(['amplifier.dat'],'r');
 fileID=fopen([name '.mu_sab.dat'],'r');
                        shortbytes=2;
                        TimeBeforeTrig=0.010;
                        offset=theseTrig(27)*nChn*shortbytes-TimeBeforeTrig*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
                        ftell(fileID)
                        fseek(fileID,offset,'bof');
                        ftell(fileID)
                        vblankmu = fread(fileID,[nChn, (0.100*FS+TimeBeforeTrig*FS)],'short')./10; %plots one second from trigger and 250ms brefore
                        figure 
                        plot (-TimeBeforeTrig*1000:1000/FS:0.100*1000-1000/FS,vblankmu(61,:))
%% find time-locked activity
temp=dir('*.sp.mat');
load(temp.name)
nChn=128;
spiketimes=nan(1*10^5,128);
for chn=1:nChn
    spiketimes(1:length(sp{chn}(:,1)),chn)=sp{chn}(:,1);
end
s=spiketimes(:);
n=length(s);
[~,IA,~] = unique(s);
out = unique(s(setdiff((1:n),IA)));

