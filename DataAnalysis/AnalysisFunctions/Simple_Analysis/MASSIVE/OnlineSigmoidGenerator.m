function [Spk_array,spk_all,p]=OnlineSigmoidGenerator(varargin)
%chns input 'A-001,B-013'
%nChn is the total number of channels recorded simulataneously
%% Online sigmoid generation
tic
filepath = pwd;
dName='amplifier';
timestart=2;
timeend=90;
%user input channels of choice?? 
if nargin==0
    chns=input('Please input channels to analyse like so "A-001,B-016":\n','s');
    nChn=input('Please input the total number of channels recorded from all headstages\n');
elseif nargin==1
    fprintf('Not enough input parameters\n')
    chns=input('Please input channels to analyse like so "A-001,B-016":\n','s');
    nChn=input('Please input the total number of channels recorded from all headstages\n');
else
    chns=varargin(1);
    nChn=varargin{2};
end

if strcmp(chns,'all')
    chnnum=Depth(1);
else
    chnsplit=split(chns,",");
    chnnum=zeros(size(chnsplit,1),1);
    for i=1:size(chnsplit,1)
        if strcmp(chnsplit{i}(1),'A')
            chnnum(i)=str2double(chnsplit{i}(end-2:end))+1;
        elseif strcmp(chnsplit{i}(1),'B')
            chnnum(i)=str2double(chnsplit{i}(end-2:end))+33;
        elseif strcmp(chnsplit{i}(1),'C')
            chnnum(i)=str2double(chnsplit{i}(end-2:end))+65;
        elseif strcmp(chnsplit{i}(1),'D')
            chnnum(i)=str2double(chnsplit{i}(end-2:end))+97;
        end
    end
end
temp=dir('IDstruct_online.mat');
if ~isempty(temp)
load('IDstruct_online.mat','Spike_trialstruct','baslinespike_trialstruct')
else
 [Spike_trialstruct,baslinespike_trialstruct]=OnlineAllExtract(chnnum,nChn,timestart,timeend);
end
%% plot sigmoid -  this can only do one channel at the moment
%rejig so record electrode just loops through the outside
loadAMP_all;
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
chn=unique(trialinfo(:,2));
chn(chn==0)=[];
count=0;
allsingletrialIDs=[];
IDs=cell(length(chn),1);

for trialloop=1:numsim_elect:length(trialinfo(:,2))
    if sum(trialinfo(trialloop:trialloop+numsim_elect-1,17)~=-1)==1 || sum(trialinfo(trialloop:trialloop+numsim_elect-1,2)==0)==numsim_elect-1
        count=count+1;
        for chnloop=1:length(chn)
            if any(trialinfo(trialloop:trialloop+numsim_elect-1,2)==chn(chnloop))
                IDs{chnloop}(count,1:3)=[trialinfo(trialloop,17) (trialloop+numsim_elect-1)./numsim_elect trialinfo(trialloop,2)];
            end
        end
        %allsingletrialIDs(count)=trialloop;%index   trialinfo(trialloop,1);
    end
end

%%
% for stimchn=1:length(chn)
%         IDs{stimchn}=sortrows(IDs{stimchn});
%         IDs{stimchn}=IDs{stimchn}(any(IDs{stimchn},2),:);
% IDs{stimchn}=IDs{stimchn}((IDs{stimchn}(:,1)~=-1),:);
% end
%%
        tparams = dir('*_exp_datafile_*.mat');
        tparams = tparams.name;
        load(tparams,'simultaneous_stim');
spk_all=[];
p=[];
for recordchn=1:length(chnnum)
    for stimchn=1:length(chn)
        IDs{stimchn}=sortrows(IDs{stimchn});
        IDs{stimchn}=IDs{stimchn}(any(IDs{stimchn},2),:);

        for ID=1:size(IDs{stimchn},1)
            t=['ID_' num2str((IDs{stimchn}(ID,2)))];
            spkcount=zeros(40,1);
            spstim=zeros(40,1);
            bspk=zeros(40,1);
            for j=1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t),'AsArray',true),2)
                t2=['Trial_' num2str(j)];
                spkcount(j)=size(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1)-size(baslinespike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1);
                spstim(j)=size(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1);
                bspk(j)=size(baslinespike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1);
            end
            spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))])(ID)=mean(spkcount(1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t),'AsArray',true),2)))/((timeend-timestart)/1000);
            p.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))])(ID)=ranksum(spstim, bspk);
        end
    end
end

%toc
if length(chnnum)>16
    %%
    chnpos=Depth(1);
    figure
    tparams = dir('*_exp_datafile_*.mat');
    if isempty(tparams)
        TrialParams = [];
        return
    end
tparams = tparams.name;
load(tparams,'StimParams');
AMP_all = unique(cell2mat(StimParams(2:end,16)));
    %%plot sigmoid
    save_minp=100;
    countsigelect=zeros(size(AMP_all));
    sz=ceil(length(chnnum)/16);
    for stimchn=1:length(chn)
        figure
        for recordchn=1:length(chnnum)
            subplot(sz,16,(recordchn))
            hold on
            
            plot(IDs{stimchn}(:,1),spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]))
            if sum(p.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))])<0.05)==0
                minp=0;
            elseif sum(p.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))])<0.05)~=length(IDs{stimchn}(:,1))
                minp=length(p.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]))...
                    -(find(fliplr(p.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))])<0.05)==0,1,'first')-1)+1;% finds minimum current level that significance occurs
            else
                minp=1;
            end
            if minp~=0
                if size(IDs{stimchn},1)<minp
                    minp=size(IDs{stimchn},1);
                end
                countsigelect(AMP_all==IDs{stimchn}(minp,1))=countsigelect(AMP_all==IDs{stimchn}(minp,1))+1;
                if minp<save_minp
                    save_minp=minp;
                end
            end
             %% model fit
%             X=IDs{stimchn}(:,1);
%             y=spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]);
%             modelfun = @(b,X) b(1)./(1+exp((b(2)-X)/b(3))); %SIGMOIDAL  b1 is asymptote, b2=x50 b3=slope  This is boltzmann sigmoid
%             beta0 = [100 1 1]; %SIGMOIDAL%% fit model
%             try
%                 mdl = fitnlm(X,y',modelfun,beta0);
%                 beta=mdl.Coefficients.Estimate;
%                 if ~any(beta<0)
%                     YFIT = beta(1)./(1+exp((beta(2)-X)/beta(3)));
%                     hold on
%                     plot(X,YFIT,'--k')
%                 else 
%                     warning('model not fit')
%                 end
%             catch
%                 warning('model not fit')
%             end
        end
        xlabel('Current (\muA)')
        ylabel('Sp/s')
    end
    %%
    countsigcumlative=0;
    figure
    hold on
    AMP_all(AMP_all==-1)=0;
    for i=1:length(AMP_all)
        countsigcumlative=countsigelect(i)+countsigcumlative;
        scatter(AMP_all(i),countsigcumlative,'filled','k');

    end
    ylabel('# significant electrodes')
    xlabel('Current (\muA)')
    
    %% collapsed sigmoid

    for stimchn=1:length(chn)
        stimchn_spkarray=zeros(nChn,length(IDs{stimchn}(:,1)));
        for recordchn=1:length(chnnum)
            stimchn_spkarray(recordchn,1:length(IDs{stimchn}(:,1)))=spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]);%IDs{stimchn}(:,1)
        end
           figure
            subplot(1,2,1)
            hold on
        plot(IDs{stimchn}(:,1),mean(stimchn_spkarray(1:64,:),1))
        xlabel('Current (\muA)')
        ylabel('Sp/s')
        title('V2 array')
                    subplot(1,2,2)
            hold on
        plot(IDs{stimchn}(:,1),mean(stimchn_spkarray(65:128,:),1))
        xlabel('Current (\muA)')
        ylabel('Sp/s')
        title('V1 array')
    end

else
    %%plot sigmoid
    for recordchn=1:length(chnnum)
     %figure
        hold on
        for stimchn=1:length(chn)
            subplot(4,16,chn(stimchn))
            plot(IDs{stimchn}(:,1),spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]))
            %% model fit
%             X=IDs{stimchn}(:,1);
%             y=spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]);
%             modelfun = @(b,X) b(1)./(1+exp((b(2)-X)/b(3))); %SIGMOIDAL  b1 is asymptote, b2=x50 b3=slope  This is boltzmann sigmoid
%             beta0 = [100 1 1]; %SIGMOIDAL%% fit model
%             try
%                 mdl = fitnlm(X,y',modelfun,beta0);
%                 beta=mdl.Coefficients.Estimate;
%                 if ~any(beta<0)
%                     YFIT = beta(1)./(1+exp((beta(2)-X)/beta(3)));
%                     hold on
%                     plot(X,YFIT,'--k')
%                 else 
%                     warning('model not fit')
%                 end
%             catch
%                 warning('model not fit')
%             end
      xlabel('Current (\muA)')
        ylabel('Sp/s')
        end
%         xlabel('Current (\muA)')
%         ylabel('Sp/s')
%         title(chnsplit{recordchn})
    end
end
toc

%% raster
% ID=34; % there will be an increase at 100ms due to offset
% for recordchn=1:length(chnnum)
% Online_raster(ID,Spk_array{recordchn})
% end




end