function TimeChangeinSpiking_SM(TimeBegin,TimeEnd, TimeStep, varargin)
% FOR SORTING AND PLOTTING
%   Plots the response curve according to trial change over time. 
%   Analyses spiking according to trial ID and Time. Creates a 3
%   dimensional plot with three axes. Intpus:
% startpointseconds - How long after the trigger do you want skip spike analysis(ms)? 
% secondstoanalyse - How long after the trigger do you want to analyse spikes for(ms)? 
%TimeStep -Timestep of spike analysis between TimeBegin and TimeEnd(ms)? 
%Eletrode of interest number
%Start trial, jump trial(usually = cond)
%End trial
%%
trig = loadTrig(0);
theseTrig = trig;
trialinfo=TrialInfo_sab;
trialinfo(1,:)=[];
cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition 
TrialParams = loadTrialParams;
maxtid=max(cell2mat(TrialParams(:,2)));
loadNREP;
try
    loadoriginalEND;
    if ((originalEND-1)/(n_REP_true*2))<maxtid
        maxtid=((originalEND-1)/(n_REP_true*2));
    end
catch
    %%code was made before loadoriginalEND
end
Loopnum=floor((TimeEnd-TimeBegin)/TimeStep);
loadStimChn;
if nargin==3 || nargin ==4
    x=zeros(Loopnum,maxtid);
    y=zeros(Loopnum,maxtid);
    z=zeros(Loopnum,maxtid);
    starttrial=1;
    trialjump=1;
    endtrial=maxtid;
elseif nargin==6
    electrode=cell2mat(varargin(1));
    starttrial=cell2mat(varargin(2));
    trialjump=cell2mat(varargin(3));
    endtrial=maxtid;
    x=zeros(Loopnum,((maxtid)/trialjump));
    y=zeros(Loopnum,((maxtid)/trialjump));
    z=zeros(Loopnum,((maxtid)/trialjump));
elseif nargin ==7
    electrode=cell2mat(varargin(1));
    starttrial=cell2mat(varargin(2));
    trialjump=cell2mat(varargin(3));
    endtrial=cell2mat(varargin(4));
    x=zeros(Loopnum,ceil((endtrial-starttrial)/trialjump));
    y=zeros(Loopnum,ceil((endtrial-starttrial)/trialjump));
    z=zeros(Loopnum,ceil((endtrial-starttrial)/trialjump));
end
    for counter=1:Loopnum
        [IDstruct,baslinespikestruct]=sortTrials_SM(TimeBegin+1+TimeStep*(counter-1),TimeStep*(counter)+TimeBegin,trig,0,starttrial,trialjump,endtrial);
        [avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct,baslinespikestruct);
        if nargin==3 %averages all trials for all electrodes
            averageallelectrodes=mean(avgnospT(:,starttrial:endtrial),1);
            z(counter,:)=(TimeStep*(counter)+TimeBegin).*ones(1,endtrial);
            x(counter,:)=1:endtrial;
            y(counter,:)=averageallelectrodes;
            x_namelabel='Maximum Current (uA)';
            title_namelabel='Average time varying changes according to trial (All electrodes)';
        elseif nargin==4 %%averages all trials for single electrode
            electrode=cell2mat(varargin(1));
            z(counter,:)=(TimeStep*(counter)+TimeBegin).*ones(1,length(avgnospT(electrode,:)));
            x(counter,:)=1:length(avgnospT(electrode,:));
            y(counter,:)=avgnospT(electrode,:);
            x_namelabel='Trials';
            title_namelabel=['Time varying changes according to trial for chn ' num2str(electrode)];
        elseif nargin==6 %allows user to select beginning trial and amount to jump by i.e. changes according to current strength for one condition 25/75
            jumparray=avgnospT(electrode,starttrial:trialjump:end);
            z(counter,:)=(TimeStep*(counter)+TimeBegin).*ones(1,length(jumparray));
            x(counter,:)=1:length(jumparray);
            y(counter,:)=jumparray;
            x_namelabel='Stimulation Amplitue(uA)';
            if (cond~=starttrial) && (starttrial~=1)
                title_namelabel=['Time varying changes for chn: ' num2str(electrode) '. Stimchn: ' num2str(stimChn(1)) ' ' num2str(stimChn(2)) ' @ ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))*100/(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))+cell2mat(trialinfo((starttrial+trialjump)*2,18)))) '/' num2str(cell2mat(trialinfo((starttrial+trialjump)*2,18))*100/((cell2mat(trialinfo((starttrial+trialjump)*2-1,18)))+cell2mat(trialinfo((starttrial+trialjump)*2,18))))];
            else
                title_namelabel=['Time varying changes for chn: ' num2str(electrode) 'Stimchn: ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,2)))];
            end
        elseif nargin==7 %allows user to select beginning trial, end trial and amount to jump by i.e. changes according to current strength for one condition 25/75
            endtrial=cell2mat(varargin(4));
            jumparray=avgnospT(electrode,starttrial:trialjump:endtrial);
            z(counter,:)=(TimeStep*(counter)+TimeBegin).*ones(1,length(jumparray));
            x(counter,:)=1:length(jumparray);
            y(counter,:)=jumparray;
            x_namelabel='Stimulation Amplitue(uA)';
            if (cond~=starttrial) && (starttrial~=1)
                title_namelabel=['Time varying changes for chn: ' num2str(electrode) '. Stimchn: ' num2str(stimChn(1)) ' ' num2str(stimChn(2)) ' @ ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))*100/(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))+cell2mat(trialinfo((starttrial+trialjump)*2,18)))) '/' num2str(cell2mat(trialinfo((starttrial+trialjump)*2,18))*100/((cell2mat(trialinfo((starttrial+trialjump)*2-1,18)))+cell2mat(trialinfo((starttrial+trialjump)*2,18))))];
            else
                title_namelabel=['Time varying changes for chn: ' num2str(electrode) '. Stimchn: ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,2)))];
            end
        end
    end
    if nargin>3
        if cond~=starttrial
            xticks_names=cell(1,floor((size(trialinfo,1)-starttrial)/(trialjump*2)));
        else
            xticks_names=cell(1,((size(trialinfo,1))/(trialjump*2)));
        end
        for trialnum=1:length(jumparray)
            if cond~=starttrial
                xticks_names(1,trialnum)={num2str(cell2mat(trialinfo(trialjump*trialnum*2+floor(starttrial/cond)*cond,18)))};
            else
                xticks_names(1,trialnum)={num2str(cell2mat(trialinfo(trialjump*trialnum*2,18)))};
            end
            
            if strcmp(cell2mat(xticks_names(1,trialnum)),'-1')
                xticks_names(1,trialnum)={num2str(0)};
            end
        end
    else
        xticks_names=cell(1,size(x,2));
        counter=0;
        for trialnum=1:size(x,2)
            if rem(trialnum,cond)==0
                counter=counter+1;
            end
            if (trialnum)==(ceil(cond/2))+cond*(counter)
                xticks_names(1,trialnum)={num2str(cell2mat(trialinfo(cond*(counter+1)*2,18)))};
                if strcmp(cell2mat(xticks_names(1,trialnum)),'-1')
                    xticks_names(1,trialnum)={num2str(0)};
                end
                %xticks_names(1,trialnum)={' '};
                %xticks_names(1,trialnum)=[];
            end 
        end
    end 
    figure
    hold on
    surf(x,z,y)
    xlabel(x_namelabel)
    zlabel('Average spiking')
    ylabel('Time (ms)')
    title(title_namelabel)
    if nargin>3
        xticks(1:size(y(1,:),2))
        %xticklabels(xticks_names(1,:))
        xlim([1 size(y(1,:),2)]) 
    else 
        xticks(1:size(y(1,:),2))
        xlim([1 (ceil((endtrial-starttrial)/trialjump))])
        for xline_number=1:cond:maxtid
            xline(xline_number,'Color','w','LineWidth',2)
        end
    end
    xticklabels(xticks_names(1,1:size(y(1,:),2)))
    colorbar
    caxis([0, max(y,[],'all')]);
end
    
%     for loopnumcount=1:maxtid/cond
%         hold on
%         xline((cond*loopnumcount)+1);
%     end
%     figure
%     hold on
%     surf(x,y,z)
%     xlabel(x_namelabel)
%     ylabel('Average spiking')
%     zlabel('Time (ms)')
%     title(title_namelabel)
%     xticks(1:size(xticks_names(1,:),2))
%     xticklabels(xticks_names(1,:))
%     xlim([1 size(xticks_names(1,:),2)])
%    