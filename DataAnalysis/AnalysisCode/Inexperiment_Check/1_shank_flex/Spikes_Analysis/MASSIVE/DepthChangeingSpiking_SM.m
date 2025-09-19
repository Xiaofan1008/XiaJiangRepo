function DepthChangeingSpiking_SM(avgnospT, chosen_trials,fignum,depthdriven, varargin)

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
loadAMP_all;
AMP=AMP_all';
trig = loadTrig(0);
theseTrig = trig;
trialinfo=loadTrialInfo;
trialinfo(1,:)=[];
cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
TrialParams = loadTrialParams;
maxtid=max(cell2mat(TrialParams(:,2)));
loadStimChn;
filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=loadMapNum;
    if E_Mapnumber>0
        %nChn=64;
        nChn=16; %per shank
    else
        nChn=32;
    end
end
if nargin ==4
    endtrial=cell2mat(varargin(1));
else
    endtrial=maxtid;
end

x=depthdriven:-50:depthdriven-(nChn-1)*50;%ones(nChn,length(chosen_trials));
x=x';
x=repelem(x,1,length(chosen_trials));
if nargin==6
    y=cell2mat(varargin(2));
else
    y=AMP(1:length(chosen_trials));
end
if ~isrow(y)
    y=y';
end
%y=cell2mat(trialinfo(chosen_trials.*2,18))';
y(y==-1)=0;
y=repelem(y,nChn,1);
z=avgnospT(:,1:length(chosen_trials));



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
    figure(fignum)
    hold on
    surf(y,x,z,'FaceColor','interp','EdgeColor','interp')
    xlabel('Amplitude (uA)')
    zlabel('Average spiking')
    ylabel('Estimated depth (um)')
    if cell2mat(trialinfo(chosen_trials(2)*2,2))~=0%rem(starttrial,cond)&&(rem((starttrial-1),cond))
        title_namelabel=['Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))*100/(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '/' num2str(cell2mat(trialinfo((chosen_trials(2)*2),18))*100/((cell2mat(trialinfo((chosen_trials(2)*2)-1,18)))+cell2mat(trialinfo((chosen_trials(2)*2),18))))];
    else
        title_namelabel=['Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
    end
    if nChn==32
    yline((depthdriven-50*(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))-1)),'Color','r','Linewidth',cell2mat(trialinfo((chosen_trials(2)*2)-1,18))*100/(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))+cell2mat(trialinfo((chosen_trials(2)*2),18)))/35,'Alpha',1)
    if cell2mat(trialinfo(chosen_trials(2)*2,2))~=0
        yline((depthdriven-50*(cell2mat(trialinfo(chosen_trials(2)*2,2))-1)),'Color','r','Linewidth',cell2mat(trialinfo((chosen_trials(2)*2),18))*100/((cell2mat(trialinfo((chosen_trials(2)*2)-1,18)))+cell2mat(trialinfo((chosen_trials(2)*2),18)))/35,'Alpha',1)
    end
    end
    title(title_namelabel)
    set(gca, 'YDir','reverse')
    ylim([depthdriven-50*(nChn-1) depthdriven])
    hcb.Title.String = "Sp/s";
    hcb.Title.Rotation=270;
    hcb.Title.Position= [40,130];
    hcb.Title.FontSize=11;
    caxis([0, 400])%max(z,[],'all')])

%     for loopnumcount=1:maxtid/cond
%         hold on
%         xline((cond*loopnumcount)+1);
%     end
 
  
    
    %%
end








% function DepthChangeingSpiking_SM(avgnospT, starttrial, trialjump, varargin)
% 
% %   Plots the response curve according to trial change over time. 
% %   Analyses spiking according to trial ID and Time. Creates a 3
% %   dimensional plot with three axes. Intpus:
% % startpointseconds - How long after the trigger do you want skip spike analysis(ms)? 
% % secondstoanalyse - How long after the trigger do you want to analyse spikes for(ms)? 
% %TimeStep -Timestep of spike analysis between TimeBegin and TimeEnd(ms)? 
% %Eletrode of interest number
% %Start trial, jump trial(usually = cond)
% %End trial
% %%
% trig = loadTrig(0);
% theseTrig = trig;
% trialinfo=loadTrialInfo;
% trialinfo(1,:)=[];
% cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
% TrialParams = loadTrialParams;
% maxtid=max(cell2mat(TrialParams(:,2)));
% loadStimChn;
% filepath = pwd;
% fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
% fileinfo = dir([filepath filesep 'info.rhs']);
% if (datetime(fileinfo.date) < fourShank_cutoff)
%     nChn=32;
%     E_Mapnumber=0;
% else
%     E_Mapnumber=loadMapNum;
%     if E_Mapnumber>0
%         nChn=64;
%     else
%         nChn=32;
%     end
% end
% if nargin ==4
%     endtrial=cell2mat(varargin(1));
% else
%     endtrial=maxtid;
% end
%     if rem((endtrial-starttrial),trialjump)&&(starttrial>trialjump)
%         x=zeros(nChn,ceil((endtrial-starttrial)/trialjump));
%         y=zeros(nChn,ceil((endtrial-starttrial)/trialjump));
%         z=zeros(nChn,ceil((endtrial-starttrial)/trialjump));
%         for trial=1:ceil((endtrial-starttrial)/trialjump)
%             z(:,trial)=avgnospT(:,starttrial+(trial-1)*trialjump);
%             y(:,trial)=cell2mat(trialinfo(((trialjump*(trial))*2),18));
%             if y(:,trial)==-1
%                 y(:,trial)=0;
%             end
%             x(:,trial)=1:nChn;
%         end
%     elseif (starttrial<trialjump)
%         x=zeros(nChn,floor((endtrial)/trialjump));
%         y=zeros(nChn,floor((endtrial)/trialjump));
%         z=zeros(nChn,floor((endtrial)/trialjump));
%         for trial=1:floor((endtrial)/trialjump)
%             z(:,trial)=avgnospT(:,starttrial+(trial-1)*trialjump);
%             y(:,trial)=cell2mat(trialinfo(((trialjump*(trial))*2),18));
%             if y(:,trial)==-1
%                 y(:,trial)=0;
%             end
%             x(:,trial)=1:nChn;
%         end
%     else
%         x=zeros(nChn,floor((endtrial-starttrial)/trialjump)+1);
%         y=zeros(nChn,floor((endtrial-starttrial)/trialjump)+1);
%         z=zeros(nChn,floor((endtrial-starttrial)/trialjump)+1);
%         for trial=1:floor((endtrial-starttrial)/trialjump)+1
%             z(:,trial)=avgnospT(:,starttrial+(trial-1)*trialjump);
%             y(:,trial)=cell2mat(trialinfo(((trialjump*(trial))*2),18));
%             if y(:,trial)==-1
%                 y(:,trial)=0;
%             end
%             x(:,trial)=1:nChn;
%         end
%     end
%    
% 
% %     figure
% %     hold on
% %     surf(x,y,z)
% %     xlabel(x_namelabel)
% %     ylabel('Average spiking')
% %     zlabel('Time (ms)')
% %     title(title_namelabel)
% %     xticks(1:size(xticks_names(1,:),2))
% %     xticklabels(xticks_names(1,:))
% %     xlim([1 size(xticks_names(1,:),2)])
% %     
%     figure
%     hold on
%     surf(y,x,z)
%     xlabel('Amplitude (uA)')
%     zlabel('Average spiking')
%     ylabel('Channel number')
%     if cell2mat(trialinfo((starttrial+trialjump)*2,2))~=0%rem(starttrial,cond)&&(rem((starttrial-1),cond))
%         title_namelabel=['Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,2))) ' ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2,2))) ' @ ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))*100/(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))+cell2mat(trialinfo((starttrial+trialjump)*2,18)))) '/' num2str(cell2mat(trialinfo((starttrial+trialjump)*2,18))*100/((cell2mat(trialinfo((starttrial+trialjump)*2-1,18)))+cell2mat(trialinfo((starttrial+trialjump)*2,18))))];
%     else
%         title_namelabel=['Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,2)))];
%     end
%     yline(cell2mat(trialinfo((starttrial+trialjump)*2-1,2)),'Color','r','Linewidth',cell2mat(trialinfo((starttrial+trialjump)*2-1,18))*100/(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))+cell2mat(trialinfo((starttrial+trialjump)*2,18)))/35,'Alpha',1)
%     yline(cell2mat(trialinfo((starttrial+trialjump)*2,2)),'Color','r','Linewidth',cell2mat(trialinfo((starttrial+trialjump)*2,18))*100/((cell2mat(trialinfo((starttrial+trialjump)*2-1,18)))+cell2mat(trialinfo((starttrial+trialjump)*2,18)))/35,'Alpha',1)
%     title(title_namelabel)
%     ylim([1 nChn])
%     hcb=colorbar;
%     hcb.Title.String = "Average spikes per trial";
%     hcb.Title.Rotation=270;
%     hcb.Title.Position= [40,130];
%     hcb.Title.FontSize=11;
%     caxis([0, 2.5]);%max(z,[],'all')])
% 
% %     for loopnumcount=1:maxtid/cond
% %         hold on
% %         xline((cond*loopnumcount)+1);
% %     end
%  
%   
%     
%     %%
% end