function Combine_trials(Subdir1, Subdir2, varargin)
%%need to combine sp so that we can see V1 vs V2 activation. also run the
%%visual stimuli

SubDir={Subdir1;Subdir2};
for loopdir=1:nargin-2
    SubDir{loopdir+2}=varargin{loopdir};
end
skip=1;
if skip~=1
concat_Trials=[];

for Trialdir=1:size(SubDir,1)
    cd(SubDir{Trialdir})
    load('IDstruct.mat','IDstruct','baslinespikestruct')

    if Trialdir==1
        concat_Trials=IDstruct;
        concat_Basline=baslinespikestruct;
    else
        for trial=1:length(fieldnames(IDstruct))
            tid=['T' num2str(trial)];
            concat_Trials.(tid)=[concat_Trials.(tid),IDstruct.(tid)];
            concat_Basline.(tid)=[concat_Basline.(tid),baslinespikestruct.(tid)];
        end
    end
end
save('Combined_trials.mat',"concat_Basline","concat_Trials")
[avgnospT,stderrspktrial,trialinfo]=AverageTrialResponse_SM(concat_Trials,concat_Basline);
end

%%
%% sigmoid

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
            if any(trialinfo(trialloop:trialloop+numsim_elect-1,17)~=-1 & trialinfo(trialloop:trialloop+numsim_elect-1,2)==chn(chnloop))
                cur=trialinfo(trialloop:trialloop+numsim_elect-1,17);
                cur=cur(cur~=-1);
                IDs{chnloop}(count,1:2)=[cur (trialloop+numsim_elect-1)./numsim_elect];
            end
        end
        %allsingletrialIDs(count)=trialloop;%index   trialinfo(trialloop,1);
    end
end


for stimchn=1:length(chn)
    IDs{stimchn}=sortrows(IDs{stimchn});
    IDs{stimchn}=IDs{stimchn}(any(IDs{stimchn},2),:);
    IDs{stimchn}=IDs{stimchn}((IDs{stimchn}(:,1)~=-1),:);
end
if skip~=1
    for stimchn=1:length(chn)
        figure
        
        for recordchn=1:size(avgnospT,1)
            subplot(8,16,recordchn)
            hold on
            plot(IDs{stimchn}(:,1),avgnospT(recordchn,IDs{stimchn}(:,2))*1000/6);
            xlabel('Current (\muA)')
            ylabel('Sp/s')
        end
        
    end
end
%% raster single
IDchosen=IDs{5}(:,2);
recordchn=[33,36];
order=Depth(0);
[pos,~]=find(order==recordchn);



for chn=1:length(recordchn)
    chnnum=['C' num2str(chn)];
    for ID=1:length(IDchosen)
        xdata = [];
        ydata = [];
        for Trialdir=1:size(SubDir,1)
            cd(SubDir{Trialdir})
            sp=loadSpikes;
            Spike_array=sp{pos(chn)};
            trig = loadTrig(0);
            trialinfo=loadTrialInfo;
            trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
            numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
            TP = loadTrialParams;
            tID = find(cell2mat(TP(numsim_elect:numsim_elect:end,2)) == IDchosen(ID));
            ID_trial=['T' num2str(IDchosen(ID))];
            theseTrig = trig(tID)./30;
            nT=length(theseTrig);
            %% Set up the raster data structure
            BIN = [-98 98];
            MAX = 400;
            
            for tr = 1:nT
                theseSp = (Spike_array(Spike_array(:,1) > theseTrig(tr)+BIN(1) & Spike_array(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
                [chninfo,ampinfo]=read_Intan_RHS2000_file;
                FS=ampinfo.amplifier_sample_rate;
                nChn=length(chninfo);
                flag=checkTrials(trig(tID(tr)),abs(BIN(1))/1000,BIN(2)/1000,nChn,FS,pos(chn), theseSp);
                if flag==1
                    continue
                end
                
                
                for i = 1:length(theseSp)
                    xdata = [xdata, (theseSp(i))];
                    ydata = [ydata, tr*(MAX/nT)];
                end
            end
            XDATA_combine.(chnnum).(ID_trial)=xdata;
            YDATA_combine.(chnnum).(ID_trial)=ydata;
        end
        figure; hold on; axM = gca;
        SMOOTHING = 2; MAX = 400;
        yScale = MAX./nT;
        for i = 1:size(xdata,2)
            line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3 ydata(i)+yScale/2.3],'Color','k','LineWidth',2);
        end
        line([0 0],[0 MAX],'Color','b');
        
        % Add the convolved spikerate
        Z = hist(xdata,BIN(1):2:BIN(2)); %#ok<HIST>
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        rate = (1000/nT)*conv(Z,window);
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(axM,BIN(1):2:BIN(2),rate,'b','LineWidth',2);
        %ylim(axM,[0 MAX]);
        % Tidy Up
        xlabel(axM,'Time (ms)');
        ylabel(axM,'Firing Rate (Sp/s)');
        %xticks(axM,150:50:300);
        %xticklabels(axM,-50:50:100);
        beautifyPlot(30,axM);
    end
end

    
        
%% raster combine all
skip=1;
if skip~=1
            xdata = [];
            ydata = [];
for Trialdir=1:size(SubDir,1)
    cd(SubDir{Trialdir})
    [chninfo,~]=read_Intan_RHS2000_file;
    sp=loadSpikes;
    nChn=length(chninfo);
    trig = loadTrig(0);
    trialinfo=loadTrialInfo;
    trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
    numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
    TP = loadTrialParams;
    for chn=1:nChn
        chnnum=['C' num2str(chn)];
        Spike_array=sp{chn};
        for ID=1:max(cell2mat(TP(:,2)))
            tID = find(cell2mat(TP(numsim_elect:numsim_elect:end,2)) == ID);
            ID_trial=['T' num2str(ID)];
            theseTrig = trig(tID)./30;
            nT=length(theseTrig);
            %% Set up the raster data structure
            BIN = [-98 98];
           MAX = 400;

            for tr = 1:nT
                theseSp = (Spike_array(Spike_array(:,1) > theseTrig(tr)+BIN(1) & Spike_array(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
                for i = 1:length(theseSp)
                    xdata = [xdata, (theseSp(i))]; 
                    ydata = [ydata, tr*(MAX/nT)];
                end
            end
            XDATA_combine.(chnnum).(ID_trial)=xdata;
            YDATA_combine.(chnnum).(ID_trial)=ydata;
        end
    end
end
save('spikeinfo.mat',"YDATA_combine","XDATA_combine")
end
end