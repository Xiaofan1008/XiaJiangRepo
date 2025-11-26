function outputMaxHist=hypothplot_sab(chn,outputMaxHist)
dbstop if error
%% Load in data from current directory
loadVarAmp;
loadNREP;
nT=n_REP_true;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir('info.rhs');
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
d = Depth(E_Mapnumber); chn = d(chn);
sp = loadSpikes(chn);
sp = denoiseSpikes(sp,chn);
TP = loadTrialParams;
trialinfo=TrialInfo_sab;
trialinfo(1,:)=[];
cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
totalnumbercond=length(trialinfo(:,2))/(2*cond); %gives the total number of condition changes i.e. 100/0 75/25 but with same electrode are classed as one condition change
loadNORECORDELECT;
loadStimChn;
loadVarAmp;
loadSingleDual;
checkmaxsize=0;
check=0;
SMOOTHING = 2;
for count=1:length(CHN)
    for recordcount=1:length(NORECORDELECT)
        clear relatedtrials
        desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)); %finds trials with desired initial electrode
        %single stim
        if (desiredchanneltrial(end)==length(trialinfo)) %need to remove the last trial in case there is a stim electrode pair resulting in 1 as the last position
            desiredchanneltrial(end)=[];
        end
        desiredchanneltrial_plus=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(CHN(count)+NORECORDELECT(recordcount)+1)),1);% finds triials with matching recording electrode spacing
        if (VarAmp==1)||(singanddual==1)
            modified_desiredchanneltrial_plus= desiredchanneltrial_plus;
            for condition=1:totalnumbercond
                if any(modified_desiredchanneltrial_plus<(cond*condition*2)) && (check>0)
                    check=check+1;
                    relatedtrials(:,check)=modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)));
                elseif any(modified_desiredchanneltrial_plus<(cond*condition*2))
                    relatedtrials=modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)));
                    check=1;
                end
                modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)))=[];
                
                if  isempty(modified_desiredchanneltrial_plus)
                    break;
                end
            end
            if exist('relatedtrials', 'var')
                relatedtrials( :, all(~relatedtrials,1) ) = [];
                desiredchannel__singleamp=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                if recordcount<2
                    desiredchannel__singleamp(desiredchannel__singleamp<checkmaxsize)=[];
                end
                desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)+NORECORDELECT(recordcount)+1); %finds trials with desired initial electrode
                if desiredchanneltrial(end)==length(trialinfo)
                    desiredchanneltrial(end)=[];
                end
                desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                if VarAmp == 1
                    desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<max(relatedtrials,[],'all'));
                end
                if length(desiredchannel__singleampmath)~=size(relatedtrials,2)
                    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                    desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<(cond*2*length(relatedtrials)));
                end
                relatedtrials_all=[desiredchannel__singleampmath';relatedtrials; desiredchannel__singleamp'];
                relatedtrials_all=(relatedtrials_all+1)./2;
            end
        else
            relatedtrials_all=(desiredchanneltrial_plus+1)./2;
        end
        
        if VarAmp~=0
            checkmaxsize=max(relatedtrials,[],'all');
        end
        maxhist=zeros(size(relatedtrials_all,1),size(relatedtrials_all,2));
        %%%%%%%%%%%%%raster input
        for counter_trials=1:size(relatedtrials_all,2)
            for counter_trials_row=1:size(relatedtrials_all,1)
                ID=relatedtrials_all(counter_trials_row,counter_trials);
                tID = find((cell2mat(TP(:,2)) == ID));
                tID(1:2:end)=[];
                tID=tID./2;
                trig = loadTrig(0);
                
                trig = trig(tID);
                trig(trig==-500)=[];
                theseTrig = trig./30;
                
                BIN = [-200 200];
                xdata = [];
                
                for tr = 1:length(trig)
                    theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
                    for i = 1:length(theseSp)
                        xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
                    end
                end
                
                
                % Add the convolved spikerate
                Z = hist(xdata,0:(abs(BIN(1))+abs(BIN(2)))); %#ok<HIST> sorts into bins of 1ms
                
                window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
                rate = (1000/nT)*conv(Z,window);%1000ms with nT number of trials
                rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                
                maxhist(counter_trials_row,counter_trials)=max(rate(200:250));
            end
        end
        struct_id=['Chn_', num2str(find(d==chn)), '_StimChn_', num2str(cell2mat(trialinfo(relatedtrials_all(2,1)*2,2))), '_', num2str(cell2mat(trialinfo(relatedtrials_all(2,1)*2-1,2)))];
        outputMaxHist.(struct_id)=maxhist;
    end
end