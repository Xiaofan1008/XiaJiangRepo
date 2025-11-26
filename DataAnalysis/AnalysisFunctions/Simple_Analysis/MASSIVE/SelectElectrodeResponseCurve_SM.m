function SelectElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,electrodeChosen)
% For plotting stim electrode response curve and electrodes between stimulating
% electrodes response curves

%INPUT - the trial information array, the average number of spikes per
%trial array, and the standard error per trial array


cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition 
totalnumbercond=length(trialinfo(:,2))/(2*cond); %gives the total number of condition changes i.e. 100/0 75/25 but with same electrode are classed as one condition change
loadNORECORDELECT;
loadStimChn;
AMP=loadAMP;
loadVarAmp;
checkmaxsize=0;
ax=zeros(1,length(AMP));
check=0;
for count=1:length(CHN)
    for recordcount=1:length(NORECORDELECT)
        clear relatedtrials
        clear result_matrix
        clear errortrial_matrix
        desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)); %finds trials with desired initial electrode
        
        figure
        hold on
        %single stim
        if (desiredchanneltrial(end)==length(trialinfo)) %need to remove the last trial in case there is a stim electrode pair resulting in 1 as the last position
            desiredchanneltrial(end)=[];
        end
        desiredchanneltrial_plus=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(CHN(count)+NORECORDELECT(recordcount)+1)),1);% finds triials with matching recording electrode spacing
        if VarAmp==1
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
                desiredchannel__singleamp(desiredchannel__singleamp<checkmaxsize)=[];
                desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)+NORECORDELECT(recordcount)+1); %finds trials with desired initial electrode
                if desiredchanneltrial(end)==length(trialinfo)
                    desiredchanneltrial(end)=[];
                end
                desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<max(relatedtrials,[],'all'));
                if length(desiredchannel__singleampmath)~=size(relatedtrials,2)
                    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                    desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<(cond*2*length(relatedtrials)));
                end
                relatedtrials_all=[desiredchannel__singleampmath';relatedtrials; desiredchannel__singleamp'];
                relatedtrials_all=(relatedtrials_all+1)./2;
                
                for plotcount=1:size(relatedtrials_all,2)
                    
                    for chncount=1:NORECORDELECT(recordcount)+2
                        if chncount==electrodeChosen-CHN(count)+1
                        for trialcount=1:size(relatedtrials_all,1)
                            result_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),(trialcount))=avgnospT((CHN(count)+chncount-1),relatedtrials_all(trialcount,plotcount));
                            errortrial_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),(trialcount))=stderrspktrial((CHN(count)+chncount-1),relatedtrials_all(trialcount,plotcount));
                        end
                        ax(plotcount)=subplot(1,size(relatedtrials_all,2),plotcount);
                        hold on
                        errorbar(1:size(relatedtrials_all,1),result_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),:),errortrial_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),:))
                       end
                    end
                    ylabel('Average number of spikes')
                    xlabel('Stimulating current(ratio)')
                    title(['Electrodes ', (num2str(CHN(count)+NORECORDELECT(recordcount)+1)),' and ', (num2str(CHN(count))), ' @ ' ,num2str(AMP(plotcount)), 'uA'])
                    xticks(1:size(relatedtrials_all,1))
                    xticklabels({'100/0','75/25','50/50','25/75','0/100'})
                    %legendentry=string(num2cell(electrodeChosen));
                    %lgd=legend(legendentry);
                    %title(lgd,'Electrode number')
                    if plotcount~=1
                        linkaxes([ax(plotcount-1) ax(plotcount)],'xy') %create same axes
                    end
                end
            end
        else
            relatedtrials_all=(desiredchanneltrial_plus+1)./2;
                for plotcount=1:size(relatedtrials_all,2)
                    
                    for chncount=1:NORECORDELECT(recordcount)+2
                        for trialcount=1:size(relatedtrials_all,1)
                            result_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),(trialcount))=avgnospT((CHN(count)+chncount-1),relatedtrials_all(trialcount,plotcount));
                            errortrial_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),(trialcount))=stderrspktrial((CHN(count)+chncount-1),relatedtrials_all(trialcount,plotcount));
                        end
                        hold on
                        errorbar(AMP,result_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),:),errortrial_matrix(((plotcount-1)*size(relatedtrials_all,1)+chncount),:))
                    end
                    ylabel('Average number of spikes')
                    xlabel('Stimulating current(uA)')
                    title(['Electrodes ', (num2str(CHN(count))),' and ', (num2str(CHN(count)+NORECORDELECT(recordcount)+1))])
                    legendentry=string(num2cell(CHN(count):CHN(count)+NORECORDELECT(recordcount)+1));
                    lgd=legend(legendentry);
                    title(lgd,'Electrode number')
                end
        end
    end
    if VarAmp~=0
        checkmaxsize=max(relatedtrials,[],'all');
    end
end
end

