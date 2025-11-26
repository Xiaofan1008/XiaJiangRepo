function StimElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial)
% Response curve at stimulating electrode

%INPUT - the trial information array, the average number of spikes per
%trial array, and the standard error per trial array


loadStimChn;
AMP=loadAMP;
ax=zeros(1,length(stimChn));
figure
hold on

for count=1:length(stimChn)
    
    desiredchanneltrial=find(cell2mat(trialinfo(:,2))==stimChn(count)); %finds trials with desired initial electrode
    if desiredchanneltrial(end)==size(trialinfo(:,2),1)
        desiredchanneltrial(end)=[];
    end
    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds trials with matching recording electrode spacing
    Desired_trial=cell2mat(trialinfo(desiredchannel__singleampmath,1));
    result_array=avgnospT((stimChn(count)),Desired_trial);
    error_array=stderrspktrial((stimChn(count)),Desired_trial);
    if length(result_array)>length(AMP)
        loopnum=length(result_array)/length(AMP);
        temp_resultarray=zeros(1,length(AMP));
        temp_error_array=zeros(1,length(AMP));
        for avgarray=1:loopnum
            temp_resultarray=temp_resultarray+result_array(1+length(AMP)*(avgarray-1):length(AMP)+length(AMP)*(avgarray-1));
            temp_error_array=temp_error_array+error_array(1+length(AMP)*(avgarray-1):length(AMP)+length(AMP)*(avgarray-1));
        end
        temp_resultarray=temp_resultarray./loopnum;
        result_array=temp_resultarray;
        temp_error_array=temp_error_array./loopnum;
        error_array=temp_error_array;
    end
    
    ax(count) = subplot(1,length(stimChn),count);
    hold on
    errorbar(AMP,result_array,error_array)
    title(['Channel ' num2str(stimChn(count))])
    ylabel('Average number of spikes')
    xlabel('Current uA')
    if count~=1
        linkaxes([ax(count-1) ax(count)],'xy') %create same axes
    end
end

end

