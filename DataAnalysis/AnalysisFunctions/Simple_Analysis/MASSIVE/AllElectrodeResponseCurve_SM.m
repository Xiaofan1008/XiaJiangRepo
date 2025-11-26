function [PlatChnArray, AmpPlatVal]=AllElectrodeResponseCurve_SM(AMPInterestSingleLinePlotINDEXDUAL,avgnospT,stderrspktrial,varargin)
% Response curve at stimulating electrode

%INPUT - the trial information array, the average number of spikes per
%trial array, and the standard error per trial array. An array of ones and
%zeros denoting whether a channel had a significant increase in spiking
AMP=loadAMP;
filepath = pwd;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir([filepath, filesep 'info.rhs']);
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
if nargin==5
   h=cell2mat(varargin(2));
   printFigures=cell2mat(varargin(1));
elseif nargin==4
    h=ones(nChn,1);
    printFigures=cell2mat(varargin(1));
else 
    h=ones(nChn,1);
    printFigures=1;
end
E_MAP = Depth(E_Mapnumber);
loadStimChn;
loadAMP_all;
AMP=AMP_all';
trialinfo=loadTrialInfo;
trialinfo(1,:)=[];
counter=0;

for count=1:length(stimChn) %loop to find the correct trials
    desiredchanneltrial=find(cell2mat(trialinfo(:,2))==stimChn(count)); %finds trials with desired initial electrode
    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
    Desired_trial=cell2mat(trialinfo(desiredchannel__singleampmath,1));
    counter=counter+1;
    if length(Desired_trial)>length(AMP)
        loopnum=length(Desired_trial)/length(AMP);
        temp_resultarray=zeros(loopnum,length(AMP));
        for reaarange=1:loopnum
            temp_resultarray(reaarange,1:length(AMP))=Desired_trial(1:length(AMP));
            Desired_trial(1:length(AMP))=[];
        end
        result_trial_array(counter:counter+1,1:length(AMP))=temp_resultarray;
        counter=counter+1;
    else
        for Amploop=1:length(AMP)
            try
                result_trial_array(count,Amploop)=(desiredchannel__singleampmath(cell2mat(trialinfo(desiredchannel__singleampmath,18))==AMP(Amploop))+1)/2;%array of mathcning trial number
            catch
                result_trial_array(Amploop)=0;
            end
        end
    end

end
PlatChnArray=zeros(nChn,1);
AmpPlatVal=zeros(nChn,1);
if printFigures==1
    figure
end
if sum(h)==nChn
    hold on
    result_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    error_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    for Chncount=1:nChn%loop for plotting
        if h(Chncount)==1
            for i=1:size(result_trial_array,1)
                for j=1:size(result_trial_array,2)
                    result_array(i,j)=avgnospT((Chncount),result_trial_array(i,j));
                    error_array(i,j)=stderrspktrial((Chncount),result_trial_array(i,j));
                end
            end
            mean_result_array=mean(result_array,1);
            mean_error_array=mean(error_array,1);
            
            if printFigures==1
            hold on
            errorbar([0 AMP(2:end)],mean_result_array,mean_error_array)
            title('Significant Channel Comparison')
            ylabel('Average number of spikes')
            xlabel('Current uA')
            end
            for i=2:length(mean_result_array)
                if mean_result_array(i)>(mean_result_array(i-1)+mean_result_array(i-1)*0.05)
                    PlatChnArray(Chncount)=AMP(i);
                    AmpPlatVal(Chncount)=mean_result_array(i);%%for plateau
                    %AmpPlatVal(Chncount)=mean_resultopfour;
                end
            end
        end
    end
    if printFigures==1
    legend_names=1:nChn;
    legend(strsplit(num2str(legend_names)))
    end
else
    if printFigures==1
    ax1=subplot(1,2,1);
    hold on
    end
    result_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    error_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    for Chncount=1:nChn%loop for plotting
        if h(Chncount)==1
            for i=1:size(result_trial_array,1)
                for j=1:size(result_trial_array,2)
                    result_array(i,j)=avgnospT((Chncount),result_trial_array(i,j));
                    error_array(i,j)=stderrspktrial((Chncount),result_trial_array(i,j));
                end
            end
            mean_result_array=mean(result_array,1);
            mean_error_array=mean(error_array,1);
            if printFigures==1
            hold on
            errorbar([0 AMP(2:end)],mean_result_array,mean_error_array)
            title('Significant Channel Comparison')
            ylabel('Average number of spikes')
            xlabel('Current uA')
            end
        end
    end
    if printFigures==1
    legend_names=1:nChn;
    legend_names=legend_names(logical(h));
    legend(strsplit(num2str(legend_names)))
    
    ax2=subplot(1,2,2);
    hold on
    end
    result_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    error_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    for Chncount=1:nChn%loop for plotting
        if h(Chncount)==0
            for i=1:size(result_trial_array,1)
                for j=1:size(result_trial_array,2)
                    result_array(i,j)=avgnospT((Chncount),result_trial_array(i,j));
                    error_array(i,j)=stderrspktrial((Chncount),result_trial_array(i,j));
                end
            end
            mean_result_array=mean(result_array,1);
            mean_error_array=mean(error_array,1);
            if printFigures==1
            hold on
            errorbar([0 AMP(2:end)],mean_result_array,mean_error_array)
            title('Not Significant Channel Comparison')
            ylabel('Average number of spikes')
            xlabel('Current uA')
            end
        end
    end
    h=abs(h-1);
    if printFigures==1
        linkaxes([ax1,ax2],'xy');

        legend_names=1:nChn;
        legend_names=legend_names(logical(h));
        legend(strsplit(num2str(legend_names)))
    end
end
end

