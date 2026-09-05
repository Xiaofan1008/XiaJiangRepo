function singleElectrodeResponseCurveAMPLITUDE_SM(trialinfo,avgnospT,stderrspktrial,Chan,starttrial,trialjump)
% Response curve at single chosen electrode based on chosen trials

%INPUT - the trial information array, the average number of spikes per
%trial array, and the standard error per trial array. Then Chosen channels
%and the stimulation channel of interest

trials=cell2mat(trialinfo(starttrial*2:trialjump*2:end,1));
figure
hold on
result_array=zeros(size(Chan,1),size(length(trials),2));
error_array=zeros(size(Chan,1),size(length(trials),2));
xticks_names=cell(size(Chan,1),size(length(trials),2));
lgend_names=cell(size(Chan,1),1);
for Chncount=1:length(Chan)%loop for plotting
    chosenchannel=Chan(Chncount);
    for trialnum=1:length(trials)
        chosentrial=trials(trialnum);
        result_array(Chncount,trialnum)=avgnospT((chosenchannel),chosentrial);
        error_array(Chncount,trialnum)=stderrspktrial((chosenchannel),chosentrial);
        if cell2mat(trialinfo(trialjump*trialnum*2,18))==-1
            xticks_names(Chncount,trialnum)={num2str(0)};
        else
            xticks_names(Chncount,trialnum)={num2str(cell2mat(trialinfo(trialjump*trialnum*2,18)))};
        end
    end
    
    hold on
    EB=errorbar(1:size(result_array,2),result_array(Chncount,:),error_array(Chncount,:));
    p=polyfit(1:size(result_array,2),result_array(Chncount,:),3);
    x1 = linspace(1,size(result_array,2));
    f1 = polyval(p,x1);
    EB.Color=[0 (Chncount-0.5)/length(Chan) 0.5];
    LI=plot(x1,f1);
    LI.Color=[0 ((Chncount)*0.75)/length(Chan) 0.8];
    lgend_names(Chncount*2-1,1)={num2str(chosenchannel)};
    lgend_names(Chncount*2,1)={['Fitted ' num2str(chosenchannel)]};
end
    xticks(1:size(result_array,2))
    xticklabels(xticks_names(1,:))
    loadStimChn;
    loadNORECORDELECT;
    title(['Stimchn w/ ' num2str(stimChn(1)) ' ' num2str(stimChn(2)) ' @ ' num2str(cell2mat(trialinfo(trials(3)*2-1,18))*100/(cell2mat(trialinfo(trials(3)*2-1,18))+cell2mat(trialinfo(trials(3)*2,18)))) '/' num2str(cell2mat(trialinfo(trials(3)*2,18))*100/((cell2mat(trialinfo(trials(3)*2-1,18)))+cell2mat(trialinfo(trials(3)*2,18))))])
    lgd=legend(lgend_names);
    title(lgd,'Electrode number')
    ylabel('Average number of spikes')
    xlabel('Current uA')
    ylim([0 3])
    xlim([1 size(result_array,2)])

end

