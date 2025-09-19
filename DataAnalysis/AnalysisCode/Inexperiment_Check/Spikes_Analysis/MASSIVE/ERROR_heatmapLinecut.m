function ERROR_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,depthdriven)
%Plots true data and prediction, then plots error
loadVarAmp;
AMP=loadAMP;
loadAMP_all;
loadStimChn;
TrialParams=loadTrialParams;
trialinfo=loadTrialInfo;
loadNORECORDELECT;
trialinfo(1,:)=[];
filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
    originalEND=size(TrialParams,1);
else
    E_Mapnumber=loadMapNum;
    loadoriginalEND;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end

[~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
truedatastruct=TrueData_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,depthdriven);
%[predictdatastruct,ratio_all]=AdditivePrediction_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT, stderrspktrial,startpointseconds, secondstoanalyse, depthdriven);
[predictdatastruct]=ModelFit(AMPInterestSingleLinePlot, depthdriven); 
trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;%trials to jump over
 endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions 
 if isempty(endtrialelect)
     endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
 end
 if endtrialelect>(originalTrialEND)
     endtrialelect=(originalTrialEND);
 end

ampindex50=(AMP==AMP_all.*2);
ampindex50=logical(sum(ampindex50,2));
ampindex50(1,1)=true;
figure
fignum=gcf;
fignum=fignum.Number;
for chosenstimchn=1:length(CHN)
    if VarAmp==1
        check=['T75_25_' num2str(CHN(chosenstimchn))];
        error75=truedatastruct.(check)-predictdatastruct.(check);
        trueampint75=truedatastruct.(check)(:,AMPInterestSingleLinePlotINDEXDUAL);
        check=['T25_75_'  num2str(CHN(chosenstimchn))];
        error25=truedatastruct.(check)-predictdatastruct.(check);
        trueampint25=truedatastruct.(check)(:,AMPInterestSingleLinePlotINDEXDUAL);
        check=['T50_50_'  num2str(CHN(chosenstimchn))];
        error50=truedatastruct.(check)-predictdatastruct.(check);%(:,ampindex50');
        trueampint50=truedatastruct.(check)(:,AMPInterestSingleLinePlotINDEXDUAL);
        
        totalRMSE=sqrt(sum(sum(error75.^2+error50.^2+error25.^2)));
        subplot(1,3,3)
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(4+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error75, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 75/25' ])
        caxis([-50 100])
        subplot(1,3,2)
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(3+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error50, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 50/50' ])
        caxis([-50 100])
        subplot(1,3,1)
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(2+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error25, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 25/75' ])
        caxis([-50 100])       
        ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
        ax.Label.String='Sp/s';
        ax.Label.Rotation=270;
        
        
        figure %plotting estimation line plot
        ax = gca;
        ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
        hold on
        SMOOTHING=1;
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        
        rate = conv(error25(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        
        rate = conv(error50(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        
        rate = conv(error75(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        
        ylabel('Sp/s')
        xlabel('Channel number')
        xlim([1 nChn])
        xline((CHN(chosenstimchn)+NORECORDELECT(1)+1),'-.k')
        xline((CHN(chosenstimchn)),'k')
        legend('75/25', '50/50','25/75', 'Stim E1','Stim E2')
        title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' ERROR'])

        %normalise
%         trueampint25norm=trueampint25./((ratio_all(1)/(ratio_all(1)+1))*(1/4)+1/(ratio_all(1)+1)*(3/4))*(ratio_all(1)+1);
%         trueampint50norm=trueampint50./((ratio_all(2)/(ratio_all(2)+1))*(2/4)+1/(ratio_all(2)+1)*(2/4))*(ratio_all(2)+1);
%         trueampint75norm=trueampint75./((ratio_all(3)/(ratio_all(3)+1))*(3/4)+1/(ratio_all(3)+1)*(1/4))*(ratio_all(3)+1);
%         figure
%         ax = gca;
%         ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
%         hold on
%         SMOOTHING=1;
%         window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
%         
%         rate = conv(trueampint25,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
%         rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
%         plot(rate)
%         rate = conv(trueampint50,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
%         rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
%         plot(rate)
%         rate = conv(trueampint75,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
%         rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
%         plot(rate)
%         ylabel('Sp/s normalised')
%         xlabel('Channel number')
%         xlim([1 nChn])
%         xline((CHN(chosenstimchn)+NORECORDELECT(1)+1),'-.k')
%         xline((CHN(chosenstimchn)),'k')
%         legend('75/25', '50/50','25/75', 'Stim E1','Stim E2')
%         title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' Normalised'])

    else
        check=['T50_50_' num2str(CHN(chosenstimchn))];
        error50=truedatastruct.(check)-predictdatastruct.(check)(:,ampindex50');
        figure
        fignum=gcf;
        fignum=fignum.Number;
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(3+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error50, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 50/50' ])
        caxis([-50 100])
        
         
        figure %plotting estimation line plot
        ax = gca;
        ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
        hold on
        SMOOTHING=1;
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        rate = conv(error50(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        ylabel('Sp/s')
        xlabel('Channel number')
        xlim([1 nChn])
        xline((CHN(chosenstimchn)+NORECORDELECT(1)+1),'-.k')
        xline((CHN(chosenstimchn)),'k')

        legend('50/50','Stim E1','Stim E2')
        title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' ERROR'])

    end
end
% loadNORECORDELECT;
% elect=8;
% check=['T100_' num2str(CHN(chosenstimchn))];
% electrodeampint100=truedatastruct.(check)(elect,:);
% check=['T75_25_' num2str(CHN(chosenstimchn))];
% electrodeampint75=truedatastruct.(check)(elect,:);
% check=['T25_75_'  num2str(CHN(chosenstimchn))];
% electrodeampint25=truedatastruct.(check)(elect,:);
% check=['T50_50_'  num2str(CHN(chosenstimchn))];
% electrodeampint50=truedatastruct.(check)(elect,:);
% check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
% electrodeampint0=truedatastruct.(check)(elect,:);
% AMP=loadAMP;
% AMP_orig=[0 AMP(2:end)];
% AMP=[0 AMP_all(2:end)'];
% AMP_double=[0 AMP(2:end).*2];
% AMP_half=AMP_double./2;
% AMP_quater=[0 AMP(2:end).*(2/4)];
% AMP_threequaters=[0 AMP(2:end).*(2*3/4)];
% [level,pos]=intersect(AMP_double, AMP_orig);
% y=[electrodeampint100 electrodeampint75 electrodeampint50 electrodeampint25 electrodeampint0];
% x1=[AMP AMP_threequaters(pos) AMP_half(pos) AMP_quater(pos) zeros(1,length(electrodeampint0))];
% x2=[zeros(1,length(electrodeampint0)) AMP_quater(pos) AMP_half(pos) AMP_threequaters(pos) AMP];
% X=[x1' x2'];
% %T = [0 0 0;1 0 0;0 1 0;1 1 0]; for multiplicative
% modelfun = @(b,X) b(1)./(1+exp(b(2)+b(3)*X(:,1)+b(4)*X(:,2)));
% %modelfun = @(b,X) b(1)+b(2)*exp(b(3)*X(:,1))+b(4)*exp(b(5)*X(:,2));
% beta0 = [350 1 1 1]; %works with top equation
% mdl = fitnlm(X,y',modelfun,beta0);
% %mdl = fitlm(X',y');
% beta=mdl.Coefficients.Estimate;
% 
% x1fit = min(x1):0.2:max(x1);
% x2fit = min(x2):0.2:max(x2);
% [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
% YFIT = beta(1)./(1+exp(beta(2)+beta(3)*X1FIT + beta(4)*X2FIT));
% figure
% mesh(X1FIT,X2FIT,YFIT)
% hold on
% scatter3(x1,x2,y,'filled')
% xlabel(['Channel ' num2str(CHN(chosenstimchn)) ' Current level (uA)'])
% ylabel(['Channel ' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1) ' Current level (uA)'])
% zlabel('Sp/s')
% title(['Spike rates from electrode ' num2str(elect)])
% 
% 
% y=[electrodeampint100 electrodeampint0];
% x1=[AMP zeros(1,length(electrodeampint0))];
% x2=[zeros(1,length(electrodeampint0)) AMP];
% X=[x1' x2'];
% modelfun = @(b,X) b(1)./(1+exp(b(2)+b(3)*X(:,1)+b(4)*X(:,2)));
% %modelfun = @(b,X) b(1)+b(2)*exp(b(3)*X(:,1))+b(4)*exp(b(5)*X(:,2));
% beta0 = [350 1 1 1]; %works with top equation
% mdl = fitnlm(X,y',modelfun,beta0);
% 
% %mdl = fitlm(X',y');linear
% beta=mdl.Coefficients.Estimate;
% 
% x1fit = min(x1):0.2:max(x1);
% x2fit = min(x2):0.2:max(x2);
% [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
% %YFIT = beta(1) + beta(2)*X1FIT + beta(3)*X2FIT; %linear fit
% YFIT = beta(1)./(1+exp(beta(2)+beta(3)*X1FIT + beta(4)*X2FIT));
% figure
% 
% s=mesh(X1FIT,X2FIT,YFIT);
% s.EdgeColor='r';
% hold on
% scatter3(x1,x2,y,'filled')
% xlabel(['Channel ' num2str(CHN(chosenstimchn)) ' Current level (uA)'])
% ylabel(['Channel ' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1) ' Current level (uA)'])
% zlabel('Sp/s')
% title(['Single stim spike rates from electrode ' num2str(elect)])



end

