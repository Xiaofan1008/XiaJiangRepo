function [electfit]=BoltzmannSigmoidFit
warning('on','MATLAB:nearlySingularMatrix')
warning('on','stats:nlinfit:IllConditionedJacobian')
warning('on','stats:nlinfit:IterationLimitExceeded')

PlotNow=0;
includeResult=1;
 trialinfo=loadTrialInfo;
 trialinfo(1,:)=[];
 loadCHN;
 loadStimChn;
chosenstimchn=1;
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
loadNORECORDELECT;
if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
    laminar=1;
else
    laminar=0;
end





%perform model fit

AMP=loadAMP;
loadAMP_all;
loadCHN;
trialinfo=loadTrialInfo;
AMP_orig=[0 AMP(2:end)];
AMP=[0 AMP_all(2:end)'];
AMP_double=[0 AMP(2:end).*2];
AMP_half=AMP_double./2;
AMP_quater=[0 AMP(2:end).*(2/4)];
AMP_threequaters=[0 AMP(2:end).*(2*3/4)];
load('truedatastruct.mat','truedatastruct')
MSEtotal=0;
electbeta=zeros(nChn,4);
electfit=zeros(nChn,6);
step=1;
for elect=1:nChn
    check=['T100_' num2str(CHN(chosenstimchn))];
    electrodeampint100=truedatastruct.(check)(elect,1:step:end);
    check=['T75_25_' num2str(CHN(chosenstimchn))];
    electrodeampint75=truedatastruct.(check)(elect,1:step:end);
    check=['T25_75_'  num2str(CHN(chosenstimchn))];
    electrodeampint25=truedatastruct.(check)(elect,1:step:end);
    check=['T50_50_'  num2str(CHN(chosenstimchn))];
    electrodeampint50=truedatastruct.(check)(elect,1:step:end);
    if laminar==1
        check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
    else
        check=['T100_' num2str(NORECORDELECT(1))];
    end
    electrodeampint0=truedatastruct.(check)(elect,1:step:end);
    
    k1=max(electrodeampint100);
    k2=max(electrodeampint0);
    [~,pos]=intersect(AMP_double, AMP_orig);
    AMP_threequatersp=AMP_threequaters(pos);
    AMP_quaterp=AMP_quater(pos);
    AMP_halfp=AMP_half(pos);
    
    if includeResult==1
        y=[electrodeampint100 electrodeampint75 electrodeampint50 electrodeampint25 electrodeampint0];
        x1=[AMP(1:step:end) AMP_threequatersp(1:step:end) AMP_halfp(1:step:end) AMP_quaterp(1:step:end) zeros(1,length(electrodeampint0))];
        x2=[zeros(1,length(electrodeampint0)) AMP_quaterp(1:step:end) AMP_halfp(1:step:end) AMP_threequatersp(1:step:end) AMP(1:step:end)];
        X=[x1' x2'];
    else
        y=[electrodeampint100 electrodeampint0];
        x1=[AMP zeros(1,length(electrodeampint0))];
        x2=[zeros(1,length(electrodeampint0)) AMP];
        X=[x1' x2'];
    end
    %modelfun = @(b,X) b(1)+b(2)*X(:,1)+b(3)*X(:,2)+b(4)*X(:,1).^2+b(5)*X(:,2).^2;%quadratic
    modelfun = @(b,X) b(1)./(1+exp(b(2)+b(3).*X(:,1)+b(4).*X(:,2))); %SIGMOIDAL  b1 is asymptote, b2=x50/slope b3=-1/slope b2=-1/slope  This is boltzmann sigmoid
    beta0 = [350 1 1 1]; %SIGMOIDAL
    %beta0 = [1 1 1 1 1]; %quad
    try
        mdl = fitnlm(X,y',modelfun,beta0);
    catch
        continue
    end
    fprintf([num2str(elect) '\n']);
    beta=mdl.Coefficients.Estimate;
    electbeta(elect,:)=beta;
    pval = coefTest(mdl);
    rsq = mdl.Rsquared.Adjusted;
    
    if PlotNow==1
        x1fit = min(x1):0.2:max(x1);
        x2fit = min(x2):0.2:max(x2);
        [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
        YFIT = beta(1)./(1+exp(beta(2)+beta(3).*X1FIT + beta(4).*X2FIT));
        figure
        mesh(X1FIT,X2FIT,YFIT);
        hold on
        scatter3(x1,x2,y,'filled')
        %scatter3(x1,x2,log((y./beta(1))./(1-(y./beta(1)))),'filled')
        xlabel(['Channel ' num2str(CHN(chosenstimchn)) ' Current level (uA)'])
        %ylabel(['Channel ' num2str(CHN(chosenstimchn+1)) ' Current level (uA)'])
        ylabel(['Channel ' num2str(stimChn(chosenstimchn+1)) ' Current level (uA)'])
        zlabel('Sp/s')
        title(['Spike rates from electrode ' num2str(elect)])
        
        
        figure
        hold on
        x1fit = AMP_orig;
        x2fit = AMP_orig;
        YFIT = beta(1)./(1+exp(beta(2)+beta(3).*x1fit + beta(4).*x2fit));
        YFIT1 = beta(1)./(1+exp(beta(2)+beta(3).*x1fit));
        YFIT2 = beta(1)./(1+exp(beta(2) + beta(4).*x2fit));
        plot(AMP_orig,YFIT, 'r')
        plot(AMP_orig,YFIT1,'k-.')
        plot(AMP_orig,YFIT2,'k--')
        plot(AMP_orig,YFIT2+YFIT1,'k')
                xlabel('current (uA)')
        ylabel('Sp/s')
        title(['Model spike rates from electrode ' num2str(elect)])
        %plot(AMP(1:step:end),electrodeampint100,'k:')
        plot(AMP_halfp(1:step:end) ,electrodeampint50, 'r:')
        plot(AMP_halfp(1:step:end),electrodeampint0(pos)+electrodeampint100(pos),'k:')
        
        %scatter(AMP(1:step:end),electrodeampint100,10,'k')
        scatter(AMP_halfp(1:step:end) ,electrodeampint50,10, 'r')
        scatter(AMP_halfp(1:step:end),electrodeampint0(pos)+electrodeampint100(pos),10,'k')

        legend('Dual model', ['E' num2str(stimChn(chosenstimchn)) ' model'],['E' num2str(stimChn(chosenstimchn+1)) 'model'],['E' num2str(stimChn(chosenstimchn)) '+E' num2str(stimChn(chosenstimchn+1)) ' model'],'Dual real',['E' num2str(stimChn(chosenstimchn)) '+E' num2str(stimChn(chosenstimchn+1)) ' real'])


    end
    
    %%not included points in the model
    check=['T100_' num2str(CHN(chosenstimchn))];
    electrodeampint100=truedatastruct.(check)(elect,step:step:end);
    check=['T75_25_' num2str(CHN(chosenstimchn))];
    electrodeampint75=truedatastruct.(check)(elect,step:step:end);
    check=['T25_75_'  num2str(CHN(chosenstimchn))];
    electrodeampint25=truedatastruct.(check)(elect,step:step:end);
    check=['T50_50_'  num2str(CHN(chosenstimchn))];
    electrodeampint50=truedatastruct.(check)(elect,step:step:end);
    if laminar==1
        check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
    else
        check=['T100_' num2str(NORECORDELECT(1))];
    end
    electrodeampint0=truedatastruct.(check)(elect,step:step:end);
    y=[electrodeampint100 electrodeampint75 electrodeampint50 electrodeampint25 electrodeampint0];
    x1=[AMP(step:step:end) AMP_threequatersp(step:step:end) AMP_halfp(step:step:end) AMP_quaterp(step:step:end) zeros(1,length(electrodeampint0))];
    x2=[zeros(1,length(electrodeampint0)) AMP_quaterp(step:step:end) AMP_halfp(step:step:end) AMP_threequatersp(step:step:end) AMP(step:step:end)];
    
    YFIT = beta(1)./(1+exp(beta(2)+beta(3).*x1 + beta(4).*x2));
%     if k1>k2
%         k=k1;
%     else
%         k=k2;
%     end
%     YFIT(YFIT>k)=k;
%     YFIT(YFIT<0)=0;
    errory=y-YFIT;
    ymu=y-mean(y);
    ssd=1-sum(errory.^2)/sum(ymu.^2);

      [warnmsg, msgid] = lastwarn;

    
    
    
    if  strcmp(msgid,'')  && (rsq>0.5) && (beta(1)<10000) && ((-beta(2)/beta(3))>0) && ((-beta(2)/beta(4))>0)  && ((beta(3)/beta(4))>(1/5) && (beta(3)/beta(4))<5)% (rsq>0.7)  &&&& (pval<0.05) %&& ((max(electrodeampint0)/max(electrodeampint100))<5 && (max(electrodeampint0)/max(electrodeampint100))>0.2) %&& ((beta(3)/beta(4))>0.5 && (beta(3)/beta(4))<2) % ((max(electrodeampint0)/max(electrodeampint100))<3 && (max(electrodeampint0)/max(electrodeampint100))>0.333)
        electfit(elect,1)=1;
        fprintf('Good fit \n')
%                 figure
%         hold on
%         x1fit = AMP_orig;
%         x2fit = AMP_orig;
%         YFIT = beta(1)./(1+exp(beta(2)+beta(3).*x1fit + beta(4).*x2fit));
%         YFIT1 = beta(1)./(1+exp(beta(2)+beta(3).*x1fit));
%         YFIT2 = beta(1)./(1+exp(beta(2) + beta(4).*x2fit));
%         plot(AMP_orig,YFIT, 'r')
%         plot(AMP_orig,YFIT1,'k-.')
%         plot(AMP_orig,YFIT2,'k--')
%         plot(AMP_orig,YFIT2+YFIT1,'k')
%                 xlabel('current (uA)')
%         ylabel('Sp/s')
%         title(['Model spike rates from electrode ' num2str(elect)])
%         %plot(AMP(1:step:end),electrodeampint100,'k:')
%         plot(AMP_halfp(1:step:end) ,electrodeampint50, 'r:')
%         plot(AMP_halfp(1:step:end),electrodeampint0(pos)+electrodeampint100(pos),'k:')
%         
%         %scatter(AMP(1:step:end),electrodeampint100,10,'k')
%         scatter(AMP_halfp(1:step:end) ,electrodeampint50,10, 'r')
%         scatter(AMP_halfp(1:step:end),electrodeampint0(pos)+electrodeampint100(pos),10,'k')
% 
%         legend('Dual model', ['E' num2str(stimChn(chosenstimchn)) ' model'],['E' num2str(stimChn(chosenstimchn+1)) 'model'],['E' num2str(stimChn(chosenstimchn)) '+E' num2str(stimChn(chosenstimchn+1)) ' model'],'Dual real',['E' num2str(stimChn(chosenstimchn)) '+E' num2str(stimChn(chosenstimchn+1)) ' real'])

    end
    lastwarn(''); % reset warning state
    electfit(elect,2)=rsq;
    electfit(elect,3:6)=beta;
    MSEtotal=MSEtotal+ssd;%summed R^2
end

save('Boltzsig.mat','electfit')
end

