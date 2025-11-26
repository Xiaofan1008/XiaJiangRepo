function [electfit,MSEtotal]=SigmoidfunctionNoDual(SubDir_Path)
folder = fileparts(which('SigmoidfunctionNoDual.m')); % Determines filepath to folder containing your .m file.
addpath(genpath(folder)); % Add that folder plus all subfolders to the path.
%fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd(SubDir_Path)%change working directory to where data is stored - currently manually input
%%
load('truedatastruct.mat','truedatastruct')
PlotNow=0;
includeResult=0;
loadNORECORDELECT;
AMP=loadAMP;
loadAMP_all;
loadCHN;
chosenstimchn=1;
trialinfo=loadTrialInfo;
AMP_orig=[0 AMP(2:end)];
AMP=[0 AMP_all(2:end)'];
AMP_double=[0 AMP(2:end).*2];
AMP_half=AMP_double./2;
AMP_quater=[0 AMP(2:end).*(2/4)];
AMP_threequaters=[0 AMP(2:end).*(2*3/4)];
MSEtotal=0;
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
if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
    laminar=1;
else
    laminar=0;
end
electbeta=zeros(nChn,4);
electfit=zeros(nChn,6);
for elect=1:nChn
    check=['T100_' num2str(CHN(chosenstimchn))];
    electrodeampint100=truedatastruct.(check)(elect,1:end);
    check=['T75_25_' num2str(CHN(chosenstimchn))];
    electrodeampint75=truedatastruct.(check)(elect,1:end);
    check=['T25_75_'  num2str(CHN(chosenstimchn))];
    electrodeampint25=truedatastruct.(check)(elect,1:end);
    check=['T50_50_'  num2str(CHN(chosenstimchn))];
    electrodeampint50=truedatastruct.(check)(elect,1:end);
    if laminar==1
        check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
    else
        check=['T100_' num2str(NORECORDELECT(1))];
    end
    electrodeampint0=truedatastruct.(check)(elect,1:end);
    
    k1=max(electrodeampint100);
    k2=max(electrodeampint0);
    [~,pos]=intersect(AMP_double, AMP_orig);
    AMP_threequatersp=AMP_threequaters(pos);
    AMP_quaterp=AMP_quater(pos);
    AMP_halfp=AMP_half(pos);
    
    if includeResult==1
        y=[electrodeampint100 electrodeampint75 electrodeampint50 electrodeampint25 electrodeampint0];
        x1=[AMP(1:2:end) AMP_threequatersp(1:2:end) AMP_halfp(1:2:end) AMP_quaterp(1:2:end) zeros(1,length(electrodeampint0))];
        x2=[zeros(1,length(electrodeampint0)) AMP_quaterp(1:2:end) AMP_halfp(1:2:end) AMP_threequatersp(1:2:end) AMP(1:2:end)];
        X=[x1' x2'];
    else
        y=[electrodeampint100 electrodeampint0];
        x1=[AMP zeros(1,length(electrodeampint0))];
        x2=[zeros(1,length(electrodeampint0)) AMP];
        X=[x1' x2'];
    end
    modelfun = @(b,X) b(1)./(1+exp(b(2)+b(3).*X(:,1)+b(4).*X(:,2))); %SIGMOIDAL  b1 is asymptote, b2=x50/slope b3=-1/slope b2=-1/slope  This is boltzmann sigmoid
    beta0 = [350 1 1 1]; %SIGMOIDAL
    try
        mdl = fitnlm(X,y',modelfun,beta0);
    catch
        continue
    end
    fprintf([num2str(elect) '\n']);
    beta=mdl.Coefficients.Estimate;
    electbeta(elect,:)=beta;
    
    
    if PlotNow==1
        x1fit = min(x1):0.2:max(x1);
        x2fit = min(x2):0.2:max(x2);
        [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
        YFIT = beta(1)./(1+exp(beta(2)+beta(3).*X1FIT + beta(4).*X2FIT));
        
%         if k1>k2
%             k=k1;
%         else
%             k=k2;
%         end
%         YFIT(YFIT>k)=k;
%         YFIT(YFIT<0)=0;
        figure
        mesh(X1FIT,X2FIT,YFIT);
        hold on
        scatter3(x1,x2,y,'filled')
        xlabel(['Channel ' num2str(CHN(chosenstimchn)) ' Current level (uA)'])
        ylabel(['Channel ' num2str(NORECORDELECT(1)) ' Current level (uA)'])
        zlabel('Sp/s')
        title(['Spike rates from electrode ' num2str(elect)])
    end
    
    %%not included points in the model
    
    check=['T75_25_' num2str(CHN(chosenstimchn))];
    electrodeampint75=truedatastruct.(check)(elect,1:end);
    check=['T25_75_'  num2str(CHN(chosenstimchn))];
    electrodeampint25=truedatastruct.(check)(elect,1:end);
    check=['T50_50_'  num2str(CHN(chosenstimchn))];
    electrodeampint50=truedatastruct.(check)(elect,1:end);

    y=[electrodeampint75 electrodeampint50 electrodeampint25 ];
    x1=[AMP_threequatersp AMP_halfp AMP_quaterp];
    x2=[AMP_quaterp AMP_halfp AMP_threequatersp];
    
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
    if (ssd>0.75) && (beta(1)<10000) && ((-beta(2)/beta(3))>0) && ((-beta(2)/beta(4))>0)
        electfit(elect,1)=1;
        fprintf('Good fit \n')
    end
    electfit(elect,2)=ssd;
    electfit(elect,3:6)=beta;
    MSEtotal=MSEtotal+errory;%summed error if positive then measured values are greater than that expected from single electrode stim
end

save('Boltzsig.mat','electfit')
end

