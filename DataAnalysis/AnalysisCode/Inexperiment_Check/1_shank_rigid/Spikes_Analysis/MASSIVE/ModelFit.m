function [predictdatastruct]=ModelFit(AMPInterestSingleLinePlot, depthdriven)


load('truedatastruct.mat','truedatastruct')
trialinfo=loadTrialInfo(0);
PlotNow=0;
includeResult=0;
loadNORECORDELECT;
AMP=loadAMP;
loadAMP_all;
loadCHN;
chosenstimchn=1;
ampinterest=AMPInterestSingleLinePlot;
AMP_orig=[0 AMP(2:end)];
AMP=[0 AMP_all(2:end)'];
AMP_double=[0 AMP(2:end).*2];
AMP_half=AMP_double./2;
AMP_quater=[0 AMP(2:end).*(2/4)];
AMP_threequaters=[0 AMP(2:end).*(2*3/4)];
[~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot/2)); %x is a throw away variable but it will te ll you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
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
electbeta=zeros(nChn,5);
Spiking50_50=zeros(nChn,length(AMP_orig));
Spiking25_75=zeros(nChn,length(AMP_orig));
Spiking75_25=zeros(nChn,length(AMP_orig));
for ampiterate=1:length(AMP_orig)
for elect=1:nChn
    AMPInterestSingleLinePlot=AMP_orig(ampiterate);
    [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot/2)); %x is a throw away variable but it will te ll you how close the value is to your chosen val
    AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
    
    check=['T100_' num2str(CHN(chosenstimchn))];
    electrodeampint100=truedatastruct.(check)(elect,:);
    check=['T75_25_' num2str(CHN(chosenstimchn))];
    electrodeampint75=truedatastruct.(check)(elect,:);
    check=['T25_75_'  num2str(CHN(chosenstimchn))];
    electrodeampint25=truedatastruct.(check)(elect,:);
    check=['T50_50_'  num2str(CHN(chosenstimchn))];
    electrodeampint50=truedatastruct.(check)(elect,:);
    check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
    %check=['T100_' num2str(NORECORDELECT(1))];
    electrodeampint0=truedatastruct.(check)(elect,:);
    [level,pos]=intersect(AMP_double, AMP_orig);
    if includeResult==1
    y=[electrodeampint100 electrodeampint75 electrodeampint50 electrodeampint25 electrodeampint0];
    x1=[AMP AMP_threequaters(pos) AMP_half(pos) AMP_quater(pos) zeros(1,length(electrodeampint0))];
    x2=[zeros(1,length(electrodeampint0)) AMP_quater(pos) AMP_half(pos) AMP_threequaters(pos) AMP];
    X=[x1' x2'];
    else 
        y=[electrodeampint100 electrodeampint0];
        x1=[AMP zeros(1,length(electrodeampint0))];
        x2=[zeros(1,length(electrodeampint0)) AMP];
        X=[x1' x2'];
    end
    
    
        %T = [0 0 0;1 0 0;0 1 0;1 1 0]; for multiplicative
    %modelfun = @(b,X) b(1)*X(:,1)+b(2)*X(:,2)+b(3); %linear
    %modelfun = @(b,X) b(1)./(1+exp(b(2)+b(3)*X(:,1)+b(4)*X(:,2))); %SIGMOIDAL
    %modelfun = @(b,X) b(1)+b(2)*X(:,1)+b(3)*X(:,2)+b(4)*X(:,1).^2+b(5)*X(:,2).^2;%quadratic

    %modelfun = @(b,X) b(1)+b(2)*exp(b(3)*X(:,1))+b(4)*exp(b(5)*X(:,2));
    %beta0 = [1 1 1]; %linear
    %beta0 = [350 1 1 1]; %SIGMOIDAL
    %beta0 = [1 1 1 1 1]; %quadratic
    %mdl = fitnlm(X,y',modelfun,beta0);
    %mdl = fitlm(X',y');
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
%     Spiking50_50(elect,ampiterate)=beta(1)./(1+exp(beta(2)+beta(3)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL) + beta(4)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)));
%     Spiking75_25(elect,ampiterate)=beta(1)./(1+exp(beta(2)+beta(3)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL) + beta(4)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)));
%     Spiking25_75(elect,ampiterate)=beta(1)./(1+exp(beta(2)+beta(3)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL) + beta(4)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)));
    Spiking50_50(elect,ampiterate)=beta(1)+beta(2)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)+(beta(3)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL))+beta(4)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL).^2+beta(5)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL).^2;
    Spiking75_25(elect,ampiterate)=beta(1)+beta(2)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)+(beta(3)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL))+beta(4)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL).^2+beta(5)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL).^2;
    Spiking25_75(elect,ampiterate)=beta(1)+beta(2)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)+(beta(3)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL))+beta(4)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL).^2+beta(5)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL).^2;
%     Spiking50_50(elect,ampiterate)=beta(1)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)+beta(2)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)+beta(3);
%     Spiking75_25(elect,ampiterate)=beta(1)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)+beta(2)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)+beta(3);
%     Spiking25_75(elect,ampiterate)=beta(1)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)+beta(2)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)+beta(3);


    if PlotNow==1
        x1fit = min(x1):0.2:max(x1);
        x2fit = min(x2):0.2:max(x2);
        [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
        %YFIT = beta(1)./(1+exp(beta(2)+beta(3)*X1FIT + beta(4)*X2FIT));
        %YFIT=beta(1)+beta(2)*X1FIT+(beta(3)*X2FIT)+beta(4)*X1FIT.^2+beta(5)*X2FIT.^2;
        YFIT=beta(1)*X1FIT+beta(2)*X2FIT+beta(3);
        figure
        m=mesh(X1FIT,X2FIT,YFIT);
        %s.EdgeColor='r';
        hold on
        scatter3(x1,x2,y,'filled')
        xlabel(['Channel ' num2str(CHN(chosenstimchn)) ' Current level (uA)'])
        ylabel(['Channel ' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1) ' Current level (uA)'])
        zlabel('Sp/s')
        title(['Spike rates from electrode ' num2str(elect)])
    end
    MSEtotal=MSEtotal+mdl.MSE;
end

end

%HEATMAPS
%25/75
        figure
        fignum=gcf;
        fignum=fignum.Number;
        subplot(1,3,1)
        DepthChangeingSpiking_SM(Spiking75_25, 1:size(Spiking75_25,2), fignum,depthdriven,size(Spiking75_25,2),AMP_orig);%plots heatmap of estimated plots
        title(['ESTIMATION Stimchn: ' num2str(CHN(chosenstimchn)) ' ' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1) ' @ 25/75' ])
        yline((depthdriven-50*(CHN(chosenstimchn)+NORECORDELECT(1)+1-1)),'Color','r','Linewidth',75/25+25/35,'Alpha',1)
        %50/50
        subplot(1,3,2)
        DepthChangeingSpiking_SM(Spiking50_50, 1:size(Spiking50_50,2), fignum,depthdriven,size(Spiking50_50,2),AMP_orig);%plots heatmap of estimated plots
        title(['ESTIMATION Stimchn: ' num2str(CHN(chosenstimchn)) ' ' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1) ' @ 50/50' ])

        %25E2/75E1
        subplot(1,3,3)
        DepthChangeingSpiking_SM(Spiking25_75, 1:size(Spiking25_75,2), fignum,depthdriven,size(Spiking25_75,2),AMP_orig);%plots heatmap of estimated plots
        title(['ESTIMATION Stimchn: ' num2str(CHN(chosenstimchn)) ' ' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1) ' @ 75/25' ])
        yline((depthdriven-50*(CHN(chosenstimchn)-1)),'Color','r','Linewidth',75/25+25/35,'Alpha',1)

        ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
        ax.Label.String='Sp/s';
        ax.Label.Rotation=270;

AMPInterestSingleLinePlot=ampinterest;
[~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP_orig(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will te ll you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)

%LINECUT
figure %plotting estimation line plot
ax = gca;
ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
hold on
SMOOTHING=1;
window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);

rate = conv(Spiking75_25(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
%rate = Additivespkquater_Eone;
plot(rate)
rate = conv(Spiking50_50(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
%rate = Additivespkquater_Eone;
plot(rate)
rate = conv(Spiking25_75(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
%rate = Additivespkquater_Eone;
plot(rate)

xline((CHN(chosenstimchn)+NORECORDELECT(1)+1),'-.k')
xline((CHN(chosenstimchn)),'k')
legend('75/25', '50/50','25/75', 'Stim E1','Stim E2')
if includeResult==1
  title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' allresults ESTIMATION'])
else
  title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' SingleChn ESTIMATION'])
end

        check=['T50_50_' num2str(CHN(chosenstimchn))];
        predictdatastruct.(check) = Spiking50_50;
        check=['T75_25_' num2str(CHN(chosenstimchn))];
        predictdatastruct.(check) = Spiking25_75;
        check=['T25_75_' num2str(CHN(chosenstimchn))];
        predictdatastruct.(check) = Spiking75_25;
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