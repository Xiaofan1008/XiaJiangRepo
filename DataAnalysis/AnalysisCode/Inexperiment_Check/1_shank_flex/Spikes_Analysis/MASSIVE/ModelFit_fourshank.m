function [predictdatastruct]=ModelFit_fourshank(AMPInterestSingleLinePlot)


load('truedatastruct.mat','truedatastruct')
trialinfo=loadTrialInfo(0);
PlotNow=0;
includeResult=1;
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

if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
    laminar=1;
else
    laminar=0;
end
electbeta=zeros(nChn,4);
Spiking50_50=zeros(nChn,length(AMP_orig));
Spiking25_75=zeros(nChn,length(AMP_orig));
Spiking75_25=zeros(nChn,length(AMP_orig));
electfit=zeros(nChn,1);
%warning('off','all')
for ampiterate=1:length(AMP_orig)
for elect=1:nChn
%     AMPInterestSingleLinePlot=AMP_orig(ampiterate);
%     [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot/2)); %x is a throw away variable but it will te ll you how close the value is to your chosen val
%     AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
%     
if includeResult==1
    step=2;
else
    step=1;
end
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

    [level,pos]=intersect(AMP_double, AMP_orig);
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
%     y1plat=max(electrodeampint100);
%     y2plat=max(electrodeampint0);
%     y0=0;
%     a1=zeros(length(AMP(1:2:end)),1);
%     a2=zeros(length(AMP(1:2:end)),1);
%     y1=electrodeampint100;
%     y2=electrodeampint0;
%         for i=1:length(AMP(1:2:end))
%             amid=AMP(i*2-1);
%             a1(i)=log(-((y1(AMP(1:2:end)==amid)-y0)/(y1plat-y0)) +1)/(-amid);
%             a2(i)=log(-((y2(AMP(1:2:end)==amid)-y0)/(y2plat-y0)) +1)/(-amid);
%         end
        
    ampplat=AMP(1:end);

    ampplat1=ampplat(find(electrodeampint100==max(electrodeampint100),1,'first'));
    ampplat2=ampplat(find(electrodeampint0==max(electrodeampint0),1,'first'));

    a1=zeros(length(AMP(1:step:end)),1);
    a2=zeros(length(AMP(1:step:end)),1);
    for i=1:length(AMP(1:step:end))
        amid=AMP(i*2-1);
        a1(i)=(electrodeampint100((find(AMP(1:step:end)==amid,1,'first')))-max(electrodeampint100))/((amid-ampplat1)^2);
        a2(i)=(electrodeampint0((find(AMP(1:step:end)==amid,1,'first')))-max(electrodeampint0))/((amid-ampplat2)^2);
    end
%      a1(isnan(a1))=[];
%      a2(isnan(a2))=[];
%      a1(isinf(a1))=[];
%      a2(isinf(a2))=[];
%     a1=sum(a1)/length(a1);
%     a2=sum(a2)/length(a2);
%     ampplat=AMP(1:2:end);
%     ampplat2=ampplat(find(electrodeampint0==max(electrodeampint0),1,'first'));
%         a22=-max(electrodeampint0)/((-ampplat2)^2);
%     a2=(a2+a22)/2;
    k1=max(electrodeampint100);
    k2=max(electrodeampint0);
    
        %T = [0 0 0;1 0 0;0 1 0;1 1 0]; for multiplicative
    %modelfun = @(b,X) b(1)*X(:,1)+b(2)*X(:,2)+b(3); %linear
     modelfun = @(b,X) b(1)./(1+exp(b(2)+b(3).*X(:,1)+b(4).*X(:,2))); %SIGMOIDAL  b1 is asymptote, b2=x50/slope b3=-1/slope b2=-1/slope  This is boltzmann sigmoid
    beta0 = [350 1 1 1]; %SIGMOIDAL
    %modelfun= @(b,X) b(1).*(X(:,1)-b(2)).^2+b(3) + b(4).*(X(:,2)-b(5)).^2 + b(6);
    %modelfun = @(b,X) b(1)+b(2)*X(:,1)+b(3)*X(:,2)+b(4)*X(:,1).^2+b(5)*X(:,2).^2;%quadratic
   %     modelfun = @(b,x) (b(1).*X(:,1).^b(2))./((X(:,1).^b(2))+ b(3)^b(2))  + (b(4).*X(:,2).^b(5))./((X(:,2).^b(5))+ b(6)^b(5)) ;
    %modelfun = @(b,X) b(1)+b(1)*exp(b(2)*X(:,1))+b(3)+b(3)*exp(b(4)*X(:,2));
    %beta0 = [1 1 1]; %linear
   % beta0 = [1 1 1 1 1]; %quadratic
    %beta0 = [1 10 350 1 10 350];
   % beta0 = [350 1 10 350 1 100];
   % beta0 = [0 350 100 1 0 350 100 1];
    %try
        mdl = fitnlm(X,y',modelfun,beta0);
    %catch
      %  continue
    %end
    fprintf([num2str(elect) '\n']);
    %mdl = fitlm(X',y');
    beta=mdl.Coefficients.Estimate;
    electbeta(elect,:)=beta;
    Spiking50_50(elect,ampiterate)=beta(1)./(1+exp(beta(2)+beta(3)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL) + beta(4)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)));
    Spiking75_25(elect,ampiterate)=beta(1)./(1+exp(beta(2)+beta(3)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL) + beta(4)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)));
    Spiking25_75(elect,ampiterate)=beta(1)./(1+exp(beta(2)+beta(3)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL) + beta(4)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)));
%     Spiking50_50(elect,ampiterate)=beta(1)+beta(2)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)+(beta(3)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL))+beta(4)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL).^2+beta(5)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL).^2;
%     Spiking75_25(elect,ampiterate)=beta(1)+beta(2)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)+(beta(3)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL))+beta(4)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL).^2+beta(5)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL).^2;
%     Spiking25_75(elect,ampiterate)=beta(1)+beta(2)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)+(beta(3)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL))+beta(4)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL).^2+beta(5)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL).^2;
% %     Spiking50_50(elect,ampiterate)=beta(1)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)+beta(2)*AMP_half(AMPInterestSingleLinePlotINDEXDUAL)+beta(3);
%     Spiking75_25(elect,ampiterate)=beta(1)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)+beta(2)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)+beta(3);
%     Spiking25_75(elect,ampiterate)=beta(1)*AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL)+beta(2)*AMP_quater(AMPInterestSingleLinePlotINDEXDUAL)+beta(3);


    if PlotNow==1
        x1fit = min(x1):0.2:max(x1);
        x2fit = min(x2):0.2:max(x2);
        [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
        YFIT = beta(1)./(1+exp(beta(2)+beta(3).*X1FIT + beta(4).*X2FIT));
        %beta=[-1.883 10 188.3 -2.46 10 246.7];
        %beta=[a1 ampplat1 k1 a2 ampplat2 k2];
        %YFIT = y0+(y1plat-y0).*(1-exp(-a1*X1FIT))+y0+(y2plat-y0).*(1-exp(-a2*X2FIT));
        %YFIT=beta(1).*(X1FIT-beta(2)).^2+beta(3) + beta(4).*(X2FIT-beta(5)).^2 + beta(6);
        %YFIT = beta(1)+(beta(2)-beta(1))./(1+10^(log(beta(3)-X1FIT))) + beta(4)+(beta(5)-beta(4))./(1+10^(log(beta(6)-X2FIT)));
        %YFIT = (beta(1).*X1FIT.^beta(2))./((X1FIT.^beta(2))+ beta(3)^beta(2))  + (beta(4).*X2FIT.^beta(5))./((X2FIT.^beta(5))+ beta(6)^beta(5)) ;
        %YFIT = beta(1)+(beta(2)-beta(1))./(1+exp((beta(3)+X1FIT)./beta(4))) + beta(5)+(beta(6)-beta(5))./(1+exp((beta(7)+X2FIT)./beta(8))); %SIGMOIDAL
        if k1>k2
            k=k1;
        else
            k=k2;
        end
        YFIT(YFIT>k)=k;
        YFIT(YFIT<0)=0;
        beta=[ 350 3 -0.5 -0.5];
        X1FIT=1:10;
        X2FIT=1:10;
        %YFIT=beta(1)+beta(2)*X1FIT+(beta(3)*X2FIT)+beta(4)*X1FIT.^2+beta(5)*X2FIT.^2;
        %YFIT=beta(1)*X1FIT+beta(2)*X2FIT+beta(3);
        %YFIT=beta(1)+beta(1)*exp(beta(2)*X1FIT)+beta(3)+beta(3)*exp(beta(4)*X2FIT)
        figure
        m=mesh(X1FIT,X2FIT,YFIT);
        %s.EdgeColor='r';
        hold on
        scatter3(x1,x2,y,'filled')
        xlabel(['Channel ' num2str(CHN(chosenstimchn)) ' Current level (uA)'])
        ylabel(['Channel ' num2str(NORECORDELECT(1)) ' Current level (uA)'])
        zlabel('Sp/s')
        title(['Spike rates from electrode ' num2str(elect)])
    end
    %MSEtotal=MSEtotal+mdl.MSE;
    
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
    X=[x1' x2'];
    %beta=[a1 max(AMP) k1 a2 max(AMP) k2];
    %YFIT=beta(1)+beta(1)*exp(beta(2)*x1)+beta(3)+beta(3)*exp(beta(4)*x2);
    %YFIT=beta(1).*(x1-beta(2)).^2+beta(3) + beta(4).*(x2-beta(5)).^2 + beta(6);
    %YFIT = (beta(1).*x1.^beta(2))./((x1.^beta(2))+ beta(3)^beta(2))  + (beta(4).*x2.^beta(5))./((x2.^beta(5))+ beta(6)^beta(5)) ;
    %YFIT = beta(1)+(beta(2)-beta(1))./(1+exp((beta(3)+x1+x2)./beta(4))); %SIGMOIDAL
    YFIT = beta(1)./(1+exp(beta(2)+beta(3).*x1 + beta(4).*x2));


    %YFIT = beta(1)+(beta(2)-beta(1))./(1+10^(log(beta(3)-x1))) + beta(4)+(beta(5)-beta(4))./(1+10^(log(beta(6)-x2)));
    if k1>k2
        k=k1;
    else
        k=k2;
    end
    YFIT(YFIT>k)=k;
    YFIT(YFIT<0)=0;
    %YFIT = y0+(y1plat-y0).*(1-exp(-a1*x1))+y0+(y2plat-y0).*(1-exp(-a2*x2));
    %YFIT=beta(1)+beta(2)*x1+(beta(3)*x2)+beta(4)*x1.^2+beta(5)*x2.^2;
    errory=y-YFIT;
    ymu=y-mean(y);
    ssd=1-sum(errory.^2)/sum(ymu.^2);
    if ssd>0.75
        electfit(elect)=1;
        fprintf('Good fit \n')
    end
    MSEtotal=MSEtotal+ssd;%mdl.MSE;
end
end



warning('on','all')
% %HEATMAPS
% %25/75
%         figure
%         fignum=gcf;
%         fignum=fignum.Number;
% %         subplot(1,3,1)
% %         DepthChangeingSpiking_SM(Spiking75_25, 1:size(Spiking75_25,2), fignum,depthdriven,size(Spiking75_25,2),AMP_orig);%plots heatmap of estimated plots
% %         title(['ESTIMATION Stimchn: ' num2str(CHN(chosenstimchn)) ' ' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1) ' @ 25/75' ])
% %         yline((depthdriven-50*(CHN(chosenstimchn)+NORECORDELECT(1)+1-1)),'Color','r','Linewidth',75/25+25/35,'Alpha',1)
% %         
%         subplot(1,4,1)
%         DepthChangeingSpiking_SM(Spiking75_25(1:16,:), 1:size(Spiking75_25,2),fignum,depthdriven,size(Spiking75_25,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<17) && (CHN(chosenstimchn)>0)
%             yline(depthdriven-(16-CHN(chosenstimchn))*50,'r')
%         end
%         if (NORECORDELECT(1)<17) && (NORECORDELECT(1)>0)
%             yline(depthdriven-(16-(NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 25/75 Shank 1')
% 
%         subplot(1,4,4)
%         DepthChangeingSpiking_SM(Spiking75_25(17:32,:), 1:size(Spiking75_25,2),fignum,depthdriven,size(Spiking75_25,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<33) && (CHN(chosenstimchn)>16)
%             yline(depthdriven-((32-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<32) && (NORECORDELECT(1)>16)
%             yline(depthdriven-((32-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 25/75 Shank 4')
%         
%         subplot(1,4,2)
%         DepthChangeingSpiking_SM(Spiking75_25(33:48,:), 1:size(Spiking75_25,2),fignum,depthdriven,size(Spiking75_25,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<49) && (CHN(chosenstimchn)>32)
%             yline(depthdriven-((48-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<49) && (NORECORDELECT(1)>32)
%             yline(depthdriven-((48-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 25/75 Shank 2')
%         
%         subplot(1,4,3)
%         DepthChangeingSpiking_SM(Spiking75_25(49:64,:), 1:size(Spiking75_25,2),fignum,depthdriven,size(Spiking75_25,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<65) && (CHN(chosenstimchn)>48)
%             yline(depthdriven-((64-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<65) && (NORECORDELECT(1)>48)
%             yline(depthdriven-((64-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 25/75 Shank 3')
%                 ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
%         ax.Label.String='Sp/s';
%         ax.Label.Rotation=270;
%         
%         
%         %50/50
%         figure
%         fignum=gcf;
%         fignum=fignum.Number;
% 
%         subplot(1,4,1)
%         DepthChangeingSpiking_SM(Spiking50_50(1:16,:), 1:size(Spiking50_50,2),fignum,depthdriven,size(Spiking50_50,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<17) && (CHN(chosenstimchn)>0)
%             yline(depthdriven-(16-CHN(chosenstimchn))*50,'r')
%         end
%         if (NORECORDELECT(1)<17) && (NORECORDELECT(1)>0)
%             yline(depthdriven-(16-(NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 50/50 Shank 1')
% 
%         subplot(1,4,4)
%         DepthChangeingSpiking_SM(Spiking50_50(17:32,:), 1:size(Spiking50_50,2),fignum,depthdriven,size(Spiking50_50,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<33) && (CHN(chosenstimchn)>16)
%             yline(depthdriven-((32-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<32) && (NORECORDELECT(1)>16)
%             yline(depthdriven-((32-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 50/50 Shank 4')
%         
%         subplot(1,4,2)
%         DepthChangeingSpiking_SM(Spiking50_50(33:48,:), 1:size(Spiking50_50,2),fignum,depthdriven,size(Spiking50_50,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<49) && (CHN(chosenstimchn)>32)
%             yline(depthdriven-((48-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<49) && (NORECORDELECT(1)>32)
%             yline(depthdriven-((48-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 50/50 Shank 2')
%         
%         subplot(1,4,3)
%         DepthChangeingSpiking_SM(Spiking50_50(49:64,:), 1:size(Spiking50_50,2),fignum,depthdriven,size(Spiking50_50,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<65) && (CHN(chosenstimchn)>48)
%             yline(depthdriven-((64-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<65) && (NORECORDELECT(1)>48)
%             yline(depthdriven-((64-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 50/50 Shank 3')
%                 ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
%         ax.Label.String='Sp/s';
%         ax.Label.Rotation=270;
%         
%         
%         %25E2/75E1
%         figure
%         fignum=gcf;
%         fignum=fignum.Number;
% 
%         subplot(1,4,1)
%         DepthChangeingSpiking_SM(Spiking25_75(1:16,:), 1:size(Spiking25_75,2),fignum,depthdriven,size(Spiking25_75,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<17) && (CHN(chosenstimchn)>0)
%             yline(depthdriven-(16-CHN(chosenstimchn))*50,'r')
%         end
%         if (NORECORDELECT(1)<17) && (NORECORDELECT(1)>0)
%             yline(depthdriven-(16-(NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 75/25 Shank 1')
% 
%         subplot(1,4,4)
%         DepthChangeingSpiking_SM(Spiking25_75(17:32,:), 1:size(Spiking25_75,2),fignum,depthdriven,size(Spiking25_75,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<33) && (CHN(chosenstimchn)>17)
%             yline(depthdriven-((32-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<32) && (NORECORDELECT(1)>16)
%             yline(depthdriven-((32-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 75/25 Shank 4')
%         
%         subplot(1,4,2)
%         DepthChangeingSpiking_SM(Spiking25_75(33:48,:), 1:size(Spiking25_75,2),fignum,depthdriven,size(Spiking25_75,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<49) && (CHN(chosenstimchn)>32)
%             yline(depthdriven-((48-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<17) && (NORECORDELECT(1)>0)
%             yline(depthdriven-((48-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 75/25 Shank 2')
%         
%         subplot(1,4,3)
%         DepthChangeingSpiking_SM(Spiking25_75(49:64,:), 1:size(Spiking25_75,2),fignum,depthdriven,size(Spiking25_75,2),AMP_orig); %plots heat maps
%         if (CHN(chosenstimchn)<65) && (CHN(chosenstimchn)>48)
%             yline(depthdriven-((64-CHN(chosenstimchn)))*50,'r')
%         end
%         if (NORECORDELECT(1)<65) && (NORECORDELECT(1)>48)
%             yline(depthdriven-((64-NORECORDELECT(1)))*50,'r')
%         end
%         title('ESTIMATION 75/25 Shank 3')
%         ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
%         ax.Label.String='Sp/s';
%         ax.Label.Rotation=270;

AMPInterestSingleLinePlot=ampinterest;
[~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP_orig(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will te ll you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)

%LINECUT

    figure %plotting estimation line plot
hold on
SMOOTHING=1;
window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
for shankplot=1:4
    
for typeplot=1:3
        if typeplot==1
            spiking=Spiking75_25;
        elseif typeplot==2
            spiking=Spiking50_50;
        else
            spiking=Spiking25_75;
        end
    if shankplot==1
        subplot(1,4,1)
        acrossshankplot=spiking(1:16,AMPInterestSingleLinePlotINDEXDUAL);
        title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(1)])
        if typeplot==1
        if (cell2mat(trialinfo(2*2,2)))<17 && (cell2mat(trialinfo(2*2,2))>0)
            yline((cell2mat(trialinfo(2*2,2))),'k')
        end
        if (cell2mat(trialinfo(2*2-1,2)))<17 && (cell2mat(trialinfo(2*2-1,2))>0)
            yline((cell2mat(trialinfo(2*2-1,2))),'k')
        end
        end
    elseif shankplot==2
        subplot(1,4,4)
        acrossshankplot=spiking(17:32,AMPInterestSingleLinePlotINDEXDUAL);
        title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(4)])
        if typeplot==1
        if ((cell2mat(trialinfo(2*2,2)))<33) && cell2mat(trialinfo(2*2,2))>16
            yline((cell2mat(trialinfo(2*2,2))-16),'k')
        end
        if ((cell2mat(trialinfo(2*2-1,2)))<33) && (cell2mat(trialinfo(2*2-1,2))>16)
            yline((cell2mat(trialinfo(2*2-1,2))-16),'k')
        end
        end
    elseif shankplot==3
        subplot(1,4,shankplot-1)
        acrossshankplot=spiking(33:48,AMPInterestSingleLinePlotINDEXDUAL);
        title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(2)])
        if typeplot==1
        if ((cell2mat(trialinfo(2*2,2)))<49) && cell2mat(trialinfo(2*2,2))>32
            yline((cell2mat(trialinfo(2*2,2))-32),'k')
        end
        if ((cell2mat(trialinfo(2*2-1,2)))<49) && (cell2mat(trialinfo(2*2-1,2))>32)
            yline((cell2mat(trialinfo(2*2-1,2))-32),'k')
        end
        end
        
    else
        subplot(1,4,shankplot-1)
        acrossshankplot=spiking(49:64,AMPInterestSingleLinePlotINDEXDUAL);
        title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA Shank: ' num2str(3)])
        if typeplot==1
        if ((cell2mat(trialinfo(2*2,2)))<65) && (cell2mat(trialinfo(2*2,2))>48)
            yline((cell2mat(trialinfo(2*2,2))-48),'k')
        end
        if ((cell2mat(trialinfo(2*2-1,2)))<65) && (cell2mat(trialinfo(2*2-1,2))>48)
            yline((cell2mat(trialinfo(2*2-1,2))-48),'k')
        end
        end
    end
    
     if typeplot==1
         ax = gca;
         ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
     end
        rate = conv(acrossshankplot,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        %acrossshankplot(p,shankplot)=rate(7);
        hold on
        plot(rate,1:16);
    
    if shankplot==2
        lgd = legend;
        lgd.Position=[0.922,0.7036,0.06,0.149];
        %lgd.String={'Stim Elect','0/100','25/75','50/50','75/25','100/0'};
        lgd.String={'Stim Elect','25/75','50/50','75/25'};
    end
    %set(gca, 'YDir','reverse')
    xlabel('Sp/s')
    ylabel('Electrode')
    ylim([1 16])
end
end






% rate = conv(Spiking75_25(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
% rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
% %rate = Additivespkquater_Eone;
% plot(rate)
% rate = conv(Spiking50_50(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
% rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
% %rate = Additivespkquater_Eone;
% plot(rate)
% rate = conv(Spiking25_75(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
% rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
% %rate = Additivespkquater_Eone;
% plot(rate)
% 
% xline((CHN(chosenstimchn)+NORECORDELECT(1)+1),'-.k')
% xline((CHN(chosenstimchn)),'k')
% legend('75/25', '50/50','25/75', 'Stim E1','Stim E2')
% if includeResult==1
%   title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' allresults ESTIMATION'])
% else
%   title([num2str(AMP_orig(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' SingleChn ESTIMATION'])
% end

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