function [p,MI,rSquare,cM,nT,fR,dp,thresh,phaseOut] = spikeWaveforms(PLOT,chan,tID,spikeRefBand,chkThresh,numBS)
dbstop if error
if nargin < 1; PLOT = zeros(6,1); end
if nargin < 2; chan = 1; end
if nargin < 3; tID = 1; end
if nargin < 4; spikeRefBand = [6, -251]; end
if nargin < 5; chkThresh = 0; thresh = 0; end
if nargin < 6; numBS = 2; end
if length(PLOT)<6
    x = size(PLOT,1); y = size(PLOT,2);
    if x>y; PLOT = [PLOT; zeros(6,1)];
    else; PLOT = [PLOT zeros(1,6)]; end
end
AMP = loadAMP + 1; nAmp = length(AMP);
if any(AMP == 1)
    ZERO = 1;
else
    ZERO = 0;
end
%% Variables
plotSW = PLOT(4); % Plot spike waveforms
plotPh = PLOT(3); % Plot preferred phase
plotRs = PLOT(1); % Plot raster
plotHs = PLOT(2); % Plot histogram
selPhase = [1,2,3,4,5,6]; % Phasebin (of 6) to look at
C = define_Colormap;
loadStimChn; nSci = length(stimChn); DUR = loadDUR; nDur = length(DUR); %#ok<NASGU>
p = zeros(size(tID,2),1); MI = zeros(size(tID,2),1); rSquare = zeros(size(tID,2),1); cM = zeros(size(tID,2),1);
if (ZERO)
    try
        p = zeros(size(tID,2)-1,1); MI = zeros(size(tID,2)-1,1); rSquare = zeros(size(tID,2)-1,1); cM = zeros(size(tID,2)-1,1);
        tID(2) = [];
        AMP(AMP == 1) = []; nAmp = length(AMP);
    catch
        p = zeros(size(tID,2),1); MI = zeros(size(tID,2),1); rSquare = zeros(size(tID,2),1); cM = zeros(size(tID,2),1);
    end
end
%% Loading in data
d = Depth; trig = loadTrig(0); sp = loadSpikes; sp = denoiseSpikes(sp{d(chan)},d(chan)); loadNREP;
output = trig_helper(spikeRefBand,chan,tID);
x = 0:0.1:1440;
y = cosd(x); thRate = zeros(nAmp,size(selPhase,2));
SPIKES = cell(length(tID),n_REP,size(selPhase,2)); SPIKER = cell(length(tID),n_REP,size(selPhase,2));
for ttID = 1:length(tID)
    tt = tID(ttID);
    spWN = [2.25 11]; % Spike selection window in ms
    psthWN = [-30 20]; % PSTH selection window in ms
    tcM = output(output(:,2) == tt,4);
    cM(ttID) = tcM(1);
    xl = -60;
    xr = 0;
    ind = output(output(:,2)==tt,:); nTrig = size(ind,1);    
    if (plotPh || plotRs || plotHs)
        R = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    end
    if (plotSW)
        S = figure; hold on; %#ok<*UNRCH>
    end
    for selP = selPhase
        tind = ind(ind(:,5)==selP,1);
        ttrig = trig(tind);
        nTrig(selP) = length(ttrig);
        %nTrig = 15;
        %% Initialise data matrices
        spWaves = []; threshSp = zeros(nTrig(selP),1);
        for n = 1:nTrig(selP)
            spWaves = [spWaves; sp(sp(:,1) >= ttrig(n)/30 + spWN(1) & sp(:,1) <= ttrig(n)/30 + spWN(2),2:end)];
            % For PSTHs
            SPIKES{ttID,n,selP} = sp(sp(:,1) >= ttrig(n)/30 + psthWN(1) & sp(:,1) <= ttrig(n)/30 + psthWN(2),1) - ttrig(n)/30;
            SPIKER{ttID,n,selP} = sp(sp(:,1) >= ttrig(n)/30 + spWN(1) & sp(:,1) <= ttrig(n)/30 + spWN(2),1) - ttrig(n)/30;
            threshSp(n) = sum(sp(:,1) >= ttrig(n)/30 + spWN(1) & sp(:,1) <= ttrig(n)/30 + spWN(2));
        end
        thRate(ttID,selP) = nanmean(threshSp(:)) ./ (diff(spWN)/1e3);
        thRateSEM(ttID,selP) = nanstd(threshSp(:)) ./ ((diff(spWN)/1e3)*sqrt(nTrig(selP)));
        if isnan(thRateSEM(ttID,selP))
            thRateSEM(ttID,selP) = 0;
        end
        %% Plot spike waveforms
        if (plotSW)
            figure(S);
            subplot(2,3,selP); hold on;
            plot(spWaves');
            meanWave = mean(spWaves);
            plot(meanWave,'Color','k','LineWidth',2);
            beautifyPlot;
            ylim([-200 200]);
            xlim([0 50]);
            yticks('');
            if selP == 1 || selP == 4
                yticks([-200 -100 0 100 200]);
                ylabel('Voltage (\muV)');
            end
            xticks('');
            if selP == 4 || selP == 5 || selP == 6
                xticks([0 15 30 45]);
                xticklabels({'0' '0.5' '1' '1.5'});
                xlabel('Time (ms)');
            end
            title(['Phase Bin: ' num2str(selP) ' | Spikes: ' num2str(size(spWaves,1))]);
        end
        %% Plot
        if (plotPh || plotRs || plotHs)
            figure(R);
            if plotPh
                xl = xl + 60; xr = xr + 60;
                patch([xl xr xr xl],[-1 -1 1 1],C{selP},'EdgeColor',C{selP},'FaceColor',C{selP});
                text(30+(selP-1)*60,1.5,num2str(selP),'FontSize',24,'HorizontalAlignment','center');                
                plot(x(1:3600),y(3601 + (cM(ttID)-30)*10:7200 + (cM(ttID)-30)*10),'Color','k','LineWidth',3);            
            end
        end
    end
    if (plotPh || plotRs || plotHs)
        figure(R);
        xticks(0:60:360);
        XTL = cell(1,7);
        for i = 1:7
            XTL{i} = cM(ttID)-30 + (i-1)*60;
            if XTL{i} < 360; XTL{i} = XTL{i} + 360; end
            if XTL{i} > 360; XTL{i} = XTL{i} - 360; end
        end
        for i = 1:7
            XTL{i} = [num2str(XTL{i}) '\circ'];
        end
        xticklabels(XTL);
        if plotRs
            % Now we need to plot rasters
            nTrials = max(nTrig); yWN = [-3.1 -1.1];
            xsc = diff(psthWN); xpix = 60/xsc;
            ysc = diff(yWN); ypix = ysc/nTrials;
            for selP = selPhase
                line([xpix*abs(psthWN(1))+(selP-1)*60 xpix*abs(psthWN(1))+(selP-1)*60],[yWN(1) yWN(2)],'Color',C{selP},'LineWidth',2);
                line([xpix*abs(xsc)+(selP-1)*60 xpix*abs(xsc)+(selP-1)*60],[yWN(1) yWN(2)],'Color',[0.5 0.5 0.5],'LineWidth',2);
                pxl = xpix*(abs(psthWN(1))+spWN(1))+(selP-1)*60; pxr = xpix*(abs(psthWN(1))+spWN(2))+(selP-1)*60;
                pyb = yWN(1); pyt = yWN(2);
                patch([pxl pxr pxr pxl],[pyb pyb pyt pyt],[0.9 0.9 0.9],'EdgeColor','none');
                endoT = 0;
                for n = 1:nTrials
                    spt = xpix*(SPIKES{ttID,n,selP} + abs(psthWN(1)));
                    if n > nTrig(selP)
                        if endoT == 0
                            line([0+(selP-1)*60 xpix*xsc+(selP-1)*60],[yWN(2)-(n-0.5)*ypix yWN(2)-(n-0.5)*ypix],'Color','r');
                            endoT = 1;
                        end
                        continue;
                    end
                    for i = 1:length(spt)
                        if spt(i) >= 56
                            continue;
                        end
                        if spt(i) <= 1.5
                            continue;
                        end
                        line([spt(i)+(selP-1)*60 spt(i)+(selP-1)*60],[yWN(2)-n*ypix yWN(2)-(n-1)*ypix],'Color','k','LineWidth',2);
                    end
                end
                % X-axis for rasters
                line([0 + (selP-1)*60 xpix*abs(xsc-1) + (selP-1)*60],[yWN(1)-0.1 yWN(1)-0.1],'Color','k','LineWidth',2);
                line([xpix*abs(psthWN(1)+25) + (selP-1)*60 xpix*abs(psthWN(1)+25) + (selP-1)*60],[yWN(1)-0.15 yWN(1)-0.1],'Color','k','LineWidth',2);
                text(xpix*abs(psthWN(1)+25) + (selP-1)*60,yWN(1)-0.2,'-25 ms','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','top');
                line([xpix*abs(psthWN(1)) + (selP-1)*60 xpix*abs(psthWN(1)) + (selP-1)*60],[yWN(1)-0.15 yWN(1)-0.1],'Color','k','LineWidth',2);
                text(xpix*abs(psthWN(1)) + (selP-1)*60,yWN(1)-0.2,'0','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','top');
                line([xpix*(abs(psthWN(1))+5) + (selP-1)*60 xpix*(abs(psthWN(1))+5) + (selP-1)*60],[yWN(1)-0.15 yWN(1)-0.1],'Color','k','LineWidth',2);
                text(xpix*(abs(psthWN(1))+5) + (selP-1)*60,yWN(1)-0.2,'5','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','top');
                line([xpix*(abs(psthWN(1))+10) + (selP-1)*60 xpix*(abs(psthWN(1))+10) + (selP-1)*60],[yWN(1)-0.15 yWN(1)-0.1],'Color','k','LineWidth',2);
                text(xpix*(abs(psthWN(1))+10) + (selP-1)*60,yWN(1)-0.2,'10','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','top');
                if (selP == 1)
                    line([0 0],yWN,'Color','k','LineWidth',2);
                    line([-5 0],[yWN(2)-0.5*ypix yWN(2)-0.5*ypix],'Color','k','LineWidth',2);
                    text(-6,yWN(2)-0.5*ypix,'Trial 1','FontSize',16,'HorizontalAlignment','right');
                    line([-5 0],[yWN(2)-9.5*ypix yWN(2)-9.5*ypix],'Color','k','LineWidth',2);
                    text(-6,yWN(2)-9.5*ypix,'Trial 10','FontSize',16,'HorizontalAlignment','right');
                    line([-5 0],[yWN(2)-19.5*ypix yWN(2)-19.5*ypix],'Color','k','LineWidth',2);
                    text(-6,yWN(2)-19.5*ypix,'Trial 20','FontSize',16,'HorizontalAlignment','right');
                    line([-5 0],[yWN(2)-(nTrig(selP)-0.5)*ypix yWN(2)-(nTrig(selP)-0.5)*ypix],'Color','k','LineWidth',2);
                    text(-6,yWN(2)-(nTrig(selP)-0.5)*ypix,['Trial ' num2str(round(nTrig(selP)))],'FontSize',16,'HorizontalAlignment','right');
                end
            end
        end
    end
    meanSp = zeros(1,selP); steSp = zeros(1,selP);
    for selP = selPhase
        spt = zeros(nTrig(selP),1);
        for n = 1:nTrig(selP)
            spt(n) = length(SPIKER{ttID,n,selP});
        end
        meanSp(selP) = nanmean(spt);
        steSp(selP) = nanstd(spt) ./ sqrt(nTrig(selP));
        if isnan(meanSp(selP)); meanSp(selP) = 0; end
        if isnan(steSp(selP)); steSp(selP) = 0; end
    end
    rateSp = 1000.* meanSp ./ diff(spWN);
    period = length(selPhase); xval = pi/3:pi/3:(3*length(selPhase))*pi/3;
    fO = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf (2*pi)/period -inf],'Upper',[inf (2*pi)/period inf]);
    data = rateSp - mean(rateSp);
    data = [data data data]; %#ok<*AGROW>
    [f, fS] = fit(xval',data','sin1',fO);
    rSquare(ttID) = fS.rsquare;
    fy = f(pi/3:0.01:4*pi+pi/3) + mean(rateSp);
    MI(ttID) = diff([min(fy) max(fy)]);
    if (plotPh || plotRs || plotHs)
        figure(R);
        if plotHs
            barWN = [-5 -4]; barsc = diff(barWN);
            steSp = barsc .* (steSp ./ max(meanSp));
            meanSp = barsc .* (meanSp ./ max(meanSp)) + barWN(1);
            for selP = selPhase
                bar(30+(selP-1)*60,meanSp(selP),55,'FaceColor',[C{selP}],'EdgeColor',C{selP},'BaseValue',barWN(1));
                errorbar(30+(selP-1)*60,meanSp(selP),steSp(selP),'Color',C{selP},'LineWidth',2);
                text(30+(selP-1)*60,-3.5,[num2str(rateSp(selP),'%3.0f') ' Sp/s'],'FontSize',18,'HorizontalAlignment','center');
            end
            beautifyPlot;
            yticks('');
        end
    end
    if plotHs; YLIM(2) = -3; end; if plotRs; YLIM(2) = -1; end; if plotPh; YLIM(2) = 2; end
    if plotPh; YLIM(1) = -1; end; if plotRs; YLIM(1) = -4; end; if plotHs; YLIM(1) = -5.5; end
    if (plotPh || plotRs || plotHs); xlim([-60 420]); ylim(YLIM); end
    % Perform a rank-sum test
    B = find(meanSp == max(meanSp),1);
    W = (B - 3); if (W <= 0); W = W + 6; end
    dif = meanSp(B) - meanSp(W);
    Wt = find(meanSp == min(meanSp),1);
    Bt = (Wt - 3); if (Bt <= 0); Bt = Bt + 6; end
    dift = meanSp(Bt) - meanSp(Wt);
    if dift > dif
        B = Bt; W = Wt;
    end
    rsX = zeros(nTrig(B),1);
    for i = 1:nTrig(B)
        rsX(i,1) = length(SPIKER{ttID,i,B});
    end
    rsY = zeros(nTrig(W),1);
    for i = 1:nTrig(W)
        rsY(i,1) = length(SPIKER{ttID,i,W});
    end
    if sum(~isnan(rsX)) == 0 || sum(~isnan(rsY)) == 0
        p(ttID) = NaN;
    else
        p(ttID) = ranksum(rsX,rsY,'tail','right');
    end
    if ~isempty(B) && ~isempty(W) && (plotHs)
        line([(30+(B-1)*60) (30+(W-1)*60)],[-5.1 -5.1],'Color','k','LineWidth',2);
        text((30+((B+W)/2-1)*60),-5.2,['p = ' num2str(p(ttID))],'FontSize',24,'HorizontalAlignment','center');
    end
    % Calculate d'
    dprime = (mean(rsX)-mean(rsY))./sqrt(0.5.*((std(rsX)^2)+(std(rsY)^2)));
    if isnan(dprime) || isinf(dprime); dprime = -1; end
    nT(ttID,:) = nTrig;
    fR(ttID) = max(rateSp);
    dp(ttID) = dprime;
end
if (chkThresh)
    AMP(2:end) = AMP(2:end)-1;
%     data = cell(8,6);
%     for nA = 1:8
%         for nP = 1:6
%             for i = 1:nT(nA,nP)
%                 data{nA,nP} = [data{nA,nP}; length(SPIKER{nA,i,nP})];
%             end
%         end
%     end
%     [stats] = permutationTest(data,numBS,AMP,thRate);
    numBootstrap = numBS; perc_data = 0.8;
    for bs = 1:numBootstrap
        for bstID = 1:length(tID)
            for sp = 1:length(selPhase)
                index = datasample(1:nT(bstID,sp),ceil(perc_data*nT(bstID,sp)));
                bsRate(bstID,sp,bs) = (length(cell2mat(SPIKER(bstID,index,sp)')) ./ ceil(perc_data*nT(bstID,sp))) ./ (diff(spWN)/1e3);                
            end
        end
    end    
    % Calculate the Naka-Rushton parameters.
    curveFit = @(prm,c)(prm(1).*c.^prm(2)./(c.^(prm(2).*prm(3))+prm(4).^(prm(2).*prm(3)))+prm(5));
    opts = optimset('Display','off');
    thRate(isnan(thRate)) = 0;
    prm0 = lsqcurvefit(curveFit,[0 1 1 1 mean(thRate(1,:))],AMP(1:end),mean(thRate(1:length(tID),:),2)',[max(mean(thRate(1:length(tID),:),2)) 0 1 0 mean(thRate(1,:))],[max(mean(thRate(1:length(tID),:),2)) 10 1 AMP(end) mean(thRate(1,:))],opts);
    for selP = 1:length(selPhase)
        % Calculate the threshold for this trial ID
        prm1 = lsqcurvefit(curveFit,[prm0(1), prm0(2), prm0(3), prm0(4), prm0(5)],AMP(1:end),thRate(1:length(tID),selP)',[prm0(1), prm0(2), prm0(3), 0, prm0(5)],[prm0(1), prm0(2), prm0(3), AMP(end), prm0(5)],opts);
        for l = 0.0001:0.0001:AMP(end)+100
            CURVE(selP,floor(l*1000)+1) = curveFit(prm1,l);
        end
        % Calculate r2
        sse = sum((thRate(1:length(tID),selP)'-CURVE(selP,1000*AMP+1)).^2);
        sst = sum((thRate(1:length(tID),selP)-mean(CURVE(selP,:))).^2);
        r2(selP) = 1 - sse/sst;
    end
    % Calculate error
    Error = zeros(numBootstrap,length(selPhase),1);
    for bs = 1:numBootstrap
        for selP = 1:length(selPhase)            
            XDATA = bsRate(1:length(tID),selP,bs);
            XDATA(isnan(XDATA)) = 0;
            prm2 = lsqcurvefit(curveFit,[prm0(1), prm0(2), prm0(3), prm0(4), prm0(5)],AMP(1:end),XDATA',[prm0(1), prm0(2), prm0(3), 0, prm0(5)],[prm0(1), prm0(2), prm0(3), AMP(end), prm0(5)],opts);
            for l = 0.0001:0.0001:AMP(end)
                Error(bs,selP,floor(l*1000)+1) = curveFit(prm2,l);
            end
        end
    end
    MIN = mean(CURVE(:,1));
    MAX = mean(CURVE(:,end));
    HG = max(CURVE,[],'all');
    distThresh = [];
    for selP = selPhase
        try
            thresh{selP} = find(CURVE(selP,:) > (MIN + 0.5*diff([MIN MAX])),1) ./ 1000;
            if thresh{selP} <= 0 || isempty(thresh{selP})
                thresh{selP} = NaN;
            end
        catch
            thresh{selP} = NaN;
        end
        for bs = 1:numBootstrap
            try
                distThresh{bs,selP} = find(Error(bs,selP,:) > (MIN + 0.5*diff([MIN MAX])),1) ./ 1000;
                if distThresh{bs,selP} <= 0 || isempty(distThresh{bs,selP})
                    distThresh{bs,selP} = NaN;
                end
            catch
                distThresh{bs,selP} = NaN;
            end
        end
        errorThresh{selP} = std(cell2mat(distThresh(:,selP)));
    end
    % Bootstrap Test Statistics
    myStat = @(x1,x2) mean(x1)-mean(x2); % Difference in means
    bsStat = zeros(6,6,numBootstrap);
    alpha = 0.05; %#ok<NASGU>
    for i = 1:6
        for j = 1:6
            for n = 1:numBootstrap
                sampleX1 = distThresh{ceil(rand(numBootstrap,1)*numBootstrap),i};
                sampleX2 = distThresh{ceil(rand(numBootstrap,1)*numBootstrap),j};
                bsStat(i,j,n) = myStat(sampleX1,sampleX2);
            end
            % Calculate the confidence interval
            CI(i,j,:) = prctile(bsStat(i,j,:),[100*0.05/2,100*(1-0.05/2)]);
            % Test the hypothesis
            H(i,j) = CI(i,j,1)>0 | CI(i,j,2)<0;
        end
    end
    % p-value
    bsP = []; pI = []; pJ = [];
    [i,j] = find(H == 1);
    for n = 1:size(i,1)
        if H(i(n),j(n)) ~= H(j(n),i(n))
            bsP(n) = 1;
            pI(n) = i(n);
            pJ = j(n);
            continue;
        end
        pI(n) = i(n);
        pJ(n) = j(n);
        %if CI(i(n),j(n),1)>0
        tmp = squeeze(sort(bsStat(i(n),j(n),:)));
        x = length(tmp(tmp <= 0));
        y = length(tmp(tmp > 0));
        bsP(n) = x ./ (x+y);
        %else
        %             tmp = squeeze(sort(bsStat(i(n),j(n),:)));
        %             x = length(tmp(tmp >= 0));
        %             y = length(tmp(tmp < 0));
        %             bsP(n) = x ./ (x+y);
        %end
    end
    [~,ii] = min(bsP);
    pI = pI(ii);
    pJ = pJ(ii);
    bsP = bsP(ii);
    phaseOut = [pI, pJ, bsP];
    if (PLOT(6))
        if isempty(pI)
            [~,pI] = min(cell2mat(thresh));
        end
        if isempty(pJ)
            [~,pJ] = max(cell2mat(thresh));
        end
        figure; hold on;
        subplot(2,1,1); hold on;
        for i = 1:6
            scatter(mean(cell2mat(distThresh(:,i))),i,200,[C{i} 'd'],'filled');
            errorbar(mean(cell2mat(distThresh(:,i))),i,std(cell2mat(distThresh(:,1))),'Horizontal','Color',C{i});
        end
        ylabel('Phasebins');
        yticks(1:6);
        yticklabels({'0\circ','60\circ','120\circ','180\circ','240\circ','300\circ',});
        xlabel('Threshold (\muA)');
        title('Mean \pm STD | Bootstrapped Phase Thresholds');
        ylim([0 7]);
        beautifyPlot(18);
        subplot(2,1,2); hold on;
        histogram((bsStat(pI,pJ,:)),'FaceColor',[0.5 0.5 0.5]);        
        YLIM = ylim;        
        ylabel('# of Bootstrapped Samples');
        xlabel('Difference in Mean Threshold (\muA)');        
        title(['Bootstrapped Difference in Mean Threshold Distribution | ' num2str(numBootstrap) ' Samples | Phasebin ' num2str(pI) ' - Phasebin ' num2str(pJ)]);
        line([0 0],[YLIM(1) YLIM(2)],'Color','k','LineStyle','--','LineWidth',2);
        text(-0.1,YLIM(2)*0.8,['p = ' num2str(bsP)],'FontSize',18,'HorizontalAlignment','right');
        actualdThresh = thresh{pI} - thresh{pJ};
        line([actualdThresh actualdThresh],[YLIM(1) YLIM(2)],'Color','r','LineWidth',2,'LineStyle','--');
        text(2,90,'Real difference in threshold','FontSize',18);
        beautifyPlot(18);
    end
    if (PLOT(5))
        figure; subplot(1,2,1); hold on;
        for selP = selPhase            
            plot(CURVE(selP,:),'Color',C{selP},'LineWidth',2);
            line([(thresh{selP})*1000 (thresh{selP})*1000],[0 HG*1.2],'Color',C{selP},'LineStyle','--','HandleVisibility','off','LineWidth',2);
            %if selP == 1 || selP == 4
            %    scatter(1000*AMP,thRate(1:length(tID),selP),300,C{selP},'x','HandleVisibility','off');
            %end
        end        
        ylim([0 HG*1.2]);
        line([0 15000],[(MIN + 0.5*diff([MIN MAX])) (MIN + 0.5*diff([MIN MAX]))],'Color','r','LineStyle','--');
        legend({['r^2 ' num2str(r2(1),'%1.2f')],['r^2 ' num2str(r2(2),'%1.2f')],['r^2 ' num2str(r2(3),'%1.2f')],['r^2 ' num2str(r2(4),'%1.2f')],['r^2 ' num2str(r2(5),'%1.2f')],['r^2 ' num2str(r2(6),'%1.2f')],'Threshold'},'Location','southeast');
        beautifyPlot(24); ylabel('Spike Rate (Sp/s)');
        xticks((0:14)*1000); xlim([-1000 (AMP(end)+1)*1000]);
        xticklabels({'n.s.','1','2','3','4','5','6','7','8','9','10','11','12','13','14'});
        xlabel('Stimulus Current (\muA)');
        subplot(1,2,2); hold on; beautifyPlot(24);
        for selP = selPhase
            bar(selP,mean(cell2mat(distThresh(:,selP))),'EdgeColor','none','FaceColor',C{selP});
            errorbar(selP,mean(cell2mat(distThresh(:,selP))),errorThresh{selP},'Color',C{selP},'HandleVisibility','off');
        end
        YLIM = [0 ceil(max(mean(cell2mat(distThresh),1))+max(cell2mat(errorThresh)))+1];
        ylim(YLIM);
        if ~isempty(bsP)
            line([pI pJ],[0.85*YLIM(2) 0.85*YLIM(2)],'Color','k','LineWidth',2);
            text(min([pI,pJ]),0.9*YLIM(2),['Bootstrap, p = ' num2str(bsP)],'FontSize',24,'HorizontalAlignment','left');
        end
        ylabel('Mean Threshold \pm STE (\muA)');
        xlabel('Phase Bins (\circ) - Rotated');
        xticks(selPhase);
%         XTL = cell(1,7);
%         for i = 1:7
%             XTL{i} = cM(3)-30 + (i-1)*60;
%             if XTL{i} < 360; XTL{i} = XTL{i} + 360; end
%             if XTL{i} > 360; XTL{i} = XTL{i} - 360; end
%         end
%         for i = 1:7
%             XTL{i} = [num2str(XTL{i}) '\circ'];
%         end
        xticklabels({'0\circ','60\circ','120\circ','180\circ','240\circ','300\circ'});
        %sgtitle('RAT0026 | Stimulation Site 6 | Recording Site 8','FontSize',24);
    end
end
end