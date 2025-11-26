%% Single Channel Response Rate
function scResponse(CHANNEL,PLOT)
dbstop if error
n_REP = []; % Basic initialisation
if nargin < 2
    PLOT = 1;
end
%% Variables
stimSt = 2.25;
baselineBIN = [-175, -21];
stimulationBIN = [stimSt, 50];
dead_channels = 0;
if ~(nargin) || CHANNEL == 0
    CHANNEL = 1:32; % The Channel of Interest
end
mChn = 1; stimChn = 0; loadStimChn;
if size(CHANNEL,2) > 1
    mChn = 0;
    nChn = length(CHANNEL);
    if nChn == 32
        A = zeros(6,6);
        for i = 1:32
            A(i+3) = i;
        end
        for i = 1:28
            A(i+2) = i;
        end
        A(1,6) = 0;
        for i = 1:4
            A(i+1) = i;
        end
        A(6,1) = 0;
    end
    XLABEL = [25,30,32,33,34,35];
    YLABEL = [7,19,32];
end
responseFig = figure('units','normalized','outerposition',[0 0 1 1]);
if ~(mChn)
    s_axes = tight_subplot(6,6,[.025 .025],[.075 .080],[.05 .01]);
    delete(s_axes([1,6,31,36]));
end
name = split(pwd,'\');
n1 = name{end-1};
n2 = split(name{end},'_');
n3 = n2{3};
n2 = n2{2};
for thisChn = CHANNEL
    %% Loading
    TrialParams = loadTrialParams; trig = loadTrig(0); sp = loadSpikes; loadNREP; AMP = loadAMP; d = Depth; loadCHN;
    %% Initialisation of data structures
    nTr = length(AMP);
    nT = zeros(n_REP,nTr*length(CHN));
    TrialParams = cell2mat(TrialParams);
    COLOR = cell(1,nTr);COLOR{1} = [0.7 0.7 0.7];COLOR{2} = [0.6 0.6 0.6];COLOR{3} = [0.5 0.5 0.5];COLOR{4} = [0.4 0.4 0.4];
    COLOR{5} = [0.3 0.3 0.3];COLOR{6} = [0.2 0.2 0.2];COLOR{7} = [0.1 0.1 0.1];COLOR{8} = [0 0 0];COLOR{9}=[0.8 0.8 0.8];COLOR{10}=[0.9 0.9 0.9];
    for i = 1:nTr*length(CHN)
        tmp = TrialParams(TrialParams(:,2) == i,1);
        nT(:,i) = tmp(1:n_REP);
    end
    preStim = zeros(n_REP,nTr*length(CHN)); postStim = zeros(n_REP,nTr*length(CHN));
    chn = d(thisChn);
    sp = denoiseSpikes(sp{chn},chn);
    SPIKES = cell(1,nTr); spWaves = cell(1,nTr); psthBIN = [-200 200]; rate = zeros(nTr,diff(psthBIN)+1);
    TICKLABELS = cell(1,length(AMP)); LEGEND = cell(1,length(AMP));
    for i = 1:length(AMP)
        TICKLABELS{i} = AMP(i);
        LEGEND{i} = [num2str(AMP(i)) ' uA'];
    end
    TICKLABELS{1} = 'N.S.';
    LEGEND{1} = 'N.S.';
    %% PSTH Methodology for window adjustment
    if (mChn) && (PLOT)
        rateFig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;
        % Line at time point 0
        line([200 200],[0 1000*nTr],'Color','r','HandleVisibility','off');
        if (mChn)
            title([n1 ' | ' n2 ' | ' n3 ' | Chn ' num2str(thisChn)]);
        else
            title(['Chn: ' num2str(thisChn)]);
        end
        ylabel('Firing Rate (Sp/s) + Offset 100 S/s/Condition');
        xlabel('Time (ms)');
        xticks(0:10:400);
        xticklabels(-200:10:200);
        box on; grid on;
        spikeFig = figure('units','normalized','outerposition',[0 0 1 1]);
        spike_axes = tight_subplot(1,nTr,[.01 .015],[.075 .115],[.075 .02]);
    end
    for r = 1:nTr*length(CHN)
        meanWave = zeros(1,49); count = 0;
        for t = 1:n_REP
            try
                thisSp = sp(sp(:,1) >= (trig(nT(t,r))/30+psthBIN(1)) & sp(:,1) <= (trig(nT(t,r))/30+psthBIN(2)),:);
                SPIKES{t} = thisSp(:,1) - (trig(nT(t,r))/30+psthBIN(1));
            catch
                SPIKES{t} = 0;
            end
        end
        rate(r,:) = psth(SPIKES(:),psthBIN,1,500);
        baseRateMean = mean(rate(r,21:abs(psthBIN(1))-21));
        stStart = find(rate(r,abs(psthBIN(1)):end)>=baseRateMean*1.2,1) - 1;
        stStop = find(rate(r,abs(psthBIN(1))+stStart:end)<=baseRateMean*1.2,1) + stStart - 1;
        if stStop <= stimSt; stStop = 11; end
        if stStop >= 20; stStop = 11; end
        if isempty(stStop); stStop = 11; end
        stimulationBIN(r,1:2) = [stimSt, stStop];
        if r == 1
            stimulationBIN(r,1:2) = [stimSt,100];
        end
        stimulationBIN(r,1:2) = [2.25 11];
        if (mChn) && (PLOT)
            for t = 1:size(TrialParams,1)/size(nT,2)
                try
                    thisSp = sp(sp(:,1) >= (trig(nT(t,r))/30+stimulationBIN(r,1)) & sp(:,1) <= (trig(nT(t,r))/30+stimulationBIN(r,2)),:);
                catch
                    fprintf('It appears that some expected trigger lines were not recorded.\n');
                end
                spWaves{t} = thisSp(:,2:end);
            end
            figure(rateFig);
            tmp = rate(r,196:204);
            rate(r,197:203) = NaN;
            plot(rate(r,:) + (r-1)*100,...
                'Color',COLOR{r});
            plot(196:204,tmp + (r-1)*100,...
                'Color','b','HandleVisibility','off');
            figure(spikeFig);
            ax = spike_axes(r); hold(ax,'on');
            for t = 1:size(TrialParams,1)/size(nT,2)
                plot(ax,spWaves{t}');
                meanWave = meanWave + sum(spWaves{t},1);
                count = count + size(spWaves{t},1);
            end
            xticks(ax,0:15:49);
            xticklabels(ax,0:0.5:49/30);
            xlabel(ax,'Time (ms)');
            meanWave = meanWave ./ count;
            plot(ax,meanWave,'Color','k','LineWidth',4);
            text(ax,15,200,[num2str(count) ' spikes'],'FontSize',12);
            text(ax,15,180,[num2str(n_REP) ' trials'],'FontSize',12);
            text(ax,15,160,[num2str(diff(stimulationBIN(r,1:2))) ' ms window']);
            title(ax,[num2str(TICKLABELS{r}) ' uA']);
            if (r == 1)
                ylabel(ax,'Voltage (uV)');
            else
                yticks(ax,'');
            end
            ylim(ax,[-300 250]);
            xlim(ax,[0 49]);
            set(ax,'FontSize',14);
            beautifyPlot(14,ax);
        end
        sgtitle([n1 ' | ' n2 ' | ' n3 ' | Chn ' num2str(thisChn)],'FontSize',30);
    end
    if (mChn) && (PLOT)
        figure(rateFig);
        axis([-50-psthBIN(1) 50-psthBIN(1)+(nTr*(10)) -10 max(max(rate))+(nTr*(100))]);
        yticks(0:50:max(max(rate))+(nTr*(100)));
        set(gca,'FontSize',18);
        legend(LEGEND,'Location','northwest');
        beautifyPlot;
    end
    %% Logic
    for r = 1:n_REP %#ok<*BDSCI>
        for t = 1:nTr*length(CHN)
            try
                thisTrig = trig(nT(r,t))/30;
                preStim(r,t) = sum(sp(:,1) >= (thisTrig + baselineBIN(1)) & sp(:,1) <= (thisTrig + baselineBIN(2)));
                postStim(r,t) = sum(sp(:,1) >= (thisTrig + stimulationBIN(t,1)) & sp(:,1) <= (thisTrig + stimulationBIN(t,2)));
            catch
                preStim(r,t) = NaN;
                postStim(r,t) = NaN;
            end
        end
    end
    if (mChn) && (PLOT)
        close;
    end
    %% Statistics
    preMean = nanmean(preStim,1) ./ (diff(baselineBIN)/1e3);
    postMean = nanmean(postStim,1) ./ (diff(stimulationBIN,1,2)/1e3)';
    preSte = (nanstd(preStim,1) ./ (diff(baselineBIN)/1e3)) ./ sqrt(n_REP);
    postSte = (nanstd(postStim,1) ./ (diff(stimulationBIN,1,2)/1e3)') ./ sqrt(n_REP);
    %statTest = [(preStim ./ (diff(baselineBIN)/1e3)); (postStim ./ (diff(stimulationBIN,1,2))'); (suppStim ./ (diff(suppressionBIN,1,2))')];
    [preP,~,stats] = anova1(preStim ./ (diff(baselineBIN)/1e3),[],'off'); %#ok<*ASGLU>
    [st,m] = multcompare(stats,'alpha',0.05,'Display','off');
    [postP,~] = anova1(postStim,[],'off');
    % Confidence intervals
    ci(1) = mean(m(:,1)) + mean(st(:,5));
    ci(2) = mean(m(:,1)) + mean(st(:,3));
    %% Curve Fitting
    if length(CHN) == 1
        curveFit = @(prm,c)(prm(1).*c.^prm(2)./(c.^(prm(2).*prm(3))+prm(4).^(prm(2).*prm(3)))+prm(5));
        opts = optimset('Display','off');
        [bestFit,resnorm] = lsqcurvefit(curveFit,[0 1 1 1 postMean(1)],AMP(1:end)+2,postMean(1:end),[-inf 0 1 0 postMean(1)],[inf 10 1 AMP(end)+2 postMean(1)],opts);
    end
    %% Plotting
    figure(responseFig); hold on;
    xdata = zeros(1,length(AMP)*length(CHN));
    for iC = 1:length(CHN)
        for iA = 1:length(AMP)
            xdata(1,iA+(iC-1)*length(AMP)) = AMP(iA) + (iC-1)*(AMP(end)+10);
        end
    end
    if (mChn)
        dotSize = 100;
        ax = gca;
        hold(ax,'on');
    else
        INDEX = find(A' == thisChn);
        ax = s_axes(INDEX);
        hold(ax,'on');
        dotSize = 40;
    end
    if any(dead_channels == thisChn)
        xlim(ax,[0 1]);
        xticks(ax,'');
        ylim(ax,[0 1]);
        yticks(ax,'');
        line(ax,[0 1],[0 1],'Color','k');
        line(ax,[0 1],[1 0],'Color','k');
        beautifyPlot(30,ax);
        drawnow;
    else
        scatter(ax,xdata,preMean,dotSize,[0,0,0]);
        errorbar(ax,xdata,preMean,preSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
        scatter(ax,xdata,postMean,dotSize,[0,0,0],'filled','r');
        errorbar(ax,xdata,postMean,postSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
        %line(ax,[-2 nTr/length(CHN)+2],[ci(1) ci(1)],'Color','r','LineStyle','--','LineWidth',1);
        %line(ax,[-2 nTr/length(CHN)+2],[ci(2) ci(2)],'Color','r','LineStyle','--','LineWidth',1);
        if length(CHN) == 1
            XData = AMP(1)+2:0.0001:AMP(end)+2;
            CURVE = curveFit(bestFit,XData);
            plot(ax,XData(1:end)-2,CURVE,'Color','k','LineWidth',2);
        else
            %nAmp = length(AMP);
            %for iii = 1:length(CHN)-1
                %line(ax,[xdata(nAmp+(iii-1)*nAmp)+5 xdata(nAmp+(iii-1)*nAmp)+5],[0 max([preMean,postMean])],'Color','b');
            %end
        end
        if (mChn)
            xlabel(ax,'Current (uA)');
            xlim(ax,[-2 xdata(end)+1]);
            xticks(ax,xdata);
            xticklabels(ax,TICKLABELS);
            ylabel(ax,'Response (Sp/s)');
        else
            if find(XLABEL == INDEX,1)
                xlabel(ax,'Current (uA)');
                xlim(ax,[-2 xdata(end)+1]);
                xticks(ax,xdata);
                xticklabels(ax,TICKLABELS);
            else
                xlim(ax,[-2 xdata(end)+1]);
                xticks(ax,xdata);
                xticklabels(ax,'');
            end
            if find(YLABEL == INDEX,1)
                ylabel(ax,'Delta Response (dSp/s)');
            end
            title(ax,num2str(thisChn));
        end
        if max(postMean) == 0
            ylim(ax,[min(postMean)-3*max(postSte) 20]);
            yticks(ax,0:10:20);
        else
            ylim(ax,[min(postMean)-3*max(postSte) max(postMean)+3*max(postSte)]);
            y = ceil(max(postMean)+3*max(postSte)) / 3;
            if (y < 100)
                if (y < 10)
                    y = round(y);
                else
                    y = round(y,-1);
                end
            else
                y = round(y,-2);
            end
            if y == 0
                y = 1;
            end
            yticks(ax,round(min(postMean)-3*max(postSte)):y:round(max(postMean)+3*max(postSte)));
            yticklabels(ax,round(min(postMean)-3*max(postSte)):y:round(max(postMean)+3*max(postSte)));
        end
        if ~(mChn)
            if thisChn == stimChn
                set(ax,'XColor','r');
            end
        else
            ylim(ax,[0 100]);
            yticks(ax,0:25:100);
        end
        beautifyPlot(14,ax);
        drawnow;
    end
end
if ~(mChn)
    sgtitle([n1 ' | ' n2 ' | ' n3],'FontSize',24);
end
end