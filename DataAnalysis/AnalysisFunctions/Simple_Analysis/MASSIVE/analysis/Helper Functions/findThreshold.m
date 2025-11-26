function [thresh, pF, r2, aBase, mStim, tStim, sat] = findThreshold(chn,PLOT,ID,nS)
dbstop if error
if nargin<2
    PLOT = 0; ID = 0;
elseif nargin<3
    ID = 1;
elseif nargin<4
    nS = 0;    
end
%% Basic Variables
baseWN = [-175, -21]; % Baseline spike extraction window
stimWN = [2.25 11]; % Window for calculating mean spike rates
% Load in variablesd = Depth; chn = d(chn);
AMP = loadAMP; nAmp = length(AMP);
DUR = loadDUR; nDur = length(DUR);
loadStimChn; nSti = length(stimChn);
sp = loadSpikes(chn); sp = denoiseSpikes(sp,chn);
tParams = loadTrialParams;
trig = loadTrig(0); nRep = length(trig)/(nAmp*nDur*nSti);
baseSp = zeros(nAmp,nRep); stimSp = zeros(nAmp,nRep);
stimMe = zeros(1,nAmp); baseMe = stimMe;
stimSte = zeros(1,nAmp); baseSte = stimSte;
CURVE = zeros(1,20001);
thresh = zeros(nDur,nSti);
%% Logic
if nS ~= 0
    stimSTART = nS; stimSTOP = nS;
else
    stimSTART = 1; stimSTOP = nSti;
end
for nS = stimSTART:stimSTOP
    AMP = loadAMP; nAmp = length(AMP); loadNREP;
    for nD = 1:nDur
        nTrig = size(trig,2)./nAmp./length(stimChn);
        if nTrig > n_REP
            nTrig = n_REP;
        end
        fParams = zeros(1,nTrig*nAmp);        
        for nA = 1:nAmp
            % Produces an ordered vector of trials
            tmp = cell2mat(tParams(cell2mat(tParams(:,2)) == (nA + nAmp*(nD-1) + nAmp*nDur*(nS-1)),1));
            fParams(1+(nA-1)*nTrig:(nA)*nTrig) = tmp(1:nTrig);
        end
        fTrig = trig(fParams)'; % An ordered list of trigger lines
        for nA = 1:nAmp
            if AMP(nA) == 0
                baseSp(nA,:) = NaN;
                stimSp(nA,:) = NaN;
                continue;
            end
            thisTrig = fTrig(1+(nA-1)*nTrig:(nA)*nTrig)./30;
            for t = 1:size(fTrig,1)/nAmp
                baseSp(nA,t) = sum(sp(:,1) >= (thisTrig(t) + baseWN(1)) & sp(:,1) <= (thisTrig(t) + baseWN(2)));
                stimSp(nA,t) = sum(sp(:,1) >= (thisTrig(t) + stimWN(1)) & sp(:,1) <= (thisTrig(t) + stimWN(2)));
            end
            baseMe(nA) = (nanmean(baseSp(nA,:)) ./ (diff(baseWN)/1e3));
            stimMe(nA) = (nanmean(stimSp(nA,:)) ./ (diff(stimWN)/1e3));
            baseSte(nA) = (nanstd(baseSp(nA,:)) ./ (diff(baseWN)/1e3)) ./ sqrt(t);
            stimSte(nA) = (nanstd(stimSp(nA,:)) ./ (diff(stimWN)/1e3)) ./ sqrt(t);
        end
        if any(AMP == 0)
            baseMe(AMP == 0) = [];
            stimMe(AMP == 0) = [];
            AMP(AMP == 0) = [];
            AMP(1) = 0;
            nAmp = length(AMP);
        elseif any(AMP == -1)
            AMP(AMP == -1) = 0;
        end
        % Curve Fitting
        curveFit = @(prm,c)(prm(1).*c.^prm(2)./(c.^(prm(2).*prm(3))+prm(4).^(prm(2).*prm(3)))+prm(5));
        opts = optimset('Display','off');
        prm0 = lsqcurvefit(curveFit,[0 1 1 1 stimMe(1)],AMP(1:end)+1,stimMe(1:end),[-inf 0 1 0 stimMe(1)],[inf 10 1 AMP(end)+1 stimMe(1)],opts);
        for l = 0.0001:0.0001:AMP(end)
            CURVE(floor(l*1000)+1) = curveFit(prm0,l+1);
        end
        if length(CURVE) > floor(l*1000)+1
            CURVE(floor(l*1000)+2:end) = [];
        end
        DoF_constant = length(stimMe(1:nAmp)) - 1;
        DoF_sigmoid = length(stimMe(1:nAmp)) - 3;
        SS_constant = sum((mean(stimMe(1:nAmp)) - stimMe(1:nAmp)).^2);
        SS_sigmoid = sum((curveFit(prm0,AMP(1:end)+1) - stimMe(1:nAmp)).^2);
        F = ((SS_constant - SS_sigmoid) / (DoF_constant - DoF_sigmoid)) / (SS_sigmoid / DoF_constant);
        pF(nD,nS) = 1 - fcdf(F, DoF_constant - DoF_sigmoid, DoF_sigmoid); %#ok<*AGROW>
        sse = sum((stimMe(1:end)-CURVE(1000*(AMP)+1)).^2);
        sst = sum((stimMe(1:end)-mean(CURVE(:))).^2);
        r2(nD,nS) = 1 - sse/sst;
        try
            thresh(nD,nS) = find(CURVE > (CURVE(1) + 0.5*diff([CURVE(1) max(CURVE)])),1) ./ 1000;
            if thresh(nD,nS) <= 0
                thresh(nD,nS) = NaN;
            end
        catch
            thresh(nD,nS) = NaN;
        end
        aBase(nD,nS) = nanmean(baseMe);
        mStim(nD,nS) = max(stimMe);
        if isnan(thresh(nD,nS))
            tStim(nD,nS) = NaN;
        else
            tStim(nD,nS) = CURVE(round(thresh(nD,nS)*1000));
        end
        sat(nD,nS) = 0;
        if (stimMe(end) <= stimMe(end-1)*1.25)
            sat(nD,nS) = 1;
        end
        if (PLOT)
            if (nS == ID)
                figure; hold on;                
                scatter(AMP*1000,stimMe,100,[0,0,0],'filled','r');
                errorbar(AMP*1000,stimMe,stimSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
                scatter(AMP*1000,baseMe,100,[0,0,0]);
                errorbar(AMP*1000,baseMe,baseSte,'Color',[0,0,0],'LineStyle','none','LineWidth',1,'HandleVisibility','off');
                YLIM = ylim;    
                plot(CURVE,'Color','k','LineWidth',2);
                if ~isnan(thresh(nD,nS))
                    line([thresh(nD,nS)*1000 thresh(nD,nS)*1000],[0 1000],'Color','r','LineWidth',2);
                    line([(AMP(1)-2)*1000 (AMP(end)+2)*1000],[CURVE(round(thresh(nD,nS)*1000)) CURVE(round(thresh(nD,nS)*1000))],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','--','HandleVisibility','off');
                end
                ylim(YLIM*1.1); xlim([-1000 (AMP(end)+1)*1000]);
                xticks(0:1000:AMP(end)*1000);
                xticklabels(0:1:AMP(end));
                beautifyPlot(14);
                legend({'Response to Stimulation','Baseline Firing Rate',['r^2 = ' num2str(r2(nD,nS),3)],['Threshold = ' num2str(thresh(nD,nS),3) ' \muA']},'Location','northwest');
                ylabel('Firing Rate (Sp/s)');
                xlabel('Stimulation Current (\muA)')
            end
        end
    end
end
end