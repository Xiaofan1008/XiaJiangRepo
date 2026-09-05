%% This function will generate a list of stim-recording pairs with the
% following attributes:
% stimChn
% recChn
% signed distance
% cortical depth
% response threshold
%% Initial Variables
base = 'X:\Tim';
lOD = {'\RAT0007\RAT0007_datafile_008_181116_141638',1,1750;...
    '\RAT0007\RAT0007_datafile_013_181116_183226',2,1800;...
    '\RAT0015\estim_pen4_001_190325_171317',3,800;...
    '\RAT0016\estim_pen2_003_190328_151025',4,850;...
    '\RAT0017\estim_pen1_009_190401_193112',5,1300;...
    '\RAT0017\estim_pen1_010_190401_201325',5,1300;...
    '\RAT0017\estim_pen1_011_190401_204831',5,1300;...
    };
baseWN = [-175, -21]; % Baseline spike extraction window
psthWN = [-200 200]; % Window for calculating mean spike rates
nChn = 32; % Number of channels to look at
%% Initial Logic
stimT = nan(size(lOD,1),32,nChn); stimTrate = nan(size(lOD,1),nChn); stimMe = zeros(8,size(lOD,1)); baseMe = stimMe; baseTrate = stimTrate;
pF = nan(nChn,size(lOD,1));
CURVE = zeros(size(lOD,1),nChn,20001);
recData = cell(length(lOD(:,1)),5);
cchn = 0; stimSt = 2.25;
%% Extended Loops
for r = 1:size(lOD,1)
    cd([base lOD{r,1}]);
    %% Loading Advanced Variables
    sp = loadSpikes; tParams = loadTrialParams; trig = loadTrig(0); d = Depth; loadStimChn; DUR = loadDUR; AMP = loadAMP;
    nAmp = length(AMP);
    for sChn = 1:length(stimChn)
        nDUR = length(DUR);
        for dDUR = 1:nDUR
            nTrig = size(trig,2)./nAmp./length(stimChn);
            fParams = zeros(1,nTrig*nAmp);
            for a = 1:nAmp
                fParams(1+(a-1)*nTrig:(a)*nTrig) = cell2mat(tParams(cell2mat(tParams(:,2)) == (a + nAmp*(dDUR-1) + nAmp*nDUR*(sChn-1)),1));
            end
            fTrig = trig(fParams)';
            % trig is now ordered such that the first nTrigs are all no stim, et cetera
            baseSp = zeros(nAmp,nTrig); stimSp = zeros(nAmp,nTrig,nChn); stimWN = zeros(nAmp,2);
            for iChn = 1:nChn                
                baseWaveforms = cell(nAmp,nTrig); stimWaveforms = cell(nAmp,nTrig);
                if (r > 2 && iChn == 2) || (r > 2 && iChn == 9) || iChn == stimChn(sChn)
                    stimT(r,sChn,iChn) = NaN;
                    continue;
                end
                tChn = d(iChn);
                spt = denoiseSpikes(sp{tChn},tChn);
                %% Metric One
                % Threshold is defined as the lowest stimulation current that elicits a
                % mean response 50% greater than the no stim response
                SPIKES = cell(1,nTrig); rate = zeros(nAmp,diff(psthWN)+1);
                for a = 1:nAmp
                    thisTrig = fTrig(1+(a-1)*nTrig:(a)*nTrig) ./ 30;
                    for t = 1:size(fTrig,1)/nAmp
                        psthSp = spt(spt(:,1) >= (thisTrig(t) + psthWN(1)) & spt(:,1) <= (thisTrig(t) + psthWN(2)));
                        SPIKES{t} = psthSp - (thisTrig(t)+psthWN(1));
                    end
                    rate(a,:) = psth(SPIKES(:),psthWN,1,500,[],0);
                    baseRateMean = mean(rate(a,21:abs(psthWN(1))-21));
                    stStart = find(rate(a,abs(psthWN(1)):end)>=baseRateMean*1.2,1) - 1;
                    stStop = find(rate(a,abs(psthWN(1))+stStart:end)<=baseRateMean*1.2,1) + stStart - 1;
                    stimWN(a,1:2) = [stimSt, min([max([stStop,7]),10])];
                    if a == 1
                        stimWN(a,1:2) = [stimSt,8];
                    end
                    for t = 1:size(fTrig,1)/nAmp
                        baseSp(a,t,iChn) = sum(spt(:,1) >= (thisTrig(t) + baseWN(1)) & spt(:,1) <= (thisTrig(t) + baseWN(2)));
                        stimSp(a,t,iChn) = sum(spt(:,1) >= (thisTrig(t) + stimWN(a,1)) & spt(:,1) <= (thisTrig(t) + stimWN(a,2)));
                        baseWaveforms{a} = [baseWaveforms{a}; spt(spt(:,1) >= (thisTrig(t) + baseWN(1)) & spt(:,1) <= (thisTrig(t) + baseWN(2)),2:end)];
                        stimWaveforms{a} = [stimWaveforms{a}; spt(spt(:,1) >= (thisTrig(t) + stimWN(a,1)) & spt(:,1) <= (thisTrig(t) + stimWN(a,2)),2:end)];
                    end                    
                    stimMe(a,r) = (mean(stimSp(a,:,iChn)) ./ (diff(stimWN(a,:))/1e3));
                    baseMe(a,r) = (mean(baseSp(a,:,iChn)) ./ (diff(baseWN)/1e3));
                end
                %% Curve Fitting
                %AMP(2) = []; nAmp = length(AMP);
                %stimMe(2,:) = []; stimWN(2,:) = [];
                curveFit = @(prm,c)(prm(1).*c.^prm(2)./(c.^(prm(2).*prm(3))+prm(4).^(prm(2).*prm(3)))+prm(5));
                prm0 = lsqcurvefit(curveFit,[0 1 1 1 stimMe(1,r)],AMP(1:end)+2,stimMe(1:nAmp,r)',[-inf 0 1 0 stimMe(1,r)],[inf 10 1 AMP(end)+2 stimMe(1,r)]);
                DoF_constant = length(stimMe(1:nAmp,r)) - 1;
                DoF_sigmoid = length(stimMe(1:nAmp,r)) - 3;
                SS_constant = sum((mean(stimMe(1:nAmp,r)) - stimMe(1:nAmp,r)).^2);
                SS_sigmoid = sum((curveFit(prm0,AMP(1:end)+2) - stimMe(1:nAmp,r)').^2);
                F = ((SS_constant - SS_sigmoid) / (DoF_constant - DoF_sigmoid)) / (SS_sigmoid / DoF_constant);
                pF(iChn,r) = 1 - fcdf(F, DoF_constant - DoF_sigmoid, DoF_sigmoid);
                bestFit = prm0;
                for l = 0:0.0001:AMP(end)+2
                    CURVE(r,iChn,floor(l*1000)+1) = curveFit(bestFit,l);
                end
                %postP = anova1(squeeze(stimSp(1:nAmp,:,iChn))' ./ (diff(stimWN(1:end,:),1,2)/1e3)',[],'off'); %#ok<*ASGLU>
                %[~,~,stats] = anova1(baseSp(1:nAmp,:,iChn)' ./ (diff(baseWN(1:end,:),1,2)/1e3)',[],'off'); %#ok<*ASGLU>
                %[st,~] = multcompare(stats,'alpha',0.001,'Display','off');
                %ci = mean(st(:,5)); && (stimMe(nAmp,r) > ci) && (postP < 0.05)
                %% Logic
                if (max(CURVE(r,iChn,:)) > 0.5*max(stimMe(2:nAmp,r))) && (max(CURVE(r,iChn,:)) > 10)  && (prm0(1) > 0) && (pF(iChn,r) < 0.05)                                        
                    stimT(r,sChn,iChn) = find(CURVE(r,iChn,2e3:end) > (CURVE(r,iChn,1) + 0.5*diff([CURVE(r,iChn,1) max(CURVE(r,iChn,2e3:end))])),1) ./ 1000;
                    if stimT(r,sChn,iChn) < 0.1
                        %AMP = loadAMP; nAmp = length(AMP);
                        continue;
                    end
                    stimTrate(r,iChn) = CURVE(r,iChn,floor((stimT(r,sChn,iChn)) * 1000));
                    baseTrate(r,iChn) = stimMe(1,r);
                    recData{r,1} = [recData{r,1}; stimChn(sChn)];
                    recData{r,2} = [recData{r,2}; iChn];
                    recData{r,3} = [recData{r,3}; (stimChn(sChn) - iChn)*50];
                    recData{r,4} = [recData{r,4}; lOD{r,3} - (stimChn(sChn)-1)*50];
                    recData{r,5} = [recData{r,5}; stimT(r,sChn,iChn)];
                end
                %AMP = loadAMP; nAmp = length(AMP);
            end
        end
    end    
end
%% Saving Output
clear srPair
for r = 1:size(lOD,1)
    cd([base lOD{r,1}]);
    for i = 1:5
        tmp = recData{r,i};
        srPair(:,i) = tmp;
    end
    save('srPair.mat','srPair');
    clear srPair
end