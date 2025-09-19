function prettyPhase
% Generates nice polar plots for use in explaining phase analysis
dbstop if error
%% Setup
% Base directory
base = 'X:\Tim';
% Recordings included in the global analysis
lOD = {
    '\RAT0017\estim_pen1_009_190401_193112',0,-1;...
    '\RAT0017\estim_pen1_010_190401_201325',0,-1;...
    '\RAT0017\estim_pen1_011_190401_204831',0,-1;...
    };
listOfChannels = {3:8,10:17,18:24};
%% Basic Variables
BAND = [4,5,6,10,20,80]; % The frequency bands to include in analysis
phaseOFFSET = [-251, -201, -251, -101, -51, -26]; % The temporal offset from stimulation to measure phase
phaseBins = [360, 6]; % Number of phase-bins to use
SW = [2.25 11]; % Spike extraction window !! IMPORTANT !!
%peakS = 1; % Normalise for peak firing rate
fakeTrials = 10000;
%% Set up channel loops
depth = Depth;
nChn = size(listOfChannels,2);
for nR = 1:size(lOD,1)
    cd([base lOD{nR}]);
    %% Calculate the phase array
    loadStimChn; trig = loadTrig(0); AMP = loadAMP; nAmp = length(AMP); TrialParams = loadTrialParams; sp = loadSpikes;
    lfp = loadLFP(depth(stimChn));
    phase = generatePhaseVector(lfp,BAND);
    for n = 1:nChn        
        % Generate a response array for each channel here
        cChn = size(listOfChannels{n},2);
        thisMean = zeros(length(BAND),phaseBins(1),nAmp);
        thisSEM = zeros(length(BAND),phaseBins(1),nAmp);
        fakeMean = zeros(length(BAND),phaseBins(1),nAmp);
        fakeSEM = zeros(length(BAND),phaseBins(1),nAmp);
        for c = 1:cChn
            iChn = depth(listOfChannels{n}(c));            
            spt = denoiseSpikes(sp{iChn},iChn);
            spt = spt(:,1);
            tphase = cell(nAmp,length(BAND));
            binCount = zeros(length(trig)/nAmp,length(BAND),nAmp);
            for a = 1:nAmp
                thisTrial = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == a));
                tTrig = cast(trig(thisTrial),'int64');
                for b = 1:length(BAND)
                    tphase{a,b}(:,1) = phase{b}(cast(trig(:)./30+phaseOFFSET(b),'int64'));
                    tphase{a,b}(tphase{a,b}(:,1) < 0,1) = tphase{a,b}(tphase{a,b}(:,1) < 0,1) + 360;
                    binCount(:,b,a) = ceil((tphase{a,b}(thisTrial,1))/(360/phaseBins(1)));
                end
                t = zeros(length(tTrig),1);
                fT = zeros(length(tTrig),fakeTrials);
                for tr = 1:length(tTrig)
                    % Real Spikes
                    chk = spt >= tTrig(tr)/30+SW(1) & spt <= tTrig(tr)/30+SW(2);
                    t(tr,1) = sum(chk);
                end
                for fTrials = 1:fakeTrials
                    fTrig = randperm(length(tTrig));
                    % Fake Spikes
                    fT(:,fTrials) = t(fTrig,1);
                end
                fakeData = nanmean(fT,2);
                for b = 1:length(BAND)
                    for p = 1:phaseBins(1)
                        if isempty(t(squeeze(binCount(:,b,a)) == p)); else
                            thisMean(b,p,a) = thisMean(b,p,a) + nanmean(t(squeeze(binCount(:,b,a)) == p));
                            thisSEM(b,p,a) = thisSEM(b,p,a) + nanstd(t(squeeze(binCount(:,b,a)) == p)) ./ sqrt(length(trig)/nAmp);
                            fakeMean(b,p,a) = fakeMean(b,p,a) + nanmean(fakeData(squeeze(binCount(:,b,a)) == p));
                            fakeSEM(b,p,a) = fakeSEM(b,p,a) + nanstd(fakeData(squeeze(binCount(:,b,a)) == p)) ./ sqrt(length(trig)/nAmp);                           
                        end
                    end
                end
            end
        end
        % Finish averaging thisMean together
        thisMean = thisMean ./ cChn;
        thisSEM = thisSEM ./ cChn;
        fakeMean = fakeMean ./ cChn;
        fakeSEM = fakeSEM ./ cChn;  
        for a = 1:nAmp
            for b = 1:length(BAND)                
                for p = 1:phaseBins(1)
                    if isnan(thisMean(b,p,a))
                        thisMean(b,p,a) = 0;
                    end
                    if isnan(thisSEM(b,p,a))
                        thisSEM(b,p,a) = 0;
                    end
                    if isnan(fakeMean(b,p,a))
                        fakeMean(b,p,a) = 0;
                    end
                    if isnan(fakeSEM(b,p,a))
                        fakeSEM(b,p,a) = 0;
                    end
                    if (sum(squeeze(binCount(:,b,a)) == p) == 0)
                        thisMean(b,p,a) = NaN;
                        thisSEM(b,p,a) = NaN;
                        fakeMean(b,p,a) = NaN;
                        fakeSEM(b,p,a) = NaN;
                    end
                end
            end
        end
        finalMean = zeros(length(BAND),phaseBins(2),nAmp); finalSEM = zeros(length(BAND),phaseBins(2),nAmp);
        fakeFinalMean = zeros(length(BAND),phaseBins(2),nAmp); fakeFinalSEM = zeros(length(BAND),phaseBins(2),nAmp);
        centreM = zeros(nAmp,length(BAND));
        % Calculate a reference point for phasebins
        tmp = 1:360;
        for a = 1:nAmp
            for b = 1:length(BAND)
%                 wa(1,:) = tmp(~isnan(thisMean(b,:,a)));
%                 wa(2,:) = thisMean(b,~isnan(thisMean(b,:,a)),a) ./ min(thisMean(b,(thisMean(b,:,a)>0),a));
%                 wa(2,:) = cast(wa(2,:) .* 6,'int32');
%                 lwa = size(wa,2);
%                 for l = 1:lwa
%                     if wa(2,l) ~= 0
%                         wa = [wa, repmat([wa(1,l);wa(2,l)],1,wa(2,l)-1)];
%                     end
%                 end
%                 circ_plot(circ_ang2rad(wa(1,:)'),'hist',[],180,true,true,'linewidth',3,'color','k');
%                 H = findall(gca,'Type','text');
%                 delete(H);
%                 clear wa
                centreM(a,b) = circ_rad2ang(circ_mean(circ_ang2rad(tmp(~isnan(thisMean(b,:,a))))',thisMean(b,~isnan(thisMean(b,:,a)),a)')) + 360;
            end
        end
        for a = 1:nAmp
            for b = 1:length(BAND)                                
                m = round(centreM(a,b));
                for p = 1:phaseBins(2)
                    data = repmat(thisMean(b,:,a),[1,3]);
                    finalMean(b,p,a) = nanmean(data(m-(180/phaseBins(2))+(p-1)*(360/phaseBins(2)):m+(180/phaseBins(2))+(p-1)*(360/phaseBins(2))));
                    data = repmat(thisSEM(b,:,a),[1,3]);
                    finalSEM(b,p,a) = nanmean(data(m-(180/phaseBins(2))+(p-1)*(360/phaseBins(2)):m+(180/phaseBins(2))+(p-1)*(360/phaseBins(2))));
                    data = repmat(fakeMean(b,:,a),[1,3]);
                    fakeFinalMean(b,p,a) = nanmean(data(m-(180/phaseBins(2))+(p-1)*(360/phaseBins(2)):m+(180/phaseBins(2))+(p-1)*(360/phaseBins(2))));
                    data = repmat(fakeSEM(b,:,a),[1,3]);
                    fakeFinalSEM(b,p,a) = nanmean(data(m-(180/phaseBins(2))+(p-1)*(360/phaseBins(2)):m+(180/phaseBins(2))+(p-1)*(360/phaseBins(2))));
                end
            end
        end
        % At this point, we have enough to plot!
        for b = 1:length(BAND)
            for a = 1:nAmp
                MAX = max(finalMean(b,:,a)) + max(finalSEM(b,:,a));
                figure; hold on;
                for i = 1:phaseBins(2)
                    bar(i,squeeze(finalMean(b,i,a)),'FaceColor','w','EdgeColor','k','LineWidth',2);
                    errorbar(i,squeeze(finalMean(b,i,a)),squeeze(finalSEM(b,i,a)),'.','Color','k');                    
                end
                ylim([0 MAX*1.2]);
                xticks(1:phaseBins(2));
                if (a==1); ylabel([num2str(BAND(b)) ' Hz']); end
                if (b==length(BAND)); xlabel('Phasebins'); end
                beautifyPlot;
                line(1:phaseBins(2),squeeze(fakeFinalMean(b,:,a)),'Color','r','LineWidth',2);
                errorbar(1:phaseBins(2),squeeze(fakeFinalMean(b,:,a)),squeeze(fakeFinalSEM(b,:,a)),'Color','r','LineWidth',2)                               
            end
            drawnow;
        end
    end
end
end