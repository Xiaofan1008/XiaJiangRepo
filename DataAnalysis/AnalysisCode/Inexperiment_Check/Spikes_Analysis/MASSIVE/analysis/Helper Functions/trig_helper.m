function output = trig_helper(BAND,loCHN,trials)
%% This function creates a data-structure that holds a variety of useful information about trigger lines
trig = loadTrig(0); TrialParams = loadTrialParams; loadNTrials; depth = Depth; loadStimChn; sp = loadSpikes; loadNREP;
nTrial = length(trials);
% Initial Structures
listOfChannels = loCHN;
BN = [2.25 11]; phaseBins = [360 6];
lfp = loadLFP(depth(stimChn));
nTrig = length(trig);
if nTrig > n_Trials
    nTrig = n_Trials;
end
phase = generatePhaseVector(lfp,BAND(1));
tphase = phase{1}(cast(trig(:)./30+BAND(2),'int64'));
tphase(tphase <= 0) = tphase(tphase <= 0) + 360;
% Generate a response array for each channel here
cChn = size(listOfChannels,2);
thisSpikes = zeros(phaseBins(1),nTrial);
binCount = zeros(n_REP,nTrial);
for c = 1:cChn
    iChn = depth(listOfChannels(c));
    spt = denoiseSpikes(sp{iChn},iChn);
    spt = spt(:,1);
    for a = 1:nTrial
        thisTrial = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == trials(a)));
        thisTrial = thisTrial(1:n_REP);
        tTrig = trig(thisTrial)./30;        
        binCount(:,a) = ceil((tphase(thisTrial))/(360/phaseBins(1)));
        t = zeros(length(tTrig),1);
        for tr = 1:n_REP
            % Real Spikes
            chk = spt >= tTrig(tr)+BN(1) & spt <= tTrig(tr)+BN(2);
            t(tr,1) = sum(chk);
            %% Each phase bin
            %phSPIKES{ph(tr,a),a} = [phSPIKES{ph(tr,a),a}; spt(spt >= sTrig(tr)+BN(1) & spt <= sTrig(tr)+BN(2)) - sTrig(tr)];
        end
        for p = 1:phaseBins(1)
            thisSpikes(p,a) = thisSpikes(p,a) + nanmean(t(binCount(:,a) == p));
        end
    end
end
tmp = 1:phaseBins(1); centreM = zeros(1,nTrial);
for a = 1:nTrial
    for i = 1:phaseBins(1)
        thisSpikes(i,a) = thisSpikes(i,a) ./ cChn;
        if sum(binCount(:,a) == i) == 0
            thisSpikes(i,a) = NaN;
        end
    end
    centreM(a) = round(circ_rad2ang(circ_mean(circ_ang2rad(tmp(~isnan(thisSpikes(:,a)))'),thisSpikes(~isnan(thisSpikes(:,a)),a))));
    if centreM(a) < 0; centreM(a) = centreM(a) + 360; end     
end
output = zeros(nTrig,5);
for i = 1:nTrig
    output(i,1) = i;
    output(i,2) = TrialParams{i,2};
    output(i,3) = tphase(i);
    for a = 1:nTrial
        output(output(:,2) == trials(a),4) = centreM(a);
    end
end
bn = cell(phaseBins(2),max(output(:,2)));
for p = 1:phaseBins(2)
    for a = 1:max(output(:,2))
        bn{p,a} = unique(output(output(:,2)==a,4));
        if p == 1; bn{p,a} = [bn{p,a}-30, bn{p,a}+30];
        else; bn{p,a} = [bn{p-1,a}(2) bn{p-1,a}(2)+60];
        end
        if bn{p,a}(1) < 0; bn{p,a}(1) = bn{p,a}(1) + 360;
        elseif bn{p,a}(1) > 360; bn{p,a}(1) = bn{p,a}(1) - 360; end
        if bn{p,a}(2) < 0; bn{p,a}(2) = bn{p,a}(2) + 360; 
        elseif bn{p,a}(2) > 360; bn{p,a}(2) = bn{p,a}(2) - 360; end
    end
end
for n = 1:nTrig
    ind = cell2mat(bn(:,output(n,2)));
    for p = 1:phaseBins(2)
        if output(n,3) > ind(p,1) && output(n,3) < ind(p,2)
            output(n,5) = p;
            break;
        end
        if output(n,3) > ind(p,1) && output(n,3) > ind(p,2) && ind(p,1) > ind(p,2)
            output(n,5) = p;
            break;
        end
        if output(n,3) < ind(p,1) && output(n,3) < ind(p,2) && ind(p,1) > ind(p,2)
            output(n,5) = p;
            break;
        end
    end
end
end