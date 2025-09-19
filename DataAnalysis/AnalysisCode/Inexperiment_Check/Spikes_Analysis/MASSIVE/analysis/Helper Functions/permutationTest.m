function [stats] = permutationTest(data,n,AMP,thRate)
% Runs a simple permutation test over the data structure
% data should be a cell array where each cell contains a vector of trials
numAmp = size(data,1);
numPh = size(data,2);
minTrialNum = min(cellfun('size',data,1),[],2);
curveFit = @(prm,c)(prm(1).*c.^prm(2)./(c.^(prm(2).*prm(3))+prm(4).^(prm(2).*prm(3)))+prm(5));
opts = optimset('Display','off');
prm0 = lsqcurvefit(curveFit,[0 1 1 1 mean(thRate(1,:))],AMP(1:end),mean(thRate(1:numAmp,:),2)',[max(mean(thRate(1:numAmp,:),2)) 0 1 0 mean(thRate(1,:))],[max(mean(thRate(1:numAmp,:),2)) 10 1 AMP(end) mean(thRate(1,:))],opts);
realCurve = zeros(numPh,floor(AMP(end)*1000)+1);
realThresh = zeros(1,numPh);
for nP = 1:numPh
    % Calculate the threshold for this trial ID
    prm1 = lsqcurvefit(curveFit,[prm0(1), prm0(2), prm0(3), prm0(4), prm0(5)],AMP(1:end),thRate(1:numAmp,nP)',[prm0(1), prm0(2), prm0(3), 0, prm0(5)],[prm0(1), prm0(2), prm0(3), AMP(end), prm0(5)],opts);
    for l = 0.0001:0.0001:AMP(end)
        realCurve(nP,floor(l*1000)+1) = curveFit(prm1,l);
    end    
end
MIN = mean(realCurve(:,1)); MAX = mean(realCurve(:,end));
for nP = 1:numPh
    try
        realThresh(nP) = find(realCurve(nP,:) > (MIN + 0.5*diff([MIN MAX])),1) ./ 1000;
        if realThresh(nP) <= 0 || isempty(realThresh(nP))
            realThresh(nP) = NaN;
        end
    catch
        realThresh(nP) = NaN;
    end
end
% Calculate the modulation index
fO = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf (2*pi)/numPh -inf],'Upper',[inf (2*pi)/numPh inf]);
xval = pi/3:pi/3:(3*numPh)*pi/3;
fitRealThresh = realThresh - mean(realThresh);
fitRealThresh = repmat(fitRealThresh,[1,3]);
[f, ~] = fit(xval',fitRealThresh','sin1',fO);
fy = f(pi/3:0.01:4*pi+pi/3) + mean(realThresh);
realMI = diff([min(fy) max(fy)]);
% Perform the permutation test
thresh = zeros(n,numPh); pMI = zeros(n,1);
for nT = 1:n
    pData = cell(numAmp,numPh);
    for nA = 1:numAmp
        % Permute within each amplitude to create a phase-invariant
        % distribution
        thisData = [];
        for nP = 1:numPh
            thisData = [thisData; data{nA,nP}(1:ceil(minTrialNum(nA)*0.85))]; %#ok<AGROW>
        end
        for nP = 1:numPh
            pData{nA,nP} = sum(thisData(randperm(size(thisData,1),ceil(minTrialNum(nA)*0.85)))) ./ (ceil(minTrialNum(nA)*0.85).*(diff([2.25,11])/1e3));            
        end        
    end
    CURVE = zeros(numPh,floor(AMP(end)*1000)+1);
    for nP = 1:numPh
        % Convert the pData structure into curves
        prm1 = lsqcurvefit(curveFit,[prm0(1), prm0(2), prm0(3), prm0(4), prm0(5)],AMP(1:end),cell2mat(pData(:,nP))',[prm0(1), prm0(2), prm0(3), 0, prm0(5)],[prm0(1), prm0(2), prm0(3), AMP(end), prm0(5)],opts);        
        for l = 0.0001:0.0001:AMP(end)
            CURVE(nP,floor(l*1000)+1) = curveFit(prm1,l);
        end        
    end
    MIN = mean(CURVE(:,1)); MAX = mean(CURVE(:,end));
    % Convert each CURVE into a threshold
    for nP = 1:numPh                
        try
            thresh(nT,nP) = find(CURVE(nP,:) > (MIN + 0.5*diff([MIN MAX])),1) ./ 1000;
            if thresh(nT,nP) <= 0 || isempty(thresh(nT,nP))
                thresh(nT,nP) = NaN;
            end
        catch
            thresh(nT,nP) = NaN;
        end
    end
    % Convert each threshold into a sinusoid    
    thisThresh = thresh(nT,:) - mean(thresh(nT,:));
    thisThresh = repmat(thisThresh,[1,3]);    
    [f, ~] = fit(xval',thisThresh','sin1',fO);
    fy = f(pi/3:0.01:4*pi+pi/3) + mean(thresh(nT,:));
    pMI(nT) = diff([min(fy) max(fy)]);
end
% Statistics on the threshold distribution
alpha = 0.05; % Significance level
pMI = sort(pMI);
figure; hold on;
histogram(pMI,0:0.1:4,'FaceColor',[0.5 0.5 0.5]);
YLIM = ylim;
xlim([0 max(pMI)+0.1]);
text(pMI(n-ceil(n*alpha))-0.01,YLIM(2)-1,[num2str((1-alpha)*100) '% CI'],'FontSize',18,'HorizontalAlignment','right');
line([pMI(n-ceil(n*alpha)) pMI(n-ceil(n*alpha))],[YLIM(1) YLIM(2)],'Color',[0.5 0.5 0.5],'LineWidth',2','LineStyle','--');
text(realMI+0.01,YLIM(2)-1,'Real Data','FontSize',18,'Color','r');
line([realMI realMI],[YLIM(1) YLIM(2)],'Color','r','LineStyle','--','LineWidth',2);
x = n - find(pMI >= realMI,1);
p = x/n;
if p < alpha
    H = 1;
    title('We reject the null hypothesis that there is no phase effect on threshold');
else
    H = 0;
    title('We accept the null hypothesis that there is no phase effect on threshold');
end
ylabel('Number of Permutations (#)');
xlabel('Modulation of Threshold due to Phase (\muA)');
beautifyPlot(18);
stats = [H,p];
figure; hold on;
subplot(1,2,1); hold on;
for nP = 1:numPh
    bar(nP,mean(thresh(:,nP)));
    errorbar(nP,mean(thresh(:,nP)),ste(thresh(:,nP)),'.');
end
ylim([0 8]);
xlabel('Phasebins');
ylabel('Threshold (\muA)');
title('Permutation of Spiking Response');
beautifyPlot(18);
subplot(1,2,2); hold on;
for nP = 1:numPh
    bar(nP,realThresh(:,nP));    
end
ylim([0 8]);
xlabel('Phasebins');
ylabel('Threshold (\muA)');
title('Actual Spiking Response');
beautifyPlot(18);
end