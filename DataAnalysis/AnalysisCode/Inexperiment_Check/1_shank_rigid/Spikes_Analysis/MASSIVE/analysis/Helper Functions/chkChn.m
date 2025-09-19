%% This script defines a set of selection criteria for excluding channels from analysis
%% Rule Three
% The stimulation mean must have a positive gradient
lobf = interpolate(postMean,1);
chk = diff([lobf(1),lobf(end)]);
if chk <= 0
    return;
end
%% Rule Four
% The maximum stimulation mean must exceed threshold
chk = max(postMean);
if chk <= ci(1)
    return;
end
%% Rule Five
% At least one stimulation point must be significantly above the
% baseline
chk = find(H == 1);
chk = sum(postMean(chk)>=ci(1));
if chk <= 0
    return;
end
%% Made it through the rules? Include the channel
collapse(thisChn,2) = 1;
collapse(thisChn,3) = min(max(find(postMean > ci(1),1),3),8);
collapse(thisChn,4) = postMean(collapse(thisChn,3)) .* diff(stimulationBIN(collapse(thisChn,3),1:2))/1e3;