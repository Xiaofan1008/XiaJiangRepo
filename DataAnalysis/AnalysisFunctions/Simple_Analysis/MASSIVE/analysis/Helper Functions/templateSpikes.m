%% Calculate a zero-condition spike template and r2t
function templateSpikes
sp = loadSpikes; nChn = size(sp,2);
trig = loadTrig(0);
TrialParams = loadTrialParams;
TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == 1));
trig = trig(TrialParams);
nTrig = length(trig);
spWN = [-500 500]; template = cell(1,nChn); r2t = cell(1,nChn); spt = cell(1,nChn);
for c = 1:nChn
    for n = 1:nTrig
        spt{c} = [spt{c}; sp{c}(sp{c}(:,1) > trig(n)/30+spWN(1) & sp{c}(:,1) < trig(n)/30+spWN(2),2:end)];
    end
    template{c} = mean(cell2mat(spt(c)),1);
    r2t{c} = mean((cell2mat(spt(c)) - template{c}).^2);
end
end