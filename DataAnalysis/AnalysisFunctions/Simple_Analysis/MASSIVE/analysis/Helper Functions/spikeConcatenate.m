function spikeConcatenate
%% Based off stimConcatenate - draws stim-associated spikes from each channel
%% Variables
spWIN = [4 16];
nChn = 1:32;
tID = 3;
%% Loading in data
d = Depth(3); trig = []; TrialParams = []; sp = [];
loadTrig;
loadSpikes;
loadTrialParams;
if ~isempty(TrialParams)
    TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) >= tID));
else
    TrialParams = 1:length(trig);
end
trig = trig(TrialParams);
nTrig = length(trig);
trig = cast(trig,'int32');
setupGraphsBasic;
n = zeros(1,length(nChn));
%% Initialise data matrices
figure; hold on;
for t = 1:nTrig
    for c = nChn
        INDEX = find(A' == c);
        subplot(ROW,COL,INDEX); hold on;
        title(['CHN: ' num2str(c)]);
        thisSp = sp{d(c)};
        thisSp = thisSp(thisSp(:,1) >= trig(t)+spWIN(1) & thisSp(:,1) <= trig(t)+spWIN(2),:);
        thisSp = denoiseSpikes(thisSp);
        n(c) = n(c) + size(thisSp,1);
        plot(thisSp(:,2:end)');
    end
end
end