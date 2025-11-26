% Load in the Intan data
generateIntanData;
FS = 30000;
% Load in the digital lines
generateTimeStamps;
% Load in the datafile and select the Trial ID for analysis
StimParams = [];AMP = [];DUR = [];CHN = [];
load([filepath, datafile],'AMP','DUR','n_REP','CHN','StimParams','n_Trials','TrialParams');
if (n_Trials ~= nStamps)
    disp('Warning: Digital Lines do not match expected number of trials');
    %return;
end
uniqueTrials = length(AMP)*length(DUR)*length(CHN);
% Build out the filters
if ~(NN_order)
    depth = Depth(1)+1;
else
    depth = Depth(2);
end
if (NN_I)
    depth = Depth(3);
end
spExist = dir([filepath '\*.sp.mat']);
if isempty(spExist)
    disp('Calculating spike thresholds . . .');
    generateFiltersandThresholds;
end
X = BIN(1):BIN(2);
