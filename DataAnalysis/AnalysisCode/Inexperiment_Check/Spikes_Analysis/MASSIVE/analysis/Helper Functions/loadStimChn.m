tparams = dir('*_exp_datafile_*.mat');
if isempty(tparams)
    TrialParams = [];
    return
end
tparams = tparams.name;
load(tparams,'TrialParams');
stimChn = unique(cell2mat(TrialParams(2:end,3)));
stimChn=stimChn(stimChn~=0);
load(tparams,'CHN');
CHN=CHN;

%%original code
% stimChn = dir('*_exp_datafile_*.mat');
% stimChn = stimChn.name;
% load(stimChn,'CHN');
% stimChn = CHN;