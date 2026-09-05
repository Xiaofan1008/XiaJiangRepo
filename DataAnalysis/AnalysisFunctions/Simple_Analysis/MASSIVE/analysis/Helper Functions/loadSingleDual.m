tparams = dir('*_exp_datafile_*.mat');
if isempty(tparams)
    TrialParams = [];
    return
end
tparams = tparams.name;
load(tparams,'singanddual');


%%original code
% stimChn = dir('*_exp_datafile_*.mat');
% stimChn = stimChn.name;
% load(stimChn,'CHN');
% stimChn = CHN;