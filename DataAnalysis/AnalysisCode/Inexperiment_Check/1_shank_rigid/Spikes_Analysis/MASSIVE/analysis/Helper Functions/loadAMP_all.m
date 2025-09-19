tparams = dir('*_exp_datafile_*.mat');
if isempty(tparams)
    TrialParams = [];
    return
end
tparams = tparams.name;
load(tparams,'StimParams');
AMP_all = unique(cell2mat(StimParams(2:end,16)));
%%original code
% stimChn = dir('*_exp_datafile_*.mat');
% stimChn = stimChn.name;
% load(stimChn,'CHN');
% stimChn = CHN;