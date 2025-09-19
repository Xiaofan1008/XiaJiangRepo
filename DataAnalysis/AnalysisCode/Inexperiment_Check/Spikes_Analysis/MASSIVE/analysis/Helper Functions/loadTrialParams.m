function TrialParams = loadTrialParams
tparams = dir('*_exp_datafile_*.mat');
if isempty(tparams)
    TrialParams = [];
    return
end
tparams = tparams.name;
load(tparams,'TrialParams');
TrialParams = TrialParams(2:end,:);
end