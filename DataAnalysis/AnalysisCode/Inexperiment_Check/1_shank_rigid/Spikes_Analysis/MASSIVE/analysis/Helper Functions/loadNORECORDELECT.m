tparams = dir('*_exp_datafile_*.mat');
if isempty(tparams)
    TrialParams = [];
    return
end
tparams = tparams.name;
try
load(tparams,'NORECORDELECT');
NORECORDELECT = NORECORDELECT(NORECORDELECT>0);
catch
    NORECORDELECT=0;
end

