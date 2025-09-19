function E_Mapnumber=loadMapNum
tparams = dir('*_exp_datafile_*.mat');
if isempty(tparams)
    TrialParams = [];
    return
end
tparams = tparams.name;
load(tparams,'E_Mapnumber');
end

