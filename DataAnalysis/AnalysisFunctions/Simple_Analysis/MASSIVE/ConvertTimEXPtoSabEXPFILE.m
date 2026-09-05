Data=load('RAT0025_exp_datafile_002.mat');
StimParams=Data.StimParams;
TrialParams=Data.TrialParams;

StimParams=repelem(StimParams,2,1);
StimParams(1,:)=[];
TrialParams=repelem(TrialParams,2,1);
TrialParams(1,:)=[];
TrialParams(3:2:end,3)=0;

