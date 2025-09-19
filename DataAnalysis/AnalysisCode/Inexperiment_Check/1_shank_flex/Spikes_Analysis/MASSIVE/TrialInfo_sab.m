function [trialinfo] = TrialInfo_sab
% calculates template of trials

%OUTPUT Example of each trial ID in order - array contains all trial info 
%TrialInfo similar structure to stim params

%INTPUT none

trialinfo={};
TrialParams=loadTrialParams;
maxtid=max(cell2mat(TrialParams(:,2)));
loadStimParams;
trialinfo(1,:)=[{'ID'}, {'StimChn'}, StimParams(1,:) ];
for ID=1:maxtid
    num=find(cell2mat(TrialParams(1:end,2)) == ID);
    trialinfo(ID+(ID-1)+1,:)=[{ID},TrialParams(num(1),3), StimParams(num(1)+1,:)];
    trialinfo((ID+1)+(ID-1)+1,:)=[{ID},TrialParams(num(2),3), StimParams(num(2)+1,:)];
end
end

