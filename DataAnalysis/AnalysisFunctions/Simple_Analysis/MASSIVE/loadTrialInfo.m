function [trialinfo] = loadTrialInfo(varargin)
% calculates template of trials

%OUTPUT Example of each trial ID in order - array contains all trial info 
%TrialInfo similar structure to stim params

%INTPUT none

trialinfo={};
TrialParams=loadTrialParams;
maxtid=max(cell2mat(TrialParams(:,2)));
loadStimParams;
trialinfo(1,:)=[{'ID'}, {'StimChn'}, StimParams(1,:) ];
count=1;
for ID=1:maxtid
    num=find(cell2mat(TrialParams(1:end,2)) == ID);
    for i=1:min(diff(find(diff(num)~=1)))
        count=count+1;
        trialinfo(count,:)=[{ID},TrialParams(num(i),3), StimParams(num(i)+1,:)];
    end
end
if nargin==1
    trialinfo(1,:)=[];
end

end

