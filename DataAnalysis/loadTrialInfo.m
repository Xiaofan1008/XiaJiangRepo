function [trialinfo] = loadTrialInfo(varargin)
% calculates template of trials

%OUTPUT Example of each trial ID in order - array contains all trial info 
%TrialInfo similar structure to stim params

%INTPUT none
global DATA_FOLDER;

TrialParams = loadTrialParams();% 1: Trial Number' 2: Trial ID; 3: Channel
StimParams = loadStimParams();
trialinfo = {};
maxtid = max(cell2mat(TrialParams(:,2)));% Find the total numebr of trials (maximum number of trial ID)
trialinfo(1,:) = [{'ID'}, {'StimChn'}, StimParams(1,:)];
for ID = 1:maxtid
    num = find(cell2mat(TrialParams(:,2)) == ID);
    trialinfo(ID+(ID-1)+1,:) = [{ID}, TrialParams(num(1),3), StimParams(num(1)+1,:)];
    trialinfo((ID+1)+(ID-1)+1,:) = [{ID}, TrialParams(num(2),3), StimParams(num(2)+1,:)];
end

if nargin==1
    trialinfo(1,:)=[];
end
end

% trialinfo={};
% TrialParams=loadTrialParams; % 1: Trial Number' 2: Trial ID; 3: Channel
% maxtid=max(cell2mat(TrialParams(:,2))); % Find the total numebr of trials (maximum number of trial ID)
% loadStimParams;
% StimParams = loadStimParams('/Users/xiaofan/Desktop/PhD Study/Data/Sabrina/S3E2_9elect_001_210511_104603');
% trialinfo(1,:)=[{'ID'}, {'StimChn'}, StimParams(1,:) ];
% for ID=1:maxtid
%     num=find(cell2mat(TrialParams(1:end,2)) == ID);
%     trialinfo(ID+(ID-1)+1,:)=[{ID},TrialParams(num(1),3), StimParams(num(1)+1,:)];
%     trialinfo((ID+1)+(ID-1)+1,:)=[{ID},TrialParams(num(2),3), StimParams(num(2)+1,:)];
% end
% if nargin==1
%     trialinfo(1,:)=[];
% end
% 
% end