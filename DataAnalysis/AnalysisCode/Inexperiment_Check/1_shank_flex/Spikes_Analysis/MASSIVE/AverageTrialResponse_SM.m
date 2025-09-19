function [avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct, baslinespikestruct)
% calculates template of trials and spiking responses

%OUTPUT matrix of average spiking avgnospT(electrodes, trials), array of
%standard error stderrspktrial(electrodes, trials),Example of each trial ID
%in order - array contains all trial info TrialInfo similar structure to
%stim params

%INTPUT Cell array structure of trials grouped into trial IDs i.e.
%cell{1}=array(electrode,trial) for trialID 1 where each coluumn contains 
%number of spikes in this trial for the given electrode

%Elecrode properties
filepath = pwd;
[amplifier_channels,~]=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);
if nChn>32
    E_Mapnumber=1;
else
    E_Mapnumber=0;
end
E_MAP = Depth(E_Mapnumber);
avgnospT=[];
stderrspktrial=[];
trialinfo={};
TrialParams=loadTrialParams;
maxtid=max(cell2mat(TrialParams(:,2)));
loadNREP;
% try%%%%%%use for timechangeinspiking
%     loadoriginalEND;
%     if ((originalEND-1)/(n_REP_true*2))<maxtid
%         maxtid=((originalEND-1)/(n_REP_true*2));
%     end
% catch
%     %%code was made before loadoriginalEND
% end
loadStimParams;
for ID=1:maxtid
    check=['T' num2str(ID)];
    IDstructnobsp.(check)=IDstruct.(check)-baslinespikestruct.(check);
    avgnospT=meanstruct(avgnospT, ID, IDstructnobsp);
    stderrspktrial=stderrorstruct(stderrspktrial, ID, IDstructnobsp);
    num=find(cell2mat(TrialParams(1:end,2)) == ID);
    trialinfo(ID+(ID-1),:)=[{ID},TrialParams(num(1),3), StimParams(num(1)+1,:)];
    trialinfo((ID+1)+(ID-1),:)=[{ID},TrialParams(num(2),3), StimParams(num(2)+1,:)];
end
tempavg=avgnospT;
avgnospT=avgnospT(E_MAP,:);
stderrspktrial=stderrspktrial(E_MAP,:);
end

