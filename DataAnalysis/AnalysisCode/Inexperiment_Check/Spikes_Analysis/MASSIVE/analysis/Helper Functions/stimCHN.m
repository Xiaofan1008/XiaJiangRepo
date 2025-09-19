% Calculate the stimulation channel from the datafile
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*exp_datafile_*.mat']);
if ~isempty(tmp)
    load(tmp.name,'StimParams');
    load(tmp.name,'TrialParams');
    if ~exist('TRIAL','var')        
        TRIAL = str2double(input('Please select a trial identifier\n','s'));        
    end
    stimROW = find(cell2mat(TrialParams(2:end,2)) == TRIAL,1);
    stimCHANNEL = StimParams{stimROW+1,1};
    stimCHANNEL = str2double(stimCHANNEL(end-2:end))+1;
    clear StimParams TrialParams tmp stimROW;
end