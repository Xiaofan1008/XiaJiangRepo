function tmp = newTrialParams(n_Trials)
%% Initialises a new cell array of Stimulus Parameters
tmp = cell(n_Trials+1,3);
tmp{1,1} = 'Trial Number';
tmp{1,2} = 'Trial ID';
tmp{1,3} = 'CHANNEL';

end