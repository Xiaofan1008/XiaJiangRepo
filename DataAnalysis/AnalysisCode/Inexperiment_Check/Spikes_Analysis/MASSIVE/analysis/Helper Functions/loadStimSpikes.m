% Helper script for loading in stimspikes data from the current filepath
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.stimsp.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    sp = load(tmp,'stim_sp');
    sp = sp.stim_sp;
end