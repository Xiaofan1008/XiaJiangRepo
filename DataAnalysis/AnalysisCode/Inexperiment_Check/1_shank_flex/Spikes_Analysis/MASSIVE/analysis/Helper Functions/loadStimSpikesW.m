% Helper script for loading in stimspikes data from the current filepath
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.stimspw.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    sp = load(tmp,'stim_spw');
    sp = sp.stim_spw;
end