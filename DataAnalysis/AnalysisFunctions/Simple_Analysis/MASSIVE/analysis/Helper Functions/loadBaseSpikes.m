% Helper script for loading in lfp data from the current filepath
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.basesp.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    sp = load(tmp,'base_sp');
    sp = sp.base_sp;
end