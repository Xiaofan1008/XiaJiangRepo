% Helper script for loading in lfp data from the current filepath
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.basespw.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    sp = load(tmp,'base_spw');
    sp = sp.base_spw;
end