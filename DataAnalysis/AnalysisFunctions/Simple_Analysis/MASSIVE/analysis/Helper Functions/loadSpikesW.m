% Helper script for loading in spikes data from the current filepath
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.spw.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    sp = load(tmp,'sp_win');
    sp = sp.sp_win;
end