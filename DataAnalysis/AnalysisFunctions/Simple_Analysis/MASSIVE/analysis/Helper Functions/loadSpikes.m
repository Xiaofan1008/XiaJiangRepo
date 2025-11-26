% Helper script for loading in spikes data from the current filepath
function sp = loadSpikes(chn)
if ~nargin
    chnchk = 0;
else
    chnchk = 1;
end
filepath = pwd;
tmp = dir([filepath filesep '*.sp.mat']);
if ~isempty(tmp)
    tmp = [filepath filesep tmp.name];
    sp = load(tmp,'sp');
     sp = sp.sp;
end
if chnchk
    sp = sp{chn};
end
end
