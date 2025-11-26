% Helper script for loading in spikes data from the current filepath
function sp = loadSpikesBR(chn)
if ~nargin
    chnchk = 0;
else
    chnchk = 1;
end
filepath = pwd;
tmp = dir([filepath '\*.spBR.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    sp = load(tmp,'sp');
    sp = sp.sp;
end
if chnchk
    sp = sp{chn};
end
end
