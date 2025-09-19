% Helper script for loading in lfp data from the current filepath
function lfp = loadLFP(chn)
lfp = [];
if ~nargin
    chnchk = 0;
else
    chnchk = 1;
end
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.lfp.dat']);
chk = strsplit(tmp.name,'_');
numChn = 32;
for i = 1:length(chk)
    if strcmp(chk{i},'analog')
        numChn = 1;
    end
end
if ~isempty(tmp)
    sz = tmp.bytes/2;
    tmp = [filepath '\' tmp.name];
    l_fid = fopen(tmp,'r');
    lfp = fread(l_fid,[numChn, sz/numChn],'float');
    fclose(l_fid);
end
if chnchk
    lfp = lfp(chn,:);
end
end
