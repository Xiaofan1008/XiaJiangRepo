% Helper script for loading in lfp data from the current filepath
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.lfp2.dat']);
if ~isempty(tmp)
    sz = tmp.bytes/2;
    tmp = [filepath '\' tmp.name];
    l_fid = fopen(tmp,'r');
    lfp2 = fread(l_fid,[32, sz/32],'float');
    fclose(l_fid);
end