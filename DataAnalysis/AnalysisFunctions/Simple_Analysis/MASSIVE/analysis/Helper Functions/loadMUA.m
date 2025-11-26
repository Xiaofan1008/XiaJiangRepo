function mu = loadMUA(trig)
% Helper script for loading in mua data from the current filepath
if ~nargin
    trig = [];
end
if ~exist('filepath','var')
    filepath = pwd;
end
tmp = dir([filepath '\*.mu.dat']);
chk = strsplit(tmp.name,'_');
nChn = 32; FS = 30000;
for i = 1:length(chk)
    if strcmp(chk{i},'analog')
        nChn = 1;
    end
end
sz = tmp.bytes/2;
sz = min(sz,30*2880000);
if ~isempty(tmp) && ~isempty(trig)
    tmp = [filepath '\' tmp.name];
    m_fid = fopen(tmp,'r');
    OFFSET = cast(nChn*2*trig,'int64');
    fseek(m_fid,OFFSET,'bof');
    mu = fread(m_fid,[nChn, sz/nChn],'int16') .* 0.195;    
    fclose(m_fid);
else
    tmp = [filepath '\' tmp.name];
    m_fid = fopen(tmp,'r');
    mu = fread(m_fid,[nChn, sz/nChn],'short');
    mu = mu ./ 10;
    fclose(m_fid);
end
end