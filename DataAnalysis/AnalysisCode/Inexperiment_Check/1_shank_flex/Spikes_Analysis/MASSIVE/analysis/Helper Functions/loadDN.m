function v = loadDN(trig)
% Helper script for loading in denoised data from the current filepath
if ~nargin
    trig = [];
end
if ~exist('filepath','var')
    filepath = pwd;
end
nChn = 32; FS = 30000;
tmp = dir([filepath '\amplifier_dn.dat']);
if ~isempty(tmp) && ~isempty(trig)
    tmp = [filepath '\' tmp.name];
    v_fid = fopen(tmp,'r');
    OFFSET = cast(nChn*2*trig,'int64');
    fseek(v_fid,OFFSET,'bof');
    v = fread(v_fid,[nChn, (FS)*30],'int16') .* 0.195;    
    fclose(v_fid);
else
    tmp = [filepath '\' tmp.name];
    v_fid = fopen(tmp,'r');    
    fseek(v_fid,0,'bof');
    v = fread(v_fid,[nChn, (FS)*30],'int16') .* 0.195;    
    fclose(v_fid);
end
