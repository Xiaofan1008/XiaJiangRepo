% Helper script for loading in lfp data from the current filepath
function trig = loadTrig(chk)
if ~nargin
    chk = 0;
end
filepath = pwd;
tmp = dir([filepath filesep '*.trig.dat']);
if ~isempty(tmp)
    sz = tmp.bytes/2;
    tmp = [filepath filesep tmp.name];
    t_fid = fopen(tmp,'r');
 trig = fread(t_fid,[1, sz],'double');% something wrong with trig size  - told to read a 1 by 884 and only got 221
    fclose(t_fid);
else    
    return;
end

if isempty(trig)
    return
end
if ~(chk)
    return;
elseif (chk==1)
    trig = trig ./ (FS/1000);
    trig = cast(trig,'int64');
end
end