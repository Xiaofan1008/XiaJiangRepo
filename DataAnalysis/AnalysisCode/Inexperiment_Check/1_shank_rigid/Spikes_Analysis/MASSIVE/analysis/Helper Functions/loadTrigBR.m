% Helper script for loading in lfp data from the current filepath
function trig = loadTrigBR(chk)
filepath = pwd;
tmp = dir([filepath '\*trigBR.dat']);
if ~isempty(tmp)
    sz = tmp.bytes/2;
    tmp = [filepath '\' tmp.name];
    t_fid = fopen(tmp,'r');
    trig = fread(t_fid,[1, sz],'double');
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
    trig = trig ./ 30;
    trig = cast(trig,'int64');
end
end