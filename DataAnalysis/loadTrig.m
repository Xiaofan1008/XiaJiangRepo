% Helper script for loading trigger times
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
 trig = fread(t_fid,[1, sz],'double'); % byte rate = 18560; 'double' = 64bits = 8 bytes; 
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