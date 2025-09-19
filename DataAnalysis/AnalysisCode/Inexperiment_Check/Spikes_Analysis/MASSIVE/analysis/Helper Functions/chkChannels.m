% Helper script for checking for bad channels on an array
if ~exist('filepath','var')
    filepath = pwd;
end
AMP_LIMIT = 1e3;
loadRAW;
loadTrig;
trig = cast(trig(trig < length(v)./30),'int32');
for t = 1:length(trig)
    v(:,trig(t)*30 - 30:trig(t)*30 + 120) = 0;
end
chk = ones(32,1);
for c = 1:size(v,1)    
    MAX = max(v(c,:));
    MIN = min(v(c,:));
    AMPLITUDE = MAX - MIN;
    if (AMPLITUDE > AMP_LIMIT)
        chk(c) = 0;
    end
end