function v = removeMean(v,stimAdjust,N,T,trig)
bBIN = [-30, 120]; % Amount to blank in samples
missed_trig = trig(trig > ((N-1)*T*1000)-200 & trig < ((N-1)*T*1000)+200);
trig = trig(trig >= ((N-1)*T*1000)+200 & trig <= (N*T*1000)-200);
if ~isempty(missed_trig)
    trig = [missed_trig trig];
end
if (N == 1)
    trig = trig - ((N-1)*T*1000);
else
    trig = trig - ((N-1)*T*1000) + 1000;
end
trig = trig * 30;
for t = 1:length(trig)
    for c = 1:size(data,1)
        v(c,trig(t)+bBIN(1):trig(t)+bBIN(2)) = v(c,trig(t)+bBIN(1):trig(t)+bBIN(2)) - stimAdjust;
    end
end
end