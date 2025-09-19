function v = final_blank(v)
%% Final pass blanking, where we give up and hit it with a huge blank
dbstop on error
%% Variables
ARTMAX = 1.5e3; % Artifact threshold
BW = [-1 4];  % Blanking window (msec)
nData = size(v,2);
nChn = size(v,1);
%% Logic Loop
chk = zeros(nChn,nData);
for c = 1:nChn
    chk(c,:) = (abs(v(c,:)) > ARTMAX);
    lines = flip(find(chk(c,:) == 1));
    dt = diff(lines);
    kill = dt == -1;
    lines(kill) = [];
    if isempty(lines)
        continue;
    end
    hold = lines(end);
    lines = lines(diff(lines) < -150);
    lines = [lines hold];
    for l = 1:length(lines)
        if lines(l) < abs(BW(1)*30)+1 || lines(l) > length(v(c,:)) - abs(BW(2)*30) - 1
            continue;
        end
        XSTART = lines(l) + (BW(1)*30);
        XSTOP = lines(l) + (BW(2)*30);
        X = [XSTART,XSTOP];
        [CF, S, mu] = polyfit(X,[v(c,XSTART),v(c,XSTOP)],1);
        X = XSTART:1:XSTOP;
        Y = polyval(CF,X,S,mu);
        v(c,XSTART:XSTOP) = Y;
    end
end
end