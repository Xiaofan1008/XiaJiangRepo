function wave = blank(wave)
dbstop on error
% This is another blanking function
%% Variables
thresh = 1.5e3; % Threshold for artefact detection
BIN = [-100 100];
nData = size(wave,2);
nChn = size(wave,1);
%% First pass threshold crossing detection
chk = zeros(nChn,nData);
for c = 1:nChn
    chk(c,:) = (abs(wave(c,:)) > thresh);
    % chk contains a logical vector of every point exceeding threshold
    % We'll use the first point as the start of the blanking window
    lines = flip(find(chk(c,:) == 1));
    dt = diff(lines);
    kill = dt == -1;
    lines(kill) = [];    
    % Discard lines that are too close together.
    % We'll run this script several times to prevent missing spots.
    if isempty(lines)
        continue;
    end
    hold = lines(end);
    lines = lines(diff(lines) < -300);
    lines = [lines hold];
    for l = 1:length(lines)
        if lines(l) < abs(BIN(1)*30)+1 || lines(l) > length(wave(c,:)) - abs(BIN(2)*30) - 1
            continue;
        end        
        while (true) % We will continue to blank the same segment until there is no change            
            a = wave(c,lines(l)+(BIN(1)*30):lines(l)+(BIN(2)*30));
            [XSTART,XSTOP] = BlankHere(wave(c,lines(l)+(BIN(1)*30):lines(l)+(BIN(2)*30)),BIN);            
            % Perform the blanking
            wave(c,lines(l)+(BIN(1)*30):lines(l)+(BIN(2)*30)) = BlankThis(wave(c,lines(l)+(BIN(1)*30):lines(l)+(BIN(2)*30)),XSTART,XSTOP);
            b = wave(c,lines(l)+(BIN(1)*30):lines(l)+(BIN(2)*30));
            if abs(sum(a - b)) <= 1e-3 % The difference between the two is really small
                break;
            end
        end
    end
end
end

function [XSTART,XSTOP] = BlankHere(wave,BIN)
% This function calculates the x positions to blank between
B_LLIMIT = 30;
B_ULIMIT = 120;
% We now know that the centre point of this wave is the first peak of
% the blanking period

%% Deal with the easy case
if std(wave(1:abs(BIN(1)*30)/2)) < 20
    XSTART = find(abs(wave(1:abs(BIN(1)*30))) < mean(abs(wave(1:abs(BIN(1)*30)/2))),1,'last');
    if XSTART < abs(BIN(1)*30)-30
        XSTART = abs(BIN(1)*30);
        for i = 1:B_LLIMIT
            if std(wave(XSTART-19:XSTART+1)) > std(wave(1:abs(BIN(1)*30)))
                XSTART = XSTART - 1;
            else
                break;
            end
        end
    end
    XSTART = XSTART - 20;
    XSTOP = XSTART + B_LLIMIT;
    for i = 1:(B_ULIMIT - B_LLIMIT)
        if std(wave(XSTOP-10:XSTOP+10)) > std(wave(abs(BIN(1)*30)+30:end))
            XSTOP = XSTOP + 1;
        else
            break
        end
    end
    return;
end
 
%% Peak-to-peak window method
% Assume the first thirty points are fine
ptp = max(wave(1:1+30)) - min(wave(1:1+30));
XSTART = []; XSTOP = [];
for n = 1:length(wave)-30
    PTP = max(wave(n:n+30)) - min(wave(n:n+30));
    if PTP > ptp * 6.0
        XSTART = n+29;
        break
    end
end

if isempty(XSTART)
    return
end

XSTART = XSTART - 15;
XSTOP = XSTART + B_LLIMIT;

if XSTOP > (diff(BIN)*30+1) - B_ULIMIT
    XSTART = [];
    XSTOP = [];
    return
end

for i = 1:(B_ULIMIT - B_LLIMIT)
    if std(wave(XSTOP-1:XSTOP+19)) > std(wave(abs(BIN(1))*30:end))
        XSTOP = XSTOP + 1;
    else
        break
    end    
end
return

end

function wave = BlankThis(wave,XSTART,XSTOP)
% This function does the actual blanking
if isempty(XSTART) && isempty(XSTOP)    
    return;
end
X = [XSTART,XSTOP];
[CF, S, mu] = polyfit(X,[wave(XSTART),wave(XSTOP)],1);
X = XSTART:1:XSTOP;
Y = polyval(CF,X,S,mu);
wave(XSTART:XSTOP) = Y;
end