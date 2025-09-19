%% Function will look at the waveform stored in v_tb and
% locate XSTART and XSTOP blanking values
function [XSTART,XSTOP] = BlankingHelper(v_tb,BIN)
% Assume the artefact crosses 500 uV.
fMAX = find(v_tb(abs(BIN(1)*30)-30:abs(BIN(1)*30)+30) >= 500,1) + (abs(BIN(1)*30)-30);
fMIN = find(v_tb(abs(BIN(1)*30)-30:abs(BIN(1)*30)+30) <= -500,1) + (abs(BIN(1)*30)-30);

if isempty(fMAX) && isempty(fMIN)
    fMAX = find(v_tb(abs(BIN(1)*30)-30:abs(BIN(1)*30)+30) >= 200,1) + (abs(BIN(1)*30)-30);
    fMIN = find(v_tb(abs(BIN(1)*30)-30:abs(BIN(1)*30)+30) <= -200,1) + (abs(BIN(1)*30)-30);
    if isempty(fMAX) && isempty(fMIN)
        % Could be running an empty trial
        % Pass XSTART/XSTOP as empty
        XSTART = [];
        XSTOP = [];
        return;
    end
end

if isempty(fMAX)
    % Only fMAX is empty
    fMAX = fMIN + 1;
end
if isempty(fMIN)
    % Only fMIN is empty
    fMIN = fMAX + 1;
end

vMax = max(v_tb);
vMin = min(v_tb);
fMAX = find(v_tb == vMax,1);
fMIN = find(v_tb == vMin,1);

% Work out whether we deflect positive or negative initially
if (fMAX > fMIN)
    % Rewind back from the first crossing
    tMIN = find(v_tb((abs(BIN(1))-4)*30:(abs(BIN(1))+4)*30) == v_tb(fMIN),1) + (abs(BIN(1))-4)*30;
    tMAX = 0;
    XSTART = find(v_tb(1:tMIN) >= mean(v_tb) - std(v_tb),1,'last');    
else
    tMAX = find(v_tb((abs(BIN(1))-4)*30:(abs(BIN(1))+4)*30) == v_tb(fMAX),1) + (abs(BIN(1))-4)*30;
    tMIN = 0;
    XSTART = find(v_tb(1:tMAX) <= mean(v_tb) + std(v_tb),1,'last');
end

if isempty(tMAX) || isempty(tMIN)
    XSTART = abs(BIN(1)*30)-15;
end

if isempty(XSTART)
    XSTART = (tMAX + tMIN);
end

if (XSTART > abs(BIN(1)*30)+(4*30))
    XSTART = abs(BIN(1)*30);
end

if (XSTART < abs(BIN(1)*30)-15)
    XSTART = abs(BIN(1)*30)-15;
end

% Don't blank more than 4 milliseconds
B_ULIMIT = XSTART + 4*30;

% Don't blank less than 1 milliseconds
B_LLIMIT = XSTART + 1*30;

% Locate the optimal position within the limit window
WINDOW = [1 5];
i = 1;
while (true)
    SD = std(v_tb(WINDOW(1)+B_LLIMIT+i:WINDOW(2)+B_LLIMIT+i));
    if ((SD < 20) && abs(v_tb(i+B_LLIMIT)) < 500)
        break;
    end
    if ((i+B_LLIMIT) >= B_ULIMIT)
        break;
    end
    i = i + 1;
end
XSTOP = i+B_LLIMIT;

if (XSTOP < abs(BIN(1)*30)-30) || (XSTART > abs(BIN(1)*30)+(4*30))
    XSTART = abs(BIN(1)*30)-10;
    XSTOP = abs(BIN(1)*30)+45;
end

if XSTART > abs(BIN(1)*30) - 6
    XSTART = abs(BIN(1)*30) - 6;
end
if XSTOP < abs(BIN(1)*30) + 45
    XSTOP = abs(BIN(1)*30) + 45;
end

% Convert XSTART and XSTOP into millisecond numbers
XSTART = XSTART / 30;
XSTOP = XSTOP / 30;