%% Spike Windowing function
function [WAVES,Sp,REJECTED] = windowSpikesOut(WAVES,Sp,threshsign)
%% Function takes in a series of spikes and applies user-defined windows to remove waveforms that aren't spikes
%% User-defined windows
%% CASE 1
% Near 0 before the spike
W1X = [2 2 12 12 2];
W1Y = [-25 25 25 -25 -25];
% Where the spike threshold should be
W2X = [9 9 17 17 9];
W2Y = [-200 -80 -80 -200 -200];
% Rising above zero after the spike threshold
W3X = [12 12 22 22 12];
W3Y = [25 75 75 25 25];
% Returning to zero before the end of the waveform
W4X = [30 30 44 44 30];
W4Y = [-10 10 10 -10 -10];
% All together now
W1 = [W1X, W1Y; W2X, W2Y; W3X, W3Y; W4X, W4Y];
%% CASE 2
% Near 0 before the spike
W1X = [7 7 9 9 7];
W1Y = [25 150 150 25 25];
% Where the spike threshold should be
W2X = [9 9 17 17 9];
W2Y = [-200 -80 -80 -200 -200];
% Rising above zero after the spike threshold
W3X = [22 22 27 27 22];
W3Y = [20 -40 -40 20 20];
% Returning to zero before the end of the waveform
W4X = [35 35 44 44 35];
W4Y = [-30 30 30 -30 -30];
% All together now
W2 = [W1X, W1Y; W2X, W2Y; W3X, W3Y; W4X, W4Y];
%% CASE 3
% Near 0 before the spike
W1X = [5 5 10 10 5];
W1Y = [50 250 250 50 50];
% Where the spike threshold should be
W2X = [9 9 17 17 9];
W2Y = [-300 -200 -200 -300 -300];
% Rising above zero after the spike threshold
W3X = [22 22 27 27 22];
W3Y = [-50 -150 -150 -50 -50];
% Returning to zero before the end of the waveform
W4X = [35 35 44 44 35];
W4Y = [-30 30 30 -30 -30];
% All together now
W3 = [W1X, W1Y; W2X, W2Y; W3X, W3Y; W4X, W4Y];
%% CASE 4
% Near 0 before the spike
W1X = [2 2 12 12 2];
W1Y = [-35 25 25 -35 -35];
% Where the spike threshold should be
W2X = [9 9 17 17 9];
W2Y = [-150 -35 -35 -150 -150];
% Rising above zero after the spike threshold
W3X = [22 22 27 27 22];
W3Y = [-30 25 25 -30 -30];
% Returning to zero before the end of the waveform
W4X = [30 30 44 44 30];
W4Y = [-30 30 30 -30 -30];
% All together now
W4 = [W1X, W1Y; W2X, W2Y; W3X, W3Y; W4X, W4Y];
%% Invert based on threshsign
W1(:,6:end) = W1(:,6:end) .* threshsign;
W2(:,6:end) = W2(:,6:end) .* threshsign;
W3(:,6:end) = W3(:,6:end) .* threshsign;
W4(:,6:end) = W4(:,6:end) .* threshsign;
%% Calculate which spike waveforms fall outside the user-defined windows
ARTMAX = 1e3;
Xs = 1:1:49;
for i = 1:length(Sp)
    chk1 = 0;
    chk2 = 0;
    chk3 = 0;
    chk4 = 0;
    [~,m] = min(WAVES(i,:),[],2);
    if (m ~= 13)
        Sp(i) = NaN;
        continue;
    end
    if max(abs(WAVES(i,:))) > ARTMAX
        Sp(i) = NaN;
        continue;
    end
    for j = 1:4
        [xi,~] = polyxpoly(Xs,WAVES(i,:),W1(j,1:5),W1(j,6:10));
        if isempty(xi)
            chk1 = 1;
            break;
        end
    end
    for j = 1:4
        [xi,~] = polyxpoly(Xs,WAVES(i,:),W2(j,1:5),W2(j,6:10));
        if isempty(xi)
            chk2 = 1;
            break;
        end
    end
    for j = 1:4
        [xi,~] = polyxpoly(Xs,WAVES(i,:),W3(j,1:5),W3(j,6:10));
        if isempty(xi)
            chk3 = 1;
            break;
        end
    end
    for j = 1:4
        [xi,~] = polyxpoly(Xs,WAVES(i,:),W4(j,1:5),W4(j,6:10));
        if isempty(xi)
            chk4 = 1;
            break;
        end
    end
    if chk1 == 1 && chk2 == 1 && chk3 == 1 && chk4 == 1
        Sp(i) = NaN;
    else
        continue;
    end    
end
REJECTED = WAVES(isnan(Sp),:);
WAVES(isnan(Sp),:) = [];
Sp(isnan(Sp)) = [];

%% How to see the rectangles being used
%mapshow(W1X,W1Y,'DisplayType','polygon','LineStyle','none');
%mapshow(W2X,W2Y,'DisplayType','polygon','LineStyle','none');
%mapshow(W3X,W3Y,'DisplayType','polygon','LineStyle','none');
%mapshow(W4X,W4Y,'DisplayType','polygon','LineStyle','none');

end