%% Spike Windowing function
function [WAVES,Sp,REJECTED] = windowSpikes(WAVES,Sp,threshsign)
%% Function takes in a series of spikes and applies user-defined windows to remove waveforms that aren't spikes
%% User-defined windows
% Large depolarisation substantially after the expected spike
W1X = [25 25 45 45 25];
W1Y = [-600 -100 -100 -600 -600];
% Large hyperpolarisation substantially after the spike
W2X = [35 35 45 45 35];
W2Y = [100 600 600 100 100];
% All together now
W1 = [W1X, W1Y; W2X, W2Y];
%% Invert based on threshsign
W1(:,6:end) = W1(:,6:end) .* threshsign;
%% Calculate which spike waveforms fall outside the user-defined windows
Xs = 1:1:49;
for i = 1:length(Sp)
    chk1 = 0;
    [~,m] = min(WAVES(i,:),[],2);
    if (m ~= 13)
        Sp(i) = NaN;
        continue;
    end
    for j = 1:size(W1,1)
        [xi,~] = polyxpoly(Xs,WAVES(i,:),W1(j,1:5),W1(j,6:10));
        if ~isempty(xi)
            chk1 = 1;
        end        
    end    
    if chk1 == 1
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