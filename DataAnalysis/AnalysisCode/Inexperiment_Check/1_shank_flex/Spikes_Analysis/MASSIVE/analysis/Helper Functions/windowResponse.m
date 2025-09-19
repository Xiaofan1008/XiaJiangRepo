%% Script takes the directory after running openIntan
% and calculates spiking rates after a defined WINDOW
%% User-defined variables
WINDOWS = zeros(10,length(time_stamps));
WINDOW = 10;
%CHANNEL = 18;
ARTIFACT_CUTOFF = 300;
threshsign = 1;
%% Main Logic
depth = Depth(2)+1;
chan = depth(CHANNEL);
[WAVES,Sp] = spikeextract(MU(chan,:),thresh(threshsign,chan),FS);
% Window the spikes
[WAVES,Sp] = windowSpikes(WAVES,Sp,threshsign);
% Sort spikes to remove artefact
if sum(Sp) ~= 0
    for l = 1:length(Sp)
        % Chuck away any spikes that exceed the threshold
        if min(WAVES(l,:)) < -ARTIFACT_CUTOFF || max(WAVES(l,:)) > ARTIFACT_CUTOFF
            Sp(l) = NaN;
        end
    end
    WAVES(isnan(Sp),:) = [];
    Sp(isnan(Sp)) = [];
else
    Sp(Sp == 0) = [];
end
% Sort spikes into windowed bins
for l = 1:length(Sp)
    COUNT = 1;
    while (COUNT < 10)
        for t = 1:length(time_stamps)
            if Sp(l) >= time_stamps(t) + ((COUNT-1)*WINDOW) && Sp(l) <= time_stamps(t) + (COUNT)*WINDOW
                WINDOWS(COUNT,t) = WINDOWS(COUNT,t) + 1;
                break;
            end            
            if Sp(l) < time_stamps(t)
                break;
            end
        end
        COUNT = COUNT + 1;
    end
end
disp(['There are ' num2str(length(Sp)) ' spikes across ' num2str(length(MU(1,:))/FS) ' seconds.']);
disp('In each window, there are: ');
for i = 1:size(WINDOWS,1)
    disp([num2str(sum(WINDOWS(i,:))) ' spikes in ' num2str(WINDOW/1000) ' seconds between ' ...
        num2str((i-1)*WINDOW/1000) ' and ' num2str(i*WINDOW/1000) ' seconds after visual stimulation.']);
end
disp(['This represents ' num2str(sum(sum(WINDOWS))/length(Sp)) ' of all spikes.']);