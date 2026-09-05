function sp = denoiseSpikes_trig(sp,trig)
if isempty(sp)
    return
end
% This function assumes that denoiseSpikes has alreay been run.
% This function builds a template waveform from a window of baseline
% activity, and matches post-triggerline spikes to that template
nTrig = length(trig);
MAX = zeros(1,nTrig);
MIN = zeros(1,nTrig);
for n = 1:nTrig
    baseSp = sp(sp(:,1) >= trig(n)-200 & sp(:,1) <= trig(n)-100,:);
    if isempty(baseSp)
        continue;
    end
    MAX(n) = max(max(baseSp(:,2:end)));
    MIN(n) = min(min(baseSp(:,2:end)));    
end
MAX = mean(MAX);
MIN = mean(MIN);
for n = 1:nTrig
    stimSp = sp(sp(:,1) >= trig(n) & sp(:,1) <= trig(n)+10,:);
    if isempty(stimSp)
        continue;
    end
    chk = sum(stimSp(:,2:end) > MAX*1.5) | sum(stimSp(:,2:end) < MIN*1.5);
    if (chk)
        sp(sp(:,1) >= trig(n) & sp(:,1) <= trig(n)+10,:) = [];
    end
end
end