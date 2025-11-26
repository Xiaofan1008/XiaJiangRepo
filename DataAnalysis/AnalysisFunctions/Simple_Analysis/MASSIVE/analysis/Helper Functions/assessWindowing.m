function assessWindowing(filepath)
%% Function used to check spike detection against the MUA waveform
%% Deal with input error
if nargin < 1
    return;
end
%% Pre-set data structures
sp = []; TrialParams = []; trig = []; mu = []; tr = 1;
%% User-defined variables
window = [-100 500];
file = split(filepath,'\');
if strcmp(file{end}(5:7),'007')
    CHANNEL = 23;
    depth = Depth(2)+1;
    trialID = 60;
else
    CHANNEL = 4;
    depth = Depth(3);
    trialID = 6;
end
%% Load in data
tic;
thisCHN = depth(CHANNEL);
loadTrig;
loadTrialParams;
TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == trialID));
trig = trig(TrialParams);
loadMUA;
MUA = mu(thisCHN,:);
clear mu;
if ~exist('sp','var') || ~exist('trig','var')
    disp('Error. Could not read in required data.');
    return;
else
    disp('Success. Required data loaded!')
end
toc;
%% Plot the data
figure;
ylim([-1000 1000]);
X = 1/30:1/30:(diff(window) + 1/30);
while (tr)
    tr = input('Which trial would you like to see?\n','s');
    tr = str2double(tr);
    if (tr == 0)
        break;
    end
    plot(X,MUA((trig(tr)+window(1))*30:(trig(tr)+window(2))*30));
    wave = MUA((trig(tr)+window(1))*30:(trig(tr)+window(2))*30);
    [sp1,spk1] = spikeextract(wave,-52.9544/2,30000);
    [sp2,spk2] = spikeextract(wave,-52.9544,30000);
    sp = [sp1; sp2];
    spk = [spk1; spk2];
    [sp,spk,~] = windowSpikesOut(sp,spk,1);
    [spk, index] = unique(spk,'rows');
    sp = sp(index,:);
    [sp, index] = unique(sp(:,1:end),'rows','stable');
    spk = spk(index);    
    for i = 1:length(spk)
        line([spk(i) spk(i)],[150 1000],'Color','r');
    end
end
end
