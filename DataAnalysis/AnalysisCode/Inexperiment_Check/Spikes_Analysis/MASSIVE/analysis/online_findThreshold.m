function online_findThreshold
dbstop if error
% Navigate to the right place
exp = dir('*datafile_*.mat');
while isempty(exp)
    % No experimental datafile here
    [~, path, ~] = ...
        uigetfile('*.mat', 'Select a .mat experimental data file', 'MultiSelect', 'off');
    cd(path);
    exp = dir('*datafile_*.mat');
end
% Check whether the datafile has been processed and has the right number of
% trigger lines
trig = [];
while isempty(trig)
    try
        loadNTrials; trig = loadTrig(0); nTrig = length(trig);
        if nTrig < n_Trials
            fprintf('Incorrect number of trigger lines detected. Please assess trigger line extraction.\n');
            keyboard;
        end
        if n_Trials < nTrig
            trig = trig(1:n_Trials);
        end
    catch
        cleanTrig;
    end
    lfp = dir('*.lfp.dat');
    if isempty(lfp)
        fprintf('Has this directory been processed?\n');
        fprintf('We will run processIntan(1,1) in this directory.\n');
        processIntan(1,1);
    end
end
%% Everything seems fine. Loop through thresholds
thresh = zeros(1,32);
for i = 1:32
    tmp = findThreshold(i,0);
    if isempty(tmp(~isnan(tmp)))
        thresh(i) = NaN;
        continue;
    else
        thresh(i) = tmp(3);
    end    
end
fprintf(['The mean ' char(177) ' ste threshold is %2.1f ' char(177) ' %2.1f\n'],nanmean(thresh),nanstd(thresh)./sqrt(32));
figure; hold on;
histogram(thresh,0:1:14);
end