function whatisN
dbstop if error
%% Load population data
[all, sig, nsig] = loadPopStats;
%% Apply Rules
% Rule 1: Baseline firing rate must exceed 5 spikes/s
bR = 1;
% Rule 2: r^2 must exceed 0.9 for thresholding fit
r2lim = 0.9;
% Rule 3: Channel must saturate for threshold to be accurate
sat = 0; % This tells applyThresholdRules to eliminate channels with saturation values of 0.
% Rule 4: Channel must reach 50 spikes/s with stimulation to be viable
maxfR = 50;
[all,~,~] = applyThresholdRules(all,sig,nsig,bR,r2lim,sat,maxfR);
% Right now, we only want one threshold value per stimulating channel.
thresh = []; %stimDepth = []; distSR = [];
recordings = unique(all.recording);
% Each recordings contributes unique stimulating channels
for nR = recordings
    nRAll.thresh = all.thresh(all.recording == nR);
    nRAll.durID = all.durID(all.recording == nR);
    nRAll.stimChn = all.stimChn(all.recording == nR);
    %nRAll.stimDepth = all.depthStim(all.recording == nR);
    %nRAll.distSR = all.distSR(all.recording == nR);
    durID = unique(nRAll.durID);
    % Each duration ID contributes unique stimulating conditions
    for nD = durID
        nDAll.thresh = nRAll.thresh(nRAll.durID == nD);
        nDAll.stimChn = nRAll.stimChn(nRAll.durID == nD);
        %nDAll.stimDepth = nRAll.stimDepth(nRAll.durID == nD);
        %nDAll.distSR = nRAll.distSR(nRAll.durID == nD);
        stimChn = unique(nDAll.stimChn);
        % For each unique stimChannel, derive a threshold
        for nS = stimChn
            nSAllthresh = nDAll.thresh(nDAll.stimChn == nS);
            %nSAllstimDepth = nDAll.stimDepth(nDAll.stimChn == nS);
            %nSAlldistSR = nDAll.distSR(nDAll.stimChn == nS);
            [~,ind] = min(nSAllthresh);
            thresh = [thresh nSAllthresh(ind)]; %#ok<*AGROW>
            %stimDepth = [stimDepth nSAllstimDepth(ind)];
            %distSR = [distSR nSAlldistSR(ind)];
        end
    end
end
fprintf('n is %3.0f\n',length(thresh));
end