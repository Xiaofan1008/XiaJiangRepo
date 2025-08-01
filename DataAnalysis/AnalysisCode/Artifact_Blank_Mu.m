function Artifact_Blank_Mu(data,N,T,trig,FS)
% This function simply blank artifact for MUA analysis. 
% Adaptive atrifact window used


% -- load trial ID --
TrialParams = loadTrialParams; 
TrialParams = cell2mat(TrialParams(:,2));
num_elect=min(diff(find(diff(TrialParams)~=0))); % number of electrodes per trial
TrialParams=TrialParams(1:num_elect:end,:);

% -- load 'num_pulse per train' and 'Pulse train period' 
loadStimParams;
StimParams_pulseinfo=cell2mat(StimParams(2:num_elect:end,8:9)); 
data = v(iChn,:);
nChn = size(data,1);

% -- trials cross two chunks --
% search the triggers around [-0.1s, 0.1s] of the
% boundary of prev and current chunk
missed_trials = TrialParams(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
missed_trig = trig(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
% remove missing triggers
missed_trials((missed_trig==-500))=[];
missed_trig((missed_trig==-500))=[];

% -- triggers login --
% find trials within current data chunk
index=trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000;
trials = TrialParams(index);
pulsetrain_info=StimParams_pulseinfo(index,:);
trig = trig(index);

if ~isempty(missed_trig) % trials cross two chunks 
    missed_pulseinfo=cell2mat(StimParams(temp+1,8:9));    missed_pulseinfo((missed_trig==-500))=[];
    trials = [missed_trials trials];
    trig = [missed_trig trig];
    pulsetrain_info = [missed_pulseinfo; pulsetrain_info];
end

% -- realign time of triggers --
if (N>1) % compensate 1s data overlap in the following chunks
    trig = trig - ((N-1)*T*FS) + FS;
end

% -- time window for voltage jump subtractation -- 
shifttime=min(diff(trig))-0.101*FS;

 % -- Artifact remove --
for t = 1:length(trig) % Loop through trials
    % Artifacts will be removed post hoc by
    % interpolating the signal from 1ms before the
    % stimulation pulse until the post-stimulation
    % signal returned to +- 150µV of the
    % pre-stimulation and remained there for 1ms. 
    thisTrig = trig(t);
    % -- Find the adaptive range --
    range_win = [1 30]; % start window range ~ 1.5 ms
    range_max = 150; % 150µV threshold
    % loop through data 
    while (range(data((thisTrig+Bin_ms): (thisTrig + range_win(1)))) > range_max) || (range(data((thisTrig+range_win(1)): (thisTrig + range_win(2))))> range_max)
        range_win = range_win + 1; % move window forward 1 sample
        if range_win(1) > win_max % maximum at 6ms
            range_win(1) = win_max;
            break
        end
    end

    %% Interpolate
    % 1. Artifact rejection - interpolate between -1ms and
    % 1st sample fall within +-150µV
    % 2. Remove residual artifact by subtracting
    % difference beween pre- sample and post- samples. 
    % 3. Reinterpolate for smoothing 
    data((thisTrig+Bin_ms) : (thisTrig+range_win(1))) = interpolate(data((thisTrig+Bin_ms) : (thisTrig+range_win(1))),1); % interpolate for artifact
    shift = diff([data(thisTrig+Bin_ms), data(thisTrig+range_win(1))]); % voltage jump
    data((thisTrig+range_win(1)) : (thisTrig+3000)) = data((thisTrig+range_win(1)) : (thisTrig+3000)) - shift; % correct data connect smoothy for artifact  
    data((thisTrig+Bin_ms) : (thisTrig+range_win(1)+1)) = interpolate(data((thisTrig+Bin_ms) : (thisTrig+range_win(1)+1)),1);
end

end