function data = simpleBlank(data,N,T,trig,mode,FS)
TrialParams = loadTrialParams; TrialParams = cell2mat(TrialParams(:,2))';
num_elect=min(diff(find(diff(TrialParams)~=0))); % number of electrodes used per trial
TrialParams=TrialParams(1:num_elect:end,:);
loadStimParams;StimParams_pulseinfo=cell2mat(StimParams(2:num_elect:end,8:9));

nChn = size(data,1); 
if mode == 1 % for HPF
    % how many samples before the trigger time should be include in the blanking window
    BIN_b = -FS/1000; % Amount to blank in samples; at 30kHz, about 1ms
    
    % find the trials might be missed if they located close to the boundary
    % of each chunk (trigger at the end of last chunk, but data in the
    % current chunk)
    missed_trials = TrialParams(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000); % gives the trial label
    missed_trig = trig(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000); % gives the sample index
    temp=find(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
    missed_trials((missed_trig==-500))=[];
    missed_trig((missed_trig==-500))=[];
    index=trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000; % triggers inside the chunk (not close to the edge) 
    trials = TrialParams(index); % trial ID 
    pulsetrain_info=StimParams_pulseinfo(index,:); % pulse information for included trials in the chunk
    trig = trig(index);
    if ~isempty(missed_trig)
        missed_pulseinfo=cell2mat(StimParams(temp+1,8:9));    missed_pulseinfo((missed_trig==-500))=[];
        trials = [missed_trials trials];
        trig = [missed_trig trig];
        pulsetrain_info = [missed_pulseinfo; pulsetrain_info];
    end
    if (N == 1)
        trig = trig - ((N-1)*T*FS);
    else
        trig = trig - ((N-1)*T*FS) + FS;
    end
    shifttime=min(diff(trig))-0.101*FS;
    wd=dir;
    if (strcmp(wd(3).folder(end-12:end-7),'220907') || strcmp(wd(3).folder(end-12:end-7),'220908')) && strcmp(wd(3).folder(end-25:end-14),'PEN2_SIGMOID')
        pulsetrain_info(:,1)=2;
        pulsetrain_info(:,2)=10000/3;%fix bug
    end
    for t = 1:length(trig) %used to go from 1
        % Don't blank zero trials
%         if trials(t) == 1
%             continue;
%         end


        if trials(t) == 5
            pause(0.01);
        end
        thisTrig = trig(t);
        for c = 1:nChn
            %for pulsetrain=1:pulsetrain_info(t,1)
            %thisTrig=round(trig(t)+((pulsetrain-1)*(pulsetrain_info(t,2)*(FS/10^6))));
            % Find the appropriate range
            if pulsetrain_info(t,1)<2 % single pulse case
            ra = [1 46]; % window of 46 samples; around 1.5ms
            range_max = 150; % Artifact considered "gone" if within +- 150mV
            % if signal range (max - min) > 150 µV, still artifact
            while range(data(c,thisTrig+ra(1):thisTrig+ra(2))) > range_max
                ra = ra + 1; % shift the window right if still noisy
                if ra(1) > 180 % stop after 180 samples
                    ra(1) = 180;
                    break;
                end
            end
            ra(1) = ra(1) + 10; % add buffer
            if ra(1) > 180
                ra(1) = 180;
            end
            else % Pulse train case
                thisTrigtemp=round(trig(t)+((pulsetrain_info(t,1)-1)*(pulsetrain_info(t,2)*(FS/10^6))));
                ra = [1 46];
                range_max = 150;
                while range(data(c,thisTrigtemp+ra(1):thisTrigtemp+ra(2))) > range_max
                    ra = ra + 1;
                    if ra(1) > 180
                        ra(1) = 180;
                        break;
                    end
                end
                ra(1) = ra(1) + 10;
                if ra(1) > 180
                    ra(1) = 180;
                end
                ra=ra+(pulsetrain_info(t,1)-1)*(pulsetrain_info(t,2)*(FS/10^6));
            end
            % Shift the entire data array to minimize artefact. Convenient
            
            %original code
            % place to reintroduce artefact is +100 msec before next trig
            if thisTrig<length(data)-FS
                data(c,thisTrig+BIN_b:thisTrig+ra(1)) = interpolate(data(c,thisTrig+BIN_b:thisTrig+ra(1)),1); % interpolate to remove the artifact
                shift = diff([data(c,thisTrig+BIN_b),data(c,thisTrig+ra(1))]); % the jump discontinuity between pre- and post-artifact values
                data(c,thisTrig+ra(1):thisTrig+shifttime) = data(c,thisTrig+ra(1):thisTrig+shifttime) - shift; % remove the step jump cause by artifact
                data(c,thisTrig+BIN_b:thisTrig+ra(1)+1) = interpolate(data(c,thisTrig+BIN_b:thisTrig+ra(1)+1),1); % keep the smooth transition (after remove the shift correction)
            end
            %if StimParams{
            %data(c,thisTrig+BIN_b) = detrend(data(c,thisTrig+BIN_b),1);
            

        end
    end
elseif mode == 2 % for LPF
    bBIN = [5, 250*30]; % Amount to blank in samples
    % Check for missed trigger lines
    missed_trig = trig(trig > ((N-1)*T*1000)-550 & trig < ((N-1)*T*1000)+550);
    trig = trig(trig >= ((N-1)*T*1000)+550 & trig <= (N*T*1000)-550);
    if ~isempty(missed_trig)
        trig = [missed_trig trig];
    end
    if (N == 1)
        trig = trig - ((N-1)*T*1000);
    else
        trig = trig - ((N-1)*T*1000) + 1000;
    end
    trig = cast(trig,'int64');
    for t = 1:length(trig)
        for c = 1:size(data,1)
            %distant_channel = 13 + c; % Get as far away as possible
            %if distant_channel > 26
            %    distant_channel = distant_channel - 26;
            %end
            %if distant_channel == 2 || distant_channel == 9
            %    distant_channel = distant_channel + 1; % Deal with dead channels
            %end
            trend1 = interpolate(data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)),1);
            wave = (data(c,trig(t)*30 + bBIN(2):trig(t)*30 + bBIN(2) + diff(bBIN)));
            trend2 = interpolate(wave,1);
            wave = wave - trend2 + trend1;
            data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)) = wave;
            %data(c,trig(t)*30 - 5:trig(t)*30 + 5) = interpolate(data(c,trig(t)*30 - 5:trig(t)*30 + 5),1);
        end
    end
end
end