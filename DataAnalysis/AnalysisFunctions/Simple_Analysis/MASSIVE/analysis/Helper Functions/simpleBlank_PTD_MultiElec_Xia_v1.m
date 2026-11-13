function data = simpleBlank_PTD_MultiElec_Xia_v1(data,N,T,trig,mode,FS)

% ================= XIA MODIFICATION =================
% Keep the raw TrialParams first, because the new mixed-prefix files contain
% extra metadata columns, including ActiveElectrodeCount in column 4.
TrialParams_raw = loadTrialParams; 
TrialParams = cell2mat(TrialParams_raw(:,2));
% =====================================================

num_elect=min(diff(find(diff(TrialParams)~=0))); % number of electrodes used per trial

% ================= XIA MODIFICATION =================
% In old files, num_elect = number of active electrodes.
% In new mixed-prefix files, num_elect = number of saved rows per trial
%                              = maximum electrode slots.
% The true number of active electrodes is saved in TrialParams column 4.
if size(TrialParams_raw,2) >= 4
    try
        activeCount_allRows = cell2mat(TrialParams_raw(:,4))';
        activeCount_perTrial = activeCount_allRows(1:num_elect:end);
    catch
        % Backward compatibility for old files without metadata.
        activeCount_perTrial = num_elect * ones(1,length(TrialParams(1:num_elect:end)));
    end
else
    % Backward compatibility for old files without metadata.
    activeCount_perTrial = num_elect * ones(1,length(TrialParams(1:num_elect:end)));
end
% =====================================================

TrialParams=TrialParams(1:num_elect:end,:);
loadStimParams;
StimParams_pulseinfo=cell2mat(StimParams(2:num_elect:end,8:9));

% if num_elect >= 2
%     StimParams_PulseDelay_info=cell2mat(StimParams(3:num_elect:end,6));
% else 
%     StimParams_PulseDelay_info=cell2mat(StimParams(2:num_elect:end,6));
% end

% ================= XIA MODIFICATION =================
% Old logic used num_elect to decide which electrode row contains the last
% artifact. This is not correct for mixed-prefix files because num_elect is
% always the maximum row number, e.g. 5, while the active electrode count can
% be 1, 2, 3, 4, or 5.
%
% New logic:
% For each trial, use activeCount_perTrial to find the PTD of the last
% active electrode. This gives the artifact that we need to blank until.
%
% Example:
%   prefix 3 in a 5-electrode file:
%   PTD = [0 5000 10000 0 0]
%   activeCount = 3
%   last active PTD = row 3 = 10000 us
StimParams_PulseDelay_info = zeros(1,length(TrialParams));

for trial_i = 1:length(TrialParams)
    activeCount_thisTrial = activeCount_perTrial(trial_i);

    if isnan(activeCount_thisTrial) || activeCount_thisTrial < 1
        % zero-control trial or invalid metadata
        StimParams_PulseDelay_info(trial_i) = 0;
    else
        activeCount_thisTrial = min(round(activeCount_thisTrial), num_elect);

        % StimParams has a header row, so data row starts at row 2.
        % For trial_i, electrode slot activeCount_thisTrial is:
        % 1 + (trial_i-1)*num_elect + activeCount_thisTrial
        stimRow = 1 + (trial_i-1)*num_elect + activeCount_thisTrial;
        StimParams_PulseDelay_info(trial_i) = StimParams{stimRow,6};
    end
end
% =====================================================

% if num_elect == 2
%     StimParams_PulseDelay_info=cell2mat(StimParams(3:num_elect:end,6)); % 2 electrode stimulation 
% elseif num_elect == 3
%     StimParams_PulseDelay_info=cell2mat(StimParams(4:num_elect:end,6)); % 3 electrode stimulation
% elseif num_elect == 4
%     StimParams_PulseDelay_info=cell2mat(StimParams(5:num_elect:end,6)); % 4 electrode stimulation
% elseif num_elect == 5
%     StimParams_PulseDelay_info=cell2mat(StimParams(6:num_elect:end,6)); % 5 electrode stimulation
% else 
%     StimParams_PulseDelay_info=cell2mat(StimParams(2:num_elect:end,6)); % Single electrode stim
% end

nChn = size(data,1); 

%% High-Pass -> Spike extraction
if mode == 1 % for HPF
    % how many samples before the trigger time should be include in the blanking window
    BIN_b = -FS/1000; % Amount to blank in samples; at 30kHz, about 1ms
    
    % find the trials might be missed if they located close to the boundary
    % of each chunk (trigger at the end of last chunk, but data in the
    % current chunk)
    % --idnetify trials inside 

    % ================= XIA MODIFICATION =================
    % Store the missed-trial mask so TrialParams, PTD, active count, and
    % pulsetrain info can stay aligned.
    missed_index = trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000;
    missed_trials = TrialParams(missed_index); % gives the trial label
    missed_trig = trig(missed_index); % gives the sample index
    temp=find(missed_index);

    valid_missed = missed_trig ~= -500;
    missed_trials = missed_trials(valid_missed);
    missed_trig = missed_trig(valid_missed);
    temp = temp(valid_missed);
    % =====================================================

    % missed_trials = TrialParams(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000); % gives the trial label
    % missed_trig = trig(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000); % gives the sample index
    % temp=find(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
    % missed_trials((missed_trig==-500))=[];
    % missed_trig((missed_trig==-500))=[];

    index=trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000; % triggers inside the chunk (not close to the edge) 
    
    trials = TrialParams(index); % trial ID 
    pulsetrain_info=StimParams_pulseinfo(index,:); % pulse information for included trials in the chunk

    % ================= XIA MODIFICATION =================
    % Keep PTD and active-electrode-count aligned with the filtered triggers.
    pulseDelay_info = StimParams_PulseDelay_info(index);
    activeCount_info = activeCount_perTrial(index);
    % =====================================================

    trig = trig(index);
    
    if ~isempty(missed_trig)

        % ================= XIA MODIFICATION =================
        % The old code used StimParams(temp+1,8:9), which is not safe for
        % multi-row trials. Use the already trial-level StimParams_pulseinfo.
        missed_pulseinfo = StimParams_pulseinfo(temp,:);
        missed_pulseDelay_info = StimParams_PulseDelay_info(temp);
        missed_activeCount_info = activeCount_perTrial(temp);

        trials = [missed_trials trials];
        trig = [missed_trig trig];
        pulsetrain_info = [missed_pulseinfo; pulsetrain_info];
        pulseDelay_info = [missed_pulseDelay_info pulseDelay_info];
        activeCount_info = [missed_activeCount_info activeCount_info];
        % =====================================================

        % missed_pulseinfo=cell2mat(StimParams(temp+1,8:9));    missed_pulseinfo((missed_trig==-500))=[];
        % trials = [missed_trials trials];
        % trig = [missed_trig trig];
        % pulsetrain_info = [missed_pulseinfo; pulsetrain_info];
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


    % ---- Main Processing Loop ---- %
    for t = 1:length(trig) %used to go from 1
        % Don't blank zero trials
%         if trials(t) == 1
%             continue;
%         end

        % ================= XIA MODIFICATION =================
        % In the new mixed-prefix file, zero-control trials have
        % active electrode count = 0. These trials should not be blanked.
        if activeCount_info(t) == 0
            continue;
        end
        % =====================================================

        if trials(t) == 5
            pause(0.01);
        end
        thisTrig = trig(t);
        for c = 1:nChn
            %for pulsetrain=1:pulsetrain_info(t,1)
            %thisTrig=round(trig(t)+((pulsetrain-1)*(pulsetrain_info(t,2)*(FS/10^6))));
            % ----- （1）Find the appropriate range ----- 

            % ================= XIA MODIFICATION =================
            % Use pulseDelay_info(t), which is the PTD of the last active
            % electrode for this trial, after chunk filtering.
            if pulsetrain_info(t,1)<2 && pulseDelay_info(t) == 0 % single pulse case without PTD
            % =====================================================

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
            elseif pulsetrain_info(t,1) >= 2 % Pulse train case
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
            
            % ================= XIA MODIFICATION =================
            elseif pulseDelay_info(t) > 0 && pulsetrain_info(t,1)<2
                thisTrigtemp = round(trig(t) + pulseDelay_info(t) * (FS/1e6));
            % =====================================================

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

                % ================= XIA MODIFICATION =================
                ra=ra+(pulseDelay_info(t) * (FS/1e6));
                % =====================================================
            end
            % Shift the entire data array to minimize artefact. Convenient
            
            % ----- (2) Perform artifact interpolation + Shift correction 
            % original code
            % place to reintroduce artefact is +100 msec before next trig
            if thisTrig<length(data)-FS
                data(c,thisTrig+BIN_b:thisTrig+ra(1)) = interpolate(data(c,thisTrig+BIN_b:thisTrig+ra(1)),1); % interpolate to remove the artifact
                shift = diff([data(c,thisTrig+BIN_b),data(c,thisTrig+ra(1))]); % the jump discontinuity between pre- and post-artifact values
                % shift = diff([data(c,thisTrig+BIN_b),data(c,round(thisTrig+ra(1)))]); % the jump discontinuity between pre- and post-artifact values
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