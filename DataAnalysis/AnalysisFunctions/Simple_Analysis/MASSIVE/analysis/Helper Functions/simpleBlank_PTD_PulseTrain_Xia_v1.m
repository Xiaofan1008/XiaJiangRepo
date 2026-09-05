function data = simpleBlank_PTD_PulseTrain_Xia_v1(data,N,T,trig,mode,FS)

TrialParams_raw = loadTrialParams;

% Original logic:
% TrialParams column 2 contains the condition ID.
% We still use this to estimate how many stimulation rows belong to each trial.
TrialParams_condition_all = cell2mat(TrialParams_raw(:,2))';
num_elect = min(diff(find(diff(TrialParams_condition_all)~=0))); % number of electrode rows used per trial

% One condition ID per trial, in the randomized trial order.
% This keeps the old variable name "TrialParams" so the rest of the code
% behaves similarly to the old version.
TrialParams = TrialParams_condition_all(1:num_elect:end);

loadStimParams;

% =========================================================================
% XIA MODIFICATION: build actual artifact event times for each trial
% =========================================================================
% Old code only read one pulse count and one PTD value per trial:
%
%   StimParams_pulseinfo = cell2mat(StimParams(2:num_elect:end,8:9));
%   StimParams_PulseDelay_info = ...
%
% That is not enough for the new event-level pulse-train design, because
% different electrode rows can have different pulse counts, e.g.
%
%   [1 0], [1 1], [2 1], [2 2], [3 2], [3 3]
%
% Also, inactive electrode rows are disabled by amplitude = -1, but may
% still have a dummy numberOfStimPulses = 1. Therefore, we must use:
%
%   firstPhaseAmplitude > 0
%
% to identify real active stimulation rows.
%
% For each active electrode row:
%
%   event times = PTD + (0:numberOfStimPulses-1) * pulseTrainPeriod
%
% For each trial, EventTimes_perTrial{k} stores all actual stimulation
% artifact times in us.
% =========================================================================

nTrialsFromStimParams = floor((size(StimParams,1)-1) / num_elect);
EventTimes_perTrial = cell(1, nTrialsFromStimParams);
IsZeroTrial_perTrial = false(1, nTrialsFromStimParams);

for tr = 1:nTrialsFromStimParams

    % StimParams has a header row, so trial tr starts at row:
    %   2 + (tr-1)*num_elect
    rows_this_trial = (2 + (tr-1)*num_elect) : (1 + tr*num_elect);

    amp_this_trial = cell2mat(StimParams(rows_this_trial,16)); % firstPhaseAmplitude
    active_rows = amp_this_trial > 0;                          % active stimulation only

    if ~any(active_rows)
        % Zero-current / baseline control trial.
        % No artifact should be blanked.
        EventTimes_perTrial{tr} = [];
        IsZeroTrial_perTrial(tr) = true;
        continue;
    end

    PTD_this_trial = cell2mat(StimParams(rows_this_trial,6));       % post-trigger delay, us
    nPulse_this_trial = cell2mat(StimParams(rows_this_trial,8));    % number of pulses
    period_this_trial = cell2mat(StimParams(rows_this_trial,9));    % pulse train period, us

    eventTimes_us = [];

    for r = find(active_rows(:))'

        thisPTD = PTD_this_trial(r);
        thisNPulse = nPulse_this_trial(r);
        thisPeriod = period_this_trial(r);

        if thisNPulse > 0
            eventTimes_us = [eventTimes_us, thisPTD + (0:thisNPulse-1) * thisPeriod];
        end
    end

    EventTimes_perTrial{tr} = unique(sort(eventTimes_us));
end

% -------------------------------------------------------------------------
% Old scalar pulse-train / PTD extraction is kept here for reference only.
% It is no longer used by the new event-level artifact blanking logic.
% -------------------------------------------------------------------------
% StimParams_pulseinfo=cell2mat(StimParams(2:num_elect:end,8:9));
%
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

    % how many samples before the first stimulation event should be included
    % in the blanking window
    BIN_b = -FS/1000; % Amount to blank in samples; at 30kHz, about 1ms

    % find the trials might be missed if they located close to the boundary
    % of each chunk (trigger at the end of last chunk, but data in the
    % current chunk)
    missed_idx = trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000;
    missed_trials = TrialParams(missed_idx); % condition ID for missed trials
    missed_trig = trig(missed_idx);          % sample index
    temp = find(missed_idx);                 % trial indices in randomized order

    % Remove invalid missed trigger placeholders
    invalid_missed = (missed_trig == -500);
    missed_trials(invalid_missed) = [];
    missed_trig(invalid_missed) = [];
    temp(invalid_missed) = [];

    % triggers inside the chunk, not close to the edge
    index = trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000;

    trials = TrialParams(index);       % condition ID
    trial_indices = find(index);       % trial indices in randomized order
    trig = trig(index);

    % XIA MODIFICATION:
    % Instead of getting one pulse-train row per trial, get the actual
    % reconstructed event times for each trial.
    trialEventTimes = EventTimes_perTrial(trial_indices);
    isZeroTrial = IsZeroTrial_perTrial(trial_indices);

    if ~isempty(missed_trig)
        trials = [missed_trials trials];
        trig = [missed_trig trig];

        missedEventTimes = EventTimes_perTrial(temp);
        missedIsZeroTrial = IsZeroTrial_perTrial(temp);

        trialEventTimes = [missedEventTimes trialEventTimes];
        isZeroTrial = [missedIsZeroTrial isZeroTrial];
    end

    if (N == 1)
        trig = trig - ((N-1)*T*FS);
    else
        trig = trig - ((N-1)*T*FS) + FS;
    end

    shifttime = min(diff(trig)) - 0.101*FS;

    % ---------------------------------------------------------------------
    % Old bug fix block kept for reference only.
    % It modified pulsetrain_info, which is no longer used in this version.
    % ---------------------------------------------------------------------
    % wd=dir;
    % if (strcmp(wd(3).folder(end-12:end-7),'220907') || strcmp(wd(3).folder(end-12:end-7),'220908')) && strcmp(wd(3).folder(end-25:end-14),'PEN2_SIGMOID')
    %     pulsetrain_info(:,1)=2;
    %     pulsetrain_info(:,2)=10000/3;%fix bug
    % end

    % ---- Main Processing Loop ---- %
    for t = 1:length(trig)

        % Don't blank zero-current/baseline trials.
        % In the new design, these are identified by all active amplitudes
        % being <= 0, so EventTimes_perTrial is empty.
        if isZeroTrial(t) || isempty(trialEventTimes{t})
            continue;
        end

        if trials(t) == 5
            pause(0.01);
        end

        thisTrig = trig(t);

        % XIA MODIFICATION:
        % Get all real stimulation event times for this trial.
        %
        % Example:
        %   PTD = [0 5] ms, event-level train 3 = [2 1]
        %   eventTimes_us = [0 5000 10000]
        %
        % We follow the old blanking style:
        %   - find adaptive recovery after the LAST event
        %   - blank/interpolate from the FIRST event to that recovery point
        eventTimes_us = trialEventTimes{t};
        firstEvent_us = min(eventTimes_us);
        lastEvent_us  = max(eventTimes_us);

        firstEvent_sample_offset = round(firstEvent_us * FS / 1e6);
        lastEvent_sample_offset  = round(lastEvent_us  * FS / 1e6);

        for c = 1:nChn

            % -----------------------------------------------------------------
            % Old scalar branch logic is kept here for reference only.
            % It is no longer used because it cannot handle event-level trains.
            % -----------------------------------------------------------------
            % if pulsetrain_info(t,1)<2 && StimParams_PulseDelay_info(t) == 0
            %     ...
            % elseif pulsetrain_info(t,1) >= 2
            %     ...
            % elseif StimParams_PulseDelay_info(t) > 0 && pulsetrain_info(t,1)<2
            %     ...
            % end

            % ----- (1) Find adaptive artifact recovery after the LAST event -----
            %
            % This keeps the old adaptive method:
            %   start after the last stimulation event
            %   slide a 46-sample window until signal range is <= 150 uV
            %
            % Then ra(1) is converted back to be relative to the original trigger.
            thisTrigtemp = round(thisTrig + lastEvent_sample_offset);

            ra = [1 46];      % window of 46 samples; around 1.5 ms at 30 kHz
            range_max = 150;  % artifact considered recovered if range <= 150 uV

            while range(data(c,thisTrigtemp+ra(1):thisTrigtemp+ra(2))) > range_max
                ra = ra + 1;

                if ra(1) > 180
                    ra(1) = 180;
                    break;
                end
            end

            ra(1) = ra(1) + 10; % add buffer, same as old code

            if ra(1) > 180
                ra(1) = 180;
            end

            % Convert ra from "relative to last event" to "relative to trigger".
            ra = ra + lastEvent_sample_offset;

            % The blanking starts before the first actual event.
            % For most of your designs firstEvent_us = 0, so this is identical
            % to the old:
            %   thisTrig + BIN_b
            blankStart = round(thisTrig + firstEvent_sample_offset + BIN_b);
            blankEnd   = round(thisTrig + ra(1));

            % ----- (2) Perform artifact interpolation + shift correction -----
            %
            % This keeps the old style:
            %   interpolate from first event to recovery after last event.
            %
            % No spike recovery is performed here.
            if thisTrig < length(data)-FS

                data(c,blankStart:blankEnd) = interpolate(data(c,blankStart:blankEnd),1);

                shift = diff([data(c,blankStart), data(c,blankEnd)]);

                data(c,blankEnd:thisTrig+shifttime) = ...
                    data(c,blankEnd:thisTrig+shifttime) - shift;

                data(c,blankStart:blankEnd+1) = ...
                    interpolate(data(c,blankStart:blankEnd+1),1);
            end

            % if StimParams{
            % data(c,thisTrig+BIN_b) = detrend(data(c,thisTrig+BIN_b),1);

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

            % distant_channel = 13 + c; % Get as far away as possible
            % if distant_channel > 26
            %     distant_channel = distant_channel - 26;
            % end
            % if distant_channel == 2 || distant_channel == 9
            %     distant_channel = distant_channel + 1; % Deal with dead channels
            % end

            trend1 = interpolate(data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)),1);
            wave = (data(c,trig(t)*30 + bBIN(2):trig(t)*30 + bBIN(2) + diff(bBIN)));
            trend2 = interpolate(wave,1);
            wave = wave - trend2 + trend1;
            data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)) = wave;

            % data(c,trig(t)*30 - 5:trig(t)*30 + 5) = interpolate(data(c,trig(t)*30 - 5:trig(t)*30 + 5),1);
        end
    end
end

end