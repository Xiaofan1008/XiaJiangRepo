function data = simpleBlank_PTD_MultiElec_Xia_v2(data,N,T,trig,mode,FS)

% ============================================================
% simpleBlank_PTD_MultiElec_Xia_v2
%
% Purpose:
%   Artifact blanking for mixed-prefix multi-electrode stimulation files.
%
% Important for new files:
%   TrialParams_raw has one row per stimulation slot.
%   Example:
%       990 trials × 5 electrode slots = 4950 rows
%
%   But trig is trial-level:
%       990 triggers
%
%   Therefore, TrialParams, activeCount_perTrial, pulseDelay_info,
%   and pulsetrain_info must all be reduced to trial-level length.
%
% Main modification:
%   Use active electrode count per trial to determine the final active
%   artifact timing for prefix 1/2/3/4/5 trials.
%
% ============================================================

%% ================= LOAD TRIAL PARAMS =================
% Keep row-level TrialParams first.
TrialParams_raw = loadTrialParams;

TrialParams_allRows = cell2mat(TrialParams_raw(:,2));
TrialParams_allRows = TrialParams_allRows(:)';   % force row vector

%% ================= DETERMINE ROWS/SLOTS PER TRIAL =================
% num_elect means "number of saved stimulation rows per trial".
% In a 5-electrode prefix file, this should be 5.
change_idx = find(diff(TrialParams_allRows) ~= 0);

if isempty(change_idx)
    num_elect = 1;
else
    num_elect = min(diff(change_idx));
end

if isempty(num_elect) || isnan(num_elect) || num_elect < 1
    num_elect = 1;
end

%% ================= REDUCE TO TRIAL-LEVEL TRIALPARAMS =================
% TrialParams must have the same length as trig.
TrialParams = TrialParams_allRows(1:num_elect:end);
TrialParams = TrialParams(:)';   % force row vector

nTrials_fromTrialParams = length(TrialParams);

%% ================= ACTIVE ELECTRODE COUNT PER TRIAL =================
% New mixed-prefix TrialParams column 4 stores the true active electrode count.
%
% Example:
%   zero control -> 0
%   prefix 1     -> 1
%   prefix 2     -> 2
%   ...
%   prefix 5     -> 5
if size(TrialParams_raw,2) >= 4
    try
        activeCount_allRows = cell2mat(TrialParams_raw(:,4));
        activeCount_allRows = activeCount_allRows(:)';  % force row vector

        activeCount_perTrial = activeCount_allRows(1:num_elect:end);
        activeCount_perTrial = activeCount_perTrial(:)'; % force row vector
    catch
        activeCount_perTrial = num_elect * ones(1, nTrials_fromTrialParams);
    end
else
    activeCount_perTrial = num_elect * ones(1, nTrials_fromTrialParams);
end

% Safety check.
if length(activeCount_perTrial) ~= length(TrialParams)
    error(['activeCount_perTrial length (%d) does not match TrialParams length (%d). ' ...
           'Check TrialParams_raw and num_elect.'], ...
           length(activeCount_perTrial), length(TrialParams));
end

% fprintf('\nArtifact blanking metadata check:\n');
% fprintf('  TrialParams_raw rows        = %d\n', size(TrialParams_raw,1));
% fprintf('  num_elect / rows per trial  = %d\n', num_elect);
% fprintf('  trial-level TrialParams     = %d\n', length(TrialParams));
% fprintf('  activeCount_perTrial        = %d\n', length(activeCount_perTrial));
% fprintf('  unique active counts        = ');
% disp(unique(activeCount_perTrial));

%% ================= LOAD STIM PARAMS =================
loadStimParams;

% This should now be trial-level pulse train info.
% Column 8 = pulse train number
% Column 9 = pulse train period, in us
%
% We use the first row of each trial. This is okay if pulse train info is
% the same across rows, which is usually true in your current files.
StimParams_pulseinfo = cell2mat(StimParams(2:num_elect:end,8:9));

% Make sure pulseinfo is trial-level.
if size(StimParams_pulseinfo,1) ~= length(TrialParams)
    error('StimParams_pulseinfo length (%d) does not match TrialParams length (%d).', ...
          size(StimParams_pulseinfo,1), length(TrialParams));
end

%% ================= FINAL ACTIVE ARTIFACT TIME PER TRIAL =================
% StimParams_PulseDelay_info is trial-level.
%
% It stores the final active artifact time for each trial, in us.
%
% For current prefix files with PulseNum = 1:
%   this is the last active PTD.
%
% For future pulse-train files:
%   final artifact time = PTD + (PulseNum - 1) * PulsePeriod

StimParams_PulseDelay_info = zeros(1, length(TrialParams));

for trial_i = 1:length(TrialParams)

    activeCount_thisTrial = activeCount_perTrial(trial_i);

    if isnan(activeCount_thisTrial) || activeCount_thisTrial < 1
        % zero-control trial or invalid metadata
        StimParams_PulseDelay_info(trial_i) = 0;
    else
        activeCount_thisTrial = min(round(activeCount_thisTrial), num_elect);

        % StimParams has one header row, so trial_i starts at:
        %   row 2 + (trial_i-1)*num_elect
        firstStimRow = 2 + (trial_i-1)*num_elect;
        activeRows = firstStimRow : (firstStimRow + activeCount_thisTrial - 1);

        % Column 6: post-trigger delay in us
        ptd_us = cell2mat(StimParams(activeRows,6));

        % Column 8: pulse train number
        pulseNum = cell2mat(StimParams(activeRows,8));

        % Column 9: pulse train period in us
        pulsePeriod_us = cell2mat(StimParams(activeRows,9));

        pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
        pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

        rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

        StimParams_PulseDelay_info(trial_i) = max(rowFinalArtifact_us);
    end
end

fprintf('  unique final artifact times us = ');
disp(unique(StimParams_PulseDelay_info));

nChn = size(data,1);

%% ============================================================
% HIGH-PASS / SPIKE EXTRACTION MODE
% ============================================================
if mode == 1

    % Blank from 1 ms before the trigger.
    BIN_b = -FS/1000;   % at 30 kHz, -30 samples = -1 ms

    %% ================= FIND MISSED TRIGGERS NEAR CHUNK BOUNDARY =================
    % missed_index and index are trial-level masks.
    % Therefore TrialParams, StimParams_pulseinfo, pulseDelay, and activeCount
    % must all be trial-level vectors.

    missed_index = trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000;
    missed_trials = TrialParams(missed_index);
    missed_trig = trig(missed_index);
    temp = find(missed_index);

    valid_missed = missed_trig ~= -500;
    missed_trials = missed_trials(valid_missed);
    missed_trig = missed_trig(valid_missed);
    temp = temp(valid_missed);

    %% ================= FIND TRIGGERS INSIDE CURRENT CHUNK =================
    index = trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000;

    trials = TrialParams(index);
    pulsetrain_info = StimParams_pulseinfo(index,:);

    % Keep PTD/final-artifact and active count aligned with filtered triggers.
    pulseDelay_info = StimParams_PulseDelay_info(index);
    activeCount_info = activeCount_perTrial(index);

    trig = trig(index);

    %% ================= ADD MISSED TRIGGERS BACK =================
    if ~isempty(missed_trig)

        missed_pulseinfo = StimParams_pulseinfo(temp,:);
        missed_pulseDelay_info = StimParams_PulseDelay_info(temp);
        missed_activeCount_info = activeCount_perTrial(temp);

        trials = [missed_trials trials];
        trig = [missed_trig trig];
        pulsetrain_info = [missed_pulseinfo; pulsetrain_info];
        pulseDelay_info = [missed_pulseDelay_info pulseDelay_info];
        activeCount_info = [missed_activeCount_info activeCount_info];
    end

    %% ================= SHIFT TRIGGERS INTO CHUNK COORDINATES =================
    if N == 1
        trig = trig - ((N-1)*T*FS);
    else
        trig = trig - ((N-1)*T*FS) + FS;
    end

    %% ================= SAFE SHIFT TIME =================
    % shifttime controls how far the post-artifact step correction extends.
    % Original code used:
    %   shifttime = min(diff(trig)) - 0.101*FS;
    % but this can fail if only one trigger exists in a chunk.

    if length(trig) >= 2
        shifttime = min(diff(trig)) - 0.101*FS;
    else
        % Conservative default if only one trigger exists in this chunk.
        shifttime = 0.05 * FS;   % 50 ms
    end

    % Make sure shifttime is at least 10 ms and integer.
    shifttime = max(round(shifttime), round(0.01*FS));

    %% ================= OLD SPECIAL BUG FIX =================
    wd = dir;
    if (strcmp(wd(3).folder(end-12:end-7),'220907') || strcmp(wd(3).folder(end-12:end-7),'220908')) && ...
            strcmp(wd(3).folder(end-25:end-14),'PEN2_SIGMOID')
        pulsetrain_info(:,1) = 2;
        pulsetrain_info(:,2) = 10000/3;
    end

    %% ================= MAIN BLANKING LOOP =================
    for t = 1:length(trig)

        % In the new mixed-prefix file, zero-control trials have active count = 0.
        % These trials should not be blanked.
        if activeCount_info(t) == 0
            continue;
        end

        thisTrig = round(trig(t));

        % Skip invalid trigger positions.
        if thisTrig <= abs(BIN_b) + 2 || thisTrig >= size(data,2)
            continue;
        end

        for c = 1:nChn

            %% ================= DETERMINE BLANKING END RANGE =================
            % ra is in samples relative to thisTrig.
            % It is shifted to the final artifact time when needed.

            if pulsetrain_info(t,1) < 2 && pulseDelay_info(t) == 0
                %% Single pulse without PTD
                artifactTrig = thisTrig;

            else
                %% Delayed multi-electrode or pulse-train case
                % pulseDelay_info is already the final active artifact time
                % in us, including pulse train if applicable.
                artifactTrig = round(thisTrig + pulseDelay_info(t) * (FS/1e6));
            end

            % Check that artifactTrig is inside data.
            if artifactTrig <= 0 || artifactTrig >= size(data,2)
                continue;
            end

            %% ----- Find when artifact returns to baseline range -----
            ra = [1 46];          % around 1.5 ms window at 30 kHz
            range_max = 150;      % artifact considered gone if range <= 150 uV

            % Prevent indexing beyond data length.
            while artifactTrig + ra(2) <= size(data,2) && ...
                    range(data(c, artifactTrig+ra(1):artifactTrig+ra(2))) > range_max

                ra = ra + 1;

                if ra(1) > 180
                    ra(1) = 180;
                    break;
                end
            end

            % Add small buffer.
            ra(1) = ra(1) + 10;

            if ra(1) > 180
                ra(1) = 180;
            end

            % Convert ra from artifact-relative to trigger-relative.
            artifactOffset_samples = artifactTrig - thisTrig;
            ra = ra + artifactOffset_samples;

            %% ================= APPLY INTERPOLATION + STEP CORRECTION =================
            blank_start = round(thisTrig + BIN_b);
            blank_end   = round(thisTrig + ra(1));

            % Boundary protection.
            if blank_start < 1
                blank_start = 1;
            end

            if blank_end > size(data,2)
                blank_end = size(data,2);
            end

            if blank_end <= blank_start + 2
                continue;
            end

            shift_end = round(thisTrig + shifttime);

            if shift_end > size(data,2)
                shift_end = size(data,2);
            end

            if shift_end <= blank_end
                shift_end = min(size(data,2), blank_end + round(0.01*FS));
            end

            if thisTrig < length(data)-FS

                % First interpolation over artifact window.
                data(c, blank_start:blank_end) = interpolate(data(c, blank_start:blank_end), 1);

                % Estimate jump discontinuity between pre-artifact and post-artifact signal.
                shift = diff([data(c, blank_start), data(c, blank_end)]);

                % Remove step shift after artifact.
                data(c, blank_end:shift_end) = data(c, blank_end:shift_end) - shift;

                % Re-interpolate the blanked region after shift correction.
                reinterp_end = min(size(data,2), blank_end + 1);
                data(c, blank_start:reinterp_end) = interpolate(data(c, blank_start:reinterp_end), 1);
            end
        end
    end

%% ============================================================
% LOW-PASS MODE
% ============================================================
elseif mode == 2

    bBIN = [5, 250*30]; % Amount to blank in samples

    % This LPF branch is mostly unchanged from the original code.
    missed_trig = trig(trig > ((N-1)*T*1000)-550 & trig < ((N-1)*T*1000)+550);
    trig = trig(trig >= ((N-1)*T*1000)+550 & trig <= (N*T*1000)-550);

    if ~isempty(missed_trig)
        trig = [missed_trig trig];
    end

    if N == 1
        trig = trig - ((N-1)*T*1000);
    else
        trig = trig - ((N-1)*T*1000) + 1000;
    end

    trig = cast(trig,'int64');

    for t = 1:length(trig)
        for c = 1:size(data,1)

            trend1 = interpolate(data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)),1);

            wave = data(c,trig(t)*30 + bBIN(2):trig(t)*30 + bBIN(2) + diff(bBIN));

            trend2 = interpolate(wave,1);

            wave = wave - trend2 + trend1;

            data(c,trig(t)*30 + bBIN(1):trig(t)*30 + bBIN(2)) = wave;
        end
    end
end

end