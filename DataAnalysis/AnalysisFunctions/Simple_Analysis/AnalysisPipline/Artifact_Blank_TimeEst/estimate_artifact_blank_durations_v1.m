% ESTIMATE_ARTIFACT_BLANK_DURATIONS_SCRIPT
%   Recomputes the REAL adaptive artifact-blanking duration for every
%   (trial, channel) pair, using the exact same recovery-detection logic
%   as simpleBlank_PTD_PulseTrain_Xia_v1.m -- but without re-running the
%   full streaming blank/filter/spike pipeline.
%
%   Trigger + trial-parameter files (already computed) are read from the
%   "Results" folder. Only small windows of the RAW amplifier.dat around
%   each trigger are read from the "DATA" folder (which may live on a
%   network share) -- the whole multi-GB file is never loaded.
%
%   Just edit the "USER SETTINGS" block below and press Run.
%
% Output:
%   summaryTable in the workspace, one row per (trial, channel):
%       AnimalID, SessionName, Trial, Channel, ConditionID,
%       nActiveElectrodes, MaxAmplitude_uA, nEvents,
%       RecoveryDuration_ms, TotalBlankSpan_ms, EventSpan_ms
%   Also saved to outputFolder (filename auto-generated per animal/session).
%
% ASSUMPTIONS (check these against your data if results look off):
%   - amplifier.dat is int16, channels interleaved sample-by-sample
%   - raw samples are scaled by x0.195 to convert to uV (same as
%     denoiseIntan_sab.m)
%   - StimParams columns: 6=PTD(us), 8=nPulses, 9=period(us), 16=firstPhaseAmplitude
%   - trig(t) is an ABSOLUTE sample index into the full continuous
%     recording (true for loadTrig(0) output as used elsewhere in this
%     pipeline)
%   - resultsFolder path looks like .../Results/<AnimalID>/<SessionName>
%     (used for auto-detecting AnimalID/SessionName -- override manually
%     below if your folder structure differs)

clear; clc;

%% ===================== USER SETTINGS =====================
resultsFolder = '/Volumes/Shared/Xia/Results/DX014/Xia_Seq_Sim2';
rawDataFolder = '/Volumes/Shared/Xia/DATA/DX014/Xia_Seq_Sim2_251202_130158';
nChn          = 64;              % 32 or 64
FS            = 30000;           % sampling rate, Hz

% Where to save the results (kept local/small -- this is just a table,
% not the raw data). Just specify a FOLDER -- the filename is generated
% automatically. Leave empty to save in MATLAB's current folder.
outputFolder = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Artifact_Time_Est';

% Animal ID, e.g. 'DX013'. Leave empty to auto-detect from resultsFolder
% (assumes .../Results/<AnimalID>/<SessionName>, i.e. the animal ID is
% the parent folder of resultsFolder). Set manually if that doesn't hold.
animalID = '';
%% ===========================================================

if isempty(outputFolder)
    outputFolder = pwd;
end
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
[~, sessionName] = fileparts(resultsFolder);

if isempty(animalID)
    [parentDir, ~] = fileparts(resultsFolder);
    [~, animalID] = fileparts(parentDir);
    fprintf('Auto-detected AnimalID = %s (from resultsFolder path)\n', animalID);
end
outputMatFile = fullfile(outputFolder, [animalID '_' sessionName '_blankDurations.mat']);

bytesPerSample = 2;   % int16
scaleFactor    = 0.195;

%% ---- Clean up macOS AppleDouble junk files (._filename) ----
% These hidden resource-fork files are auto-created by macOS when files
% are copied to/from network shares, and can make dir('*.ext')-based
% loaders (loadTrig, loadTrialParams, etc.) match more than one file,
% causing cryptic errors like "Too many input arguments" on a simple "/".
removeAppleDoubleFiles(resultsFolder);
removeAppleDoubleFiles(rawDataFolder);

%% ---- Load trigger + trial parameters from the RESULTS folder ----
origDir = pwd;
dirCleanup = onCleanup(@() cd(origDir)); %#ok<NASGU>
cd(resultsFolder);

trig = loadTrig(0);
TrialParams_raw = loadTrialParams;
loadStimParams; % script -- creates StimParams in this workspace

cd(origDir);

%% ---- Reconstruct per-trial event times (same logic as blanking fn) ----
TrialParams_condition_all = cell2mat(TrialParams_raw(:,2))';
num_elect = min(diff(find(diff(TrialParams_condition_all)~=0)));
TrialParams_cond = TrialParams_condition_all(1:num_elect:end);

nTrialsFromStimParams = floor((size(StimParams,1)-1) / num_elect);

if nTrialsFromStimParams ~= numel(TrialParams_cond)
    warning('estimate_artifact_blank_durations:trialMismatch', ...
        ['nTrialsFromStimParams (%d) does not match TrialParams trial ' ...
         'count (%d). Using the smaller of the two -- check your ' ...
         'StimParams/TrialParams files if this is unexpected.'], ...
        nTrialsFromStimParams, numel(TrialParams_cond));
end
nTrials = min(nTrialsFromStimParams, numel(TrialParams_cond));

if numel(trig) ~= nTrials
    warning('estimate_artifact_blank_durations:trigMismatch', ...
        ['Number of triggers (%d) does not match number of trials (%d). ' ...
         'Using the smaller of the two, in order.'], numel(trig), nTrials);
end
nTrials = min(nTrials, numel(trig));

EventTimes_perTrial = cell(1, nTrials);
IsZeroTrial   = false(1, nTrials);
nActiveElec   = zeros(1, nTrials);
maxAmplitude  = nan(1, nTrials);
conditionID   = TrialParams_cond(1:nTrials);

for tr = 1:nTrials
    rows_this_trial = (2 + (tr-1)*num_elect) : (1 + tr*num_elect);
    amp_this_trial = cell2mat(StimParams(rows_this_trial,16));
    active_rows = amp_this_trial > 0;

    if ~any(active_rows)
        EventTimes_perTrial{tr} = [];
        IsZeroTrial(tr) = true;
        continue;
    end

    nActiveElec(tr)  = sum(active_rows);
    maxAmplitude(tr) = max(amp_this_trial(active_rows));

    PTD_this_trial    = cell2mat(StimParams(rows_this_trial,6));
    nPulse_this_trial = cell2mat(StimParams(rows_this_trial,8));
    period_this_trial = cell2mat(StimParams(rows_this_trial,9));

    eventTimes_us = [];
    for r = find(active_rows(:))'
        thisPTD    = PTD_this_trial(r);
        thisNPulse = nPulse_this_trial(r);
        thisPeriod = period_this_trial(r);
        if thisNPulse > 0
            eventTimes_us = [eventTimes_us, thisPTD + (0:thisNPulse-1) * thisPeriod]; %#ok<AGROW>
        end
    end
    EventTimes_perTrial{tr} = unique(sort(eventTimes_us));
end

%% ---- Open RAW amplifier.dat on the shared drive ----
rawFile = fullfile(rawDataFolder, 'amplifier.dat');
fid = fopen(rawFile, 'r');
if fid == -1
    error('estimate_artifact_blank_durations:noFile', ...
        'Could not open raw file: %s', rawFile);
end
fidCleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fileInfo = dir(rawFile);
nTotalSamples = fileInfo.bytes / (bytesPerSample * nChn);

%% ---- Per-trial, per-channel blank-duration estimation ----
BIN_b            = -FS/1000;  % same pre-event margin as blanking code (~1ms)
range_max        = 150;       % uV, same recovery threshold
searchWindow      = 46;       % samples
maxSearchSamples  = 180;      % same cap as blanking code (~6ms search cap)
postBuffer        = 10;       % same buffer added after search
marginSamples     = 400;      % extra safety margin for the raw read window

maxRows = nTrials * nChn;
trialCol = nan(maxRows,1); chanCol = nan(maxRows,1); condCol = nan(maxRows,1);
nActiveCol = nan(maxRows,1); ampCol = nan(maxRows,1); nEventsCol = nan(maxRows,1);
recCol = nan(maxRows,1); totalSpanCol = nan(maxRows,1); eventSpanCol = nan(maxRows,1);
idx = 0;

progressStep = max(1, floor(nTrials/20));
tic;

for t = 1:nTrials
    if IsZeroTrial(t) || isempty(EventTimes_perTrial{t})
        continue;
    end

    if mod(t, progressStep) == 0
        fprintf('Trial %d / %d  (%.1f s elapsed)\n', t, nTrials, toc);
    end

    thisTrig = double(trig(t));
    eventTimes_us = EventTimes_perTrial{t};
    firstEvent_us = min(eventTimes_us);
    lastEvent_us  = max(eventTimes_us);
    firstEvent_offset = round(firstEvent_us * FS / 1e6);
    lastEvent_offset  = round(lastEvent_us  * FS / 1e6);

    winStartSample = thisTrig + firstEvent_offset + BIN_b;
    winEndSample   = thisTrig + lastEvent_offset + maxSearchSamples + searchWindow + postBuffer + marginSamples;

    if winStartSample < 1 || winEndSample > nTotalSamples
        warning('estimate_artifact_blank_durations:outOfBounds', ...
            'Trial %d: raw-data window out of file bounds, skipping.', t);
        continue;
    end

    nSamplesToRead = winEndSample - winStartSample + 1;
    byteOffset = (winStartSample - 1) * nChn * bytesPerSample;

    fseek(fid, byteOffset, 'bof');
    rawChunk = fread(fid, [nChn, nSamplesToRead], 'int16') .* scaleFactor;

    if size(rawChunk,2) < nSamplesToRead
        warning('estimate_artifact_blank_durations:shortRead', ...
            'Trial %d: short read near end of file, skipping.', t);
        continue;
    end

    localLastEvent  = thisTrig + lastEvent_offset  - winStartSample + 1;
    localFirstEvent = thisTrig + firstEvent_offset - winStartSample + 1;
    localBlankStart = localFirstEvent + BIN_b;

    for c = 1:nChn
        try
            ra = [1 searchWindow];
            while range(rawChunk(c, localLastEvent+ra(1):localLastEvent+ra(2))) > range_max
                ra = ra + 1;
                if ra(1) > maxSearchSamples
                    ra(1) = maxSearchSamples;
                    break;
                end
            end
            ra(1) = ra(1) + postBuffer;
            if ra(1) > maxSearchSamples
                ra(1) = maxSearchSamples;
            end

            localBlankEnd = localLastEvent + ra(1);

            % Recovery-only duration: measured from the LAST stimulation
            % event to the point signal recovered. This is the true
            % adaptive artifact-recovery quantity (what should scale with
            % current amplitude), independent of PTD/electrode spacing.
            recoveryDurationMs = (localBlankEnd - localLastEvent) / FS * 1000;

            % Total blanked span: from just before the FIRST event to
            % recovery after the LAST event. For sequential-stim trials
            % this also includes the fixed inter-electrode PTD spacing,
            % so it's NOT a pure artifact-recovery measure -- kept here
            % for reference (e.g. total signal loss / data completeness).
            totalBlankSpanMs = (localBlankEnd - localBlankStart) / FS * 1000;

            % Gap between first and last stimulation event in this trial
            % (driven by PTD/pulse-train design, not by the artifact
            % itself) -- useful to sanity-check where the inflation in
            % the old single metric was coming from.
            eventSpanMs = (lastEvent_us - firstEvent_us) / 1000;

            idx = idx + 1;
            trialCol(idx)      = t;
            chanCol(idx)       = c;
            condCol(idx)       = conditionID(t);
            nActiveCol(idx)    = nActiveElec(t);
            ampCol(idx)        = maxAmplitude(t);
            nEventsCol(idx)    = numel(eventTimes_us);
            recCol(idx)        = recoveryDurationMs;
            totalSpanCol(idx)  = totalBlankSpanMs;
            eventSpanCol(idx)  = eventSpanMs;
        catch ME
            warning('estimate_artifact_blank_durations:channelError', ...
                'Trial %d, channel %d failed: %s', t, c, ME.message);
        end
    end
end

keep = 1:idx;
nRows = numel(keep);
animalCol  = repmat({animalID}, nRows, 1);
sessionCol = repmat({sessionName}, nRows, 1);

summaryTable = table(animalCol, sessionCol, trialCol(keep), chanCol(keep), condCol(keep), ...
    nActiveCol(keep), ampCol(keep), nEventsCol(keep), ...
    recCol(keep), totalSpanCol(keep), eventSpanCol(keep), ...
    'VariableNames', {'AnimalID','SessionName','Trial','Channel','ConditionID','nActiveElectrodes', ...
                       'MaxAmplitude_uA','nEvents', ...
                       'RecoveryDuration_ms','TotalBlankSpan_ms','EventSpan_ms'});

save(outputMatFile, 'summaryTable');
fprintf('Done. %d rows saved to: %s\n', height(summaryTable), outputMatFile);

% Quick sanity-check summary, grouped by condition (edit grouping as
% needed once you know which ConditionIDs are simultaneous vs sequential)
% RecoveryDuration_ms is the metric that should be reported for the
% reviewer (true adaptive artifact-recovery time, from last stim event
% to signal recovery). TotalBlankSpan_ms/EventSpan_ms are kept for
% reference so you can see how much of the old inflated number was just
% PTD/electrode spacing.
try
    g = groupsummary(summaryTable, 'ConditionID', {'mean','max'}, ...
        {'RecoveryDuration_ms','TotalBlankSpan_ms','EventSpan_ms'});
    disp(g);
catch
    % groupsummary not available in older MATLAB -- skip
end

%% ===================== Local functions =====================
function removeAppleDoubleFiles(folderPath)
% Deletes hidden "._*" AppleDouble resource-fork files that macOS
% creates when copying to non-native filesystems (e.g. SMB shares).
% Safe to delete -- they hold no scientific data, only macOS metadata.
staleFiles = dir(fullfile(folderPath, '._*'));
for k = 1:numel(staleFiles)
    try
        delete(fullfile(staleFiles(k).folder, staleFiles(k).name));
    catch ME
        warning('Could not delete %s: %s', staleFiles(k).name, ME.message);
    end
end
end