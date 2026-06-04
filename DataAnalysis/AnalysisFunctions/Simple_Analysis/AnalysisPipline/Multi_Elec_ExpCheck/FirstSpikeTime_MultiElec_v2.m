% ============================================================
%  Extract First Post-Blank Spike Time per Trial
%  For mixed-prefix stimulation files
%
%  Purpose:
%    For each recording channel and each sequential prefix trial,
%    find the first detected spike after the final active artifact.
%
%  Important:
%    - Prefix 1 has PTD = 0, but it is NOT simultaneous.
%    - Simultaneous trials are identified by conditionType_trial == 2.
%    - Sequential prefix trials are conditionType_trial == 1.
%    - Metadata are read directly from randomized StimParams columns 26–30
%      to avoid trial-order mismatch.
%
%  Output:
%    firstSpikeTimes{ch}(trial)
%    hasSpike(ch, trial)
%    lastActivePTD_ms(trial)
%    activeCount_trial
%    prefixLength_trial
%    isi_ms_trial
%    conditionType_trial
%    conditionSetID_trial
%    trialAmps
%
%  Saved file:
%    FirstSpikeTimes_Prefix.mat
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER INPUT ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

win_ms = 15;       % search window length after final artifact, in ms
FS = 30000;        % sampling rate

% If true, print one example trial per prefix to confirm metadata alignment.
debug_print_trial_content = true;

%% ====================== LOAD FOLDER ======================
if ~isfolder(data_folder)
    error('Data folder does not exist.');
end
cd(data_folder);

fprintf('\nData folder:\n%s\n', data_folder);

%% ====================== LOAD SPIKES ======================
% This should contain sp_clipped.
sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No *sp_xia.mat file found.');

load(sp_file(1).name, 'sp_clipped');
fprintf('Loaded spike file: %s\n', sp_file(1).name);

%% ====================== LOAD TRIGGERS ======================
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD EXPERIMENT PARAMETERS ======================
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat file found.');

S = load(param_file(1).name);

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;   % rows/slots per trial
n_Trials          = S.n_Trials;

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== REMOVE HEADER ROW ======================
StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;
if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

%% ====================== TRIAL METADATA FROM STIMPARAMS ======================
% Important:
%   Do NOT use the separate metadata arrays directly here, because they may
%   not be randomized in the same order as StimParams.
%
% Instead, derive metadata directly from randomized StimParams columns:
%   26 = ActiveElectrodeCount
%   27 = PrefixLength
%   28 = ISI_ms
%   29 = ConditionType
%   30 = ConditionSetID

if size(StimParams_data,2) < 30
    error('StimParams does not contain columns 26–30. Cannot use mixed-prefix metadata.');
end

firstRow_eachTrial = 1:simultaneous_stim:size(StimParams_data,1);

activeCount_trial    = cell2mat(StimParams_data(firstRow_eachTrial,26));
prefixLength_trial   = cell2mat(StimParams_data(firstRow_eachTrial,27));
isi_ms_trial         = cell2mat(StimParams_data(firstRow_eachTrial,28));
conditionType_trial  = cell2mat(StimParams_data(firstRow_eachTrial,29));
conditionSetID_trial = cell2mat(StimParams_data(firstRow_eachTrial,30));

% Force column vectors.
activeCount_trial    = activeCount_trial(:);
prefixLength_trial   = prefixLength_trial(:);
isi_ms_trial         = isi_ms_trial(:);
conditionType_trial  = conditionType_trial(:);
conditionSetID_trial = conditionSetID_trial(:);

fprintf('\nUsing trial metadata from StimParams columns 26–30.\n');

%% ====================== EXTRACT AMPLITUDE PER TRIAL ======================
trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

%% ====================== CALCULATE FINAL ACTIVE ARTIFACT TIME ======================
% lastActivePTD_ms(tr) = final active artifact time for trial tr.
%
% For each active row:
%   final artifact time =
%       PTD_us + (PulseNum - 1) * PulsePeriod_us
%
% Then for the trial:
%   lastActivePTD_us(tr) = max(final artifact time across active rows)
%
% For current prefix files with PulseNum = 1:
%   Prefix 1, ISI 5 ms -> 0 ms
%   Prefix 3, ISI 5 ms -> 10 ms
%   Prefix 5, ISI 5 ms -> 20 ms

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials

    if conditionType_trial(tr) == 0
        % Zero-control trial.
        lastActivePTD_us(tr) = 0;
        continue;
    end

    if conditionType_trial(tr) == 1
        % Sequential prefix trial.
        nActive_this = prefixLength_trial(tr);
    elseif conditionType_trial(tr) == 2
        % Full simultaneous trial.
        nActive_this = activeCount_trial(tr);
    else
        nActive_this = activeCount_trial(tr);
    end

    if isnan(nActive_this) || nActive_this < 1
        lastActivePTD_us(tr) = 0;
        continue;
    end

    nActive_this = min(round(nActive_this), simultaneous_stim);

    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);
    activeRows = rr(1:nActive_this);

    % Column 6: post-trigger delay, us.
    ptd_us = cell2mat(StimParams_data(activeRows,6));

    % Column 8: pulse train number.
    pulseNum = cell2mat(StimParams_data(activeRows,8));

    % Column 9: pulse train period, us.
    pulsePeriod_us = cell2mat(StimParams_data(activeRows,9));

    pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
    pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

    rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowFinalArtifact_us);
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

%% ====================== IDENTIFY TRIAL TYPES ======================
% Sequential prefix trials are used for first-spike detection.
isPrefixTrial = conditionType_trial == 1;

% Simultaneous and zero-control trials are saved but not searched.
isSimultaneous = conditionType_trial == 2;
isZeroControl  = conditionType_trial == 0;

fprintf('\nDetected prefixes: ');
disp(unique(prefixLength_trial(isPrefixTrial))');

fprintf('Detected ISIs (ms): ');
disp(unique(isi_ms_trial(isPrefixTrial))');

fprintf('Detected condition set IDs: ');
disp(unique(conditionSetID_trial(conditionSetID_trial > 0))');

fprintf('Detected amplitudes: ');
disp(unique(trialAmps)');

fprintf('Detected condition types: ');
disp(unique(conditionType_trial)');

fprintf('Detected final active artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================
% This helps confirm that Prefix 1/2/3/4/5 are aligned correctly with
% randomized StimParams.

if debug_print_trial_content

    prefix_all = unique(prefixLength_trial(isPrefixTrial & prefixLength_trial > 0));
    set_all = unique(conditionSetID_trial(conditionSetID_trial > 0));
    amp_all = unique(trialAmps(isPrefixTrial));
    amp_all = amp_all(amp_all > 0);
    isi_all = unique(isi_ms_trial(isPrefixTrial));

    if ~isempty(prefix_all) && ~isempty(set_all) && ~isempty(amp_all) && ~isempty(isi_all)

        debug_set = set_all(1);
        debug_amp = amp_all(end);   % usually higher amp is easier to check
        debug_isi = isi_all(1);

        fprintf('\n================ DEBUG TRIAL CONTENT CHECK ================\n');
        fprintf('Debug Set = %d | Amp = %.1f uA | ISI = %.1f ms\n', ...
            debug_set, debug_amp, debug_isi);

        for ip = 1:numel(prefix_all)

            prefix_val = prefix_all(ip);

            tr_debug = find(conditionSetID_trial == debug_set & ...
                            conditionType_trial == 1 & ...
                            prefixLength_trial == prefix_val & ...
                            isi_ms_trial == debug_isi & ...
                            trialAmps == debug_amp, ...
                            1, 'first');

            if isempty(tr_debug)
                fprintf('\nPrefix %d: no matching trial found.\n', prefix_val);
                continue;
            end

            rr = (tr_debug-1)*simultaneous_stim + (1:simultaneous_stim);

            stimNames_debug = StimParams_data(rr,1);
            ptd_debug = cell2mat(StimParams_data(rr,6));
            pulseNum_debug = cell2mat(StimParams_data(rr,8));
            pulsePeriod_debug = cell2mat(StimParams_data(rr,9));
            amp_debug = cell2mat(StimParams_data(rr,16));
            activeCount_debug = cell2mat(StimParams_data(rr,26));
            prefix_debug = cell2mat(StimParams_data(rr,27));
            isi_debug = cell2mat(StimParams_data(rr,28));
            condType_debug = cell2mat(StimParams_data(rr,29));
            setID_debug = cell2mat(StimParams_data(rr,30));

            fprintf('\nPrefix %d | Trial %d\n', prefix_val, tr_debug);
            fprintf('  conditionType = %d, setID = %d, ISI = %.1f ms, amp = %.1f uA\n', ...
                conditionType_trial(tr_debug), conditionSetID_trial(tr_debug), ...
                isi_ms_trial(tr_debug), trialAmps(tr_debug));

            fprintf('  StimNames:       ');
            disp(stimNames_debug');

            fprintf('  PTD_us:          ');
            disp(ptd_debug');

            fprintf('  PulseNum:        ');
            disp(pulseNum_debug');

            fprintf('  PulsePeriod_us:  ');
            disp(pulsePeriod_debug');

            fprintf('  Amp_col16:       ');
            disp(amp_debug');

            fprintf('  ActiveCount_col26: ');
            disp(activeCount_debug');

            fprintf('  Prefix_col27:      ');
            disp(prefix_debug');

            fprintf('  ISI_col28:         ');
            disp(isi_debug');

            fprintf('  CondType_col29:    ');
            disp(condType_debug');

            fprintf('  SetID_col30:       ');
            disp(setID_debug');

            fprintf('  FinalArtifact_us = %.1f\n', lastActivePTD_us(tr_debug));
        end

        fprintf('===========================================================\n');
    end
end

%% ====================== INITIALIZE OUTPUTS ======================
nChn = numel(sp_clipped);

firstSpikeTimes = cell(nChn,1);
hasSpike = zeros(nChn, n_Trials);

fprintf('\nExtracting first post-blank spike times...\n');

%% ====================== MAIN LOOP ======================
for ch = 1:nChn

    S_ch = sp_clipped{ch};
    firstSpikeTimes{ch} = nan(n_Trials,1);

    if isempty(S_ch)
        continue;
    end

    spike_times_ms = S_ch(:,1);

    for tr = 1:nTrials_use

        % Only detect first post-blank spike for sequential prefix trials.
        % Do not use PTD == 0 to identify simultaneous trials, because
        % Prefix 1 also has PTD = 0.
        if ~isPrefixTrial(tr)
            continue;
        end

        % Trigger time in ms.
        t0 = trig(tr) / FS * 1000;

        % Search after the final active artifact.
        win_start = lastActivePTD_ms(tr);
        win_end   = lastActivePTD_ms(tr) + win_ms;

        rel_t = spike_times_ms - t0;

        spk_in_win = rel_t(rel_t >= win_start & rel_t <= win_end);

        if ~isempty(spk_in_win)
            firstSpikeTimes{ch}(tr) = spk_in_win(1);
            hasSpike(ch,tr) = 1;
        end
    end

    fprintf('Ch %2d: first spikes found in %d/%d prefix trials.\n', ...
        ch, sum(hasSpike(ch,isPrefixTrial)), sum(isPrefixTrial(1:nTrials_use)));
end

fprintf('==============================================================\n');

%% ====================== SAVE ======================
save_name = 'FirstSpikeTimes_Prefix.mat';

save(save_name, ...
     'firstSpikeTimes', ...
     'hasSpike', ...
     'lastActivePTD_ms', ...
     'lastActivePTD_us', ...
     'activeCount_trial', ...
     'prefixLength_trial', ...
     'isi_ms_trial', ...
     'conditionType_trial', ...
     'conditionSetID_trial', ...
     'trialAmps', ...
     'isPrefixTrial', ...
     'isSimultaneous', ...
     'isZeroControl', ...
     'win_ms', ...
     'FS', ...
     'n_Trials', ...
     'nTrials_use', ...
     '-v7.3');

fprintf('Saved to: %s\n', fullfile(data_folder, save_name));