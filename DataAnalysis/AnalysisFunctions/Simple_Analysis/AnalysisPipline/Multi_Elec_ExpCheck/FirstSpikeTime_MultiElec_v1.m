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
data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

win_ms = 10;       % search window length after final artifact, in ms
FS = 30000;        % sampling rate

%% ====================== LOAD FOLDER ======================
cd(data_folder);

%% ====================== LOAD SPIKES ======================
% This should contain sp_clipped.
sp_file = dir('*sp_xia.mat');
assert(~isempty(sp_file), 'No *sp_xia.mat file found.');

load(sp_file(1).name, 'sp_clipped');
fprintf('Loaded spike file: %s\n', sp_file(1).name);

%% ====================== LOAD TRIGGERS ======================
trig = loadTrig(0);

%% ====================== LOAD EXPERIMENT PARAMETERS ======================
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat file found.');

% Load full file so new metadata variables are available.
S = load(param_file(1).name);

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;   % rows/slots per trial
n_Trials          = S.n_Trials;

%% ====================== EXTRACT AMPLITUDE PER TRIAL ======================
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
trialAmps(trialAmps == -1) = 0;

%% ====================== LOAD MIXED-PREFIX METADATA ======================
% New design file should contain these variables.
% If not found, stop because recovery should only be done on new files.
requiredVars = {'active_electrode_count_by_trial', ...
                'prefix_length_by_trial', ...
                'isi_ms_by_trial', ...
                'condition_type_by_trial', ...
                'condition_set_id_by_trial'};

for i = 1:numel(requiredVars)
    assert(isfield(S, requiredVars{i}), ...
        'Missing metadata variable "%s". This code requires the new mixed-prefix file.', ...
        requiredVars{i});
end

activeCount_trial    = S.active_electrode_count_by_trial(:);
prefixLength_trial   = S.prefix_length_by_trial(:);
isi_ms_trial         = S.isi_ms_by_trial(:);
conditionType_trial  = S.condition_type_by_trial(:);
conditionSetID_trial = S.condition_set_id_by_trial(:);

%% ====================== CALCULATE LAST ACTIVE PTD ======================
% lastActivePTD_ms(tr) = time of the last active stimulation artifact
% for trial tr.
%
% Example:
%   Prefix 1, ISI 5 ms -> lastActivePTD_ms = 0
%   Prefix 3, ISI 5 ms -> lastActivePTD_ms = 10
%   Prefix 5, ISI 5 ms -> lastActivePTD_ms = 20

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials
    activeCount_this = activeCount_trial(tr);

    if isnan(activeCount_this) || activeCount_this < 1
        lastActivePTD_us(tr) = 0;
    else
        activeCount_this = min(round(activeCount_this), simultaneous_stim);

        % StimParams has a header row.
        % Data row for trial tr, active electrode slot activeCount_this:
        stimRow = 1 + (tr-1)*simultaneous_stim + activeCount_this;

        lastActivePTD_us(tr) = StimParams{stimRow,6};
    end
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

%% ====================== IDENTIFY TRIAL TYPES ======================
% Sequential prefix trials are used for first-spike detection.
isPrefixTrial = conditionType_trial == 1;

% Prefix 1 is normally a source condition and does not need recovery,
% but we still allow first-spike detection for all prefix trials.
% The recovery code will decide which trials are used.
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

fprintf('Detected last active PTDs (ms): ');
disp(unique(lastActivePTD_ms)');

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

    for tr = 1:n_Trials

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
        ch, sum(hasSpike(ch,isPrefixTrial)), sum(isPrefixTrial));
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
     '-v7.3');

fprintf('Saved to: %s\n', fullfile(data_folder, save_name));