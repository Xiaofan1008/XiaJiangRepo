%% ============================================================
% SpikeFiltering_PrefixRecovery_Cleanup.m
%
% Purpose:
%   Clean spike waveforms after prefix-based spike recovery.
%
% Input:
%   *.sp_xia_PrefixRecovery.mat
%   variable: sp_seq
%
% Output:
%   *.sp_xia_PrefixRecovery_SSD.mat
%   final cleaned variable: sp_corr
%
% Filtering steps:
%   1) Start-slope / rapid swing filter
%   2) Morphology width filter
%   3) Zero-crossing check
%   4) Template correlation filter
%
% Important for new mixed-prefix file:
%   The evoked template window is based on the final active artifact time:
%
%       template window = lastActivePTD_ms + evoked_template_after_PTD_ms
%
%   Example for ISI = 5 ms:
%       Prefix 1: last PTD = 0 ms  -> template window 1–10 ms
%       Prefix 3: last PTD = 10 ms -> template window 11–20 ms
%       Prefix 5: last PTD = 20 ms -> template window 21–30 ms
%
% ============================================================

clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';
FS = 30000;

%% ================= FILTERING PARAMETERS =================

% 1. Start-slope filter:
% Removes waveforms with an unrealistically fast voltage swing near the
% beginning of the waveform.
Start_Slope_Thresh_uV = 100;     % max allowed voltage swing in sliding 0.2 ms window

% 2. Morphology width filter:
% Removes waveforms with impossible trough-to-peak timing.
min_trough_peak_ms = 0.15;       % too fast = likely artifact
max_trough_peak_ms = 1.0;        % too slow = likely noise

% 3. Template correlation filter:
corr_thresh = 0.6;               % lower = more permissive, higher = stricter
do_corr_filter = 1;

%% ================= WINDOWS =================

% Baseline window: Used as fallback template if there are not enough evoked spikes.
baseline_window_ms = [-100 -5];

% Evoked template window after final active PTD:
% For each trial:
%   evoked template window =
%       [lastActivePTD_ms + 1, lastActivePTD_ms + 10]
%
% This avoids using spikes inside the artifact/recovery region for long
% prefixes.
evoked_template_after_PTD_ms = [1 30];

% Cleanup window:
% This is the time range where the correlation filter is applied
% Example:[0 40] means only remove low-correlation spikes from 0–40 ms after trigger.
cleanup_window_ms = [0 60];

%% ================= FOLDER & BASE NAME =================
if ~isfolder(data_folder)
    error('Folder does not exist.');
end
cd(data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');

if numel(underscores) > 4
    base_name = last_folder(1 : underscores(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Data folder:\n%s\n', data_folder);
fprintf('Base name: %s\n', base_name);

%% ================= LOAD PREFIX-RECOVERED SPIKES =================
% New input file from prefix recovery code.
fname_sp = [base_name '.sp_xia_PrefixRecovery.mat'];

assert(isfile(fname_sp), 'Cannot find %s. Please run prefix recovery first.', fname_sp);

S_in = load(fname_sp);

if isfield(S_in, 'sp_seq')
    sp_in = S_in.sp_seq;
else
    error('Variable "sp_seq" not found in %s.', fname_sp);
end

nCh = numel(sp_in);
fprintf('\nLoaded recovered spike file: %s\n', fname_sp);
fprintf('Number of spike channels: %d\n', nCh);

%% ================= LOAD TRIGGERS =================
if isempty(dir('*.trig.dat'))
    cur = pwd;
    cleanTrig_sabquick;
    cd(cur);
end

trig = loadTrig(0);
nTrials = numel(trig);

fprintf('Number of triggers: %d\n', nTrials);

%% ================= LOAD EXPERIMENT PARAMETERS =================
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat found.');

% Load the full experiment file so that the new mixed-prefix metadata are available.
S_exp = load(param_file(1).name);

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;   % rows/slots per trial
n_Trials          = S_exp.n_Trials;

% Check trial number consistency.
if n_Trials ~= nTrials
    warning('n_Trials in exp file (%d) does not match trigger number (%d). Using nTrials = min of both.', ...
        n_Trials, nTrials);
end

nTrials_use = min(n_Trials, nTrials);

%% ================= LOAD MIXED-PREFIX METADATA =================
% These variables are required for the new prefix stimulation file.
requiredVars = {'active_electrode_count_by_trial', ...
                'prefix_length_by_trial', ...
                'isi_ms_by_trial', ...
                'condition_type_by_trial', ...
                'condition_set_id_by_trial'};

for i = 1:numel(requiredVars)
    assert(isfield(S_exp, requiredVars{i}), ...
        'Missing metadata variable "%s". This code requires the new mixed-prefix file.', ...
        requiredVars{i});
end

activeCount_trial    = S_exp.active_electrode_count_by_trial(:);
prefixLength_trial   = S_exp.prefix_length_by_trial(:);
isi_ms_trial         = S_exp.isi_ms_by_trial(:);
conditionType_trial  = S_exp.condition_type_by_trial(:);
conditionSetID_trial = S_exp.condition_set_id_by_trial(:);

%% ================= EXTRACT AMPLITUDE PER TRIAL =================
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;

%% ================= CALCULATE LAST ACTIVE PTD =================
% lastActivePTD_ms(tr) = timing of final active stimulation artifact
% relative to the trigger for trial tr.
%
% This is used to build the dynamic evoked template window.

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials
    activeCount_this = activeCount_trial(tr);

    if isnan(activeCount_this) || activeCount_this < 1
        lastActivePTD_us(tr) = 0;
    else
        activeCount_this = min(round(activeCount_this), simultaneous_stim);

        % StimParams has one header row.
        stimRow = 1 + (tr-1)*simultaneous_stim + activeCount_this;

        lastActivePTD_us(tr) = StimParams{stimRow,6};
    end
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

fprintf('\nDetected prefixes: ');
disp(unique(prefixLength_trial(conditionType_trial == 1))');

fprintf('Detected ISIs (ms): ');
disp(unique(isi_ms_trial(conditionType_trial == 1))');

fprintf('Detected condition types: ');
disp(unique(conditionType_trial)');

fprintf('Detected last active PTDs (ms): ');
disp(unique(lastActivePTD_ms)');

%% ================= 1) START SLOPE FILTER =================
% This removes spikes with very fast waveform swings near waveform onset.
% These are likely residual blanking/stimulation artifacts.

sp_slope = sp_in;

check_dur_ms = 0.8;       % scan first 0.8 ms of waveform
slide_win_ms = 0.2;       % use 0.2 ms sliding window

n_check = round(check_dur_ms * FS / 1000);
n_slide = round(slide_win_ms * FS / 1000);

fprintf('\nStart Slope Filtering (Swing > %d uV in %.1f ms window, scanning first %.1f ms)...\n', ...
    Start_Slope_Thresh_uV, slide_win_ms, check_dur_ms);

for ch = 1:nCh

    if isempty(sp_in{ch})
        continue;
    end

    wfs = sp_in{ch}(:,2:end);

    limit_idx = min(n_check, size(wfs, 2));

    if limit_idx > n_slide

        early_wfs = wfs(:, 1:limit_idx);
        max_swing = zeros(size(early_wfs,1), 1);

        for k = 1 : (limit_idx - n_slide + 1)
            win_data = early_wfs(:, k : k+n_slide-1);
            swing = max(win_data, [], 2) - min(win_data, [], 2);
            max_swing = max(max_swing, swing);
        end

        valid_mask = max_swing < Start_Slope_Thresh_uV;
        sp_slope{ch} = sp_in{ch}(valid_mask, :);

        removed = sum(~valid_mask);
        if removed > 0
            fprintf('Ch %2d: Removed %d artifact spikes by start-slope filter.\n', ch, removed);
        end

    else
        sp_slope{ch} = sp_in{ch};
    end
end

%% ================= 2) MORPHOLOGY FILTER =================
% This removes waveforms with unrealistic trough-to-peak width.

sp_morph = sp_slope;

fprintf('\nMorphology Filtering (trough-to-peak width %.2f to %.2f ms)...\n', ...
    min_trough_peak_ms, max_trough_peak_ms);

for ch = 1:nCh

    if isempty(sp_slope{ch})
        continue;
    end

    wfs = sp_slope{ch}(:,2:end);

    % Find absolute minimum and maximum positions.
    [~, min_idx] = min(wfs, [], 2);
    [~, max_idx] = max(wfs, [], 2);

    % Absolute time difference between trough and peak.
    width_ms = abs(max_idx - min_idx) / (FS / 1000);

    valid_idx = (width_ms >= min_trough_peak_ms) & ...
                (width_ms <= max_trough_peak_ms);

    sp_morph{ch} = sp_slope{ch}(valid_idx, :);

    removed = sum(~valid_idx);
    if removed > 0
        fprintf('Ch %2d: Removed %d spikes by morphology width filter.\n', ch, removed);
    end
end

%% ================= 3) ZERO-CROSSING CHECK =================
% This removes purely negative waveforms without a positive phase.

sp_zcross = sp_morph;

fprintf('\nZero-Crossing Check (removing waveforms without positive phase)...\n');

for ch = 1:nCh

    if isempty(sp_zcross{ch})
        continue;
    end

    wfs = sp_zcross{ch}(:,2:end);

    has_positive_phase = max(wfs, [], 2) > 0;

    removed = sum(~has_positive_phase);

    if removed > 0
        sp_zcross{ch} = sp_zcross{ch}(has_positive_phase, :);
        fprintf('Ch %2d: Removed %d spikes by zero-crossing check.\n', ch, removed);
    end
end

%% ================= 4) TEMPLATE CORRELATION FILTER =================
% This builds a waveform template for each channel.
%
% Template hierarchy:
%   1) Use evoked spikes first.
%   2) If not enough evoked spikes, use baseline spikes as fallback.
%
% For the new prefix file, the evoked template window is dynamic:
%
%   [lastActivePTD_ms + evoked_template_after_PTD_ms(1),
%    lastActivePTD_ms + evoked_template_after_PTD_ms(2)]
%
% This makes the template window safe for Prefix 1, 2, 3, 4, and 5.

sp_corr = sp_zcross;

if do_corr_filter

    fprintf('\nCorrelation Filtering (Evoked Template First, Baseline Fallback)...\n');

    for ch = 1:nCh

        if isempty(sp_corr{ch})
            continue;
        end

        spt = sp_corr{ch}(:,1);      % spike times in ms
        wfs = sp_corr{ch}(:,2:end);  % spike waveforms
        nSpikes = size(wfs,1);

        %% ----- A. Build dynamic evoked and baseline masks -----
        mask_evoked = false(nSpikes,1);
        mask_base   = false(nSpikes,1);

        for tr = 1:nTrials_use

            t0 = trig(tr) / FS * 1000;

            %% Evoked template window
            % Only use prefix/sequential and simultaneous stimulation trials,
            % not zero-control trials.
            %
            % conditionType:
            %   0 = zero-control
            %   1 = sequential prefix
            %   2 = simultaneous
            if conditionType_trial(tr) ~= 0

                evoked_start = lastActivePTD_ms(tr) + evoked_template_after_PTD_ms(1);
                evoked_end   = lastActivePTD_ms(tr) + evoked_template_after_PTD_ms(2);

                in_evoked = (spt >= t0 + evoked_start) & ...
                             (spt <= t0 + evoked_end);

                mask_evoked = mask_evoked | in_evoked;
            end

            %% Baseline template window
            in_base = (spt >= t0 + baseline_window_ms(1)) & ...
                      (spt <= t0 + baseline_window_ms(2));

            mask_base = mask_base | in_base;
        end

        wfs_evoked = wfs(mask_evoked,:);
        wfs_base   = wfs(mask_base,:);

        %% ----- B. Choose template -----
        % Keep your preferred hierarchy:
        %   evoked first, baseline fallback.
        if size(wfs_evoked, 1) >= 5
            template = mean(wfs_evoked, 1);
            source_str = 'Dynamic_Evoked';
        elseif size(wfs_base, 1) >= 5
            template = mean(wfs_base, 1);
            source_str = 'Baseline_Fallback';
        else
            fprintf('Ch %2d: Not enough spikes for template. Skipping correlation filter.\n', ch);
            continue;
        end

        %% ----- C. Correlate each spike with the template -----
        corr_vals = zeros(nSpikes,1);

        for i = 1:nSpikes
            corr_vals(i) = corr(template(:), wfs(i,:)');
        end

        %% ----- D. Apply correlation cleanup only inside cleanup_window_ms -----
        % This means low-correlation spikes are only removed if they occur in
        % the early evoked/artifact-related time window.
        bad_mask = false(nSpikes,1);

        for tr = 1:nTrials_use

            t0 = trig(tr) / FS * 1000;

            in_cleanup_window = (spt >= t0 + cleanup_window_ms(1)) & ...
                                (spt <= t0 + cleanup_window_ms(2));

            if any(in_cleanup_window)
                bad_mask(in_cleanup_window & (corr_vals < corr_thresh)) = true;
            end
        end

        removed_count = sum(bad_mask);
        sp_corr{ch} = sp_corr{ch}(~bad_mask, :);

        fprintf('Ch %2d (%s template): Removed %d spikes by correlation filter.\n', ...
            ch, source_str, removed_count);
    end
end

%% ================= SAVE OUTPUT =================
QC_params = struct();

QC_params.input_file = fname_sp;
QC_params.output_type = 'PrefixRecovery_SSD';

QC_params.Start_Slope_Thresh_uV = Start_Slope_Thresh_uV;
QC_params.check_dur_ms = check_dur_ms;
QC_params.slide_win_ms = slide_win_ms;

QC_params.min_trough_peak_ms = min_trough_peak_ms;
QC_params.max_trough_peak_ms = max_trough_peak_ms;

QC_params.baseline_window_ms = baseline_window_ms;
QC_params.evoked_template_after_PTD_ms = evoked_template_after_PTD_ms;
QC_params.cleanup_window_ms = cleanup_window_ms;

QC_params.corr_thresh = corr_thresh;
QC_params.do_corr_filter = do_corr_filter;

QC_params.FS = FS;

out_name = [base_name '.sp_xia_PrefixRecovery_SSD.mat'];

if isfile(out_name)
    delete(out_name);
end

save(out_name, ...
    'sp_in', ...
    'sp_slope', ...
    'sp_morph', ...
    'sp_zcross', ...
    'sp_corr', ...
    'QC_params', ...
    'lastActivePTD_ms', ...
    'activeCount_trial', ...
    'prefixLength_trial', ...
    'isi_ms_trial', ...
    'conditionType_trial', ...
    'conditionSetID_trial', ...
    'trialAmps', ...
    '-v7.3');

fprintf('\nSaved cleaned prefix-recovered spikes to:\n%s\n', fullfile(data_folder, out_name));