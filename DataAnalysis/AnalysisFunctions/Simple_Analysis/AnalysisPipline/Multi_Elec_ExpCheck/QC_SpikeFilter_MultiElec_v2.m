%% ============================================================
% SpikeFiltering_PrefixRecovery_Cleanup.m
% Purpose:
%   Clean spike waveforms after prefix-based spike recovery.
% Input:
%   *.sp_xia_PrefixRecovery.mat
%   variable: sp_seq
% Output:
%   *.sp_xia_PrefixRecovery_SSD.mat
%   final cleaned variable: sp_corr
% Filtering steps:
%   1) Start-slope / rapid swing filter
%   2) Morphology width filter
%   3) Zero-crossing check
%   4) Channel-level template correlation filter
% ============================================================

clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

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

% Minimum number of candidate-clean spikes needed to build channel template.
min_template_spikes = 10;

%% ================= WINDOWS =================

% Cleanup window:
% Correlation filter is applied only inside this time range relative to each trigger.
%
% Important:
%   Keep this the same for all active-count conditions.
%   This makes later Prefix 1/2/3/4/5 comparisons fairer.
cleanup_window_ms = [0 60];

%% ================= FOLDER & BASE NAME =================
if ~isfolder(data_folder)
    error('Folder does not exist.');
end
cd(data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');

% More robust than >4 because some folders have exactly 4 underscores.
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Data folder:\n%s\n', data_folder);
fprintf('Base name: %s\n', base_name);

%% ================= LOAD PREFIX-RECOVERED SPIKES =================
% Safer than manually constructing the filename.
rec_file = dir('*sp_xia_PrefixRecovery.mat');
assert(~isempty(rec_file), 'Cannot find *sp_xia_PrefixRecovery.mat. Please run prefix recovery first.');

fname_sp = rec_file(1).name;
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
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ================= LOAD EXPERIMENT PARAMETERS =================
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat found.');

S_exp = load(param_file(1).name);

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;   % rows/slots per trial
n_Trials          = S_exp.n_Trials;

if n_Trials ~= nTrig
    warning('n_Trials in exp file (%d) does not match trigger number (%d). Using nTrials_use = min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ================= REMOVE HEADER ROW =================
StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;
if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

%% ================= TRIAL METADATA FROM STIMPARAMS =================
% Important:
%   Do NOT use separate metadata arrays directly.
%   They may not be randomized in the same order as StimParams/triggers.
%
% StimParams columns:
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

activeCount_trial    = activeCount_trial(:);
prefixLength_trial   = prefixLength_trial(:);
isi_ms_trial         = isi_ms_trial(:);
conditionType_trial  = conditionType_trial(:);
conditionSetID_trial = conditionSetID_trial(:);

fprintf('\nUsing trial metadata from StimParams columns 26–30.\n');

%% ================= EXTRACT AMPLITUDE PER TRIAL =================
trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

%% ================= CALCULATE FINAL ACTIVE ARTIFACT TIME =================
% lastActivePTD_ms(tr) = final active artifact time for trial tr.
%
% For each active row:
%   final artifact time =
%       PTD_us + (PulseNum - 1) * PulsePeriod_us
%
% Then for each trial:
%   lastActivePTD_us(tr) = max(final artifact time across active rows)

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials

    if conditionType_trial(tr) == 0
        lastActivePTD_us(tr) = 0;
        continue;
    end

    if conditionType_trial(tr) == 1
        nActive_this = prefixLength_trial(tr);
    elseif conditionType_trial(tr) == 2
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

    % Column 6: PTD in us.
    ptd_us = cell2mat(StimParams_data(activeRows,6));

    % Column 8: pulse train number.
    pulseNum = cell2mat(StimParams_data(activeRows,8));

    % Column 9: pulse train period in us.
    pulsePeriod_us = cell2mat(StimParams_data(activeRows,9));

    pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
    pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

    rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowFinalArtifact_us);
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

fprintf('\nDetected prefixes: ');
disp(unique(prefixLength_trial(conditionType_trial == 1))');

fprintf('Detected ISIs (ms): ');
disp(unique(isi_ms_trial(conditionType_trial == 1))');

fprintf('Detected condition types: ');
disp(unique(conditionType_trial)');

fprintf('Detected final active artifact times / last PTDs (ms): ');
disp(unique(lastActivePTD_ms)');

%% ================= 1) START SLOPE FILTER =================
sp_slope = sp_in;

check_dur_ms = 0.8;       % scan first 0.8 ms of waveform
slide_win_ms = 0.2;       % use 0.2 ms sliding window

n_check = round(check_dur_ms * FS / 1000);
n_slide = round(slide_win_ms * FS / 1000);

fprintf('\nStart Slope Filtering (Swing > %d uV in %.1f ms window, scanning first %.1f ms)...\n', ...
    Start_Slope_Thresh_uV, slide_win_ms, check_dur_ms);

removed_slope = zeros(nCh,1);

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
        removed_slope(ch) = removed;

        if removed > 0
            fprintf('Ch %2d: Removed %d artifact spikes by start-slope filter.\n', ch, removed);
        end

    else
        sp_slope{ch} = sp_in{ch};
    end
end

%% ================= 2) MORPHOLOGY FILTER =================
sp_morph = sp_slope;

fprintf('\nMorphology Filtering (trough-to-peak width %.2f to %.2f ms)...\n', ...
    min_trough_peak_ms, max_trough_peak_ms);

removed_morph = zeros(nCh,1);

for ch = 1:nCh

    if isempty(sp_slope{ch})
        continue;
    end

    wfs = sp_slope{ch}(:,2:end);

    [~, min_idx] = min(wfs, [], 2);
    [~, max_idx] = max(wfs, [], 2);

    width_ms = abs(max_idx - min_idx) / (FS / 1000);

    valid_idx = (width_ms >= min_trough_peak_ms) & ...
                (width_ms <= max_trough_peak_ms);

    sp_morph{ch} = sp_slope{ch}(valid_idx, :);

    removed = sum(~valid_idx);
    removed_morph(ch) = removed;

    if removed > 0
        fprintf('Ch %2d: Removed %d spikes by morphology width filter.\n', ch, removed);
    end
end

%% ================= 3) ZERO-CROSSING CHECK =================
sp_zcross = sp_morph;

fprintf('\nZero-Crossing Check (removing waveforms without positive phase)...\n');

removed_zcross = zeros(nCh,1);

for ch = 1:nCh

    if isempty(sp_zcross{ch})
        continue;
    end

    wfs = sp_zcross{ch}(:,2:end);

    has_positive_phase = max(wfs, [], 2) > 0;

    removed = sum(~has_positive_phase);
    removed_zcross(ch) = removed;

    if removed > 0
        sp_zcross{ch} = sp_zcross{ch}(has_positive_phase, :);
        fprintf('Ch %2d: Removed %d spikes by zero-crossing check.\n', ch, removed);
    end
end

%% ================= 4) CHANNEL-LEVEL TEMPLATE CORRELATION FILTER =================
% This builds one channel-level waveform template from all remaining
% candidate-clean spikes after:
%   slope -> morphology -> zero-crossing
%
% The template does not care about:
%   prefix length
%   active electrode count
%   condition type
%   stimulation set/order
%   spike timing
%
% Then the same correlation threshold is applied to spikes inside
% cleanup_window_ms around each trigger.

sp_corr = sp_zcross;

corr_template = cell(nCh,1);
corr_source_n = zeros(nCh,1);
removed_corr = zeros(nCh,1);

if do_corr_filter

    fprintf('\nChannel-Level Correlation Filtering...\n');
    fprintf('Template source: all candidate-clean spikes after slope/morphology/zero-crossing.\n');
    fprintf('Cleanup window: %.1f to %.1f ms relative to each trigger.\n', ...
        cleanup_window_ms(1), cleanup_window_ms(2));

    for ch = 1:nCh

        if isempty(sp_corr{ch})
            continue;
        end

        spt = sp_corr{ch}(:,1);      % spike times in ms
        wfs = sp_corr{ch}(:,2:end);  % spike waveforms
        nSpikes = size(wfs,1);

        %% ----- A. Build one channel-level template -----
        if nSpikes < min_template_spikes
            fprintf('Ch %2d: Only %d candidate-clean spikes. Skipping correlation filter.\n', ...
                ch, nSpikes);
            continue;
        end

        template = mean(wfs, 1);
        corr_template{ch} = template;
        corr_source_n(ch) = nSpikes;

        %% ----- B. Correlate each spike with the template -----
        corr_vals = zeros(nSpikes,1);

        for i = 1:nSpikes
            corr_vals(i) = corr(template(:), wfs(i,:)');
        end

        % Protect against flat waveforms or zero-variance template/waveform.
        corr_vals(isnan(corr_vals)) = -1;

        %% ----- C. Apply correlation cleanup only inside cleanup_window_ms -----
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
        removed_corr(ch) = removed_count;

        sp_corr{ch} = sp_corr{ch}(~bad_mask, :);

        fprintf('Ch %2d: Template N = %d | Removed %d spikes by correlation filter.\n', ...
            ch, nSpikes, removed_count);
    end
end

%% ================= REMOVAL SUMMARY BY CONDITION =================
% This summary is important for checking whether filtering removes more
% spikes from one prefix than another.
%
% Summary is based on final correlation stage:
%   before correlation = sp_zcross
%   after correlation  = sp_corr

SummaryRows = {};
row_i = 0;

condition_labels = {'ZeroControl', 'Prefix', 'Simultaneous'};

for condType = 0:2

    if condType == 0
        prefix_list = 0;
    elseif condType == 1
        prefix_list = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));
        prefix_list = sort(prefix_list(:))';
    elseif condType == 2
        prefix_list = 0;
    end

    for pp = 1:numel(prefix_list)

        prefix_val = prefix_list(pp);

        if condType == 0
            trial_list = find(conditionType_trial == 0);
            label = 'ZeroControl';
        elseif condType == 1
            trial_list = find(conditionType_trial == 1 & prefixLength_trial == prefix_val);
            label = sprintf('Prefix_%d', prefix_val);
        elseif condType == 2
            trial_list = find(conditionType_trial == 2);
            label = 'Simultaneous';
        end

        trial_list = trial_list(trial_list <= nTrials_use);

        before_count = 0;
        after_count  = 0;

        for ch = 1:nCh

            if isempty(sp_zcross{ch}) && isempty(sp_corr{ch})
                continue;
            end

            if ~isempty(sp_zcross{ch})
                spt_before = sp_zcross{ch}(:,1);
            else
                spt_before = [];
            end

            if ~isempty(sp_corr{ch})
                spt_after = sp_corr{ch}(:,1);
            else
                spt_after = [];
            end

            for tt = 1:numel(trial_list)

                tr = trial_list(tt);
                t0 = trig(tr) / FS * 1000;

                before_count = before_count + sum(spt_before >= t0 + cleanup_window_ms(1) & ...
                                                  spt_before <= t0 + cleanup_window_ms(2));

                after_count = after_count + sum(spt_after >= t0 + cleanup_window_ms(1) & ...
                                                spt_after <= t0 + cleanup_window_ms(2));
            end
        end

        removed_count = before_count - after_count;

        if before_count > 0
            removed_percent = removed_count / before_count * 100;
        else
            removed_percent = NaN;
        end

        row_i = row_i + 1;
        SummaryRows(row_i,1:6) = {label, condType, prefix_val, before_count, after_count, removed_percent};
    end
end

RemovalSummary = cell2table(SummaryRows, ...
    'VariableNames', {'ConditionLabel', 'ConditionType', 'PrefixLength', ...
                      'BeforeCorr', 'AfterCorr', 'RemovedPercent'});

fprintf('\n================ Removal Summary inside cleanup window ================\n');
disp(RemovalSummary);
fprintf('======================================================================\n');

%% ================= SAVE OUTPUT =================
QC_params = struct();

QC_params.input_file = fname_sp;
QC_params.output_type = 'PrefixRecovery_SSD';

QC_params.Start_Slope_Thresh_uV = Start_Slope_Thresh_uV;
QC_params.check_dur_ms = check_dur_ms;
QC_params.slide_win_ms = slide_win_ms;

QC_params.min_trough_peak_ms = min_trough_peak_ms;
QC_params.max_trough_peak_ms = max_trough_peak_ms;

QC_params.cleanup_window_ms = cleanup_window_ms;

QC_params.corr_thresh = corr_thresh;
QC_params.do_corr_filter = do_corr_filter;
QC_params.min_template_spikes = min_template_spikes;
QC_params.template_strategy = 'channel_level_all_candidate_clean_spikes';

QC_params.FS = FS;

out_name = [base_name '.sp_xia_SSD.mat'];

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
    'corr_template', ...
    'corr_source_n', ...
    'removed_slope', ...
    'removed_morph', ...
    'removed_zcross', ...
    'removed_corr', ...
    'RemovalSummary', ...
    'lastActivePTD_ms', ...
    'lastActivePTD_us', ...
    'activeCount_trial', ...
    'prefixLength_trial', ...
    'isi_ms_trial', ...
    'conditionType_trial', ...
    'conditionSetID_trial', ...
    'trialAmps', ...
    '-v7.3');

fprintf('\nSaved cleaned prefix-recovered spikes to:\n%s\n', fullfile(data_folder, out_name));