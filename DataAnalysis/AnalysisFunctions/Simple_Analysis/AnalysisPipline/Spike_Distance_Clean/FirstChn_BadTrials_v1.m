%% ========================================================================
% BAD-TRIAL ENTRY FOR Ch1 SINGLE DATASETS ONLY - v1
%
% PURPOSE
%   Record the final bad trials of the Ch1 single-pulse condition(s) for all
%   single datasets of one animal in ONE place (distance_bad_trials below),
%   using the relative trial numbers printed / plotted by
%   Clean_Ch1_Single_v1.m. Only the Ch1 set is touched; every other set /
%   condition already saved in a dataset's _BadTrials.mat is left alone.
%
% HOW TO FILL distance_bad_trials
%   One entry per single dataset (one per distance block, e.g. 100-400 um,
%   500-700 um). For each entry:
%     matFilePath : the dataset's <base_name>_TrialSpikeCounts.mat
%     bad         : one row per amplitude, {Amp_uA, [relative bad trial numbers]}
%   Example:
%     distance_bad_trials(1).matFilePath = '.../Linearity_S1_100_400um_Single1_TrialSpikeCounts.mat';
%     distance_bad_trials(1).bad = {
%         3,  [2 9];
%         5,  [4 11 12];
%         10, [];          % [] = Ch1 at 10 uA has NO bad trials
%     };
%   Each listed amplitude REPLACES whatever was saved before for Ch1 at that
%   amplitude (it is the final list, not a patch). Amplitudes you do not list
%   are left untouched.
%
% OUTPUT
%   <base_name>_BadTrials.mat in each dataset folder, same format as
%   BadTrial_ManualConfirm_v1 (BadConfirmed, AllBadAbsoluteTrialIDs,
%   SourceTrialCountFile), so SpikeTRatio and the other scripts read it as
%   before. A timestamped backup of the previous file is kept.
% ========================================================================

clear; clc;

%% ============================ USER SETTINGS ===========================

Target_Channel = 1;              % this script is meant for Ch1 only
Condition_Tolerance = 1e-4;
Create_Backup_If_Output_Exists = true;

distance_bad_trials = struct('matFilePath',{},'bad',{});

distance_bad_trials(1).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1/Linearity_S1_100_400um_Single1_TrialSpikeCounts.mat';
distance_bad_trials(1).bad = {
  % Amp(uA)  relative bad trials
     3,      [1,2,6:8,16,19,20,22,24,26,27,29];
     5,      [1,2,6,7,10,18,26,27,28];
    10,      [7,8,18,20,22:24,26];
};

distance_bad_trials(2).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1/Linearity_S1_500_700um_Single1_TrialSpikeCounts.mat';
distance_bad_trials(2).bad = {
  % Amp(uA)  relative bad trials
     3,      [4,18,26:28,16:18];
     5,      [1,3,9,10,16,17];
    10,      [1:4,23,26,30,11,15];
};

%% =========================== APPLY EACH DATASET ========================

if isempty(distance_bad_trials)
    error('distance_bad_trials is empty.');
end

fprintf('\n============================================================\n');
fprintf('Ch%d SINGLE BAD TRIALS (v1) - %d dataset(s)\n',Target_Channel,numel(distance_bad_trials));
fprintf('============================================================\n');

for d = 1:numel(distance_bad_trials)

    matFilePath = distance_bad_trials(d).matFilePath;
    bad_rows    = distance_bad_trials(d).bad;

    fprintf('\n--- Dataset %d ---\n%s\n',d,matFilePath);

    if ~isfile(matFilePath)
        error('Dataset %d: file does not exist:\n%s',d,matFilePath);
    end
    if isempty(bad_rows) || size(bad_rows,2) ~= 2
        error('Dataset %d: bad must have 2 columns {Amp_uA, [trials]} and at least one row.',d);
    end

    LoadData = load(matFilePath,'SpikeCounts','TrialSummary','sim_stim');
    SpikeCounts  = LoadData.SpikeCounts;
    TrialSummary = LoadData.TrialSummary;
    sim_stim     = double(LoadData.sim_stim);

    if sim_stim ~= 1
        error('Dataset %d is not a single-pulse dataset (simultaneous_stim = %d).',d,sim_stim);
    end

    % ---- find the Ch1 set ----
    si = [];
    for ss = 1:numel(SpikeCounts.set)
        sc = SpikeCounts.set(ss).stim_channels;
        sc = sc(sc > 0);
        if numel(sc) == 1 && sc(1) == Target_Channel
            si = ss;
            break;
        end
    end
    if isempty(si)
        error('Dataset %d has no single-pulse set for Ch%d.',d,Target_Channel);
    end

    % ---- existing saved file ----
    [dataset_folder,count_file_name] = fileparts(matFilePath);
    base_name = regexprep(count_file_name,'_TrialSpikeCounts$','');
    output_path = fullfile(dataset_folder,sprintf('%s_BadTrials.mat',base_name));

    BadConfirmed = empty_confirmed_struct();
    if isfile(output_path)
        ExistingLoad = load(output_path,'BadConfirmed');
        if isfield(ExistingLoad,'BadConfirmed')
            BadConfirmed = ExistingLoad.BadConfirmed;
            fprintf('Existing saved file found (%d condition(s) already confirmed).\n',numel(BadConfirmed));
        else
            warning('Existing output file has no BadConfirmed - starting fresh.');
        end
    else
        fprintf('No existing saved file - this run creates it.\n');
    end

    % ---- apply each amplitude row ----
    amp_values = [SpikeCounts.set(si).amp.amp_value];

    for r = 1:size(bad_rows,1)

        amp_value = bad_rows{r,1};
        ai = find(abs(amp_values-amp_value) < Condition_Tolerance,1);
        if isempty(ai)
            error('Dataset %d, row %d: amp %g uA not found for Ch%d (available: %s).', ...
                d,r,amp_value,Target_Channel,num2str(amp_values));
        end
        pi = 1;     % single-pulse dataset: one PTD
        cond = SpikeCounts.set(si).amp(ai).ptd(pi);
        total_trials = cond.total_trials;

        final_trials = validate_trial_list(bad_rows{r,2},total_trials,d,r);

        condition_mask = TrialSummary.SetIndex == si & ...
            TrialSummary.AmplitudeIndex == ai & ...
            TrialSummary.PTDIndex == pi;

        absolute_trial_ids = zeros(size(final_trials));
        for k = 1:numel(final_trials)
            trial_mask = condition_mask & TrialSummary.ConditionTrialIndex == final_trials(k);
            if sum(trial_mask) ~= 1
                error('Dataset %d, row %d: condition trial %d could not be uniquely matched.',d,r,final_trials(k));
            end
            absolute_trial_ids(k) = TrialSummary.AbsoluteTrialID(trial_mask);
        end
        absolute_trial_ids = sort(unique(absolute_trial_ids));

        R = empty_confirmed_struct();
        R(1).set = si;
        R(1).stimChannels = SpikeCounts.set(si).stim_channels;
        R(1).stimOrderLabel = SpikeCounts.set(si).set_name;
        R(1).stimType = 'single';
        R(1).amp = amp_values(ai);
        R(1).ptd_ms = cond.PTD_ms;
        R(1).distance_um = SpikeCounts.set(si).distance_um;
        R(1).distance_note = SpikeCounts.set(si).distance_note;
        R(1).n_total_trials = total_trials;
        R(1).condition_trial_indices = final_trials;
        R(1).absolute_trial_ids = absolute_trial_ids;
        R(1).good_condition_trial_indices = setdiff(1:total_trials,final_trials);
        R(1).n_bad = numel(final_trials);
        R(1).confirmed_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

        match_idx = [];
        for k = 1:numel(BadConfirmed)
            if BadConfirmed(k).set == si && ...
                    abs(BadConfirmed(k).amp-amp_values(ai)) < Condition_Tolerance && ...
                    abs(BadConfirmed(k).ptd_ms-cond.PTD_ms) < Condition_Tolerance
                match_idx = k;
                break;
            end
        end
        if isempty(match_idx)
            BadConfirmed(end+1) = R; %#ok<SAGROW>
        else
            BadConfirmed(match_idx) = R;
        end

        fprintf('  Ch%d | Amp %g uA | bad %d/%d, relative: %s\n', ...
            Target_Channel,amp_values(ai),numel(final_trials),total_trials,number_list(final_trials));
    end

    % ---- all-bad list + save ----
    AllBadAbsoluteTrialIDs = [];
    for k = 1:numel(BadConfirmed)
        AllBadAbsoluteTrialIDs = union(AllBadAbsoluteTrialIDs,BadConfirmed(k).absolute_trial_ids(:).');
    end
    AllBadAbsoluteTrialIDs = sort(unique(AllBadAbsoluteTrialIDs(:).'));

    if isfile(output_path) && Create_Backup_If_Output_Exists
        backup_path = fullfile(dataset_folder,sprintf('%s_BadTrials_BACKUP_%s.mat', ...
            base_name,datestr(now,'yyyymmdd_HHMMSS')));
        copyfile(output_path,backup_path);
        fprintf('Backup of previous file: %s\n',backup_path);
    end

    SourceTrialCountFile = matFilePath;
    save(output_path,'BadConfirmed','AllBadAbsoluteTrialIDs','SourceTrialCountFile','-v7.3');

    fprintf('Saved: %s\n',output_path);
    fprintf('Conditions confirmed in this file: %d | total bad trials (all sets): %d\n', ...
        numel(BadConfirmed),numel(AllBadAbsoluteTrialIDs));
end

fprintf('\nDone. Re-run Clean_Ch1_Single_v1 to check what is left.\n');

%% =========================== LOCAL FUNCTIONS ==========================

function trials = validate_trial_list(trials,total_trials,dataset_index,row_index)
if isempty(trials)
    trials = [];
    return;
end
trials = unique(double(trials(:).'));
invalid = ~isfinite(trials) | trials < 1 | trials > total_trials | fix(trials) ~= trials;
if any(invalid)
    error('Dataset %d, row %d: invalid relative trial numbers (must be 1-%d): %s', ...
        dataset_index,row_index,total_trials,num2str(trials(invalid)));
end
end

function R = empty_confirmed_struct()
R = struct('set',{},'stimChannels',{},'stimOrderLabel',{},'stimType',{}, ...
    'amp',{},'ptd_ms',{},'distance_um',{},'distance_note',{},'n_total_trials',{}, ...
    'condition_trial_indices',{},'absolute_trial_ids',{},'good_condition_trial_indices',{}, ...
    'n_bad',{},'confirmed_at',{});
end

function output = number_list(values)
if isempty(values)
    output = '(none)';
else
    output = strtrim(sprintf('%g ',values));
end
end