%% ========================================================================
% BAD-TRIAL MANUAL CONFIRMATION
%
% PURPOSE
%   You've looked at BadTrial_ReferenceCriteria_v1.m's printed flagged-
%   trial list (and/or your own inspection) for one or more conditions,
%   and decided on the FINAL bad-trial list for each. This script records
%   that decision. Bad trials are GLOBAL (excluded from the whole dataset,
%   not per channel), matching the responding-channel confirmation tool's
%   Overrides pattern.
%
%   Edit Overrides below (one row per condition you're confirming this
%   run - you don't have to do every condition at once), then run.
%
% OVERRIDES TABLE
%   Overrides = {
%     % Set  Amp  PTD(ms)  ResetAllFirst  AddBad              RemoveBad
%       1,   5,   0,       true,          [2,7,12:14,19],     [];
%       1,   5,   3,       false,         [5],                [23];
%   };
%
%   AddBad / RemoveBad use the RELATIVE trial ID within that condition
%   (ConditionTrialIndex - i.e. "trial 7 of this condition's 30 trials"),
%   matching what BadTrial_ReferenceCriteria_v1.m prints.
%
%   ResetAllFirst = true  : ignore whatever is already saved for this
%                           condition. The new final bad-trial list =
%                           AddBad.
%   ResetAllFirst = false : start from whatever is already saved for this
%                           condition (empty, with a warning, if nothing
%                           was saved yet), then add AddBad and remove
%                           RemoveBad - a patch, so you can fix one or two
%                           trials without retyping the whole list.
%   AddBad and RemoveBad must not share a trial on the same row - that is
%   treated as a contradiction and raises an error.
%
%   PTD is ignored (leave as 0) for a single-pulse dataset
%   (simultaneous_stim == 1).
%
% OUTPUT
%   <base_name>_BadTrials.mat, in the same folder as matFilePath
%   Main saved variable: BadConfirmed (struct array, one row per confirmed
%   condition). A matching condition from a previous run is REPLACED (not
%   duplicated); conditions not mentioned in Overrides this run are left
%   untouched. Also saves AllBadAbsoluteTrialIDs, the union of every
%   confirmed condition's bad AbsoluteTrialID - what you exclude later.
%   A timestamped backup of the previous file is kept.
% ========================================================================

clear;

%% ============================ USER SETTINGS ===========================

% Full path to the saved <base_name>_TrialSpikeCounts.mat file
matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_SimSeq1/Linearity_S1_100_400um_SimSeq1_TrialSpikeCounts.mat';

Condition_Tolerance = 1e-4;

%% --------------------- CONDITIONS BEING CONFIRMED NOW ------------------

Overrides = {
  % Set  Amp  PTD(ms)  ResetAllFirst  AddBad              RemoveBad
  %  100_400um SimSeq
    % 3uA 
    % Set 1->100, 2-> 300, 3->100, 4->400, 5->200, 6->300, 7->400, 8->200
    % 1,   3,   0, true, [9,13,5,7,9,15,29],  [];
    % 1,   3,   5, true, [2,4,6,10,12,14,19,20,24,28],  [];
    % 2,   3,   0, true, [10,12,13,23],  [];
    % 2,   3,   5, true, [4,7,8,10,14,16,25:30],  [];
    % % 3,   3,   0, true, [],  [];
    % 3,   3,   5, true, [3,8,21,19,25],  [];
    % 4,   3,   0, true, [3,5,6,7,12,17,21,30],  [];
    % 4,   3,   5, true, [9,10,11,15,18,23,27,30],  [];
    % 5,   3,   0, true, [5,6,13,14,21,22],  [];
    % 5,   3,   5, true, [7,8,15,18,21,22,24],  [];
    % % 6,   3,   0, true, [5,6,13,14,21,22],  [];
    % 6,   3,   5, true, [1,3,7,8,10,12,14,20,23,26,28,30],  [];
    % % 7,   3,   0, true, [5,6,13,14,21,22],  [];
    % 7,   3,   5, true, [1,4,15,16,18,20,23,25],  [];
    % % 8,   3,   0, true, [5,6,13,14,21,22],  [];
    % 8,   3,   5, true, [2:4,7,11,19,20,24],  [];

    % 5uA 
    % 1,   5,   0, true, [],  [];
    % 1,   5,   5, true, [], [];    

    % 10uA
    1,   10,  0, true, [],  []; % 100um 
    1,   10,  5, true, [2,8,28], [];
    2,   10,  0, true, [3,4,7,8],  []; % 300um 
    2,   10,  5, true, [3,4,5,13,20], [];
    3,   10,  0, true, [],  []; % 100um 
    3,   10,  5, true, [8,15,17,19,27], [];
    4,   10,  0, true, [4  11  22  26  30],  []; % 400um 
    4,   10,  5, true, [17,27], [];
    5,   10,  0, true, [1,5,6,8,11,16:18,24],  []; % 200um 
    5,   10,  5, true, [3  11 18  22  23  29], [];
    6,   10,  5, true, [10,20,25,28],  []; % 300um 
    7,   10,  5, true, [7,14,17,23,24],  []; % 400um 
    8,   10,  5, true, [9,10,18,19,21],  []; % 200um 
    
    
    
   
    % 100_400um Single
    % 1,  3, 0,   true, [3,4,5,7,10,12,13,19,20,22],[];
    % 1,  5, 0,   true, [9,10,30],[];
    % 1,  10,0,   true, [11,20:22,26],[];
    % 2,  3, 0,   true, [1,2,6:8,16,19,20,22,24,26,27,29],[];
    % 2,  5, 0,   true  [2,6,7,10,18,26:28],[];
    % 2,  10, 0,   true  [2,7,8,18,20,21:24,26],[];
    % 
    % 3,  3, 0,   true, [2,9,11,15,22,23,25,26,27],[];
    % 3,  5, 0,   true  [13,23,27],[];
    % 3,  10, 0,   true  [2,9,15,16,24],[];
    % 
    % 4,  3, 0,   true, [1,7,8,10,12:15,17,18,20,21,24:26,30],[];
    % 4,  5, 0,   true  [6,23,24],[];
    % 4,  10, 0,   true  [8,13,20,26],[];
    % 
    % 5,  3, 0,   true, [1,4,10,12,14,16,19,22,23],[];
    % 5,  5, 0,   true  [4,9],[];
    % 5,  10, 0,   true  [2,10,20],[];

    % 500_700 um SimSeq
    % 1, 3, 0, true, [7,27],[];
    % 1, 5, 0, true, [2,19,26],[];
    % 1, 10, 0, true, [15,28,30],[];
    % 1, 3, 5, true, [4],[];
    % 1, 5, 5, true, [10,28],[];
    % 1, 10, 5, true, [2,18],[];
    % 
    % 2, 3, 0, true, [8],[];
    % 2, 5, 0, true, [2   9  13  14],[];
    % 2, 10, 0, true, [2 3],[];
    % 2, 3, 5, true, [3   7  17  18],[];
    % 2, 5, 5, true, [8],[];
    % 2, 10, 5, true, [14],[];
    % 
    % 3, 3, 5, true, [2   4   7   9  10  11  13  29],[];
    % 3, 5, 5, true, [30],[];
    % 3, 10, 5, true, [2],[];
    % 
    % 4, 3, 0, true, [2   5   6  10  13  16  24  26],[];
    % 4, 5, 0, true, [3  12],[];
    % 4, 10, 0, true, [10  11  12],[];
    % 4, 3, 5, true, [5   8  13  16],[];
    % 4, 5, 5, true, [14],[];
    % 4, 10, 5, true, [],[];
    % 
    % 5, 3, 5, true, [2   3   4  10  11  15  17  18],[];
    % 5, 5, 5, true, [13],[];
    % 5, 10, 5, true, [17],[];
    % 
    % 6, 3, 5, true, [7   9  11  30],[];
    % 6, 5, 5, true, [],[];
    % 6, 10, 5, true, [13],[];

    % 500_700um Single 
    % 1, 3, 0, true, [1,7,10,30], [];
    % 1, 5, 0, true, [6,14,30], [];
    % 1, 10,0, true, [9], [];
    % 
    % 2, 3, 0, true, [9], [];
    % 2, 5, 0, true, [8,12,14,23], [];
    % 2, 10,0, true, [], [];
    % 
    % 3, 3, 0, true, [1,7,9,23,24], [];
    % 3, 5, 0, true, [7], [];
    % 3, 10,0, true, [2,6], [];
    % 
    % 4, 3, 0, true, [6,11,30], [];
    % 4, 5, 0, true, [20,23,30], [];
    % 4, 10,0, true, [2], [];
    

};

Create_Backup_If_Output_Exists = true;

%% =========================== INITIAL CHECKS ===========================

if ~isfile(matFilePath)
    error('File does not exist:\n%s',matFilePath);
end
if isempty(Overrides)
    error('Overrides is empty - nothing to confirm.');
end
if size(Overrides,2) ~= 6
    error('Overrides must have 6 columns: Set, Amp, PTD(ms), ResetAllFirst, AddBad, RemoveBad.');
end

fprintf('\n============================================================\n');
fprintf('BAD-TRIAL MANUAL CONFIRMATION\n');
fprintf('============================================================\n');
fprintf('Source file: %s\n',matFilePath);
fprintf('Rows in Overrides: %d\n',size(Overrides,1));

%% ===================== LOAD THE SPIKE-COUNT TABLE =====================

LoadData = load(matFilePath,'SpikeCounts','TrialSummary','sim_stim');

SpikeCounts  = LoadData.SpikeCounts;
TrialSummary = LoadData.TrialSummary;
sim_stim     = double(LoadData.sim_stim);

nSets = numel(SpikeCounts.set);

%% ================ LOAD EXISTING SAVED FILE (IF ANY) ====================

[dataset_folder,count_file_name] = fileparts(matFilePath);
base_name = regexprep(count_file_name,'_TrialSpikeCounts$','');
output_name = sprintf('%s_BadTrials.mat',base_name);
output_path = fullfile(dataset_folder,output_name);

if isfile(output_path)
    ExistingLoad = load(output_path,'BadConfirmed');
    if isfield(ExistingLoad,'BadConfirmed')
        BadConfirmed = ExistingLoad.BadConfirmed;
        fprintf('\nExisting saved file found (%d condition(s) already confirmed):\n%s\n', ...
            numel(BadConfirmed),output_path);
    else
        warning('Existing output file does not contain BadConfirmed - starting fresh.');
        BadConfirmed = empty_confirmed_struct();
    end
else
    fprintf('\nNo existing saved file - this run creates it.\n');
    BadConfirmed = empty_confirmed_struct();
end

%% ===================== VALIDATE + APPLY EACH OVERRIDE ROW ==============

nRows = size(Overrides,1);

for r = 1:nRows

    si          = Overrides{r,1};
    amp_value   = Overrides{r,2};
    ptd_value   = Overrides{r,3};
    resetFirst  = Overrides{r,4};
    addBad      = Overrides{r,5};
    removeBad   = Overrides{r,6};

    if isempty(si) || si < 1 || si > nSets || fix(si) ~= si
        error('Row %d: set %s is not a valid stimulation set (1-%d).',r,mat2str(si),nSets);
    end

    amp_values = [SpikeCounts.set(si).amp.amp_value];
    ai = find(abs(amp_values-amp_value) < Condition_Tolerance,1);
    if isempty(ai)
        error('Row %d: amp %g uA does not match any decoded amplitude for Set %d (%s).', ...
            r,amp_value,si,num2str(amp_values));
    end

    ptd_values_ms = [SpikeCounts.set(si).amp(ai).ptd.PTD_ms];
    pi = find(abs(ptd_values_ms-ptd_value) < Condition_Tolerance,1);
    if isempty(pi)
        error('Row %d: ptd %g ms does not match any decoded PTD for Set %d, Amp %g uA (%s).', ...
            r,ptd_value,si,amp_value,num2str(ptd_values_ms));
    end

    cond = SpikeCounts.set(si).amp(ai).ptd(pi);
    total_trials = cond.total_trials;
    current_ptd_ms = cond.PTD_ms;

    addBad    = validate_trial_list(addBad,total_trials,r,'AddBad');
    removeBad = validate_trial_list(removeBad,total_trials,r,'RemoveBad');

    overlap = intersect(addBad,removeBad);
    if ~isempty(overlap)
        error(['Row %d: trial(s) %s appear in BOTH AddBad and RemoveBad - ' ...
            'contradictory input.'],r,num2str(overlap));
    end

    %% ---------------- FIND EXISTING SAVED ENTRY (IF ANY) ---------------

    match_idx = [];
    for k = 1:numel(BadConfirmed)
        if BadConfirmed(k).set == si && ...
                abs(BadConfirmed(k).amp-amp_values(ai)) < Condition_Tolerance && ...
                abs(BadConfirmed(k).ptd_ms-current_ptd_ms) < Condition_Tolerance
            match_idx = k;
            break;
        end
    end

    if resetFirst || isempty(match_idx)
        if ~resetFirst && isempty(match_idx)
            warning(['Row %d: ResetAllFirst is false but no existing saved entry was found ' ...
                'for this condition - starting from an empty list.'],r);
        end
        baseline_trials = [];
    else
        baseline_trials = BadConfirmed(match_idx).condition_trial_indices;
    end

    final_trials = union(baseline_trials(:).',addBad);
    final_trials = setdiff(final_trials,removeBad);
    final_trials = sort(unique(final_trials),'ascend');

    %% ---------------- CONVERT TO ABSOLUTE TRIAL IDs ---------------------

    condition_mask = TrialSummary.SetIndex == si & ...
        TrialSummary.AmplitudeIndex == ai & ...
        TrialSummary.PTDIndex == pi;

    absolute_trial_ids = zeros(size(final_trials));
    for trial_position = 1:numel(final_trials)
        trial_mask = condition_mask & TrialSummary.ConditionTrialIndex == final_trials(trial_position);
        if sum(trial_mask) ~= 1
            error('Row %d: condition trial %d could not be uniquely matched in TrialSummary.', ...
                r,final_trials(trial_position));
        end
        absolute_trial_ids(trial_position) = TrialSummary.AbsoluteTrialID(trial_mask);
    end
    absolute_trial_ids = sort(unique(absolute_trial_ids));

    good_trial_indices = setdiff(1:total_trials,final_trials);

    %% ---------------- BUILD RECORD --------------------------------------

    stim_channels = SpikeCounts.set(si).stim_channels;

    if sim_stim == 1
        stim_type = 'single';
    elseif abs(current_ptd_ms) < Condition_Tolerance
        stim_type = 'simultaneous';
    else
        stim_type = 'sequential';
    end

    R = empty_confirmed_struct();
    R(1).set = si;
    R(1).stimChannels = stim_channels;
    R(1).stimOrderLabel = SpikeCounts.set(si).set_name;
    R(1).stimType = stim_type;
    R(1).amp = amp_values(ai);
    R(1).ptd_ms = current_ptd_ms;
    R(1).distance_um = SpikeCounts.set(si).distance_um;
    R(1).distance_note = SpikeCounts.set(si).distance_note;
    R(1).n_total_trials = total_trials;
    R(1).condition_trial_indices = final_trials;
    R(1).absolute_trial_ids = absolute_trial_ids;
    R(1).good_condition_trial_indices = good_trial_indices;
    R(1).n_bad = numel(final_trials);
    R(1).confirmed_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

    if isempty(match_idx)
        BadConfirmed(end+1) = R; %#ok<SAGROW>
    else
        BadConfirmed(match_idx) = R;
    end

    fprintf('\nRow %d: Set %d (%s) | Amp %g uA | PTD %g ms | %s | ResetAllFirst=%d\n', ...
        r,si,R(1).stimOrderLabel,R(1).amp,R(1).ptd_ms,stim_type,resetFirst);
    fprintf('  Bad trials (%d/%d), relative: %s\n', ...
        R(1).n_bad,total_trials,number_list(final_trials));
end

%% =========================== BUILD ALL-BAD LIST ========================

AllBadAbsoluteTrialIDs = [];
for k = 1:numel(BadConfirmed)
    AllBadAbsoluteTrialIDs = union(AllBadAbsoluteTrialIDs,BadConfirmed(k).absolute_trial_ids(:).');
end
AllBadAbsoluteTrialIDs = sort(unique(AllBadAbsoluteTrialIDs(:).'));

%% =========================== SAVE RESULT ==============================

if isfile(output_path) && Create_Backup_If_Output_Exists
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    backup_name = sprintf('%s_BadTrials_BACKUP_%s.mat',base_name,timestamp);
    backup_path = fullfile(dataset_folder,backup_name);
    copyfile(output_path,backup_path);
    fprintf('\nExisting output backed up to:\n%s\n',backup_path);
end

SourceTrialCountFile = matFilePath;

save(output_path,'BadConfirmed','AllBadAbsoluteTrialIDs','SourceTrialCountFile','-v7.3');

fprintf('\n============================================================\n');
fprintf('BAD-TRIAL MANUAL CONFIRMATION COMPLETE\n');
fprintf('Saved to: %s\n',output_path);
fprintf('Total conditions confirmed (all-time, this file): %d\n',numel(BadConfirmed));
fprintf('Total bad trials across all confirmed conditions: %d\n',numel(AllBadAbsoluteTrialIDs));
fprintf('============================================================\n');

fprintf('\nAll confirmed conditions in this file:\n');
for k = 1:numel(BadConfirmed)
    R = BadConfirmed(k);
    fprintf('  Set %d (%s) | Amp %g uA | PTD %g ms | %s | Bad %d/%d\n', ...
        R.set,R.stimOrderLabel,R.amp,R.ptd_ms,R.stimType,R.n_bad,R.n_total_trials);
end

%% =========================== LOCAL FUNCTIONS ==========================

function trials = validate_trial_list(trials,total_trials,row_index,field_name)
if isempty(trials)
    trials = [];
    return;
end
trials = unique(double(trials(:).'));
invalid = ~isfinite(trials) | trials < 1 | trials > total_trials | fix(trials) ~= trials;
if any(invalid)
    error('Row %d: %s contains invalid relative trial numbers (must be 1-%d): %s', ...
        row_index,field_name,total_trials,num2str(trials(invalid)));
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