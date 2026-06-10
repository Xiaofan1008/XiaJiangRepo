%% ========================================================================
%  MANUAL BAD TRIAL ENTRY for Mixed-Prefix Multi-Electrode Stimulation
%  Purpose:
%    Manually flag bad trials globally for ALL recording channels.
%  Manual entry format:
%    Bad_Groups(c).Set = setID;
%    Bad_Groups(c).Trials = {
%        % Type,   Amp,  PrefixOrActiveCount,  ISI_ms,  Bad relative trial IDs
%          'seq',  5,    5,                    5,       [1 4 8];
%          'seq',  10,   5,                    5,       [];
%          'sim',  5,    5,                    NaN,     [3 9];
%    };
%  Meaning:
%    For 'seq':
%       PrefixOrActiveCount = PrefixLength
%       ISI_ms              = ISI
%    For 'sim':
%       PrefixOrActiveCount = ActiveElectrodeCount
%       ISI_ms              = ignored, use NaN
%  Output:
%       BadTrials{ich}
%       GoodTrials{ich}
%       BadTrialTable
%       BadTrialInfo
% ========================================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis'));

%% ====================== USER SETTINGS ========================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

Electrode_Type = 2;   % 0 rigid, 1 single-shank flex, 2 four-shank flex

% If true:
%   load existing bad-trial file and add new entries.
%
% If false:
%   start from empty BadTrials and overwrite the file.
append_to_existing = true;

% Tolerances for matching amplitude and ISI.
Amp_Tol_uA = 0.01;
ISI_Tol_ms = 0.01;

%% ========================================================================
%  MANUAL BAD TRIALS
%
%  Each row:
%       Type, Amp, PrefixOrActiveCount, ISI_ms, Bad relative trial IDs
%  Type:
%       'seq' = sequential prefix/recovery condition
%       'sim' = simultaneous condition
%  For 'seq':
%       PrefixOrActiveCount = PrefixLength
%       ISI_ms = ISI
%  For 'sim':
%       PrefixOrActiveCount = ActiveElectrodeCount
%       ISI_ms = NaN
%  Example:
%       'seq', 10, 5, 5, [1 3 7]
%       means Set X, 10 uA, Prefix 5, ISI 5 ms,
%       relative trials 1, 3, and 7 are bad.
%
%       'sim', 10, 5, NaN, [2 8]
%       means Set X, 10 uA, simultaneous 5-electrode condition,
%       relative trials 2 and 8 are bad.
% ========================================================================

Bad_Groups = [];

% ================= Set 1 =================
c = 1;
Bad_Groups(c).Set = 1;
Bad_Groups(c).Trials = {
    % Type,   Amp,  Prefix/ActiveCount,  ISI_ms,  Bad relative trial IDs
      'seq',  5,    5,                    5,       [5:9,14:15,17:19,21:23,25:28,30];
      'seq',  10,   5,                    5,       [6,7,14,17,18,20,21,26,28:29];
      'sim',  5,    5,                    NaN,     [5,10,11,18,20,22:25,29:30];
      'sim',  10,   5,                    NaN,     [2,15,19,20:24,29];
};

% ================= Set 2 =================
c = 2;
Bad_Groups(c).Set = 2;
Bad_Groups(c).Trials = {
    % Type,   Amp,  Prefix/ActiveCount,  ISI_ms,  Bad relative trial IDs
      'seq',  5,    5,                    5,       [1,2,9:10,13,15:16,24,29];
      'seq',  10,   5,                    5,       [4:5,10,14,16,23,26:27,30];
};

% ================= Set 3 =================
c = 3;
Bad_Groups(c).Set = 3;
Bad_Groups(c).Trials = {
    % Type,   Amp,  Prefix/ActiveCount,  ISI_ms,  Bad relative trial IDs
      'seq',  5,    5,                    5,       [2:3,5,9,16:21,23,29];
      'seq',  10,   5,                    5,       [1:5,8,18:19,27,28,30];
};

%% ====================== INITIAL SETUP ========================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end

cd(data_folder);
fprintf('\nRunning manual bad-trial entry in:\n%s\n\n', data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');

if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ====================== ELECTRODE MAP ========================

d = Depth_s(Electrode_Type);
nCh_Total = length(d);

%% ====================== LOAD EXPERIMENT FILE =================

f_exp = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));

if isempty(f_exp)
    error('No *_exp_datafile_*.mat found.');
end

S_exp = load(fullfile(data_folder, f_exp(1).name));

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

%% ====================== REMOVE HEADER ROW ====================

StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(StimParams_data,1), expected_rows);
end

%% ====================== TRIAL METADATA =======================

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

%% ====================== AMPLITUDE PER TRIAL ==================

trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

%% ====================== LOAD OR CREATE BAD TRIALS FILE ========

bad_filename = sprintf('%s.MixedPrefixBadTrials.mat', base_name);
save_file_path = fullfile(data_folder, bad_filename);

if append_to_existing && isfile(save_file_path)

    fprintf('\nFound existing file:\n%s\nLoading previous bad trials...\n', save_file_path);

    tmp = load(save_file_path);

    if isfield(tmp, 'BadTrials')
        BadTrials = tmp.BadTrials;
    else
        warning('Existing file found but BadTrials missing. Creating empty BadTrials.');
        BadTrials = cell(nCh_Total, 1);
    end

    if isfield(tmp, 'BadTrialTable')
        ExistingBadTrialTable = tmp.BadTrialTable;
    else
        ExistingBadTrialTable = table();
    end

else

    if append_to_existing
        fprintf('\nNo existing bad-trial file found. Creating new file:\n%s\n', save_file_path);
    else
        fprintf('\nappend_to_existing = false. Starting fresh and overwriting:\n%s\n', save_file_path);
    end

    BadTrials = cell(nCh_Total, 1);
    ExistingBadTrialTable = table();
end

if numel(BadTrials) < nCh_Total
    BadTrials{nCh_Total} = [];
end

%% ====================== PROCESS MANUAL ENTRIES ===============

NewRows = {};
rowCounter = 0;
count_added = 0;
count_skipped_empty = 0;
count_skipped_notfound = 0;
count_skipped_outofrange = 0;
count_duplicate = 0;

if isempty(Bad_Groups)

    fprintf('\nNo Bad_Groups entries. Nothing to add.\n');

else

    fprintf('\nProcessing manual bad-trial entries...\n');

    for c = 1:numel(Bad_Groups)

        tgt_set = Bad_Groups(c).Set;

        if ~isfield(Bad_Groups(c), 'Trials') || isempty(Bad_Groups(c).Trials)
            fprintf('  [SKIP] Set %d has no Trials table.\n', tgt_set);
            continue;
        end

        TrialsTable = Bad_Groups(c).Trials;

        for r = 1:size(TrialsTable, 1)

            tgt_type = lower(string(TrialsTable{r,1}));
            tgt_amp  = TrialsTable{r,2};
            tgt_count_or_prefix = TrialsTable{r,3};
            tgt_isi  = TrialsTable{r,4};
            bad_relative_list = TrialsTable{r,5};

            if isempty(bad_relative_list)
                count_skipped_empty = count_skipped_empty + 1;
                continue;
            end

            %% ----- Find matching absolute trial list -----
            if tgt_type == "seq" || tgt_type == "prefix"

                condition_label = 'seq';
                tgt_prefix = tgt_count_or_prefix;
                tgt_active_count = NaN;

                matching_trials = find(conditionType_trial == 1 & ...
                                       conditionSetID_trial == tgt_set & ...
                                       abs(trialAmps - tgt_amp) < Amp_Tol_uA & ...
                                       prefixLength_trial == tgt_prefix & ...
                                       abs(isi_ms_trial - tgt_isi) < ISI_Tol_ms);

            elseif tgt_type == "sim" || tgt_type == "simultaneous"

                condition_label = 'sim';
                tgt_prefix = NaN;
                tgt_active_count = tgt_count_or_prefix;

                matching_trials = find(conditionType_trial == 2 & ...
                                       conditionSetID_trial == tgt_set & ...
                                       abs(trialAmps - tgt_amp) < Amp_Tol_uA & ...
                                       activeCount_trial == tgt_active_count);

            else

                fprintf('  [SKIP] Unknown type "%s" in Set %d row %d. Use ''seq'' or ''sim''.\n', ...
                    tgt_type, tgt_set, r);
                continue;
            end

            if isempty(matching_trials)

                fprintf('  [SKIP] No matching trials: Set %d | Type %s | Amp %.1f | Count/Prefix %.1f | ISI %.1f\n', ...
                    tgt_set, condition_label, tgt_amp, tgt_count_or_prefix, tgt_isi);

                count_skipped_notfound = count_skipped_notfound + 1;
                continue;
            end

            %% ----- Convert relative trial IDs to absolute trial IDs -----
            for b = 1:numel(bad_relative_list)

                rel_trial = bad_relative_list(b);

                if rel_trial < 1 || rel_trial > numel(matching_trials)

                    fprintf('  [SKIP] Set %d | Type %s | Amp %.1f | rel trial %d out of range. Only %d trials exist.\n', ...
                        tgt_set, condition_label, tgt_amp, rel_trial, numel(matching_trials));

                    count_skipped_outofrange = count_skipped_outofrange + 1;
                    continue;
                end

                abs_trial = matching_trials(rel_trial);

                %% ----- Apply globally to all channels -----
                is_new_global = false;

                for ich = 1:nCh_Total

                    current_list = BadTrials{ich};

                    if isempty(current_list)
                        current_list = [];
                    end

                    new_list = unique([current_list(:); abs_trial(:)]);

                    if numel(new_list) > numel(current_list)
                        is_new_global = true;
                    end

                    BadTrials{ich} = new_list;
                end

                if is_new_global

                    count_added = count_added + 1;

                    fprintf('  -> Added GLOBAL bad trial: Set %d | %s | Amp %.1f | Prefix %.1f | ActiveCount %.1f | ISI %.1f | Rel %d -> Abs %d\n', ...
                        tgt_set, condition_label, tgt_amp, tgt_prefix, tgt_active_count, tgt_isi, rel_trial, abs_trial);

                    rowCounter = rowCounter + 1;

                    NewRows(rowCounter,:) = { ...
                        tgt_set, ...
                        condition_label, ...
                        tgt_amp, ...
                        tgt_prefix, ...
                        tgt_active_count, ...
                        tgt_isi, ...
                        rel_trial, ...
                        abs_trial};

                else

                    count_duplicate = count_duplicate + 1;

                    fprintf('  [DUPLICATE] Already marked: Set %d | %s | Amp %.1f | Rel %d -> Abs %d\n', ...
                        tgt_set, condition_label, tgt_amp, rel_trial, abs_trial);
                end
            end
        end
    end
end

fprintf('\nManual entry processing finished.\n');
fprintf('  Added new global bad trials : %d\n', count_added);
fprintf('  Duplicate entries skipped   : %d\n', count_duplicate);
fprintf('  Empty rows skipped          : %d\n', count_skipped_empty);
fprintf('  Not-found rows skipped      : %d\n', count_skipped_notfound);
fprintf('  Out-of-range trials skipped : %d\n', count_skipped_outofrange);

%% ====================== CREATE / UPDATE BAD TRIAL TABLE =======

if isempty(NewRows)

    NewBadTrialTable = table();

else

    NewBadTrialTable = cell2table(NewRows, ...
        'VariableNames', { ...
        'SetID', ...
        'ConditionType', ...
        'Amp_uA', ...
        'PrefixLength', ...
        'ActiveCount', ...
        'ISI_ms', ...
        'RelativeTrialID', ...
        'AbsoluteTrialID'});

end

if append_to_existing && ~isempty(ExistingBadTrialTable)

    if ~isempty(NewBadTrialTable)
        BadTrialTable = [ExistingBadTrialTable; NewBadTrialTable];
        BadTrialTable = unique(BadTrialTable, 'rows', 'stable');
    else
        BadTrialTable = ExistingBadTrialTable;
    end

else

    BadTrialTable = NewBadTrialTable;
end

%% ====================== GOOD TRIALS ===========================

GoodTrials = cell(nCh_Total, 1);

allTrialIDs = (1:n_Trials)';

for ich = 1:nCh_Total

    if isempty(BadTrials{ich})
        BadTrials{ich} = [];
    end

    BadTrials{ich} = unique(BadTrials{ich}(:));
    GoodTrials{ich} = setdiff(allTrialIDs, BadTrials{ich});
end

%% ====================== SAVE ================================

BadTrialInfo = struct();
BadTrialInfo.data_folder = data_folder;
BadTrialInfo.base_name = base_name;
BadTrialInfo.Electrode_Type = Electrode_Type;
BadTrialInfo.append_to_existing = append_to_existing;
BadTrialInfo.Amp_Tol_uA = Amp_Tol_uA;
BadTrialInfo.ISI_Tol_ms = ISI_Tol_ms;
BadTrialInfo.n_Trials = n_Trials;
BadTrialInfo.nCh_Total = nCh_Total;
BadTrialInfo.manual_entry_format = ...
    'Bad_Groups(c).Trials = {Type, Amp, PrefixOrActiveCount, ISI_ms, BadRelativeTrials}';
BadTrialInfo.created_or_updated_on = datestr(now);

fprintf('\nSaving bad trials to:\n%s\n', save_file_path);

save(save_file_path, ...
    'BadTrials', ...
    'GoodTrials', ...
    'BadTrialTable', ...
    'BadTrialInfo', ...
    'Bad_Groups', ...
    'activeCount_trial', ...
    'prefixLength_trial', ...
    'isi_ms_trial', ...
    'conditionType_trial', ...
    'conditionSetID_trial', ...
    'trialAmps', ...
    '-v7.3');

fprintf('Save complete.\n');

%% ====================== QUICK CHECK ===========================

% total_bad_unique = numel(unique(BadTrials{1}));
% 
% fprintf('\nQuick check:\n');
% fprintf('  Number of globally bad absolute trials: %d\n', total_bad_unique);
% 
% if ~isempty(BadTrialTable)
%     disp(BadTrialTable);
% else
%     fprintf('  BadTrialTable is empty.\n');
% end