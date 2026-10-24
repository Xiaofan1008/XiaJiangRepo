%% ========================================================================
% BAD-TRIAL ENTRY FOR THE OTHER (ChB) SINGLE ELECTRODES - v1
%
% PURPOSE
%   Record the final bad trials of every single-pulse electrode EXCEPT the
%   anchor (Ch1, done with BadTrial_Ch1_Single_v1) for all single datasets of
%   one animal, using the relative trial numbers shown on the figures of
%   Clean_ChB_Single_v1.m. One row per electrode, one column per amplitude.
%
% HOW TO FILL distance_bad_trials
%   Amps = [3 5 10];                       % the amplitude columns, in this order
%   distance_bad_trials(1).matFilePath = '.../<base_name>_TrialSpikeCounts.mat';
%   distance_bad_trials(1).bad = {
%     % StimCh   3 uA        5 uA        10 uA
%        3,      [2 9],      [],         [4];
%        5,      [],         [7],        [];
%   };
%   - a row = one electrode; its columns follow Amps
%   - [] = you looked at that condition and it has NO bad trials
%   - a (electrode, amplitude) you do not list is left UNTOUCHED (not reviewed)
%   - a listed condition REPLACES whatever was saved for it (final list, not a patch)
%   - only the electrodes / amplitudes you list are changed; everything else in the
%     dataset's _BadTrials.mat (including Ch1) stays as it is
%
% OUTPUT
%   <base_name>_BadTrials.mat in each dataset folder, same format as
%   BadTrial_ManualConfirm_v1 (BadConfirmed, AllBadAbsoluteTrialIDs,
%   SourceTrialCountFile). A timestamped backup of the previous file is kept.
%   After saving, per dataset it prints
%     - a summary of the bad trials per electrode and amplitude
%     - a REVIEW CHECKLIST: single-pulse conditions with no saved entry at all,
%       so a forgotten condition cannot be mistaken for "no bad trials"
% ========================================================================

clear; clc;

%% ============================ USER SETTINGS ===========================

Anchor_Channel = 1;              % handled by BadTrial_Ch1_Single_v1 (not allowed here)
Condition_Tolerance = 1e-4;
Create_Backup_If_Output_Exists = false;

Amps = [3 5 10];                 % amplitude columns of every bad-trial table below (uA)

distance_bad_trials = struct('matFilePath',{},'bad',{});

distance_bad_trials(1).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1/Linearity_S1_100_400um_Single1_TrialSpikeCounts.mat';
distance_bad_trials(1).bad = {
  % StimCh   3 uA     5 uA     10 uA
     3,      [3,4,7,12,13,16,19,20,21],      [1,2,3,7,8,9,10,26,28],      [];
     5,      [5,11,13,14,25:27,29],      [4,8,13,14,21:23,27],      [2,8,9,15,16,24];
     7,      [3,4,10,22],      [3,4,9,13,15,19,22,27],      [1,2,3,10,20,30];   
     9,      [1,2,5,6,7,10:13,14,16,20,21,27,28], [6,13,16,23,24,29,30], [1,8,11,15,28];

};

distance_bad_trials(2).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1/Linearity_S1_500_700um_Single1_TrialSpikeCounts.mat';
distance_bad_trials(2).bad = {
  % StimCh   3 uA     5 uA     10 uA
     11,     [1,10,8,10,12,13,14,16,26], [1,3,8,9,17,19,22,30], [2,3,5,8:10,12,14,16,18,24,25,28,29];
     13,     [1,11,27,28], [2,6,13:15,26], [7,10,14,15,20,24,21,2,4];
     15,     [1,10,13,14,16,20,21,25], [1,7,10,13,16,20,26], [3,4,7,13,19,22,23,29];
};

%% =========================== APPLY EACH DATASET ========================

if isempty(distance_bad_trials), error('distance_bad_trials is empty.'); end
Amps = double(Amps(:).');
if isempty(Amps), error('Amps is empty.'); end

fprintf('\n============================================================\n');
fprintf('ChB SINGLE BAD TRIALS (v1) - %d dataset(s) | amplitude columns: %s uA\n', ...
    numel(distance_bad_trials),num2str(Amps));
fprintf('============================================================\n');

for d = 1:numel(distance_bad_trials)

    matFilePath = distance_bad_trials(d).matFilePath;
    bad_rows    = distance_bad_trials(d).bad;

    fprintf('\n--- Dataset %d ---\n%s\n',d,matFilePath);

    if ~isfile(matFilePath)
        error('Dataset %d: file does not exist:\n%s',d,matFilePath);
    end
    if isempty(bad_rows) || size(bad_rows,2) ~= 1+numel(Amps)
        error('Dataset %d: bad needs %d columns (StimCh + one per amplitude in Amps = [%s]) and at least one row.', ...
            d,1+numel(Amps),num2str(Amps));
    end

    LoadData = load(matFilePath,'SpikeCounts','TrialSummary','sim_stim');
    SpikeCounts  = LoadData.SpikeCounts;
    TrialSummary = LoadData.TrialSummary;
    sim_stim     = double(LoadData.sim_stim);

    if sim_stim ~= 1
        error('Dataset %d is not a single-pulse dataset (simultaneous_stim = %d).',d,sim_stim);
    end

    % ---- single-pulse set of every electrode ----
    nSets = numel(SpikeCounts.set);
    set_stim = nan(nSets,1);
    for ss = 1:nSets
        sc = SpikeCounts.set(ss).stim_channels;
        sc = sc(sc > 0);
        if numel(sc) == 1, set_stim(ss) = sc(1); end
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

    % ---- apply each electrode row, each amplitude column ----
    listed_stims = [];
    for r = 1:size(bad_rows,1)

        stim_ch = bad_rows{r,1};
        if ~isnumeric(stim_ch) || ~isscalar(stim_ch)
            error('Dataset %d, row %d: first column must be a single stim channel number.',d,r);
        end
        if stim_ch == Anchor_Channel
            error('Dataset %d, row %d: Ch%d is the anchor - enter it with BadTrial_Ch1_Single_v1.',d,r,Anchor_Channel);
        end
        if any(listed_stims == stim_ch)
            error('Dataset %d: Ch%d appears in more than one row.',d,stim_ch);
        end
        listed_stims(end+1) = stim_ch; %#ok<SAGROW>

        si = find(set_stim == stim_ch,1);
        if isempty(si)
            error('Dataset %d, row %d: no single-pulse set for Ch%d (available: %s).', ...
                d,r,stim_ch,num2str(set_stim(~isnan(set_stim)).'));
        end

        amp_values = [SpikeCounts.set(si).amp.amp_value];

        for k = 1:numel(Amps)

            trials_in = bad_rows{r,1+k};
            if ~isnumeric(trials_in)
                error('Dataset %d, Ch%d, %g uA: the trial list must be numeric, e.g. [2 9] or [].',d,stim_ch,Amps(k));
            end

            ai = find(abs(amp_values-Amps(k)) < Condition_Tolerance,1);
            if isempty(ai)
                error('Dataset %d, Ch%d: amp %g uA not found (available: %s).',d,stim_ch,Amps(k),num2str(amp_values));
            end
            pi = 1;     % single-pulse dataset: one PTD
            cond = SpikeCounts.set(si).amp(ai).ptd(pi);
            total_trials = cond.total_trials;

            final_trials = validate_trial_list(trials_in,total_trials,d,stim_ch,Amps(k));

            condition_mask = TrialSummary.SetIndex == si & ...
                TrialSummary.AmplitudeIndex == ai & ...
                TrialSummary.PTDIndex == pi;

            absolute_trial_ids = zeros(size(final_trials));
            for t = 1:numel(final_trials)
                trial_mask = condition_mask & TrialSummary.ConditionTrialIndex == final_trials(t);
                if sum(trial_mask) ~= 1
                    error('Dataset %d, Ch%d, %g uA: condition trial %d could not be uniquely matched.', ...
                        d,stim_ch,Amps(k),final_trials(t));
                end
                absolute_trial_ids(t) = TrialSummary.AbsoluteTrialID(trial_mask);
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
            for m = 1:numel(BadConfirmed)
                if BadConfirmed(m).set == si && ...
                        abs(BadConfirmed(m).amp-amp_values(ai)) < Condition_Tolerance && ...
                        abs(BadConfirmed(m).ptd_ms-cond.PTD_ms) < Condition_Tolerance
                    match_idx = m;
                    break;
                end
            end
            if isempty(match_idx)
                BadConfirmed(end+1) = R; %#ok<SAGROW>
            else
                BadConfirmed(match_idx) = R;
            end
        end
    end

    % ---- all-bad list + save ----
    AllBadAbsoluteTrialIDs = [];
    for m = 1:numel(BadConfirmed)
        AllBadAbsoluteTrialIDs = union(AllBadAbsoluteTrialIDs,BadConfirmed(m).absolute_trial_ids(:).');
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

    % ---- summary: bad trials per electrode x amplitude (all single-pulse sets) ----
    fprintf('\nSummary - bad trials per electrode and amplitude (bad/total; "-" = no saved entry):\n');
    fprintf('%-6s %-9s', 'Stim','Dist(um)');
    for k = 1:numel(Amps), fprintf(' %12s',sprintf('%g uA',Amps(k))); end
    fprintf('\n');
    missing = {};
    for ss = 1:nSets
        if isnan(set_stim(ss)), continue; end
        d_um = SpikeCounts.set(ss).distance_um;
        if isnan(d_um), dtxt = 'n/a'; else, dtxt = sprintf('%g',d_um); end
        fprintf('Ch%-4d %-9s',set_stim(ss),dtxt);
        for k = 1:numel(Amps)
            hit = find([BadConfirmed.set] == ss & abs([BadConfirmed.amp]-Amps(k)) < Condition_Tolerance,1);
            if isempty(hit)
                fprintf(' %12s','-');
                missing{end+1} = sprintf('Ch%d @ %g uA',set_stim(ss),Amps(k)); %#ok<SAGROW>
            else
                fprintf(' %12s',sprintf('%d/%d',BadConfirmed(hit).n_bad,BadConfirmed(hit).n_total_trials));
            end
        end
        fprintf('\n');
    end

    fprintf('\nREVIEW CHECKLIST (single-pulse conditions with NO saved entry):\n');
    if isempty(missing)
        fprintf('  none - every electrode and amplitude has been confirmed.\n');
    else
        fprintf('  %s\n',strjoin(missing,' | '));
        fprintf('  (these are NOT "no bad trials" - they have simply not been reviewed yet)\n');
    end
    fprintf('Total bad trials in this file (all sets): %d\n',numel(AllBadAbsoluteTrialIDs));
end

fprintf('\nDone. Re-run Clean_ChB_Single_v1 to check what is left.\n');

%% =========================== LOCAL FUNCTIONS ==========================

function trials = validate_trial_list(trials,total_trials,dataset_index,stim_ch,amp_value)
if isempty(trials)
    trials = [];
    return;
end
trials = unique(double(trials(:).'));
invalid = ~isfinite(trials) | trials < 1 | trials > total_trials | fix(trials) ~= trials;
if any(invalid)
    error('Dataset %d, Ch%d, %g uA: invalid relative trial numbers (must be 1-%d): %s', ...
        dataset_index,stim_ch,amp_value,total_trials,num2str(trials(invalid)));
end
end

function R = empty_confirmed_struct()
R = struct('set',{},'stimChannels',{},'stimOrderLabel',{},'stimType',{}, ...
    'amp',{},'ptd_ms',{},'distance_um',{},'distance_note',{},'n_total_trials',{}, ...
    'condition_trial_indices',{},'absolute_trial_ids',{},'good_condition_trial_indices',{}, ...
    'n_bad',{},'confirmed_at',{});
end