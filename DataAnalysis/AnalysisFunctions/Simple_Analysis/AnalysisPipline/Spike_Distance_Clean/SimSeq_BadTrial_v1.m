%% ========================================================================
% BAD-TRIAL ENTRY FOR THE SIM / SEQ (PAIRED) DATASETS - v1
%
% PURPOSE
%   Record the final bad trials of the paired datasets (all distances of one
%   dataset in one place) from the template printed by Clean_Paired_v1.m
%   (section 5). The template is already filled in with the NEW flagged trials
%   as suggestions: paste it below and edit only what is inside the brackets.
%
% HOW A ROW READS
%   Amps = [3 5 10];                          % the amplitude columns, in this order
%   distance_bad_trials(1).matFilePath = '.../<base_name>_TrialSpikeCounts.mat';
%   distance_bad_trials(1).bad = {
%     % Dist  Order  PTD    3 uA       5 uA     10 uA
%       100,   1,     0,    [2 9],     [],      [4];      % Ch1 -> Ch3 (circle)  Sim
%       100,   1,     5,    [],        [7],     [];       % Ch1 -> Ch3 (circle)  Seq
%       100,   2,     5,    [1 5],     [3],     [];       % Ch3 -> Ch1 (diamond) Seq
%   };
%   Dist  = distance in um (as in the figure title)
%   Order = 1 circle (the lower channel number is stimulated first, e.g. Ch1 -> ChB)
%           2 diamond (the higher channel number first, e.g. ChB -> Ch1)
%   PTD   = 0 for simultaneous, otherwise the PTD in ms (e.g. 5)
%   A PTD may exist for only one stim order (e.g. Sim); list only the rows the
%   template shows - a row that does not exist in the dataset stops with an error.
%
%   [] = you looked at that condition and it has NO bad trials
%   a condition you do not list is left UNTOUCHED (not reviewed)
%   a listed condition REPLACES whatever was saved for it (final list, not a patch)
%
% OUTPUT
%   <base_name>_BadTrials.mat in each dataset folder, same format as
%   BadTrial_ManualConfirm_v1 (BadConfirmed, AllBadAbsoluteTrialIDs,
%   SourceTrialCountFile), so SpikeTRatio and Clean_Paired_v1 read it as before.
%   A timestamped backup of the previous file is kept. The bad trials of a
%   paired file only change the ACTUAL; the prediction is unaffected.
%   After saving, per dataset it prints
%     - the bad / total trials of every condition, with the exclusion rate
%     - the exclusion rate per distance and PTD (does one distance lose far
%       more than the others?) and a warning above Warn_Exclusion_Fraction
%     - a REVIEW CHECKLIST of conditions with no saved entry
% ========================================================================

clear; clc;

%% ============================ USER SETTINGS ===========================

Condition_Tolerance = 1e-4;
Distance_Tol_um = 1;
Create_Backup_If_Output_Exists = true;
Warn_Exclusion_Fraction = 0.30;   % warn when a condition loses more than this share of its trials

Amps = [3 5 10];                  % amplitude columns of every bad-trial table below (uA)

distance_bad_trials = struct('matFilePath',{},'bad',{});

distance_bad_trials(1).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_SimSeq1/Linearity_S1_100_400um_SimSeq1_TrialSpikeCounts.mat';
distance_bad_trials(1).bad = {
  % Dist  Order  PTD   3 uA    5 uA    10 uA
    100,   1,   0,   [10,22,23,26,24],   [1 4 7 9 14 19 20 25 29 27],   [3 5 6 7 10 12 17:19 23 27];   % Ch1 -> Ch3 (circle)  Sim
    100,   1,   5,   [1 2 9 15:17 21:23 26:30],   [1:3 6 10:13 17 22 23 25 30],   [7 10:11 13 16 17 20:26];   % Ch1 -> Ch3 (circle)  Seq 5 ms
    100,   2,   5,   [2 4 7 9 10 12 14 15:17 25 28 30],   [1:2 5 6 10 12 18 21 26 28],   [1:7 11 18 22:24 27:30];   % Ch3 -> Ch1 (diamond)  Seq 5 ms

    200,   1,   0,   [11 28 30],   [10 25],   [3 9 20];   % Ch1 -> Ch5 (circle)  Sim
    200,   1,   5,   [6 17],   [8 17],   [2 5 20];   % Ch1 -> Ch5 (circle)  Seq 5 ms
    200,   2,   5,   [8 13 27],   [11 16 20],   [3 17];   % Ch5 -> Ch1 (diamond)  Seq 5 ms

    300,   1,   0,   [1 2 4 5 6 7],   [5 16 23 29 30],   [2 11 23];   % Ch1 -> Ch7 (circle)  Sim
    300,   1,   5,   [5 9 15 17 22 24],   [1],   [6 8 9 12 24 26];   % Ch1 -> Ch7 (circle)  Seq 5 ms
    300,   2,   5,   [4 18 19 22],   [2 6 11 14 19 27],   [12 16 22];   % Ch7 -> Ch1 (diamond)  Seq 5 ms

    400,   1,   0,   [8 10 16 19 23 24 28],   [21 24 25],   [5 8 10 17 23 24 25 27 28];   % Ch1 -> Ch9 (circle)  Sim
    400,   1,   5,   [4 20],   [7 10 21 28],   [10];   % Ch1 -> Ch9 (circle)  Seq 5 ms
    400,   2,   5,   [8 24],   [5],   [12 13];   % Ch9 -> Ch1 (diamond)  Seq 5 ms
};

distance_bad_trials(2).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_SimSeq1/Linearity_S1_500_700um_SimSeq1_TrialSpikeCounts.mat';
distance_bad_trials(2).bad = {
  % Dist  Order  PTD    3 uA     5 uA     10 uA
    500,   1,     0,    [],      [],      [];
    500,   1,     5,    [],      [],      [];
    500,   2,     5,    [],      [],      [];
};

%% =========================== APPLY EACH DATASET ========================

if isempty(distance_bad_trials), error('distance_bad_trials is empty.'); end
Amps = double(Amps(:).');
if isempty(Amps), error('Amps is empty.'); end

fprintf('\n============================================================\n');
fprintf('PAIRED BAD TRIALS (v1) - %d dataset(s) | amplitude columns: %s uA\n', ...
    numel(distance_bad_trials),num2str(Amps));
fprintf('============================================================\n');

for d = 1:numel(distance_bad_trials)

    matFilePath = distance_bad_trials(d).matFilePath;
    bad_rows    = distance_bad_trials(d).bad;

    fprintf('\n--- Dataset %d ---\n%s\n',d,matFilePath);

    if ~isfile(matFilePath)
        error('Dataset %d: file does not exist:\n%s',d,matFilePath);
    end
    if isempty(bad_rows) || size(bad_rows,2) ~= 3+numel(Amps)
        error('Dataset %d: bad needs %d columns (Dist, Order, PTD + one per amplitude in Amps = [%s]) and at least one row.', ...
            d,3+numel(Amps),num2str(Amps));
    end

    LoadData = load(matFilePath,'SpikeCounts','TrialSummary','sim_stim');
    SpikeCounts  = LoadData.SpikeCounts;
    TrialSummary = LoadData.TrialSummary;
    sim_stim     = double(LoadData.sim_stim);

    if sim_stim < 2
        error('Dataset %d is not a paired dataset (simultaneous_stim = %d).',d,sim_stim);
    end

    % ---- describe every two-channel set: distance and order (1 = lower channel first) ----
    nSets = numel(SpikeCounts.set);
    set_dist = nan(nSets,1);
    set_ord  = nan(nSets,1);
    for ss = 1:nSets
        sc = SpikeCounts.set(ss).stim_channels;
        sc = sc(sc > 0);
        if numel(sc) ~= 2, continue; end
        set_dist(ss) = SpikeCounts.set(ss).distance_um;
        if sc(1) < sc(2), set_ord(ss) = 1; else, set_ord(ss) = 2; end
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

    % ---- apply each row (distance, order, PTD) and each amplitude column ----
    used_keys = zeros(0,3);
    for r = 1:size(bad_rows,1)

        dist_v = bad_rows{r,1}; ord_v = bad_rows{r,2}; ptd_v = bad_rows{r,3};
        if ~(isnumeric(dist_v) && isscalar(dist_v) && isnumeric(ord_v) && isscalar(ord_v) && isnumeric(ptd_v) && isscalar(ptd_v))
            error('Dataset %d, row %d: Dist, Order and PTD must be single numbers.',d,r);
        end
        if ~any(ord_v == [1 2])
            error('Dataset %d, row %d: Order must be 1 (circle) or 2 (diamond).',d,r);
        end
        if any(all(abs(used_keys-[dist_v ord_v ptd_v]) < 1e-9,2))
            error('Dataset %d, row %d: the same Dist / Order / PTD appears in more than one row.',d,r);
        end
        used_keys(end+1,:) = [dist_v ord_v ptd_v]; %#ok<SAGROW>

        cand = find(abs(set_dist-dist_v) <= Distance_Tol_um & set_ord == ord_v);
        if isempty(cand)
            error(['Dataset %d, row %d: no stim set at %g um with Order %d.\n' ...
                'Sets in this dataset: %s'],d,r,dist_v,ord_v,list_sets(SpikeCounts,set_dist,set_ord));
        elseif numel(cand) > 1
            error('Dataset %d, row %d: %g um / Order %d matches more than one set (%s) - cannot tell them apart.', ...
                d,r,dist_v,ord_v,num2str(cand.'));
        end
        si = cand;
        amp_values = [SpikeCounts.set(si).amp.amp_value];

        for k = 1:numel(Amps)

            trials_in = bad_rows{r,3+k};
            if ~isnumeric(trials_in)
                error('Dataset %d, row %d, %g uA: the trial list must be numeric, e.g. [2 9] or [].',d,r,Amps(k));
            end

            ai = find(abs(amp_values-Amps(k)) < Condition_Tolerance,1);
            if isempty(ai)
                error('Dataset %d, row %d: amp %g uA not found (available: %s).',d,r,Amps(k),num2str(amp_values));
            end
            ptd_values = [SpikeCounts.set(si).amp(ai).ptd.PTD_ms];
            pi = find(abs(ptd_values-ptd_v) < Condition_Tolerance,1);
            if isempty(pi)
                error('Dataset %d, row %d (%g um, Order %d, %g uA): PTD %g ms does not exist here (available: %s).', ...
                    d,r,dist_v,ord_v,Amps(k),ptd_v,num2str(ptd_values));
            end
            cond = SpikeCounts.set(si).amp(ai).ptd(pi);
            total_trials = cond.total_trials;

            final_trials = validate_trial_list(trials_in,total_trials,d,r,Amps(k));

            condition_mask = TrialSummary.SetIndex == si & ...
                TrialSummary.AmplitudeIndex == ai & ...
                TrialSummary.PTDIndex == pi;

            absolute_trial_ids = zeros(size(final_trials));
            for t = 1:numel(final_trials)
                trial_mask = condition_mask & TrialSummary.ConditionTrialIndex == final_trials(t);
                if sum(trial_mask) ~= 1
                    error('Dataset %d, row %d, %g uA: condition trial %d could not be uniquely matched.', ...
                        d,r,Amps(k),final_trials(t));
                end
                absolute_trial_ids(t) = TrialSummary.AbsoluteTrialID(trial_mask);
            end
            absolute_trial_ids = sort(unique(absolute_trial_ids));

            if abs(cond.PTD_ms) < Condition_Tolerance, stim_type = 'simultaneous'; else, stim_type = 'sequential'; end

            R = empty_confirmed_struct();
            R(1).set = si;
            R(1).stimChannels = SpikeCounts.set(si).stim_channels;
            R(1).stimOrderLabel = SpikeCounts.set(si).set_name;
            R(1).stimType = stim_type;
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

    % ---- summary per condition ----
    fprintf('\nSummary - bad/total trials per condition ("-" = no saved entry):\n');
    fprintf('%-9s %-6s %-14s %-9s %-4s', 'Dist(um)','Order','Stim order','PTD(ms)','Set');
    for k = 1:numel(Amps), fprintf(' %14s',sprintf('%g uA',Amps(k))); end
    fprintf(' | %8s\n','excl %');

    [~,set_order_print] = sortrows([set_dist,set_ord]);
    missing = {};
    warn_list = {};
    dist_tot = containers.Map('KeyType','char','ValueType','any');   % key 'dist|ptd' -> [bad total]

    for q = 1:numel(set_order_print)
        ss = set_order_print(q);
        if isnan(set_dist(ss)), continue; end
        ptd_all = unique(cell2mat(arrayfun(@(a) [a.ptd.PTD_ms],SpikeCounts.set(ss).amp,'UniformOutput',false)));
        for pt = ptd_all(:).'
            fprintf('%-9g %-6d %-14s %-9g %-4d',set_dist(ss),set_ord(ss),SpikeCounts.set(ss).set_name,pt,ss);
            row_bad = 0; row_tot = 0;
            for k = 1:numel(Amps)
                hit = find([BadConfirmed.set] == ss & abs([BadConfirmed.amp]-Amps(k)) < Condition_Tolerance & ...
                    abs([BadConfirmed.ptd_ms]-pt) < Condition_Tolerance,1);
                if isempty(hit)
                    fprintf(' %14s','-');
                    ai_chk = find(abs([SpikeCounts.set(ss).amp.amp_value]-Amps(k)) < Condition_Tolerance,1);
                    if ~isempty(ai_chk)
                        pv = [SpikeCounts.set(ss).amp(ai_chk).ptd.PTD_ms];
                        if any(abs(pv-pt) < Condition_Tolerance)
                            missing{end+1} = sprintf('%g um, Order %d, PTD %g, %g uA',set_dist(ss),set_ord(ss),pt,Amps(k)); %#ok<SAGROW>
                        end
                    end
                else
                    nb = BadConfirmed(hit).n_bad; nt = BadConfirmed(hit).n_total_trials;
                    fprintf(' %14s',sprintf('%d/%d (%.0f%%)',nb,nt,100*nb/nt));
                    row_bad = row_bad+nb; row_tot = row_tot+nt;
                    if nt > 0 && nb/nt > Warn_Exclusion_Fraction
                        warn_list{end+1} = sprintf('%g um, Order %d, PTD %g, %g uA: %d/%d (%.0f%%)', ...
                            set_dist(ss),set_ord(ss),pt,Amps(k),nb,nt,100*nb/nt); %#ok<SAGROW>
                    end
                end
            end
            if row_tot > 0, fprintf(' | %7.0f%%\n',100*row_bad/row_tot); else, fprintf(' | %8s\n','-'); end
            key = sprintf('%g|%g',set_dist(ss),pt);
            if isKey(dist_tot,key), v = dist_tot(key); else, v = [0 0]; end
            dist_tot(key) = v+[row_bad row_tot];
        end
    end

    % ---- exclusion rate per distance and PTD ----
    fprintf('\nExclusion rate per distance and PTD (all stim orders and amplitudes with a saved entry):\n');
    keys_all = keys(dist_tot);
    parsed = nan(numel(keys_all),2);
    for q = 1:numel(keys_all)
        tk = sscanf(keys_all{q},'%f|%f');
        parsed(q,:) = tk(:).';
    end
    [~,oo] = sortrows(parsed);
    fprintf('%-9s %-9s %12s\n','Dist(um)','PTD(ms)','excluded');
    for q = oo(:).'
        v = dist_tot(keys_all{q});
        if v(2) > 0
            fprintf('%-9g %-9g %12s\n',parsed(q,1),parsed(q,2),sprintf('%d/%d (%.0f%%)',v(1),v(2),100*v(1)/v(2)));
        end
    end

    if ~isempty(warn_list)
        fprintf('\nWARNING - conditions that lost more than %.0f%% of their trials:\n',100*Warn_Exclusion_Fraction);
        fprintf('  %s\n',strjoin(warn_list,'\n  '));
    end

    fprintf('\nREVIEW CHECKLIST (paired conditions with NO saved entry):\n');
    if isempty(missing)
        fprintf('  none - every condition has been confirmed.\n');
    else
        fprintf('  %s\n',strjoin(missing,'\n  '));
        fprintf('  (these are NOT "no bad trials" - they have simply not been reviewed yet)\n');
    end
    fprintf('Total bad trials in this file: %d\n',numel(AllBadAbsoluteTrialIDs));
end

fprintf('\nDone. Re-run Clean_Paired_v1 to check what is left.\n');

%% =========================== LOCAL FUNCTIONS ==========================

function trials = validate_trial_list(trials,total_trials,dataset_index,row_index,amp_value)
if isempty(trials)
    trials = [];
    return;
end
trials = unique(double(trials(:).'));
invalid = ~isfinite(trials) | trials < 1 | trials > total_trials | fix(trials) ~= trials;
if any(invalid)
    error('Dataset %d, row %d, %g uA: invalid relative trial numbers (must be 1-%d): %s', ...
        dataset_index,row_index,amp_value,total_trials,num2str(trials(invalid)));
end
end

function txt = list_sets(SpikeCounts,set_dist,set_ord)
parts = {};
for ss = 1:numel(SpikeCounts.set)
    if isnan(set_dist(ss)), continue; end
    parts{end+1} = sprintf('%g um / Order %d (%s)',set_dist(ss),set_ord(ss),SpikeCounts.set(ss).set_name); %#ok<AGROW>
end
txt = strjoin(parts,'; ');
end

function R = empty_confirmed_struct()
R = struct('set',{},'stimChannels',{},'stimOrderLabel',{},'stimType',{}, ...
    'amp',{},'ptd_ms',{},'distance_um',{},'distance_note',{},'n_total_trials',{}, ...
    'condition_trial_indices',{},'absolute_trial_ids',{},'good_condition_trial_indices',{}, ...
    'n_bad',{},'confirmed_at',{});
end