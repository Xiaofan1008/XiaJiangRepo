%% ========================================================================
% BAD-TRIAL ENTRY FOR THE SIM / SEQ (PAIRED) DATASETS - ONE AMPLITUDE - v2
%
% PURPOSE
%   Record the final bad trials of ONE amplitude of the paired datasets, from the
%   template printed by Clean_Paired_v2.m (section 5). The result goes to its own
%   file, tagged with the amplitude, e.g.  <base_name>_10uA_BadTrials.mat
%   Singles are NOT touched (they keep the untagged _BadTrials.mat).
%   Run it once per amplitude (change Amp only), always with the same rules.
%
% HOW A ROW READS
%   Amp = 10;                                  % the amplitude of this run (uA)
%   distance_bad_trials(1).matFilePath = '.../<base_name>_TrialSpikeCounts.mat';
%   distance_bad_trials(1).bad = {
%     % Dist  Order  PTD   trials (10 uA)
%       100,   1,     0,    [4];        % Ch1 -> Ch3 (circle)  Sim
%       100,   1,     5,    [];         % Ch1 -> Ch3 (circle)  Seq
%       100,   2,     5,    [1 5];      % Ch3 -> Ch1 (diamond) Seq
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
% TAG
%   Tag = 'auto' -> '<Amp>uA' (10 -> '10uA', 2.5 -> '2.5uA'). Use the same tag in
%   Clean_Paired_v2 (BadTrial_Tag) so the check reads this file. A trial of one
%   amplitude can never affect another amplitude, so the files stay independent.
%
% OUTPUT
%   <base_name>_<tag>_BadTrials.mat in each dataset folder, same format as
%   BadTrial_ManualConfirm_v1 (BadConfirmed, AllBadAbsoluteTrialIDs,
%   SourceTrialCountFile). A timestamped backup of the previous file is kept.
%   The bad trials of a paired file only change the ACTUAL; the prediction is
%   unaffected. After saving, per dataset it prints
%     - bad / total trials of every condition at this amplitude, with exclusion rate
%     - the exclusion rate per distance and PTD (does one distance lose far
%       more than the others?) and a warning above Warn_Exclusion_Fraction
%     - a REVIEW CHECKLIST of conditions at this amplitude with no saved entry
% ========================================================================

clear; clc;

%% ============================ USER SETTINGS ===========================

Condition_Tolerance = 1e-4;
Distance_Tol_um = 1;
Create_Backup_If_Output_Exists = false;
Warn_Exclusion_Fraction = 0.30;   % warn when a condition loses more than this share of its trials

Amp = 5;                         % the ONE amplitude of this run (uA)
Tag = 'auto';                     % 'auto' -> '<Amp>uA'; or type your own, e.g. '10uA'

distance_bad_trials = struct('matFilePath',{},'bad',{});

distance_bad_trials(1).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_SimSeq1/Linearity_S1_100_400um_SimSeq1_TrialSpikeCounts.mat';
distance_bad_trials(1).bad = {
  % Dist  Order  PTD   trials (5 uA)
    100,   1,   0,   [1 3 4 7 9 11 16:19 20 25 27];   % Ch1 -> Ch3 (circle)  Sim  [0 saved + 7 new]
    100,   1,   5,   [2:4 11 13 16:18 22 25 26 27 29 30];   % Ch1 -> Ch3 (circle)  Seq 5 ms  [0 saved + 7 new]
    100,   2,   5,   [1 3 4 5 6 8 9 10 11 13 15 17 18 20 21 24 25 26 28 30];   % Ch3 -> Ch1 (diamond)  Seq 5 ms  [0 saved + 5 new]

    200,   1,   0,   [1 3 4 7:9 11:13 24:26 28 30];   % Ch1 -> Ch5 (circle)  Sim  [0 saved + 9 new]
    200,   1,   5,   [6 7 8 10 11 12 13 17 19 21 22 23 27 29];   % Ch1 -> Ch5 (circle)  Seq 5 ms  [0 saved + 4 new]
    200,   2,   5,   [6 8 10 11 13 14 16 18 20 21 23 24 29 30];   % Ch5 -> Ch1 (diamond)  Seq 5 ms  [0 saved + 5 new]

    300,   1,   0,   [1 4 5 6 10 11 14 15 16 18 20 23 25 28 29 30];   % Ch1 -> Ch7 (circle)  Sim  [0 saved + 6 new]
    300,   1,   5,   [1 11 14 16 24 28];   % Ch1 -> Ch7 (circle)  Seq 5 ms  [0 saved + 4 new]
    300,   2,   5,   [2 6 11 12 14 15 17 19 21 27];   % Ch7 -> Ch1 (diamond)  Seq 5 ms  [0 saved + 6 new]

    400,   1,   0,   [4 7 9 14:16 23:26 30];   % Ch1 -> Ch9 (circle)  Sim  [0 saved + 5 new]
    400,   1,   5,   [2 5 6 7 8 10 14 15 17 19 21 23 24 25 26 27 28 30];   % Ch1 -> Ch9 (circle)  Seq 5 ms  [0 saved + 4 new]
    400,   2,   5,   [5 7 8 11 13 14 19 30];   % Ch9 -> Ch1 (diamond)  Seq 5 ms  [0 saved + 3 new]
};

distance_bad_trials(2).matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_SimSeq1/Linearity_S1_500_700um_SimSeq1_TrialSpikeCounts.mat';
distance_bad_trials(2).bad = {
  % Dist  Order  PTD   trials (5 uA)
    500,   1,   0,   [2 6 9 11 12 13 14 17 18 19 21 22 24 26 27 30];   % Ch1 -> Ch11 (circle)  Sim  [0 saved + 9 new]
    500,   1,   5,   [2 5 6 8 10 11 16];   % Ch1 -> Ch11 (circle)  Seq 5 ms  [0 saved + 4 new]
    500,   2,   5,   [1 2 3 4 9 13 14 18 20 22 30];   % Ch11 -> Ch1 (diamond)  Seq 5 ms  [0 saved + 8 new]

    600,   1,   0,   [2 4 5 6 7 13 14 15 16 18 19 20 22 23:25 26 27 28];   % Ch1 -> Ch13 (circle)  Sim  [0 saved + 7 new]
    600,   1,   5,   [1 2 4 5 7 9 10 12 17 30];   % Ch1 -> Ch13 (circle)  Seq 5 ms  [0 saved + 6 new]
    600,   2,   5,   [1 2 3 4 6 8 9 10 12 13 15 16 17 19 21 25 26 30];   % Ch13 -> Ch1 (diamond)  Seq 5 ms  [0 saved + 9 new]

    700,   1,   0,   [1 2 3 10 12 14 15 17 19 20 21 22 24 26:29 30];   % Ch1 -> Ch15 (circle)  Sim  [0 saved + 7 new]
    700,   1,   5,   [1 2 3 4 5 6 12 14 15 17 19 20 23 30];   % Ch1 -> Ch15 (circle)  Seq 5 ms  [0 saved + 8 new]
    700,   2,   5,   [1 2 4 6 13 14 15 16 17 18 19 20 23 24 26 29];   % Ch15 -> Ch1 (diamond)  Seq 5 ms  [0 saved + 4 new]
};

%% =========================== APPLY EACH DATASET ========================

if isempty(distance_bad_trials), error('distance_bad_trials is empty.'); end
if ~(isnumeric(Amp) && isscalar(Amp) && Amp > 0), error('Amp must be one positive number (uA).'); end
Amp = double(Amp);
if ischar(Tag) && strcmpi(Tag,'auto'), Tag = sprintf('%guA',Amp); end
if ~(ischar(Tag) && ~isempty(Tag)), error('Tag must be a non-empty text (or ''auto'').'); end

fprintf('\n============================================================\n');
fprintf('PAIRED BAD TRIALS (one amplitude, v2) - %d dataset(s) | %g uA | tag ''%s''\n', ...
    numel(distance_bad_trials),Amp,Tag);
fprintf('============================================================\n');

for d = 1:numel(distance_bad_trials)

    matFilePath = distance_bad_trials(d).matFilePath;
    bad_rows    = distance_bad_trials(d).bad;

    fprintf('\n--- Dataset %d ---\n%s\n',d,matFilePath);

    if ~isfile(matFilePath)
        error('Dataset %d: file does not exist:\n%s',d,matFilePath);
    end
    if isempty(bad_rows) || size(bad_rows,2) ~= 4
        error('Dataset %d: bad needs 4 columns {Dist, Order, PTD, [trials]} and at least one row (this script handles one amplitude: Amp = %g).', ...
            d,Amp);
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
    output_path = fullfile(dataset_folder,sprintf('%s_%s_BadTrials.mat',base_name,Tag));

    BadConfirmed = empty_confirmed_struct();
    if isfile(output_path)
        ExistingLoad = load(output_path,'BadConfirmed');
        if isfield(ExistingLoad,'BadConfirmed')
            BadConfirmed = ExistingLoad.BadConfirmed;
            fprintf('Existing %s file found (%d condition(s) already confirmed).\n',Tag,numel(BadConfirmed));
        else
            warning('Existing output file has no BadConfirmed - starting fresh.');
        end
    else
        fprintf('No existing %s file - this run creates it.\n',Tag);
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

        for k = 1   % single amplitude (loop kept so the block below is unchanged)

            trials_in = bad_rows{r,4};
            if ~isnumeric(trials_in)
                error('Dataset %d, row %d, %g uA: the trial list must be numeric, e.g. [2 9] or [].',d,r,Amp);
            end

            ai = find(abs(amp_values-Amp) < Condition_Tolerance,1);
            if isempty(ai)
                error('Dataset %d, row %d: amp %g uA not found (available: %s).',d,r,Amp,num2str(amp_values));
            end
            ptd_values = [SpikeCounts.set(si).amp(ai).ptd.PTD_ms];
            pi = find(abs(ptd_values-ptd_v) < Condition_Tolerance,1);
            if isempty(pi)
                error('Dataset %d, row %d (%g um, Order %d, %g uA): PTD %g ms does not exist here (available: %s).', ...
                    d,r,dist_v,ord_v,Amp,ptd_v,num2str(ptd_values));
            end
            cond = SpikeCounts.set(si).amp(ai).ptd(pi);
            total_trials = cond.total_trials;

            final_trials = validate_trial_list(trials_in,total_trials,d,r,Amp);

            condition_mask = TrialSummary.SetIndex == si & ...
                TrialSummary.AmplitudeIndex == ai & ...
                TrialSummary.PTDIndex == pi;

            absolute_trial_ids = zeros(size(final_trials));
            for t = 1:numel(final_trials)
                trial_mask = condition_mask & TrialSummary.ConditionTrialIndex == final_trials(t);
                if sum(trial_mask) ~= 1
                    error('Dataset %d, row %d, %g uA: condition trial %d could not be uniquely matched.', ...
                        d,r,Amp,final_trials(t));
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
        backup_path = fullfile(dataset_folder,sprintf('%s_%s_BadTrials_BACKUP_%s.mat', ...
            base_name,Tag,datestr(now,'yyyymmdd_HHMMSS')));
        copyfile(output_path,backup_path);
        fprintf('Backup of previous file: %s\n',backup_path);
    end

    SourceTrialCountFile = matFilePath;
    save(output_path,'BadConfirmed','AllBadAbsoluteTrialIDs','SourceTrialCountFile','-v7.3');
    fprintf('Saved: %s\n',output_path);

    % ---- summary per condition ----
    fprintf('\nSummary at %g uA - bad/total trials per condition ("-" = no saved entry):\n',Amp);
    fprintf('%-9s %-6s %-14s %-9s %-4s %14s\n','Dist(um)','Order','Stim order','PTD(ms)','Set',sprintf('%g uA',Amp));

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
            for k = 1
                hit = find([BadConfirmed.set] == ss & abs([BadConfirmed.amp]-Amp) < Condition_Tolerance & ...
                    abs([BadConfirmed.ptd_ms]-pt) < Condition_Tolerance,1);
                if isempty(hit)
                    fprintf(' %14s','-');
                    ai_chk = find(abs([SpikeCounts.set(ss).amp.amp_value]-Amp) < Condition_Tolerance,1);
                    if ~isempty(ai_chk)
                        pv = [SpikeCounts.set(ss).amp(ai_chk).ptd.PTD_ms];
                        if any(abs(pv-pt) < Condition_Tolerance)
                            missing{end+1} = sprintf('%g um, Order %d, PTD %g, %g uA',set_dist(ss),set_ord(ss),pt,Amp); %#ok<SAGROW>
                        end
                    end
                else
                    nb = BadConfirmed(hit).n_bad; nt = BadConfirmed(hit).n_total_trials;
                    fprintf(' %14s',sprintf('%d/%d (%.0f%%)',nb,nt,100*nb/nt));
                    row_bad = row_bad+nb; row_tot = row_tot+nt;
                    if nt > 0 && nb/nt > Warn_Exclusion_Fraction
                        warn_list{end+1} = sprintf('%g um, Order %d, PTD %g, %g uA: %d/%d (%.0f%%)', ...
                            set_dist(ss),set_ord(ss),pt,Amp,nb,nt,100*nb/nt); %#ok<SAGROW>
                    end
                end
            end
            fprintf('\n');
            key = sprintf('%g|%g',set_dist(ss),pt);
            if isKey(dist_tot,key), v = dist_tot(key); else, v = [0 0]; end
            dist_tot(key) = v+[row_bad row_tot];
        end
    end

    % ---- exclusion rate per distance and PTD ----
    fprintf('\nExclusion rate per distance and PTD (all stim orders, %g uA, conditions with a saved entry):\n',Amp);
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

    fprintf('\nREVIEW CHECKLIST (paired conditions at %g uA with NO saved entry):\n',Amp);
    if isempty(missing)
        fprintf('  none - every condition has been confirmed.\n');
    else
        fprintf('  %s\n',strjoin(missing,'\n  '));
        fprintf('  (these are NOT "no bad trials" - they have simply not been reviewed yet)\n');
    end
    fprintf('Total bad trials in this file: %d\n',numel(AllBadAbsoluteTrialIDs));
end

fprintf('\nDone. Re-run Clean_Paired_v2 (Amp_To_Clean = %g) to check what is left.\n',Amp);

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