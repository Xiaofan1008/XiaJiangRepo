%% =============================================================
%  Manual Random Fraction Spike Deletion for Final Pulse-Train Conditions
%
%  Purpose:
%    Practice/exploratory tool to manually remove spikes from selected
%    time windows in final recovered pulse-train conditions.
%
%  Rule format:
%    ManualDeleteRules = {
%    % DepthChannels      Family      SetID   Amp_uA   WinStart_ms   WinEnd_ms   DeleteFraction
%      [25 26 27],        'AutoSim',  1,      10,      19.8,         20.4,       1.0
%      [25 26 27],        'AutoSim',  1,      10,      24.8,         25.4,       0.5
%      [17:22 25 31],     'Seq',      2,      10,      14.8,         15.4,       1.0
%    };
%
%  Family:
%    'AutoSim' = final AutoSim/simultaneous condition only
%    'Seq'     = final sequential condition only
%    'All'     = both final AutoSim and final Seq conditions
%
%  NaN option:
%    SetID = NaN means all selected sets
%    Amp_uA = NaN means all selected amplitudes
%
%  Output if dry_run = true:
%    ManualDelete_Practice_DryRunSummary.mat
%
%  Output if dry_run = false:
%    *_ManualCleanPractice.mat
%
%  Saved variables when dry_run = false:
%    sp_manualclean
%    DeletionSummary
%    ManualDeleteRules
%    ManualCleanParams
% =============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ========================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Spike source:
%   'auto'         = prefer PrefixRecovery_SSD, then PrefixRecovery, then base
%   'recovery_ssd' = *_sp_xia_PrefixRecovery_SSD.mat, variable sp_corr
%   'recovery'     = *_sp_xia_PrefixRecovery.mat, variable sp_seq
%   'base'         = *_sp_xia.mat, variable sp_clipped or sp
spike_source = 'auto';

% Amplitudes to include when searching for final conditions.
% [] means all non-zero amplitudes.
Plot_Amps = [10];

% Stimulation set IDs to include.
% [] means all detected sets.
Plot_SetIDs = [];

% Families to include when building final condition list.
plot_auto_sim   = true;
plot_sequential = true;

% Include zero-current control?
include_zero_control = false;

% Dry run:
%   true  = only report candidate/deleted spike counts, do not save cleaned spikes
%   false = actually delete and save a new spike file
dry_run = false;

% Random seed for reproducible fractional deletion.
% If DeleteFraction = 0.5, this makes the random selection repeatable.
random_seed = 1;

%% ----------- Manual delete rules ------------ %%
% Each row is one complete rule.
%
% Columns:
%   1. DepthChannels
%   2. Family: 'AutoSim', 'Seq', or 'All'
%   3. SetID
%   4. Amp_uA
%   5. WinStart_ms
%   6. WinEnd_ms
%   7. DeleteFraction
%
ManualDeleteRules = {
% DepthChannels          Family      SetID   Amp_uA   WinStart_ms   WinEnd_ms   DeleteFraction
  [25],                  'AutoSim',  1,      10,      19.8,         20.4,       1.0
  [25],                  'AutoSim',  1,      10,      24.8,         25.4,       0.5
};






% Output names.
dryrun_summary_name = 'ManualDelete_Practice_DryRunSummary.mat';

%% ====================== CHECK FOLDER =========================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end

cd(data_folder);
fprintf('Manual deletion practice tool in folder:\n%s\n\n', data_folder);

rng(random_seed);

%% ====================== LOAD SPIKES ==========================

sp = [];
spike_file_used = '';
spike_variable_used = '';

switch lower(spike_source)

    case 'auto'

        ssd_recovery_file = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));
        prefix_file       = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
        prefix_file       = prefix_file(~contains({prefix_file.name}, 'SSD'));
        base_file         = dir(fullfile(data_folder, '*sp_xia.mat'));

        if ~isempty(ssd_recovery_file)

            spike_file_used = ssd_recovery_file(1).name;
            S_sp = load(fullfile(data_folder, spike_file_used));

            if isfield(S_sp, 'sp_corr')
                sp = S_sp.sp_corr;
                spike_variable_used = 'sp_corr';
            else
                error('PrefixRecovery_SSD file found but variable sp_corr was not found.');
            end

        elseif ~isempty(prefix_file)

            spike_file_used = prefix_file(1).name;
            S_sp = load(fullfile(data_folder, spike_file_used));

            if isfield(S_sp, 'sp_seq')
                sp = S_sp.sp_seq;
                spike_variable_used = 'sp_seq';
            else
                error('PrefixRecovery file found but variable sp_seq was not found.');
            end

        elseif ~isempty(base_file)

            spike_file_used = base_file(1).name;
            S_sp = load(fullfile(data_folder, spike_file_used));

            if isfield(S_sp, 'sp_clipped')
                sp = S_sp.sp_clipped;
                spike_variable_used = 'sp_clipped';
            elseif isfield(S_sp, 'sp')
                sp = S_sp.sp;
                spike_variable_used = 'sp';
            else
                error('Base spike file found but no usable spike variable was found.');
            end

        else
            error('No usable spike file found.');
        end

    case 'recovery_ssd'

        sp_file = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));
        assert(~isempty(sp_file), 'No *sp_xia_PrefixRecovery_SSD.mat file found.');

        spike_file_used = sp_file(1).name;
        S_sp = load(fullfile(data_folder, spike_file_used));

        if isfield(S_sp, 'sp_corr')
            sp = S_sp.sp_corr;
            spike_variable_used = 'sp_corr';
        else
            error('sp_corr not found in %s.', spike_file_used);
        end

    case 'recovery'

        sp_file = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
        sp_file = sp_file(~contains({sp_file.name}, 'SSD'));
        assert(~isempty(sp_file), 'No *sp_xia_PrefixRecovery.mat file found.');

        spike_file_used = sp_file(1).name;
        S_sp = load(fullfile(data_folder, spike_file_used));

        if isfield(S_sp, 'sp_seq')
            sp = S_sp.sp_seq;
            spike_variable_used = 'sp_seq';
        else
            error('sp_seq not found in %s.', spike_file_used);
        end

    case 'base'

        sp_file = dir(fullfile(data_folder, '*sp_xia.mat'));
        assert(~isempty(sp_file), 'No *sp_xia.mat file found.');

        spike_file_used = sp_file(1).name;
        S_sp = load(fullfile(data_folder, spike_file_used));

        if isfield(S_sp, 'sp_clipped')
            sp = S_sp.sp_clipped;
            spike_variable_used = 'sp_clipped';
        elseif isfield(S_sp, 'sp')
            sp = S_sp.sp;
            spike_variable_used = 'sp';
        else
            error('No usable spike variable found in %s.', spike_file_used);
        end

    otherwise
        error('Unknown spike_source: %s', spike_source);
end

nCh = numel(sp);

fprintf('Loaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);
fprintf('Loaded %d spike channels.\n', nCh);

% This is the spike variable that will be modified.
sp_manualclean = sp;

%% ====================== LOAD TRIGGERS ========================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD EXPERIMENT PARAMETERS ===========

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
TrialParams       = S_exp.TrialParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;

if isfield(S_exp, 'E_MAP')
    E_MAP = S_exp.E_MAP;
else
    E_MAP = [];
    warning('E_MAP not found. Stim labels may be less clean.');
end

if isfield(S_exp, 'StimMeta')
    StimMeta = S_exp.StimMeta;
else
    error('StimMeta was not found. This pulse-train delete tool requires StimMeta.');
end

fprintf('Loaded exp datafile: %s\n', fileDIR(1).name);
fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match nTrig (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== SAMPLING RATE ========================

try
    [~, freq_params] = read_Intan_RHS2000_file;
    FS = freq_params.amplifier_sample_rate;
catch
    FS = 30000;
    warning('Could not read info.rhs. Using FS = 30000 Hz.');
end

fprintf('Sampling rate: %.1f Hz\n', FS);

%% ====================== REMOVE HEADER ROWS ===================

StimParams_data  = StimParams(2:end,:);
TrialParams_data = TrialParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(StimParams_data,1), expected_rows);
end

if size(TrialParams_data,1) ~= expected_rows
    warning('TrialParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(TrialParams_data,1), expected_rows);
end

%% ====================== TRIAL CONDITION IDS ==================

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of trial-level condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA =================

stimSet_trial        = NaN(n_Trials,1);
trainLevel_trial     = NaN(n_Trials,1);
totalLevels_trial    = NaN(n_Trials,1); %#ok<NASGU>
trialAmps            = NaN(n_Trials,1);
isAutoSim_trial      = false(n_Trials,1);
isZeroControl        = false(n_Trials,1);
eventEnd_ms_trial    = NaN(n_Trials,1);

eventTimes_ms_trial  = cell(n_Trials,1);

for tr = 1:n_Trials

    condID = conditionID_trial(tr);

    if condID < 1 || condID > numel(StimMeta)
        warning('Trial %d has invalid condition ID %d.', tr, condID);
        continue;
    end

    meta = StimMeta(condID);

    if isfield(meta, 'StimSetIndex')
        stimSet_trial(tr) = meta.StimSetIndex;
    end

    if isfield(meta, 'TrainLevel')
        trainLevel_trial(tr) = meta.TrainLevel;
    end

    if isfield(meta, 'TotalTrainLevels')
        totalLevels_trial(tr) = meta.TotalTrainLevels; %#ok<NASGU>
    end

    if isfield(meta, 'IsAutoSimultaneous')
        isAutoSim_trial(tr) = logical(meta.IsAutoSimultaneous);
    end

    if isfield(meta, 'IsZeroCurrentControl')
        isZeroControl(tr) = logical(meta.IsZeroCurrentControl);
    end

    if isfield(meta, 'EventTimesIncluded_ms')
        eventTimes_ms_trial{tr} = meta.EventTimesIncluded_ms;
    else
        eventTimes_ms_trial{tr} = [];
    end

    if isfield(meta, 'EventEndTime_ms')
        eventEnd_ms_trial(tr) = meta.EventEndTime_ms;
    elseif ~isempty(eventTimes_ms_trial{tr})
        eventEnd_ms_trial(tr) = max(eventTimes_ms_trial{tr});
    else
        eventEnd_ms_trial(tr) = 0;
    end

    % Use actual randomized StimParams rows to get amplitude.
    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);

    if max(rr) <= size(StimParams_data,1)

        amp_vec = cell2mat(StimParams_data(rr,16));
        amp_vec = amp_vec(:)';

        if all(amp_vec <= 0)
            trialAmps(tr) = 0;
        else
            trialAmps(tr) = max(amp_vec(amp_vec > 0));
        end
    end
end

%% ====================== APPLY USER FILTERS ===================

% Amplitudes.
all_amps = unique(trialAmps(~isnan(trialAmps)));
all_amps = all_amps(all_amps > 0);
all_amps = sort(all_amps(:))';

if isempty(Plot_Amps)
    Plot_Amps_selected = all_amps;
else
    Plot_Amps_selected = intersect(all_amps, Plot_Amps);
end

% Sets.
all_sets = unique(stimSet_trial(~isnan(stimSet_trial)));
all_sets = all_sets(all_sets > 0);
all_sets = sort(all_sets(:))';

if isempty(Plot_SetIDs)
    SetIDs_selected = all_sets;
else
    SetIDs_selected = intersect(all_sets, Plot_SetIDs);
end

fprintf('\nSelected amplitudes: ');
disp(Plot_Amps_selected);

fprintf('Selected set IDs: ');
disp(SetIDs_selected);

fprintf('Use AutoSim final conditions:    %d\n', plot_auto_sim);
fprintf('Use sequential final conditions: %d\n', plot_sequential);

%% ====================== ELECTRODE MAP ========================

cur_dir = pwd;

try
    cd(data_folder);
    d = Depth_s(Electrode_Type);
    cd(cur_dir);
catch ME
    cd(cur_dir);
    warning('Depth_s failed. Falling back to direct mapping.');
    warning('%s', ME.message);

    % Fallback large enough for all requested channels in delete rules.
    maxDepthNeeded = getMaxDepthChannelFromRules(ManualDeleteRules);
    d = 1:maxDepthNeeded;
end

%% ====================== BUILD FINAL CONDITION LIST ============
% Each FinalCond is one final recovered condition:
%   AutoSim: highest TrainLevel for that amp/set/family
%   Seq:     highest TrainLevel for that amp/set/family

FinalConds = struct( ...
    'ConditionName', {}, ...
    'ConditionFamily', {}, ...
    'StimSetID', {}, ...
    'Amplitude_uA', {}, ...
    'FinalTrainLevel', {}, ...
    'TrialList', {}, ...
    'ConditionLabel', {}, ...
    'EventTimes_ms', {}, ...
    'EventEnd_ms', {});

for ai = 1:numel(Plot_Amps_selected)

    ampVal = Plot_Amps_selected(ai);

    %% ---------------- AutoSim final conditions ----------------
    if plot_auto_sim

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trials_base = find(abs(trialAmps - ampVal) < 1e-9 & ...
                               stimSet_trial == setID & ...
                               isAutoSim_trial == true);

            if ~include_zero_control
                trials_base = trials_base(isZeroControl(trials_base) == false);
            end

            trials_base = trials_base(trials_base <= nTrials_use);

            if isempty(trials_base)
                continue;
            end

            finalLevel = max(trainLevel_trial(trials_base));
            trials_final = trials_base(trainLevel_trial(trials_base) == finalLevel);

            if isempty(trials_final)
                continue;
            end

            tr_rep = trials_final(1);

            conditionLabel = getStimLabelForTrial( ...
                tr_rep, StimParams_data, simultaneous_stim, E_MAP, true);

            conditionName = sprintf('AutoSim | %s | Set %g | %.1f uA | final level %g', ...
                conditionLabel, setID, ampVal, finalLevel);

            C.ConditionName = conditionName;
            C.ConditionFamily = 'AutoSim';
            C.StimSetID = setID;
            C.Amplitude_uA = ampVal;
            C.FinalTrainLevel = finalLevel;
            C.TrialList = trials_final(:)';
            C.ConditionLabel = conditionLabel;
            C.EventTimes_ms = eventTimes_ms_trial{tr_rep};
            C.EventEnd_ms = eventEnd_ms_trial(tr_rep);

            FinalConds(end+1) = C; %#ok<SAGROW>
        end
    end

    %% ---------------- Sequential final conditions --------------
    if plot_sequential

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trials_base = find(abs(trialAmps - ampVal) < 1e-9 & ...
                               stimSet_trial == setID & ...
                               isAutoSim_trial == false);

            if ~include_zero_control
                trials_base = trials_base(isZeroControl(trials_base) == false);
            end

            trials_base = trials_base(trials_base <= nTrials_use);

            if isempty(trials_base)
                continue;
            end

            finalLevel = max(trainLevel_trial(trials_base));
            trials_final = trials_base(trainLevel_trial(trials_base) == finalLevel);

            if isempty(trials_final)
                continue;
            end

            tr_rep = trials_final(1);

            conditionLabel = getStimLabelForTrial( ...
                tr_rep, StimParams_data, simultaneous_stim, E_MAP, false);

            conditionName = sprintf('Seq | %s | Set %g | %.1f uA | final level %g', ...
                conditionLabel, setID, ampVal, finalLevel);

            C.ConditionName = conditionName;
            C.ConditionFamily = 'Seq';
            C.StimSetID = setID;
            C.Amplitude_uA = ampVal;
            C.FinalTrainLevel = finalLevel;
            C.TrialList = trials_final(:)';
            C.ConditionLabel = conditionLabel;
            C.EventTimes_ms = eventTimes_ms_trial{tr_rep};
            C.EventEnd_ms = eventEnd_ms_trial(tr_rep);

            FinalConds(end+1) = C; %#ok<SAGROW>
        end
    end
end

if isempty(FinalConds)
    error('No final conditions found. Check Plot_Amps, Plot_SetIDs, and metadata.');
end

fprintf('\n================ FINAL CONDITIONS AVAILABLE ================\n');
for ci = 1:numel(FinalConds)
    fprintf('%2d | %s | Trials %d\n', ...
        ci, FinalConds(ci).ConditionName, numel(FinalConds(ci).TrialList));
end
fprintf('============================================================\n');

%% ====================== VALIDATE RULES =======================

validateManualDeleteRules(ManualDeleteRules);

%% ====================== APPLY MANUAL DELETE RULES =============

Summary_RuleIndex = [];
Summary_DepthChannel = [];
Summary_ConditionName = {};
Summary_ConditionFamily = {};
Summary_ConditionLabel = {};
Summary_StimSetID = [];
Summary_Amplitude_uA = [];
Summary_FinalTrainLevel = [];
Summary_WindowStart_ms = [];
Summary_WindowEnd_ms = [];
Summary_DeleteFraction = [];
Summary_NumberOfTrials = [];
Summary_CandidateSpikeCount = [];
Summary_DeletedSpikeCount = [];
Summary_PredictedRemainingSpikeCount = [];
Summary_DryRun = [];

fprintf('\nApplying manual delete rules...\n');

for ri = 1:size(ManualDeleteRules, 1)

    ruleDepthChannels = ManualDeleteRules{ri, 1};
    ruleFamily        = ManualDeleteRules{ri, 2};
    ruleSetID         = ManualDeleteRules{ri, 3};
    ruleAmp           = ManualDeleteRules{ri, 4};
    winStart_ms       = ManualDeleteRules{ri, 5};
    winEnd_ms         = ManualDeleteRules{ri, 6};
    deleteFraction    = ManualDeleteRules{ri, 7};

    ruleFamilyNorm = lower(strtrim(ruleFamily));

    fprintf('\nRule %d:\n', ri);
    fprintf('  Depth channels: ');
    disp(ruleDepthChannels);
    fprintf('  Family: %s | SetID: %g | Amp: %g | Window: %.3f to %.3f ms | DeleteFraction: %.3f\n', ...
        ruleFamily, ruleSetID, ruleAmp, winStart_ms, winEnd_ms, deleteFraction);

    %% ----- Find matching final conditions -----
    matchedCondIdx = [];

    for ci = 1:numel(FinalConds)

        C = FinalConds(ci);

        familyMatch = false;

        switch ruleFamilyNorm
            case 'all'
                familyMatch = true;
            case 'seq'
                familyMatch = strcmpi(C.ConditionFamily, 'Seq');
            case 'autosim'
                familyMatch = strcmpi(C.ConditionFamily, 'AutoSim');
            otherwise
                error('Unknown family "%s" in rule %d. Use AutoSim, Seq, or All.', ...
                    ruleFamily, ri);
        end

        if ~familyMatch
            continue;
        end

        if ~isnan(ruleSetID) && C.StimSetID ~= ruleSetID
            continue;
        end

        if ~isnan(ruleAmp) && abs(C.Amplitude_uA - ruleAmp) > 1e-9
            continue;
        end

        matchedCondIdx(end+1) = ci; %#ok<AGROW>
    end

    if isempty(matchedCondIdx)
        warning('Rule %d did not match any final conditions.', ri);
        continue;
    end

    %% ----- Apply to each channel and matched condition -----
    for dd = 1:numel(ruleDepthChannels)

        depthCh = ruleDepthChannels(dd);

        if depthCh < 1 || depthCh > numel(d)
            warning('Rule %d: DepthChannel %d is outside Depth_s map. Skipped.', ...
                ri, depthCh);
            continue;
        end

        spikeCh = d(depthCh);

        if spikeCh < 1 || spikeCh > nCh
            warning('Rule %d: DepthChannel %d maps to spike channel %d, outside spike file range. Skipped.', ...
                ri, depthCh, spikeCh);
            continue;
        end

        if isempty(sp_manualclean{spikeCh})
            sp_times_current = [];
        else
            sp_times_current = sp_manualclean{spikeCh}(:,1);
        end

        for mm = 1:numel(matchedCondIdx)

            ci = matchedCondIdx(mm);
            C = FinalConds(ci);

            trials_this = C.TrialList;
            trials_this = trials_this(trials_this <= nTrials_use);

            if isempty(trials_this)
                continue;
            end

            %% ----- Find candidate spikes in this channel/condition/window -----
            candidateIdx = [];

            if ~isempty(sp_times_current)

                for ti = 1:numel(trials_this)

                    tr = trials_this(ti);

                    if tr > numel(trig)
                        continue;
                    end

                    t0_ms = trig(tr) / FS * 1000;

                    absStart = t0_ms + winStart_ms;
                    absEnd   = t0_ms + winEnd_ms;

                    idx_this = find(sp_times_current >= absStart & ...
                                    sp_times_current <  absEnd);

                    candidateIdx = [candidateIdx; idx_this(:)]; %#ok<AGROW>
                end
            end

            candidateIdx = unique(candidateIdx);
            nCandidate = numel(candidateIdx);

            if nCandidate == 0
                nDelete = 0;
                deleteIdx = [];
            elseif deleteFraction >= 1
                nDelete = nCandidate;
                deleteIdx = candidateIdx;
            elseif deleteFraction <= 0
                nDelete = 0;
                deleteIdx = [];
            else
                nDelete = round(deleteFraction * nCandidate);
                nDelete = min(max(nDelete, 0), nCandidate);

                if nDelete > 0
                    randOrder = randperm(nCandidate, nDelete);
                    deleteIdx = candidateIdx(randOrder);
                else
                    deleteIdx = [];
                end
            end

            predictedRemaining = nCandidate - nDelete;

            %% ----- Actually delete if not dry run -----
            if ~dry_run && nDelete > 0

                deleteMask = false(size(sp_manualclean{spikeCh}, 1), 1);
                deleteMask(deleteIdx) = true;

                sp_manualclean{spikeCh}(deleteMask, :) = [];

                % Update current spike times for this channel, because later
                % matched conditions/rules should operate on the already
                % modified spike list.
                if isempty(sp_manualclean{spikeCh})
                    sp_times_current = [];
                else
                    sp_times_current = sp_manualclean{spikeCh}(:,1);
                end
            end

            %% ----- Record summary row -----
            Summary_RuleIndex = [Summary_RuleIndex; ri]; %#ok<AGROW>
            Summary_DepthChannel = [Summary_DepthChannel; depthCh]; %#ok<AGROW>
            Summary_ConditionName = [Summary_ConditionName; {C.ConditionName}]; %#ok<AGROW>
            Summary_ConditionFamily = [Summary_ConditionFamily; {C.ConditionFamily}]; %#ok<AGROW>
            Summary_ConditionLabel = [Summary_ConditionLabel; {C.ConditionLabel}]; %#ok<AGROW>
            Summary_StimSetID = [Summary_StimSetID; C.StimSetID]; %#ok<AGROW>
            Summary_Amplitude_uA = [Summary_Amplitude_uA; C.Amplitude_uA]; %#ok<AGROW>
            Summary_FinalTrainLevel = [Summary_FinalTrainLevel; C.FinalTrainLevel]; %#ok<AGROW>
            Summary_WindowStart_ms = [Summary_WindowStart_ms; winStart_ms]; %#ok<AGROW>
            Summary_WindowEnd_ms = [Summary_WindowEnd_ms; winEnd_ms]; %#ok<AGROW>
            Summary_DeleteFraction = [Summary_DeleteFraction; deleteFraction]; %#ok<AGROW>
            Summary_NumberOfTrials = [Summary_NumberOfTrials; numel(trials_this)]; %#ok<AGROW>
            Summary_CandidateSpikeCount = [Summary_CandidateSpikeCount; nCandidate]; %#ok<AGROW>
            Summary_DeletedSpikeCount = [Summary_DeletedSpikeCount; nDelete]; %#ok<AGROW>
            Summary_PredictedRemainingSpikeCount = [Summary_PredictedRemainingSpikeCount; predictedRemaining]; %#ok<AGROW>
            Summary_DryRun = [Summary_DryRun; dry_run]; %#ok<AGROW>

            fprintf('  Ch %d | %s | candidates %d | delete %d | remain %d\n', ...
                depthCh, C.ConditionName, nCandidate, nDelete, predictedRemaining);
        end
    end
end

%% ====================== CREATE DELETION SUMMARY TABLE =========

DeletionSummary = table( ...
    Summary_RuleIndex, ...
    Summary_DepthChannel, ...
    Summary_ConditionName, ...
    Summary_ConditionFamily, ...
    Summary_ConditionLabel, ...
    Summary_StimSetID, ...
    Summary_Amplitude_uA, ...
    Summary_FinalTrainLevel, ...
    Summary_WindowStart_ms, ...
    Summary_WindowEnd_ms, ...
    Summary_DeleteFraction, ...
    Summary_NumberOfTrials, ...
    Summary_CandidateSpikeCount, ...
    Summary_DeletedSpikeCount, ...
    Summary_PredictedRemainingSpikeCount, ...
    Summary_DryRun, ...
    'VariableNames', { ...
        'RuleIndex', ...
        'DepthChannel', ...
        'ConditionName', ...
        'ConditionFamily', ...
        'ConditionLabel', ...
        'StimSetID', ...
        'Amplitude_uA', ...
        'FinalTrainLevel', ...
        'WindowStart_ms', ...
        'WindowEnd_ms', ...
        'DeleteFraction', ...
        'NumberOfTrials', ...
        'CandidateSpikeCount', ...
        'DeletedSpikeCount', ...
        'PredictedRemainingSpikeCount', ...
        'DryRun'});

%% ====================== PARAMS STRUCT =========================

ManualCleanParams = struct();

ManualCleanParams.DataFolder = data_folder;
ManualCleanParams.InputSpikeFile = spike_file_used;
ManualCleanParams.InputSpikeVariable = spike_variable_used;
ManualCleanParams.ExpDataFileUsed = fileDIR(1).name;

ManualCleanParams.Electrode_Type = Electrode_Type;
ManualCleanParams.DepthToSpikeChannelMap = d;

ManualCleanParams.Plot_Amps_Input = Plot_Amps;
ManualCleanParams.Plot_Amps_Selected = Plot_Amps_selected;
ManualCleanParams.Plot_SetIDs_Input = Plot_SetIDs;
ManualCleanParams.SetIDs_Selected = SetIDs_selected;

ManualCleanParams.PlotAutoSim = plot_auto_sim;
ManualCleanParams.PlotSequential = plot_sequential;
ManualCleanParams.IncludeZeroControl = include_zero_control;
ManualCleanParams.OnlyFinalTrainLevel = true;

ManualCleanParams.DryRun = dry_run;
ManualCleanParams.RandomSeed = random_seed;
ManualCleanParams.RandomDeletionMode = 'random_candidate_spikes_with_fixed_seed';
ManualCleanParams.Note = ['Practice/manual artifact cleaning tool. Fractional deletion ', ...
                          'randomly removes candidate spikes in user-defined windows. ', ...
                          'Original spike file is not overwritten.'];

ManualCleanParams.FS = FS;
ManualCleanParams.nTrials_use = nTrials_use;
ManualCleanParams.DateGenerated = datestr(now);

%% ====================== SAVE OUTPUT ===========================

if dry_run

    dryrun_path = fullfile(data_folder, dryrun_summary_name);

    save(dryrun_path, ...
        'ManualDeleteRules', ...
        'DeletionSummary', ...
        'ManualCleanParams', ...
        'FinalConds', ...
        '-v7.3');

    fprintf('\nDRY RUN ONLY. No spikes were deleted.\n');
    fprintf('Saved dry-run summary to:\n%s\n', dryrun_path);

else

    [~, inBase, ~] = fileparts(spike_file_used);
    out_file_name = [inBase '_Clean.mat'];
    out_file_path = fullfile(data_folder, out_file_name);

    % Save the cleaned spike variable.
    sp_corr = sp_manualclean;
    
    save(out_file_path, ...
        'sp_corr', ...
        'ManualDeleteRules', ...
        'DeletionSummary', ...
        'ManualCleanParams', ...
        'FinalConds', ...
        '-v7.3');

    fprintf('\nManual practice cleaning complete.\n');
    fprintf('Original spike file was NOT overwritten.\n');
    fprintf('Saved manually cleaned spike file to:\n%s\n', out_file_path);
end

fprintf('\nDeletion summary rows: %d\n', height(DeletionSummary));
fprintf('Total candidate spikes: %d\n', sum(DeletionSummary.CandidateSpikeCount));
fprintf('Total deleted spikes:   %d\n', sum(DeletionSummary.DeletedSpikeCount));

fprintf('\nFinished manual random-fraction deletion tool.\n');

%% =============================================================
%  LOCAL HELPER FUNCTIONS
%% =============================================================

function validateManualDeleteRules(ManualDeleteRules)

    if isempty(ManualDeleteRules)
        error('ManualDeleteRules is empty. Please define at least one rule.');
    end

    if size(ManualDeleteRules, 2) ~= 7
        error(['ManualDeleteRules must have 7 columns: ', ...
               'DepthChannels, Family, SetID, Amp_uA, WinStart_ms, WinEnd_ms, DeleteFraction.']);
    end

    for ri = 1:size(ManualDeleteRules, 1)

        depthChs = ManualDeleteRules{ri, 1};
        family   = ManualDeleteRules{ri, 2};
        winStart = ManualDeleteRules{ri, 5};
        winEnd   = ManualDeleteRules{ri, 6};
        frac     = ManualDeleteRules{ri, 7};

        if isempty(depthChs) || ~isnumeric(depthChs)
            error('Rule %d: DepthChannels must be a numeric vector.', ri);
        end

        if ~(ischar(family) || isstring(family))
            error('Rule %d: Family must be AutoSim, Seq, or All.', ri);
        end

        familyNorm = lower(strtrim(char(family)));

        if ~ismember(familyNorm, {'autosim', 'seq', 'all'})
            error('Rule %d: Family must be AutoSim, Seq, or All.', ri);
        end

        if ~isnumeric(winStart) || ~isnumeric(winEnd) || winEnd <= winStart
            error('Rule %d: invalid window. WinEnd_ms must be larger than WinStart_ms.', ri);
        end

        if ~isnumeric(frac) || frac < 0 || frac > 1
            error('Rule %d: DeleteFraction must be between 0 and 1.', ri);
        end
    end
end

function maxDepthNeeded = getMaxDepthChannelFromRules(ManualDeleteRules)

    maxDepthNeeded = 1;

    if isempty(ManualDeleteRules)
        return;
    end

    for ri = 1:size(ManualDeleteRules, 1)
        chs = ManualDeleteRules{ri, 1};

        if isnumeric(chs) && ~isempty(chs)
            maxDepthNeeded = max(maxDepthNeeded, max(chs(:)));
        end
    end
end

function label = getStimLabelForTrial(tr, StimParams_data, simultaneous_stim, E_MAP, isAutoSim)
    % Build stimulation label for one representative trial.

    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);
    rr = rr(rr <= size(StimParams_data,1));

    if isempty(rr)
        label = sprintf('Trial%d', tr);
        return;
    end

    stimNames = StimParams_data(rr,1);

    try
        ampVec = cell2mat(StimParams_data(rr,16));
    catch
        ampVec = ones(numel(rr),1);
    end

    activeRows = ampVec > 0;
    stimNames_active = stimNames(activeRows);

    label = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim);
end

function setLabel = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim)
    % Build stimulation label for condition identification.
    %
    % AutoSim:
    %   Ch35+Ch39
    %
    % Sequential:
    %   Ch35→Ch39

    if isempty(stimNames_active)
        setLabel = 'NoActiveStim';
        return;
    end

    stimNames_active = unique(stimNames_active, 'stable');

    labelParts = cell(1, numel(stimNames_active));

    for i = 1:numel(stimNames_active)

        stimName = stimNames_active{i};

        chNum = convertStimNameUsingEMap(stimName, E_MAP);

        if isnan(chNum)
            labelParts{i} = sprintf('%s', stimName);
        else
            labelParts{i} = sprintf('Ch%d', chNum);
        end
    end

    if isAutoSim
        setLabel = strjoin(labelParts, '+');
    else
        setLabel = strjoin(labelParts, '→');
    end
end

function chNum = convertStimNameUsingEMap(stimName, E_MAP)
    % Convert Intan stimulation label to channel number using E_MAP.
    %
    % E_MAP convention:
    %   E_MAP{1,1} = header / map name
    %   E_MAP{2,1} = channel 1
    %   E_MAP{3,1} = channel 2
    %
    % Therefore:
    %   channel number = row index - 1

    chNum = NaN;

    if isempty(stimName)
        return;
    end

    if isstring(stimName)
        stimName = char(stimName);
    end

    stimName = strtrim(stimName);

    if isempty(E_MAP)
        chNum = parseNumberFromStimName(stimName);
        return;
    end

    if iscell(E_MAP)

        for r = 1:size(E_MAP,1)

            thisName = E_MAP{r,1};

            if isstring(thisName)
                thisName = char(thisName);
            end

            if ischar(thisName)

                thisName = strtrim(thisName);

                if strcmp(thisName, stimName)
                    chNum = r - 1;
                    return;
                end
            end
        end

    elseif isstring(E_MAP)

        E_list = cellstr(E_MAP);

        for r = 1:numel(E_list)
            if strcmp(strtrim(E_list{r}), stimName)
                chNum = r - 1;
                return;
            end
        end

    elseif ischar(E_MAP)

        E_list = cellstr(E_MAP);

        for r = 1:numel(E_list)
            if strcmp(strtrim(E_list{r}), stimName)
                chNum = r - 1;
                return;
            end
        end
    end

    chNum = parseNumberFromStimName(stimName);
end

function chNum = parseNumberFromStimName(stimName)
    % Fallback if E_MAP lookup fails.

    chNum = NaN;

    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
    end
end