%% Raster + PSTH Viewer for Filtered Pulse-Train Recovery Results
% ------------------------------------------------------------
% Purpose:
%   Plot raster + PSTH for pulse-train / event-level stimulation files.
%
% Main use here:
%   Check filtered recovered spike result:
%       *.sp_xia_PrefixRecovery_SSD.mat
%       variable: sp_corr
% ------------------------------------------------------------

clear all
% close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== CHOOSE FOLDER ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Single2';

%% ====================== USER SETTINGS ======================

raster_chn_start = 33;
raster_chn_end   = 64;   % Depth_s index
Electrode_Type   = 2;   % 0: rigid; 1: single-shank flex; 2: four-shank flex

% Stimulation set selection.
% [] means all detected sets.
SetIDs_to_plot = [];

% Amplitudes to plot.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% ============================================================
% IMPORTANT LEVEL SETTINGS
%
% Sequential/interleaved and AutoSim/simultaneous have different final levels.
%
% Sequential example:
%   Level 1 = [0]
%   Level 2 = [0 5]
%   ...
%   Level 6 = [0 5 10 15 20 25]
%
% AutoSim example:
%   Level 1 = [0]
%   Level 2 = [0 10]
%   Level 3 = [0 10 20]
%
% Therefore:
%   use 6 for final sequential result
%   use 3 for final AutoSim result
% ============================================================

% Sequential/interleaved train levels to plot.
% [] means all detected sequential levels.
seq_train_levels_to_plot = 3;

% AutoSim/simultaneous train levels to plot.
% [] means all detected AutoSim levels.
autosim_train_levels_to_plot = 3;

% Plot sequential / interleaved condition?
plot_sequential = true;

% Plot AutoSim / simultaneous condition?
plot_auto_sim = true;

% Include zero-current control trials?
% Usually 0 for this recovery-check viewer.
include_zero_control = 0;

% Raster/PSTH parameters.
ras_win        = [-60 80];   % ms
bin_ms_raster = 1;            % PSTH bin size, ms
smooth_ms      = 5;           % PSTH smoothing width

% Spike file source:
%   'raw'          = load *.sp.mat, variable sp
%   'xia'          = load *.sp_xia.mat, variable sp_clipped
%   'ssd'          = load *.sp_xia_SSD.mat, variable sp_corr
%   'recovery'     = load *.sp_xia_PrefixRecovery.mat, variable sp_seq
%   'recovery_ssd' = load *.sp_xia_PrefixRecovery_SSD.mat, variable sp_corr
%   'auto'         = priority: recovery_ssd > recovery > ssd > xia > raw
%
% For checking the filtered recovered spike result, use:
spike_source = 'recovery_ssd';

% This is only used when spike_source = 'raw'.
% For recovery_ssd checking, it does nothing.
save_sp_xia = false;

% If true, overwrite existing *.sp_xia.mat when spike_source = 'raw'.
overwrite_sp_xia = false;

% If true, print one example trial for each selected condition.
debug_print_trial_content = true;

%% ====================== CHECK FOLDER ======================

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end

cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ====================== BASE NAME ======================

parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');

if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ====================== SAMPLING RATE ======================

FS = 30000;

%% ====================== LOAD EXPERIMENT DATAFILE ======================

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
TrialParams       = S_exp.TrialParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;
E_MAP             = S_exp.E_MAP;

if isfield(S_exp, 'StimMeta')
    StimMeta = S_exp.StimMeta;
else
    error('StimMeta was not found. This pulse-train raster viewer requires StimMeta.');
end

if isfield(S_exp, 'trainMode')
    trainMode = S_exp.trainMode;
else
    trainMode = NaN;
end

fprintf('Loaded exp datafile: %s\n', fileDIR(1).name);
fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);
fprintf('trainMode: %g\n', trainMode);

%% ====================== LOAD SPIKES ======================
% All loaded spike data are assigned to sp_clipped internally.
% This keeps the rest of the plotting code unchanged.
%
% For example:
%   recovery_ssd loads sp_corr, then:
%       sp_clipped = sp_corr;

sp_clipped = [];

switch lower(spike_source)

    case 'raw'
        sp_files = dir(fullfile(data_folder, '*.sp.mat'));
        sp_files = sp_files(~contains({sp_files.name}, 'sp_xia'));

        assert(~isempty(sp_files), 'No raw .sp.mat file found in the current folder.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading raw spike file: %s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp')
            sp = S_sp.sp;
        else
            error('Variable "sp" not found in %s.', sp_filename);
        end

        sp_clipped = sp;

        if save_sp_xia
            out_xia = [base_name '.sp_xia.mat'];

            if exist(out_xia, 'file') && ~overwrite_sp_xia
                fprintf('Existing %s found. Not overwritten because overwrite_sp_xia = false.\n', out_xia);
            else
                save(out_xia, 'sp_clipped', '-v7.3');
                fprintf('Saved raw spikes as %s\n', out_xia);
            end
        end

    case 'xia'
        sp_files = dir(fullfile(data_folder, '*.sp_xia.mat'));
        sp_files = sp_files(~contains({sp_files.name}, 'PrefixRecovery'));
        sp_files = sp_files(~contains({sp_files.name}, 'SSD'));

        assert(~isempty(sp_files), 'No ordinary .sp_xia.mat file found.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading Xia spike file: %s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp_clipped')
            sp_clipped = S_sp.sp_clipped;
        else
            error('Variable "sp_clipped" not found in %s.', sp_filename);
        end

    case 'ssd'
        sp_files = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
        sp_files = sp_files(~contains({sp_files.name}, 'PrefixRecovery'));

        assert(~isempty(sp_files), 'No ordinary .sp_xia_SSD.mat file found.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading ordinary SSD spike file: %s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp_corr')
            sp_clipped = S_sp.sp_corr;
        else
            error('Variable "sp_corr" not found in %s.', sp_filename);
        end

    case 'recovery'
        sp_files = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
        sp_files = sp_files(~contains({sp_files.name}, 'SSD'));

        assert(~isempty(sp_files), 'No *sp_xia_PrefixRecovery.mat file found.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading recovered spike file: %s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp_seq')
            sp_clipped = S_sp.sp_seq;
        else
            error('Variable "sp_seq" not found in %s.', sp_filename);
        end

    case 'recovery_ssd'
        sp_files = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));

        assert(~isempty(sp_files), 'No *sp_xia_PrefixRecovery_SSD.mat file found.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading filtered recovered spike file: %s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp_corr')
            sp_clipped = S_sp.sp_corr;
        else
            error('Variable "sp_corr" not found in %s.', sp_filename);
        end

    case 'auto'
        sp_files_recovery_ssd = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));

        sp_files_recovery = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery.mat'));
        sp_files_recovery = sp_files_recovery(~contains({sp_files_recovery.name}, 'SSD'));

        sp_files_ssd = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
        sp_files_ssd = sp_files_ssd(~contains({sp_files_ssd.name}, 'PrefixRecovery'));

        sp_files_xia = dir(fullfile(data_folder, '*.sp_xia.mat'));
        sp_files_xia = sp_files_xia(~contains({sp_files_xia.name}, 'PrefixRecovery'));
        sp_files_xia = sp_files_xia(~contains({sp_files_xia.name}, 'SSD'));

        sp_files_raw = dir(fullfile(data_folder, '*.sp.mat'));
        sp_files_raw = sp_files_raw(~contains({sp_files_raw.name}, 'sp_xia'));

        if ~isempty(sp_files_recovery_ssd)
            sp_filename = fullfile(data_folder, sp_files_recovery_ssd(1).name);
            fprintf('Auto source: loading filtered recovered spike file: %s\n', sp_filename);
            S_sp = load(sp_filename);

            if isfield(S_sp, 'sp_corr')
                sp_clipped = S_sp.sp_corr;
            else
                error('Variable "sp_corr" not found in %s.', sp_filename);
            end

        elseif ~isempty(sp_files_recovery)
            sp_filename = fullfile(data_folder, sp_files_recovery(1).name);
            fprintf('Auto source: loading recovered spike file: %s\n', sp_filename);
            S_sp = load(sp_filename);

            if isfield(S_sp, 'sp_seq')
                sp_clipped = S_sp.sp_seq;
            else
                error('Variable "sp_seq" not found in %s.', sp_filename);
            end

        elseif ~isempty(sp_files_ssd)
            sp_filename = fullfile(data_folder, sp_files_ssd(1).name);
            fprintf('Auto source: loading ordinary SSD spike file: %s\n', sp_filename);
            S_sp = load(sp_filename);

            if isfield(S_sp, 'sp_corr')
                sp_clipped = S_sp.sp_corr;
            else
                error('Variable "sp_corr" not found in %s.', sp_filename);
            end

        elseif ~isempty(sp_files_xia)
            sp_filename = fullfile(data_folder, sp_files_xia(1).name);
            fprintf('Auto source: loading Xia spike file: %s\n', sp_filename);
            S_sp = load(sp_filename);

            if isfield(S_sp, 'sp_clipped')
                sp_clipped = S_sp.sp_clipped;
            else
                error('Variable "sp_clipped" not found in %s.', sp_filename);
            end

        elseif ~isempty(sp_files_raw)
            sp_filename = fullfile(data_folder, sp_files_raw(1).name);
            fprintf('Auto source: loading raw spike file: %s\n', sp_filename);
            S_sp = load(sp_filename);

            if isfield(S_sp, 'sp')
                sp = S_sp.sp;
                sp_clipped = sp;

                if save_sp_xia
                    out_xia = [base_name '.sp_xia.mat'];
                    save(out_xia, 'sp_clipped', '-v7.3');
                    fprintf('Saved raw spikes as %s\n', out_xia);
                end
            else
                error('Variable "sp" not found in %s.', sp_filename);
            end

        else
            error('No suitable spike file found for spike_source = auto.');
        end

    otherwise
        error('Unknown spike_source: %s. Use raw, xia, ssd, recovery, recovery_ssd, or auto.', spike_source);
end

fprintf('Spike source used for plotting: %s\n', spike_source);

%% ====================== LOAD TRIGGERS ======================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match trigger number (%d). Using min of both for plotting.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== REMOVE HEADER ROWS ======================

StimParams_data  = StimParams(2:end,:);
TrialParams_data = TrialParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

if size(TrialParams_data,1) ~= expected_rows
    warning('TrialParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(TrialParams_data,1), expected_rows);
end

%% ====================== TRIAL-LEVEL CONDITION IDS ======================
% TrialParams columns:
%   1 = trial number
%   2 = condition ID
%   3 = internal electrode ID
%
% We use the first row of each trial to get the trial-level condition ID.

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA FROM StimMeta ======================

stimSet_trial      = NaN(n_Trials,1);
trainLevel_trial   = NaN(n_Trials,1);
totalLevels_trial  = NaN(n_Trials,1);
amp_trial          = NaN(n_Trials,1);
isAutoSim_trial    = false(n_Trials,1);
isZero_trial       = false(n_Trials,1);
eventEnd_ms_trial  = NaN(n_Trials,1);

eventTimes_ms_trial = cell(n_Trials,1);
pulseCount_trial    = cell(n_Trials,1);

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
        totalLevels_trial(tr) = meta.TotalTrainLevels;
    end

    if isfield(meta, 'IsAutoSimultaneous')
        isAutoSim_trial(tr) = logical(meta.IsAutoSimultaneous);
    end

    if isfield(meta, 'IsZeroCurrentControl')
        isZero_trial(tr) = logical(meta.IsZeroCurrentControl);
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

    if isfield(meta, 'PulseCountPerElectrode')
        pulseCount_trial{tr} = meta.PulseCountPerElectrode;
    else
        pulseCount_trial{tr} = [];
    end

    % Use actual randomized StimParams rows for amplitude.
    % This is safer than relying only on StimMeta.
    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);

    if max(rr) <= size(StimParams_data,1)
        amp_vec = cell2mat(StimParams_data(rr,16));
        amp_vec = amp_vec(:)';

        if all(amp_vec <= 0)
            amp_trial(tr) = 0;
        else
            amp_trial(tr) = max(amp_vec(amp_vec > 0));
        end
    end
end

%% ====================== APPLY USER FILTERS ======================

% Amplitudes.
all_amps = unique(amp_trial(~isnan(amp_trial)));

if include_zero_control == 0
    all_amps = all_amps(all_amps > 0);
end

if isempty(Amps_to_plot)
    Amps_selected = all_amps;
else
    Amps_selected = intersect(all_amps, Amps_to_plot);
end

% Set IDs.
all_sets = unique(stimSet_trial(~isnan(stimSet_trial)));

if include_zero_control == 0
    all_sets = all_sets(all_sets > 0);
end

if isempty(SetIDs_to_plot)
    SetIDs_selected = all_sets;
else
    SetIDs_selected = intersect(all_sets, SetIDs_to_plot);
end

%% ====================== FAMILY-SPECIFIC TRAIN LEVEL FILTERS ======================
% Sequential and AutoSim can have different maximum train levels.
%
% Example:
%   Sequential final level = 6
%   AutoSim final level    = 3
%
% Therefore, we keep separate selected level lists.

all_seq_levels = unique(trainLevel_trial(~isnan(trainLevel_trial) & ...
                                         isAutoSim_trial == 0 & ...
                                         isZero_trial == 0));
all_seq_levels = sort(all_seq_levels(:))';

all_autosim_levels = unique(trainLevel_trial(~isnan(trainLevel_trial) & ...
                                             isAutoSim_trial == 1 & ...
                                             isZero_trial == 0));
all_autosim_levels = sort(all_autosim_levels(:))';

if isempty(seq_train_levels_to_plot)
    SeqLevels_selected = all_seq_levels;
else
    SeqLevels_selected = intersect(all_seq_levels, seq_train_levels_to_plot);
    SeqLevels_selected = sort(SeqLevels_selected(:))';
end

if isempty(autosim_train_levels_to_plot)
    AutoSimLevels_selected = all_autosim_levels;
else
    AutoSimLevels_selected = intersect(all_autosim_levels, autosim_train_levels_to_plot);
    AutoSimLevels_selected = sort(AutoSimLevels_selected(:))';
end

% This combined list is only used for debug printing.
TrainLevels_debug = unique([SeqLevels_selected, AutoSimLevels_selected]);

fprintf('\nDetected amplitudes: ');
disp(unique(amp_trial(~isnan(amp_trial)))');

fprintf('Selected amplitudes: ');
disp(Amps_selected');

fprintf('\nDetected StimSetIndex values: ');
disp(unique(stimSet_trial(~isnan(stimSet_trial)))');

fprintf('Selected StimSetIndex values: ');
disp(SetIDs_selected');

fprintf('\nDetected Seq TrainLevels: ');
disp(all_seq_levels);

fprintf('Selected Seq TrainLevels: ');
disp(SeqLevels_selected);

fprintf('\nDetected AutoSim TrainLevels: ');
disp(all_autosim_levels);

fprintf('Selected AutoSim TrainLevels: ');
disp(AutoSimLevels_selected);

fprintf('\nDetected event end times (ms): ');
disp(unique(eventEnd_ms_trial(~isnan(eventEnd_ms_trial)))');

fprintf('\nPlot sequential: %d\n', plot_sequential);
fprintf('Plot AutoSim:    %d\n', plot_auto_sim);

%% ====================== DEBUG TRIAL CONTENT CHECK ======================

if debug_print_trial_content

    fprintf('\n================ DEBUG PULSE-TRAIN TRIAL CONTENT CHECK ================\n');

    for autoFlag = [0 1]

        if autoFlag == 0
            debugFamily = 'Seq';
            DebugLevels_thisFamily = SeqLevels_selected;
        else
            debugFamily = 'AutoSim';
            DebugLevels_thisFamily = AutoSimLevels_selected;
        end

        if isempty(DebugLevels_thisFamily)
            fprintf('\n--- Debug family: %s has no selected levels. Skipped. ---\n', debugFamily);
            continue;
        end

        fprintf('\n--- Debug family: %s ---\n', debugFamily);

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            for aa = 1:numel(Amps_selected)

                amp_val = Amps_selected(aa);

                for li = 1:numel(DebugLevels_thisFamily)

                    level_val = DebugLevels_thisFamily(li);

                    tr_debug = find(stimSet_trial == setID & ...
                                    amp_trial == amp_val & ...
                                    trainLevel_trial == level_val & ...
                                    isAutoSim_trial == logical(autoFlag) & ...
                                    isZero_trial == 0, ...
                                    1, 'first');

                    if isempty(tr_debug)
                        continue;
                    end

                    rr = (tr_debug-1)*simultaneous_stim + (1:simultaneous_stim);

                    stimNames_debug   = StimParams_data(rr,1);
                    ptd_debug         = cell2mat(StimParams_data(rr,6));
                    pulseNum_debug    = cell2mat(StimParams_data(rr,8));
                    pulsePeriod_debug = cell2mat(StimParams_data(rr,9));
                    amp_debug         = cell2mat(StimParams_data(rr,16));

                    activeRows_debug = amp_debug > 0;
                    stimNames_active = stimNames_debug(activeRows_debug);
                    stimLabel_debug  = buildStimLabelFromStimNames(stimNames_active, E_MAP, logical(autoFlag));

                    fprintf('\n%s | Set %g | Amp %.1f uA | TrainLevel %g | Trial %d | CondID %d\n', ...
                        debugFamily, setID, amp_val, level_val, tr_debug, conditionID_trial(tr_debug));

                    fprintf('  IsAutoSimultaneous:        %d\n', isAutoSim_trial(tr_debug));
                    fprintf('  Stim label from E_MAP:     %s\n', stimLabel_debug);

                    fprintf('  Actual active StimParams CHANNEL: ');
                    disp(stimNames_active');

                    fprintf('  EventTimesIncluded_ms:     ');
                    disp(eventTimes_ms_trial{tr_debug});

                    fprintf('  EventEndTime_ms:           %.3f\n', eventEnd_ms_trial(tr_debug));

                    fprintf('  PulseCountPerElectrode:    ');
                    disp(pulseCount_trial{tr_debug});

                    fprintf('  StimParams PTD_us:         ');
                    disp(ptd_debug');

                    fprintf('  StimParams PulseNum:       ');
                    disp(pulseNum_debug');

                    fprintf('  StimParams PulsePeriod:    ');
                    disp(pulsePeriod_debug');

                    fprintf('  StimParams Amp_col16:      ');
                    disp(amp_debug');
                end
            end
        end
    end

    fprintf('=======================================================================\n');
end

%% ====================== ELECTRODE MAP ======================

d = Depth_s(Electrode_Type);

%% ====================== RASTER / PSTH SETUP ======================

edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

% Color by amplitude.
[Amps_all_for_color, ~] = unique(amp_trial(:));
cmap = lines(numel(Amps_all_for_color));

%% ====================== MAIN RASTER + PSTH LOOP ======================
% Plot structure:
%
%   Channel
%     Set
%       Family: Seq or AutoSim
%         TrainLevel
%           all selected amplitudes plotted together
%
% Therefore:
%   one figure = recording channel × set × family × train/event level
%
% This avoids mixing Seq and AutoSim trials.

for ich = raster_chn_start:raster_chn_end

    ch = d(ich);

    if ch > numel(sp_clipped)
        continue;
    end

    if isempty(sp_clipped{ch})
        continue;
    end

    for si = 1:numel(SetIDs_selected)

        setID = SetIDs_selected(si);

        %% ----- Build list of families to plot -----
        familyList = {};

        if plot_sequential
            familyList{end+1} = 'seq';
        end

        if plot_auto_sim
            familyList{end+1} = 'autosim';
        end

        for iFamily = 1:numel(familyList)

            familyName = familyList{iFamily};

            if strcmp(familyName, 'seq')
                autoFlag = false;
                familyTitle = 'Seq';
                Levels_this_family = SeqLevels_selected;

            elseif strcmp(familyName, 'autosim')
                autoFlag = true;
                familyTitle = 'AutoSim';
                Levels_this_family = AutoSimLevels_selected;

            else
                continue;
            end

            % If this family has no selected level, skip it.
            if isempty(Levels_this_family)
                continue;
            end

            for li = 1:numel(Levels_this_family)

                level_val = Levels_this_family(li);

                %% ----- Check whether this set/family/level has any selected amplitude -----
                hasAnyAmp = false;

                for aa = 1:numel(Amps_selected)

                    amp_val = Amps_selected(aa);

                    tlist_check = find(stimSet_trial == setID & ...
                                       amp_trial == amp_val & ...
                                       trainLevel_trial == level_val & ...
                                       isAutoSim_trial == autoFlag);

                    if include_zero_control == 0
                        tlist_check = tlist_check(isZero_trial(tlist_check) == 0);
                    end

                    tlist_check = tlist_check(tlist_check <= nTrials_use);

                    if ~isempty(tlist_check)
                        hasAnyAmp = true;
                        break;
                    end
                end

                if ~hasAnyAmp
                    continue;
                end

                %% ----- Find one example trial for event metadata -----
                tr_example = [];

                for aa = 1:numel(Amps_selected)

                    amp_val = Amps_selected(aa);

                    tlist_example = find(stimSet_trial == setID & ...
                                         amp_trial == amp_val & ...
                                         trainLevel_trial == level_val & ...
                                         isAutoSim_trial == autoFlag);

                    if include_zero_control == 0
                        tlist_example = tlist_example(isZero_trial(tlist_example) == 0);
                    end

                    tlist_example = tlist_example(tlist_example <= nTrials_use);

                    if ~isempty(tlist_example)
                        tr_example = tlist_example(1);
                        break;
                    end
                end

                if isempty(tr_example)
                    continue;
                end

                %% ----- Build set label across all levels from this family -----
                % Important:
                %   Level 1 in sequential mode may only activate one electrode.
                %   If we build the label using only Level 1, the title may
                %   show only one channel.
                %
                %   Therefore, collect stimulation channels across selected
                %   levels for the same set and family.

                setLabel = buildSetLabelAcrossLevels( ...
                    setID, autoFlag, Amps_selected, Levels_this_family, ...
                    stimSet_trial, amp_trial, trainLevel_trial, ...
                    isAutoSim_trial, isZero_trial, ...
                    StimParams_data, simultaneous_stim, E_MAP);

                %% ----- Event metadata from example trial -----
                eventTimes_ms = eventTimes_ms_trial{tr_example};
                eventEnd_ms   = eventEnd_ms_trial(tr_example);
                pulseCounts   = pulseCount_trial{tr_example};
                totalLevels   = totalLevels_trial(tr_example);

                if isempty(eventTimes_ms)
                    eventText = '[]';
                else
                    eventText = num2str(eventTimes_ms);
                end

                if isempty(pulseCounts)
                    pulseText = '[]';
                else
                    pulseText = num2str(pulseCounts);
                end

                %% ----- Create figure -----
                figName = sprintf('RasterPSTH_%s_%s_Ch%d_Set%d_Level%g_AllAmps', ...
                                  spike_source, familyTitle, ich, setID, level_val);

                figure('Color','w','Name',figName);
                tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                ax1 = nexttile([3 1]);
                hold(ax1,'on'); box(ax1,'off');

                title(ax1, sprintf(['%s | Rec Ch %d | Set%d: %s | ' ...
                                    'Level %.0f/%.0f | pulses [%s] | events [%s] ms'], ...
                                    familyTitle, ich, setID, setLabel, ...
                                    level_val, totalLevels, pulseText, eventText), ...
                                    'Interpreter','none');

                ax2 = nexttile;
                hold(ax2,'on'); box(ax2,'off');

                %% ----- Storage for raster amplitude grouping -----
                y_cursor = 0;
                ytick_vals = [];
                ytick_labels = {};

                maxRate = 0;
                legend_labels = {};

                %% ----- Loop through amplitudes inside the same figure -----
                for aa = 1:numel(Amps_selected)

                    amp_val = Amps_selected(aa);

                    %% ----- Find trials for this amplitude and family -----
                    amp_trials = find(stimSet_trial == setID & ...
                                      amp_trial == amp_val & ...
                                      trainLevel_trial == level_val & ...
                                      isAutoSim_trial == autoFlag);

                    if include_zero_control == 0
                        amp_trials = amp_trials(isZero_trial(amp_trials) == 0);
                    end

                    amp_trials = amp_trials(amp_trials <= nTrials_use);

                    nTr = numel(amp_trials);

                    if nTr == 0
                        continue;
                    end

                    %% ----- Choose colour for this amplitude -----
                    amp_color_idx = find(Amps_all_for_color == amp_val, 1, 'first');

                    if isempty(amp_color_idx)
                        color = [0 0 0];
                    else
                        color = cmap(amp_color_idx,:);
                    end

                    %% ----- Raster and PSTH counts -----
                    counts = zeros(1, numel(edges)-1);

                    for t = 1:nTr

                        tr = amp_trials(t);
                        t0 = trig(tr) / FS * 1000;

                        tt = sp_clipped{ch}(:,1);
                        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;

                        % Raster.
                        % Trials are stacked by amplitude.
                        y0 = y_cursor + t;

                        for spike_t = tt'
                            plot(ax1, [spike_t spike_t], [y0-0.4 y0+0.4], ...
                                 'Color', color, 'LineWidth', 1.1);
                        end

                        % PSTH count.
                        counts = counts + histcounts(tt, edges);
                    end

                    %% ----- PSTH for this amplitude -----
                    rate = filter(g, 1, counts/(nTr*bin_s));

                    plot(ax2, ctrs, rate, ...
                         'Color', color, ...
                         'LineWidth', 1.6);

                    maxRate = max(maxRate, max(rate));

                    %% ----- Raster y-axis group label -----
                    ytick_vals(end+1) = y_cursor + nTr/2; %#ok<SAGROW>
                    ytick_labels{end+1} = sprintf('%g uA', amp_val); %#ok<SAGROW>

                    y_cursor = y_cursor + nTr;

                    legend_labels{end+1} = sprintf('%g uA', amp_val); %#ok<SAGROW>
                end % amplitude loop

                %% ----- Event lines on raster -----
                drawEventLines(ax1, eventTimes_ms, eventEnd_ms);

                %% ----- Finalize raster axis -----
                xlim(ax1, ras_win);
                ylim(ax1, [0 max(1,y_cursor+1)]);
                yticks(ax1, ytick_vals);
                yticklabels(ax1, ytick_labels);
                ylabel(ax1, 'Amplitude');

                %% ----- Event lines on PSTH -----
                drawEventLines(ax2, eventTimes_ms, eventEnd_ms);

                %% ----- Finalize PSTH axis -----
                xlim(ax2, ras_win);

                if maxRate <= 0
                    ylim(ax2, [0 50]);
                else
                    ylim(ax2, [0 max(50, ceil(maxRate*1.1/10)*10)]);
                end

                xlabel(ax2, 'Time (ms)');
                ylabel(ax2, 'Rate (sp/s)');

                if ~isempty(legend_labels)
                    legend(ax2, legend_labels, 'Box','off','Location','northeast');
                end

            end % train level
        end % family
    end % set
end % channel

fprintf('\nFinished pulse-train raster + PSTH plotting.\n');

%% ====================== LOCAL HELPER FUNCTIONS ======================

function drawEventLines(ax, eventTimes_ms, eventEnd_ms)
    % Draw event timing lines:
    %   0 ms: red dashed
    %   other events: black dashed
    %   final event: blue dotted

    xline(ax, 0, 'r--', 'LineWidth', 1);

    if ~isempty(eventTimes_ms)

        for ee = 1:numel(eventTimes_ms)

            thisEvent = eventTimes_ms(ee);

            % Avoid drawing duplicate black line at 0 ms.
            if abs(thisEvent) < 1e-9
                continue;
            end

            xline(ax, thisEvent, 'k--', 'LineWidth', 1);
        end
    end

    if isfinite(eventEnd_ms) && eventEnd_ms > 0
        xline(ax, eventEnd_ms, 'b:', 'LineWidth', 1.5);
    end
end

function setLabel = buildSetLabelAcrossLevels( ...
    setID, autoFlag, Amps_selected, Levels_this_family, ...
    stimSet_trial, amp_trial, trainLevel_trial, ...
    isAutoSim_trial, isZero_trial, ...
    StimParams_data, simultaneous_stim, E_MAP)
    % Build a stimulation set label using selected levels.
    %
    % Why this is needed:
    %   In sequential event-level stimulation, Level 1 may activate only
    %   one electrode, e.g. [1 0].
    %
    %   If we build the label from Level 1 only, the title may incorrectly
    %   show only one channel.
    %
    %   Therefore, this function collects all active stimulation channel
    %   names across selected levels and amplitudes for the same set/family.
    %
    % Note:
    %   If you only plot the final level, e.g. Seq Level 6 or AutoSim Level 3,
    %   this still works because both channels are active in the final level.

    stimNames_all = {};

    for aa = 1:numel(Amps_selected)

        amp_val = Amps_selected(aa);

        for li = 1:numel(Levels_this_family)

            level_val = Levels_this_family(li);

            tr_label = find(stimSet_trial == setID & ...
                            amp_trial == amp_val & ...
                            trainLevel_trial == level_val & ...
                            isAutoSim_trial == autoFlag & ...
                            isZero_trial == 0, ...
                            1, 'first');

            if isempty(tr_label)
                continue;
            end

            rr_label = (tr_label-1)*simultaneous_stim + (1:simultaneous_stim);

            amp_label = cell2mat(StimParams_data(rr_label,16));
            activeRows = amp_label > 0;

            stimNames_this = StimParams_data(rr_label(activeRows),1);

            stimNames_all = [stimNames_all; stimNames_this(:)]; %#ok<AGROW>
        end
    end

    if isempty(stimNames_all)
        setLabel = sprintf('Set%d', setID);
        return;
    end

    % Remove duplicates while keeping first appearance order.
    stimNames_all = unique(stimNames_all, 'stable');

    setLabel = buildStimLabelFromStimNames(stimNames_all, E_MAP, autoFlag);
end

function setLabel = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim)
    % Build stimulation label for figure title.
    %
    % The label only shows mapped channel numbers.
    %
    % Example:
    %   Sequential:
    %       Ch22→Ch18
    %
    %   AutoSim:
    %       Ch22+Ch18

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
    % Convert Intan stimulation label such as 'A-017' to the corresponding
    % channel number using E_MAP.
    %
    % In your E_MAP format:
    %   E_MAP{1,1} is the array/map name.
    %   E_MAP{2,1} is channel 1.
    %   E_MAP{3,1} is channel 2.
    %   ...
    %
    % Therefore:
    %   if E_MAP{row,1} matches stimName,
    %   channel number = row - 1.

    chNum = NaN;

    if isempty(stimName)
        return;
    end

    if isstring(stimName)
        stimName = char(stimName);
    end

    stimName = strtrim(stimName);

    %% ----- Main case: E_MAP is a cell array -----
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
    end

    %% ----- If E_MAP is a string array -----
    if isstring(E_MAP)

        hit = find(strcmp(strtrim(E_MAP), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    %% ----- If E_MAP is a char matrix -----
    if ischar(E_MAP)

        nameList = cellstr(E_MAP);
        hit = find(strcmp(strtrim(nameList), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    %% ----- Fallback only if E_MAP matching fails -----
    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
        warning('Stim channel %s was not found in E_MAP. Falling back to parsed number %d.', ...
            stimName, chNum);
    else
        warning('Stim channel %s was not found in E_MAP and could not be parsed.', stimName);
    end
end