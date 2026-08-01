%% Raster + PSTH Viewer for Pulse-Train / Event-Level Stimulation Files
% ------------------------------------------------------------
% Purpose:
%   1) Load spike file.
%   2) Optionally save raw *.sp.mat as *.sp_xia.mat with variable sp_clipped.
%   3) Plot raster + PSTH grouped by:
%        recording channel
%        stimulation set
%        amplitude
%        train/event level
% ------------------------------------------------------------

clear all
% close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== CHOOSE FOLDER ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX030/Xia_Elec2_Train1';

%% ====================== USER SETTINGS ======================

raster_chn_start = 1;
raster_chn_end   = 64;   % Depth_s index
Electrode_Type   = 2;    % 0: rigid; 1: single-shank flex; 2: four-shank flex

% Stimulation set selection.
% [] means all detected sets.
SetIDs_to_plot = [1];

% Amplitudes to plot.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% Train/event levels to plot.
% [] means all detected train levels.
% Example:
%   [1 3 6] plots only level 1, level 3, and level 6.
train_levels_to_plot = [];

% Include automatically added simultaneous condition?
% 1 = include, 0 = exclude.
include_auto_sim = 1;

% Include zero-current control trials?
% 1 = include, 0 = exclude.
include_zero_control = 0;

% Raster/PSTH parameters.
ras_win       = [-20 100];   % ms
bin_ms_raster = 1;           % PSTH bin size, ms
smooth_ms     = 2;           % PSTH smoothing width

% Spike file source:
%   'raw'  = load *.sp.mat, variable sp
%   'xia'  = load *.sp_xia.mat, variable sp_clipped
%   'ssd'  = load *.sp_xia_SSD.mat, variable sp_corr
%   'auto' = priority: SSD > xia > raw
spike_source = 'raw';

% If true, save raw *.sp.mat as *.sp_xia.mat with variable sp_clipped.
% This is kept from your original code because it is important for your pipeline.
save_sp_xia = true;

% If true, overwrite existing *.sp_xia.mat.
overwrite_sp_xia = true;

% If true, print one example trial for each selected condition.
debug_print_trial_content = false;

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
% This section keeps your original raw .sp.mat -> .sp_xia.mat saving.
% It also adds a choice of spike source for later checking.

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

        % Keep original saving step.
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
        assert(~isempty(sp_files), 'No .sp_xia.mat file found.');

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
        assert(~isempty(sp_files), 'No .sp_xia_SSD.mat file found.');

        sp_filename = fullfile(data_folder, sp_files(1).name);
        fprintf('Loading SSD spike file: %s\n', sp_filename);

        S_sp = load(sp_filename);

        if isfield(S_sp, 'sp_corr')
            sp_clipped = S_sp.sp_corr;
        else
            error('Variable "sp_corr" not found in %s.', sp_filename);
        end

    case 'auto'
        sp_files_ssd = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
        sp_files_xia = dir(fullfile(data_folder, '*.sp_xia.mat'));
        sp_files_raw = dir(fullfile(data_folder, '*.sp.mat'));
        sp_files_raw = sp_files_raw(~contains({sp_files_raw.name}, 'sp_xia'));

        if ~isempty(sp_files_ssd)
            sp_filename = fullfile(data_folder, sp_files_ssd(1).name);
            fprintf('Auto source: loading SSD spike file: %s\n', sp_filename);
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
        error('Unknown spike_source: %s. Use raw, xia, ssd, or auto.', spike_source);
end

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

    if isfield(meta, 'Amp_uA')
        amp_vec = meta.Amp_uA;
        amp_vec = amp_vec(:)';

        if all(amp_vec <= 0)
            amp_trial(tr) = 0;
        else
            amp_trial(tr) = max(amp_vec(amp_vec > 0));
        end
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
    end

    if isfield(meta, 'PulseCountPerElectrode')
        pulseCount_trial{tr} = meta.PulseCountPerElectrode;
    else
        pulseCount_trial{tr} = [];
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

% Train levels.
all_train_levels = unique(trainLevel_trial(~isnan(trainLevel_trial)));

if isempty(train_levels_to_plot)
    TrainLevels_selected = all_train_levels;
else
    TrainLevels_selected = intersect(all_train_levels, train_levels_to_plot);
end

fprintf('\nDetected amplitudes: ');
disp(unique(amp_trial(~isnan(amp_trial)))');

fprintf('Selected amplitudes: ');
disp(Amps_selected');

fprintf('\nDetected StimSetIndex values: ');
disp(unique(stimSet_trial(~isnan(stimSet_trial)))');

fprintf('Selected StimSetIndex values: ');
disp(SetIDs_selected');

fprintf('\nDetected TrainLevels: ');
disp(all_train_levels');

fprintf('Selected TrainLevels: ');
disp(TrainLevels_selected');

fprintf('\nDetected event end times (ms): ');
disp(unique(eventEnd_ms_trial(~isnan(eventEnd_ms_trial)))');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================

if debug_print_trial_content

    fprintf('\n================ DEBUG PULSE-TRAIN TRIAL CONTENT CHECK ================\n');

    for si = 1:numel(SetIDs_selected)

        setID = SetIDs_selected(si);

        for aa = 1:numel(Amps_selected)

            amp_val = Amps_selected(aa);

            for li = 1:numel(TrainLevels_selected)

                level_val = TrainLevels_selected(li);

                tr_debug = find(stimSet_trial == setID & ...
                                amp_trial == amp_val & ...
                                trainLevel_trial == level_val & ...
                                isZero_trial == 0, ...
                                1, 'first');

                if include_auto_sim == 0 && ~isempty(tr_debug) && isAutoSim_trial(tr_debug)
                    tr_debug = [];
                end

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
                stimLabel_debug  = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim_trial(tr_debug));

                fprintf('\nSet %g | Amp %.1f uA | TrainLevel %g | Trial %d | CondID %d\n', ...
                    setID, amp_val, level_val, tr_debug, conditionID_trial(tr_debug));

                fprintf('  IsAutoSimultaneous:        %d\n', isAutoSim_trial(tr_debug));
                fprintf('  Stim label from E_MAP:     %s\n', stimLabel_debug);

                fprintf('  Actual active StimParams CHANNEL: ');
                disp(stimNames_active');

                fprintf('  EventTimesIncluded_ms:     ');
                disp(eventTimes_ms_trial{tr_debug});

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

                % Recompute event times from active StimParams rows.
                eventTimes_fromRows_us = [];

                for rridx = find(activeRows_debug(:))'
                    eventTimes_fromRows_us = [eventTimes_fromRows_us, ...
                        ptd_debug(rridx) + (0:pulseNum_debug(rridx)-1) * pulsePeriod_debug(rridx)];
                end

                eventTimes_fromRows_ms = unique(sort(eventTimes_fromRows_us)) ./ 1000;

                fprintf('  Recomputed events from active rows (ms): ');
                disp(eventTimes_fromRows_ms);

            end
        end
    end

    fprintf('=======================================================================\n');
end

%% ====================== PULSE TRAIN PERIOD ======================
% Kept for compatibility/debugging. Not used as a main grouping variable.

pulseTrain_all = cell2mat(StimParams_data(:,9));
pulseTrain = pulseTrain_all(firstRow_eachTrial);
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:)); %#ok<ASGLU>
n_PULSE = numel(PulsePeriods); %#ok<NASGU>

%% ====================== ELECTRODE MAP ======================

d = Depth_s(Electrode_Type);

%% ====================== RASTER / PSTH SETUP ======================

edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

% Color by amplitude.
[Amps_all_for_color, ~, ampIdx_all] = unique(amp_trial(:));
cmap = lines(numel(Amps_all_for_color));

%% ====================== MAIN RASTER + PSTH LOOP ======================
% New pulse-train plot structure:
%
%   Channel
%     Set
%       TrainLevel
%         all selected amplitudes plotted together



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

        for li = 1:numel(TrainLevels_selected)

            level_val = TrainLevels_selected(li);

            %% ----- Check whether this set/level has any selected amplitude -----
            hasAnyAmp = false;

            for aa = 1:numel(Amps_selected)

                amp_val = Amps_selected(aa);

                tlist_check = find(stimSet_trial == setID & ...
                                   amp_trial == amp_val & ...
                                   trainLevel_trial == level_val);

                if include_zero_control == 0
                    tlist_check = tlist_check(isZero_trial(tlist_check) == 0);
                end

                if include_auto_sim == 0
                    tlist_check = tlist_check(isAutoSim_trial(tlist_check) == 0);
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

            %% ----- Find one example trial for title/event metadata -----
            % Prefer a non-zero stimulation trial if possible.
            tr_example = [];

            for aa = 1:numel(Amps_selected)

                amp_val = Amps_selected(aa);

                tlist_example = find(stimSet_trial == setID & ...
                                     amp_trial == amp_val & ...
                                     trainLevel_trial == level_val);

                if include_zero_control == 0
                    tlist_example = tlist_example(isZero_trial(tlist_example) == 0);
                end

                if include_auto_sim == 0
                    tlist_example = tlist_example(isAutoSim_trial(tlist_example) == 0);
                end

                tlist_example = tlist_example(tlist_example <= nTrials_use);

                % Prefer active stimulation trials for the set label.
                tlist_nonzero = tlist_example(isZero_trial(tlist_example) == 0);

                if ~isempty(tlist_nonzero)
                    tr_example = tlist_nonzero(1);
                    break;
                elseif isempty(tr_example) && ~isempty(tlist_example)
                    tr_example = tlist_example(1);
                end
            end

            if isempty(tr_example)
                continue;
            end

            %% ----- Build display metadata from the example trial -----
            rr_label = (tr_example-1)*simultaneous_stim + (1:simultaneous_stim);

            stimNames_label = StimParams_data(rr_label,1);
            amp_label       = cell2mat(StimParams_data(rr_label,16));
            activeRows      = amp_label > 0;

            stimNames_active = stimNames_label(activeRows);

            % This helper uses StimParams channel names, e.g. A-017,
            % and converts them through E_MAP for the title label.
            setLabel = buildStimLabelFromStimNames( ...
                stimNames_active, E_MAP, isAutoSim_trial(tr_example));

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

            if isAutoSim_trial(tr_example)
                stimTypeText = 'AutoSim';
            elseif isZero_trial(tr_example)
                stimTypeText = 'ZeroControl';
            else
                stimTypeText = 'PulseTrain';
            end

            %% ----- Create figure -----
            figName = sprintf('RasterPSTH_Ch%d_Set%d_Level%g_AllAmps', ...
                              ich, setID, level_val);

            figure('Color','w','Name',figName);
            tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

            ax1 = nexttile([3 1]);
            hold(ax1,'on'); box(ax1,'off');

            title(ax1, sprintf(['Rec Ch %d | Set%d: %s | %s | ' ...
                                'Level %.0f/%.0f | pulses [%s] | events [%s] ms'], ...
                                ich, setID, setLabel, stimTypeText, ...
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

                %% ----- Find trials for this amplitude -----
                amp_trials = find(stimSet_trial == setID & ...
                                  amp_trial == amp_val & ...
                                  trainLevel_trial == level_val);

                if include_zero_control == 0
                    amp_trials = amp_trials(isZero_trial(amp_trials) == 0);
                end

                if include_auto_sim == 0
                    amp_trials = amp_trials(isAutoSim_trial(amp_trials) == 0);
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
                ytick_vals(end+1) = y_cursor + nTr/2;
                ytick_labels{end+1} = sprintf('%g uA', amp_val);

                y_cursor = y_cursor + nTr;

                legend_labels{end+1} = sprintf('%g uA', amp_val);
            end % amplitude loop

            %% ----- Event lines on raster -----
            xline(ax1, 0, 'r--');

            if ~isempty(eventTimes_ms)
                for ee = 1:numel(eventTimes_ms)

                    thisEvent = eventTimes_ms(ee);

                    % Avoid drawing duplicate black line at 0 ms.
                    if abs(thisEvent) < 1e-9
                        continue;
                    end

                    xline(ax1, thisEvent, 'k--', 'LineWidth', 1);
                end
            end

            % Final event line.
            % This is useful for seeing where the last stimulation event occurred.
            if isfinite(eventEnd_ms) && eventEnd_ms > 0
                xline(ax1, eventEnd_ms, 'b:', 'LineWidth', 1.5);
            end

            %% ----- Finalize raster axis -----
            xlim(ax1, ras_win);
            ylim(ax1, [0 max(1,y_cursor+1)]);
            yticks(ax1, ytick_vals);
            yticklabels(ax1, ytick_labels);
            ylabel(ax1, 'Amplitude');

            %% ----- Event lines on PSTH -----
            xline(ax2, 0, 'r--');

            if ~isempty(eventTimes_ms)
                for ee = 1:numel(eventTimes_ms)

                    thisEvent = eventTimes_ms(ee);

                    if abs(thisEvent) < 1e-9
                        continue;
                    end

                    xline(ax2, thisEvent, 'k--', 'LineWidth', 1);
                end
            end

            if isfinite(eventEnd_ms) && eventEnd_ms > 0
                xline(ax2, eventEnd_ms, 'b:', 'LineWidth', 1.5);
            end

            %% ----- Finalize PSTH axis -----
            xlim(ax2, ras_win);
            ylim(ax2, [0 max(50, ceil(maxRate*1.1/10)*10)]);
            xlabel(ax2, 'Time (ms)');
            ylabel(ax2, 'Rate (sp/s)');

            if ~isempty(legend_labels)
                legend(ax2, legend_labels, 'Box','off','Location','northeast');
            end

        end % train level
    end % set
end % channel

fprintf('\nFinished pulse-train raster + PSTH plotting.\n');

%% ====================== LOCAL HELPER FUNCTIONS ======================
% These helper functions are only used for title labels.
% They do not change the data analysis.

function setLabel = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim)

    if isempty(stimNames_active)
        setLabel = 'NoActiveStim';
        return;
    end

    % Remove duplicates while keeping order.
    stimNames_active = unique(stimNames_active, 'stable');

    labelParts = cell(1, numel(stimNames_active));

    for i = 1:numel(stimNames_active)

        stimName = stimNames_active{i};

        chNum = convertStimNameUsingEMap(stimName, E_MAP);

        if isnan(chNum)
            % If conversion fails, keep original Intan name.
            labelParts{i} = sprintf('%s', stimName);
        else
            % Show both channel number and original Intan name.
            labelParts{i} = sprintf('Ch%d(%s)', chNum, stimName);
        end
    end

    if isAutoSim
        setLabel = ['AutoSim: ' strjoin(labelParts, '+')];
    else
        setLabel = strjoin(labelParts, '→');
    end
end

function chNum = convertStimNameUsingEMap(stimName, E_MAP)
    % Convert Intan stimulation label such as 'A-017' to the corresponding
    % channel number using E_MAP.
    %
    % In your E_MAP format:
    %   E_MAP is usually a 129x1 cell array.
    %   E_MAP{1,1} is the array/map name, e.g. '4SHANK ...'
    %   E_MAP{2,1} is channel 1
    %   E_MAP{3,1} is channel 2
    %   ...
    %
    % Therefore:
    %   if E_MAP{row,1} matches stimName,
    %   then channel number = row - 1.
    %
    % Example:
    %   E_MAP{19,1} = 'A-016'
    %   channel number = 19 - 1 = 18

    chNum = NaN;

    if isempty(stimName)
        return;
    end

    % Make sure stimName is char.
    if isstring(stimName)
        stimName = char(stimName);
    end

    % Remove spaces just in case.
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

                    % Row 1 is header, so channel number = row - 1.
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
    % This fallback is not the preferred mapping.
    % It only prevents the code from crashing if a channel name cannot be
    % found in E_MAP.
    %
    % Example:
    %   'A-017' -> 17
    %
    % But this may not match your actual electrode/channel map.
    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
        warning('Stim channel %s was not found in E_MAP. Falling back to parsed number %d.', ...
            stimName, chNum);
    else
        warning('Stim channel %s was not found in E_MAP and could not be parsed.', stimName);
    end
end
         