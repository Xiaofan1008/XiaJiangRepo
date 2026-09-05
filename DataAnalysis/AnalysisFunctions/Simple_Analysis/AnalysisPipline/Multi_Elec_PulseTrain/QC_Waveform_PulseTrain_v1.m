%% ============================================================
%  Waveform Viewer for SSD-Filtered Pulse-Train Recovery Files
%
%  Purpose:
%    Plot spike waveforms by recording channel so that noisy / bad channels
%    can be identified after pulse-train spike recovery and SSD filtering.
%
%  Input spike file:
%    *.sp_xia_PrefixRecovery_SSD.mat
%
%  Required spike variable:
%    sp_corr
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train3';

% Recording channels to plot, using Depth_s index.
spike_chn_start = 1;
spike_chn_end   = 64;

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Stimulation set IDs to plot.
% [] means all detected sets.
SetIDs_to_plot = [2];

% Amplitudes to plot.
% [] means all non-zero amplitudes.
Amps_to_plot = []; 

% Sequential/interleaved train levels to plot.
% [] means all detected sequential levels.
seq_train_levels_to_plot = 6;

% AutoSim/simultaneous train levels to plot.
% [] means all detected AutoSim levels.
autosim_train_levels_to_plot = 3;

% Plot sequential / interleaved family?
plot_sequential = true;

% Plot AutoSim / simultaneous family?
plot_auto_sim = true;

% Include zero-current control trials?
% Usually false for checking recovered stimulation responses.
include_zero_control = false;

% Time window after trigger used to collect waveforms.
waveform_time_window_ms = [-10 40];

% Bin width for grouping waveforms by spike latency.
% [0 40] with 2 ms bins gives 20 bins, matching 4 x 5 layout.
bin_ms = 2;

% Optional waveform amplitude threshold for plotting.
% Waveforms with any point larger than this absolute value are not plotted.
% Set Inf to disable.
amp_threshold = 500;

% If true, align each waveform to its trough before plotting.
% This helps inspect waveform shape independent of small timing jitter.
align_waveforms_to_trough = true;

% Plot individual waveforms and/or mean waveform.
plot_individual_waveforms = true;
plot_mean_waveform = true;

% To avoid very slow/crowded plots, randomly limit the number of waveforms
% per bin and amplitude.
max_waveforms_per_bin_amp = 200;

% Figure layout.
layout_row = [];
layout_col = [];

% Sampling frequency.
FS = 30000;

% Print example trial information.
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

%% ====================== LOAD FILTERED RECOVERED SPIKES ======================
% For the current pulse-train pipeline, we specifically load:
%
%   *.sp_xia_PrefixRecovery_SSD.mat
%
% This avoids accidentally loading the ordinary *.sp_xia_SSD.mat file.

ssd_files = dir(fullfile(data_folder, '*sp_xia_PrefixRecovery_SSD.mat'));
assert(~isempty(ssd_files), 'No *sp_xia_PrefixRecovery_SSD.mat file found in the current folder.');

ssd_filename = fullfile(data_folder, ssd_files(1).name);
fprintf('Loading filtered recovered spike file:\n%s\n', ssd_filename);

S_sp = load(ssd_filename);

if isfield(S_sp, 'sp_corr')
    sp_use = S_sp.sp_corr;
else
    error('Variable "sp_corr" not found in %s.', ssd_filename);
end

nCh = numel(sp_use);

fprintf('Loaded filtered recovered spikes from: %s\n', ssd_files(1).name);
fprintf('Number of spike channels: %d\n', nCh);

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

%% ====================== LOAD EXPERIMENT PARAMETERS ======================

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
    error('StimMeta was not found. This pulse-train waveform viewer requires StimMeta.');
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

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
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
% Each trial has simultaneous_stim rows.
% We use the first row of each trial to get the condition ID.

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA FROM StimMeta ======================
% This replaces the old mixed-prefix metadata columns 26–30.
%
% Old variables:
%   conditionSetID_trial -> stimSet_trial
%   prefixLength_trial   -> trainLevel_trial
%   conditionType_trial  -> isAutoSim_trial / isZero_trial
%   lastActivePTD_ms     -> eventEnd_ms_trial
%   isi_ms_trial         -> removed

stimSet_trial      = NaN(n_Trials,1);
trainLevel_trial   = NaN(n_Trials,1);
totalLevels_trial  = NaN(n_Trials,1);
trialAmps          = NaN(n_Trials,1);
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
            trialAmps(tr) = 0;
        else
            trialAmps(tr) = max(amp_vec(amp_vec > 0));
        end
    end
end

% Compatibility name.
% In pulse-train data, lastActivePTD_ms is simply the final event time.
lastActivePTD_ms = eventEnd_ms_trial;

%% ====================== APPLY USER FILTERS ======================

% Set IDs.
SetIDs_all = unique(stimSet_trial(~isnan(stimSet_trial)));

if include_zero_control == 0
    SetIDs_all = SetIDs_all(SetIDs_all > 0);
end

if isempty(SetIDs_to_plot)
    SetIDs_selected = SetIDs_all;
else
    SetIDs_selected = intersect(SetIDs_all, SetIDs_to_plot);
end

% Amplitudes.
Amps_all = unique(trialAmps(~isnan(trialAmps)));

if include_zero_control == 0
    Amps_all = Amps_all(Amps_all > 0);
end

if isempty(Amps_to_plot)
    Amps_selected = Amps_all;
else
    Amps_selected = intersect(Amps_all, Amps_to_plot);
end

% Family-specific train levels.
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

% Color by amplitude.
[Amps_for_color, ~] = unique(trialAmps(:));
cmap = lines(numel(Amps_for_color));

fprintf('\nDetected amplitudes: ');
disp(unique(trialAmps(~isnan(trialAmps)))');

fprintf('Selected amplitudes: ');
disp(Amps_selected');

fprintf('\nDetected set IDs: ');
disp(unique(stimSet_trial(~isnan(stimSet_trial)))');

fprintf('Selected set IDs: ');
disp(SetIDs_selected');

fprintf('\nDetected Seq train levels: ');
disp(all_seq_levels);

fprintf('Selected Seq train levels: ');
disp(SeqLevels_selected);

fprintf('\nDetected AutoSim train levels: ');
disp(all_autosim_levels);

fprintf('Selected AutoSim train levels: ');
disp(AutoSimLevels_selected);

fprintf('\nDetected final event times (ms): ');
disp(unique(eventEnd_ms_trial(~isnan(eventEnd_ms_trial)))');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================

if debug_print_trial_content

    fprintf('\n================ DEBUG PULSE-TRAIN WAVEFORM CHECK ================\n');

    for autoFlag = [0 1]

        if autoFlag == 0
            debugFamily = 'Seq';
            DebugLevels_thisFamily = SeqLevels_selected;
        else
            debugFamily = 'AutoSim';
            DebugLevels_thisFamily = AutoSimLevels_selected;
        end

        if isempty(DebugLevels_thisFamily)
            fprintf('\n--- %s has no selected levels. Skipped. ---\n', debugFamily);
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
                                    trialAmps == amp_val & ...
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

                    fprintf('\n%s | Set %g | Amp %.1f uA | Level %g | Trial %d | CondID %d\n', ...
                        debugFamily, setID, amp_val, level_val, tr_debug, conditionID_trial(tr_debug));

                    fprintf('  Stim label:              %s\n', stimLabel_debug);
                    fprintf('  IsAutoSimultaneous:      %d\n', isAutoSim_trial(tr_debug));

                    fprintf('  EventTimesIncluded_ms:   ');
                    disp(eventTimes_ms_trial{tr_debug});

                    fprintf('  EventEndTime_ms:         %.3f\n', eventEnd_ms_trial(tr_debug));

                    fprintf('  PulseCountPerElectrode:  ');
                    disp(pulseCount_trial{tr_debug});

                    fprintf('  StimParams CHANNEL:      ');
                    disp(stimNames_active');

                    fprintf('  StimParams PTD_us:       ');
                    disp(ptd_debug');

                    fprintf('  StimParams PulseNum:     ');
                    disp(pulseNum_debug');

                    fprintf('  StimParams Period_us:    ');
                    disp(pulsePeriod_debug');

                    fprintf('  StimParams Amp_col16:    ');
                    disp(amp_debug');
                end
            end
        end
    end

    fprintf('===================================================================\n');
end

%% ====================== ELECTRODE MAP ======================

d = Depth_s(Electrode_Type);

%% ====================== WAVEFORM TIME BINS ======================

bin_edges = waveform_time_window_ms(1):bin_ms:waveform_time_window_ms(2);
nBins = numel(bin_edges) - 1;

% Automatically choose subplot layout if not manually specified.
if isempty(layout_row) || isempty(layout_col)

    % Aim for a reasonably wide figure.
    layout_col = ceil(sqrt(nBins));
    layout_row = ceil(nBins / layout_col);

    fprintf('Auto subplot layout: %d rows x %d columns for %d bins.\n', ...
        layout_row, layout_col, nBins);

else
    layout_capacity = layout_row * layout_col;

    if nBins > layout_capacity
        warning(['Number of time bins (%d) is larger than subplot layout capacity (%d). ', ...
                 'Automatically increasing layout size.'], ...
                 nBins, layout_capacity);

        layout_col = ceil(sqrt(nBins));
        layout_row = ceil(nBins / layout_col);

        fprintf('Updated subplot layout: %d rows x %d columns for %d bins.\n', ...
            layout_row, layout_col, nBins);
    end
end
%% ====================== WAVEFORM X AXIS ======================

example_ch = find(~cellfun(@isempty, sp_use), 1, 'first');

if isempty(example_ch)
    error('All channels are empty in sp_use.');
end

wf_len = size(sp_use{example_ch}, 2) - 1;
t_wave = (0:wf_len-1) / FS * 1000;

%% ====================== MAIN WAVEFORM PLOTTING LOOP ======================
% Figure structure:
%
%   one figure = recording channel × set × family × train level
%
% Each subplot = spike latency bin.
% Within each subplot, waveform colors represent amplitudes.

for ich = spike_chn_start:spike_chn_end

    ch = d(ich);

    if ch > nCh
        continue;
    end

    if isempty(sp_use{ch})
        continue;
    end

    sp_times_all = sp_use{ch}(:,1);
    sp_wave_all  = sp_use{ch}(:,2:end);

    % Optional amplitude threshold for plotting.
    % This only controls plotting, not the saved spike file.
    valid_idx = all(abs(sp_wave_all) <= amp_threshold, 2);

    sp_times_all = sp_times_all(valid_idx);
    sp_wave_all  = sp_wave_all(valid_idx,:);

    if isempty(sp_times_all)
        continue;
    end

    for si = 1:numel(SetIDs_selected)

        setID = SetIDs_selected(si);

        %% ----- Build family list -----
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

            if isempty(Levels_this_family)
                continue;
            end

            for li = 1:numel(Levels_this_family)

                level_val = Levels_this_family(li);

                %% ----- Trial group -----
                group_trials = find(stimSet_trial == setID & ...
                                    trainLevel_trial == level_val & ...
                                    isAutoSim_trial == autoFlag & ...
                                    ismember(trialAmps, Amps_selected));

                if include_zero_control == 0
                    group_trials = group_trials(isZero_trial(group_trials) == 0);
                end

                group_trials = group_trials(group_trials <= nTrials_use);

                if isempty(group_trials)
                    continue;
                end

                %% ----- Build set label across selected levels -----
                setLabel = buildSetLabelAcrossLevels( ...
                    setID, autoFlag, Amps_selected, Levels_this_family, ...
                    stimSet_trial, trialAmps, trainLevel_trial, ...
                    isAutoSim_trial, isZero_trial, ...
                    StimParams_data, simultaneous_stim, E_MAP);

                %% ----- Representative metadata -----
                tr_example = group_trials(1);

                eventTimes_this = eventTimes_ms_trial{tr_example};
                eventEnd_this   = eventEnd_ms_trial(tr_example);
                pulseCounts     = pulseCount_trial{tr_example};
                totalLevels     = totalLevels_trial(tr_example);

                if isempty(eventTimes_this)
                    eventText = '[]';
                else
                    eventText = num2str(eventTimes_this);
                end

                if isempty(pulseCounts)
                    pulseText = '[]';
                else
                    pulseText = num2str(pulseCounts);
                end

                finalEvent_ms = max(eventEnd_ms_trial(group_trials));

                %% ----- Collect waveforms by latency bin and amplitude -----
                all_spikes_by_bin_amp = cell(nBins, numel(Amps_selected));

                for tt = 1:numel(group_trials)

                    tr = group_trials(tt);
                    t0_ms = trig(tr) / FS * 1000;

                    rel_times_all = sp_times_all - t0_ms;

                    in_win = rel_times_all >= waveform_time_window_ms(1) & ...
                             rel_times_all <  waveform_time_window_ms(2);

                    if ~any(in_win)
                        continue;
                    end

                    rel_times = rel_times_all(in_win);
                    waveforms = sp_wave_all(in_win,:);

                    amp_val_trial = trialAmps(tr);
                    amp_pos = find(Amps_selected == amp_val_trial, 1, 'first');

                    if isempty(amp_pos)
                        continue;
                    end

                    for j = 1:numel(rel_times)

                        bin_idx = find(rel_times(j) >= bin_edges(1:end-1) & ...
                                       rel_times(j) <  bin_edges(2:end), ...
                                       1, 'first');

                        if isempty(bin_idx)
                            continue;
                        end

                        all_spikes_by_bin_amp{bin_idx, amp_pos}(end+1,:) = waveforms(j,:);
                    end
                end

                %% ----- Skip if no waveforms -----
                all_waves = cell2mat(all_spikes_by_bin_amp(:));

                if isempty(all_waves)
                    continue;
                end

                %% ----- Y limit -----
                y_max = max(abs(all_waves(:)));

                if isempty(y_max) || y_max == 0 || isnan(y_max)
                    y_lim = [-100 100];
                else
                    y_lim = [-1 1] * ceil(y_max/50)*50;
                end

                %% ----- Figure -----
                figName = sprintf('SSDRecoveredWaveforms_%s_Ch%d_Set%d_Level%g', ...
                    familyTitle, ich, setID, level_val);

                figure('Name', figName, ...
                       'Color', 'w', ...
                       'Position', [100 100 1400 800]);

                tiledlayout(layout_row, layout_col, ...
                            'Padding', 'compact', ...
                            'TileSpacing', 'compact');

                sgtitle(sprintf(['SSD recovered waveforms | Rec Ch %d | %s | Set %d: %s | ' ...
                                 'Level %.0f/%.0f | pulses [%s] | events [%s] ms | Final event %.1f ms'], ...
                    ich, familyTitle, setID, setLabel, ...
                    level_val, totalLevels, pulseText, eventText, finalEvent_ms), ...
                    'Interpreter', 'none');

                %% ----- Plot bins -----
                for b = 1:nBins

                    ax = nexttile;
                    hold(ax, 'on');
                    box(ax, 'off');

                    spike_count = 0;

                    bin_start = bin_edges(b);
                    bin_end   = bin_edges(b+1);

                    for aa = 1:numel(Amps_selected)

                        waves = all_spikes_by_bin_amp{b,aa};

                        if isempty(waves)
                            continue;
                        end

                        % Randomly limit waveforms if too many.
                        if size(waves,1) > max_waveforms_per_bin_amp
                            rand_idx = randperm(size(waves,1), max_waveforms_per_bin_amp);
                            waves_plot = waves(rand_idx,:);
                        else
                            waves_plot = waves;
                        end

                        spike_count = spike_count + size(waves,1);

                        %% ----- Optional trough alignment -----
                        if align_waveforms_to_trough

                            aligned_waves = zeros(size(waves_plot));

                            for k = 1:size(waves_plot,1)
                                [~, min_idx] = min(waves_plot(k,:));
                                shift = ceil(size(waves_plot,2)/2) - min_idx;
                                aligned_waves(k,:) = circshift(waves_plot(k,:), shift, 2);
                            end

                            waves_to_plot = aligned_waves;

                        else
                            waves_to_plot = waves_plot;
                        end

                        %% ----- Color for this amplitude -----
                        amp_global_idx = find(Amps_for_color == Amps_selected(aa), 1, 'first');

                        if isempty(amp_global_idx)
                            amp_color = [0 0 0];
                        else
                            amp_color = cmap(amp_global_idx,:);
                        end

                        %% ----- Plot individual waveforms -----
                        if plot_individual_waveforms
                            plot(ax, t_wave, waves_to_plot', ...
                                 'Color', [amp_color 0.25]);
                        end

                        %% ----- Plot mean waveform -----
                        if plot_mean_waveform
                            mean_wf = mean(waves_to_plot, 1);
                            plot(ax, t_wave, mean_wf, ...
                                 'Color', amp_color, ...
                                 'LineWidth', 1.8);
                        end
                    end

                    title(ax, sprintf('%.0f–%.0f ms | %d spikes', ...
                          bin_start, bin_end, spike_count), ...
                          'Interpreter', 'none');

                    xlabel(ax, 'Waveform time (ms)');
                    ylabel(ax, 'uV');

                    ylim(ax, y_lim);
                    yticks(ax, linspace(y_lim(1), y_lim(2), 3));
                    xticks(ax, round(linspace(t_wave(1), t_wave(end), 3), 2));

                    grid(ax, 'off');
                    axis(ax, 'square');

                    % Mark if this latency bin contains the final stimulation event.
                    if isfinite(finalEvent_ms) && ...
                       finalEvent_ms >= bin_start && finalEvent_ms < bin_end

                        text(ax, 0.05, 0.90, 'Final event bin', ...
                             'Units', 'normalized', ...
                             'FontSize', 8, ...
                             'FontWeight', 'bold');
                    end
                end

                %% ----- Legend -----
                legend_handles = gobjects(numel(Amps_selected),1);
                legend_labels  = cell(numel(Amps_selected),1);

                for aa = 1:numel(Amps_selected)

                    amp_global_idx = find(Amps_for_color == Amps_selected(aa), 1, 'first');

                    if isempty(amp_global_idx)
                        amp_color = [0 0 0];
                    else
                        amp_color = cmap(amp_global_idx,:);
                    end

                    legend_handles(aa) = plot(nan, nan, '-', ...
                        'Color', amp_color, ...
                        'LineWidth', 1.8);

                    legend_labels{aa} = sprintf('%g uA', Amps_selected(aa));
                end

                legend(legend_handles, legend_labels, ...
                       'Location', 'northeastoutside', ...
                       'Box', 'off');

            end % train level
        end % family
    end % set
end % channel

fprintf('\nFinished pulse-train SSD-recovered waveform plotting.\n');

%% ====================== LOCAL HELPER FUNCTIONS ======================

function setLabel = buildSetLabelAcrossLevels( ...
    setID, autoFlag, Amps_selected, Levels_this_family, ...
    stimSet_trial, trialAmps, trainLevel_trial, ...
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

    stimNames_all = {};

    for aa = 1:numel(Amps_selected)

        amp_val = Amps_selected(aa);

        for li = 1:numel(Levels_this_family)

            level_val = Levels_this_family(li);

            tr_label = find(stimSet_trial == setID & ...
                            trialAmps == amp_val & ...
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