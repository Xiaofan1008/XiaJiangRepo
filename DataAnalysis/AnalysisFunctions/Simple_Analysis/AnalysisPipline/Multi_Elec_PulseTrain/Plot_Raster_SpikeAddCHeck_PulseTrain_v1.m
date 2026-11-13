%% ============================================================
%  Pulse-Train Spike Recovery Check Plot
%
%  Purpose:
%    Plot recovered spike rasters and PSTHs to check whether pulse-train
%    spike recovery looks correct.
%
%  Input spike file:
%    *.sp_xia_PrefixRecovery.mat
%
%  Required spike variable:
%    sp_seq
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

% Recording channels to plot, using Depth_s index.
channels_to_plot = 20:25;

Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Amplitudes to plot.
% [] means all non-zero amplitudes.
amps_to_plot = [10];

% Train/event levels to plot.
% [] means all detected levels.
% Example:
%   [1 2 3 4 5 6] for sequential pulse train.
%   AutoSim may only have [1 2 3], and the code will automatically skip
%   levels that do not exist.
train_levels_to_plot = [];

% Stimulation set IDs to plot.
% [] means all detected sets.
sets_to_plot = [];

% Plot sequential/interleaved condition?
plot_sequential = true;

% Plot AutoSim/simultaneous condition?
% AutoSim is plotted in a separate figure from sequential.
plot_auto_sim = true;

% Include zero-current trials?
% Usually false for recovery checking.
include_zero_control = false;

% Number of trials to plot in each raster.
nTrials_to_plot = 30;

% Plotting time window around trigger.
ras_win = [-5 40];      % ms

% PSTH settings.
bin_ms = 1;             % bin size for PSTH
smooth_ms = 5;          % smoothing width

% Raster line settings.
raster_line_width = 1.1;

% Print example trial information to confirm metadata alignment.
debug_print_trial_content = true;

%% ====================== CHECK FOLDER ======================

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end

cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ====================== LOAD RECOVERED SPIKES ======================

rec_file = dir('*sp_xia_PrefixRecovery.mat');
assert(~isempty(rec_file), 'No *sp_xia_PrefixRecovery.mat file found. Run the recovery code first.');

load(rec_file(1).name, 'sp_seq');
fprintf('Loaded recovered spike file: %s\n', rec_file(1).name);

nChn_spike = numel(sp_seq);

%% ====================== LOAD TRIGGERS ======================

if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD EXPERIMENT PARAMETERS ======================

fDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fDIR), 'No *_exp_datafile_*.mat found.');

S = load(fDIR(1).name);

StimParams        = S.StimParams;
TrialParams       = S.TrialParams;
simultaneous_stim = S.simultaneous_stim;   % rows/slots per trial
n_Trials          = S.n_Trials;
E_MAP             = S.E_MAP;

if isfield(S, 'StimMeta')
    StimMeta = S.StimMeta;
else
    error('StimMeta was not found. This pulse-train check code requires StimMeta.');
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== SAMPLING RATE ======================

try
    [~, freq_params] = read_Intan_RHS2000_file;
    FS = freq_params.amplifier_sample_rate;
catch
    FS = 30000;
    warning('Could not read info.rhs. Using FS = 30000 Hz.');
end

fprintf('Sampling rate: %.1f Hz\n', FS);

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

%% ====================== TRIAL CONDITION IDS ======================
% TrialParams columns:
%   column 1 = trial number
%   column 2 = condition ID
%   column 3 = internal electrode ID
%
% Each trial has simultaneous_stim rows.
% We use the first row of each trial to get the trial-level condition ID.

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of trial-level condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA FROM StimMeta ======================
% This replaces the old prefix metadata columns 26–30.
%
% New pulse-train metadata:
%   stimSet_trial       = stimulation set/order ID
%   trainLevel_trial    = event/train level
%   trialAmps           = current amplitude
%   isAutoSim_trial     = AutoSim/simultaneous condition flag
%   isZeroControl       = zero-current control flag
%   eventTimes_ms_trial = event times in ms
%   eventEnd_ms_trial   = last event time in ms

stimSet_trial        = NaN(n_Trials,1);
trainLevel_trial     = NaN(n_Trials,1);
totalLevels_trial    = NaN(n_Trials,1);
trialAmps            = NaN(n_Trials,1);
isAutoSim_trial      = false(n_Trials,1);
isZeroControl        = false(n_Trials,1);
eventEnd_ms_trial    = NaN(n_Trials,1);

eventTimes_ms_trial  = cell(n_Trials,1);
pulseCount_trial     = cell(n_Trials,1);

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

    if isfield(meta, 'PulseCountPerElectrode')
        pulseCount_trial{tr} = meta.PulseCountPerElectrode;
    else
        pulseCount_trial{tr} = [];
    end

    % Use actual randomized StimParams rows to get amplitude.
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

%% ====================== DEBUG TRIAL CONTENT CHECK ======================

if debug_print_trial_content

    fprintf('\n================ DEBUG PULSE-TRAIN RECOVERY CHECK ================\n');

    all_sets_debug = unique(stimSet_trial(~isnan(stimSet_trial)));
    all_sets_debug = all_sets_debug(all_sets_debug > 0);

    all_amps_debug = unique(trialAmps(~isnan(trialAmps)));
    all_amps_debug = all_amps_debug(all_amps_debug > 0);

    all_levels_debug = unique(trainLevel_trial(~isnan(trainLevel_trial)));

    if ~isempty(all_sets_debug) && ~isempty(all_amps_debug) && ~isempty(all_levels_debug)

        debug_set = all_sets_debug(1);
        debug_amp = all_amps_debug(end);

        for autoFlag = [0 1]

            if autoFlag == 0
                debug_family = 'Seq';
            else
                debug_family = 'AutoSim';
            end

            fprintf('\n--- Debug family: %s ---\n', debug_family);

            for li = 1:numel(all_levels_debug)

                level_val = all_levels_debug(li);

                tr_debug = find(stimSet_trial == debug_set & ...
                                trialAmps == debug_amp & ...
                                trainLevel_trial == level_val & ...
                                isAutoSim_trial == logical(autoFlag) & ...
                                isZeroControl == 0, ...
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
                    debug_family, debug_set, debug_amp, level_val, tr_debug, conditionID_trial(tr_debug));

                fprintf('  Stim label:              %s\n', stimLabel_debug);

                fprintf('  Actual active channels:  ');
                disp(stimNames_active');

                fprintf('  EventTimesIncluded_ms:   ');
                disp(eventTimes_ms_trial{tr_debug});

                fprintf('  EventEndTime_ms:         %.3f\n', eventEnd_ms_trial(tr_debug));

                fprintf('  PulseCountPerElectrode:  ');
                disp(pulseCount_trial{tr_debug});

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

    fprintf('===================================================================\n');
end

%% ====================== APPLY USER FILTERS ======================

% Amplitudes.
all_amps = unique(trialAmps(~isnan(trialAmps)));
all_amps = all_amps(all_amps > 0);   % exclude zero-current

if isempty(amps_to_plot)
    amps_sel = all_amps;
else
    amps_sel = intersect(all_amps, amps_to_plot);
end

% Sets.
all_sets = unique(stimSet_trial(~isnan(stimSet_trial)));
all_sets = all_sets(all_sets > 0);

if isempty(sets_to_plot)
    set_sel = all_sets;
else
    set_sel = intersect(all_sets, sets_to_plot);
end

% Train/event levels.
all_levels = unique(trainLevel_trial(~isnan(trainLevel_trial)));
all_levels = sort(all_levels(:))';

if isempty(train_levels_to_plot)
    level_sel = all_levels;
else
    level_sel = intersect(all_levels, train_levels_to_plot);
    level_sel = sort(level_sel(:))';
end

fprintf('\nSelected amplitudes: ');
disp(amps_sel');

fprintf('Selected sets: ');
disp(set_sel');

fprintf('Selected train/event levels: ');
disp(level_sel);

fprintf('Plot sequential: %d\n', plot_sequential);
fprintf('Plot AutoSim:    %d\n', plot_auto_sim);

%% ====================== PSTH SETTINGS ======================

edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms / 1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

%% ====================== CHANNEL MAP ======================

d = Depth_s(Electrode_Type);

%% ====================== MAIN PLOTTING LOOP ======================
% One figure:
%   recording channel × set × amplitude × condition family
%
% Condition family:
%   Seq     = interleaved/sequential pulse train
%   AutoSim = simultaneous pulse train
%
% Raster subplots:
%   one subplot per train/event level
%
% Bottom subplot:
%   PSTH overlay across train/event levels

for ich = 1:length(channels_to_plot)

    ch_plot = channels_to_plot(ich);
    ch_spike = d(ch_plot);

    if ch_spike > nChn_spike
        warning('Channel %d maps to spike channel %d, but sp_seq only has %d channels. Skipped.', ...
                ch_plot, ch_spike, nChn_spike);
        continue;
    end

    if isempty(sp_seq{ch_spike})
        fprintf('Channel %d has no spikes. Skipped.\n', ch_plot);
        continue;
    end

    spike_times_ch = sp_seq{ch_spike}(:,1);

    for i_set = 1:numel(set_sel)

        set_id = set_sel(i_set);

        for i_amp = 1:numel(amps_sel)

            amp_val = amps_sel(i_amp);

            % Plot sequential and AutoSim as separate figure families.
            familyList = {};

            if plot_sequential
                familyList{end+1} = 'seq';
            end

            if plot_auto_sim
                familyList{end+1} = 'autosim';
            end

            for i_family = 1:numel(familyList)

                familyName = familyList{i_family};

                if strcmp(familyName, 'seq')
                    autoFlag = false;
                    familyTitle = 'Seq';
                    stimJoinMode = false;   % false means use arrow: Ch22→Ch18
                elseif strcmp(familyName, 'autosim')
                    autoFlag = true;
                    familyTitle = 'AutoSim';
                    stimJoinMode = true;    % true means use plus: Ch22+Ch18
                else
                    continue;
                end

                %% ----- Build level condition list for this family -----
                plotLevels = [];

                for li = 1:numel(level_sel)

                    level_val = level_sel(li);

                    tlist_check = find(stimSet_trial == set_id & ...
                                       trialAmps == amp_val & ...
                                       trainLevel_trial == level_val & ...
                                       isAutoSim_trial == autoFlag);

                    if ~include_zero_control
                        tlist_check = tlist_check(isZeroControl(tlist_check) == 0);
                    end

                    tlist_check = tlist_check(tlist_check <= nTrials_use);

                    if ~isempty(tlist_check)
                        plotLevels(end+1) = level_val; %#ok<AGROW>
                    end
                end

                if isempty(plotLevels)
                    continue;
                end

                nRaster = numel(plotLevels);

                %% ----- Build stimulation set label using all levels in this set/family -----
                
                stimNames_all_for_label = {};
                
                for li_label = 1:numel(plotLevels)
                
                    level_label = plotLevels(li_label);
                
                    tr_label = find(stimSet_trial == set_id & ...
                                    trialAmps == amp_val & ...
                                    trainLevel_trial == level_label & ...
                                    isAutoSim_trial == autoFlag & ...
                                    isZeroControl == 0, ...
                                    1, 'first');
                
                    if isempty(tr_label)
                        continue;
                    end
                
                    rr_label = (tr_label-1)*simultaneous_stim + (1:simultaneous_stim);
                
                    amp_label = cell2mat(StimParams_data(rr_label,16));
                    activeRows = amp_label > 0;
                
                    stimNames_this = StimParams_data(rr_label(activeRows),1);
                
                    stimNames_all_for_label = [stimNames_all_for_label; stimNames_this(:)]; %#ok<AGROW>
                end
                
                if isempty(stimNames_all_for_label)
                    set_label = sprintf('Set%d', set_id);
                else
                    % Remove duplicates while keeping order.
                    stimNames_all_for_label = unique(stimNames_all_for_label, 'stable');
                
                    set_label = buildStimLabelFromStimNames(stimNames_all_for_label, E_MAP, autoFlag);
                end

                %% ----- Create figure -----
                figName = sprintf('PulseTrainRecoveryCheck_%s_Ch%d_Set%d_Amp%g', ...
                                  familyTitle, ch_plot, set_id, amp_val);

                figure('Name', figName, ...
                       'Color', 'w', ...
                       'Position', [100 100 900 900]);

                tiledlayout(nRaster + 1, 1, ...
                            'TileSpacing', 'compact', ...
                            'Padding', 'compact');

                % Colors for each train level.
                cond_colors = lines(nRaster);

                % Store PSTH curves for bottom plot.
                psth_all = cell(nRaster,1);
                psth_labels = cell(nRaster,1);
                maxRate = 0;

                %% ====================== RASTER SUBPLOTS ======================
                for cc = 1:nRaster

                    level_val = plotLevels(cc);
                    thisColor = cond_colors(cc,:);

                    tlist = find(stimSet_trial == set_id & ...
                                 trialAmps == amp_val & ...
                                 trainLevel_trial == level_val & ...
                                 isAutoSim_trial == autoFlag);

                    if ~include_zero_control
                        tlist = tlist(isZeroControl(tlist) == 0);
                    end

                    tlist = tlist(tlist <= nTrials_use);

                    if ~isempty(tlist)
                        tlist = tlist(1:min(nTrials_to_plot, numel(tlist)));
                    end

                    axR = nexttile;
                    hold(axR, 'on');
                    box(axR, 'off');

                    counts = zeros(1, numel(edges)-1);

                    if isempty(tlist)

                        title(axR, sprintf('Level %g | No trials', level_val), ...
                              'Interpreter', 'none');

                        psth_all{cc} = zeros(1, numel(ctrs));
                        psth_labels{cc} = sprintf('Level %g', level_val);

                    else

                        %% ----- Representative event times for this level -----
                        tr_rep = tlist(1);

                        eventTimes_this = eventTimes_ms_trial{tr_rep};
                        eventEnd_this   = eventEnd_ms_trial(tr_rep);

                        if isempty(eventTimes_this)
                            eventText = '[]';
                        else
                            eventText = num2str(eventTimes_this);
                        end

                        %% ----- Plot raster trials -----
                        for k = 1:numel(tlist)

                            tr = tlist(k);

                            if tr > length(trig)
                                continue;
                            end

                            t0 = trig(tr) / FS * 1000;

                            rel_t = spike_times_ch - t0;
                            rel_t = rel_t(rel_t >= ras_win(1) & rel_t <= ras_win(2));

                            % Raster: one horizontal row per trial.
                            for s = 1:numel(rel_t)
                                plot(axR, [rel_t(s) rel_t(s)], [k-0.4 k+0.4], ...
                                     'Color', thisColor, ...
                                     'LineWidth', raster_line_width);
                            end

                            % PSTH counts.
                            counts = counts + histcounts(rel_t, edges);
                        end

                        %% ----- Compute PSTH -----
                        rate = counts / (numel(tlist) * bin_s);
                        rate = filter(g, 1, rate);

                        psth_all{cc} = rate;
                        psth_labels{cc} = sprintf('Level %g', level_val);

                        maxRate = max(maxRate, max(rate));

                        %% ----- Raster title -----
                        title(axR, sprintf('Level %g | events [%s] ms | final %.1f ms | %d trials', ...
                                           level_val, eventText, eventEnd_this, numel(tlist)), ...
                                           'Interpreter', 'none');

                        %% ----- Event lines on raster -----
                        drawEventLines(axR, eventTimes_this, eventEnd_this);
                    end

                    %% ----- Raster axis formatting -----
                    xlim(axR, ras_win);
                    ylim(axR, [0 max(1, numel(tlist)+1)]);
                    ylabel(axR, 'Trial');

                    if cc < nRaster
                        set(axR, 'XTickLabel', []);
                    else
                        xlabel(axR, 'Time from trigger (ms)');
                    end

                end % level raster loop

                %% ====================== BOTTOM PSTH SUBPLOT ======================
                axP = nexttile;
                hold(axP, 'on');
                box(axP, 'off');

                for cc = 1:nRaster

                    if isempty(psth_all{cc})
                        continue;
                    end

                    plot(axP, ctrs, psth_all{cc}, ...
                         'Color', cond_colors(cc,:), ...
                         'LineWidth', 1.8);
                end

                %% ----- Event lines on bottom PSTH -----
                % Because levels are nested, using the union of all event
                % times gives the full event structure for this family.
                allEventTimes_forPSTH = [];

                for cc = 1:nRaster

                    level_val = plotLevels(cc);

                    tlist_tmp = find(stimSet_trial == set_id & ...
                                     trialAmps == amp_val & ...
                                     trainLevel_trial == level_val & ...
                                     isAutoSim_trial == autoFlag);

                    if ~include_zero_control
                        tlist_tmp = tlist_tmp(isZeroControl(tlist_tmp) == 0);
                    end

                    tlist_tmp = tlist_tmp(tlist_tmp <= nTrials_use);

                    if isempty(tlist_tmp)
                        continue;
                    end

                    allEventTimes_forPSTH = [allEventTimes_forPSTH, eventTimes_ms_trial{tlist_tmp(1)}]; %#ok<AGROW>
                end

                allEventTimes_forPSTH = unique(sort(allEventTimes_forPSTH));

                if isempty(allEventTimes_forPSTH)
                    finalEvent_forPSTH = NaN;
                else
                    finalEvent_forPSTH = max(allEventTimes_forPSTH);
                end

                drawEventLines(axP, allEventTimes_forPSTH, finalEvent_forPSTH);

                xlim(axP, ras_win);

                if maxRate <= 0
                    ylim(axP, [0 50]);
                else
                    ylim(axP, [0 max(50, ceil(maxRate*1.1/10)*10)]);
                end

                xlabel(axP, 'Time from trigger (ms)');
                ylabel(axP, 'Rate (sp/s)');
                title(axP, 'PSTH overlay across train levels', 'Interpreter', 'none');

                legend(axP, psth_labels, 'Box', 'off', 'Location', 'northeast');

                %% ====================== FIGURE TITLE ======================
                sgtitle(sprintf('%s recovery check | Rec Ch %d | Set %d: %s | %g uA', ...
                    familyTitle, ch_plot, set_id, set_label, amp_val), ...
                    'Interpreter', 'none');

            end % family
        end % amplitude
    end % set
end % channel

fprintf('\nFinished pulse-train recovery check plotting.\n');

%% ====================== LOCAL HELPER FUNCTIONS ======================

function drawEventLines(ax, eventTimes_ms, eventEnd_ms)
    % Draw event timing lines:
    %   0 ms: red dashed
    %   other events: black dashed
    %   final event: blue dotted
    %
    % This function is used for both raster and PSTH axes.

    axes(ax); %#ok<LAXES>

    % Trigger / first event line.
    xline(ax, 0, 'r--', 'LineWidth', 1);

    if ~isempty(eventTimes_ms)

        for ee = 1:numel(eventTimes_ms)

            thisEvent = eventTimes_ms(ee);

            % Avoid drawing another black line at 0 ms.
            if abs(thisEvent) < 1e-9
                continue;
            end

            xline(ax, thisEvent, 'k--', 'LineWidth', 1);
        end
    end

    % Final event line.
    if isfinite(eventEnd_ms) && eventEnd_ms > 0
        xline(ax, eventEnd_ms, 'b:', 'LineWidth', 1.5);
    end
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

    % Remove duplicates while keeping original order.
    stimNames_active = unique(stimNames_active, 'stable');

    labelParts = cell(1, numel(stimNames_active));

    for i = 1:numel(stimNames_active)

        stimName = stimNames_active{i};

        chNum = convertStimNameUsingEMap(stimName, E_MAP);

        if isnan(chNum)
            % If conversion fails, keep original name as a fallback.
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
    %   then channel number = row - 1.

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
    % This fallback is not preferred because the parsed number may not equal
    % the actual mapped channel number. It only prevents the code from
    % crashing if one name is missing from E_MAP.

    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
        warning('Stim channel %s was not found in E_MAP. Falling back to parsed number %d.', ...
            stimName, chNum);
    else
        warning('Stim channel %s was not found in E_MAP and could not be parsed.', stimName);
    end
end