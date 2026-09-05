%% =============================================================
%  Pulse-Train Raster + PSTH All-Channel Viewer
%
%  Purpose:
%    Plot all recording channels in one figure for each pulse-train
%    condition/stimulation style.
% =============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ========================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

Electrode_Type   = 2;       % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Recording channels to plot, using Depth_s index.
raster_chn_start = 1;
raster_chn_end   = 64;

% Plotting windows.
ras_win        = [-20 60];   % ms, time relative to trigger
bin_ms_raster  = 1;         % ms, PSTH bin size
smooth_ms      = 3;         % ms, PSTH smoothing width

% Which amplitudes to plot.
% [] means all non-zero amplitudes.
Plot_Amps = [10];

% Which stimulation set/order IDs to plot.
% [] means all detected sets.
Plot_SetIDs = [];

% Plot AutoSim / simultaneous conditions?
plot_auto_sim = true;

% Plot sequential conditions?
plot_sequential = true;

% Include zero-current trials?
include_zero_control = false;

% Number of trials to plot in each raster.
nTrials_to_plot = 30;

% Fixed figure size.
% Format: [left bottom width height]
fig_position = [50 50 1300 850];

% Save figures?
save_figures = false;
fig_out_folder = fullfile(data_folder, 'AllChannel_RasterPSTH_PulseTrain');

% Print condition summary before plotting.
debug_print_condition_summary = true;

% Plot raster using dots or vertical tick marks.
% 'dot'  = compact
% 'tick' = like spike-event raster
raster_style = 'dot';

% Raster marker / line settings.
raster_marker_size = 4;
raster_line_width  = 0.7;

%% ====================== COLOUR SETTINGS ======================

Color.AutoSim = [0.4660 0.6740 0.1880];    % green
Color.Seq1    = [0 0.4470 0.7410];         % blue
Color.Seq2    = [0.8500 0.3250 0.0980];    % orange
Color.Seq3    = [0.4940 0.1840 0.5560];    % purple
Color.Seq4    = [0.3010 0.7450 0.9330];    % cyan

Color.Event0     = [0.85 0 0];             % red
Color.EventOther = [0.15 0.15 0.15];       % black/grey
Color.EventFinal = [0 0.2 0.85];           % blue

%% ===================== INITIAL SETUP =========================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end

if save_figures && ~isfolder(fig_out_folder)
    mkdir(fig_out_folder);
end

cd(data_folder);
fprintf('Raster+PSTH all-channel plotting in folder:\n%s\n\n', data_folder);

%% ===================== LOAD SPIKES ===========================

sp = [];

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

nCh = numel(sp);

fprintf('Loaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);
fprintf('Loaded %d spike channels.\n', nCh);

%% ===================== LOAD TRIGGERS =========================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ===================== LOAD EXPERIMENT PARAMETERS =============

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
    warning('E_MAP not found. Stim channel labels may be less clean.');
end

if isfield(S_exp, 'StimMeta')
    StimMeta = S_exp.StimMeta;
else
    error('StimMeta was not found. This pulse-train viewer requires StimMeta.');
end

fprintf('Loaded exp datafile: %s\n', fileDIR(1).name);
fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match nTrig (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ===================== SAMPLING RATE =========================

try
    [~, freq_params] = read_Intan_RHS2000_file;
    FS = freq_params.amplifier_sample_rate;
catch
    FS = 30000;
    warning('Could not read info.rhs. Using FS = 30000 Hz.');
end

fprintf('Sampling rate: %.1f Hz\n', FS);

%% ===================== REMOVE HEADER ROWS ====================

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

%% ===================== TRIAL CONDITION IDS ===================
% TrialParams columns:
%   column 1 = trial number
%   column 2 = condition ID
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

%% ===================== BUILD TRIAL METADATA FROM StimMeta =====
% New pulse-train metadata:
%   stimSet_trial       = stimulation set/order ID
%   trainLevel_trial    = event/train level
%   trialAmps           = current amplitude
%   isAutoSim_trial     = AutoSim/simultaneous flag
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

%% ===================== APPLY USER FILTERS =====================

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

fprintf('Plot AutoSim:    %d\n', plot_auto_sim);
fprintf('Plot sequential: %d\n', plot_sequential);

%% ===================== ELECTRODE MAP =========================

cur_dir = pwd;

try
    cd(data_folder);
    d = Depth_s(Electrode_Type);
    cd(cur_dir);
catch ME
    cd(cur_dir);
    warning('Depth_s failed. Falling back to direct mapping.');
    warning('%s', ME.message);
    d = 1:raster_chn_end;
end

depth_range = raster_chn_start : min(raster_chn_end, numel(d));
nChPlot = numel(depth_range);

fprintf('Plotting Depth_s channels %d to %d (%d channels).\n', ...
    depth_range(1), depth_range(end), nChPlot);

%% ===================== PSTH KERNEL ===========================

edges = ras_win(1) : bin_ms_raster : ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

if smooth_ms <= 0
    g = 1;
else
    xg = -round(3*smooth_ms):round(3*smooth_ms);
    g = exp(-0.5*(xg / smooth_ms).^2);
    g = g / sum(g);
end

%% ===================== BUILD FINAL CONDITION LIST =============
% Each FinalCond is one figure to plot:
%   AutoSim Set X final level
%   Seq Set Y final level
%   Seq Set Z final level
%
% The final train level is automatically chosen as the highest available
% TrainLevel for that amplitude, set, and family.

FinalConds = struct( ...
    'Family', {}, ...
    'SetID', {}, ...
    'Amp', {}, ...
    'FinalLevel', {}, ...
    'TrialList', {}, ...
    'Label', {}, ...
    'EventTimes', {}, ...
    'EventEnd', {}, ...
    'Color', {});

for ai = 1:numel(Plot_Amps_selected)

    ampVal = Plot_Amps_selected(ai);

    %% ---------------- AutoSim final conditions ----------------
    if plot_auto_sim

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trials_base = find(trialAmps == ampVal & ...
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

            label = getStimLabelForTrial( ...
                tr_rep, StimParams_data, simultaneous_stim, E_MAP, true);

            C.Family     = 'AutoSim';
            C.SetID      = setID;
            C.Amp        = ampVal;
            C.FinalLevel = finalLevel;
            C.TrialList  = trials_final(:)';
            C.Label      = label;
            C.EventTimes = eventTimes_ms_trial{tr_rep};
            C.EventEnd   = eventEnd_ms_trial(tr_rep);
            C.Color      = Color.AutoSim;

            FinalConds(end+1) = C; %#ok<SAGROW>
        end
    end

    %% ---------------- Sequential final conditions ----------------
    if plot_sequential

        seqCounter = 0;

        for si = 1:numel(SetIDs_selected)

            setID = SetIDs_selected(si);

            trials_base = find(trialAmps == ampVal & ...
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

            label = getStimLabelForTrial( ...
                tr_rep, StimParams_data, simultaneous_stim, E_MAP, false);

            seqCounter = seqCounter + 1;

            switch seqCounter
                case 1
                    thisColor = Color.Seq1;
                case 2
                    thisColor = Color.Seq2;
                case 3
                    thisColor = Color.Seq3;
                otherwise
                    thisColor = Color.Seq4;
            end

            C.Family     = 'Seq';
            C.SetID      = setID;
            C.Amp        = ampVal;
            C.FinalLevel = finalLevel;
            C.TrialList  = trials_final(:)';
            C.Label      = label;
            C.EventTimes = eventTimes_ms_trial{tr_rep};
            C.EventEnd   = eventEnd_ms_trial(tr_rep);
            C.Color      = thisColor;

            FinalConds(end+1) = C; %#ok<SAGROW>
        end
    end
end

if isempty(FinalConds)
    error('No final conditions were found. Check Plot_Amps, Plot_SetIDs, and StimMeta.');
end

%% ===================== DEBUG CONDITION SUMMARY ================

if debug_print_condition_summary

    fprintf('\n================ FINAL CONDITIONS TO PLOT ================\n');

    for i = 1:numel(FinalConds)

        fprintf('%2d | %s | Set %g | Amp %.1f uA | FinalLevel %g | Trials %d | Label %s\n', ...
            i, FinalConds(i).Family, FinalConds(i).SetID, FinalConds(i).Amp,  FinalConds(i).FinalLevel, ...
            numel(FinalConds(i).TrialList), FinalConds(i).Label);

        fprintf('     Events: ');
        disp(FinalConds(i).EventTimes);
    end

    fprintf('==========================================================\n');
end

%% ========================================================================
%  MAIN CONDITION LOOP
%% ========================================================================

for ci = 1:numel(FinalConds)

    C = FinalConds(ci);

    trials_this_all = C.TrialList;
    trials_this_all = trials_this_all(trials_this_all <= nTrials_use);

    if isempty(trials_this_all)
        continue;
    end

    trials_plot = trials_this_all(1:min(nTrials_to_plot, numel(trials_this_all)));

    figTitle = sprintf('%s | Set %g: %s | Final level %g | Amp %.1f uA | nTrials=%d',C.Family, C.SetID, C.Label, C.FinalLevel, C.Amp, numel(trials_this_all));

    figName = sprintf('AllCh_%s_Set%g_Level%g_Amp%g',  C.Family, C.SetID, C.FinalLevel, C.Amp);

    figure('Color','w', 'Name', figName, 'Position', fig_position);

    tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');

    sgtitle(figTitle, 'FontSize', 14, 'shouhFontWeight','bold', 'Interpreter','none');

    %% ===================== CHANNEL LOOP =====================

    for idxDepth = 1:nChPlot

        ich = depth_range(idxDepth);   % Depth_s index
        ch  = d(ich);                  % spike channel index

        ax = nexttile;
        hold(ax, 'on');
        box(ax, 'off');

        if ch < 1 || ch > nCh || isempty(sp{ch})
            axis(ax, 'off');
            continue;
        end

        sp_times = sp{ch}(:,1);

        [allTrialSpikes, rate_s, yMaxPSTH] = get_raster_psth_for_channel(sp_times, trig, FS, trials_plot, ras_win, edges, bin_s, ctrs, g);

        nTr = numel(trials_plot);

        %% ----- Plot PSTH -----
        yyaxis(ax, 'left');

        if any(rate_s)
            % plot(ax, ctrs, rate_s, 'Color', C.Color, 'LineWidth', 1.4);
            plot(ax, ctrs, rate_s, 'Color', [0 0 0], 'LineWidth', 1.4);
        end

        xlim(ax, ras_win);
        ylim(ax, [0 yMaxPSTH]);
        ylabel(ax, 'Rate');

        ax.YAxis(1).Color = C.Color;

        %% ----- Plot raster -----
        yyaxis(ax, 'right');

        for ti = 1:nTr

            tt = allTrialSpikes{ti};

            if isempty(tt)
                continue;
            end

            switch lower(raster_style)

                case 'dot'
                    % plot(ax, tt, ti*ones(size(tt)), '.','Color', C.Color,'MarkerSize', raster_marker_size);
                    plot(ax, tt, ti*ones(size(tt)), '.','Color', [0.35 0.35 0.35],'MarkerSize', raster_marker_size);

                case 'tick'
                    for ss = 1:numel(tt)
                        plot(ax, [tt(ss) tt(ss)], [ti-0.35 ti+0.35],'Color', C.Color,'LineWidth', raster_line_width);
                    end

                otherwise
                    error('Unknown raster_style: %s', raster_style);
            end
        end

        ylim(ax, [0 nTr+1]);
        set(ax, 'YTick', []);
        ax.YAxis(2).Color = [0.2 0.2 0.2];

        %% ----- Event lines -----
        % drawEventLines(ax, C.EventTimes, C.EventEnd, Color);

        xlim(ax, ras_win);

        %% ----- Title -----
        title(ax, sprintf('Ch %d', ich), ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');

        %% ----- Axis formatting -----
        ax.FontSize = 9;
        ax.LineWidth = 0.8;
        ax.TickDir = 'out';
        ax.Box = 'off';

        % No grid lines.

        if idxDepth > nChPlot - ceil(sqrt(nChPlot))
            xlabel(ax, 'Time (ms)');
        end
    end

    if save_figures
        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, fullfile(fig_out_folder, safeFigName(figName, 'png')));
    end
end

fprintf('\nFinished pulse-train all-channel raster + PSTH plotting.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function [allTrialSpikes, rate_s, yMaxPSTH] = get_raster_psth_for_channel( ...
    sp_times, trig, FS, trials_this, ras_win, edges, bin_s, ctrs, g)

    allTrialSpikes = cell(numel(trials_this),1);
    counts = zeros(1, numel(edges)-1);

    for ti = 1:numel(trials_this)

        tr = trials_this(ti);

        if tr > numel(trig)
            continue;
        end

        t0 = trig(tr) / FS * 1000;

        tt = sp_times;
        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;

        allTrialSpikes{ti} = tt(:)';
        counts = counts + histcounts(tt, edges);
    end

    if isempty(trials_this) || ~any(counts)
        rate_s = zeros(size(ctrs));
    else
        rate = counts / (numel(trials_this) * bin_s);
        rate_s = conv(rate, g, 'same');
    end

    maxRate = max(rate_s);

    if maxRate <= 0
        yMaxPSTH = 50;
    else
        yMaxPSTH = max(50, ceil(maxRate * 1.15 / 10) * 10);
    end
end

function drawEventLines(ax, eventTimes_ms, eventEnd_ms, Color)
    % Draw event timing lines:
    %   0 ms: red dashed
    %   other events: dark dashed
    %   final event: blue dotted

    yl = ylim(ax);

    xline(ax, 0, '--', ...
        'Color', Color.Event0, ...
        'LineWidth', 1.0, ...
        'Alpha', 0.65, ...
        'HandleVisibility', 'off');

    if ~isempty(eventTimes_ms)

        for ee = 1:numel(eventTimes_ms)

            thisEvent = eventTimes_ms(ee);

            if abs(thisEvent) < 1e-9
                continue;
            end

            xline(ax, thisEvent, '--', ...
                'Color', Color.EventOther, ...
                'LineWidth', 0.8, ...
                'Alpha', 0.45, ...
                'HandleVisibility', 'off');
        end
    end

    if isfinite(eventEnd_ms) && eventEnd_ms > 0

        xline(ax, eventEnd_ms, ':', ...
            'Color', Color.EventFinal, ...
            'LineWidth', 1.2, ...
            'Alpha', 0.75, ...
            'HandleVisibility', 'off');
    end

    ylim(ax, yl);
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
    % Build stimulation label for figure title.
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

function outName = safeFigName(figTitle, ext)
    % Make safe file name from figure title.

    outName = regexprep(figTitle, '[^\w\d\-_. ]', '_');
    outName = strrep(outName, ' ', '_');
    outName = [outName '.' ext];
end