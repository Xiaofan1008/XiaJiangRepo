%% ============================================================
%  Pulse-Train Simultaneous vs Sequential Report Plot
%
%  Purpose:
%    Generate supervisor-style figures comparing actual AutoSim and
%    sequential pulse-train responses.
%
%  Figure structure:
%    One figure = one recording channel x one amplitude
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train3';

% Recording channels to plot, using Depth_s index.
channels_to_plot = [1:2];

Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Spike source:
%   'recovery_ssd' = *_sp_xia_PrefixRecovery_SSD.mat, variable sp_corr
%   'recovery'     = *_sp_xia_PrefixRecovery.mat, variable sp_seq
spike_source = 'recovery_ssd';

% Amplitudes to plot.
% [] means all non-zero amplitudes.
amps_to_plot = [10];

% Optional: only plot selected stimulation set IDs.
% [] means all detected final AutoSim/Seq sets.
sets_to_plot = [];

% Plot sequential conditions?
plot_sequential = true;

% Plot AutoSim/simultaneous conditions?
plot_auto_sim = true;

% Include zero-current control?
include_zero_control = false;

% Number of trials to plot in each raster.
nTrials_to_plot = 30;

% Time window around trigger.
ras_win = [-10 60];      % ms

% PSTH settings.
bin_ms = 1;             % bin size, ms
smooth_ms = 2;          % smoothing width in bins/ms approximately

% Raster marker settings.
raster_marker_size = 6;

% Save figures?
save_figures = false;
fig_out_folder = fullfile(data_folder, 'ReportPlot_SimVsSeq_RasterPSTH');

% Print debug information?
debug_print_condition_summary = true;

%% ====================== COLOUR SETTINGS ======================

Color.AutoSim = [0.4660 0.6740 0.1880];    % green
Color.Seq1    = [0 0.4470 0.7410];         % blue
Color.Seq2    = [0.8500 0.3250 0.0980];    % orange
Color.Seq3    = [0.4940 0.1840 0.5560];    % purple
Color.Seq4    = [0.3010 0.7450 0.9330];    % cyan

Color.Event0      = [0.85 0 0];            % red
Color.EventOther  = [0.15 0.15 0.15];      % dark grey/black
Color.EventFinal  = [0 0.2 0.85];          % blue

%% ====================== CHECK FOLDER ======================

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end

if save_figures && ~isfolder(fig_out_folder)
    mkdir(fig_out_folder);
end

cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ====================== LOAD SPIKES ======================

switch lower(spike_source)

    case 'recovery_ssd'

        sp_file = dir('*sp_xia_PrefixRecovery_SSD.mat');
        assert(~isempty(sp_file), ...
            'No *sp_xia_PrefixRecovery_SSD.mat file found.');

        S_sp = load(sp_file(1).name);

        if isfield(S_sp, 'sp_corr')
            sp_use = S_sp.sp_corr;
        else
            error('sp_corr not found in %s.', sp_file(1).name);
        end

        fprintf('Loaded SSD recovered spike file: %s\n', sp_file(1).name);

    case 'recovery'

        sp_file = dir('*sp_xia_PrefixRecovery.mat');
        sp_file = sp_file(~contains({sp_file.name}, 'SSD'));

        assert(~isempty(sp_file), ...
            'No *sp_xia_PrefixRecovery.mat file found.');

        S_sp = load(sp_file(1).name);

        if isfield(S_sp, 'sp_seq')
            sp_use = S_sp.sp_seq;
        else
            error('sp_seq not found in %s.', sp_file(1).name);
        end

        fprintf('Loaded recovered spike file: %s\n', sp_file(1).name);

    otherwise
        error('Unknown spike_source: %s', spike_source);
end

nChn_spike = numel(sp_use);

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
simultaneous_stim = S.simultaneous_stim;
n_Trials          = S.n_Trials;

if isfield(S, 'E_MAP')
    E_MAP = S.E_MAP;
else
    E_MAP = [];
    warning('E_MAP not found. Stim labels may be less clean.');
end

if isfield(S, 'StimMeta')
    StimMeta = S.StimMeta;
else
    error('StimMeta was not found. This code requires StimMeta.');
end

fprintf('Loaded exp datafile: %s\n', fDIR(1).name);
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

firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

trialNumber_trial = cell2mat(TrialParams_data(firstRow_eachTrial,1)); %#ok<NASGU>
conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
conditionID_trial = conditionID_trial(:);

if numel(conditionID_trial) ~= n_Trials
    warning('Number of trial-level condition IDs does not match n_Trials.');
end

%% ====================== BUILD TRIAL METADATA ======================

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

    % Get actual amplitude from randomized StimParams rows.
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

%% ====================== APPLY USER FILTERS ======================

all_amps = unique(trialAmps(~isnan(trialAmps)));
all_amps = all_amps(all_amps > 0);

if isempty(amps_to_plot)
    amps_sel = all_amps;
else
    amps_sel = intersect(all_amps, amps_to_plot);
end

all_sets = unique(stimSet_trial(~isnan(stimSet_trial)));
all_sets = all_sets(all_sets > 0);

if isempty(sets_to_plot)
    set_sel = all_sets;
else
    set_sel = intersect(all_sets, sets_to_plot);
end

fprintf('\nSelected amplitudes: ');
disp(amps_sel');

fprintf('Selected sets: ');
disp(set_sel');

fprintf('Plot sequential: %d\n', plot_sequential);
fprintf('Plot AutoSim:    %d\n', plot_auto_sim);

%% ====================== PSTH SETTINGS ======================

edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms / 1000;

% Smooth kernel.
if smooth_ms <= 0
    g = 1;
else
    xg = -round(3*smooth_ms):round(3*smooth_ms);
    g = exp(-0.5*(xg / smooth_ms).^2);
    g = g / sum(g);
end

%% ====================== CHANNEL MAP ======================

cur_dir = pwd;

try
    cd(data_folder);
    d = Depth_s(Electrode_Type);
    cd(cur_dir);
catch ME
    cd(cur_dir);
    warning('Depth_s failed. Falling back to direct mapping.');
    warning('%s', ME.message);
    d = 1:max(channels_to_plot);
end

%% ====================== BUILD FINAL CONDITION TABLE ======================
% Each row in FinalCondTable is one condition to plot:
%   Family = 'AutoSim' or 'Seq'
%   SetID
%   Amp
%   FinalLevel = automatically detected max train level for that family/set/amp
%   TrialList
%   Label

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

for i_amp = 1:numel(amps_sel)

    amp_val = amps_sel(i_amp);

    %% ----- AutoSim final conditions -----
    if plot_auto_sim

        candidate_sets_auto = unique(stimSet_trial( ...
            trialAmps == amp_val & ...
            isAutoSim_trial == true & ...
            isZeroControl == false));

        candidate_sets_auto = candidate_sets_auto(~isnan(candidate_sets_auto));

        if ~isempty(sets_to_plot)
            candidate_sets_auto = intersect(candidate_sets_auto, set_sel);
        end

        for iset = 1:numel(candidate_sets_auto)

            set_id = candidate_sets_auto(iset);

            trial_base = find(trialAmps == amp_val & ...
                              stimSet_trial == set_id & ...
                              isAutoSim_trial == true);

            if ~include_zero_control
                trial_base = trial_base(isZeroControl(trial_base) == false);
            end

            trial_base = trial_base(trial_base <= nTrials_use);

            if isempty(trial_base)
                continue;
            end

            final_level = max(trainLevel_trial(trial_base));

            trial_final = trial_base(trainLevel_trial(trial_base) == final_level);

            if isempty(trial_final)
                continue;
            end

            tr_rep = trial_final(1);

            label = getStimLabelForTrial(tr_rep, StimParams_data, simultaneous_stim, E_MAP, true);

            newCond.Family     = 'AutoSim';
            newCond.SetID      = set_id;
            newCond.Amp        = amp_val;
            newCond.FinalLevel = final_level;
            newCond.TrialList  = trial_final(:)';
            newCond.Label      = label;
            newCond.EventTimes = eventTimes_ms_trial{tr_rep};
            newCond.EventEnd   = eventEnd_ms_trial(tr_rep);
            newCond.Color      = Color.AutoSim;

            FinalConds(end+1) = newCond; %#ok<SAGROW>
        end
    end

    %% ----- Sequential final conditions -----
    if plot_sequential

        candidate_sets_seq = unique(stimSet_trial( ...
            trialAmps == amp_val & ...
            isAutoSim_trial == false & ...
            isZeroControl == false));

        candidate_sets_seq = candidate_sets_seq(~isnan(candidate_sets_seq));

        if ~isempty(sets_to_plot)
            candidate_sets_seq = intersect(candidate_sets_seq, set_sel);
        end

        seq_counter = 0;

        for iset = 1:numel(candidate_sets_seq)

            set_id = candidate_sets_seq(iset);

            trial_base = find(trialAmps == amp_val & ...
                              stimSet_trial == set_id & ...
                              isAutoSim_trial == false);

            if ~include_zero_control
                trial_base = trial_base(isZeroControl(trial_base) == false);
            end

            trial_base = trial_base(trial_base <= nTrials_use);

            if isempty(trial_base)
                continue;
            end

            final_level = max(trainLevel_trial(trial_base));

            trial_final = trial_base(trainLevel_trial(trial_base) == final_level);

            if isempty(trial_final)
                continue;
            end

            tr_rep = trial_final(1);

            label = getStimLabelForTrial(tr_rep, StimParams_data, simultaneous_stim, E_MAP, false);

            seq_counter = seq_counter + 1;

            switch seq_counter
                case 1
                    thisColor = Color.Seq1;
                case 2
                    thisColor = Color.Seq2;
                case 3
                    thisColor = Color.Seq3;
                otherwise
                    thisColor = Color.Seq4;
            end

            newCond.Family     = 'Seq';
            newCond.SetID      = set_id;
            newCond.Amp        = amp_val;
            newCond.FinalLevel = final_level;
            newCond.TrialList  = trial_final(:)';
            newCond.Label      = label;
            newCond.EventTimes = eventTimes_ms_trial{tr_rep};
            newCond.EventEnd   = eventEnd_ms_trial(tr_rep);
            newCond.Color      = thisColor;

            FinalConds(end+1) = newCond; %#ok<SAGROW>
        end
    end
end

if isempty(FinalConds)
    error('No final conditions found. Check amplitude, set, and metadata filters.');
end

%% ====================== DEBUG CONDITION SUMMARY ======================

if debug_print_condition_summary

    fprintf('\n================ FINAL CONDITIONS TO PLOT ================\n');

    for i = 1:numel(FinalConds)

        fprintf('%2d | %s | Set %g | Amp %.1f uA | FinalLevel %g | Trials %d | Label %s\n', ...
            i, ...
            FinalConds(i).Family, ...
            FinalConds(i).SetID, ...
            FinalConds(i).Amp, ...
            FinalConds(i).FinalLevel, ...
            numel(FinalConds(i).TrialList), ...
            FinalConds(i).Label);

        fprintf('     Events: ');
        disp(FinalConds(i).EventTimes);
    end

    fprintf('==========================================================\n');
end

%% ====================== MAIN REPORT PLOTTING LOOP ======================

for ich = 1:length(channels_to_plot)

    ch_plot = channels_to_plot(ich);
    ch_spike = d(ch_plot);

    if ch_spike > nChn_spike
        warning('Channel %d maps to spike channel %d, but spike file only has %d channels. Skipped.', ...
                ch_plot, ch_spike, nChn_spike);
        continue;
    end

    if isempty(sp_use{ch_spike})
        fprintf('Channel %d has no spikes. Skipped.\n', ch_plot);
        continue;
    end

    spike_times_ch = sp_use{ch_spike}(:,1);

    for i_amp = 1:numel(amps_sel)

        amp_val = amps_sel(i_amp);

        cond_idx_this_amp = find([FinalConds.Amp] == amp_val);

        if isempty(cond_idx_this_amp)
            continue;
        end

        CondsPlot = FinalConds(cond_idx_this_amp);

        nRaster = numel(CondsPlot);

        figName = sprintf('Report_SimVsSeq_RasterPSTH_Ch%d_Amp%g', ...
                          ch_plot, amp_val);

        figure('Name', figName, ...
               'Color', 'w', ...
               'Position', [100 80 700 850]);

        tiledlayout(nRaster + 1, 1, ...
                    'TileSpacing', 'compact', ...
                    'Padding', 'compact');

        psth_all = cell(nRaster,1);
        psth_labels = cell(nRaster,1);
        maxRate = 0;

        %% ====================== RASTER PANELS ======================
        for ic = 1:nRaster

            C = CondsPlot(ic);

            tlist = C.TrialList;
            tlist = tlist(tlist <= nTrials_use);

            if ~isempty(tlist)
                tlist_plot = tlist(1:min(nTrials_to_plot, numel(tlist)));
            else
                tlist_plot = [];
            end

            axR = nexttile;
            hold(axR, 'on');
            box(axR, 'off');

            counts = zeros(1, numel(edges)-1);

            if isempty(tlist_plot)

                title(axR, sprintf('%s | Set %g | Level %g | No trials', ...
                    C.Family, C.SetID, C.FinalLevel), ...
                    'Interpreter', 'none');

                psth_all{ic} = zeros(1, numel(ctrs));

            else

                %% ----- Raster trials -----
                for k = 1:numel(tlist_plot)

                    tr = tlist_plot(k);

                    if tr > length(trig)
                        continue;
                    end

                    t0 = trig(tr) / FS * 1000;

                    rel_t = spike_times_ch - t0;
                    rel_t = rel_t(rel_t >= ras_win(1) & rel_t <= ras_win(2));

                    y = k .* ones(size(rel_t));

                    plot(axR, rel_t, y, '.', ...
                        'Color', C.Color, ...
                        'MarkerSize', raster_marker_size);

                    counts = counts + histcounts(rel_t, edges);
                end

                %% ----- PSTH -----
                rate = counts / (numel(tlist_plot) * bin_s);
                rate = conv(rate, g, 'same');

                psth_all{ic} = rate;
                maxRate = max(maxRate, max(rate));

                %% ----- Raster title -----
                eventText = num2str(C.EventTimes);

                title(axR, sprintf('%s | Set %g: %s | final level %g | events [%s] ms | %d trials', ...
                    C.Family, C.SetID, C.Label, C.FinalLevel, eventText, numel(tlist_plot)), ...
                    'Interpreter', 'none');

                %% ----- Event lines -----
                drawEventLines(axR, C.EventTimes, C.EventEnd, Color);
            end

            xlim(axR, ras_win);
            ylim(axR, [0 max(1, numel(tlist_plot)+1)]);
            ylabel(axR, 'Trial');

            if ic < nRaster
                set(axR, 'XTickLabel', []);
            else
                xlabel(axR, 'Time from trigger (ms)');
            end

            improveAxes(axR);

            psth_labels{ic} = sprintf('%s %s', C.Family, C.Label);
        end

        %% ====================== PSTH OVERLAY PANEL ======================
        axP = nexttile;
        hold(axP, 'on');
        box(axP, 'off');

        for ic = 1:nRaster

            if isempty(psth_all{ic})
                continue;
            end

            C = CondsPlot(ic);

            plot(axP, ctrs, psth_all{ic}, ...
                'Color', C.Color, ...
                'LineWidth', 2.0);
        end

        %% ----- Draw union of event lines -----
        allEventTimes = [];
        allEventEnds = [];

        for ic = 1:nRaster
            allEventTimes = [allEventTimes, CondsPlot(ic).EventTimes]; %#ok<AGROW>
            allEventEnds = [allEventEnds, CondsPlot(ic).EventEnd]; %#ok<AGROW>
        end

        allEventTimes = unique(sort(allEventTimes));
        finalEvent = max(allEventEnds);

        drawEventLines(axP, allEventTimes, finalEvent, Color);

        xlim(axP, ras_win);

        if maxRate <= 0
            ylim(axP, [0 50]);
        else
            ylim(axP, [0 max(50, ceil(maxRate*1.1/10)*10)]);
        end

        xlabel(axP, 'Time from trigger (ms)');
        ylabel(axP, 'Rate (sp/s)');
        title(axP, 'PSTH overlay: AutoSim vs Sequential final train levels', ...
            'Interpreter', 'none');

        legend(axP, psth_labels, ...
            'Box', 'off', ...
            'Location', 'northeast');

        improveAxes(axP);

        %% ====================== FIGURE TITLE ======================
        sgtitle(sprintf('Sim vs Seq pulse-train response | Rec Ch %d | %.1f uA', ...
            ch_plot, amp_val), ...
            'Interpreter', 'none');

        if save_figures
            saveas(gcf, fullfile(fig_out_folder, safeFigName(figName, 'png')));
        end

    end % amplitude
end % channel

fprintf('\nFinished Sim vs Seq report plotting.\n');

%% ============================================================
%  LOCAL HELPER FUNCTIONS
% ============================================================

function drawEventLines(ax, eventTimes_ms, eventEnd_ms, Color)
    % Draw event timing lines.
    %
    % 0 ms: red dashed
    % other events: black dashed
    % final event: blue dotted

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
                'LineWidth', 0.9, ...
                'Alpha', 0.45, ...
                'HandleVisibility', 'off');
        end
    end

    if isfinite(eventEnd_ms) && eventEnd_ms > 0

        xline(ax, eventEnd_ms, ':', ...
            'Color', Color.EventFinal, ...
            'LineWidth', 1.4, ...
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
    % Sequential:
    %   Ch22→Ch18
    %
    % AutoSim:
    %   Ch22+Ch18

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
    % Your E_MAP format:
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

function improveAxes(ax)

    ax.FontSize = 10;
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.Box = 'off';

    ax.GridAlpha = 0.15;
end

function outName = safeFigName(figTitle, ext)
    % Make safe file name from figure title.

    outName = regexprep(figTitle, '[^\w\d\-_. ]', '_');
    outName = strrep(outName, ' ', '_');
    outName = [outName '.' ext];
end