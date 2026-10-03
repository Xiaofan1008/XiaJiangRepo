%% ========================================================================
% RASTER + PSTH, ALL CHANNELS, ONE FIGURE PER (CONDITION x CHANNEL BLOCK)
%
% PURPOSE
%   Visual QC / responding-channel-by-eye tool. For every selected
%   stimulation condition (set x amplitude x PTD), plot every recording
%   channel's raster + PSTH so you can judge responsiveness and spot bad
%   trials/channels before confirming the responding-channel list.
%
% INPUT CONVENTION (matches SpikeFiltering_PCA_ReasonSave_v3.m /
% SpikeRecovery_SequentialStitch_v1.m)
%   spike_file : one of
%       *.sp_xia_filtered.mat            (sp_clipped,  times only)
%       *.sp_xia_waveforms_filtered.mat  (sp_waveforms, times + shape)
%       *.sp_xia_stitched.mat            (sp_clipped,  stitched sequential)
%       *.sp_xia_waveforms_stitched.mat  (sp_waveforms, stitched sequential)
%   Either sp_clipped or sp_waveforms is accepted; only the time column is
%   used, so pointing at the (smaller) times-only file is enough.
%   Cell array indexed by PHYSICAL recording channel (not map position).
%
%   dataset_folder : the folder containing this recording's own
%       *_exp_datafile_*.mat (StimParams, simultaneous_stim, E_MAP,
%       n_Trials) and trigger file - loaded via loadTrig(0), exactly as in
%       your existing scripts. For a stitched sequential file this MUST be
%       the sequential dataset's own folder (its own triggers/StimParams),
%       not the single-pulse dataset's.
%
%   Electrode_Type / ChnMap(Electrode_Type, nCh) maps map position ->
%   physical channel, same convention as SpikeFiltering_PCA_ReasonSave_v3.m.
%   Probe layout differs by animal, so this is a required setting, not
%   auto-detected.
%
% CONDITIONS
%   Decoded the same way as RespondingChn_identify_SimSeq_v2.m: unique
%   stimulation-channel orders are "sets" (A->B and B->A are different
%   sets), crossed with amplitude and, when simultaneous_stim >= 2, PTD.
%   simultaneous_stim == 1 (single-pulse files) is also supported: no PTD
%   axis in that case.
%
% LAYOUT
%   Channels are plotted in plain channel-number (map-position) order -
%   NOT arranged by shank/probe geometry. One figure holds at most
%   Max_Channels_Per_Figure channels; a condition with more channels than
%   that is split into further figures (e.g. 96 -> 64 + 32), each clearly
%   named with the condition and its channel range.
%
% NOT YET INCLUDED (future scripts in this pipeline)
%   - Responding-channel highlighting/overlay (needs the auto-reference +
%     manual-confirm responding-channel tools, not built yet)
%   - Bad-trial exclusion (needs the dedicated condition-aware bad-trial
%     tool, not built yet). Every trial in the condition is used as-is.
% ========================================================================

clear;
% close all;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_stitched.mat';
dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1';

% Electrode type (must match this dataset - probe differs by animal):
%   0 = rigid single-shank probe
%   1 = flexible single-shank probe
%   2 = four-shank flexible probe
%   3 = 64chn+32chn hybrid
Electrode_Type = 3;

FS_fallback = 30000;   % used only if fs/FS is absent from spike_file

% Empty = plot every map-position channel, in channel-number order
Plot_Channels = [];

% Empty means all available sets/amplitudes/PTDs
Plot_Sets = [1];
Plot_Amps = [10];
Plot_PTDs = [];

Condition_Tolerance = 1e-4;

% Channels per figure. A condition with more channels than this is split
% into further figures (64, then the remainder), never one giant figure.
Max_Channels_Per_Figure = 64;

%% --------------------------- FIGURE OPTIONS ---------------------------

% 'docked' : every figure is a tab in one MATLAB figure container
% 'normal' : every figure is its own window
Figure_Window_Style = 'docked';
fig_position = [50 50 1600 900];   % only used when Figure_Window_Style = 'normal'
Build_Figures_Invisibly = true;    % build off-screen, then reveal (faster)

%% --------------------------- PLOTTING OPTIONS -------------------------

ras_win = [-50 80];     % ms, relative to the first stimulation pulse (t=0)
bin_ms_raster = 1;
smooth_ms = 5;

Raster_Marker_Size = 4;
PSTH_Line_Width = 1.4;
Minimum_PSTH_YMax = 50;

%% =========================== INITIAL CHECKS ===========================

if ~isfile(spike_file)
    error('Spike file does not exist:\n%s',spike_file);
end
if ~isfolder(dataset_folder)
    error('Dataset folder does not exist:\n%s',dataset_folder);
end
if ~isnumeric(ras_win) || numel(ras_win) ~= 2 || any(~isfinite(ras_win)) || ras_win(2) <= ras_win(1)
    error('ras_win must be [start end], with end > start.');
end
if ~isscalar(bin_ms_raster) || ~isfinite(bin_ms_raster) || bin_ms_raster <= 0
    error('bin_ms_raster must be positive.');
end
if ~isscalar(smooth_ms) || ~isfinite(smooth_ms) || smooth_ms <= 0
    error('smooth_ms must be positive.');
end
if ~isscalar(Max_Channels_Per_Figure) || Max_Channels_Per_Figure < 1
    error('Max_Channels_Per_Figure must be a positive integer.');
end
valid_window_styles = {'docked','normal'};
if ~any(strcmpi(Figure_Window_Style,valid_window_styles))
    error('Figure_Window_Style must be ''docked'' or ''normal''.');
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('RASTER + PSTH, ALL CHANNELS\n');
fprintf('============================================================\n');
fprintf('Spike file: %s\n',spike_file);
fprintf('Dataset folder: %s\n',dataset_folder);
fprintf('Figure style: %s\n',Figure_Window_Style);
fprintf('Raster window: [%g,%g) ms\n',ras_win);
fprintf('Max channels per figure: %d\n',Max_Channels_Per_Figure);

%% ========================= LOAD SPIKE TIMES ============================
% Accept either sp_clipped (times only) or sp_waveforms (times + shape);
% only column 1 (time, ms) is needed for raster/PSTH.

spike_vars_present = who('-file',spike_file);

if ismember('sp_clipped',spike_vars_present)
    SpikeLoad = load(spike_file,'sp_clipped');
    sp_raw = SpikeLoad.sp_clipped;
    spike_source_var = 'sp_clipped';
elseif ismember('sp_waveforms',spike_vars_present)
    SpikeLoad = load(spike_file,'sp_waveforms');
    sp_raw = SpikeLoad.sp_waveforms;
    spike_source_var = 'sp_waveforms';
else
    error('Spike file contains neither sp_clipped nor sp_waveforms:\n%s',spike_file);
end

nCh = numel(sp_raw);

sp = cell(nCh,1);
for ch = 1:nCh
    cellData = sp_raw{ch};
    if isempty(cellData)
        sp{ch} = [];
    else
        spike_times = double(cellData(:,1));
        spike_times = spike_times(isfinite(spike_times));
        if ~issorted(spike_times)
            spike_times = sort(spike_times);
        end
        sp{ch} = spike_times;
    end
end
clear sp_raw SpikeLoad;

if ismember('fs',spike_vars_present)
    FsLoad = load(spike_file,'fs');
    FS = double(FsLoad.fs);
elseif ismember('FS',spike_vars_present)
    FsLoad = load(spike_file,'FS');
    FS = double(FsLoad.FS);
else
    warning('fs/FS not found in spike file - using fallback FS = %g Hz.',FS_fallback);
    FS = FS_fallback;
end

fprintf('\nSpike variable used: %s\n',spike_source_var);
fprintf('Physical channels in spike file: %d\n',nCh);
fprintf('Sampling rate: %g Hz\n',FS);

%% ========================= ELECTRODE MAPPING ===========================

map_nums_plus = double(ChnMap(Electrode_Type,nCh));
map_nums_plus = map_nums_plus(:);

if isempty(map_nums_plus)
    error('ChnMap(%d,%d) returned an empty channel map.',Electrode_Type,nCh);
end

nMapPos = numel(map_nums_plus);
fprintf('Map positions: %d\n',nMapPos);

%% =========================== LOAD TRIGGERS ============================

cd(dataset_folder);

if isempty(dir('*.trig.dat'))
    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);
end

trigger_files = dir('*.trig.dat');
if isempty(trigger_files)
    error('No *.trig.dat file could be found or created in:\n%s',dataset_folder);
end

trig = double(loadTrig(0));
trig = trig(:);
nTrig = numel(trig);
trig_ms = trig/FS*1000;

fprintf('Triggers loaded: %d\n',nTrig);

%% ==================== LOAD EXPERIMENT PARAMETERS =====================

experiment_files = dir('*_exp_datafile_*.mat');
experiment_files = remove_metadata_and_backups(experiment_files);

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found in:\n%s',dataset_folder);
end
if numel(experiment_files) > 1
    error('Multiple current experiment files were found in:\n%s',dataset_folder);
end

experiment_file = experiment_files(1).name;

ExpLoad = load(experiment_file,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required_variables = {'StimParams','simultaneous_stim','E_MAP','n_Trials'};
for variable_index = 1:numel(required_variables)
    if ~isfield(ExpLoad,required_variables{variable_index})
        error('Experiment file is missing variable: %s',required_variables{variable_index});
    end
end

StimParams = ExpLoad.StimParams;
sim_stim = double(ExpLoad.simultaneous_stim);
E_MAP = ExpLoad.E_MAP;
n_Trials = double(ExpLoad.n_Trials);

if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.',nTrig,n_Trials);
elseif nTrig > n_Trials
    warning('%d triggers loaded for %d trials; extras are ignored.',nTrig,n_Trials);
    trig_ms = trig_ms(1:n_Trials);
end

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Stimulation events per trial (simultaneous_stim): %d\n',sim_stim);
fprintf('Trials used: %d\n',n_Trials);

%% =========================== DECODE AMPLITUDES ========================

amplitudes_all = cell2mat(StimParams(2:end,16));
amplitudes_all = double(amplitudes_all(:));

first_event_rows = 1:sim_stim:numel(amplitudes_all);
trialAmps = amplitudes_all(first_event_rows);
trialAmps = trialAmps(1:n_Trials);
trialAmps(trialAmps == -1) = 0;

if sim_stim >= 2
    second_event_rows = 2:sim_stim:numel(amplitudes_all);
    secondPulseAmps = amplitudes_all(second_event_rows);
    secondPulseAmps = secondPulseAmps(1:n_Trials);
    secondPulseAmps(secondPulseAmps == -1) = 0;

    if any(abs(trialAmps-secondPulseAmps) > 1e-6)
        warning(['Some trials contain different first- and second-pulse ' ...
            'amplitudes. Conditions use the first-pulse amplitude.']);
    end
end

[Amps,~,ampIdx] = unique(trialAmps(:));
nAMP = numel(Amps);

%% ============================== DECODE PTDs =============================

if sim_stim >= 2
    PTD_all_us = cell2mat(StimParams(2:end,6));
    PTD_all_us = double(PTD_all_us(:));

    second_event_rows = 2:sim_stim:numel(PTD_all_us);
    trialPTD_us = PTD_all_us(second_event_rows);
    trialPTD_us = trialPTD_us(1:n_Trials);

    [PTDs,~,ptdIdx] = unique(trialPTD_us(:));
    PTDs_ms = PTDs/1000;
    nPTD = numel(PTDs);
else
    PTDs = 0;
    PTDs_ms = 0;
    nPTD = 1;
    ptdIdx = ones(n_Trials,1);
end

%% ==================== DECODE STIMULATION CHANNEL ORDERS =================

stimNames = StimParams(2:end,1);
stimNames = stimNames(1:n_Trials*sim_stim);

[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));
if any(~isMapped)
    error('%d stimulation entries could not be mapped through E_MAP.',sum(~isMapped));
end

stimSeq = zeros(n_Trials,sim_stim);
for trial_id = 1:n_Trials
    rows_this_trial = (trial_id-1)*sim_stim+(1:sim_stim);
    mapped_channels = idx_all(rows_this_trial);
    mapped_channels = mapped_channels(mapped_channels > 0);
    stimSeq(trial_id,1:numel(mapped_channels)) = mapped_channels(:).';
end

[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

fprintf('\nAmplitudes: %s uA\n',num2str(Amps(:).'));
if sim_stim >= 2
    fprintf('PTDs: %s ms\n',num2str(PTDs_ms(:).'));
end
fprintf('Ordered stimulation sets: %d\n',nSets);
for si = 1:nSets
    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);
    fprintf('  Set %d: %s\n',si,format_order(stim_channels));
end

%% ======================= SELECT CHANNELS TO PLOT =======================

if isempty(Plot_Channels)
    plot_channels = 1:nMapPos;
else
    plot_channels = unique(double(Plot_Channels(:).'),'stable');
    valid_plot_channels = isfinite(plot_channels) & plot_channels >= 1 & ...
        plot_channels <= nMapPos & fix(plot_channels) == plot_channels;
    if any(~valid_plot_channels)
        warning('Invalid Plot_Channels were ignored: %s', ...
            num2str(plot_channels(~valid_plot_channels)));
        plot_channels = plot_channels(valid_plot_channels);
    end
end

plot_channels = sort(plot_channels,'ascend');   % plain channel-number order

if isempty(plot_channels)
    error('No valid Plot_Channels remain.');
end

fprintf('Plot channels: %d (%s ... %s)\n',numel(plot_channels), ...
    num2str(plot_channels(1)),num2str(plot_channels(end)));

channel_blocks = {};
for start_idx = 1:Max_Channels_Per_Figure:numel(plot_channels)
    end_idx = min(start_idx+Max_Channels_Per_Figure-1,numel(plot_channels));
    channel_blocks{end+1} = plot_channels(start_idx:end_idx); %#ok<SAGROW>
end
fprintf('Channel blocks per condition: %d\n',numel(channel_blocks));

%% ===================== CACHE SELECTED SPIKE TIMES =====================
% Map position -> physical channel -> spike-time vector

spike_time_cache = cell(nMapPos,1);
valid_spike_cache = false(nMapPos,1);

for channel_position = 1:numel(plot_channels)
    ich = plot_channels(channel_position);
    phys_ch = map_nums_plus(ich);

    valid_mapping = isfinite(phys_ch) && phys_ch >= 1 && phys_ch <= nCh && fix(phys_ch) == phys_ch;

    if ~valid_mapping || isempty(sp{phys_ch})
        continue;
    end

    spike_time_cache{ich} = sp{phys_ch};
    valid_spike_cache(ich) = true;
end

fprintf('Channels containing spike data: %d/%d\n', ...
    sum(valid_spike_cache(plot_channels)),numel(plot_channels));

clear sp;

%% ======================= SELECT CONDITIONS TO PLOT ======================

if isempty(Plot_Sets)
    selected_sets = 1:nSets;
else
    selected_sets = unique(double(Plot_Sets(:).'),'stable');
    valid_sets = selected_sets >= 1 & selected_sets <= nSets & fix(selected_sets) == selected_sets;
    if any(~valid_sets)
        warning('Invalid Plot_Sets were ignored.');
        selected_sets = selected_sets(valid_sets);
    end
end
if isempty(selected_sets)
    error('No valid stimulation sets were selected.');
end

if isempty(Plot_Amps)
    selected_amps = Amps(:).';
else
    selected_amps = double(Plot_Amps(:).');
end

if isempty(Plot_PTDs)
    selected_ptds = PTDs_ms(:).';
else
    selected_ptds = double(Plot_PTDs(:).');
end

fprintf('Selected sets: %s\n',number_list(selected_sets));
fprintf('Selected amplitudes: %s uA\n',number_list(selected_amps));
if sim_stim >= 2
    fprintf('Selected PTDs: %s ms\n',number_list(selected_ptds));
end

%% =================== PRECOMPUTE CONDITION TRIALS =====================

condition_trials = cell(nSets,nAMP,nPTD);
for si = 1:nSets
    for ai = 1:nAMP
        for pi = 1:nPTD
            condition_trials{si,ai,pi} = find(combClass == si & ampIdx == ai & ptdIdx == pi);
        end
    end
end

%% =========================== PSTH SETTINGS ============================

edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs = edges(1:end-1)+diff(edges)/2;
bin_s = bin_ms_raster/1000;

smooth_samples = max(1,round(smooth_ms/bin_ms_raster));
g = exp(-0.5 * ((0:smooth_samples-1)/(smooth_samples/2)).^2);
g = g/sum(g);

%% =====================================================================
% MAIN CONDITION LOOP
% ======================================================================

figures_created = 0;

for si = selected_sets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    % Physical distance between the stim channels of this set (NaN, with
    % a note, if undefined - e.g. channels on different probes)
    dist_label = '';
    if numel(stim_channels) == 2
        [dist_um,dist_note] = ChnPairDistance(stim_channels(1),stim_channels(2),Electrode_Type);
        if isnan(dist_um)
            dist_label = sprintf(' | Dist n/a (%s)',dist_note);
        else
            dist_label = sprintf(' | Dist %g um',dist_um);
        end
    end

    for amp_value = selected_amps

        ai = find(abs(Amps-amp_value) < Condition_Tolerance,1);
        if isempty(ai), continue; end

        for ptd_value = selected_ptds

            pi = find(abs(PTDs_ms-ptd_value) < Condition_Tolerance,1);
            if isempty(pi), continue; end

            trials_this = condition_trials{si,ai,pi};
            if isempty(trials_this), continue; end

            current_ptd_ms = PTDs_ms(pi);
            nConditionTrials = numel(trials_this);

            %% ---------------- CONDITION / STIM-TYPE LABEL --------------

            if sim_stim == 1
                stimulation_label = format_order(stim_channels);
                stimulation_mode = 'Single';
            elseif abs(current_ptd_ms) < Condition_Tolerance
                stimulation_label = format_simultaneous(stim_channels);
                stimulation_mode = 'Simultaneous';
            else
                stimulation_label = format_order(stim_channels);
                stimulation_mode = 'Sequential';
            end

            if sim_stim >= 2
                condition_core = sprintf('Set %d (%s) | Amp %g uA | PTD %g ms | nTrials=%d | %s%s', ...
                    si,stimulation_label,Amps(ai),current_ptd_ms,nConditionTrials,stimulation_mode,dist_label);
            else
                condition_core = sprintf('Set %d (%s) | Amp %g uA | nTrials=%d | %s%s', ...
                    si,stimulation_label,Amps(ai),nConditionTrials,stimulation_mode,dist_label);
            end

            %% ---------------- CHANNEL-BLOCK LOOP (<=64 PER FIGURE) ------

            for block_idx = 1:numel(channel_blocks)

                channels_this_block = channel_blocks{block_idx};
                nChPlot = numel(channels_this_block);

                if numel(channel_blocks) > 1
                    figTitle = sprintf('%s | Chn %d-%d', ...
                        condition_core,channels_this_block(1),channels_this_block(end));
                else
                    figTitle = condition_core;
                end

                %% ---------------- CREATE FIGURE -------------------------

                if Build_Figures_Invisibly
                    initial_visibility = 'off';
                else
                    initial_visibility = 'on';
                end

                if strcmpi(Figure_Window_Style,'docked')
                    fig = figure('Color','w','Name',figTitle,'NumberTitle','off', ...
                        'WindowStyle','docked','Visible',initial_visibility);
                else
                    fig = figure('Color','w','Name',figTitle,'NumberTitle','off', ...
                        'WindowStyle','normal','Position',fig_position,'Visible',initial_visibility);
                end

                layout = tiledlayout(fig,'flow','TileSpacing','compact','Padding','compact');
                title(layout,figTitle,'FontSize',13,'FontWeight','bold','Interpreter','none');

                figures_created = figures_created+1;

                %% ---------------- CHANNEL LOOP (PLAIN NUMBER ORDER) -----

                for channel_position = 1:nChPlot

                    ich = channels_this_block(channel_position);
                    ax = nexttile(layout);
                    hold(ax,'on');

                    if ~valid_spike_cache(ich)
                        title(ax,sprintf('Ch %d',ich),'FontSize',11,'FontWeight','bold');
                        axis(ax,'off');
                        continue;
                    end

                    spike_times = spike_time_cache{ich};

                    %% ------------ COLLECT ALL RASTER POINTS -------------

                    raster_x_cells = cell(nConditionTrials,1);
                    raster_y_cells = cell(nConditionTrials,1);

                    for trial_position = 1:nConditionTrials

                        trial_id = trials_this(trial_position);
                        trigger_time_ms = trig_ms(trial_id);

                        absolute_start = trigger_time_ms+ras_win(1);
                        absolute_end = trigger_time_ms+ras_win(2);

                        first_index = first_index_geq(spike_times,absolute_start);
                        end_index = first_index_geq(spike_times,absolute_end);

                        if first_index >= end_index, continue; end

                        relative_spikes = spike_times(first_index:end_index-1)-trigger_time_ms;
                        relative_spikes = relative_spikes(:);

                        raster_x_cells{trial_position} = relative_spikes;
                        raster_y_cells{trial_position} = repmat(trial_position,numel(relative_spikes),1);
                    end

                    if nConditionTrials == 0
                        all_raster_x = [];
                        all_raster_y = [];
                    else
                        all_raster_x = vertcat(raster_x_cells{:});
                        all_raster_y = vertcat(raster_y_cells{:});
                    end

                    %% ------------ PSTH -----------------------------------

                    if nConditionTrials == 0
                        rate_s = zeros(size(ctrs));
                    else
                        counts = histcounts(all_raster_x,edges);
                        rate = counts/(nConditionTrials*bin_s);
                        rate_s = filter(g,1,rate);
                    end

                    maxRate = max(rate_s);
                    yMaxPSTH = max(Minimum_PSTH_YMax,ceil(maxRate*1.1/10)*10);

                    %% ------------ LEFT Y-AXIS: PSTH -----------------------

                    yyaxis(ax,'left');
                    if any(rate_s)
                        plot(ax,ctrs,rate_s,'LineWidth',PSTH_Line_Width);
                    end
                    xlim(ax,ras_win);
                    ylim(ax,[0 yMaxPSTH]);
                    ylabel(ax,'Rate (sp/s)');

                    %% ------------ RIGHT Y-AXIS: RASTER --------------------

                    yyaxis(ax,'right');
                    if ~isempty(all_raster_x)
                        plot(ax,all_raster_x,all_raster_y,'.','LineStyle','none', ...
                            'Color',[0 0 0],'MarkerSize',Raster_Marker_Size);
                    end

                    if nConditionTrials > 0
                        ylim(ax,[0 nConditionTrials+1]);
                    else
                        ylim(ax,[0 1]);
                    end
                    set(ax,'YTick',[]);

                    %% ------------ STIMULATION MARKERS ---------------------

                    xline(ax,0,'r--','LineWidth',1);
                    if sim_stim >= 2 && current_ptd_ms > Condition_Tolerance
                        xline(ax,current_ptd_ms,'k:','LineWidth',1);
                    end
                    xlim(ax,ras_win);

                    %% ------------ CHANNEL TITLE ---------------------------

                    title(ax,sprintf('Ch %d',ich),'FontSize',11,'FontWeight','bold', ...
                        'Interpreter','none');

                    if channel_position > nChPlot-ceil(sqrt(nChPlot))
                        xlabel(ax,'Time (ms)');
                    end
                end

                %% ---------------- DISPLAY COMPLETED FIGURE ----------------

                if Build_Figures_Invisibly
                    fig.Visible = 'on';
                end
                drawnow limitrate;
            end
        end
    end
end

fprintf('\n============================================================\n');
fprintf('RASTER + PSTH PLOTTING COMPLETE\n');
fprintf('Figure tabs created: %d\n',figures_created);
fprintf('Figures saved: NO\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function index = first_index_geq(sorted_values,target)
nValues = numel(sorted_values);
low = 1;
high = nValues+1;
while low < high
    middle = floor((low+high)/2);
    if middle <= nValues && sorted_values(middle) < target
        low = middle+1;
    else
        high = middle;
    end
end
index = low;
end

function files = remove_metadata_and_backups(files)
if isempty(files), return; end
keep = true(size(files));
for file_index = 1:numel(files)
    file_name = files(file_index).name;
    if startsWith(file_name,'._') || contains(file_name,'BACKUP','IgnoreCase',true)
        keep(file_index) = false;
    end
end
files = files(keep);
end

function label = format_order(stim_channels)
stim_channels = stim_channels(stim_channels > 0);
if isempty(stim_channels)
    label = '(none)';
elseif numel(stim_channels) == 1
    label = sprintf('Ch%d',stim_channels(1));
else
    labels = arrayfun(@(c) sprintf('Ch%d',c),stim_channels,'UniformOutput',false);
    label = strjoin(labels,' -> ');
end
end

function label = format_simultaneous(stim_channels)
stim_channels = stim_channels(stim_channels > 0);
if isempty(stim_channels)
    label = '(none)';
elseif numel(stim_channels) == 1
    label = sprintf('Ch%d',stim_channels(1));
else
    labels = arrayfun(@(c) sprintf('Ch%d',c),stim_channels,'UniformOutput',false);
    label = strjoin(labels,' + ');
end
end

function output = number_list(values)
if isempty(values)
    output = '(none)';
else
    output = strtrim(sprintf('%g ',values));
end
end