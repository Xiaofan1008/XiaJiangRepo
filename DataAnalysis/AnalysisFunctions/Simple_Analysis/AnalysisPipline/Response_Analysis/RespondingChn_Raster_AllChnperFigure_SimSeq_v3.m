%% ========================================================================
% FAST MULTI-ISI RASTER + PSTH: ALL SELECTED CHANNELS
%
% SPEED IMPROVEMENTS
%   1. Figures are created as docked tabs.
%   2. Figures remain invisible until all tiles are complete.
%   3. Only spike-time columns from sp_corr are retained in memory.
%   4. Fast binary searches extract spikes from each trial window.
%   5. One raster graphics object is created per channel, rather than one
%      object for every trial.
%   6. PSTH counts are calculated once from all collected relative spikes.
%
% OUTPUT
%   One docked figure tab per:
%
%       selected set × amplitude × PTD
%
%   Figures are displayed but not saved.
% ========================================================================

clear;
% close all;

addpath(genpath( ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';

FS = 30000;

% Electrode type:
%   0 = rigid single-shank probe
%   1 = flexible single-shank probe
%   2 = four-shank flexible probe
Electrode_Type = 2;

%% ------------------------- CONDITION SELECTION ------------------------

% Channel Selectiomn
Plot_Channels = 1:64;
% Empty means all available sets
Plot_Sets = [1];
% Empty means all available amplitudes
Plot_Amps = [5];
% Empty means all available PTDs
% Example: Plot_PTDs = [0 3 5 10 20];
Plot_PTDs = [];

Condition_Tolerance = 1e-4;

%% ----------------------------- QC OPTIONS -----------------------------

% false = show all trials
% true  = load and exclude bad trials separately for each channel
Use_Bad_Trials = true;

% false = ignore bad-channel files
% true  = load and omit bad channels from each set
Use_Bad_Channels = false;

%% --------------------------- FIGURE OPTIONS ---------------------------

% 'docked':
%   Display every condition as a tab in one MATLAB figure container.
%
% 'normal':
%   Display every condition in a separate window.
Figure_Window_Style = 'docked';

% Only used when Figure_Window_Style = 'normal'
fig_position = [50 50 1600 900];

% Build each figure invisibly before displaying it
% This avoids repeated rendering while tiles are being added.
Build_Figures_Invisibly = true;

%% --------------------------- PLOTTING OPTIONS -------------------------

ras_win = [-50 80];

bin_ms_raster = 1;
smooth_ms = 5;

Raster_Marker_Size = 4;
PSTH_Line_Width = 1.4;
Minimum_PSTH_YMax = 50;

%% =========================== INITIAL CHECKS ===========================

if ~isfolder(data_folder)
    error('Dataset folder does not exist:\n%s',data_folder);
end

if ~isnumeric(Plot_Channels) || isempty(Plot_Channels)
    error('Plot_Channels must contain at least one channel index.');
end

if ~isnumeric(ras_win) || ...
        numel(ras_win) ~= 2 || ...
        any(~isfinite(ras_win)) || ...
        ras_win(2) <= ras_win(1)

    error('ras_win must be [start end], with end > start.');
end

if ~isscalar(bin_ms_raster) || ...
        ~isfinite(bin_ms_raster) || ...
        bin_ms_raster <= 0

    error('bin_ms_raster must be positive.');
end

if ~isscalar(smooth_ms) || ...
        ~isfinite(smooth_ms) || ...
        smooth_ms <= 0

    error('smooth_ms must be positive.');
end

if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('FS must be a positive sampling rate.');
end

valid_window_styles = {'docked','normal'};

if ~any(strcmpi(Figure_Window_Style,valid_window_styles))
    error('Figure_Window_Style must be ''docked'' or ''normal''.');
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

cd(data_folder);

fprintf('\n============================================================\n');
fprintf('FAST MULTI-ISI RASTER + PSTH\n');
fprintf('============================================================\n');
fprintf('Dataset: %s\n',data_folder);
fprintf('Figure style: %s\n',Figure_Window_Style);
fprintf('Raster window: [%g,%g) ms\n',ras_win);
fprintf('Use bad trials: %s\n',logical_text(Use_Bad_Trials));
fprintf('Use bad channels: %s\n',logical_text(Use_Bad_Channels));

%% ============================ LOAD sp_corr ============================

ssd_files = dir('*.sp_xia_SSD.mat');
ssd_files = remove_metadata_and_backups(ssd_files);

if isempty(ssd_files)
    error('No *.sp_xia_SSD.mat file was found.');
end

if numel(ssd_files) > 1
    error('Multiple current *.sp_xia_SSD.mat files were found.');
end

ssd_file = ssd_files(1).name;
base_name = erase(ssd_file,'.sp_xia_SSD.mat');

if ~ismember('sp_corr',who('-file',ssd_file))
    error('The SSD file does not contain sp_corr:\n%s',ssd_file);
end

SpikeLoad = load(ssd_file,'sp_corr');
sp = SpikeLoad.sp_corr;
nSpChannels = numel(sp);

if ~iscell(sp)
    error('sp_corr must be a cell array.');
end

fprintf('\nSpike file: %s\n',ssd_file);
fprintf('Spike-data cells: %d\n',nSpChannels);

%% ===================== LOAD RESPONDING CHANNELS =======================

Resp = [];
hasResp = false;

responding_file = sprintf( ...
    '%s_MultiISI_RespondingChannels.mat',base_name);

if isfile(responding_file)

    RespondingLoad = load(responding_file,'Responding');

    if isfield(RespondingLoad,'Responding')
        Resp = RespondingLoad.Responding;
        hasResp = true;

        fprintf('Responding file: %s\n',responding_file);
    else
        warning('Responding file does not contain Responding.');
    end

else
    warning(['Responding file not found. Channels will not be ' ...
        'highlighted:\n%s'],responding_file);
end

%% ==================== OPTIONALLY LOAD BAD CHANNELS ====================

BadCh_perSet = {};
BadCh_global = [];
bad_channel_file = '';

if Use_Bad_Channels

    bad_channel_patterns = { ...
        '*_MultiISIsBadChannels.mat', ...
        '*.MultiISIsBadChannels.mat', ...
        '*_MultiISIBadChannels.mat', ...
        '*.MultiISIBadChannels.mat', ...
        '*.BadChannels.mat'};

    bad_channel_file = find_first_qc_file( ...
        bad_channel_patterns);

    if isempty(bad_channel_file)

        warning(['Use_Bad_Channels is true, but no bad-channel ' ...
            'file was found.']);

    else
        BadChannelLoad = load(bad_channel_file);

        if isfield(BadChannelLoad,'BadCh_perSet')
            BadCh_perSet = BadChannelLoad.BadCh_perSet;

        elseif isfield(BadChannelLoad,'BadCh')
            BadCh_global = BadChannelLoad.BadCh;

        else
            warning(['Bad-channel file contains neither ' ...
                'BadCh_perSet nor BadCh.']);

            bad_channel_file = '';
        end
    end
end

if isempty(bad_channel_file)
    fprintf('Bad-channel file: not applied\n');
else
    fprintf('Bad-channel file: %s\n',bad_channel_file);
end

%% ===================== OPTIONALLY LOAD BAD TRIALS =====================

BadTrials = [];
bad_trial_file = '';

if Use_Bad_Trials

    bad_trial_patterns = { ...
        '*.MultiISIsBadTrials.mat', ...
        '*_MultiISIsBadTrials.mat', ...
        '*.MultiISIBadTrials.mat', ...
        '*_MultiISIBadTrials.mat', ...
        '*.SimSeqBadTrials.mat', ...
        '*_SimSeqBadTrials.mat', ...
        '*.BadTrials.mat'};

    bad_trial_file = find_first_qc_file( ...
        bad_trial_patterns);

    if isempty(bad_trial_file)

        warning(['Use_Bad_Trials is true, but no bad-trial ' ...
            'file was found.']);

    else
        BadTrialLoad = load(bad_trial_file);

        if isfield(BadTrialLoad,'BadTrials')
            BadTrials = BadTrialLoad.BadTrials;
        else
            warning('Bad-trial file does not contain BadTrials.');
            bad_trial_file = '';
        end
    end
end

if isempty(bad_trial_file)
    fprintf('Bad-trial file: not applied\n');
else
    fprintf('Bad-trial file: %s\n',bad_trial_file);
end

%% =========================== LOAD TRIGGERS ============================

if isempty(dir('*.trig.dat'))
    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);
end

trigger_files = dir('*.trig.dat');

if isempty(trigger_files)
    error('No *.trig.dat file could be found or created.');
end

trig = double(loadTrig(0));
trig = trig(:);
nTrig = numel(trig);

% Calculate trigger times once
trig_ms = trig/FS*1000;

fprintf('Triggers loaded: %d\n',nTrig);

%% ==================== LOAD EXPERIMENT PARAMETERS =====================

experiment_files = dir('*_exp_datafile_*.mat');
experiment_files = remove_metadata_and_backups(experiment_files);

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found.');
end

if numel(experiment_files) > 1
    error('Multiple current experiment files were found.');
end

experiment_file = experiment_files(1).name;

ExpLoad = load(experiment_file, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');

required_variables = { ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials'};

for variable_index = 1:numel(required_variables)

    variable_name = required_variables{variable_index};

    if ~isfield(ExpLoad,variable_name)
        error('Experiment file is missing variable: %s', ...
            variable_name);
    end
end

StimParams = ExpLoad.StimParams;
sim_stim = double(ExpLoad.simultaneous_stim);
E_MAP = ExpLoad.E_MAP;
n_Trials = double(ExpLoad.n_Trials);

if sim_stim ~= 2
    error(['This script requires paired stimulation, but ' ...
        'simultaneous_stim = %d.'],sim_stim);
end

if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.', ...
        nTrig,n_Trials);

elseif nTrig > n_Trials
    warning('%d triggers loaded for %d trials; extras are ignored.', ...
        nTrig,n_Trials);

    trig_ms = trig_ms(1:n_Trials);
end

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Trials used: %d\n',n_Trials);

%% =========================== DECODE AMPLITUDES ========================

amplitudes_all = cell2mat(StimParams(2:end,16));
amplitudes_all = double(amplitudes_all(:));

first_event_rows = 1:sim_stim:numel(amplitudes_all);
second_event_rows = 2:sim_stim:numel(amplitudes_all);

trialAmps = amplitudes_all(first_event_rows);
secondPulseAmps = amplitudes_all(second_event_rows);

trialAmps = trialAmps(1:n_Trials);
secondPulseAmps = secondPulseAmps(1:n_Trials);

trialAmps(trialAmps == -1) = 0;
secondPulseAmps(secondPulseAmps == -1) = 0;

if any(abs(trialAmps-secondPulseAmps) > 1e-6)
    warning(['Some trials contain different first- and second-pulse ' ...
        'amplitudes. Conditions use the first-pulse amplitude.']);
end

[Amps,~,ampIdx] = unique(trialAmps(:));
nAMP = numel(Amps);

%% ============================== DECODE PTDs ===========================

PTD_all_us = cell2mat(StimParams(2:end,6));
PTD_all_us = double(PTD_all_us(:));

trialPTD_us = PTD_all_us(second_event_rows);
trialPTD_us = trialPTD_us(1:n_Trials);

[PTDs,~,ptdIdx] = unique(trialPTD_us(:));
PTDs_ms = PTDs/1000;
nPTD = numel(PTDs);

%% ==================== DECODE STIMULATION ORDERS ======================

stimNames = StimParams(2:end,1);
stimNames = stimNames(1:n_Trials*sim_stim);

[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));

if any(~isMapped)
    error('%d stimulation entries could not be mapped.', ...
        sum(~isMapped));
end

stimSeq = zeros(n_Trials,sim_stim);

for trial_id = 1:n_Trials

    rows_this_trial = ...
        (trial_id-1)*sim_stim+(1:sim_stim);

    mapped_channels = idx_all(rows_this_trial);
    mapped_channels = mapped_channels(mapped_channels > 0);

    stimSeq(trial_id,1:numel(mapped_channels)) = ...
        mapped_channels(:).';
end

[uniqueComb,~,combClass] = unique( ...
    stimSeq,'rows','stable');

nSets = size(uniqueComb,1);

fprintf('\nAmplitudes: %s uA\n',num2str(Amps(:).'));
fprintf('PTDs: %s ms\n',num2str(PTDs_ms(:).'));
fprintf('Ordered sets: %d\n',nSets);

for si = 1:nSets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    fprintf('  Set %d: %s\n', ...
        si,format_order(stim_channels));
end

%% ========================= ELECTRODE MAPPING ==========================

d = double(Depth_s(Electrode_Type));
d = d(:);

if isempty(d)
    error('Depth_s returned an empty channel map.');
end

plot_channels = unique( ...
    double(Plot_Channels(:).'),'stable');

valid_plot_channels = ...
    isfinite(plot_channels) & ...
    plot_channels >= 1 & ...
    plot_channels <= numel(d) & ...
    fix(plot_channels) == plot_channels;

if any(~valid_plot_channels)

    warning('Invalid Plot_Channels were ignored: %s', ...
        num2str(plot_channels(~valid_plot_channels)));

    plot_channels = ...
        plot_channels(valid_plot_channels);
end

if isempty(plot_channels)
    error('No valid Plot_Channels remain.');
end

fprintf('Plot channels: %s\n',number_list(plot_channels));

%% ===================== CACHE SELECTED SPIKE TIMES =====================

% Retain only the first column of sp_corr for selected channels.
% This releases the much larger waveform matrices before plotting.

spike_time_cache = cell(numel(d),1);
valid_spike_cache = false(numel(d),1);

fprintf('\nCaching selected spike-time columns...\n');

for channel_position = 1:numel(plot_channels)

    ich = plot_channels(channel_position);
    spike_channel = d(ich);

    valid_mapping = ...
        isfinite(spike_channel) && ...
        spike_channel >= 1 && ...
        spike_channel <= nSpChannels && ...
        fix(spike_channel) == spike_channel;

    if ~valid_mapping || isempty(sp{spike_channel})
        continue;
    end

    spike_times = double(sp{spike_channel}(:,1));
    spike_times = spike_times(isfinite(spike_times));

    % Binary window searches require sorted spike times
    if ~issorted(spike_times)
        spike_times = sort(spike_times);
    end

    spike_time_cache{ich} = spike_times(:);
    valid_spike_cache(ich) = true;
end

fprintf('Channels containing spike data: %d/%d\n', ...
    sum(valid_spike_cache(plot_channels)), ...
    numel(plot_channels));

% Release waveform data before creating figures
clear sp SpikeLoad;

%% ======================== SELECT CONDITIONS ==========================

if isempty(Plot_Sets)

    selected_sets = 1:nSets;

else
    selected_sets = unique( ...
        double(Plot_Sets(:).'),'stable');

    valid_sets = ...
        selected_sets >= 1 & ...
        selected_sets <= nSets & ...
        fix(selected_sets) == selected_sets;

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
fprintf('Selected PTDs: %s ms\n',number_list(selected_ptds));

%% =================== PRECOMPUTE CONDITION TRIALS =====================

condition_trials = cell(nSets,nAMP,nPTD);

for si = 1:nSets
    for ai = 1:nAMP
        for pi = 1:nPTD

            condition_trials{si,ai,pi} = find( ...
                combClass == si & ...
                ampIdx == ai & ...
                ptdIdx == pi);
        end
    end
end

%% =========================== PSTH SETTINGS ============================

edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs = edges(1:end-1)+diff(edges)/2;
bin_s = bin_ms_raster/1000;

smooth_samples = max( ...
    1,round(smooth_ms/bin_ms_raster));

g = exp(-0.5 * ...
    ((0:smooth_samples-1)/(smooth_samples/2)).^2);

g = g/sum(g);

%% ======================= COUNT EXPECTED FIGURES =======================

expected_figures = 0;

for si = selected_sets
    for amp_value = selected_amps
        for ptd_value = selected_ptds

            ai = find(abs(Amps-amp_value) < ...
                Condition_Tolerance,1);

            pi = find(abs(PTDs_ms-ptd_value) < ...
                Condition_Tolerance,1);

            if isempty(ai) || isempty(pi)
                continue;
            end

            if ~isempty(condition_trials{si,ai,pi})
                expected_figures = expected_figures+1;
            end
        end
    end
end

fprintf('Expected figure tabs: %d\n',expected_figures);

%% =====================================================================
% MAIN CONDITION LOOPS
% ======================================================================

figures_created = 0;

for si = selected_sets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    %% ---------------- BAD CHANNELS FOR THIS SET -----------------------

    if Use_Bad_Channels && ~isempty(bad_channel_file)

        bad_channels_this_set = ...
            get_bad_channels_for_set( ...
            BadCh_perSet,BadCh_global,si);

    else
        bad_channels_this_set = [];
    end

    channels_this_set = plot_channels;

    if Use_Bad_Channels

        channels_this_set = setdiff( ...
            channels_this_set, ...
            bad_channels_this_set, ...
            'stable');
    end

    if isempty(channels_this_set)

        warning(['No channels remain for Set %d after applying ' ...
            'bad-channel exclusion.'],si);

        continue;
    end

    for amp_value = selected_amps

        ai = find(abs(Amps-amp_value) < ...
            Condition_Tolerance,1);

        if isempty(ai)
            continue;
        end

        for ptd_value = selected_ptds

            pi = find(abs(PTDs_ms-ptd_value) < ...
                Condition_Tolerance,1);

            if isempty(pi)
                continue;
            end

            trials_this = condition_trials{si,ai,pi};

            if isempty(trials_this)
                continue;
            end

            current_ptd_ms = PTDs_ms(pi);
            nConditionTrials = numel(trials_this);

            %% ---------------- CONDITION TITLE -------------------------

            if abs(current_ptd_ms) < Condition_Tolerance

                stimulation_label = ...
                    format_simultaneous(stim_channels);

                stimulation_mode = 'Simultaneous';

            else
                stimulation_label = ...
                    format_order(stim_channels);

                stimulation_mode = 'Sequential';
            end

            figTitle = sprintf( ...
                ['Set %d (%s) | Amp %g uA | PTD %g ms | ' ...
                'nTrials=%d | %s'], ...
                si, ...
                stimulation_label, ...
                Amps(ai), ...
                current_ptd_ms, ...
                nConditionTrials, ...
                stimulation_mode);

            %% ---------------- CREATE DOCKED FIGURE ---------------------

            if Build_Figures_Invisibly
                initial_visibility = 'off';
            else
                initial_visibility = 'on';
            end

            if strcmpi(Figure_Window_Style,'docked')

                % Do not specify Position for a docked figure.
                % Setting Position would undock it.
                fig = figure( ...
                    'Color','w', ...
                    'Name',figTitle, ...
                    'NumberTitle','off', ...
                    'WindowStyle','docked', ...
                    'Visible',initial_visibility);

            else
                fig = figure( ...
                    'Color','w', ...
                    'Name',figTitle, ...
                    'NumberTitle','off', ...
                    'WindowStyle','normal', ...
                    'Position',fig_position, ...
                    'Visible',initial_visibility);
            end

            layout = tiledlayout( ...
                fig, ...
                'flow', ...
                'TileSpacing','compact', ...
                'Padding','compact');

            title( ...
                layout, ...
                figTitle, ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'Interpreter','none');

            figures_created = figures_created+1;
            nChPlot = numel(channels_this_set);

            %% ---------------- CHANNEL LOOP ----------------------------

            for channel_position = 1:nChPlot

                ich = channels_this_set(channel_position);

                ax = nexttile(layout);
                hold(ax,'on');

                if ~valid_spike_cache(ich)

                    title(ax,sprintf('Ch %d',ich), ...
                        'FontSize',11, ...
                        'FontWeight','bold');

                    axis(ax,'off');
                    continue;
                end

                spike_times = spike_time_cache{ich};

                %% ------------ CHANNEL-SPECIFIC VALID TRIALS -----------

                if Use_Bad_Trials && ~isempty(bad_trial_file)

                    bad_trials_this_channel = ...
                        get_bad_trials_for_channel( ...
                        BadTrials,ich);

                else
                    bad_trials_this_channel = [];
                end

                valid_trials = setdiff( ...
                    trials_this, ...
                    bad_trials_this_channel, ...
                    'stable');

                nValidTrials = numel(valid_trials);

                %% ------------ COLLECT ALL RASTER POINTS ---------------

                raster_x_cells = cell(nValidTrials,1);
                raster_y_cells = cell(nValidTrials,1);

                for trial_position = 1:nValidTrials

                    trial_id = valid_trials(trial_position);
                    trigger_time_ms = trig_ms(trial_id);

                    absolute_start = ...
                        trigger_time_ms+ras_win(1);

                    absolute_end = ...
                        trigger_time_ms+ras_win(2);

                    % Use binary searches rather than scanning the entire
                    % spike vector with a logical comparison.
                    first_index = first_index_geq( ...
                        spike_times,absolute_start);

                    end_index = first_index_geq( ...
                        spike_times,absolute_end);

                    if first_index >= end_index
                        continue;
                    end

                    relative_spikes = ...
                        spike_times(first_index:end_index-1) ...
                        -trigger_time_ms;

                    relative_spikes = relative_spikes(:);

                    raster_x_cells{trial_position} = ...
                        relative_spikes;

                    raster_y_cells{trial_position} = ...
                        repmat( ...
                        trial_position, ...
                        numel(relative_spikes), ...
                        1);
                end

                if nValidTrials == 0

                    all_raster_x = [];
                    all_raster_y = [];

                else
                    all_raster_x = ...
                        vertcat(raster_x_cells{:});

                    all_raster_y = ...
                        vertcat(raster_y_cells{:});
                end

                %% ------------ PSTH CALCULATION ------------------------

                if nValidTrials == 0

                    rate_s = zeros(size(ctrs));

                else
                    counts = histcounts(all_raster_x,edges);
                    rate = counts/(nValidTrials*bin_s);
                    rate_s = filter(g,1,rate);
                end

                maxRate = max(rate_s);

                yMaxPSTH = max( ...
                    Minimum_PSTH_YMax, ...
                    ceil(maxRate*1.1/10)*10);

                %% ------------ RESPONDING LABEL ------------------------

                isResp = false;

                if hasResp && ...
                        si <= numel(Resp.set) && ...
                        ai <= numel(Resp.set(si).amp) && ...
                        isfield(Resp.set(si).amp(ai),'ptd') && ...
                        pi <= numel(Resp.set(si).amp(ai).ptd) && ...
                        isfield( ...
                        Resp.set(si).amp(ai).ptd(pi),'channel') && ...
                        ich <= numel( ...
                        Resp.set(si).amp(ai).ptd(pi).channel)

                    R = Resp.set(si).amp(ai).ptd(pi).channel(ich);

                    if isfield(R,'is_responsive') && ...
                            R.is_responsive

                        isResp = true;
                    end
                end

                %% ------------ CHANNEL APPEARANCE ----------------------

                if isResp
                    ax.Color = [1.0 0.88 0.88];
                    ax.LineWidth = 1.4;
                end

                %% ------------ LEFT Y-AXIS: PSTH -----------------------

                yyaxis(ax,'left');

                if any(rate_s)

                    plot(ax,ctrs,rate_s, ...
                        'LineWidth',PSTH_Line_Width);
                end

                xlim(ax,ras_win);
                ylim(ax,[0 yMaxPSTH]);
                ylabel(ax,'Rate (sp/s)');

                %% ------------ RIGHT Y-AXIS: RASTER --------------------

                yyaxis(ax,'right');

                % One graphics object contains all raster points for this
                % channel. This is much faster than one plot per trial.
                if ~isempty(all_raster_x)

                    plot(ax, ...
                        all_raster_x, ...
                        all_raster_y, ...
                        '.', ...
                        'LineStyle','none', ...
                        'Color',[0 0 0], ...
                        'MarkerSize',Raster_Marker_Size);
                end

                if nValidTrials > 0
                    ylim(ax,[0 nValidTrials+1]);
                else
                    ylim(ax,[0 1]);
                end

                % Retain the clean presentation setting
                set(ax,'YTick',[]);

                %% ------------ STIMULATION MARKERS ---------------------

                xline(ax,0, ...
                    'r--', ...
                    'LineWidth',1);

                if current_ptd_ms > Condition_Tolerance

                    xline(ax,current_ptd_ms, ...
                        'k:', ...
                        'LineWidth',1);
                end

                xlim(ax,ras_win);

                %% ------------ CHANNEL TITLE ---------------------------

                if isResp
                    channel_title = sprintf('Ch %d RESP',ich);
                else
                    channel_title = sprintf('Ch %d',ich);
                end

                title_handle = title( ...
                    ax, ...
                    channel_title, ...
                    'FontSize',11, ...
                    'FontWeight','bold', ...
                    'Interpreter','none');

                if isResp
                    title_handle.Color = [0.7 0 0];
                end

                if channel_position > ...
                        nChPlot-ceil(sqrt(nChPlot))

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

fprintf('\n============================================================\n');
fprintf('FAST RASTER + PSTH PLOTTING COMPLETE\n');
fprintf('Figure tabs created: %d\n',figures_created);
fprintf('Figures saved: NO\n');
fprintf('Original experiment files modified: NO\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function index = first_index_geq(sorted_values,target)
% Return the index of the first value greater than or equal to target.
%
% If every value is smaller than target, this returns:
%       numel(sorted_values)+1
%
% This binary search avoids repeatedly scanning a full spike vector.

nValues = numel(sorted_values);

low = 1;
high = nValues+1;

while low < high

    middle = floor((low+high)/2);

    if middle <= nValues && ...
            sorted_values(middle) < target

        low = middle+1;
    else
        high = middle;
    end
end

index = low;
end

function files = remove_metadata_and_backups(files)
% Remove macOS metadata files and timestamped backup files.

if isempty(files)
    return;
end

keep = true(size(files));

for file_index = 1:numel(files)

    file_name = files(file_index).name;

    if startsWith(file_name,'._') || ...
            contains(file_name,'BACKUP','IgnoreCase',true)

        keep(file_index) = false;
    end
end

files = files(keep);
end

function file_name = find_first_qc_file(patterns)
% Find the first QC file while ignoring metadata and backup files.

file_name = '';

for pattern_index = 1:numel(patterns)

    candidates = dir(patterns{pattern_index});
    candidates = remove_metadata_and_backups(candidates);

    if ~isempty(candidates)
        file_name = candidates(1).name;
        return;
    end
end
end

function bad_channels = get_bad_channels_for_set( ...
    BadCh_perSet,BadCh_global,set_index)
% Return bad displayed-channel indices for one stimulation set.

bad_channels = [];

if ~isempty(BadCh_perSet)

    if iscell(BadCh_perSet)

        if set_index <= numel(BadCh_perSet)
            bad_channels = BadCh_perSet{set_index};
        end

    elseif isnumeric(BadCh_perSet)
        bad_channels = BadCh_perSet;
    end

elseif ~isempty(BadCh_global)

    if iscell(BadCh_global)

        if set_index <= numel(BadCh_global)
            bad_channels = BadCh_global{set_index};
        end

    elseif isnumeric(BadCh_global)
        bad_channels = BadCh_global;
    end
end

if isempty(bad_channels)
    bad_channels = [];
else
    bad_channels = unique( ...
        double(bad_channels(:).'),'stable');

    bad_channels = bad_channels(isfinite(bad_channels));
end
end

function bad_trials = get_bad_trials_for_channel( ...
    BadTrials,channel_index)
% Support numeric global lists or channel-specific cell arrays.

bad_trials = [];

if isempty(BadTrials)
    return;
end

if isnumeric(BadTrials)

    bad_trials = BadTrials;

elseif iscell(BadTrials) && ...
        channel_index <= numel(BadTrials)

    bad_trials = BadTrials{channel_index};
end

if isempty(bad_trials)
    bad_trials = [];
else
    bad_trials = unique( ...
        double(bad_trials(:).'),'stable');

    bad_trials = bad_trials( ...
        isfinite(bad_trials) & ...
        bad_trials >= 1 & ...
        fix(bad_trials) == bad_trials);
end
end

function label = format_order(stim_channels)
% Format sequential stimulation as ChA -> ChB.

stim_channels = stim_channels(stim_channels > 0);

if isempty(stim_channels)

    label = '(none)';

elseif numel(stim_channels) == 1

    label = sprintf('Ch%d',stim_channels(1));

else
    labels = arrayfun(@(channel_number) ...
        sprintf('Ch%d',channel_number), ...
        stim_channels, ...
        'UniformOutput',false);

    label = strjoin(labels,' -> ');
end
end

function label = format_simultaneous(stim_channels)
% Format simultaneous stimulation as ChA + ChB.

stim_channels = stim_channels(stim_channels > 0);

if isempty(stim_channels)

    label = '(none)';

elseif numel(stim_channels) == 1

    label = sprintf('Ch%d',stim_channels(1));

else
    labels = arrayfun(@(channel_number) ...
        sprintf('Ch%d',channel_number), ...
        stim_channels, ...
        'UniformOutput',false);

    label = strjoin(labels,' + ');
end
end

function output = number_list(values)
% Convert a numeric vector into command-window text.

if isempty(values)
    output = '(none)';
else
    output = strtrim(sprintf('%g ',values));
end
end

function output = logical_text(value)
% Convert a logical value to ON or OFF.

if value
    output = 'ON';
else
    output = 'OFF';
end
end