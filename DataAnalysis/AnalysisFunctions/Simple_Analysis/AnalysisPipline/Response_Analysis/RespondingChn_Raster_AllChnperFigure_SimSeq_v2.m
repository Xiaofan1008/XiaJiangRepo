%% ========================================================================
% MULTI-ISI RASTER + PSTH: ALL SELECTED CHANNELS
%
% PURPOSE
%   Plot all selected recording channels in the same figure for every
%   selected:
%
%       ordered stimulation set × amplitude × PTD
%
% SPIKE DATA
%   Strictly requires:
%
%       *.sp_xia_SSD.mat
%       variable: sp_corr
%
% RESPONDING CHANNELS
%   Loads:
%
%       <base_name>_MultiISI_RespondingChannels.mat
%
%   Responding channels are highlighted with a pale-red background.
%
% OPTIONAL QC FILES
%   Use_Bad_Trials:
%       false = show every trial
%       true  = load bad trials and exclude them from each channel
%
%   Use_Bad_Channels:
%       false = ignore bad-channel files
%       true  = load bad channels and omit them from the figure
%
% CHANNEL SELECTION
%   Plot_Channels is entered manually using displayed channel indices.
%
% EXAMPLES
%       Plot_Channels = 1:64;
%       Plot_Channels = 35:64;
%       Plot_Channels = [35:40 42:48 50:64];
%
% FIGURES
%   One figure is produced for each selected condition.
%   Figures are displayed but are not saved.
% ========================================================================

clear;
close all;

addpath(genpath( ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================
data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';
% Sampling rate in Hz
FS = 30000;

%% -------------------------- CONDITION SELECTION -----------------------
% Electrode type:
%   0 = rigid single-shank probe
%   1 = flexible single-shank probe
%   2 = four-shank flexible probe
Electrode_Type = 2;

% Channel selection
Plot_Channels = 1:64;

% Empty means plot every available set.
% Example: Plot_Sets = [1 2];
Plot_Sets = [];

% Empty means plot every available amplitude.
% Example: Plot_Amps = [5 10];
Plot_Amps = [];

% Empty means plot every available PTD.
% PTD is entered in milliseconds.
% Example: Plot_PTDs = [0 5 10 20];
Plot_PTDs = [];

% Numerical tolerance for amplitude and PTD matching
Condition_Tolerance = 1e-4;

%% ----------------------------- QC OPTIONS -----------------------------

% false:
%   Do not load or apply bad-trial information.
%
% true:
%   Load a bad-trial file and exclude the bad trials separately for each
%   recording channel.
Use_Bad_Trials = false;

% false:
%   Do not load or apply bad-channel information.
%
% true:
%   Load BadCh_perSet or BadCh and omit those channels from each figure.
Use_Bad_Channels = false;

%% --------------------------- PLOTTING OPTIONS -------------------------

% Raster window relative to the first pulse
ras_win = [-50 80];

% PSTH bin width in milliseconds
bin_ms_raster = 1;

% Existing one-sided smoothing-kernel duration
smooth_ms = 5;

% Figure position and size
fig_position = [50 50 1600 900];

% Raster marker size
Raster_Marker_Size = 4;

% PSTH line width
PSTH_Line_Width = 1.4;

% Minimum PSTH y-axis limit
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

    error('bin_ms_raster must be a positive number.');
end

if ~isscalar(smooth_ms) || ...
        ~isfinite(smooth_ms) || ...
        smooth_ms <= 0

    error('smooth_ms must be a positive number.');
end

if ~isscalar(FS) || ~isfinite(FS) || FS <= 0
    error('FS must be a positive sampling rate in Hz.');
end

% Remember and restore the original MATLAB folder
starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

% Move automatically to the selected dataset
cd(data_folder);

fprintf('\n============================================================\n');
fprintf('MULTI-ISI RASTER + PSTH\n');
fprintf('============================================================\n');
fprintf('Dataset: %s\n',data_folder);
fprintf('Electrode type: %d\n',Electrode_Type);
fprintf('Sampling rate: %g Hz\n',FS);
fprintf('Raster window: [%g,%g) ms\n',ras_win);
fprintf('Use bad trials: %s\n',logical_text(Use_Bad_Trials));
fprintf('Use bad channels: %s\n',logical_text(Use_Bad_Channels));

%% ============================ LOAD sp_corr ============================

ssd_files = dir('*.sp_xia_SSD.mat');

if isempty(ssd_files)
    error('No *.sp_xia_SSD.mat file was found.');
end

if numel(ssd_files) > 1
    warning('Multiple SSD files found. Using: %s',ssd_files(1).name);
end

ssd_file = ssd_files(1).name;
base_name = erase(ssd_file,'.sp_xia_SSD.mat');

if ~ismember('sp_corr',who('-file',ssd_file))
    error(['The SSD file does not contain sp_corr.\n' ...
        'Required file: %s'],ssd_file);
end

SpikeLoad = load(ssd_file,'sp_corr');
sp = SpikeLoad.sp_corr;
nSpChannels = numel(sp);

if ~iscell(sp)
    error('sp_corr must be a cell array containing one cell per channel.');
end

fprintf('\nSpike file: %s\n',ssd_file);
fprintf('Spike variable: sp_corr\n');
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

        fprintf('Responding-channel file: %s\n',responding_file);
    else
        warning(['The responding-channel file exists but does not ' ...
            'contain Responding.']);
    end

else
    warning(['No MultiISI responding-channel file was found.\n' ...
        'Expected file: %s\n' ...
        'Channels will not be highlighted as responding.'], ...
        responding_file);
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

        warning(['Use_Bad_Channels is true, but no suitable ' ...
            'bad-channel file was found.']);

    else
        BadChannelLoad = load(bad_channel_file);

        if isfield(BadChannelLoad,'BadCh_perSet')
            BadCh_perSet = BadChannelLoad.BadCh_perSet;

        elseif isfield(BadChannelLoad,'BadCh')
            BadCh_global = BadChannelLoad.BadCh;

        else
            warning(['Bad-channel file does not contain BadCh_perSet ' ...
                'or BadCh. Bad channels will not be applied.']);

            bad_channel_file = '';
        end
    end
end

if Use_Bad_Channels && ~isempty(bad_channel_file)
    fprintf('Bad-channel file: %s\n',bad_channel_file);

elseif ~Use_Bad_Channels
    fprintf('Bad-channel file: not used\n');
end

%% ===================== OPTIONALLY LOAD BAD TRIALS =====================

BadTrials = [];
bad_trial_file = '';

if Use_Bad_Trials

    bad_trial_patterns = { ...
        '*_MultiISIsBadTrials.mat', ...
        '*.MultiISIsBadTrials.mat', ...
        '*_MultiISIBadTrials.mat', ...
        '*.MultiISIBadTrials.mat', ...
        '*.SimSeqBadTrials.mat', ...
        '*_SimSeqBadTrials.mat', ...
        '*.BadTrials.mat'};

    bad_trial_file = find_first_qc_file( ...
        bad_trial_patterns);

    if isempty(bad_trial_file)

        warning(['Use_Bad_Trials is true, but no suitable ' ...
            'bad-trial file was found.']);

    else
        BadTrialLoad = load(bad_trial_file);

        if isfield(BadTrialLoad,'BadTrials')
            BadTrials = BadTrialLoad.BadTrials;
        else
            warning(['Bad-trial file does not contain BadTrials. ' ...
                'Bad trials will not be applied.']);

            bad_trial_file = '';
        end
    end
end

if Use_Bad_Trials && ~isempty(bad_trial_file)
    fprintf('Bad-trial file: %s\n',bad_trial_file);

elseif ~Use_Bad_Trials
    fprintf('Bad-trial file: not used\n');
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

fprintf('Trigger file: %s\n',trigger_files(1).name);
fprintf('Triggers loaded: %d\n',nTrig);

%% ==================== LOAD EXPERIMENT PARAMETERS =====================

experiment_files = dir('*_exp_datafile_*.mat');

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found.');
end

if numel(experiment_files) > 1
    warning('Multiple experiment files found. Using: %s', ...
        experiment_files(1).name);
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
sim_stim   = double(ExpLoad.simultaneous_stim);
E_MAP      = ExpLoad.E_MAP;
n_Trials   = double(ExpLoad.n_Trials);

if sim_stim ~= 2
    error(['This script requires two stimulation events per trial, but ' ...
        'simultaneous_stim = %d.'],sim_stim);
end

if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.', ...
        nTrig,n_Trials);

elseif nTrig > n_Trials
    warning('%d triggers were loaded for %d trials; extras are ignored.', ...
        nTrig,n_Trials);

    trig = trig(1:n_Trials);
end

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Trials/triggers used: %d/%d\n',n_Trials,n_Trials);

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

amplitude_mismatch = ...
    abs(trialAmps-secondPulseAmps) > 1e-6;

if any(amplitude_mismatch)
    warning(['The two pulses have different amplitudes in %d trials. ' ...
        'Conditions will use the first-pulse amplitude.'], ...
        sum(amplitude_mismatch));
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
    error('%d stimulation entries could not be mapped through E_MAP.', ...
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

% Preserve stimulation order
[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

fprintf('\nDetected amplitudes: %s uA\n', ...
    num2str(Amps(:).'));

fprintf('Detected PTDs: %s ms\n', ...
    num2str(PTDs_ms(:).'));

fprintf('Detected ordered stimulation sets: %d\n',nSets);

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

% Plot_Channels is manually selected by the user.
plot_channels = unique(double(Plot_Channels(:).'),'stable');

valid_plot_channels = ...
    isfinite(plot_channels) & ...
    plot_channels >= 1 & ...
    plot_channels <= numel(d) & ...
    fix(plot_channels) == plot_channels;

if any(~valid_plot_channels)

    warning(['The following Plot_Channels are invalid and will be ' ...
        'ignored: %s'], ...
        num2str(plot_channels(~valid_plot_channels)));

    plot_channels = plot_channels(valid_plot_channels);
end

if isempty(plot_channels)
    error('No valid Plot_Channels remain after channel validation.');
end

fprintf('\nUser-selected plot channels: %s\n', ...
    compact_number_list(plot_channels));

%% ======================== SELECT CONDITIONS ==========================

if isempty(Plot_Sets)
    selected_sets = 1:nSets;
else
    selected_sets = unique(double(Plot_Sets(:).'),'stable');

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

fprintf('Selected sets: %s\n',num2str(selected_sets));
fprintf('Selected amplitudes: %s uA\n',num2str(selected_amps));
fprintf('Selected PTDs: %s ms\n',num2str(selected_ptds));

%% =========================== PSTH SETTINGS ============================

edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs = edges(1:end-1)+diff(edges)/2;
bin_s = bin_ms_raster/1000;

% Retain the existing smoothing method
smooth_samples = max(1,round(smooth_ms/bin_ms_raster));

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

            trials_this = find( ...
                combClass == si & ...
                ampIdx == ai & ...
                ptdIdx == pi);

            if ~isempty(trials_this)
                expected_figures = expected_figures+1;
            end
        end
    end
end

fprintf('Expected figures: %d\n',expected_figures);

%% =====================================================================
% MAIN CONDITION LOOPS
% ======================================================================

figures_created = 0;

for si = selected_sets

    stim_channels = uniqueComb(si,uniqueComb(si,:) > 0);

    % Read the bad channels for this stimulation set
    if Use_Bad_Channels && ~isempty(bad_channel_file)

        bad_channels_this_set = get_bad_channels_for_set( ...
            BadCh_perSet,BadCh_global,si);

    else
        bad_channels_this_set = [];
    end

    % If bad-channel exclusion is enabled, omit those channel tiles
    channels_this_set = plot_channels;

    if Use_Bad_Channels
        channels_this_set = setdiff( ...
            channels_this_set, ...
            bad_channels_this_set, ...
            'stable');
    end

    if isempty(channels_this_set)

        warning(['No plot channels remain for Set %d after applying ' ...
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

            current_ptd_ms = PTDs_ms(pi);

            trials_this = find( ...
                combClass == si & ...
                ampIdx == ai & ...
                ptdIdx == pi);

            if isempty(trials_this)
                continue;
            end

            nConditionTrials = numel(trials_this);

            %% ---------------- CONDITION LABEL -------------------------

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

            %% ---------------- CREATE FIGURE ---------------------------

            figure( ...
                'Color','w', ...
                'Name',figTitle, ...
                'Position',fig_position);

            tiledlayout( ...
                'flow', ...
                'TileSpacing','compact', ...
                'Padding','compact');

            sgtitle( ...
                figTitle, ...
                'FontSize',14, ...
                'FontWeight','bold', ...
                'Interpreter','none');

            figures_created = figures_created+1;

            %% ---------------- CHANNEL LOOP ----------------------------

            nChPlot = numel(channels_this_set);

            for channel_position = 1:nChPlot

                ich = channels_this_set(channel_position);
                spike_channel = d(ich);

                ax = nexttile;
                hold(ax,'on');

                %% ------------ CHECK CHANNEL MAPPING -------------------

                valid_spike_channel = ...
                    isfinite(spike_channel) && ...
                    spike_channel >= 1 && ...
                    spike_channel <= nSpChannels && ...
                    fix(spike_channel) == spike_channel;

                if ~valid_spike_channel || isempty(sp{spike_channel})

                    title(ax,sprintf('Ch %d',ich), ...
                        'FontSize',11, ...
                        'FontWeight','bold');

                    axis(ax,'off');
                    continue;
                end

                spike_times = double(sp{spike_channel}(:,1));

                %% ------------ APPLY CHANNEL-SPECIFIC BAD TRIALS -------

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

                %% ------------ COLLECT RASTER AND PSTH -----------------

                allTrialSpikes = cell(nValidTrials,1);
                counts = zeros(1,numel(edges)-1);

                for trial_position = 1:nValidTrials

                    trial_id = valid_trials(trial_position);
                    trigger_ms = trig(trial_id)/FS*1000;

                    relative_spikes = spike_times-trigger_ms;

                    relative_spikes = relative_spikes( ...
                        relative_spikes >= ras_win(1) & ...
                        relative_spikes <  ras_win(2));

                    allTrialSpikes{trial_position} = ...
                        relative_spikes;

                    counts = counts+histcounts( ...
                        relative_spikes,edges);
                end

                if nValidTrials == 0

                    rate_s = zeros(size(ctrs));

                else
                    rate = counts/(nValidTrials*bin_s);
                    rate_s = filter(g,1,rate);
                end

                maxRate = max(rate_s);

                yMaxPSTH = max( ...
                    Minimum_PSTH_YMax, ...
                    ceil(maxRate*1.1/10)*10);

                %% ------------ CHECK RESPONDING LABEL -----------------

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

                %% ------------ RESPONDING-CHANNEL APPEARANCE ----------

                if isResp
                    set(ax,'Color',[1.0 0.88 0.88]);
                    ax.XColor = [0.6 0 0];
                    ax.YColor = [0.6 0 0];
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

                for trial_position = 1:nValidTrials

                    relative_spikes = ...
                        allTrialSpikes{trial_position};

                    if isempty(relative_spikes)
                        continue;
                    end

                    plot(ax, ...
                        relative_spikes, ...
                        trial_position*ones(size(relative_spikes)), ...
                        '.', ...
                        'Color',[0 0 0], ...
                        'MarkerSize',Raster_Marker_Size);
                end

                if nValidTrials > 0
                    ylim(ax,[0 nValidTrials+1]);
                else
                    ylim(ax,[0 1]);
                end

                % Retain the existing presentation setting:
                % do not display raster trial-number ticks.
                set(ax,'YTick',[]);

                %% ------------ STIMULATION-TIME MARKERS ----------------

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
                    set(title_handle,'Color',[0.7 0 0]);
                end

                %% ------------ X-AXIS LABELS ---------------------------

                if channel_position > ...
                        nChPlot-ceil(sqrt(nChPlot))

                    xlabel(ax,'Time (ms)');
                end
            end
        end
    end
end

fprintf('\n============================================================\n');
fprintf('RASTER + PSTH PLOTTING COMPLETE\n');
fprintf('Figures created: %d\n',figures_created);
fprintf('Figures saved: NO\n');
fprintf('Original experiment files modified: NO\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function file_name = find_first_qc_file(patterns)
% Find the first suitable QC file while ignoring backups and macOS
% metadata files beginning with ._.

file_name = '';

for pattern_index = 1:numel(patterns)

    candidates = dir(patterns{pattern_index});

    for candidate_index = 1:numel(candidates)

        candidate_name = candidates(candidate_index).name;

        if startsWith(candidate_name,'._')
            continue;
        end

        if contains(candidate_name,'BACKUP', ...
                'IgnoreCase',true)
            continue;
        end

        file_name = candidate_name;
        return;
    end
end
end

function bad_channels = get_bad_channels_for_set( ...
    BadCh_perSet,BadCh_global,set_index)
% Return displayed bad-channel indices for one stimulation set.

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
% Support either:
%   1. One numeric global bad-trial list
%   2. A cell array containing one bad-trial list per channel

bad_trials = [];

if isempty(BadTrials)
    return;
end

if isnumeric(BadTrials)

    bad_trials = BadTrials;

elseif iscell(BadTrials)

    if channel_index <= numel(BadTrials)
        bad_trials = BadTrials{channel_index};
    end
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

function output = compact_number_list(values)
% Produce a compact readable list for the command window.

if isempty(values)
    output = '(none)';
else
    output = strtrim(sprintf('%g ',values));
end
end

function output = logical_text(value)
% Convert a logical setting to ON or OFF.

if value
    output = 'ON';
else
    output = 'OFF';
end
end