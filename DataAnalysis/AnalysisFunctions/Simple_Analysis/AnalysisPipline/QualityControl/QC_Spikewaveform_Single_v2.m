%% =============================================================
% SINGLE-STIMULATION SPIKE-WAVEFORM REVIEW
%
% PURPOSE
%   Plot post-stimulation spike waveforms for manual identification of
%   channels that remain biased by stimulation artifacts.
%
% DATA HANDLING
%   - Prioritizes sp_corr from the SSD file.
%   - Does not use responding-channel results.
%   - Does not load or remove bad trials.
%   - Does not modify or save experimental data.
%
% FIGURE ORGANIZATION
%   One figure is produced for each:
%       recording channel × stimulation set
%
%   Each tile represents one post-stimulation time bin.
%   Different stimulation amplitudes are displayed with different colors.
% =============================================================

clear;
close all;

addpath(genpath( ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ========================= USER SETTINGS ==============================

% Single-stimulation dataset
data_folder = ...
    '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_Single1';

% Probe type:
%   0 = rigid single shank
%   1 = flexible single shank
%   2 = flexible four shank
Electrode_Type = 2;

% Default recording-channel range, using Depth_s channel indices
spike_chn_start = 1;
spike_chn_end   = 64;

% Optional explicit recording-channel selection.
% Empty means use spike_chn_start:spike_chn_end.
% Examples:
%   Channels_To_Plot = [];
%   Channels_To_Plot = 31:64;
%   Channels_To_Plot = [35 38 42 47];
Channels_To_Plot = [];

% Optional stimulation-set selection.
% Empty means plot every stimulation set.
% Examples:
%   Sets_To_Plot = [];
%   Sets_To_Plot = 1;
%   Sets_To_Plot = [1 2];
Sets_To_Plot = [];

% Sampling frequency
FS = 30000;

% Post-stimulation interval included in the waveform review
waveform_window_ms = [0 30];

% Width of each post-stimulation time bin
bin_ms = 2;

% Set Inf to retain all waveform amplitudes.
% A finite value removes waveforms exceeding that absolute voltage, which
% can hide artifacts and is therefore not recommended for initial QC.
amp_threshold = Inf;

% false = plot waveforms without realignment; recommended for artifact QC
% true  = align waveform minima to the center using non-circular shifting
Align_Waveforms = false;

% Figure layout
layout_row = 3;
layout_col = 5;

figure_position = [100 100 1400 800];

%% ======================= VALIDATE SETTINGS ============================

if ~isfolder(data_folder)
    error('Folder does not exist:\n%s',data_folder);
end

if numel(waveform_window_ms) ~= 2 || ...
        waveform_window_ms(2) <= waveform_window_ms(1)
    error('waveform_window_ms must be [start end], with end > start.');
end

window_duration_ms = diff(waveform_window_ms);
nBins_exact = window_duration_ms/bin_ms;

if abs(nBins_exact-round(nBins_exact)) > 1e-9
    error(['The waveform-window duration must be exactly divisible ' ...
           'by bin_ms.']);
end

nBins = round(nBins_exact);

if layout_row*layout_col < nBins
    error(['The selected layout has %d tiles but %d waveform bins ' ...
           'are required.'],layout_row*layout_col,nBins);
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder));

cd(data_folder);

fprintf('\n');
fprintf('============================================================\n');
fprintf('SINGLE-STIMULATION SPIKE-WAVEFORM REVIEW\n');
fprintf('============================================================\n');
fprintf('Dataset:\n%s\n',data_folder);
fprintf('Waveform interval: [%g,%g) ms\n',waveform_window_ms);
fprintf('Bin width: %g ms\n',bin_ms);
fprintf('Number of bins: %d\n',nBins);
fprintf('Amplitude plotting threshold: %g uV\n',amp_threshold);
fprintf('Waveform alignment: %s\n', ...
    logical_to_text(Align_Waveforms));

%% ========================= FIND SPIKE FILE ============================

ssd_files  = dir('*.sp_xia_SSD.mat');
base_files = dir('*.sp_xia.mat');

if ~isempty(ssd_files)

    if numel(ssd_files) > 1
        warning('Multiple SSD files found. Using: %s', ...
            ssd_files(1).name);
    end

    spike_file = ssd_files(1).name;

elseif ~isempty(base_files)

    if numel(base_files) > 1
        warning('Multiple base spike files found. Using: %s', ...
            base_files(1).name);
    end

    spike_file = base_files(1).name;

else
    error('No *.sp_xia_SSD.mat or *.sp_xia.mat file was found.');
end

%% ========================== LOAD SPIKES ===============================

available_variables = who('-file',spike_file);

if ismember('sp_corr',available_variables)

    spike_variable = 'sp_corr';

elseif ismember('sp_SSD',available_variables)

    spike_variable = 'sp_SSD';
    warning('sp_corr is missing. Using sp_SSD.');

elseif ismember('sp_in',available_variables)

    spike_variable = 'sp_in';
    warning('sp_corr is missing. Using sp_in.');

elseif ismember('sp_clipped',available_variables)

    spike_variable = 'sp_clipped';
    warning('sp_corr is missing. Using sp_clipped.');

elseif ismember('sp',available_variables)

    spike_variable = 'sp';
    warning('sp_corr is missing. Using sp.');

else
    error('No usable spike variable was found in:\n%s',spike_file);
end

SpikeLoad = load(spike_file,spike_variable);
sp_use = SpikeLoad.(spike_variable);

nCh = numel(sp_use);

fprintf('Spike file: %s\n',spike_file);
fprintf('Spike variable: %s\n',spike_variable);
fprintf('Spike channels: %d\n',nCh);

%% ========================== LOAD TRIGGERS =============================

if isempty(dir('*.trig.dat'))

    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);
end

if isempty(dir('*.trig.dat'))
    error('No *.trig.dat file could be found or created.');
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Triggers loaded: %d\n',nTrig);

%% ================= LOAD EXPERIMENT PARAMETERS =========================

experiment_files = dir('*_exp_datafile_*.mat');

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found.');
end

if numel(experiment_files) > 1
    warning('Multiple experiment files found. Using: %s', ...
        experiment_files(1).name);
end

ExpLoad = load(experiment_files(1).name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');

required_fields = {
    'StimParams'
    'simultaneous_stim'
    'E_MAP'
    'n_Trials'
};

for k = 1:numel(required_fields)

    if ~isfield(ExpLoad,required_fields{k})
        error('Experiment file is missing variable: %s', ...
            required_fields{k});
    end
end

StimParams        = ExpLoad.StimParams;
simultaneous_stim = ExpLoad.simultaneous_stim;
E_MAP             = ExpLoad.E_MAP;
n_Trials          = ExpLoad.n_Trials;

fprintf('Experiment file: %s\n',experiment_files(1).name);
fprintf('Experiment trials: %d\n',n_Trials);
fprintf('Stimulation events per trial: %d\n',simultaneous_stim);

if simultaneous_stim ~= 1
    warning(['simultaneous_stim is %d rather than 1. This may not be ' ...
             'a single-stimulation dataset.'],simultaneous_stim);
end

if nTrig < n_Trials

    error(['Only %d triggers were loaded for %d trials. Check the ' ...
           'trigger file before waveform review.'],nTrig,n_Trials);

elseif nTrig > n_Trials

    warning(['There are %d triggers but %d experiment trials. Only ' ...
             'experiment trial IDs will be used.'],nTrig,n_Trials);
end

%% ========================= DECODE AMPLITUDES ==========================

trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);

if numel(trialAmps) ~= n_Trials
    error(['Decoded %d trial amplitudes for %d experiment trials. ' ...
           'Check StimParams and simultaneous_stim.'], ...
           numel(trialAmps),n_Trials);
end

% Convert experiment control value -1 to 0 uA before unique indexing
trialAmps(trialAmps == -1) = 0;

[Amps,~,ampIdx] = unique(trialAmps(:));

n_AMP = numel(Amps);
cmap = lines(n_AMP);

fprintf('Available amplitudes: %s uA\n',num2str(Amps(:).'));

%% ==================== DECODE STIMULATION SETS =========================

E_NAME = E_MAP(2:end);
stimNames = StimParams(2:end,1);

[isMapped,idx_all] = ismember(stimNames,E_NAME);

if any(~isMapped)
    warning('%d stimulation entries could not be mapped through E_MAP.', ...
        sum(~isMapped));
end

stimChPerTrial_all = cell(n_Trials,1);

for trial_id = 1:n_Trials

    rows_this_trial = ...
        (trial_id-1)*simultaneous_stim + ...
        (1:simultaneous_stim);

    stimulation_indices = ...
        unique(idx_all(rows_this_trial),'stable');

    stimulation_indices = ...
        stimulation_indices(stimulation_indices > 0).';

    stimChPerTrial_all{trial_id} = ...
        stimulation_indices;
end

comb = zeros(n_Trials,simultaneous_stim);

for trial_id = 1:n_Trials

    stimulation_indices = ...
        stimChPerTrial_all{trial_id};

    comb(trial_id,1:numel(stimulation_indices)) = ...
        stimulation_indices;
end

[uniqueComb,~,combClass] = ...
    unique(comb,'rows','stable');

nSets = size(uniqueComb,1);

fprintf('Detected stimulation sets: %d\n',nSets);

for set_id = 1:nSets

    stimIdx = uniqueComb(set_id,:);
    stimIdx = stimIdx(stimIdx > 0);

    stimLabel = strjoin( ...
        arrayfun(@(channel_index) ...
        sprintf('Ch%d',channel_index), ...
        stimIdx, ...
        'UniformOutput',false), ...
        ', ');

    fprintf('  Set %d: %s\n',set_id,stimLabel);
end

%% ========================= ELECTRODE MAP ==============================

% Depth_s maps:
%   depth/channel index -> spike-data channel index
d = Depth_s(Electrode_Type);

if isempty(d)
    error('Depth_s returned an empty channel map.');
end

%% ================= SELECT RECORDING CHANNELS ==========================

if isempty(Channels_To_Plot)

    selected_channels = ...
        spike_chn_start:spike_chn_end;

else

    selected_channels = ...
        unique(double(Channels_To_Plot(:).'),'stable');
end

valid_channel_mask = ...
    selected_channels >= 1 & ...
    selected_channels <= numel(d);

if any(~valid_channel_mask)

    warning('Ignoring recording-channel indices outside 1:%d.', ...
        numel(d));

    selected_channels = ...
        selected_channels(valid_channel_mask);
end

if isempty(selected_channels)
    error('No valid recording channels were selected.');
end

%% ================= SELECT STIMULATION SETS ============================

if isempty(Sets_To_Plot)

    selected_sets = 1:nSets;

else

    selected_sets = ...
        unique(double(Sets_To_Plot(:).'),'stable');

    valid_set_mask = ...
        selected_sets >= 1 & ...
        selected_sets <= nSets & ...
        fix(selected_sets) == selected_sets;

    if any(~valid_set_mask)
        warning('Ignoring unavailable stimulation-set indices.');
    end

    selected_sets = selected_sets(valid_set_mask);
end

if isempty(selected_sets)
    error('No valid stimulation sets were selected.');
end

maximum_figures = ...
    numel(selected_channels)*numel(selected_sets);

fprintf('Recording channels selected: %s\n', ...
    num2str(selected_channels));

fprintf('Stimulation sets selected: %s\n', ...
    num2str(selected_sets));

fprintf('Maximum number of figures: %d\n',maximum_figures);

%% ================= DETERMINE WAVEFORM LENGTH ==========================

example_ch = find(~cellfun(@isempty,sp_use),1,'first');

if isempty(example_ch)
    error('All channels are empty in the selected spike variable.');
end

if size(sp_use{example_ch},2) < 2
    error('Spike data do not contain waveform columns.');
end

wf_len = size(sp_use{example_ch},2)-1;

% Time within the stored waveform snippet
t_wave = (0:wf_len-1)/FS*1000;

fprintf('Stored waveform samples: %d\n',wf_len);
fprintf('Stored waveform duration: %.3f ms\n',t_wave(end));

%% ================= SPIKE-WAVEFORM PLOTTING ============================

figures_created = 0;

for ich = selected_channels

    % Convert the selected depth/channel index to the spike-data cell index
    ch = d(ich);

    if ch < 1 || ch > nCh

        warning('Ch %d maps outside the spike-data channel range.',ich);
        continue;
    end

    if isempty(sp_use{ch})

        fprintf('Ch %d has no spikes. Skipping.\n',ich);
        continue;
    end

    if size(sp_use{ch},2)-1 ~= wf_len

        warning(['Ch %d has a different waveform length. ' ...
                 'Skipping this channel.'],ich);
        continue;
    end

    sp_times = double(sp_use{ch}(:,1));
    sp_wave  = double(sp_use{ch}(:,2:end));

    %% ---------------- OPTIONAL AMPLITUDE LIMIT ------------------------

    if isfinite(amp_threshold)

        valid_waveform = ...
            all(abs(sp_wave) <= amp_threshold,2);

        removed_waveforms = ...
            sum(~valid_waveform);

        if removed_waveforms > 0
            fprintf(['Ch %d: %d waveforms exceed the plotting threshold ' ...
                     'and will not be shown.\n'], ...
                     ich,removed_waveforms);
        end

        sp_times = sp_times(valid_waveform);
        sp_wave  = sp_wave(valid_waveform,:);
    end

    if isempty(sp_times)
        fprintf('Ch %d has no waveforms after filtering. Skipping.\n',ich);
        continue;
    end

    %% ---------------- STIMULATION-SET LOOP ----------------------------

    for set_id = selected_sets

        trial_ids = find(combClass == set_id);

        if isempty(trial_ids)
            continue;
        end

        all_spikes_by_bin_amp = cell(nBins,n_AMP);

        %% ------------ COLLECT WAVEFORMS BY BIN AND AMP ---------------

        for trial_position = 1:numel(trial_ids)

            trial_id = trial_ids(trial_position);
            trigger_ms = trig(trial_id)/FS*1000;
            amp_id = ampIdx(trial_id);

            absolute_window = ...
                trigger_ms+waveform_window_ms;

            spike_mask = ...
                sp_times >= absolute_window(1) & ...
                sp_times <  absolute_window(2);

            if ~any(spike_mask)
                continue;
            end

            relative_times = ...
                sp_times(spike_mask)-trigger_ms;

            waveforms = ...
                sp_wave(spike_mask,:);

            for spike_id = 1:numel(relative_times)

                bin_index = floor( ...
                    (relative_times(spike_id)-waveform_window_ms(1)) ...
                    /bin_ms)+1;

                if bin_index < 1 || bin_index > nBins
                    continue;
                end

                all_spikes_by_bin_amp{bin_index,amp_id}(end+1,:) = ...
                    waveforms(spike_id,:); %#ok<AGROW>
            end
        end

        %% ---------------- CHECK AVAILABLE WAVEFORMS ------------------

        nonempty_cells = ...
            ~cellfun(@isempty,all_spikes_by_bin_amp);

        if ~any(nonempty_cells,'all')
            continue;
        end

        all_waves = ...
            cell2mat(all_spikes_by_bin_amp(nonempty_cells));

        if isempty(all_waves)
            continue;
        end

        y_max = max(abs(all_waves(:)),[],'omitnan');

        if ~isfinite(y_max) || y_max <= 0
            y_max = 50;
        end

        y_limit_value = max(50,ceil(y_max/50)*50);
        y_lim = [-y_limit_value y_limit_value];

        %% ---------------- STIMULATION LABEL ---------------------------

        % Keep the original straightforward channel index labels
        stimIdx = uniqueComb(set_id,:);
        stimIdx = stimIdx(stimIdx > 0);

        stimLabel = strjoin( ...
            arrayfun(@(channel_index) ...
            sprintf('Ch%d',channel_index), ...
            stimIdx, ...
            'UniformOutput',false), ...
            ', ');

        %% ---------------- CREATE FIGURE -------------------------------

        figure_name = sprintf( ...
            'Recording Ch %d | StimSet %d (%s) | Single Pulse', ...
            ich,set_id,stimLabel);

        figure( ...
            'Name',figure_name, ...
            'NumberTitle','off', ...
            'Color','w', ...
            'Position',figure_position);

        tiledlayout( ...
            layout_row,layout_col, ...
            'Padding','compact', ...
            'TileSpacing','compact');

        sgtitle(figure_name, ...
            'FontWeight','bold', ...
            'Interpreter','none');

        figures_created = figures_created+1;

        %% ---------------- PLOT EACH TIME BIN --------------------------

        for bin_index = 1:nBins

            ax = nexttile;
            hold(ax,'on');

            for amp_id = 1:n_AMP

                waves = ...
                    all_spikes_by_bin_amp{bin_index,amp_id};

                if isempty(waves)
                    continue;
                end

                if Align_Waveforms

                    waves_to_plot = ...
                        align_waveforms_without_wrap(waves);

                else

                    waves_to_plot = waves;
                end

                % Blend the amplitude color with white instead of using
                % an unsupported four-element RGBA line color.
                waveform_color = ...
                    0.55*cmap(amp_id,:)+0.45*[1 1 1];

                plot(ax, ...
                    t_wave, ...
                    waves_to_plot', ...
                    'Color',waveform_color, ...
                    'LineWidth',0.5);
            end

            %% ---------------- TILE ANNOTATION -------------------------

            spike_count = 0;

            for amp_id = 1:n_AMP
                spike_count = spike_count+size( ...
                    all_spikes_by_bin_amp{bin_index,amp_id},1);
            end

            bin_start_ms = waveform_window_ms(1)+ ...
                (bin_index-1)*bin_ms;

            bin_end_ms = bin_start_ms+bin_ms;

            title(ax, ...
                sprintf('%g-%g ms (%d spikes)', ...
                bin_start_ms,bin_end_ms,spike_count), ...
                'Interpreter','none');

            xlabel(ax,'Waveform time (ms)');
            ylabel(ax,'Voltage (uV)');

            ylim(ax,y_lim);

            yticks(ax,linspace(y_lim(1),y_lim(2),3));
            xticks(ax,round(linspace(t_wave(1),t_wave(end),3),3));

            axis(ax,'square');
            grid(ax,'on');
            box(ax,'off');
        end

        %% ---------------- AMPLITUDE LEGEND ----------------------------

        legend_handles = gobjects(n_AMP,1);
        legend_labels  = cell(n_AMP,1);

        for amp_id = 1:n_AMP

            legend_handles(amp_id) = plot( ...
                nan,nan,'-', ...
                'Color',cmap(amp_id,:), ...
                'LineWidth',1.5);

            legend_labels{amp_id} = ...
                sprintf('%g uA',Amps(amp_id));
        end

        legend( ...
            legend_handles, ...
            legend_labels, ...
            'Location','northeastoutside');
    end
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('WAVEFORM REVIEW COMPLETE\n');
fprintf('Figures created: %d\n',figures_created);
fprintf('No bad trials were loaded or removed.\n');
fprintf('No responding-channel labels were used.\n');
fprintf('No experimental files were modified.\n');
fprintf('============================================================\n');

%% ========================== LOCAL FUNCTIONS ===========================

function aligned = align_waveforms_without_wrap(waves)
% Align each waveform minimum to the center without circular wrapping.
% Vacated samples are filled with NaN and are not connected in the plot.

nWaveforms = size(waves,1);
nSamples   = size(waves,2);

target_index = ceil(nSamples/2);

aligned = nan(size(waves));

for waveform_id = 1:nWaveforms

    waveform = waves(waveform_id,:);

    [~,minimum_index] = min(waveform);

    sample_shift = target_index-minimum_index;

    if sample_shift == 0

        aligned(waveform_id,:) = waveform;

    elseif sample_shift > 0

        destination = ...
            (1+sample_shift):nSamples;

        source = ...
            1:(nSamples-sample_shift);

        aligned(waveform_id,destination) = ...
            waveform(source);

    else

        left_shift = abs(sample_shift);

        destination = ...
            1:(nSamples-left_shift);

        source = ...
            (1+left_shift):nSamples;

        aligned(waveform_id,destination) = ...
            waveform(source);
    end
end

end

function output = logical_to_text(value)

if value
    output = 'enabled';
else
    output = 'disabled';
end

end