clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder
data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

%% Choice
% ================= XIA MODIFICATION =================
% Spike filtering is removed from this quick plotting script.
% The raw spikes from .sp.mat will still be saved as sp_clipped in .sp_xia.mat
% so later analysis scripts can load .sp_xia.mat normally.

raster_chn_start = 1;
raster_chn_end = 32; %nChn
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% New plotting choices for mixed-prefix stimulation files.
%
% ActiveCounts_to_plot:
%   Which prefix lengths / active electrode counts to plot.
%   Example:
%       [1 2 3 4 5] plots all prefixes.
%       [5] only plots full 5-electrode prefix trials.
%       [] plots all detected active counts.
%
% ISI_to_plot_ms:
%   Which ISI values to plot, in ms.
%   Example:
%       [5] only plots 5 ms ISI.
%       [3 5 10] plots selected ISIs.
%       [] plots all detected ISIs.
%
% ConditionTypes_to_plot:
%   1 = sequential prefix / recovery condition
%   2 = full simultaneous condition
%   [] = plot all detected condition types except zero-control trials.
%
% SetIDs_to_plot:
%   Which entered stimulation orders to plot.
%   Example:
%       [1 2 3] plots three entered orders.
%       [] plots all detected sets.
ActiveCounts_to_plot = [1 2 3 4 5];  % [] for all prefixes
ISI_to_plot_ms = [5];                % [] for all ISIs
ConditionTypes_to_plot = [1];        % 1 = prefix/recovery, 2 = simultaneous
SetIDs_to_plot = [];                 % [] for all entered orders
% =====================================================

%% Spike Amplitude Filtering Parameters
% These parameters are kept here for compatibility / future use.
% They are not used in the current quick plotting script.
pos_limit = 100;    % upper bound (µV)
neg_limit = -100;  % lower bound (µV)

baseline_window_ms = [-60, -5];        % Baseline window (ms)
response_window_ms = [2, 25];          % Response window (ms)

%% Waveform templete filtering parameters
% These parameters are kept here for compatibility / future use.
% They are not used in the current quick plotting script.
template_window_ms = [5 15];   % use 0–5 ms spikes after trigger to build template
baseline_window_ms = [-60 0];    % window to filter
corr_thresh        = 0.70;       % threshold (increase to be stricter)

%% Load folders
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

% Extract file name
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);  % 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback if no underscores
end

%% Pre Set
FS=30000; % Sampling frequency

% Load .sp.mat file
% sp_files = dir('*.sp.mat');
sp_files = dir(fullfile(data_folder, '*.sp.mat'));
assert(~isempty(sp_files), 'No .sp.mat file found in the current folder.');

% sp_filename = sp_files(1).name;
sp_filename = fullfile(data_folder, sp_files(1).name);
fprintf('Loading spike file: %s\n', sp_filename);

S = load(sp_filename);
if isfield(S, 'sp')
    sp = S.sp;
else
    error('Variable "sp" not found in %s.', sp_filename);
end

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0); 

%% Load StimParams and decode amplitudes, stimulation sets, ISI
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

% ================= XIA MODIFICATION =================
% Load all variables from the experiment datafile.
% This allows us to use the new mixed-prefix metadata if it exists.
S = load(fullfile(data_folder, fileDIR(1).name));

StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
CHN                = S.CHN;
E_MAP              = S.E_MAP;
n_Trials           = S.n_Trials;
% =====================================================

%% Trial amplitude list
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);

% Convert disabled / zero-control amplitude from -1 to 0 for plotting label.
trialAmps(trialAmps == -1) = 0;

[Amps, ~, ampIdx] = unique(trialAmps(:));
n_AMP = numel(Amps);
cmap = lines(n_AMP);  % color map for amplitudes

% ================= XIA MODIFICATION =================
% Read mixed-prefix metadata if available.
%
% New mixed-prefix files contain:
%   active_electrode_count_by_trial
%   prefix_length_by_trial
%   isi_ms_by_trial
%   condition_type_by_trial
%   condition_set_id_by_trial
%
% If these variables do not exist, the code falls back to old-style decoding.
if isfield(S, 'active_electrode_count_by_trial')
    activeCount_trial = S.active_electrode_count_by_trial(:);
    prefixLength_trial = S.prefix_length_by_trial(:);
    isi_ms_trial = S.isi_ms_by_trial(:);
    conditionType_trial = S.condition_type_by_trial(:);
    conditionSetID_trial = S.condition_set_id_by_trial(:);
else
    fprintf('\nNo mixed-prefix metadata found. Using old-style decoding.\n');

    activeCount_trial = simultaneous_stim * ones(n_Trials,1);
    prefixLength_trial = simultaneous_stim * ones(n_Trials,1);
    conditionType_trial = ones(n_Trials,1);       % treat old files as normal stimulation condition
    conditionSetID_trial = ones(n_Trials,1);      % temporary label

    % For old files, estimate ISI/PTD from the last stored stimulation row.
    if simultaneous_stim > 4
        old_ptd_all = cell2mat(StimParams(6:simultaneous_stim:end, 6));
    elseif simultaneous_stim == 4
        old_ptd_all = cell2mat(StimParams(5:simultaneous_stim:end, 6));
    elseif simultaneous_stim == 3
        old_ptd_all = cell2mat(StimParams(4:simultaneous_stim:end, 6));
    elseif simultaneous_stim == 2
        old_ptd_all = cell2mat(StimParams(3:simultaneous_stim:end, 6));
    else
        old_ptd_all = zeros(n_Trials,1);
    end

    isi_ms_trial = old_ptd_all(:) ./ 1000;
end

% Calculate last active PTD for each trial.
% This is not used as the main grouping variable anymore, but it is useful
% for checking and for possible future spike blanking from 0 to final artifact.
lastActivePTD_us = zeros(n_Trials,1);
for tr = 1:n_Trials
    activeCount_this = activeCount_trial(tr);

    if isnan(activeCount_this) || activeCount_this < 1
        lastActivePTD_us(tr) = 0;
    else
        activeCount_this = min(round(activeCount_this), simultaneous_stim);
        stimRow = 1 + (tr-1)*simultaneous_stim + activeCount_this;
        lastActivePTD_us(tr) = StimParams{stimRow,6};
    end
end

PTD_us = lastActivePTD_us(:);
[PTD_values, ~, ptdIdx] = unique(PTD_us);

fprintf('\nDetected active electrode counts:'); disp(unique(activeCount_trial)');
fprintf('Detected prefix lengths:'); disp(unique(prefixLength_trial)');
fprintf('Detected ISIs (ms):'); disp(unique(isi_ms_trial)');
fprintf('Detected condition types:'); disp(unique(conditionType_trial)');
fprintf('Detected condition set IDs:'); disp(unique(conditionSetID_trial)');
fprintf('Detected last active PTDs (us):'); disp(unique(PTD_us)');
% =====================================================

%% Stimulation set selection
% ================= XIA MODIFICATION =================
% In the new mixed-prefix file, ConditionSetID directly indicates the
% entered stimulation order.
SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));

if isempty(SetIDs_to_plot)
    SetIDs_selected = SetIDs_all;
else
    SetIDs_selected = intersect(SetIDs_all, SetIDs_to_plot);
end

if isempty(ActiveCounts_to_plot)
    ActiveCounts_selected = unique(activeCount_trial(activeCount_trial > 0));
else
    ActiveCounts_selected = intersect(unique(activeCount_trial(activeCount_trial > 0)), ActiveCounts_to_plot);
end

if isempty(ISI_to_plot_ms)
    ISIs_selected = unique(isi_ms_trial(conditionType_trial == 1));
else
    ISIs_selected = intersect(unique(isi_ms_trial(conditionType_trial == 1)), ISI_to_plot_ms);
end

if isempty(ConditionTypes_to_plot)
    ConditionTypes_selected = unique(conditionType_trial(conditionType_trial > 0));
else
    ConditionTypes_selected = intersect(unique(conditionType_trial(conditionType_trial > 0)), ConditionTypes_to_plot);
end

fprintf('\nSet IDs selected for plotting:'); disp(SetIDs_selected');
fprintf('Prefixes selected for plotting:'); disp(ActiveCounts_selected');
fprintf('ISIs selected for plotting (ms):'); disp(ISIs_selected');
fprintf('Condition types selected for plotting:'); disp(ConditionTypes_selected');
% =====================================================

%% Pulse Train Period
pulseTrain_all = cell2mat(StimParams(2:end,9));  % Column 9: Pulse Train Period
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);  % take 1 per trial
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

%% Electrode Map
d = Depth_s(Electrode_Type); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

%% Save .sp.mat as .sp_xia.mat for later analysis
% ================= XIA MODIFICATION =================
% No spike filtering is applied in this quick plotting script.
%
% The original .sp.mat contains variable:
%   sp
%
% Later analysis scripts expect .sp_xia.mat to contain:
%   sp_clipped
%
% So here we simply rename sp -> sp_clipped and save it.
% The spike contents are unchanged.
sp_clipped = sp;
save([base_name '.sp_xia.mat'], 'sp_clipped', '-v7.3');
fprintf('\nSaved unfiltered spikes as %s.sp_xia.mat\n', base_name);
% =====================================================

%% Raster Plot Parameters
ras_win         = [-20 100];   % ms
bin_ms_raster   = 1;           % bin size
smooth_ms       = 2;           % smoothing window
% raster_chn_start = 1;
% raster_chn_end = 32; %nChn

%% === Initialize structure to store first-spike times ===
firstSpikeTimes = cell(raster_chn_end, 1); % each cell: vector of first-spike times per trial (ms)
fprintf('\nComputing First Spike Times per Trial\n');
post_spike_window_ms = [5,8];

%% Raster Plot 
 %% ========================= RASTER + PSTH ========================= %%
edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2); 
g = g / sum(g);

% ================= XIA MODIFICATION =================
% New plotting logic:
%   Channel
%     Set/order
%       Condition type
%         Prefix
%           ISI
%             Amplitude
%
% ConditionType:
%   1 = sequential prefix / recovery
%   2 = full simultaneous
%
% For ConditionType == 2, ISI is ignored because simultaneous stimulation
% has ISI_ms = 0 by definition.
% =====================================================

for ich = raster_chn_start:raster_chn_end
    ch = d(ich);

    if ch > numel(sp_clipped)
        continue;
    end

    if isempty(sp_clipped{ch}), continue; end

    for si = 1:numel(SetIDs_selected)
        setID = SetIDs_selected(si);

        % Build a readable label from CHN.
        % In the mixed-prefix design, condition set ID corresponds to the
        % entered electrode order.
        if setID <= size(CHN,1)
            stimVec = CHN(setID,:);
            stimVec = stimVec(stimVec > 0);
            setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimVec, 'UniformOutput', false), '→');
        else
            setLabel = sprintf('Set%d', setID);
        end

        for ci = 1:numel(ConditionTypes_selected)
            condType = ConditionTypes_selected(ci);

            % =====================================================
            % Sequential prefix/recovery condition
            % =====================================================
            if condType == 1

                for pi = 1:n_PULSE
                    pulse_val = PulsePeriods(pi);

                    for aci = 1:numel(ActiveCounts_selected)
                        activeCount = ActiveCounts_selected(aci);

                        for ii = 1:numel(ISIs_selected)
                            isi_val = ISIs_selected(ii);

                            % Find trials for this set, prefix, ISI, and pulse period.
                            trials_this_group = find(conditionSetID_trial == setID & ...
                                                     conditionType_trial == condType & ...
                                                     activeCount_trial == activeCount & ...
                                                     isi_ms_trial == isi_val & ...
                                                     pulseIdx == pi);

                            if isempty(trials_this_group), continue; end

                            % --------- FIGURE ---------
                            figName = sprintf('Ch %d | Set %d %s | Prefix %d | ISI %.1f ms | Pulse %d µs', ...
                                              ich, setID, setLabel, activeCount, isi_val, pulse_val);

                            figure('Color','w','Name',figName);
                            tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                            ax1 = nexttile([3 1]);
                            hold(ax1,'on'); box(ax1,'off');
                            title(ax1, sprintf('Raster — Ch %d | Set %d %s | Prefix %d | ISI %.1f ms', ...
                                               ich, setID, setLabel, activeCount, isi_val), ...
                                               'Interpreter','none');

                            ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');

                            % ===== PSTH storage =====
                            psth_curves = cell(1, n_AMP);
                            maxRate = 0;
                            y_cursor = 0;
                            ytick_vals = [];
                            ytick_labels = {};

                            % ================= LOOP AMPLITUDES ==================
                            for ai = 1:n_AMP
                                amp_val = Amps(ai);
                                color = cmap(ai,:);

                                % Do not include zero-control trials in these
                                % condition-specific plots. Zero trials have
                                % ConditionType = 0 and ActiveCount = 0.
                                amp_trials = find(ampIdx == ai & ...
                                                  pulseIdx == pi & ...
                                                  conditionSetID_trial == setID & ...
                                                  conditionType_trial == condType & ...
                                                  activeCount_trial == activeCount & ...
                                                  isi_ms_trial == isi_val);

                                nTr = numel(amp_trials);
                                if nTr == 0
                                    psth_curves{ai} = zeros(1, numel(ctrs));
                                    continue;
                                end

                                % ======== RASTER + PSTH COUNTING ========
                                counts = zeros(1, numel(edges)-1);
                                for t = 1:nTr
                                    tr = amp_trials(t);

                                    if tr > length(trig)
                                        continue;
                                    end

                                    t0 = trig(tr)/FS*1000;

                                    tt = sp_clipped{ch}(:,1);
                                    tt = tt(tt >= t0+ras_win(1) & tt <= t0+ras_win(2)) - t0;

                                    % ---- RASTER ----
                                    y0 = y_cursor + t;
                                    for spike_t = tt'
                                        plot(ax1, [spike_t spike_t], [y0-0.4 y0+0.4], ...
                                             'Color', color, 'LineWidth', 1.1);
                                    end

                                    % ---- PSTH ----
                                    counts = counts + histcounts(tt, edges);
                                end

                                % y-axis structure
                                ytick_vals(end+1) = y_cursor + nTr/2;
                                ytick_labels{end+1} = sprintf('%d µA', amp_val);
                                y_cursor = y_cursor + nTr;

                                % ---- Compute PSTH ----
                                rate = filter(g, 1, counts/(nTr*bin_s));
                                psth_curves{ai} = rate;
                                maxRate = max(maxRate, max(rate));
                            end

                            % Finalize raster axis
                            xline(ax1, 0, 'r--');
                            xlim(ax1, ras_win);
                            ylim(ax1, [0 max(1,y_cursor)]);
                            yticks(ax1, ytick_vals);
                            yticklabels(ax1, ytick_labels);
                            ylabel(ax1, 'Amplitude');

                            % Finalize PSTH
                            for ai = 1:n_AMP
                                plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.6);
                            end
                            xline(ax2, 0, 'r--');
                            xlim(ax2, ras_win);
                            ylim(ax2, [0 max(50,ceil(maxRate*1.1/10)*10)]);
                            xlabel(ax2, 'Time (ms)');
                            ylabel(ax2, 'Rate (sp/s)');
                            legend(ax2, arrayfun(@(a) sprintf('%.0f µA',a), Amps, 'UniformOutput', false), ...
                                   'Box','off','Location','northeast');

                        end % ISI
                    end % Active count / Prefix
                end % Pulse

            % =====================================================
            % Full simultaneous condition
            % =====================================================
            elseif condType == 2

                for pi = 1:n_PULSE
                    pulse_val = PulsePeriods(pi);

                    % In the mixed-prefix file, full simultaneous normally
                    % has ActiveElectrodeCount = simultaneous_stim and ISI_ms = 0.
                    activeCounts_sim = intersect(ActiveCounts_selected, unique(activeCount_trial(conditionType_trial == 2)));

                    if isempty(activeCounts_sim)
                        activeCounts_sim = unique(activeCount_trial(conditionType_trial == 2));
                    end

                    for aci = 1:numel(activeCounts_sim)
                        activeCount = activeCounts_sim(aci);

                        trials_this_group = find(conditionSetID_trial == setID & ...
                                                 conditionType_trial == condType & ...
                                                 activeCount_trial == activeCount & ...
                                                 pulseIdx == pi);

                        if isempty(trials_this_group), continue; end

                        % --------- FIGURE ---------
                        figName = sprintf('Ch %d | Set %d %s | Simultaneous | Pulse %d µs', ...
                                          ich, setID, setLabel, pulse_val);

                        figure('Color','w','Name',figName);
                        tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                        ax1 = nexttile([3 1]);
                        hold(ax1,'on'); box(ax1,'off');
                        title(ax1, sprintf('Raster — Ch %d | Set %d %s | Simultaneous', ...
                                           ich, setID, setLabel), ...
                                           'Interpreter','none');

                        ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');

                        % ===== PSTH storage =====
                        psth_curves = cell(1, n_AMP);
                        maxRate = 0;
                        y_cursor = 0;
                        ytick_vals = [];
                        ytick_labels = {};

                        % ================= LOOP AMPLITUDES ==================
                        for ai = 1:n_AMP
                            amp_val = Amps(ai);
                            color = cmap(ai,:);

                            amp_trials = find(ampIdx == ai & ...
                                              pulseIdx == pi & ...
                                              conditionSetID_trial == setID & ...
                                              conditionType_trial == condType & ...
                                              activeCount_trial == activeCount);

                            nTr = numel(amp_trials);
                            if nTr == 0
                                psth_curves{ai} = zeros(1, numel(ctrs));
                                continue;
                            end

                            % ======== RASTER + PSTH COUNTING ========
                            counts = zeros(1, numel(edges)-1);
                            for t = 1:nTr
                                tr = amp_trials(t);

                                if tr > length(trig)
                                    continue;
                                end

                                t0 = trig(tr)/FS*1000;

                                tt = sp_clipped{ch}(:,1);
                                tt = tt(tt >= t0+ras_win(1) & tt <= t0+ras_win(2)) - t0;

                                % ---- RASTER ----
                                y0 = y_cursor + t;
                                for spike_t = tt'
                                    plot(ax1, [spike_t spike_t], [y0-0.4 y0+0.4], ...
                                         'Color', color, 'LineWidth', 1.1);
                                end

                                % ---- PSTH ----
                                counts = counts + histcounts(tt, edges);
                            end

                            % y-axis structure
                            ytick_vals(end+1) = y_cursor + nTr/2;
                            ytick_labels{end+1} = sprintf('%d µA', amp_val);
                            y_cursor = y_cursor + nTr;

                            % ---- Compute PSTH ----
                            rate = filter(g, 1, counts/(nTr*bin_s));
                            psth_curves{ai} = rate;
                            maxRate = max(maxRate, max(rate));
                        end

                        % Finalize raster axis
                        xline(ax1, 0, 'r--');
                        xlim(ax1, ras_win);
                        ylim(ax1, [0 max(1,y_cursor)]);
                        yticks(ax1, ytick_vals);
                        yticklabels(ax1, ytick_labels);
                        ylabel(ax1, 'Amplitude');

                        % Finalize PSTH
                        for ai = 1:n_AMP
                            plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.6);
                        end
                        xline(ax2, 0, 'r--');
                        xlim(ax2, ras_win);
                        ylim(ax2, [0 max(50,ceil(maxRate*1.1/10)*10)]);
                        xlabel(ax2, 'Time (ms)');
                        ylabel(ax2, 'Rate (sp/s)');
                        legend(ax2, arrayfun(@(a) sprintf('%.0f µA',a), Amps, 'UniformOutput', false), ...
                               'Box','off','Location','northeast');

                    end % active count
                end % pulse
            end % condition type
        end % condition type
    end % Set
end % Channel