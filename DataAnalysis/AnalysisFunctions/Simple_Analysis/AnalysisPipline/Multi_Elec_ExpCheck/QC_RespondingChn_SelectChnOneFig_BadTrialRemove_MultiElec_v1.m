%% =============================================================
%  Cleaned Raster + PSTH Selected-Channel Viewer
%  for Mixed-Prefix Multi-Electrode Stimulation
%
%  Purpose:
%    Plot selected recording channels in one figure for each condition.
%    The selected channels can be continuous or discontinuous.
%    Bad trials are removed before plotting.
%
%  Bad trial file:
%    *_MixedPrefixBadTrials.mat
%
%  Responding channel file:
%    *_RespondingChannels_FullSeqAndSim.mat
%    or
%    *_RespondingChannels_AllPrefixesAndSim.mat
%
%  ConditionType:
%    0 = zero-control
%    1 = sequential prefix/recovery
%    2 = full simultaneous
% =============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ========================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

Electrode_Type = 2;       % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Recording channels to plot, using Depth_s index.
% This can be continuous or discontinuous.

raster_channels_to_plot = [8:10];

FS = 30000;

% Plotting windows.
ras_win       = [-50 80];   % ms, time relative to trigger
bin_ms_raster = 1;          % ms, PSTH bin size
smooth_ms     = 5;          % ms, smoothing width

% Which amplitudes to plot.
% [] means all non-zero amplitudes.
Plot_Amps = [10];

% Which set/order IDs to plot.
% [] means all detected sets.
Plot_SetIDs = [];

% Which condition types to plot:
%   1 = sequential prefix/recovery
%   2 = simultaneous
% [] means both detected non-zero condition types.
Plot_ConditionTypes = [1 2];

% Which prefixes to plot for sequential trials.
% [] means all detected prefixes.
% For full-sequence + simultaneous overview, use [5].
Plot_Prefixes = [5];

% Which ISIs to plot for sequential trials.
% [] means all detected ISIs.
Plot_ISIs_ms = [5];

% Full simultaneous active count.
% For 5-electrode experiment, usually 5.
full_prefix_length = 5;

% If true, simultaneous trials are pooled across all set/order IDs.
% If false, simultaneous trials must match the current setID.
pool_simultaneous_across_sets = false;

%% ====================== BAD TRIAL SETTINGS ====================

% If true, load *_MixedPrefixBadTrials.mat and remove bad trials.
remove_bad_trials = true;

% If true, skip plotting if all trials in a condition are removed.
skip_if_no_clean_trials = true;

%% ====================== RESPONDING CHANNEL SETTINGS ===========

% Which responding-channel file to load:
%   'FullSeqAndSim'
%   'AllPrefixesAndSim'
RespondingFileMode = 'FullSeqAndSim';

% Highlight responding channels.
highlight_responding_channels = true;

% Background colour for responding channels.
resp_bg_color = [1.0 0.88 0.88];

%% ====================== FIGURE SETTINGS =======================

fig_position = [50 50 1600 900];

%% ===================== INITIAL SETUP =========================

if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end

cd(data_folder);
fprintf('Cleaned Raster+PSTH plotting in folder:\n%s\n\n', data_folder);

%% ===================== BASE NAME =============================

parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');

if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ===================== LOAD SPIKES ===========================

sp = [];

ssd_file    = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
prefix_file = dir(fullfile(data_folder, '*.sp_xia_PrefixRecovery.mat'));
base_file   = dir(fullfile(data_folder, '*.sp_xia.mat'));

if ~isempty(ssd_file)

    spike_file_used = ssd_file(1).name;
    S_sp = load(fullfile(data_folder, spike_file_used));

    if isfield(S_sp, 'sp_corr')
        sp = S_sp.sp_corr;
        spike_variable_used = 'sp_corr';
    elseif isfield(S_sp, 'sp_SSD')
        sp = S_sp.sp_SSD;
        spike_variable_used = 'sp_SSD';
    elseif isfield(S_sp, 'sp_in')
        sp = S_sp.sp_in;
        spike_variable_used = 'sp_in';
    elseif isfield(S_sp, 'sp_pca')
        sp = S_sp.sp_pca;
        spike_variable_used = 'sp_pca';
    else
        error('SSD file found but no usable spike variable was found.');
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
    error('No spike file found.');
end

nCh = numel(sp);

fprintf('Loaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);
fprintf('Loaded %d spike channels.\n', nCh);

%% ===================== LOAD RESPONDING CHANNELS ==============

Resp = [];
hasResp = false;

resp_file = sprintf('%s_RespondingChannels_%s.mat', base_name, RespondingFileMode);

if isfile(resp_file)

    tmp = load(resp_file);

    if isfield(tmp, 'Responding')
        Resp = tmp.Responding;
        hasResp = true;
        fprintf('Loaded responding-channel data from %s\n', resp_file);
    else
        warning('Responding file found but Responding struct missing.');
    end

else
    fprintf('No responding-channel file found: %s\n', resp_file);
    fprintf('Channels will not be highlighted as RESP.\n');
end

%% ===================== LOAD BAD TRIALS ========================

BadTrials = {};
BadTrials_global = [];
hasBadTrials = false;

bad_trial_file = sprintf('%s.MixedPrefixBadTrials.mat', base_name);
bad_trial_path = fullfile(data_folder, bad_trial_file);

if remove_bad_trials

    if isfile(bad_trial_path)

        tmp = load(bad_trial_path);

        if isfield(tmp, 'BadTrials')

            BadTrials = tmp.BadTrials;
            BadTrials_global = collect_global_bad_trials(BadTrials);
            hasBadTrials = true;

            fprintf('Loaded bad-trial file: %s\n', bad_trial_file);
            fprintf('Number of unique bad absolute trials: %d\n', numel(BadTrials_global));

        else
            warning('Bad-trial file found, but BadTrials variable is missing. No trials removed.');
        end

    else
        fprintf('No bad-trial file found: %s\n', bad_trial_file);
        fprintf('No bad trials will be removed.\n');
    end
else
    fprintf('remove_bad_trials = false. Plotting raw trials.\n');
end

%% ===================== LOAD TRIGGERS =========================

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ===================== LOAD STIM PARAMS ======================

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;
n_Trials          = S_exp.n_Trials;
E_MAP             = S_exp.E_MAP; %#ok<NASGU>

if isfield(S_exp, 'CHN')
    CHN = S_exp.CHN;
else
    CHN = [];
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match nTrig (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ===================== REMOVE HEADER ROW =====================

StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;

if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
        size(StimParams_data,1), expected_rows);
end

%% ===================== TRIAL METADATA ========================

if size(StimParams_data,2) < 30
    error('StimParams does not contain columns 26–30.');
end

firstRow_eachTrial = 1:simultaneous_stim:size(StimParams_data,1);

activeCount_trial    = cell2mat(StimParams_data(firstRow_eachTrial,26));
prefixLength_trial   = cell2mat(StimParams_data(firstRow_eachTrial,27));
isi_ms_trial         = cell2mat(StimParams_data(firstRow_eachTrial,28));
conditionType_trial  = cell2mat(StimParams_data(firstRow_eachTrial,29));
conditionSetID_trial = cell2mat(StimParams_data(firstRow_eachTrial,30));

activeCount_trial    = activeCount_trial(:);
prefixLength_trial   = prefixLength_trial(:);
isi_ms_trial         = isi_ms_trial(:);
conditionType_trial  = conditionType_trial(:);
conditionSetID_trial = conditionSetID_trial(:);

fprintf('\nUsing trial metadata from StimParams columns 26–30.\n');

%% ===================== AMPLITUDES ============================

trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

Amps = unique(trialAmps(:));
Amps = sort(Amps(:))';

all_nonzero_amps = Amps(Amps > 0);

if isempty(Plot_Amps)
    Plot_Amps_selected = all_nonzero_amps;
else
    Plot_Amps_selected = intersect(all_nonzero_amps, Plot_Amps);
end

%% ===================== FINAL ACTIVE ARTIFACT TIME =============

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials

    if conditionType_trial(tr) == 0
        lastActivePTD_us(tr) = 0;
        continue;
    end

    if conditionType_trial(tr) == 1
        nActive_this = prefixLength_trial(tr);
    elseif conditionType_trial(tr) == 2
        nActive_this = activeCount_trial(tr);
    else
        nActive_this = activeCount_trial(tr);
    end

    if isnan(nActive_this) || nActive_this < 1
        lastActivePTD_us(tr) = 0;
        continue;
    end

    nActive_this = min(round(nActive_this), simultaneous_stim);

    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);
    activeRows = rr(1:nActive_this);

    ptd_us = cell2mat(StimParams_data(activeRows,6));
    pulseNum = cell2mat(StimParams_data(activeRows,8));
    pulsePeriod_us = cell2mat(StimParams_data(activeRows,9));

    pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
    pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

    rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowFinalArtifact_us);
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

fprintf('Detected final active artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ===================== APPLY USER FILTERS =====================

% Set/order IDs.
SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));

if isempty(Plot_SetIDs)
    SetIDs_selected = SetIDs_all;
else
    SetIDs_selected = intersect(SetIDs_all, Plot_SetIDs);
end

% Condition types.
ConditionTypes_all = unique(conditionType_trial(conditionType_trial > 0));

if isempty(Plot_ConditionTypes)
    ConditionTypes_selected = ConditionTypes_all;
else
    ConditionTypes_selected = intersect(ConditionTypes_all, Plot_ConditionTypes);
end

% Prefixes.
Prefix_all = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));
Prefix_all = sort(Prefix_all(:))';

if isempty(Plot_Prefixes)
    Prefixes_selected = Prefix_all;
else
    Prefixes_selected = intersect(Prefix_all, Plot_Prefixes);
end

% ISIs.
ISI_all = unique(isi_ms_trial(conditionType_trial == 1));

if isempty(Plot_ISIs_ms)
    ISIs_selected = ISI_all;
else
    ISIs_selected = intersect(ISI_all, Plot_ISIs_ms);
end

fprintf('\nSelected set IDs: ');
disp(SetIDs_selected');

fprintf('Selected amplitudes: ');
disp(Plot_Amps_selected');

fprintf('Selected condition types: ');
disp(ConditionTypes_selected');

fprintf('Selected prefixes: ');
disp(Prefixes_selected);

fprintf('Selected ISIs (ms): ');
disp(ISIs_selected');

fprintf('Selected raster channels Depth_s index: ');
disp(raster_channels_to_plot);

fprintf('Responding file mode: %s\n', RespondingFileMode);
fprintf('Pool simultaneous across sets: %d\n', pool_simultaneous_across_sets);
fprintf('Remove bad trials: %d\n', remove_bad_trials);

%% ===================== ELECTRODE MAP =========================

d = Depth_s(Electrode_Type);

% Use manually selected channels instead of continuous start:end range.
depth_range = unique(raster_channels_to_plot(:)', 'stable');

% Remove channels outside Depth_s range.
depth_range = depth_range(depth_range >= 1 & depth_range <= numel(d));

if isempty(depth_range)
    error('No valid channels in raster_channels_to_plot. Check channel list.');
end

nChPlot = numel(depth_range);

fprintf('Valid selected channels after range check: ');
disp(depth_range);

%% ===================== PSTH KERNEL ===========================

edges = ras_win(1) : bin_ms_raster : ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
g = g / sum(g);

%% ========================================================================
%  MAIN CONDITION LOOPS
% ========================================================================

for si = 1:numel(SetIDs_selected)

    setID = SetIDs_selected(si);

    % ----- Build set label -----
    if ~isempty(CHN) && setID <= size(CHN,1)

        stimVec = CHN(setID,:);
        stimVec = stimVec(stimVec > 0);

        setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), ...
            stimVec, 'UniformOutput', false), '→');
    else
        setLabel = sprintf('Set%d', setID);
    end

    for ai = 1:numel(Plot_Amps_selected)

        ampVal = Plot_Amps_selected(ai);

        %% =============================================================
        %  Sequential prefix/recovery condition
        % =============================================================

        if ismember(1, ConditionTypes_selected)

            for pi = 1:numel(Prefixes_selected)

                prefixVal = Prefixes_selected(pi);

                for ii = 1:numel(ISIs_selected)

                    isiVal = ISIs_selected(ii);

                    trials_raw = find(conditionType_trial == 1 & ...
                                      conditionSetID_trial == setID & ...
                                      prefixLength_trial == prefixVal & ...
                                      abs(isi_ms_trial - isiVal) < 1e-6 & ...
                                      abs(trialAmps - ampVal) < 1e-6);

                    trials_raw = trials_raw(trials_raw <= nTrials_use);

                    if isempty(trials_raw)
                        continue;
                    end

                    trials_this = remove_bad_trials_from_condition( ...
                        trials_raw, BadTrials_global, remove_bad_trials, hasBadTrials);

                    nRaw = numel(trials_raw);
                    nClean = numel(trials_this);
                    nRemoved = nRaw - nClean;

                    if isempty(trials_this) && skip_if_no_clean_trials
                        fprintf('Skipping Set %d | Prefix %d | ISI %.1f | Amp %.1f: all %d trials removed.\n', ...
                            setID, prefixVal, isiVal, ampVal, nRaw);
                        continue;
                    end

                    finalArtifact_ms = max(lastActivePTD_ms(trials_raw));

                    % figTitle = sprintf('Set %d (%s) | Prefix %d | ISI %.1f ms | Amp %.1f uA | nRaw=%d | nClean=%d | removed=%d', ...
                    %     setID, setLabel, prefixVal, isiVal, ampVal, nRaw, nClean, nRemoved);

                    figTitle = sprintf('Set %d (%s) | Prefix %d | ISI %.1f ms | Amp %.1f uA', setID, setLabel, prefixVal, isiVal, ampVal);

                    figure('Color','w', ...
                           'Name', figTitle, ...
                           'Position', fig_position);

                    tiledlayout('flow', ...
                                'TileSpacing','compact', ...
                                'Padding','compact');

                    sgtitle(figTitle, ...
                            'FontSize', 14, ...
                            'FontWeight','bold', ...
                            'Interpreter','none');

                    % ----- Channel loop -----
                    for idxDepth = 1:nChPlot

                        ich = depth_range(idxDepth);
                        ch  = d(ich);

                        ax = nexttile;
                        hold(ax,'on');

                        if ch < 1 || ch > nCh || isempty(sp{ch}) || isempty(trials_this)
                            axis(ax, 'off');
                            continue;
                        end

                        sp_times = sp{ch}(:,1);

                        [allTrialSpikes, rate_s, yMaxPSTH] = ...
                            get_raster_psth_for_channel( ...
                                sp_times, trig, FS, trials_this, ras_win, ...
                                edges, bin_s, ctrs, g);

                        nTr = numel(trials_this);

                        % ----- Is this channel responsive? -----
                        isResp = false;

                        if highlight_responding_channels && hasResp
                            isResp = get_resp_status_prefix_by_value( ...
                                Resp, setID, ampVal, prefixVal, isiVal, ich);
                        end

                        if isResp
                            set(ax, 'Color', resp_bg_color);
                            ax.XColor = [0.6 0 0];
                            ax.YColor = [0.6 0 0];
                            ax.LineWidth = 1.4;
                        end

                        % ----- Plot PSTH -----
                        yyaxis(ax,'left');

                        if any(rate_s)
                            plot(ax, ctrs, rate_s, 'LineWidth', 1.4);
                        end

                        xlim(ax, ras_win);
                        ylim(ax, [0 yMaxPSTH]);
                        ylabel(ax, 'Rate');

                        % ----- Plot raster -----
                        yyaxis(ax,'right');

                        for ti = 1:nTr
                            tt = allTrialSpikes{ti};

                            if isempty(tt)
                                continue;
                            end

                            plot(ax, tt, ti*ones(size(tt)), '.', ...
                                'Color', [0 0 0], ...
                                'MarkerSize', 4);
                        end

                        ylim(ax, [0 nTr+1]);
                        set(ax, 'YTick', []);

                        % ----- Event lines -----
                        xline(ax, 0, 'r--', 'LineWidth', 1);

                        if isfinite(finalArtifact_ms) && finalArtifact_ms > 0
                            xline(ax, finalArtifact_ms, 'k--', 'LineWidth', 1);
                        end

                        xlim(ax, ras_win);

                        % ----- Title -----
                        if isResp
                            chLabel = sprintf('Ch %d RESP', ich);
                        else
                            chLabel = sprintf('Ch %d', ich);
                        end

                        ht = title(ax, chLabel, ...
                            'FontSize', 10, ...
                            'FontWeight', 'bold', ...
                            'Interpreter', 'none');

                        if isResp
                            set(ht, 'Color', [0.7 0 0]);
                        end

                        if idxDepth > nChPlot - ceil(sqrt(nChPlot))
                            xlabel(ax, 'Time (ms)');
                        end
                    end
                end
            end
        end

        %% =============================================================
        %  Simultaneous condition
        % =============================================================

        if ismember(2, ConditionTypes_selected)

            if pool_simultaneous_across_sets

                trials_raw = find(conditionType_trial == 2 & ...
                                  activeCount_trial == full_prefix_length & ...
                                  abs(trialAmps - ampVal) < 1e-6);

                simLabel = 'Simultaneous pooled';

                % For pooled simultaneous trials, highlighting is ambiguous.
                % The code will check each current setID, but if simultaneous
                % only exists in one set, use pool_simultaneous_across_sets = false
                % for cleaner interpretation.

            else

                trials_raw = find(conditionType_trial == 2 & ...
                                  conditionSetID_trial == setID & ...
                                  activeCount_trial == full_prefix_length & ...
                                  abs(trialAmps - ampVal) < 1e-6);

                simLabel = 'Simultaneous';
            end

            trials_raw = trials_raw(trials_raw <= nTrials_use);

            if isempty(trials_raw)
                continue;
            end

            trials_this = remove_bad_trials_from_condition( ...
                trials_raw, BadTrials_global, remove_bad_trials, hasBadTrials);

            nRaw = numel(trials_raw);
            nClean = numel(trials_this);
            nRemoved = nRaw - nClean;

            if isempty(trials_this) && skip_if_no_clean_trials
                fprintf('Skipping Set %d | %s | Amp %.1f: all %d trials removed.\n', ...
                    setID, simLabel, ampVal, nRaw);
                continue;
            end

            finalArtifact_ms = max(lastActivePTD_ms(trials_raw));

            % figTitle = sprintf('Set %d (%s) | %s | Amp %.1f uA | nRaw=%d | nClean=%d | removed=%d', ...
            %     setID, setLabel, simLabel, ampVal, nRaw, nClean, nRemoved);

            figTitle = sprintf('Set %d (%s) | %s | Amp %.1f uA', setID, setLabel, simLabel, ampVal);

            figure('Color','w', ...
                   'Name', figTitle, ...
                   'Position', fig_position);

            tiledlayout('flow', ...
                        'TileSpacing','compact', ...
                        'Padding','compact');

            sgtitle(figTitle, ...
                    'FontSize', 14, ...
                    'FontWeight','bold', ...
                    'Interpreter','none');

            % ----- Channel loop -----
            for idxDepth = 1:nChPlot

                ich = depth_range(idxDepth);
                ch  = d(ich);

                ax = nexttile;
                hold(ax,'on');

                if ch < 1 || ch > nCh || isempty(sp{ch}) || isempty(trials_this)
                    axis(ax, 'off');
                    continue;
                end

                sp_times = sp{ch}(:,1);

                [allTrialSpikes, rate_s, yMaxPSTH] = ...
                    get_raster_psth_for_channel( ...
                        sp_times, trig, FS, trials_this, ras_win, ...
                        edges, bin_s, ctrs, g);

                nTr = numel(trials_this);

                % ----- Is this channel responsive? -----
                isResp = false;

                if highlight_responding_channels && hasResp

                    if pool_simultaneous_across_sets
                        % Pooled simultaneous highlighting is less clean.
                        % This checks the current setID only.
                        respSetForSim = setID;
                    else
                        respSetForSim = setID;
                    end

                    isResp = get_resp_status_sim_by_value( ...
                        Resp, respSetForSim, ampVal, full_prefix_length, ich);
                end

                if isResp
                    set(ax, 'Color', resp_bg_color);
                    ax.XColor = [0.6 0 0];
                    ax.YColor = [0.6 0 0];
                    ax.LineWidth = 1.4;
                end

                % ----- Plot PSTH -----
                yyaxis(ax,'left');

                if any(rate_s)
                    plot(ax, ctrs, rate_s, 'LineWidth', 1.4);
                end

                xlim(ax, ras_win);
                ylim(ax, [0 yMaxPSTH]);
                ylabel(ax, 'Rate');

                % ----- Plot raster -----
                yyaxis(ax,'right');

                for ti = 1:nTr
                    tt = allTrialSpikes{ti};

                    if isempty(tt)
                        continue;
                    end

                    plot(ax, tt, ti*ones(size(tt)), '.', ...
                        'Color', [0 0 0], ...
                        'MarkerSize', 4);
                end

                ylim(ax, [0 nTr+1]);
                set(ax, 'YTick', []);

                % ----- Event lines -----
                xline(ax, 0, 'r--', 'LineWidth', 1);

                if isfinite(finalArtifact_ms) && finalArtifact_ms > 0
                    xline(ax, finalArtifact_ms, 'k--', 'LineWidth', 1);
                end

                xlim(ax, ras_win);

                % ----- Title -----
                if isResp
                    chLabel = sprintf('Ch %d RESP', ich);
                else
                    chLabel = sprintf('Ch %d', ich);
                end

                ht = title(ax, chLabel, ...
                    'FontSize', 10, ...
                    'FontWeight', 'bold', ...
                    'Interpreter', 'none');

                if isResp
                    set(ht, 'Color', [0.7 0 0]);
                end

                if idxDepth > nChPlot - ceil(sqrt(nChPlot))
                    xlabel(ax, 'Time (ms)');
                end
            end
        end
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function BadTrials_global = collect_global_bad_trials(BadTrials)

    BadTrials_global = [];

    if isempty(BadTrials)
        return;
    end

    for ii = 1:numel(BadTrials)

        if isempty(BadTrials{ii})
            continue;
        end

        BadTrials_global = [BadTrials_global; BadTrials{ii}(:)]; %#ok<AGROW>
    end

    BadTrials_global = unique(BadTrials_global(:));
end

function trials_clean = remove_bad_trials_from_condition( ...
    trials_raw, BadTrials_global, remove_bad_trials, hasBadTrials)

    if remove_bad_trials && hasBadTrials && ~isempty(BadTrials_global)
        trials_clean = setdiff(trials_raw(:), BadTrials_global(:), 'stable');
    else
        trials_clean = trials_raw(:);
    end

    trials_clean = trials_clean(:)';
end

function [allTrialSpikes, rate_s, yMaxPSTH] = get_raster_psth_for_channel( ...
    sp_times, trig, FS, trials_this, ras_win, edges, bin_s, ctrs, g)

    allTrialSpikes = cell(numel(trials_this),1);
    counts = zeros(1, numel(edges)-1);

    for ti = 1:numel(trials_this)

        tr = trials_this(ti);
        t0 = trig(tr) / FS * 1000;

        tt = sp_times;
        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;

        allTrialSpikes{ti} = tt;
        counts = counts + histcounts(tt, edges);
    end

    if ~any(counts)
        rate_s = zeros(size(ctrs));
    else
        rate = counts / (numel(trials_this) * bin_s);
        rate_s = filter(g, 1, rate);
    end

    maxRate = max(rate_s);
    yMaxPSTH = max(50, ceil(maxRate * 1.1 / 10) * 10);
end

%% ===================== VALUE-BASED RESPONDING LOOKUP =====================

function isResp = get_resp_status_prefix_by_value( ...
    Resp, targetSetID, targetAmp, targetPrefix, targetISI, ich)

    isResp = false;

    if ~isfield(Resp, 'set') || isempty(Resp.set)
        return;
    end

    %% ----- Find correct set by setID -----
    si = find_resp_set_index(Resp, targetSetID);

    if si == 0
        return;
    end

    %% ----- Find correct amplitude by amp_value -----
    ai = find_resp_amp_index(Resp.set(si), targetAmp);

    if ai == 0
        return;
    end

    if ~isfield(Resp.set(si).amp(ai), 'prefix') || ...
       isempty(Resp.set(si).amp(ai).prefix)
        return;
    end

    %% ----- Find correct prefix by prefix length and ISI -----
    pi = find_resp_prefix_index( ...
        Resp.set(si).amp(ai).prefix, targetPrefix, targetISI);

    if pi == 0
        return;
    end

    if ~isfield(Resp.set(si).amp(ai).prefix(pi), 'channel')
        return;
    end

    if ich > numel(Resp.set(si).amp(ai).prefix(pi).channel)
        return;
    end

    R = Resp.set(si).amp(ai).prefix(pi).channel(ich);

    if isfield(R, 'is_responsive') && R.is_responsive
        isResp = true;
    end
end

function isResp = get_resp_status_sim_by_value( ...
    Resp, targetSetID, targetAmp, targetActiveCount, ich)

    isResp = false;

    if ~isfield(Resp, 'set') || isempty(Resp.set)
        return;
    end

    %% ----- Find correct set by setID -----
    si = find_resp_set_index(Resp, targetSetID);

    if si == 0
        return;
    end

    %% ----- Find correct amplitude by amp_value -----
    ai = find_resp_amp_index(Resp.set(si), targetAmp);

    if ai == 0
        return;
    end

    if ~isfield(Resp.set(si).amp(ai), 'sim') || ...
       ~isfield(Resp.set(si).amp(ai).sim, 'channel')
        return;
    end

    S = Resp.set(si).amp(ai).sim;

    %% ----- Check active count if available -----
    if isfield(S, 'active_count')
        if abs(S.active_count - targetActiveCount) > 1e-6
            return;
        end
    elseif isfield(S, 'activeCount')
        if abs(S.activeCount - targetActiveCount) > 1e-6
            return;
        end
    end

    if ich > numel(S.channel)
        return;
    end

    R = S.channel(ich);

    if isfield(R, 'is_responsive') && R.is_responsive
        isResp = true;
    end
end

function si = find_resp_set_index(Resp, targetSetID)

    si = 0;

    if ~isfield(Resp, 'set') || isempty(Resp.set)
        return;
    end

    for k = 1:numel(Resp.set)

        if isfield(Resp.set(k), 'setID')
            if Resp.set(k).setID == targetSetID
                si = k;
                return;
            end
        else
            % Fallback for older structures.
            if k == targetSetID
                si = k;
                return;
            end
        end
    end
end

function ai = find_resp_amp_index(RespSet, targetAmp)

    ai = 0;

    if ~isfield(RespSet, 'amp') || isempty(RespSet.amp)
        return;
    end

    for k = 1:numel(RespSet.amp)

        if isfield(RespSet.amp(k), 'amp_value')
            ampVal = RespSet.amp(k).amp_value;
        elseif isfield(RespSet.amp(k), 'ampVal')
            ampVal = RespSet.amp(k).ampVal;
        elseif isfield(RespSet.amp(k), 'amplitude')
            ampVal = RespSet.amp(k).amplitude;
        else
            continue;
        end

        if abs(ampVal - targetAmp) < 1e-6
            ai = k;
            return;
        end
    end
end

function pi = find_resp_prefix_index(prefixStruct, targetPrefix, targetISI)

    pi = 0;

    for k = 1:numel(prefixStruct)

        P = prefixStruct(k);

        %% ----- Prefix length -----
        if isfield(P, 'prefix_length')
            prefixVal = P.prefix_length;
        elseif isfield(P, 'prefixLength')
            prefixVal = P.prefixLength;
        else
            prefixVal = NaN;
        end

        if isnan(prefixVal) || abs(prefixVal - targetPrefix) > 1e-6
            continue;
        end

        %% ----- ISI -----
        if isfield(P, 'isi_ms')
            isiVal = P.isi_ms;
        elseif isfield(P, 'ISI_ms')
            isiVal = P.ISI_ms;
        else
            % If ISI was not stored in Responding, accept the prefix match.
            isiVal = targetISI;
        end

        if abs(isiVal - targetISI) > 1e-6
            continue;
        end

        pi = k;
        return;
    end
end