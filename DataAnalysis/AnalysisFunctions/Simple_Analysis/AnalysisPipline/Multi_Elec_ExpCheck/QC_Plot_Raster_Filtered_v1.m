%% ============================================================
%  Raster + PSTH Viewer for SSD-Filtered Mixed-Prefix Stimulation Files
%  Purpose:
%    Plot raster + PSTH using SSD-filtered spike data.
%  Input spike file:
%    *.sp_xia_SSD.mat
%  Required spike variable:
%    sp_corr
%  ConditionType:
%    0 = zero-control
%    1 = sequential prefix/recovery
%    2 = full simultaneous
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== CHOOSE FOLDER ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

%% ====================== USER SETTINGS ======================

raster_chn_start = 1;
raster_chn_end   = 32;   % Depth_s index
Electrode_Type   = 2;    % 0: rigid; 1: single-shank flex; 2: four-shank flex

% Which prefixes to plot.
% [] means all detected prefixes.
ActiveCounts_to_plot = [5];

% Which ISIs to plot, in ms.
% [] means all detected ISIs.
ISI_to_plot_ms = [5];

% Which condition types to plot.
%   1 = sequential prefix/recovery
%   2 = full simultaneous
% [] means all detected condition types except zero-control.
ConditionTypes_to_plot = [];

% Which stimulation set/order IDs to plot.
% [] means all detected sets.
SetIDs_to_plot = [1];

% Amplitudes to plot.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% Raster/PSTH parameters.
ras_win        = [-20 100];   % ms
bin_ms_raster = 1;           % PSTH bin size, ms
smooth_ms     = 2;           % PSTH smoothing width

% If true, print one example trial for each prefix.
debug_print_trial_content = true;

% If true, simultaneous trials are pooled across all set/order IDs.
% Use this only if simultaneous stimulation does not depend on electrode order.
pool_simultaneous_across_sets = true;

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

%% ====================== SAMPLING RATE ======================
FS = 30000;

%% ====================== LOAD SSD-FILTERED SPIKES ======================
% This code expects the SSD-filtered file to contain variable sp_corr.

ssd_files = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
assert(~isempty(ssd_files), 'No *.sp_xia_SSD.mat file found in the current folder.');

ssd_filename = fullfile(data_folder, ssd_files(1).name);
fprintf('Loading SSD-filtered spike file: %s\n', ssd_filename);

S_sp = load(ssd_filename);

if isfield(S_sp, 'sp_corr')
    sp_plot = S_sp.sp_corr;
else
    error('Variable "sp_corr" not found in %s.', ssd_filename);
end

fprintf('Loaded SSD-filtered spikes from: %s\n', ssd_files(1).name);

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

%% ====================== LOAD STIM PARAMS ======================
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;   % rows/slots per trial
n_Trials          = S_exp.n_Trials;
E_MAP             = S_exp.E_MAP;

if isfield(S_exp, 'CHN')
    CHN = S_exp.CHN;
else
    CHN = [];
end

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match trigger number (%d). Using min of both for plotting.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== REMOVE HEADER ROW ======================
StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;
if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

%% ====================== TRIAL METADATA FROM STIMPARAMS ======================
% Important:
%   Read trial metadata directly from randomized StimParams columns 26–30.
%
% StimParams columns:
%   26 = ActiveElectrodeCount
%   27 = PrefixLength
%   28 = ISI_ms
%   29 = ConditionType
%   30 = ConditionSetID

if size(StimParams_data,2) < 30
    error('StimParams does not contain columns 26–30. Cannot use mixed-prefix metadata.');
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

%% ====================== AMPLITUDES ======================
trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0 before unique().
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

[Amps, ~, ampIdx] = unique(trialAmps(:));
n_AMP = numel(Amps);
cmap = lines(n_AMP);

%% ====================== PULSE TRAIN PERIOD ======================
% Kept mainly for compatibility.
pulseTrain_all = cell2mat(StimParams_data(:,9));
pulseTrain = pulseTrain_all(firstRow_eachTrial);
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

%% ====================== LAST ACTIVE ARTIFACT TIME ======================
% Calculate final artifact time for each trial.
%
% For each active row:
%   final artifact time =
%       PTD_us + (PulseNum - 1) * PulsePeriod_us
%
% Then for the whole trial:
%   lastActivePTD_us(tr) = max(final artifact time across active rows)

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

    rowLastArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;
    lastActivePTD_us(tr) = max(rowLastArtifact_us);
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

%% ====================== APPLY USER FILTERS ======================

% Set/order IDs.
SetIDs_all = unique(conditionSetID_trial(conditionSetID_trial > 0));

if isempty(SetIDs_to_plot)
    SetIDs_selected = SetIDs_all;
else
    SetIDs_selected = intersect(SetIDs_all, SetIDs_to_plot);
end

% Prefixes.
prefix_all = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));

if isempty(ActiveCounts_to_plot)
    Prefixes_selected = prefix_all;
else
    Prefixes_selected = intersect(prefix_all, ActiveCounts_to_plot);
end

% ISIs for sequential prefix trials.
isi_all = unique(isi_ms_trial(conditionType_trial == 1));

if isempty(ISI_to_plot_ms)
    ISIs_selected = isi_all;
else
    ISIs_selected = intersect(isi_all, ISI_to_plot_ms);
end

% Condition types.
conditionTypes_all = unique(conditionType_trial(conditionType_trial > 0));

if isempty(ConditionTypes_to_plot)
    ConditionTypes_selected = conditionTypes_all;
else
    ConditionTypes_selected = intersect(conditionTypes_all, ConditionTypes_to_plot);
end

% Amplitudes.
all_nonzero_amps = Amps(Amps > 0);

if isempty(Amps_to_plot)
    Amps_selected = all_nonzero_amps;
else
    Amps_selected = intersect(all_nonzero_amps, Amps_to_plot);
end

fprintf('\nDetected amplitudes: ');
disp(Amps');

fprintf('Selected amplitudes: ');
disp(Amps_selected');

fprintf('\nDetected set IDs: ');
disp(SetIDs_all');

fprintf('Selected set IDs: ');
disp(SetIDs_selected');

fprintf('\nDetected prefixes: ');
disp(prefix_all');

fprintf('Selected prefixes: ');
disp(Prefixes_selected');

fprintf('\nDetected ISIs (ms): ');
disp(isi_all');

fprintf('Selected ISIs (ms): ');
disp(ISIs_selected');

fprintf('\nDetected condition types: ');
disp(unique(conditionType_trial)');

fprintf('Selected condition types: ');
disp(ConditionTypes_selected');

fprintf('\nDetected final active artifact times (us): ');
disp(unique(lastActivePTD_us)');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================
if debug_print_trial_content && ~isempty(SetIDs_selected) && ~isempty(Amps_selected) && ~isempty(ISIs_selected)

    fprintf('\n================ DEBUG TRIAL CONTENT CHECK ================\n');

    debug_set = SetIDs_selected(1);
    debug_amp = Amps_selected(1);
    debug_isi = ISIs_selected(1);

    fprintf('Debug Set = %d | Amp = %.1f uA | ISI = %.1f ms\n', ...
        debug_set, debug_amp, debug_isi);

    for ip = 1:numel(Prefixes_selected)

        prefix_val = Prefixes_selected(ip);

        tr_debug = find(conditionSetID_trial == debug_set & ...
                        conditionType_trial == 1 & ...
                        prefixLength_trial == prefix_val & ...
                        isi_ms_trial == debug_isi & ...
                        trialAmps == debug_amp, ...
                        1, 'first');

        if isempty(tr_debug)
            fprintf('\nPrefix %d: no matching trial found.\n', prefix_val);
            continue;
        end

        rr = (tr_debug-1)*simultaneous_stim + (1:simultaneous_stim);

        stimNames_debug = StimParams_data(rr,1);
        ptd_debug = cell2mat(StimParams_data(rr,6));
        pulseNum_debug = cell2mat(StimParams_data(rr,8));
        pulsePeriod_debug = cell2mat(StimParams_data(rr,9));
        amp_debug = cell2mat(StimParams_data(rr,16));
        activeCount_debug = cell2mat(StimParams_data(rr,26));
        prefix_debug = cell2mat(StimParams_data(rr,27));
        isi_debug = cell2mat(StimParams_data(rr,28));
        condType_debug = cell2mat(StimParams_data(rr,29));
        setID_debug = cell2mat(StimParams_data(rr,30));

        fprintf('\nPrefix %d | Trial %d\n', prefix_val, tr_debug);
        fprintf('  conditionType = %d, setID = %d, ISI = %.1f ms, amp = %.1f uA\n', ...
            conditionType_trial(tr_debug), conditionSetID_trial(tr_debug), ...
            isi_ms_trial(tr_debug), trialAmps(tr_debug));

        fprintf('  StimNames:       ');
        disp(stimNames_debug');

        fprintf('  PTD_us:          ');
        disp(ptd_debug');

        fprintf('  PulseNum:        ');
        disp(pulseNum_debug');

        fprintf('  PulsePeriod_us:  ');
        disp(pulsePeriod_debug');

        fprintf('  Amp_col16:       ');
        disp(amp_debug');

        fprintf('  ActiveCount_col26: ');
        disp(activeCount_debug');

        fprintf('  Prefix_col27:      ');
        disp(prefix_debug');

        fprintf('  ISI_col28:         ');
        disp(isi_debug');

        fprintf('  CondType_col29:    ');
        disp(condType_debug');

        fprintf('  SetID_col30:       ');
        disp(setID_debug');

        fprintf('  FinalArtifact_us = %.1f\n', lastActivePTD_us(tr_debug));
    end

    fprintf('===========================================================\n');
end

%% ====================== ELECTRODE MAP ======================
d = Depth_s(Electrode_Type);

%% ====================== RASTER / PSTH SETUP ======================
edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

%% ====================== MAIN RASTER + PSTH LOOP ======================
% Plot structure:
%   Channel
%     Set/order
%       Condition type
%         Prefix
%           ISI
%             Amplitude

for ich = raster_chn_start:raster_chn_end

    ch = d(ich);

    if ch > numel(sp_plot)
        continue;
    end

    if isempty(sp_plot{ch})
        continue;
    end

    for si = 1:numel(SetIDs_selected)

        setID = SetIDs_selected(si);

        %% ----- Build set label -----
        if ~isempty(CHN) && setID <= size(CHN,1)
            stimVec = CHN(setID,:);
            stimVec = stimVec(stimVec > 0);
            setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimVec, 'UniformOutput', false), '→');
        else
            setLabel = sprintf('Set%d', setID);
        end

        for ci = 1:numel(ConditionTypes_selected)

            condType = ConditionTypes_selected(ci);

            %% =====================================================
            % Sequential prefix/recovery condition
            % ======================================================
            if condType == 1

                for pi = 1:n_PULSE

                    pulse_val = PulsePeriods(pi);

                    for aci = 1:numel(Prefixes_selected)

                        prefixVal = Prefixes_selected(aci);

                        for ii = 1:numel(ISIs_selected)

                            isi_val = ISIs_selected(ii);

                            %% ----- Check whether group has data -----
                            group_trials = find(conditionSetID_trial == setID & ...
                                                conditionType_trial == 1 & ...
                                                prefixLength_trial == prefixVal & ...
                                                isi_ms_trial == isi_val & ...
                                                pulseIdx == pi & ...
                                                ismember(trialAmps, Amps_selected));

                            group_trials = group_trials(group_trials <= nTrials_use);

                            if isempty(group_trials)
                                continue;
                            end

                            %% ----- Figure -----
                            figName = sprintf('SSD | Ch %d | Set %d %s | Prefix %d | ISI %.1f ms | Pulse %d us', ...
                                              ich, setID, setLabel, prefixVal, isi_val, pulse_val);

                            figure('Color','w','Name',figName);
                            tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                            ax1 = nexttile([3 1]);
                            hold(ax1,'on'); box(ax1,'off');

                            title(ax1, sprintf('SSD Raster — Ch %d | Set %d %s | Prefix %d | ISI %.1f ms', ...
                                               ich, setID, setLabel, prefixVal, isi_val), ...
                                               'Interpreter','none');

                            ax2 = nexttile;
                            hold(ax2,'on'); box(ax2,'off');

                            %% ----- PSTH storage -----
                            psth_curves = cell(1, numel(Amps_selected));
                            maxRate = 0;
                            y_cursor = 0;
                            ytick_vals = [];
                            ytick_labels = {};
                            legend_labels = cell(1, numel(Amps_selected));

                            %% ----- Loop amplitudes -----
                            for aa = 1:numel(Amps_selected)

                                amp_val = Amps_selected(aa);
                                amp_idx = find(Amps == amp_val, 1, 'first');

                                if isempty(amp_idx)
                                    continue;
                                end

                                color = cmap(amp_idx,:);

                                amp_trials = find(ampIdx == amp_idx & ...
                                                  pulseIdx == pi & ...
                                                  conditionSetID_trial == setID & ...
                                                  conditionType_trial == 1 & ...
                                                  prefixLength_trial == prefixVal & ...
                                                  isi_ms_trial == isi_val);

                                amp_trials = amp_trials(amp_trials <= nTrials_use);

                                nTr = numel(amp_trials);

                                if nTr == 0
                                    psth_curves{aa} = zeros(1, numel(ctrs));
                                    legend_labels{aa} = sprintf('%g uA', amp_val);
                                    continue;
                                end

                                counts = zeros(1, numel(edges)-1);

                                %% ----- Trial loop -----
                                for t = 1:nTr

                                    tr = amp_trials(t);
                                    t0 = trig(tr) / FS * 1000;

                                    tt = sp_plot{ch}(:,1);
                                    tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;

                                    % Raster.
                                    y0 = y_cursor + t;
                                    for spike_t = tt'
                                        plot(ax1, [spike_t spike_t], [y0-0.4 y0+0.4], ...
                                             'Color', color, 'LineWidth', 1.1);
                                    end

                                    % PSTH.
                                    counts = counts + histcounts(tt, edges);
                                end

                                ytick_vals(end+1) = y_cursor + nTr/2;
                                ytick_labels{end+1} = sprintf('%g uA', amp_val);
                                y_cursor = y_cursor + nTr;

                                rate = filter(g, 1, counts/(nTr*bin_s));
                                psth_curves{aa} = rate;
                                maxRate = max(maxRate, max(rate));

                                legend_labels{aa} = sprintf('%g uA', amp_val);
                            end

                            %% ----- Final artifact line -----
                            finalArtifact_trials = lastActivePTD_ms(group_trials);
                            finalArtifact_ms = max(finalArtifact_trials);

                            %% ----- Finalize raster axis -----
                            xline(ax1, 0, 'r--');

                            if isfinite(finalArtifact_ms) && finalArtifact_ms > 0
                                xline(ax1, finalArtifact_ms, 'k--', 'LineWidth', 1);
                            end

                            xlim(ax1, ras_win);
                            ylim(ax1, [0 max(1,y_cursor)]);
                            yticks(ax1, ytick_vals);
                            yticklabels(ax1, ytick_labels);
                            ylabel(ax1, 'Amplitude');

                            %% ----- Finalize PSTH -----
                            for aa = 1:numel(Amps_selected)

                                amp_val = Amps_selected(aa);
                                amp_idx = find(Amps == amp_val, 1, 'first');

                                if isempty(amp_idx) || isempty(psth_curves{aa})
                                    continue;
                                end

                                plot(ax2, ctrs, psth_curves{aa}, ...
                                     'Color', cmap(amp_idx,:), ...
                                     'LineWidth', 1.6);
                            end

                            xline(ax2, 0, 'r--');

                            if isfinite(finalArtifact_ms) && finalArtifact_ms > 0
                                xline(ax2, finalArtifact_ms, 'k--', 'LineWidth', 1);
                            end

                            xlim(ax2, ras_win);
                            ylim(ax2, [0 max(50, ceil(maxRate*1.1/10)*10)]);
                            xlabel(ax2, 'Time (ms)');
                            ylabel(ax2, 'Rate (sp/s)');
                            legend(ax2, legend_labels, 'Box','off','Location','northeast');

                        end % ISI
                    end % prefix
                end % pulse

            %% =====================================================
            % Full simultaneous condition
            % ======================================================
            elseif condType == 2

                for pi = 1:n_PULSE

                    pulse_val = PulsePeriods(pi);

                    activeCounts_sim = unique(activeCount_trial(conditionType_trial == 2));
                    activeCounts_sim = activeCounts_sim(activeCounts_sim > 0);

                    for aci = 1:numel(activeCounts_sim)

                        activeCount = activeCounts_sim(aci);

                        if pool_simultaneous_across_sets
                            group_trials = find(conditionType_trial == 2 & ...
                                                activeCount_trial == activeCount & ...
                                                pulseIdx == pi & ...
                                                ismember(trialAmps, Amps_selected));
                        else
                            group_trials = find(conditionSetID_trial == setID & ...
                                                conditionType_trial == 2 & ...
                                                activeCount_trial == activeCount & ...
                                                pulseIdx == pi & ...
                                                ismember(trialAmps, Amps_selected));
                        end

                        group_trials = group_trials(group_trials <= nTrials_use);

                        if isempty(group_trials)
                            continue;
                        end

                        %% ----- Figure -----
                        if pool_simultaneous_across_sets
                            figName = sprintf('SSD | Ch %d | Set %d %s | Simultaneous pooled | Pulse %d us', ...
                                              ich, setID, setLabel, pulse_val);
                        else
                            figName = sprintf('SSD | Ch %d | Set %d %s | Simultaneous | Pulse %d us', ...
                                              ich, setID, setLabel, pulse_val);
                        end

                        figure('Color','w','Name',figName);
                        tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                        ax1 = nexttile([3 1]);
                        hold(ax1,'on'); box(ax1,'off');

                        if pool_simultaneous_across_sets
                            title(ax1, sprintf('SSD Raster — Ch %d | Set %d %s | Simultaneous pooled', ...
                                               ich, setID, setLabel), ...
                                               'Interpreter','none');
                        else
                            title(ax1, sprintf('SSD Raster — Ch %d | Set %d %s | Simultaneous', ...
                                               ich, setID, setLabel), ...
                                               'Interpreter','none');
                        end

                        ax2 = nexttile;
                        hold(ax2,'on'); box(ax2,'off');

                        %% ----- PSTH storage -----
                        psth_curves = cell(1, numel(Amps_selected));
                        maxRate = 0;
                        y_cursor = 0;
                        ytick_vals = [];
                        ytick_labels = {};
                        legend_labels = cell(1, numel(Amps_selected));

                        %% ----- Loop amplitudes -----
                        for aa = 1:numel(Amps_selected)

                            amp_val = Amps_selected(aa);
                            amp_idx = find(Amps == amp_val, 1, 'first');

                            if isempty(amp_idx)
                                continue;
                            end

                            color = cmap(amp_idx,:);

                            if pool_simultaneous_across_sets
                                amp_trials = find(ampIdx == amp_idx & ...
                                                  pulseIdx == pi & ...
                                                  conditionType_trial == 2 & ...
                                                  activeCount_trial == activeCount);
                            else
                                amp_trials = find(ampIdx == amp_idx & ...
                                                  pulseIdx == pi & ...
                                                  conditionSetID_trial == setID & ...
                                                  conditionType_trial == 2 & ...
                                                  activeCount_trial == activeCount);
                            end

                            amp_trials = amp_trials(amp_trials <= nTrials_use);

                            nTr = numel(amp_trials);

                            if nTr == 0
                                psth_curves{aa} = zeros(1, numel(ctrs));
                                legend_labels{aa} = sprintf('%g uA', amp_val);
                                continue;
                            end

                            counts = zeros(1, numel(edges)-1);

                            %% ----- Trial loop -----
                            for t = 1:nTr

                                tr = amp_trials(t);
                                t0 = trig(tr) / FS * 1000;

                                tt = sp_plot{ch}(:,1);
                                tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;

                                % Raster.
                                y0 = y_cursor + t;
                                for spike_t = tt'
                                    plot(ax1, [spike_t spike_t], [y0-0.4 y0+0.4], ...
                                         'Color', color, 'LineWidth', 1.1);
                                end

                                % PSTH.
                                counts = counts + histcounts(tt, edges);
                            end

                            ytick_vals(end+1) = y_cursor + nTr/2;
                            ytick_labels{end+1} = sprintf('%g uA', amp_val);
                            y_cursor = y_cursor + nTr;

                            rate = filter(g, 1, counts/(nTr*bin_s));
                            psth_curves{aa} = rate;
                            maxRate = max(maxRate, max(rate));

                            legend_labels{aa} = sprintf('%g uA', amp_val);
                        end

                        %% ----- Final artifact line -----
                        finalArtifact_trials = lastActivePTD_ms(group_trials);
                        finalArtifact_ms = max(finalArtifact_trials);

                        %% ----- Finalize raster -----
                        xline(ax1, 0, 'r--');

                        if isfinite(finalArtifact_ms) && finalArtifact_ms > 0
                            xline(ax1, finalArtifact_ms, 'k--', 'LineWidth', 1);
                        end

                        xlim(ax1, ras_win);
                        ylim(ax1, [0 max(1,y_cursor)]);
                        yticks(ax1, ytick_vals);
                        yticklabels(ax1, ytick_labels);
                        ylabel(ax1, 'Amplitude');

                        %% ----- Finalize PSTH -----
                        for aa = 1:numel(Amps_selected)

                            amp_val = Amps_selected(aa);
                            amp_idx = find(Amps == amp_val, 1, 'first');

                            if isempty(amp_idx) || isempty(psth_curves{aa})
                                continue;
                            end

                            plot(ax2, ctrs, psth_curves{aa}, ...
                                 'Color', cmap(amp_idx,:), ...
                                 'LineWidth', 1.6);
                        end

                        xline(ax2, 0, 'r--');

                        if isfinite(finalArtifact_ms) && finalArtifact_ms > 0
                            xline(ax2, finalArtifact_ms, 'k--', 'LineWidth', 1);
                        end

                        xlim(ax2, ras_win);
                        ylim(ax2, [0 max(50, ceil(maxRate*1.1/10)*10)]);
                        xlabel(ax2, 'Time (ms)');
                        ylabel(ax2, 'Rate (sp/s)');
                        legend(ax2, legend_labels, 'Box','off','Location','northeast');

                    end % activeCount
                end % pulse
            end % condition type
        end % condition type
    end % set
end % channel