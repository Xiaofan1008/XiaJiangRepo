%% ============================================================
%  Waveform Viewer for SSD-Filtered Mixed-Prefix Stimulation Files
%
%  Purpose:
%    Plot spike waveforms by recording channel so that bad channels can be
%    identified and excluded.
%
%  Input spike file:
%    *.sp_xia_SSD.mat
%
%  Required spike variable:
%    sp_corr
%
%  Grouping:
%    Channel
%      Set/order ID
%        Condition type
%          Prefix length / active count
%          ISI
%          Amplitude
%
%  ConditionType:
%    0 = zero-control
%    1 = sequential prefix/recovery
%    2 = full simultaneous
%
%  Important:
%    Trial metadata are read directly from randomized StimParams columns
%    26–30. This avoids trial-order mismatch.
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

% Recording channels to plot, using Depth_s index.
spike_chn_start = 1;
spike_chn_end   = 32;

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Which condition types to plot:
%   1 = sequential prefix/recovery
%   2 = simultaneous
% [] means plot both detected non-zero condition types.
ConditionTypes_to_plot = [1];

% Sequential prefix lengths to plot.
% [] means all detected prefixes.
PrefixLengths_to_plot = [5];

% ISIs to plot, in ms.
% [] means all detected ISIs.
ISI_to_plot_ms = [5];

% Stimulation set/order IDs to plot.
% [] means all detected sets.
SetIDs_to_plot = [3];

% Amplitudes to plot.
% [] means all non-zero amplitudes.
Amps_to_plot = [];

% For simultaneous condition:
% If true, simultaneous trials are pooled across all set/order IDs because
% simultaneous stimulation does not depend on electrode order.
pool_simultaneous_across_sets = true;

% Time window after trigger used to collect waveforms.
waveform_time_window_ms = [0 40];

% Bin width for grouping waveforms by spike latency.
bin_ms = 2;

% Optional waveform amplitude threshold for plotting.
% Waveforms with any point larger than this absolute value are not plotted.
% Set Inf to disable.
amp_threshold = 500;

% If true, align each waveform to its trough before plotting.
% This helps inspect waveform shape independent of small timing jitter.
align_waveforms_to_trough = true;

% Plot individual waveforms and mean waveform.
plot_individual_waveforms = true;
plot_mean_waveform = true;

% To avoid very slow/crowded plots, randomly limit the number of waveforms
% per bin and amplitude.
max_waveforms_per_bin_amp = 200;

% Figure layout.
layout_row = 4;
layout_col = 5;

% Sampling frequency.
FS = 30000;

% Print one example trial per prefix.
debug_print_trial_content = true;

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

%% ====================== LOAD SSD-FILTERED SPIKES ======================

ssd_files = dir(fullfile(data_folder, '*.sp_xia_SSD.mat'));
assert(~isempty(ssd_files), 'No *.sp_xia_SSD.mat file found in the current folder.');

ssd_filename = fullfile(data_folder, ssd_files(1).name);
fprintf('Loading SSD-filtered spike file:\n%s\n', ssd_filename);

S_sp = load(ssd_filename);

if isfield(S_sp, 'sp_corr')
    sp_use = S_sp.sp_corr;
else
    error('Variable "sp_corr" not found in %s.', ssd_filename);
end

nCh = numel(sp_use);

fprintf('Loaded SSD-filtered spikes from: %s\n', ssd_files(1).name);
fprintf('Number of spike channels: %d\n', nCh);

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

%% ====================== LOAD EXPERIMENT PARAMETERS ======================

fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

S_exp = load(fullfile(data_folder, fileDIR(1).name));

StimParams        = S_exp.StimParams;
simultaneous_stim = S_exp.simultaneous_stim;
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
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
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

%% ====================== AMPLITUDE PER TRIAL ======================

trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

[Amps, ~, ampIdx] = unique(trialAmps(:));
n_AMP = numel(Amps);
cmap = lines(n_AMP);

%% ====================== FINAL ACTIVE ARTIFACT TIME ======================

% For each active row:
%   final artifact time =
%       PTD_us + (PulseNum - 1) * PulsePeriod_us
%
% Then for each trial:
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

    rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowFinalArtifact_us);
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

% Prefix lengths for sequential trials.
Prefix_all = unique(prefixLength_trial(conditionType_trial == 1 & prefixLength_trial > 0));
Prefix_all = sort(Prefix_all(:))';

if isempty(PrefixLengths_to_plot)
    Prefixes_selected = Prefix_all;
else
    Prefixes_selected = intersect(Prefix_all, PrefixLengths_to_plot);
end

% ISIs for sequential prefix trials.
ISI_all = unique(isi_ms_trial(conditionType_trial == 1));

if isempty(ISI_to_plot_ms)
    ISIs_selected = ISI_all;
else
    ISIs_selected = intersect(ISI_all, ISI_to_plot_ms);
end

% Condition types.
ConditionTypes_all = unique(conditionType_trial(conditionType_trial > 0));

if isempty(ConditionTypes_to_plot)
    ConditionTypes_selected = ConditionTypes_all;
else
    ConditionTypes_selected = intersect(ConditionTypes_all, ConditionTypes_to_plot);
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
disp(Prefix_all');

fprintf('Selected prefixes: ');
disp(Prefixes_selected');

fprintf('\nDetected ISIs (ms): ');
disp(ISI_all');

fprintf('Selected ISIs (ms): ');
disp(ISIs_selected');

fprintf('\nDetected condition types: ');
disp(unique(conditionType_trial)');

fprintf('Selected condition types: ');
disp(ConditionTypes_selected');

fprintf('\nDetected final active artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================

if debug_print_trial_content && ~isempty(SetIDs_selected) && ...
        ~isempty(Amps_selected) && ~isempty(ISIs_selected)

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

%% ====================== WAVEFORM TIME BINS ======================

bin_edges = waveform_time_window_ms(1):bin_ms:waveform_time_window_ms(2);
nBins = numel(bin_edges) - 1;

if nBins > layout_row * layout_col
    warning('Number of time bins (%d) is larger than subplot layout capacity (%d). Some bins will not fit clearly.', ...
        nBins, layout_row * layout_col);
end

%% ====================== WAVEFORM X AXIS ======================

example_ch = find(~cellfun(@isempty, sp_use), 1, 'first');

if isempty(example_ch)
    error('All channels are empty in sp_use.');
end

wf_len = size(sp_use{example_ch}, 2) - 1;
t_wave = (0:wf_len-1) / FS * 1000;

%% ====================== MAIN WAVEFORM PLOTTING LOOP ======================

for ich = spike_chn_start:spike_chn_end

    ch = d(ich);

    if ch > nCh
        continue;
    end

    if isempty(sp_use{ch})
        continue;
    end

    sp_times_all = sp_use{ch}(:,1);
    sp_wave_all  = sp_use{ch}(:,2:end);

    % Optional amplitude threshold for plotting.
    valid_idx = all(abs(sp_wave_all) <= amp_threshold, 2);

    sp_times_all = sp_times_all(valid_idx);
    sp_wave_all  = sp_wave_all(valid_idx,:);

    if isempty(sp_times_all)
        continue;
    end

    for si = 1:numel(SetIDs_selected)

        setID = SetIDs_selected(si);

        %% ----- Build set label -----
        if ~isempty(CHN) && setID <= size(CHN,1)
            stimVec = CHN(setID,:);
            stimVec = stimVec(stimVec > 0);
            setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), ...
                         stimVec, 'UniformOutput', false), '→');
        else
            setLabel = sprintf('Set%d', setID);
        end

        for ci = 1:numel(ConditionTypes_selected)

            condType = ConditionTypes_selected(ci);

            %% =====================================================
            % Sequential prefix/recovery condition
            % =====================================================
            if condType == 1

                for pi = 1:numel(Prefixes_selected)

                    prefixVal = Prefixes_selected(pi);

                    for ii = 1:numel(ISIs_selected)

                        isi_val = ISIs_selected(ii);

                        %% ----- Trial group -----
                        group_trials = find(conditionSetID_trial == setID & ...
                                            conditionType_trial == 1 & ...
                                            prefixLength_trial == prefixVal & ...
                                            isi_ms_trial == isi_val & ...
                                            ismember(trialAmps, Amps_selected));

                        group_trials = group_trials(group_trials <= nTrials_use);

                        if isempty(group_trials)
                            continue;
                        end

                        %% ----- Collect waveforms by bin and amplitude -----
                        all_spikes_by_bin_amp = cell(nBins, numel(Amps_selected));

                        for tt = 1:numel(group_trials)

                            tr = group_trials(tt);
                            t0_ms = trig(tr) / FS * 1000;

                            rel_times_all = sp_times_all - t0_ms;

                            in_win = rel_times_all >= waveform_time_window_ms(1) & ...
                                     rel_times_all <  waveform_time_window_ms(2);

                            if ~any(in_win)
                                continue;
                            end

                            rel_times = rel_times_all(in_win);
                            waveforms = sp_wave_all(in_win,:);

                            amp_val_trial = trialAmps(tr);
                            amp_pos = find(Amps_selected == amp_val_trial, 1, 'first');

                            if isempty(amp_pos)
                                continue;
                            end

                            for j = 1:numel(rel_times)

                                bin_idx = find(rel_times(j) >= bin_edges(1:end-1) & ...
                                               rel_times(j) <  bin_edges(2:end), ...
                                               1, 'first');

                                if isempty(bin_idx)
                                    continue;
                                end

                                all_spikes_by_bin_amp{bin_idx, amp_pos}(end+1,:) = waveforms(j,:);
                            end
                        end

                        %% ----- Skip if no waveforms -----
                        all_waves = cell2mat(all_spikes_by_bin_amp(:));

                        if isempty(all_waves)
                            continue;
                        end

                        %% ----- Y limit -----
                        y_max = max(abs(all_waves(:)));

                        if isempty(y_max) || y_max == 0 || isnan(y_max)
                            y_lim = [-100 100];
                        else
                            y_lim = [-1 1] * ceil(y_max/50)*50;
                        end

                        finalArtifact_trials = lastActivePTD_ms(group_trials);
                        finalArtifact_ms = max(finalArtifact_trials);

                        %% ----- Figure -----
                        figName = sprintf('SSD Waveforms | Ch%d | Set%d | Prefix%d | ISI%gms', ...
                            ich, setID, prefixVal, isi_val);

                        figure('Name', figName, ...
                               'Color', 'w', ...
                               'Position', [100 100 1400 800]);

                        tiledlayout(layout_row, layout_col, ...
                                    'Padding', 'compact', ...
                                    'TileSpacing', 'compact');

                        sgtitle(sprintf('SSD Waveforms | Rec Ch %d | Set %d: %s | Prefix %d | ISI %.1f ms | Final artifact %.1f ms', ...
                            ich, setID, setLabel, prefixVal, isi_val, finalArtifact_ms), ...
                            'Interpreter', 'none');

                        %% ----- Plot bins -----
                        for b = 1:nBins

                            ax = nexttile;
                            hold(ax, 'on');
                            box(ax, 'off');

                            spike_count = 0;

                            bin_start = bin_edges(b);
                            bin_end   = bin_edges(b+1);

                            for aa = 1:numel(Amps_selected)

                                waves = all_spikes_by_bin_amp{b,aa};

                                if isempty(waves)
                                    continue;
                                end

                                % Randomly limit waveforms if too many.
                                if size(waves,1) > max_waveforms_per_bin_amp
                                    rand_idx = randperm(size(waves,1), max_waveforms_per_bin_amp);
                                    waves_plot = waves(rand_idx,:);
                                else
                                    waves_plot = waves;
                                end

                                spike_count = spike_count + size(waves,1);

                                %% ----- Optional trough alignment -----
                                if align_waveforms_to_trough
                                    aligned_waves = zeros(size(waves_plot));

                                    for k = 1:size(waves_plot,1)
                                        [~, min_idx] = min(waves_plot(k,:));
                                        shift = ceil(size(waves_plot,2)/2) - min_idx;
                                        aligned_waves(k,:) = circshift(waves_plot(k,:), shift, 2);
                                    end

                                    waves_to_plot = aligned_waves;
                                else
                                    waves_to_plot = waves_plot;
                                end

                                %% ----- Plot individual waveforms -----
                                if plot_individual_waveforms
                                    plot(ax, t_wave, waves_to_plot', ...
                                         'Color', [cmap(find(Amps == Amps_selected(aa),1,'first'),:) 0.25]);
                                end

                                %% ----- Plot mean waveform -----
                                if plot_mean_waveform
                                    mean_wf = mean(waves_to_plot, 1);
                                    plot(ax, t_wave, mean_wf, ...
                                         'Color', cmap(find(Amps == Amps_selected(aa),1,'first'),:), ...
                                         'LineWidth', 1.8);
                                end
                            end

                            title(ax, sprintf('%.0f–%.0f ms | %d spikes', ...
                                  bin_start, bin_end, spike_count), ...
                                  'Interpreter', 'none');

                            xlabel(ax, 'Waveform time (ms)');
                            ylabel(ax, 'uV');

                            ylim(ax, y_lim);
                            yticks(ax, linspace(y_lim(1), y_lim(2), 3));
                            xticks(ax, round(linspace(t_wave(1), t_wave(end), 3), 2));

                            grid(ax, 'on');
                            axis(ax, 'square');

                            % Mark if this bin overlaps the final artifact time.
                            if isfinite(finalArtifact_ms) && ...
                               finalArtifact_ms >= bin_start && finalArtifact_ms < bin_end

                                text(ax, 0.05, 0.90, 'Final artifact bin', ...
                                     'Units', 'normalized', ...
                                     'FontSize', 8, ...
                                     'FontWeight', 'bold');
                            end
                        end

                        %% ----- Legend -----
                        legend_handles = gobjects(numel(Amps_selected),1);
                        legend_labels  = cell(numel(Amps_selected),1);

                        for aa = 1:numel(Amps_selected)
                            amp_global_idx = find(Amps == Amps_selected(aa), 1, 'first');
                            legend_handles(aa) = plot(nan, nan, '-', ...
                                'Color', cmap(amp_global_idx,:), ...
                                'LineWidth', 1.8);
                            legend_labels{aa} = sprintf('%g uA', Amps_selected(aa));
                        end

                        legend(legend_handles, legend_labels, ...
                               'Location', 'northeastoutside', ...
                               'Box', 'off');

                    end % ISI
                end % prefix

            %% =====================================================
            % Full simultaneous condition
            % =====================================================
            elseif condType == 2

                activeCounts_sim = unique(activeCount_trial(conditionType_trial == 2));
                activeCounts_sim = activeCounts_sim(activeCounts_sim > 0);

                for aci = 1:numel(activeCounts_sim)

                    activeCount = activeCounts_sim(aci);

                    if pool_simultaneous_across_sets

                        group_trials = find(conditionType_trial == 2 & ...
                                            activeCount_trial == activeCount & ...
                                            ismember(trialAmps, Amps_selected));

                    else

                        group_trials = find(conditionSetID_trial == setID & ...
                                            conditionType_trial == 2 & ...
                                            activeCount_trial == activeCount & ...
                                            ismember(trialAmps, Amps_selected));
                    end

                    group_trials = group_trials(group_trials <= nTrials_use);

                    if isempty(group_trials)
                        continue;
                    end

                    %% ----- Collect waveforms by bin and amplitude -----
                    all_spikes_by_bin_amp = cell(nBins, numel(Amps_selected));

                    for tt = 1:numel(group_trials)

                        tr = group_trials(tt);
                        t0_ms = trig(tr) / FS * 1000;

                        rel_times_all = sp_times_all - t0_ms;

                        in_win = rel_times_all >= waveform_time_window_ms(1) & ...
                                 rel_times_all <  waveform_time_window_ms(2);

                        if ~any(in_win)
                            continue;
                        end

                        rel_times = rel_times_all(in_win);
                        waveforms = sp_wave_all(in_win,:);

                        amp_val_trial = trialAmps(tr);
                        amp_pos = find(Amps_selected == amp_val_trial, 1, 'first');

                        if isempty(amp_pos)
                            continue;
                        end

                        for j = 1:numel(rel_times)

                            bin_idx = find(rel_times(j) >= bin_edges(1:end-1) & ...
                                           rel_times(j) <  bin_edges(2:end), ...
                                           1, 'first');

                            if isempty(bin_idx)
                                continue;
                            end

                            all_spikes_by_bin_amp{bin_idx, amp_pos}(end+1,:) = waveforms(j,:);
                        end
                    end

                    all_waves = cell2mat(all_spikes_by_bin_amp(:));

                    if isempty(all_waves)
                        continue;
                    end

                    y_max = max(abs(all_waves(:)));

                    if isempty(y_max) || y_max == 0 || isnan(y_max)
                        y_lim = [-100 100];
                    else
                        y_lim = [-1 1] * ceil(y_max/50)*50;
                    end

                    finalArtifact_trials = lastActivePTD_ms(group_trials);
                    finalArtifact_ms = max(finalArtifact_trials);

                    %% ----- Figure -----
                    if pool_simultaneous_across_sets
                        simLabel = 'Simultaneous pooled across sets';
                    else
                        simLabel = sprintf('Simultaneous Set %d', setID);
                    end

                    figName = sprintf('SSD Waveforms | Ch%d | %s | Active%d', ...
                        ich, simLabel, activeCount);

                    figure('Name', figName, ...
                           'Color', 'w', ...
                           'Position', [100 100 1400 800]);

                    tiledlayout(layout_row, layout_col, ...
                                'Padding', 'compact', ...
                                'TileSpacing', 'compact');

                    sgtitle(sprintf('SSD Waveforms | Rec Ch %d | %s | Active count %d | Final artifact %.1f ms', ...
                        ich, simLabel, activeCount, finalArtifact_ms), ...
                        'Interpreter', 'none');

                    %% ----- Plot bins -----
                    for b = 1:nBins

                        ax = nexttile;
                        hold(ax, 'on');
                        box(ax, 'off');

                        spike_count = 0;

                        bin_start = bin_edges(b);
                        bin_end   = bin_edges(b+1);

                        for aa = 1:numel(Amps_selected)

                            waves = all_spikes_by_bin_amp{b,aa};

                            if isempty(waves)
                                continue;
                            end

                            if size(waves,1) > max_waveforms_per_bin_amp
                                rand_idx = randperm(size(waves,1), max_waveforms_per_bin_amp);
                                waves_plot = waves(rand_idx,:);
                            else
                                waves_plot = waves;
                            end

                            spike_count = spike_count + size(waves,1);

                            if align_waveforms_to_trough
                                aligned_waves = zeros(size(waves_plot));

                                for k = 1:size(waves_plot,1)
                                    [~, min_idx] = min(waves_plot(k,:));
                                    shift = ceil(size(waves_plot,2)/2) - min_idx;
                                    aligned_waves(k,:) = circshift(waves_plot(k,:), shift, 2);
                                end

                                waves_to_plot = aligned_waves;
                            else
                                waves_to_plot = waves_plot;
                            end

                            if plot_individual_waveforms
                                plot(ax, t_wave, waves_to_plot', ...
                                     'Color', [cmap(find(Amps == Amps_selected(aa),1,'first'),:) 0.25]);
                            end

                            if plot_mean_waveform
                                mean_wf = mean(waves_to_plot, 1);
                                plot(ax, t_wave, mean_wf, ...
                                     'Color', cmap(find(Amps == Amps_selected(aa),1,'first'),:), ...
                                     'LineWidth', 1.8);
                            end
                        end

                        title(ax, sprintf('%.0f–%.0f ms | %d spikes', ...
                              bin_start, bin_end, spike_count), ...
                              'Interpreter', 'none');

                        xlabel(ax, 'Waveform time (ms)');
                        ylabel(ax, 'uV');

                        ylim(ax, y_lim);
                        yticks(ax, linspace(y_lim(1), y_lim(2), 3));
                        xticks(ax, round(linspace(t_wave(1), t_wave(end), 3), 2));

                        grid(ax, 'on');
                        axis(ax, 'square');

                        if isfinite(finalArtifact_ms) && ...
                           finalArtifact_ms >= bin_start && finalArtifact_ms < bin_end

                            text(ax, 0.05, 0.90, 'Final artifact bin', ...
                                 'Units', 'normalized', ...
                                 'FontSize', 8, ...
                                 'FontWeight', 'bold');
                        end
                    end

                    %% ----- Legend -----
                    legend_handles = gobjects(numel(Amps_selected),1);
                    legend_labels  = cell(numel(Amps_selected),1);

                    for aa = 1:numel(Amps_selected)
                        amp_global_idx = find(Amps == Amps_selected(aa), 1, 'first');
                        legend_handles(aa) = plot(nan, nan, '-', ...
                            'Color', cmap(amp_global_idx,:), ...
                            'LineWidth', 1.8);
                        legend_labels{aa} = sprintf('%g uA', Amps_selected(aa));
                    end

                    legend(legend_handles, legend_labels, ...
                           'Location', 'northeastoutside', ...
                           'Box', 'off');

                end % active count
            end % condition type
        end % condition type
    end % set
end % channel