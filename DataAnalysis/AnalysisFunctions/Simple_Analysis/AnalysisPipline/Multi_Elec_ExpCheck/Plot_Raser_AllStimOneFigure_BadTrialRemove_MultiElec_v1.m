%% ============================================================
%  Clean Raster + PSTH Plot: Prefix 5 Sequential Sets vs One Sim Reference
%  Figure structure:
%    One figure = one recording channel x one amplitude x one ISI
%    Raster subplots:
%       Seq Set 1, Prefix 5
%       Seq Set 2, Prefix 5
%       Seq Set 3, Prefix 5
%       ...
%       Simultaneous reference
%    Bottom subplot:
%       PSTH overlay for all sequential sets + one simultaneous reference
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ======================

data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

% Recording channels to plot, using Depth_s index.
channels_to_plot = 1:25;

Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Amplitudes to plot. [] means all non-zero amplitudes.
amps_to_plot = [10];

% Full sequential prefix to plot.
TargetSeqPrefix = 5;

% Sequential ISI to plot.
isi_to_plot_ms = [5];

% Stimulation set/order IDs to plot.
% [] means all detected sequential sets.
sets_to_plot = [];

% Full simultaneous active electrode count.
TargetSimActiveCount = 5;

% Simultaneous reference set.
% [] means automatically use the first set that has simultaneous trials.
% You can force it, for example:
% SimReferenceSetID = 1;
SimReferenceSetID = [];

% Bad trial removal.
remove_bad_trials = true;

% Number of trials to plot in each raster.
nTrials_to_plot = 30;

% Plotting time window around trigger.
ras_win = [-10 80];      % ms

% PSTH settings.
bin_ms = 1;
smooth_ms = 10;

% Raster line settings.
raster_line_width = 1.1;

% Print one example trial per set to confirm alignment.
debug_print_trial_content = true;

%% ====================== CHECK FOLDER ======================

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end

cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ====================== BASE NAME ======================

parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');

if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ====================== LOAD SPIKES ======================

[sp, spike_file_used, spike_variable_used] = load_spike_file_for_plot(data_folder);

fprintf('Loaded spike file: %s\n', spike_file_used);
fprintf('Using spike variable: %s\n', spike_variable_used);

nChn_spike = numel(sp);

%% ====================== LOAD BAD TRIALS ======================

BadTrials = {};
BadTrialFile = '';

if remove_bad_trials

    bad_file = sprintf('%s.MixedPrefixBadTrials.mat', base_name);
    bad_path = fullfile(data_folder, bad_file);

    if isfile(bad_path)
        tmp_bad = load(bad_path);

        if isfield(tmp_bad, 'BadTrials')
            BadTrials = tmp_bad.BadTrials;
            BadTrialFile = bad_file;
            fprintf('Loaded bad-trial file: %s\n', bad_file);
        else
            warning('Bad-trial file found but variable BadTrials missing. No bad trials removed.');
        end
    else
        fprintf('No bad-trial file found: %s\n', bad_file);
        fprintf('No bad trials will be removed.\n');
    end
else
    fprintf('remove_bad_trials = false. No bad trials removed.\n');
end

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
simultaneous_stim = S.simultaneous_stim;
n_Trials          = S.n_Trials;

if isfield(S, 'CHN')
    CHN = S.CHN;
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

%% ====================== SAMPLING RATE ======================

try
    [~, freq_params] = read_Intan_RHS2000_file;
    FS = freq_params.amplifier_sample_rate;
catch
    FS = 30000;
    warning('Could not read info.rhs. Using FS = 30000 Hz.');
end

fprintf('Sampling rate: %.1f Hz\n', FS);

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
    error('StimParams does not contain columns 26-30. Cannot use mixed-prefix metadata.');
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

fprintf('\nUsing trial metadata from StimParams columns 26-30.\n');

%% ====================== AMPLITUDE PER TRIAL ======================

trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

%% ====================== FINAL ACTIVE ARTIFACT TIME ======================

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

    % Column 6: PTD in us.
    ptd_us = cell2mat(StimParams_data(activeRows,6));

    % Column 8: pulse train number.
    pulseNum = cell2mat(StimParams_data(activeRows,8));

    % Column 9: pulse train period in us.
    pulsePeriod_us = cell2mat(StimParams_data(activeRows,9));

    pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
    pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

    rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowFinalArtifact_us);
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

%% ====================== APPLY USER FILTERS ======================

isSeqTrial = conditionType_trial == 1;
isSimTrial = conditionType_trial == 2;

% Amplitudes.
all_amps = unique(trialAmps(conditionType_trial > 0));
all_amps = all_amps(all_amps > 0);

if isempty(amps_to_plot)
    amps_sel = all_amps;
else
    amps_sel = intersect(all_amps, amps_to_plot);
end

% ISIs.
all_isis = unique(isi_ms_trial(isSeqTrial));

if isempty(isi_to_plot_ms)
    isi_sel = all_isis;
else
    isi_sel = intersect(all_isis, isi_to_plot_ms);
end

% Sequential set/order IDs only.
all_seq_sets = unique(conditionSetID_trial(isSeqTrial & ...
                                           prefixLength_trial == TargetSeqPrefix & ...
                                           conditionSetID_trial > 0));
all_seq_sets = sort(all_seq_sets(:))';

if isempty(sets_to_plot)
    set_sel = all_seq_sets;
else
    set_sel = intersect(all_seq_sets, sets_to_plot);
end

fprintf('\nSelected amplitudes: ');
disp(amps_sel');

fprintf('Selected sequential prefix: %d\n', TargetSeqPrefix);

fprintf('Selected ISIs (ms): ');
disp(isi_sel');

fprintf('Selected sequential set IDs: ');
disp(set_sel');

%% ====================== FIND SIMULTANEOUS REFERENCE SET ======================

sim_sets_available = unique(conditionSetID_trial(isSimTrial & ...
                                                 activeCount_trial == TargetSimActiveCount & ...
                                                 conditionSetID_trial > 0));
sim_sets_available = sort(sim_sets_available(:))';

if isempty(sim_sets_available)
    warning('No simultaneous trials found for ActiveCount = %d.', TargetSimActiveCount);
    SimReferenceSetID_used = NaN;
else
    if isempty(SimReferenceSetID)
        SimReferenceSetID_used = sim_sets_available(1);
        fprintf('Automatic simultaneous reference set: Set %d\n', SimReferenceSetID_used);
    else
        SimReferenceSetID_used = SimReferenceSetID;

        if ~ismember(SimReferenceSetID_used, sim_sets_available)
            warning('Requested SimReferenceSetID = %d was not found among available sim sets.', ...
                SimReferenceSetID_used);
        else
            fprintf('Using requested simultaneous reference set: Set %d\n', SimReferenceSetID_used);
        end
    end
end

fprintf('Available simultaneous set IDs: ');
disp(sim_sets_available);

fprintf('Detected final artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ====================== DEBUG TRIAL CONTENT CHECK ======================

if debug_print_trial_content && ~isempty(set_sel) && ~isempty(amps_sel) && ~isempty(isi_sel)

    debug_amp = amps_sel(end);
    debug_isi = isi_sel(1);

    fprintf('\n================ DEBUG TRIAL CONTENT CHECK ================\n');
    fprintf('Debug Amp = %.1f uA | ISI = %.1f ms | Prefix = %d\n', ...
        debug_amp, debug_isi, TargetSeqPrefix);

    for i_set = 1:numel(set_sel)

        debug_set = set_sel(i_set);

        tr_debug = find(conditionSetID_trial == debug_set & ...
                        conditionType_trial == 1 & ...
                        prefixLength_trial == TargetSeqPrefix & ...
                        isi_ms_trial == debug_isi & ...
                        trialAmps == debug_amp, ...
                        1, 'first');

        if isempty(tr_debug)
            fprintf('\nSeq Set %d: no matching trial found.\n', debug_set);
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

        fprintf('\nSeq Set %d | Trial %d\n', debug_set, tr_debug);
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

    if isfinite(SimReferenceSetID_used)

        tr_sim_debug = find(conditionSetID_trial == SimReferenceSetID_used & ...
                            conditionType_trial == 2 & ...
                            activeCount_trial == TargetSimActiveCount & ...
                            trialAmps == debug_amp, ...
                            1, 'first');

        if ~isempty(tr_sim_debug)

            rr = (tr_sim_debug-1)*simultaneous_stim + (1:simultaneous_stim);

            fprintf('\nSim Reference Set %d | Trial %d\n', SimReferenceSetID_used, tr_sim_debug);
            fprintf('  conditionType = %d, setID = %d, amp = %.1f uA\n', ...
                conditionType_trial(tr_sim_debug), ...
                conditionSetID_trial(tr_sim_debug), ...
                trialAmps(tr_sim_debug));

            fprintf('  StimNames:       ');
            disp(StimParams_data(rr,1)');

            fprintf('  PTD_us:          ');
            disp(cell2mat(StimParams_data(rr,6))');

            fprintf('  Amp_col16:       ');
            disp(cell2mat(StimParams_data(rr,16))');

            fprintf('  ActiveCount_col26: ');
            disp(cell2mat(StimParams_data(rr,26))');

            fprintf('  CondType_col29:    ');
            disp(cell2mat(StimParams_data(rr,29))');

            fprintf('  SetID_col30:       ');
            disp(cell2mat(StimParams_data(rr,30))');

            fprintf('  FinalArtifact_us = %.1f\n', lastActivePTD_us(tr_sim_debug));
        end
    end

    fprintf('===========================================================\n');
end

%% ====================== PSTH SETTINGS ======================

edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms / 1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

%% ====================== CHANNEL MAP ======================

d = Depth_s(Electrode_Type);

%% ====================== MAIN PLOTTING LOOP ======================

% One figure:
%   recording channel x amplitude x ISI
%
% Raster subplots:
%   one subplot per sequential set/order
%   one extra subplot for simultaneous reference
%
% Bottom subplot:
%   PSTH overlay for all sequential sets + one simultaneous reference

for ich = 1:length(channels_to_plot)

    ch_plot = channels_to_plot(ich);
    ch_spike = d(ch_plot);

    if ch_spike > nChn_spike
        warning('Channel %d maps to spike channel %d, but spike file only has %d channels. Skipped.', ...
                ch_plot, ch_spike, nChn_spike);
        continue;
    end

    if isempty(sp{ch_spike})
        fprintf('Channel %d has no spikes. Skipped.\n', ch_plot);
        continue;
    end

    spike_times_ch = sp{ch_spike}(:,1);

    %% ----- Bad trials for this channel -----
    bad_trs_ch = get_bad_trials_for_channel(BadTrials, ch_plot);

    for i_amp = 1:numel(amps_sel)

        amp_val = amps_sel(i_amp);

        for i_isi = 1:numel(isi_sel)

            isi_val = isi_sel(i_isi);

            %% ----- Build condition list -----
            plotConds = struct([]);
            cond_i = 0;

            % Sequential sets.
            for i_set = 1:numel(set_sel)

                set_id = set_sel(i_set);

                cond_i = cond_i + 1;
                plotConds(cond_i).type = 'seq';
                plotConds(cond_i).set_id = set_id;
                plotConds(cond_i).label = sprintf('Seq Set %d', set_id);
            end

            % One simultaneous reference.
            if isfinite(SimReferenceSetID_used)

                cond_i = cond_i + 1;
                plotConds(cond_i).type = 'sim';
                plotConds(cond_i).set_id = SimReferenceSetID_used;
                plotConds(cond_i).label = sprintf('Sim ref Set %d', SimReferenceSetID_used);
            end

            nRaster = numel(plotConds);

            if nRaster == 0
                continue;
            end

            %% ----- Check whether this figure has any clean data -----
            hasAnyData = false;

            for cc = 1:nRaster

                tlist_raw = get_condition_trials( ...
                    plotConds(cc), ...
                    conditionType_trial, conditionSetID_trial, ...
                    prefixLength_trial, isi_ms_trial, activeCount_trial, ...
                    trialAmps, amp_val, isi_val, ...
                    TargetSeqPrefix, TargetSimActiveCount, ...
                    nTrials_use);

                tlist_clean = setdiff(tlist_raw(:), bad_trs_ch(:), 'stable');

                if ~isempty(tlist_clean)
                    hasAnyData = true;
                    break;
                end
            end

            if ~hasAnyData
                continue;
            end

            %% ----- Colors -----
            nSeq = numel(set_sel);
            seq_colors = lines(max(nSeq,1));
            sim_color = [0.2 0.2 0.2];

            %% ----- Create figure -----
            figName = sprintf('CleanRasterPSTH_Ch%d_Amp%g_ISI%gms_Prefix%d_vs_Sim', ...
                              ch_plot, amp_val, isi_val, TargetSeqPrefix);

            figure('Name', figName, ...
                   'Color', 'w', ...
                   'Position', [100 100 700 700]);

            tiledlayout(nRaster + 1, 1, ...
                        'TileSpacing', 'compact', ...
                        'Padding', 'compact');

            psth_all = cell(nRaster,1);
            psth_labels = cell(nRaster,1);
            psth_colors = cell(nRaster,1);
            psth_styles = cell(nRaster,1);

            maxRate = 0;

            %% ====================== RASTER SUBPLOTS ======================
            for cc = 1:nRaster

                cond = plotConds(cc);

                %% ----- Colour and line style -----
                if strcmp(cond.type, 'seq')

                    seq_idx = find(set_sel == cond.set_id, 1);

                    if isempty(seq_idx)
                        thisColor = [0 0 0];
                    else
                        thisColor = seq_colors(seq_idx,:);
                    end

                    thisStyle = '-';

                elseif strcmp(cond.type, 'sim')

                    thisColor = sim_color;
                    thisStyle = '--';

                else

                    thisColor = [0 0 0];
                    thisStyle = '-';
                end

                %% ----- Raw and cleaned trial list -----
                tlist_raw = get_condition_trials( ...
                    cond, ...
                    conditionType_trial, conditionSetID_trial, ...
                    prefixLength_trial, isi_ms_trial, activeCount_trial, ...
                    trialAmps, amp_val, isi_val, ...
                    TargetSeqPrefix, TargetSimActiveCount, ...
                    nTrials_use);

                tlist_clean_all = setdiff(tlist_raw(:), bad_trs_ch(:), 'stable');

                nRaw = numel(tlist_raw);
                nCleanAll = numel(tlist_clean_all);
                nRemoved = nRaw - nCleanAll;

                % Limit number of raster trials after cleaning.
                tlist_plot = tlist_clean_all;

                if ~isempty(tlist_plot)
                    tlist_plot = tlist_plot(1:min(nTrials_to_plot, numel(tlist_plot)));
                end

                axR = nexttile;
                hold(axR, 'on');
                box(axR, 'off');

                counts = zeros(1, numel(edges)-1);

                if isempty(tlist_plot)

                    title(axR, sprintf('%s | No clean trials | raw %d, clean %d, removed %d', ...
                          cond.label, nRaw, nCleanAll, nRemoved), ...
                          'Interpreter', 'none');

                    psth_all{cc} = zeros(1, numel(ctrs));
                    psth_labels{cc} = cond.label;
                    psth_colors{cc} = thisColor;
                    psth_styles{cc} = thisStyle;
                    ptd_this = NaN;

                else

                    %% ----- Get representative final artifact time -----
                    ptd_vals = unique(lastActivePTD_ms(tlist_plot));
                    ptd_vals = ptd_vals(~isnan(ptd_vals));

                    if isempty(ptd_vals)
                        ptd_this = NaN;
                    else
                        ptd_this = max(ptd_vals);
                    end

                    %% ----- Plot raster trials -----
                    for k = 1:numel(tlist_plot)

                        tr = tlist_plot(k);

                        if tr > length(trig)
                            continue;
                        end

                        t0 = trig(tr) / FS * 1000;

                        rel_t = spike_times_ch - t0;
                        rel_t = rel_t(rel_t >= ras_win(1) & rel_t <= ras_win(2));

                        for s = 1:numel(rel_t)
                            plot(axR, [rel_t(s) rel_t(s)], [k-0.4 k+0.4], ...
                                 'Color', thisColor, ...
                                 'LineWidth', raster_line_width);
                        end

                        counts = counts + histcounts(rel_t, edges);
                    end

                    %% ----- Compute PSTH using plotted clean trials -----
                    rate = counts / (numel(tlist_plot) * bin_s);
                    rate = filter(g, 1, rate);

                    psth_all{cc} = rate;
                    psth_labels{cc} = cond.label;
                    psth_colors{cc} = thisColor;
                    psth_styles{cc} = thisStyle;

                    maxRate = max(maxRate, max(rate));

                    %% ----- Raster title -----
                    % title(axR, sprintf('%s | Final artifact %.1f ms | raw %d, clean %d, plotted %d, removed %d', ...
                    %                    cond.label, ptd_this, nRaw, nCleanAll, numel(tlist_plot), nRemoved), ...
                    %                    'Interpreter', 'none');
                    title(axR, sprintf('%s | Final artifact %.1f ms', cond.label, ptd_this),'Interpreter', 'none');
                end

                %% ----- Raster axis formatting -----
                xline(axR, 0, 'r--', 'LineWidth', 1);

                if isfinite(ptd_this) && ptd_this > 0
                    xline(axR, ptd_this, 'k--', 'LineWidth', 1);
                end

                xlim(axR, ras_win);
                ylim(axR, [0 max(1, numel(tlist_plot)+1)]);
                ylabel(axR, 'Trial');

                if cc < nRaster
                    set(axR, 'XTickLabel', []);
                else
                    xlabel(axR, 'Time from trigger (ms)');
                end
            end

            %% ====================== BOTTOM PSTH SUBPLOT ======================

            axP = nexttile;
            hold(axP, 'on');
            box(axP, 'off');

            for cc = 1:nRaster

                if isempty(psth_all{cc})
                    continue;
                end

                plot(axP, ctrs, psth_all{cc}, ...
                     'Color', psth_colors{cc}, ...
                     'LineStyle', psth_styles{cc}, ...
                     'LineWidth', 2);
            end

            xline(axP, 0, 'r--', 'LineWidth', 1);
            xlim(axP, ras_win);

            if maxRate <= 0
                ylim(axP, [0 50]);
            else
                ylim(axP, [0 max(50, ceil(maxRate*1.1/10)*10)]);
            end

            xlabel(axP, 'Time from trigger (ms)');
            ylabel(axP, 'Rate (sp/s)');
            title(axP, 'PSTH overlay', 'Interpreter', 'none');

            legend(axP, psth_labels, 'Box', 'off', 'Location', 'northeast');

            %% ====================== FIGURE TITLE ======================

            sgtitle(sprintf('Rec Ch %d | Prefix %d Seq Sets vs Sim Ref | %g uA | Seq ISI %.1f ms', ...
                ch_plot, TargetSeqPrefix, amp_val, isi_val), ...
                'Interpreter', 'none');

        end
    end
end

%% ========================================================================
%  HELPER FUNCTIONS
%% ========================================================================

function [sp, spike_file_used, spike_variable_used] = load_spike_file_for_plot(data_folder)

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
        error('No spike file found in %s.', data_folder);
    end
end

function bad_trs = get_bad_trials_for_channel(BadTrials, ch_plot)

    bad_trs = [];

    if isempty(BadTrials)
        return;
    end

    if iscell(BadTrials)

        if ch_plot <= numel(BadTrials) && ~isempty(BadTrials{ch_plot})
            bad_trs = BadTrials{ch_plot}(:);
        else
            bad_trs = [];
        end

    elseif isnumeric(BadTrials)

        bad_trs = BadTrials(:);

    else

        bad_trs = [];
    end
end

function tlist = get_condition_trials( ...
    cond, ...
    conditionType_trial, conditionSetID_trial, ...
    prefixLength_trial, isi_ms_trial, activeCount_trial, ...
    trialAmps, amp_val, isi_val, ...
    TargetSeqPrefix, TargetSimActiveCount, ...
    nTrials_use)

    if strcmp(cond.type, 'seq')

        tlist = find(conditionType_trial == 1 & ...
                     conditionSetID_trial == cond.set_id & ...
                     prefixLength_trial == TargetSeqPrefix & ...
                     abs(isi_ms_trial - isi_val) < 1e-6 & ...
                     abs(trialAmps - amp_val) < 1e-6);

    elseif strcmp(cond.type, 'sim')

        tlist = find(conditionType_trial == 2 & ...
                     conditionSetID_trial == cond.set_id & ...
                     activeCount_trial == TargetSimActiveCount & ...
                     abs(trialAmps - amp_val) < 1e-6);

    else

        tlist = [];
    end

    tlist = tlist(tlist <= nTrials_use);
end