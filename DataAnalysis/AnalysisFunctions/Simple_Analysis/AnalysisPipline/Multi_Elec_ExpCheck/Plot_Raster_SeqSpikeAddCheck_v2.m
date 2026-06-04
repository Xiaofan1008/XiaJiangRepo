%% ============================================================
%  Prefix Recovery Check Plot
%  Purpose:
%    Plot recovered spike rasters and PSTHs to check whether prefix-based
%    spike recovery looks correct.
%  Figure structure:
%    One figure = one recording channel × one stimulation set × one amplitude × one ISI
%  Input spike file:
%    *.sp_xia_PrefixRecovery.mat
%  Required spike variable:
%    sp_seq
%  Important modification:
%    Trial metadata are read directly from randomized StimParams columns 26–30.
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';

% Recording channels to plot, using Depth_s index.
channels_to_plot = 17:25;

% Amplitudes to plot. [] means all non-zero amplitudes.
amps_to_plot = [10];

% Prefixes to plot. [] means all detected sequential prefixes.
prefix_to_plot = [1 2 3 4 5];

% ISIs to plot, in ms. [] means all detected ISIs.
isi_to_plot_ms = [];

% Stimulation set/order IDs to plot. [] means all detected sets.
sets_to_plot = [];

% Number of trials to plot in each raster.
nTrials_to_plot = 30;

% Plotting time window around trigger.
ras_win = [-5 40];      % ms

% PSTH settings.
bin_ms = 1;             % bin size for PSTH
smooth_ms = 2;          % smoothing width

% Electrode type for Depth_s mapping.
Electrode_Type = 2;     % 0 rigid, 1 single-shank flex, 2 four-shank flex

% Raster line settings.
raster_line_width = 1.1;

% Print one example trial per prefix to confirm alignment.
debug_print_trial_content = true;

%% ====================== CHECK FOLDER ======================
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end

cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ====================== LOAD RECOVERED SPIKES ======================
rec_file = dir('*sp_xia_PrefixRecovery.mat');
assert(~isempty(rec_file), 'No *sp_xia_PrefixRecovery.mat file found. Run the recovery code first.');

load(rec_file(1).name, 'sp_seq');
fprintf('Loaded recovered spike file: %s\n', rec_file(1).name);

nChn_spike = numel(sp_seq);

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
simultaneous_stim = S.simultaneous_stim;   % rows/slots per trial
n_Trials          = S.n_Trials;
E_MAP             = S.E_MAP;

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
% Important:
%   Do NOT use separate metadata arrays directly here.
%   They may not be aligned with randomized StimParams.
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

%% ====================== AMPLITUDE PER TRIAL ======================
trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

%% ====================== FINAL ACTIVE ARTIFACT TIME ======================
% lastActivePTD_ms tells where the final artifact occurs for each trial.
%
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

%% ====================== DEBUG TRIAL CONTENT CHECK ======================
if debug_print_trial_content

    isPrefixTrial_tmp = conditionType_trial == 1;

    prefix_all_debug = unique(prefixLength_trial(isPrefixTrial_tmp & prefixLength_trial > 0));
    set_all_debug = unique(conditionSetID_trial(conditionSetID_trial > 0));
    amp_all_debug = unique(trialAmps(isPrefixTrial_tmp));
    amp_all_debug = amp_all_debug(amp_all_debug > 0);
    isi_all_debug = unique(isi_ms_trial(isPrefixTrial_tmp));

    if ~isempty(prefix_all_debug) && ~isempty(set_all_debug) && ...
       ~isempty(amp_all_debug) && ~isempty(isi_all_debug)

        debug_set = set_all_debug(1);
        debug_amp = amp_all_debug(end);
        debug_isi = isi_all_debug(1);

        fprintf('\n================ DEBUG TRIAL CONTENT CHECK ================\n');
        fprintf('Debug Set = %d | Amp = %.1f uA | ISI = %.1f ms\n', ...
            debug_set, debug_amp, debug_isi);

        for ip = 1:numel(prefix_all_debug)

            prefix_val = prefix_all_debug(ip);

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
end

%% ====================== APPLY USER FILTERS ======================
% Only sequential prefix trials are plotted.
isPrefixTrial = conditionType_trial == 1;

% Amplitudes.
all_amps = unique(trialAmps(isPrefixTrial));
all_amps = all_amps(all_amps > 0);   % exclude zero-control

if isempty(amps_to_plot)
    amps_sel = all_amps;
else
    amps_sel = intersect(all_amps, amps_to_plot);
end

% Prefixes.
all_prefixes = unique(prefixLength_trial(isPrefixTrial & prefixLength_trial > 0));
all_prefixes = sort(all_prefixes(:))';

if isempty(prefix_to_plot)
    prefix_sel = all_prefixes;
else
    prefix_sel = intersect(all_prefixes, prefix_to_plot);
end

% ISIs.
all_isis = unique(isi_ms_trial(isPrefixTrial));

if isempty(isi_to_plot_ms)
    isi_sel = all_isis;
else
    isi_sel = intersect(all_isis, isi_to_plot_ms);
end

% Set/order IDs.
all_sets = unique(conditionSetID_trial(isPrefixTrial & conditionSetID_trial > 0));

if isempty(sets_to_plot)
    set_sel = all_sets;
else
    set_sel = intersect(all_sets, sets_to_plot);
end

fprintf('\nSelected amplitudes: ');
disp(amps_sel');

fprintf('Selected prefixes: ');
disp(prefix_sel);

fprintf('Selected ISIs (ms): ');
disp(isi_sel');

fprintf('Selected set IDs: ');
disp(set_sel');

fprintf('Detected final artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ====================== PSTH SETTINGS ======================
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms / 1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

% Color map for prefixes.
prefix_colors = lines(numel(prefix_sel));

%% ====================== CHANNEL MAP ======================
d = Depth_s(Electrode_Type);

%% ====================== MAIN PLOTTING LOOP ======================
% One figure:
%   recording channel × set ID × amplitude × ISI
%
% Raster subplots:
%   one subplot per prefix
%
% Bottom subplot:
%   overlay PSTH curves for all prefixes

for ich = 1:length(channels_to_plot)

    ch_plot = channels_to_plot(ich);
    ch_spike = d(ch_plot);

    if ch_spike > nChn_spike
        warning('Channel %d maps to spike channel %d, but sp_seq only has %d channels. Skipped.', ...
                ch_plot, ch_spike, nChn_spike);
        continue;
    end

    if isempty(sp_seq{ch_spike})
        fprintf('Channel %d has no spikes. Skipped.\n', ch_plot);
        continue;
    end

    spike_times_ch = sp_seq{ch_spike}(:,1);

    for i_set = 1:numel(set_sel)

        set_id = set_sel(i_set);

        %% ----- Build set label -----
        if ~isempty(CHN) && set_id <= size(CHN,1)
            stimVec = CHN(set_id,:);
            stimVec = stimVec(stimVec > 0);
            set_label = strjoin(arrayfun(@(x) sprintf('Ch%d',x), ...
                                  stimVec,'UniformOutput',false),'→');
        else
            set_label = sprintf('Set%d', set_id);
        end

        for i_amp = 1:numel(amps_sel)

            amp_val = amps_sel(i_amp);

            for i_isi = 1:numel(isi_sel)

                isi_val = isi_sel(i_isi);

                %% ----- Check whether this figure has any data -----
                hasAnyData = false;

                for ip = 1:numel(prefix_sel)

                    prefix_val = prefix_sel(ip);

                    tlist_check = find(conditionType_trial == 1 & ...
                                       conditionSetID_trial == set_id & ...
                                       prefixLength_trial == prefix_val & ...
                                       isi_ms_trial == isi_val & ...
                                       trialAmps == amp_val);

                    tlist_check = tlist_check(tlist_check <= nTrials_use);

                    if ~isempty(tlist_check)
                        hasAnyData = true;
                        break;
                    end
                end

                if ~hasAnyData
                    continue;
                end

                %% ----- Create figure -----
                nRaster = numel(prefix_sel);

                figName = sprintf('PrefixRecoveryCheck_Ch%d_Set%d_Amp%g_ISI%gms', ...
                                  ch_plot, set_id, amp_val, isi_val);

                figure('Name', figName, ...
                       'Color', 'w', ...
                       'Position', [100 100 900 900]);

                tiledlayout(nRaster + 1, 1, ...
                            'TileSpacing', 'compact', ...
                            'Padding', 'compact');

                % Store PSTH curves for bottom plot.
                psth_all = cell(numel(prefix_sel),1);
                psth_labels = cell(numel(prefix_sel),1);
                maxRate = 0;

                %% ====================== RASTER SUBPLOTS ======================
                for ip = 1:numel(prefix_sel)

                    prefix_val = prefix_sel(ip);
                    thisColor = prefix_colors(ip,:);

                    % Find trials for this prefix condition.
                    tlist = find(conditionType_trial == 1 & ...
                                 conditionSetID_trial == set_id & ...
                                 prefixLength_trial == prefix_val & ...
                                 isi_ms_trial == isi_val & ...
                                 trialAmps == amp_val);

                    tlist = tlist(tlist <= nTrials_use);

                    % Select only first N trials for raster.
                    if ~isempty(tlist)
                        tlist = tlist(1:min(nTrials_to_plot, numel(tlist)));
                    end

                    axR = nexttile;
                    hold(axR, 'on');
                    box(axR, 'off');

                    counts = zeros(1, numel(edges)-1);

                    if isempty(tlist)

                        title(axR, sprintf('Prefix %d | No trials', prefix_val), ...
                              'Interpreter', 'none');

                        psth_all{ip} = zeros(1, numel(ctrs));
                        psth_labels{ip} = sprintf('Prefix %d', prefix_val);
                        ptd_this = NaN;

                    else

                        %% ----- Get representative final artifact time -----
                        ptd_vals = unique(lastActivePTD_ms(tlist));
                        ptd_vals = ptd_vals(~isnan(ptd_vals));

                        if isempty(ptd_vals)
                            ptd_this = NaN;
                        else
                            % Use max to avoid misleading title if anything odd exists.
                            ptd_this = max(ptd_vals);
                        end

                        %% ----- Plot raster trials -----
                        for k = 1:numel(tlist)

                            tr = tlist(k);

                            if tr > length(trig)
                                continue;
                            end

                            t0 = trig(tr) / FS * 1000;

                            rel_t = spike_times_ch - t0;
                            rel_t = rel_t(rel_t >= ras_win(1) & rel_t <= ras_win(2));

                            % Raster: one horizontal row per trial.
                            for s = 1:numel(rel_t)
                                plot(axR, [rel_t(s) rel_t(s)], [k-0.4 k+0.4], ...
                                     'Color', thisColor, ...
                                     'LineWidth', raster_line_width);
                            end

                            % PSTH counts.
                            counts = counts + histcounts(rel_t, edges);
                        end

                        %% ----- Compute PSTH for this prefix -----
                        rate = counts / (numel(tlist) * bin_s);
                        rate = filter(g, 1, rate);

                        psth_all{ip} = rate;
                        psth_labels{ip} = sprintf('Prefix %d', prefix_val);

                        maxRate = max(maxRate, max(rate));

                        %% ----- Raster title -----
                        title(axR, sprintf('Prefix %d | Final artifact %.1f ms | %d trials', ...
                                           prefix_val, ptd_this, numel(tlist)), ...
                                           'Interpreter', 'none');
                    end

                    %% ----- Raster axis formatting -----
                    xline(axR, 0, 'r--', 'LineWidth', 1);

                    % Final active artifact line.
                    if isfinite(ptd_this) && ptd_this > 0
                        xline(axR, ptd_this, 'k--', 'LineWidth', 1);
                    end

                    xlim(axR, ras_win);
                    ylim(axR, [0 max(1, numel(tlist)+1)]);
                    ylabel(axR, 'Trial');

                    if ip < numel(prefix_sel)
                        set(axR, 'XTickLabel', []);
                    else
                        xlabel(axR, 'Time from trigger (ms)');
                    end

                end % prefix raster loop

                %% ====================== BOTTOM PSTH SUBPLOT ======================
                axP = nexttile;
                hold(axP, 'on');
                box(axP, 'off');

                for ip = 1:numel(prefix_sel)

                    if isempty(psth_all{ip})
                        continue;
                    end

                    plot(axP, ctrs, psth_all{ip}, ...
                         'Color', prefix_colors(ip,:), ...
                         'LineWidth', 1.8);
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
                title(axP, 'PSTH overlay across prefixes', 'Interpreter', 'none');
                legend(axP, psth_labels, 'Box', 'off', 'Location', 'northeast');

                %% ====================== FIGURE TITLE ======================
                sgtitle(sprintf('Prefix Recovery Check | Rec Ch %d | Set %d: %s | %g uA | ISI %.1f ms', ...
                    ch_plot, set_id, set_label, amp_val, isi_val), ...
                    'Interpreter', 'none');

            end % ISI
        end % amplitude
    end % set
end % channel