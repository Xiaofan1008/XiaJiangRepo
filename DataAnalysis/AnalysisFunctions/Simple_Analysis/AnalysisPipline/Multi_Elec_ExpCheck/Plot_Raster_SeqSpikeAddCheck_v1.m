%% ============================================================
%  Prefix Recovery Check Plot
%
%  Purpose:
%    Plot recovered spike rasters and PSTHs to check whether prefix-based
%    spike recovery looks correct.
%
%  Figure structure:
%    One figure = one recording channel × one stimulation set × one amplitude × one ISI
%
%  Input spike file:
%    *.sp_xia_PrefixRecovery.mat
%
%  Required spike variable:
%    sp_seq
%
%  Important:
%    This code does NOT compare original vs recovered spikes.
%    It only plots the recovered spikes, to keep the figure clean.
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

% Recording channels to plot.
channels_to_plot = 15:20;

% Amplitudes to plot. [] means all non-zero amplitudes.
amps_to_plot = [10];

% Prefixes to plot.[] means all detected sequential prefixes.
prefix_to_plot = [1 2 3 4 5];

% ISIs to plot, in ms.[] means all detected ISIs.
isi_to_plot_ms = [5];

% Stimulation set/order IDs to plot.[] means all detected sets.
sets_to_plot = [];

% Number of trials to plot in each raster.
% If the condition has more trials than this, only the first N trials are drawn.
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

%% ====================== CHECK FOLDER ======================
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ====================== LOAD RECOVERED SPIKES ======================
% This file should be created by the prefix recovery code.
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

%% ====================== LOAD EXPERIMENT PARAMETERS ======================
fDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fDIR), 'No *_exp_datafile_*.mat found.');

% Load full file so metadata are available.
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

%% ====================== SAMPLING RATE ======================
% Spike times are in ms, but trigger samples need FS conversion.
% If info.rhs is available, read FS from Intan header.
% Otherwise use 30000 Hz.
try
    [~, freq_params] = read_Intan_RHS2000_file;
    FS = freq_params.amplifier_sample_rate;
catch
    FS = 30000;
    warning('Could not read info.rhs. Using FS = 30000 Hz.');
end

%% ====================== LOAD MIXED-PREFIX METADATA ======================
% These variables are required for the new prefix experiment.
requiredVars = {'active_electrode_count_by_trial', ...
                'prefix_length_by_trial', ...
                'isi_ms_by_trial', ...
                'condition_type_by_trial', ...
                'condition_set_id_by_trial'};

for i = 1:numel(requiredVars)
    assert(isfield(S, requiredVars{i}), ...
        'Missing metadata variable "%s". This plotting code requires the new mixed-prefix file.', ...
        requiredVars{i});
end

activeCount_trial    = S.active_electrode_count_by_trial(:);
prefixLength_trial   = S.prefix_length_by_trial(:);
isi_ms_trial         = S.isi_ms_by_trial(:);
conditionType_trial  = S.condition_type_by_trial(:);
conditionSetID_trial = S.condition_set_id_by_trial(:);

%% ====================== AMPLITUDE PER TRIAL ======================
% Use the first row/slot of each trial to get trial amplitude.
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;

%% ====================== LAST ACTIVE PTD ======================
% lastActivePTD_ms tells where the final artifact occurs for each trial.
%
% Example for ISI = 5 ms:
%   Prefix 1 -> 0 ms
%   Prefix 2 -> 5 ms
%   Prefix 3 -> 10 ms
%   Prefix 4 -> 15 ms
%   Prefix 5 -> 20 ms

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials
    activeCount_this = activeCount_trial(tr);

    if isnan(activeCount_this) || activeCount_this < 1
        lastActivePTD_us(tr) = 0;
    else
        activeCount_this = min(round(activeCount_this), simultaneous_stim);

        % StimParams has one header row.
        stimRow = 1 + (tr-1)*simultaneous_stim + activeCount_this;
        lastActivePTD_us(tr) = StimParams{stimRow,6};
    end
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

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

%% ====================== PSTH SETTINGS ======================
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms / 1000;

% Simple smoothing kernel.
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

    for i_set = 1:numel(set_sel)
        set_id = set_sel(i_set);

        %% ----- Build set label -----
        % CHN stores the entered electrode order.
        if ~isempty(CHN) && set_id <= size(CHN,1)
            stimVec = CHN(set_id,:);
            stimVec = stimVec(stimVec > 0);
            set_label = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimVec,'UniformOutput',false),'→');
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

                % Use one extra row at the bottom for PSTH.
                figName = sprintf('PrefixRecoveryCheck_Ch%d_Set%d_Amp%g_ISI%gms', ...
                                  ch_plot, set_id, amp_val, isi_val);

                figure('Name', figName, ...
                       'Color', 'w', ...
                       'Position', [100 100 1500 900]);

                tl = tiledlayout(nRaster + 1, 1, ...
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
                    else
                        %% ----- Get representative PTD for this prefix -----
                        ptd_vals = unique(lastActivePTD_ms(tlist));
                        ptd_vals = ptd_vals(~isnan(ptd_vals));

                        if isempty(ptd_vals)
                            ptd_this = NaN;
                        else
                            ptd_this = ptd_vals(1);
                        end

                        %% ----- Plot raster trials -----
                        for k = 1:numel(tlist)
                            tr = tlist(k);

                            if tr > length(trig)
                                continue;
                            end

                            t0 = trig(tr) / FS * 1000;

                            spike_times = sp_seq{ch_spike}(:,1);
                            rel_t = spike_times - t0;

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

                        %% ----- Raster axis formatting -----
                        title(axR, sprintf('Prefix %d | Last PTD %.1f ms | %d trials', ...
                                           prefix_val, ptd_this, numel(tlist)), ...
                                           'Interpreter', 'none');
                    end

                    % Trigger line.
                    xline(axR, 0, 'r--', 'LineWidth', 1);

                    % Last active PTD line.
                    % Only draw if it is later than 0 ms.
                    if exist('ptd_this', 'var') && isfinite(ptd_this) && ptd_this > 0
                        xline(axR, ptd_this, 'k--', 'LineWidth', 1);
                    end

                    xlim(axR, ras_win);
                    ylim(axR, [0 max(1, numel(tlist)+1)]);
                    ylabel(axR, 'Trial');

                    % Only the last raster subplot needs x-axis label.
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
                sgtitle(sprintf('Prefix Recovery Check | Rec Ch %d | Set %d: %s | %g µA | ISI %.1f ms', ...
                    ch_plot, set_id, set_label, amp_val, isi_val), ...
                    'Interpreter', 'none');

            end % ISI
        end % amplitude
    end % set
end % channel