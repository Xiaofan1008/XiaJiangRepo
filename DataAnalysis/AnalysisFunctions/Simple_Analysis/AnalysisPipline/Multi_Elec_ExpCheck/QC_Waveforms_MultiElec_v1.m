%% ============================================================
%  Spike Waveform Plot for Prefix-Recovery Stimulation Files
%
%  Purpose:
%    Check spike waveform shapes after prefix recovery / filtering.
%
%  Figure structure:
%    One figure = recording channel × stimulation set × prefix × ISI
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

spike_chn_start = 1;
spike_chn_end   = 64;   % Depth_s index
Electrode_Type  = 2;    % 0: rigid; 1: single-shank flex; 2: four-shank flex

% ================= XIA MODIFICATION =================
% New metadata-based selection.
%
% prefix_to_plot:
%   Which sequential prefixes to plot.
%   [] means all detected prefixes.
%
% isi_to_plot_ms:
%   Which ISI values to plot.
%   [] means all detected ISIs.
%
% sets_to_plot:
%   Which condition set/order IDs to plot.
%   [] means all detected sets.
%
% amps_to_plot:
%   Which amplitudes to plot.
%   [] means all non-zero amplitudes.
prefix_to_plot = [1 2 3 4 5];
isi_to_plot_ms = [5];
sets_to_plot   = [];
amps_to_plot   = [];
% =====================================================

FS = 30000;              % Sampling frequency

win_ms = 40;             % waveform latency window after trigger to inspect
bin_ms = 2;              % latency bin size
nBins  = win_ms / bin_ms;

amp_threshold = 300;     % waveform amplitude plotting limit; set Inf to disable

layout_row = 4;          % subplot layout
layout_col = 5;          % 4 × 5 = 20 bins for 0–40 ms at 2 ms/bin

align_waveforms = 1;     % 1 = align waveform troughs; 0 = plot original waveform timing

%% ================= CHECK FOLDER =================
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ================= BASE NAME =================
parts       = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');

if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ================= LOAD SPIKES =================
% ================= XIA MODIFICATION =================
% Load spike files in the new preferred order:
%   1) PrefixRecovery_SSD file: final filtered recovered spikes, sp_corr
%   2) PrefixRecovery file: recovered spikes before SSD filtering, sp_seq
%   3) Base sp_xia file: original clipped spikes, sp_clipped
sp_use = [];
sp_source_label = '';

fprintf('\nTrying to load spike file...\n');

prefix_ssd_file = [base_name '.sp_xia_PrefixRecovery_SSD.mat'];
prefix_rec_file = [base_name '.sp_xia_PrefixRecovery.mat'];
base_file       = [base_name '.sp_xia.mat'];

if isfile(prefix_ssd_file)
    fprintf('Found prefix-recovery SSD file: %s\n', prefix_ssd_file);
    S_sp = load(prefix_ssd_file);

    if isfield(S_sp, 'sp_corr')
        sp_use = S_sp.sp_corr;
        sp_source_label = 'PrefixRecovery SSD / sp_corr';
    elseif isfield(S_sp, 'sp_seq_filtered')
        sp_use = S_sp.sp_seq_filtered;
        sp_source_label = 'PrefixRecovery SSD / sp_seq_filtered';
    else
        error('No usable spike variable found in %s.', prefix_ssd_file);
    end

elseif isfile(prefix_rec_file)
    fprintf('Found prefix-recovery file: %s\n', prefix_rec_file);
    S_sp = load(prefix_rec_file);

    if isfield(S_sp, 'sp_seq')
        sp_use = S_sp.sp_seq;
        sp_source_label = 'PrefixRecovery / sp_seq';
    else
        error('No sp_seq variable found in %s.', prefix_rec_file);
    end

elseif isfile(base_file)
    fprintf('Falling back to base spike file: %s\n', base_file);
    S_sp = load(base_file);

    if isfield(S_sp, 'sp_clipped')
        sp_use = S_sp.sp_clipped;
        sp_source_label = 'Base / sp_clipped';
    elseif isfield(S_sp, 'sp')
        sp_use = S_sp.sp;
        sp_source_label = 'Base / sp';
    else
        error('No usable spike variable found in %s.', base_file);
    end

else
    error('No usable spike file found.');
end

nCh = numel(sp_use);
fprintf('Using spike source: %s\n', sp_source_label);
fprintf('Number of spike channels: %d\n', nCh);
% =====================================================

%% ================= LOAD TRIGGERS =================
if isempty(dir('*.trig.dat'))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ================= LOAD StimParams & METADATA =================
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

% Load full file so all metadata are available.
S = load(fullfile(data_folder, fileDIR(1).name));

StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;   % rows/slots per trial
n_Trials           = S.n_Trials;
E_MAP              = S.E_MAP;

if isfield(S, 'CHN')
    CHN = S.CHN;
else
    CHN = [];
end

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ================= LOAD NEW PREFIX METADATA =================
% These variables are required for the new mixed-prefix stimulation file.
requiredVars = {'active_electrode_count_by_trial', ...
                'prefix_length_by_trial', ...
                'isi_ms_by_trial', ...
                'condition_type_by_trial', ...
                'condition_set_id_by_trial'};

for i = 1:numel(requiredVars)
    assert(isfield(S, requiredVars{i}), ...
        'Missing metadata variable "%s". This code requires the new mixed-prefix file.', ...
        requiredVars{i});
end

activeCount_trial    = S.active_electrode_count_by_trial(:);
prefixLength_trial   = S.prefix_length_by_trial(:);
isi_ms_trial         = S.isi_ms_by_trial(:);
conditionType_trial  = S.condition_type_by_trial(:);
conditionSetID_trial = S.condition_set_id_by_trial(:);

%% ================= AMPLITUDES =================
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:simultaneous_stim:end);

% Important: convert -1 to 0 before unique(), so ampIdx stays correct.
trialAmps(trialAmps == -1) = 0;

[Amps, ~, ampIdx] = unique(trialAmps(:));
n_AMP = numel(Amps);
cmap  = lines(n_AMP);

%% ================= LAST ACTIVE PTD =================
% Used only for title/interpretation.
% For ISI = 5 ms:
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

        % StimParams has a header row.
        stimRow = 1 + (tr-1)*simultaneous_stim + activeCount_this;
        lastActivePTD_us(tr) = StimParams{stimRow,6};
    end
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

%% ================= APPLY USER FILTERS =================
% Only sequential prefix trials are included in this waveform QC plot.
isPrefixTrial = conditionType_trial == 1;

% Prefixes
all_prefixes = unique(prefixLength_trial(isPrefixTrial & prefixLength_trial > 0));
all_prefixes = sort(all_prefixes(:))';

if isempty(prefix_to_plot)
    prefix_sel = all_prefixes;
else
    prefix_sel = intersect(all_prefixes, prefix_to_plot);
end

% ISIs
all_isis = unique(isi_ms_trial(isPrefixTrial));

if isempty(isi_to_plot_ms)
    isi_sel = all_isis;
else
    isi_sel = intersect(all_isis, isi_to_plot_ms);
end

% Sets
all_sets = unique(conditionSetID_trial(isPrefixTrial & conditionSetID_trial > 0));

if isempty(sets_to_plot)
    set_sel = all_sets;
else
    set_sel = intersect(all_sets, sets_to_plot);
end

% Amplitudes
all_amps = unique(trialAmps(isPrefixTrial));
all_amps = all_amps(all_amps > 0);   % exclude zero-control

if isempty(amps_to_plot)
    amps_sel = all_amps;
else
    amps_sel = intersect(all_amps, amps_to_plot);
end

fprintf('\nSelected prefixes: ');
disp(prefix_sel);

fprintf('Selected ISIs (ms): ');
disp(isi_sel');

fprintf('Selected set IDs: ');
disp(set_sel');

fprintf('Selected amplitudes: ');
disp(amps_sel');

fprintf('Detected last active PTDs (ms): ');
disp(unique(lastActivePTD_ms)');

%% ================= ELECTRODE MAP =================
d = Depth_s(Electrode_Type);

%% ================= WAVEFORM TIME AXIS =================
example_ch = find(~cellfun(@isempty, sp_use), 1, 'first');

if isempty(example_ch)
    error('All channels are empty in sp_use.');
end

wf_len = size(sp_use{example_ch}, 2) - 1;
t_wave = (0:wf_len-1) / FS * 1000;

%% ================= SPIKE WAVEFORM PLOTTING =================
% Main loop:
%   channel -> set -> prefix -> ISI
%
% In each figure:
%   subplots are latency bins.
%   waveform colors show amplitude.

for ich = spike_chn_start:spike_chn_end

    ch = d(ich);

    if ch > nCh || isempty(sp_use{ch})
        continue;
    end

    %% ----- Get spike times and waveforms for this channel -----
    sp_times = sp_use{ch}(:,1);
    sp_wave  = sp_use{ch}(:,2:end);

    % Optional waveform amplitude threshold for plotting.
    if isfinite(amp_threshold)
        valid_idx = all(abs(sp_wave) <= amp_threshold, 2);
        sp_times  = sp_times(valid_idx);
        sp_wave   = sp_wave(valid_idx,:);
    end

    if isempty(sp_times)
        continue;
    end

    for i_set = 1:numel(set_sel)

        set_id = set_sel(i_set);

        %% ----- Build stimulation set label -----
        % CHN stores the entered sequence/order.
        if ~isempty(CHN) && set_id <= size(CHN,1)
            stimIdx = CHN(set_id,:);
            stimIdx = stimIdx(stimIdx > 0);
            stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), ...
                                  stimIdx, 'UniformOutput', false), ', ');
        else
            stimLabel = sprintf('Set%d', set_id);
        end

        for ip = 1:numel(prefix_sel)

            prefix_val = prefix_sel(ip);

            for ii = 1:numel(isi_sel)

                isi_val = isi_sel(ii);

                %% ----- Select trials for this figure -----
                trial_mask = conditionType_trial == 1 & ...
                             conditionSetID_trial == set_id & ...
                             prefixLength_trial == prefix_val & ...
                             isi_ms_trial == isi_val;

                % Limit to valid trial range if trigger number differs.
                trial_mask((nTrials_use+1):end) = false;

                if ~any(trial_mask)
                    continue;
                end

                trial_ids = find(trial_mask);

                %% ----- Representative PTD for this prefix -----
                ptd_vals = unique(lastActivePTD_ms(trial_ids));
                ptd_vals = ptd_vals(~isnan(ptd_vals));

                if isempty(ptd_vals)
                    ptd_text = 'NA';
                else
                    ptd_text = sprintf('%.1f ms', ptd_vals(1));
                end

                %% ----- Collect waveforms by latency bin and amplitude -----
                all_spikes_by_bin_amp = cell(nBins, n_AMP);

                for idx = 1:numel(trial_ids)

                    tr = trial_ids(idx);

                    if tr > nTrials_use
                        continue;
                    end

                    t0_ms = trig(tr) / FS * 1000;
                    amp_id = ampIdx(tr);

                    % Only plot amplitudes selected by the user.
                    amp_val_this = Amps(amp_id);
                    if ~ismember(amp_val_this, amps_sel)
                        continue;
                    end

                    % Spikes after trigger within window.
                    mask_sp = sp_times >= t0_ms & sp_times < (t0_ms + win_ms);

                    if ~any(mask_sp)
                        continue;
                    end

                    rel_times = sp_times(mask_sp) - t0_ms;
                    waveforms = sp_wave(mask_sp,:);

                    for j = 1:numel(rel_times)
                        bin_idx = floor(rel_times(j) / bin_ms) + 1;

                        if bin_idx >= 1 && bin_idx <= nBins
                            all_spikes_by_bin_amp{bin_idx, amp_id}(end+1,:) = waveforms(j,:);
                        end
                    end
                end

                %% ----- Skip empty figures -----
                all_waves = cell2mat(all_spikes_by_bin_amp(:));

                if isempty(all_waves)
                    continue;
                end

                %% ----- Determine y limits -----
                y_max = max(abs(all_waves(:)));

                if y_max == 0
                    y_lim = [-50 50];
                else
                    y_lim = [-1 1] * ceil(y_max/50)*50;
                end

                %% ----- Create figure -----
                figName = sprintf('Waveforms_Ch%d_Set%d_Prefix%d_ISI%gms', ...
                                  ich, set_id, prefix_val, isi_val);

                figure('Name', figName, ...
                       'Color','w', ...
                       'Position', [100 100 1400 800]);

                tiledlayout(layout_row, layout_col, ...
                            'Padding','compact', ...
                            'TileSpacing','compact');

                %% ----- Plot waveform bins -----
                for b = 1:nBins

                    nexttile;
                    hold on;

                    spike_count = 0;

                    for a = 1:n_AMP

                        amp_val = Amps(a);

                        % Only plot selected amplitudes.
                        if ~ismember(amp_val, amps_sel)
                            continue;
                        end

                        waves = all_spikes_by_bin_amp{b,a};

                        if isempty(waves)
                            continue;
                        end

                        spike_count = spike_count + size(waves, 1);

                        %% ----- Optional waveform alignment -----
                        if align_waveforms == 1
                            aligned_waves = zeros(size(waves));

                            for k = 1:size(waves,1)
                                [~, min_idx] = min(waves(k,:));
                                shift = ceil(size(waves,2)/2) - min_idx;
                                aligned_waves(k,:) = circshift(waves(k,:), shift, 2);
                            end

                            waves_to_plot = aligned_waves;
                        else
                            waves_to_plot = waves;
                        end

                        % Plot each waveform for this amplitude.
                        plot(t_wave, waves_to_plot', ...
                             'Color', [cmap(a,:) 0.3]);
                    end

                    title(sprintf('%d–%d ms (%d spikes)', ...
                          (b-1)*bin_ms, b*bin_ms, spike_count));

                    xlabel('Waveform time (ms)');
                    ylabel('uV');

                    ylim(y_lim);
                    yticks(linspace(y_lim(1), y_lim(2), 3));
                    xticks(round(linspace(t_wave(1), t_wave(end), 3)));

                    axis square;
                    grid on;
                end

                %% ----- Legend -----
                legend_handles = gobjects(numel(amps_sel),1);
                legend_labels  = cell(numel(amps_sel),1);

                for aa = 1:numel(amps_sel)
                    amp_val = amps_sel(aa);

                    amp_idx_for_color = find(Amps == amp_val, 1, 'first');

                    legend_handles(aa) = plot(nan,nan,'-', ...
                        'Color', cmap(amp_idx_for_color,:), ...
                        'LineWidth',1.5);

                    legend_labels{aa} = sprintf('%g µA', amp_val);
                end

                legend(legend_handles, legend_labels, ...
                       'Location','northeastoutside');

                %% ----- Figure title -----
                sgtitle(sprintf(['Waveform QC | %s | Rec Ch %d | Set %d: %s | ' ...
                                 'Prefix %d | ISI %.1f ms | Last PTD %s'], ...
                    sp_source_label, ich, set_id, stimLabel, ...
                    prefix_val, isi_val, ptd_text), ...
                    'Interpreter','none');

            end % ISI
        end % prefix
    end % set
end % channel