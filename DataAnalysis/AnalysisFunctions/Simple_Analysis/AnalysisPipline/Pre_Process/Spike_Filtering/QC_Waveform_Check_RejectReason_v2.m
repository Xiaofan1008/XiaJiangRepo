%% ============================================================
% QC_FilterCheck_Waveforms.m
% Sanity-check for SpikeFiltering_Evoked_Cleanup.m: for each channel,
% overlay the waveforms that were REJECTED against the ones that were
% KEPT, so you can see by eye whether the amplitude/zero-crossing/PCA
% filter is throwing out real spikes or keeping garbage.
%
% Enter the RAW waveforms file and the FILTERED waveforms file directly
% (no folder entry, no guessing) -- same convention as
% SpikeFiltering_Evoked_Cleanup.m. "Rejected" = any spike present in the
% raw file's time column but not in the filtered file's time column
% (matched by exact spike time, since filtering only ever removes rows,
% never changes the surviving ones).
% ============================================================
clear all;clc;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
raw_wf_file      = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_waveforms.mat';
filtered_wf_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_waveforms_filtered.mat';

% Optional: per-spike reject-reason log saved by SpikeFiltering_Evoked_Cleanup.m
% (<base>.sp_xia_rejectlog.mat). When given, rejected waveforms are
% color-split by which stage removed them (amplitude/zero-crossing/PCA)
% instead of a single generic "rejected" group. Leave as '' to skip this
% and fall back to the old lumped-together behaviour (e.g. for filtered
% files produced before this field existed).
rejectlog_file   = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_rejectlog.mat';

chns_to_plot = [10 24 44 51 52 58 60 62];
                    % channels to produce DETAIL figures for -- any list,
                    % not necessarily contiguous or in order, e.g.
                    % [10 24 44], or 1:64 for every channel, or
                    % [48:52 60 62] to mix a range with specific channels.
                    % (the summary bar chart above always covers every
                    % channel regardless of this setting)

max_wave_per_group = 300;  % if a channel has more than this many kept or
                            % rejected spikes, randomly subsample down to
                            % this many for plotting only (keeps figures
                            % readable/fast; counts/title still show the
                            % true totals)

FS = 30000;         % sampling frequency (overridden below if 'fs' is
                    % saved inside the files)

%% ================= LOAD RAW + FILTERED =================
assert(isfile(raw_wf_file), 'Cannot find %s.', raw_wf_file);
S_raw = load(raw_wf_file);
assert(isfield(S_raw, 'sp_waveforms'), 'Variable "sp_waveforms" not found in %s.', raw_wf_file);
wf_raw = S_raw.sp_waveforms;

assert(isfile(filtered_wf_file), 'Cannot find %s.', filtered_wf_file);
S_filt = load(filtered_wf_file);
assert(isfield(S_filt, 'sp_waveforms'), 'Variable "sp_waveforms" not found in %s.', filtered_wf_file);
wf_filt = S_filt.sp_waveforms;

if isfield(S_filt, 'fs')
    FS = S_filt.fs;
elseif isfield(S_raw, 'fs')
    FS = S_raw.fs;
end

if isfield(S_filt, 'QC_params')
    QC_params = S_filt.QC_params;
    fprintf('Loaded QC_params from filtered file (amp_ceiling=%g uV, pca_alpha=%s).\n', ...
        QC_params.amp_ceiling_uV, mat2str(QC_params.pca_alpha_probe));
else
    QC_params = [];
end

nCh = numel(wf_raw);
assert(numel(wf_filt) == nCh, ...
    'Raw (%d channels) and filtered (%d channels) waveform files disagree -- check these are the same dataset.', ...
    nCh, numel(wf_filt));

% -------- Optional reject-reason log --------
have_reject_log = ~isempty(rejectlog_file);
if have_reject_log
    assert(isfile(rejectlog_file), 'Cannot find %s.', rejectlog_file);
    S_log = load(rejectlog_file);
    assert(isfield(S_log, 'reject_log'), 'Variable "reject_log" not found in %s.', rejectlog_file);
    reject_log = S_log.reject_log;
    assert(numel(reject_log) == nCh, ...
        'Reject log (%d channels) does not match raw waveforms (%d channels) -- check these are the same dataset.', ...
        numel(reject_log), nCh);
    fprintf('Loaded per-spike reject log from %s -- rejected waveforms will be split by reason.\n', rejectlog_file);
end

% Reason groups plotted for rejected spikes, in draw order (first drawn
% = furthest back) and their display colors. 'kept' never appears here.
reason_names  = {'amplitude', 'zero_crossing', 'pca'};
reason_labels = {'amplitude ceiling', 'zero-crossing', 'PCA/Mahalanobis'};
reason_colors = [0.85 0.55 0.10;   % amplitude   - orange
                 0.55 0.55 0.55;   % zero_crossing - gray
                 0.60 0.10 0.70];  % pca         - purple

%% ================= PER-CHANNEL KEPT/REJECTED SPLIT =================
nKept     = zeros(nCh,1);
nRejected = zeros(nCh,1);
kept_wf     = cell(nCh,1);
rejected_wf = cell(nCh,1);              % used only when no reject log
rejected_wf_by_reason = cell(nCh, numel(reason_names));  % used when reject log present

for ch = 1:nCh
    raw_ch  = wf_raw{ch};
    filt_ch = wf_filt{ch};
    if isempty(raw_ch), continue; end

    raw_times = raw_ch(:,1);
    if isempty(filt_ch)
        is_kept = false(size(raw_times));
    else
        filt_times = filt_ch(:,1);
        is_kept = ismember(raw_times, filt_times);
    end

    kept_wf{ch} = raw_ch(is_kept, 2:end);
    nKept(ch)   = sum(is_kept);

    if have_reject_log
        log_ch = reject_log{ch};
        assert(isequal(log_ch.time, raw_times), ...
            'Ch %d: reject log spike times do not match the raw waveforms file -- not from the same run?', ch);
        reason_ch = log_ch.reason;
        for r = 1:numel(reason_names)
            is_r = strcmp(reason_ch, reason_names{r});
            rejected_wf_by_reason{ch,r} = raw_ch(is_r, 2:end);
        end
        nRejected(ch) = sum(~is_kept);
    else
        rejected_wf{ch} = raw_ch(~is_kept, 2:end);
        nRejected(ch)   = sum(~is_kept);
    end
end

pctRejected = 100 * nRejected ./ max(nKept + nRejected, 1);

%% ================= SUMMARY: % REJECTED PER CHANNEL =================
figure('Name','QC Filter Summary -- % spikes rejected per channel','Color','w', ...
       'Position',[100 100 1200 400]);
bar(1:nCh, pctRejected, 'FaceColor',[0.85 0.33 0.1]);
xlabel('Channel'); ylabel('% rejected');
title('Fraction of spikes removed by SpikeFiltering\_Evoked\_Cleanup, per channel', 'Interpreter','none');
xlim([0 nCh+1]); grid on;
for ch = 1:nCh
    if nKept(ch)+nRejected(ch) == 0, continue; end
    text(ch, pctRejected(ch)+1, sprintf('%d', nRejected(ch)), ...
        'HorizontalAlignment','center', 'FontSize', 7);
end
fprintf('\nChannel summary (kept / rejected / %% rejected):\n');
for ch = 1:nCh
    if nKept(ch)+nRejected(ch) == 0, continue; end
    fprintf('Ch %3d: %5d kept, %5d rejected (%.1f%%)\n', ch, nKept(ch), nRejected(ch), pctRejected(ch));
end

%% ================= DETAIL FIGURES: KEPT vs REJECTED WAVEFORMS =================
wf_len = [];
for ch = 1:nCh
    if ~isempty(kept_wf{ch}), wf_len = size(kept_wf{ch},2); break; end
    if have_reject_log
        for r = 1:numel(reason_names)
            if ~isempty(rejected_wf_by_reason{ch,r}), wf_len = size(rejected_wf_by_reason{ch,r},2); break; end
        end
    elseif ~isempty(rejected_wf{ch})
        wf_len = size(rejected_wf{ch},2);
    end
    if ~isempty(wf_len), break; end
end
assert(~isempty(wf_len), 'No waveforms found in either file -- nothing to plot.');
t_wave = (0:wf_len-1) / FS * 1000;

for ch = chns_to_plot
    if ch > nCh, continue; end
    if nKept(ch) == 0 && nRejected(ch) == 0, continue; end

    figure('Name', sprintf('Ch %d -- kept vs rejected waveforms (%d kept, %d rejected)', ...
                            ch, nKept(ch), nRejected(ch)), ...
           'Color','w','Position',[100 100 900 500]);
    hold on;

    mean_handles = [];
    mean_labels  = {};

    if have_reject_log
        % -------- Rejected, split by reason (plot first, so kept draws on top) --------
        reason_counts = zeros(1, numel(reason_names));
        for r = 1:numel(reason_names)
            rw = rejected_wf_by_reason{ch,r};
            reason_counts(r) = size(rw,1);
            if isempty(rw), continue; end
            rw_plot = rw;
            if size(rw_plot,1) > max_wave_per_group
                samp = randperm(size(rw_plot,1), max_wave_per_group);
                rw_plot = rw_plot(samp,:);
            end
            plot(t_wave, rw_plot', 'Color', [reason_colors(r,:) 0.25], 'LineWidth', 0.5);
        end
    else
        % -------- Rejected, lumped together (no reject log available) --------
        rw = rejected_wf{ch};
        if ~isempty(rw)
            if size(rw,1) > max_wave_per_group
                samp = randperm(size(rw,1), max_wave_per_group);
                rw = rw(samp,:);
            end
            plot(t_wave, rw', 'Color', [0.6 0.6 0.6 0.25], 'LineWidth', 0.5);
        end
    end

    % -------- Kept --------
    kw = kept_wf{ch};
    if ~isempty(kw)
        if size(kw,1) > max_wave_per_group
            samp = randperm(size(kw,1), max_wave_per_group);
            kw = kw(samp,:);
        end
        plot(t_wave, kw', 'Color', [0.1 0.4 0.8 0.15], 'LineWidth', 0.5);
    end

    % -------- Mean traces on top --------
    if have_reject_log
        for r = 1:numel(reason_names)
            rw = rejected_wf_by_reason{ch,r};
            if isempty(rw), continue; end
            h = plot(t_wave, mean(rw,1), 'Color', reason_colors(r,:), 'LineWidth', 2);
            mean_handles(end+1) = h; %#ok<SAGROW>
            mean_labels{end+1}  = sprintf('mean rejected: %s (n=%d)', reason_labels{r}, reason_counts(r)); %#ok<SAGROW>
        end
    elseif ~isempty(rejected_wf{ch})
        h = plot(t_wave, mean(rejected_wf{ch},1), 'Color', [0.6 0 0], 'LineWidth', 2);
        mean_handles(end+1) = h;
        mean_labels{end+1}  = 'mean rejected';
    end
    if ~isempty(kept_wf{ch})
        h = plot(t_wave, mean(kept_wf{ch},1), 'Color', [0 0.2 0.6], 'LineWidth', 2);
        mean_handles = [h, mean_handles]; %#ok<AGROW>
        mean_labels  = [{'mean kept'}, mean_labels]; %#ok<AGROW>
    end

    % -------- Amplitude ceiling reference line, if known --------
    if ~isempty(QC_params) && isfield(QC_params,'amp_ceiling_uV')
        yline(QC_params.amp_ceiling_uV,  'k--');
        yline(-QC_params.amp_ceiling_uV, 'k--');
    end

    xlabel('Time (ms)'); ylabel('\muV');
    if have_reject_log
        title(sprintf('Ch %d: %d kept (blue) vs %d rejected, by reason', ch, nKept(ch), nRejected(ch)));
    else
        title(sprintf('Ch %d: %d kept (blue), %d rejected (gray/red mean)', ch, nKept(ch), nRejected(ch)));
    end
    if ~isempty(mean_handles)
        legend(mean_handles, mean_labels, 'Location','northeastoutside');
    end
    grid on; box off;
end