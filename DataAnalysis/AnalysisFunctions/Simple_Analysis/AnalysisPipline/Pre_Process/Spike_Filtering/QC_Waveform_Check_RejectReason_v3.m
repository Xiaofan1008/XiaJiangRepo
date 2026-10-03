%% ============================================================
% QC_FilterCheck_Waveforms.m
% Sanity-check for SpikeFiltering_Evoked_Cleanup.m: for each channel,
% overlay the waveforms that were REJECTED (split by which stage caught
% them -- amplitude ceiling / zero-crossing / PCA) against the ones that
% were KEPT, so you can see by eye whether the filter is throwing out
% real spikes or keeping garbage.
%
% Reads ONLY the compact "<base>.sp_xia_rejectlog.mat" file saved by
% SpikeFiltering_Evoked_Cleanup.m -- NOT the raw or filtered waveforms
% files. That file already holds, per channel: the exact counts/
% percentages for every group (kept/amplitude/zero_crossing/pca) and a
% capped random subsample of example waveforms per group (enough to
% redraw these plots), so there's no need to load the full waveform
% files here at all.
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
rejectlog_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_rejectlog.mat';

chns_to_plot = [10 24 44 51 52 58 60 62];
                    % channels to produce DETAIL figures for -- any list,
                    % not necessarily contiguous or in order, e.g.
                    % [10 24 44], or 1:64 for every channel, or
                    % [48:52 60 62] to mix a range with specific channels.
                    % (the summary bar chart below always covers every
                    % channel regardless of this setting)

FS = 30000;         % sampling frequency (overridden below if 'fs' is
                    % saved inside the file)

%% ================= LOAD COMPACT REJECT LOG =================
assert(isfile(rejectlog_file), 'Cannot find %s.', rejectlog_file);
S = load(rejectlog_file);
assert(isfield(S, 'reject_summary') && isfield(S, 'reject_examples') && isfield(S, 'group_names'), ...
    ['%s does not look like a compact reject log (missing reject_summary/' ...
     'reject_examples/group_names) -- re-run the current SpikeFiltering_Evoked_Cleanup.m ' ...
     'to produce one.'], rejectlog_file);

reject_summary  = S.reject_summary;
reject_examples = S.reject_examples;
group_names     = S.group_names;      % e.g. {'kept','amplitude','zero_crossing','pca'}
nCh = numel(reject_summary);

if isfield(S, 'fs'), FS = S.fs; end

if isfield(S, 'QC_params')
    QC_params = S.QC_params;
    fprintf('Loaded QC_params (amp_ceiling=%g uV, pca_alpha=%s).\n', ...
        QC_params.amp_ceiling_uV, mat2str(QC_params.pca_alpha_probe));
else
    QC_params = [];
end

% Index of 'kept' vs the rejected reason groups within group_names.
kept_g = find(strcmp(group_names, 'kept'));
assert(~isempty(kept_g), '"kept" group not found in group_names.');
reason_g = setdiff(1:numel(group_names), kept_g, 'stable');

% Display label/color for each rejected-reason group, matched by name so
% this doesn't depend on column order in group_names.
reason_label_map = struct('amplitude','amplitude ceiling', ...
                           'zero_crossing','zero-crossing', ...
                           'pca','PCA/Mahalanobis');
reason_color_map = struct('amplitude',   [0.85 0.55 0.10], ...   % orange
                           'zero_crossing', [0.55 0.55 0.55], ... % gray
                           'pca',         [0.60 0.10 0.70]);      % purple

%% ================= SUMMARY: % REJECTED PER CHANNEL =================
nKept     = arrayfun(@(s) s.n_kept, reject_summary)';
nRaw      = arrayfun(@(s) s.n_raw,  reject_summary)';
nRejected = nRaw - nKept;
pctRejected = 100 * nRejected ./ max(nRaw, 1);

figure('Name','QC Filter Summary -- % spikes rejected per channel','Color','w', ...
       'Position',[100 100 1200 400]);
bar(1:nCh, pctRejected, 'FaceColor',[0.85 0.33 0.1]);
xlabel('Channel'); ylabel('% rejected');
title('Fraction of spikes removed by SpikeFiltering\_Evoked\_Cleanup, per channel', 'Interpreter','none');
xlim([0 nCh+1]); grid on;
for ch = 1:nCh
    if nRaw(ch) == 0, continue; end
    text(ch, pctRejected(ch)+1, sprintf('%d', nRejected(ch)), ...
        'HorizontalAlignment','center', 'FontSize', 7);
end

fprintf('\nChannel summary (kept / rejected / %% rejected, by reason):\n');
for ch = 1:nCh
    if nRaw(ch) == 0, continue; end
    reason_str = '';
    for r = reason_g
        gname = group_names{r};
        n_g = reject_summary(ch).(['n_' gname]);
        if n_g == 0, continue; end
        reason_str = [reason_str, sprintf('%s=%d ', gname, n_g)]; %#ok<AGROW>
    end
    fprintf('Ch %3d: %5d kept, %5d rejected (%.1f%%) [%s]\n', ...
        ch, nKept(ch), nRejected(ch), pctRejected(ch), strtrim(reason_str));
end

%% ================= DETAIL FIGURES: KEPT vs REJECTED WAVEFORMS, BY REASON =================
wf_len = [];
for ch = 1:nCh
    for g = 1:numel(group_names)
        if ~isempty(reject_examples{ch,g})
            wf_len = size(reject_examples{ch,g}, 2);
            break;
        end
    end
    if ~isempty(wf_len), break; end
end
assert(~isempty(wf_len), 'No example waveforms found in the reject log -- nothing to plot.');
t_wave = (0:wf_len-1) / FS * 1000;

for ch = chns_to_plot
    if ch > nCh, continue; end
    if nRaw(ch) == 0, continue; end

    figure('Name', sprintf('Ch %d -- kept vs rejected waveforms (%d kept, %d rejected)', ...
                            ch, nKept(ch), nRejected(ch)), ...
           'Color','w','Position',[100 100 900 500]);
    hold on;

    % -------- Rejected examples, split by reason (plot first, so kept draws on top) --------
    for r = reason_g
        gname = group_names{r};
        rw = reject_examples{ch, r};
        if isempty(rw), continue; end
        color = reason_color_map.(gname);
        plot(t_wave, rw', 'Color', [color 0.25], 'LineWidth', 0.5);
    end

    % -------- Kept examples --------
    kw = reject_examples{ch, kept_g};
    if ~isempty(kw)
        plot(t_wave, kw', 'Color', [0.1 0.4 0.8 0.15], 'LineWidth', 0.5);
    end

    % -------- Mean traces on top (computed from the saved EXAMPLES; true
    % counts for the legend/title come from reject_summary, not from how
    % many examples happened to be saved) --------
    mean_handles = [];
    mean_labels  = {};
    for r = reason_g
        gname = group_names{r};
        rw = reject_examples{ch, r};
        if isempty(rw), continue; end
        n_g = reject_summary(ch).(['n_' gname]);
        h = plot(t_wave, mean(rw,1), 'Color', reason_color_map.(gname), 'LineWidth', 2);
        mean_handles(end+1) = h; %#ok<SAGROW>
        mean_labels{end+1}  = sprintf('mean rejected: %s (n=%d)', reason_label_map.(gname), n_g); %#ok<SAGROW>
    end
    if ~isempty(kw)
        h = plot(t_wave, mean(kw,1), 'Color', [0 0.2 0.6], 'LineWidth', 2);
        mean_handles = [h, mean_handles]; %#ok<AGROW>
        mean_labels  = [{sprintf('mean kept (n=%d)', nKept(ch))}, mean_labels]; %#ok<AGROW>
    end

    % -------- Amplitude ceiling reference line, if known --------
    if ~isempty(QC_params) && isfield(QC_params,'amp_ceiling_uV')
        yline(QC_params.amp_ceiling_uV,  'k--');
        yline(-QC_params.amp_ceiling_uV, 'k--');
    end

    xlabel('Time (ms)'); ylabel('\muV');
    title(sprintf('Ch %d: %d kept (blue) vs %d rejected, by reason', ch, nKept(ch), nRejected(ch)));
    if ~isempty(mean_handles)
        legend(mean_handles, mean_labels, 'Location','northeastoutside');
    end
    grid on; box off;
end