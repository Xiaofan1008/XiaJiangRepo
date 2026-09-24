%% Spike detection across all electrodes
%  + separate threshold settings per probe (and per channel if needed)
%  + rejection of bad trials, judged separately for each probe
close all; clear all; clc;
filename = '/Volumes/MACData/Data/Data_Xia/DX035/Xia_150um_SimSeq2_260924_174921/Xia_150um_SimSeq2.mu_sab.dat';
saveFile = '/Volumes/MACData/Data/Data_Xia/DX035/Xia_150um_SimSeq2_260924_174921/Xia_150um_SimSeq2.sp.mat';

% -------- PARAMETERS --------
useAdaptiveThresh = false;  % false = global, true = adaptive
strobe = 0;                 % set 1 if using stimulation trigger blanking
amplifier_gain = 0.195;     % microvolts per bit
snippet_ms = 3;             % total snippet duration
fs = 30000;                 % sampling rate (Hz)
win_ms = 5;                 % PSTH window
step_ms = 1;                % PSTH step
pre_ms = 50;                % pre-trigger window (ms)
post_ms = 300;              % post-trigger window (ms)
stim_dur_ms = 0.5;          % stimulation duration (ms)
artefact_cut_ms = 10;       % (used if strobe=1)
half_snip = round((snippet_ms/2) * fs / 1000);
Electrode_Type = 3; % 0:single shank rigid; 1:single shank flex; 2:four shank flex; 3:64chn+32chn hybrid

% -------- PROBES (map positions) --------
probe_names = {'4-shank (Ports A+B)', 'Single shank (Port C)'};
probe_elecs = {1:64, 65:96};

% -------- THRESHOLD SETTINGS PER PROBE --------
% threshold_method: 'mad' = coeff x robust noise (median|x|/0.6745), not
%                   inflated by artifacts or strong firing
%                   'std' = mean - coeff x SD (original method)
threshold_method = 'mad';
coeff_probe      = [3.5, 3.5];    % threshold multiplier   [4-shank, single shank]
ptp_min_probe    = [20,  20];     % snippet peak-to-peak minimum (uV)
ptp_max_probe    = [400, 400];    % snippet peak-to-peak maximum (uV)
abs_max_probe    = [500, 500];    % reject snippets exceeding this anywhere (uV)

% Per-channel overrides of the multiplier (map position, coeff), applied
% after the probe values. Example: coeff_override = [70 4.5; 71 4.5];
coeff_override = [];
bad_channels   = [];              % map positions to skip entirely (dead / noisy electrodes)

% -------- BAD-TRIAL REJECTION --------
% Every rule compares a trial with the OTHER trials of the same channel
% (robust z-score = distance from the median in robust SDs). The stimulus
% response and the stimulus artifact are in every trial, so they no longer
% count against a trial - only trials that stand out are rejected.
reject_bad_trials = true;
reject_scope      = 'probe';  % 'probe' = each probe keeps its own good trials
                              % 'all'   = a trial bad on either probe is dropped for both
art_uV            = 500;      % artifact: a channel's largest value in the trial must exceed art_uV
art_k             = 6;        %   AND be > art_k robust SDs above that channel's usual trial maximum;
art_frac          = 0.25;     %   trial rejected if this fraction of the probe's channels are flagged
art_ignore_ms     = [0 10];   % ignore this window after the trigger in the artifact check ([] = none)
noise_k           = 4;        % noise: pre-trigger noise > noise_k robust SDs above usual (median over channels)
spike_outlier_k   = 6;        % spike count: probe total > spike_outlier_k robust SDs above usual
rule_max_frac     = 0.3;      % safety: a rule that would reject more than 30% of trials is ignored
show_trial_figure = true;

% -------- Load supporting data --------
trig = loadTrig(0);

% Physical recording channel for each map position (electrode 1..96 in
% probe order), via the shared ProbeMAP/Depth_s convention already used
% throughout the rest of the pipeline -- not a separately hand-built
% shift, so this stays correct automatically if ProbeMAP ever changes.
map_nums_plus = Depth_s(Electrode_Type);   % physical channel per map position

nChannels = numel(map_nums_plus);
nProbes   = numel(probe_elecs);
probe_of  = zeros(1, nChannels);
for p = 1:nProbes
    probe_of(probe_elecs{p}) = p;
end
if any(probe_of == 0)
    error('Every map position must belong to a probe (check probe_elecs).');
end

% Per-channel settings from the probe values + overrides
coeff_ch   = coeff_probe(probe_of);
ptp_min_ch = ptp_min_probe(probe_of);
ptp_max_ch = ptp_max_probe(probe_of);
abs_max_ch = abs_max_probe(probe_of);
for i = 1:size(coeff_override, 1)
    coeff_ch(coeff_override(i,1)) = coeff_override(i,2);
end

% Load all data once (int16 = 4x less memory than double; converted per channel)
fid = fopen(filename, 'r');
data = fread(fid, [nChannels, Inf], 'int16=>int16');
fclose(fid);

pre_samps  = round(pre_ms * fs / 1000);
post_samps = round(post_ms * fs / 1000);
nSamples   = size(data,2);
nTrials    = size(trig,2);

% Trials whose window lies inside the recording
tr_ok  = (trig(1,:) - pre_samps >= 1) & (trig(1,:) + post_samps - 1 <= nSamples);
tr_idx = find(tr_ok);
nV     = numel(tr_idx);
if nV < nTrials
    fprintf('%d trial(s) fall outside the recording and are skipped.\n', nTrials - nV);
end

% Samples of each trial window (same window as the alignment below)
r_offs = -pre_samps : post_samps - 1;
idx_tr = trig(1, tr_idx)' + r_offs;                  % trials x samples
t_rel_ms  = r_offs / fs * 1000;
cols_art  = true(size(r_offs));                     % artifact check window
if strobe
    cols_art = cols_art & ~(t_rel_ms >= 0 & t_rel_ms < artefact_cut_ms);
end
if ~isempty(art_ignore_ms)
    cols_art = cols_art & ~(t_rel_ms >= art_ignore_ms(1) & t_rel_ms < art_ignore_ms(2));
end
cols_base = t_rel_ms < 0;                           % noise measured before the trigger only

% -------- Storage --------
locs_all     = cell(1, nChannels);
thr_ch       = NaN(1, nChannels);
noise_ch     = NaN(1, nChannels);
trial_maxabs = NaN(nChannels, nV);
trial_noise  = NaN(nChannels, nV);

% -------- Loop through electrodes: detection + trial metrics --------
for e_number = 1:nChannels
    if ismember(e_number, bad_channels)
        fprintf('Electrode %d skipped (bad channel)\n', e_number);
        continue
    end
    p  = probe_of(e_number);
    ch = map_nums_plus(e_number);
    channelData = double(data(ch,:)) * amplifier_gain;     % uV
    channelData = channelData - mean(channelData);           % = detrend 'constant'
    robust_noise = median(abs(channelData)) / 0.6745;

    % ----- Threshold for this channel -----
    if strcmpi(threshold_method, 'mad')
        noise  = robust_noise;
        thresh = -coeff_ch(e_number) * noise;
    else
        noise  = std(channelData);
        thresh = mean(channelData) - coeff_ch(e_number) * noise;
    end

    % ----- Spike detection -----
    if ~useAdaptiveThresh
        [~, locs] = findpeaks(-channelData, 'MinPeakHeight', -thresh);
        locs = locs(:);
    else
        window_ms  = 100;
        win_samps  = round(fs * window_ms / 1000);
        step_samps = round(win_samps/2);
        starts     = 1:step_samps:nSamples-win_samps;
        seg_locs   = cell(numel(starts), 1);
        warnState  = warning('off','signal:findpeaks:largeMinPeakHeight');
        for w = 1:numel(starts)
            seg = channelData(starts(w):starts(w)+win_samps-1);
            if strcmpi(threshold_method, 'mad')
                thr = -coeff_ch(e_number) * median(abs(seg)) / 0.6745;
            else
                thr = mean(seg) - coeff_ch(e_number) * std(seg);
            end
            [~, sl] = findpeaks(-seg, 'MinPeakHeight', -thr);
            seg_locs{w} = sl(:) + (starts(w)-1);
        end
        warning(warnState);
        locs = unique(vertcat(seg_locs{:}));
    end

    % Remove edge spikes
    locs = locs(locs > half_snip & locs < (nSamples - half_snip));

    % ----- Amplitude filtering (probe-specific limits, vectorised) -----
    keep  = false(size(locs));
    offs  = -half_snip:half_snip;
    chunk = 50000;
    for i0 = 1:chunk:numel(locs)
        i1    = min(numel(locs), i0 + chunk - 1);
        snips = channelData(locs(i0:i1) + offs);
        ptp   = max(snips, [], 2) - min(snips, [], 2);
        keep(i0:i1) = ptp >= ptp_min_ch(e_number) & ptp <= ptp_max_ch(e_number) & ...
                      all(abs(snips) <= abs_max_ch(e_number), 2);
    end
    locs = locs(keep);

    % ----- Per-trial metrics for bad-trial rejection -----
    seg = channelData(idx_tr);
    trial_maxabs(e_number, :) = max(abs(seg(:, cols_art)), [], 2)';
    trial_noise(e_number, :)  = (median(abs(seg(:, cols_base)), 2) / 0.6745)';
    clear seg

    locs_all{e_number} = locs;
    thr_ch(e_number)   = thresh;
    noise_ch(e_number) = noise;
    fprintf('E%-3d probe %d | coeff %.1f | noise %5.1f uV | thr %6.1f uV | %6d spikes\n', ...
        e_number, p, coeff_ch(e_number), noise, thresh, numel(locs));
end
clear data channelData

% -------- Spike count per trial and channel --------
spk_count = zeros(nChannels, nV);
ws_all    = trig(1, tr_idx) - pre_samps;
we_all    = trig(1, tr_idx) + post_samps - 1;
edges     = reshape([ws_all; we_all + 1], 1, []);
use_hist  = all(diff(edges) > 0);          % windows in time order and not overlapping
for e_number = 1:nChannels
    L = locs_all{e_number};
    if isempty(L), continue; end
    if use_hist
        c = histcounts(L, edges);
        spk_count(e_number, :) = c(1:2:end);
    else
        for ti = 1:nV
            spk_count(e_number, ti) = sum(L >= ws_all(ti) & L <= we_all(ti));
        end
    end
end

% -------- Bad-trial rejection, per probe --------
% robust z-score of each channel's trial values relative to its other trials
robust_z = @(M) (M - median(M, 2, 'omitnan')) ./ ...
                max(1.4826 * median(abs(M - median(M, 2, 'omitnan')), 2, 'omitnan'), eps);

trial_bad_probe = false(nProbes, nV);
tb_art = false(nProbes, nV); tb_noise = false(nProbes, nV); tb_spk = false(nProbes, nV);
noise_z_probe = cell(1, nProbes);
spk_z_probe   = cell(1, nProbes);
fprintf('\n');
for p = 1:nProbes
    chans = setdiff(probe_elecs{p}, bad_channels);
    if isempty(chans), continue; end

    % artifact: large AND unusual for that channel, on many channels at once
    A       = trial_maxabs(chans, :);
    flagged = A > art_uV & robust_z(A) > art_k;
    tb_art(p,:) = sum(flagged, 1) >= max(1, ceil(art_frac * numel(chans)));

    % noise before the trigger, unusual across the probe
    noise_z_probe{p} = median(robust_z(trial_noise(chans, :)), 1, 'omitnan');
    tb_noise(p,:)    = noise_z_probe{p} > noise_k;

    % total spike count on the probe, unusual
    spk_z_probe{p} = robust_z(sum(spk_count(chans, :), 1));
    tb_spk(p,:)    = spk_z_probe{p} > spike_outlier_k;

    % safety: a rule that flags a large share of trials is catching something
    % present in (almost) every trial, not bad trials
    names = {'artifact', 'noise', 'spike count'};
    rules = {tb_art(p,:), tb_noise(p,:), tb_spk(p,:)};
    for r = 1:3
        if mean(rules{r}) > rule_max_frac
            warning('%s: the %s rule would reject %.0f%% of trials - ignored. Loosen it or check the data.', ...
                probe_names{p}, names{r}, 100*mean(rules{r}));
            rules{r}(:) = false;
        end
    end
    [tb_art(p,:), tb_noise(p,:), tb_spk(p,:)] = rules{:};

    if reject_bad_trials
        trial_bad_probe(p,:) = tb_art(p,:) | tb_noise(p,:) | tb_spk(p,:);
    end
    fprintf('%s: artifact %d | noise %d | spike count %d  =>  kept %d of %d trials\n', ...
        probe_names{p}, sum(tb_art(p,:)), sum(tb_noise(p,:)), sum(tb_spk(p,:)), ...
        sum(~trial_bad_probe(p,:)), nV);
end
if strcmpi(reject_scope, 'all')
    trial_bad_probe = repmat(any(trial_bad_probe, 1), nProbes, 1);
    fprintf('reject_scope = ''all'': kept %d of %d trials on both probes\n', sum(~trial_bad_probe(1,:)), nV);
end
good_trials_probe = cell(1, nProbes);
for p = 1:nProbes
    good_trials_probe{p} = tr_idx(~trial_bad_probe(p,:));
end

% -------- Peri-trigger alignment (good trials only) --------
AllSpikes = struct('rel_spike_times', [], 'trial_numbers', [], 'channel', [], ...
                   'probe', [], 'coeff', [], 'threshold_uV', [], 'noise_uV', [], ...
                   'good_trials', [], 'n_good_trials', []);

for e_number = 1:nChannels
    p    = probe_of(e_number);
    good = good_trials_probe{p};
    locs = locs_all{e_number};

    rel_c = cell(numel(good), 1);
    trl_c = cell(numel(good), 1);
    for ti = 1:numel(good)
        t = good(ti);
        win_start = trig(1,t) - pre_samps;
        win_end   = trig(1,t) + post_samps - 1;
        spikes_in_window = locs(locs >= win_start & locs <= win_end);
        spikes_rel_ms = (spikes_in_window - trig(1,t)) / fs * 1000;

        if strobe
            cut_start = 0; cut_end = artefact_cut_ms;
            keep_mask = ~(spikes_rel_ms >= cut_start & spikes_rel_ms < cut_end);
            spikes_rel_ms = spikes_rel_ms(keep_mask);
            shift_mask = spikes_rel_ms >= cut_end;
            spikes_rel_ms(shift_mask) = spikes_rel_ms(shift_mask) - (cut_end - cut_start);
        end

        rel_c{ti} = spikes_rel_ms(:);
        trl_c{ti} = t * ones(numel(spikes_rel_ms), 1);
    end

    AllSpikes(e_number).rel_spike_times = vertcat(rel_c{:});
    AllSpikes(e_number).trial_numbers   = vertcat(trl_c{:});   % original trial numbers
    AllSpikes(e_number).channel         = e_number;
    AllSpikes(e_number).probe           = p;
    AllSpikes(e_number).coeff           = coeff_ch(e_number);
    AllSpikes(e_number).threshold_uV    = thr_ch(e_number);
    AllSpikes(e_number).noise_uV        = noise_ch(e_number);
    AllSpikes(e_number).good_trials     = good;
    AllSpikes(e_number).n_good_trials   = numel(good);
end

% -------- Trial figure --------
if show_trial_figure
    figure('Name', 'Bad-trial rejection', 'Color', 'w', 'Position', [80 150 700*nProbes 450]);
    for p = 1:nProbes
        if isempty(spk_z_probe{p}), continue; end
        ax  = subplot(1, nProbes, p);
        bad = trial_bad_probe(p,:);
        plot(ax, tr_idx, spk_z_probe{p}, 'k.'); hold(ax, 'on');
        plot(ax, tr_idx, noise_z_probe{p}, '.', 'Color', [0.2 0.4 0.9]);
        plot(ax, tr_idx(tb_art(p,:)), zeros(1, sum(tb_art(p,:))), 'm^', 'MarkerFaceColor', 'm');
        plot(ax, tr_idx(bad), spk_z_probe{p}(bad), 'ro');
        yline(ax, spike_outlier_k, 'k--'); yline(ax, noise_k, '--', 'Color', [0.2 0.4 0.9]);
        hold(ax, 'off');
        xlabel(ax, 'Trial'); ylabel(ax, 'Robust z (vs other trials)');
        legend(ax, {'Spike count', 'Pre-trigger noise', 'Artifact', 'Rejected'}, 'Location', 'northwest');
        title(ax, sprintf('%s: %d of %d kept', probe_names{p}, sum(~bad), nV));
    end
end

%% -------- Standard-format outputs for the rest of the pipeline --------
% sp: cell array indexed by PHYSICAL recording channel (the same
% convention used everywhere else), spike times in ms, ALL trials (not
% only the good ones) -- this is what Plot_Raster_Seq_PTDVeried.m and
% the other scripts look for when they load a *.sp.mat file and check
% for a variable called 'sp'.
% NOTE: this detector doesn't keep spike waveform snippets, so sp{ch} is
% Nx1 (spike times only). Anything reading sp{ch}(:,1) works fine;
% anything reading sp{ch}(:,2:end) for waveform-based filtering will
% find nothing there.
sp = cell(1, nChannels);
for e_number = 1:nChannels
    ch = map_nums_plus(e_number);
    sp{ch} = locs_all{e_number}(:) / fs * 1000;   % samples -> ms
end

% BadTrials: cell array indexed by PHYSICAL recording channel, listing
% ORIGINAL trial numbers considered bad -- same filename pattern
% (*.BadTrials.mat) and variable name already used across the rest of
% the pipeline. Rejection here was decided per PROBE (not per channel),
% so every channel on the same probe currently shares the same list.
BadTrials = cell(1, nChannels);
for e_number = 1:nChannels
    ch = map_nums_plus(e_number);
    p  = probe_of(e_number);
    BadTrials{ch} = tr_idx(trial_bad_probe(p,:));
end

% -------- Save once at end --------
% saveFile (*.sp.mat) now also contains 'sp' in the standard format, so
% it can be loaded directly by the rest of the pipeline -- e.g.
% Plot_Raster_Seq_PTDVeried.m's dir('*.sp.mat') + isfield(S,'sp') check
% -- with no changes needed to those other scripts. AllSpikes (with the
% richer per-probe QC info) stays in the same file too.
save(saveFile, 'sp', 'AllSpikes', 'fs', 'pre_ms', 'post_ms', 'stim_dur_ms', ...
     'good_trials_probe', 'trial_bad_probe', 'tb_art', 'tb_noise', 'tb_spk', 'tr_idx', ...
     'probe_names', 'probe_elecs', 'coeff_ch', 'thr_ch', 'threshold_method', '-v7.3');
fprintf('All spike data (including standard-format "sp") saved to %s\n', saveFile);

% Bad trials saved separately, in the *.BadTrials.mat convention already
% used elsewhere in the pipeline.
badTrialsFile = strrep(saveFile, '.sp.mat', '.BadTrials.mat');
save(badTrialsFile, 'BadTrials', '-v7.3');
fprintf('Bad-trial list (standard format) saved to %s\n', badTrialsFile);