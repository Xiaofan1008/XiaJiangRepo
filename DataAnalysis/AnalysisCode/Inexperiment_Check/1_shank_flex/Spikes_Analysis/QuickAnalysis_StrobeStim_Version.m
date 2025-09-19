%% Quick Analysis Code (Strobe Stimulation Version)
clear all
close all
addpath(genpath('./MASSIVE'));
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));

%% ------------------ Parameters ------------------ %%
filepath = pwd;
[amplifier_channels, frequency_parameters] = read_Intan_RHS2000_file;
nChn = size(amplifier_channels, 2);
FS = frequency_parameters.amplifier_sample_rate;

dName = 'amplifier';
vFID = fopen([filepath filesep dName '.dat'], 'r');

window_ms = [-100, 100];
window_samps = round(window_ms / 1000 * FS);
total_sample = diff(window_samps) + 1;
trialLength = total_sample;

if isempty(dir('*.trig.dat'))
    cleanTrig_Strobequick;
end
trig = loadTrig(0);
nTrials = length(trig);
Trial_start = 1;

t_axis = (window_samps(1):window_samps(2)) / FS * 1000;
d = Depth_s(1); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks FlexMUA_all = zeros(nChn, trialLength, nTrials);
Spikes = cell(nChn, nTrials);

preS  = round(10.0e-3 * FS);
postS = round(10.0e-3 * FS);
alignRad = max(1, round(3e-3 * FS));
pp_min_uV   = 30;
pp_max_uV   = 500;
tp_min_ms   = 0.30;
tp_max_ms   = 1.20;
tp_min_samp = round(tp_min_ms/1000 * FS);
tp_max_samp = round(tp_max_ms/1000 * FS);

% Filter
firOrder = 512;
Fpass = [300 3000] / (FS/2);
bpFIR = fir1(firOrder, Fpass, 'bandpass', hamming(firOrder+1), 'scale');
gd = firOrder / 2;

%% ---------- Thresholds ---------- %%
thr = zeros(nChn, 1);
num_samples_threshold = min(floor(10 * FS), floor(dir([filepath filesep dName '.dat']).bytes / (nChn * 2)));
fseek(vFID, 0, 'bof');
data_thres_uV = fread(vFID, [nChn, num_samples_threshold], 'int16') * 0.195;

for ch_idx = 1:nChn
    thisChn = d(ch_idx);
    x = double(data_thres_uV(thisChn, :));
    x_fir = filter(bpFIR, 1, x);
    if gd > 0
        x_fir = [x_fir(gd+1:end), zeros(1, gd)];
    end
    sigma = median(abs(x_fir)) / 0.6745;
    thr(ch_idx) = -4.5 * sigma;
end

%% ---------- Process each trial ---------- %%
for tr = 1:nTrials
    fprintf('Reading trial %d/%d\n', tr, nTrials);
    start_index = trig(tr) + window_samps(1);
    start_byte = start_index * nChn * 2;
    fseek(vFID, start_byte, 'bof');
    data = fread(vFID, [nChn, total_sample], 'int16') * 0.195;

    stim_idx = -window_samps(1) + 1;

    for ch = 1:nChn
        thisChn = d(ch);
        raw = data(thisChn, :);
        sig_filt = filter(bpFIR, 1, raw);
        if gd > 0
            signal = [sig_filt(gd+1:end), zeros(1, gd)];
        else
            signal = sig_filt;
        end
        MUA_all(ch,:,tr) = signal;

        cross_idx = find(signal(2:end) < thr(ch) & signal(1:end-1) >= thr(ch)) + 1;
        refrac_samp = max(1, round(1e-3 * FS));
        if numel(cross_idx) > 1
            keep = [true, diff(cross_idx) > refrac_samp];
            cross_idx = cross_idx(keep);
        end

        trough_idx = [];
        for k = 1:numel(cross_idx)
            i0 = cross_idx(k);
            a = max(1, i0 - preS);
            b = min(total_sample, i0 + postS);
            if (a + preS > b - postS)
                continue
            end
            [~, rel] = min(signal(a:b));
            ip = a + rel - 1;
            trough_idx(end+1,1) = ip;
        end

        if isempty(trough_idx)
            Spikes{thisChn, tr} = struct('idx', [], 't_ms', [], 'wf', zeros(0, preS+postS+1));
        else
            [trough_idx, order] = sort(trough_idx(:));
            vals = signal(trough_idx);
            keep_idx = true(size(trough_idx));
            last_keep = 1;
            for m = 2:numel(trough_idx)
                if trough_idx(m) - trough_idx(last_keep) <= refrac_samp
                    if vals(m) < vals(last_keep)
                        keep_idx(last_keep) = false;
                        last_keep = m;
                    else
                        keep_idx(m) = false;
                    end
                else
                    last_keep = m;
                end
            end
            trough_idx = trough_idx(keep_idx);

            wf_list  = {};
            idx_list = [];
            for k = 1:numel(trough_idx)
                ip = trough_idx(k);
                if ip - preS < 1 || ip + postS > total_sample
                    continue
                end
                wf = signal(ip - preS : ip + postS);
                if max(abs(wf)) > 500
                    continue
                end
                [~, i_trough] = min(wf);
                [~, i_peak_rel] = max(wf(i_trough:end));
                i_peak = i_trough + i_peak_rel - 1;
                is_biphasic = (min(wf) < 0) && (max(wf(i_trough:end)) > 0);
                if ~is_biphasic, continue; end
                pp_amp = max(wf) - min(wf);
                if ~(pp_amp >= pp_min_uV && pp_amp <= pp_max_uV), continue; end
                tp_width = i_peak - i_trough;
                if ~(tp_width >= tp_min_samp && tp_width <= tp_max_samp), continue; end
                wf_list{end+1}  = wf;
                idx_list(end+1) = ip;
            end
            if ~isempty(idx_list)
                wf_mat = cat(1, wf_list{:});
            else
                wf_mat = zeros(0, preS+postS+1);
            end
            spk_times_ms = (idx_list - stim_idx) / FS * 1000;
            Spikes{thisChn, tr} = struct('idx', idx_list, 't_ms', spk_times_ms, 'wf', wf_mat);
        end
    end
end

fclose(vFID);

%% ------------------ Raster + PSTH in Same Subplot per Channel ------------------ %%
figure('Name','Raster + PSTH per Channel','Color','w');

% Define window for plotting
plotWindow = [-5, 10];  % ms

% PSTH settings
binWidth_ms = 1;
bins = plotWindow(1):binWidth_ms:plotWindow(2);
binCenters = (bins(1:end-1) + bins(2:end)) / 2;
nBins = numel(binCenters);

nCols = 8;
nRows = ceil(nChn / nCols);

% Precompute PSTHs and global max firing rate
psth_all = nan(nChn, nBins);
max_FR = 0;

for ch = 1:nChn
    all_spike_times = [];
    for tr = 1:nTrials
        if isempty(Spikes{ch,tr}), continue; end
        spk = Spikes{ch,tr}.t_ms;
        spk = spk(spk >= plotWindow(1) & spk <= plotWindow(2));
        all_spike_times = [all_spike_times, spk];
    end
    if isempty(all_spike_times), continue; end
    counts = histcounts(all_spike_times, bins);
    psth = counts / (nTrials * binWidth_ms * 1e-3);  % Hz
    psth_all(ch,:) = psth;
    max_FR = max(max_FR, max(psth));
end

% Plot each channel
for ch = 1:nChn
    subplot(nRows, nCols, ch); hold on;

    % --- Raster (spikes per trial, same as before) ---
    for tr = 1:nTrials
        if isempty(Spikes{ch,tr}), continue; end
        spk_times = Spikes{ch,tr}.t_ms;
        spk_times = spk_times(spk_times >= plotWindow(1) & spk_times <= plotWindow(2));
        y = tr * ones(size(spk_times));  % plot at trial number
        plot(spk_times, y, 'k.', 'MarkerSize', 2);
    end

    % --- PSTH overlay (firing rate curve) ---
    psth = psth_all(ch,:);
    plot(binCenters, psth, 'Color', [0 0 0.6], 'LineWidth', 2);
    % --- Time zero line ---
    xline(0, 'Color', [0.6 0 0], 'LineWidth', 1.2);
    % Formatting
    title(sprintf('Ch %d', ch));
    xlim(plotWindow);
    ylim([0, max_FR * 1.1]);
    if ch > (nChn - nCols)
        xlabel('Time (ms)');
    end
    if mod(ch-1,nCols)==0
        ylabel('Firing Rate (sps)');
    end
end

sgtitle('Raster + PSTH per Channel');

%% ------------------ Plot Average Spike Waveform Per Channel ------------------ %%
figure('Name','Average Spike Waveforms Per Channel','Color','w');
nCols = 8; % Adjust based on layout
nRows = ceil(nChn / nCols);

% Preallocate and compute average waveforms
AvgWF = nan(nChn, preS + postS + 1);
for ch = 1:nChn
    wf_all = [];
    for tr = 1:nTrials
        if isempty(Spikes{ch,tr}), continue; end
        wf = Spikes{ch,tr}.wf;
        if ~isempty(wf)
            wf_all = [wf_all; wf];  % concatenate spikes
        end
    end
    if ~isempty(wf_all)
        AvgWF(ch,:) = mean(wf_all, 1, 'omitnan');  % average per channel
    end
end

% Determine global Y-limits
ymin = min(AvgWF(:), [], 'omitnan')*1.1;
ymax = max(AvgWF(:), [], 'omitnan')*1.1;

% Time axis
t_wf = (-preS:postS) / FS * 1000;

% Plot
for ch = 1:nChn
    subplot(nRows, nCols, ch);
    if ~isnan(mean(AvgWF(ch,:)))
        plot(t_wf, AvgWF(ch,:), 'LineWidth', 1.2);
        title(sprintf('Ch %d', ch));
        xlabel('Time (ms)');
        ylabel('uV');
        xlim([t_wf(1), t_wf(end)]);
        ylim([ymin ymax]);  % Same scale for all
    else
        title(sprintf('Ch %d (no spikes)', ch));
        ylim([ymin ymax]);  % Keep axis consistent even for empty
    end
end

sgtitle('Average Spike Waveforms by Channel');