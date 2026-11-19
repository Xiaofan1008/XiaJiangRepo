%% Spike detection across all electrodes
close all;clear all;clc;
filename = 'C:\Users\dshu0011\Downloads\Analysis\DX035\exp2_260924_152504\exp2.mu_sab.dat';
saveFile = 'C:\Users\dshu0011\Downloads\Analysis\DX035\exp2_260924_152504\Spikes_mua.mat';

% -------- PARAMETERS --------
useAdaptiveThresh = false;  % false = global, true = adaptive
strobe = 0;                 % set 1 if using stimulation trigger blanking
coeff = 3.5;                % threshold multiplier
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


% -------- Load supporting data --------
trig = loadTrig(0);
E_MAP = ProbeMAP;
map = E_MAP(2:97,7); %E_MAP(2:33,5); % 32 channels %map = E_MAP(2:97,7); bilateral %map = E_MAP(2:33,5) 1shank flexible
convert_label_to_num = @(label) str2double(label(3:end));
map_nums = cellfun(convert_label_to_num, map);
map_nums_plus = map_nums + 1; % shift for MATLAB 1-indexing
if numel(map_nums_plus) > 32
    map_nums_plus(33:end) = map_nums_plus(33:end) + 32;
end
if numel(map_nums_plus) > 64
    map_nums_plus(65:end) = map_nums_plus(65:end) + 32;
end
if numel(map_nums_plus) > 96
    map_nums_plus(97:end) = map_nums_plus(97:end) + 32;
end


% Load all data once
fid = fopen(filename, 'r');
data = fread(fid, [length(map), Inf], 'int16=>double');
fclose(fid);

nChannels = length(map);
pre_samps = round(pre_ms * fs / 1000);
post_samps = round(post_ms * fs / 1000);
nSamples = size(data,2);
nTrials = size(trig,2);

% -------- Storage structure --------
AllSpikes = struct('rel_spike_times', [], 'trial_numbers', [], 'channel', []);

% -------- Loop through electrodes --------
for e_number = 1:nChannels
    fprintf('Processing electrode %d / %d ...\n', e_number, nChannels);

    ch = map_nums_plus(e_number);
    channelData = detrend(data(ch,:), 'constant');

    % ----- Spike detection -----
    if ~useAdaptiveThresh
        thresh = mean(channelData) - coeff * std(channelData);
        [~, locs] = findpeaks(-channelData, ...
            'MinPeakHeight', -thresh);
    else
        window_ms = 100;
        win_samps = round(fs * window_ms / 1000);
        step_samps = round(win_samps/2);
        locs = [];
        for w = 1:step_samps:nSamples-win_samps
            seg = channelData(w:w+win_samps-1);
            thr = mean(seg) - coeff * std(seg);
            warnState = warning('off','signal:findpeaks:largeMinPeakHeight');
            [~, seg_locs] = findpeaks(-seg, ...
                'MinPeakHeight', -thr);
            seg_locs = seg_locs + (w-1);
            locs = [locs; seg_locs(:)];
        end
        locs = unique(locs);
    end

    % Remove edge spikes
    valid = locs > half_snip & locs < (nSamples - half_snip);
    locs = locs(valid);

    % ----- Amplitude filtering -----
    keep = false(size(locs));
    for i = 1:length(locs)
        snip = channelData(locs(i)-half_snip : locs(i)+half_snip);
        snip_uV = snip * amplifier_gain;
        ptp_amp = max(snip_uV) - min(snip_uV);
        if ptp_amp >= 20 && ptp_amp <= 400 && all(abs(snip_uV) <= 500)
            keep(i) = true;
        end
    end
    locs = locs(keep);

    % ----- Peri-trigger alignment -----
    rel_spike_times = [];
    trial_numbers   = [];

    for t = 1:nTrials
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

        rel_spike_times = [rel_spike_times; spikes_rel_ms(:)];
        trial_numbers   = [trial_numbers; t * ones(length(spikes_rel_ms), 1)];
    end

    % Store in struct
    AllSpikes(e_number).rel_spike_times = rel_spike_times;
    AllSpikes(e_number).trial_numbers   = trial_numbers;
    AllSpikes(e_number).channel         = e_number;
end

% -------- Save once at end --------
save(saveFile, 'AllSpikes', 'fs', 'pre_ms', 'post_ms', 'stim_dur_ms', '-v7.3');
fprintf('✅ All spike data saved to %s\n', saveFile);