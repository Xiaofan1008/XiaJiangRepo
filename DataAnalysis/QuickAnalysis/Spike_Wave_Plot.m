%% Plot Parameters
win_ms = 300;         % total time window after each trigger (ms)
bin_ms = 10;          % width of each time bin (ms)
nBins  = win_ms / bin_ms;
nChn   = numel(sp);   % total number of channels
amp_threshold = 150;  % max absolute amplitude allowed (µV)


%% Experiment Parameter Load %%
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
fileDIR = fileDIR(1).name;
S = load(fileDIR,'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
CHN = S.CHN;
E_MAP = S.E_MAP;
n_Trials = S.n_Trials;

%% Plot 

for ch = 1:nChn
    if isempty(sp{ch}), continue; end

    fprintf('Processing Channel %d...\n', ch);

    % Spike data
    sp_times = sp{ch}(:,1);        % Spike times (ms)
    sp_wave  = sp{ch}(:,2:end);    % Spike waveforms

    % === Remove spikes exceeding amplitude threshold ===
    valid_idx = all(abs(sp_wave) <= amp_threshold, 2);
    sp_times  = sp_times(valid_idx);
    sp_wave   = sp_wave(valid_idx,:);

    if isempty(sp_times)
        fprintf('All spikes exceeded threshold for channel %d — skipping.\n', ch);
        continue;
    end

    % Initialize: each bin collects spikes across all trials
    all_spikes_by_bin = cell(nBins, 1);

    % Loop over all triggers
    for i = 1:length(trig)
        t0_ms = trig(i) / FS * 1000;              % Convert to ms
        t_window_idx = sp_times >= t0_ms & sp_times < (t0_ms + win_ms);
        if ~any(t_window_idx), continue; end

        rel_times = sp_times(t_window_idx) - t0_ms;   % Time relative to trigger
        waveforms = sp_wave(t_window_idx, :);

        % Assign each spike to a time bin
        for j = 1:length(rel_times)
            bin_idx = floor(rel_times(j) / bin_ms) + 1;
            if bin_idx >= 1 && bin_idx <= nBins
                all_spikes_by_bin{bin_idx}(end+1,:) = waveforms(j,:);
            end
        end
    end

    % === Determine common Y limit for this channel ===
    all_waves_cat = cell2mat(all_spikes_by_bin(~cellfun('isempty',all_spikes_by_bin)));
    if isempty(all_waves_cat)
        fprintf('No valid spikes after filtering for channel %d\n', ch);
        continue;
    end
    y_max = max(abs(all_waves_cat(:)));
    y_lim = [-1 1] * ceil(y_max/50)*50;   % round to nearest 50 µV

    % === Plot all bins for this channel ===
    figure('Name', sprintf('Channel %d — Spikes by 10 ms Bins', ch), ...
           'Color','w','Position', [100 100 1400 800]);
    tiledlayout(6, 5, 'Padding','compact', 'TileSpacing','compact');

    for b = 1:nBins
        nexttile; hold on;
        bin_spikes = all_spikes_by_bin{b};

        if isempty(bin_spikes)
            title(sprintf('%d–%d ms (0 spikes)', (b-1)*bin_ms, b*bin_ms));
            axis off;
            continue;
        end

        % Align spikes by their min value
        aligned_waves = zeros(size(bin_spikes));
        for k = 1:size(bin_spikes,1)
            [~, min_idx] = min(bin_spikes(k,:));
            shift = ceil(size(bin_spikes,2)/2) - min_idx;
            aligned_waves(k,:) = circshift(bin_spikes(k,:), shift, 2);
        end

        % Time axis for waveform (in ms)
        t_wave = (0:size(sp_wave,2)-1) / FS * 1000;

        % Plot aligned waveforms
        plot(t_wave, aligned_waves', 'Color', [0.2 0.2 0.2 0.3]);
        title(sprintf('%d–%d ms (%d spikes)', ...
                      (b-1)*bin_ms, b*bin_ms, size(aligned_waves,1)));

        xlabel('Time (ms)');
        ylabel('Voltage (µV)');
        axis square;
        grid on;

        % Consistent axis limits & ticks for this channel
        ylim(y_lim);
        yticks(linspace(y_lim(1), y_lim(2), 3));  % e.g., -max, 0, +max
        xticks(round(linspace(t_wave(1), t_wave(end), 3)));
    end
end