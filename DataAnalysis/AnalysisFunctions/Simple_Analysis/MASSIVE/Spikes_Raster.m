%% Load all spikes and plot a few
load('/Volumes/MACData/Data/Data_Xia/DX035/Xia_100um_Single1_260924_160608/Xia_100um_Single1.sp.mat');

strobe = 0;              % set to 1 to blank first 10 ms
cut_start = 0;           % start of blank window (ms)
cut_end   = 10;          % end of blank window (ms)

for e = 33:64
    rel_spike_times = AllSpikes(e).rel_spike_times(:);  % force column vector
    trial_numbers   = AllSpikes(e).trial_numbers(:);
    nTrials = max(trial_numbers);

    % --- Blank first 10 ms if requested ---
    if strobe
        keep_idx = ~(rel_spike_times >= cut_start & rel_spike_times < cut_end);
        rel_spike_times = rel_spike_times(keep_idx);
        trial_numbers   = trial_numbers(keep_idx);
    end

    % --- Defensive checks ---
    if isempty(rel_spike_times)
        fprintf('No spikes for E%d, skipping\n', e);
        continue;
    end
    if isempty(nTrials) || isnan(nTrials)
        nTrials = 1;
    end

    % --- Compute PSTH ---
    win_ms = 2;
    step_ms = 1;
    edges_start = -pre_ms : step_ms : (post_ms - win_ms);
    fr_sliding = zeros(1, numel(edges_start));

    for i = 1:numel(edges_start)
        window_start = edges_start(i);
        window_end = window_start + win_ms;

        % Count spikes (scalar)
        nSpikes = sum(rel_spike_times >= window_start & rel_spike_times < window_end);
        if numel(nSpikes) > 1
            nSpikes = nSpikes(1);  % just in case
        end

        fr_sliding(i) = nSpikes / (win_ms/1000) / nTrials;
    end
    
    % --- Plot ---
    figure('Name', sprintf('E%d', e), 'Color', 'w');
    subplot(2,1,1);
    plot(rel_spike_times, trial_numbers, 'k.', 'MarkerSize', 4);
    xline(0, 'r--');
    rectangle('Position', [0, 0, stim_dur_ms, nTrials+1], ...
              'FaceColor', [1 0 0 0.2], 'EdgeColor', 'none');
    ylabel('Trial');
    title(sprintf('Raster – E%d', e));
    set(gca,'YDir','normal'); 
    xlim([-50 100]);
    ylim([-5 nTrials+5]);
    nTrials_raster = numel(unique(trial_numbers));
    text(0.98, 0.98, sprintf('%d trials', nTrials_raster), ...
        'Units','normalized', 'HorizontalAlignment','right', ...
        'VerticalAlignment','top', 'FontSize',10, 'Color','k');

    subplot(2,1,2);
    plot(edges_start + win_ms/2, fr_sliding, 'k', 'LineWidth', 1.5);
    xlabel('Time (ms)');
    ylabel('Firing rate (spikes/s)');
    xlim([-50 100]);
end
