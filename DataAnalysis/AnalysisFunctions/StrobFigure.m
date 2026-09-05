%% Load spikes
% load('\\storage.erc.monash.edu\shares\R-MNHS-Syncitium\Shared\Dasha\Analysis\DX012\Strobe_new_251125_211617\Spikes_full.mat');
load('/Volumes/Shared/Dasha/Analysis/DX012/Strobe_new_251125_211617/Spikes_full.mat');

strobe = 0;              
cut_start = 0;           
cut_end   = 10;           

% Parameters
pre_ms  = 50;            
post_ms = 200;           
stim_dur_ms = 0;         
win_ms = 5;              
step_ms = 1;             

for e = 22
    rel_spike_times = AllSpikes(e).rel_spike_times(:);  
    trial_numbers   = AllSpikes(e).trial_numbers(:);
    nTrials = max(trial_numbers);

    if strobe
        keep_idx = ~(rel_spike_times >= cut_start & rel_spike_times < cut_end);
        rel_spike_times = rel_spike_times(keep_idx);
        trial_numbers   = trial_numbers(keep_idx);
    end

    if isempty(rel_spike_times)
        fprintf('No spikes for E%d, skipping\n', e);
        continue;
    end
    if isempty(nTrials) || isnan(nTrials)
        nTrials = 1;
    end

    % --- Compute PSTH ---
    edges_start = -pre_ms : step_ms : (post_ms - win_ms);
    nBins = numel(edges_start);
    fr_sliding = zeros(1, nBins);

    for i = 1:nBins
        window_start = edges_start(i);
        window_end   = window_start + win_ms;
        nSpikes = sum(rel_spike_times >= window_start & rel_spike_times < window_end);
        fr_sliding(i) = nSpikes / (win_ms/1000) / nTrials;  % spikes/s
    end

    % --- Plot figure ---
    figure('Color','w','Position',[100 100 400 300]); hold on;

    % Overlay spikes as dots with small jitter
    spike_jitter_height = 0.1;
    for t = 1:nTrials
        spikes_trial = rel_spike_times(trial_numbers == t);
        jitter_y = t + (rand(size(spikes_trial))-0.5)*spike_jitter_height;  % tiny jitter
        plot(spikes_trial, jitter_y, '.', 'Color', [0 0 0]+0.1, 'MarkerSize', 3);
    end

    % PSTH overlay (scaled visually, Y-axis shows real firing rate)
    yMax_spikes = nTrials + 1;  % height for dots
    psth_scale_factor = yMax_spikes / max(fr_sliding); 
    plot(edges_start + win_ms/2, fr_sliding * psth_scale_factor, 'k', 'LineWidth', 1.6);

    % Adjust Y-ticks to reflect actual firing rate
    yTicks = get(gca,'YTick');
    set(gca,'YTickLabel', round(yTicks / psth_scale_factor));
    ylabel('Firing rate (spikes/s)');

    % Stimulus rectangle (optional)
    if stim_dur_ms>0
        yL = ylim;
        rectangle('Position',[0 yL(1) stim_dur_ms diff(yL)], ...
                  'FaceColor',[0.5 0.5 0.5 0.1],'EdgeColor','none');
    end

    % Labels and aesthetics
    xlabel('Time (ms)');
    %title(sprintf('Raster + PSTH – E%d', e));
    xlim([-pre_ms post_ms]);
    ylim([0 yMax_spikes]);
    set(gca, 'FontName','Arial','FontSize',11,'LineWidth',1,'TickDir','out','Box','off');

    hold off;
end
