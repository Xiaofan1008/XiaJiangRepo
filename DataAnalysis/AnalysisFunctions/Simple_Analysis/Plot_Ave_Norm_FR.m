% === Load Data ===
load('FR_heatmap_norm_single.mat');
load('FR_heatmap_norm_sim.mat');
load('FR_heatmap_norm_seq.mat');

Amps = [0, 1, 2, 3, 4, 5];

%% === SINGLE STIMULATION ===
[nChn, n_AMP, nSets_single] = size(FR_heatmap_norm_single);
FR_summary_single = struct();

for si = 1:nSets_single
    mean_FR_single = zeros(1, n_AMP);
    sem_FR_single  = zeros(1, n_AMP);

    for ai = 1:n_AMP
        data = FR_heatmap_norm_single(:, ai, si);  % use all data (no exclusion)
        if ~isempty(data)
            mean_FR_single(ai) = mean(data, 'omitnan');
            sem_FR_single(ai) = std(data, 'omitnan') / sqrt(sum(~isnan(data)));
        end
    end

    % Save result
    set_name = sprintf('set%d', si);
    FR_summary_single.(set_name).mean = mean_FR_single;
    FR_summary_single.(set_name).sem = sem_FR_single;

    % Plot
    figure('Color','w', 'Name', sprintf('Single Stim Set %d', si));
    errorbar(Amps, mean_FR_single, sem_FR_single, '-o', 'LineWidth', 1.5);
    title(sprintf('Single Stimulation Set %d', si));
    xlabel('Amplitude (µA)');
    ylabel('Normalized Firing Rate');
    legend('Single', 'Location', 'northwest');
    ylim([0, max(mean_FR_single + sem_FR_single) + 0.1]);
    box off; 
end


%% === SIMULTANEOUS STIMULATION ===
[~, ~, nSets_sim] = size(FR_heatmap_norm_sim);
FR_summary_sim = struct();

for si = 1:nSets_sim
    mean_FR_sim = zeros(1, n_AMP);
    sem_FR_sim  = zeros(1, n_AMP);

    for ai = 1:n_AMP
        data = FR_heatmap_norm_sim(:, ai, si);
        if ~isempty(data)
            mean_FR_sim(ai) = mean(data, 'omitnan');
            sem_FR_sim(ai) = std(data, 'omitnan') / sqrt(sum(~isnan(data)));
        end
    end

    set_name = sprintf('set%d', si);
    FR_summary_sim.(set_name).mean = mean_FR_sim;
    FR_summary_sim.(set_name).sem = sem_FR_sim;

    figure('Color','w', 'Name', sprintf('Simultaneous Stim Set %d', si));
    errorbar(Amps, mean_FR_sim, sem_FR_sim, '-s', 'LineWidth', 1.5);
    title(sprintf('Simultaneous Stimulation Set %d', si));
    xlabel('Amplitude (µA)');
    ylabel('Normalized Firing Rate');
    legend('Simultaneous', 'Location', 'northwest');
    ylim([0, max(mean_FR_sim + sem_FR_sim) + 0.1]);
    box off; 
end


%% === SEQUENTIAL STIMULATION ===
[~, ~, nSets_seq] = size(FR_heatmap_norm_seq);
FR_summary_seq = struct();

for si = 1:nSets_seq
    mean_FR_seq = zeros(1, n_AMP);
    sem_FR_seq  = zeros(1, n_AMP);

    for ai = 1:n_AMP
        data = FR_heatmap_norm_seq(:, ai, si);
        if ~isempty(data)
            mean_FR_seq(ai) = mean(data, 'omitnan');
            sem_FR_seq(ai) = std(data, 'omitnan') / sqrt(sum(~isnan(data)));
        end
    end

    set_name = sprintf('set%d', si);
    FR_summary_seq.(set_name).mean = mean_FR_seq;
    FR_summary_seq.(set_name).sem = sem_FR_seq;

    figure('Color','w', 'Name', sprintf('Sequential Stim Set %d', si));
    errorbar(Amps, mean_FR_seq, sem_FR_seq, '-^', 'LineWidth', 1.5);
    title(sprintf('Sequential Stimulation Set %d', si));
    xlabel('Amplitude (µA)');
    ylabel('Normalized Firing Rate');
    legend('Sequential', 'Location', 'northwest');
    ylim([0, max(mean_FR_seq + sem_FR_seq) + 0.1]);
    box off;
end