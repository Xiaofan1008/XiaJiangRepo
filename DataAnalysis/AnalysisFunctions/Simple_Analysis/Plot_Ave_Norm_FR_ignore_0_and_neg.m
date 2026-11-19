% Plot the average normalized firing rate 

load('FR_heatmap_norm_sim.mat')
load('FR_heatmap_norm_seq.mat')
load('FR_heatmap_norm_single.mat')

%% === SINGLE STIMULATION ===
[nChn, n_AMP, nSets_single] = size(FR_heatmap_norm_single);
FR_summary_single = struct();
Amps = [0,1,2,3,4,5];
for si = 1:nSets_single
    mean_FR = zeros(1, n_AMP);
    sem_FR  = zeros(1, n_AMP);

    for ai = 1:n_AMP
        data = FR_heatmap_norm_single(:, ai, si);
        valid = data(data > 0);  % exclude zero and negative
        if ~isempty(valid)
            mean_FR(ai) = mean(valid);
            sem_FR(ai) = std(valid) / sqrt(length(valid));
        end
    end

    % Save result
    set_name = sprintf('set%d', si);
    FR_summary_single.(set_name).mean = mean_FR;
    FR_summary_single.(set_name).sem = sem_FR;

    % Plot
    figure('Color','w', 'Name', sprintf('Single Stim Set %d', si));
    errorbar(Amps, mean_FR, sem_FR, '-o', 'LineWidth', 2);
    title(sprintf('Single Stimulation Set %d', si));
    xlabel('Amplitude (µA)');
    ylabel('Normalized Firing Rate');
    legend('Single', 'Location', 'northwest');
    ylim([0, 1.1 * max(mean_FR + sem_FR)]);
    box off;
end


%% === SIMULTANEOUS STIMULATION ===
[~, ~, nSets_sim] = size(FR_heatmap_norm_sim);
FR_summary_sim = struct();

for si = 1:nSets_sim
    mean_FR = zeros(1, n_AMP);
    sem_FR  = zeros(1, n_AMP);

    for ai = 1:n_AMP
        data = FR_heatmap_norm_sim(:, ai, si);
        valid = data(data > 0);
        if ~isempty(valid)
            mean_FR(ai) = mean(valid);
            sem_FR(ai) = std(valid) / sqrt(length(valid));
        end
    end

    set_name = sprintf('set%d', si);
    FR_summary_sim.(set_name).mean = mean_FR;
    FR_summary_sim.(set_name).sem = sem_FR;

    figure('Color','w', 'Name', sprintf('Simultaneous Stim Set %d', si));
    errorbar(Amps, mean_FR, sem_FR, '-s', 'LineWidth', 2);
    title(sprintf('Simultaneous Stimulation Set %d', si));
    xlabel('Amplitude (µA)');
    ylabel('Normalized Firing Rate');
    legend('Simultaneous', 'Location', 'northwest');
    ylim([0, 1.1 * max(mean_FR + sem_FR)]);
    box off;
end


%% === SEQUENTIAL STIMULATION ===
[~, ~, nSets_seq] = size(FR_heatmap_norm_seq);
FR_summary_seq = struct();

for si = 1:nSets_seq
    mean_FR = zeros(1, n_AMP);
    sem_FR  = zeros(1, n_AMP);

    for ai = 1:n_AMP
        data = FR_heatmap_norm_seq(:, ai, si);
        valid = data(data > 0);
        if ~isempty(valid)
            mean_FR(ai) = mean(valid);
            sem_FR(ai) = std(valid) / sqrt(length(valid));
        end
    end

    set_name = sprintf('set%d', si);
    FR_summary_seq.(set_name).mean = mean_FR;
    FR_summary_seq.(set_name).sem = sem_FR;

    figure('Color','w', 'Name', sprintf('Sequential Stim Set %d', si));
    errorbar(Amps, mean_FR, sem_FR, '-^', 'LineWidth', 2);
    title(sprintf('Sequential Stimulation Set %d', si));
    xlabel('Amplitude (µA)');
    ylabel('Normalized Firing Rate');
    legend('Sequential', 'Location', 'northwest');
    ylim([0, 1.1 * max(mean_FR + sem_FR)]);
    box off;
end