load('FR_heatmap_norm_sim.mat')
load('FR_heatmap_norm_seq.mat')
load('FR_heatmap_norm_single.mat')

Amps = [0,1,2,3,4,5];
n_AMP = length(Amps);

%% === STEP 1: Compute mean FR for each stim set ===

% --- Single ---
[~, ~, nSets_single] = size(FR_heatmap_norm_single);
mean_FR_single_all = zeros(nSets_single, n_AMP);
sem_FR_single_all  = zeros(nSets_single, n_AMP);
for si = 1:nSets_single
    for ai = 1:n_AMP
        data = FR_heatmap_norm_single(:, ai, si);
        mean_FR_single_all(si, ai) = mean(data, 'omitnan');
        sem_FR_single_all(si, ai)  = std(data, 'omitnan') / sqrt(sum(~isnan(data)));
    end
end

% --- Simultaneous ---
[~, ~, nSets_sim] = size(FR_heatmap_norm_sim);
mean_FR_sim_all = zeros(nSets_sim, n_AMP);
sem_FR_sim_all  = zeros(nSets_sim, n_AMP);
for si = 1:nSets_sim
    for ai = 1:n_AMP
        data = FR_heatmap_norm_sim(:, ai, si);
        mean_FR_sim_all(si, ai) = mean(data, 'omitnan');
        sem_FR_sim_all(si, ai)  = std(data, 'omitnan') / sqrt(sum(~isnan(data)));
    end
end

% --- Sequential ---
[~, ~, nSets_seq] = size(FR_heatmap_norm_seq);
mean_FR_seq_all = zeros(nSets_seq, n_AMP);
sem_FR_seq_all  = zeros(nSets_seq, n_AMP);
for si = 1:nSets_seq
    for ai = 1:n_AMP
        data = FR_heatmap_norm_seq(:, ai, si);
        mean_FR_seq_all(si, ai) = mean(data, 'omitnan');
        sem_FR_seq_all(si, ai)  = std(data, 'omitnan') / sqrt(sum(~isnan(data)));
    end
end


%% === STEP 2: Combine sets ===

% --- Single: (Set1+2) and (Set3+4), then average the two curves
avg_12_mean = (mean_FR_single_all(1,:) + mean_FR_single_all(2,:)) / 2;
avg_34_mean = (mean_FR_single_all(3,:) + mean_FR_single_all(4,:)) / 2;
mean_FR_single_final = (avg_12_mean + avg_34_mean) / 2;

avg_12_sem = sqrt((sem_FR_single_all(1,:).^2 + sem_FR_single_all(2,:).^2) / 4);
avg_34_sem = sqrt((sem_FR_single_all(3,:).^2 + sem_FR_single_all(4,:).^2) / 4);
sem_FR_single_final = sqrt((avg_12_sem.^2 + avg_34_sem.^2) / 4);

% --- Simultaneous: average across sets
mean_FR_sim_final = mean(mean_FR_sim_all, 1, 'omitnan');
sem_FR_sim_final  = sqrt(mean(sem_FR_sim_all.^2, 1, 'omitnan'));

% --- Sequential: average across sets
mean_FR_seq_final = mean(mean_FR_seq_all, 1, 'omitnan');
sem_FR_seq_final  = sqrt(mean(sem_FR_seq_all.^2, 1, 'omitnan'));


%% === STEP 3: Plot ===
figure('Color','w', 'Name', 'Average FR Across Stimulation Types');
hold on;

errorbar(Amps, mean_FR_single_final, sem_FR_single_final, '-o', ...
    'LineWidth', 2, 'DisplayName', 'Single (avg of 1&2 + 3&4)');
errorbar(Amps, mean_FR_sim_final, sem_FR_sim_final, '-s', ...
    'LineWidth', 2, 'DisplayName', 'Simultaneous (avg across sets)');
errorbar(Amps, mean_FR_seq_final, sem_FR_seq_final, '-^', ...
    'LineWidth', 2, 'DisplayName', 'Sequential (avg across sets)');

xlabel('Amplitude (µA)');
ylabel('Normalized Firing Rate');
title('Average Normalized Firing Rate (Across Sets)');
legend('Location','northwest');
ylim([0, 1.1 * max([mean_FR_single_final + sem_FR_single_final, ...
                   mean_FR_sim_final + sem_FR_sim_final, ...
                   mean_FR_seq_final + sem_FR_seq_final])]);
box off;