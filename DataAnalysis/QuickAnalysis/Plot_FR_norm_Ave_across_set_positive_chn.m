load('FR_heatmap_norm_sim.mat')
load('FR_heatmap_norm_seq.mat')
load('FR_heatmap_norm_single.mat')

Amps = [0,1,2,3,4,5];
n_AMP = length(Amps);

%% === STEP 1: Compute mean & SEM for each stim set ===

% --- Single ---
[~, ~, nSets_single] = size(FR_heatmap_norm_single);
mean_FR_single_all = zeros(nSets_single, n_AMP);
sem_FR_single_all  = zeros(nSets_single, n_AMP);

for si = 1:nSets_single
    for ai = 1:n_AMP
        data = FR_heatmap_norm_single(:, ai, si);
        data = data(data > 0);  % ignore zero & negative
        if isempty(data)
            mean_FR_single_all(si, ai) = NaN;
            sem_FR_single_all(si, ai)  = NaN;
        else
            mean_FR_single_all(si, ai) = mean(data);
            sem_FR_single_all(si, ai)  = std(data) / sqrt(length(data));
        end
    end
end

% --- Simultaneous ---
[~, ~, nSets_sim] = size(FR_heatmap_norm_sim);
mean_FR_sim_all = zeros(nSets_sim, n_AMP);
sem_FR_sim_all  = zeros(nSets_sim, n_AMP);

for si = 1:nSets_sim
    for ai = 1:n_AMP
        data = FR_heatmap_norm_sim(:, ai, si);
        data = data(data > 0);
        if isempty(data)
            mean_FR_sim_all(si, ai) = NaN;
            sem_FR_sim_all(si, ai)  = NaN;
        else
            mean_FR_sim_all(si, ai) = mean(data);
            sem_FR_sim_all(si, ai)  = std(data) / sqrt(length(data));
        end
    end
end

% --- Sequential ---
[~, ~, nSets_seq] = size(FR_heatmap_norm_seq);
mean_FR_seq_all = zeros(nSets_seq, n_AMP);
sem_FR_seq_all  = zeros(nSets_seq, n_AMP);

for si = 1:nSets_seq
    for ai = 1:n_AMP
        data = FR_heatmap_norm_seq(:, ai, si);
        data = data(data > 0);
        if isempty(data)
            mean_FR_seq_all(si, ai) = NaN;
            sem_FR_seq_all(si, ai)  = NaN;
        else
            mean_FR_seq_all(si, ai) = mean(data);
            sem_FR_seq_all(si, ai)  = std(data) / sqrt(length(data));
        end
    end
end


%% === STEP 2: Combine across sets ===

% --- Single: (Set1+2) and (Set3+4), then average the two curves
avg_12_mean = nanmean([mean_FR_single_all(1,:); mean_FR_single_all(2,:)], 1);
avg_34_mean = nanmean([mean_FR_single_all(3,:); mean_FR_single_all(4,:)], 1);
mean_FR_single_final = nanmean([avg_12_mean; avg_34_mean], 1);

avg_12_sem = sqrt(nanmean([sem_FR_single_all(1,:).^2; sem_FR_single_all(2,:).^2], 1));
avg_34_sem = sqrt(nanmean([sem_FR_single_all(3,:).^2; sem_FR_single_all(4,:).^2], 1));
sem_FR_single_final = sqrt(nanmean([avg_12_sem.^2; avg_34_sem.^2], 1));

% --- Simultaneous: average across all sets
mean_FR_sim_final = nanmean(mean_FR_sim_all, 1);
sem_FR_sim_final  = sqrt(nanmean(sem_FR_sim_all.^2, 1));

% --- Sequential: average across all sets
mean_FR_seq_final = nanmean(mean_FR_seq_all, 1);
sem_FR_seq_final  = sqrt(nanmean(sem_FR_seq_all.^2, 1));


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
title('Average Normalized Firing Rate (Above Baseline FR)');
legend('Location','northwest');
ylim([0, 1.1 * max([mean_FR_single_final + sem_FR_single_final, ...
                   mean_FR_sim_final + sem_FR_sim_final, ...
                   mean_FR_seq_final + sem_FR_seq_final], [], 'all', 'omitnan')]);
box off;