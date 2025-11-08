% ============================================================
% Compare and Plot Firing Rate Curves Across Stimulation Types
% ============================================================
% This script:
% 1. Loads firing rate .mat files for Single, Simultaneous, and Sequential stimulations
% 2. Plots mean ± SEM firing rate curves for a chosen channel
% 3. Optionally plots the summed firing rate (two single sets),
%    aligning amplitudes before summation
% ============================================================

clear all; close all; clc;

%% === User Input === %%
% Folder paths
folder_single = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single3';
folder_simul  = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim3';
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq3';

% Channel to plot
ch_to_plot = 18;

% Whether to plot the summed single-stim sets
plot_sum = true;      % set false to disable
sum_sets = [1, 2];    % two sets from the single stimulation to sum

%% === Load Firing Rate Data === %%
fprintf('Loading firing rate files...\n');

file_single = dir(fullfile(folder_single, '*FR_AllCh_AllSets.mat'));
file_simul  = dir(fullfile(folder_simul,  '*FR_AllCh_AllSets.mat'));
file_seq    = dir(fullfile(folder_seq,    '*FR_AllCh_AllSets.mat'));

assert(~isempty(file_single), 'Single-stim file not found!');
assert(~isempty(file_simul),  'Simultaneous-stim file not found!');
assert(~isempty(file_seq),    'Sequential-stim file not found!');

S_single = load(fullfile(folder_single, file_single(1).name));
S_simul  = load(fullfile(folder_simul,  file_simul(1).name));
S_seq    = load(fullfile(folder_seq,    file_seq(1).name));

FR_single = S_single.FR_summary;
FR_simul  = S_simul.FR_summary;
FR_seq    = S_seq.FR_summary;

fprintf('Files loaded successfully.\n');

%% === Plot Setup === %%
figure; hold on; box off;
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
cmap = lines(6);
lw = 2;

fprintf('Plotting firing rate for channel %d...\n', ch_to_plot);

%% === 1. Plot Single Stimulation === %%
for s = 1:numel(FR_single)
    amps = FR_single(s).amps(:);
    mean_FR = FR_single(s).mean(ch_to_plot,:);
    sem_FR  = FR_single(s).sem(ch_to_plot,:);
    stimCh  = FR_single(s).stimCh;
    if all(isnan(mean_FR)), continue; end

    label_text = sprintf('Single Set %d [Ch %s]', s, num2str(stimCh));
    errorbar(amps, mean_FR, sem_FR, 'o-', 'LineWidth', lw, ...
        'MarkerSize', 6, 'Color', cmap(1,:), 'DisplayName', label_text);
end

%% === 2. Plot Simultaneous Stimulation === %%
for s = 1:numel(FR_simul)
    amps = FR_simul(s).amps(:);
    mean_FR = FR_simul(s).mean(ch_to_plot,:);
    sem_FR  = FR_simul(s).sem(ch_to_plot,:);
    stimCh  = FR_simul(s).stimCh;
    if all(isnan(mean_FR)), continue; end

    label_text = sprintf('Simultaneous Set %d [Ch %s]', s, num2str(stimCh));
    errorbar(amps, mean_FR, sem_FR, '^-', 'LineWidth', lw, ...
        'MarkerSize', 6, 'Color', cmap(2,:), 'DisplayName', label_text);
end

%% === 3. Plot Sequential Stimulation === %%
for s = 1:numel(FR_seq)
    amps = FR_seq(s).amps(:);
    mean_FR = FR_seq(s).mean(ch_to_plot,:);
    sem_FR  = FR_seq(s).sem(ch_to_plot,:);
    stimCh  = FR_seq(s).stimCh;
    if all(isnan(mean_FR)), continue; end

    label_text = sprintf('Sequential Set %d [Ch %s]', s, num2str(stimCh));
    errorbar(amps, mean_FR, sem_FR, 's-', 'LineWidth', lw, ...
        'MarkerSize', 6, 'Color', cmap(3,:), 'DisplayName', label_text);
end

%% === 4. Optional: Sum of Two Single-Stim Sets (Dashed Line) === %%
if plot_sum
    fprintf('Computing summed firing rate from single sets %s...\n', num2str(sum_sets));

    % Extract data for both sets
    data1 = FR_single(sum_sets(1));
    data2 = FR_single(sum_sets(2));

    amps1 = data1.amps(:);
    amps2 = data2.amps(:);

    % --- Align amplitudes between the two sets ---
    [amps_common, idx1, idx2] = intersect(amps1, amps2);

    if isempty(amps_common)
        warning('No matching amplitudes found between Set %d and Set %d!', sum_sets(1), sum_sets(2));
    else
        % Extract matching mean ± SEM
        mean1 = data1.mean(ch_to_plot, idx1);
        sem1  = data1.sem(ch_to_plot, idx1);
        mean2 = data2.mean(ch_to_plot, idx2);
        sem2  = data2.sem(ch_to_plot, idx2);

        % --- Compute summed FR and propagated SEM ---
        mean_sum = mean1 + mean2;
        sem_sum  = sqrt(sem1.^2 + sem2.^2);

        % Stim channel info
        stimCh1 = num2str(data1.stimCh);
        stimCh2 = num2str(data2.stimCh);
        label_sum = sprintf('Sum Single Sets %d+%d [Ch %s + %s]', ...
            sum_sets(1), sum_sets(2), stimCh1, stimCh2);

        % --- Plot dashed line ---
        errorbar(amps_common, mean_sum, sem_sum, 's--', 'LineWidth', lw, ...
            'MarkerSize', 6, 'Color', [0.3 0.3 0.3], 'DisplayName', label_sum);
    end
end

%% === Finalize Plot === %%
xlabel('Stimulation Amplitude (\muA)', 'FontWeight', 'bold');
ylabel('Baseline-corrected Firing Rate (spikes/s)', 'FontWeight', 'bold');
title(sprintf('Channel %d | Firing Rate Comparison Across Stimulation Types', ch_to_plot), 'FontWeight', 'bold');
legend('show', 'Location', 'best');
grid on;

fprintf('Done.\n');