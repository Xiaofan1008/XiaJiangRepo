%% === Plot FR Summary Across Conditions: Simultaneous vs Sequential vs Separate ===

data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Stim_Set1_FR_correct'; 

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Using data from folder:\n%s\n\n', data_folder);

% Load all *_FR_perChn_perSet.mat files
files = dir(fullfile(data_folder, '*_FR_perChn_perSet.mat'));
if numel(files) < 3
    error('Need at least 3 *_FR_perChn_perSet.mat files: simultaneous, sequential, and separate.');
end

% Sort and load files
for i = 1:length(files)
    fname = files(i).name;
    if contains(lower(fname), 'sim') || contains(lower(fname), 'simultaneous')
        S_sim = load(fullfile(data_folder, fname));
    elseif contains(lower(fname), 'seq') || contains(lower(fname), 'sequential')
        S_seq = load(fullfile(data_folder, fname));
    elseif contains(lower(fname), 'separate') || contains(lower(fname), 'single')
        S_sep = load(fullfile(data_folder, fname));
    else
        warning('Unknown file type: %s', fname);
    end
end

% Validate
if ~exist('S_sim', 'var') || ~exist('S_seq', 'var') || ~exist('S_sep', 'var')
    error('Missing one of the required files (sim, seq, or separate).');
end

% Extract from loaded data
Amps = S_sim.Amps;  % Assume all have same amplitudes
nAmp = numel(Amps);
nChn = size(S_sim.FR_data(1).FR_perCh, 1);

% === Sum of separate stim sets ===
sep_FR1 = S_sep.FR_data(1).FR_perCh;  % [ch × amp]
sep_FR2 = S_sep.FR_data(2).FR_perCh;
FR_sep_sum = sep_FR1 + sep_FR2;

% === Get average and SEM across channels ===
mean_FR_sim = mean(S_sim.FR_data(1).FR_perCh, 1, 'omitnan');
sem_FR_sim  = std(S_sim.FR_data(1).FR_perCh, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(S_sim.FR_data(1).FR_perCh)));

mean_FR_seq = mean(S_seq.FR_data(1).FR_perCh, 1, 'omitnan');
sem_FR_seq  = std(S_seq.FR_data(1).FR_perCh, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(S_seq.FR_data(1).FR_perCh)));

mean_FR_sep = mean(FR_sep_sum, 1, 'omitnan');
sem_FR_sep  = std(FR_sep_sum, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(FR_sep_sum)));



%% === Plot Combined Firing Rate Curves ===
figure('Color','w', 'Name', 'FR Curve Comparison');

% Get stim channel labels
label_sim = strjoin(arrayfun(@(x) sprintf('Ch%d', x), S_sim.FR_data(1).stimChannels, 'UniformOutput', false), '+');
label_seq = strjoin(arrayfun(@(x) sprintf('Ch%d', x), S_seq.FR_data(1).stimChannels, 'UniformOutput', false), '+');
label_sep1 = strjoin(arrayfun(@(x) sprintf('Ch%d', x), S_sep.FR_data(1).stimChannels, 'UniformOutput', false), '+');
label_sep2 = strjoin(arrayfun(@(x) sprintf('Ch%d', x), S_sep.FR_data(2).stimChannels, 'UniformOutput', false), '+');

errorbar(Amps, mean_FR_sim, sem_FR_sim, '-o', 'LineWidth', 2); hold on;
errorbar(Amps, mean_FR_seq, sem_FR_seq, '-s', 'LineWidth', 2);
errorbar(Amps, mean_FR_sep(2:6), sem_FR_sep(2:6), '-^', 'LineWidth', 2);

legend({ ...
    sprintf('Simultaneous (%s)', label_sim), ...
    sprintf('Sequential (%s)', label_seq), ...
    sprintf('Separate Sum (%s)', label_seq)}, ...
    'Location','northwest');

xlabel('Amplitude (µA)');
ylabel('Average Firing Rate (spikes/s)');
title('Comparison of Firing Rates Across Conditions');
ylim([0, 1.1 * max([mean_FR_sim, mean_FR_seq, mean_FR_sep(2:6)] + ...
    [sem_FR_sim, sem_FR_seq, sem_FR_sep(2:6)], [], 'omitnan')]);
box off;

%% === NEW SECTION: Plot Firing Rate Curve for a Specific Channel ===
% ---------------------------------------------------------------
ch_target = 30;  % desired channel number
fprintf('\nPlotting FR curve for Channel %d...\n', ch_target);

% Extract FR for this channel from each condition
FR_sim_ch = S_sim.FR_data(1).FR_perCh(ch_target, :);
FR_seq_ch = S_seq.FR_data(1).FR_perCh(ch_target, :);
FR_sep1_ch = S_sep.FR_data(1).FR_perCh(ch_target, :);
FR_sep2_ch = S_sep.FR_data(2).FR_perCh(ch_target, :);

% Plot the channel-specific FR curve
figure('Color','w', 'Name', sprintf('FR Curve – Ch%d', ch_target));
hold on;

p1 = plot(Amps, FR_sim_ch, '-o', 'LineWidth', 2);
p2 = plot(Amps, FR_seq_ch, '-s', 'LineWidth', 2);
p3 = plot(Amps, FR_sep1_ch(2:6), '-^', 'LineWidth', 2);
p4 = plot(Amps, FR_sep2_ch(2:6), '-v', 'LineWidth', 2);

xlabel('Amplitude (µA)');
ylabel(sprintf('Firing Rate (spikes/s) – Ch%d', ch_target));
title(sprintf('Firing Rate Curve – Ch%d (Per Stim Type)', ch_target));
legend([p1, p2, p3, p4], { ...
    sprintf('Simultaneous (%s)', label_sim), ...
    sprintf('Sequential (%s)', label_seq), ...
    sprintf('Separate (%s)', label_sep1), ...
    sprintf('Separate (%s)', label_sep2)}, ...
    'Location', 'northwest');
box off;
ylim([0, 1.1 * max([FR_sim_ch, FR_seq_ch, FR_sep1_ch(2:6), FR_sep2_ch(2:6)], [], 'omitnan')]);
fprintf('✓ Firing rate curve plotted for Channel %d.\n', ch_target);