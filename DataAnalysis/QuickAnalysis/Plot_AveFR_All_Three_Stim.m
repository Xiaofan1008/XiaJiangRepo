%% === Plot FR Summary Across Conditions: Simultaneous vs Sequential vs Separate ===

data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Stim_Set2_FR_correct'; 

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