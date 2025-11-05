clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder

% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1'; 
data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim2';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Seq5_New_251014_194221';

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

% Extract file name
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);  % 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback if no underscores
end

%% Choice
Spike_filtering =0;


%% Pre Set
FS=30000; % Sampling frequency
% Load .sp.mat file
% sp_files = dir('*.sp.mat');
sp_files = dir(fullfile(data_folder, '*.sp.mat'));
assert(~isempty(sp_files), 'No .sp.mat file found in the current folder.');
% sp_filename = sp_files(1).name;
sp_filename = fullfile(data_folder, sp_files(1).name);
fprintf('Loading spike file: %s\n', sp_filename);
S = load(sp_filename);
if isfield(S, 'sp')
    sp = S.sp;
else
    error('Variable "sp" not found in %s.', sp_filename);
end

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0); 

%% Spike Amplitude Filtering Parameters
pos_limit = 100;    % upper bound (µV)
neg_limit = -100;  % lower bound (µV)

baseline_window_ms = [-60, -5];        % Baseline window (ms)
response_window_ms = [2, 10];          % Response window (ms)

%% Load StimParams and decode amplitudes, stimulation sets, ISI
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
% S = load(fileDIR(1).name,'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
S = load(fullfile(data_folder, fileDIR(1).name), 'StimParams', 'simultaneous_stim', 'CHN', 'E_MAP', 'n_Trials');
StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
CHN                = S.CHN;
E_MAP              = S.E_MAP;
n_Trials           = S.n_Trials;

% Trial amplitude list
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps, ~, ampIdx] = unique(trialAmps(:));
if any(Amps == -1), Amps(Amps == -1) = 0; end
n_AMP = numel(Amps);
cmap = lines(n_AMP);  % color map for amplitudes

% Stimulation set decoding
E_NAME = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_NAME);
stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = unique(idx_all(rr)); v = v(v > 0).';
    stimChPerTrial_all{t} = v;
end
comb = zeros(n_Trials, simultaneous_stim);
for i = 1:n_Trials
    v = stimChPerTrial_all{i};  
    comb(i,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb, 'rows');
combClass_win = combClass;
nSets = size(uniqueComb,1);

% Pulse Train Period (inter-pulse interval)
pulseTrain_all = cell2mat(StimParams(2:end,9));  % Column 9: Pulse Train Period
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);  % take 1 per trial

[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

% Electrode Map
d = Depth_s(1); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

%% Spike Amplitude Filtering (Before Plotting)
% % sp_clipped = sp;   % copy original spike structure

if Spike_filtering == 1

    fprintf('\n===== Spike Amplitude Filtering =====\n');
    fprintf('Criteria: max ≤ %.1f µV and min ≥ %.1f µV\n', pos_limit, neg_limit);
    
    for ch = 1:numel(sp)
        if isempty(sp{ch}), continue; end
    
        waveforms = sp{ch}(:, 2:end);
    
        % Amplitude-based filtering
        max_vals = max(waveforms, [], 2);
        min_vals = min(waveforms, [], 2);
        in_range = (max_vals <= pos_limit) & (min_vals >= neg_limit);
    
        % Zero-crossing check (must include both positive and negative values)
        has_positive = any(waveforms > 0, 2);
        has_negative = any(waveforms < 0, 2);
        crosses_zero = has_positive & has_negative;

        % Keep spikes within amplitude bounds
        % valid_idx = (max_vals <= pos_limit) & (min_vals >= neg_limit);
        % Final combined condition
        valid_idx = in_range & crosses_zero;

        % Apply filter
        sp_clipped{ch} = sp{ch}(valid_idx, :);
    
        % Print summary
        n_total = size(sp{ch}, 1);
        n_keep  = sum(valid_idx);
        fprintf('Channel %2d: kept %4d / %4d spikes (%.1f%%)\n', ...
            ch, n_keep, n_total, 100*n_keep/n_total);
    end
    
    fprintf('=====================================\n\n');
    save([base_name '.sp_xia.mat'],'sp_clipped');
else
    load([base_name '.sp_xia.mat']);

    % load([base_name '.sp.mat']);
    % sp_clipped = sp;

    % load([base_name '.sp_xia_FirstPulse.mat']);
    % sp_clipped = sp_seq;
end

%% === Baseline-Corrected Firing Rate per Channel === %%
fprintf('\n• Computing Baseline-Corrected FR per Channel \n');

nChn = numel(sp_clipped);
nTrials = length(trialAmps);
FR_baseline = nan(nChn, nTrials);
FR_response = nan(nChn, nTrials);
FR_corrected = nan(nChn, nTrials);

win_baseline = baseline_window_ms;  % e.g. [-60, -5]
win_response = response_window_ms;  % e.g. [2, 20]

baseline_s = diff(win_baseline) / 1000;
response_s = diff(win_response) / 1000;

for ch = 1:nChn
    spk = sp_clipped{d(ch)};
    if isempty(spk), continue; end

    for tr = 1:nTrials
        t0 = trig(tr) / FS * 1000;
        rel_spk = spk(:,1) - t0;

        FR_baseline(ch,tr) = sum(rel_spk >= win_baseline(1) & rel_spk < win_baseline(2)) / baseline_s;
        FR_response(ch,tr) = sum(rel_spk >= win_response(1) & rel_spk < win_response(2)) / response_s;
        FR_corrected(ch,tr) = FR_response(ch,tr) - FR_baseline(ch,tr);
        % FR_corrected(ch,tr) = FR_response(ch,tr);
    end
end

%% === Average FR per Channel per Amp per Stim Set ===
FR_heatmap = nan(nChn, n_AMP, nSets);
% FR_heatmap = nan(nChn, 6, nSets);
% FR_heatmap = nan(nChn, 5, nSets);
for si = 1:nSets
    trial_set = find(combClass_win == si);
    for ai = 1:n_AMP
        trial_amp = find(trialAmps == Amps(ai));
        trials_here = intersect(trial_set, trial_amp);
        if isempty(trials_here), continue; end

        FR_heatmap(:, ai, si) = mean(FR_corrected(:, trials_here), 2);
    end
end

%% === Normalize per-channel within each stim set ===
% max_FR_single = nan(nSets, nChn);
% load('max_FR_per_channel_singleStim.mat');  % loads max_FR_single
FR_heatmap_norm = nan(size(FR_heatmap));
for si = 1:nSets
    for ch = 1:nChn
        max_val = max(FR_heatmap(ch,:,si), [], 'omitnan');
        % max_FR_single(si,ch) = max_val;
        if max_val > 0            
            % FR_heatmap_norm(ch,:,si) = FR_heatmap(ch,:,si) / max_val;
                % FR_heatmap_norm(ch,:,si) = FR_heatmap(ch,:,si) / max_FR_single(1,ch);
                FR_heatmap_norm(ch,:,si) = FR_heatmap(ch,:,si);            
        else
            FR_heatmap_norm(ch,:,si) = 0;
        end
    end
end
% save('max_FR_per_channel_singleStim.mat', 'max_FR_single');
%% === Plot Heatmaps: One per Stimulation Set ===
for si = 1:nSets
    figure('Color','w', 'Name', sprintf('Stim Set %d', si));
    imagesc(FR_heatmap_norm(:,:,si));
    colormap(parula);
    colorbar;

    stimChs = uniqueComb(si, uniqueComb(si,:) > 0);
    stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimChs, 'UniformOutput', false), '+');

    title(sprintf('%s (Simultaneous)', stimLabel), 'Interpreter','none');
    xlabel('Amplitude (µA)');
    xticks(1:n_AMP);
    xticklabels(arrayfun(@(x) num2str(x), Amps, 'UniformOutput', false));
    ylabel('Channel');
    % caxis([0 1]);  % normalized
end



%% === Plot Average FR per Stim Set with Error Bars ===
fprintf('\n• Plotting Average Firing Rate Curves per Stim Set \n');

mean_FR_set = cell(1, nSets);
sem_FR_set  = cell(1, nSets);

for si = 1:nSets
    data = FR_heatmap_norm(:,:,si);  % [channels × amplitudes]
    
    % Exclude 0 and negative FR values
    % data(data <= 0) = NaN;
    
    % Compute average and SEM across channels (omit NaNs)
    mean_FR = mean(data, 1, 'omitnan');
    sem_FR  = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));

    % Save to variables for future access
    mean_FR_set{si} = mean_FR;
    sem_FR_set{si}  = sem_FR;

    % Plot with error bars
    figure('Color','w', 'Name', sprintf('Avg FR – Stim Set %d', si));
    errorbar(Amps, mean_FR, sem_FR, '-o', 'LineWidth', 1.5);
    hold on;

    stimChs = uniqueComb(si, uniqueComb(si,:) > 0);
    stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimChs, 'UniformOutput', false), '+');

    title(sprintf('Avg Firing Rate – %s (Simultaneous)', stimLabel), 'Interpreter','none');
    xlabel('Amplitude (µA)');
    ylabel('Mean Firing Rate (spikes/s)');
    ylim([0, 1.1 * max(mean_FR + sem_FR, [], 'omitnan')]);
    legend(sprintf('Stim Set %d', si), 'Location', 'northwest');
    box off;
end

%% === Save Firing Rate per Channel for Each Stim Set ===
fprintf('\n• Saving Firing Rate per Channel for Each Stim Set\n');

FR_data = struct();

for si = 1:nSets
    stimChs = uniqueComb(si, uniqueComb(si,:) > 0);
    stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimChs, 'UniformOutput', false), '+');

    FR_data(si).stimSet      = si;
    FR_data(si).stimChannels = stimChs;
    FR_data(si).stimLabel    = stimLabel;
    FR_data(si).Amps         = Amps;
    FR_data(si).FR_perCh     = FR_heatmap_norm(:,:,si);  % [channels × amplitudes]
    FR_data(si).meanFR       = mean_FR_set{si};          % [1 × amplitudes]
    FR_data(si).semFR        = sem_FR_set{si};           % [1 × amplitudes]
end

output_filename = sprintf('%s_FR_perChn_perSet.mat', base_name);
save(output_filename, 'FR_data', 'Amps', 'uniqueComb', 'FR_heatmap_norm');
fprintf('• Saved firing rate summary to: %s\n', output_filename);