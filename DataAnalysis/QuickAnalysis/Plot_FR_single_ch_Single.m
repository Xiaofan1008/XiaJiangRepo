clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder

data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single2'; 
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim2';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq2';

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
ch_to_plot = 18; 
stimSet_to_plot = 1;
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
fprintf('\n• Computing Baseline-Corrected FR per Channel\n');

nChn = numel(sp_clipped);
nTrials = length(trialAmps);
FR_baseline = nan(nChn, nTrials);
FR_response = nan(nChn, nTrials);
FR_corrected = nan(nChn, nTrials);

win_baseline = baseline_window_ms;  
win_response = response_window_ms;  

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

%% === Plot Firing Rate for a Specific Channel === %%  
fprintf('\nPlotting firing rate for Channel %d\n', ch_to_plot);

stim_mask = (combClass_win == stimSet_to_plot);

% Extract amplitudes and FRs for the chosen set
amps_this = trialAmps(stim_mask);
FR_this   = FR_corrected(ch_to_plot, stim_mask);

% Group by amplitude
[amps_unique, ~, amp_idx] = unique(amps_this);
mean_FR = nan(size(amps_unique));
sem_FR  = nan(size(amps_unique));

for i = 1:numel(amps_unique)
    vals = FR_this(amp_idx == i);
    vals = vals(~isnan(vals));
    mean_FR(i) = mean(vals);
    sem_FR(i)  = std(vals) / sqrt(numel(vals));
end

% Plot with error bars
figure; hold on;
errorbar(amps_unique, mean_FR, sem_FR, 'o-', 'LineWidth', 2, ...
    'MarkerSize', 6, 'Color', cmap(mod(stimSet_to_plot-1, size(cmap,1))+1, :));
xlabel('Stimulation Amplitude (\muA)');
ylabel('Baseline-corrected Firing Rate (spikes/s)');
box off

% Mark titles/legend
stim_electrodes = uniqueComb(stimSet_to_plot, :);
stim_electrodes = stim_electrodes(stim_electrodes > 0);
stim_str = sprintf('Set %d: [Ch %s]', stimSet_to_plot, num2str(stim_electrodes));

title(sprintf('Channel %d | %s', ch_to_plot, stim_str), 'FontWeight', 'bold');
legend({sprintf('Stim Set %d', stimSet_to_plot)}, 'Location', 'best');