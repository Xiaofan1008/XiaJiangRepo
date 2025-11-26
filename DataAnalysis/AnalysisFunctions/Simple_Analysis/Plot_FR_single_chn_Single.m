clear;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single2';  % <--- CHANGE THIS
cd(data_folder);

channel_to_plot = 30;

% Extract file name
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);  % 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback if no underscores
end


% Load spike data
% S = load(fullfile(data_folder, '*.sp_xia.mat'));
S = load([base_name '.sp_xia.mat']);
sp_clipped = S.sp_clipped;

% Load triggers
trig = loadTrig(0); 

% Load experiment info
expfile = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
Sexp = load(fullfile(data_folder, expfile(1).name), 'StimParams', 'E_MAP', 'simultaneous_stim', 'n_Trials');
StimParams = Sexp.StimParams;
E_MAP = Sexp.E_MAP;
simultaneous_stim = Sexp.simultaneous_stim;
n_Trials = Sexp.n_Trials;
FS = 30000;

%% === Decode Stim Params ===
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps, ~, ampIdx] = unique(trialAmps(:));
n_AMP = numel(Amps);

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
nSets = size(uniqueComb,1);

% Electrode Map
d = Depth_s(1); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

%% === Firing Rate per Trial ===
baseline_window_ms = [-60, -5];
response_window_ms = [2, 10];
baseline_s = diff(baseline_window_ms) / 1000;
response_s = diff(response_window_ms) / 1000;

nChn = numel(sp_clipped);
FR_Trials = cell(nSets, 1);

for si = 1:nSets
    stimChs = uniqueComb(si, uniqueComb(si,:) > 0);
    stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimChs, 'UniformOutput', false), '+');
    trial_idxs = find(combClass == si);

    thisFR = nan(nChn, n_AMP, numel(trial_idxs));  % [chn × amp × trials]

    for t_idx = 1:numel(trial_idxs)
        tr = trial_idxs(t_idx);
        amp = trialAmps(tr);
        amp_idx = find(Amps == amp);

        t0 = trig(tr) / FS * 1000;

        for ch = 1:nChn
            spk = sp_clipped{d(ch)};
            if isempty(spk), continue; end
            rel_spk = spk(:,1) - t0;

            n_baseline = sum(rel_spk >= baseline_window_ms(1) & rel_spk < baseline_window_ms(2));
            n_response = sum(rel_spk >= response_window_ms(1) & rel_spk < response_window_ms(2));
            fr = (n_response / response_s) - (n_baseline / baseline_s);

            thisFR(ch, amp_idx, t_idx) = fr;
        end
    end

    FR_Trials{si}.stimSet = si;
    FR_Trials{si}.stimChannels = stimChs;
    FR_Trials{si}.stimLabel = stimLabel;
    FR_Trials{si}.FR_trials = thisFR;  % [chn × amp × trials]
end

%% === Save ===
fprintf('Saving FR_Trials structure...\n');
save('FR_perTrial.mat', 'FR_Trials', 'Amps', 'uniqueComb');
fprintf('Done. File saved as: FR_perTrial.mat\n');


%% === Plot Firing Rate Curve for a Specific Channel ===

fprintf('\nPlotting firing rate for Channel %d...\n', channel_to_plot);

colors = lines(nSets);
figure('Color','w', 'Name', sprintf('FR per Trial – Channel %d', channel_to_plot)); hold on;

for si = 1:nSets
    FR = FR_Trials{si}.FR_trials;  % [chn × amp × trial]
    stimLabel = FR_Trials{si}.stimLabel;

    if channel_to_plot > size(FR,1)
        warning('Channel %d not found in stim set %d, skipping.', channel_to_plot, si);
        continue;
    end

    fr_ch = squeeze(FR(channel_to_plot,:,:));  % [amp × trials]
    mean_fr = mean(fr_ch, 2, 'omitnan');
    sem_fr  = std(fr_ch, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(fr_ch), 2));

    errorbar(Amps, mean_fr, sem_fr, '-o', 'LineWidth', 1.5, 'Color', colors(si,:), ...
        'DisplayName', sprintf('Set %d: %s', si, stimLabel));
end

xlabel('Amplitude (µA)');
ylabel('Firing Rate (spikes/s)');
title(sprintf('Firing Rate vs Amplitude – Channel %d', channel_to_plot));
legend('show', 'Location', 'northwest');
box off;