% ============================================================
% Compute mean ± SEM firing rate per channel and per stimulation set
%
% This code calculates baseline-corrected firing rates for every channel
% and stimulation set, performs a permutation test (p < 0.05) to detect
% significantly responsive channels, prints the results (channel numbers),
% and saves all outputs to a .mat file.
% ============================================================

clear all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% === Folder Selection ===
data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq1';  % change as needed
assert(isfolder(data_folder), 'Folder not found!');
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

% Extract base name
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);
else
    base_name = last_folder;
end
isSequential = contains(lower(data_folder), 'seq'); % detect sequential mode

%% === Parameters ===
Spike_filtering = 0;
FS = 30000; % Sampling rate
baseline_window_ms = [-60, -5];
response_window_ms = [2, 15];
pos_limit = 100;
neg_limit = -100;
nPerm = 5000;  % number of permutations for the permutation test

%% === Load Spike Data ===
sp_files = dir(fullfile(data_folder, '*.sp.mat'));
assert(~isempty(sp_files), 'No .sp.mat file found.');
sp_filename = fullfile(data_folder, sp_files(1).name);
fprintf('Loading spike file: %s\n', sp_filename);
S = load(sp_filename);
sp = S.sp;

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);

%% === Load StimParams and Decode Sets ===
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
S = load(fullfile(data_folder, fileDIR(1).name), ...
    'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
StimParams = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP = S.E_MAP;
n_Trials = S.n_Trials;

% Trial amplitudes
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps,~,ampIdx] = unique(trialAmps(:));
if any(Amps==-1), Amps(Amps==-1)=0; end
n_AMP = numel(Amps);

% Stimulation set decoding
E_NAME = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~,idx_all] = ismember(stimNames,E_NAME);
stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = unique(idx_all(rr)); v = v(v>0).';
    stimChPerTrial_all{t} = v;
end
comb = zeros(n_Trials, simultaneous_stim);
for i = 1:n_Trials
    v = stimChPerTrial_all{i};
    comb(i,1:numel(v)) = v;
end
[uniqueComb,~,combClass] = unique(comb,'rows');
combClass_win = combClass;
nSets = size(uniqueComb,1);

% Pulse train period
pulseTrain_all = cell2mat(StimParams(2:end,9));
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);
[PulsePeriods,~,pulseIdx] = unique(pulseTrain(:));

d = Depth_s(1);

%% === Load Filtered Spikes ===
if Spike_filtering == 1
    fprintf('\nApplying spike amplitude filtering...\n');
    for ch = 1:numel(sp)
        if isempty(sp{ch}), continue; end
        waveforms = sp{ch}(:,2:end);
        max_vals = max(waveforms,[],2);
        min_vals = min(waveforms,[],2);
        in_range = (max_vals <= pos_limit) & (min_vals >= neg_limit);
        has_positive = any(waveforms>0,2);
        has_negative = any(waveforms<0,2);
        valid_idx = in_range & has_positive & has_negative;
        sp_clipped{ch} = sp{ch}(valid_idx,:);
    end
    save([base_name '.sp_xia.mat'],'sp_clipped');
else
    if isSequential
        load([base_name '.sp_xia_FirstPulse.mat']);
        sp_clipped = sp_seq;
    else
        load([base_name '.sp_xia.mat']);
    end
end

%% === Compute Firing Rates ===
fprintf('\nComputing firing rates for all channels and stim sets...\n');
nChn = numel(sp_clipped);
nTrials = length(trialAmps);
FR_baseline = nan(nChn,nTrials);
FR_response = nan(nChn,nTrials);
FR_corrected = nan(nChn,nTrials);

baseline_s = diff(baseline_window_ms)/1000;
response_s = diff(response_window_ms)/1000;

for ch = 1:nChn
    spk = sp_clipped{d(ch)};
    if isempty(spk), continue; end
    for tr = 1:nTrials
        t0 = trig(tr)/FS*1000;
        rel_spk = spk(:,1)-t0;
        FR_baseline(ch,tr) = sum(rel_spk >= baseline_window_ms(1) & rel_spk < baseline_window_ms(2)) / baseline_s;
        FR_response(ch,tr) = sum(rel_spk >= response_window_ms(1) & rel_spk < response_window_ms(2)) / response_s;
        FR_corrected(ch,tr) = FR_response(ch,tr) - FR_baseline(ch,tr);
    end
end

%% === Permutation Test for Responsiveness ===
FR_summary = struct;
responsive_counts = zeros(1, nSets);
responsive_channels = cell(1, nSets);

fprintf('\nPerforming permutation test for responsiveness (nPerm = %d)...\n', nPerm);

for s = 1:nSets
    stim_mask = (combClass_win == s);
    amps_this = trialAmps(stim_mask);
    [amps_unique,~,amp_idx] = unique(amps_this);
    nAmp = numel(amps_unique);

    FR_summary(s).stimSet = s;
    FR_summary(s).amps = amps_unique;
    FR_summary(s).mean = nan(nChn,nAmp);
    FR_summary(s).sem  = nan(nChn,nAmp);
    FR_summary(s).raw  = cell(nChn,nAmp);
    FR_summary(s).stimCh = uniqueComb(s,uniqueComb(s,:)>0);
    FR_summary(s).pval  = nan(nChn,1);
    FR_summary(s).sig   = false(nChn,1);

    for ch = 1:nChn
        FR_base = FR_baseline(ch, stim_mask);
        FR_resp = FR_response(ch, stim_mask);
        if all(isnan(FR_base)) || all(isnan(FR_resp)), continue; end

        % --- Permutation test ---
        diff_obs = mean(FR_resp - FR_base, 'omitnan');
        valid_idx = ~isnan(FR_resp) & ~isnan(FR_base);
        FR_diff = FR_resp(valid_idx) - FR_base(valid_idx);

        if numel(FR_diff) < 5
            p = NaN;  % too few trials
        else
            nTrials_valid = numel(FR_diff);
            null_diffs = zeros(1,nPerm);
            for p_i = 1:nPerm
                signs = (rand(1,nTrials_valid) > 0.5)*2 - 1; % ±1 random flip
                null_diffs(p_i) = mean(signs .* FR_diff);
            end
            p = sum(null_diffs >= diff_obs) / nPerm; % one-tailed test
        end

        FR_summary(s).pval(ch) = p;
        FR_summary(s).sig(ch)  = (p < 0.05);

        % --- Mean ± SEM per amplitude ---
        for i = 1:nAmp
            vals = FR_corrected(ch, stim_mask);
            vals = vals(amp_idx == i);
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            FR_summary(s).mean(ch,i) = mean(vals);
            FR_summary(s).sem(ch,i)  = std(vals)/sqrt(numel(vals));
            FR_summary(s).raw{ch,i}  = vals;
        end
    end

    % Count responsive channels
    responsive_idx = find(FR_summary(s).sig);
    responsive_counts(s) = numel(responsive_idx);
    responsive_channels{s} = responsive_idx;

    if ~isempty(responsive_idx)
        fprintf('Stim Set %d [Ch %s]: %d responsive channels (p<0.05)\n', ...
            s, num2str(FR_summary(s).stimCh), responsive_counts(s));
        fprintf('→ Responding channels: %s\n', num2str(responsive_idx(:)'));
    else
        fprintf('Stim Set %d [Ch %s]: 0 responsive channels.\n', ...
            s, num2str(FR_summary(s).stimCh));
    end
end

%% === Save Results ===
% save_name = sprintf('%s_FR_AllCh_AllSets_permTest.mat', base_name);
% save(fullfile(data_folder, save_name), 'FR_summary', 'Amps', 'E_MAP', ...
%     'uniqueComb', 'combClass_win', 'responsive_counts', 'responsive_channels', '-v7.3');
% 
% fprintf('\nPermutation test results saved to:\n%s\n', fullfile(data_folder, save_name));
% fprintf('\nTotal responsive channels: %d\n', sum(responsive_counts));