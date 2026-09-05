% ============================================================
% Compute mean ± SEM firing rate per channel and per stimulation set
% ============================================================
% This code calculates baseline-corrected firing rates for every channel
% and stimulation set, performs a signed-rank test (p < 0.05) to detect
% significantly responsive channels, prints the results (channel numbers),
% and saves all outputs to a .mat file.
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

% Choose Folder
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim1';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq1';

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
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

%% Parameters
Spike_filtering = 0;
FS = 30000; % Sampling rate
baseline_window_ms = [-60, -5];
response_window_ms = [2, 15];
pos_limit = 100;
neg_limit = -100;

%% Load spike data
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

%% Load StimParams and decode sets
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
S = load(fullfile(data_folder, fileDIR(1).name), ...
    'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
StimParams = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP = S.E_MAP;
n_Trials = S.n_Trials;

% Decode trial structure
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps,~,ampIdx] = unique(trialAmps(:));
if any(Amps==-1), Amps(Amps==-1)=0; end
n_AMP = numel(Amps);

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

pulseTrain_all = cell2mat(StimParams(2:end,9));
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);
[PulsePeriods,~,pulseIdx] = unique(pulseTrain(:));

d = Depth_s(1); % channel depth mapping

%% Spike Amplitude Filtering (optional)
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

%% === Compute Firing Rate for Each Channel & Set === %%
fprintf('\nComputing firing rates for all channels and all stim sets...\n');

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

%% === Compute Mean ± SEM per Channel & Set, Perform Signed-Rank Test === %%
% FR_summary = struct;
% responsive_counts = zeros(1, nSets);
% responsive_channels = cell(1, nSets);
% 
% fprintf('\nPerforming signed-rank test for responsiveness...\n');
% 
% for s = 1:nSets
%     stim_mask = (combClass_win == s);
%     amps_this = trialAmps(stim_mask);
%     [amps_unique,~,amp_idx] = unique(amps_this);
%     nAmp = numel(amps_unique);
% 
%     FR_summary(s).stimSet = s;
%     FR_summary(s).amps = amps_unique;
%     FR_summary(s).mean = nan(nChn,nAmp);
%     FR_summary(s).sem  = nan(nChn,nAmp);
%     FR_summary(s).raw  = cell(nChn,nAmp);
%     FR_summary(s).stimCh = uniqueComb(s,uniqueComb(s,:)>0);
%     FR_summary(s).pval  = nan(nChn,1);
%     FR_summary(s).sig   = false(nChn,1);
% 
%     for ch = 1:nChn
%         FR_base = FR_baseline(ch, stim_mask);
%         FR_resp = FR_response(ch, stim_mask);
%         if all(isnan(FR_base)) || all(isnan(FR_resp)), continue; end
% 
%         % --- Wilcoxon signed-rank test (right-tailed) ---
%         try
%             p = signrank(FR_resp, FR_base, 'tail', 'right');
%         catch
%             p = NaN;
%         end
%         FR_summary(s).pval(ch) = p;
%         FR_summary(s).sig(ch)  = (p < 0.05);
% 
%         % Compute mean ± SEM
%         for i = 1:nAmp
%             vals = FR_corrected(ch, stim_mask);
%             vals = vals(amp_idx == i);
%             vals = vals(~isnan(vals));
%             if isempty(vals), continue; end
%             FR_summary(s).mean(ch,i) = mean(vals);
%             FR_summary(s).sem(ch,i)  = std(vals)/sqrt(numel(vals));
%             FR_summary(s).raw{ch,i}  = vals;
%         end
%     end
% 
%     % --- Count and list responsive channels ---
%     responsive_idx = find(FR_summary(s).sig);
%     responsive_counts(s) = numel(responsive_idx);
%     responsive_channels{s} = responsive_idx;
% 
%     if ~isempty(responsive_idx)
%         fprintf('Stim Set %d [Ch %s]: %d responsive channels (p<0.05)\n', ...
%             s, num2str(FR_summary(s).stimCh), responsive_counts(s));
%         fprintf('→ Responding channel numbers: %s\n', num2str(responsive_idx(:)'));
%     else
%         fprintf('Stim Set %d [Ch %s]: 0 responsive channels.\n', ...
%             s, num2str(FR_summary(s).stimCh));
%     end
% end
% 
% %% === Save Results === %%
% save_name = sprintf('%s_FR_AllCh_AllSets.mat', base_name);
% save(fullfile(data_folder, save_name), 'FR_summary', 'Amps', 'E_MAP', ...
%     'uniqueComb', 'combClass_win', 'responsive_counts', 'responsive_channels', '-v7.3');
% 
% fprintf('\nFiring rate summary and test results saved to:\n%s\n', ...
%     fullfile(data_folder, save_name));
% fprintf('\nTotal responsive channels across all sets: %d\n', sum(responsive_counts));


%% === Compute Mean ± SEM per Channel & Set, Perform Signed-Rank Test & Response freaction threshold === %%
FR_summary = struct;
responsive_counts = zeros(1, nSets);
responsive_channels = cell(1, nSets);

response_fraction_thresh = 0.1; 

fprintf('\nPerforming signed-rank test with response fraction check (%.0f%% threshold)...\n', ...
    response_fraction_thresh*100);

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
    FR_summary(s).respFrac = nan(nChn,1); % store response fraction

    for ch = 1:nChn
        FR_base = FR_baseline(ch, stim_mask);
        FR_resp = FR_response(ch, stim_mask);

        if all(isnan(FR_base)) || all(isnan(FR_resp)), continue; end
        valid_idx = ~isnan(FR_base) & ~isnan(FR_resp);
        FR_base = FR_base(valid_idx);
        FR_resp = FR_resp(valid_idx);

        % --- Step 1: Compute response fraction ---
        nTrials_valid = numel(FR_base);
        if nTrials_valid < 3
            continue; % skip channels with too few trials
        end
        n_higher = sum(FR_resp > FR_base);
        resp_fraction = n_higher / nTrials_valid;
        FR_summary(s).respFrac(ch) = resp_fraction;

        % --- Step 2: Apply fraction threshold ---
        if resp_fraction < response_fraction_thresh
            FR_summary(s).pval(ch) = NaN;
            FR_summary(s).sig(ch)  = false;
            continue;
        end

        % --- Step 3: Wilcoxon signed-rank test (right-tailed) ---
        try
            p = signrank(FR_resp, FR_base, 'tail', 'right');
        catch
            p = NaN;
        end
        FR_summary(s).pval(ch) = p;
        FR_summary(s).sig(ch)  = (p < 0.05);

        % --- Step 4: Compute mean ± SEM per amplitude ---
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

    % --- Step 5: Report responsive channels ---
    responsive_idx = find(FR_summary(s).sig);
    responsive_counts(s) = numel(responsive_idx);
    responsive_channels{s} = responsive_idx;

    if ~isempty(responsive_idx)
        fprintf('Stim Set %d [Ch %s]: %d responsive channels (p<0.05, ≥%.0f%% resp trials)\n', ...
            s, num2str(FR_summary(s).stimCh), responsive_counts(s), response_fraction_thresh*100);
        fprintf('→ Responding channels: %s\n', num2str(responsive_idx(:)'));
    else
        fprintf('Stim Set %d [Ch %s]: 0 responsive channels.\n', ...
            s, num2str(FR_summary(s).stimCh));
    end
end

%% === Save Results === %%
% save_name = sprintf('%s_FR_AllCh_AllSets_FRcheck.mat', base_name);
% save(fullfile(data_folder, save_name), 'FR_summary', 'Amps', 'E_MAP', ...
%     'uniqueComb', 'combClass_win', 'responsive_counts', 'responsive_channels', '-v7.3');

fprintf('\nResults saved to:\n%s\n', fullfile(data_folder, save_name));
fprintf('\nTotal responsive channels (≥%.0f%% rule): %d\n', ...
    response_fraction_thresh*100, sum(responsive_counts));