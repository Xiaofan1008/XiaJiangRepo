% ============================================================
% Compute mean ± SEM firing rate per channel and per stimulation set
% 
% This code calculates baseline-corrected firing rates for every channel
% and stimulation set, performs a signed-rank test (p < 0.05) to detect
% significantly responsive channels, prints the results (channel numbers),
% and saves all outputs to a .mat file.
% ============================================================

clear
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ========= Enter Parameters ==========
% Choose Folder
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Single2_251106_114943';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim2_251106_120305';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq2_5ms_251106_120811';

selectedAmps = [4,5];
response_fraction_thresh = 0.4;  
if isempty(selectedAmps)
    error('No amplitude entered. Please specify at least one.');
end
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


%% === Compute Mean ± SEM per Channel & Set, Perform Signed-Rank Test (Amplitude-Specific) === %%
% FR_summary = struct;
% responsive_counts = zeros(1, nSets);
% responsive_channels = cell(1, nSets);
% 
% fprintf('\nPerforming signed-rank test with response fraction check (≥%.0f%% trials)...\n', ...
%     response_fraction_thresh*100);
% 
% fprintf('\nAvailable amplitudes: %s µA\n', num2str(unique(trialAmps)));
% 
% for s = 1:nSets
%     stim_mask_all = (combClass_win == s);
%     amps_this_all = trialAmps(stim_mask_all);
%     [amps_unique,~,amp_idx_all] = unique(amps_this_all);
%     nAmp = numel(amps_unique);
% 
%     FR_summary(s).stimSet = s;
%     FR_summary(s).amps = amps_unique;
%     FR_summary(s).mean = nan(nChn,nAmp);
%     FR_summary(s).sem  = nan(nChn,nAmp);
%     FR_summary(s).raw  = cell(nChn,nAmp);
%     FR_summary(s).stimCh = uniqueComb(s,uniqueComb(s,:)>0);
%     FR_summary(s).pval  = nan(nChn,nAmp);
%     FR_summary(s).sig   = false(nChn,nAmp);
%     FR_summary(s).respFrac = nan(nChn,nAmp);
% 
%     for amp_i = 1:nAmp
%         amp_val = amps_unique(amp_i);
%         if ~ismember(amp_val, selectedAmps)
%             continue; % skip amplitudes not selected
%         end
% 
%         stim_mask = stim_mask_all & (trialAmps == amp_val);
%         if sum(stim_mask) < 3, continue; end  % skip amplitudes with too few trials
% 
%         fprintf('\n--- Stim Set %d [Ch %s] | %.0f µA ---\n', ...
%             s, num2str(FR_summary(s).stimCh), amp_val);
% 
%         for ch = 1:nChn
%             FR_base = FR_baseline(ch, stim_mask);
%             FR_resp = FR_response(ch, stim_mask);
%             if all(isnan(FR_base)) || all(isnan(FR_resp)), continue; end
%             valid_idx = ~isnan(FR_base) & ~isnan(FR_resp);
%             FR_base = FR_base(valid_idx);
%             FR_resp = FR_resp(valid_idx);
% 
%             % --- Step 1: Compute response fraction ---
%             nTrials_valid = numel(FR_base);
%             if nTrials_valid < 3, continue; end
%             n_higher = sum(FR_resp > FR_base);
%             resp_fraction = n_higher / nTrials_valid;
%             FR_summary(s).respFrac(ch,amp_i) = resp_fraction;
% 
%             % --- Step 2: Compute firing rate mean ± SEM (always computed) ---
%             vals = FR_corrected(ch, stim_mask);
%             vals = vals(~isnan(vals));
%             if ~isempty(vals)
%                 FR_summary(s).mean(ch,amp_i) = mean(vals);
%                 FR_summary(s).sem(ch,amp_i)  = std(vals)/sqrt(numel(vals));
%                 FR_summary(s).raw{ch,amp_i}  = vals;
%             end
% 
%             % --- Step 3: Skip significance test if response fraction too low ---
%             if resp_fraction < response_fraction_thresh
%                 FR_summary(s).pval(ch,amp_i) = NaN;
%                 FR_summary(s).sig(ch,amp_i)  = false;
%                 continue;
%             end
% 
%             % --- Step 4: Signed-rank test (right-tailed) ---
%             try
%                 p = signrank(FR_resp, FR_base, 'tail', 'right');
%             catch
%                 p = NaN;
%             end
%             FR_summary(s).pval(ch,amp_i) = p;
%             FR_summary(s).sig(ch,amp_i)  = (p < 0.05);
% 
% 
%             % --- ZERO-AWARE SIGN TEST ---
%             % diffvals = FR_resp - FR_base;
%             % 
%             % n_pos  = sum(diffvals > 0);
%             % n_neg  = sum(diffvals < 0);
%             % n_zero = sum(diffvals == 0);
%             % 
%             % N = n_pos + n_neg + n_zero;
%             % 
%             % % If all equal, no responsiveness
%             % if N == 0 || n_pos == 0
%             %     p = 1;
%             % else
%             %     p = 1 - binocdf(n_pos - 1, N, 0.5);
%             % end
%             % 
%             % FR_summary(s).pval(ch,amp_i) = p;
%             % FR_summary(s).sig(ch,amp_i)  = (p < 0.05);
% 
%         end
% 
%         % --- Step 5: Report significant channels for this amplitude ---
%         responsive_idx = find(FR_summary(s).sig(:,amp_i));
%         responsive_counts(s) = responsive_counts(s) + numel(responsive_idx);
%         responsive_channels{s, amp_i} = responsive_idx;
% 
%         if ~isempty(responsive_idx)
%             fprintf('→ %.0f µA: %d responsive channels (≥%.0f%% resp trials)\n', ...
%                 amp_val, numel(responsive_idx), response_fraction_thresh*100);
%             fprintf('   Channels: %s\n', num2str(responsive_idx(:)'));
%         else
%             fprintf('→ %.0f µA: 0 responsive channels.\n', amp_val);
%         end
%     end
% end

%% === Compute Mean ± SEM per Channel & Set, Perform Paired t-test (Amplitude-Specific) === %%
FR_summary = struct;
responsive_counts = zeros(1, nSets);
responsive_channels = cell(1, nSets);

fprintf('\nPerforming paired t-test with response fraction check (≥%.0f%% trials)...\n', ...
    response_fraction_thresh*100);

fprintf('\nAvailable amplitudes: %s µA\n', num2str(unique(trialAmps)));

for s = 1:nSets
    stim_mask_all = (combClass_win == s);
    amps_this_all = trialAmps(stim_mask_all);
    [amps_unique,~,amp_idx_all] = unique(amps_this_all);
    nAmp = numel(amps_unique);

    FR_summary(s).stimSet = s;
    FR_summary(s).amps = amps_unique;
    FR_summary(s).mean = nan(nChn,nAmp);
    FR_summary(s).sem  = nan(nChn,nAmp);
    FR_summary(s).raw  = cell(nChn,nAmp);
    FR_summary(s).stimCh = uniqueComb(s,uniqueComb(s,:)>0);
    FR_summary(s).pval  = nan(nChn,nAmp);
    FR_summary(s).sig   = false(nChn,nAmp);
    FR_summary(s).respFrac = nan(nChn,nAmp);

    for amp_i = 1:nAmp
        amp_val = amps_unique(amp_i);

        if ~ismember(amp_val, selectedAmps)
            continue;
        end
        stim_mask = stim_mask_all & (trialAmps == amp_val);
        if sum(stim_mask) < 3
            continue;
        end
        fprintf('\n--- Stim Set %d [Ch %s] | %.0f µA ---\n', ...
            s, num2str(FR_summary(s).stimCh), amp_val);
        for ch = 1:nChn
            
            FR_base = FR_baseline(ch, stim_mask);
            FR_resp = FR_response(ch, stim_mask);

            if all(isnan(FR_base)) || all(isnan(FR_resp)), continue; end

            valid_idx = ~isnan(FR_base) & ~isnan(FR_resp);
            FR_base = FR_base(valid_idx);
            FR_resp = FR_resp(valid_idx);

            if numel(FR_base) < 3
                continue;
            end

            %% --- Step 1: Compute response fraction ---
            nTrials_valid = numel(FR_base);
            n_higher = sum(FR_resp > FR_base);
            resp_fraction = n_higher / nTrials_valid;
            FR_summary(s).respFrac(ch,amp_i) = resp_fraction;

            %% --- Step 2: Always compute FR mean ± SEM ---
            vals = FR_corrected(ch, stim_mask);
            vals = vals(~isnan(vals));
            if ~isempty(vals)
                FR_summary(s).mean(ch,amp_i) = mean(vals);
                FR_summary(s).sem(ch,amp_i)  = std(vals)/sqrt(numel(vals));
                FR_summary(s).raw{ch,amp_i}  = vals;
            end

            %% --- Step 3: Require a minimum response fraction ---
            if resp_fraction < response_fraction_thresh
                FR_summary(s).pval(ch,amp_i) = NaN;
                FR_summary(s).sig(ch,amp_i)  = false;
                continue;
            end

            % ---------------------------------------------------------
            % --- Step 4: Paired t-test: (resp > baseline)?
            %     H0: FR_resp - FR_base <= 0
            %     H1: FR_resp - FR_base > 0
            % ---------------------------------------------------------
            try
                [~,p] = ttest(FR_resp, FR_base, 'Tail', 'right');
            catch
                p = NaN;
            end

            FR_summary(s).pval(ch,amp_i) = p;
            FR_summary(s).sig(ch,amp_i)  = (p < 0.05);

        end

        %% --- Step 5: Report significant channels ---
        responsive_idx = find(FR_summary(s).sig(:,amp_i));
        responsive_counts(s) = responsive_counts(s) + numel(responsive_idx);
        responsive_channels{s, amp_i} = responsive_idx;

        if ~isempty(responsive_idx)
            fprintf('→ %.0f µA: %d responsive channels (≥%.0f%% resp trials)\n', ...
                amp_val, numel(responsive_idx), response_fraction_thresh*100);
            fprintf('   Channels: %s\n', num2str(responsive_idx(:)'));
        else
            fprintf('→ %.0f µA: 0 responsive channels.\n', amp_val);
        end
    end
end

%% === Compute Population Firing Rate for Selected Amplitudes & Plot Tuning Curves ===
fprintf('\n===== Computing Population FR (Significant Channels Only) =====\n');

for s = 1:nSets

    amps_unique = FR_summary(s).amps;
    nAmp        = numel(amps_unique);

    % containers
    avgFR_sig   = nan(1, nAmp);   % average FR across sig channels
    semFR_sig   = nan(1, nAmp);   % SEM across sig channels
    nSig_amp    = zeros(1,nAmp);  % number of sig channels

    for amp_i = 1:nAmp

        amp_val = amps_unique(amp_i);

        % skip un-selected amplitudes
        if ~ismember(amp_val, selectedAmps)
            continue;
        end

        sig_channels = find(FR_summary(s).sig(:,amp_i));

        if isempty(sig_channels)
            fprintf('Set %d | %.0f µA → No sig channels, skipping.\n', s, amp_val);
            continue;
        end

        nSig_amp(amp_i) = numel(sig_channels);

        vals = FR_summary(s).mean(sig_channels, amp_i);
        vals = vals(~isnan(vals));

        avgFR_sig(amp_i) = mean(vals);
        semFR_sig(amp_i) = std(vals) / sqrt(numel(vals));

        fprintf('Set %d | %.0f µA → AVG FR = %.3f sp/s (SEM=%.3f, n=%d)\n', ...
            s, amp_val, avgFR_sig(amp_i), semFR_sig(amp_i), nSig_amp(amp_i));
    end

    % ---- Save into FR_summary ----
    FR_summary(s).avgSigFR  = avgFR_sig;
    FR_summary(s).avgSigSEM = semFR_sig;
    FR_summary(s).nSigCh    = nSig_amp;

    % % ---- Plot tuning curve for this stim set ----
    % figure('Color','w'); hold on;
    % amp_vals_plot = amps_unique(ismember(amps_unique, selectedAmps));
    % avg_vals_plot = avgFR_sig(ismember(amps_unique, selectedAmps));
    % sem_vals_plot = semFR_sig(ismember(amps_unique, selectedAmps));
    % 
    % errorbar(amp_vals_plot, avg_vals_plot, sem_vals_plot, ...
    %     '-o', 'LineWidth', 2, 'MarkerSize', 7);
    % 
    % xlabel('Amplitude (µA)', 'FontSize', 9);
    % ylabel('Mean Firing Rate (sp/s)', 'FontSize', 9);
    % title(sprintf('Population FR Tuning Curve — Set %d (Stim Ch %s)', ...
    %     s, num2str(FR_summary(s).stimCh)));
    % box off;

end
fprintf('===== Population FR finished. =====\n');

%% === Save Results ===
save_name = sprintf('%s_FR_SigCh_ByAmp.mat', base_name);
save(fullfile(data_folder, save_name), 'FR_summary', 'Amps', 'E_MAP', ...
    'uniqueComb', 'combClass_win', 'responsive_channels', 'selectedAmps', ...
    '-v7.3');

fprintf('\nResults saved to:\n%s\n', fullfile(data_folder, save_name));