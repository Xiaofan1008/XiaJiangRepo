%% ============================================================
%   SCRIPT A: SpikeFiltering_Save.m
%   Performs spike filtering ONLY, then saves sp_clipped
%   Xia, 2025
%% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% =============== USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq4_5ms_251125_154235';
Electrode_Type = 1;
FS = 30000;

% Filtering parameters
pos_limit = 100;
neg_limit = -100;
template_window_ms  = [5 15];
baseline_window_ms  = [-60 0];
corr_thresh         = 0.70;
SSD_threshold_factor = 16;

%% =============== Load Data =====================
cd(data_folder);
fprintf('Directory: %s\n', data_folder);

% --- Load spike file (.sp.mat) ---
sp_files = dir('*.sp.mat');
assert(~isempty(sp_files), 'No .sp.mat file found.');
S = load(sp_files(1).name);
sp = S.sp;

% --- Load triggers ---
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

% --- Load StimParams (PTD extraction only) ---
fileDIR = dir('*_exp_datafile_*.mat');
SS = load(fileDIR(1).name, 'StimParams','simultaneous_stim','n_Trials');
StimParams        = SS.StimParams;
simultaneous_stim = SS.simultaneous_stim;
n_Trials          = SS.n_Trials;

% PTD (µs)
if simultaneous_stim > 1
    PTD_us = cell2mat(StimParams(3:simultaneous_stim:end,6));
else
    PTD_us = zeros(n_Trials,1);
end

%% ============================================================
%   STEP 1 — BLANK SPIKES BETWEEN 0 → PTD
%% ============================================================
fprintf('\n=== STEP 1: Blanking spikes between 0 → PTD ===\n');
sp_blank = sp;

for ch = 1:numel(sp_blank)
    if isempty(sp_blank{ch}), continue; end
    spike_times = sp_blank{ch}(:,1);
    keep = true(size(spike_times));

    for tr = 1:length(trig)
        t0_ms  = trig(tr)/FS*1000;
        PTD_ms = PTD_us(tr)/1000;
        kill = spike_times >= t0_ms & spike_times < (t0_ms + PTD_ms);
        keep(kill) = false;
    end
    sp_blank{ch} = sp_blank{ch}(keep,:);
end

%% ============================================================
%   STEP 2 — SSD Waveform Consistency Filtering
%% ============================================================

fprintf('\n=== STEP 2: SSD waveform filtering ===\n');
sp_ssd = sp_blank;

for ch = 1:numel(sp_ssd)
    if isempty(sp_ssd{ch}), continue; end

    waves = sp_ssd{ch}(:,2:end);
    if size(waves,1) < 5
        continue;
    end

    mean_wave = mean(waves,1);
    SSD = sum((waves - mean_wave).^2,2);
    meanSSD = mean(SSD);

    keep = SSD <= SSD_threshold_factor * meanSSD;
    sp_ssd{ch} = sp_ssd{ch}(keep,:);
end

%% ============================================================
%   STEP 3 — TEMPLATE CORRELATION FILTERING (baseline only)
%% ============================================================

fprintf('\n=== STEP 3: Template correlation filtering ===\n');
sp_clipped = sp_ssd;

for ch = 1:numel(sp_clipped)
    if isempty(sp_clipped{ch}), continue; end

    waves = sp_clipped{ch}(:,2:end);
    times = sp_clipped{ch}(:,1);

    % Build template from evoked spikes
    evoked = false(size(times));
    for tr = 1:length(trig)
        t0 = trig(tr)/FS*1000;
        evoked = evoked | (times >= t0+template_window_ms(1) & times <= t0+template_window_ms(2));
    end

    if sum(evoked) < 5, continue; end
    template = mean(waves(evoked,:),1);

    % Correlation for all spikes
    corr_vals = zeros(size(times));
    for i = 1:numel(times)
        corr_vals(i) = corr(template(:), waves(i,:)');
    end

    keep = true(size(times));

    % Baseline trial-by-trial filtering
    for tr = 1:length(trig)
        t0 = trig(tr)/FS*1000;
        baseline_mask = (times >= t0 + baseline_window_ms(1)) & (times < t0 + baseline_window_ms(2));
        kill_idx = baseline_mask & (corr_vals < corr_thresh);
        keep(kill_idx) = false;
    end

    sp_clipped{ch} = sp_clipped{ch}(keep,:);
end

%% Save results
parts = split(data_folder, filesep);
folderName = parts{end};
base_name = folderName;
save_name = [base_name '.sp_xia.mat'];
save(save_name, 'sp_clipped');
fprintf('\nSaved filtered spikes to:\n%s\n', save_name);