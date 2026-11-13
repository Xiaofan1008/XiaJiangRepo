% ============================================================
% Compute mean ± SEM firing rate per channel and per stimulation set
% ============================================================
% This code calculates baseline-corrected firing rates for every channel
% and every stimulation set, stores mean, SEM, and raw trial values,
% and saves the results into a .mat file.
% The resulting .mat file can later be loaded and compared across
% stimulation types (Single, Simultaneous, Sequential).
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder
data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim4';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq4_5ms';

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
isSequential = contains(lower(data_folder), 'seq'); % if it is the sequential stimulation mode

%% Parameters
Spike_filtering = 0;
FS = 30000; % Sampling rate
baseline_window_ms = [-60, -5];
response_window_ms = [2, 10];
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

d = Depth_s(1); % channel depth mapping: 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

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

%% === Compute Mean ± SEM per Channel & Set === %%
FR_summary = struct;
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
    
    for ch = 1:nChn
        FR_thisCh = FR_corrected(ch,stim_mask);
        for i = 1:nAmp
            vals = FR_thisCh(amp_idx==i);
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            FR_summary(s).mean(ch,i) = mean(vals);
            FR_summary(s).sem(ch,i)  = std(vals)/sqrt(numel(vals));
            FR_summary(s).raw{ch,i}  = vals;
        end
    end
end

%% === Save Results === %%
save_name = sprintf('%s_FR_AllCh_AllSets.mat', base_name);
save(fullfile(data_folder, save_name),'FR_summary','Amps','E_MAP','uniqueComb','combClass_win','-v7.3');
fprintf('Firing rate summary saved to:\n%s\n', fullfile(data_folder, save_name));