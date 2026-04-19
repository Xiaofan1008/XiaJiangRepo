%% ========================================================================
%  Spike Count Summary Extractor (Trial-Centric)
%  Extracts baseline and post-stimulus spike counts per trial,
%  grouped strictly by Set -> Amp -> PTD -> Trial.
%  Shows all 64 channels in a single table per trial.
% ========================================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ========================
data_folder      = '/Volumes/MACData/Data/Data_Xia/DX021/Xia_ISI_SimSeq2';
Electrode_Type   = 1;          % 0 rigid, 1 single-shank flex, 2 four-shank flex

raster_chn_start = 1;          % Depth_s index start
raster_chn_end   = 32;         % Depth_s index end

% ---- Counting Windows (Relative to 1st Pulse Trigger) ----
baseline_win_ms  = [-50 0];    % Window to count baseline spikes
post_win_ms      = [2 40];   % Window to count evoked/response spikes

FS = 30000;                    % Sampling rate

%% ===================== INITIAL SETUP =========================
if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end
cd(data_folder);
fprintf('Extracting spike counts in folder:\n%s\n\n', data_folder);

% ---- Extract base_name ----
parts       = split(data_folder, filesep);
last_folder = parts{end};
u           = strfind(last_folder, '_');
if numel(u) > 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

%% ===================== LOAD SPIKES ===========================
ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];

if isfile(ssd_file)
    S = load(ssd_file);
    if     isfield(S,'sp_pca'),  sp = S.sp_pca;
    elseif isfield(S,'sp_corr'), sp = S.sp_corr;
    elseif isfield(S,'sp_SSD'),  sp = S.sp_SSD;
    else, error('SSD file missing usable variable.');
    end
    fprintf('Loaded SSD filtered spikes from %s\n', ssd_file);
elseif isfile(base_file)
    S = load(base_file);
    if     isfield(S,'sp_clipped'), sp = S.sp_clipped;
    elseif isfield(S,'sp'),         sp = S.sp;
    else, error('Base spike file missing usable variable.');
    end
    fprintf('Loaded base spikes from %s\n', base_file);
else
    error('No spike file %s or %s found.', ssd_file, base_file);
end
nCh = numel(sp);

%% ===================== LOAD TRIGGERS =========================
if isempty(dir('*.trig.dat'))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);   
nTrig = numel(trig);

%% ===================== LOAD StimParams & DECODE ==============
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat file found.');
S = load(fileDIR(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');

StimParams = S.StimParams;
sim_stim   = S.simultaneous_stim;
E_MAP      = S.E_MAP;
n_Trials   = S.n_Trials;

% ---- Amplitudes ----
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:sim_stim:end);
[Amps,~,ampIdx] = unique(trialAmps);
Amps(Amps == -1) = 0;
nAMP = numel(Amps);

% ---- PTDs ----
if sim_stim > 1
    PTD_all = cell2mat(StimParams(3:sim_stim:end,6));   % µs
else
    PTD_all = zeros(n_Trials,1);    
end
[PTDs,~,ptdIdx] = unique(PTD_all);
nPTD = numel(PTDs);

% ---- Stimulation sets ----
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));
stimSeq = zeros(n_Trials, sim_stim);
for t = 1:n_Trials
    rr = (t-1)*sim_stim + (1:sim_stim);
    v  = idx_all(rr);
    v  = v(v>0);
    stimSeq(t,1:numel(v)) = v;
end
[uniqueComb,~,combClass] = unique(stimSeq, 'rows','stable');
nSets = size(uniqueComb,1);

% ---- Electrode Map ----
d = Depth_s(Electrode_Type);   
depth_range = raster_chn_start : min(raster_chn_end, numel(d));

%% ===================== MAIN COUNTING LOOP ====================
SpikeCounts = struct(); % Initialize the master output structure
fprintf('\nStarting spike counting loop...\n');

for si = 1:nSets
    % [METADATA LEVEL 1: SET]
    stim_channels = uniqueComb(si, uniqueComb(si,:)>0);
    set_label = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stim_channels, 'UniformOutput', false), ' + ');
    
    SpikeCounts.set(si).set_name      = set_label;
    SpikeCounts.set(si).stim_channels = stim_channels;
    
    for ai = 1:nAMP
        % [METADATA LEVEL 2: AMP]
        amp_val = Amps(ai);
        
        SpikeCounts.set(si).amp(ai).amp_label = sprintf('%.1f uA', amp_val);
        SpikeCounts.set(si).amp(ai).amp_value = amp_val;
        
        for pi = 1:nPTD
            % [METADATA LEVEL 3: PTD]
            ptd_ms_val = PTDs(pi) / 1000;
            trials_this_condition = find(combClass==si & ampIdx==ai & ptdIdx==pi);
            
            cond_str = sprintf('Set: %s | Amp: %.1f uA | PTD: %.1f ms', set_label, amp_val, ptd_ms_val);
            
            SpikeCounts.set(si).amp(ai).ptd(pi).condition_summary = cond_str;
            SpikeCounts.set(si).amp(ai).ptd(pi).PTD_ms            = ptd_ms_val;
            SpikeCounts.set(si).amp(ai).ptd(pi).total_trials      = numel(trials_this_condition);
            
            % Skip if this condition combination never happened
            if isempty(trials_this_condition)
                continue;
            end
            
            nTr_cond = numel(trials_this_condition);
            nChns    = numel(depth_range);
            
            % [MODIFIED] SWAPPED LOOP ORDER: Loop by Trial first
            for k = 1:nTr_cond
                
                tr = trials_this_condition(k); % Absolute trial ID
                t0 = trig(tr) / FS * 1000;     % Trigger time of 1st pulse in ms
                
                % Define exact windows anchored to the 1st pulse
                win_base_start = t0 + baseline_win_ms(1);
                win_base_end   = t0 + baseline_win_ms(2);
                
                win_post_start = t0 + post_win_ms(1);
                win_post_end   = t0 + post_win_ms(2);
                
                % [METADATA LEVEL 4: TRIAL]
                SpikeCounts.set(si).amp(ai).ptd(pi).trial(k).relative_trial_ID = k;
                SpikeCounts.set(si).amp(ai).ptd(pi).trial(k).absolute_trial_ID = tr;
                
                % Pre-allocate 64-row arrays for this specific trial's channels
                % Using zeros() ensures dead channels will automatically have 0 spikes
                depth_idx_arr = zeros(nChns, 1);
                intan_ch_arr  = zeros(nChns, 1);
                base_spikes   = zeros(nChns, 1);
                post_spikes   = zeros(nChns, 1);
                
                % [MODIFIED] Loop through all Channels for this specific trial
                for idxDepth = 1:nChns
                    
                    ich = depth_range(idxDepth);   % Depth index (1-64)
                    ch  = d(ich);                  % Intan hardware channel
                    
                    depth_idx_arr(idxDepth) = ich;
                    intan_ch_arr(idxDepth)  = ch;
                    
                    % If channel is dead/empty, it skips counting and leaves the array at 0
                    if ch < 1 || ch > nCh || isempty(sp{ch})
                        continue;
                    end
                    
                    sp_times = sp{ch}(:,1); % All spike times for this channel in ms
                    
                    % Count spikes falling in the baseline window
                    base_spikes(idxDepth) = sum(sp_times >= win_base_start & sp_times < win_base_end);
                    
                    % Count spikes falling in the post-stimulus window
                    post_spikes(idxDepth) = sum(sp_times >= win_post_start & sp_times < win_post_end);
                    
                end % End Channel Loop
                
                % Combine the 64 channel rows into a single table for this trial
                Channel_Results = table(depth_idx_arr, intan_ch_arr, base_spikes, post_spikes, ...
                    'VariableNames', {'Depth_Index', 'Intan_Channel', 'Baseline_Spikes', 'Post_Spikes'});
                
                SpikeCounts.set(si).amp(ai).ptd(pi).trial(k).Channel_Results = Channel_Results;
                
            end % End Trial Loop
        end % End PTD Loop
    end % End Amp Loop
end % End Set Loop

%% ===================== SAVE RESULTS ==========================
save_name = sprintf('%s_ISI_TrialSpikeCounts_PerTrial.mat', base_name);
full_save_path = fullfile(data_folder, save_name);

save(full_save_path, 'SpikeCounts', 'baseline_win_ms', 'post_win_ms', ...
     'Amps', 'PTDs', 'uniqueComb');

fprintf('\nFinished. Extracted spike counts for all trials.\n');
fprintf('Data saved to: %s\n', full_save_path);