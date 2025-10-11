%% Quick Analysis Code
clear all
close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
% addpath(genpath
% ('/Volumes/MACData/Data/Data_Sabrina/Experimental_Design'));


%% ================================================ %%
%%                 Data Processing                  %%
%% ================================================ %%

%  ------------------ Parameters ------------------  %
% Intan parameters
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

dName='amplifier';
vFID = fopen([filepath filesep dName '.dat'],'r'); % read data file

% Quick Analysis parameters 
window_ms = [-100,100]; % time window for data extract
window_samps = round(window_ms/1000 * FS);
total_sample = diff(window_samps) + 1;
trialLength = diff(window_samps)+1;

% ---- Artifact blank parameters ---- %
% Enable/disable second stimulation artifact blanking
enable_second_artifact = 1;

% 1st stimulation artifact (typically 0 ms)
preBlank1_ms   = 1;    % ms to blank BEFORE 1st stim
postBlank1_ms  = 2;    % ms to blank AFTER  1st stim
preBlank1_samp  = round(preBlank1_ms/1000 * FS);
postBlank1_samp = round(postBlank1_ms/1000 * FS);

% 2nd stimulation artifact (e.g. 5.5 ms)
secondStim_ms      = 5.5;  % onset time relative to stim trigger (ms)
preBlank2_ms        = 0.5; % ms to blank BEFORE 2nd stim
postBlank2_ms       = 2.0; % ms to blank AFTER  2nd stim
preBlank2_samp      = round(preBlank2_ms/1000 * FS);
postBlank2_samp     = round(postBlank2_ms/1000 * FS);
secondStim_samp     = round(secondStim_ms/1000 * FS);

% Cross-fade tapers at the edges of the filled region
taper_pre_ms  = 0.8;   % ms, PRE side
taper_post_ms = 0.8;   % ms, POST side

% Detection guard around the filled segment
guard_pre_ms  = max(taper_pre_ms, 0.5);
guard_post_ms = max(taper_post_ms, 0.5);
guard_pre_samp  = max(1, round(guard_pre_ms/1000  * FS));
guard_post_samp = max(1, round(guard_post_ms/1000 * FS));

% Trigger and parameters
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);
trig_all = trig; 
d = Depth_s(1); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex
% Extract trigger
nTrials = 440; % number of trials used for analysis
Trial_start = 1;
if (nTrials<=length(trig)); trig = trig(Trial_start:Trial_start+nTrials-1); end
nTrig = length(trig); % length of trigger

% load trial parameters
TrialParams = loadTrialParams;
trialIDs = cell2mat(TrialParams(:,2));
trialIDs = trialIDs(1:nTrials);

% load StimParameters
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
fileDIR = fileDIR(1).name;
S = load(fileDIR,'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
CHN = S.CHN;
E_MAP = S.E_MAP;
n_Trials = S.n_Trials;

% extract amplitude
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps_all = trialAmps_all(1:simultaneous_stim:end);
trialAmps = trialAmps_all(Trial_start : Trial_start + nTrials - 1);

% Group by amplitude
[uniqAmps,~,ampIdx] = unique(trialAmps(:));
Amps = uniqAmps;
if any(Amps == -1)        % treat -1 as 0 µA, if present
    Amps(Amps == -1) = 0;
end
n_AMP = numel(Amps);

% Time vector for plotting
t_axis = (window_samps(1):window_samps(2)) / FS * 1000;

% Preallocate
MUA_all = zeros(nChn, trialLength, nTrials);
Spikes = cell(nChn, nTrials); 

% Spike waveform window ~3 ms (e.g., 1 ms pre, 2 ms post)
preS  = round(1.0e-3 * FS);
postS = round(3.0e-3 * FS);

% For better alignment, local negative peak window
alignRad = max(1, round(3e-3 * FS));

% Spike detection parameters
pp_min_uV   = 30;           % min peak-to-peak 
pp_max_uV   = 500;          % max peak-to-peak
tp_min_ms   = 0.30;         % trough-to-peak min width
tp_max_ms   = 1.20;         % trough-to-peak max width
tp_min_samp = round(tp_min_ms/1000 * FS);
tp_max_samp = round(tp_max_ms/1000 * FS);

%% ------------------ Filter design ------------------ %% 
centerFreq = 4000;   % Hz
halfBandwidth = 3700; % Hz
f1 = centerFreq - halfBandwidth;  % lower cutoff
f2 = centerFreq + halfBandwidth;  % upper cutoff
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency1', 300, ...
    'HalfPowerFrequency2', 3000, ...
    'SampleRate', FS);

% ------------------ Filter design (causal FIR) ------------------ %
% Linear-phase FIR band-pass; causal -> no pre-0 leakage
hp = 300;                        % Hz
lp = 3000;                       % Hz
Fpass = [hp lp] / (FS/2);        % normalized
firOrder = 512;                  % >= ~FS/hp
bpFIR = fir1(firOrder, Fpass, 'bandpass', hamming(firOrder+1), 'scale');

% Group delay (samples) of linear-phase FIR = firOrder/2 (constant)
gd = firOrder/2;                 % integer samples

%% ------------------ Threshold calculation ------------------ %% 
thr = zeros(nChn, 1);
thres_duration_sec   = 10;   % Length of data to use for threshold calc
num_samples_threshold = min( ...
    floor(thres_duration_sec * FS), ...
    floor(dir([filepath filesep dName '.dat']).bytes / (nChn * 2)) ...
); % number of threshold samples
assert(num_samples_threshold > 5*FS, ...
    'Recording shorter than ~5 s; not enough for thresholds.'); % warrning
% read data 
fseek(vFID, 0, 'bof');
data_thres_uV = fread(vFID, [nChn, num_samples_threshold], 'int16') * 0.195;
% Loop over channels
for ch_idx = 1:nChn
    thisChn = d(ch_idx);   % Depth-mapped channel index
    x = double(data_thres_uV(thisChn, :));
    
    x_fir = filter(bpFIR, 1, x); % filtering
    if gd > 0 % compensate filter delay
        x_fir = [x_fir(gd+1:end), zeros(1, gd)];
    end

    sigma = median(abs(x_fir)) / 0.6745;
    thr(ch_idx) = -4 * sigma; % Compute threshold as -4.5 × standard deviation    
end

%% ------------------ Data Processing ------------------ %%
for tr = 1:nTrials
    fprintf('Reading trial %d/%d\n', tr, nTrials);
    start_index = trig(tr) + window_samps(1);
    start_byte = start_index * nChn * 2; % 2 bytes per int16
    fseek(vFID, start_byte, 'bof');
    data = fread(vFID, [nChn, total_sample], 'int16') * 0.195; % in uV

    stim_idx = -window_samps(1) + 1;
    % First artifact blanking window
    blank1 = max(1, stim_idx - preBlank1_samp) : ...
             min(total_sample, stim_idx + postBlank1_samp);    
    % Second artifact blanking window (optional)
    if enable_second_artifact
        stim2_idx = stim_idx + secondStim_samp;
        blank2 = max(1, stim2_idx - preBlank2_samp) : ...
                 min(total_sample, stim2_idx + postBlank2_samp);
    else
        blank2 = [];
    end    
    % Merge blanking ranges
    blank_range = unique([blank1, blank2]);
    keep_range  = setdiff(1:total_sample, blank_range);
    
    for ch = 1:nChn       
        thisChn = d(ch);
        raw = data(thisChn, :);
        %% ---- Artifact Blank (Multi-window) ---- %%
        stimSample = -window_samps(1) + 1;  % stim index (1st stim) in trial window
        
        % Define artifact blanking windows
        blank1 = max(1, stimSample - preBlank1_samp) : ...
                 min(total_sample, stimSample + postBlank1_samp);
        
        if enable_second_artifact
            stim2_idx = stimSample + secondStim_samp;
            blank2 = max(1, stim2_idx - preBlank2_samp) : ...
                     min(total_sample, stim2_idx + postBlank2_samp);
        else
            blank2 = [];
        end
        
        % Merge and interpolate over both blank regions
        blank_range = unique([blank1, blank2]);
        keep_range  = setdiff(1:total_sample, blank_range);
        
        raw_interp = raw;
        raw_interp(blank_range) = interp1(keep_range, raw(keep_range), blank_range, 'linear', 'extrap');
        
        % Filter & rectify
        sig_filt_plot = filter(bpFIR, 1, raw_interp);
        if gd > 0
            sig_filt_plot = [sig_filt_plot(gd+1:end), zeros(1, gd)];
        end
        
        smoothed = movmean(sig_filt_plot, 15);   % Optional: smoothing for plot
        MUA_all(ch, :, tr) = sig_filt_plot;      % Store filtered trace
        %% ---- Baseline firing rate ---- %%       
        sig_filt = filter(bpFIR, 1, raw_interp);
        % Compensate constant group delay
        if gd > 0
            signal = [sig_filt(gd+1:end), zeros(1, gd)];
        else
            signal = sig_filt;
        end

        %% ---- Spike detection ---- %%       
        % Threshold crossings
        cross_idx = find(signal(2:end) < thr(ch) & signal(1:end-1) >= thr(ch)) + 1;

        % Guard region around 1st artifact
        guard1_start = max(1, stim_idx - preBlank1_samp - guard_pre_samp);
        guard1_end   = min(total_sample, stim_idx + postBlank1_samp + guard_post_samp);
        
        % Guard region around 2nd artifact (optional)
        if enable_second_artifact
            stim2_idx = stim_idx + secondStim_samp;
            guard2_start = max(1, stim2_idx - preBlank2_samp - guard_pre_samp);
            guard2_end   = min(total_sample, stim2_idx + postBlank2_samp + guard_post_samp);
        else
            guard2_start = [];
            guard2_end   = [];
        end
        
        % Combine both guard masks
        mask_keep = true(size(cross_idx));
        if ~isempty(guard2_start)
            mask_keep = ~(...
                (cross_idx >= guard1_start & cross_idx <= guard1_end) | ...
                (cross_idx >= guard2_start & cross_idx <= guard2_end) );
        else
            mask_keep = ~(cross_idx >= guard1_start & cross_idx <= guard1_end);
        end
        cross_idx = cross_idx(mask_keep);


        % Refractory period 1ms
        refrac_samp = max(1, round(1e-3 * FS));
        % Space crossings: drop any within refractory of the previous kept one
        if numel(cross_idx) > 1
            keep = [true, diff(cross_idx) > refrac_samp];
            cross_idx = cross_idx(keep);
        end

        % --- Locate true troughs without local-window constraint ---
        trough_idx = [];    % absolute sample indices of troughs
        for k = 1:numel(cross_idx)
            i0 = cross_idx(k);
        
            % Full extraction window around the crossing
            a = max(1, i0 - preS);
            b = min(total_sample, i0 + postS);
        
            % Ensure full data window
            if (a + preS > b - postS)
                continue
            end
        
            % Most negative sample anywhere in that full window
            [~, rel] = min(signal(a:b));
            ip = a + rel - 1;   % absolute index
            trough_idx(end+1,1) = ip;
        end

        % If nothing was found, save an empty struct and move on
        if isempty(trough_idx)
            Spikes{thisChn, tr} = struct('idx', [], 't_ms', [], 'wf', zeros(0, preS+postS+1));
        else
            [trough_idx, order] = sort(trough_idx(:)); 
            vals = signal(trough_idx);
            keep_idx = true(size(trough_idx));
            last_keep = 1;

            for m = 2:numel(trough_idx)
                if trough_idx(m) - trough_idx(last_keep) <= refrac_samp
                    if vals(m) < vals(last_keep)
                         keep_idx(last_keep) = false;  % drop previous
                         last_keep = m; 
                    else 
                        keep_idx(m) = false;
                    end
                else 
                    last_keep = m;
                end
            end
            trough_idx = trough_idx(keep_idx);
        % --- Extract waveforms centered ---
            wf_list  = {};
            idx_list = [];

            for k = 1:numel(trough_idx)
                ip = trough_idx(k);
        
                % Ensure full window fits inside trial
                if ip - preS < 1 || ip + postS > total_sample
                    continue
                end        
                wf = signal(ip - preS : ip + postS);
        
                % Hard artifact rejection (absolute µV limit)
                if max(abs(wf)) > 500
                    continue
                end
        
                % Spike-shape features
                [~, i_trough] = min(wf);          % trough location within wf
                % restrict positive peak search to AFTER trough to enforce biphasic sequence
                [~, i_peak_rel] = max(wf(i_trough:end));
                i_peak = i_trough + i_peak_rel - 1;
        
                % (1) biphasic: negative first, then positive later
                is_biphasic = (min(wf) < 0) && (max(wf(i_trough:end)) > 0);
                % if ~is_biphasic, continue; end
        
                % (2) peak-to-peak amplitude (µV)
                pp_amp = max(wf) - min(wf);
                % if ~(pp_amp >= pp_min_uV && pp_amp <= pp_max_uV), continue; end
        
                % (3) trough-to-peak lag (samples)
                tp_width = i_peak - i_trough;
                % if ~(tp_width >= tp_min_samp && tp_width <= tp_max_samp), continue; end
        
                % Keep this spike
                wf_list{end+1}  = wf;              
                idx_list(end+1) = ip; 
            end
        
            % Save into Spikes{ch,tr}
            if ~isempty(idx_list)
                wf_mat = cat(1, wf_list{:});
            else
                wf_mat = zeros(0, preS+postS+1);
            end
            spk_times_ms = (idx_list - stim_idx) / FS * 1000;
        
            Spikes{thisChn, tr} = struct( ...
                'idx', idx_list, ...
                't_ms', spk_times_ms, ...
                'wf',  wf_mat ...
            );
        end
    end
end

% fclose(vFID);

% ------- Average spike waveform per amplitude per channel ------ %
wf_len = preS + postS + 1;
AvgWF_amp = nan(nChn, wf_len, n_AMP);   % mean waveform per (ch, amp)
Nspk_amp  = zeros(nChn, n_AMP);         % num of spikes in each mean

for ch = 1:nChn
    for k = 1:n_AMP
        rows = [];
        for tr = 1:nTrials
            if ampIdx(tr) ~= k, continue; end
            if isempty(Spikes{ch,tr}), continue; end
            wfs = Spikes{ch,tr}.wf;   % [Nspk x wf_len]
            if ~isempty(wfs)
                rows = [rows; wfs];
            end
        end
        if ~isempty(rows)
            AvgWF_amp(ch,:,k) = mean(rows,1,'omitnan');
            Nspk_amp(ch,k)    = size(rows,1);
        end
    end
end

% ------ Stimulation Channels ------ %
E_NAME = E_MAP(2:end); % Channel name e.g. "A-001"
if exist('StimParams','var') && ~isempty(StimParams)
    % Extract first column for channel 
    stimNames = StimParams(2:end,1);
    [tf, idx_all] = ismember(stimNames, E_NAME);
end
% group channels per trial
stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v  = unique(idx_all(rr));                                  % dedupe within trial
    v  = v(v>=0).';                                             % drop zeros, row vector
    stimChPerTrial_all{t} = v;
end
% find unique channels 
uniqueCh = unique([stimChPerTrial_all{:}]); 
% find unique channel set for each trial
comb = zeros(n_Trials, simultaneous_stim);
for i = 1:n_Trials
    v = stimChPerTrial_all{i};  
    comb(i,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb, 'rows');
combClass_win = combClass(Trial_start : Trial_start + nTrials - 1);

%% ================================== %%
%%                Plot                %%
%% ================================== %%
% ------ (1) Waveform plots ----- %
QuickAnalysis.plot_all_spike_waveforms(preS, postS, FS, nChn, nTrials, Spikes, uniqueCh);
QuickAnalysis.plot_mean_spike_waveforms_per_amp_per_stim_set(preS, postS, FS, nChn, n_AMP, Amps, Spikes, ampIdx, combClass_win, uniqueComb);
ch_spike = 22;
Spike_window = 5; % time window for the spike subplot (ms)
QuickAnalysis.plot_spikes_by_time_grid(ch_spike, Spikes, FS, preS, postS, Amps, ampIdx, nTrials, Spike_window, 20)

ch_filter_plot = 22;  % change to your desired channel
filtered_window = [-5,15];
QuickAnalysis.plot_filtered_signal_by_amplitude(MUA_all, ampIdx, Amps, t_axis, ch_filter_plot, filtered_window);


% ------ (2) Raster & PSTH plots  ------ %
% QuickAnalysis.raster_across_all_sets(nChn, nTrials, Spikes, Trial_start);
% QuickAnalysis.all_channels_psth_different_amp_stim_set(nChn, n_AMP, Amps, ampIdx, Spikes, combClass_win, uniqueComb);
% ch = 11; 
% QuickAnalysis.raster_different_stimset_and_amp(ch, n_AMP, Amps, ampIdx, Spikes, combClass_win, uniqueComb);
% QuickAnalysis.raster_psth_all_chn_per_stimset(nChn, Spikes, combClass_win, uniqueComb);

% QuickAnalysis.raster_psth_stacked_per_set_amp(nChn, Amps, ampIdx, Spikes, combClass_win, uniqueComb);
QuickAnalysis.raster_psth_stacked_per_set_amp_same_color(nChn, Amps, ampIdx, Spikes, combClass_win, uniqueComb);

% ------ (3) Heatmap ------ %
% QuickAnalysis.heatmap_fr_per_chn_per_amp_per_stimset(nChn, n_AMP, Amps, ampIdx, Spikes, combClass_win, uniqueComb, E_NAME);
QuickAnalysis.heatmap_per_chn_per_amp_per_stimset_same_color(nChn, n_AMP, Amps, ampIdx, Spikes, combClass_win, uniqueComb, E_NAME)
% ------ (4) Tuning curve ------ %
ch_tune = 19; resp_ms = [0 20];
% QuickAnalysis.tuning_fr_vs_amp_different_sets(ch_tune, resp_ms, Amps, n_AMP, ampIdx, Spikes, combClass_win, uniqueComb);
QuickAnalysis.tuning_fr_amp_different_sets(ch_tune, resp_ms, Amps, n_AMP, ampIdx, Spikes, combClass_win, uniqueComb);