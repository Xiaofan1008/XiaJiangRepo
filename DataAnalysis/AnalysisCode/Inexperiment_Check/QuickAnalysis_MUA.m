%% Quick Analysis Code
clear all
close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Sabrina/Experimental_Design'));



%% Load Intan parameters
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

dName='amplifier';
% info = dir([filepath filesep dName '.dat']);
% info = info.bytes/2;
% nL = (ceil(info / (nChn*FS*double(T)))+1);
vFID = fopen([filepath filesep dName '.dat'],'r'); % read data file
%% Quick Analysis parameters
nTrials = 400; % number of trials used for analysis
window_ms = [-100,100]; % time window for data extract
window_samps = round(window_ms/1000 * FS);
total_sample = diff(window_samps) + 1;
trialLength = diff(window_samps)+1;
blank_ms = 2;%[-1,1]; % artifact blanking window ±1 ms
blank_samps = round(blank_ms/1000 * FS);

%% Load trigger and parameters
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);
trig = trig(1:nTrials);
nTrig = length(trig); % length of trigger
d = Depth;
% load trial parameters
TrialParams = loadTrialParams;
trialIDs = cell2mat(TrialParams(:,2));
trialIDs = trialIDs(1:nTrials);
%% Bandpass filter design for MUA (300-3000 Hz)
centerFreq = 4000;   % Hz
halfBandwidth = 3700; % Hz
f1 = centerFreq - halfBandwidth;  % lower cutoff
f2 = centerFreq + halfBandwidth;  % upper cutoff
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency1', f1, ...
    'HalfPowerFrequency2', f2, ...
    'SampleRate', FS);

% Time vector for plotting
t_axis = (window_samps(1):window_samps(2)) / FS * 1000;

% Preallocate
MUA_all = zeros(nChn, trialLength, nTrials);
Spikes = cell(nChn, nTrials); 

% baseline firing rate
bl_ms   = [-60, -5]; % baseline firing rate time 
bl_samp = round(bl_ms/1000 * FS);

% Spike waveform window ~3 ms (e.g., 1 ms pre, 2 ms post)

preS  = round(1.0e-3 * FS);
postS = round(2.0e-3 * FS);

% For better alignment, local negative peak within ±0.5 ms of crossing
alignRad = max(1, round(0.5e-3 * FS));

for tr = 1:nTrials
    fprintf('Reading trial %d/%d\n', tr, nTrials);
    start_index = trig(tr) + window_samps(1);
    start_byte = start_index * nChn * 2; % 2 bytes per int16
    fseek(vFID, start_byte, 'bof');
    data = fread(vFID, [nChn, total_sample], 'int16') * 0.195; % in uV

    stim_idx = -window_samps(1) + 1;
    blank_range = max(1, stim_idx - blank_samps) : min(total_sample, stim_idx + blank_samps);
    keep_range = setdiff(1:total_sample, blank_range);
    
    for ch = 1:nChn
        
        thisChn = d(ch);
        raw = data(thisChn, :);
        %% ---- Artifact Blank ---- %%
        % Interpolate blank
        raw_interp = raw;
        raw_interp(blank_range) = interp1(keep_range, raw(keep_range), blank_range, 'linear', 'extrap');
        % raw(blank_range(1):blank_range(end)) = interpolate(raw(blank_range(1):blank_range(end)), 1);
        % Filter & rectify
        filtered = filtfilt(bpFilt, raw_interp);
        % mua = abs(filtered);
        % smoothed = movmean(mua, 15);
        smoothed = movmean(filtered, 15);
        % Store
        MUA_all(ch, :, tr) = smoothed;
        MUA_unsmooth(ch, :, tr) = filtered;

        %% ---- baseline firing rate ---- %%
        bl_start = stim_idx + bl_samp(1);
        bl_end   = stim_idx + bl_samp(2) - 1; 

        signal = filtered;     % data extract
        seg = signal(bl_start:bl_end);  % data extract

        % noise estimate and negative threshold
        sigma = median(abs(seg)) / 0.6745;        % MAD -> sigma
        thr   = -4.5 * sigma;                     % threshold

        %% ---- Spike detection ----
        cross_idx = find(signal(2:end) < thr & signal(1:end-1) >= thr) + 1;
        cross_idx = cross_idx(:); 
        % Refractory (e.g., 0.6 ms)
        refrac_samp = max(1, round(0.6e-3 * FS));
        if numel(cross_idx) > 1
            keep = [true; diff(cross_idx) > refrac_samp];   % (n x 1) logical
        else
            keep = true(size(cross_idx));                   % 0 or 1 crossing
        end        
        cross_idx = cross_idx(keep);

        peaks = [];  % sample indices at peak

        if ~isempty(cross_idx)
            for k = 1:numel(cross_idx)
                i0 = cross_idx(k);
                a  = max(1, i0 - alignRad);
                b  = min(total_sample, i0 + alignRad);
                [~, rel] = min(signal(a:b));     % negative peak
                ip = a + rel - 1;           % absolute index in this trial window
                % Check waveform bounds
                if ip - preS < 1 || ip + postS > total_sample
                    continue
                end
                peaks(end+1) = ip;
            end
        end

        % Extract waveforms; reject artifacts with |wf| > 500 µV
        wf_list = {}; % waveform 
        idx_list = [];
        for k = 1:numel(peaks)
            ip = peaks(k);
            wf = signal(ip - preS : ip + postS);
            if max(abs(wf)) > 500      % artifact rejection
                continue
            end
            wf_list{end+1} = wf; 
            idx_list(end+1) = ip; 
        end
        % Save into Spikes{ch,tr}
        if ~isempty(idx_list)
            wf_mat = cat(1, wf_list{:});              % [Nspk x (preS+postS+1)]
        else
            wf_mat = zeros(0, preS+postS+1);
        end
        spk_times_ms = (idx_list - stim_idx) / FS * 1000;  % relative to stim (ms)
        
        Spikes{thisChn, tr} = struct( ...
            'idx', idx_list, ...
            't_ms', spk_times_ms, ...
            'wf',  wf_mat ...
        );           
    end
end

fclose(vFID);


%% plot all trials 
for ch = 1:nChn
    subplot(4, 8, ch);
    MUA_ch = MUA_all(ch,:,:);
    hold on;
    ymax = max(abs(MUA_ch),[],"all");
    for tr = 1:nTrials
        plot(t_axis, MUA_all(ch, :,tr),'k');
    end
    hold off;
    title(sprintf('Ch %d', ch));
    xlabel('Time (ms)');
    ylabel('µV');
    ylim([-ymax,ymax]);
    xlim([-2,10])
end
sgtitle('Spike waveform across trials');


wf_t  = (-preS:postS) / FS * 1000;   % ms

% --- Find global y-limit ---
ymax_global = 0;
for ch = 1:nChn
    for tr = 1:nTrials
        if ~isempty(Spikes{ch,tr}) && isfield(Spikes{ch,tr},'wf')
            wfs = Spikes{ch,tr}.wf;
            if ~isempty(wfs)
                ymax_global = max(ymax_global, max(abs(wfs),[],'all'));
            end
        end
    end
end

% --- Plot ---
figure;
for ch = 1:nChn
    subplot(4, 8, ch);
    hold on;

    for tr = 1:nTrials
        if ~isempty(Spikes{ch,tr}) && isfield(Spikes{ch,tr},'wf')
            wfs = Spikes{ch,tr}.wf;
            if ~isempty(wfs)
                plot(wf_t, wfs', 'k');   % plot each spike
            end
        end
    end

    title(sprintf('Ch %d', ch));
    xlabel('Time (ms)');
    ylabel('µV');
    ylim([-ymax_global, ymax_global]);   % same for all
    xlim([wf_t(1), wf_t(end)]);
    hold off;
end
sgtitle('Spike waveforms across all trials');
