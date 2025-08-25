%% Quick Analysis Code
clear all
close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Sabrina/Experimental_Design'));


%% ———————————————————————————————————————————————— %%
%% --------------- Data Processing ---------------- %%
%% ———————————————————————————————————————————————— %%

%% ------------------ Parameters ------------------ %%
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
blank_ms = 2;%[-1,1]; % artifact blanking window ±1 ms
blank_samps = round(blank_ms/1000 * FS);

% Trigger and parameters
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);
trig_all = trig; 
d = Depth_s(0);
% Extract trigger
nTrials = 1000; % number of trials used for analysis
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
postS = round(2.0e-3 * FS);

% For better alignment, local negative peak within ±0.5 ms of crossing
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
    
    % Extract the first 20 s for this channel
    signal = data_thres_uV(thisChn, :);
    
    % Apply band-pass filter
    signal_filt = filtfilt(bpFilt, double(signal));
    
    sigma = median(abs(signal_filt)) / 0.6745;
    % Compute threshold as -4.5 × standard deviation
    thr(ch_idx) = -4.5 * sigma;
end

%% ------ Baseline firing rate ------ %% 
maxBlankSamps = floor(10 * FS);  % up to 10s
firstTrig = trig_all(1);
blankSamps    = min(maxBlankSamps, max(0, firstTrig - 1)); % data for calculation

if blankSamps < 1
    warning('No pre-stim blank data available before the first trigger.');
else 
    fseek(vFID, 0, 'bof');
    blankBlock_uV = fread(vFID, [nChn, blankSamps], 'int16') * 0.195;  % [chn x samp] µV

    refrac_samp = round(1e-3*FS);
    for ch = 1:nChn
        thisChn = d(ch);
        x = double(blankBlock_uV(thisChn,:));
        x = filtfilt(bpFilt, x);                % filter for spike analysis
        thr_ch = thr(ch_idx);


    end



end


%% ------------------ Data Processing ------------------ %%

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
        % Interpolate artifact
        raw_interp = raw;
        raw_interp(blank_range) = interp1(keep_range, raw(keep_range), blank_range, 'linear', 'extrap');
        % raw(blank_range(1):blank_range(end)) = interpolate(raw(blank_range(1):blank_range(end)), 1);
        % Filter & rectify
        filtered = filtfilt(bpFilt, raw_interp);
        smoothed = movmean(filtered, 15);
        % Store
        MUA_all(ch, :, tr) = smoothed;

        %% ---- Baseline firing rate ---- %%
        % bl_start = stim_idx + bl_samp(1);
        % bl_end   = stim_idx + bl_samp(2) - 1; 

        signal = filtered;     % data extract

        %% ---- Spike detection ---- %%
        cross_idx = find(signal(2:end) < thr(ch) & signal(1:end-1) >= thr(ch)) + 1;
        cross_idx = cross_idx(:); 
        % Refractory period
        refrac_samp = max(1, round(1e-3 * FS));
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

            % ------ Addition Criteria ------ %
            [~, i_trough] = min(wf);                     % trough index
            [~, i_peak] = max(wf);  
            % 1) biphasic check 
            is_biphasic = (min(wf) < 0) && (max(wf(i_trough:end)) > 0);
            if ~is_biphasic, continue; end
            % % 2) peak-to-peak amplitude constraint
            pp_amp = max(wf) - min(wf);
            pp_ok  = (pp_amp >= pp_min_uV) && (pp_amp <= pp_max_uV);
            if ~pp_ok, continue; end
            % 3) through -> peak width constraint
            tp_width = i_peak - i_trough; % samples
            tp_ok = tp_width >= tp_min_samp && tp_width <= tp_max_samp;
            if ~tp_ok, continue; end
          

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





%% —————————————————————————————————————————— %%
%% ****************** Plot ****************** %% 
%% —————————————————————————————————————————— %%
%% Spike Waveform Plot 
wf_t  = (-preS:postS) / FS * 1000;   % waveform time ms

% ------ (1) Spikes waveform plot all ------ %
% Find global y-limit
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
cmap = lines(nTrials);   % nTrials distinct color

cmap_seq = turbo(nTrials);                     % deterministic, perceptually uniform
trialColor = @(tr) cmap_seq(tr,:); 
% Plot
figure;
for ch = 1:nChn
    subplot(4, 8, ch);
    hold on;

    % Soft red background if this is a stimulation channel
    if ismember(ch, uniqueCh)
        set(gca, 'Color', [1 0.95 0.95]);    % light pink
        baseLW = 1.8;                       % slightly thicker for stim chans
    else
        baseLW = 1.2;
    end

    for tr = 1:nTrials
        if ~isempty(Spikes{ch,tr}) && isfield(Spikes{ch,tr},'wf')
            wfs = Spikes{ch,tr}.wf;
            if ~isempty(wfs)
                % plot(wf_t, wfs', 'Color',cmap(tr,:),'LineWidth', 0.5);   % plot each spike
                plot(wf_t, wfs', 'Color',trialColor(tr),'LineWidth', 0.5);   % plot each spike
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
sgtitle(sprintf('Spike waveforms across %d trials', nTrials),'FontWeight','bold');


% ------ (2) Avg Spike waveform per AMP per CHH ------ % 
cmap = turbo(n_AMP);
% y-axis
abs_max = max(abs(AvgWF_amp), [], 'all', 'omitnan');
y_lim = ceil(abs_max/5)*5;

figure('Name','Spike avg waveform','NumberTitle','off','Color','w');
hLines = gobjects(1,n_AMP);

for ch = 1:nChn
     subplot(4,8,ch); hold on
        
     % Soft red background if this is a stimulation channel
    if ismember(ch, uniqueCh)
        set(gca, 'Color', [1 0.95 0.95]);    % light pink
        baseLW = 1.8;                       % slightly thicker for stim chans
    else
        baseLW = 1.2;
    end

     for k = 1:n_AMP
        if ~all(isnan(AvgWF_amp(ch,:,k)))
            hLines(k) = plot(wf_t, squeeze(AvgWF_amp(ch,:,k)), ...
                'LineWidth', 1.5, 'Color', cmap(k,:));
        end
     end

    xline(0,'r-');
    title(sprintf('Ch %d', ch), 'FontSize', 8);
    xlim([wf_t(1) wf_t(end)]);
    ylim([-y_lim y_lim]);    
    xlabel('Time (ms)','FontSize',8);
    ylabel('\muV','FontSize',8);
    box off
end
sgtitle(sprintf('Average Spike waveforms across Amplitudes (%d trials)', nTrials),'FontWeight','bold');
% legend
legLabels = arrayfun(@(a) sprintf('%g \\muA', a), Amps, 'UniformOutput', false);
lgd = legend(hLines, legLabels, 'Orientation','horizontal', ...
             'Box','off','FontSize',10);
lgd.Position = [0.1, 0.01, 0.8, 0.05];   % center at bottom


% ------ (3) Spike waveform per AMP one CHN ------ % 
ch_view = 8;
absmax = 0;
for tr = 1:nTrials
    S = Spikes{ch_view,tr};
    if ~isempty(S) && isfield(S,'wf') && ~isempty(S.wf)
        absmax = max(absmax, max(abs(S.wf),[],'all'));
    end
end
y_lim = ceil(absmax/5)*5;

nRows = ceil(sqrt(n_AMP));
nCols = ceil(n_AMP / nRows);

figure('Name',sprintf('Ch %d: spikes by amplitude', ch_view), ...
       'NumberTitle','off','Color','w');

for k = 1:n_AMP
    subplot(nRows, nCols, k); hold on
    rows = [];
    for tr = 1:nTrials
        if ampIdx(tr) ~= k, continue; end
        S = Spikes{ch_view,tr};
        if ~isempty(S) && isfield(S,'wf') && ~isempty(S.wf)
            rows = [rows; S.wf]; % [Nspk x wf_len]
        end
    end
    if ~isempty(rows)
        nSpk = size(rows,1);
        rowColors = lines(nSpk);
        for r = 1: nSpk
            plot(wf_t, rows(r,:), 'Color', rowColors(r,:), 'LineWidth', 0.5);
        end
        % plot all spikes (thin, faint)
        % plot(wf_t, rows', 'Color', 'k', 'LineWidth', 0.5); % alpha supported R2022a+; if not, drop the 0.25
        % mean waveform (bold)
        plot(wf_t, mean(rows,1,'omitnan'), 'Color', 'b', 'LineWidth', 2.0);
        nSpk = size(rows,1);
    else
        nSpk = 0;
    end
    xline(0,'r-');
    title(sprintf('%g \\muA  (n=%d)', Amps(k), nSpk), 'FontSize', 10);
    xlim([wf_t(1) wf_t(end)]);
    ylim([-y_lim y_lim]);
    xlabel('Time (ms)'); ylabel('\muV');
    axis square 
    box off
end
sgtitle(sprintf('Channel %d: spikes grouped by amplitude', ch_view), 'FontWeight','bold');

% ------ (4) Ave Spike per Amp one Chn ------ % 
ch_view = 23;   % <--- change this to whichever channel you want
cmap = turbo(n_AMP);

abs_max = max(abs(AvgWF_amp(ch_view,:,:)), [], 'all', 'omitnan');
y_lim   = ceil(abs_max/5)*5;

figure('Name',sprintf('Ch %d: Avg spikes by amplitude', ch_view), ...
       'NumberTitle','off','Color','w'); 
hold on

hLines = gobjects(1,n_AMP);

for k = 1:n_AMP
    y = squeeze(AvgWF_amp(ch_view,:,k));  % [1 x wf_len]
    if ~all(isnan(y))
        hLines(k) = plot(wf_t, y, 'LineWidth', 1.8, 'Color', cmap(k,:));
    end
end

% formatting
xline(0,'r-','LineWidth',2);
title(sprintf('Channel %d: Avg spike waveforms', ch_view), 'FontSize', 12);
xlim([wf_t(1) wf_t(end)]);
ylim([-y_lim y_lim]);
xlabel('Time (ms)');
ylabel('\muV');
box off
axis square

% legend
legLabels = arrayfun(@(a) sprintf('%g \\muA', a), Amps, 'UniformOutput', false);
valid = isgraphics(hLines);

lgd = legend(hLines(valid), legLabels(valid), ...
       'Orientation','horizontal','Box','off','FontSize',10, ...
       'NumColumns', ceil(numel(Amps)/2), ...
       'Location','southoutside');

% ------- (5) Spike Waveforms by Stim Set ------ % 
classes_here = unique(combClass_win(:)');   % combos appear in the window
for cc = 1:length(classes_here)
    tr_idx = find(combClass_win == cc); % trials belong to this combo
    stimIdx = uniqueComb(cc, :);

    % global y-limit
    ymax_global = 0;
    for ch = 1:nChn
        for t = 1:length(tr_idx)
            S = Spikes{ch,tr_idx(t)};
            if ~isempty(S) && isfield(S,'wf') && ~isempty(S.wf)
                ymax_global = max(ymax_global, max(abs(S.wf),[],'all'));
            end
        end
    end
    y_lim = ceil(ymax_global/50)*50;  % round up to nearest 50 µV
    comboLabel = strjoin(E_NAME(stimIdx), ' + ');
    cmap = jet(numel(tr_idx));

    figure('Color','w','Name',sprintf('Spikes | %s', comboLabel),'NumberTitle','off');
    
    for ch = 1:nChn
        subplot(4,8,ch); hold on;

        % stimulation channels
        if ismember(ch, stimIdx)
            set(gca, 'Color', [1 0.95 0.95]);   % light pink
        end

        nSpk_ch = 0;
        for t = 1:length(tr_idx)
            S = Spikes{ch,tr_idx(t)};
            if ~isempty(S) && isfield(S,'wf') && ~isempty(S.wf)
                plot(wf_t, S.wf','Color', cmap(t,:), 'LineWidth', 0.5);
                nSpk_ch = nSpk_ch + size(S.wf, 1);
            end
        end

        xline(0,'r-','LineWidth',1.2);
        title(sprintf('Ch %d  (n=%d )', ch, nSpk_ch), 'FontSize', 10);
        xlim([wf_t(1) wf_t(end)]);
        ylim([-y_lim y_lim]);   % same across all subplots
        xlabel('Time (ms)','FontSize',8);
        ylabel('\muV','FontSize',8);
        box off;
    end
    sgtitle(sprintf('Spike waveforms | Stim combo: %s  | Trials: %d', comboLabel, numel(tr_idx)), ...
            'FontWeight','bold');    
end

% ------ (6) Mean Spike waveforms per AMP per Stim Set ------ %
classes_here = unique(combClass_win(:)'); % combos appear
for cc = 1:length(classes_here)
    tr_idx  = find(combClass_win == cc); % trials of this combo
    stimIdx = uniqueComb(cc,:); 
    
    % Compute mean waveform per CHN per AMP
    AvgWF_amp_combo = nan(nChn, wf_len, n_AMP);
    for ch = 1:nChn
        for k = 1:n_AMP
            rows = [];
            for tt = 1:length(tr_idx(:)')
                if ampIdx(tt) ~= k, continue; end
                S = Spikes{ch, tr_idx(tt)};
                if ~isempty(S) && isfield(S,'wf') && ~isempty(S.wf)
                    rows = [rows; S.wf];
                end
            end
            if ~isempty(rows)
                AvgWF_amp_combo(ch,:,k) = mean(rows, 1, 'omitnan');
            end
        end
    end

    % y-limit
    abs_max = max(abs(AvgWF_amp_combo), [], 'all', 'omitnan');
    y_lim   = 50 * ceil(abs_max/50);

    cmap = turbo(n_AMP);
    % plot 
    comboLabel = strjoin(E_NAME(stimIdx), ' + ');
    figure('Color','w','Name',sprintf('Mean spikes waveforms| %s', comboLabel),'NumberTitle','off');
    for ch = 1:nChn
        subplot(4,8,ch); hold on
        if ismember(ch, stimIdx) % highlight stimulated channels
            set(gca,'Color',[1 0.95 0.95]);
        end

        for k = 1:n_AMP
            y = AvgWF_amp_combo(ch,:,k);
            if ~all(isnan(y))
                hL(k) = plot(wf_t, squeeze(y), 'LineWidth', 1.6, 'Color', cmap(k,:));
            end
        end

        xline(0,'r-');
        title(sprintf('Ch %d', ch), 'FontSize', 8);
        xlim([wf_t(1) wf_t(end)]);
        ylim([-y_lim y_lim]);
        xlabel('Time (ms)','FontSize',8);
        ylabel('\muV','FontSize',8);
        box off
    end

    legLabels = arrayfun(@(a) sprintf('%g \\muA', a), Amps, 'UniformOutput', false);
    valid = isgraphics(hL);   % some amplitudes may be empty for this combo
    if any(valid)
        lgd = legend(hL(valid), legLabels(valid), 'Orientation','horizontal', ...
                     'Box','off','FontSize',10);
        lgd.Position = [0.1, 0.01, 0.8, 0.05];
    end
    sgtitle(sprintf('Mean spike waveforms by amplitude | Stim combo: %s', comboLabel), ...
            'FontWeight','bold');
end



%% Raster Plot 
ch_raster = 23;     % channel for raster plot
ras_win   = [-20 100];% ms window for raster display
bin_ms    = 1;      % PSTH bin (ms)

% Collect spike times
spkTimesPerTrial = cell(nTrials,1);
for tr = 1:nTrials
    if ~isempty(Spikes{ch_raster,tr}) && isfield(Spikes{ch_raster,tr},'t_ms')
        tt = Spikes{ch_raster,tr}.t_ms;
        spkTimesPerTrial{tr} = tt(tt >= ras_win(1) & tt <= ras_win(2));
    else
        spkTimesPerTrial{tr} = [];
    end
end

% Figure with raster (top & larger) and PSTH (bottom & smaller)
figure();
if ~isempty(bin_ms) && bin_ms > 0
    tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile([3 1]);
else
    tl = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile;
end
title(tl, sprintf('Raster | Ch %d  (%d trials)', ch_raster, nTrials), 'FontWeight','bold');
% ax1 = nexttile; 
hold(ax1,'on');
for tr = 1:nTrials
    tt = spkTimesPerTrial{tr};
    if ~isempty(tt)
        % draw a vertical tick per spike
        y0 = Trial_start + tr - 1;
        for k = 1:numel(tt)
            plot(ax1, [tt(k)+1 tt(k)+1], [y0-1 y0+1], ...
                 'Color', 'k', 'LineWidth', 1.5);
        end
    end
end
xline(0, 'r', 'LineWidth', 1.5);
xlim(ax1, ras_win);
ylim(ax1, [Trial_start-0.5, Trial_start+nTrials-0.5]);
yt = round(linspace(Trial_start, Trial_start + nTrials - 1, 5));
xt = round(linspace(ras_win(1),ras_win(2),25));
yticks(ax1, yt);
xticks(ax1, xt);
ylabel(ax1, 'Trial #');
xlabel(ax1, 'Time (ms)');
box(ax1,'off');

% PSTH
if ~isempty(bin_ms) && bin_ms > 0
    ax2 = nexttile; hold(ax2,'on');

    edges  = ras_win(1):bin_ms:ras_win(2);           % bin edges (ms)
    counts = zeros(1, numel(edges)-1);               % init histogram
    ctrs = edges(1:end-1) + diff(edges)/2;

    % count spikes
    for tr = 1:nTrials
        tt = spkTimesPerTrial{tr};                   % spike times (ms) for this trial
        if isempty(tt), continue; end
        counts = counts + histcounts(tt, edges);
    end

    % convert to firing rate (spike/s)
    bin_s = bin_ms/1000;
    rate  = counts / (nTrials * bin_s);
    % smooth data
    smooth_ms = 3; % smooth window
    smooth_bins = max(1, round(smooth_ms/bin_ms));
    % rate_s = rate; % preallocate
    sigma_bins  = max(1, round(smooth_bins/2));

    k   = 0:(smooth_bins-1);                               % ONLY past & current samples
    g   = exp(-0.5*(k./sigma_bins).^2);                    % one-sided Gaussian
    g   = g ./ sum(g);                                     % normalize to sum=1

    % Causal smoothing (no leakage before 0): y[n] = sum g[k]*x[n-k]
    rate_s = filter(g, 1, rate);                           % FIR, causal

    % rate_s = smoothdata(rate, 'movmean', [smooth_bins-1 0]);

    % rate_s      = smoothdata(rate, 'gaussian', smooth_bins);
    plot(ax2, ctrs, rate_s, 'k-', 'LineWidth', 1.8);
    % bar(ax2, ctrs, rate, 1.0, 'FaceColor',[0.2 0.2 0.2], 'EdgeColor','none');

    xline(ax2, 0, 'r', 'LineWidth', 1.5);
    xlim(ax2, ras_win);
    % ylim(ax2, [0, 300]);
    ylabel(ax2, 'Rate (sp/s)');
    xlabel(ax2, 'Time (ms)');
    box(ax2,'off');
end

%% ------ Normalized Spike Analysis ------ %% 
 
% Parameters
resp_ms = [0 20]; % post stimulation window (ms)
