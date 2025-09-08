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

% ---- Artifact blank parameters ---- %
preBlank_ms   = 1;   % ms to blank BEFORE stim
postBlank_ms  = 2;   % ms to blank AFTER  stim
preBlank_samp  = round(preBlank_ms/1000 * FS);
postBlank_samp = round(postBlank_ms/1000 * FS);

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
d = Depth_s(0); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex
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
    thr(ch_idx) = -4.5 * sigma; % Compute threshold as -4.5 × standard deviation    
end

%% ------------------ Data Processing ------------------ %%
for tr = 1:nTrials
    fprintf('Reading trial %d/%d\n', tr, nTrials);
    start_index = trig(tr) + window_samps(1);
    start_byte = start_index * nChn * 2; % 2 bytes per int16
    fseek(vFID, start_byte, 'bof');
    data = fread(vFID, [nChn, total_sample], 'int16') * 0.195; % in uV

    stim_idx = -window_samps(1) + 1;
    blank_range = max(1, stim_idx - preBlank_samp) : min(total_sample, stim_idx + postBlank_samp);
    keep_range = setdiff(1:total_sample, blank_range);
    
    for ch = 1:nChn       
        thisChn = d(ch);
        raw = data(thisChn, :);
        %% ---- Artifact Blank ---- %%
        stimSample  = -window_samps(1) + 1;  % stim index in the trial window
        segStartIdx = max(1,           stimSample - preBlank_samp);
        segEndIdx   = min(total_sample, stimSample + postBlank_samp);

        % Interpolate artifact
        raw_interp = raw;
        raw_interp(blank_range) = interp1(keep_range, raw(keep_range), blank_range, 'linear', 'extrap');
        % Filter & rectify
        sig_filt_plot = filter(bpFIR, 1, raw_interp);
        if gd > 0
            sig_filt_plot = [sig_filt_plot(gd+1:end), zeros(1, gd)];
        end
        smoothed = movmean(sig_filt_plot, 15);   % if you still want this quick display
        % Store
        MUA_all(ch, :, tr) = sig_filt_plot;

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
        guardStartIdx = max(1,segStartIdx - guard_pre_samp);
        guardEndIdx   = min(total_sample, segEndIdx  + guard_post_samp);
        mask_keep = ~(cross_idx >= guardStartIdx & cross_idx <= guardEndIdx);
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
                if ~is_biphasic, continue; end
        
                % (2) peak-to-peak amplitude (µV)
                pp_amp = max(wf) - min(wf);
                if ~(pp_amp >= pp_min_uV && pp_amp <= pp_max_uV), continue; end
        
                % (3) trough-to-peak lag (samples)
                tp_width = i_peak - i_trough;
                if ~(tp_width >= tp_min_samp && tp_width <= tp_max_samp), continue; end
        
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


%% ==================== Plot ==================== %% 

%% -------- Spike Waveform Plot -------- %
wf_t  = (-preS:postS) / FS * 1000;   % waveform time ms

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
for i = 1:nChn
    % ch_view = 20;
    ch_view = i;
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
end

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
    % comboLabel = strjoin(E_NAME(stimIdx), ' + ');
    comboLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx(stimIdx>0), 'UniformOutput', false), ' + ');
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
            for tt = 1:numel(tr_idx)
                tr = tr_idx(tt);
                if ampIdx(tr) ~= k, continue; end
                S = Spikes{ch, tr};
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
    y_lim   = 10 * ceil(abs_max/10);

    cmap = turbo(n_AMP);
    % plot 
    % comboLabel = strjoin(E_NAME(stimIdx), ' + ');
     comboLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx(stimIdx>0), 'UniformOutput', false), ' + ');
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
        % ylim([-y_lim y_lim]);
        ylim([-200 200]);
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

% ------ (7) Mean Spike pre Amp per Stim Set one Chn ------ %
for cc = 1:nChn
    ch_view = cc;                                 % <-- channel to visualize
    classes_here = unique(combClass_win(:)');     % stimulation-set classes present
    nComb  = numel(classes_here);
    
    % Grid for the stim sets
    nRows = ceil(sqrt(nComb));
    nCols = ceil(nComb / nRows);
    
    % Precompute means and spike counts
    AvgWF_ch_amp_bySet = nan(wf_len, n_AMP, nComb);
    nSpk_set = zeros(1, nComb);
    
    for ii = 1:nComb
        cc = classes_here(ii);
        tr_idx = find(combClass_win == cc);
    
        for k = 1:n_AMP
            rows = [];
            for t = 1:numel(tr_idx)
                tr = tr_idx(t);
                if ampIdx(tr) ~= k, continue; end
                S = Spikes{ch_view, tr};
                if ~isempty(S) && isfield(S,'wf') && ~isempty(S.wf)
                    rows = [rows; S.wf]; %#ok<AGROW>
                    nSpk_set(ii) = nSpk_set(ii) + size(S.wf,1);
                end
            end
            if ~isempty(rows)
                AvgWF_ch_amp_bySet(:,k,ii) = mean(rows, 1, 'omitnan');
            end
        end
    end
    
    % Global y-limit
    abs_max = max(abs(AvgWF_ch_amp_bySet), [], 'all', 'omitnan');
    y_lim   = 10 * ceil(abs_max/10);
    
    % Plot
    cmap = turbo(n_AMP);
    figure('Color','w','Name',sprintf('Ch %d | Mean spikes per amplitude across stim sets', ch_view), ...
           'NumberTitle','off');
    tl = tiledlayout(nRows,nCols,'TileSpacing','compact');
    
    for ii = 1:nComb
        cc = classes_here(ii);
        stimIdx    = uniqueComb(cc,:);
        % comboLabel = strjoin(E_NAME(stimIdx), ' + ');
        comboLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx(stimIdx>0), 'UniformOutput', false), ' + ');
    
        ax = nexttile; hold(ax,'on');
        if ismember(ch_view, stimIdx)
            set(ax,'Color',[1 0.95 0.95]);
        end
    
        for k = 1:n_AMP
            y = AvgWF_ch_amp_bySet(:,k,ii);
            if ~all(isnan(y))
                plot(wf_t, y, 'LineWidth', 1.6, 'Color', cmap(k,:));
            end
        end
    
        xline(0,'r-');
        title(sprintf('%s  (n=%d)', comboLabel, nSpk_set(ii)), ...
              'FontSize', 9, 'Interpreter','none');
        xlim([wf_t(1) wf_t(end)]);
        ylim([-y_lim y_lim]);
        xlabel('Time (ms)'); ylabel('\muV');
        box off
    end
    
    % Detached legend (create dummy lines just for legend)
    legLabels = arrayfun(@(a) sprintf('%g \\muA', a), Amps, 'UniformOutput', false);
    dummyLines = gobjects(1,n_AMP);
    for k = 1:n_AMP
        dummyLines(k) = plot(nan, nan, 'Color', cmap(k,:), 'LineWidth', 1.6);
    end
    lgd = legend(dummyLines, legLabels, ...
        'Orientation','horizontal','Box','off','FontSize',10, ...
        'NumColumns', min(n_AMP, ceil(n_AMP/2)), ...
        'Location','southoutside');  % put below entire layout
    lgd.Layout.Tile = 'south';
    
    sgtitle(sprintf('Channel %d: Mean spike waveforms by amplitude across stim sets', ch_view), ...
            'FontWeight','bold');
end


%% -------- Raster Plot ---------%
% -------(1): across all sets ------- %
for cc = 1:nChn
    ch_raster = cc;     % channel for raster plot
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
                plot(ax1, [tt(k) tt(k)], [y0-1 y0+1], ...
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
        % Smoothing
        k   = 0:(smooth_bins-1);                               % ONLY past & current samples
        g   = exp(-0.5*(k./sigma_bins).^2);                    % one-sided Gaussian
        g   = g ./ sum(g);                                     % normalize to sum=1    
        % Causal smoothing (no leakage before 0): y[n] = sum g[k]*x[n-k]
        rate_s = filter(g, 1, rate);                           % FIR, causal
    
        % rate_s = smoothdata(rate, 'movmean', [smooth_bins-1 0]);
    
        % rate_s      = smoothdata(rate, 'gaussian', smooth_bins);
        
        % plot(ax2, ctrs, rate_s, 'k-', 'LineWidth', 1.8);
        % bar(ax2, ctrs, rate, 1.0, 'FaceColor',[0.2 0.2 0.2], 'EdgeColor','none');
           
        % y-lim logic
        maxVal = max(rate_s);
        if maxVal > 100
            plotColor = 'b';        % highlight in red
            plot(ax2, ctrs, rate_s, plotColor, 'LineWidth', 1.8);
            ylim(ax2,'auto');       % allow autoscale
        else
            plotColor = 'k';
            plot(ax2, ctrs, rate_s, plotColor, 'LineWidth', 1.8);
            ylim(ax2,[0 100]);      % fixed limit
        end


        xline(ax2, 0, 'r', 'LineWidth', 1.5);
        xlim(ax2, ras_win);
        % ylim(ax2, [0, 300]);
        ylabel(ax2, 'Rate (sp/s)');
        xlabel(ax2, 'Time (ms)');
        box(ax2,'off');
    end
end


% ------ (2) One channel across Stim Sets ------ %
 for cc = 1:nChn
    ch_raster   = cc;          % channel to visualize
    ras_win     = [-20 100];   % ms
    bin_ms      = 1;           % PSTH bin (ms)
    psth_thresh = 200;         % if peak > this => auto ylim + red curve, else [0,200] + black
    
    classes_here = unique(combClass_win(:)');   % stim sets present
    nSets = numel(classes_here);
    if nSets == 0
        warning('No stimulation sets found.');
    else
        % layout math
        nCols = ceil(sqrt(nSets));
        nRowsSets = ceil(nSets / nCols);
        nRows = 2 * nRowsSets;
    
        fig = figure('Color','w','Name',sprintf('Ch %d | Rasters by stimulation set', ch_raster), ...
                     'NumberTitle','off');
    
        axRaster = gobjects(1, nSets);
        axPSTH   = gobjects(1, nSets);
        tileIndex = @(r,c) (r-1)*nCols + c;
    
        for ii = 1:nSets
            cc     = classes_here(ii);
            tr_idx = find(combClass_win == cc);
            nThis  = numel(tr_idx);
            if nThis == 0, continue; end
    
            % Stim-set label
            stimIdx = uniqueComb(cc,:); stimIdx = stimIdx(stimIdx>0);
            setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');
    
            % grid placement
            col    = mod(ii-1, nCols) + 1;
            rowSet = floor((ii-1) / nCols) + 1;
            rRaster = 2*rowSet - 1;
            rPSTH   = rRaster + 1;
    
            axRaster(ii) = subplot(nRows, nCols, tileIndex(rRaster, col));
            axPSTH(ii)   = subplot(nRows, nCols, tileIndex(rPSTH,   col));
    
            % Collect spike times
            spkTimesPerTrial = cell(nThis,1);
            for k = 1:nThis
                tr = tr_idx(k);
                if ~isempty(Spikes{ch_raster,tr}) && isfield(Spikes{ch_raster,tr},'t_ms')
                    tt = Spikes{ch_raster,tr}.t_ms;
                    spkTimesPerTrial{k} = tt(tt >= ras_win(1) & tt <= ras_win(2));
                else
                    spkTimesPerTrial{k} = [];
                end
            end
    
            % Raster
            axes(axRaster(ii)); hold on;
             % Background highlight if this channel was stimulated
            if ismember(ch_raster, stimIdx)
                set(gca,'Color',[1 0.9 0.9]);  % light pink background
            end

            for k = 1:nThis
                tt = spkTimesPerTrial{k};
                if isempty(tt), continue; end
                y0 = k;
                for s = 1:numel(tt)
                    plot([tt(s) tt(s)], [y0-0.4 y0+0.4], 'k', 'LineWidth', 1.0);
                end
            end
            xline(0, 'r', 'LineWidth', 1.2);
            xlim(ras_win);
            ylim([0.5, nThis+0.5]);
            ylabel('Trial #');
            box off;
            title(sprintf('Set %d: %s (Trials: %d)', ii, setLabel, nThis), ...
                  'Interpreter','none','FontSize',9);
    
            % PSTH (no title)
            axes(axPSTH(ii)); hold on;
            edges = ras_win(1):bin_ms:ras_win(2);
            ctrs  = edges(1:end-1) + diff(edges)/2;
            counts = zeros(1, numel(edges)-1);
            for k = 1:nThis
                counts = counts + histcounts(spkTimesPerTrial{k}, edges);
            end
            bin_s = bin_ms/1000;
            rate  = counts / (nThis * bin_s);
    
            % smoothing
            smooth_ms = 3; smooth_bins = max(1, round(smooth_ms/bin_ms));
            sigma_bins = max(1, round(smooth_bins/2));
            kk = 0:(smooth_bins-1);
            g = exp(-0.5*(kk./sigma_bins).^2); g = g/sum(g);
            rate_s = filter(g,1,rate);
    
            peakVal = max(rate_s);

             % Background highlight if this channel was stimulated
            if ismember(ch_raster, stimIdx)
                set(gca,'Color',[1 0.9 0.9]);  % light pink background
            end

            if peakVal > psth_thresh
                plot(ctrs, rate_s, 'b-', 'LineWidth', 1.4);
                ylim('auto');
            else
                plot(ctrs, rate_s, 'k-', 'LineWidth', 1.4);
                ylim([0 psth_thresh]);
            end
            xline(0, 'r', 'LineWidth', 1.2);
            xlim(ras_win);
            ylabel('Rate (sp/s)'); xlabel('Time (ms)');
            box off;
        end
    
        % Resize raster vs PSTH
        for ii = 1:nSets
            pR = get(axRaster(ii),'Position');
            pP = get(axPSTH(ii),  'Position');
            left  = pR(1); width = pR(3);
            totalTop   = pR(2)+pR(4);
            bottomPair = pP(2);
            totalHeight = totalTop-bottomPair;
            hRaster = 0.72*totalHeight;
            hGap    = 0.06*totalHeight;
            hPSTH   = 0.22*totalHeight;
            bottomPSTH   = bottomPair;
            bottomRaster = bottomPSTH+hPSTH+hGap;
            set(axRaster(ii),'Position',[left bottomRaster width hRaster]);
            set(axPSTH(ii),  'Position',[left bottomPSTH  width hPSTH ]);
        end
    
        sgtitle(sprintf('Ch %d', ch_raster),'FontWeight','bold');
    end
 end

 % ------ (3) Raster & PSTH shared plot ------ %%
ch_raster   = 32;          % channel to visualize
ras_win     = [-20 100];   % ms
bin_ms      = 1;           % PSTH bin (ms)
psth_thresh = 200;         % threshold for red line

classes_here = unique(combClass_win(:)');   % stim sets present
nSets = numel(classes_here);

nCols = ceil(sqrt(nSets));
nRows = ceil(nSets / nCols);

% === Pre-pass: find global trial count & PSTH max ===
globalYmaxTrials = 0;
globalYmaxRate   = 0;

for ii = 1:nSets
    tr_idx = find(combClass_win == classes_here(ii));
    nThis  = numel(tr_idx);
    if nThis == 0, continue; end
    globalYmaxTrials = max(globalYmaxTrials, nThis);

    % Collect spike times
    spkTimesPerTrial = cell(nThis,1);
    for k = 1:nThis
        tr = tr_idx(k);
        if ~isempty(Spikes{ch_raster,tr}) && isfield(Spikes{ch_raster,tr},'t_ms')
            tt = Spikes{ch_raster,tr}.t_ms;
            spkTimesPerTrial{k} = tt(tt >= ras_win(1) & tt <= ras_win(2));
        else
            spkTimesPerTrial{k} = [];
        end
    end

    % PSTH
    edges = ras_win(1):bin_ms:ras_win(2);
    counts = zeros(1, numel(edges)-1);
    for k = 1:nThis
        counts = counts + histcounts(spkTimesPerTrial{k}, edges);
    end
    bin_s = bin_ms/1000;
    rate  = counts / (nThis * bin_s);

    % smoothing
    smooth_ms = 3; smooth_bins = max(1, round(smooth_ms/bin_ms));
    sigma_bins = max(1, round(smooth_bins/2));
    kk = 0:(smooth_bins-1);
    g = exp(-0.5*(kk./sigma_bins).^2); g = g/sum(g);
    rate_s = filter(g,1,rate);

    globalYmaxRate = max(globalYmaxRate, max(rate_s));
end

%% === Plot ===
fig = figure('Color','w','Name',sprintf('Ch %d | Raster+PSTH by stim set', ch_raster), ...
             'NumberTitle','off');

for ii = 1:nSets
    cc     = classes_here(ii);
    tr_idx = find(combClass_win == cc);
    nThis  = numel(tr_idx);
    if nThis == 0, continue; end

    % Stim-set label
    stimIdx = uniqueComb(cc,:); stimIdx = stimIdx(stimIdx>0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');

    % Subplot
    ax = subplot(nRows,nCols,ii); hold on;

    % Background if stim channel
    if ismember(ch_raster, stimIdx)
        set(ax,'Color',[1 0.9 0.9]); % light pink
    end

    % Collect spike times
    spkTimesPerTrial = cell(nThis,1);
    for k = 1:nThis
        tr = tr_idx(k);
        if ~isempty(Spikes{ch_raster,tr}) && isfield(Spikes{ch_raster,tr},'t_ms')
            tt = Spikes{ch_raster,tr}.t_ms;
            spkTimesPerTrial{k} = tt(tt >= ras_win(1) & tt <= ras_win(2));
        else
            spkTimesPerTrial{k} = [];
        end
    end

    % Raster
    for k = 1:nThis
        tt = spkTimesPerTrial{k};
        if isempty(tt), continue; end
        y0 = nThis - k + 1;
        for s = 1:numel(tt)
            plot([tt(s) tt(s)], [y0-0.4 y0+0.4], 'k', 'LineWidth', 0.8);
        end
    end

    % PSTH (no scaling, true spikes/s)
    edges = ras_win(1):bin_ms:ras_win(2);
    ctrs  = edges(1:end-1) + diff(edges)/2;
    counts = zeros(1, numel(edges)-1);
    for k = 1:nThis
        counts = counts + histcounts(spkTimesPerTrial{k}, edges);
    end
    bin_s = bin_ms/1000;
    rate  = counts / (nThis * bin_s);

    % smoothing
    smooth_ms = 3; smooth_bins = max(1, round(smooth_ms/bin_ms));
    sigma_bins = max(1, round(smooth_bins/2));
    kk = 0:(smooth_bins-1);
    g = exp(-0.5*(kk./sigma_bins).^2); g = g/sum(g);
    rate_s = filter(g,1,rate);

    if max(rate_s) > psth_thresh
        plot(ctrs, rate_s, 'r-', 'LineWidth', 1.5);
    else
        plot(ctrs, rate_s, 'b-', 'LineWidth', 1.5);
    end

    % Formatting
    xline(0,'r','LineWidth',1.2);
    xlim(ras_win);
    ylim([0 globalYmaxRate*1.1]);   % unified PSTH scale
    title(sprintf('Set %d: %s (Trials: %d)', ii, setLabel, nThis), ...
          'Interpreter','none','FontSize',9);
    xlabel('Time (ms)'); ylabel('Firing rate (sp/s)');
end

sgtitle(sprintf('Channel %d — Raster + PSTH by stimulation set', ch_raster), ...
        'FontWeight','bold');



% ------ (4) All 32 channels, Different Stim Set ------
for ss = 1: length(classes_here)
    stim_set_to_plot = classes_here(ss);   % pick stim set
    ras_win     = [-20 50];              % ms
    bin_ms      = 1;                      % PSTH bin (ms)
    psth_thresh = 200;                    % sp/s threshold
    
    tr_idx = find(combClass_win == stim_set_to_plot);
    nThis  = numel(tr_idx);
    stimIdx = uniqueComb(stim_set_to_plot,:); stimIdx = stimIdx(stimIdx>0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');
    
    
    edges = ras_win(1):bin_ms:ras_win(2);
    ctrs  = edges(1:end-1) + diff(edges)/2;
    bin_s = bin_ms/1000;
    
    figure('Color','w','Name',sprintf('Raster Plot %s', setLabel),'NumberTitle','off');
    for ch = 1:nChn
        ax = subplot(4,8,ch); hold(ax,'on');
    
        % pink background if stim channel
        if ismember(ch, stimIdx)
            set(ax,'Color',[1 0.92 0.92]);
        end
    
        % ---- Collect spike times ----
        spkTimesPerTrial = cell(nThis,1);
        totalSpikesThisChannel = 0;
        for k = 1:nThis
            tr = tr_idx(k);
            if ~isempty(Spikes{ch,tr}) && isfield(Spikes{ch,tr},'t_ms')
                tt = Spikes{ch,tr}.t_ms;
                spkTimesPerTrial{k} = tt(tt >= ras_win(1) & tt <= ras_win(2));
                totalSpikesThisChannel = totalSpikesThisChannel + numel(spkTimesPerTrial{k}); 
            else
                spkTimesPerTrial{k} = [];
            end
        end
    
        % ---- PSTH ----
        counts = zeros(1, numel(edges)-1);
        for k = 1:nThis
            counts = counts + histcounts(spkTimesPerTrial{k}, edges);
        end
        rate  = counts / (nThis * bin_s);
    
        % smoothing
        smooth_ms   = 3;
        smooth_bins = max(1, round(smooth_ms/bin_ms));
        sigma_bins  = max(1, round(smooth_bins/2));
        kk = 0:(smooth_bins-1);
        g  = exp(-0.5*(kk./sigma_bins).^2); g = g/sum(g);
        rate_s = filter(g,1,rate);
    
        % y-limits + PSTH color
        ymax_rate = max(rate_s);
        if ymax_rate > psth_thresh
            psth_color = [0.85 0 0];    % red
            yLimUse = [0, ceil(ymax_rate/10)*10];
        else
            psth_color = [0 0 0];       % black
            yLimUse = [0 psth_thresh];
        end
        plot(ax, ctrs, rate_s, '-', 'Color', psth_color, 'LineWidth', 2.5);
    
        % ---- Raster overlay ----
        yTop = max(yLimUse(2), 1);
        yRasterLow  = 5;
        yRasterHigh = 0.90 * yTop;
        for k = 1:nThis
            tt = spkTimesPerTrial{k};
            if isempty(tt), continue; end
            y0 = yRasterLow + (k-0.5)/nThis * (yRasterHigh - yRasterLow);
            for s = 1:numel(tt)
                plot(ax, [tt(s) tt(s)], [y0-0.005*yTop y0+0.005*yTop], 'k', 'Color', [0.2 0.2 0.2],'LineWidth', 0.5);
            end
        end
    
        % formatting
        xline(ax, 0, 'r', 'LineWidth', 1.2);
        xlim(ax, ras_win);
        ylim(ax, yLimUse);
    
        % mark stim channels with *
        if ismember(ch, stimIdx)
            title(ax, sprintf('Ch %d* | n=%d', ch, totalSpikesThisChannel), 'FontSize', 9);
        else
            title(ax, sprintf('Ch %d | n=%d', ch, totalSpikesThisChannel), 'FontSize', 9);
        end
    
        if ch > 24
            xlabel(ax, 'Time (ms)');
        end
        if mod(ch-1,8)==0
            ylabel(ax, 'Rate (sp/s)');
        end
        box(ax,'off');
    end
    sgtitle(sprintf('Stim set: %s ', setLabel), ...
            'FontWeight','bold');
end

% ------ (5) One channel, different amplitude ------ %
ch_raster   = 11;         % channel to visualize
ras_win     = [-20 60];  % ms window
bin_ms      = 1;          % PSTH bin size (ms)
psth_thresh = 200;        % highlight PSTH red if max > threshold
only_nonempty = true;     % if true, skip amplitudes with 0 trials
limitToStimSet = [];      % [] = all stim sets, or set to a class id

amps_to_plot = 1:n_AMP;
if only_nonempty
    has_trials = arrayfun(@(k) any(ampIdx == k), amps_to_plot);
    amps_to_plot = amps_to_plot(has_trials);
end
nA = numel(amps_to_plot);

% y-limit
globalMaxRate = 0;
for a = 1:nA
    kAmp = amps_to_plot(a);
    tr_idx = find(ampIdx == kAmp);
    if ~isempty(limitToStimSet)
        tr_idx = tr_idx(combClass_win(tr_idx) == limitToStimSet);
    end
    nThis = numel(tr_idx);
    if nThis == 0, continue; end

    % collect PSTH counts
    edges = ras_win(1):bin_ms:ras_win(2);
    counts = zeros(1, numel(edges)-1);
    for t = 1:nThis
        if isempty(Spikes{ch_raster,tr_idx(t)}), continue; end
        tt = Spikes{ch_raster,tr_idx(t)}.t_ms;
        counts = counts + histcounts(tt, edges);
    end
    bin_s = bin_ms/1000;
    rate = counts / (nThis * bin_s);

    % smooth
    smooth_ms = 3; smooth_bins = max(1, round(smooth_ms/bin_ms));
    sigma_bins = max(1, round(smooth_bins/2));
    kk = 0:(smooth_bins-1);
    g  = exp(-0.5*(kk./sigma_bins).^2); g = g/sum(g);
    rate_s = filter(g,1,rate);

    globalMaxRate = max(globalMaxRate, max(rate_s));
end
yLimUse = [0, ceil(globalMaxRate*1.1/10)*10]; % round up

% Plot all amplitudes
nCols = ceil(sqrt(nA));
nRows = ceil(nA / nCols);

fig = figure('Color','w','Name',sprintf('Ch %d | Raster+PSTH by amplitude', ch_raster), ...
             'NumberTitle','off');

for a = 1:nA
    kAmp = amps_to_plot(a);
    amp_uA = Amps(kAmp);
    tr_idx = find(ampIdx == kAmp);
    if ~isempty(limitToStimSet)
        tr_idx = tr_idx(combClass_win(tr_idx) == limitToStimSet);
    end
    nThis = numel(tr_idx);
    if nThis == 0, continue; end

    % Collect spike times
    spkTimesPerTrial = cell(nThis,1);
    for t = 1:nThis
        tr = tr_idx(t);
        if ~isempty(Spikes{ch_raster,tr}) && isfield(Spikes{ch_raster,tr},'t_ms')
            tt = Spikes{ch_raster,tr}.t_ms;
            spkTimesPerTrial{t} = tt(tt >= ras_win(1) & tt <= ras_win(2));
        else
            spkTimesPerTrial{t} = [];
        end
    end

    % Subplot
    ax = subplot(nRows, nCols, a); hold(ax,'on');

    % PSTH
    edges = ras_win(1):bin_ms:ras_win(2);
    ctrs  = edges(1:end-1) + diff(edges)/2;
    counts = zeros(1, numel(edges)-1);
    for t = 1:nThis
        counts = counts + histcounts(spkTimesPerTrial{t}, edges);
    end
    bin_s = bin_ms/1000;
    rate  = counts / (nThis * bin_s);

    % smooth
    smooth_ms = 3; smooth_bins = max(1, round(smooth_ms/bin_ms));
    sigma_bins = max(1, round(smooth_bins/2));
    kk = 0:(smooth_bins-1);
    g  = exp(-0.5*(kk./sigma_bins).^2); g = g/sum(g);
    rate_s = filter(g,1,rate);

    % PSTH line
    if max(rate_s) > psth_thresh
        plot(ax, ctrs, rate_s, 'r-', 'LineWidth', 2.5);
    else
        plot(ax, ctrs, rate_s, 'k-', 'LineWidth', 2.5);
    end

    % Raster overlay (lighter)
    yTop        = yLimUse(2);
    yRasterLow  = 5;
    yRasterHigh = 0.90 * yTop;
    tickHalf    = 0.005 * yTop;
    for t = 1:nThis
        tt = spkTimesPerTrial{t};
        if isempty(tt), continue; end
        y0 = yRasterLow + (t-0.5)/nThis * (yRasterHigh - yRasterLow);
        for s = 1:numel(tt)
            plot(ax, [tt(s) tt(s)], [y0 - tickHalf, y0 + tickHalf], ...
                 'Color', [0.2 0.2 0.2], 'LineWidth', 0.5);
        end
    end

    % Stim line
    xline(ax, 0, 'r-', 'LineWidth', 1.0);

    % Cosmetics
    xlim(ax, ras_win);
    ylim(ax, yLimUse);
    title(ax, sprintf('%g \\muA  (Trials: %d)', amp_uA, nThis), 'FontSize', 9);
    xlabel(ax, 'Time (ms)'); ylabel(ax, 'Rate (sp/s)');
    box(ax,'off');
end

sgtitle(sprintf('Channel %d — Raster + PSTH by amplitude', ch_raster), 'FontWeight','bold');


% ------ (6) Different Stim Set & Different AMP
ch_raster    = 9;         % channel to visualize
ras_win      = [-20 80];  % ms window
bin_ms       = 1;          % PSTH bin size (ms)
psth_thresh  = 300;        % highlight PSTH red if max > threshold

% Optional filters (leave [] to include all)
sets_to_plot = [];         % e.g., [1 3 5] for specific stim-set classes
amps_to_plot = 1:n_AMP;    % indices into Amps
only_nonempty = true;      % skip empty set×amp cells

% ---- pick which stim sets to show ----
all_sets = unique(combClass_win(:))';
if ~isempty(sets_to_plot)
    sets_here = intersect(all_sets, sets_to_plot);
else
    sets_here = all_sets;
end

% If requested, drop amplitudes that have no trials anywhere
if only_nonempty
    has_trials_any = false(size(amps_to_plot));
    for jj = 1:numel(amps_to_plot)
        kAmp = amps_to_plot(jj);
        has_trials_any(jj) = any(ampIdx == kAmp);
    end
    amps_here = amps_to_plot(has_trials_any);
else
    amps_here = amps_to_plot;
end

nSets = numel(sets_here);
nAmps = numel(amps_here);
if nSets == 0 || nAmps == 0
    warning('Nothing to plot for the chosen filters.'); 
    return
end

%  global y-limit (max firing rate) 
edges = ras_win(1):bin_ms:ras_win(2);
bin_s = bin_ms/1000;
globalMaxRate = 0;

for ii = 1:nSets
    set_id = sets_here(ii);
    for jj = 1:nAmps
        kAmp = amps_here(jj);
        tr_idx = find(ampIdx == kAmp & combClass_win == set_id);
        nThis  = numel(tr_idx);
        if nThis == 0, continue; end

        counts = zeros(1, numel(edges)-1);
        for t = 1:nThis
            S = Spikes{ch_raster, tr_idx(t)};
            if isempty(S), continue; end
            counts = counts + histcounts(S.t_ms, edges);
        end
        rate = counts / (nThis * bin_s);

        % causal Gaussian-like smoothing (same as before)
        smooth_ms  = 3;
        smooth_bins = max(1, round(smooth_ms/bin_ms));
        sigma_bins  = max(1, round(smooth_bins/2));
        k = 0:(smooth_bins-1);
        g = exp(-0.5*(k./sigma_bins).^2); g = g/sum(g);
        rate_s = filter(g,1,rate);

        globalMaxRate = max(globalMaxRate, max(rate_s));
    end
end

% ---------- Figure layout ----------
figName = sprintf('Ch %d | Raster+PSTH by stim set × amplitude', ch_raster);
fig = figure('Color','w','Name',figName,'NumberTitle','off');
tl  = tiledlayout(nSets, nAmps, 'TileSpacing','compact','Padding','compact');
sgtitle(tl, figName, 'FontWeight','bold');

ctrs = edges(1:end-1) + diff(edges)/2;

for ii = 1:nSets
    set_id  = sets_here(ii);
    stimIdx = uniqueComb(set_id,:); stimIdx = stimIdx(stimIdx>0);

    % Label for stim set
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');
    

    for jj = 1:nAmps
        kAmp = amps_here(jj);
        amp_uA = Amps(kAmp);

        ax = nexttile(tl); hold(ax,'on');

        % background if ch_raster is stimulated in this set
        if ismember(ch_raster, stimIdx)
            set(ax, 'Color', [1 0.92 0.92]);  % light red
        end

        % Trials for this (set, amp)
        tr_idx = find(ampIdx == kAmp & combClass_win == set_id);
        nThis  = numel(tr_idx);
        if nThis == 0
            % empty cell -> annotate and move on
            title(ax, sprintf('%s | %g \\muA (n=0)', setLabel, amp_uA), ...
                  'Interpreter','none','FontSize',8);
            xlim(ax, ras_win); ylim(ax, yLimUse);
            box(ax,'off'); continue;
        end

        % Collect spike times restricted to window
        spkPerTrial = cell(nThis,1);
        for t = 1:nThis
            S = Spikes{ch_raster, tr_idx(t)};
            if isempty(S) || ~isfield(S,'t_ms'), spkPerTrial{t} = []; continue; end
            tt = S.t_ms;
            spkPerTrial{t} = tt(tt >= ras_win(1) & tt <= ras_win(2));
        end

        % ---- PSTH ----
        counts = zeros(1, numel(edges)-1);
        for t = 1:nThis
            counts = counts + histcounts(spkPerTrial{t}, edges);
        end
        rate  = counts / (nThis * bin_s);

        smooth_ms  = 3;
        smooth_bins = max(1, round(smooth_ms/bin_ms));
        sigma_bins  = max(1, round(smooth_bins/2));
        k = 0:(smooth_bins-1);
        g = exp(-0.5*(k./sigma_bins).^2); g = g/sum(g);
        rate_s = filter(g,1,rate);
        yLimUse = [0, max(10, 10*ceil(1.1*globalMaxRate/10))];
        if max(rate_s) > psth_thresh
            % yLimUse = [0, max(10, 10*ceil(1.1*globalMaxRate/10))];
            plot(ax, ctrs, rate_s, '-', 'Color', [0.6 0 0], 'LineWidth', 2.5);
        else
            % yLimUse = [0 psth_thresh];
            plot(ax, ctrs, rate_s, 'k-', 'LineWidth', 2.5);
        end

        % ---- Raster overlay ----
        yTop        = yLimUse(2);
        yRasterLow  = 5;
        yRasterHigh = 0.90 * yTop;
        tickHalf    = 0.004 * yTop;       % thin ticks
        for t = 1:nThis
            tt = spkPerTrial{t};
            if isempty(tt), continue; end
            y0 = yRasterLow + (t-0.5)/nThis * (yRasterHigh - yRasterLow);
            for s = 1:numel(tt)
                plot(ax, [tt(s) tt(s)], [y0 - tickHalf, y0 + tickHalf], ...
                     'Color', [0 0 0 0.35], 'LineWidth', 0.6);
            end
        end
        xline(ax, 0, 'r-', 'LineWidth', 1);
        xlim(ax, ras_win); ylim(ax, yLimUse);
        if ii == nSets, xlabel(ax, 'Time (ms)'); end
        if jj == 1,     ylabel(ax, 'Rate (sp/s)'); end
        title(ax, sprintf('%s | %g µA (n=%d)', setLabel, amp_uA, nThis), ...
              'Interpreter','none','FontSize',11);
        box(ax,'off');
    end
end

% ------- (7) One channel, PSTH of all AMP -------- %
ch_raster   = 9;         % channel to visualize
ras_win     = [-20 60];   % ms window
bin_ms      = 1;          % PSTH bin size (ms)
psth_thresh = 200;        % threshold for red highlight
only_nonempty = true;
limitToStimSet = [];      % [] = all stim sets, or set to a class id

amps_to_plot = 1:n_AMP;
if only_nonempty
    has_trials = arrayfun(@(k) any(ampIdx == k), amps_to_plot);
    amps_to_plot = amps_to_plot(has_trials);
end
nA = numel(amps_to_plot);

% Prepare figure
fig = figure('Color','w','Name',sprintf('Ch %d | PSTH overlay by amplitude', ch_raster), ...
             'NumberTitle','off');
hold on;

cmap = lines(nA);  % distinct colors per amplitude
globalMaxRate = 0;

for a = 1:nA
    kAmp = amps_to_plot(a);
    amp_uA = Amps(kAmp);

    % Trials for this amplitude (and stim set filter if requested)
    tr_idx = find(ampIdx == kAmp);
    if ~isempty(limitToStimSet)
        tr_idx = tr_idx(combClass_win(tr_idx) == limitToStimSet);
    end
    nThis = numel(tr_idx);
    if nThis == 0, continue; end

    % Collect spike times
    edges = ras_win(1):bin_ms:ras_win(2);
    counts = zeros(1, numel(edges)-1);
    for t = 1:nThis
        tr = tr_idx(t);
        if isempty(Spikes{ch_raster,tr}), continue; end
        tt = Spikes{ch_raster,tr}.t_ms;
        counts = counts + histcounts(tt, edges);
    end
    bin_s = bin_ms/1000;
    rate  = counts / (nThis * bin_s);

    % Smooth
    smooth_ms = 3; smooth_bins = max(1, round(smooth_ms/bin_ms));
    sigma_bins = max(1, round(smooth_bins/2));
    kk = 0:(smooth_bins-1);
    g  = exp(-0.5*(kk./sigma_bins).^2); g = g/sum(g);
    rate_s = filter(g,1,rate);

    % Update max
    globalMaxRate = max(globalMaxRate, max(rate_s));

    % Plot PSTH
    plot(edges(1:end-1) + diff(edges)/2, rate_s, ...
         'LineWidth', 2, 'Color', cmap(a,:), ...
         'DisplayName', sprintf('%g µA (n=%d)', amp_uA, nThis));
end

xline(0,'r--','LineWidth',2);
xlim(ras_win);
ylim([0, ceil(globalMaxRate*1.1/10)*10]);
xlabel('Time (ms)');
ylabel('Firing rate (sp/s)');
title(sprintf('Channel %d — PSTH overlay by amplitude', ch_raster));
legend('show','Location','northeast');
box off;

% -------- (8) One CHN, PSTH of different AMP, different Stim Set ------- %
ch_raster   = 11;           % channel to visualize
ras_win     = [-20 60];    % ms window
bin_ms      = 1;           % PSTH bin size (ms)
smooth_ms   = 3;           % Gaussian FIR smoothing (causal)
only_nonempty = true;      % keep only amps that appear anywhere
classes_here = unique(combClass_win(:)');   % stim sets present
nSets        = numel(classes_here);

% Choose amplitudes that appear at least once (globally)
amps_to_plot = 1:n_AMP;
if only_nonempty
    has_trials_any = arrayfun(@(k) any(ampIdx == k), amps_to_plot);
    amps_to_plot = amps_to_plot(has_trials_any);
end
nA = numel(amps_to_plot);

% Time binning
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;

% Smoothing kernel (causal Gaussian FIR)
smooth_bins = max(1, round(smooth_ms/bin_ms));
sigma_bins  = max(1, round(smooth_bins/2));
kk = 0:(smooth_bins-1);
g  = exp(-0.5*(kk./sigma_bins).^2); g = g / sum(g);

% ---------- Pass 1: compute global Y limit across ALL sets & amps ----------
globalMaxRate = 0;
hasData_set_amp = false(nSets, nA);   % track where data exists

for ii = 1:nSets
    set_id = classes_here(ii);
    tr_set = find(combClass_win == set_id);

    for ai = 1:nA
        kAmp = amps_to_plot(ai);
        tr_idx = intersect(tr_set, find(ampIdx == kAmp));
        nThis  = numel(tr_idx);
        if nThis == 0, continue; end

        counts = zeros(1, numel(edges)-1);
        for t = 1:nThis
            tr = tr_idx(t);
            S  = Spikes{ch_raster, tr};
            if isempty(S) || ~isfield(S,'t_ms') || isempty(S.t_ms), continue; end
            counts = counts + histcounts(S.t_ms, edges);
        end
        rate   = counts / (nThis * bin_s);
        rate_s = filter(g, 1, rate);

        hasData_set_amp(ii, ai) = true;
        globalMaxRate = max(globalMaxRate, max(rate_s));
    end
end
yLimUse = [0, max(1, ceil(globalMaxRate*1.1/10)*10)];  % shared y-limit

% ---------- Pass 2: plot one panel per stim set with PSTH overlays ----------
nCols = ceil(sqrt(nSets));
nRows = ceil(nSets / nCols);

fig = figure('Color','w','Name',sprintf('Ch %d | PSTH by amp × stim set', ch_raster), ...
             'NumberTitle','off');
tl  = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

cmap      = lines(nA);   % color per amplitude
ampLabels = arrayfun(@(k) sprintf('%.0f \\muA', Amps(k)), amps_to_plot, 'UniformOutput', false);
legendHandles = gobjects(1, nA);  % collect one handle per amp for global legend
haveLegendHandle = false(1, nA);

for ii = 1:nSets
    set_id = classes_here(ii);
    % Stim electrodes for this set -> label and panel highlight
    stimIdx = uniqueComb(set_id,:); stimIdx = stimIdx(stimIdx > 0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');

    ax = nexttile(tl); hold(ax,'on');

    % light red background if this plotted channel is a stim electrode
    if ismember(ch_raster, stimIdx)
        set(ax,'Color',[1 0.92 0.92]);
    end

    tr_set = find(combClass_win == set_id);

    % overlay one PSTH per amplitude (only those with data in this set)
    for ai = 1:nA
        if ~hasData_set_amp(ii, ai), continue; end
        kAmp  = amps_to_plot(ai);
        tr_idx = intersect(tr_set, find(ampIdx == kAmp));
        nThis  = numel(tr_idx);
        if nThis == 0, continue; end

        counts = zeros(1, numel(edges)-1);
        for t = 1:nThis
            tr = tr_idx(t);
            S  = Spikes{ch_raster, tr};
            if isempty(S) || ~isfield(S,'t_ms') || isempty(S.t_ms), continue; end
            counts = counts + histcounts(S.t_ms, edges);
        end
        rate   = counts / (nThis * bin_s);
        rate_s = filter(g, 1, rate);

        % plot line, and capture first handle per amplitude for the global legend
        if ~haveLegendHandle(ai)
            legendHandles(ai)   = plot(ax, ctrs, rate_s, 'LineWidth', 2, ...
                                       'Color', cmap(ai,:), 'DisplayName', ampLabels{ai});
            haveLegendHandle(ai) = true;
        else
            plot(ax, ctrs, rate_s, 'LineWidth', 2, 'Color', cmap(ai,:));
        end
    end

    xline(ax, 0, 'r-', 'LineWidth', 1.0);
    xlim(ax, ras_win);
    ylim(ax, yLimUse);
    title(ax, sprintf('Set %d: %s', ii, setLabel), 'FontSize', 9, 'Interpreter','none');
    if ii > (nRows-1)*nCols, xlabel(ax,'Time (ms)'); end
    if mod(ii-1, nCols) == 0, ylabel(ax,'Rate (sp/s)'); end
    box(ax,'off');
end

% Global legend (centered bottom) using only amplitudes that actually appeared anywhere
valid = isgraphics(legendHandles);
lgd = legend(legendHandles(valid), ampLabels(valid), ...
             'Orientation','horizontal', 'Location','southoutside');
lgd.Box = 'off';

% Global title
sgtitle(tl, sprintf('Channel %d — PSTH overlays by amplitude within each stimulation set', ch_raster), ...
        'FontWeight','bold');

% -------- (9) All CHNs, PSTH of different AMP ------- %
ras_win        = [-20 60];   % ms
bin_ms         = 1;          % PSTH bin size (ms)
smooth_ms      = 3;          % Gaussian FIR smoothing (causal)
limitToStimSet = [];         % [] -> all stim sets; or set to a class id

% Pick amplitudes that actually exist in the dataset (global, not per ch)
amps_to_plot = 1:n_AMP;
has_trials_global = arrayfun(@(k) any(ampIdx == k), amps_to_plot);
amps_to_plot = amps_to_plot(has_trials_global);
nA = numel(amps_to_plot);

% Time binning
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;

% Smoothing kernel (causal Gaussian FIR)
smooth_bins = max(1, round(smooth_ms/bin_ms));
sigma_bins  = max(1, round(smooth_bins/2));
kk = 0:(smooth_bins-1);
g  = exp(-0.5*(kk./sigma_bins).^2); g = g / sum(g);

% Precompute PSTHs and global max for y-limits
PSTH = cell(nChn, nA);        % each cell: 1 x (num bins) smoothed rate
hasData = false(nChn, nA);
globalMaxRate = 0;

for ch = 1:nChn
    for ai = 1:nA
        kAmp   = amps_to_plot(ai);
        tr_idx = find(ampIdx == kAmp);
        if ~isempty(limitToStimSet)
            tr_idx = tr_idx(combClass_win(tr_idx) == limitToStimSet);
        end
        nThis = numel(tr_idx);
        if nThis == 0, continue; end

        counts = zeros(1, numel(edges)-1);
        for t = 1:nThis
            tr = tr_idx(t);
            S  = Spikes{ch, tr};
            if isempty(S) || ~isfield(S,'t_ms') || isempty(S.t_ms), continue; end
            counts = counts + histcounts(S.t_ms, edges);
        end

        rate = counts / (nThis * bin_s);
        rate_s = filter(g, 1, rate);            % causal smoothing
        PSTH{ch, ai} = rate_s;
        hasData(ch, ai) = true;
        globalMaxRate = max(globalMaxRate, max(rate_s));
    end
end

yLimUse = [0, max(1, ceil(globalMaxRate*1.1/10)*10)];  % shared y-limit

% Figure & layout
nCols = 4;                        % 8x4 grid typical for 32 channels
nRows = ceil(nChn / nCols);
fig = figure('Color','w', 'Name','PSTH overlay by amplitude — all channels', ...
             'NumberTitle','off');
cmap = lines(nA);
ampLabels = arrayfun(@(k) sprintf('%.0f \\muA', Amps(k)), amps_to_plot, 'UniformOutput', false);
lineHandles = gobjects(1, nA);    % for a clean global legend, capture first use

for ch = 1:nChn
    ax = nexttile; hold(ax,'on');
    firstPlotted = false(1, nA);

    for ai = 1:nA
        if ~hasData(ch, ai), continue; end
        y = PSTH{ch, ai};
        if ~firstPlotted(ai)
            % keep a handle for legend only once (first channel where it appears)
            lineHandles(ai) = plot(ax, ctrs, y, 'LineWidth', 2, 'Color', cmap(ai,:), ...
                                   'DisplayName', ampLabels{ai});
            firstPlotted(ai) = true;
        else
            plot(ax, ctrs, y, 'LineWidth', 2, 'Color', cmap(ai,:));
        end
    end

    xline(ax, 0, 'r-', 'LineWidth', 1.0);
    xlim(ax, ras_win);
    ylim(ax, yLimUse);
    title(ax, sprintf('Ch %d', ch), 'FontSize', 9);
    if ch > (nRows-1)*nCols, xlabel(ax, 'Time (ms)'); end
    if mod(ch-1, nCols) == 0, ylabel(ax, 'Rate (sp/s)'); end
    box(ax,'off');
    % axis(ax,'square');
end

% Single overall legend (bottom, outside)
validHandles = isgraphics(lineHandles);
lgd = legend(lineHandles(validHandles), ampLabels(validHandles), ...
             'Orientation','horizontal', 'Location','southoutside');
lgd.Box = 'off';
lgd.Layout.Tile = 'south';   % attach legend to bottom row of layout


sgtitle('PSTH by amplitude, all channels', 'FontWeight','bold');

% -------- (10) All CHNs, PSTH of different AMP, Stim Set ------- %
classes_here = unique(combClass_win(:)');   % stimulation sets present

ras_win        = [-20 60];   % ms
bin_ms         = 1;          % PSTH bin size (ms)
smooth_ms      = 3;          % Gaussian FIR smoothing (causal)

% Pick amplitudes that actually exist (global)
amps_to_plot = 1:n_AMP;
has_trials_global = arrayfun(@(k) any(ampIdx == k), amps_to_plot);
amps_to_plot = amps_to_plot(has_trials_global);
nA = numel(amps_to_plot);

% Time binning
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;

% Smoothing kernel (causal Gaussian FIR)
smooth_bins = max(1, round(smooth_ms/bin_ms));
sigma_bins  = max(1, round(smooth_bins/2));
kk = 0:(smooth_bins-1);
g  = exp(-0.5*(kk./sigma_bins).^2); g = g / sum(g);

cmap = lines(nA);
ampLabels = arrayfun(@(k) sprintf('%.0f \\muA', Amps(k)), amps_to_plot, 'UniformOutput', false);

for ss = 1:numel(classes_here)
    set_id = classes_here(ss);

    % Stim set label
    stimIdx = uniqueComb(set_id,:); stimIdx = stimIdx(stimIdx>0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');

    % ---- Compute PSTHs (same as before) ----
    PSTH = cell(nChn, nA);
    hasData = false(nChn, nA);
    globalMaxRate = 0;

    for ch = 1:nChn
        for ai = 1:nA
            kAmp = amps_to_plot(ai);
            tr_idx = find(ampIdx == kAmp & combClass_win == set_id);
            nThis = numel(tr_idx);
            if nThis == 0, continue; end

            counts = zeros(1, numel(edges)-1);
            for t = 1:nThis
                tr = tr_idx(t);
                S  = Spikes{ch, tr};
                if isempty(S) || ~isfield(S,'t_ms') || isempty(S.t_ms), continue; end
                counts = counts + histcounts(S.t_ms, edges);
            end

            rate   = counts / (nThis * bin_s);
            rate_s = filter(g, 1, rate);
            PSTH{ch, ai} = rate_s;
            hasData(ch, ai) = true;
            globalMaxRate   = max(globalMaxRate, max(rate_s));
        end
    end

    yLimUse = [0, max(1, ceil(globalMaxRate*1.1/10)*10)];

    % ---- Figure ----
    nCols = 4; nRows = ceil(nChn/nCols);
    fig = figure('Color','w','Name',sprintf('PSTH overlay | Set %s', setLabel),...
                 'NumberTitle','off');
    cmap = lines(nA);
    ampLabels = arrayfun(@(k) sprintf('%.0f \\muA', Amps(k)), amps_to_plot,'UniformOutput',false);
    lineHandles = gobjects(1,nA);
    firstSeen   = false(1,nA);

    % ---- Subplots ----
    for ch = 1:nChn
        ax = subplot(nRows,nCols,ch); hold(ax,'on');

         if ismember(ch, stimIdx)
            set(ax, 'Color', [1 0.92 0.92]);  % light red
         end

        for ai = 1:nA
            if ~hasData(ch,ai), continue; end
            y = PSTH{ch,ai};

            if ~firstSeen(ai)
                lineHandles(ai) = plot(ax, ctrs, y, 'LineWidth', 2, 'Color', cmap(ai,:), ...
                                       'DisplayName', ampLabels{ai});
                firstSeen(ai) = true;
            else
                plot(ax, ctrs, y, 'LineWidth', 2, 'Color', cmap(ai,:));
            end
        end

        xline(ax,0,'r-','LineWidth',1.0);
        xlim(ax, ras_win); ylim(ax, yLimUse);
        title(ax, sprintf('Ch %d', ch), 'FontSize', 9);
        box(ax,'off');

    end
    valid = isgraphics(lineHandles);
    lgd = legend(lineHandles(valid), ampLabels(valid), ...
        'Orientation','horizontal','NumColumns',min(nA,6),'Box','off');
    lgd.Position = [0.3 0.01 0.4 0.05];   % [x y w h] normalized units
    sgtitle(sprintf('PSTH by amplitude — all channels | Stim set: %s', setLabel), ...
            'FontWeight','bold');
end
%% ------ Tuning Curve ------ %%

% -------- (1) Firing rate vs. AMP (All sets) -------- %
ch_tune  = 22;                 % channel 
resp_ms  = [0 20];            % response window (ms)
resp_s   = diff(resp_ms)/1000;

classes_here = unique(combClass_win(:)'); 
n_AMP   = numel(Amps);

meanFR_amp = nan(1, n_AMP);
semFR_amp  = nan(1, n_AMP);

for k = 1:n_AMP
    trials_k = find(ampIdx == k);
    % count spikes in window for each trial
    counts   = zeros(numel(trials_k),1);
    for ii = 1:numel(trials_k)
        tr = trials_k(ii);
        S  = Spikes{ch_tune, tr};
        if ~isempty(S) && isfield(S,'t_ms') && ~isempty(S.t_ms)
            counts(ii) = sum(S.t_ms >= resp_ms(1) & S.t_ms <= resp_ms(2));
        end
    end
    fr = counts / resp_s;             % spikes/s per trial
    meanFR_amp(k) = mean(fr,'omitnan');
    semFR_amp(k)  = std(fr,'omitnan') / sqrt(max(1,numel(fr)));
end

figure('Color','w','Name',sprintf('Ch %d | FR vs Amplitude', ch_tune),'NumberTitle','off');
errorbar(Amps, meanFR_amp, semFR_amp, '-o','LineWidth',1.8);
xlabel('Amplitude (\muA)'); ylabel('Firing rate (sp/s)');
title(sprintf('Channel %d — %d–%d ms window', ch_tune, resp_ms(1), resp_ms(2)));
box off;

% -------- (2) Firing rate vs. AMP (DIfferent sets) -------- %
ch_tune = 9;                 % channel
resp_ms = [0 20];            % response window (ms, relative to stim)
resp_s  = diff(resp_ms)/1000;

classes_here = unique(combClass_win(:)');   % stim sets present
nSets        = numel(classes_here);
n_AMP        = numel(Amps);

% Prepare legend labels for each stim set
setLabels = strings(1, nSets);
for ii = 1:nSets
    set_id  = classes_here(ii);
    stimIdx = uniqueComb(set_id, :); 
    stimIdx = stimIdx(stimIdx > 0);
    setLabels(ii) = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx, 'UniformOutput', false), ' + ');
end

% Compute mean FR (and SEM) per amplitude within each stim set
meanFR_set_amp = nan(nSets, n_AMP);
semFR_set_amp  = nan(nSets, n_AMP);

for ii = 1:nSets
    set_id   = classes_here(ii);
    trials_S = find(combClass_win == set_id);  % trials belonging to this stim set

    for k = 1:n_AMP
        trials_kS = intersect(trials_S, find(ampIdx == k));   % trials with this amp AND this set
        if isempty(trials_kS), continue; end

        % spike counts per trial in window
        counts = zeros(numel(trials_kS),1);
        for jj = 1:numel(trials_kS)
            tr = trials_kS(jj);
            S  = Spikes{ch_tune, tr};
            if ~isempty(S) && isfield(S,'t_ms') && ~isempty(S.t_ms)
                counts(jj) = sum(S.t_ms >= resp_ms(1) & S.t_ms <= resp_ms(2));
            end
        end
        fr = counts / resp_s;  % sp/s for each trial
        meanFR_set_amp(ii,k) = mean(fr, 'omitnan');
        semFR_set_amp(ii,k)  = std(fr,  'omitnan') / sqrt(max(1, numel(fr)));
        stdFR_set_amp(ii,k) = std(fr, 'omitnan');   % standard deviation
    end
end

% Plot: overlay one line per stim set across amplitude
figure('Color','w','Name',sprintf('Ch %d | FR vs Amplitude by Stim Set', ch_tune),'NumberTitle','off');
hold on;
cmap = lines(nSets);   % distinct colors per set
pltHandles = []; 
for ii = 1:nSets
    y = meanFR_set_amp(ii,:);
    % e = semFR_set_amp(ii,:);
    e = stdFR_set_amp(ii,:);
    % Only plot where we actually have data for this set
    hasData = ~isnan(y);
    if any(hasData)
        % Keep Amps aligned with available data
       h = plot(Amps(hasData), y(hasData), '-o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', cmap(ii,:));
       hErr =  errorbar(Amps(hasData), y(hasData), e(hasData), 'LineStyle', 'none', 'Color', cmap(ii,:), 'LineWidth', 1.5);
        
       set(hErr, 'HandleVisibility', 'off'); 
       pltHandles(end+1) = h;
    end
end

box off;
xlabel('Amplitude (µA)');
ylabel('Firing rate (sp/s)');
title(sprintf('Channel %d — %d to %d ms', ch_tune, resp_ms(1), resp_ms(2)));
% leg = legend(setLabels,'Location', 'SouthEast');
% leg.Box = 'off';

if ~isempty(pltHandles)
    leg = legend(pltHandles,setLabels, 'Location', 'SouthEast', 'Interpreter','none');
    leg.Box = 'off';
end

% Optional: enforce same x-limits across all sets, and nice y-limits
xlim([min(Amps) max(Amps)]);
yl = ylim; ylim([0 max(yl(2), 1)]);


%% ------ Heat Map ------ %% 
 % -------- (1) FR per CHN per AMP per Stim Set -------- % 
resp_ms  = [0 20];                % response window (ms)
resp_s   = diff(resp_ms)/1000;
classes_here = unique(combClass_win(:)');   % stim sets present
nSets    = numel(classes_here);

% loop over stim sets
for si = 1:nSets
    set_id   = classes_here(si);
    stimIdx  = uniqueComb(set_id, :); stimIdx = stimIdx(stimIdx > 0);
    if exist('E_NAME','var') && ~isempty(E_NAME) && max(stimIdx) <= numel(E_NAME)
        setLabel = strjoin(E_NAME(stimIdx), ' + ');
    else
        setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimIdx,'UniformOutput',false),' + ');
    end
    
    % trials of this stim set
    trials_S = find(combClass_win == set_id);
    
    % Preallocate [channels × amplitudes]
    FR_mean = nan(nChn, n_AMP);
    maxFR   = 0;
    
    % compute FR for each channel × amplitude
    for ch = 1:nChn
        for k = 1:n_AMP
            trials_kS = intersect(trials_S, find(ampIdx == k));
            if isempty(trials_kS), continue; end
            
            counts = zeros(numel(trials_kS),1);
            for jj = 1:numel(trials_kS)
                tr = trials_kS(jj);
                S  = Spikes{ch, tr};
                if ~isempty(S) && isfield(S,'t_ms') && ~isempty(S.t_ms)
                    counts(jj) = sum(S.t_ms >= resp_ms(1) & S.t_ms <= resp_ms(2));
                end
            end
            fr_trials = counts / resp_s;  % sp/s
            FR_mean(ch, k) = mean(fr_trials, 'omitnan');
            maxFR = max(maxFR, FR_mean(ch,k));
        end
    end
    
    % ---- Plot heatmap for this stim set ----
    figure('Color','w','Name',sprintf('Stim set %s | FR heatmap', setLabel), ...
           'NumberTitle','off');
    imagesc(FR_mean);
    colormap(jet);
    caxis([0, maxFR]);     % normalize to this set
    cb = colorbar; ylabel(cb,'Firing rate (sp/s)');
    set(gca,'YDir','normal');
    
    % axis labels
    yticks(1:nChn);
    yticklabels(arrayfun(@(c) sprintf('Ch%d',c), 1:nChn, 'UniformOutput',false));
    xticks(1:n_AMP);
    xticklabels(arrayfun(@(a) sprintf('%g µA', a), Amps, 'UniformOutput',false));
    xlabel('Amplitude');
    ylabel('Channel');
    title(sprintf('Firing rate heatmap | Stim set: %s', setLabel), 'Interpreter','none');
end