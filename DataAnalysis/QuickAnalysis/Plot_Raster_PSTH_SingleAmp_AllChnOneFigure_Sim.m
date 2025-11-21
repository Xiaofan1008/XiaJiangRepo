clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder

data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim2'; 
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Single3_251106_131759';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq1';

if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

% Extract file name
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);  % 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback if no underscores
end

%% Choice
Spike_filtering = 0;
raster_chn_start = 1;
raster_chn_end = 32; %nChn
Electrode_Type = 1; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% Parameters for plotting
ras_win       = [-10 20];   % ms, time window to plot (you can change this)
bin_ms_raster = 1;          % bin size for PSTH (ms)
smooth_ms     = 3;          % smoothing window (ms)
Plot_Amps     = [5];        % amplitudes (µA) to plot; one figure per amp & stim set

%% Pre Set
FS=30000; % Sampling frequency
% Load .sp.mat file
% sp_files = dir('*.sp.mat');
sp_files = dir(fullfile(data_folder, '*.sp.mat'));
assert(~isempty(sp_files), 'No .sp.mat file found in the current folder.');
% sp_filename = sp_files(1).name;
sp_filename = fullfile(data_folder, sp_files(1).name);
fprintf('Loading spike file: %s\n', sp_filename);
S = load(sp_filename);
if isfield(S, 'sp')
    sp = S.sp;
else
    error('Variable "sp" not found in %s.', sp_filename);
end

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0); 

%% Spike Amplitude Filtering Parameters
pos_limit = 100;    % upper bound (µV)
neg_limit = -100;  % lower bound (µV)

baseline_window_ms = [-60, -5];        % Baseline window (ms)
response_window_ms = [2, 15];          % Response window (ms)

%% Waveform templete filtering parameters
template_window_ms = [0 10];   % use 0–5 ms spikes after trigger to build template
baseline_window_ms = [-60 0];    % window to filter
corr_thresh        = 0.70;       % threshold (increase to be stricter)

%% Load StimParams and decode amplitudes, stimulation sets, ISI
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
% S = load(fileDIR(1).name,'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
S = load(fullfile(data_folder, fileDIR(1).name), 'StimParams', 'simultaneous_stim', 'CHN', 'E_MAP', 'n_Trials');
StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
CHN                = S.CHN;
E_MAP              = S.E_MAP;
n_Trials           = S.n_Trials;

% Trial amplitude list
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps, ~, ampIdx] = unique(trialAmps(:));
if any(Amps == -1), Amps(Amps == -1) = 0; end
n_AMP = numel(Amps);
cmap = lines(n_AMP);  % color map for amplitudes

% Stimulation set decoding
E_NAME = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_NAME);
stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = unique(idx_all(rr)); v = v(v > 0).';
    stimChPerTrial_all{t} = v;
end
comb = zeros(n_Trials, simultaneous_stim);
for i = 1:n_Trials
    v = stimChPerTrial_all{i};  
    comb(i,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb, 'rows');
combClass_win = combClass;
nSets = size(uniqueComb,1);

% Pulse Train Period (inter-pulse interval)
pulseTrain_all = cell2mat(StimParams(2:end,9));  % Column 9: Pulse Train Period
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);  % take 1 per trial

[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

%% Load Significant Channel File
% File should be named:  <base_name>_FR_SigCh_ByAmp.mat
sig_file = fullfile(data_folder, sprintf('%s_FR_SigCh_ByAmp.mat', base_name));
hasSigFile = false;
if exist(sig_file, 'file')
    fprintf('Loading significant channel file:\n%s\n', sig_file);
    load(sig_file, 'FR_summary', 'Amps', 'responsive_channels');
    hasSigFile = true;
else
    fprintf('No significant-channel file found.\n');
end

%% Electrode Map
d = Depth_s(Electrode_Type); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

%% Spike Amplitude Filtering (Before Plotting)
% % sp_clipped = sp;   % copy original spike structure

% if Spike_filtering == 1
% 
%     fprintf('\n===== Spike Amplitude Filtering =====\n');
%     fprintf('Criteria: max ≤ %.1f µV and min ≥ %.1f µV\n', pos_limit, neg_limit);
% 
%     for ch = 1:numel(sp)
%         if isempty(sp{ch}), continue; end
% 
%         waveforms = sp{ch}(:, 2:end);
% 
%         % Amplitude-based filtering
%         max_vals = max(waveforms, [], 2);
%         min_vals = min(waveforms, [], 2);
%         in_range = (max_vals <= pos_limit) & (min_vals >= neg_limit);
% 
%         % Zero-crossing check (must include both positive and negative values)
%         has_positive = any(waveforms > 0, 2);
%         has_negative = any(waveforms < 0, 2);
%         crosses_zero = has_positive & has_negative;
% 
%         % Keep spikes within amplitude bounds
%         % valid_idx = (max_vals <= pos_limit) & (min_vals >= neg_limit);
%         % Final combined condition
%         valid_idx = in_range & crosses_zero;
% 
%         % Apply filter
%         sp_clipped{ch} = sp{ch}(valid_idx, :);
% 
%         % Print summary
%         n_total = size(sp{ch}, 1);
%         n_keep  = sum(valid_idx);
%         fprintf('Channel %2d: kept %4d / %4d spikes (%.1f%%)\n', ...
%             ch, n_keep, n_total, 100*n_keep/n_total);
%     end
% 
%     fprintf('=====================================\n\n');
%     save([base_name '.sp_xia.mat'],'sp_clipped');
% else
%     load([base_name '.sp_xia.mat']);
% 
%     % load([base_name '.sp.mat']);
%     % sp_clipped = sp;
% 
%     % load([base_name '.sp_xia_FirstPulse.mat']);
%     % sp_clipped = sp_seq;
% end

if Spike_filtering == 1
    fprintf('\nSpike Waveform Consistency Filtering (SSD-based)\n'); 

    SSD_threshold_factor = 10;  % from Allison-Walker (2022)
    t_axis = (0:48) / FS * 1000;

    for ch = 1:numel(sp)
        if isempty(sp{ch}), continue; end

        waveforms = sp{ch}(:, 2:end);   % exclude timestamps
        nSpikes = size(waveforms, 1);
        if nSpikes < 5
            sp_clipped{ch} = sp{ch};
            continue;
        end

        % --- Compute mean spike waveform ---
        mean_wave = mean(waveforms, 1);

        % --- Compute sum of squared differences (SSD) for each spike ---
        diff_mat = waveforms - mean_wave;      % deviation matrix
        SSD = sum(diff_mat.^2, 2);             % one SSD value per spike
        mean_SSD = mean(SSD);

        % --- Reject spikes whose SSD exceeds threshold ---
        valid_idx = SSD <= (SSD_threshold_factor * mean_SSD);

        % --- Apply filter ---
        sp_clipped{ch} = sp{ch}(valid_idx, :);

        % --- Summary ---
        n_total = size(sp{ch}, 1);
        n_keep  = sum(valid_idx);
        fprintf('Channel %2d: kept %4d / %4d spikes (%.1f%%)\n', ...
            ch, n_keep, n_total, 100*n_keep/n_total);
    end
    
   %% Waveform Correlation Template + Trial-by-Trial Baseline Filtering
    fprintf('\nWaveform Correlation Filtering\n')

    for ch = 1:numel(sp_clipped)   
        if isempty(sp_clipped{ch}), continue; end   
        waveforms = sp_clipped{ch}(:,2:end);
        sp_times  = sp_clipped{ch}(:,1);
        nSpikes   = size(waveforms,1);    
        if nSpikes < 10
            continue;
        end    
        % --- 1. Extract evoked spikes to build template ---
        evoked_mask = false(nSpikes,1);
    
        for tr = 1:length(trig)
            t0 = trig(tr)/FS*1000;
            evoked_mask = evoked_mask | ...
                (sp_times >= t0 + template_window_ms(1) & ...
                 sp_times <= t0 + template_window_ms(2));
        end    
        evoked_waves = waveforms(evoked_mask,:);
        if size(evoked_waves,1) < 5
            fprintf('Ch %d: not enough evoked spikes for template\n', ch);
            continue;
        end    
        template = mean(evoked_waves,1);    
        % --- 2. Compute correlation of every spike vs template ---
        corr_vals = zeros(nSpikes,1);
        for i = 1:nSpikes
            corr_vals(i) = corr(template(:), waveforms(i,:)');
        end
    
        % --- 3. TRIAL-BY-TRIAL FILTERING (baseline only) ---
        final_keep = true(nSpikes,1);   % keep everything unless filtered
    
        for tr = 1:length(trig)
            t0 = trig(tr)/FS*1000;
            % spikes in baseline window for this trial
            baseline_mask = (sp_times >= t0 + baseline_window_ms(1)) & ...
                            (sp_times <  t0 + baseline_window_ms(2));    
            idx_trial = find(baseline_mask);    
            if isempty(idx_trial), continue; end  
            % evaluate correlation only for these spikes
            bad_idx = idx_trial(corr_vals(idx_trial) < corr_thresh);  
            % remove them
            final_keep(bad_idx) = false;
        end
    
        %% --- 4. Keep evoked spikes + good baseline spikes ---
        beforeN = nSpikes;
        sp_clipped{ch} = sp_clipped{ch}(final_keep,:);
        afterN  = sum(final_keep);
    
        % fprintf('Ch %2d: kept %4d / %4d spikes (baseline filtered, corr >= %.2f)\n', ...
        %     ch, afterN, beforeN, corr_thresh);
    end    
    fprintf('Correlation Filtering Complete\n');
    save([base_name '.sp_xia.mat'], 'sp_clipped');
else
    load([base_name '.sp_xia.mat']);

    % load([base_name '.sp.mat']);
    % sp_clipped = sp;

    % load([base_name '.sp_xia_FirstPulse.mat']);
    % sp_clipped = sp_seq;
end

% if Spike_filtering == 1
%     fprintf('\n===== Spike Amplitude + Zero-Crossing + SSD Filtering =====\n');
% 
%     SSD_threshold_factor = 16;   % from Allison-Walker (2022)
%     t_axis = (0:48) / FS * 1000; % waveform sample timing (ms)
% 
%     for ch = 1:numel(sp)
%         if isempty(sp{ch}), continue; end
% 
%         waveforms = sp{ch}(:, 2:end);   % exclude timestamps
%         nSpikes = size(waveforms, 1);
%         if nSpikes < 5
%             sp_clipped{ch} = sp{ch};
%             continue;
%         end
% 
%         %% === 1. Amplitude-based filtering ===
%         max_vals = max(waveforms, [], 2);
%         min_vals = min(waveforms, [], 2);
%         in_range = (max_vals <= pos_limit) & (min_vals >= neg_limit);
% 
%         %% === 2. Zero-crossing check ===
%         has_positive = any(waveforms > 0, 2);
%         has_negative = any(waveforms < 0, 2);
%         crosses_zero = has_positive & has_negative;
% 
%         %% === 3. Combine amplitude + zero-crossing filters ===
%         valid_idx_lvl1 = in_range & crosses_zero;
%         waveforms_lvl1 = waveforms(valid_idx_lvl1, :);
% 
%         %% === 4. SSD-based waveform consistency filter ===
%         if size(waveforms_lvl1, 1) > 5
%             % Compute mean waveform
%             mean_wave = mean(waveforms_lvl1, 1);
% 
%             % Compute sum of squared differences for each spike
%             diff_mat = waveforms_lvl1 - mean_wave;
%             SSD = sum(diff_mat.^2, 2);
%             mean_SSD = mean(SSD);
% 
%             % Keep spikes whose SSD ≤ 16×mean SSD
%             valid_idx_lvl2 = SSD <= (SSD_threshold_factor * mean_SSD);
% 
%             % Combine both levels
%             final_valid = valid_idx_lvl1;
%             final_valid(valid_idx_lvl1) = valid_idx_lvl2;
%         else
%             % Too few spikes for SSD calculation
%             final_valid = valid_idx_lvl1;
%         end
% 
%         %% === 5. Apply combined filters ===
%         sp_clipped{ch} = sp{ch}(final_valid, :);
% 
%         %% === 6. Summary printout ===
%         n_total = size(sp{ch}, 1);
%         n_keep  = sum(final_valid);
%         fprintf('Channel %2d: kept %4d / %4d spikes (%.1f%%)\n', ...
%             ch, n_keep, n_total, 100 * n_keep / n_total);
%     end
% 
%     fprintf('=====================================\n\n');
%     save([base_name '.sp_xia.mat'], 'sp_clipped');
% 
% else
%     load([base_name '.sp_xia.mat']);
% 
%     % load([base_name '.sp.mat']);
%     % sp_clipped = sp;
% 
%     % load([base_name '.sp_xia_FirstPulse.mat']);
%     % sp_clipped = sp_seq;
% end
 %% ================== Raster + PSTH (All Channels, Grid) ==================

% --------- PSTH kernel and bin edges ---------
edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
g = g / sum(g);

% --------- Channel list for plotting ---------
ch_list = raster_chn_start:raster_chn_end;
nChPlot = numel(ch_list);

% ===== Loop over stimulation sets =====
for set_id = 1:nSets
    
    stimIdx = uniqueComb(set_id, :);
    stimIdx = stimIdx(stimIdx > 0);
    setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx, ...
                        'UniformOutput', false), ' + ');
    
    % ===== Loop over amplitudes to plot =====
    for ai = 1:numel(Plot_Amps)
        
        amp_val = Plot_Amps(ai);
        amp_idx_match = find(abs(Amps - amp_val) < 1e-6);
        if isempty(amp_idx_match)
            fprintf('%.1f µA not found in this dataset — skipping.\n', amp_val);
            continue;
        end
        
        % Trials for this stim set & amplitude (all pulse periods collapsed)
        trial_mask = (combClass_win == set_id) & (ampIdx == amp_idx_match);
        trial_ids  = find(trial_mask);
        nTr        = numel(trial_ids);
        
        if nTr == 0
            fprintf('Set %d (%s) | %.1f µA: no trials → skipping.\n', ...
                    set_id, setLabel, amp_val);
            continue;
        end

        % -------- Determine significant channels for this stim set & amplitude --------
        sigChThis = [];
        if hasSigFile
            try
                % --- use amplitude list inside FR_summary for this set ---
                amps_set = FR_summary(set_id).amps;       
                % find which column corresponds to this amplitude
                ai_sig = find(abs(amps_set - amp_val) < 1e-6);        
                if ~isempty(ai_sig)
                    % get significant channels for THIS set & THIS amplitude
                    sigChThis = find(FR_summary(set_id).sig(:, ai_sig));
                end
            catch
                sigChThis = [];
            end
        end
        
        fprintf('Set %d (%s) | %.0f µA: %d significant channels detected.\n', ...
                set_id, setLabel, amp_val, numel(sigChThis));
        
        %% ===== Create figure for (set, amplitude) =====
        % Grid for channels: roughly square
        nCols = ceil(sqrt(nChPlot));
        nRows = ceil(nChPlot / nCols);
        
        figTitle = sprintf('Stim %s | %.0f µA | Simultaneous Stimualtion', setLabel, amp_val);
        figure('Color','w', 'Name', figTitle,'Position',[50 50 1600 900]);
        sgtitle(figTitle, 'FontSize', 16, 'FontWeight', 'bold');
        
        tl = tiledlayout(nRows, nCols, ...
            'TileSpacing','compact', 'Padding','compact');
        title(tl, sprintf('Stim %s | %.0f µA | Simultaneous Stimualtion',setLabel, amp_val), 'FontSize', 12, 'FontWeight', 'bold');
        
        %% Loop over channels in this figure
        for ci = 1:nChPlot         
            recCh = ch_list(ci);   % recording channel index
            ch    = d(recCh);      % index into sp_clipped            
            if ch < 1 || ch > numel(sp_clipped)
                continue;
            end
            if isempty(sp_clipped{ch})
                nexttile; axis off;
                continue;
            end            
            S_ch = sp_clipped{ch};           
            % --- Collect spikes per trial ---
            allTrialSpikes = cell(nTr,1);
            counts = zeros(1, numel(edges)-1);
            
            for ti = 1:nTr
                tr = trial_ids(ti);
                t0 = trig(tr) / FS * 1000;   % ms
                tt = S_ch(:,1);
                % spikes in plotting window relative to trigger
                tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                allTrialSpikes{ti} = tt;
                counts = counts + histcounts(tt, edges);
            end         
            % --- Compute PSTH ---
            if nTr > 0
                rate   = counts / (nTr * bin_s);   % sp/s
                rate_s = filter(g, 1, rate);       % smoothed
            else
                rate_s = zeros(1, numel(ctrs));
            end
            maxRate = max(rate_s);
            % Make sure PSTH axis always goes at least up to 50 sp/s
            yMaxPSTH = max(100, ceil(maxRate * 1.1 / 10)*10);
            
            % ===== Plot in a single subplot: PSTH (left) + Raster (right) =====
            ax = nexttile(tl);
            if ismember(recCh, sigChThis)
                ax.LineWidth = 2;
                ax.XColor = [0.8 0 0];
                ax.YColor = [0.8 0 0];              
            end
            hold(ax, 'on');           
            % --- PSTH on LEFT y-axis ---
            yyaxis(ax, 'left');
            if any(rate_s)
                plot(ax, ctrs, rate_s, 'LineWidth', 1.6);
            end
            xlim(ax, ras_win);
            ylim(ax, [0 yMaxPSTH]);
            ylabel(ax, 'Firing rate (sp/s)');          
            % --- Raster on RIGHT y-axis ---
            yyaxis(ax, 'right');
            for ti = 1:nTr
                tt = allTrialSpikes{ti};
                if isempty(tt), continue; end
            
                % --- small dot raster ---
                plot(ax, tt, ti * ones(size(tt)), '.', ...
                    'Color', [0 0 0], 'MarkerSize', 4);  
            end
            ylim(ax, [0 nTr+1]);
            set(ax, 'YTick', []);           
            % --- Shared x-axis and trigger line ---
            xline(ax, 0, 'r--', 'LineWidth', 1, 'HandleVisibility','off');
            xlim(ax, ras_win);          
            % --- Labels / title ---
            if ismember(recCh, sigChThis)
                chTitle = sprintf('Ch %d  (SIG)', recCh);
            else
                chTitle = sprintf('Ch %d', recCh);
            end
            title(ax, chTitle, 'FontSize', 12, 'FontWeight','bold');
            % title(ax, sprintf('Ch %d', recCh), 'FontSize', 12, 'FontWeight','bold');
            if ci > (nRows-1)*nCols
                xlabel(ax, 'Time (ms)');
            end
            
        end % channel loop
        
    end % amplitude loop
end % stim set loop

fprintf('\nRaster + PSTH plotting finished.\n');