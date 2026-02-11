clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder

% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Single5_251014_184742'; 
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Sim5_251014_183532';
data_folder = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_1';

%% Choice
Spike_filtering = 1;
raster_chn_start = 35;
raster_chn_end = 40; %nChn
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex
PTD_to_plot = [0 5 10 12 15 17];   % e.g., [500 1000], empty for all PTD
PTD_to_plot = PTD_to_plot.*1000;
%% Spike Amplitude Filtering Parameters
pos_limit = 100;    % upper bound (µV)
neg_limit = -100;  % lower bound (µV)

baseline_window_ms = [-60, -5];        % Baseline window (ms)
response_window_ms = [2, 25];          % Response window (ms)

%% Waveform templete filtering parameters
template_window_ms = [5 15];   % use 0–5 ms spikes after trigger to build template
baseline_window_ms = [-60 0];    % window to filter
corr_thresh        = 0.70;       % threshold (increase to be stricter)

%% Load folders
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

% %Stim channel without order
% stimChPerTrial_all = cell(n_Trials,1);
% for t = 1:n_Trials
%     rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
%     v = unique(idx_all(rr)); v = v(v > 0).';
%     stimChPerTrial_all{t} = v;
% end
% comb = zeros(n_Trials, simultaneous_stim);
% for i = 1:n_Trials
%     v = stimChPerTrial_all{i};  
%     comb(i,1:numel(v)) = v;
% end
% [uniqueComb, ~, combClass] = unique(comb, 'rows');
% combClass_win = combClass;
% nSets = size(uniqueComb,1);

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    % KEEP ORDER — do NOT use unique()
    v = idx_all(rr);
    % Remove invalid zeros only, keep order
    v = v(v > 0);
    stimChPerTrial_all{t} = v(:)';   % row vector
end
% Build combination matrix preserving order
comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t, 1:numel(v)) = v;
end
% Now each distinct ordered sequence becomes a different set
[uniqueComb, ~, combClass] = unique(comb, 'rows', 'stable');
nSets = size(uniqueComb,1);
combClass_win = combClass;

% Pulse Train Period (inter-pulse interval)
pulseTrain_all = cell2mat(StimParams(2:end,9));  % Column 9: Pulse Train Period
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);  % take 1 per trial
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

% Extract POST-TRIGGER DELAY (PTD)
% Column 6 stores PTD for the 2nd pulse of each trial block
% For sequential stimulation: row 2 = delayed pulse
if simultaneous_stim > 1
    ptd_all = cell2mat(StimParams(3:simultaneous_stim:end, 6));  % µs
else
    ptd_all = zeros(n_Trials,1); % Single-pulse case, no PTD
end

PTD_us = ptd_all(:);
[PTD_values, ~, ptdIdx] = unique(PTD_us);

%% ============================================================
%   BLANK SPIKES BETWEEN 0 ms AND PTD FOR EACH TRIAL
%   This removes all spikes occurring between:
%       trigger_time  -->  trigger_time + PTD (in ms)
if Spike_filtering ==1
    fprintf('\nBlanking spikes from 0 → PTD for each trial...\n');
    
    sp_blank = sp;   % start from raw spike timestamps
    
    for ch = 1:numel(sp_blank)
        if isempty(sp_blank{ch}), continue; end
    
        spike_times = sp_blank{ch}(:,1);    % in ms
    
        % build a logical mask of spikes to keep
        keep_mask = true(size(spike_times));
    
        % loop through trials
        for tr = 1:length(trig)
            t0_ms   = trig(tr) / FS * 1000;     % trigger time in ms
            PTD_ms  = PTD_us(tr) / 1000;        % convert µs → ms
    
            % blank all spikes between [t0, t0 + PTD]
            blank_mask = spike_times >= t0_ms & spike_times < (t0_ms + PTD_ms);
            keep_mask(blank_mask) = false;
        end
    
        % apply blanking
        sp_blank{ch} = sp_blank{ch}(keep_mask,:);
    end
    
    fprintf('Spike blanking complete.\n');
    
    % from here on, the pipeline uses sp_blank instead of sp or sp_clipped
    sp = sp_blank;
end 
if isempty(PTD_to_plot)
    PTD_selected = PTD_values;
else
    PTD_selected = intersect(PTD_values, PTD_to_plot);
end
fprintf('\nDetected PTDs (µs):'); disp(PTD_values');
fprintf('PTDs selected for plotting:'); disp(PTD_selected');
n_PTD = numel(PTD_selected);
% Electrode Map
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
%     % load([base_name '.sp_xia.mat']);
% 
%     % load([base_name '.sp.mat']);
%     % sp_clipped = sp;
% 
%     load([base_name '.sp_xia_FirstPulse.mat']);
%     sp_clipped = sp_seq;
% end

% load([base_name '.sp_xia_FirstPulse.mat']);
% sp = sp_seq;
if Spike_filtering == 1
    fprintf('\nSpike Waveform Consistency Filtering (SSD-based)\n'); 

    SSD_threshold_factor = 16;  % from Allison-Walker (2022)
    t_axis = (0:48) / FS * 1000;

    for ch = 1:numel(sp)
        if isempty(sp{ch}), continue; end

        waveforms = sp{ch}(:, 2:end);   % exclude timestamps
        nSpikes = size(waveforms, 1);
        if nSpikes < 5
            sp_clipped{ch} = sp{ch};
            continue;
        end

        % --- Compute mean spike waveform for this channel ---
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
        %%  Waveform Correlation Template + Trial-by-Trial Baseline Filtering
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
        
            % --- 4. Keep evoked spikes + good baseline spikes ---
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
%     % load([base_name '.sp_xia.mat']);
% 
%     % load([base_name '.sp.mat']);
%     % sp_clipped = sp;
% 
%     load([base_name '.sp_xia_FirstPulse.mat']);
%     sp_clipped = sp_seq;
% end


%% Raster Plot Parameters
ras_win         = [-20 100];   % ms
bin_ms_raster   = 1;           % bin size
smooth_ms       = 2;           % smoothing window
% raster_chn_start = 1;
% raster_chn_end = 32; %nChn

%% === Initialize structure to store first-spike times ===
firstSpikeTimes = cell(raster_chn_end, 1); % each cell: vector of first-spike times per trial (ms)
fprintf('\nComputing First Spike Times per Trial\n');
post_spike_window_ms = [5,8];

%% Raster Plot 
 %% ========================= RASTER + PSTH ========================= %%
edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2); 
g = g / sum(g);

for ich = raster_chn_start:raster_chn_end
    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end

    for si = 1:nSets
        stimVec = uniqueComb(si, :);
        stimVec = stimVec(stimVec > 0);
        setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimVec, 'UniformOutput', false), '→');

        for pi = 1:n_PULSE
            pulse_val = PulsePeriods(pi);

            for i_PTD = 1:n_PTD
                ptd_val = PTD_selected(i_PTD);

                % --------- FIXED: restrict to THIS PTD ---------
                trials_this_period = find( combClass_win == si & ...
                                           pulseIdx == pi & ...
                                           ampIdx >= 1 & ...   % any amp
                                           PTD_us == ptd_val );  % << FIXED
                if isempty(trials_this_period), continue; end

                % --------- FIGURE ---------
                figName = sprintf('Ch %d | Set %s | Pulse %d µs | PTD %d µs', ...
                                  ich, setLabel, pulse_val, ptd_val);
                figure('Color','w','Name',figName);
                tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                ax1 = nexttile([3 1]);
                hold(ax1,'on'); box(ax1,'off');
                title(ax1, sprintf('Raster — Ch %d | Set %s | PTD %d µs', ...
                                   ich, setLabel, ptd_val), 'Interpreter','none');

                ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');

                % ===== PSTH storage =====
                psth_curves = cell(1, n_AMP);
                maxRate = 0;
                y_cursor = 0;
                ytick_vals = [];
                ytick_labels = {};

                % ================= LOOP AMPLITUDES ==================
                for ai = 1:n_AMP
                    amp_val = Amps(ai);
                    color = cmap(ai,:);

                    % -------- FIXED TRIAL FILTERING --------
                    amp_trials = find( ampIdx == ai & ...
                                       pulseIdx == pi & ...
                                       combClass_win == si & ...
                                       PTD_us == ptd_val );   % << FIXED

                    nTr = numel(amp_trials);
                    if nTr == 0
                        psth_curves{ai} = zeros(1, numel(ctrs));
                        continue;
                    end

                    % ======== RASTER + PSTH COUNTING ========
                    counts = zeros(1, numel(edges)-1);
                    for t = 1:nTr
                        tr = amp_trials(t);
                        t0 = trig(tr)/FS*1000;

                        tt = sp_clipped{ch}(:,1);
                        tt = tt(tt >= t0+ras_win(1) & tt <= t0+ras_win(2)) - t0;

                        % ---- RASTER ----
                        y0 = y_cursor + t;
                        for spike_t = tt'
                            plot(ax1, [spike_t spike_t], [y0-0.4 y0+0.4], 'Color', color, 'LineWidth', 1.1);
                        end

                        % ---- PSTH ----
                        counts = counts + histcounts(tt, edges);
                    end

                    % y-axis structure
                    ytick_vals(end+1) = y_cursor + nTr/2;
                    ytick_labels{end+1} = sprintf('%d µA', amp_val);
                    y_cursor = y_cursor + nTr;

                    % ---- Compute PSTH ----
                    rate = filter(g, 1, counts/(nTr*bin_s));
                    psth_curves{ai} = rate;
                    maxRate = max(maxRate, max(rate));
                end

                % Finalize raster axis
                xline(ax1, 0, 'r--');
                xlim(ax1, ras_win);
                ylim(ax1, [0 y_cursor]);
                yticks(ax1, ytick_vals);
                yticklabels(ax1, ytick_labels);
                ylabel(ax1, 'Amplitude');

                % Finalize PSTH
                for ai = 1:n_AMP
                    plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.6);
                end
                xline(ax2, 0, 'r--');
                xlim(ax2, ras_win);
                ylim(ax2, [0 max(50,ceil(maxRate*1.1/10)*10)]);
                xlabel(ax2, 'Time (ms)');
                ylabel(ax2, 'Rate (sp/s)');
                legend(ax2, arrayfun(@(a) sprintf('%.0f µA',a), Amps, 'UniformOutput', false), ...
                       'Box','off','Location','northeast');

            end % PTD
        end % Pulse
    end % Set
end % Channel