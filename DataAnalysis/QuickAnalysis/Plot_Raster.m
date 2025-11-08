clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder

data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim6'; 
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single2_251104_122817';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Seq5_New_251014_194221';

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
Spike_filtering =1;


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

% Electrode Map
d = Depth_s(1); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

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

% if Spike_filtering == 1
%     fprintf('\n===== Spike Waveform Consistency Filtering (SSD-based) =====\n'); 
% 
%     SSD_threshold_factor = 16;  % from Allison-Walker (2022)
%     t_axis = (0:48) / FS * 1000;
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
%         % --- Compute mean spike waveform for this channel ---
%         mean_wave = mean(waveforms, 1);
% 
%         % --- Compute sum of squared differences (SSD) for each spike ---
%         diff_mat = waveforms - mean_wave;      % deviation matrix
%         SSD = sum(diff_mat.^2, 2);             % one SSD value per spike
%         mean_SSD = mean(SSD);
% 
%         % --- Reject spikes whose SSD exceeds threshold ---
%         valid_idx = SSD <= (SSD_threshold_factor * mean_SSD);
% 
%         % --- Apply filter ---
%         sp_clipped{ch} = sp{ch}(valid_idx, :);
% 
%         % --- Summary ---
%         n_total = size(sp{ch}, 1);
%         n_keep  = sum(valid_idx);
%         fprintf('Channel %2d: kept %4d / %4d spikes (%.1f%%)\n', ...
%             ch, n_keep, n_total, 100*n_keep/n_total);
%     end
% 
%     fprintf('=====================================\n\n');
%     save([base_name '.sp_xia.mat'], 'sp_clipped');
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
    fprintf('\n===== Spike Amplitude + Zero-Crossing + SSD Filtering =====\n');

    SSD_threshold_factor = 16;   % from Allison-Walker (2022)
    t_axis = (0:48) / FS * 1000; % waveform sample timing (ms)

    for ch = 1:numel(sp)
        if isempty(sp{ch}), continue; end

        waveforms = sp{ch}(:, 2:end);   % exclude timestamps
        nSpikes = size(waveforms, 1);
        if nSpikes < 5
            sp_clipped{ch} = sp{ch};
            continue;
        end

        %% === 1. Amplitude-based filtering ===
        max_vals = max(waveforms, [], 2);
        min_vals = min(waveforms, [], 2);
        in_range = (max_vals <= pos_limit) & (min_vals >= neg_limit);

        %% === 2. Zero-crossing check ===
        has_positive = any(waveforms > 0, 2);
        has_negative = any(waveforms < 0, 2);
        crosses_zero = has_positive & has_negative;

        %% === 3. Combine amplitude + zero-crossing filters ===
        valid_idx_lvl1 = in_range & crosses_zero;
        waveforms_lvl1 = waveforms(valid_idx_lvl1, :);

        %% === 4. SSD-based waveform consistency filter ===
        if size(waveforms_lvl1, 1) > 5
            % Compute mean waveform
            mean_wave = mean(waveforms_lvl1, 1);

            % Compute sum of squared differences for each spike
            diff_mat = waveforms_lvl1 - mean_wave;
            SSD = sum(diff_mat.^2, 2);
            mean_SSD = mean(SSD);

            % Keep spikes whose SSD ≤ 16×mean SSD
            valid_idx_lvl2 = SSD <= (SSD_threshold_factor * mean_SSD);

            % Combine both levels
            final_valid = valid_idx_lvl1;
            final_valid(valid_idx_lvl1) = valid_idx_lvl2;
        else
            % Too few spikes for SSD calculation
            final_valid = valid_idx_lvl1;
        end

        %% === 5. Apply combined filters ===
        sp_clipped{ch} = sp{ch}(final_valid, :);

        %% === 6. Summary printout ===
        n_total = size(sp{ch}, 1);
        n_keep  = sum(final_valid);
        fprintf('Channel %2d: kept %4d / %4d spikes (%.1f%%)\n', ...
            ch, n_keep, n_total, 100 * n_keep / n_total);
    end

    fprintf('=====================================\n\n');
    save([base_name '.sp_xia.mat'], 'sp_clipped');

else
    load([base_name '.sp_xia.mat']);

    % load([base_name '.sp.mat']);
    % sp_clipped = sp;

    % load([base_name '.sp_xia_FirstPulse.mat']);
    % sp_clipped = sp_seq;
end

%% Raster Plot Parameters
ras_win         = [-20 100];   % ms
bin_ms_raster   = 2;           % bin size
smooth_ms       = 3;           % smoothing window
raster_chn_start = 1;
raster_chn_end = 32; %nChn


%% Raster Plot 
 edges = ras_win(1):bin_ms_raster:ras_win(2);
    ctrs  = edges(1:end-1) + diff(edges)/2;
    bin_s = bin_ms_raster / 1000;
    g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
    g = g / sum(g);

    for ich = raster_chn_start:raster_chn_end
        ch = d(ich);
        if isempty(sp_clipped{ch}), continue; end

        for si = 1:nSets
            set_id = si;
            stimIdx = uniqueComb(set_id, :); stimIdx = stimIdx(stimIdx > 0);
            setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx, 'UniformOutput', false), ' + ');

            for pi = 1:n_PULSE
                pulse_val = PulsePeriods(pi);
                trials_this_period = find(pulseIdx == pi & combClass_win == set_id);
                if isempty(trials_this_period), continue; end

                figure('Color','w','Name',sprintf('Ch %d | Set %s | %d µs', ich, setLabel, pulse_val));
                tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
                ax1 = nexttile([3 1]); hold(ax1,'on'); box(ax1,'off');
                ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');

                psth_curves = cell(1, n_AMP);
                maxRate = 0;
                y_cursor = 0;
                ytick_vals = [];
                ytick_labels = {};
                total_spikes = 0;

                for ai = 1:n_AMP
                    amp_val = Amps(ai);
                    color = cmap(ai,:);
                    amp_trials = find(ampIdx == ai & pulseIdx == pi & combClass_win == set_id);
                    nTr = numel(amp_trials);
                    spkTimesPerTrial = cell(nTr,1);
                    counts = zeros(1, numel(edges)-1);
                    total_spikes_amp = 0;

                    for t = 1:nTr
                        tr = amp_trials(t);
                        S_ch = sp_clipped{ch};
                        t0 = trig(tr)/FS*1000;
                        tt = S_ch(:,1);
                        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                        spkTimesPerTrial{t} = tt;
                        counts = counts + histcounts(tt, edges);
                        total_spikes_amp = total_spikes_amp + numel(tt);
                    end
                    total_spikes = total_spikes + total_spikes_amp;

                    for t = 1:nTr
                        tt = spkTimesPerTrial{t};
                        y0 = y_cursor + t;
                        for k = 1:numel(tt)
                            plot(ax1, [tt(k) tt(k)], [y0-0.4 y0+0.4], 'Color', color, 'LineWidth', 1.2);
                        end
                    end

                    ytick_vals(end+1) = y_cursor + max(nTr, 1)/2;
                    ytick_labels{end+1} = sprintf('%.0f µA', amp_val);

                    % Horizontal divider
                    if nTr > 0
                        plot(ax1, ras_win, [y_cursor+nTr y_cursor+nTr], 'Color', [0.8 0.8 0.8], 'LineStyle','--');
                    end

                    % ---- Add spike count label ----
                    if nTr > 0
                        text(ax1, ras_win(2)+2, y_cursor + nTr/2, ...
                             sprintf('%d spikes', total_spikes_amp), ...
                             'Color', color, 'FontSize', 9, 'HorizontalAlignment','left');
                    end

                    y_cursor = y_cursor + max(nTr,1);

                    if nTr > 0
                        rate = counts / (nTr * bin_s);
                        rate_s = filter(g, 1, rate);
                        psth_curves{ai} = rate_s;
                        maxRate = max(maxRate, max(rate_s));
                    else
                        psth_curves{ai} = zeros(1, numel(ctrs));
                    end

                    % fprintf('Ch %2d | Set %s | %4d µs | %3d µA — Total spikes: %d\n', ...
                    %     ich, setLabel, pulse_val, amp_val, total_spikes_amp);
                end
                fprintf('Ch %2d | Set %s | %4d µs — Total spikes: %d\n', ...
                        ich, setLabel, pulse_val, total_spikes);

                % Finalize raster plot
                xline(ax1, 0, 'r', 'LineWidth', 1.5);
                xlim(ax1, [ras_win(1), ras_win(2) + 15]);  % extra space for spike count label
                ylim(ax1, [0, y_cursor]);
                yticks(ax1, ytick_vals);
                yticklabels(ax1, ytick_labels);
                ylabel(ax1, 'Amplitude');
                title(ax1, sprintf('Raster — Ch %d | Set: %s | %d µs', ich, setLabel, pulse_val), 'Interpreter','none');

                % PSTH
                for ai = 1:n_AMP
                    plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.5);
                end
                xline(ax2, 0, 'r', 'LineWidth', 1.5);
                xlim(ax2, ras_win);
                if maxRate > 0
                    ylim(ax2, [0 ceil(maxRate*1.1/10)*10]);
                else
                    ylim(ax2, [0 1]);
                end
                xlabel(ax2, 'Time (ms)');
                ylabel(ax2, 'Rate (sp/s)');
                legend(ax2, arrayfun(@(a) sprintf('%.0f µA', a), Amps, 'UniformOutput', false), ...
                    'Box','off','Location','northeast');
            end
        end
    end