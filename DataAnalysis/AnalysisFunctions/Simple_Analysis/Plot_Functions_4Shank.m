% clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder
% data_folder = uigetdir(pwd, 'Select the folder containing your data');
% if isequal(data_folder, 0)
%     error('No folder selected. Script aborted.');
% end
% cd(data_folder);  % Change to the selected folder

data_folder = '/Volumes/MACData/Data/Data_Xia/CJ268/Electrode_Test_Xia_251022_172613'; 
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Seq5_New_251014_194221';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Single5_251014_184742';
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

% Extract file name
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 2
    % Keep only the part before the last two underscore segments (assuming timestamp has two parts)
    base_name = last_folder(1 : underscores(end-1) - 1);  % 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback if no underscores
end

%% Choice
Spike_filtering = 0;
spike_plot_per_amp_StimSet = 1;
raster_psth_per_amp_StimSet = 1;
firing_rate_curve_per_amp_StimSet = 0;
Shank_to_plot =4;
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
% Load Trigger
% if isempty(dir('*.trig.dat'))
%     cleanTrig_sabquick;
% end
% trig = loadTrig(0);

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);  % if loadTrig handles current folder, no need to modify


%% Spike Amplitude Filtering Parameters
pos_limit = 80;    % upper bound (µV)
neg_limit = -100;  % lower bound (µV)

%% Spike Waveform Parameters
win_ms = 300;         % total time window after each trigger (ms)
bin_ms = 5;          % bin size (ms)
% nBins  = win_ms / bin_ms;
nBins  = 100 / bin_ms; % plot 100ms 
nChn   = numel(sp);
amp_threshold = 100;  % max allowed amplitude (µV)
spike_chn_start = 20;
spike_chn_end = 20; %nChn
% switch Shank_to_plot
%     case 1
%         spike_chn_start = 1;
%         spike_chn_end = 16;
%     case 2
%         spike_chn_start = 33;
%         spike_chn_end = 48;
%     case 3
%         spike_chn_start = 49;
%         spike_chn_end = 64;
%     case 4
%         spike_chn_start = 17;
%         spike_chn_end = 32;
% end

%% Raster Plot Parameters
ras_win         = [-20 300];   % ms
bin_ms_raster   = 2;           % bin size
smooth_ms       = 3;           % smoothing window
raster_chn_start = spike_chn_start;
raster_chn_end = spike_chn_end; %nChn

%% Firing Rate Curve Plot Parameters
channels_to_plot = [1,2,3,18,19,20,35,38,50,53] ;          % target channel
time_window_ms  = [2,20];    % Firing rate window in ms

%% Load StimParams and decode amplitudes & stimulation sets
% fileDIR = dir('*_exp_datafile_*.mat');
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

% Electrode Map
d = Depth_s(1); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

%% Spike Amplitude Filtering (Before Plotting)
sp_clipped = sp;   % copy original spike structure

if Spike_filtering == 1

    fprintf('\n===== Spike Amplitude Filtering =====\n');
    fprintf('Criteria: max ≤ %.1f µV and min ≥ %.1f µV\n', pos_limit, neg_limit);
    
    for ch = 1:numel(sp)
        if isempty(sp{ch}), continue; end
    
        waveforms = sp{ch}(:, 2:end);
    
        % Amplitude-based filtering
        max_vals = max(waveforms, [], 2);
        min_vals = min(waveforms, [], 2);
        in_range = (max_vals <= pos_limit) & (min_vals >= neg_limit);
    
        % Zero-crossing check (must include both positive and negative values)
        has_positive = any(waveforms > 0, 2);
        has_negative = any(waveforms < 0, 2);
        crosses_zero = has_positive & has_negative;
    
        % Keep spikes within amplitude bounds
        % valid_idx = (max_vals <= pos_limit) & (min_vals >= neg_limit);
        % Final combined condition
        valid_idx = in_range & crosses_zero;
        % Apply filter
        sp_clipped{ch} = sp{ch}(valid_idx, :);
    
        % Print summary
        n_total = size(sp{ch}, 1);
        n_keep  = sum(valid_idx);
        fprintf('Channel %2d: kept %4d / %4d spikes (%.1f%%)\n', ...
            ch, n_keep, n_total, 100*n_keep/n_total);
    end
    
    fprintf('=====================================\n\n');
    save([base_name '.sp_xia.mat'],'sp_clipped');
else
    load([base_name '.sp_xia.mat']);
end


if spike_plot_per_amp_StimSet == 1
    %% Spike plotting loop
    for ich = spike_chn_start:spike_chn_end %1:nChn
        ch = d(ich); % channel map to Intan
        if isempty(sp_clipped{ch}), continue; end
        fprintf('Processing Channel %d...\n', ich);
    
        % Spike data
        sp_times = sp_clipped{ch}(:,1);
        sp_wave  = sp_clipped{ch}(:,2:end);
        valid_idx = all(abs(sp_wave) <= amp_threshold, 2);
        sp_times = sp_times(valid_idx);
        sp_wave  = sp_wave(valid_idx,:);
        if isempty(sp_times), continue; end
    
        t_wave = (0:size(sp_wave,2)-1) / FS * 1000;
    
        for s = 1:nSets
            set_id = s;
            trial_mask = (combClass_win == set_id);
            if ~any(trial_mask), continue; end
    
            trial_ids = find(trial_mask);
            all_spikes_by_bin_amp = cell(nBins, n_AMP);
    
            for idx = 1:length(trial_ids)
                tr = trial_ids(idx);
                t0_ms = trig(tr) / FS * 1000;
                amp_id = ampIdx(tr);
                spk_mask = sp_times >= t0_ms & sp_times < (t0_ms + win_ms);
                if ~any(spk_mask), continue; end
    
                rel_times = sp_times(spk_mask) - t0_ms;
                waveforms = sp_wave(spk_mask,:);
    
                for j = 1:length(rel_times)
                    bin_idx = floor(rel_times(j) / bin_ms) + 1;
                    if bin_idx >= 1 && bin_idx <= nBins
                        all_spikes_by_bin_amp{bin_idx, amp_id}(end+1,:) = waveforms(j,:);
                    end
                end
            end
    
            % Determine shared Y-limit
            all_waves = cell2mat(all_spikes_by_bin_amp(:));
            if isempty(all_waves), continue; end
            y_max = max(abs(all_waves(:)));
            y_lim = [-1 1] * ceil(y_max/50)*50;
    
            % ==== PLOTTING ====
            stimIdx = uniqueComb(set_id,:); stimIdx = stimIdx(stimIdx > 0);
            stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx, 'UniformOutput', false), ', ');
            figure('Name', sprintf('Ch %d | StimSet %d (%s)', ich, set_id, stimLabel), ...
                   'Color','w','Position', [100 100 700 400]); %'Position', [100 100 1400 800]);
            tiledlayout(4, 5, 'Padding','compact', 'TileSpacing','compact');
    
            for b = 1:nBins
                nexttile; hold on;
                for a = 1:n_AMP
                    this_bin = all_spikes_by_bin_amp{b,a};
                    if isempty(this_bin), continue; end
    
                    % Align by min
                    aligned_waves = zeros(size(this_bin));
                    for k = 1:size(this_bin,1)
                        [~, min_idx] = min(this_bin(k,:));
                        shift = ceil(size(this_bin,2)/2) - min_idx;
                        aligned_waves(k,:) = circshift(this_bin(k,:), shift, 2);
                    end
    
                    plot(t_wave, aligned_waves', 'Color', [cmap(a,:) 0.3]);
                end
    
                % Count total spikes in this bin (across all amplitudes)
                spike_count = 0;
                for a = 1:n_AMP
                    spike_count = spike_count + size(all_spikes_by_bin_amp{b,a}, 1);
                end
                
                % Plot title includes spike count
                title(sprintf('%d–%d ms (%d spikes)', (b-1)*bin_ms, b*bin_ms, spike_count));
                xlabel('Time (ms)');
                ylabel('Voltage (µV)');
                ylim(y_lim);
                yticks(linspace(y_lim(1), y_lim(2), 3));
                xticks(round(linspace(t_wave(1), t_wave(end), 3)));
                axis square; grid on;
            end
    
            % Legend (global, separate)
            legend_entries = arrayfun(@(a) sprintf('%g µA', Amps(a)), 1:n_AMP, 'UniformOutput', false);
            legend_colors = arrayfun(@(a) plot(nan,nan,'-', 'Color', cmap(a,:), 'LineWidth',1.5), 1:n_AMP);
            legend(legend_colors, legend_entries, 'Location','northeastoutside');
        end
    end
end


 if raster_psth_per_amp_StimSet == 1    
    % Time axis
    edges = ras_win(1):bin_ms_raster:ras_win(2);
    ctrs  = edges(1:end-1) + diff(edges)/2;
    bin_s = bin_ms_raster / 1000;
    
    % Smoothing kernel
    g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
    g = g / sum(g);
    
    %% Main Loop
    for si = 1:nSets
        set_id = si;
        stimIdx = uniqueComb(set_id, :); stimIdx = stimIdx(stimIdx > 0);
        setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx, 'UniformOutput', false), ' + ');
    
        for ich = raster_chn_start:raster_chn_end %1:nChn
            ch = d(ich);
            if isempty(sp_clipped{ch}), continue; end
    
            figure('Color','w','Name',sprintf('Ch %d | Set %s', ich, setLabel));
            tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
            ax1 = nexttile([3 1]); hold(ax1,'on'); box(ax1,'off');
            ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');
    
            psth_curves = cell(1, n_AMP);
            maxRate = 0;
            y_cursor = 0;
            ytick_vals = [];
            ytick_labels = {};
    
            for ai = 1:n_AMP
                amp_val = Amps(ai);
                if amp_val == 0, continue; end   % SKIP 0 µA trials
                color = cmap(ai,:);
                amp_trials = find(ampIdx == ai & combClass_win == set_id);
                nTr = numel(amp_trials);
                spkTimesPerTrial = cell(nTr,1);
                counts = zeros(1, numel(edges)-1);
    
                for t = 1:nTr
                    tr = amp_trials(t);
                    S_ch = sp_clipped{ch};
                    t0 = trig(tr)/FS*1000;
                    tt = S_ch(:,1);
                    tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
                    spkTimesPerTrial{t} = tt;
                    counts = counts + histcounts(tt, edges);
                end
    
                for t = 1:nTr
                    tt = spkTimesPerTrial{t};
                    y0 = y_cursor + t;
                    for k = 1:numel(tt)
                        plot(ax1, [tt(k) tt(k)], [y0-0.4 y0+0.4], 'Color', color, 'LineWidth', 1.2);
                    end
                end
    
                ytick_vals(end+1) = y_cursor + max(nTr, 1)/2;
                ytick_labels{end+1} = sprintf('%.0f µA', amp_val);
                if nTr > 0
                    plot(ax1, ras_win, [y_cursor+nTr y_cursor+nTr], 'Color', [0.8 0.8 0.8], 'LineStyle','--');
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
            end
    
            xline(ax1, 0, 'r', 'LineWidth', 1.5);
            % xlim(ax1, ras_win);
            xlim(ax1, [-20,100]);
            ylim(ax1, [0, y_cursor]);
            yticks(ax1, ytick_vals);
            yticklabels(ax1, ytick_labels);
            ylabel(ax1, 'Amplitude');
            title(ax1, sprintf('Raster — Ch %d | Set: %s', ich, setLabel), 'Interpreter','none');
    
            for ai = 1:n_AMP
                y = psth_curves{ai};
                if isempty(y) || numel(y) ~= numel(ctrs), continue; end 
                plot(ax2, ctrs, y, 'Color', cmap(ai,:), 'LineWidth', 1.5);
                % plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.5);
            end
            xline(ax2, 0, 'r', 'LineWidth', 1.5);
            % xlim(ax2, ras_win);
            xlim(ax2, [-20,100]);
            if maxRate > 0
                ylim(ax2, [0 ceil(maxRate*1.1/10)*10]);
            else
                ylim(ax2, [0 1]);
            end
            xlabel(ax2, 'Time (ms)');
            ylabel(ax2, 'Rate (sp/s)');
            legend(ax2, arrayfun(@(a) sprintf('%.0f µA', a), Amps, 'UniformOutput', false), 'Box','off','Location','northeast');
        end
    end
 end


 if firing_rate_curve_per_amp_StimSet == 1
    fprintf('\n===== Firing Rate Plot for Channel %d =====\n', channels_to_plot);

    nSets  = size(uniqueComb,1);
    n_AMP  = numel(Amps);
    cmap   = lines(nSets);  % one color per stim set
    win_s  = diff(time_window_ms) / 1000;

    firing_rates = nan(nSets, n_AMP);        % mean firing rate
    firing_sems  = nan(nSets, n_AMP);        % SEM

    for s = 1:nSets
        trial_ids = find(combClass_win == s);

        for a = 1:n_AMP
            trials_this_amp = trial_ids(ampIdx(trial_ids) == a);
            if isempty(trials_this_amp), continue; end

            all_FR = [];  % collect firing rates from all channels

            for ch_i = 1:length(channels_to_plot)
                ch = d(channels_to_plot(ch_i));
                if isempty(sp_clipped{ch}), continue; end

                sp_times = sp_clipped{ch}(:,1);  % spike times in ms
                fr_per_trial = nan(numel(trials_this_amp), 1);

                for ti = 1:numel(trials_this_amp)
                    tr = trials_this_amp(ti);
                    t0_ms = trig(tr) / FS * 1000;
                    rel_sp_times = sp_times - t0_ms;
                    spk_mask = rel_sp_times >= time_window_ms(1) & rel_sp_times < time_window_ms(2);
                    fr_per_trial(ti) = sum(spk_mask) / win_s;
                end

                all_FR = [all_FR; fr_per_trial];  % concatenate across channels
            end

            firing_rates(s,a) = mean(all_FR);
            firing_sems(s,a)  = std(all_FR) / sqrt(numel(all_FR));  % SEM
        end
    end

    % Plotting 
    figure('Color','w'); hold on;
    for s = 1:nSets
        label = sprintf('Set %s', strjoin(arrayfun(@(x) sprintf('Ch%d', x), ...
            uniqueComb(s, uniqueComb(s,:) > 0), 'UniformOutput', false), '+'));

        errorbar(Amps, firing_rates(s,:), firing_sems(s,:), '-o', ...
            'LineWidth', 1.5, 'Color', cmap(s,:), 'DisplayName', label, ...
            'CapSize', 5, 'MarkerFaceColor', cmap(s,:));
    end

    xlabel('Amplitude (µA)');
    ylabel(sprintf('Firing Rate (sp/s)\n[%d–%d ms]', time_window_ms(1), time_window_ms(2)));
    title(sprintf('Mean Firing Rate — Channels %s', mat2str(channels_to_plot)));
    legend('Location','best', 'Interpreter','none');
    box off;

 end
