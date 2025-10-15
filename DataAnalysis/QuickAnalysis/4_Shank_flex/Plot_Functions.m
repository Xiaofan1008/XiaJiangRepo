%% Plot Choice
spike_plot_per_amp_StimSet = 1;
raster_psth_per_amp_StimSet = 1;





%% Pre Set
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
FS=30000;
% Load .sp.mat file
sp_files = dir('*.sp.mat');
assert(~isempty(sp_files), 'No .sp.mat file found in the current folder.');
sp_filename = sp_files(1).name;
fprintf('Loading spike file: %s\n', sp_filename);
S = load(sp_filename);
if isfield(S, 'sp')
    sp = S.sp;
else
    error('Variable "sp" not found in %s.', sp_filename);
end
% Load Trigger
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);

%% Parameters
win_ms = 300;         % total time window after each trigger (ms)
bin_ms = 10;          % bin size (ms)
nBins  = win_ms / bin_ms;
nChn   = numel(sp);
amp_threshold = 150;  % max allowed amplitude (µV)

%% Load StimParams and decode amplitudes & stimulation sets
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
S = load(fileDIR(1).name,'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');
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
d = Depth_s(2); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

if spike_plot_per_amp_StimSet == 1
    %% Spike plotting loop
    for ich = 1:nChn
        ch = d(ich); % channel map to Intan
        if isempty(sp{ch}), continue; end
        fprintf('Processing Channel %d...\n', ich);
    
        % Spike data
        sp_times = sp{ch}(:,1);
        sp_wave  = sp{ch}(:,2:end);
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
                   'Color','w','Position', [100 100 1400 800]);
            tiledlayout(6, 5, 'Padding','compact', 'TileSpacing','compact');
    
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
    %% Parameters
    ras_win   = [-20 300];   % ms
    bin_ms    = 1;           % bin size
    smooth_ms = 3;           % smoothing window
    
    % Time axis
    edges = ras_win(1):bin_ms:ras_win(2);
    ctrs  = edges(1:end-1) + diff(edges)/2;
    bin_s = bin_ms / 1000;
    
    % Smoothing kernel
    g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
    g = g / sum(g);
    
    %% Main Loop
    for si = 1:nSets
        set_id = si;
        stimIdx = uniqueComb(set_id, :); stimIdx = stimIdx(stimIdx > 0);
        setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx, 'UniformOutput', false), ' + ');
    
        for ch = 1:nChn
            ich = ch;
            ch = d(ich);
            if isempty(sp{ch}), continue; end
    
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
                color = cmap(ai,:);
                amp_trials = find(ampIdx == ai & combClass_win == set_id);
                nTr = numel(amp_trials);
                spkTimesPerTrial = cell(nTr,1);
                counts = zeros(1, numel(edges)-1);
    
                for t = 1:nTr
                    tr = amp_trials(t);
                    S_ch = sp{ch};
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
                ytick_labels{end+1} = sprintf('%.0f \x03bcA', amp_val);
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
            xlim(ax1, ras_win);
            ylim(ax1, [0, y_cursor]);
            yticks(ax1, ytick_vals);
            yticklabels(ax1, ytick_labels);
            ylabel(ax1, 'Amplitude');
            title(ax1, sprintf('Raster — Ch %d | Set: %s', ich, setLabel), 'Interpreter','none');
    
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
            legend(ax2, arrayfun(@(a) sprintf('%.0f µA', a), Amps, 'UniformOutput', false), 'Box','off','Location','northeast');
        end
    end
 end
