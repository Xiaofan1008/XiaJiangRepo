clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));


%% Spike Filtering Parameters
Spike_filtering = 0; 
pos_limit = 100;    % upper bound (µV)
neg_limit = -100;  % lower bound (µV)
width_min_ms = 0.01;
width_max_ms = 0.4;

%% Spike waveform channel 
spike_chn_start = 22;
spike_chn_end = 22; %nChn

selectedAmps = [5]; 
Electrode_Type = 1; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

%% Choose Folder

data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim1'; 
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
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

%% Pre Set
FS=30000; % Sampling frequency
sp_files = dir(fullfile(data_folder, '*.sp.mat'));
assert(~isempty(sp_files), 'No .sp.mat file found in the current folder.');
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

%% Spike Amplitude Filtering (Before Plotting)
% % sp_clipped = sp;   % copy original spike structure

% if Spike_filtering == 1
%     fprintf('\n===== Spike Amplitude Filtering =====\n'); 
% 
%     t_axis = (0:48) / FS * 1000;
%     for ch = 1:numel(sp)
%         if isempty(sp{ch}), continue; end
% 
%         waveforms = sp{ch}(:, 2:end);
%         nSpikes = size(waveforms, 1);
% 
%         % Amplitude-based filtering
%         max_vals = max(waveforms, [], 2);
%         min_vals = min(waveforms, [], 2);
%         in_range = (max_vals <= pos_limit) & (min_vals >= neg_limit);
% 
%         % Zero-crossing check
%         has_positive = any(waveforms > 0, 2);
%         has_negative = any(waveforms < 0, 2);
%         crosses_zero = has_positive & has_negative;
% 
%         % Spike width calculation (min to next max)
%         [~, min_idx] = min(waveforms, [], 2);  % index of neg peak
%         spike_widths = nan(nSpikes,1);
% 
%         for i = 1:nSpikes
%             mn = min_idx(i);
%             if mn >= size(waveforms, 2), continue; end
%             [~, rel_max_idx] = max(waveforms(i, mn:end));
%             spike_widths(i) = t_axis(mn + rel_max_idx - 1) - t_axis(mn);
%         end
%         width_valid = spike_widths >= width_min_ms & spike_widths <= width_max_ms;
% 
%         % Keep spikes within amplitude bounds
%         % valid_idx = (max_vals <= pos_limit) & (min_vals >= neg_limit);
%         % Final combined condition
%         valid_idx = in_range & crosses_zero;
%         % valid_idx = in_range & crosses_zero & width_valid;
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

%% Load StimParams and decode amplitudes, stimulation sets, ISI
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
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
d = Depth_s(Electrode_Type); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex


%% Spike Waveform Parameters
% win_ms = 300;         % total time window after each trigger (ms)
% bin_ms = 5;          % bin size (ms)

pre_ms  = 5;     % time BEFORE trigger to include (ms)
post_ms = 20;    % time AFTER trigger to include (ms)

bin_ms = 5;       % bin size (unchanged)
% nBins  = win_ms / bin_ms;
win_ms   = pre_ms + post_ms; % total window width

% New bin edges spanning BEFORE & AFTER trigger
edges = -pre_ms : bin_ms : post_ms;
nBins = numel(edges) - 1;
amp_threshold = 100;  % max allowed amplitude (µV)

layout_row = 1; % plot time window layout in rows
layout_col = 5; % plot time window layout in columns


%% Spike plotting loop
for ich = spike_chn_start:spike_chn_end
    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end

    % Spike data
    sp_times = sp_clipped{ch}(:,1);
    sp_wave  = sp_clipped{ch}(:,2:end);
    valid_idx = all(abs(sp_wave) <= amp_threshold, 2);
    sp_times = sp_times(valid_idx);
    sp_wave  = sp_wave(valid_idx,:);
    if isempty(sp_times), continue; end
    t_wave   = (0:size(sp_wave,2)-1) / FS * 1000;

    for s = 1:nSets
        trial_mask = (combClass_win == s);
        trial_ids = find(trial_mask);
        if isempty(trial_ids), continue; end

        %% Select amplitudes
        ampValList = Amps(:)';
        ampIndexMap = containers.Map(ampValList, 1:n_AMP);

        valid_amp_idx = [];
        for a = selectedAmps
            if isKey(ampIndexMap, a)
                valid_amp_idx(end+1) = ampIndexMap(a);
            end
        end
        if isempty(valid_amp_idx), continue; end

        %% Allocate bins (REFINED)
        all_spikes_by_bin_amp = cell(nBins, numel(valid_amp_idx));

        %% --- Loop through trials and collect spikes ---
        for idx = 1:length(trial_ids)
            tr = trial_ids(idx);
            amp_id_all = ampIdx(tr);
            if ~ismember(amp_id_all, valid_amp_idx)
                continue
            end

            amp_id_local = find(valid_amp_idx == amp_id_all);

            t0_ms = trig(tr) / FS * 1000;

            % New window: include pre-trigger + post-trigger
            spk_mask = sp_times >= (t0_ms - pre_ms) & sp_times <= (t0_ms + post_ms);
            if ~any(spk_mask), continue; end

            rel_times = sp_times(spk_mask) - t0_ms;    % now includes negative times
            waveforms = sp_wave(spk_mask,:);

            % Assign to bins (with negative times allowed)
            bin_ids = discretize(rel_times, edges);

            for j = 1:numel(rel_times)
                b = bin_ids(j);
                if ~isnan(b)
                    all_spikes_by_bin_amp{b, amp_id_local}(end+1,:) = waveforms(j,:);
                end
            end
        end

        %% Determine y-limits
        all_waves = cell2mat(all_spikes_by_bin_amp(:));
        if isempty(all_waves), continue; end
        y_lim = [-1 1] * ceil(max(abs(all_waves(:))) / 50) * 50;

        %% --- PLOT ---
        stimIdx = uniqueComb(s,:); stimIdx = stimIdx(stimIdx > 0);
        stimLabel = strjoin(arrayfun(@(x)sprintf('Ch%d',x), stimIdx,'UniformOutput',false), ', ');

        figure('Name', sprintf('Ch %d | StimSet %d (%s)', ich, s, stimLabel), ...
               'Color','w','Position',[100 100 1400 800]);
        tiledlayout(layout_row, layout_col, 'Padding','compact','TileSpacing','compact');

        for b = 1:nBins
            nexttile; hold on;

            total_spikes = 0;
            for a_local = 1:numel(valid_amp_idx)
                amp_global_id = valid_amp_idx(a_local);
                this_bin = all_spikes_by_bin_amp{b, a_local};
                if isempty(this_bin), continue; end

                % Align by minimum
                aligned_waves = zeros(size(this_bin));
                for k = 1:size(this_bin,1)
                    [~, min_idx] = min(this_bin(k,:));
                    shift = ceil(size(this_bin,2)/2) - min_idx;
                    aligned_waves(k,:) = circshift(this_bin(k,:), shift, 2);
                end

                plot(t_wave, aligned_waves', 'Color', [cmap(amp_global_id,:) 0.4]);
                total_spikes = total_spikes + size(this_bin,1);
            end

            title(sprintf('%d–%d ms (%d spikes)', edges(b), edges(b+1), total_spikes));
            xlabel('Time (ms)');
            ylabel('Voltage (µV)');
            ylim(y_lim);
            axis square;
            box off;
        end

        %% Legend
        legend_handles = arrayfun(@(id) plot(nan,nan,'-','Color',cmap(id,:),'LineWidth',2), valid_amp_idx);
        legend_labels  = arrayfun(@(id) sprintf('%g µA', Amps(id)), valid_amp_idx,'UniformOutput',false);
        legend(legend_handles, legend_labels, 'Location','northeastoutside');
    end
end