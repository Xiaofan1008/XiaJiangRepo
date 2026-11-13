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
spike_chn_start = 1;
spike_chn_end = 64; %nChn
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

%% Choose Folder

% data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim5'; 
data_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Single1_251203_121428';
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
win_ms = 300;         % total time window after each trigger (ms)
bin_ms = 5;          % bin size (ms)
% nBins  = win_ms / bin_ms;
nBins  = 100 / bin_ms;
nChn   = numel(sp);
amp_threshold = 100;  % max allowed amplitude (µV)

layout_row = 4; % plot time window layout in rows
layout_col = 5; % plot time window layout in columns

%% Spike plotting loop
for ich = spike_chn_start:spike_chn_end %1:nChn
    ch = d(ich); % channel map to Intan
    if isempty(sp_clipped{ch}), continue; end
    % fprintf('Processing Channel %d...\n', ich);

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

        % Determine Y-limit
        all_waves = cell2mat(all_spikes_by_bin_amp(:));
        if isempty(all_waves), continue; end
        y_max = max(abs(all_waves(:)));
        y_lim = [-1 1] * ceil(y_max/50)*50;

        % Plot
        stimIdx = uniqueComb(set_id,:); stimIdx = stimIdx(stimIdx > 0);
        stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx, 'UniformOutput', false), ', ');
        figure('Name', sprintf('Ch %d | StimSet %d (%s)', ich, set_id, stimLabel), ...
               'Color','w','Position', [100 100 1400 800]);        
        % tiledlayout(4, 5, 'Padding','compact', 'TileSpacing','compact');
        tiledlayout(layout_row, layout_col, 'Padding','compact', 'TileSpacing','compact');


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