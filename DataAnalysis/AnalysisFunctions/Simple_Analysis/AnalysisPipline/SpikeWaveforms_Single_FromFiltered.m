clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
spike_chn_start = 33;
spike_chn_end   = 64;   % nChn
Electrode_Type  = 2;    % 0: rigid; 1: single-shank flex; 2: four-shank flex

data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Single2_251128_125722';

FS = 30000;            % Sampling frequency
win_ms       = 300;    % total window after trigger (ms)
bin_ms       = 5;      % bin size (ms)
nBins        = 100/bin_ms;   % we only use first 100 ms
amp_threshold = 100;   % optional plotting limit (set Inf to disable)

layout_row = 4;        % subplot layout
layout_col = 5;

%% ================= CHECK FOLDER =================
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% ================= BASE NAME (your preferred way) =================
parts       = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);  % e.g. 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback
end

%% ================= LOAD FILTERED SPIKES =================
fprintf('Loading filtered spikes from sp_xia.mat...\n');
filtered_file = [base_name '.sp_xia.mat'];
assert(isfile(filtered_file), ...
    'Filtered spike file %s not found. Run your spike filtering script first.', filtered_file);

load(filtered_file, 'sp_clipped');   % expects cell array {ch} with [time, waveform...]

%% ================= LOAD TRIGGERS =================
if isempty(dir('*.trig.dat'))
    cur_dir = pwd; 
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);

%% ================= LOAD StimParams & DECODE =================
fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');
S = load(fullfile(data_folder, fileDIR(1).name), ...
         'StimParams', 'simultaneous_stim', 'CHN', 'E_MAP', 'n_Trials');

StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
CHN                = S.CHN;
E_MAP              = S.E_MAP;
n_Trials           = S.n_Trials;

%% ---- Amplitudes ----
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:simultaneous_stim:end);
[Amps, ~, ampIdx] = unique(trialAmps(:));
Amps(Amps == -1) = 0;
n_AMP = numel(Amps);
cmap  = lines(n_AMP);

%% ---- Stimulation sets (order-insensitive, as in your original) ----
E_NAME    = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_NAME);

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v  = unique(idx_all(rr)); 
    v  = v(v > 0).';
    stimChPerTrial_all{t} = v;
end

comb = zeros(n_Trials, simultaneous_stim);
for i = 1:n_Trials
    v = stimChPerTrial_all{i};
    comb(i,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb, 'rows');
combClass_win = combClass;
nSets         = size(uniqueComb,1);

%% ---- Pulse Train Period (for completeness; not directly used in plotting) ----
pulseTrain_all = cell2mat(StimParams(2:end,9));    % col 9: pulse train period
pulseTrain     = pulseTrain_all(1:simultaneous_stim:end);
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:)); %#ok<ASGLU>
n_PULSE = numel(PulsePeriods); %#ok<NASGU>

%% ================= ELECTRODE MAP =================
d = Depth_s(Electrode_Type);   % returns Intan channel index for each Depth_s index

%% ================= SPIKE WAVEFORM PLOTTING =================
t_wave = (0:48) / FS * 1000;   % assumes 49-sample waveforms, adjust if needed

for ich = spike_chn_start:spike_chn_end
    ch = d(ich);   % map to Intan channel index
    if ch > numel(sp_clipped) || isempty(sp_clipped{ch}), continue; end

    sp_times = sp_clipped{ch}(:,1);
    sp_wave  = sp_clipped{ch}(:,2:end);

    % Optional plotting-only amplitude limit
    valid_idx = all(abs(sp_wave) <= amp_threshold, 2);
    sp_times  = sp_times(valid_idx);
    sp_wave   = sp_wave(valid_idx,:);
    if isempty(sp_times), continue; end

    for s = 1:nSets
        set_id    = s;
        trial_mask = (combClass_win == set_id);
        if ~any(trial_mask), continue; end

        trial_ids = find(trial_mask);
        all_spikes_by_bin_amp = cell(nBins, n_AMP);

        % -------- collect waveforms per bin × amplitude --------
        for idx = 1:length(trial_ids)
            tr    = trial_ids(idx);
            t0_ms = trig(tr) / FS * 1000;
            amp_id = ampIdx(tr);

            mask_sp = sp_times >= t0_ms & sp_times < (t0_ms + win_ms);
            if ~any(mask_sp), continue; end

            rel_times = sp_times(mask_sp) - t0_ms;
            waveforms = sp_wave(mask_sp,:);

            for j = 1:length(rel_times)
                bin_idx = floor(rel_times(j) / bin_ms) + 1;
                if bin_idx >= 1 && bin_idx <= nBins
                    all_spikes_by_bin_amp{bin_idx, amp_id}(end+1,:) = waveforms(j,:);
                end
            end
        end

        % Determine y-limits
        all_waves = cell2mat(all_spikes_by_bin_amp(:));
        if isempty(all_waves), continue; end
        y_max = max(abs(all_waves(:)));
        y_lim = [-1 1] * ceil(y_max/50)*50;

        % Figure title
        stimIdx   = uniqueComb(set_id,:); 
        stimIdx   = stimIdx(stimIdx > 0);
        stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimIdx, 'UniformOutput', false), ', ');

        figure('Name', sprintf('Ch %d | StimSet %d (%s)', ich, set_id, stimLabel), ...
               'Color','w','Position', [100 100 1400 800]);        
        tiledlayout(layout_row, layout_col, 'Padding','compact', 'TileSpacing','compact');

        % -------- per-bin waveform panels --------
        for b = 1:nBins
            nexttile; hold on;

            for a = 1:n_AMP
                waves = all_spikes_by_bin_amp{b,a};
                if isempty(waves), continue; end

                % Align each waveform to its minimum
                aligned_waves = zeros(size(waves));
                for k = 1:size(waves,1)
                    [~, min_idx] = min(waves(k,:));
                    shift = ceil(size(waves,2)/2) - min_idx;
                    aligned_waves(k,:) = circshift(waves(k,:), shift, 2);
                end

                plot(t_wave, aligned_waves', 'Color', [cmap(a,:) 0.3]);
            end

            % spike count in this bin (all amplitudes)
            spike_count = 0;
            for a = 1:n_AMP
                spike_count = spike_count + size(all_spikes_by_bin_amp{b,a}, 1);
            end

            title(sprintf('%d–%d ms (%d spikes)', ...
                  (b-1)*bin_ms, b*bin_ms, spike_count));
            xlabel('Time (ms)');
            ylabel('Voltage (µV)');
            ylim(y_lim);
            yticks(linspace(y_lim(1), y_lim(2), 3));
            xticks(round(linspace(t_wave(1), t_wave(end), 3)));
            axis square; grid on;
        end

        % -------- global legend for amplitudes --------
        legend_handles = gobjects(n_AMP,1);
        legend_labels  = cell(n_AMP,1);
        for a = 1:n_AMP
            legend_handles(a) = plot(nan,nan,'-','Color',cmap(a,:), 'LineWidth',1.5);
            legend_labels{a}  = sprintf('%g µA', Amps(a));
        end
        legend(legend_handles, legend_labels, 'Location','northeastoutside');
    end
end