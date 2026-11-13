clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
spike_chn_start = 1;
spike_chn_end   = 32;   % nChn (Depth_s index)
Electrode_Type  = 1;    % 0: rigid; 1: single-shank flex; 2: four-shank flex

data_folder = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Exp1_Seq4';

FS = 30000;            % Sampling frequency
win_ms       = 100;    % total window after trigger (ms)
bin_ms       = 2;      % bin size (ms)
nBins        = 40/bin_ms;   % we only use first 100 ms
amp_threshold = 300;   % optional plotting limit (set Inf to disable)

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

%% ================= LOAD QC / FILTERED SPIKES =================
sp_use = [];
fprintf('Trying to load QC / filtered spikes...\n');

qc_file   = [base_name '.sp_xia_QC.mat'];
ssd_file  = [base_name '.sp_xia_SSD.mat'];
base_file = [base_name '.sp_xia.mat'];

if isfile(qc_file)
    fprintf('Found QC file: %s\n', qc_file);
    S = load(qc_file);
    if isfield(S,'sp_qc')
        sp_use = S.sp_qc;
    elseif isfield(S,'sp_corr')
        sp_use = S.sp_corr;
    elseif isfield(S,'sp_clipped')
        sp_use = S.sp_clipped;
    else
        error('No usable spike variable (sp_qc/sp_corr/sp_clipped) in %s', qc_file);
    end

elseif isfile(ssd_file)
    fprintf('Found SSD file: %s\n', ssd_file);
    S = load(ssd_file);
    if isfield(S,'sp_pca')
        sp_use = S.sp_pca;
    elseif isfield(S,'sp_corr')
        sp_use = S.sp_corr;
    elseif isfield(S,'sp_SSD')
        sp_use = S.sp_SSD;
    elseif isfield(S,'sp_in')
        sp_use = S.sp_in;
    else
        error('No usable spike variable (sp_corr/sp_SSD/sp_in) in %s', ssd_file);
    end

elseif isfile(base_file)
    fprintf('Falling back to base spike file: %s\n', base_file);
    S = load(base_file);
    if isfield(S,'sp_clipped')
        sp_use = S.sp_clipped;
    elseif isfield(S,'sp')
        sp_use = S.sp;
    else
        error('No usable spike variable (sp_clipped/sp) in %s', base_file);
    end
else
    error('No spike file found: %s, %s, or %s', qc_file, ssd_file, base_file);
end

nCh = numel(sp_use);

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

%% ---- Stimulation sets (order-sensitive) ----
E_NAME    = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_NAME);

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    % KEEP ORDER — NOT unique()
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

%% ---- Pulse Train Period (for completeness; not directly used in plotting) ----
pulseTrain_all = cell2mat(StimParams(2:end,9));    % col 9: pulse train period
pulseTrain     = pulseTrain_all(1:simultaneous_stim:end);
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods); 

%% ================= ELECTRODE MAP =================
d = Depth_s(Electrode_Type);   % returns Intan channel index for each Depth_s index

%% ================= SPIKE WAVEFORM PLOTTING =================

% We will set t_wave based on the actual waveform length in the first non-empty channel
example_ch = find(~cellfun(@isempty, sp_use), 1, 'first');
if isempty(example_ch)
    error('All channels are empty in sp_use. Nothing to plot.');
end
wf_len = size(sp_use{example_ch}, 2) - 1;  % subtract timestamp column
t_wave = (0:wf_len-1) / FS * 1000;        % ms

for ich = spike_chn_start:spike_chn_end
    ch = d(ich);   % map depth index to Intan channel index
    if ch > nCh || isempty(sp_use{ch}), continue; end

    sp_times = sp_use{ch}(:,1);
    sp_wave  = sp_use{ch}(:,2:end);

    % Optional amplitude limit (for plotting only)
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

        figure('Name', sprintf('Ch %d | StimSet %d (%s)| Sequential Pulse', ich, set_id, stimLabel), ...
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
            legend_handles(a) = plot(nan,nan,'-', 'Color', cmap(a,:), 'LineWidth',1.5);
            legend_labels{a}  = sprintf('%g µA', Amps(a));
        end
        legend(legend_handles, legend_labels, 'Location','northeastoutside');
    end
end