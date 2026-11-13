%% ============================================================
%   Plot_SpikeWaveforms_FromFiltered.m
%   Loads ONLY filtered spikes (sp_xia.mat)
%   Plots spike waveforms by time-bin × amplitude × stim set
%   Xia, 2025
%% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
spike_chn_start = 28;
spike_chn_end   = 30;

plot_amps       = [6];     % amplitudes to include
Electrode_Type  = 2;

data_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq4_5ms_251125_154235';

win_ms      = 300;     % window after stimulus to search spikes
bin_ms      = 5;       % time bin for grouping waveforms
nBins       = floor(100/bin_ms);  
amp_threshold = 100;   % for plotting only (optional)
layout_row  = 4;
layout_col  = 5;

FS = 30000;

%% ================= CHECK FOLDER =================
if ~isfolder(data_folder), error('Invalid folder'); end
cd(data_folder);

%% ================= EXTRACT BASE NAME =================
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);  % 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback if no underscores
end

%% ================= LOAD FILTERED SPIKES =================
fprintf('Loading filtered spikes from sp_xia.mat...\n');
filtered_file = [base_name '.sp_xia.mat'];
assert(isfile(filtered_file), ...
    'Filtered spike file sp_xia.mat not found. Run SpikeFiltering_Save.m first.');
load(filtered_file, 'sp_clipped');

%% ================= LOAD TRIGGERS =================
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

%% ================= LOAD STIM PARAMS =================
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR),'No *_exp_datafile_*.mat found.');

S = load(fileDIR(1).name, ...
         'StimParams','simultaneous_stim','E_MAP','n_Trials');

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP             = S.E_MAP;
n_Trials          = S.n_Trials;

%% ================= DECODE AMPLITUDES =================
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:simultaneous_stim:end);

[Amps,~,ampIdx] = unique(trialAmps(:));
Amps(Amps==-1) = 0;
n_AMP = numel(Amps);
cmap  = lines(n_AMP);

ampMask = ismember(trialAmps, plot_amps);

%% ================= DECODE STIM SETS =================
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

stimChPerTrial_all = cell(n_Trials,1);
for t=1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = idx_all(rr);
    v = v(v>0);
    stimChPerTrial_all{t} = v(:)';
end

comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t,1:numel(v)) = v;
end

[uniqueComb,~,combClass_win] = unique(comb,'rows','stable');
nSets = size(uniqueComb,1);

%% ================= DECODE PTD =================
postTrigDelay_all = cell2mat(StimParams(2:end,6));
postTrigDelay     = postTrigDelay_all(2:simultaneous_stim:end);

[uniqueDelays,~,delayIdx] = unique(postTrigDelay);
n_DELAYS = numel(uniqueDelays);

%% ================= SETUP ELECTRODE MAP =================
d = Depth_s(Electrode_Type);

%% ============================================================
%                MAIN SPIKE WAVEFORM PLOTTING
%% ============================================================

for ich = spike_chn_start:spike_chn_end

    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end

    sp_times = sp_clipped{ch}(:,1);
    sp_wave  = sp_clipped{ch}(:,2:end);

    % optional amplitude limit for display
    ok = all(abs(sp_wave)<=amp_threshold,2);
    sp_times = sp_times(ok);
    sp_wave  = sp_wave(ok,:);
    if isempty(sp_times), continue; end

    t_wave = (0:size(sp_wave,2)-1)/FS*1000;

    % ---- Iterate PTDs ----
    for delay_i = 1:n_DELAYS
        delay_val = uniqueDelays(delay_i);

        % ---- Iterate stimulation sets ----
        for s = 1:nSets
            set_id = s;

            % sequence-sensitive set + PTD + amplitude choice
            trial_mask = (combClass_win == set_id) & ...
                         (delayIdx == delay_i) & ...
                         ampMask;

            if ~any(trial_mask), continue; end

            trial_ids = find(trial_mask);

            stimIdx = uniqueComb(set_id,:);
            stimIdx = stimIdx(stimIdx>0);

            figName = sprintf('Ch %d | Stim %s | PTD %d ms | Amps %s', ...
                ich, strjoin(string(stimIdx),','), delay_val/1000, mat2str(plot_amps));

            figure('Name',figName,'Color','w','Position',[100 100 1400 800]);
            tiledlayout(layout_row,layout_col,'Padding','compact','TileSpacing','compact');

            all_spikes_by_bin_amp = cell(nBins,n_AMP);

            % ---- Extract waveforms per trial ----
            for tr = trial_ids(:)'
                t0_ms = trig(tr)/FS*1000;
                amp_id = ampIdx(tr);

                mask_sp = sp_times>=t0_ms & sp_times<(t0_ms+win_ms);
                if ~any(mask_sp), continue; end

                rel_t = sp_times(mask_sp) - t0_ms;
                wavef = sp_wave(mask_sp,:);

                for k = 1:numel(rel_t)
                    b = floor(rel_t(k)/bin_ms)+1;
                    if b>=1 && b<=nBins
                        all_spikes_by_bin_amp{b,amp_id}(end+1,:) = wavef(k,:);
                    end
                end
            end

            % ---- Plot per bin ----
            allW = cell2mat(all_spikes_by_bin_amp(:));
            if isempty(allW), continue; end

            y_max = max(abs(allW(:)));
            y_lim = [-1 1]*ceil(y_max/50)*50;

            for b = 1:nBins
                nexttile; hold on

                for a = 1:n_AMP
                    if ~ismember(Amps(a), plot_amps), continue; end
                    waves = all_spikes_by_bin_amp{b,a};
                    if isempty(waves), continue; end

                    % align each waveform to its minimum
                    aligned = zeros(size(waves));
                    for k = 1:size(waves,1)
                        [~,mid] = min(waves(k,:));
                        shift = ceil(size(waves,2)/2)-mid;
                        aligned(k,:) = circshift(waves(k,:),shift,2);
                    end

                    plot(t_wave, aligned', 'Color', [cmap(a,:) 0.30]);
                end

                % count spikes in this bin
                cnt = 0;
                for a = 1:n_AMP
                    if ismember(Amps(a), plot_amps)
                        cnt = cnt + size(all_spikes_by_bin_amp{b,a},1);
                    end
                end

                title(sprintf('%d-%d ms (%d spikes)', ...
                      (b-1)*bin_ms, b*bin_ms, cnt));
                xlabel('Time (ms)');
                ylabel('µV');
                ylim(y_lim);
                axis square;
            end

            % ---- Legend ----
            keepA = find(ismember(Amps, plot_amps));
            h_leg = gobjects(numel(keepA),1);
            for ii = 1:numel(keepA)
                a = keepA(ii);
                h_leg(ii) = plot(nan,nan,'-','Color',cmap(a,:),'LineWidth',1.5);
            end

            legend(h_leg, arrayfun(@(x) sprintf('%g µA',x), ...
                   Amps(keepA), 'UniformOutput',false), ...
                   'Location','northeastoutside');

        end
    end
end