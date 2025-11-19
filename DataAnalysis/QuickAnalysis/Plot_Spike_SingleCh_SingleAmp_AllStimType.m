clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ===================== USER INPUTS ===================== %%
folder_single = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
folder_sim    = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Sim1';
folder_seq    = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq1';

target_channels = [19];      % channels to plot (recording channel index)
plot_amps       = [5];       % amplitudes (µA) → one figure per amp
Electrode_Type  = 1;         % 0: rigid, 1: single shank flex, 2: four shank flex

FS  = 30000;                 % Hz
pre_ms  = 5;                 % time BEFORE trigger (ms) for waveform binning
post_ms = 20;                % time AFTER trigger (ms)
bin_ms  = 5;                 % bin width for time windows (ms)

% Spike waveform amplitude threshold (to exclude crazy artefacts)
amp_threshold = 100;         % µV (abs value)

% Spike-count summary window (for printed per-condition stats)
spikeCount_window_ms = [0 20];   % e.g. [0 20] ms after trigger

%% ===================== DEFINE TIME BINS ===================== %%
edges = -pre_ms : bin_ms : post_ms;    % e.g. -5:5:20 → bins: [-5 0), [0 5), [5 10)...
nBins = numel(edges) - 1;
t_wave = (0:48) / FS * 1000;           % waveform sample times (ms), assuming 49 samples

%% LOAD DATA: SINGLE
cd(folder_single);
fprintf('\nLoading SINGLE dataset from:\n%s\n', folder_single);

% spikes
f_sp = dir('*sp_xia.mat');
assert(~isempty(f_sp), 'No *sp_xia.mat found in SINGLE folder.');
load(f_sp(1).name, 'sp_clipped');
D_single.sp = sp_clipped;

% triggers
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_single.trig = loadTrig(0);

% Stim params
S = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams;
D_single.simN   = S.simultaneous_stim;
D_single.E_MAP  = S.E_MAP;
D_single.nTrials = S.n_Trials;

% amplitudes
amps_all = cell2mat(Stim(2:end,16));
D_single.trialAmps = amps_all(1:D_single.simN:end);

% stim sets
E_NAME   = D_single.E_MAP(2:end);
stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, E_NAME);
stimChPerTrial = cell(D_single.nTrials,1);
for t = 1:D_single.nTrials
    rr = (t-1)*D_single.simN + (1:D_single.simN);
    v = unique(idx(rr)); v = v(v>0).';
    stimChPerTrial{t} = v;
end
comb = zeros(D_single.nTrials, D_single.simN);
for t = 1:D_single.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end
[D_single.uniqueComb,~,D_single.combClass] = unique(comb,'rows');

%%  LOAD DATA: SIMULTANEOUS
cd(folder_sim);
fprintf('\nLoading SIMULTANEOUS dataset from:\n%s\n', folder_sim);

f_sp = dir('*sp_xia.mat');
assert(~isempty(f_sp), 'No *sp_xia.mat found in SIM folder.');
load(f_sp(1).name, 'sp_clipped');
D_sim.sp = sp_clipped;

if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_sim.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams;
D_sim.simN    = S.simultaneous_stim;
D_sim.E_MAP   = S.E_MAP;
D_sim.nTrials = S.n_Trials;

amps_all = cell2mat(Stim(2:end,16));
D_sim.trialAmps = amps_all(1:D_sim.simN:end);

E_NAME   = D_sim.E_MAP(2:end);
stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, E_NAME);
stimChPerTrial = cell(D_sim.nTrials,1);
for t = 1:D_sim.nTrials
    rr = (t-1)*D_sim.simN + (1:D_sim.simN);
    v = unique(idx(rr)); v = v(v>0).';
    stimChPerTrial{t} = v;
end
comb = zeros(D_sim.nTrials, D_sim.simN);
for t = 1:D_sim.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end
[D_sim.uniqueComb,~,D_sim.combClass] = unique(comb,'rows');

%% LOAD DATA: SEQUENTIAL (First Pulse)
cd(folder_seq);
fprintf('\nLoading SEQUENTIAL (FirstPulse) dataset from:\n%s\n', folder_seq);

f_sp = dir('*sp_xia_FirstPulse.mat');
assert(~isempty(f_sp), 'No *sp_xia_FirstPulse.mat found in SEQ folder.');
load(f_sp(1).name, 'sp_seq');
D_seq.sp = sp_seq;    % first-pulse spikes only

if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
D_seq.trig = loadTrig(0);

S = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim = S.StimParams;
D_seq.simN    = S.simultaneous_stim;
D_seq.E_MAP   = S.E_MAP;
D_seq.nTrials = S.n_Trials;

amps_all = cell2mat(Stim(2:end,16));
D_seq.trialAmps = amps_all(1:D_seq.simN:end);

E_NAME   = D_seq.E_MAP(2:end);
stimNames = Stim(2:end,1);
[~, idx] = ismember(stimNames, E_NAME);
stimChPerTrial = cell(D_seq.nTrials,1);
for t = 1:D_seq.nTrials
    rr = (t-1)*D_seq.simN + (1:D_seq.simN);
    v = unique(idx(rr)); v = v(v>0).';
    stimChPerTrial{t} = v;
end
comb = zeros(D_seq.nTrials, D_seq.simN);
for t = 1:D_seq.nTrials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end
[D_seq.uniqueComb,~,D_seq.combClass] = unique(comb,'rows');

%% ===================== ELECTRODE MAP ===================== %%
d = Depth_s(Electrode_Type);   % maps recording channel index → sp{} index

%% ============================================================
%  MAIN LOOP: CHANNELS × AMPLITUDES
%  One figure per (channel, amplitude):
%    Row1: Single set 1
%    Row2: Single set 2 (if exists)
%    Row3: Sim set 1
%    Row4: Seq set 1
%  Columns: time bins (edges)
% ============================================================

stimRow_info = struct( ...
    'label',   {'Single Set 1','Single Set 2','Simultaneous','Sequential'}, ...
    'dataset', {'single','single','sim','seq'}, ...
    'set_idx', {1,2,1,1} );   % we use first 2 sets in SINGLE, first set in SIM/SEQ

for ch_rec = target_channels
    ch_idx = d(ch_rec);   % index into sp cell array

    for amp_val = plot_amps

        fprintf('\n==============================================\n');
        fprintf('Channel %d  |  Amplitude %.1f µA\n', ch_rec, amp_val);
        fprintf('==============================================\n');

        figure('Color','w','Position',[100 100 1600 900]);
        tl = tiledlayout(4, nBins, 'TileSpacing','compact','Padding','compact');
        title(tl, sprintf('Spike Waveforms — Ch %d — %.1f µA', ch_rec, amp_val), ...
              'FontSize',16);

        % ---- Loop over the 4 rows (stim types) ----
        for row = 1:4

            ds_name = stimRow_info(row).dataset;
            set_id  = stimRow_info(row).set_idx;
            rowLabel = stimRow_info(row).label;

            % Select dataset struct
            switch ds_name
                case 'single', D = D_single;
                case 'sim',    D = D_sim;
                case 'seq',    D = D_seq;
            end

            % Check that this stimulation set exists
            if set_id > size(D.uniqueComb,1)
                % No such set → leave row blank but label it
                if nBins >= 1
                    ax = nexttile( (row-1)*nBins + 1 );
                    text(ax,0.5,0.5,sprintf('%s (No set %d)',rowLabel,set_id), ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                    axis(ax,'off');
                end
                continue;
            end

            % Stim channels for this set
            stimCh = D.uniqueComb(set_id, D.uniqueComb(set_id,:)>0);
            stimCh_str = strjoin(arrayfun(@(x)sprintf('Ch%d',x), stimCh,'UniformOutput',false), ', ');

            % Spike data for this dataset and channel
            S_ch = D.sp{ch_idx};
            if isempty(S_ch)
                % No spikes at all; just annotate row
                if nBins >= 1
                    ax = nexttile( (row-1)*nBins + 1 );
                    text(ax,0.5,0.5,sprintf('%s — %s\n(no spikes)',rowLabel,stimCh_str), ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                    axis(ax,'off');
                end
                continue;
            end

            sp_times = S_ch(:,1);        % ms
            sp_wave  = S_ch(:,2:end);    % waveform

            % Basic amplitude sanity filter
            valid_idx = all(abs(sp_wave) <= amp_threshold, 2);
            sp_times  = sp_times(valid_idx);
            sp_wave   = sp_wave(valid_idx,:);
            if isempty(sp_times)
                continue;
            end

            % Trials for this set + amplitude
            trial_mask = (D.combClass == set_id) & (abs(D.trialAmps - amp_val) < 1e-6);
            trial_ids  = find(trial_mask);
            nTrials    = numel(trial_ids);

            if nTrials == 0
                if nBins >= 1
                    ax = nexttile( (row-1)*nBins + 1 );
                    text(ax,0.5,0.5,sprintf('%s — %s\n(no trials at %.1f µA)', ...
                        rowLabel,stimCh_str,amp_val), ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                    axis(ax,'off');
                end
                continue;
            end

            % -------- Collect waveforms per time bin --------
            all_spikes_by_bin = cell(nBins, 1);

            for tr = trial_ids'
                t0_ms = D.trig(tr)/FS*1000;

                % spikes around this trigger
                mask = sp_times >= (t0_ms - pre_ms) & sp_times <= (t0_ms + post_ms);
                if ~any(mask), continue; end

                rel_times = sp_times(mask) - t0_ms;  % ms, can be <0
                waves     = sp_wave(mask,:);

                bin_ids = discretize(rel_times, edges);
                for j = 1:numel(rel_times)
                    b = bin_ids(j);
                    if ~isnan(b)
                        all_spikes_by_bin{b}(end+1,:) = waves(j,:);
                    end
                end
            end

            % Determine y-limit from all waveforms in this row
            all_waves_row = cell2mat(all_spikes_by_bin(:));
            if isempty(all_waves_row)
                if nBins >= 1
                    ax = nexttile( (row-1)*nBins + 1 );
                    text(ax,0.5,0.5,sprintf('%s — %s\n(no spikes in window)', ...
                        rowLabel,stimCh_str), ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                    axis(ax,'off');
                end
                continue;
            end
            y_lim = [-1 1] * ceil(max(abs(all_waves_row(:))) / 50) * 50;

            % -------- Per-condition spike counts (user-defined window) --------
            totalSpikes_cond = 0;
            perTrial_counts  = zeros(nTrials,1);
            for i_tr = 1:nTrials
                tr = trial_ids(i_tr);
                t0_ms = D.trig(tr)/FS*1000;
                mask = sp_times >= (t0_ms + spikeCount_window_ms(1)) & ...
                       sp_times <= (t0_ms + spikeCount_window_ms(2));
                n_sp = sum(mask);
                perTrial_counts(i_tr) = n_sp;
                totalSpikes_cond = totalSpikes_cond + n_sp;
            end
            mean_perTrial_cond = totalSpikes_cond / nTrials;

            fprintf('%s — StimCh [%s] | Ch %d | %.1f µA | window [%g %g] ms: ', ...
                rowLabel, stimCh_str, ch_rec, amp_val, ...
                spikeCount_window_ms(1), spikeCount_window_ms(2));
            fprintf('total = %d, mean = %.2f spikes/trial\n', ...
                totalSpikes_cond, mean_perTrial_cond);

            % -------- Plot per-bin waveforms in this row (4 × nBins grid) --------
            for b = 1:nBins
                tile_idx = (row-1)*nBins + b;
                ax = nexttile(tile_idx);
                hold(ax,'on'); box(ax,'off');

                bin_waves = all_spikes_by_bin{b};
                total_spikes_bin = size(bin_waves,1);

                if total_spikes_bin > 0
                    % Align each waveform by its negative peak
                    aligned_waves = zeros(size(bin_waves));
                    for k = 1:total_spikes_bin
                        [~, min_idx] = min(bin_waves(k,:));
                        shift = ceil(size(bin_waves,2)/2) - min_idx;
                        aligned_waves(k,:) = circshift(bin_waves(k,:), shift, 2);
                    end

                    plot(ax, t_wave, aligned_waves', 'Color', [0 0 0 0.25]);
                end

                % Bin time range
                t1 = edges(b);
                t2 = edges(b+1);
                mean_perTrial_bin = total_spikes_bin / nTrials;

                title(ax, sprintf('%d–%d ms (total=%d, mean=%.2f/trial)', ...
                    t1, t2, total_spikes_bin, mean_perTrial_bin), ...
                    'FontSize',9);
                ylim(ax, y_lim);
                xlim(ax, [min(t_wave) max(t_wave)]);
                axis square;
                box off;
                if b == 1
                    % Bold stimulation type + channels (left-side labels)
                    ylabel(ax, { ...
                        sprintf('\\bf%s\\rm', rowLabel), ...     % e.g. "Single Set 1"
                        sprintf('Stim Ch: %s', stimCh_str) ...   % e.g. "Stim Ch: Ch5"
                        }, 'FontSize', 12, 'Interpreter','tex');
                end
                if row == 4
                    xlabel(ax, 'Waveform time (ms)');
                end
            end

        end % row loop

    end % amp loop
end % channel loop

fprintf('\nDone.\n');