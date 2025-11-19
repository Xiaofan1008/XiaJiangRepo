%% ============================================================
% Raster + PSTH Plot for Specific Amplitudes (Sequential Stim)
% ============================================================
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder
data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Seq1';

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
    base_name = last_folder(1 : underscores(end-1) - 1);
else
    base_name = last_folder;
end

%% Choice
Spike_filtering = 0;
raster_chn_start = 30;
raster_chn_end   = 30; % nChn
Plot_Amps        = [5]; % µA amplitudes to plot
Electrode_Type = 1; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

%% Pre Set
FS = 30000; % Sampling frequency

sp_files = dir(fullfile(data_folder, '*.sp_xia_FirstPulse.mat'));
if isempty(sp_files)
    sp_files = dir(fullfile(data_folder, '*.sp_xia.mat'));
end
assert(~isempty(sp_files), 'No spike file found.');
load(sp_files(1).name, 'sp_seq', 'sp_clipped');
if exist('sp_seq','var')
    sp_clipped = sp_seq;
end
fprintf('Loaded spike data: %s\n', sp_files(1).name);

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);

%% Load StimParams and Decode
S = load(dir('*_exp_datafile_*.mat').name, ...
    'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP             = S.E_MAP;
n_Trials          = S.n_Trials;

% --- Amplitudes ---
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
[Amps, ~, ampIdx] = unique(trialAmps(:));
if any(Amps == -1), Amps(Amps == -1) = 0; end
n_AMP = numel(Amps);
cmap = lines(n_AMP);

% --- Stim Sets ---
E_NAME = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_NAME);
stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = unique(idx_all(rr)); v = v(v>0).';
    stimChPerTrial_all{t} = v;
end
comb = zeros(n_Trials, simultaneous_stim);
for i = 1:n_Trials
    v = stimChPerTrial_all{i};
    comb(i,1:numel(v)) = v;
end
[uniqueComb, ~, combClass] = unique(comb,'rows');
combClass_win = combClass;
nSets = size(uniqueComb,1);

% --- Pulse Periods ---
pulseTrain_all = cell2mat(StimParams(2:end,9)); 
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

% --- Electrode Map ---
d = Depth_s(Electrode_Type);

%% Raster Plot Parameters
ras_win       = [-20 100]; % ms
bin_ms_raster = 2;
smooth_ms     = 3;

edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;
g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
g = g / sum(g);

%% Plot Raster + PSTH
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

            figure('Color','w','Name',sprintf('Ch %d | Set %s (Seq)', ich, setLabel));
            tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
            ax1 = nexttile([3 1]); hold(ax1,'on'); box(ax1,'off');
            ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');

            maxRate = 0; 
            y_cursor = 0;
            ytick_vals = [];
            ytick_labels = {};
            total_spikes = 0;

            % === Loop through specific amplitudes only ===
            for ai = 1:length(Plot_Amps)
                amp_val = Plot_Amps(ai);
                amp_idx_match = find(abs(Amps - amp_val) < 1e-6);
                if isempty(amp_idx_match)
                    fprintf('⚠️ %.1f µA not found — skipping.\n', amp_val);
                    continue;
                end

                color = cmap(amp_idx_match,:);
                amp_trials = find(ampIdx == amp_idx_match & pulseIdx == pi & combClass_win == set_id);
                nTr = numel(amp_trials);
                if nTr == 0, continue; end

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

                % --- Raster ---
                for t = 1:nTr
                    tt = spkTimesPerTrial{t};
                    y0 = y_cursor + t;
                    for k = 1:numel(tt)
                        plot(ax1, [tt(k) tt(k)], [y0-0.4 y0+0.4], ...
                             'Color', color, 'LineWidth', 1.2);
                    end
                end

                ytick_vals(end+1) = y_cursor + nTr/2;
                ytick_labels{end+1} = sprintf('%.0f µA', amp_val);
                if nTr > 0
                    plot(ax1, ras_win, [y_cursor+nTr y_cursor+nTr], ...
                        'Color', [0.85 0.85 0.85], 'LineStyle','--');
                end
                y_cursor = y_cursor + nTr;

                % --- PSTH ---
                rate = counts / (nTr * bin_s);
                rate_s = filter(g, 1, rate);
                maxRate = max(maxRate, max(rate_s));
                plot(ax2, ctrs, rate_s, 'Color', color, 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('%.0f µA (n=%d)', amp_val, nTr));
            end

            % === Finalize Plots ===
            % (red 0 ms line ignored in legend)
            xline(ax1, 0, 'r', 'LineWidth', 1.5, 'HandleVisibility','off');
            xlim(ax1, [ras_win(1), ras_win(2)+10]);
            ylim(ax1, [0, y_cursor]);
            yticks(ax1, ytick_vals);
            yticklabels(ax1, ytick_labels);
            ylabel(ax1, 'Amplitude');
            title(ax1, sprintf('Raster — Ch %d | Set: %s (Seq)', ich, setLabel), 'Interpreter','none');

            xline(ax2, 0, 'r', 'LineWidth', 1.5, 'HandleVisibility','off');
            xlim(ax2, ras_win);
            ylim(ax2, [0 ceil(maxRate*1.1/10)*10]);
            xlabel(ax2, 'Time (ms)');
            ylabel(ax2, 'Rate (sp/s)');
            legend(ax2, 'show', 'Box','off','Location','northeast');

            fprintf('Ch %2d | Set %s — Total spikes: %d\n', ich, setLabel, total_spikes);
        end
    end
end