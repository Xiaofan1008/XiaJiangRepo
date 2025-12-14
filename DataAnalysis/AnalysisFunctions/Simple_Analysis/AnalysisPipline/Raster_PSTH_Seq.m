%% ============================================================
%   SCRIPT: RasterPlot filtered data
%   Loads pre-filtered spikes from sp_xia.mat
%   Plots RASTER + PSTH
%   Xia, 2025
%% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq4_5ms_251125_154235';
Electrode_Type = 1;  % 0 rigid, 1 flex, 2 4-shank flex
FS = 30000;

raster_chn_start = 1;
raster_chn_end   = 32;

ras_win       = [-20 100];
bin_ms        = 1;
smooth_ms     = 3;

%% ================= LOAD FILTERED SPIKES =================
cd(data_folder);
fprintf('Directory: %s\n', data_folder);

% Extract file name
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) >= 4
    base_name = last_folder(1 : underscores(end-1) - 1);  % 'Xia_Exp1_Seq'
else
    base_name = last_folder;  % fallback if no underscores
end

filtered_file = [base_name '.sp_xia.mat'];
assert(isfile(filtered_file), ...
    'Filtered spike file sp_xia.mat not found. Run SpikeFiltering_Save.m first.');

load(filtered_file, 'sp_clipped');
fprintf('Loaded filtered spikes from sp_xia.mat\n');

%% ================= LOAD TRIGGERS =================
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig = loadTrig(0);

%% ================= LOAD StimParams =================
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found.');

SS = load(fileDIR(1).name, ...
          'StimParams','simultaneous_stim','CHN','E_MAP','n_Trials');

StimParams        = SS.StimParams;
simultaneous_stim = SS.simultaneous_stim;
CHN               = SS.CHN;
E_MAP             = SS.E_MAP;
n_Trials          = SS.n_Trials;

%% ============================================================
%    DECODING BLOCK (AMPLITUDES, SETS, PULSE PERIOD, PTD)
%% ============================================================

% ---- Amplitude decoding ----
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);

[Amps, ~, ampIdx] = unique(trialAmps(:));
Amps(Amps == -1) = 0;
n_AMP = numel(Amps);
cmap = lines(n_AMP);

% ---- Decode stimulation sets (order preserved) ----
E_NAME = E_MAP(2:end);
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_NAME);

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = idx_all(rr);
    v = v(v > 0); 
    stimChPerTrial_all{t} = v(:)';
end

comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t,1:numel(v)) = v;
end

[uniqueComb,~,combClass_win] = unique(comb,'rows','stable');
nSets = size(uniqueComb,1);

% ---- Pulse Train Period (µs) ----
pulseTrain_all = cell2mat(StimParams(2:end,9));
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);

[PulsePeriods,~,pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

% ---- POST-TRIGGER DELAY (PTD, µs) ----
if simultaneous_stim > 1
    ptd_all = cell2mat(StimParams(3:simultaneous_stim:end,6));
else
    ptd_all = zeros(n_Trials,1);
end

PTD_us = ptd_all(:);
[PTD_values,~,ptdIdx] = unique(PTD_us);

%% ============================================================
%   Raster Plot Pre-Setup
%% ============================================================

edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

d = Depth_s(Electrode_Type);

%% ============================================================
%                   RASTER + PSTH PLOTTING
%% ============================================================

for ich = raster_chn_start:raster_chn_end
    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end

    for si = 1:nSets
        stimVec = uniqueComb(si,:);
        stimVec = stimVec(stimVec > 0);
        setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), ...
                           stimVec, 'UniformOutput', false), '→');

        for pi = 1:n_PULSE
            pulse_val = PulsePeriods(pi);

            for i_PTD = 1:numel(PTD_values)
                ptd_val = PTD_values(i_PTD);

                trials_sel = find( combClass_win == si & ...
                                   pulseIdx == pi & ...
                                   PTD_us == ptd_val );

                if isempty(trials_sel), continue; end

                figName = sprintf('Ch %d | Set %s | Pulse %d µs | PTD %d µs', ...
                                   ich, setLabel, pulse_val, ptd_val);
                figure('Color','w','Name',figName);
                tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                ax1 = nexttile([3 1]); hold(ax1,'on'); box(ax1,'off');
                title(ax1, sprintf('Raster — Ch %d | Set %s | PTD %d µs', ...
                                   ich, setLabel, ptd_val));

                ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');

                psth_curves = cell(1, n_AMP);
                y_cursor = 0;
                ytick_vals = [];
                ytick_labels = {};
                maxRate = 0;

                %% ---- Loop over amplitudes ----
                for ai = 1:n_AMP
                    amp_val = Amps(ai);
                    color   = cmap(ai,:);

                    amp_trials = trials_sel(ampIdx(trials_sel) == ai);
                    if isempty(amp_trials)
                        psth_curves{ai} = zeros(1, numel(ctrs));
                        continue;
                    end

                    counts = zeros(1, numel(edges)-1);

                    for tr = amp_trials(:)'
                        t0 = trig(tr)/FS*1000;
                        tt = sp_clipped{ch}(:,1);
                        tt = tt(tt >= t0+ras_win(1) & tt <= t0+ras_win(2)) - t0;

                        % raster
                        y0 = y_cursor + 1;
                        for spike_t = tt'
                            plot(ax1, [spike_t spike_t], ...
                                      [y0-0.4 y0+0.4], ...
                                      'Color', color, 'LineWidth', 1.1);
                        end

                        y_cursor = y_cursor + 1;

                        % PSTH
                        counts = counts + histcounts(tt, edges);
                    end

                    rate = filter(g, 1, counts / (numel(amp_trials)*bin_s));
                    psth_curves{ai} = rate;
                    maxRate = max(maxRate, max(rate));

                    ytick_vals(end+1) = y_cursor - numel(amp_trials)/2;
                    ytick_labels{end+1} = sprintf('%d µA', amp_val);
                end

                % Finalize raster
                xline(ax1,0,'r--');
                xlim(ax1, ras_win);
                ylim(ax1, [0 y_cursor]);
                yticks(ax1, ytick_vals);
                yticklabels(ax1, ytick_labels);
                ylabel(ax1, 'Amplitude');

                % PSTH
                for ai = 1:n_AMP
                    plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.6);
                end
                xline(ax2,0,'r--');
                xlim(ax2, ras_win);
                ylim(ax2, [0 max(50, ceil(maxRate*1.1/10)*10)]);
                xlabel(ax2, 'Time (ms)');
                ylabel(ax2, 'Rate (sp/s)');
                legend(ax2, arrayfun(@(a) sprintf('%.0f µA',a), Amps,'UniformOutput',false), ...
                       'Location','northeast','Box','off');
            end
        end
    end
end