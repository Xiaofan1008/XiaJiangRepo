clear all
% close all
% addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% Choose Folder

% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Single5_251014_184742';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX009/Xia_Exp1_Sim5_251014_183532';
data_folder = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_600_700um_SimSeq1';

%% Choice
raster_chn_start = 48;
raster_chn_end = 64; %nChn
Electrode_Type = 3; % 0:single shank rigid; 1:single shank flex; 2:four shank flex
PTD_to_plot = [];   % e.g., [500 1000], empty for all PTD
PTD_to_plot = PTD_to_plot.*1000;

%% Load folder
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

%% Pre Set
FS=30000; % Sampling frequency
% Load .sp_xia.mat file (produced directly by the spike detection script)
sp_files = dir(fullfile(data_folder, '*.sp_xia.mat'));
assert(~isempty(sp_files), 'No .sp_xia.mat file found in the current folder.');
sp_filename = fullfile(data_folder, sp_files(1).name);
fprintf('Loading spike file: %s\n', sp_filename);
S = load(sp_filename);
if isfield(S, 'sp_clipped')
    sp_clipped = S.sp_clipped;
else
    error('Variable "sp_clipped" not found in %s.', sp_filename);
end

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);

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
    % KEEP ORDER — do NOT use unique()
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

% Pulse Train Period (inter-pulse interval)
pulseTrain_all = cell2mat(StimParams(2:end,9));  % Column 9: Pulse Train Period
pulseTrain = pulseTrain_all(1:simultaneous_stim:end);  % take 1 per trial
[PulsePeriods, ~, pulseIdx] = unique(pulseTrain(:));
n_PULSE = numel(PulsePeriods);

% Extract POST-TRIGGER DELAY (PTD)
% Column 6 stores PTD for the 2nd pulse of each trial block
% For sequential stimulation: row 2 = delayed pulse
if simultaneous_stim > 4
    ptd_all = cell2mat(StimParams(6:simultaneous_stim:end, 6));
elseif simultaneous_stim == 4
    ptd_all = cell2mat(StimParams(5:simultaneous_stim:end, 6));
elseif simultaneous_stim == 3
    ptd_all = cell2mat(StimParams(4:simultaneous_stim:end, 6));
elseif simultaneous_stim == 2
    ptd_all = cell2mat(StimParams(3:simultaneous_stim:end, 6));  % µs
else
    ptd_all = zeros(n_Trials,1); % Single-pulse case, no PTD
end

PTD_us = ptd_all(:);
[PTD_values, ~, ptdIdx] = unique(PTD_us);

if isempty(PTD_to_plot)
    PTD_selected = PTD_values;
else
    PTD_selected = intersect(PTD_values, PTD_to_plot);
end
fprintf('\nDetected PTDs (µs):'); disp(PTD_values');
fprintf('PTDs selected for plotting:'); disp(PTD_selected');
n_PTD = numel(PTD_selected);
% Electrode Map
d = Depth_s(Electrode_Type); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex

%% Raster Plot Parameters
ras_win         = [-20 100];   % ms
bin_ms_raster   = 1;           % bin size
smooth_ms       = 2;           % smoothing window
% raster_chn_start = 1;
% raster_chn_end = 32; %nChn

%% === Initialize structure to store first-spike times ===
firstSpikeTimes = cell(raster_chn_end, 1); % each cell: vector of first-spike times per trial (ms)
fprintf('\nComputing First Spike Times per Trial\n');
post_spike_window_ms = [5,8];

%% Raster Plot
 %% ========================= RASTER + PSTH ========================= %%
edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

for ich = raster_chn_start:raster_chn_end
    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end

    for si = 1:nSets
        stimVec = uniqueComb(si, :);
        stimVec = stimVec(stimVec > 0);
        setLabel = strjoin(arrayfun(@(x) sprintf('Ch%d', x), stimVec, 'UniformOutput', false), '→');

        for pi = 1:n_PULSE
            pulse_val = PulsePeriods(pi);

            for i_PTD = 1:n_PTD
                ptd_val = PTD_selected(i_PTD);

                % --------- FIXED: restrict to THIS PTD ---------
                trials_this_period = find( combClass_win == si & ...
                                           pulseIdx == pi & ...
                                           ampIdx >= 1 & ...   % any amp
                                           PTD_us == ptd_val );  % << FIXED
                if isempty(trials_this_period), continue; end

                % --------- FIGURE ---------
                figName = sprintf('Ch %d | Set %s | Pulse %d µs | PTD %d µs', ...
                                  ich, setLabel, pulse_val, ptd_val);
                figure('Color','w','Name',figName);
                tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

                ax1 = nexttile([3 1]);
                hold(ax1,'on'); box(ax1,'off');
                title(ax1, sprintf('Raster — Ch %d | Set %s | PTD %d µs', ...
                                   ich, setLabel, ptd_val), 'Interpreter','none');

                ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');

                % ===== PSTH storage =====
                psth_curves = cell(1, n_AMP);
                maxRate = 0;
                y_cursor = 0;
                ytick_vals = [];
                ytick_labels = {};

                % ================= LOOP AMPLITUDES ==================
                for ai = 1:n_AMP
                    amp_val = Amps(ai);
                    color = cmap(ai,:);

                    % -------- FIXED TRIAL FILTERING --------
                    amp_trials = find( ampIdx == ai & ...
                                       pulseIdx == pi & ...
                                       combClass_win == si & ...
                                       PTD_us == ptd_val );   % << FIXED

                    nTr = numel(amp_trials);
                    if nTr == 0
                        psth_curves{ai} = zeros(1, numel(ctrs));
                        continue;
                    end

                    % ======== RASTER + PSTH COUNTING ========
                    counts = zeros(1, numel(edges)-1);
                    for t = 1:nTr
                        tr = amp_trials(t);
                        t0 = trig(tr)/FS*1000;

                        tt = sp_clipped{ch}(:,1);
                        tt = tt(tt >= t0+ras_win(1) & tt <= t0+ras_win(2)) - t0;

                        % ---- RASTER ----
                        y0 = y_cursor + t;
                        for spike_t = tt'
                            plot(ax1, [spike_t spike_t], [y0-0.4 y0+0.4], 'Color', color, 'LineWidth', 1.1);
                        end

                        % ---- PSTH ----
                        counts = counts + histcounts(tt, edges);
                    end

                    % y-axis structure
                    ytick_vals(end+1) = y_cursor + nTr/2;
                    ytick_labels{end+1} = sprintf('%d µA', amp_val);
                    y_cursor = y_cursor + nTr;

                    % ---- Compute PSTH ----
                    rate = filter(g, 1, counts/(nTr*bin_s));
                    psth_curves{ai} = rate;
                    maxRate = max(maxRate, max(rate));
                end

                % Finalize raster axis
                xline(ax1, 0, 'r--');
                xlim(ax1, ras_win);
                ylim(ax1, [0 y_cursor]);
                yticks(ax1, ytick_vals);
                yticklabels(ax1, ytick_labels);
                ylabel(ax1, 'Amplitude');

                % Finalize PSTH
                for ai = 1:n_AMP
                    plot(ax2, ctrs, psth_curves{ai}, 'Color', cmap(ai,:), 'LineWidth', 1.6);
                end
                xline(ax2, 0, 'r--');
                xlim(ax2, ras_win);
                ylim(ax2, [0 max(50,ceil(maxRate*1.1/10)*10)]);
                xlabel(ax2, 'Time (ms)');
                ylabel(ax2, 'Rate (sp/s)');
                legend(ax2, arrayfun(@(a) sprintf('%.0f µA',a), Amps, 'UniformOutput', false), ...
                       'Box','off','Location','northeast');

            end % PTD
        end % Pulse
    end % Set
end % Channel