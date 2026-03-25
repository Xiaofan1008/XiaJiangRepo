%% Strobe (Flash) Raster & PSTH Viewer
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX018/strobe_LFP_260319_102436'; 

raster_chn_start = 1;
raster_chn_end   = 64; % nChn
Electrode_Type   = 2;  % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% [MODIFIED] Artifact Blanking Settings
artifact_centers = [2, 34, 68, 102, 135, 170]; % Target timestamps to erase (ms)
blanking_window  = 1;                        % Tolerance (+/- ms around center)

% Raster Plot Parameters (Expanded for visual cortex response)
ras_win         = [-50 200];   % ms (You can change this window length here)
bin_ms_raster   = 2;           % bin size
smooth_ms       = 10;          % smoothing window

%% ====================== LOAD DATA ======================
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
fprintf('Changed directory to:\n%s\n', data_folder);

FS = 30000; % Sampling frequency

% Load .sp.mat file
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

% Load Triggers
if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0); 
nTr = length(trig); % The total number of strobe flashes is your trial count

% Electrode Map
d = Depth_s(Electrode_Type); 
sp_clipped = sp; 

%% ====================== RASTER PLOT SETUP ======================
edges = ras_win(1):bin_ms_raster:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms_raster / 1000;

% Gaussian smoothing kernel
g = exp(-0.5 * ((0:smooth_ms-1) / (smooth_ms/2)).^2);
g = g / sum(g);

%% ====================== MAIN PLOTTING LOOP ======================
for ich = raster_chn_start:raster_chn_end
    ch = d(ich);
    
    % Skip if channel has no data
    if isempty(sp_clipped{ch}), continue; end
    
    % Create Figure
    figure('Color','w','Name',sprintf('Strobe Response - Ch %d', ich), 'Position', [100 100 800 600]);
    tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
    
    ax1 = nexttile([3 1]); hold(ax1,'on'); box(ax1,'off');
    ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');
    
    counts = zeros(1, numel(edges)-1);
    total_spikes = 0;
    
    % Loop through every strobe flash (Trial)
    for t = 1:nTr
        S_ch = sp_clipped{ch};
        t0 = trig(t) / FS * 1000;
        
        % Extract spikes within the chosen window around this trigger
        tt = S_ch(:,1);
        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
        
        % [MODIFIED] Scrub out the specific artifact timestamps
        valid_spikes_mask = true(size(tt)); % Assume all are real spikes initially
        
        for a = 1:length(artifact_centers)
            center = artifact_centers(a);
            % Find any spikes hiding inside the blanking window for this center
            is_artifact = (tt >= (center - blanking_window)) & (tt <= (center + blanking_window));
            valid_spikes_mask(is_artifact) = false; % Mark them for deletion
        end
        
        % Keep only the true spikes
        tt = tt(valid_spikes_mask);
        
        % Add to histogram counts
        counts = counts + histcounts(tt, edges);
        total_spikes = total_spikes + numel(tt);
        
        % Plot raster ticks for this trial
        for k = 1:numel(tt)
            plot(ax1, [tt(k) tt(k)], [t-0.4 t+0.4], 'Color', 'k', 'LineWidth', 1.2);
        end
    end
    
    % --- Finalize Raster Plot (Top) ---
    xline(ax1, 0, 'r', 'LineWidth', 1.5);
    xlim(ax1, ras_win);
    ylim(ax1, [0, nTr + 1]);
    ylabel(ax1, 'Trial Number', 'FontWeight', 'bold');
    title(ax1, sprintf('Strobe Raster — Ch %d | Total Flashes: %d', ich, nTr), 'Interpreter','none', 'FontWeight', 'bold');
    set(ax1, 'XTickLabel', []); % Hide X-labels to prevent crowding the PSTH
    
    % --- Finalize PSTH Plot (Bottom) ---
    rate = counts / (nTr * bin_s);
    rate_s = filter(g, 1, rate); % Apply smoothing
    
    plot(ax2, ctrs, rate_s, 'Color', 'k', 'LineWidth', 1.5);
    
    xline(ax2, 0, 'r', 'LineWidth', 1.5);
    xlim(ax2, ras_win);
    
    % Dynamically set Y-limits so it doesn't crash if the channel is completely dead
    maxRate = max(rate_s);
    if maxRate > 0
        ylim(ax2, [0, ceil(maxRate*1.1/10)*10 + 1]); 
    else
        ylim(ax2, [0 1]);
    end
    
    xlabel(ax2, 'Time (ms)', 'FontWeight', 'bold');
    ylabel(ax2, 'Rate (sp/s)', 'FontWeight', 'bold');
end