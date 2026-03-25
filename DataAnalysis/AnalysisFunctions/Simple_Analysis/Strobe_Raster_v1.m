%% Strobe (Flash) Raster & PSTH Viewer
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX018/strobe_LFP_2_260319_103140'; 
raster_chn_start = 32;
raster_chn_end   = 64; % nChn
Electrode_Type   = 2;  %j 0:single shank rigid; 1:single shank flex; 2:four shank flex

% Time-Blanking Artifact Settings
artifact_centers = [-34, 0, 2, 34, 68, 102, 136, 170]; % Target timestamps to erase (ms)
blanking_window  = 1;                        % Tolerance (+/- ms around center)

% Waveform Amplitude Filtering Settings
amp_min = 20;   % Minimum peak-to-peak amplitude (µV) to be considered a real spike
amp_max = 200;  % Maximum peak-to-peak amplitude (µV) (removes giant artifacts)
abs_max = 200;  % Absolute maximum voltage limit (+/- µV)

% Raster Plot Parameters (Expanded for visual cortex response)
ras_win         = [-50 200];   % ms (You can change this window length here)

% Sliding PSTH Parameters
win_ms  = 4;    % Sliding window size (ms)
step_ms = 1;    % Sliding step size (ms)

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

%% ====================== MAIN PLOTTING LOOP ======================
for ich = raster_chn_start:raster_chn_end
    ch = d(ich);
    
    % Skip if channel has no data
    if isempty(sp_clipped{ch}), continue; end
    
    % [MODIFIED] Idea 1: Apply Amplitude Filter Upfront
    waveforms = sp_clipped{ch}(:, 2:end); % Get just the waveform voltages
    ptp_amp   = max(waveforms, [], 2) - min(waveforms, [], 2); % Peak-to-Peak
    
    % Find spikes that pass the criteria
    valid_spikes = (ptp_amp >= amp_min) & (ptp_amp <= amp_max) & ...
                   all(waveforms <= abs_max, 2) & all(waveforms >= -abs_max, 2);
               
    S_ch_filtered = sp_clipped{ch}(valid_spikes, :); % Keep only the clean spikes
    
    % Skip if filtering removed everything
    if isempty(S_ch_filtered), continue; end
    
    % Create Figure
    figure('Color','w','Name',sprintf('Strobe Response - Ch %d', ich), 'Position', [100 100 800 600]);
    tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
    
    ax1 = nexttile([3 1]); hold(ax1,'on'); box(ax1,'off');
    ax2 = nexttile; hold(ax2,'on'); box(ax2,'off');
    
    % [MODIFIED] Idea 3: Pre-allocate arrays for optimized fast-plotting
    rel_spike_times = [];
    trial_numbers   = [];
    
    % Loop through every strobe flash (Trial)
    for t = 1:nTr
        t0 = trig(t) / FS * 1000;
        
        % Extract spikes within the chosen window around this trigger
        tt = S_ch_filtered(:,1);
        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
        
        % [MODIFIED] Time-Blanking: Scrub out the specific artifact timestamps
        valid_spikes_mask = true(size(tt)); % Assume all are real spikes initially
        
        for a = 1:length(artifact_centers)
            center = artifact_centers(a);
            % Find any spikes hiding inside the blanking window for this center
            is_artifact = (tt >= (center - blanking_window)) & (tt <= (center + blanking_window));
            valid_spikes_mask(is_artifact) = false; % Mark them for deletion
        end
        
        % Keep only the true spikes
        tt = tt(valid_spikes_mask);
        
        % [MODIFIED] Idea 3: Collect data for the scatter plot instead of drawing lines
        rel_spike_times = [rel_spike_times; tt(:)];
        trial_numbers   = [trial_numbers; t * ones(numel(tt), 1)];
    end
    
    % --- Finalize Raster Plot (Top) ---
    % [MODIFIED] Idea 3: Draw all spikes instantly as dots
    plot(ax1, rel_spike_times, trial_numbers, 'k.', 'MarkerSize', 4);
    
    xline(ax1, 0, 'r', 'LineWidth', 1.5);
    xlim(ax1, ras_win);
    ylim(ax1, [0, nTr + 1]);
    ylabel(ax1, 'Trial Number', 'FontWeight', 'bold');
    title(ax1, sprintf('Strobe Raster — Ch %d | Total Flashes: %d', ich, nTr), 'Interpreter','none', 'FontWeight', 'bold');
    set(ax1, 'XTickLabel', []); % Hide X-labels to prevent crowding the PSTH
    
    % --- Finalize PSTH Plot (Bottom) ---
    % [MODIFIED] Idea 2: Sliding Window PSTH Math
    edges_start = ras_win(1) : step_ms : (ras_win(2) - win_ms);
    nBins = length(edges_start);
    fr_sliding = zeros(1, nBins);
    
    for b = 1:nBins
        w_start = edges_start(b);
        w_end   = w_start + win_ms;
        % Count spikes in this specific window across all trials
        nSpikes_in_win = sum(rel_spike_times >= w_start & rel_spike_times < w_end);
        % Convert to Rate (spikes per second per trial)
        fr_sliding(b) = nSpikes_in_win / (win_ms / 1000) / nTr;
    end
    
    % Plot the sliding PSTH
    plot(ax2, edges_start + win_ms/2, fr_sliding, 'k', 'LineWidth', 1.5);
    
    xline(ax2, 0, 'r', 'LineWidth', 1.5);
    xlim(ax2, ras_win);
    
    % Dynamically set Y-limits so it doesn't crash if the channel is completely dead
    maxRate = max(fr_sliding);
    if maxRate > 0
        ylim(ax2, [0, ceil(maxRate*1.1/10)*10 + 1]); 
    else
        ylim(ax2, [0 1]);
    end
    
    xlabel(ax2, 'Time (ms)', 'FontWeight', 'bold');
    ylabel(ax2, 'Rate (sp/s)', 'FontWeight', 'bold');
end