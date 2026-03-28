%% Strobe (Flash) Raster & PSTH Viewer (Overlay Format for IEEE)
clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER SETTINGS ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX019/strobe_260320_110254'; 
raster_chn_start = 52;
raster_chn_end   = 52; 
Electrode_Type   = 2;  

% Time-Blanking Artifact Settings
artifact_centers = [-170, -136, -102, -68, -34, 0, 2, 34, 68, 102, 136, 170]; 
blanking_window  = 1;                        

% Waveform Amplitude Filtering Settings
amp_min = 20;   
amp_max = 200;  
abs_max = 200;  

% Raster Plot Parameters 
ras_win = [-200 200];   % ms 
bin_size_ms  = 1;      
smooth_sigma = 5;      

% Saving and Visibility Settings
show_title = true;
save_figs = false;      
show_figs = true;     
save_dir  = '/Volumes/MACData/Data/Data_Xia/Figures/IEEE_TBME_Paper/Visul_Stim_Raster/DX019/Strobe'; 

% [MODIFIED 1] Manual Trial Rejection (Blacklist)
% Type the specific trial numbers you want to delete inside the brackets for any given channel.
bad_trials_list = cell(64, 1);
% Example: bad_trials_list{24} = [12, 45, 108]; 
bad_trials_list{52} = [1:60,440:540]; % Put bad trial numbers for Ch 24 here
bad_trials_list{31} = []; % Put bad trial numbers for Ch 31 here

%% ====================== LOAD DATA ======================
if ~isfolder(data_folder)
    error('The specified folder does not exist. Please check the path.');
end
cd(data_folder);
FS = 30000; 
sp_files = dir(fullfile(data_folder, '*.sp.mat'));
sp_filename = fullfile(data_folder, sp_files(1).name);
S = load(sp_filename);
sp = S.sp;

if isempty(dir(fullfile(data_folder, '*.trig.dat')))
    cur_dir = pwd; cd(data_folder);
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0); 
nTr = length(trig); 
d = Depth_s(Electrode_Type); 
sp_clipped = sp; 

if save_figs && ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% ====================== MAIN PLOTTING LOOP ======================
for ich = raster_chn_start:raster_chn_end
    ch = d(ich);
    if isempty(sp_clipped{ch}), continue; end
    
    waveforms = sp_clipped{ch}(:, 2:end); 
    ptp_amp   = max(waveforms, [], 2) - min(waveforms, [], 2); 
    
    valid_spikes = (ptp_amp >= amp_min) & (ptp_amp <= amp_max) & ...
                   all(waveforms <= abs_max, 2) & all(waveforms >= -abs_max, 2);
    S_ch_filtered = sp_clipped{ch}(valid_spikes, :); 
    
    if isempty(S_ch_filtered), continue; end
    
    % [MODIFIED 2] Apply the Blacklist for this specific channel
    bad_trials = bad_trials_list{ich};
    valid_trials = setdiff(1:nTr, bad_trials); % Creates a list of only the good trials
    N_valid = length(valid_trials); % The new total trial count
    
    % Skip plotting if you accidentally blacklisted every single trial
    if N_valid == 0, continue; end 
    
    if show_figs, fig_vis = 'on'; else, fig_vis = 'off'; end
    
    fig = figure('Color','w','Name',sprintf('Strobe - Ch %d', ich), 'Visible', fig_vis);
    set(fig, 'Units', 'centimeters');
    set(fig, 'Position', [2, 2, 8.89, 8.89]); 
    set(fig, 'PaperPositionMode', 'auto'); 
             
    ax = axes(fig); 
    hold(ax, 'on');
    
    rel_spike_times = [];
    trial_numbers   = [];
    
    % [MODIFIED 3] Loop only through valid trials and remap the Y-axis so there are no blank gaps
    for t_idx = 1:N_valid
        t = valid_trials(t_idx); % Get the actual original trial number to pull the right trigger
        
        t0 = trig(t) / FS * 1000;
        tt = S_ch_filtered(:,1);
        tt = tt(tt >= t0 + ras_win(1) & tt <= t0 + ras_win(2)) - t0;
        
        valid_spikes_mask = true(size(tt)); 
        for a = 1:length(artifact_centers)
            center = artifact_centers(a);
            is_artifact = (tt >= (center - blanking_window)) & (tt <= (center + blanking_window));
            valid_spikes_mask(is_artifact) = false; 
        end
        tt = tt(valid_spikes_mask);
        
        rel_spike_times = [rel_spike_times; tt(:)];
        trial_numbers   = [trial_numbers; t_idx * ones(numel(tt), 1)]; % Use t_idx to prevent vertical gaps
    end
    
    % --- Calculate Smoothed PSTH ---
    edges = ras_win(1) : bin_size_ms : ras_win(2);
    bin_centers = edges(1:end-1) + bin_size_ms/2;
    spike_counts = histcounts(rel_spike_times, edges);
    
    % [MODIFIED 4] Update PSTH Math to divide by N_valid instead of nTr
    fr_raw = spike_counts / (bin_size_ms / 1000) / N_valid;
    fr_smoothed = smoothdata(fr_raw, 'gaussian', smooth_sigma * 5);
    
    yyaxis right
    
    N_spikes = length(rel_spike_times);
    dash_w = 0.8; 
    x_v = zeros(3 * N_spikes, 1); y_v = zeros(3 * N_spikes, 1);
    x_v(1:3:end) = rel_spike_times; x_v(2:3:end) = rel_spike_times; x_v(3:3:end) = NaN;
    y_v(1:3:end) = trial_numbers - dash_w/2; y_v(2:3:end) = trial_numbers + dash_w/2; y_v(3:3:end) = NaN;
    
    plot(x_v, y_v, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    
    % [MODIFIED 5] Update Y-limit and N text label to use N_valid
    ylim([0, N_valid + 1]); 
    set(gca, 'YColor', 'none'); 
    
    text(ras_win(2)*0.95, N_valid*0.95, sprintf('N = %d', N_valid), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'FontName', 'Arial', 'FontWeight', 'bold', 'Color', 'k');

    yyaxis left
    
    plot(bin_centers, fr_smoothed, 'k-', 'LineWidth', 1.5);
    xline(0, '--k', 'LineWidth', 1.5);
    
    maxRate = max(fr_smoothed);
    if maxRate > 0
        ylim([0, ceil(maxRate*1.2/10)*10]); 
    else
        ylim([0 10]);
    end
    xlim(ras_win);
    
    ylabel('Firing rate (sp/s)', 'FontWeight', 'bold', 'Color', 'k');
    xlabel('Time (ms)', 'FontWeight', 'bold', 'Color', 'k');
    
    set(gca, 'YColor', 'k'); 
    set(gca, 'XColor', 'k');
    
    axis square;
    box off;
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
    if show_title
        title(sprintf('Ch %d Response', ich), 'FontWeight', 'normal');
    end
    
    if save_figs
        save_name = fullfile(save_dir, sprintf('Strobe_raster_Ch%d.tif', ich));
        print(fig, save_name, '-dtiff', '-r300'); 
        fprintf('Saved: %s\n', save_name);
    end
    if ~show_figs, close(fig); end
end
fprintf('Done processing channels.\n');