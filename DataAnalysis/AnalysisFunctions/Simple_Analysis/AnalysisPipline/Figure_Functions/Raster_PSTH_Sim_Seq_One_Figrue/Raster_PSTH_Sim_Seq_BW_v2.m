%% ============================================================
%   Simultaneous vs Sequential: Final Conference Figure
%   - FIXED: Curve visibility near zero (removed clipping/white patch)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

% Plot Settings
save_figure = false;
save_dir = '/Users/xiaofan/Desktop/PhD Study/Conference/IEEE EMBC/Figures/2_Raster_PSTH_Sample';
fig_name    = 'GrandAverage_Recruitment_Curve_BarSEM.tiff';

%% =============== USER SETTINGS ==============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim4';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq4_5ms_new';
Electrode_Type   = 1;
target_channels  = [21];
plot_amp         = 5;     % µA
plot_PTD_ms      = 5;     % ms (Time of second pulse)
stim_set_id_Sim = 1;      
stim_set_id_Seq = 1;

% Figure Settings
ras_win   = [-30 30];
bin_ms    = 1;
smooth_ms = 5; 
FS        = 30000;
blank_artifact_win = [0 0]; 
stim_width_us = 100;   % the width of the stimulation line
stim_width_ms = stim_width_us / 1000; 
raster_color  = [0.2 0.2 0.2]; 
dash_width_ms = 0.5;   % length of the raster points        

%% ============================================================
%                 HELPER & DATA LOADING
% ============================================================
function sp = load_ssd_spikes(folder)
    cd(folder);
    f = dir('*sp_xia_SSD.mat');
    if isempty(f), error('No sp_xia_SSD.mat in %s', folder); end
    S = load(f(1).name);
    if isfield(S,'sp_corr'), sp=S.sp_corr;
    elseif isfield(S,'sp_SSD'), sp=S.sp_SSD;
    elseif isfield(S,'sp_in'), sp=S.sp_in;
    else, error('SSD file missing variables.'); end
end

% --- SIM ---
cd(folder_sim);
sp_sim = load_ssd_spikes(folder_sim);
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_sim = loadTrig(0);
Ssim = load(dir('*_exp_datafile_*.mat').name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim; E_MAP_sim = Ssim.E_MAP; nTr_sim = Ssim.n_Trials;
amps_all_sim = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;
stimNames = Stim_sim(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_sim(2:end));
comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim, rr = (t-1)*simN_sim + (1:simN_sim); v = idx_all(rr); v = v(v>0); comb_sim(t,1:numel(v)) = v(:).'; end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable'); nSets_sim = size(uniqueComb_sim,1);

% --- SEQ ---
cd(folder_seq);
sp_seq = load_ssd_spikes(folder_seq);
if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_seq = loadTrig(0);
Sseq = load(dir('*_exp_datafile_*.mat').name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP; nTr_seq = Sseq.n_Trials;
amps_all_seq = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000;
[PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);
stimNames = Stim_seq(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v = idx_all(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

% --- Condition Selection ---
ai_sim = find(Amps_sim == plot_amp,1);
ai_seq = find(Amps_seq == plot_amp,1);
pi_seq = find(PTDs_ms == plot_PTD_ms,1);
trials_sim = find(combClass_sim==stim_set_id_Sim & ampIdx_sim==ai_sim);
trials_seq = find(combClass_seq==stim_set_id_Seq & ampIdx_seq==ai_seq & ptdIdx_seq==pi_seq);

% --- Kernels ---
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2); g = g / sum(g);
win_sec = 5; b_sec = ones(1, win_sec) / win_sec; a_sec = 1;
d = Depth_s(Electrode_Type);

%% ============================================================
%                PLOTTING LOOP
% ============================================================
for ch_idx = 1:length(target_channels)
    target_channel = target_channels(ch_idx);
    recCh = d(target_channel);
    
    S_ch_sim = sp_sim{recCh};
    S_ch_seq = sp_seq{recCh};
    
    % Compute Rates
    allSpikes_sim = cell(numel(trials_sim),1); counts_sim = zeros(1,length(edges)-1);
    for i = 1:numel(trials_sim)
        t0 = trig_sim(trials_sim(i))/FS*1000; tt = S_ch_sim(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        allSpikes_sim{i} = tt; counts_sim = counts_sim + histcounts(tt,edges);
    end
    rate_sim_raw = filter(g, 1, counts_sim/(max(1,numel(trials_sim))*bin_s));
    rate_sim_sec = filter(b_sec, a_sec, rate_sim_raw);
    
    allSpikes_seq = cell(numel(trials_seq),1); counts_seq = zeros(1,length(edges)-1);
    for i = 1:numel(trials_seq)
        t0 = trig_seq(trials_seq(i))/FS*1000; tt = S_ch_seq(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        allSpikes_seq{i} = tt; counts_seq = counts_seq + histcounts(tt,edges);
    end
    rate_seq_raw = filter(g, 1, counts_seq/(max(1,numel(trials_seq))*bin_s));
    rate_seq_sec = filter(b_sec, a_sec, rate_seq_raw);
    
    mask_blank = ctrs >= blank_artifact_win(1) & ctrs <= blank_artifact_win(2);
    rate_sim_sec(mask_blank) = NaN; rate_seq_sec(mask_blank) = NaN;
    
    plot_overlay_figure(1, 'Final Result', rate_sim_sec, rate_seq_sec, ...
        allSpikes_sim, allSpikes_seq, ctrs, ras_win, stim_width_ms, ...
        plot_PTD_ms, raster_color, dash_width_ms, target_channel, plot_amp);
end

%% ============================================================
%   LOCAL PLOTTING FUNCTION
% ============================================================
function plot_overlay_figure(fig_num, type_str, r_sim, r_seq, sp_sim, sp_seq, ...
                             ctrs, ras_win, stim_width, ptd_ms, r_col, dash_w, ch, amp)
    
    figure('Color','w','Position',[100, 200, 700, 500]); 
    
    % Determine Max Rate
    max_rate = max([max(r_sim), max(r_seq)]) * 1.2; 
    if max_rate == 0 || isnan(max_rate), max_rate = 10; end
    
    % FIXED: Set Lower Limit slightly below zero so thick lines aren't clipped
    % y_lims = [-max_rate*0.02, max_rate];
    y_lims = [0, max_rate];

    % -----------------------------------------------------------
    % 1. LEFT AXIS SETUP (Backgrounds + PSTH)
    % -----------------------------------------------------------
    yyaxis left
    ylim(y_lims);
    hold on;
    
    % A. Draw Backgrounds (Scaled to rate limits)
    % patch([ras_win(1) ras_win(2) ras_win(2) ras_win(1)], ...
    %       [max_rate*0.5 max_rate*0.5 max_rate max_rate], ...
    %       [0.92 0.94 0.96], 'EdgeColor', 'none', 'HandleVisibility', 'off'); % Light Blue
    
    patch([ras_win(1) ras_win(2) ras_win(2) ras_win(1)], ...
          [max_rate*0.5 max_rate*0.5 max_rate max_rate], ...
          [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off'); % Light Grey
    
    % B. Draw Stimulation Bars
    patch([0, 0+stim_width, 0+stim_width, 0], [y_lims(1), y_lims(1), y_lims(2), y_lims(2)], ...
          [0.6 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    patch([ptd_ms, ptd_ms+stim_width, ptd_ms+stim_width, ptd_ms], ...
          [y_lims(1), y_lims(1), y_lims(2), y_lims(2)], ...
          [0.6 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');


    % C. Plot Curves
    p1 = plot(ctrs, r_sim, 'k--', 'LineWidth', 3.0, 'DisplayName', 'Simultaneous');
    p2 = plot(ctrs, r_seq, 'k-',  'LineWidth', 3.0, 'DisplayName', 'Sequential');
    
    ylabel('Firing rate (Sp/s)', 'FontWeight','bold','FontSize',18,'Color','k', 'FontName', 'Times New Roman');
    xlabel('Time (ms)', 'FontWeight','bold','FontSize',18,'Color','k', 'FontName', 'Times New Roman');
    xlim(ras_win);
    xticks(ras_win(1):5:ras_win(2)); % Sets ticks every 5ms across the whole window
    % -----------------------------------------------------------
    % 2. RIGHT AXIS SETUP (Rasters)
    % -----------------------------------------------------------
    yyaxis right
    
    % CRITICAL: Make Right Axis Transparent so Left Axis shows through
    ax_r = gca; 
    set(ax_r, 'Color', 'none'); 
    set(ax_r, 'YColor', 'none'); 
    set(ax_r, 'XColor', 'k');    
    
    ylim([0 1]); 
    ylabel(''); 
    hold on;
    
    % Raster Plotting
    y_seq_min = 0.05; y_seq_max = 0.45;
    y_sim_min = 0.55; y_sim_max = 0.95;

    n_seq = numel(sp_seq);
    if n_seq > 0
        h_tick = (y_seq_max - y_seq_min) / n_seq;
        num_spikes = sum(cellfun(@length, sp_seq));
        x_v = zeros(1, num_spikes*3); y_v = zeros(1, num_spikes*3);
        idx = 1;
        for i = 1:n_seq
            tt = sp_seq{i}; if isempty(tt), continue; end; tt = tt(:)'; 
            y_pos = y_seq_min + (i-1)*h_tick;
            for t_val = tt
                 x_v(idx:idx+2) = [t_val-dash_w/2, t_val+dash_w/2, NaN];
                 y_v(idx:idx+2) = [y_pos, y_pos, NaN];
                 idx = idx + 3;
            end
        end
        x_v = x_v(1:idx-1); y_v = y_v(1:idx-1);
        plot(x_v, y_v, '-', 'Color', r_col, 'LineWidth', 0.8); 
    end
    
    n_sim = numel(sp_sim);
    if n_sim > 0
        h_tick = (y_sim_max - y_sim_min) / n_sim;
        num_spikes = sum(cellfun(@length, sp_sim));
        x_v = zeros(1, num_spikes*3); y_v = zeros(1, num_spikes*3);
        idx = 1;
        for i = 1:n_sim
            tt = sp_sim{i}; if isempty(tt), continue; end; tt = tt(:)';
            y_pos = y_sim_min + (i-1)*h_tick;
            for t_val = tt
                 x_v(idx:idx+2) = [t_val-dash_w/2, t_val+dash_w/2, NaN];
                 y_v(idx:idx+2) = [y_pos, y_pos, NaN];
                 idx = idx + 3;
            end
        end
        x_v = x_v(1:idx-1); y_v = y_v(1:idx-1);
        plot(x_v, y_v, '-', 'Color', r_col, 'LineWidth', 0.8);
    end
    
    % Text Labels
    text(ras_win(1)+2, 0.95, 'Simultaneous', 'VerticalAlignment','top','FontSize',18,'FontWeight','bold','Color','k', 'FontName', 'Times New Roman');
    text(ras_win(1)+2, 0.45, 'Sequential', 'VerticalAlignment','top','FontSize',18,'FontWeight','bold','Color','k', 'FontName', 'Times New Roman');

    % -----------------------------------------------------------
    % 3. FINAL VISIBILITY RESTORATION
    % -----------------------------------------------------------
    yyaxis left
    set(gca, 'Box', 'off');     
    set(gca, 'Color', 'none');  
    set(gca, 'Layer', 'top');   
    set(gca, 'YColor', 'k');    
    set(gca, 'XColor', 'k');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 16, 'LineWidth', 3);
    % set(gca, 'Units', 'normalized', 'Position', [0.13 0.13 0.85 0.82]);
    % title(sprintf('Ch %d | %d µA', ch, amp), 'Color', 'k', 'FontSize', 14);
    legend([p1, p2], 'Location', 'northeast', 'Box', 'off', 'FontSize', 18);
    axis square;
    hold off;
end

%% Save Figure
if save_figure
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    exportgraphics(gcf, fullfile(save_dir, fig_name),'ContentType', 'vector');
    fprintf('\nFigure saved to: %s\n', fullfile(save_dir, fig_name));
end