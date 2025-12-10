%% ============================================================
%   Simultaneous vs Sequential Raster + PSTH (multiple channels)
% ============================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% =============== USER SETTINGS ==============================

folder_sim = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim1_251125_112055';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq1_5ms_251125_112735';

Electrode_Type   = 1;
target_channels  = [1:32];
plot_amp         = 10;     % µA
plot_PTD_ms      = 5;      % ms
stim_set_id_Sim = 1;      % stimulation set
stim_set_id_Seq = 1;

ras_win   = [-20 20];
bin_ms    = 1;
smooth_ms = 2;
FS        = 30000;

% Softer colours for paper
sim_col = [0.25 0.45 0.75];
seq_col = [0.75 0.35 0.35];

%% ============================================================
%                 HELPER FUNCTION 
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

%% ============================================================
%              LOAD SIMULTANEOUS DATA (SSD)
% ============================================================
cd(folder_sim);
sp_sim = load_ssd_spikes(folder_sim);

if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_sim = loadTrig(0);

Ssim = load(dir('*_exp_datafile_*.mat').name, ...
            'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim_sim = Ssim.StimParams;
simN_sim  = Ssim.simultaneous_stim;
E_MAP_sim = Ssim.E_MAP;
nTr_sim   = Ssim.n_Trials;

amps_all_sim  = cell2mat(Stim_sim(2:end,16));
trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim);
Amps_sim(Amps_sim==-1) = 0;

% Order-sensitive stim sets
stimNames = Stim_sim(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP_sim(2:end));
stimChPerTrial = cell(nTr_sim,1);
for t = 1:nTr_sim
    rr = (t-1)*simN_sim + (1:simN_sim);
    v  = idx_all(rr); v = v(v>0);
    stimChPerTrial{t} = v(:).';
end
comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim
    v = stimChPerTrial{t};
    comb_sim(t,1:numel(v)) = v;
end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable');
nSets_sim = size(uniqueComb_sim,1);

%% ============================================================
%              LOAD SEQUENTIAL DATA (SSD)
% ============================================================
cd(folder_seq);
sp_seq = load_ssd_spikes(folder_seq);

if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
trig_seq = loadTrig(0);

Sseq = load(dir('*_exp_datafile_*.mat').name, ...
            'StimParams','simultaneous_stim','E_MAP','n_Trials');
Stim_seq = Sseq.StimParams;
simN_seq  = Sseq.simultaneous_stim;  
E_MAP_seq = Sseq.E_MAP;
nTr_seq   = Sseq.n_Trials;

amps_all_seq  = cell2mat(Stim_seq(2:end,16));
trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq);
Amps_seq(Amps_seq==-1) = 0;

% PTD
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6));
PTD_all_ms = PTD_all_us / 1000;
[PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

% Ordered sets
stimNames = Stim_seq(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP_seq(2:end));
stimChPerTrial = cell(nTr_seq,1);
for t = 1:nTr_seq
    rr = (t-1)*simN_seq + (1:simN_seq);
    v  = idx_all(rr); v = v(v>0);
    stimChPerTrial{t} = v(:).';
end
comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq
    v = stimChPerTrial{t};
    comb_seq(t,1:numel(v)) = v;
end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable');
nSets_seq = size(uniqueComb_seq,1);

%% ============================================================
%          CHECK REQUESTED CONDITION
% ============================================================
if stim_set_id_Sim > nSets_sim || stim_set_id_Seq > nSets_seq
    error('stim_set_id exceeds available sets.');
end

ai_sim = find(Amps_sim == plot_amp,1);
ai_seq = find(Amps_seq == plot_amp,1);
pi_seq = find(PTDs_ms == plot_PTD_ms,1);

trials_sim = find(combClass_sim==stim_set_id_Sim & ampIdx_sim==ai_sim);
trials_seq = find(combClass_seq==stim_set_id_Seq & ampIdx_seq==ai_seq & ptdIdx_seq==pi_seq);

%% ============================================================
%            PSTH Kernel
% ============================================================
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;

g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g / sum(g);

d = Depth_s(Electrode_Type);

%% ============================================================
%                PLOTTING LOOP
% ============================================================

for ch_idx = 1:length(target_channels)

    target_channel = target_channels(ch_idx);
    recCh = d(target_channel);

    S_ch_sim = sp_sim{recCh};
    S_ch_seq = sp_seq{recCh};

    %% --- Compute SIM PSTH ---
    allSpikes_sim = cell(numel(trials_sim),1);
    counts_sim = zeros(1,length(edges)-1);
    for i = 1:numel(trials_sim)
        t0 = trig_sim(trials_sim(i))/FS*1000;
        tt = S_ch_sim(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        allSpikes_sim{i} = tt;
        counts_sim = counts_sim + histcounts(tt,edges);
    end
    rate_sim = filter(g,1, counts_sim/(max(1,numel(trials_sim))*bin_s));

    %% --- Compute SEQ PSTH ---
    allSpikes_seq = cell(numel(trials_seq),1);
    counts_seq = zeros(1,length(edges)-1);
    for i = 1:numel(trials_seq)
        t0 = trig_seq(trials_seq(i))/FS*1000;
        tt = S_ch_seq(:,1);
        tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
        allSpikes_seq{i} = tt;
        counts_seq = counts_seq + histcounts(tt,edges);
    end
    rate_seq = filter(g,1, counts_seq/(max(1,numel(trials_seq))*bin_s));

    %% ============================================================
    %                   FIGURE FOR THIS CHANNEL
    %% ============================================================

    stimCh_sim = uniqueComb_sim(stim_set_id_Sim, uniqueComb_sim(stim_set_id_Sim,:)>0);
    stimCh_seq = uniqueComb_seq(stim_set_id_Seq, uniqueComb_seq(stim_set_id_Seq,:)>0);

    fig = figure('Color','w','Position',[200 200 520 680]);
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    title(tl, sprintf('Ch %d | %.0f µA | PTD %.0f ms | Set %d', ...
        target_channel, plot_amp, plot_PTD_ms, stim_set_id_Seq), ...
        'FontSize',12,'FontWeight','bold');

   % --- RASTER SIM ---
        ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'off');
        
        nTrials = numel(trials_sim);
        
        for t = 1:nTrials
            tt = allSpikes_sim{t};
            y  = t;   % TRUE TRIAL NUMBER
            if ~isempty(tt)
                plot(ax1, tt, y*ones(size(tt)), '.', 'Color', sim_col, 'MarkerSize',5);
            end
        end
        
        xline(ax1,0,'k-');
        ylabel(ax1,'Trials');
        xlim(ax1,ras_win);
        ylim(ax1,[0 nTrials+1]);  
        title(ax1, sprintf('Simultaneous — Stim Ch: %s', mat2str(stimCh_sim)));
        
        % Reduce height
        pos = get(ax1,'Position');
        pos(4) = pos(4)*0.55;
        set(ax1,'Position',pos);

   % --- RASTER SEQ ---
    ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'off');    
    nTrials = numel(trials_seq);    
    for t = 1:nTrials
        tt = allSpikes_seq{t};
        y  = t;  
        if ~isempty(tt)
            plot(ax2, tt, y*ones(size(tt)), '.', 'Color', seq_col, 'MarkerSize',5);
        end
    end    
    xline(ax2,0,'k-');
    ylabel(ax2,'Trials');
    xlim(ax2,ras_win);
    ylim(ax2,[0 nTrials+1]);
    title(ax2, sprintf('Sequential (Inter-Stimulus Interval %.0f ms) — Stim Ch: %s', ...
        plot_PTD_ms, mat2str(stimCh_seq)));    
    % Reduce height
    pos = get(ax2,'Position');
    pos(4) = pos(4)*0.55;
    set(ax2,'Position',pos);

   % --- PSTH ---
    ax3 = nexttile(tl); hold(ax3,'on'); box(ax3,'off');
    
    plot(ax3, ctrs, rate_sim, 'Color', sim_col, 'LineWidth',1.8);
    plot(ax3, ctrs, rate_seq, 'Color', seq_col, 'LineWidth',1.8);
    
    xline(ax3,0,'k-','HandleVisibility','off');    
    ylabel(ax3,'Rate (sp/s)');
    xlabel(ax3,'Time (ms)');
    xlim(ax3,ras_win);
    title(ax3,'PSTH');    
    legend(ax3, {'Simultaneous','Sequential'}, ...
           'Location','northeast','Box','off','FontSize',9);    
    % Reduce height
    pos = get(ax3,'Position');
    pos(4) = pos(4)*0.55;
    set(ax3,'Position',pos);

   %% ============================================================
    %   FIGURE STYLE: Split-Zone Raster (Background) + Overlayed PSTH
    % ============================================================
    
    % 1. Setup Figure
    figure('Color','w','Position',[200 200 700 500]); 
    
    % --- RIGHT AXIS: Background Patches & Raster Dots ---
    yyaxis right
    ylim([0 1]);
    yticks([]); % Hide the ticks on the right side for a cleaner look
    ylabel(''); % No label needed for raster position
    hold on;
    
    % Define zones
    seq_y_range = [0 0.5];   % Bottom half
    sim_y_range = [0.5 1.0]; % Top half
    
    % A. Plot Background Patches
    % Sequential Zone (Bottom)
    patch([ras_win(1) ras_win(2) ras_win(2) ras_win(1)], ...
          [seq_y_range(1) seq_y_range(1) seq_y_range(2) seq_y_range(2)], ...
          seq_col, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Simultaneous Zone (Top)
    patch([ras_win(1) ras_win(2) ras_win(2) ras_win(1)], ...
          [sim_y_range(1) sim_y_range(1) sim_y_range(2) sim_y_range(2)], ...
          sim_col, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % B. Plot Raster Dots (Vectorized for speed and "Cloud" look)
    % --- Sequential Raster ---
    % --- Sequential Raster ---
    all_seq_x = [];
    all_seq_y = [];
    for i = 1:numel(allSpikes_seq)
        tt = allSpikes_seq{i};
        if isempty(tt), continue; end
        
        % FORCE ROW VECTOR: tt(:)' ensures it is always 1xN
        tt = tt(:)';      
        % Random jitter within the bottom zone
        y_vals = seq_y_range(1) + 0.05 + rand(size(tt)) * (diff(seq_y_range)-0.1);
        
        % Concatenate safely
        all_seq_x = [all_seq_x, tt]; 
        all_seq_y = [all_seq_y, y_vals];
    end
    dot_size = 15;
    dot_color = 'k';
    scatter(all_seq_x, all_seq_y, dot_size, dot_color, 'filled', ...
    'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');

    % --- Simultaneous Raster ---
    all_sim_x = [];
    all_sim_y = [];
    for i = 1:numel(allSpikes_sim)
        tt = allSpikes_sim{i};
        if isempty(tt), continue; end
        % FORCE ROW VECTOR
        tt = tt(:)'; 
        % Random jitter within the top zone
        y_vals = sim_y_range(1) + 0.05 + rand(size(tt)) * (diff(sim_y_range)-0.1);     
        all_sim_x = [all_sim_x, tt];
        all_sim_y = [all_sim_y, y_vals];
    end
    scatter(all_sim_x, all_sim_y, dot_size, dot_color, 'filled', ...
    'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');    
    
    %% Gaussian Smoothing
    sigma_bins = 2; % Start with a small sigma, like 1 or 2, and adjust.

    % 1. Create the Gaussian Kernel
    kernel_size = 2 * ceil(2 * sigma_bins) + 1; % Ensures a reasonable kernel length
    gauss_kernel = gausswin(kernel_size, kernel_size / (2*sigma_bins));
    gauss_kernel = gauss_kernel / sum(gauss_kernel); % Normalize area to 1
    
    % 2. Apply Convolution (Smoothing)
    rate_seq_smooth_gauss = conv(rate_seq, gauss_kernel, 'same');
    rate_sim_smooth_gauss = conv(rate_sim, gauss_kernel, 'same');
    
    % --- LEFT AXIS: PSTH Curves ---
    yyaxis left
    hold on;
    
    % Plot Curves
    % p1 = plot(ctrs, rate_sim,'-','Color', sim_col, 'LineWidth', 2.5, 'DisplayName', 'Simultaneous');
    % p2 = plot(ctrs, rate_seq,'-', 'Color', seq_col, 'LineWidth', 2.5, 'DisplayName', 'Sequential');
    
    p1 = plot(ctrs, rate_sim_smooth_gauss, '-', 'Color', sim_col, 'LineWidth', 2.5, 'DisplayName', 'Simultaneous');
    p2 = plot(ctrs, rate_seq_smooth_gauss, '-', 'Color', seq_col, 'LineWidth', 2.5, 'DisplayName', 'Sequential');
    % Formatting Left Axis
    ylabel('Firing rate (Sp/s)', 'FontWeight', 'bold');
    xlim(ras_win);
    set(gca, 'YColor', 'k'); % Make left axis text black
    
    % Stimulus Line (Zero line)
    xline(0, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % --- FINAL TOUCHES ---
    title(sprintf('Ch %d | %.0f µA | Inter-Stimulus Interval %.0f ms', ...
        target_channel, plot_amp, plot_PTD_ms));
    xlabel('Time (ms)');
    legend([p1, p2], 'Location', 'northeast', 'Box', 'off');
    
    % Ensure the background is transparent so patches show through
    set(gca, 'Color', 'none');
    hold off;
    axis square;
end

