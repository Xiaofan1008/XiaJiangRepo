%% ============================================================
%   SINGLE-CHANNEL VERIFICATION: Raster & PSTH (Single Folder)
%   - Purpose: Visualizes raw spikes to verify FWHM Duration metric.
%   - Input: ONE data folder containing BOTH Sim and Seq trials.
%   - Output: 1 Figure per Channel, per Stimulation Set. 
%   - Layout: Top Row = Combined Rasters | Bottom Row = PSTH & FWHM Lines.
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= 1. USER SETTINGS =========================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim2'; 
Electrode_Type = 2;

% --- Visual & Channel Settings ---
channels_to_plot  = [1:64];       % List the channels you want to inspect
amps_to_plot      = [5, 10];   % E.g., one Low, one Mid, one High amplitude
remove_bad_trials = true;         % Set to true to drop QC.BadTrials
Seq_PTD_Target    = 5;            % The Sequential delay (ms)

% --- PSTH & Duration Math Settings ---
plot_win_ms     = [-20 40];    % X-axis limits for the plot
baseline_win_ms = [-50 -10];   % Background noise window
analysis_win_ms = [2 20];      % Search window for response
psth_bin_ms     = 1;          
FS              = 30000;

%% ================= 2. LOAD & PARSE DATA =====================
fprintf('Loading Dataset...\n');
[R, sp, trig, S, QC] = load_experiment_data(data_folder);

Stim = S.StimParams; 
simN = S.simultaneous_stim; 
if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end

amps_all  = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
[Amps,~,ampIdx] = unique(trialAmps); Amps(Amps==-1) = 0;

if simN > 1, PTD_all_us = cell2mat(Stim(3:simN:end,6)); else, PTD_all_us = zeros(nTr,1); end
PTD_all_ms = PTD_all_us / 1000; 
[PTDs_ms,~,ptdIdx] = unique(PTD_all_ms);

stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, S.E_MAP(2:end));
comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[uniqueComb,~,combClass] = unique(comb,'rows','stable'); 

ptd_sim_target = find(abs(PTDs_ms - 0) < 0.001);
ptd_seq_target = find(abs(PTDs_ms - Seq_PTD_Target) < 0.001);

d = Depth_s(Electrode_Type); 
full_edges = (baseline_win_ms(1)-10) : psth_bin_ms : (plot_win_ms(2)+10);
centers = full_edges(1:end-1) + psth_bin_ms/2;

% [MODIFIED] Calculate total number of stimulation Sets (Orders)
nSets = size(uniqueComb, 1);

%% ================= 3. PROCESS & PLOT PER CHANNEL ============
for c = 1:length(channels_to_plot)
    ch = channels_to_plot(c);
    ch_idx = find(d == ch, 1);
    if isempty(ch_idx), continue; end
    
    % [MODIFIED] Added an outer loop to generate plots for EVERY stimulation order
    for ss = 1:nSets
        
        nAmps = length(amps_to_plot);
        figTitle = sprintf('Ch %02d | Set %d | Raster & PSTH Verification', ch, ss);
        figure('Color','w', 'Position', [100, 100, 300*nAmps, 500], 'Name', figTitle); % Adjusted height for 2 rows
        
        for a = 1:nAmps
            target_amp = amps_to_plot(a);
            [~, ai] = min(abs(Amps - target_amp)); 
            
            tr_sim = find(ampIdx == ai & combClass == ss & ptdIdx == ptd_sim_target);
            tr_seq = find(ampIdx == ai & combClass == ss & ptdIdx == ptd_seq_target);
            
            % Skip plotting this column if there is literally no data for this Set/Amp combination
            if isempty(tr_sim) && isempty(tr_seq), continue; end
            
            if remove_bad_trials && ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
                tr_sim = setdiff(tr_sim, QC.BadTrials{ch_idx}); 
                tr_seq = setdiff(tr_seq, QC.BadTrials{ch_idx}); 
            end
            
            spikes_sim = extract_trials(tr_sim, trig, sp{ch_idx}, full_edges, FS);
            spikes_seq = extract_trials(tr_seq, trig, sp{ch_idx}, full_edges, FS);
            
            % --- ROW 1: COMBINED RASTERS ---
            % [MODIFIED] Changed layout to 2 rows and plotted both rasters together
            subplot(2, nAmps, a); hold on;
            
            % Plot Sim (Bottom)
            plot_raster(spikes_sim, tr_sim, 'k', 0);
            nSimTrials = max(1, length(tr_sim)); % Use max(1) to avoid math collapsing
            
            % Plot Seq (Top, offset by nSimTrials)
            plot_raster(spikes_seq, tr_seq, [0.5 0.5 0.5], nSimTrials);
            nSeqTrials = max(1, length(tr_seq));
            
            % Draw a subtle horizontal divider line
            yline(nSimTrials + 0.5, ':', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
            
            title(sprintf('%.1f \\muA', Amps(ai)), 'FontWeight', 'bold', 'FontSize', 12);
            if a == 1, ylabel('Trials (Sim \downarrow | Seq \uparrow)', 'FontWeight', 'bold'); end
            xlim(plot_win_ms); ylim([0.5, nSimTrials + nSeqTrials + 0.5]);
            set(gca, 'TickDir', 'out'); box off;
            
            % --- ROW 2: OVERLAID PSTH & METRIC ---
            subplot(2, nAmps, a + nAmps); hold on;
            
            % [MODIFIED] Function now returns the precise coordinates for the FWHM lines
            [rate_sim, mask_sim, lines_sim] = calc_psth_fwhm(spikes_sim, length(tr_sim), full_edges, centers, baseline_win_ms, analysis_win_ms);
            [rate_seq, mask_seq, lines_seq] = calc_psth_fwhm(spikes_seq, length(tr_seq), full_edges, centers, baseline_win_ms, analysis_win_ms);
            
            % [MODIFIED] Swapped Line Styles (Sim=Dashed, Seq=Solid)
            plot(centers, rate_sim, '-k', 'LineWidth', 2);
            plot(centers, rate_seq, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
            
            % [MODIFIED] Draw sleek horizontal/vertical lines instead of solid block shading
            for i = 1:length(lines_sim)
                plot(lines_sim(i).x, lines_sim(i).y, '-k', 'LineWidth', 1.5);
                plot([lines_sim(i).x(1), lines_sim(i).x(1)], [0, lines_sim(i).y(1)], ':k');
                plot([lines_sim(i).x(2), lines_sim(i).x(2)], [0, lines_sim(i).y(2)], ':k');
            end
            for i = 1:length(lines_seq)
                plot(lines_seq(i).x, lines_seq(i).y, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
                plot([lines_seq(i).x(1), lines_seq(i).x(1)], [0, lines_seq(i).y(1)], ':', 'Color', [0.5 0.5 0.5]);
                plot([lines_seq(i).x(2), lines_seq(i).x(2)], [0, lines_seq(i).y(2)], ':', 'Color', [0.5 0.5 0.5]);
            end
            
            xlim(plot_win_ms); xlabel('Time (ms)'); 
            if a == 1, ylabel('Firing Rate (Hz)', 'FontWeight', 'bold'); end
            set(gca, 'TickDir', 'out'); box off;
            
            max_y = max([rate_sim, rate_seq, 10]); ylim([0 max_y*1.15]);
            
            dur_sim = sum(mask_sim) * psth_bin_ms; dur_seq = sum(mask_seq) * psth_bin_ms;
            text(plot_win_ms(1)+2, max_y*1.05, sprintf('Sim: %d ms', dur_sim), 'Color', 'k', 'FontSize', 9);
            text(plot_win_ms(1)+2, max_y*0.90, sprintf('Seq: %d ms', dur_seq), 'Color', [0.4 0.4 0.4], 'FontSize', 9);
        end
    end
end

%% ==================== HELPER FUNCTIONS =========================
function spikes = extract_trials(tr_ids, trig, sp_data, edges, FS)
    spikes = []; win_pad = [edges(1), edges(end)];
    for k = 1:length(tr_ids)
        t0 = trig(tr_ids(k))/FS*1000; tt = sp_data(:,1) - t0;
        mask = tt >= win_pad(1) & tt <= win_pad(2);
        spikes = [spikes; [tt(mask), repmat(k, sum(mask), 1)]]; %#ok<AGROW>
    end
end

% [MODIFIED] Added y_offset to stack rasters on top of each other
function plot_raster(spikes, tr_ids, clr, y_offset)
    if isempty(spikes), return; end
    line([spikes(:,1)'; spikes(:,1)'], [(spikes(:,2)' + y_offset)-0.4; (spikes(:,2)' + y_offset)+0.4], 'Color', clr);
    yline(0, 'k', 'LineWidth', 1);
end

% [MODIFIED] Returns exact x/y coordinates for the width lines instead of a solid mask block
function [smooth_rate, valid_mask, fwhm_lines] = calc_psth_fwhm(spikes, nTr, edges, centers, base_win, anal_win)
    fwhm_lines = []; % Initialize struct for drawing lines
    if isempty(spikes) || nTr == 0
        smooth_rate = zeros(size(centers)); valid_mask = false(size(centers)); return;
    end
    
    bin_size = edges(2) - edges(1); counts = histcounts(spikes(:,1), edges);
    rate_hz = (counts / nTr) * (1000 / bin_size); smooth_rate = smoothdata(rate_hz, 'gaussian', 5);
    
    base_mask = centers >= base_win(1) & centers <= base_win(2);
    mu = mean(smooth_rate(base_mask)); sigma = std(smooth_rate(base_mask)); if sigma == 0, sigma = 1; end 
    threshold = max(mu + 3 * sigma, 10);
    
    anal_mask = centers >= anal_win(1) & centers <= anal_win(2); 
    anal_rate = smooth_rate; anal_rate(~anal_mask) = 0; 
    
    [pks, locs] = findpeaks(anal_rate, 'MinPeakHeight', threshold);
    valid_mask = false(size(smooth_rate));
    
    for i = 1:length(pks)
        half_height = pks(i)/2;
        
        % Trace backwards and forwards to find exact 50% crossings
        left_idx = locs(i);
        while left_idx > 1 && anal_rate(left_idx-1) >= half_height, left_idx = left_idx - 1; end
        
        right_idx = locs(i);
        while right_idx < length(anal_rate) && anal_rate(right_idx+1) >= half_height, right_idx = right_idx + 1; end
        
        valid_mask(left_idx:right_idx) = true;
        
        % Save the horizontal line coordinates for plotting
        fwhm_lines(end+1).x = [centers(left_idx), centers(right_idx)]; %#ok<AGROW>
        fwhm_lines(end).y   = [half_height, half_height];
    end
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
    f_bt = dir('*.BadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end