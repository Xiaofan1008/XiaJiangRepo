%% ============================================================
%   SPATIAL ANALYSIS: Single Dataset (Curves & Last-Crossing R50)
%   - Logic: 
%       1. Set Matching: Explicitly pairs Sim/Seq by Electrode IDs.
%       2. Profile Plots: Global, Inner, and Outer Recruitment Curves (100um bins).
%       3. Normalized Plots: Density relative to total shank (16 ch).
%       4. R50 Math: "Last Crossing" logic to handle noisy gaps.
%   - Style: IEEE Publication Ready (Black & White Markers)
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim8';
Electrode_Type = 2; % 2 = 4-Shank
dist_bin_edges = 0 : 50 : 800;  % High-res for R50 Math
plot_bin_edges = 0 : 100 : 800; % Smooth-res for Figures
Seq_PTD_Target = 5; 

save_figures = false;
save_dir     = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank';

%% =================== 1. LOAD DATA ====================
fprintf('Loading Data from: %s\n', data_folder);
[R, S, QC] = load_experiment_data(data_folder);
[Amps, ptd_sim_idx, ptd_seq_idx, uniqueComb] = parse_stim_params(S, Seq_PTD_Target);
nSets = size(uniqueComb, 1);
nCh_Total = 64; 

%% ================= 2. PROCESS SPATIAL DATA =================
SpatialResults = struct();
bin_centers_50  = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;
bin_centers_100 = plot_bin_edges(1:end-1) + diff(plot_bin_edges)/2;

for ss = 1:nSets
    % Identify Stimulation Electrodes for this Set
    stimCh_seq = uniqueComb(ss,:); stimCh_seq = stimCh_seq(stimCh_seq > 0);
    [~, shank_id] = get_electrode_coords(stimCh_seq(1), Electrode_Type);
    
    % Find matching Sim Set (0ms) by Electrode ID
    sim_ss_idx = [];
    for idx = 1:nSets
        stimCh_sim = uniqueComb(idx,:); stimCh_sim = stimCh_sim(stimCh_sim > 0);
        if isequal(sort(stimCh_seq), sort(stimCh_sim))
            sim_ss_idx = idx; break;
        end
    end
    if isempty(sim_ss_idx), continue; end
    
    % Get recording channels on the same shank
    active_ch = [];
    for ch = 1:nCh_Total
        [~, s_id] = get_electrode_coords(ch, Electrode_Type);
        if s_id == shank_id, active_ch = [active_ch; ch]; end
    end
    
    % Distance and Zone Classification
    dist_to_nearest = nan(length(active_ch), 1);
    is_inner = false(length(active_ch), 1);
    stim_y = (mod(stimCh_seq-1, 16)) * 50;
    y_min = min(stim_y); y_max = max(stim_y);
    
    for i = 1:length(active_ch)
        ch_y = (mod(active_ch(i)-1, 16)) * 50;
        dist_to_nearest(i) = min(abs(ch_y - stim_y));
        if ch_y > y_min && ch_y < y_max, is_inner(i) = true; end
    end
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai); if curr_amp == 0, continue; end
        
        [sim_resp, sim_valid] = get_resp_vector(R, sim_ss_idx, ai, ptd_sim_idx, active_ch, QC);
        [seq_resp, seq_valid] = get_resp_vector(R, ss, ai, ptd_seq_idx, active_ch, QC);
        
        SpatialResults.Set(ss).Amp(ai).Val = curr_amp;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Sim = (sum(sim_resp & sim_valid) / sum(sim_valid)) * 100;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Seq = (sum(seq_resp & seq_valid) / sum(seq_valid)) * 100;
        
        % A. Profile Probabilities (100um bins)
        SpatialResults.Set(ss).Amp(ai).Prob_Global_Sim = calc_binned_prob(dist_to_nearest, sim_resp, sim_valid, plot_bin_edges, 'prob');
        SpatialResults.Set(ss).Amp(ai).Prob_Global_Seq = calc_binned_prob(dist_to_nearest, seq_resp, seq_valid, plot_bin_edges, 'prob');
        SpatialResults.Set(ss).Amp(ai).Prob_Inner_Sim = calc_binned_prob(dist_to_nearest(is_inner), sim_resp(is_inner), sim_valid(is_inner), plot_bin_edges, 'prob');
        SpatialResults.Set(ss).Amp(ai).Prob_Inner_Seq = calc_binned_prob(dist_to_nearest(is_inner), seq_resp(is_inner), seq_valid(is_inner), plot_bin_edges, 'prob');
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Sim = calc_binned_prob(dist_to_nearest(~is_inner), sim_resp(~is_inner), sim_valid(~is_inner), plot_bin_edges, 'prob');
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Seq = calc_binned_prob(dist_to_nearest(~is_inner), seq_resp(~is_inner), seq_valid(~is_inner), plot_bin_edges, 'prob');
        
        % B. Profile Normalized Density (Total Shank Fraction)
        SpatialResults.Set(ss).Amp(ai).Norm_Global_Sim = calc_binned_prob(dist_to_nearest, sim_resp, sim_valid, plot_bin_edges, 'norm');
        SpatialResults.Set(ss).Amp(ai).Norm_Global_Seq = calc_binned_prob(dist_to_nearest, seq_resp, seq_valid, plot_bin_edges, 'norm');
        
        % C. R50 Math (Using Last-Crossing Logic)
        prob_50_sim = calc_binned_prob(dist_to_nearest, sim_resp, sim_valid, dist_bin_edges, 'prob');
        prob_50_seq = calc_binned_prob(dist_to_nearest, seq_resp, seq_valid, dist_bin_edges, 'prob');
        SpatialResults.Set(ss).Amp(ai).R50_Sim = interpolate_r50(bin_centers_50, prob_50_sim);
        SpatialResults.Set(ss).Amp(ai).R50_Seq = interpolate_r50(bin_centers_50, prob_50_seq);
    end
end

%% ================= 3. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     SPATIAL SUMMARY: ACTIVE %% AND R50 (RADIUS OF ACTIVATION)       \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-6s | %-12s | %-12s | %-10s | %-10s\n', 'Set', 'Amp', 'Active% Sim', 'Active% Seq', 'R50 Sim', 'R50 Seq');
fprintf('--------------------------------------------------------------------\n');
for ss = 1:nSets
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
        fprintf('Set %02d | %3.1fuA | %5.1f%%       | %5.1f%%       | %5.1fum   | %5.1fum\n', ...
            ss, D.Val, D.TotalPerc_Sim, D.TotalPerc_Seq, D.R50_Sim, D.R50_Seq);
    end
end

%% ================= 4. PLOTTING: RECRUITMENT CURVES (3-PANEL) =================
for ss = 1:nSets
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
        figure('Units', 'centimeters', 'Position', [2, 2, 24, 8], 'Color', 'w', 'Name', 'Probability Profiles');
        t = tiledlayout(1,3, 'TileSpacing', 'compact');
        zones = {'Global', 'Inner', 'Outer'};
        data_sim = {D.Prob_Global_Sim, D.Prob_Inner_Sim, D.Prob_Outer_Sim};
        data_seq = {D.Prob_Global_Seq, D.Prob_Inner_Seq, D.Prob_Outer_Seq};
        for p = 1:3
            nexttile; hold on;
            plot(bin_centers_100, data_sim{p}, '--ok', 'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
            plot(bin_centers_100, data_seq{p}, '-sk', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
            title(zones{p}); ylim([0 1.1]); xlim([0 800]); xticks(0:200:800); grid on;
        end
        xlabel(t, 'Distance to Nearest Stimulation Site (µm)', 'FontSize', 10);
        sgtitle(sprintf('Probability Profiles: Set %d | %.1f µA', ss, D.Val));
    end
end

%% ================= 5. PLOTTING: SHANK-NORMALIZED DENSITY =================
for ss = 1:nSets
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
        figure('Units', 'centimeters', 'Position', [5, 10, 10, 8], 'Color', 'w', 'Name', 'Normalized Density');
        hold on;
        plot(bin_centers_100, D.Norm_Global_Sim, '--ok', 'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
        plot(bin_centers_100, D.Norm_Global_Seq, '-sk', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
        xlabel('Distance (µm)'); ylabel('Fraction of Total Shank (N=16)');
        title(sprintf('Global Density: Set %d | %.1f µA', ss, D.Val));
        xlim([0 800]); ylim([0 0.5]); grid on; axis square;
    end
end

%% ================= 6. PLOT: SUMMARY CURVES =================
for ss = 1:nSets
    figure('Units', 'centimeters', 'Position', [5, 5, 18, 8], 'Color', 'w', 'Name', 'Summary Plots');
    t = tiledlayout(1,2);
    nexttile; hold on; % Active %
    amps = [SpatialResults.Set(ss).Amp.Val];
    plot(amps, [SpatialResults.Set(ss).Amp.TotalPerc_Sim], '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
    plot(amps, [SpatialResults.Set(ss).Amp.TotalPerc_Seq], '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
    xlabel('Amplitude (\muA)'); ylabel('Active Channel %'); title('Active Percentage'); grid on; axis square;
    nexttile; hold on; % R50
    plot(amps, [SpatialResults.Set(ss).Amp.R50_Sim], '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    plot(amps, [SpatialResults.Set(ss).Amp.R50_Seq], '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Amplitude (\muA)'); ylabel('R50 Distance (\mum)'); title('R_{50} (Last Crossing)');
    grid on; axis square;
end

%% ================= 7. SAVE RESULTS =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
out_name = fullfile(save_dir, ['Result_Spatial_SameShank_DX013_' exp_id '.mat']);
save(out_name, 'SpatialResults', 'Amps', 'dist_bin_edges', 'Electrode_Type');
fprintf('\n>>> Spatial Data Saved: %s\n', out_name);

%% ================= HELPER FUNCTIONS =================
function r50 = interpolate_r50(x, y)
    r50 = NaN; valid_mask = ~isnan(y);
    if sum(valid_mask) < 2, return; end
    xv = x(valid_mask); yv = y(valid_mask);
    
    % LAST CROSSING LOGIC: Scan from the outside in for the last point >= 50%
    idx = find(yv >= 0.5, 1, 'last'); 
    
    if isempty(idx)
        r50 = xv(1); % If no point >= 50%, set to minimal recorded distance
    elseif idx == length(yv)
        r50 = xv(end); % If last point is >= 50%, radius extends to end of recording
    else
        % Linear interpolation between the last point >= 50% and the next point below
        x1 = xv(idx);   y1 = yv(idx);
        x2 = xv(idx+1); y2 = yv(idx+1);
        r50 = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1);
    end
end

function metric = calc_binned_prob(dist, resp, valid, edges, mode)
    metric = nan(1, length(edges)-1);
    for b = 1:length(edges)-1
        mask = dist >= edges(b) & dist < edges(b+1) & valid;
        if sum(mask) > 0
            if strcmp(mode, 'prob'), metric(b) = sum(resp(mask)) / sum(mask);
            else, metric(b) = sum(resp(mask)) / 16; end
        elseif strcmp(mode, 'norm'), metric(b) = 0; end
    end
end

function [resp, valid] = get_resp_vector(R, s_idx, a_idx, p_idx, ch_subset, QC)
    resp = zeros(length(ch_subset), 1); valid = true(length(ch_subset), 1);
    try
        chn_list = R.set(s_idx).amp(a_idx).ptd(p_idx).channel;
        for i = 1:length(ch_subset)
            c = ch_subset(i);
            if ~isempty(QC.BadCh) && s_idx <= length(QC.BadCh) && ismember(c, QC.BadCh{s_idx}), valid(i) = false; continue; end
            if c <= length(chn_list) && isfield(chn_list(c), 'is_responsive') && chn_list(c).is_responsive, resp(i) = 1; end
        end
    catch; end
end

function [Amps, ptd_sim_idx, ptd_seq_idx, uniqueComb] = parse_stim_params(S, target_ptd)
    Stim = S.StimParams; simN = S.simultaneous_stim;
    if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end
    amps_all = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
    [Amps,~,~] = unique(trialAmps); Amps(Amps==-1) = 0;
    PTD_all_us = cell2mat(Stim(3:simN:end,6)); PTDs_ms = PTD_all_us / 1000;
    unique_ptds = unique(PTDs_ms);
    ptd_sim_idx = find(abs(unique_ptds - 0) < 0.001);
    ptd_seq_idx = find(abs(unique_ptds - target_ptd) < 0.001);
    stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, S.E_MAP(2:end));
    comb = zeros(nTr, simN);
    for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
    [uniqueComb,~,~] = unique(comb,'rows','stable');
end

function [coords, shank_id] = get_electrode_coords(ch_idx, type)
    if type == 2 % 4-Shank
        if ch_idx <= 16, shank_id = 1; x = 0; 
        elseif ch_idx <= 32, shank_id = 4; x = 600;
        elseif ch_idx <= 48, shank_id = 2; x = 200;
        else, shank_id = 3; x = 400; end
        y = (mod(ch_idx-1, 16)) * 50;
    else, shank_id = 1; x = 0; y = (ch_idx-1) * 50; end
    coords = [x, y];
end

function [R, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*RespondingChannels.mat'); R = load(f(1).name).Responding;
    f_exp = dir('*_exp_datafile_*.mat'); S = load(f_exp(1).name);
    QC.BadCh = []; f_bc = dir('*.BadChannels.mat'); 
    if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
end