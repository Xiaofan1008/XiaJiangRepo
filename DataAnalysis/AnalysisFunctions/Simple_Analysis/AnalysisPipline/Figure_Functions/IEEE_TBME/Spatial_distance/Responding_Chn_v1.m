%% ============================================================
%   SPATIAL ANALYSIS: Single Dataset (Curves & R50)
%   - Logic: 
%       1. Same-Shank Analysis (Active Shank only).
%       2. Profile Plots: Global, Inner, and Outer Recruitment Curves.
%       3. Summary Plots: Active Chn % and R50 vs. Amplitude.
%   - Style: IEEE Publication Ready (Black & White Markers)
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Seq_Sim5';

Electrode_Type = 2; % 2 = 4-Shank
dist_bin_edges = 0 : 50 : 800; 
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
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;

for ss = 1:nSets
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh > 0);
    [~, shank_id] = get_electrode_coords(stimCh(1), Electrode_Type);
    
    % Identify recording channels on the same shank as stimulation
    active_ch = [];
    for ch = 1:nCh_Total
        [~, s_id] = get_electrode_coords(ch, Electrode_Type);
        if s_id == shank_id, active_ch = [active_ch; ch]; end
    end
    
    % Identify Inner vs. Outer spatial zones
    dist_to_nearest = nan(length(active_ch), 1);
    is_inner = false(length(active_ch), 1);
    stim_y = (mod(stimCh-1, 16)) * 50;
    y_min = min(stim_y); y_max = max(stim_y);
    
    for i = 1:length(active_ch)
        ch_y = (mod(active_ch(i)-1, 16)) * 50;
        dist_to_nearest(i) = min(abs(ch_y - stim_y));
        if ch_y > y_min && ch_y < y_max, is_inner(i) = true; end
    end
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        if curr_amp == 0, continue; end
        
        [sim_resp, sim_valid] = get_resp_vector(R, ss, ai, ptd_sim_idx, active_ch, QC);
        [seq_resp, seq_valid] = get_resp_vector(R, ss, ai, ptd_seq_idx, active_ch, QC);
        
        SpatialResults.Set(ss).Amp(ai).Val = curr_amp;
        
        % A. Active % calculation (Relative to total valid channels on the shank)
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Sim = (sum(sim_resp & sim_valid) / sum(sim_valid)) * 100;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Seq = (sum(seq_resp & seq_valid) / sum(seq_valid)) * 100;
        
        % B. Recruitment Curve Probabilities (0-1)
        SpatialResults.Set(ss).Amp(ai).Prob_Global_Sim = calc_binned_prob(dist_to_nearest, sim_resp, sim_valid, dist_bin_edges);
        SpatialResults.Set(ss).Amp(ai).Prob_Global_Seq = calc_binned_prob(dist_to_nearest, seq_resp, seq_valid, dist_bin_edges);
        
        SpatialResults.Set(ss).Amp(ai).Prob_Inner_Sim = calc_binned_prob(dist_to_nearest(is_inner), sim_resp(is_inner), sim_valid(is_inner), dist_bin_edges);
        SpatialResults.Set(ss).Amp(ai).Prob_Inner_Seq = calc_binned_prob(dist_to_nearest(is_inner), seq_resp(is_inner), seq_valid(is_inner), dist_bin_edges);
        
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Sim = calc_binned_prob(dist_to_nearest(~is_inner), sim_resp(~is_inner), sim_valid(~is_inner), dist_bin_edges);
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Seq = calc_binned_prob(dist_to_nearest(~is_inner), seq_resp(~is_inner), seq_valid(~is_inner), dist_bin_edges);
        
        % C. Effective Reach (R50)
        SpatialResults.Set(ss).Amp(ai).R50_Sim = interpolate_r50(bin_centers, SpatialResults.Set(ss).Amp(ai).Prob_Global_Sim);
        SpatialResults.Set(ss).Amp(ai).R50_Seq = interpolate_r50(bin_centers, SpatialResults.Set(ss).Amp(ai).Prob_Global_Seq);
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
% These replace the bar plots to show spatial decay profiles clearly
for ss = 1:nSets
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
        
        figure('Units', 'centimeters', 'Position', [2, 2, 24, 8], 'Color', 'w');
        t = tiledlayout(1,3, 'TileSpacing', 'compact');
        
        % Plotting Style: Simultaneous (Dashed/White Circle) vs Sequential (Solid/Black Square)
        zones = {'Global', 'Inner', 'Outer'};
        data_sim = {D.Prob_Global_Sim, D.Prob_Inner_Sim, D.Prob_Outer_Sim};
        data_seq = {D.Prob_Global_Seq, D.Prob_Inner_Seq, D.Prob_Outer_Seq};
        
        for p = 1:3
            nexttile; hold on;
            plot(bin_centers, data_sim{p}, '--ok', 'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
            plot(bin_centers, data_seq{p}, '-sk', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
            title(zones{p}); ylim([0 1.1]); grid on;
            if p == 1, ylabel('Recruitment Probability'); end
        end
        
        xlabel(t, 'Distance to Nearest Stimulation Site (µm)', 'FontName', 'Arial', 'FontSize', 10);
        sgtitle(sprintf('Spatial Profiles: Set %d | %.1f µA', ss, D.Val));
    end
end

%% ================= 5. PLOT: SUMMARY CURVES =================
for ss = 1:nSets
    figure('Units', 'centimeters', 'Position', [5, 5, 18, 8], 'Color', 'w');
    t = tiledlayout(1,2);
    
    % Left Panel: Total Active Channel Percentage
    nexttile; hold on;
    amps = [SpatialResults.Set(ss).Amp.Val];
    p_sim = [SpatialResults.Set(ss).Amp.TotalPerc_Sim];
    p_seq = [SpatialResults.Set(ss).Amp.TotalPerc_Seq];
    plot(amps, p_sim, '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
    plot(amps, p_seq, '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
    xlabel('Amplitude (\muA)'); ylabel('Active Channel %'); title('Active Percentage');
    legend('Location', 'northwest', 'Box', 'off'); axis square;
    
    % Right Panel: R50 Distance
    nexttile; hold on;
    r_sim = [SpatialResults.Set(ss).Amp.R50_Sim];
    r_seq = [SpatialResults.Set(ss).Amp.R50_Seq];
    plot(amps, r_sim, '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    plot(amps, r_seq, '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
    xlabel('Amplitude (\muA)'); ylabel('R50 Distance (\mum)'); title('R_50');
    axis square;
end

%% ================= 6. SAVE RESULTS =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
out_name = fullfile(save_dir, ['Result_Spatial_SameShank_R50_DX015_' exp_id '.mat']);
save(out_name, 'SpatialResults', 'Amps', 'dist_bin_edges', 'Electrode_Type');
fprintf('\n>>> Spatial Data Saved: %s\n', out_name);

%% ================= HELPER FUNCTIONS =================
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

function prob = calc_binned_prob(dist, resp, valid, edges)
    prob = nan(1, length(edges)-1);
    for b = 1:length(edges)-1
        mask = dist >= edges(b) & dist < edges(b+1) & valid;
        if sum(mask) > 0, prob(b) = sum(resp(mask)) / sum(mask); end
    end
end

function r50 = interpolate_r50(x, y)
    r50 = NaN; valid_mask = ~isnan(y);
    if sum(valid_mask) < 2, return; end
    xv = x(valid_mask); yv = y(valid_mask);
    idx = find(yv < 0.5, 1, 'first');
    if isempty(idx) && all(yv >= 0.5), r50 = max(xv); return; end
    if idx == 1, r50 = min(xv); return; end
    if ~isempty(idx)
        x1 = xv(idx-1); y1 = yv(idx-1); x2 = xv(idx); y2 = yv(idx);
        r50 = x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1);
    end
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