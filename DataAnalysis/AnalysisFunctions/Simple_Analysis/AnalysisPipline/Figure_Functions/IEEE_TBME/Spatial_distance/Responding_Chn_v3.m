%% ============================================================
%   SPATIAL ANALYSIS: Single Dataset (Curves & R80 Cumulative)
%   - Logic: 
%       1. Probe Detection: 16-ch (4-shank) vs 32-ch (1-shank).
%       2. Neighbor Filter: Removes isolated active channels (noise).
%       3. R80 Reach: Distance containing 80% of total cleaned activity.
%       4. Summary Plots: Global potency and reach vs amplitude.
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim8';
Seq_PTD_Target = 5; 
save_figures = false;
save_dir     = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Resp_SameShank_v2';

% High-Resolution 50um Binning
dist_bin_edges = 0 : 50 : 800;  
plot_bin_edges = 0 : 50 : 800; 

%% =================== 1. LOAD DATA ====================
fprintf('Loading Data from: %s\n', data_folder);
[R, S, QC] = load_experiment_data(data_folder);
[Amps, ptd_sim_idx, ptd_seq_idx, uniqueComb] = parse_stim_params(S, Seq_PTD_Target);
nSets = size(uniqueComb, 1);

% --- DETECT PROBE TYPE ---
nCh_Raw = length(fieldnames(R.set(1).amp(1).ptd(1).channel)); 
if nCh_Raw <= 32
    active_N = 32; Electrode_Type = 1; nCh_Total = 32;
else
    active_N = 16; Electrode_Type = 2; nCh_Total = 64;
end
fprintf('Detected Probe: %d channels. Denominator: %d\n', nCh_Raw, active_N);

%% ================= 2. PROCESS SPATIAL DATA =================
SpatialResults = struct();
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;

for ss = 1:nSets
    stimCh_seq = uniqueComb(ss,:); stimCh_seq = stimCh_seq(stimCh_seq > 0);
    [~, shank_id] = get_electrode_coords(stimCh_seq(1), Electrode_Type);
    
    sim_ss_idx = [];
    for idx = 1:nSets
        stimCh_sim = uniqueComb(idx,:); stimCh_sim = stimCh_sim(stimCh_sim > 0);
        if isequal(sort(stimCh_seq), sort(stimCh_sim)), sim_ss_idx = idx; break; end
    end
    if isempty(sim_ss_idx), continue; end
    
    active_ch = [];
    for ch = 1:nCh_Total
        [~, s_id] = get_electrode_coords(ch, Electrode_Type);
        if s_id == shank_id, active_ch = [active_ch; ch]; end
    end
    
    y_coords = (mod(active_ch-1, 16)) * 50; 
    stim_y = (mod(stimCh_seq-1, 16)) * 50;
    y_min = min(stim_y); y_max = max(stim_y);
    
    dist_to_nearest = zeros(length(active_ch), 1);
    is_inner = false(length(active_ch), 1);
    for i = 1:length(active_ch)
        dist_to_nearest(i) = min(abs(y_coords(i) - stim_y));
        if y_coords(i) > y_min && y_coords(i) < y_max, is_inner(i) = true; end
    end
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai); if curr_amp == 0, continue; end
        [sim_resp, sim_valid] = get_resp_vector(R, sim_ss_idx, ai, ptd_sim_idx, active_ch, QC);
        [seq_resp, seq_valid] = get_resp_vector(R, ss, ai, ptd_seq_idx, active_ch, QC);
        
        % --- NEIGHBOR FILTER ---
        clean_sim = filter_isolated_channels(sim_resp, y_coords, 100);
        clean_seq = filter_isolated_channels(seq_resp, y_coords, 100);
        
        SpatialResults.Set(ss).Amp(ai).Val = curr_amp;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Sim = (sum(sim_resp & sim_valid) / sum(sim_valid)) * 100;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Seq = (sum(seq_resp & seq_valid) / sum(seq_valid)) * 100;
        
        % A. Profile Probabilities (Raw data for curves)
        SpatialResults.Set(ss).Amp(ai).Prob_Global_Sim = calc_binned_prob(dist_to_nearest, sim_resp, sim_valid, plot_bin_edges, 'prob', active_N);
        SpatialResults.Set(ss).Amp(ai).Prob_Global_Seq = calc_binned_prob(dist_to_nearest, seq_resp, seq_valid, plot_bin_edges, 'prob', active_N);
        SpatialResults.Set(ss).Amp(ai).Prob_Inner_Sim  = calc_binned_prob(dist_to_nearest(is_inner), sim_resp(is_inner), sim_valid(is_inner), plot_bin_edges, 'prob', active_N);
        SpatialResults.Set(ss).Amp(ai).Prob_Inner_Seq  = calc_binned_prob(dist_to_nearest(is_inner), seq_resp(is_inner), seq_valid(is_inner), plot_bin_edges, 'prob', active_N);
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Sim  = calc_binned_prob(dist_to_nearest(~is_inner), sim_resp(~is_inner), sim_valid(~is_inner), plot_bin_edges, 'prob', active_N);
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Seq  = calc_binned_prob(dist_to_nearest(~is_inner), seq_resp(~is_inner), seq_valid(~is_inner), plot_bin_edges, 'prob', active_N);

        % B. R80 Reach Calculation (Cleaned data)
        SpatialResults.Set(ss).Amp(ai).R80_Sim = calculate_R_cumulative(dist_to_nearest, clean_sim, 0.80);
        SpatialResults.Set(ss).Amp(ai).R80_Seq = calculate_R_cumulative(dist_to_nearest, clean_seq, 0.80);
        
        % Store density
        SpatialResults.Set(ss).Amp(ai).Norm_Global_Sim = calc_binned_prob(dist_to_nearest, sim_resp, sim_valid, dist_bin_edges, 'norm', active_N);
        SpatialResults.Set(ss).Amp(ai).Norm_Global_Seq = calc_binned_prob(dist_to_nearest, seq_resp, seq_valid, dist_bin_edges, 'norm', active_N);
    end
end

%% ================= 3. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     SPATIAL SUMMARY: ACTIVE %% AND R80 (80%% CUMULATIVE REACH)      \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-6s | %-12s | %-12s | %-10s | %-10s\n', 'Set', 'Amp', 'Active% Sim', 'Active% Seq', 'R80 Sim', 'R80 Seq');
fprintf('--------------------------------------------------------------------\n');
for ss = 1:nSets
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
        fprintf('Set %02d | %3.1fuA | %5.1f%%       | %5.1f%%       | %5.1fum   | %5.1fum\n', ...
            ss, D.Val, D.TotalPerc_Sim, D.TotalPerc_Seq, D.R80_Sim, D.R80_Seq);
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
            plot(bin_centers, data_sim{p}, '--ok', 'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
            plot(bin_centers, data_seq{p}, '-sk', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
            title(zones{p}); ylim([0 1.1]); xlim([0 800]); grid on;
            if p == 3
                legend('Location', 'northeast', 'Box', 'off'); 
            end
        end
        xlabel(t, 'Distance to Nearest Stimulation Site (µm)', 'FontSize', 10);
        sgtitle(sprintf('Set %d | %.1f µA (HD 50µm Bins)', ss, D.Val));
    end
end

%% ================= 5. PLOTTING: SHANK-NORMALIZED DENSITY =================
for ss = 1:nSets
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
        figure('Units', 'centimeters', 'Position', [5, 10, 10, 8], 'Color', 'w', 'Name', 'Density');
        hold on;
        plot(bin_centers, D.Norm_Global_Sim, '--ok', 'LineWidth', 1.2, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
        plot(bin_centers, D.Norm_Global_Seq, '-sk', 'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
        xlabel('Distance (µm)'); ylabel(sprintf('Fraction of Shank (N=%d)', active_N));
        title(sprintf('Total Density: Set %d | %.1f µA', ss, D.Val));
        xlim([0 800]); ylim([0 0.2]);  axis square; grid on;
        % Add this right before the end of the Section 5 loops
        legend('Simultaneous', 'Sequential', 'Location', 'northeast', 'Box', 'off');
    end
end

%% ================= 6. SUMMARY PLOTS (SET TRENDS) =================
figure('Units', 'centimeters', 'Position', [5, 5, 20, 10], 'Color', 'w', 'Name', 'Dataset Summary');
t = tiledlayout(1,2);
nexttile; hold on; % Potency
for ss = 1:nSets
    vals = [SpatialResults.Set(ss).Amp.Val];
    plot(vals, [SpatialResults.Set(ss).Amp.TotalPerc_Sim], '--ok', 'MarkerFaceColor', 'w');
    plot(vals, [SpatialResults.Set(ss).Amp.TotalPerc_Seq], '-sk', 'MarkerFaceColor', 'k');
end
xlabel('Amplitude (\muA)'); ylabel('Active Channel %'); title('Active Percentage'); axis square;
legend('Sim', 'Seq', 'Location', 'northwest', 'Box', 'off');

nexttile; hold on; % Reach
for ss = 1:nSets
    vals = [SpatialResults.Set(ss).Amp.Val];
    plot(vals, [SpatialResults.Set(ss).Amp.R80_Sim], '--ok', 'MarkerFaceColor', 'w');
    plot(vals, [SpatialResults.Set(ss).Amp.R80_Seq], '-sk', 'MarkerFaceColor', 'k');
end
xlabel('Amplitude (\muA)'); ylabel('R_{80} Reach (\mum)'); title('Effective Reach (80%)'); axis square;
legend('Sim', 'Seq', 'Location', 'northwest', 'Box', 'off');
%% ================= 7. SAVE RESULTS =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
out_name = fullfile(save_dir, ['Result_Spatial_R80_DX013_' exp_id '.mat']);
save(out_name, 'SpatialResults', 'Amps', 'dist_bin_edges', 'active_N');
fprintf('\n>>> Spatial Data Saved: %s\n', out_name);

%% ================= HELPER FUNCTIONS =================
function r_cum = calculate_R_cumulative(dists, resp, threshold)
    r_cum = 0; active_idx = find(resp > 0);
    if isempty(active_idx), return; end
    [sorted_dist, ~] = sort(dists(active_idx));
    n_total = length(active_idx);
    cutoff_idx = ceil(threshold * n_total);
    r_cum = sorted_dist(cutoff_idx);
end

function cleaned = filter_isolated_channels(resp, y, radius)
    cleaned = resp; active_idx = find(resp > 0);
    for i = 1:length(active_idx)
        idx = active_idx(i); other_active = setdiff(active_idx, idx);
        dist_to_others = abs(y(idx) - y(other_active));
        if isempty(dist_to_others) || min(dist_to_others) > radius
            cleaned(idx) = 0; 
        end
    end
end

function metric = calc_binned_prob(dist, resp, valid, edges, mode, denom)
    metric = nan(1, length(edges)-1);
    for b = 1:length(edges)-1
        mask = dist >= edges(b) & dist < edges(b+1) & valid;
        if sum(mask) > 0
            if strcmp(mode, 'prob'), metric(b) = sum(resp(mask)) / sum(mask);
            else, metric(b) = sum(resp(mask)) / denom; end
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