%% ============================================================
%   SPATIAL ANALYSIS: Separate Datasets (Outer Leakage & D90 Containment)
%   - Logic: 
%       1. Probe Detection: 16-ch (4-shank) vs 32-ch (1-shank). Isolates to same shank.
%       2. "Big Point" Distance: Target zone = 0um. Outer channels measured from boundary edge.
%       3. Neighbor Filter: Removes isolated active channels (noise) with 100um tolerance.
%       4. D90 Spread: Distance required to contain 90% of ALL cleaned activity.
%       5. Outer Density: Fraction of TOTAL response that leaked to specific outer distances.
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
% --- MODIFIED: Specify separate folders for Sim and Seq data ---
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Sim'; % Update this path
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Seq'; 

Seq_PTD_Target = 5.5; 
% save_figures = false;
save_dir     = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Outer_Spatial_Percentage';

% Bins for Outer Distance (Starts at 0 boundary, goes up to 800um away)
dist_bin_edges = 0 : 50 : 1600;  
plot_bin_edges = 0 : 50 : 1600; 

%% =================== 1. LOAD DATA ====================
fprintf('Loading Sim Data from: %s\n', folder_sim);
[R_sim, S_sim, QC_sim] = load_experiment_data(folder_sim);
[Amps_sim, ptd_sim_idx, ~, uniqueComb_sim] = parse_stim_params(S_sim, 0); % Target 0 for Sim

fprintf('Loading Seq Data from: %s\n', folder_seq);
[R_seq, S_seq, QC_seq] = load_experiment_data(folder_seq);
[Amps_seq, ~, ptd_seq_idx, uniqueComb_seq] = parse_stim_params(S_seq, Seq_PTD_Target);

% --- MODIFIED: Find common amplitudes tested in both files ---
Amps = intersect(Amps_sim, Amps_seq);
nSets_seq = size(uniqueComb_seq, 1);

% --- DETECT PROBE TYPE ---
nCh_Raw = length(R_seq.set(1).amp(1).ptd(1).channel); 
if nCh_Raw <= 32
    active_N = 32; Electrode_Type = 1; nCh_Total = 32;
else
    active_N = 16; Electrode_Type = 2; nCh_Total = 64;
end
fprintf('Detected Probe: %d channels. Denominator: %d\n', nCh_Raw, active_N);

%% ================= 2. PROCESS SPATIAL DATA =================
SpatialResults = struct();
bin_centers = dist_bin_edges(1:end-1) + diff(dist_bin_edges)/2;

for ss = 1:nSets_seq
    % Find Seq Set Indices
    stimCh_seq = uniqueComb_seq(ss,:); stimCh_seq = stimCh_seq(stimCh_seq > 0);
    [~, shank_id] = get_electrode_coords(stimCh_seq(1), Electrode_Type);
    
    % --- MODIFIED: Search for matching Sim set in the SIM uniqueComb ---
    sim_ss_idx = [];
    for idx = 1:size(uniqueComb_sim, 1)
        stimCh_sim = uniqueComb_sim(idx,:); stimCh_sim = stimCh_sim(stimCh_sim > 0);
        if isequal(sort(stimCh_seq), sort(stimCh_sim)), sim_ss_idx = idx; break; end
    end
    if isempty(sim_ss_idx), continue; end % Skip if this combination wasn't tested in Sim
    
    % Isolate recording channels to the SAME SHANK
    active_ch = [];
    for ch = 1:nCh_Total
        [~, s_id] = get_electrode_coords(ch, Electrode_Type);
        if s_id == shank_id, active_ch = [active_ch; ch]; end
    end
    
    % Get physical Y coordinates
    y_coords = (mod(active_ch-1, 16)) * 50; 
    stim_y = (mod(stimCh_seq-1, 16)) * 50;
    y_min = min(stim_y); y_max = max(stim_y);
    
    % "BIG POINT" DISTANCE CALCULATION
    dist_to_boundary = zeros(length(active_ch), 1);
    is_outer = false(length(active_ch), 1);
    
    for i = 1:length(active_ch)
        if y_coords(i) < y_min
            dist_to_boundary(i) = y_min - y_coords(i); 
            is_outer(i) = true;
        elseif y_coords(i) > y_max
            dist_to_boundary(i) = y_coords(i) - y_max; 
            is_outer(i) = true;
        else
            dist_to_boundary(i) = 0; 
            is_outer(i) = false;
        end
    end
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai); if curr_amp == 0, continue; end
        
        % --- MODIFIED: Find exact amplitude index for each file ---
        ai_sim = find(Amps_sim == curr_amp);
        ai_seq = find(Amps_seq == curr_amp);
        
        % --- MODIFIED: Pull from separated R and QC structs ---
        [sim_resp, sim_valid] = get_resp_vector(R_sim, sim_ss_idx, ai_sim, ptd_sim_idx, active_ch, QC_sim);
        [seq_resp, seq_valid] = get_resp_vector(R_seq, ss, ai_seq, ptd_seq_idx, active_ch, QC_seq);
        
        % --- NEIGHBOR FILTER (100um Tolerance) ---
        clean_sim = filter_isolated_channels(sim_resp, y_coords, 100);
        clean_seq = filter_isolated_channels(seq_resp, y_coords, 100);
        
        SpatialResults.Set(ss).Amp(ai).Val = curr_amp;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Sim = (sum(clean_sim & sim_valid) / sum(sim_valid)) * 100;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Seq = (sum(clean_seq & seq_valid) / sum(seq_valid)) * 100;
        
        % A. PANEL A DATA: Probability of Outer Leakage 
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Sim = calc_binned_prob(dist_to_boundary(is_outer), clean_sim(is_outer), sim_valid(is_outer), plot_bin_edges, 'prob', 0);
        SpatialResults.Set(ss).Amp(ai).Prob_Outer_Seq = calc_binned_prob(dist_to_boundary(is_outer), clean_seq(is_outer), seq_valid(is_outer), plot_bin_edges, 'prob', 0);
        
        % B. PANEL B DATA: Density / Spatial Footprint 
        total_active_sim = sum(clean_sim);
        total_active_seq = sum(clean_seq);
        SpatialResults.Set(ss).Amp(ai).Density_Outer_Sim = calc_binned_prob(dist_to_boundary(is_outer), clean_sim(is_outer), sim_valid(is_outer), plot_bin_edges, 'density', total_active_sim);
        SpatialResults.Set(ss).Amp(ai).Density_Outer_Seq = calc_binned_prob(dist_to_boundary(is_outer), clean_seq(is_outer), seq_valid(is_outer), plot_bin_edges, 'density', total_active_seq);
        
        % C. PANEL C DATA: D90 Total Containment Spread
        SpatialResults.Set(ss).Amp(ai).D90_Sim = calculate_R_cumulative(dist_to_boundary, clean_sim, 0.90);
        SpatialResults.Set(ss).Amp(ai).D90_Seq = calculate_R_cumulative(dist_to_boundary, clean_seq, 0.90);
    end
end

%% ================= 3. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     SPATIAL SUMMARY: OUTER LEAKAGE & D90 CONTAINMENT EXTENT        \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-6s | %-12s | %-12s | %-10s | %-10s\n', 'Set', 'Amp', 'Active% Sim', 'Active% Seq', 'D90 Sim', 'D90 Seq');
fprintf('--------------------------------------------------------------------\n');
for ss = 1:nSets_seq
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
        fprintf('Set %02d | %3.1fuA | %5.1f%%       | %5.1f%%       | %5.1fum   | %5.1fum\n', ...
            ss, D.Val, D.TotalPerc_Sim, D.TotalPerc_Seq, D.D90_Sim, D.D90_Seq);
    end
end

%% ================= 4. PLOTTING: 3-PANEL PER AMPLITUDE =================
for ss = 1:nSets_seq
    for ai = 1:length(Amps)
        D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val < 2, continue; end
        
        figure('Units', 'centimeters', 'Position', [2, 2, 26, 8], 'Color', 'w', 'Name', sprintf('Spatial Profile - Set %d - %.1fuA', ss, D.Val));
        t = tiledlayout(1,3, 'TileSpacing', 'compact');
        
        % Panel 1: Probability Decay (Lines)
        nexttile; hold on;
        plot(bin_centers, D.Prob_Outer_Sim, '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'DisplayName', 'Sim');
        plot(bin_centers, D.Prob_Outer_Seq, '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', 'DisplayName', 'Seq');
        title('A. Outer Leakage Probability'); 
        xlabel('Distance from Boundary (\mum)'); ylabel('Probability of Activation');
        ylim([0 1.1]); xlim([0 800]); grid off; box off;
        legend('Location', 'northeast', 'Box', 'off'); 

        % Panel 2: Density Footprint (Grouped Bars)
        nexttile; hold on;
        bar_data = [D.Density_Outer_Sim(:), D.Density_Outer_Seq(:)];
        b = bar(bin_centers, bar_data, 'grouped');
        b(1).FaceColor = 'w'; b(1).EdgeColor = 'k'; b(1).LineWidth = 1.2;
        b(2).FaceColor = 'k'; b(2).EdgeColor = 'k'; b(2).LineWidth = 1.2;
        title('B. Spatial Density'); 
        xlabel('Distance from Boundary (\mum)'); ylabel('Fraction of Total Response');
        ylim([0 0.3]); xlim([0 800]); grid off; box off;
        
        % Panel 3: D90 Placeholder text
        nexttile; hold on;
        text(0.5, 0.8, sprintf('D_{90} Containment (Sim):\n%.1f \\mum', D.D90_Sim), 'HorizontalAlignment', 'center', 'FontSize', 12);
        text(0.5, 0.2, sprintf('D_{90} Containment (Seq):\n%.1f \\mum', D.D90_Seq), 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
        title('C. Max Spread (D_{90})'); 
        axis off;
        
        sgtitle(sprintf('Set %d | %.1f \\muA', ss, D.Val), 'FontWeight', 'bold');
    end
end

%% ================= 5. SUMMARY PLOT: D90 TUNING CURVE =================
figure('Units', 'centimeters', 'Position', [5, 5, 12, 10], 'Color', 'w', 'Name', 'Dataset Summary D90');
hold on; 
for ss = 1:nSets_seq
    vals = [SpatialResults.Set(ss).Amp.Val];
    plot(vals, [SpatialResults.Set(ss).Amp.D90_Sim], '--ok', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
    plot(vals, [SpatialResults.Set(ss).Amp.D90_Seq], '-sk', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
end
xlabel('Amplitude (\muA)'); ylabel('Outward Spread D_{90} (\mum)'); 
title('Effective Containment Boundary'); axis square;
set(gca, 'TickDir', 'out', 'Box', 'off');
legend('Sim', 'Seq', 'Location', 'northwest', 'Box', 'off');

%% ================= 6. SAVE RESULTS =================
% --- MODIFIED: Extracts exp_id based on the sequential folder name ---
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_seq, filesep); exp_id = parts{end}; 
out_name = fullfile(save_dir, ['Result_Spatial_D90_DX005_' exp_id '.mat']);
save(out_name, 'SpatialResults', 'Amps', 'dist_bin_edges', 'active_N');
fprintf('\n>>> Spatial Data Saved: %s\n', out_name);

%% ================= HELPER FUNCTIONS =================
function d_cum = calculate_R_cumulative(dists, resp, threshold)
    d_cum = 0; active_idx = find(resp > 0);
    if isempty(active_idx), return; end
    
    [sorted_dist, ~] = sort(dists(active_idx));
    n_total = length(active_idx);
    cutoff_idx = ceil(threshold * n_total);
    d_cum = sorted_dist(cutoff_idx);
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

function metric = calc_binned_prob(dist, resp, valid, edges, mode, total_active_on_shank)
    metric = nan(1, length(edges)-1);
    for b = 1:length(edges)-1
        mask = dist >= edges(b) & dist < edges(b+1) & valid;
        if sum(mask) > 0
            if strcmp(mode, 'prob')
                metric(b) = sum(resp(mask)) / sum(mask); 
            elseif strcmp(mode, 'density')
                if total_active_on_shank > 0
                    metric(b) = sum(resp(mask)) / total_active_on_shank; 
                else
                    metric(b) = 0;
                end
            end
        elseif strcmp(mode, 'density')
            metric(b) = 0; 
        end
    end
end

function [resp, valid] = get_resp_vector(R, s_idx, a_idx, p_idx, ch_subset, QC)
    resp = zeros(length(ch_subset), 1); valid = true(length(ch_subset), 1);
    try
        chn_list = R.set(s_idx).amp(a_idx).ptd(p_idx).channel;
        for i = 1:length(ch_subset)
            c = ch_subset(i);
            if ~isempty(QC.BadCh) && s_idx <= length(QC.BadCh) && ismember(c, QC.BadCh{s_idx})
                valid(i) = false; continue; 
            end
            if c <= length(chn_list) && isfield(chn_list(c), 'is_responsive') && chn_list(c).is_responsive
                resp(i) = 1; 
            end
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