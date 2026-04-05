%% ============================================================
%   SPATIAL ANALYSIS: Single Dataset (R_EXT RAW ARCHIVER) - SPLIT FILES
%   - Logic: 
%       1. Handles Sim and Seq data stored in separate recording folders.
%       2. Probe Detection: Isolates active channels to the stimulated shank.
%       3. Distance: Measures absolute distance from the outer edge of the target zone.
%       4. Neighbor Filter: Removes isolated active channels (noise, 100um tolerance).
%       5. RADIUS PREVIEW: Calculates R_ext (Maximal Responding Distance).
%       6. RAW ARCHIVE: Saves exact channel locations for population calculation.
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
% Define the two separate folders here
data_folder_seq = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Seq'; 
data_folder_sim = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Sim';

Seq_PTD_Target = 5.5; 
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Radius';

%% =================== 1. LOAD DATA ====================
fprintf('Loading Sequential Data from: %s\n', data_folder_seq);
[R_seq, S_seq, QC_seq] = load_experiment_data(data_folder_seq);
[Amps_seq, ~, ptd_seq_idx, uniqueComb_seq] = parse_stim_params(S_seq, Seq_PTD_Target);

fprintf('Loading Simultaneous Data from: %s\n', data_folder_sim);
[R_sim, S_sim, QC_sim] = load_experiment_data(data_folder_sim);
[Amps_sim, ptd_sim_idx, ~, uniqueComb_sim] = parse_stim_params(S_sim, Seq_PTD_Target);

% Safety Check: Only process amplitudes that exist in BOTH files
Amps = intersect(Amps_seq, Amps_sim);
nSets_Seq = size(uniqueComb_seq, 1);

% --- DETECT PROBE TYPE (Using Seq as reference) ---
nCh_Raw = length(R_seq.set(1).amp(1).ptd(1).channel); 
if nCh_Raw <= 32
    active_N = 32; Electrode_Type = 1; nCh_Total = 32;
else
    active_N = 16; Electrode_Type = 2; nCh_Total = 64;
end
fprintf('Detected Probe: %d channels. Denominator: %d\n', nCh_Raw, active_N);

%% ================= 2. PROCESS SPATIAL DATA =================
SpatialResults = struct();
for ss = 1:nSets_Seq
    % Extract Sequential Set
    stimCh_seq = uniqueComb_seq(ss,:); stimCh_seq = stimCh_seq(stimCh_seq > 0);
    [~, shank_id] = get_electrode_coords(stimCh_seq(1), Electrode_Type);
    
    % Find matching Simultaneous Set Index in the Sim file
    sim_ss_idx = [];
    for idx = 1:size(uniqueComb_sim, 1)
        stimCh_sim = uniqueComb_sim(idx,:); stimCh_sim = stimCh_sim(stimCh_sim > 0);
        if isequal(sort(stimCh_seq), sort(stimCh_sim))
            sim_ss_idx = idx; 
            break; 
        end
    end
    if isempty(sim_ss_idx), continue; end
    
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
    
    for i = 1:length(active_ch)
        if y_coords(i) < y_min
            dist_to_boundary(i) = y_min - y_coords(i); 
        elseif y_coords(i) > y_max
            dist_to_boundary(i) = y_coords(i) - y_max; 
        else
            dist_to_boundary(i) = 0; 
        end
    end
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai); if curr_amp == 0, continue; end
        
        % Find the correct amplitude index for both files individually
        ai_seq = find(Amps_seq == curr_amp);
        ai_sim = find(Amps_sim == curr_amp);
        
        [sim_resp, sim_valid] = get_resp_vector(R_sim, sim_ss_idx, ai_sim, ptd_sim_idx, active_ch, QC_sim);
        [seq_resp, seq_valid] = get_resp_vector(R_seq, ss, ai_seq, ptd_seq_idx, active_ch, QC_seq);
        
        % --- NEIGHBOR FILTER (100um Tolerance) ---
        clean_sim = filter_isolated_channels(sim_resp, y_coords, 100);
        clean_seq = filter_isolated_channels(seq_resp, y_coords, 100);
        
        % Total Potency
        SpatialResults.Set(ss).Amp(ai).Val = curr_amp;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Sim = (sum(clean_sim & sim_valid) / sum(sim_valid)) * 100;
        SpatialResults.Set(ss).Amp(ai).TotalPerc_Seq = (sum(clean_seq & seq_valid) / sum(seq_valid)) * 100;
        
        % --- QUICK PREVIEW CALCULATION (R_ext) ---
        SpatialResults.Set(ss).Amp(ai).R_ext_Sim = calculate_R_ext(dist_to_boundary, clean_sim);
        SpatialResults.Set(ss).Amp(ai).R_ext_Seq = calculate_R_ext(dist_to_boundary, clean_seq);
        
        % --- THE RAW ARCHIVE ---
        SpatialResults.Set(ss).Amp(ai).Dist_to_Boundary = dist_to_boundary;
        SpatialResults.Set(ss).Amp(ai).Clean_Sim        = clean_sim;
        SpatialResults.Set(ss).Amp(ai).Clean_Seq        = clean_seq;
    end
end

%% ================= 3. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     SPATIAL SUMMARY: R_ext (Maximal Responding Distance) PREVIEW   \n');
fprintf('====================================================================\n');
fprintf('%-6s | %-6s | %-12s | %-12s\n', 'Set', 'Amp', 'R_ext Sim', 'R_ext Seq');
fprintf('--------------------------------------------------------------------\n');
for ss = 1:nSets_Seq
    for ai = 1:length(Amps)
        try
            D = SpatialResults.Set(ss).Amp(ai); if isempty(D.Val) || D.Val == 0, continue; end
            fprintf('Set %02d | %3.1fuA | %5.1fum     | %5.1fum\n', ...
                ss, D.Val, D.R_ext_Sim, D.R_ext_Seq);
        catch
            continue; % Skips safely if amplitude was skipped due to mismatch
        end
    end
end

%% ================= 4. SAVE RESULTS =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder_seq, filesep); exp_id = parts{end};
% Adding '_Split' to the filename to denote it was processed from two files
out_name = fullfile(save_dir, ['Result_Spatial_RawRadius_DX005_' exp_id '.mat']);
save(out_name, 'SpatialResults', 'Amps', 'active_N');
fprintf('\n>>> Raw Spatial Data Saved: %s\n', out_name);

%% ================= HELPER FUNCTIONS =================
function r_max = calculate_R_ext(dists, resp)
    active_idx = find(resp > 0);
    if isempty(active_idx)
        r_max = 0; 
    else
        r_max = max(dists(active_idx));
    end
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
    
    % Handles finding the correct index even if one file doesn't have PTD=5
    sim_match = find(abs(unique_ptds - 0) < 0.001);
    seq_match = find(abs(unique_ptds - target_ptd) < 0.001);
    ptd_sim_idx = 1; if ~isempty(sim_match), ptd_sim_idx = sim_match(1); end
    ptd_seq_idx = 1; if ~isempty(seq_match), ptd_seq_idx = seq_match(1); end
    
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