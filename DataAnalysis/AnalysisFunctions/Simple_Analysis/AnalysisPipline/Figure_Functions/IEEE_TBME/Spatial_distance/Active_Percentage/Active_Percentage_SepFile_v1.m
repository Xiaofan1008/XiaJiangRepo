%% ============================================================
%   POTENCY ANALYSIS: Dual-File (Dual-Shank Cluster)
%   - Logic: 
%       1. Cross-Loading: Loads separate Sim and Seq folders.
%       2. Shank Pairing: Pairs shanks {1,4} or {2,3} for 32-ch cluster.
%       3. Denominator: Fixed N=32 for consistent comparison.
%   - Style: IEEE TBME Publication Style
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS =================
sim_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim5'; 
seq_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq5_5ms';

Seq_PTD_Target = 5; 
save_dir     = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Spatial_Active_Percentage';

%% =================== 1. LOAD DATA (DUAL) ====================
fprintf('Loading SIM Data: %s\n', sim_folder);
[R_sim, S_sim, QC_sim] = load_experiment_data(sim_folder);
[Amps_sim, ptd_sim_idx, ~, uniqueComb_sim] = parse_stim_params(S_sim, 0); 

fprintf('Loading SEQ Data: %s\n', seq_folder);
[R_seq, S_seq, QC_seq] = load_experiment_data(seq_folder);
[Amps_seq, ~, ptd_seq_idx, uniqueComb_seq] = parse_stim_params(S_seq, Seq_PTD_Target);

% Common ground for looping
nSets = size(uniqueComb_seq, 1);
Amps = intersect(Amps_sim, Amps_seq);

% --- CORRECTED PROBE DETECTION ---
nCh_Raw = length(R_seq.set(1).amp(1).ptd(1).channel); % Corrected from fieldnames
if nCh_Raw <= 32
    active_N = 32; Electrode_Type = 1; nCh_Total = 32;
else
    active_N = 32; Electrode_Type = 2; nCh_Total = 64; % Forced to 32 for dual-shank cluster
end
fprintf('Detected Probe: %d channels. Potency Denominator: %d\n', nCh_Raw, active_N);

%% ================= 2. PROCESS POTENCY DATA =================
PotencyResults = struct();

for ss = 1:nSets
    % Identify Stimulation Electrodes from SEQ file
    stimCh_seq = uniqueComb_seq(ss,:); stimCh_seq = stimCh_seq(stimCh_seq > 0);
    [~, stim_shank_id] = get_electrode_coords(stimCh_seq(1), Electrode_Type);
    
    % --- DUAL SHANK PAIRING LOGIC ---
    if Electrode_Type == 2 
        if ismember(stim_shank_id, [1, 4]), target_shanks = [1, 4];
        else, target_shanks = [2, 3]; end
    else
        target_shanks = 1;
    end
    
    % --- CROSS-FILE MATCHING ---
    sim_ss_idx = [];
    for idx = 1:size(uniqueComb_sim, 1)
        stimCh_sim = uniqueComb_sim(idx,:); stimCh_sim = stimCh_sim(stimCh_sim > 0);
        if isequal(sort(stimCh_seq), sort(stimCh_sim))
            sim_ss_idx = idx; break;
        end
    end
    if isempty(sim_ss_idx), continue; end
    
    % Get all channels in the 32-channel cluster
    active_ch = [];
    for ch = 1:nCh_Total
        [~, s_id] = get_electrode_coords(ch, Electrode_Type);
        if ismember(s_id, target_shanks), active_ch = [active_ch; ch]; end
    end
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        ai_sim = find(Amps_sim == curr_amp);
        ai_seq = find(Amps_seq == curr_amp);
        
        % Fetch response vectors for the full 32-channel cluster
        [sim_resp, sim_valid] = get_resp_vector(R_sim, sim_ss_idx, ai_sim, ptd_sim_idx, active_ch, QC_sim);
        [seq_resp, seq_valid] = get_resp_vector(R_seq, ss, ai_seq, ptd_seq_idx, active_ch, QC_seq);
        
        PotencyResults.Set(ss).Amp(ai).Val = curr_amp;
        PotencyResults.Set(ss).Amp(ai).TotalPerc_Sim = (sum(sim_resp & sim_valid) / sum(sim_valid)) * 100;
        PotencyResults.Set(ss).Amp(ai).TotalPerc_Seq = (sum(seq_resp & seq_valid) / sum(seq_valid)) * 100;
        PotencyResults.Set(ss).Amp(ai).N_Valid = sum(sim_valid);
    end
end

%% ================= 3. COMMAND WINDOW SUMMARY =================
fprintf('\n====================================================================\n');
fprintf('     POTENCY SUMMARY (DUAL-FILE): 32-CH CLUSTER                      \n');
fprintf('====================================================================\n');
for ss = 1:size(PotencyResults.Set, 2)
    for ai = 1:length(Amps)
        if isempty(PotencyResults.Set(ss).Amp(ai).Val), continue; end
        D = PotencyResults.Set(ss).Amp(ai);
        fprintf('Set %02d | %3.1fuA | Sim: %5.1f%% | Seq: %5.1f%% | N=%d\n', ...
            ss, D.Val, D.TotalPerc_Sim, D.TotalPerc_Seq, D.N_Valid);
    end
end

%% ================= 6. SUMMARY PLOT =================
figure('Units', 'centimeters', 'Position', [5, 5, 12, 10], 'Color', 'w', 'Name', 'Potency_Trend');
hold on;
for ss = 1:size(PotencyResults.Set, 2)
    if isempty(PotencyResults.Set(ss).Amp), continue; end
    vals = [PotencyResults.Set(ss).Amp.Val];
    h1 = plot(vals, [PotencyResults.Set(ss).Amp.TotalPerc_Sim], '--ok', 'MarkerFaceColor', 'w');
    h2 = plot(vals, [PotencyResults.Set(ss).Amp.TotalPerc_Seq], '-sk', 'MarkerFaceColor', 'k');
end
xlabel('Amplitude (\muA)'); ylabel('Active Channel % (N=32 Cluster)');
title('Network Potency (Separate Files)'); grid on; axis square;
legend([h1, h2], {'Sim', 'Seq'}, 'Location', 'northwest', 'Box', 'off');

%% ================= 7. SAVE RESULTS =================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(seq_folder, filesep); exp_id = parts{end};
out_name = fullfile(save_dir, ['Result_ActivePercentage_DX011_' exp_id '.mat']);
save(out_name, 'PotencyResults', 'Amps', 'active_N');
fprintf('\n>>> Potency Data Saved: %s\n', out_name);

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
    if type == 2 
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