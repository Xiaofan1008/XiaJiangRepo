%% ============================================================
%   Response Magnitude: SPATIAL BINNED SPREAD (SYMMETRY MAPPING)
%   - Corrected: Handles mismatched trial counts between Sim and Seq folders.
%   - Logic: 
%       1. Sequential folder is "Master" (contains 4->8 and 8->4).
%       2. Simultaneous data is mapped as a reference to both directions.
%% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
sim_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim6';
seq_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Seq6_5ms';

Electrode_Type = 1; 
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance';

post_win_ms = [2 20]; 
pre_win_ms  = [50 0]; 
FS = 30000;         

%% =================== 1. LOAD DATA ====================
fprintf('Loading Sequential (Master) Data...\n');
[R_seq, sp_seq, trig_seq, S_seq, QC_seq] = load_experiment_data(seq_folder);

fprintf('Loading Simultaneous (Reference) Data...\n');
[R_sim, sp_sim, trig_sim, S_sim, QC_sim] = load_experiment_data(sim_folder);

% --- A. Extract Sequential Params ---
Stim_seq = S_seq.StimParams; 
simN = S_seq.simultaneous_stim; 
E_MAP = S_seq.E_MAP;

amps_seq = cell2mat(Stim_seq(2:end,16)); 
trialAmps_seq = amps_seq(1:simN:end);
[Amps,~,ampIdx_seq] = unique(trialAmps_seq); Amps(Amps==-1) = 0;

nTr_seq = (size(Stim_seq, 1) - 1) / simN;
comb_seq = zeros(nTr_seq, simN);
stimNames_seq = Stim_seq(2:end,1); [~, idx_all_seq] = ismember(stimNames_seq, E_MAP(2:end));
for t = 1:nTr_seq, rr = (t-1)*simN + (1:simN); v = idx_all_seq(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); 
nSets_seq = size(uniqueComb_seq,1);

% --- B. Extract Sim Params & Create Physical Lookup ---
Stim_sim = S_sim.StimParams;
nTr_sim = (size(Stim_sim, 1) - 1) / simN;

% ---> FIXED: Calculate ampIdx_sim for the Sim folder specifically to avoid dimension errors
amps_sim = cell2mat(Stim_sim(2:end,16)); 
trialAmps_sim = amps_sim(1:simN:end);
[~,~,ampIdx_sim] = unique(trialAmps_sim); 

comb_sim = zeros(nTr_sim, simN);
stimNames_sim = Stim_sim(2:end,1); [~, idx_all_sim] = ismember(stimNames_sim, E_MAP(2:end));
for t = 1:nTr_sim, rr = (t-1)*simN + (1:simN); v = idx_all_sim(rr); v = v(v>0); comb_sim(t,1:numel(v)) = v(:).'; end
[uniqueComb_sim,~,~] = unique(comb_sim,'rows','stable');

% Sorted keys to allow order-blind matching of the physical pairs
sim_lookup_keys = sort(uniqueComb_sim, 2); 

d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

%% =================== 2. PROCESS PER SEQUENTIAL SET ====================
Raw_Sim_All = nan(nCh_Total, length(Amps), nSets_seq);
Raw_Seq_All = nan(nCh_Total, length(Amps), nSets_seq);

post_dur_ms = post_win_ms(2) - post_win_ms(1);
pre_dur_ms  = pre_win_ms(2) - pre_win_ms(1);

for ss = 1:nSets_seq
    % --- A. Find the Matching Sim Set ---
    current_seq_pair = uniqueComb_seq(ss,:);
    sorted_seq_pair = sort(current_seq_pair);
    
    sim_match_idx = find(all(sim_lookup_keys == sorted_seq_pair, 2), 1);
    has_sim_match = ~isempty(sim_match_idx);
    
    % --- B. Identify Union Population ---
    local_resp_mask = false(nCh_Total, 1);
    for ai = 1:length(Amps)
        try
            this = R_seq.set(ss).amp(ai).ptd(1).channel;
            for ch = 1:min(length(this), nCh_Total)
                if isfield(this(ch), 'is_responsive') && this(ch).is_responsive, local_resp_mask(ch) = true; end
            end
        catch; end
    end
    if has_sim_match
        for ai = 1:length(Amps)
            try
                this = R_sim.set(sim_match_idx).amp(ai).ptd(1).channel;
                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch), 'is_responsive') && this(ch).is_responsive, local_resp_mask(ch) = true; end
                end
            catch; end
        end
    end
    
    local_resp_indices = find(local_resp_mask);
    if isempty(local_resp_indices), continue; end
    
    % --- C. Calculate Spikes ---
    for k = 1:length(local_resp_indices)
        ch_idx = local_resp_indices(k);
        recCh = d(ch_idx);
        
        % Sequential Spikes (Uses Seq folder data)
        if ~(~isempty(QC_seq.BadCh) && ss <= length(QC_seq.BadCh) && ismember(ch_idx, QC_seq.BadCh{ss}))
            bad_trs_seq = []; if ~isempty(QC_seq.BadTrials) && ch_idx <= length(QC_seq.BadTrials), bad_trs_seq = QC_seq.BadTrials{ch_idx}; end
            S_ch_seq = sp_seq{recCh};
            for ai = 1:length(Amps)
                tr = setdiff(find(combClass_seq==ss & ampIdx_seq==ai), bad_trs_seq);
                if ~isempty(tr)
                    post = get_spike_count(tr, trig_seq, S_ch_seq, post_win_ms, FS);
                    pre  = get_spike_count(tr, trig_seq, S_ch_seq, pre_win_ms, FS);
                    Raw_Seq_All(ch_idx, ai, ss) = max(0, post - (pre/pre_dur_ms)*post_dur_ms);
                end
            end
        end
        
        % Simultaneous Spikes (Uses Sim folder data)
        if has_sim_match && ~(~isempty(QC_sim.BadCh) && sim_match_idx <= length(QC_sim.BadCh) && ismember(ch_idx, QC_sim.BadCh{sim_match_idx}))
            bad_trs_sim = []; if ~isempty(QC_sim.BadTrials) && ch_idx <= length(QC_sim.BadTrials), bad_trs_sim = QC_sim.BadTrials{ch_idx}; end
            S_ch_sim = sp_sim{recCh};
            
            % Identify trials in Sim folder for this physical pair
            % ---> FIXED: Uses ampIdx_sim and sim_pair_logic to ensure compatible array sizes
            sim_pair_logic = all(comb_sim == uniqueComb_sim(sim_match_idx,:), 2);
            
            for ai = 1:length(Amps)
                tr_sim = setdiff(find(sim_pair_logic & ampIdx_sim==ai), bad_trs_sim);
                if ~isempty(tr_sim)
                    post = get_spike_count(tr_sim, trig_sim, S_ch_sim, post_win_ms, FS);
                    pre  = get_spike_count(tr_sim, trig_sim, S_ch_sim, pre_win_ms, FS);
                    Raw_Sim_All(ch_idx, ai, ss) = max(0, post - (pre/pre_dur_ms)*post_dur_ms);
                end
            end
        end
    end
end

%% ===================== 3. PLOT & SAVE ======================
Results = struct();
Results.Amps = Amps;

for ss = 1:nSets_seq
    stimCh = uniqueComb_seq(ss,:);
    shank1 = ceil(stimCh(1) / 16);
    if ceil(stimCh(2) / 16) ~= shank1, continue; end
    C_idx = (stimCh(1) + stimCh(2)) / 2;
    
    Results.SpatialSets(ss).stimCh = stimCh;
    
    label_text = sprintf('Seq: Ch%d -> Ch%d', stimCh(1), stimCh(2));
    figure('Color','w','Position',[100 100 1200 800],'Name', label_text);
    
    t = tiledlayout('flow','TileSpacing','compact');
    title(t, [label_text ' vs. Simultaneous Reference'], 'FontWeight','bold');

    xlabel(t, 'Distance from Centroid (µm)', 'FontSize', 12);
    ylabel(t, 'Mean Spike Count (Baseline Subtracted)', 'FontSize', 12);


    bin_edges = -825:50:825; bin_centers = -800:50:800;
    num_bins = length(bin_centers);

    for ai = 1:length(Amps)
        ax = nexttile; hold on;
        x_dist = []; y_sim = []; y_seq = [];
        for ch_idx = 1:nCh_Total
            if ceil(ch_idx / 16) == shank1
                dist = (ch_idx - C_idx) * 50;
                s_v = Raw_Sim_All(ch_idx, ai, ss); q_v = Raw_Seq_All(ch_idx, ai, ss);
                if ~isnan(s_v) || ~isnan(q_v)
                    x_dist(end+1) = dist; y_sim(end+1) = s_v; y_seq(end+1) = q_v;
                end
            end
        end
        
        Results.SpatialSets(ss).Amp(ai).Sim_Spikes = y_sim;
        Results.SpatialSets(ss).Amp(ai).Seq_Spikes = y_seq;
        Results.SpatialSets(ss).Amp(ai).distances = x_dist;

        if ~isempty(x_dist)
            b_sim = nan(num_bins, 1); b_seq = nan(num_bins, 1);
            idx_bin = discretize(x_dist, bin_edges);
            for b = 1:num_bins
                hit = (idx_bin == b);
                if any(hit), b_sim(b) = mean(y_sim(hit),'omitnan'); b_seq(b) = mean(y_seq(hit),'omitnan'); end
            end
            b_plot = bar(ax, bin_centers, [b_sim, b_seq], 'grouped', 'BarWidth', 1.0);
            if ~isempty(b_plot)
                b_plot(1).FaceColor = [0 0.3 0.8]; b_plot(2).FaceColor = [0.85 0.33 0.10];
                b_plot(1).EdgeColor = 'none'; b_plot(2).EdgeColor = 'none';
            end
        end
        xline(ax, -100, '--k'); xline(ax, 100, '--k');
        title(ax, sprintf('%.1f \\muA', Amps(ai))); xlim([-800 800]);
    end
    linkaxes(t.Children, 'xy');
    lg = legend(b_plot, 'Simultaneous', 'Sequential', 'Box', 'off');
    lg.Layout.Tile = 'east';
end

%% =================== 4. SAVE ====================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
[~, folder_name] = fileparts(seq_folder);
out_name = fullfile(save_dir, ['SpikeCount_vs_Distance_DX012_' folder_name '.mat']);
save(out_name, 'Results');
fprintf('\n>>> Results Saved to: %s\n', out_name);

%% ==================== HELPER FUNCTIONS =========================
function count_val = get_spike_count(tr_ids, trig, sp_data, count_win, FS)
    nTr = numel(tr_ids); total_spikes = 0;
    for k = 1:nTr
        tr = tr_ids(k); t0 = trig(tr)/FS*1000; tt = sp_data(:,1) - t0;
        mask = tt >= count_win(1) & tt <= count_win(2);
        total_spikes = total_spikes + sum(mask);
    end
    if nTr > 0, count_val = total_spikes / nTr; else, count_val = NaN; end
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*_RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;
    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end
    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);
    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.BadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
    f_bt = dir('*_BadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end