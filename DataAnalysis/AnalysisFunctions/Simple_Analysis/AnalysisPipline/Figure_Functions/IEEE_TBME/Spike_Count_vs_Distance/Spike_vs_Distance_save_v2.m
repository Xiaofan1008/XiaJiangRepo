%% ============================================================
%   Response Magnitude: SPATIAL BINNED SPREAD (BASELINE SUBTRACTED)
%   - Metric: Mean Spike Count (Evoked - Expected Baseline) vs. Signed Distance
%   - Logic: 
%       1. Loop Set-by-Set.
%       2. Extract spike counts and subtract baseline (Firing Rate * Post Window).
%       3. Generate a separate figure for each Set.
%       4. Bin data into 50um segments centered on electrodes.
%       5. SAVE comprehensive results for population analysis.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim7';
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount_vs_Distance';

% 1. Analysis Window
post_win_ms = [2 20]; 
pre_win_ms  = [50 0]; % Using 50ms for baseline rate

% 2. Plotting
FS = 30000;         

%% =================== 1. LOAD DATA ====================
[R, sp, trig, S, QC] = load_experiment_data(data_folder);

% --- Extract Stim Params ---
Stim = S.StimParams; 
simN = S.simultaneous_stim; 
E_MAP = S.E_MAP;

if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end
amps_all  = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
[Amps,~,ampIdx] = unique(trialAmps); Amps(Amps==-1) = 0;

if simN > 1
    PTD_all_us = cell2mat(Stim(3:simN:end,6)); 
else
    PTD_all_us = zeros(nTr,1);
end
PTD_all_ms = PTD_all_us / 1000; 
[PTDs_ms,~,ptdIdx] = unique(PTD_all_ms);
stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP(2:end));

comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[uniqueComb,~,combClass] = unique(comb,'rows','stable'); 
nSets = size(uniqueComb,1);

d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

%% =================== 2. PROCESS PER SET ====================
% Output Arrays: [Channels x Amps x Sets]
Raw_Sim_All = nan(nCh_Total, length(Amps), nSets);
Raw_Seq_All = nan(nCh_Total, length(Amps), nSets);

post_dur_ms = post_win_ms(2) - post_win_ms(1);
pre_dur_ms  = pre_win_ms(2) - pre_win_ms(1);

for ss = 1:nSets
    local_resp_mask = false(nCh_Total, 1);
    
    % Check Sim
    ptd0 = find(abs(PTDs_ms - 0) < 0.001);
    if ~isempty(ptd0)
        for ai = 1:length(Amps)
            try
                this = R.set(ss).amp(ai).ptd(ptd0).channel;
                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch) = true;
                    end
                end
            catch; end
        end
    end
    
    % Check Seq (5ms)
    ptd5 = find(abs(PTDs_ms - 5) < 0.001);
    if ~isempty(ptd5)
        for ai = 1:length(Amps)
            try
                this = R.set(ss).amp(ai).ptd(ptd5).channel;
                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch) = true;
                    end
                end
            catch; end
        end
    end
    
    local_resp_indices = find(local_resp_mask);
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    fprintf('  Set %d (Ch %s): %d Responding Channels\n', ss, num2str(stimCh), length(local_resp_indices));
    if isempty(local_resp_indices), continue; end
    
    for k = 1:length(local_resp_indices)
        ch_idx = local_resp_indices(k);
        recCh  = d(ch_idx);
        S_ch   = sp{recCh};
        
        bad_trs = [];
        if ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
            bad_trs = QC.BadTrials{ch_idx};
        end
        if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss}), continue; end
        
        curve_sim = nan(1, length(Amps));
        curve_seq = nan(1, length(Amps));
        
        for ai = 1:length(Amps)
            if ~isempty(ptd0)
                tr = setdiff(find(combClass==ss & ptdIdx==ptd0 & ampIdx==ai), bad_trs);
                if ~isempty(tr)
                    post_val = get_spike_count(tr, trig, S_ch, post_win_ms, FS);
                    pre_val  = get_spike_count(tr, trig, S_ch, pre_win_ms, FS);
                    expected_baseline = (pre_val / pre_dur_ms) * post_dur_ms;
                    curve_sim(ai) = max(0, post_val - expected_baseline); 
                end
            end
            if ~isempty(ptd5)
                tr = setdiff(find(combClass==ss & ptdIdx==ptd5 & ampIdx==ai), bad_trs);
                if ~isempty(tr)
                    post_val = get_spike_count(tr, trig, S_ch, post_win_ms, FS);
                    pre_val  = get_spike_count(tr, trig, S_ch, pre_win_ms, FS);
                    expected_baseline = (pre_val / pre_dur_ms) * post_dur_ms;
                    curve_seq(ai) = max(0, post_val - expected_baseline); 
                end
            end
        end
        Raw_Sim_All(ch_idx, :, ss) = curve_sim;
        Raw_Seq_All(ch_idx, :, ss) = curve_seq;
    end
end

%% ===================== 3. PLOT & PREPARE DATA STRUCT ======================
% Initialize Results struct for later population analysis
Results = struct();
Results.data_folder = data_folder;
Results.Amps = Amps;
Results.post_win = post_win_ms;
Results.pre_win = pre_win_ms;

for ss = 1:nSets
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    if length(stimCh) ~= 2, continue; end
    
    shank1 = ceil(stimCh(1) / 16);
    shank2 = ceil(stimCh(2) / 16);
    if shank1 ~= shank2, continue; end
    
    C_idx = (stimCh(1) + stimCh(2)) / 2;
    
    % ---> ADDED: Capture set-specific data for saving
    Results.SpatialSets(ss).stimCh = stimCh;
    Results.SpatialSets(ss).shankID = shank1;
    Results.SpatialSets(ss).centroidIdx = C_idx;
    
    figName = sprintf('Spatial_Set_%d_Ch%d_Ch%d', ss, stimCh(1), stimCh(2));
    figure('Color', 'w', 'Position', [100 100 1200 800], 'Name', figName);
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Spatial Spread - Set %d (Stim: Ch %d & %d)', ss, stimCh(1), stimCh(2)), 'FontWeight', 'bold');
    
    bin_edges = -825:50:825; bin_centers = -800:50:800;
    num_bins = length(bin_centers);
    
    for ai = 1:length(Amps)
        ax = nexttile; hold(ax, 'on');
        
        % Data for this amplitude
        x_dist_vals = [];
        y_sim_vals = [];
        y_seq_vals = [];
        valid_ch_ids = [];
        
        for ch_idx = 1:nCh_Total
            if ceil(ch_idx / 16) == shank1
                dist = (ch_idx - C_idx) * 50;
                val_sim = Raw_Sim_All(ch_idx, ai, ss);
                val_seq = Raw_Seq_All(ch_idx, ai, ss);
                
                if ~isnan(val_sim) || ~isnan(val_seq)
                    x_dist_vals(end+1) = dist;
                    y_sim_vals(end+1)  = val_sim;
                    y_seq_vals(end+1)  = val_seq;
                    valid_ch_ids(end+1) = ch_idx;
                end
            end
        end
        
        Results.SpatialSets(ss).Amp(ai).val = Amps(ai);
        Results.SpatialSets(ss).Amp(ai).distances = x_dist_vals;
        Results.SpatialSets(ss).Amp(ai).Sim_Spikes = y_sim_vals;
        Results.SpatialSets(ss).Amp(ai).Seq_Spikes = y_seq_vals;
        Results.SpatialSets(ss).Amp(ai).ch_indices = valid_ch_ids;

        % Binning and Plotting logic (Grouped Bars)
        if ~isempty(x_dist_vals)
            b_sim = nan(num_bins, 1); b_seq = nan(num_bins, 1);
            idx_bin = discretize(x_dist_vals, bin_edges);
            for b = 1:num_bins
                hit = (idx_bin == b);
                if any(hit)
                    b_sim(b) = mean(y_sim_vals(hit), 'omitnan');
                    b_seq(b) = mean(y_seq_vals(hit), 'omitnan');
                end
            end
            b_plot = bar(ax, bin_centers, [b_sim, b_seq], 'grouped', 'BarWidth', 1.0);
            if ~isempty(b_plot)
                b_plot(1).FaceColor = [0 0.3 0.8]; b_plot(1).EdgeColor = 'none';
                b_plot(2).FaceColor = [0.85 0.33 0.10]; b_plot(2).EdgeColor = 'none';
            end
        end
        xline(ax, -100, '--k', 'HandleVisibility', 'off');
        xline(ax, 100, '--k', 'HandleVisibility', 'off');
        title(ax, sprintf('%.1f \\muA', Amps(ai))); xlim(ax, [-800 800]);
    end
    linkaxes(t.Children, 'xy');

    % Create a global legend that represents all the bars
    lg = legend(b_plot, 'Simultaneous', 'Sequential');
    lg.Layout.Tile = 'east'; % Moves it to the right of all subplots
    lg.Box = 'off';
    lg.FontName = 'Arial';
    lg.FontSize = 10;
end

%% =================== 4. SAVE RESULTS ====================
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); 
out_name = fullfile(save_dir, ['SpikeCount_vs_Distance_DX013_' parts{end} '.mat']);
save(out_name, 'Results');
fprintf('\n>>> Comprehensive Results Saved to: %s\n', out_name);

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