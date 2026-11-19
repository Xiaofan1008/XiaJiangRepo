%% ============================================================
%   Response Magnitude: SPATIAL BINNED SPREAD (BASELINE SUBTRACTED)
%   - Metric: Mean Spike Count (Evoked - Expected Baseline) vs. Signed Distance
%   - Logic: 
%       1. Loop Set-by-Set.
%       2. Extract spike counts and subtract baseline (Firing Rate * Post Window).
%       3. Generate a separate figure for each Set.
%       4. Bin data into 50um segments centered on electrodes.
%       5. Plot subplots for all amplitudes as grouped bar charts.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_4';
Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% 1. Analysis Window
post_win_ms = [2 20]; 
pre_win_ms  = [50 0]; 

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

% ---> CHANGED: Calculate durations once so we can use them for the firing rate math
post_dur_ms = post_win_ms(2) - post_win_ms(1);
pre_dur_ms  = pre_win_ms(2) - pre_win_ms(1);

for ss = 1:nSets
    % --- A. Identify Union Population for THIS Set ---
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
    
    % --- B. Calculate Raw Spikes ---
    for k = 1:length(local_resp_indices)
        ch_idx = local_resp_indices(k);
        recCh  = d(ch_idx);
        S_ch   = sp{recCh};
        
        bad_trs = [];
        if ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
            bad_trs = QC.BadTrials{ch_idx};
        end
        if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss}), continue; end
        
        % 1. Get Raw Curves
        curve_sim = nan(1, length(Amps));
        curve_seq = nan(1, length(Amps));
        
        for ai = 1:length(Amps)
            if ~isempty(ptd0)
                tr = setdiff(find(combClass==ss & ptdIdx==ptd0 & ampIdx==ai), bad_trs);
                if ~isempty(tr)
                    % ---> CHANGED: Baseline Firing Rate Subtraction
                    post_val = get_spike_count(tr, trig, S_ch, post_win_ms, FS);
                    pre_val  = get_spike_count(tr, trig, S_ch, pre_win_ms, FS);
                    
                    baseline_rate = pre_val / pre_dur_ms; % Calculate spikes per ms
                    expected_baseline = baseline_rate * post_dur_ms; % Scale to post-win duration
                    
                    curve_sim(ai) = max(0, post_val - expected_baseline); 
                end
            end
            if ~isempty(ptd5)
                tr = setdiff(find(combClass==ss & ptdIdx==ptd5 & ampIdx==ai), bad_trs);
                if ~isempty(tr)
                    % ---> CHANGED: Baseline Firing Rate Subtraction
                    post_val = get_spike_count(tr, trig, S_ch, post_win_ms, FS);
                    pre_val  = get_spike_count(tr, trig, S_ch, pre_win_ms, FS);
                    
                    baseline_rate = pre_val / pre_dur_ms; % Calculate spikes per ms
                    expected_baseline = baseline_rate * post_dur_ms; % Scale to post-win duration
                    
                    curve_seq(ai) = max(0, post_val - expected_baseline); 
                end
            end
        end
        
        % Store Values Strictly
        Raw_Sim_All(ch_idx, :, ss) = curve_sim;
        Raw_Seq_All(ch_idx, :, ss) = curve_seq;
    end
end

%% ===================== 3. PLOT RESULTS (SPATIAL DIAGNOSTIC) ======================
for ss = 1:nSets
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    
    % We need exactly 2 stimulation channels for Simultaneous Centroid math
    if length(stimCh) ~= 2
        continue; 
    end
    
    % --- Shank Verification ---
    % Calculate which shank (1-4) the stimulators are on
    shank1 = ceil(stimCh(1) / 16);
    shank2 = ceil(stimCh(2) / 16);
    
    if shank1 ~= shank2
        fprintf('>>> Skipping Set %d: Stim channels cross multiple shanks.\n', ss);
        continue;
    end
    
    % --- Distance Math ---
    % Calculate the Centroid Index
    C_idx = (stimCh(1) + stimCh(2)) / 2;
    
    % Create Figure
    figName = sprintf('Spatial_Set_%d_Ch%d_Ch%d', ss, stimCh(1), stimCh(2));
    figure('Color', 'w', 'Position', [100 100 1200 800], 'Name', figName);
    
    % Use TiledLayout for automatic subplot grid calculation
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Spatial Spread - Set %d (Stim: Ch %d & %d)', ss, stimCh(1), stimCh(2)), 'FontWeight', 'bold');
    xlabel(t, 'Distance from Centroid (\\mum)', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel(t, 'Mean Spike Count (Baseline Subtracted)', 'FontWeight', 'bold', 'FontSize', 12);
    
    % --- Define Zero-Centered Bins ---
    % Edges shifted by 25um so the Centers land exactly on 0, 50, 100, etc.
    bin_edges = -825:50:825; 
    bin_centers = -800:50:800;
    num_bins = length(bin_centers);
    
    % Loop through Amplitudes
    for ai = 1:length(Amps)
        current_amp = Amps(ai);
        
        % Open next subplot tile
        ax = nexttile; 
        hold(ax, 'on');
        
        % Temporary arrays to hold valid plotting data
        x_vals = [];
        y_sim = [];
        y_seq = [];
        
        for ch_idx = 1:nCh_Total
            % ---> CRITICAL FILTER: Only plot if recording channel is on the same shank!
            rec_shank = ceil(ch_idx / 16);
            if rec_shank == shank1
                
                % Calculate Signed Distance from Centroid
                dist = (ch_idx - C_idx) * 50;
                
                % Grab raw data
                val_sim = Raw_Sim_All(ch_idx, ai, ss);
                val_seq = Raw_Seq_All(ch_idx, ai, ss);
                
                % Only add to plot array if data exists
                if ~isnan(val_sim) || ~isnan(val_seq)
                    x_vals(end+1) = dist;
                    y_sim(end+1)  = val_sim;
                    y_seq(end+1)  = val_seq;
                end
            end
        end
        
        % --- BIN THE DATA ---
        if ~isempty(x_vals)
            binned_sim = nan(num_bins, 1);
            binned_seq = nan(num_bins, 1);
            
            % Discretize x_vals to find which bin index (1 to num_bins) each distance falls into
            bin_indices = discretize(x_vals, bin_edges);
            
            for b = 1:num_bins
                % Find all raw data points that belong to bin 'b'
                idx_in_bin = (bin_indices == b);
                
                if any(idx_in_bin)
                    % Calculate mean spike count for this bin (ignoring NaNs)
                    binned_sim(b) = mean(y_sim(idx_in_bin), 'omitnan');
                    binned_seq(b) = mean(y_seq(idx_in_bin), 'omitnan');
                end
            end
            
            % Plot the Grouped Bar Chart
            b_plot = bar(ax, bin_centers, [binned_sim, binned_seq], 'grouped', 'BarWidth', 1.0);
            
            % Apply Colors and Names
            if ~isempty(b_plot)
                b_plot(1).FaceColor = [0 0.3 0.8];      % Sim Blue
                b_plot(1).EdgeColor = 'none';
                b_plot(1).DisplayName = 'Simultaneous';
                
                b_plot(2).FaceColor = [0.85 0.33 0.10]; % Seq Orange
                b_plot(2).EdgeColor = 'none';
                b_plot(2).DisplayName = 'Sequential';
            end
        end
        
        % Add vertical dashed lines at the exact stimulation locations (-100 and +100)
        xline(ax, -100, '--k', 'HandleVisibility', 'off');
        xline(ax, 100, '--k', 'HandleVisibility', 'off');
        
        % Formatting for the tile
        title(ax, sprintf('%.1f \\muA', current_amp));
        xlim(ax, [-800 800]); % Lock the X window to the shank distance
        
        % Only put the legend on the very first subplot to save space
        if ai == 1
            legend(ax, 'Location', 'northeast', 'Box', 'off');
        end
    end
    
    % ---> THE AXIS LOCK: Forces all subplots in this figure to share the same Y-axis scale
    linkaxes(t.Children, 'xy');
end

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