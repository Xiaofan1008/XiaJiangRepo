%% ============================================================
%   PEAK FIRING RATE LATENCY ANALYSIS
%   - Metric: Time to Peak of Smoothed PSTH (ms)
%   - Logic: 
%       1. For each Set & Condition, find RESPONDING channels.
%       2. Compute PSTH (0.5ms bins) -> Smooth (Gaussian).
%       3. Find time of maximum rate.
%   - Output: Saves 'PeakData' separated by Set and Amp.
%   - Plot: Grouped Bar Plots (Side-by-Side) for Sim vs Seq.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Seq_Sim7';
Electrode_Type = 2;

% 1. Analysis Parameters
post_win_ms   = [2 20];    % Analysis Window (Look for peak here)
psth_bin_ms   = 0.5;       % Bin size for PSTH calculation
smooth_sigma  = 5;         % Gaussian smoothing sigma (ms)

% 2. Condition Settings
Seq_PTD_Target = 5;        

% 3. Plotting
plot_histograms = true;
hist_bin_width  = 1;       % Width of bars in the final plot (ms)
fig_pos = [100 100 1200 800];

%% =================== 1. LOAD DATA ====================
[R, sp, trig, S, QC] = load_experiment_data(data_folder);

% --- Extract Stim Params ---
Stim = S.StimParams; 
simN = S.simultaneous_stim; 
if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end

% Amplitudes
amps_all  = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
[Amps,~,ampIdx] = unique(trialAmps); Amps(Amps==-1) = 0;

% PTDs
if simN > 1
    PTD_all_us = cell2mat(Stim(3:simN:end,6)); 
else
    PTD_all_us = zeros(nTr,1);
end
PTD_all_ms = PTD_all_us / 1000; 
[PTDs_ms,~,ptdIdx] = unique(PTD_all_ms);

% Parse Sets
stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, S.E_MAP(2:end));
comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[uniqueComb,~,combClass] = unique(comb,'rows','stable'); 
nSets = size(uniqueComb,1);

% Identify PTD indices
ptd_sim_idx = find(abs(PTDs_ms - 0) < 0.001);
ptd_seq_idx = find(abs(PTDs_ms - Seq_PTD_Target) < 0.001);

if isempty(ptd_sim_idx), warning('No Simultaneous (0ms) data found!'); end
if isempty(ptd_seq_idx), warning('No Sequential (%.1fms) data found!', Seq_PTD_Target); end

d = Depth_s(Electrode_Type); 
nCh_Total = length(d);
FS = 30000;

%% ================= 2. DEFINE SMOOTHING KERNEL =================
% Create Gaussian Kernel
k_width = 3 * smooth_sigma; 
t_k = -k_width : psth_bin_ms : k_width;
kernel = exp(-t_k.^2 / (2 * smooth_sigma^2));
kernel = kernel / sum(kernel); % Normalize area to 1

% Define Calculation Edges (High Resolution)
calc_edges = post_win_ms(1) : psth_bin_ms : post_win_ms(2);
calc_centers = calc_edges(1:end-1) + psth_bin_ms/2;

%% ================= 3. ANALYZE PEAK LATENCY =================
% Store results: Organized by Set -> Amp
% PeakData.Set(s).Amp(a).Sim = [...]
PeakData = struct();

fprintf('Analyzing Peak Latency (%d Sets, %d Amps)...\n', nSets, length(Amps));

for ss = 1:nSets
    % Get Active Channels Name for this Set (for label)
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    set_label = ['Set ' num2str(ss) ' (Ch ' num2str(stimCh) ')'];
    
    PeakData.Set(ss).Label = set_label;
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        
        peak_times_sim = [];
        peak_times_seq = [];
        
        % --- A. Simultaneous (PTD=0) ---
        if ~isempty(ptd_sim_idx)
            tr_sim = find(ampIdx == ai & combClass == ss & ptdIdx == ptd_sim_idx);
            
            try
                resp_list = R.set(ss).amp(ai).ptd(ptd_sim_idx).channel;
                for ch = 1:min(length(resp_list), nCh_Total)
                    if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                        
                        % Check Bad Channel
                        if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch, QC.BadCh{ss}), continue; end
                        
                        % Get Peak Time
                        peak_t = get_peak_time(tr_sim, trig, sp{d(ch)}, calc_edges, calc_centers, kernel, FS, QC.BadTrials, ch);
                        if ~isnan(peak_t)
                            peak_times_sim = [peak_times_sim; peak_t]; %#ok<*AGROW>
                        end
                    end
                end
            catch; end
        end
        
        % --- B. Sequential (PTD=Target) ---
        if ~isempty(ptd_seq_idx)
            tr_seq = find(ampIdx == ai & combClass == ss & ptdIdx == ptd_seq_idx);
            
            try
                resp_list = R.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
                for ch = 1:min(length(resp_list), nCh_Total)
                    if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                        
                        if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch, QC.BadCh{ss}), continue; end
                        
                        peak_t = get_peak_time(tr_seq, trig, sp{d(ch)}, calc_edges, calc_centers, kernel, FS, QC.BadTrials, ch);
                        if ~isnan(peak_t)
                            peak_times_seq = [peak_times_seq; peak_t];
                        end
                    end
                end
            catch; end
        end
        
        % Store Data Per Set Per Amp
        PeakData.Set(ss).Amp(ai).Val = curr_amp;
        PeakData.Set(ss).Amp(ai).Sim = peak_times_sim;
        PeakData.Set(ss).Amp(ai).Seq = peak_times_seq;
        
    end
end

%% ================= 4. PLOT HISTOGRAMS (Grouped Bars) =================
if plot_histograms
    
    % Define Edges for the Bar Plot (Visual Bins)
    plot_edges = post_win_ms(1) : hist_bin_width : post_win_ms(2);
    % Centers for 'bar' function
    plot_centers = plot_edges(1:end-1) + hist_bin_width/2;
    
    for ss = 1:nSets
        
        % Create One Figure Per Set
        figTitle = sprintf('Peak Latency - %s', PeakData.Set(ss).Label);
        figure('Color','w', 'Position', fig_pos, 'Name', figTitle);
        t = tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');
        title(t, figTitle, 'FontSize', 16, 'FontWeight','bold');
        
        for ai = 1:length(Amps)
            curr_amp = PeakData.Set(ss).Amp(ai).Val;
            sim_data = PeakData.Set(ss).Amp(ai).Sim;
            seq_data = PeakData.Set(ss).Amp(ai).Seq;
            
            % Skip if no data
            if isempty(sim_data) && isempty(seq_data), continue; end
            
            nexttile; hold on;
            
            % 1. Calculate Histogram Counts (Normalized to Probability/Fraction)
            N_sim = histcounts(sim_data, plot_edges, 'Normalization', 'probability');
            N_seq = histcounts(seq_data, plot_edges, 'Normalization', 'probability');
            
            % 2. Plot Grouped Bars (Side-by-Side)
            % Input to bar needs to be [Bins x Groups]
            b = bar(plot_centers, [N_sim; N_seq]', 'grouped');
            
            % 3. Style the Bars
            b(1).FaceColor = [0 0.4 0.8]; % Sim (Blue)
            b(1).EdgeColor = 'none';
            b(1).DisplayName = 'Simultaneous';
            
            b(2).FaceColor = [0.8 0.1 0.1]; % Seq (Red)
            b(2).EdgeColor = 'none';
            b(2).DisplayName = 'Sequential';
            
            % 4. Add Median Lines (Optional, for reference)
            if ~isempty(sim_data)
                xline(median(sim_data), '--', 'Color', [0 0.3 0.8], 'LineWidth', 1.5, 'HandleVisibility','off');
            end
            if ~isempty(seq_data)
                xline(median(seq_data), '--', 'Color', [0.7 0 0], 'LineWidth', 1.5, 'HandleVisibility','off');
            end
            
            % Format
            title(sprintf('%.1f \\muA', curr_amp), 'FontWeight','bold');
            if ai == 1
                ylabel('Fraction of Channels');
                legend('Location','best');
            end
            if ai == length(Amps)
                xlabel('Time from Pulse (ms)');
            end
            xlim(post_win_ms);
            ylim([0 1]); % Probability is 0-1
            box off;
        end
    end
end

%% ================= 5. SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX015/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_PeakLatency_Separated_' exp_id '.mat']);

save(out_filename, 'PeakData', 'Amps', 'post_win_ms', 'smooth_sigma');
fprintf('\n>>> Peak Latency Data Saved: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function peak_time = get_peak_time(tr_ids, trig, sp_data, edges, centers, kernel, FS, BadTrials, ch_idx)
    % 1. Filter Bad Trials
    if ~isempty(BadTrials) && ch_idx <= length(BadTrials)
        tr_ids = setdiff(tr_ids, BadTrials{ch_idx});
    end
    
    if isempty(tr_ids)
        peak_time = NaN; return;
    end
    
    % 2. Accumulate Spikes
    all_spikes = [];
    for k = 1:length(tr_ids)
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        
        % Get spikes in window relative to t0
        win_pad = [edges(1)-5, edges(end)+5]; 
        
        tt = sp_data(:,1) - t0;
        mask = tt >= win_pad(1) & tt <= win_pad(2);
        all_spikes = [all_spikes; tt(mask)]; %#ok<AGROW>
    end
    
    if isempty(all_spikes)
        peak_time = NaN; return;
    end
    
    % 3. Compute Raw PSTH
    counts = histcounts(all_spikes, edges);
    
    % 4. Smooth PSTH
    smooth_counts = conv(counts, kernel, 'same');
    
    % 5. Find Peak
    [~, max_idx] = max(smooth_counts);
    peak_time = centers(max_idx);
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