%% ============================================================
%   PEAK FIRING RATE LATENCY ANALYSIS (Separate Sim/Seq Files)
%   - Metric: Time to Peak of Smoothed PSTH (ms)
%   - Logic: 
%       1. Load Sim and Seq datasets separately.
%       2. For each Set & Condition, find RESPONDING channels.
%       3. Compute PSTH (0.5ms bins) -> Smooth (Gaussian).
%       4. Find time of maximum rate.
%   - Output: Saves 'PeakData' separated by Set and Amp.
%   - Plot: Grouped Bar Plots (Side-by-Side) for Sim vs Seq.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Sim';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX005/Xia_Exp1_Seq';
Electrode_Type = 2;

% 1. Analysis Parameters
post_win_ms   = [2 20];    % Analysis Window (Look for peak here)
psth_bin_ms   = 0.5;       % Bin size for PSTH calculation
smooth_sigma  = 5;         % Gaussian smoothing sigma (ms)

% 2. Condition Settings
Seq_PTD_Target = 5.5;        

% 3. Plotting
plot_histograms = true;
hist_bin_width  = 1;       % Width of bars in the final plot (ms)
fig_pos = [100 100 1200 800];

%% =================== 1. LOAD DATA ====================
% Load Sim Data
fprintf('Loading Sim Data...\n');
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);

% Load Seq Data
fprintf('Loading Seq Data...\n');
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% --- Extract Stim Params (Use Sim as Master for Amps) ---
Stim_sim = Ssim.StimParams; 
simN_sim = Ssim.simultaneous_stim; 
if isfield(Ssim, 'n_Trials'), nTr_sim = Ssim.n_Trials; else, nTr_sim = (size(Stim_sim, 1) - 1) / simN_sim; end

% Amplitudes (Sim)
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps,~,ampIdx_sim] = unique(trialAmps_sim); Amps(Amps==-1) = 0;

% Amplitudes (Seq - Just to define indices)
Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim;
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[~,~,ampIdx_seq] = unique(trialAmps_seq); 

% PTDs (Sim)
if simN_sim > 1, PTD_all_sim = cell2mat(Stim_sim(3:simN_sim:end,6)); else, PTD_all_sim = zeros(nTr_sim,1); end
PTD_ms_sim = PTD_all_sim / 1000; [PTDs_sim,~,ptdIdx_sim] = unique(PTD_ms_sim);

% PTDs (Seq)
PTD_all_seq = cell2mat(Stim_seq(3:simN_seq:end,6)); 
PTD_ms_seq = PTD_all_seq / 1000; [PTDs_seq,~,ptdIdx_seq] = unique(PTD_ms_seq);

% Parse Sets (Sim)
stimNames_sim = Stim_sim(2:end,1); [~, idx_all_sim] = ismember(stimNames_sim, Ssim.E_MAP(2:end));
comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim, rr = (t-1)*simN_sim + (1:simN_sim); v = idx_all_sim(rr); v = v(v>0); comb_sim(t,1:numel(v)) = v(:).'; end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable'); 
nSets_sim = size(uniqueComb_sim,1);

% Parse Sets (Seq)
stimNames_seq = Stim_seq(2:end,1); [~, idx_all_seq] = ismember(stimNames_seq, Sseq.E_MAP(2:end));
comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v = idx_all_seq(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); 
nSets_seq = size(uniqueComb_seq,1);

% Assume Sets are matched
nSets = nSets_sim;

% Identify PTD indices (Sim should have 0, Seq should have Target)
ptd_sim_idx = find(abs(PTDs_sim - 0) < 0.001);
ptd_seq_idx = find(abs(PTDs_seq - Seq_PTD_Target) < 0.001);

if isempty(ptd_sim_idx), warning('No Simultaneous (0ms) data found in Sim file!'); end
if isempty(ptd_seq_idx), warning('No Sequential (%.1fms) data found in Seq file!', Seq_PTD_Target); end

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
PeakData = struct();

fprintf('Analyzing Peak Latency (%d Sets, %d Amps)...\n', nSets, length(Amps));

for ss = 1:nSets
    % Get Active Channels Name for this Set (from Sim file)
    stimCh = uniqueComb_sim(ss,:); stimCh = stimCh(stimCh>0);
    set_label = ['Set ' num2str(ss) ' (Ch ' num2str(stimCh) ')'];
    
    PeakData.Set(ss).Label = set_label;
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        
        peak_times_sim = [];
        peak_times_seq = [];
        
        % --- A. Simultaneous (From Sim Data) ---
        if ~isempty(ptd_sim_idx)
            tr_sim = find(ampIdx_sim == ai & combClass_sim == ss & ptdIdx_sim == ptd_sim_idx);
            
            try
                % Check Rsim
                resp_list = Rsim.set(ss).amp(ai).ptd(ptd_sim_idx).channel;
                for ch = 1:min(length(resp_list), nCh_Total)
                    if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                        
                        % Check Bad Channel (Sim)
                        if ~isempty(QC_Sim.BadCh) && ss <= length(QC_Sim.BadCh) && ismember(ch, QC_Sim.BadCh{ss}), continue; end
                        
                        % Get Peak Time (Use Sim Spikes & Trig)
                        peak_t = get_peak_time(tr_sim, trig_sim, sp_sim{d(ch)}, calc_edges, calc_centers, kernel, FS, QC_Sim.BadTrials, ch);
                        if ~isnan(peak_t)
                            peak_times_sim = [peak_times_sim; peak_t]; %#ok<*AGROW>
                        end
                    end
                end
            catch; end
        end
        
        % --- B. Sequential (From Seq Data) ---
        if ~isempty(ptd_seq_idx)
            tr_seq = find(ampIdx_seq == ai & combClass_seq == ss & ptdIdx_seq == ptd_seq_idx);
            
            try
                % Check Rseq
                resp_list = Rseq.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
                for ch = 1:min(length(resp_list), nCh_Total)
                    if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                        
                        % Check Bad Channel (Seq)
                        if ~isempty(QC_Seq.BadCh) && ss <= length(QC_Seq.BadCh) && ismember(ch, QC_Seq.BadCh{ss}), continue; end
                        
                        % Get Peak Time (Use Seq Spikes & Trig)
                        peak_t = get_peak_time(tr_seq, trig_seq, sp_seq{d(ch)}, calc_edges, calc_centers, kernel, FS, QC_Seq.BadTrials, ch);
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
        
        legend_shown = false;

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
            
            % 4. Add Median Lines
            if ~isempty(sim_data)
                xline(median(sim_data), '--', 'Color', [0 0.3 0.8], 'LineWidth', 1.5, 'HandleVisibility','off');
            end
            if ~isempty(seq_data)
                xline(median(seq_data), '--', 'Color', [0.7 0 0], 'LineWidth', 1.5, 'HandleVisibility','off');
            end
            
            % Format
            title(sprintf('%.1f \\muA', curr_amp), 'FontWeight','bold');
            ylabel('Fraction of Channels');
            % if ai == 1
            %     ylabel('Fraction of Channels');
            %     legend('Location','best');
            % end
            if ~legend_shown
                legend('Location','best');
                legend_shown = true;
            end
            if ai == length(Amps)
                xlabel('Time from Pulse (ms)');
            end
            xlim(post_win_ms);
            ylim([0 1]); 
            box off;
        end
    end
end

%% ================= 5. SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Peak_FR_Time/DX005/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
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
        all_spikes = [all_spikes; tt(mask)]; 
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