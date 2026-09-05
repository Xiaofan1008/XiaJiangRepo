%% ============================================================
%   RESPONSE DURATION ANALYSIS (Temporal Span - 3*SD Absolute)
%   - Metric: Total Duration = Time spent above 3*SD baseline threshold
%   - Method: 
%       1. Compute PSTH (1ms bins) and smooth (10ms Gaussian).
%       2. Calculate Baseline Statistics (-50 to -10ms).
%       3. Threshold = Mean + 3*SD (min floor of 10 Hz).
%       4. Duration = Sum of bins exceeding threshold in analysis window.
%   - Structure: Separate folders for Sim and Seq datasets.
%   - Logic: Order-sensitive Sequential sets matched to Order-agnostic Sim sets.
% ============================================================
clear; close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% Separate folders for Sim and Seq data
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX09/Xia_Exp1_Sim3'; 
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX09/Xia_Exp1_Seq3'; 

Electrode_Type = 2;

% 1. Duration Analysis Parameters
baseline_win_ms = [-50 -10];   % Window to calc Mean & SD (Background noise)
analysis_win_ms = [2 20];      % Window to search for response (Onset must be here)
psth_bin_ms     = 1;           % Bin size 

% 2. Threshold Rules
std_threshold   = 3;           % Threshold = Mean + N * SD 

% 3. Condition Settings
Seq_PTD_Target  = 5;        

% 4. Plotting
fig_pos = [100 100 800 500];

%% =================== 1. LOAD & PARSE DATA ====================
fprintf('Loading Simultaneous Data...\n');
[R_sim, sp_sim, trig_sim, S_sim, QC_sim] = load_experiment_data(folder_sim);
[Amps_sim, ampIdx_sim, ptdIdx_sim, ptd_sim_idx, uniqueComb_sim, combClass_sim] = parse_stim_params(S_sim, 0);

fprintf('Loading Sequential Data...\n');
[R_seq, sp_seq, trig_seq, S_seq, QC_seq] = load_experiment_data(folder_seq);
[Amps_seq, ampIdx_seq, ptdIdx_seq, ptd_seq_idx, uniqueComb_seq, combClass_seq] = parse_stim_params(S_seq, Seq_PTD_Target);

% Use Sequential combinations as the master list (since order matters here)
nSets = size(uniqueComb_seq, 1);
Amps = Amps_seq; % Assume amplitude range is defined by the main experiment file

d = Depth_s(Electrode_Type); 
nCh_Total = length(d);
FS = 30000;

%% ================= 2. ANALYZE DURATION =================
DurationData = struct();
fprintf('Analyzing Response Duration (%d Sets, %d Amps)...\n', nSets, length(Amps));

% Define Edges for PSTH
full_edges = (baseline_win_ms(1)-10) : psth_bin_ms : (analysis_win_ms(2)+10);
full_centers = full_edges(1:end-1) + psth_bin_ms/2;

for ss = 1:nSets
    % Master is Sequential Set order
    stimCh_seq = uniqueComb_seq(ss,:); stimCh_seq = stimCh_seq(stimCh_seq>0);
    set_label = ['Seq Order (Ch ' num2str(stimCh_seq) ')'];
    DurationData.Set(ss).Label = set_label;
    
    % Find the matching Simultaneous Set (Order independent)
    sim_ss_idx = [];
    for idx = 1:size(uniqueComb_sim, 1)
        stimCh_sim = uniqueComb_sim(idx,:); stimCh_sim = stimCh_sim(stimCh_sim>0);
        if isequal(sort(stimCh_sim), sort(stimCh_seq))
            sim_ss_idx = idx; break;
        end
    end
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        durations_sim = []; durations_seq = [];
        
        % --- A. Simultaneous Data Extraction ---
        if ~isempty(ptd_sim_idx) && ~isempty(sim_ss_idx)
            % Use the dynamically matched sim_ss_idx
            tr_sim = find(ampIdx_sim == ai & combClass_sim == sim_ss_idx & ptdIdx_sim == ptd_sim_idx);
            try
                resp_list = R_sim.set(sim_ss_idx).amp(ai).ptd(ptd_sim_idx).channel;
                for ch = 1:min(length(resp_list), nCh_Total)
                    if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                        if ~isempty(QC_sim.BadCh) && sim_ss_idx <= length(QC_sim.BadCh) && ismember(ch, QC_sim.BadCh{sim_ss_idx}), continue; end
                        
                        dur = calc_duration(tr_sim, trig_sim, sp_sim{d(ch)}, full_edges, full_centers, ...
                                            baseline_win_ms, analysis_win_ms, std_threshold, FS, QC_sim.BadTrials, ch);
                        if ~isnan(dur), durations_sim = [durations_sim; dur]; end 
                    end
                end
            catch; end
        end
        
        % --- B. Sequential Data Extraction ---
        if ~isempty(ptd_seq_idx)
            tr_seq = find(ampIdx_seq == ai & combClass_seq == ss & ptdIdx_seq == ptd_seq_idx);
            try
                resp_list = R_seq.set(ss).amp(ai).ptd(ptd_seq_idx).channel;
                for ch = 1:min(length(resp_list), nCh_Total)
                    if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                        if ~isempty(QC_seq.BadCh) && ss <= length(QC_seq.BadCh) && ismember(ch, QC_seq.BadCh{ss}), continue; end
                        
                        dur = calc_duration(tr_seq, trig_seq, sp_seq{d(ch)}, full_edges, full_centers, ...
                                            baseline_win_ms, analysis_win_ms, std_threshold, FS, QC_seq.BadTrials, ch);
                        if ~isnan(dur), durations_seq = [durations_seq; dur]; end
                    end
                end
            catch; end
        end
        
        med_sim = median(durations_sim); med_seq = median(durations_seq);
        if isnan(med_sim), med_sim = 0; end
        if isnan(med_seq), med_seq = 0; end
        
        fprintf('Set %d | Amp%4.1fµA: Sim:Median = %4.1f ms (N=%2d) | Seq:Median = %4.1f ms (N=%2d)\n', ...
            ss, curr_amp, med_sim, length(durations_sim), med_seq, length(durations_seq));
        
        % Store Data
        DurationData.Set(ss).Amp(ai).Val = curr_amp;
        DurationData.Set(ss).Amp(ai).Sim = durations_sim;
        DurationData.Set(ss).Amp(ai).Seq = durations_seq;
    end
end

%% ================= 3. PLOT BOXPLOTS (Per Amplitude) =================
fprintf('Generating Box Plots...\n');
warning('off', 'all');
for ss = 1:nSets
    figTitle = sprintf('Response Duration (3*SD) - %s', DurationData.Set(ss).Label);
    figure('Color','w', 'Position', fig_pos, 'Name', figTitle);
    t = tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');
    title(t, figTitle, 'FontSize', 16, 'FontWeight','bold');
    
    for ai = 1:length(Amps)
        curr_amp = DurationData.Set(ss).Amp(ai).Val;
        sim_data = DurationData.Set(ss).Amp(ai).Sim;
        seq_data = DurationData.Set(ss).Amp(ai).Seq;
        
        if isempty(sim_data) && isempty(seq_data), continue; end
        
        ax = nexttile; axes(ax); hold on; 
        
        data_all = [sim_data; seq_data];
        groups   = [ones(size(sim_data)); 2*ones(size(seq_data))];
        has_sim = ~isempty(sim_data); has_seq = ~isempty(seq_data);
        
        current_labels = {};
        if has_sim, current_labels{end+1} = 'Sim'; end
        if has_seq, current_labels{end+1} = 'Seq'; end
        
        if ~isempty(data_all)
            boxplot(data_all, groups, 'Labels', current_labels, 'Widths', 0.5, 'Symbol','o');
            h = findobj(gca, 'Tag', 'Box');
            if has_seq && has_sim
                if length(h) >= 1, patch(get(h(1),'XData'), get(h(1),'YData'), 'k', 'FaceAlpha', 0.1); end
                if length(h) >= 2, patch(get(h(2),'XData'), get(h(2),'YData'), 'w', 'FaceAlpha', 0.1); end
            elseif has_seq && ~has_sim
                if length(h) >= 1, patch(get(h(1),'XData'), get(h(1),'YData'), 'k', 'FaceAlpha', 0.1); end
            elseif ~has_seq && has_sim
                if length(h) >= 1, patch(get(h(1),'XData'), get(h(1),'YData'), 'w', 'FaceAlpha', 0.1); end
            end
        end
        
        ylabel('Duration (ms)');
        title(sprintf('%.1f \\muA', curr_amp), 'FontWeight','bold');
        box off; ylim([0 30]); 
    end
end
warning('on', 'all');

%% ================= 4. SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_seq, filesep); 
exp_id = parts{end}; % Grab the parent folder name (e.g., DX013) to name the file
out_filename = fullfile(save_dir, ['Result_Duration_3SD_DX009_' exp_id '.mat']);

save(out_filename, 'DurationData', 'Amps', 'baseline_win_ms', 'analysis_win_ms', 'std_threshold');
fprintf('\n>>> Duration Data (3*SD) Saved: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
% Parses stimulation parameters and returns PTD indices correctly
function [Amps, ampIdx, ptdIdx, target_ptd_idx, uniqueComb, combClass] = parse_stim_params(S, target_ptd)
    Stim = S.StimParams;
    simN = S.simultaneous_stim;
    if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end
    amps_all  = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
    [Amps,~,ampIdx] = unique(trialAmps); Amps(Amps==-1) = 0;
    
    if simN > 1, PTD_all_us = cell2mat(Stim(3:simN:end,6)); else, PTD_all_us = zeros(nTr,1); end
    PTDs_ms = PTD_all_us / 1000;
    [unique_ptds,~,ptdIdx] = unique(PTDs_ms);
    target_ptd_idx = find(abs(unique_ptds - target_ptd) < 0.001);
    
    stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, S.E_MAP(2:end));
    comb = zeros(nTr, simN);
    for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
    [uniqueComb,~,combClass] = unique(comb,'rows','stable');
end

% Calculates Total Duration using smoothed absolute 3*SD threshold
function duration = calc_duration(tr_ids, trig, sp_data, edges, centers, base_win, anal_win, nSD, FS, BadTrials, ch_idx)
    if ~isempty(BadTrials) && ch_idx <= length(BadTrials)
        tr_ids = setdiff(tr_ids, BadTrials{ch_idx});
    end
    if isempty(tr_ids), duration = NaN; return; end
    
    all_spikes = [];
    for k = 1:length(tr_ids)
        tr = tr_ids(k); t0 = trig(tr)/FS*1000;
        win_pad = [edges(1), edges(end)];
        tt = sp_data(:,1) - t0;
        mask = tt >= win_pad(1) & tt <= win_pad(2);
        all_spikes = [all_spikes; tt(mask)]; %#ok<AGROW>
    end
    if isempty(all_spikes), duration = NaN; return; end
    
    bin_size_ms = edges(2) - edges(1);
    counts = histcounts(all_spikes, edges);
    rate_hz = (counts / length(tr_ids)) * (1000 / bin_size_ms);
    
    smooth_rate = smoothdata(rate_hz, 'gaussian', 10);
    
    base_mask = centers >= base_win(1) & centers <= base_win(2);
    base_data = smooth_rate(base_mask);
    
    mu = mean(base_data);
    sigma = std(base_data);
    if sigma == 0, sigma = 1; end 
    
    threshold = max(mu + nSD * sigma, 10);
    
    anal_mask = centers >= anal_win(1) & centers <= anal_win(2);
    anal_rate = smooth_rate(anal_mask);
    
    if ~any(anal_rate > threshold)
        duration = NaN; 
    else
        valid_mask = (anal_rate > threshold);
        duration = sum(valid_mask) * bin_size_ms;
    end
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