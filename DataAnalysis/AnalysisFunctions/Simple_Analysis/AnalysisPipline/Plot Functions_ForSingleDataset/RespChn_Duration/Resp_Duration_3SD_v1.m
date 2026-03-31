%% ============================================================
%   RESPONSE DURATION ANALYSIS (Temporal Span - 3*SD Absolute)
%   - Metric: Total Duration = Time spent above 3*SD baseline threshold
%   - Method: 
%       1. Compute PSTH (1ms bins) and smooth (10ms Gaussian).
%       2. Calculate Baseline Statistics (-50 to -10ms).
%       3. Threshold = Mean + 3*SD (min floor of 10 Hz).
%       4. Duration = Sum of all bins exceeding threshold in analysis window.
%   - Output: Box Plots (Sim vs Seq) per Amplitude.
%   - Saves: 'DurationData' for Group Analysis.
% ============================================================
clear; 
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim8';
Electrode_Type = 2;

% 1. Duration Analysis Parameters
baseline_win_ms = [-50 -10];   % Window to calc Mean & SD (Background noise)
analysis_win_ms = [2 20];      % Window to search for response (Onset must be here)
psth_bin_ms     = 1;           % Bin size (1ms is standard for latency/duration)

% 2. Threshold Rules
std_threshold   = 3;           % Threshold = Mean + N * SD (User choice: 2 or 3)

% 3. Condition Settings
Seq_PTD_Target  = 5;        

% 4. Plotting
fig_pos = [100 100 800 500];

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
    % Get Active Channels Name for this Set (for label)
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    set_label = ['Set ' num2str(ss) ' (Ch ' num2str(stimCh) ')'];
    
    DurationData.Set(ss).Label = set_label;
    
    for ai = 1:length(Amps)
        curr_amp = Amps(ai);
        durations_sim = []; durations_seq = [];
        
        % --- A. Simultaneous (PTD=0) ---
        if ~isempty(ptd_sim_idx)
            tr_sim = find(ampIdx == ai & combClass == ss & ptdIdx == ptd_sim_idx);
            try
                resp_list = R.set(ss).amp(ai).ptd(ptd_sim_idx).channel;
                for ch = 1:min(length(resp_list), nCh_Total)
                    if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                        if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch, QC.BadCh{ss}), continue; end
                        
                        dur = calc_duration(tr_sim, trig, sp{d(ch)}, full_edges, full_centers, ...
                                            baseline_win_ms, analysis_win_ms, std_threshold, FS, QC.BadTrials, ch);
                        if ~isnan(dur), durations_sim = [durations_sim; dur]; end %#ok<*AGROW>
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
                        
                        dur = calc_duration(tr_seq, trig, sp{d(ch)}, full_edges, full_centers, ...
                                            baseline_win_ms, analysis_win_ms, std_threshold, FS, QC.BadTrials, ch);
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
        
        % Store Data Per Set Per Amp
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
        box off; ylim([0 30]); % Increased Y-limit slightly since total duration will naturally be longer than FWHM
    end
end
warning('on', 'all');

%% ================= 4. SAVE RESULTS =================
% Save to new 3SD directory to prevent overwriting FWHM data
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_3SD';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Duration_3SD_DX013_' exp_id '.mat']);

save(out_filename, 'DurationData', 'Amps', 'baseline_win_ms', 'analysis_win_ms', 'std_threshold');
fprintf('\n>>> Duration Data (3*SD) Saved: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
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
    
    % 10ms Gaussian smoothing for clean absolute thresholding
    smooth_rate = smoothdata(rate_hz, 'gaussian', 10);
    
    base_mask = centers >= base_win(1) & centers <= base_win(2);
    base_data = smooth_rate(base_mask);
    
    mu = mean(base_data);
    sigma = std(base_data);
    if sigma == 0, sigma = 1; end 
    
    % 3*SD Threshold with a 10 Hz minimum floor
    threshold = max(mu + nSD * sigma, 10);
    
    anal_mask = centers >= anal_win(1) & centers <= anal_win(2);
    anal_rate = smooth_rate(anal_mask);
    
    % Total Duration: Sum of all bins above the 3SD threshold
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