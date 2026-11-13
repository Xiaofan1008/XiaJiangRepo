%% ============================================================
%   RESPONSE DURATION ANALYSIS (Temporal Span) - DUAL DATASET
%   - Metric: Duration = End_Time - Start_Time (ms)
%   - Method: Consecutive Bin Thresholding (Maunsell-style)
%       * Onset:  2 consecutive bins > Mean + N*SD
%       * Offset: 5 consecutive bins < Mean + N*SD (Bridging Gap)
%   - Logic: 
%       1. Loads Simultaneous Data -> Calculates Durations -> Stores in .Sim
%       2. Loads Sequential Data   -> Calculates Durations -> Stores in .Seq
%       3. Matches data based on Set Index and Amplitude.
%   - Output: Box Plots (Sim vs Seq) per Amplitude.
%   - Saves: 'DurationData' for Group Analysis.
% ============================================================
clear;close all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% --- Define Separate Data Folders ---
sim_data_folder = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Exp1_Sim4'; % Path to Sim Data
seq_data_folder = '/Volumes/MACData/Data/Data_Xia/DX006/Xia_Exp1_Seq4'; % Path to Seq Data

Electrode_Type = 1;

% 1. Duration Analysis Parameters
baseline_win_ms = [-50 -10];   % Window to calc Mean & SD (Background noise)
analysis_win_ms = [2 15];      % Window to search for response Onset
psth_bin_ms     = 1;           % Bin size (1ms is standard for latency/duration)

% 2. Threshold Rules
% Why 2*SD? A 2*SD threshold (P < 0.05) is more sensitive to the "tails" 
% of the response than 3*SD. We couple this with a "Consecutive Bin" rule
% to ensure robustness against random noise.
std_threshold   = 2;          % Threshold = Mean + N * SD
n_bins_onset    = 2;          % Require N consecutive bins > Thresh to start
n_bins_offset   = 3;          % Require N consecutive bins < Thresh to stop 
                              % (This bridges small gaps < 3ms in the firing)

% 3. Condition Settings
% Note: This target PTD is used to identify relevant trials in the Seq file
Seq_PTD_Target  = 5.5;        

% 4. Plotting
fig_pos = [100 100 800 500];

%% ================= 1. ANALYZE BOTH DATASETS =================
% We will loop twice: 1=Simultaneous Folder, 2=Sequential Folder
folders_to_process = {sim_data_folder, seq_data_folder};
modes              = {'Sim', 'Seq'};

% Initialize Main Result Structure
DurationData = struct();

fprintf('Starting Analysis of 2 Datasets...\n');

for d_idx = 1:2
    curr_folder = folders_to_process{d_idx};
    curr_mode   = modes{d_idx}; % 'Sim' or 'Seq'
    
    fprintf('\n---> Loading %s Data from: %s\n', curr_mode, curr_folder);
    
    % --- Load Data for this specific folder ---
    [R, sp, trig, S, QC] = load_experiment_data(curr_folder);
    
    % --- Extract Stim Params ---
    Stim = S.StimParams; 
    simN = S.simultaneous_stim; 
    if isfield(S, 'n_Trials'), nTr = S.n_Trials; else, nTr = (size(Stim, 1) - 1) / simN; end
    
    % Extract Amplitudes
    amps_all  = cell2mat(Stim(2:end,16)); trialAmps = amps_all(1:simN:end);
    [Amps,~,ampIdx] = unique(trialAmps); Amps(Amps==-1) = 0;
    
    % Extract PTDs
    if simN > 1, PTD_all_us = cell2mat(Stim(3:simN:end,6)); else, PTD_all_us = zeros(nTr,1); end
    PTD_all_ms = PTD_all_us / 1000; 
    [PTDs_ms,~,ptdIdx] = unique(PTD_all_ms);
    
    % Parse Sets (Active Channels)
    stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, S.E_MAP(2:end));
    comb = zeros(nTr, simN);
    for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
    [uniqueComb,~,combClass] = unique(comb,'rows','stable'); 
    nSets = size(uniqueComb,1);
    
    % Identify Target PTD Index based on current mode
    if strcmp(curr_mode, 'Sim')
        target_ptd_idx = find(abs(PTDs_ms - 0) < 0.001);
    else
        target_ptd_idx = find(abs(PTDs_ms - Seq_PTD_Target) < 0.001);
    end
    
    % Setup Depth/Channels
    d = Depth_s(Electrode_Type); 
    nCh_Total = length(d);
    FS = 30000;
    
    % Define Edges for PSTH Calculation
    full_edges = (baseline_win_ms(1)-10) : psth_bin_ms : (analysis_win_ms(2)+10);
    full_centers = full_edges(1:end-1) + psth_bin_ms/2;

    % --- Loop through Sets and Amplitudes ---
    for ss = 1:nSets
        % Label generation
        stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
        set_label = ['Set ' num2str(ss) ' (Ch ' num2str(stimCh) ')'];
        
        % Ensure Structure Exists
        if ~isfield(DurationData, 'Set') || length(DurationData.Set) < ss
             DurationData.Set(ss).Label = set_label;
        end
        
        for ai = 1:length(Amps)
            curr_amp = Amps(ai);
            durations = [];
            
            % Find Amplitude Index in the Structure
            store_idx = -1;
            
            % We must check if .Amp exists AND if it is a valid structure.
            % When MATLAB expands DurationData.Set to a new index (e.g. Set 2), 
            % it initializes .Amp as [], which causes the crash.
            if isfield(DurationData.Set(ss), 'Amp') && isstruct(DurationData.Set(ss).Amp)
                existing_amps = [DurationData.Set(ss).Amp.Val];
                idx_match = find(abs(existing_amps - curr_amp) < 0.01);
                if ~isempty(idx_match), store_idx = idx_match; end
            end
            
            % If Amp doesn't exist (or struct was empty), create new entry
            if store_idx == -1
                % Check again if we can append or need to start at 1
                if isfield(DurationData.Set(ss), 'Amp') && isstruct(DurationData.Set(ss).Amp)
                    store_idx = length(DurationData.Set(ss).Amp) + 1;
                else
                    store_idx = 1;
                end
                
                % Initialize the structure fields safely
                DurationData.Set(ss).Amp(store_idx).Val = curr_amp;
                DurationData.Set(ss).Amp(store_idx).Sim = [];
                DurationData.Set(ss).Amp(store_idx).Seq = [];
            end
            
            % --- Process Trials ---
            if ~isempty(target_ptd_idx)
                tr_list = find(ampIdx == ai & combClass == ss & ptdIdx == target_ptd_idx);
                
                try
                    resp_list = R.set(ss).amp(ai).ptd(target_ptd_idx).channel;
                    for ch = 1:min(length(resp_list), nCh_Total)
                        if isfield(resp_list(ch), 'is_responsive') && resp_list(ch).is_responsive
                            if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch, QC.BadCh{ss}), continue; end
                            
                            % --- CORE CALCULATION ---
                            dur = calc_duration(tr_list, trig, sp{d(ch)}, full_edges, full_centers, ...
                                                baseline_win_ms, analysis_win_ms, ...
                                                std_threshold, n_bins_onset, n_bins_offset, FS, QC.BadTrials, ch);
                            if ~isnan(dur)
                                durations = [durations; dur]; %#ok<*AGROW>
                            end
                        end
                    end
                catch; end
            end
            
            % --- Store in Correct Field (.Sim or .Seq) ---
            if strcmp(curr_mode, 'Sim')
                DurationData.Set(ss).Amp(store_idx).Sim = durations;
            else
                DurationData.Set(ss).Amp(store_idx).Seq = durations;
            end
        end
    end
end

%% ================= 3. PLOT BOXPLOTS (Per Amplitude) =================
fprintf('\nGenerating Box Plots...\n');

% Re-sort Amps inside the structure just in case they were added out of order
% (Optional robustness step)

for ss = 1:length(DurationData.Set)
    if isempty(DurationData.Set(ss).Label), continue; end
    
    figTitle = sprintf('Response Duration - %s', DurationData.Set(ss).Label);
    figure('Color','w', 'Position', fig_pos, 'Name', figTitle);
    t = tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');
    title(t, figTitle, 'FontSize', 16, 'FontWeight','bold');
    
    % Access Amps for this Set
    if ~isfield(DurationData.Set(ss), 'Amp'), continue; end
    SetAmps = [DurationData.Set(ss).Amp.Val];
    [SetAmps, sortIdx] = sort(SetAmps);
    DurationData.Set(ss).Amp = DurationData.Set(ss).Amp(sortIdx); % Sort struct
    
    for ai = 1:length(SetAmps)
        curr_amp = DurationData.Set(ss).Amp(ai).Val;
        sim_data = DurationData.Set(ss).Amp(ai).Sim;
        seq_data = DurationData.Set(ss).Amp(ai).Seq;
        
        % Skip if NO data at all
        if isempty(sim_data) && isempty(seq_data), continue; end
        
        % 1. Create Tile & Make Current
        ax = nexttile; 
        axes(ax); hold on; 
        
        % 2. Prepare Data & Logic
        data_all = [sim_data; seq_data];
        groups   = [ones(size(sim_data)); 2*ones(size(seq_data))];
        
        has_sim = ~isempty(sim_data);
        has_seq = ~isempty(seq_data);
        
        % 3. Define Labels Dynamically (Handles missing Sim or Seq)
        current_labels = {};
        if has_sim, current_labels{end+1} = 'Sim'; end
        if has_seq, current_labels{end+1} = 'Seq'; end
        
        % 4. Plot (Only if data exists)
        if ~isempty(data_all)
            boxplot(data_all, groups, 'Labels', current_labels, 'Widths', 0.5, 'Symbol','o');
            
            % 5. Style Adjustments (Robust Coloring)
            h = findobj(gca, 'Tag', 'Box');
            % Note: 'h' is usually returned in reverse order (Group 2 then Group 1)
            
            if has_seq && has_sim
                % Both present: h(1)=Seq, h(2)=Sim
                if length(h) >= 1, patch(get(h(1),'XData'), get(h(1),'YData'), 'k', 'FaceAlpha', 0.1); end % Seq
                if length(h) >= 2, patch(get(h(2),'XData'), get(h(2),'YData'), 'w', 'FaceAlpha', 0.1); end % Sim
            elseif has_seq && ~has_sim
                % Only Seq present
                if length(h) >= 1, patch(get(h(1),'XData'), get(h(1),'YData'), 'k', 'FaceAlpha', 0.1); end
            elseif ~has_seq && has_sim
                % Only Sim present
                if length(h) >= 1, patch(get(h(1),'XData'), get(h(1),'YData'), 'w', 'FaceAlpha', 0.1); end
            end
        end
        
        ylabel('Duration (ms)');
        title(sprintf('%.1f \\muA', curr_amp), 'FontWeight','bold');
        box off;
        ylim([0 20]); % Adjusted slightly higher to see full range
    end
end

%% ================= 4. SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Duration_Analysis/DX006/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

% Create ID from Sim folder name
parts = split(sim_data_folder, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_Duration_Separated_' exp_id '.mat']);

% Save Analysis Parameters
save(out_filename, 'DurationData', 'baseline_win_ms', 'analysis_win_ms', ...
     'std_threshold', 'n_bins_onset', 'n_bins_offset');
fprintf('\n>>> Duration Data Saved: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function duration = calc_duration(tr_ids, trig, sp_data, edges, centers, base_win, anal_win, nSD, nOn, nOff, FS, BadTrials, ch_idx)
    % Calculates Duration using Consecutive Bin Method
    
    % 1. Filter Bad Trials
    if ~isempty(BadTrials) && ch_idx <= length(BadTrials)
        tr_ids = setdiff(tr_ids, BadTrials{ch_idx});
    end
    if isempty(tr_ids), duration = NaN; return; end
    
    % 2. Accumulate Spikes (Raster)
    all_spikes = [];
    for k = 1:length(tr_ids)
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        
        win_pad = [edges(1), edges(end)];
        tt = sp_data(:,1) - t0;
        mask = tt >= win_pad(1) & tt <= win_pad(2);
        all_spikes = [all_spikes; tt(mask)]; %#ok<AGROW>
    end
    if isempty(all_spikes), duration = NaN; return; end
    
    % 3. Compute PSTH (Counts per bin)
    counts = histcounts(all_spikes, edges);
    
    % 4. Baseline Statistics (Noise Floor)
    base_mask = centers >= base_win(1) & centers <= base_win(2);
    base_data = counts(base_mask);
    
    mu = mean(base_data);
    sigma = std(base_data);
    if sigma == 0, sigma = 1; end % Safety check
    
    threshold = mu + nSD * sigma;
    
    % 5. Find Onset (Start of Response)
    % We strictly search within the 'analysis_win_ms' (e.g., 2-15ms)
    anal_indices = find(centers >= anal_win(1) & centers <= anal_win(2));
    
    onset_idx = NaN;
    for i = 1 : length(anal_indices) - nOn + 1
        idx = anal_indices(i);
        % Rule: Must exceed threshold for 'nOn' consecutive bins
        if all(counts(idx : idx + nOn - 1) > threshold)
            onset_idx = idx;
            break; 
        end
    end
    
    if isnan(onset_idx)
        duration = NaN; return; % No onset detected -> No duration
    end
    
    % 6. Find Offset (End of Response)
    % Search forward from onset.
    offset_idx = NaN;
    search_indices = (onset_idx + 1) : length(centers) - nOff;
    
    for i = 1 : length(search_indices)
        idx = search_indices(i);
        % Rule: Must fall BELOW threshold for 'nOff' consecutive bins (Silence)
        if all(counts(idx : idx + nOff - 1) < threshold)
            offset_idx = idx; % Use the start of the silence as the cutoff
            break;
        end
    end
    
    % If response never dies out, cap it at analysis end
    if isnan(offset_idx)
        t_start = centers(onset_idx);
        t_end   = anal_win(2); 
    else
        t_start = centers(onset_idx);
        t_end   = centers(offset_idx);
    end
    
    duration = t_end - t_start;
    
    if duration < 0, duration = NaN; end
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