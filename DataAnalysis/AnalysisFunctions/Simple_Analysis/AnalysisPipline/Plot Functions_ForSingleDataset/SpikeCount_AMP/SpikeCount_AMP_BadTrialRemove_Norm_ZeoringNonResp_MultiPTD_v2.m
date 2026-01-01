%% ============================================================
%   Response Magnitude: GLOBAL MAX AT REF AMP NORMALIZATION
%   - Metric: Mean Spike Count (2-20ms)
%   - Filter 1: If Channel NOT responsive -> Count = 0
%   - Filter 2: If Channel BAD -> Count = NaN 
%   - Filter 3: Exclude BAD TRIALS
%   - Normalization: Value / MAX(All_Sim@Ref, All_Seq@Ref)
%   - Refined: Combined Dataset (Sim=0ms, Seq=5ms)
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
% [MODIFIED] Single dataset folder
data_folder = '/Volumes/MACData/Data/Data_Xia/DX013/Xia_Exp1_Seq_Sim2';
Electrode_Type = 2;

% 1. Analysis Window
post_win_ms = [2 20]; 

% 2. NORMALIZATION SETTINGS
Ref_Amp = 5;            % Look at this amplitude (uA)
min_ref_response = 1;   % Floor to avoid dividing by noise

% 3. Plotting
FS = 30000;         
jitter_width = 0.2; dot_size = 20;          

%% =================== 1. LOAD DATA ====================
[R, sp, trig, S, QC] = load_experiment_data(data_folder);

% --- Extract Stim Params ---
Stim = S.StimParams; 
simN = S.simultaneous_stim; 
E_MAP = S.E_MAP;
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
stimNames = Stim(2:end,1); [~, idx_all] = ismember(stimNames, E_MAP(2:end));
comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[uniqueComb,~,combClass] = unique(comb,'rows','stable'); 
nSets = size(uniqueComb,1);

%% ================= 2. IDENTIFY UNION POPULATION =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_channels_mask = false(nCh_Total, 1);

% Scan R for PTD=0 or PTD=5 responsive channels
for si=1:numel(R.set)
    for ai=1:numel(R.set(si).amp)
        for pi=1:numel(R.set(si).amp(ai).ptd)
            curr_ptd = R.set(si).amp(ai).ptd(pi).PTD_ms;
            
            % Check if relevant PTD (0 or 5)
            if abs(curr_ptd - 0) > 0.001 && abs(curr_ptd - 5) > 0.001
                continue; 
            end
            
            this = R.set(si).amp(ai).ptd(pi).channel; 
            for ch=1:min(length(this),nCh_Total)
                if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                    resp_channels_mask(ch)=true; 
                end
            end
        end
    end
end
resp_channels = find(resp_channels_mask); 
fprintf('Analyzing Fixed Population (Union of Sim & Seq 5ms): %d Channels\n', length(resp_channels));

%% =================== 3. COMPUTE RAW COUNTS & NORMALIZE =================
% Init Output Arrays [Channels x Amps x Sets]
Norm_Sim = nan(length(resp_channels), length(Amps), nSets);
Norm_Seq = nan(length(resp_channels), length(Amps), nSets, numel(PTDs_ms));
Raw_Sim_All = nan(length(resp_channels), length(Amps), nSets);
Raw_Seq_All = nan(length(resp_channels), length(Amps), nSets, numel(PTDs_ms));

for ci = 1:length(resp_channels)
    ch_idx = resp_channels(ci); recCh = d(ch_idx);    
    S_ch = sp{recCh};
    
    % Get Bad Trials for this channel
    bad_trs = []; 
    if ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
        bad_trs = QC.BadTrials{ch_idx}; 
    end
    
    % --- A. Calculate Raw Sim (PTD=0) ---
    ptd_sim_idx = find(abs(PTDs_ms - 0) < 0.001);
    
    if ~isempty(ptd_sim_idx)
        for ss = 1:nSets
            % Check Bad Channel
            is_bad_ch = false; 
            if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
                is_bad_ch = true;
            end
            
            for ai = 1:length(Amps)
                if is_bad_ch, Raw_Sim_All(ci, ai, ss) = NaN; continue; end
                
                % Check response
                is_resp = 0; 
                try is_resp = R.set(ss).amp(ai).ptd(ptd_sim_idx).channel(ch_idx).is_responsive; catch, end
                if ~is_resp, Raw_Sim_All(ci, ai, ss) = 0; continue; end
                
                % Filter by Set, Amp, PTD=0
                tr_ids = setdiff(find(ampIdx == ai & combClass == ss & ptdIdx == ptd_sim_idx), bad_trs); 
                if ~isempty(tr_ids)
                    Raw_Sim_All(ci, ai, ss) = get_spike_count(tr_ids, trig, S_ch, post_win_ms, FS); 
                end
            end
        end
    end
    
    % --- B. Calculate Raw Seq (PTD=5ms only) ---
    for p = 1:numel(PTDs_ms)
        
        % [MODIFIED] Filter: Only process 5ms
        if abs(PTDs_ms(p) - 5) > 0.001, continue; end
        
        for ss = 1:nSets
            is_bad_ch = false; 
            if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
                is_bad_ch = true;
            end
            
            for ai = 1:length(Amps)
                if is_bad_ch, Raw_Seq_All(ci, ai, ss, p) = NaN; continue; end
                
                is_resp = 0; 
                try is_resp = R.set(ss).amp(ai).ptd(p).channel(ch_idx).is_responsive; catch, end
                if ~is_resp, Raw_Seq_All(ci, ai, ss, p) = 0; continue; end
                
                tr_ids = setdiff(find(combClass==ss & ptdIdx==p & ampIdx==ai), bad_trs);
                if ~isempty(tr_ids)
                    Raw_Seq_All(ci, ai, ss, p) = get_spike_count(tr_ids, trig, S_ch, post_win_ms, FS); 
                end
            end
        end
    end
    
    % --- C. NORMALIZE TO GLOBAL MAX AT REFERENCE AMPLITUDE ---
    norm_denom = NaN;
    
    % 1. Get Sim values at Ref Amp
    ref_idx = find(abs(Amps - Ref_Amp) < 0.1);
    val_sim_ref = [];
    if ~isempty(ref_idx) && ~isempty(ptd_sim_idx)
        slice = Raw_Sim_All(ci, ref_idx, :);
        val_sim_ref = slice(:);
    end
    
    % 2. Get Seq values at Ref Amp (Only 5ms slices)
    val_seq_ref = [];
    if ~isempty(ref_idx)
        % We need to grab only the PTD=5 slices
        ptd5_idx = find(abs(PTDs_ms - 5) < 0.001);
        if ~isempty(ptd5_idx)
            slice = Raw_Seq_All(ci, ref_idx, :, ptd5_idx); 
            val_seq_ref = slice(:); 
        end
    end
    
    % 3. Find GLOBAL MAX (Highest of Sim or Seq-5ms)
    all_ref_vals = [val_sim_ref; val_seq_ref];
    if ~isempty(all_ref_vals)
        norm_denom = max(all_ref_vals, [], 'omitnan');
    end
    
    % 4. Apply Normalization
    if ~isnan(norm_denom) && norm_denom > min_ref_response
        Norm_Sim(ci, :, :)       = Raw_Sim_All(ci, :, :) / norm_denom;
        Norm_Seq(ci, :, :, :)    = Raw_Seq_All(ci, :, :, :) / norm_denom;
    else
        norm_denom = min_ref_response;
        Norm_Sim(ci, :, :)       = Raw_Sim_All(ci, :, :) / norm_denom;
        Norm_Seq(ci, :, :, :)    = Raw_Seq_All(ci, :, :, :) / norm_denom;
    end
end 

%% ===================== 4. PLOT RESULTS ======================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% --- Plot Sim (PTD=0) ---
sim_base_col = [0 0.3 0.8];
for ss = 1:nSets
    data_set = squeeze(Norm_Sim(:, :, ss));
    % Check if empty (e.g. if this set was purely seq)
    if all(all(isnan(data_set) | data_set==0)), continue; end
    
    AvgSim = mean(data_set, 1, 'omitnan');
    SEMSim = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
    
    if nSets > 1, col = sim_base_col * (0.5 + 0.5*(ss/nSets)); else, col = sim_base_col; end
    
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    lbl = sprintf('Sim Set %d %d', stimCh(1),stimCh(2));
    
    if any(~isnan(AvgSim))
        for i = 1:length(Amps)
            valid = ~isnan(data_set(:,i));
            if any(valid)
                x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(i);
                scatter(x_jit, data_set(valid,i), dot_size, col, 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off'); 
            end
        end
        plot_shaded_error(Amps, AvgSim, SEMSim, col);
        plot(Amps, AvgSim, '-o', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, 'DisplayName', lbl);
    end
end

% --- Plot Seq (PTD=5ms) ---
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 
for p = 1:numel(PTDs_ms)
    if abs(PTDs_ms(p) - 5) > 0.001, continue; end
    
    for ss = 1:nSets
        data_set = squeeze(Norm_Seq(:, :, ss, p));
        AvgSeq = mean(data_set, 1, 'omitnan');
        SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
        col = set_colors(mod(ss-1,3)+1, :);
        
        stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
        
        if any(~isnan(AvgSeq))
            for i = 1:length(Amps)
                valid = ~isnan(data_set(:,i));
                if any(valid)
                    x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(i); 
                    scatter(x_jit, data_set(valid,i), dot_size, col, 'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off'); 
                end
            end
            plot_shaded_error(Amps, AvgSeq, SEMSeq, col);
            plot(Amps, AvgSeq, '-s', 'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', 'DisplayName', sprintf('Seq Set %d %d', stimCh(1),stimCh(2)));
        end
    end
end

yline(1.0, '--k', sprintf('Ref Max @ %.0fuA', Ref_Amp), 'HandleVisibility','off');
xlabel('Amplitude (uA)', 'FontWeight','bold'); 
ylabel(sprintf('Normalized (%.0fuA)', Ref_Amp), 'FontWeight','bold');
title(sprintf('Response Magnitude (Global Max Norm @ %.0fuA)', Ref_Amp), 'FontWeight','bold');
legend('Location','best','Box','off'); box off; 
ylim([0 3.0]);

%% ================= SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX013/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_SpikeNormGlobalRef_' num2str(Ref_Amp) 'uA_Zeroed_5ms_' exp_id '.mat']);

ResultNorm = struct();
ResultNorm.Raw_Sim = Raw_Sim_All; ResultNorm.Raw_Seq = Raw_Seq_All; 
ResultNorm.Norm_Sim = Norm_Sim; ResultNorm.Norm_Seq = Norm_Seq; ResultNorm.Amps = Amps;

save(out_filename, 'ResultNorm');
fprintf('\n>>> Results Saved to: %s\n', out_filename);

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

function plot_shaded_error(x, y, se, col)
    if numel(x) < 2, return; end
    x=x(:)'; y=y(:)'; se=se(:)'; valid = ~isnan(y); x=x(valid); y=y(valid); se=se(valid);
    if isempty(x), return; end
    fill([x fliplr(x)], [y+se fliplr(y-se)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
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