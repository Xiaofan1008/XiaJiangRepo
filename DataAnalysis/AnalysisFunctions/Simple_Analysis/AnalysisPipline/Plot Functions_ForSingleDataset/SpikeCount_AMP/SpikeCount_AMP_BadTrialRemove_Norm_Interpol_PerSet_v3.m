%% ============================================================
%   Response Magnitude: PER-SET NORMALIZATION (Separate Folders)
%   - Metric: Mean Spike Count (2-20ms)
%   - Logic: 
%       1. Loop Set-by-Set.
%       2. Identify Union of responding channels for THAT Set (using Rsim & Rseq).
%       3. Print Responding Channel Count.
%       4. IF 5uA exists -> Use Real Value.
%       5. IF 5uA missing -> Interpolate.
%       6. Normalize (Divide by Max(Sim_5uA, Seq_5uA)).
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX007/Xia_Exp1_Sim1';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX007/Xia_Exp1_Seq';
Electrode_Type = 2;

% 1. Analysis Window
post_win_ms = [2 20]; 

% 2. NORMALIZATION SETTINGS
Ref_Amp = 5;            % Target Amplitude (uA)
min_ref_response = 0.01; % Floor to avoid dividing by 0
PTD_Ref = 5.5;
% 3. Plotting
FS = 30000;         
jitter_width = 0.2; dot_size = 20;          

%% =================== 1. LOAD DATA ====================
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% --- Extract Stim Params (Sim) ---
Stim_sim = Ssim.StimParams; simN_sim = Ssim.simultaneous_stim; E_MAP_sim = Ssim.E_MAP;
if isfield(Ssim, 'n_Trials'), nTr_sim = Ssim.n_Trials; else, nTr_sim = (size(Stim_sim, 1) - 1) / simN_sim; end
amps_all_sim  = cell2mat(Stim_sim(2:end,16)); trialAmps_sim = amps_all_sim(1:simN_sim:end);
[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim); Amps_sim(Amps_sim==-1) = 0;

% Parse Sim Sets
stimNames_sim = Stim_sim(2:end,1); [~, idx_all_sim] = ismember(stimNames_sim, E_MAP_sim(2:end));
comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim, rr = (t-1)*simN_sim + (1:simN_sim); v = idx_all_sim(rr); v = v(v>0); comb_sim(t,1:numel(v)) = v(:).'; end
[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable'); nSets_sim = size(uniqueComb_sim,1);

% --- Extract Stim Params (Seq) ---
Stim_seq = Sseq.StimParams; simN_seq = Sseq.simultaneous_stim; E_MAP_seq = Sseq.E_MAP;
if isfield(Sseq, 'n_Trials'), nTr_seq = Sseq.n_Trials; else, nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq; end
amps_all_seq  = cell2mat(Stim_seq(2:end,16)); trialAmps_seq = amps_all_seq(1:simN_seq:end);
[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq); Amps_seq(Amps_seq==-1) = 0;
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6)); PTD_all_ms = PTD_all_us / 1000; [PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

stimNames_seq = Stim_seq(2:end,1); [~, idx_all_seq] = ismember(stimNames_seq, E_MAP_seq(2:end));
comb_seq = zeros(nTr_seq, simN_seq); for t = 1:nTr_seq, rr = (t-1)*simN_seq + (1:simN_seq); v = idx_all_seq(rr); v = v(v>0); comb_seq(t,1:numel(v)) = v(:).'; end
[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable'); nSets_seq = size(uniqueComb_seq,1);

% Assume Sim sets dictate the structure (usually identical)
nSets = nSets_sim;
Amps = Amps_sim; % Use Sim Amps as reference axis

d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

%% =================== 2. PROCESS PER SET ====================
% Output Arrays: [Channels x Amps x Sets]
Norm_Sim = nan(nCh_Total, length(Amps), nSets);
Norm_Seq = nan(nCh_Total, length(Amps), nSets);
Raw_Sim_All = nan(nCh_Total, length(Amps), nSets);
Raw_Seq_All = nan(nCh_Total, length(Amps), nSets);

% Check if Ref_Amp exists exactly
ref_amp_idx = find(abs(Amps - Ref_Amp) < 0.001);
has_ref_amp = ~isempty(ref_amp_idx);

fprintf('Processing %d Sets. Ref Amp (%.1fuA) found? %d\n', nSets, Ref_Amp, has_ref_amp);

for ss = 1:nSets
    
    % --- A. Identify Union Population for THIS Set (Using Rsim and Rseq) ---
    local_resp_mask = false(nCh_Total, 1);
    
    % Check Sim (Rsim, PTD=0 usually index 1)
    % Note: R structures usually have PTD=0 at index 1 for Sim files
    try
        % Loop amps to find any response
        for ai = 1:length(Amps)
             this = Rsim.set(ss).amp(ai).ptd(1).channel;
             for ch = 1:min(length(this), nCh_Total)
                 if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                     local_resp_mask(ch) = true;
                 end
             end
        end
    catch; end
    
    % Check Seq (Rseq, find PTD=5ms)
    ptd5_seq = find(abs(PTDs_ms - PTD_Ref) < 0.001);
    if ~isempty(ptd5_seq)
        try
            for ai = 1:length(Amps)
                this = Rseq.set(ss).amp(ai).ptd(ptd5_seq).channel;
                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch) = true;
                    end
                end
            end
        catch; end
    end
    
    local_resp_indices = find(local_resp_mask);
    
    % [PRINT] Responding Channels for this Set
    stimCh = uniqueComb_sim(ss,:); stimCh = stimCh(stimCh>0);
    fprintf('  Set %d (Ch %s): %d Responding Channels\n', ss, num2str(stimCh), length(local_resp_indices));
    
    if isempty(local_resp_indices), continue; end
    
    % --- B. Calculate & Normalize ---
    for k = 1:length(local_resp_indices)
        ch_idx = local_resp_indices(k);
        recCh  = d(ch_idx);
        
        % Check Bad Channels (Sim & Seq)
        is_bad = false;
        if ~isempty(QC_Sim.BadCh)
             if iscell(QC_Sim.BadCh) 
                 if ss<=length(QC_Sim.BadCh) && ismember(ch_idx, QC_Sim.BadCh{ss}), is_bad=true; end
             elseif ismember(ch_idx, QC_Sim.BadCh), is_bad=true; end
        end
        if ~isempty(QC_Seq.BadCh)
             if iscell(QC_Seq.BadCh)
                 if ss<=length(QC_Seq.BadCh) && ismember(ch_idx, QC_Seq.BadCh{ss}), is_bad=true; end
             elseif ismember(ch_idx, QC_Seq.BadCh), is_bad=true; end
        end
        if is_bad, continue; end
        
        % 1. Get Raw Curves
        curve_sim = nan(1, length(Amps));
        curve_seq = nan(1, length(Amps));
        
        % Get Bad Trials
        bt_sim = []; if ~isempty(QC_Sim.BadTrials) && ch_idx<=length(QC_Sim.BadTrials), bt_sim = QC_Sim.BadTrials{ch_idx}; end
        bt_seq = []; if ~isempty(QC_Seq.BadTrials) && ch_idx<=length(QC_Seq.BadTrials), bt_seq = QC_Seq.BadTrials{ch_idx}; end
        
        S_ch_sim = sp_sim{recCh};
        S_ch_seq = sp_seq{recCh};
        
        for ai = 1:length(Amps)
            % Sim
            tr = setdiff(find(combClass_sim==ss & ampIdx_sim==ai), bt_sim);
            if ~isempty(tr), curve_sim(ai) = get_spike_count(tr, trig_sim, S_ch_sim, post_win_ms, FS); end
            
            % Seq (PTD=5ms)
            if ~isempty(ptd5_seq)
                tr = setdiff(find(combClass_seq==ss & ptdIdx_seq==ptd5_seq & ampIdx_seq==ai), bt_seq);
                if ~isempty(tr), curve_seq(ai) = get_spike_count(tr, trig_seq, S_ch_seq, post_win_ms, FS); end
            end
        end
        
        Raw_Sim_All(ch_idx, :, ss) = curve_sim;
        Raw_Seq_All(ch_idx, :, ss) = curve_seq;
        
        % 2. Determine Denominator (Real vs Interp)
        val_sim_ref = NaN;
        val_seq_ref = NaN;
        
        if has_ref_amp
            % Use Real Data
            val_sim_ref = curve_sim(ref_amp_idx);
            val_seq_ref = curve_seq(ref_amp_idx);
        else
            % Use Interpolation
            valid_sim = ~isnan(curve_sim);
            if sum(valid_sim) >= 2
                val_sim_ref = interp1(Amps(valid_sim), curve_sim(valid_sim), Ref_Amp, 'linear', 'extrap');
            elseif sum(valid_sim) == 1
                val_sim_ref = curve_sim(valid_sim);
            end
            
            valid_seq = ~isnan(curve_seq);
            if sum(valid_seq) >= 2
                val_seq_ref = interp1(Amps(valid_seq), curve_seq(valid_seq), Ref_Amp, 'linear', 'extrap');
            elseif sum(valid_seq) == 1
                val_seq_ref = curve_seq(valid_seq);
            end
        end
        
        % 3. Normalize
        vals = [max(0, val_sim_ref), max(0, val_seq_ref)];
        if any(~isnan(vals))
            denom = max(vals, [], 'omitnan');
        else
            denom = NaN;
        end
        
        if ~isnan(denom) && denom > min_ref_response
            Norm_Sim(ch_idx, :, ss) = curve_sim / denom;
            Norm_Seq(ch_idx, :, ss) = curve_seq / denom;
        else
            % Fallback (Raw)
            Norm_Sim(ch_idx, :, ss) = curve_sim;
            Norm_Seq(ch_idx, :, ss) = curve_seq;
        end
    end
end

%% ===================== 3. PLOT RESULTS ======================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% --- Plot Sim ---
sim_base_col = [0 0.3 0.8];
for ss = 1:nSets
    data_set = squeeze(Norm_Sim(:, :, ss));
    if all(all(isnan(data_set))), continue; end
    
    AvgSim = mean(data_set, 1, 'omitnan');
    SEMSim = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
    
    if nSets > 1, col = sim_base_col * (0.5 + 0.5*(ss/nSets)); else, col = sim_base_col; end
    
    stimCh = uniqueComb_sim(ss,:); stimCh = stimCh(stimCh>0);
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

% --- Plot Seq ---
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20]; 
for ss = 1:nSets
    data_set = squeeze(Norm_Seq(:, :, ss));
    if all(all(isnan(data_set))), continue; end
    
    AvgSeq = mean(data_set, 1, 'omitnan');
    SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
    col = set_colors(mod(ss-1,3)+1, :);
    
    stimCh = uniqueComb_seq(ss,:); stimCh = stimCh(stimCh>0);
    
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

yline(1.0, '--k', sprintf('Ref Max @ %.0fuA', Ref_Amp), 'HandleVisibility','off');
xlabel('Amplitude (uA)', 'FontWeight','bold'); 
ylabel(sprintf('Normalized (%.0fuA)', Ref_Amp), 'FontWeight','bold');
title(sprintf('Response Magnitude (Per-Set Norm @ %.0fuA)', Ref_Amp), 'FontWeight','bold');
legend('Location','best','Box','off'); box off; 
ylim([0 3.0]);

%% ================= SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX007/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(folder_sim, filesep); exp_id = parts{end};
out_filename = fullfile(save_dir, ['Result_SpikeNormGlobalRef_' num2str(Ref_Amp) 'uA_Zeroed_5ms_' exp_id '.mat']);

ResultNorm = struct();
ResultNorm.Raw_Sim = Raw_Sim_All; 
ResultNorm.Raw_Seq = Raw_Seq_All; 
ResultNorm.Norm_Sim = Norm_Sim; 
ResultNorm.Norm_Seq = Norm_Seq; 
ResultNorm.Amps = Amps;

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