%% ============================================================
%   Response Magnitude: PER-SET NORMALIZATION
%   - Metric: Mean Spike Count (2-20ms)
%   - Logic: 
%       1. Loop through each Stimulation Set (e.g., Set 1, Set 2).
%       2. Identify responding channels SPECIFIC to that Set (Sim or Seq).
%       3. Calculate Raw Counts & Normalize using that Set's max.
%   - Result: Removes "dilution" from non-responding spatial areas.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX015/Xia_Seq_Sim4';
Electrode_Type = 2;

% 1. Analysis Window
post_win_ms = [2 20]; 

% 2. NORMALIZATION SETTINGS
Ref_Amp = 5;            % Target Amplitude for Normalization (uA)
min_ref_response = 0;   % Floor to avoid dividing by noise

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
Norm_Sim = nan(nCh_Total, length(Amps), nSets);
Norm_Seq = nan(nCh_Total, length(Amps), nSets); % Only for PTD=5ms
Raw_Sim_All = nan(nCh_Total, length(Amps), nSets);
Raw_Seq_All = nan(nCh_Total, length(Amps), nSets);

fprintf('Processing %d Sets...\n', nSets);

for ss = 1:nSets
    
    % --- A. Identify Union Population for THIS Set ---
    % A channel is valid if it responds to Set 'ss' in either Sim or Seq(5ms)
    local_resp_mask = false(nCh_Total, 1);
    
    % Check Sim (PTD=0)
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
    
    % Check Seq (PTD=5ms)
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
    
    if isempty(local_resp_indices)
        fprintf('  Set %d: No responsive channels found. Skipping.\n', ss);
        continue;
    end
    
    % --- B. Calculate Raw Counts & Normalize (Per Channel in this Set) ---
    for k = 1:length(local_resp_indices)
        ch_idx = local_resp_indices(k);
        recCh  = d(ch_idx);
        S_ch   = sp{recCh};
        
        bad_trs = [];
        if ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
            bad_trs = QC.BadTrials{ch_idx};
        end
        
        % Check Bad Channel for this Set
        if ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
            continue; % Leave as NaN
        end
        
        % 1. Calculate Raw Curves
        curve_sim = nan(1, length(Amps));
        curve_seq = nan(1, length(Amps));
        
        for ai = 1:length(Amps)
            % Sim
            if ~isempty(ptd0)
                tr = setdiff(find(combClass==ss & ptdIdx==ptd0 & ampIdx==ai), bad_trs);
                if ~isempty(tr), curve_sim(ai) = get_spike_count(tr, trig, S_ch, post_win_ms, FS); end
            end
            % Seq
            if ~isempty(ptd5)
                tr = setdiff(find(combClass==ss & ptdIdx==ptd5 & ampIdx==ai), bad_trs);
                if ~isempty(tr), curve_seq(ai) = get_spike_count(tr, trig, S_ch, post_win_ms, FS); end
            end
        end
        
        % Store Raw
        Raw_Sim_All(ch_idx, :, ss) = curve_sim;
        Raw_Seq_All(ch_idx, :, ss) = curve_seq;
        
        % 2. Interpolate Max at Ref_Amp
        est_sim = NaN; est_seq = NaN;
        
        valid = ~isnan(curve_sim);
        if sum(valid) >= 2, est_sim = interp1(Amps(valid), curve_sim(valid), Ref_Amp, 'linear', 'extrap');
        elseif sum(valid)==1, est_sim = curve_sim(valid); end
        
        valid = ~isnan(curve_seq);
        if sum(valid) >= 2, est_seq = interp1(Amps(valid), curve_seq(valid), Ref_Amp, 'linear', 'extrap');
        elseif sum(valid)==1, est_seq = curve_seq(valid); end
        
        % 3. Normalize (Denominator = Max of Sim or Seq at Ref_Amp)
        est_vals = [max(0, est_sim), max(0, est_seq)];
        if any(~isnan(est_vals))
            denom = max(est_vals, [], 'omitnan');
        else
            denom = NaN;
        end
        
        if ~isnan(denom) && denom > min_ref_response
            Norm_Sim(ch_idx, :, ss) = curve_sim / denom;
            Norm_Seq(ch_idx, :, ss) = curve_seq / denom;
        else
            % Fallback (Don't divide by small noise)
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
    % Check if empty
    if all(all(isnan(data_set))), continue; end
    
    AvgSim = mean(data_set, 1, 'omitnan');
    SEMSim = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));
    
    if nSets > 1, col = sim_base_col * (0.5 + 0.5*(ss/nSets)); else, col = sim_base_col; end
    
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    lbl = sprintf('Sim Set %d %d', stimCh(1),stimCh(2));
    
    if any(~isnan(AvgSim))
        % Plot individual points (jittered)
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

yline(1.0, '--k', sprintf('Ref Max @ %.0fuA (Interp)', Ref_Amp), 'HandleVisibility','off');
xlabel('Amplitude (uA)', 'FontWeight','bold'); 
ylabel(sprintf('Normalized (%.0fuA)', Ref_Amp), 'FontWeight','bold');
title(sprintf('Response Magnitude (Per-Set Norm @ %.0fuA)', Ref_Amp), 'FontWeight','bold');
legend('Location','best','Box','off'); box off; 
ylim([0 3.0]);

%% ================= SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX015/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};
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