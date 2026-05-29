%% ============================================================
%   Response Magnitude: PER-SET NORMALIZATION (Separate Folders)
%   Simultaneous vs Sequential Version with Ordered Seq Sets
%
%   - Metric: Mean Spike Count (2-20 ms)
%   - Main changes:
%       1. Loop over SEQUENTIAL ordered sets.
%       2. Match each Seq ordered set to the corresponding Sim unordered set.
%       3. Allow amplitude exclusion by Mode + SeqSetIndex + Amp.
%       4. Excluded amplitudes are set to NaN, not deleted.
%       5. Ref_Amp normalization is estimated if Ref_Amp is missing/excluded.
%       6. Plot connects across excluded/missing amplitudes using valid points.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
folder_sim = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_Exp1_Sim2';
folder_seq = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_Exp1_Seq2';
Electrode_Type = 2;

% 1. Analysis Window
post_win_ms = [2 20];

% 2. NORMALIZATION SETTINGS
Ref_Amp = 5;              % Target amplitude for normalization, uA
min_ref_response = 0.01;  % Floor to avoid unstable division
PTD_Ref = 5;              % Sequential PTD to analyze, ms

% 3. Plotting
FS = 30000;
jitter_width = 0.2;
dot_size = 20;

% ========================================================================
% [NEW] Amplitude exclusion rules
%
% Format:
%   ExcludeAmpRules = {
%       % Mode, SeqSetIndex, Amp
%       'Seq', 1, 8;
%       'Sim', 2, 10;
%   };
%
% Meaning:
%   'Seq', ss_seq, amp:
%       exclude this amplitude from the sequential ordered set ss_seq.
%
%   'Sim', ss_seq, amp:
%       exclude this amplitude from the matched simultaneous curve used for
%       sequential ordered set ss_seq.
%
% Excluded amplitudes are set to NaN, not removed from the Amps axis.
% ========================================================================
ExcludeAmpRules = {
    % Mode,   SeqSetIndex,   Amp
    % 'Seq',  1,             8;
    % 'Sim',  2,             10;
    % 'Seq',  1,             8;
    % 'Seq',  2,             8;
    % 'Sim',  1,             8;
    % 'Sim',  2,             8;
    % 
    % 'Seq',  1,             10;
    % 'Seq',  2,             10;
    % 'Sim',  1,             10;
    % 'Sim',  2,             10;
};

%% =================== 1. LOAD DATA ====================
[Rsim, sp_sim, trig_sim, Ssim, QC_Sim] = load_experiment_data(folder_sim);
[Rseq, sp_seq, trig_seq, Sseq, QC_Seq] = load_experiment_data(folder_seq);

% --- Extract Stim Params (Sim) ---
Stim_sim = Ssim.StimParams;
simN_sim = Ssim.simultaneous_stim;
E_MAP_sim = Ssim.E_MAP;

if isfield(Ssim, 'n_Trials')
    nTr_sim = Ssim.n_Trials;
else
    nTr_sim = (size(Stim_sim, 1) - 1) / simN_sim;
end

amps_all_sim  = cell2mat(Stim_sim(2:end,16));
trialAmps_sim = amps_all_sim(1:simN_sim:end);
trialAmps_sim(trialAmps_sim == -1) = 0;

[Amps_sim,~,ampIdx_sim] = unique(trialAmps_sim);
Amps_sim(Amps_sim == -1) = 0;

% Parse Sim Sets
stimNames_sim = Stim_sim(2:end,1);
[~, idx_all_sim] = ismember(stimNames_sim, E_MAP_sim(2:end));

comb_sim = zeros(nTr_sim, simN_sim);
for t = 1:nTr_sim
    rr = (t-1)*simN_sim + (1:simN_sim);
    v = idx_all_sim(rr);
    v = v(v>0);
    comb_sim(t,1:numel(v)) = v(:).';
end

[uniqueComb_sim,~,combClass_sim] = unique(comb_sim,'rows','stable');
nSets_sim = size(uniqueComb_sim,1);

% --- Extract Stim Params (Seq) ---
Stim_seq = Sseq.StimParams;
simN_seq = Sseq.simultaneous_stim;
E_MAP_seq = Sseq.E_MAP;

if isfield(Sseq, 'n_Trials')
    nTr_seq = Sseq.n_Trials;
else
    nTr_seq = (size(Stim_seq, 1) - 1) / simN_seq;
end

amps_all_seq  = cell2mat(Stim_seq(2:end,16));
trialAmps_seq = amps_all_seq(1:simN_seq:end);
trialAmps_seq(trialAmps_seq == -1) = 0;

[Amps_seq,~,ampIdx_seq] = unique(trialAmps_seq);
Amps_seq(Amps_seq == -1) = 0;

% Sequential PTD
PTD_all_us = cell2mat(Stim_seq(3:simN_seq:end,6));
PTD_all_ms = PTD_all_us / 1000;
[PTDs_ms,~,ptdIdx_seq] = unique(PTD_all_ms);

ptd5_seq = find(abs(PTDs_ms - PTD_Ref) < 0.001);
if isempty(ptd5_seq)
    warning('PTD_Ref %.3f ms was not found in sequential dataset. Seq curves will be NaN.', PTD_Ref);
elseif numel(ptd5_seq) > 1
    warning('Multiple PTD entries matched PTD_Ref %.3f ms. Using the first match.', PTD_Ref);
    ptd5_seq = ptd5_seq(1);
end

% Parse Seq Sets
stimNames_seq = Stim_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP_seq(2:end));

comb_seq = zeros(nTr_seq, simN_seq);
for t = 1:nTr_seq
    rr = (t-1)*simN_seq + (1:simN_seq);
    v = idx_all_seq(rr);
    v = v(v>0);
    comb_seq(t,1:numel(v)) = v(:).';
end

[uniqueComb_seq,~,combClass_seq] = unique(comb_seq,'rows','stable');
nSets_seq = size(uniqueComb_seq,1);

% ========================================================================
% [MODIFIED]
% Use sequential ordered sets as the main analysis structure.
%
% Each sequential ordered set gets matched to a simultaneous unordered set.
% Example:
%     Seq [A B] -> Sim [A B]
%     Seq [B A] -> Sim [A B]
% ========================================================================
nSets = nSets_seq;

% Use the union of Sim and Seq amplitude values as the shared amplitude axis
Amps = unique([Amps_sim(:); Amps_seq(:)]).';
nAMP = numel(Amps);

d = Depth_s(Electrode_Type);
nCh_Total = length(d);

% Find matched Sim set for each Seq ordered set
matchedSimSet = nan(nSets_seq, 1);

for ss_seq = 1:nSets_seq

    seq_pair = uniqueComb_seq(ss_seq,:);
    seq_pair = seq_pair(seq_pair > 0);
    seq_pair_sorted = sort(seq_pair);

    for ss_sim = 1:nSets_sim
        sim_pair = uniqueComb_sim(ss_sim,:);
        sim_pair = sim_pair(sim_pair > 0);
        sim_pair_sorted = sort(sim_pair);

        if numel(seq_pair_sorted) == numel(sim_pair_sorted) && ...
                all(seq_pair_sorted == sim_pair_sorted)
            matchedSimSet(ss_seq) = ss_sim;
            break;
        end
    end
end

fprintf('\nSequential ordered set matching:\n');
for ss_seq = 1:nSets_seq
    seq_pair = uniqueComb_seq(ss_seq,:);
    seq_pair = seq_pair(seq_pair > 0);

    if isnan(matchedSimSet(ss_seq))
        fprintf('  Seq Set %d [%s] -> No matching Sim set found\n', ...
            ss_seq, num2str(seq_pair));
    else
        sim_pair = uniqueComb_sim(matchedSimSet(ss_seq),:);
        sim_pair = sim_pair(sim_pair > 0);
        fprintf('  Seq Set %d [%s] -> Sim Set %d [%s]\n', ...
            ss_seq, num2str(seq_pair), matchedSimSet(ss_seq), num2str(sim_pair));
    end
end

%% =================== 2. PROCESS PER SEQ SET ====================
% Output Arrays: [Channels x Amps x SeqOrderedSets]
Norm_Sim = nan(nCh_Total, nAMP, nSets);
Norm_Seq = nan(nCh_Total, nAMP, nSets);
Raw_Sim_All = nan(nCh_Total, nAMP, nSets);
Raw_Seq_All = nan(nCh_Total, nAMP, nSets);

% Save reference estimate values and method labels
RefVal_Sim = nan(nCh_Total, nSets);
RefVal_Seq = nan(nCh_Total, nSets);
RefMethod_Sim = strings(nCh_Total, nSets);
RefMethod_Seq = strings(nCh_Total, nSets);

fprintf('\nProcessing %d sequential ordered sets. Ref Amp = %.1f uA\n', nSets, Ref_Amp);

for ss_seq = 1:nSets

    ss_sim = matchedSimSet(ss_seq);

    if isnan(ss_sim)
        fprintf('WARNING: Seq Set %d has no matched Sim set. Skipping this set.\n', ss_seq);
        continue;
    end

    ss_sim = double(ss_sim);

    % --- A. Identify Union Population for THIS Seq Set and matched Sim Set ---
    local_resp_mask = false(nCh_Total, 1);

    % Check Sim responding channels from matched Sim set
    try
        for aiR = 1:numel(Rsim.set(ss_sim).amp)
            if isfield(Rsim.set(ss_sim).amp(aiR), 'ptd') && ...
                    ~isempty(Rsim.set(ss_sim).amp(aiR).ptd)

                this = Rsim.set(ss_sim).amp(aiR).ptd(1).channel;

                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch) = true;
                    end
                end
            end
        end
    catch ME
        fprintf('WARNING: Could not read Rsim responding channels for Sim Set %d. %s\n', ...
            ss_sim, ME.message);
    end

    % Check Seq responding channels from current ordered Seq set and PTD_Ref
    if ~isempty(ptd5_seq)
        try
            for aiR = 1:numel(Rseq.set(ss_seq).amp)
                if isfield(Rseq.set(ss_seq).amp(aiR), 'ptd') && ...
                        numel(Rseq.set(ss_seq).amp(aiR).ptd) >= ptd5_seq

                    this = Rseq.set(ss_seq).amp(aiR).ptd(ptd5_seq).channel;

                    for ch = 1:min(length(this), nCh_Total)
                        if isfield(this(ch), 'is_responsive') && this(ch).is_responsive
                            local_resp_mask(ch) = true;
                        end
                    end
                end
            end
        catch ME
            fprintf('WARNING: Could not read Rseq responding channels for Seq Set %d. %s\n', ...
                ss_seq, ME.message);
        end
    end

    local_resp_indices = find(local_resp_mask);

    seqCh = uniqueComb_seq(ss_seq,:);
    seqCh = seqCh(seqCh > 0);

    simCh = uniqueComb_sim(ss_sim,:);
    simCh = simCh(simCh > 0);

    fprintf('  Seq Set %d [%s] matched Sim Set %d [%s]: %d Responding Channels\n', ...
        ss_seq, num2str(seqCh), ss_sim, num2str(simCh), length(local_resp_indices));

    if isempty(local_resp_indices)
        continue;
    end

    % --- B. Calculate & Normalize ---
    for k = 1:length(local_resp_indices)

        ch_idx = local_resp_indices(k);
        recCh  = d(ch_idx);

        % Check Bad Channels using matched Sim set and current Seq set
        is_bad = false;

        if ~isempty(QC_Sim.BadCh)
            if iscell(QC_Sim.BadCh)
                if ss_sim <= length(QC_Sim.BadCh) && ismember(ch_idx, QC_Sim.BadCh{ss_sim})
                    is_bad = true;
                end
            elseif ismember(ch_idx, QC_Sim.BadCh)
                is_bad = true;
            end
        end

        if ~isempty(QC_Seq.BadCh)
            if iscell(QC_Seq.BadCh)
                if ss_seq <= length(QC_Seq.BadCh) && ismember(ch_idx, QC_Seq.BadCh{ss_seq})
                    is_bad = true;
                end
            elseif ismember(ch_idx, QC_Seq.BadCh)
                is_bad = true;
            end
        end

        if is_bad
            continue;
        end

        % 1. Get Raw Curves
        curve_sim = nan(1, nAMP);
        curve_seq = nan(1, nAMP);

        % Get Bad Trials
        bt_sim = [];
        if ~isempty(QC_Sim.BadTrials) && ch_idx <= length(QC_Sim.BadTrials)
            bt_sim = QC_Sim.BadTrials{ch_idx};
        end

        bt_seq = [];
        if ~isempty(QC_Seq.BadTrials) && ch_idx <= length(QC_Seq.BadTrials)
            bt_seq = QC_Seq.BadTrials{ch_idx};
        end

        S_ch_sim = sp_sim{recCh};
        S_ch_seq = sp_seq{recCh};

        for ai = 1:nAMP

            amp_val = Amps(ai);

            % ------------------------------------------------------------
            % Sim: matched Sim set, match by amplitude VALUE
            % ------------------------------------------------------------
            tr_sim = find(combClass_sim == ss_sim & ...
                          abs(trialAmps_sim - amp_val) < 0.001);
            tr_sim = setdiff(tr_sim, bt_sim);

            if ~isempty(tr_sim)
                curve_sim(ai) = get_spike_count(tr_sim, trig_sim, S_ch_sim, post_win_ms, FS);
            end

            % ------------------------------------------------------------
            % Seq: current ordered Seq set, match by amplitude VALUE and PTD
            % ------------------------------------------------------------
            if ~isempty(ptd5_seq)
                tr_seq = find(combClass_seq == ss_seq & ...
                              ptdIdx_seq == ptd5_seq & ...
                              abs(trialAmps_seq - amp_val) < 0.001);
                tr_seq = setdiff(tr_seq, bt_seq);

                if ~isempty(tr_seq)
                    curve_seq(ai) = get_spike_count(tr_seq, trig_seq, S_ch_seq, post_win_ms, FS);
                end
            end
        end

        % =================================================================
        % [NEW] Apply mode-specific, set-specific amplitude exclusion.
        % Excluded amplitudes are set to NaN before reference estimation.
        % =================================================================
        for ai = 1:nAMP
            amp_val = Amps(ai);

            if is_excluded_amp('Sim', ss_seq, amp_val, ExcludeAmpRules)
                curve_sim(ai) = NaN;
            end

            if is_excluded_amp('Seq', ss_seq, amp_val, ExcludeAmpRules)
                curve_seq(ai) = NaN;
            end
        end

        Raw_Sim_All(ch_idx, :, ss_seq) = curve_sim;
        Raw_Seq_All(ch_idx, :, ss_seq) = curve_seq;

        % =================================================================
        % [MODIFIED] Estimate Ref_Amp response if real Ref_Amp is missing
        % or excluded.
        % =================================================================
        [val_sim_ref, method_sim] = estimate_ref_response(Amps, curve_sim, Ref_Amp);
        [val_seq_ref, method_seq] = estimate_ref_response(Amps, curve_seq, Ref_Amp);

        RefVal_Sim(ch_idx, ss_seq) = val_sim_ref;
        RefVal_Seq(ch_idx, ss_seq) = val_seq_ref;
        RefMethod_Sim(ch_idx, ss_seq) = method_sim;
        RefMethod_Seq(ch_idx, ss_seq) = method_seq;

        % Denominator uses the larger valid Ref_Amp response across Sim/Seq
        vals = [max(0, val_sim_ref), max(0, val_seq_ref)];

        if any(~isnan(vals))
            denom = max(vals, [], 'omitnan');
        else
            denom = NaN;
        end

        if ~isnan(denom) && denom > min_ref_response
            Norm_Sim(ch_idx, :, ss_seq) = curve_sim / denom;
            Norm_Seq(ch_idx, :, ss_seq) = curve_seq / denom;
        else
            Norm_Sim(ch_idx, :, ss_seq) = NaN;
            Norm_Seq(ch_idx, :, ss_seq) = NaN;
        end
    end
end

%% ===================== 3. PLOT RESULTS ======================
figure('Color','w', 'Position',[100 100 700 500]); hold on;

% --- Plot Sim ---
sim_base_col = [0 0.3 0.8];

for ss_seq = 1:nSets

    data_set = squeeze(Norm_Sim(:, :, ss_seq));

    if all(all(isnan(data_set)))
        continue;
    end

    AvgSim = mean(data_set, 1, 'omitnan');
    SEMSim = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));

    if nSets > 1
        col = sim_base_col * (0.5 + 0.5*(ss_seq/nSets));
    else
        col = sim_base_col;
    end

    ss_sim = matchedSimSet(ss_seq);
    if isnan(ss_sim)
        continue;
    end

    simCh = uniqueComb_sim(ss_sim,:);
    simCh = simCh(simCh > 0);

    lbl = sprintf('Sim Set %d %d', simCh(1), simCh(2));

    if any(~isnan(AvgSim))

        % Scatter only valid data points
        for i = 1:length(Amps)
            valid = ~isnan(data_set(:,i));
            if any(valid)
                x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(i);
                scatter(x_jit, data_set(valid,i), dot_size, col, ...
                    'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off');
            end
        end

        % Plot line using only valid mean points so NaN does not interrupt line
        valid_line = ~isnan(AvgSim);
        if any(valid_line)
            plot_shaded_error(Amps(valid_line), AvgSim(valid_line), SEMSim(valid_line), col);
            plot(Amps(valid_line), AvgSim(valid_line), '-o', ...
                'Color', col, 'LineWidth', 2, 'MarkerFaceColor', col, ...
                'DisplayName', lbl);
        end
    end
end

% --- Plot Seq ---
set_colors = [0.85 0.33 0.10; 0.60 0.20 0.60; 0.20 0.60 0.20];

for ss_seq = 1:nSets

    data_set = squeeze(Norm_Seq(:, :, ss_seq));

    if all(all(isnan(data_set)))
        continue;
    end

    AvgSeq = mean(data_set, 1, 'omitnan');
    SEMSeq = std(data_set, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data_set), 1));

    col = set_colors(mod(ss_seq-1,3)+1, :);

    seqCh = uniqueComb_seq(ss_seq,:);
    seqCh = seqCh(seqCh > 0);

    if any(~isnan(AvgSeq))

        % Scatter only valid data points
        for i = 1:length(Amps)
            valid = ~isnan(data_set(:,i));
            if any(valid)
                x_jit = (rand(sum(valid),1) - 0.5) * jitter_width + Amps(i);
                scatter(x_jit, data_set(valid,i), dot_size, col, ...
                    'filled', 'MarkerFaceAlpha', 0.2, 'HandleVisibility','off');
            end
        end

        % Plot line using only valid mean points so NaN does not interrupt line
        valid_line = ~isnan(AvgSeq);
        if any(valid_line)
            plot_shaded_error(Amps(valid_line), AvgSeq(valid_line), SEMSeq(valid_line), col);
            plot(Amps(valid_line), AvgSeq(valid_line), '-s', ...
                'Color', col, 'LineWidth', 2, 'MarkerFaceColor', 'w', ...
                'DisplayName', sprintf('Seq Set %d>%d', seqCh(1), seqCh(2)));
        end
    end
end

yline(1.0, '--k', sprintf('Ref Max @ %.0f uA', Ref_Amp), 'HandleVisibility','off');
xlabel('Amplitude (uA)', 'FontWeight','bold');
ylabel(sprintf('Normalized (estimated %.0f uA ref)', Ref_Amp), 'FontWeight','bold');
title(sprintf('Response Magnitude: Ordered Seq Sets, Per-Set Norm @ %.0f uA', Ref_Amp), ...
    'FontWeight','bold');
legend('Location','best','Box','off');
box off;
ylim([0 3.0]);

%% ================= SAVE RESULTS =================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/SpikeCount/DX018/';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

parts = split(folder_sim, filesep);
exp_id = parts{end};

out_filename = fullfile(save_dir, ...
    ['Result_SpikeNormGlobalRef_' num2str(Ref_Amp) ...
    'uA_Zeroed_5ms_' exp_id '.mat']);

ResultNorm = struct();
ResultNorm.Raw_Sim = Raw_Sim_All;
ResultNorm.Raw_Seq = Raw_Seq_All;
ResultNorm.Norm_Sim = Norm_Sim;
ResultNorm.Norm_Seq = Norm_Seq;
ResultNorm.Amps = Amps;

% Save set matching and metadata
ResultNorm.uniqueComb_sim = uniqueComb_sim;
ResultNorm.uniqueComb_seq = uniqueComb_seq;
ResultNorm.matchedSimSet = matchedSimSet;
ResultNorm.PTD_Ref = PTD_Ref;
ResultNorm.Ref_Amp = Ref_Amp;
ResultNorm.post_win_ms = post_win_ms;

% Save exclusion rules
ResultNorm.ExcludeAmpRules = ExcludeAmpRules;

% Save reference estimate metadata
ResultNorm.RefVal_Sim = RefVal_Sim;
ResultNorm.RefVal_Seq = RefVal_Seq;
ResultNorm.RefMethod_Sim = RefMethod_Sim;
ResultNorm.RefMethod_Seq = RefMethod_Seq;

save(out_filename, 'ResultNorm');
fprintf('\n>>> Results Saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================

function count_val = get_spike_count(tr_ids, trig, sp_data, count_win, FS)

    if isempty(sp_data)
        count_val = NaN;
        return;
    end

    nTr = numel(tr_ids);
    total_spikes = 0;

    for k = 1:nTr
        tr = tr_ids(k);

        if tr > numel(trig)
            continue;
        end

        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;
        mask = tt >= count_win(1) & tt <= count_win(2);
        total_spikes = total_spikes + sum(mask);
    end

    if nTr > 0
        count_val = total_spikes / nTr;
    else
        count_val = NaN;
    end
end

function plot_shaded_error(x, y, se, col)

    if numel(x) < 2
        return;
    end

    x = x(:)';
    y = y(:)';
    se = se(:)';

    valid = ~isnan(x) & ~isnan(y) & ~isnan(se);
    x = x(valid);
    y = y(valid);
    se = se(valid);

    if isempty(x)
        return;
    end

    fill([x fliplr(x)], [y+se fliplr(y-se)], col, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)

    cd(folder);

    f = dir('*RespondingChannels.mat');
    if isempty(f)
        error('No Responding file in %s', folder);
    end
    R = load(f(1).name).Responding;

    f = dir('*sp_xia_SSD.mat');
    if isempty(f)
        f = dir('*sp_FirstPulse.mat');
    end
    if isempty(f)
        f = dir('*sp_xia.mat');
    end
    if isempty(f)
        error('No Spike file in %s', folder);
    end

    S_sp = load(f(1).name);

    if isfield(S_sp,'sp_corr')
        sp = S_sp.sp_corr;
    elseif isfield(S_sp,'sp_SSD')
        sp = S_sp.sp_SSD;
    elseif isfield(S_sp,'sp_seq')
        sp = S_sp.sp_seq;
    else
        sp = S_sp.sp_clipped;
    end

    if isempty(dir('*.trig.dat'))
        cleanTrig_sabquick;
    end

    trig = loadTrig(0);

    S = load(dir('*_exp_datafile_*.mat').name);

    QC.BadCh = [];
    QC.BadTrials = [];

    f_bc = dir('*.BadChannels.mat');
    if ~isempty(f_bc)
        tmp = load(f_bc(1).name);
        if isfield(tmp, 'BadCh_perSet')
            QC.BadCh = tmp.BadCh_perSet;
        elseif isfield(tmp, 'BadCh')
            QC.BadCh = tmp.BadCh;
        end
    end

    % Prefer SimSeqBadTrials if available, otherwise use any BadTrials file
    f_bt = dir('*.SimSeqBadTrials.mat');
    if isempty(f_bt)
        f_bt = dir('*.BadTrials.mat');
    end

    if ~isempty(f_bt)
        tmp = load(f_bt(1).name);
        if isfield(tmp, 'BadTrials')
            QC.BadTrials = tmp.BadTrials;
        end
    end
end

function tf = is_excluded_amp(mode_name, seq_set_idx, amp_val, ExcludeAmpRules)

    tf = false;

    if isempty(ExcludeAmpRules)
        return;
    end

    for r = 1:size(ExcludeAmpRules, 1)

        rule_mode = ExcludeAmpRules{r, 1};
        rule_set  = ExcludeAmpRules{r, 2};
        rule_amp  = ExcludeAmpRules{r, 3};

        if strcmpi(rule_mode, mode_name) && ...
                rule_set == seq_set_idx && ...
                abs(rule_amp - amp_val) < 0.001
            tf = true;
            return;
        end
    end
end

function [ref_val, method] = estimate_ref_response(Amps, curve, Ref_Amp)

    ref_val = NaN;
    method = "none";

    Amps = Amps(:)';
    curve = curve(:)';

    valid = ~isnan(curve);

    if ~any(valid)
        return;
    end

    % Case 1: real Ref_Amp exists and is valid
    ref_idx = find(abs(Amps - Ref_Amp) < 0.001, 1);

    if ~isempty(ref_idx) && valid(ref_idx)
        ref_val = curve(ref_idx);
        method = "real";
        return;
    end

    valid_amps = Amps(valid);
    valid_vals = curve(valid);

    % Sort by amplitude for interpolation
    [valid_amps, sort_idx] = sort(valid_amps);
    valid_vals = valid_vals(sort_idx);

    % Case 2: interpolation if Ref_Amp lies inside valid amplitude range
    if numel(valid_amps) >= 2 && ...
            Ref_Amp >= min(valid_amps) && ...
            Ref_Amp <= max(valid_amps)

        ref_val = interp1(valid_amps, valid_vals, Ref_Amp, 'linear');
        method = "interp";
        return;
    end

    % Case 3: nearest valid amplitude if interpolation is not possible
    [~, nearest_idx] = min(abs(valid_amps - Ref_Amp));
    ref_val = valid_vals(nearest_idx);
    method = "nearest";
end