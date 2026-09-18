%% ============================================================
%   Paired-Stimulation Linearity Analysis vs Single-Electrode Prediction
%
%   For each electrode pair (Set) x Amplitude:
%     R_A, R_B      = single-electrode-alone responses (from a SEPARATE
%                     single-electrode dataset), using the exact same
%                     net spike-count/trial metric ([2,40]ms window,
%                     baseline-corrected) as the paired-data side, so
%                     the comparison is apples-to-apples.
%     R_predicted   = R_A + R_B  (constant across ISI)
%     R_actual(ISI) = measured paired response at each ISI (from the
%                     same metric already used in your ISI tuning-curve
%                     script)
%     Ratio(ISI)    = R_actual(ISI) / R_predicted
%                     ~1 = linear summation; <1 = sub-additive;
%                     >1 = supra-additive
%
%   Both R_A/R_B and R_actual are pooled the same way: MEAN across the
%   SAME responding-channel population (taken from the ISI dataset's
%   Responding file for that Set), matching the pooling convention
%   already used in your ISI tuning-curve plot.
%
%   Two plots per Set (one curve per amplitude in each):
%     Plot A: Ratio vs ISI, with a reference line at 1.0
%     Plot B: R_actual(ISI) vs ISI, with each amplitude's R_predicted
%             overlaid as a matching-color dashed horizontal line
%
%   ASSUMPTIONS (flagging these -- check if results look off):
%     - The single-electrode dataset was recorded with the same
%       headstage/electrode setup, so Depth_s(Electrode_Type) and each
%       file's own E_MAP both refer to the same physical electrodes as
%       the paired (ISI) dataset -- electrode identity is matched via
%       E_MAP channel numbering, not by any other key.
%     - No bad-trial file currently exists for the single-electrode
%       dataset. The code checks for one anyway (same pattern as the
%       paired side) and will just use it automatically if you add one
%       later -- no code change needed.
%     - Responding-channel population is taken ONLY from the paired
%       (ISI) dataset's Responding file, for both R_A/R_B and R_actual,
%       per your instruction to keep pooling consistent.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder        = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';   % paired ISI dataset
single_elec_folder = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_Single1';   % single-electrode-alone dataset

Electrode_Type = 2; % 0:single shank rigid; 1:single shank flex; 2:four shank flex

% Select which ISIs (PTDs) to analyze. Leave [] to auto-detect ALL tested.
target_ISIs = [];

% Amplitudes to run the comparison on. Matched against the
% single-electrode dataset's own amplitudes (tolerance below).
target_Amps = [5 10]; % uA

% Fixed macro window + baseline window (SAME for paired and single-elec)
fixed_macro_win_ms = [2, 40];
baseline_win_ms    = [-50, -10];

% Tolerance (uA) for matching amplitudes between the two datasets
amp_match_tol = 1e-4;

% Bad-trial exclusion -- applies to BOTH datasets if a BadTrials file
% exists for them (single-electrode side currently has none, so this
% is a no-op there until one is added).
exclude_bad_trials = true;

% Artifact scrubber for the PAIRED data (same as your ISI script)
% Format: [Set, Amplitude (uA), ISI (ms)]
manual_artifacts = [
];

% PSTH kernel settings (only used internally by the shared spike-count
% helper function; the PSTH traces themselves aren't plotted here)
psth_win_ms = [-50 100];
FS = 30000;
bin_ms     = 1;
sigma_bins = 3;

bin_s = bin_ms/1000;
kernel_size = 2 * ceil(2*sigma_bins) + 1;
g_sym = gausswin(kernel_size); g_sym = g_sym / sum(g_sym);
edges_psth = psth_win_ms(1):bin_ms:psth_win_ms(2);

%% ================= RESTORE WORKING DIRECTORY ON EXIT =========
origDir = pwd;
dirCleanup = onCleanup(@() cd(origDir)); %#ok<NASGU>

%% =================== 1. LOAD PAIRED (ISI) DATA ====================
removeAppleDoubleFiles(data_folder);
removeAppleDoubleFiles(single_elec_folder);

[R, sp, trig, S, QC] = load_experiment_data(data_folder);

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

stimNames = Stim(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));
comb = zeros(nTr, simN);
for t = 1:nTr, rr = (t-1)*simN + (1:simN); v = idx_all(rr); v = v(v>0); comb(t,1:numel(v)) = v(:).'; end
[uniqueComb,~,combClass] = unique(comb,'rows','stable');
nSets = size(uniqueComb,1);

if isempty(target_Amps)
    target_Amps = Amps(:)';
    fprintf('Auto-detected target_Amps: %s uA\n', num2str(target_Amps));
end

if isempty(target_ISIs)
    target_ISIs = sort(unique([0; PTDs_ms(:)]))';
    fprintf('Auto-detected target_ISIs: %s ms\n', num2str(target_ISIs));
end

%% ================= 2. IDENTIFY POPULATION PER SET (paired) =============
d = Depth_s(Electrode_Type); nCh_Total = length(d);
resp_channels_per_set = cell(nSets, 1);

[~, ai_targets] = ismember(target_Amps, Amps);
if any(ai_targets == 0)
    error('One or more Target Amplitudes not found in the paired (ISI) dataset.');
end

fprintf('Analyzing Responding Channels (UNION MASK across %s uA) for ISIs %s ms:\n', ...
    num2str(target_Amps), num2str(target_ISIs));

for ss = 1:nSets
    local_resp_mask = false(nCh_Total, 1);

    if ss <= numel(R.set)
        for a_idx = 1:length(ai_targets)
            ai_targ = ai_targets(a_idx);

            if ai_targ <= numel(R.set(ss).amp)
                for pi=1:numel(R.set(ss).amp(ai_targ).ptd)
                    curr_ptd = R.set(ss).amp(ai_targ).ptd(pi).PTD_ms;

                    if ~any(abs(target_ISIs - curr_ptd) < 0.001), continue; end

                    this = R.set(ss).amp(ai_targ).ptd(pi).channel;
                    for ch=1:min(length(this),nCh_Total)
                        if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                            local_resp_mask(ch)=true;
                        end
                    end
                end
            end
        end
    end

    resp_channels_per_set{ss} = find(local_resp_mask);
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    fprintf('  Set %d (Ch:%s) -> %d Channels\n', ss, num2str(stimCh), length(resp_channels_per_set{ss}));
end

%% =================== 3. COMPUTE SPIKE COUNTS (paired) =================
SpikeCount_sim = nan(nCh_Total, length(target_Amps), nSets);
SpikeCount_seq = nan(nCh_Total, length(target_Amps), nSets, length(target_ISIs));

for ss = 1:nSets
    current_set_channels = resp_channels_per_set{ss};

    for ci = 1:length(current_set_channels)
        ch_idx = current_set_channels(ci);
        recCh = d(ch_idx);
        S_ch = sp{recCh};

        bad_trs = [];
        if exclude_bad_trials && ~isempty(QC.BadTrials) && ch_idx <= length(QC.BadTrials)
            bad_trs = QC.BadTrials{ch_idx};
        end
        if exclude_bad_trials && ~isempty(QC.BadCh) && ss <= length(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
            continue;
        end

        % --- SIMULTANEOUS (PTD = 0) ---
        p_sim_idx = find(target_ISIs == 0);
        if ~isempty(p_sim_idx)
            ptd_sim_raw = find(abs(PTDs_ms - 0) < 0.001);
            if ~isempty(ptd_sim_raw)
                for a_idx = 1:length(target_Amps)
                    ai_raw = ai_targets(a_idx);

                    tr_ids = find(combClass == ss & ampIdx == ai_raw & ptdIdx == ptd_sim_raw);
                    tr_ids = setdiff(tr_ids, bad_trs);

                    if isempty(tr_ids), continue; end

                    [count_val, ~] = get_spike_count_macro(tr_ids, trig, S_ch, ...
                        fixed_macro_win_ms, baseline_win_ms, edges_psth, g_sym, bin_s, FS);

                    SpikeCount_sim(ch_idx, a_idx, ss) = count_val;
                end
            end
        end

        % --- SEQUENTIAL (Target ISIs > 0) ---
        for p_idx = 1:length(target_ISIs)
            curr_isi = target_ISIs(p_idx);
            if curr_isi == 0, continue; end

            p_raw = find(abs(PTDs_ms - curr_isi) < 0.001);
            if isempty(p_raw), continue; end

            for a_idx = 1:length(target_Amps)
                ai_raw = ai_targets(a_idx);

                tr_ids = find(combClass==ss & ptdIdx==p_raw & ampIdx==ai_raw);
                tr_ids = setdiff(tr_ids, bad_trs);

                if isempty(tr_ids), continue; end

                [count_val, ~] = get_spike_count_macro(tr_ids, trig, S_ch, ...
                    fixed_macro_win_ms, baseline_win_ms, edges_psth, g_sym, bin_s, FS);

                SpikeCount_seq(ch_idx, a_idx, ss, p_idx) = count_val;
            end
        end
    end
end

%% =================== 3.5 ARTIFACT SCRUBBER (paired) ===================
if ~isempty(manual_artifacts)
    fprintf('\nApplying Artifact Scrubber...\n');
    for r = 1:size(manual_artifacts, 1)
        scrub_set = manual_artifacts(r, 1);
        scrub_amp = manual_artifacts(r, 2);
        scrub_isi = manual_artifacts(r, 3);

        a_scrub_idx = find(abs(target_Amps - scrub_amp) < 0.001);
        p_scrub_idx = find(abs(target_ISIs - scrub_isi) < 0.001);

        if isempty(a_scrub_idx) || isempty(p_scrub_idx), continue; end

        if scrub_isi == 0
            if scrub_set <= nSets
                SpikeCount_sim(:, a_scrub_idx, scrub_set) = NaN;
                fprintf('  Scrubbed: Set %d | %.1f uA | %d ms\n', scrub_set, scrub_amp, scrub_isi);
            end
        else
            if scrub_set <= nSets
                SpikeCount_seq(:, a_scrub_idx, scrub_set, p_scrub_idx) = NaN;
                fprintf('  Scrubbed: Set %d | %.1f uA | %d ms\n', scrub_set, scrub_amp, scrub_isi);
            end
        end
    end
    fprintf('\n');
end

%% =================== 4. LOAD SINGLE-ELECTRODE DATA =====================
[sp_se, trig_se, S_se, QC_se] = load_single_electrode_data(single_elec_folder);

Stim_se  = S_se.StimParams;
simN_se  = S_se.simultaneous_stim;
E_MAP_se = S_se.E_MAP;
if isfield(S_se, 'n_Trials'), nTr_se = S_se.n_Trials; else, nTr_se = (size(Stim_se,1)-1)/simN_se; end

amps_se_all = cell2mat(Stim_se(2:end,16));
trialAmps_se = amps_se_all(1:simN_se:end);
[Amps_se,~,~] = unique(trialAmps_se); Amps_se(Amps_se==-1) = 0;

stimNames_se = Stim_se(2:end,1);
[~, idx_all_se] = ismember(stimNames_se, E_MAP_se(2:end));
trialElectrode_se = zeros(nTr_se,1);
for t = 1:nTr_se
    rr = (t-1)*simN_se + (1:simN_se);
    v = idx_all_se(rr); v = v(v>0);
    if ~isempty(v)
        trialElectrode_se(t) = v(1); % single-electrode stim: one entry per trial
    end
end

fprintf('\nSingle-electrode dataset: %d trials, amplitudes tested: %s uA\n', ...
    nTr_se, num2str(Amps_se(:).'));

%% =================== 5. COMPUTE R_predicted / R_actual / RATIO =========
R_predicted = nan(nSets, length(target_Amps));
R_actual    = nan(nSets, length(target_Amps), length(target_ISIs));

for ss = 1:nSets

    current_set_channels = resp_channels_per_set{ss};
    stimCh = uniqueComb(ss, uniqueComb(ss,:)>0);

    if numel(stimCh) ~= 2
        fprintf('Set %d does not have exactly 2 stimulating electrodes (found %d) -- skipping linearity for this set.\n', ...
            ss, numel(stimCh));
        continue;
    end

    if isempty(current_set_channels)
        fprintf('Set %d has no responding channels -- skipping linearity for this set.\n', ss);
        continue;
    end

    elecA = stimCh(1);
    elecB = stimCh(2);

    for a_idx = 1:length(target_Amps)
        target_amp = target_Amps(a_idx);

        ai_se = find(abs(Amps_se - target_amp) < amp_match_tol, 1);
        if isempty(ai_se)
            fprintf('WARNING: Set %d | Amp %.1f uA: no matching amplitude found in single-electrode dataset. Skipping.\n', ...
                ss, target_amp);
            continue;
        end
        matched_amp_se = Amps_se(ai_se);

        R_A = compute_pooled_single_electrode_response(elecA, matched_amp_se, ...
            current_set_channels, d, sp_se, trig_se, trialElectrode_se, trialAmps_se, ...
            QC_se, exclude_bad_trials, fixed_macro_win_ms, baseline_win_ms, ...
            edges_psth, g_sym, bin_s, FS);

        R_B = compute_pooled_single_electrode_response(elecB, matched_amp_se, ...
            current_set_channels, d, sp_se, trig_se, trialElectrode_se, trialAmps_se, ...
            QC_se, exclude_bad_trials, fixed_macro_win_ms, baseline_win_ms, ...
            edges_psth, g_sym, bin_s, FS);

        if isnan(R_A) || isnan(R_B)
            fprintf('WARNING: Set %d | Amp %.1f uA: could not compute R_A/R_B (no valid trials/channels). Skipping.\n', ...
                ss, target_amp);
            continue;
        end

        R_predicted(ss, a_idx) = R_A + R_B;

        for p_idx = 1:length(target_ISIs)
            isi_val = target_ISIs(p_idx);

            if isi_val == 0
                data_set = squeeze(SpikeCount_sim(current_set_channels, a_idx, ss));
            else
                data_set = squeeze(SpikeCount_seq(current_set_channels, a_idx, ss, p_idx));
            end

            data_set = data_set(~isnan(data_set));
            if ~isempty(data_set)
                R_actual(ss, a_idx, p_idx) = mean(data_set);
            end
        end
    end
end

Ratio = R_actual ./ R_predicted; % broadcasts R_predicted (nSets x nAmp) across the ISI dimension

%% ===================== 6. PLOTS =========================================
amp_colors = lines(length(target_Amps));

for ss = 1:nSets

    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    if numel(stimCh) ~= 2, continue; end
    if all(isnan(R_predicted(ss,:))), continue; end

    %% ---- Plot A: Ratio vs ISI ----
    figure('Color','w', 'Position',[100 100 800 600]);
    hold on;
    yline(1.0, 'k--', 'LineWidth', 1); % linear-summation reference

    for a_idx = 1:length(target_Amps)
        y = squeeze(Ratio(ss, a_idx, :));
        valid_idx = ~isnan(y);
        if ~any(valid_idx), continue; end

        col = amp_colors(a_idx, :);
        lbl = sprintf('%.1f uA', target_Amps(a_idx));

        plot(target_ISIs(valid_idx), y(valid_idx), '-o', 'Color', col, ...
            'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', lbl);
    end

    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold', 'FontSize', 12);
    ylabel('Ratio (Actual / Predicted Linear Sum)', 'FontWeight','bold', 'FontSize', 12);
    title(sprintf('Linearity Ratio -- Set %d (Stim Channels: %s)', ss, num2str(stimCh)), ...
        'FontWeight','bold', 'FontSize', 14);
    xticks(sort(target_ISIs));
    box off;
    lgd = legend('Location','best','Box','off'); title(lgd, 'Amplitudes');

    %% ---- Plot B: Raw actual vs ISI, with predicted lines overlaid ----
    figure('Color','w', 'Position',[950 100 800 600]);
    hold on;

    for a_idx = 1:length(target_Amps)
        y = squeeze(R_actual(ss, a_idx, :));
        valid_idx = ~isnan(y);

        col = amp_colors(a_idx, :);
        lbl = sprintf('%.1f uA', target_Amps(a_idx));

        if any(valid_idx)
            plot(target_ISIs(valid_idx), y(valid_idx), '-o', 'Color', col, ...
                'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', lbl);
        end

        if ~isnan(R_predicted(ss, a_idx))
            yline(R_predicted(ss, a_idx), '--', 'Color', col, 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%.1f uA predicted (A+B)', target_Amps(a_idx)));
        end
    end

    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold', 'FontSize', 12);
    ylabel('Net Mean Spikes', 'FontWeight','bold', 'FontSize', 12);
    title(sprintf('Actual vs Predicted Linear Sum -- Set %d (Stim Channels: %s)', ss, num2str(stimCh)), ...
        'FontWeight','bold', 'FontSize', 14);
    xticks(sort(target_ISIs));
    box off;
    lgd = legend('Location','best','Box','off'); title(lgd, 'Amplitudes');
end

%% ============================================================
%   7. SAVE RESULTS
% ============================================================
save_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_Analysis/DX023/';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
parts = split(data_folder, filesep); exp_id = parts{end};

if length(target_Amps) > 4
    amp_str = 'AllAmps';
else
    amp_str = strjoin(string(target_Amps), '_');
end

out_filename = fullfile(save_dir, ['Result_Linearity_' char(amp_str) 'uA_' exp_id '.mat']);

ResultLin = struct();
ResultLin.Metadata.Created = datestr(now);
ResultLin.Metadata.Metric = 'Net Mean Spike Count per Trial, mean-pooled across responding channels';
ResultLin.Metadata.FixedMacroWin = fixed_macro_win_ms;
ResultLin.Metadata.BaselineWin = baseline_win_ms;
ResultLin.Metadata.TargetAmps = target_Amps;
ResultLin.Metadata.TargetISIs = target_ISIs;
ResultLin.Metadata.BadTrialsExcluded = exclude_bad_trials;
ResultLin.Metadata.AmpMatchTolerance_uA = amp_match_tol;
ResultLin.Metadata.PairedDataFolder = data_folder;
ResultLin.Metadata.SingleElectrodeFolder = single_elec_folder;
ResultLin.Metadata.Stimulation_Sets = uniqueComb;
ResultLin.Metadata.Responding_Channels = resp_channels_per_set;

ResultLin.R_predicted = R_predicted;   % [nSets x nAmp]
ResultLin.R_actual    = R_actual;      % [nSets x nAmp x nISI]
ResultLin.Ratio       = Ratio;         % [nSets x nAmp x nISI]

save(out_filename, 'ResultLin');
fprintf('\n>>> Linearity results saved to: %s\n', out_filename);

%% ==================== HELPER FUNCTIONS =========================
function [count_val, psth_trace] = get_spike_count_macro(tr_ids, trig, sp_data, ...
    count_win, base_win, psth_edges, g_sym, bin_s, FS)

    nTr = numel(tr_ids);
    all_psth_spikes = [];
    total_net_spikes_in_window = 0;

    dur_evoked = count_win(2) - count_win(1);
    dur_base   = base_win(2) - base_win(1);

    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;

        mask_count = (tt >= count_win(1) & tt <= count_win(2));
        evoked_count = sum(mask_count);

        mask_base = (tt >= base_win(1) & tt <= base_win(2));
        base_count = sum(mask_base);

        net_count = max(0, evoked_count - (base_count * (dur_evoked / dur_base)));
        total_net_spikes_in_window = total_net_spikes_in_window + net_count;

        all_psth_spikes = [all_psth_spikes; tt(tt >= psth_edges(1) & tt <= psth_edges(end))]; %#ok<AGROW>
    end

    if nTr > 0
        count_val = total_net_spikes_in_window / nTr;
    else
        count_val = NaN;
    end

    h_psth = histcounts(all_psth_spikes, psth_edges);
    rate_psth = h_psth / (nTr * bin_s);
    psth_trace = conv(rate_psth, g_sym, 'same');
end

function [R, sp, trig, S, QC] = load_experiment_data(folder)
    cd(folder);
    f = dir('*_MultiISI_RespondingChannels.mat'); if isempty(f), error('No Responding file in %s', folder); end
    R = load(f(1).name).Responding;

    f = dir('*sp_xia_SSD.mat'); if isempty(f), f=dir('*sp_xia.mat'); end
    if isempty(f), error('No Spike file in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end

    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end; trig = loadTrig(0);
    S = load(dir('*_exp_datafile_*.mat').name);

    QC.BadCh = []; QC.BadTrials = [];
    f_bc = dir('*.MultiISIsBadChannels.mat'); if ~isempty(f_bc), tmp = load(f_bc(1).name); QC.BadCh = tmp.BadCh_perSet; end
    f_bt = dir('*.MultiISIsBadTrials.mat'); if ~isempty(f_bt), tmp = load(f_bt(1).name); QC.BadTrials = tmp.BadTrials; end
end

function [sp, trig, S, QC] = load_single_electrode_data(folder)
% Same style as load_experiment_data, but for a single-electrode-only
% dataset that has no Responding/BadChannels file of its own. Checks for
% a BadTrials file anyway (several common naming patterns) so it starts
% working automatically if one is added later -- no code change needed.

    cd(folder);

    f = dir('*sp_xia_SSD.mat'); if isempty(f), f = dir('*sp_xia.mat'); end
    if isempty(f), error('No spike file found in %s', folder); end
    S_sp = load(f(1).name);
    if isfield(S_sp,'sp_corr'), sp = S_sp.sp_corr; elseif isfield(S_sp,'sp_SSD'), sp = S_sp.sp_SSD; else, sp = S_sp.sp_in; end

    if isempty(dir('*.trig.dat')), cleanTrig_sabquick; end
    trig = loadTrig(0);

    expFile = dir('*_exp_datafile_*.mat');
    if isempty(expFile), error('No experiment datafile found in %s', folder); end
    S = load(expFile(1).name, 'StimParams','simultaneous_stim','E_MAP','n_Trials');

    QC.BadTrials = [];
    f_bt = dir('*.MultiISIsBadTrials.mat');
    if isempty(f_bt), f_bt = dir('*.BadTrials.mat'); end
    if ~isempty(f_bt)
        tmp = load(f_bt(1).name);
        if isfield(tmp,'BadTrials'), QC.BadTrials = tmp.BadTrials; end
    end
end

function R = compute_pooled_single_electrode_response(targetElectrode, targetAmp, ...
    channels_list, d, sp, trig, trialElectrode, trialAmps, QC, exclude_bad_trials, ...
    count_win, base_win, psth_edges, g_sym, bin_s, FS)
% Pools the net spike-count/trial metric across channels_list for trials
% where targetElectrode was stimulated alone at targetAmp -- same
% pooling convention (mean across channels) as the paired-data side.

    tr_ids_base = find(trialElectrode == targetElectrode & abs(trialAmps - targetAmp) < 1e-4);

    per_channel_vals = nan(numel(channels_list), 1);

    for ci = 1:numel(channels_list)
        ch_idx = channels_list(ci);

        if ch_idx > numel(d), continue; end
        recCh = d(ch_idx);
        if recCh < 1 || recCh > numel(sp) || isempty(sp{recCh}), continue; end

        bad_trs = [];
        if exclude_bad_trials && ~isempty(QC.BadTrials) && ch_idx <= numel(QC.BadTrials)
            bad_trs = QC.BadTrials{ch_idx};
        end

        tr_ids = setdiff(tr_ids_base, bad_trs);
        if isempty(tr_ids), continue; end

        [count_val, ~] = get_spike_count_macro(tr_ids, trig, sp{recCh}, ...
            count_win, base_win, psth_edges, g_sym, bin_s, FS);

        per_channel_vals(ci) = count_val;
    end

    per_channel_vals = per_channel_vals(~isnan(per_channel_vals));
    if isempty(per_channel_vals)
        R = NaN;
    else
        R = mean(per_channel_vals);
    end
end

function removeAppleDoubleFiles(folderPath)
% Deletes hidden "._*" AppleDouble resource-fork files that macOS
% creates when copying to non-native filesystems (e.g. SMB shares).
% Safe to delete -- they hold no scientific data, only macOS metadata.
% Without this, dir('*.ext')-based loaders (loadTrig, load_experiment_data,
% etc.) can match more than one file and throw cryptic errors.
staleFiles = dir(fullfile(folderPath, '._*'));
for k = 1:numel(staleFiles)
    try
        delete(fullfile(staleFiles(k).folder, staleFiles(k).name));
    catch ME
        warning('Could not delete %s: %s', staleFiles(k).name, ME.message);
    end
end
end