%% ============================================================
%   EXPORT DATA FOR COLLABORATOR (Paired ISI + Single-Electrode)
%
%   Produces a clean, shareable export -- NOT the raw .dat files, and
%   NOT your internal QC bookkeeping:
%     - Only RESPONDING channels are included (silently -- no "bad
%       channel" list is exported, non-responding channels are simply
%       absent).
%     - Bad trials are excluded (silently -- no "bad trial" list or
%       original trial count is exported).
%     - Trial numbers are renumbered fresh (1..N) per channel/condition,
%       so the original session's trial count/order isn't inferable.
%
%   Two tables are exported, as both CSV and .mat:
%     1) SpikeCounts -- one row per (condition x channel x trial):
%        the baseline-corrected net spike count in a fixed window
%        ([2,40]ms by default), a ready-to-use summary metric.
%     2) SpikeTimes -- one row per (condition x channel x trial x spike):
%        individual spike times relative to the trigger, within a wider
%        window (matches your raster/PSTH window by default), so your
%        colleague can build rasters/PSTHs or recompute counts with a
%        different window himself.
%
%   Spike WAVEFORMS are intentionally NOT exported (large, not needed
%   for rasters/PSTHs or a linearity model).
%
%   An auto-generated README.txt documents every column and setting
%   used, filled in from the actual values below -- so the docs can
%   never drift out of sync with the actual export.
%
%   ASSUMPTION (single-electrode responding channels): an electrode can
%   participate in more than one Set (pairing). For the single-electrode
%   export, "responding channels for electrode E" = the UNION of
%   responding channels across every Set that includes E in the paired
%   dataset's Responding file. Flagging this since there's no
%   single-electrode-specific Responding file to draw from directly.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder        = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';   % paired ISI dataset
single_elec_folder = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_Single1';   % single-electrode dataset

Electrode_Type = 2;

% Window for the precomputed baseline-corrected spike-count column
fixed_macro_win_ms = [2, 40];
baseline_win_ms    = [-50, -10];

% Window for the exported raw spike times (default matches typical
% raster/PSTH window -- widen/narrow as needed)
spike_time_win_ms = [-50, 100];

FS = 30000;

% Only export these conditions. Leave empty ({}) to export ALL
% conditions that have at least one responding channel.
%   Paired:  {Set, Amp(uA), ISI(ms)}
%   Single:  {Electrode, Amp(uA)}
IncludePairedConditions = {
    % 1, 5, 0;
    % 1, 5, 5;
    2,   10,      0;
    2,   10,      3;
    2,   10,      4;
    2,   10,      5;
    2,   10,      6;
    2,   10,      7;
    2,   10,      8;
    2,   10,      9;
    2,   10,      10;
    2,   10,      11;
    2,   10,      12;
    2,   10,      13;
    2,   10,      14;
    2,   10,      15;
    2,   10,      17;
    2,   10,      20;
};
IncludeSingleConditions = {
    % 5, 5;
};

output_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Luke_Data/DX023_v2';

%% ================= RESTORE WORKING DIRECTORY ON EXIT =========
origDir = pwd;
dirCleanup = onCleanup(@() cd(origDir)); %#ok<NASGU>

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% =================== 1. LOAD PAIRED (ISI) DATA ====================
removeAppleDoubleFiles(data_folder);
removeAppleDoubleFiles(single_elec_folder);

[R, sp, trig, S, QC] = load_experiment_data(data_folder);

[~, AnimalID] = fileparts(fileparts(data_folder));
[~, SessionName_paired] = fileparts(data_folder);
[~, SessionName_single] = fileparts(single_elec_folder);

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

d = Depth_s(Electrode_Type); nCh_Total = length(d);

%% ============== 2. RESPONDING CHANNELS PER SET (paired) ==============
resp_channels_per_set = cell(nSets, 1);

for ss = 1:nSets
    local_resp_mask = false(nCh_Total, 1);
    if ss <= numel(R.set)
        for ai = 1:numel(R.set(ss).amp)
            for pi = 1:numel(R.set(ss).amp(ai).ptd)
                this = R.set(ss).amp(ai).ptd(pi).channel;
                for ch = 1:min(length(this), nCh_Total)
                    if isfield(this(ch),'is_responsive') && this(ch).is_responsive
                        local_resp_mask(ch) = true;
                    end
                end
            end
        end
    end
    resp_channels_per_set{ss} = find(local_resp_mask);
end

%% =================== 3. EXPORT PAIRED DATA =====================
Sc_DatasetType = {}; Sc_Animal = {}; Sc_Session = {}; Sc_ElecA = []; Sc_ElecB = [];
Sc_Amp = []; Sc_ISI = []; Sc_Chan = []; Sc_Trial = []; Sc_Count = [];

St_DatasetType = {}; St_Animal = {}; St_Session = {}; St_ElecA = []; St_ElecB = [];
St_Amp = []; St_ISI = []; St_Chan = []; St_Trial = []; St_SpikeTime = [];

fprintf('\nExporting PAIRED (ISI) data...\n');

for ss = 1:nSets
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    if numel(stimCh) < 1, continue; end
    elecA = stimCh(1);
    elecB = NaN; if numel(stimCh) >= 2, elecB = stimCh(2); end

    channels_this_set = resp_channels_per_set{ss};
    if isempty(channels_this_set), continue; end

    for ai = 1:numel(Amps)
        target_amp = Amps(ai);

        for pi = 1:numel(PTDs_ms)
            target_isi = PTDs_ms(pi);

            if ~isempty(IncludePairedConditions)
                keepRow = false;
                for r = 1:size(IncludePairedConditions,1)
                    if IncludePairedConditions{r,1} == ss && ...
                            abs(IncludePairedConditions{r,2} - target_amp) < 1e-4 && ...
                            abs(IncludePairedConditions{r,3} - target_isi) < 1e-4
                        keepRow = true; break;
                    end
                end
                if ~keepRow, continue; end
            end

            tr_ids_base = find(combClass == ss & ampIdx == ai & ptdIdx == pi);
            if isempty(tr_ids_base), continue; end

            for ci = 1:numel(channels_this_set)
                ch_idx = channels_this_set(ci);
                recCh = d(ch_idx);
                if recCh < 1 || recCh > numel(sp) || isempty(sp{recCh}), continue; end

                bad_trs = [];
                if ~isempty(QC.BadTrials) && ch_idx <= numel(QC.BadTrials)
                    bad_trs = QC.BadTrials{ch_idx};
                end
                if ~isempty(QC.BadCh) && ss <= numel(QC.BadCh) && ismember(ch_idx, QC.BadCh{ss})
                    continue;
                end

                tr_ids = setdiff(tr_ids_base, bad_trs);
                if isempty(tr_ids), continue; end

                [countsPerTrial, spikeTimesPerTrial] = get_counts_and_times_per_trial( ...
                    tr_ids, trig, sp{recCh}, fixed_macro_win_ms, baseline_win_ms, ...
                    spike_time_win_ms, FS);

                nValid = numel(tr_ids);
                freshTrialIdx = (1:nValid)';

                Sc_DatasetType = [Sc_DatasetType; repmat({'Paired'}, nValid, 1)]; %#ok<AGROW>
                Sc_Animal      = [Sc_Animal; repmat({AnimalID}, nValid, 1)]; %#ok<AGROW>
                Sc_Session     = [Sc_Session; repmat({SessionName_paired}, nValid, 1)]; %#ok<AGROW>
                Sc_ElecA       = [Sc_ElecA; repmat(elecA, nValid, 1)]; %#ok<AGROW>
                Sc_ElecB       = [Sc_ElecB; repmat(elecB, nValid, 1)]; %#ok<AGROW>
                Sc_Amp         = [Sc_Amp; repmat(target_amp, nValid, 1)]; %#ok<AGROW>
                Sc_ISI         = [Sc_ISI; repmat(target_isi, nValid, 1)]; %#ok<AGROW>
                Sc_Chan        = [Sc_Chan; repmat(ch_idx, nValid, 1)]; %#ok<AGROW>
                Sc_Trial       = [Sc_Trial; freshTrialIdx]; %#ok<AGROW>
                Sc_Count       = [Sc_Count; countsPerTrial]; %#ok<AGROW>

                for k = 1:nValid
                    nSpk = numel(spikeTimesPerTrial{k});
                    if nSpk == 0, continue; end

                    St_DatasetType = [St_DatasetType; repmat({'Paired'}, nSpk, 1)]; %#ok<AGROW>
                    St_Animal      = [St_Animal; repmat({AnimalID}, nSpk, 1)]; %#ok<AGROW>
                    St_Session     = [St_Session; repmat({SessionName_paired}, nSpk, 1)]; %#ok<AGROW>
                    St_ElecA       = [St_ElecA; repmat(elecA, nSpk, 1)]; %#ok<AGROW>
                    St_ElecB       = [St_ElecB; repmat(elecB, nSpk, 1)]; %#ok<AGROW>
                    St_Amp         = [St_Amp; repmat(target_amp, nSpk, 1)]; %#ok<AGROW>
                    St_ISI         = [St_ISI; repmat(target_isi, nSpk, 1)]; %#ok<AGROW>
                    St_Chan        = [St_Chan; repmat(ch_idx, nSpk, 1)]; %#ok<AGROW>
                    St_Trial       = [St_Trial; repmat(freshTrialIdx(k), nSpk, 1)]; %#ok<AGROW>
                    St_SpikeTime   = [St_SpikeTime; spikeTimesPerTrial{k}]; %#ok<AGROW>
                end
            end
        end
    end
    fprintf('  Set %d (Ch:%s) done.\n', ss, num2str(stimCh));
end

%% =================== 4. LOAD + EXPORT SINGLE-ELECTRODE DATA =====================
fprintf('\nExporting SINGLE-ELECTRODE data...\n');

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
    if ~isempty(v), trialElectrode_se(t) = v(1); end
end

allElectrodes_se = unique(trialElectrode_se(trialElectrode_se > 0));

% Responding channels per electrode = union across every paired Set that
% includes that electrode (see ASSUMPTION note at top of file).
resp_channels_per_electrode = containers.Map('KeyType','double','ValueType','any');
for e = allElectrodes_se(:)'
    unionMask = false(nCh_Total,1);
    for ss = 1:nSets
        stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
        if ismember(e, stimCh)
            chs = resp_channels_per_set{ss};
            unionMask(chs) = true;
        end
    end
    resp_channels_per_electrode(e) = find(unionMask);
end

for e = allElectrodes_se(:)'
    channels_this_elec = resp_channels_per_electrode(e);
    if isempty(channels_this_elec), continue; end

    for ai = 1:numel(Amps_se)
        target_amp = Amps_se(ai);

        if ~isempty(IncludeSingleConditions)
            keepRow = false;
            for r = 1:size(IncludeSingleConditions,1)
                if IncludeSingleConditions{r,1} == e && ...
                        abs(IncludeSingleConditions{r,2} - target_amp) < 1e-4
                    keepRow = true; break;
                end
            end
            if ~keepRow, continue; end
        end

        tr_ids_base = find(trialElectrode_se == e & abs(trialAmps_se - target_amp) < 1e-4);
        if isempty(tr_ids_base), continue; end

        for ci = 1:numel(channels_this_elec)
            ch_idx = channels_this_elec(ci);
            if ch_idx > numel(d), continue; end
            recCh = d(ch_idx);
            if recCh < 1 || recCh > numel(sp_se) || isempty(sp_se{recCh}), continue; end

            bad_trs = [];
            if ~isempty(QC_se.BadTrials) && ch_idx <= numel(QC_se.BadTrials)
                bad_trs = QC_se.BadTrials{ch_idx};
            end

            tr_ids = setdiff(tr_ids_base, bad_trs);
            if isempty(tr_ids), continue; end

            [countsPerTrial, spikeTimesPerTrial] = get_counts_and_times_per_trial( ...
                tr_ids, trig_se, sp_se{recCh}, fixed_macro_win_ms, baseline_win_ms, ...
                spike_time_win_ms, FS);

            nValid = numel(tr_ids);
            freshTrialIdx = (1:nValid)';

            Sc_DatasetType = [Sc_DatasetType; repmat({'Single'}, nValid, 1)]; %#ok<AGROW>
            Sc_Animal      = [Sc_Animal; repmat({AnimalID}, nValid, 1)]; %#ok<AGROW>
            Sc_Session     = [Sc_Session; repmat({SessionName_single}, nValid, 1)]; %#ok<AGROW>
            Sc_ElecA       = [Sc_ElecA; repmat(e, nValid, 1)]; %#ok<AGROW>
            Sc_ElecB       = [Sc_ElecB; repmat(NaN, nValid, 1)]; %#ok<AGROW>
            Sc_Amp         = [Sc_Amp; repmat(target_amp, nValid, 1)]; %#ok<AGROW>
            Sc_ISI         = [Sc_ISI; repmat(NaN, nValid, 1)]; %#ok<AGROW>
            Sc_Chan        = [Sc_Chan; repmat(ch_idx, nValid, 1)]; %#ok<AGROW>
            Sc_Trial       = [Sc_Trial; freshTrialIdx]; %#ok<AGROW>
            Sc_Count       = [Sc_Count; countsPerTrial]; %#ok<AGROW>

            for k = 1:nValid
                nSpk = numel(spikeTimesPerTrial{k});
                if nSpk == 0, continue; end

                St_DatasetType = [St_DatasetType; repmat({'Single'}, nSpk, 1)]; %#ok<AGROW>
                St_Animal      = [St_Animal; repmat({AnimalID}, nSpk, 1)]; %#ok<AGROW>
                St_Session     = [St_Session; repmat({SessionName_single}, nSpk, 1)]; %#ok<AGROW>
                St_ElecA       = [St_ElecA; repmat(e, nSpk, 1)]; %#ok<AGROW>
                St_ElecB       = [St_ElecB; repmat(NaN, nSpk, 1)]; %#ok<AGROW>
                St_Amp         = [St_Amp; repmat(target_amp, nSpk, 1)]; %#ok<AGROW>
                St_ISI         = [St_ISI; repmat(NaN, nSpk, 1)]; %#ok<AGROW>
                St_Chan        = [St_Chan; repmat(ch_idx, nSpk, 1)]; %#ok<AGROW>
                St_Trial       = [St_Trial; repmat(freshTrialIdx(k), nSpk, 1)]; %#ok<AGROW>
                St_SpikeTime   = [St_SpikeTime; spikeTimesPerTrial{k}]; %#ok<AGROW>
            end
        end
    end
    fprintf('  Electrode %d done.\n', e);
end

%% =================== 5. BUILD + SAVE TABLES =====================
SpikeCounts = table(Sc_DatasetType, Sc_Animal, Sc_Session, Sc_ElecA, Sc_ElecB, ...
    Sc_Amp, Sc_ISI, Sc_Chan, Sc_Trial, Sc_Count, 'VariableNames', ...
    {'DatasetType','AnimalID','SessionName','ElectrodeA','ElectrodeB', ...
     'Amplitude_uA','ISI_ms','Channel','TrialIndex','SpikeCount'});

SpikeTimes = table(St_DatasetType, St_Animal, St_Session, St_ElecA, St_ElecB, ...
    St_Amp, St_ISI, St_Chan, St_Trial, St_SpikeTime, 'VariableNames', ...
    {'DatasetType','AnimalID','SessionName','ElectrodeA','ElectrodeB', ...
     'Amplitude_uA','ISI_ms','Channel','TrialIndex','SpikeTime_ms'});

csv_counts_path = fullfile(output_dir, [AnimalID '_SpikeCounts.csv']);
csv_times_path  = fullfile(output_dir, [AnimalID '_SpikeTimes.csv']);
mat_path        = fullfile(output_dir, [AnimalID '_ExportForCollaborator.mat']);
readme_path     = fullfile(output_dir, 'README.txt');

writetable(SpikeCounts, csv_counts_path);
writetable(SpikeTimes, csv_times_path);

ExportMetadata = struct();
ExportMetadata.AnimalID = AnimalID;
ExportMetadata.PairedSessionName = SessionName_paired;
ExportMetadata.SingleElectrodeSessionName = SessionName_single;
ExportMetadata.Created = datestr(now);
ExportMetadata.FixedMacroWindow_ms = fixed_macro_win_ms;
ExportMetadata.BaselineWindow_ms = baseline_win_ms;
ExportMetadata.SpikeTimeWindow_ms = spike_time_win_ms;
ExportMetadata.SamplingRate_Hz = FS;
ExportMetadata.Note = ['Only responding channels and non-excluded trials are ' ...
    'included. TrialIndex is a fresh 1..N renumbering, independent per ' ...
    'channel and condition -- it does NOT reflect original recording ' ...
    'trial order or count, and may not align across channels within ' ...
    'the same condition.'];

save(mat_path, 'SpikeCounts', 'SpikeTimes', 'ExportMetadata');

%% =================== 6. WRITE README =====================
fid = fopen(readme_path, 'w');
fprintf(fid, 'DATA EXPORT README\n');
fprintf(fid, '===================\n\n');
fprintf(fid, 'Animal: %s\n', AnimalID);
fprintf(fid, 'Paired (ISI) session: %s\n', SessionName_paired);
fprintf(fid, 'Single-electrode session: %s\n', SessionName_single);
fprintf(fid, 'Created: %s\n\n', datestr(now));

fprintf(fid, 'FILES\n-----\n');
fprintf(fid, '%s_SpikeCounts.csv / .mat (SpikeCounts table)\n', AnimalID);
fprintf(fid, '  One row per (condition x channel x trial): a single baseline-\n');
fprintf(fid, '  corrected net spike count per trial, in a fixed [%g, %g] ms\n', ...
    fixed_macro_win_ms(1), fixed_macro_win_ms(2));
fprintf(fid, '  window relative to the stimulation trigger, baseline-corrected\n');
fprintf(fid, '  against a [%g, %g] ms pre-trigger window.\n\n', ...
    baseline_win_ms(1), baseline_win_ms(2));

fprintf(fid, '%s_SpikeTimes.csv / .mat (SpikeTimes table)\n', AnimalID);
fprintf(fid, '  One row per individual spike, with its time (ms) relative to\n');
fprintf(fid, '  the stimulation trigger, within a [%g, %g] ms window. Use this\n', ...
    spike_time_win_ms(1), spike_time_win_ms(2));
fprintf(fid, '  to build rasters/PSTHs, or to recompute spike counts yourself\n');
fprintf(fid, '  with a different window. Spike WAVEFORMS are not included.\n\n');

fprintf(fid, 'COLUMNS (both tables)\n----------------------\n');
fprintf(fid, 'DatasetType   "Paired" (two electrodes stimulated together) or\n');
fprintf(fid, '              "Single" (one electrode stimulated alone)\n');
fprintf(fid, 'AnimalID      Animal identifier\n');
fprintf(fid, 'SessionName   Recording session identifier\n');
fprintf(fid, 'ElectrodeA    Stimulating electrode # (or the single electrode,\n');
fprintf(fid, '              for DatasetType=="Single")\n');
fprintf(fid, 'ElectrodeB    Second stimulating electrode # (NaN for "Single" rows)\n');
fprintf(fid, 'Amplitude_uA  Stimulation current amplitude, microamps\n');
fprintf(fid, 'ISI_ms        Inter-stimulus interval between ElectrodeA and\n');
fprintf(fid, '              ElectrodeB pulses, milliseconds. 0 = simultaneous.\n');
fprintf(fid, '              NaN for "Single" rows (not applicable).\n');
fprintf(fid, 'Channel       Recording channel number\n');
fprintf(fid, 'TrialIndex    A fresh 1..N renumbering of included trials, specific\n');
fprintf(fid, '              to that (Channel, condition) combination. This is NOT\n');
fprintf(fid, '              the original recording trial number, and does not\n');
fprintf(fid, '              necessarily line up across different channels within\n');
fprintf(fid, '              the same condition.\n');
fprintf(fid, 'SpikeCount    (SpikeCounts table only) net baseline-corrected spike\n');
fprintf(fid, '              count for that trial, in the fixed window above.\n');
fprintf(fid, 'SpikeTime_ms  (SpikeTimes table only) one spike''s time relative to\n');
fprintf(fid, '              the trigger, milliseconds. Negative = before trigger.\n\n');

fprintf(fid, 'WHAT''S ALREADY BEEN FILTERED OUT\n---------------------------------\n');
fprintf(fid, '- Only channels that showed a significant evoked response\n');
fprintf(fid, '  ("responding channels") for a given condition are included.\n');
fprintf(fid, '- Trials flagged as noisy/artifactual during quality control have\n');
fprintf(fid, '  already been excluded.\n');
fprintf(fid, '- Spike waveforms are not included (only spike times and counts).\n\n');

fprintf(fid, 'SAMPLING RATE: %d Hz\n', FS);
fclose(fid);

fprintf('\n>>> Export complete.\n');
fprintf('    %s\n    %s\n    %s\n    %s\n', csv_counts_path, csv_times_path, mat_path, readme_path);

%% ==================== HELPER FUNCTIONS =========================
function [countsPerTrial, spikeTimesPerTrial] = get_counts_and_times_per_trial( ...
    tr_ids, trig, sp_data, count_win, base_win, spike_win, FS)

    nTr = numel(tr_ids);
    countsPerTrial = nan(nTr,1);
    spikeTimesPerTrial = cell(nTr,1);

    dur_evoked = count_win(2) - count_win(1);
    dur_base   = base_win(2) - base_win(1);

    for k = 1:nTr
        tr = tr_ids(k);
        t0 = trig(tr)/FS*1000;
        tt = sp_data(:,1) - t0;

        evoked_count = sum(tt >= count_win(1) & tt <= count_win(2));
        base_count   = sum(tt >= base_win(1) & tt <= base_win(2));
        net_count = max(0, evoked_count - (base_count * (dur_evoked / dur_base)));

        countsPerTrial(k) = net_count;
        spikeTimesPerTrial{k} = tt(tt >= spike_win(1) & tt <= spike_win(2));
    end
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

function removeAppleDoubleFiles(folderPath)
staleFiles = dir(fullfile(folderPath, '._*'));
for k = 1:numel(staleFiles)
    try
        delete(fullfile(staleFiles(k).folder, staleFiles(k).name));
    catch ME
        warning('Could not delete %s: %s', staleFiles(k).name, ME.message);
    end
end
end