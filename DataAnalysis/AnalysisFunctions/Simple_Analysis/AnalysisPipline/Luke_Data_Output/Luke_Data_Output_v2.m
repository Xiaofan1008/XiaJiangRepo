%% ============================================================
%   EXPORT DATA FOR Modeling (Paired ISI + Single-Electrode)
%
%   Produces a clean, shareable export -- NOT the raw .dat files, and
%   NOT your internal QC bookkeeping:
%     - Only RESPONDING channels are included (silently -- non-
%       responding channels are simply absent, no list is exported).
%     - Bad trials are excluded (silently -- no list or original trial
%       count is exported).
%     - Trial position in each array IS the trial index (1..N), fresh
%       and independent per channel/condition -- the original session's
%       trial count/order isn't inferable.
%     - Spike WAVEFORMS are not exported (not needed for rasters/PSTHs
%       or a linearity model, and take a lot of space).
%
%   SELECTION: list the (Set, Amp, ISI) rows you want in
%   IncludePairedConditions below. Leave it empty to export everything.
%   The matching SINGLE-ELECTRODE conditions (both component electrodes
%   of each selected Set, at the same amplitude) are derived
%   automatically from uniqueComb -- you only maintain ONE list.
%
%   OUTPUT:
%     <name>_Export.mat   -- Paired / Single,
%       nested MATLAB structs: Paired is Set -> Amp -> ISI -> Channel;
%       Single is Set -> Amp -> Electrode -> Channel (scoped per Set, see
%       below) -- each leaf holding SpikeCount and SpikeTimes_ms per trial
%       (array/cell position = trial index; no separate trial-index
%       column needed).
%     Paired_SpikeCounts.csv, Paired_SpikeTimes.csv,
%     Single_SpikeCounts.csv, Single_SpikeTimes.csv -- flat, portable
%       equivalents of the same data (CSV can't nest, so these stay
%       tabular with an explicit TrialIndex column).
%     README.txt -- column/setting documentation, auto-filled from the
%       actual values below.
%
%   SINGLE-ELECTRODE DATA IS SCOPED PER SET: for each exported Set (a
%   pairing of two electrodes), the Single-electrode data for its two
%   component electrodes uses the SAME responding-channel population as
%   that Set's Paired data (resp_channels_per_set{ss}) -- not a global
%   per-electrode union. If the same physical electrode appears in more
%   than one Set (e.g. a "seq order swapped" condition), its data is
%   exported once per Set, each time scoped to that Set's own channels.
%   This guarantees Paired(k) and Single(k) for the same Set can be
%   combined directly (e.g. R_predicted = R_A + R_B) without a mismatched
%   channel population.
% ============================================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ================= USER SETTINGS ============================
data_folder        = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';   % paired ISI dataset
single_elec_folder = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_Single1';   % single-electrode dataset

Electrode_Type = 2;

fixed_macro_win_ms = [2, 40];      % window for the precomputed SpikeCount
baseline_win_ms    = [-50, -10];   % baseline window for that correction
spike_time_win_ms  = [-50, 100];   % window for exported raw spike times

FS = 30000;

% Only export these (Set, Amp, ISI) rows. Leave empty ({}) to export
% EVERYTHING that has at least one responding channel. The matching
% single-electrode conditions are derived automatically -- no separate
% list to maintain.
IncludePairedConditions = {
  % Set  Amp(uA)  ISI(ms)
    1,   5,      0;
    1,   5,      3;
    1,   5,      4;
    1,   5,      5;
    1,   5,      6;
    1,   5,      7;
    1,   5,      8;
    1,   5,      9;
    1,   5,      10;
    1,   5,      11;
    1,   5,      12;
    1,   5,      13;
    1,   5,      14;
    1,   5,      15;
    1,   5,      17;
    1,   5,      20;
};

output_dir = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Luke_Data/DX023_5uA_v1';
output_name = 'DX023_5uA';   % base filename for this export bundle

%% ================= RESTORE WORKING DIRECTORY ON EXIT =========
origDir = pwd;
dirCleanup = onCleanup(@() cd(origDir)); %#ok<NASGU>
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% =================== 1. LOAD PAIRED (ISI) DATA ====================
removeAppleDoubleFiles(data_folder);
removeAppleDoubleFiles(single_elec_folder);

[R, sp, trig, S, QC] = load_experiment_data(data_folder);
[~, AnimalID] = fileparts(fileparts(data_folder));

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

%% =========== 3. DERIVE THE FULL LIST OF (Set,Amp,ISI) TO EXPORT ===========
if isempty(IncludePairedConditions)
    pairedRows = zeros(0,3);
    for ss = 1:nSets
        for ai = 1:numel(Amps)
            for pi = 1:numel(PTDs_ms)
                pairedRows(end+1,:) = [ss, Amps(ai), PTDs_ms(pi)]; %#ok<AGROW>
            end
        end
    end
else
    pairedRows = cell2mat(IncludePairedConditions);
end

% The matching single-electrode conditions are simply the same (Set, Amp)
% pairs -- Single is built per-Set in section 5 below, reusing each
% Set's own resp_channels_per_set{ss}.
setAmpPairs = unique(pairedRows(:,[1 2]), 'rows');

%% =================== 4. BUILD NESTED STRUCT: PAIRED =====================
fprintf('\nExporting PAIRED (ISI) data...\n');
Paired = struct('SetNumber',{}, 'Electrodes',{}, 'Amp',{});

setsToDo = unique(pairedRows(:,1));
for ss = setsToDo(:)'
    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    channels_this_set = resp_channels_per_set{ss};
    if isempty(channels_this_set), continue; end

    setEntry.SetNumber = ss;
    setEntry.Electrodes = stimCh;
    setEntry.Amp = struct('Amplitude_uA',{}, 'ISI',{});

    ampsToDo = unique(pairedRows(pairedRows(:,1)==ss, 2));
    for target_amp = ampsToDo(:)'
        ai = find(abs(Amps - target_amp) < 1e-4, 1);
        if isempty(ai), continue; end

        ampEntry.Amplitude_uA = target_amp;
        ampEntry.ISI = struct('ISI_ms',{}, 'Channel',{});

        isisToDo = unique(pairedRows(pairedRows(:,1)==ss & pairedRows(:,2)==target_amp, 3));
        for target_isi = isisToDo(:)'
            pi = find(abs(PTDs_ms - target_isi) < 1e-4, 1);
            if isempty(pi), continue; end

            tr_ids_base = find(combClass == ss & ampIdx == ai & ptdIdx == pi);
            if isempty(tr_ids_base), continue; end

            isiEntry.ISI_ms = target_isi;
            isiEntry.Channel = struct('Channel',{}, 'SpikeCount',{}, 'SpikeTimes_ms',{});

            for ch_idx = channels_this_set(:)'
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

                chanEntry.Channel = ch_idx;
                chanEntry.SpikeCount = countsPerTrial;
                chanEntry.SpikeTimes_ms = spikeTimesPerTrial;
                isiEntry.Channel(end+1) = chanEntry;
            end

            if ~isempty(isiEntry.Channel)
                ampEntry.ISI(end+1) = isiEntry;
            end
        end

        if ~isempty(ampEntry.ISI)
            setEntry.Amp(end+1) = ampEntry;
        end
    end

    if ~isempty(setEntry.Amp)
        Paired(end+1) = setEntry;
    end
    fprintf('  Set %d (Ch:%s) done.\n', ss, num2str(stimCh));
end

%% =================== 5. LOAD + BUILD NESTED STRUCT: SINGLE (per-Set) =====================
fprintf('\nExporting SINGLE-ELECTRODE data (scoped per Set)...\n');
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

Single = struct('SetNumber',{}, 'Electrodes',{}, 'Amp',{});

for r = 1:size(setAmpPairs,1)
    ss = setAmpPairs(r,1); target_amp = setAmpPairs(r,2);
    if ss > nSets, continue; end

    stimCh = uniqueComb(ss,:); stimCh = stimCh(stimCh>0);
    % SAME channel population as this Set's Paired data -- not a union.
    channels_this_set = resp_channels_per_set{ss};
    if isempty(channels_this_set), continue; end

    ai_se = find(abs(Amps_se - target_amp) < 1e-4, 1);
    if isempty(ai_se)
        fprintf('WARNING: Set %d: no matching %.1f uA in single-electrode dataset. Skipping.\n', ss, target_amp);
        continue;
    end
    matched_amp_se = Amps_se(ai_se);

    % Find (or start) this Set's entry in Single
    existingIdx = [];
    for k = 1:numel(Single)
        if Single(k).SetNumber == ss, existingIdx = k; break; end
    end
    if isempty(existingIdx)
        clear setEntrySE
        setEntrySE.SetNumber = ss;
        setEntrySE.Electrodes = stimCh;
        setEntrySE.Amp = struct('Amplitude_uA',{}, 'Electrode',{});
        Single(end+1) = setEntrySE;
        existingIdx = numel(Single);
    end

    clear ampEntrySE
    ampEntrySE.Amplitude_uA = target_amp;
    ampEntrySE.Electrode = struct('Electrode',{}, 'Channel',{});

    for e = stimCh(:)'
        tr_ids_base = find(trialElectrode_se == e & abs(trialAmps_se - matched_amp_se) < 1e-4);
        if isempty(tr_ids_base)
            fprintf('WARNING: Set %d, Electrode %d: no trials found in single-electrode dataset. Skipping.\n', ss, e);
            continue;
        end

        clear elecEntrySE
        elecEntrySE.Electrode = e;
        elecEntrySE.Channel = struct('Channel',{}, 'SpikeCount',{}, 'SpikeTimes_ms',{});

        for ch_idx = channels_this_set(:)'
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

            clear chanEntrySE
            chanEntrySE.Channel = ch_idx;
            chanEntrySE.SpikeCount = countsPerTrial;
            chanEntrySE.SpikeTimes_ms = spikeTimesPerTrial;
            elecEntrySE.Channel(end+1) = chanEntrySE;
        end

        if ~isempty(elecEntrySE.Channel)
            ampEntrySE.Electrode(end+1) = elecEntrySE;
        end
    end

    if ~isempty(ampEntrySE.Electrode)
        Single(existingIdx).Amp(end+1) = ampEntrySE;
    end
    fprintf('  Set %d | %.1f uA done.\n', ss, target_amp);
end

%% =================== 6. METADATA =====================
ExportMetadata = struct();
ExportMetadata.AnimalID = AnimalID;
ExportMetadata.FixedMacroWindow_ms = fixed_macro_win_ms;
ExportMetadata.BaselineWindow_ms = baseline_win_ms;
ExportMetadata.SpikeTimeWindow_ms = spike_time_win_ms;
ExportMetadata.SamplingRate_Hz = FS;

%% =================== 7. SAVE .mat =====================
mat_path = fullfile(output_dir, [output_name '_Export.mat']);
save(mat_path, 'Paired', 'Single', 'ExportMetadata');

%% =================== 8. FLATTEN + SAVE CSVs =====================
[PairedCounts, PairedTimes] = flatten_paired(Paired);
[SingleCounts, SingleTimes] = flatten_single(Single);

writetable(PairedCounts, fullfile(output_dir, 'Paired_SpikeCounts.csv'));
writetable(PairedTimes,  fullfile(output_dir, 'Paired_SpikeTimes.csv'));
writetable(SingleCounts, fullfile(output_dir, 'Single_SpikeCounts.csv'));
writetable(SingleTimes,  fullfile(output_dir, 'Single_SpikeTimes.csv'));

%% =================== 9. WRITE README =====================
readme_path = fullfile(output_dir, 'README.txt');
fid = fopen(readme_path, 'w');
fprintf(fid, 'DATA EXPORT README\nAnimal: %s\n\n', AnimalID);

fprintf(fid, 'FILES\n-----\n');
fprintf(fid, '%s_Export.mat\n', output_name);
fprintf(fid, '  Paired: struct array, one entry per stimulation Set.\n');
fprintf(fid, '    .SetNumber, .Electrodes ([ElecA ElecB])\n');
fprintf(fid, '    .Amp(i).Amplitude_uA\n');
fprintf(fid, '    .Amp(i).ISI(j).ISI_ms   (0 = simultaneous)\n');
fprintf(fid, '    .Amp(i).ISI(j).Channel(k).Channel        recording channel #\n');
fprintf(fid, '    .Amp(i).ISI(j).Channel(k).SpikeCount     [Nx1] one value per trial\n');
fprintf(fid, '    .Amp(i).ISI(j).Channel(k).SpikeTimes_ms  {Nx1} one vector per trial\n');
fprintf(fid, '    SpikeCount(t) and SpikeTimes_ms{t} refer to the SAME trial t.\n');
fprintf(fid, '    Trial order is a fresh renumbering, not the original recording order.\n\n');

fprintf(fid, '  Single: struct array, one entry per Set (same Set numbers\n');
fprintf(fid, '  as Paired) -- the single-electrode-alone data for that\n');
fprintf(fid, '  Set''s two component electrodes, using the SAME responding-channel\n');
fprintf(fid, '  population as that Set''s Paired data (so Single(k) and Paired(k) for\n');
fprintf(fid, '  the same Set can be combined directly, e.g. R_predicted = R_A + R_B).\n');
fprintf(fid, '  If a physical electrode appears in more than one Set, its data is\n');
fprintf(fid, '  exported once per Set, each time scoped to that Set''s own channels.\n');
fprintf(fid, '    .SetNumber, .Electrodes ([ElecA ElecB])\n');
fprintf(fid, '    .Amp(i).Amplitude_uA\n');
fprintf(fid, '    .Amp(i).Electrode(j).Electrode          which of the pair''s electrodes\n');
fprintf(fid, '    .Amp(i).Electrode(j).Channel(k).Channel / .SpikeCount / .SpikeTimes_ms\n');
fprintf(fid, '    (same per-trial layout as Paired above)\n\n');

fprintf(fid, '  ExportMetadata: AnimalID, FixedMacroWindow_ms ([%g %g]), ', ...
    fixed_macro_win_ms(1), fixed_macro_win_ms(2));
fprintf(fid, 'BaselineWindow_ms ([%g %g]), SpikeTimeWindow_ms ([%g %g]), SamplingRate_Hz (%d).\n\n', ...
    baseline_win_ms(1), baseline_win_ms(2), spike_time_win_ms(1), spike_time_win_ms(2), FS);

fprintf(fid, 'Paired_SpikeCounts.csv / Paired_SpikeTimes.csv\n');
fprintf(fid, 'Single_SpikeCounts.csv / Single_SpikeTimes.csv\n');
fprintf(fid, '  Flat, portable equivalents of the .mat data (same numbers, just\n');
fprintf(fid, '  tabular instead of nested -- CSV can''t nest). Columns:\n');
fprintf(fid, '  Paired: Set, ElectrodeA, ElectrodeB, Amplitude_uA, ISI_ms, Channel, TrialIndex, [SpikeCount | SpikeTime_ms]\n');
fprintf(fid, '  Single: Set, Electrode, Amplitude_uA, Channel, TrialIndex, [SpikeCount | SpikeTime_ms]\n');
fprintf(fid, '  TrialIndex is the same fresh renumbering described above.\n\n');

fprintf(fid, 'SpikeCount: net baseline-corrected spike count in the fixed window above.\n');
fprintf(fid, 'SpikeTime_ms / SpikeTimes_ms: spike time(s) relative to the stimulation\n');
fprintf(fid, '  trigger, within the spike-time window above. No waveforms included.\n\n');

fprintf(fid, 'Only responding channels and QC-passed trials are included; both are\n');
fprintf(fid, 'applied silently (no channel/trial exclusion lists are exported).\n');
fclose(fid);

fprintf('\n>>> Export complete: %s\n', output_dir);

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

function [countsTable, timesTable] = flatten_paired(PairedStruct)
% Flattens the nested Paired struct into two flat tables for CSV export.
c_Set=[]; c_EA=[]; c_EB=[]; c_Amp=[]; c_ISI=[]; c_Chan=[]; c_Trial=[]; c_Count=[];
t_Set=[]; t_EA=[]; t_EB=[]; t_Amp=[]; t_ISI=[]; t_Chan=[]; t_Trial=[]; t_Time=[];

for si = 1:numel(PairedStruct)
    setEntry = PairedStruct(si);
    elecs = setEntry.Electrodes;
    ea = elecs(1); eb = NaN; if numel(elecs)>=2, eb = elecs(2); end

    for ai = 1:numel(setEntry.Amp)
        ampEntry = setEntry.Amp(ai);
        for pi = 1:numel(ampEntry.ISI)
            isiEntry = ampEntry.ISI(pi);
            for ci = 1:numel(isiEntry.Channel)
                chanEntry = isiEntry.Channel(ci);
                nTr = numel(chanEntry.SpikeCount);
                trialIdx = (1:nTr)';

                c_Set=[c_Set; repmat(setEntry.SetNumber,nTr,1)]; %#ok<AGROW>
                c_EA=[c_EA; repmat(ea,nTr,1)]; %#ok<AGROW>
                c_EB=[c_EB; repmat(eb,nTr,1)]; %#ok<AGROW>
                c_Amp=[c_Amp; repmat(ampEntry.Amplitude_uA,nTr,1)]; %#ok<AGROW>
                c_ISI=[c_ISI; repmat(isiEntry.ISI_ms,nTr,1)]; %#ok<AGROW>
                c_Chan=[c_Chan; repmat(chanEntry.Channel,nTr,1)]; %#ok<AGROW>
                c_Trial=[c_Trial; trialIdx]; %#ok<AGROW>
                c_Count=[c_Count; chanEntry.SpikeCount]; %#ok<AGROW>

                for k = 1:nTr
                    nSpk = numel(chanEntry.SpikeTimes_ms{k});
                    if nSpk==0, continue; end
                    t_Set=[t_Set; repmat(setEntry.SetNumber,nSpk,1)]; %#ok<AGROW>
                    t_EA=[t_EA; repmat(ea,nSpk,1)]; %#ok<AGROW>
                    t_EB=[t_EB; repmat(eb,nSpk,1)]; %#ok<AGROW>
                    t_Amp=[t_Amp; repmat(ampEntry.Amplitude_uA,nSpk,1)]; %#ok<AGROW>
                    t_ISI=[t_ISI; repmat(isiEntry.ISI_ms,nSpk,1)]; %#ok<AGROW>
                    t_Chan=[t_Chan; repmat(chanEntry.Channel,nSpk,1)]; %#ok<AGROW>
                    t_Trial=[t_Trial; repmat(k,nSpk,1)]; %#ok<AGROW>
                    t_Time=[t_Time; chanEntry.SpikeTimes_ms{k}]; %#ok<AGROW>
                end
            end
        end
    end
end

countsTable = table(c_Set,c_EA,c_EB,c_Amp,c_ISI,c_Chan,c_Trial,c_Count, 'VariableNames', ...
    {'Set','ElectrodeA','ElectrodeB','Amplitude_uA','ISI_ms','Channel','TrialIndex','SpikeCount'});
timesTable = table(t_Set,t_EA,t_EB,t_Amp,t_ISI,t_Chan,t_Trial,t_Time, 'VariableNames', ...
    {'Set','ElectrodeA','ElectrodeB','Amplitude_uA','ISI_ms','Channel','TrialIndex','SpikeTime_ms'});
end

function [countsTable, timesTable] = flatten_single(SingleStruct)
% Flattens the nested Single struct into two flat tables for CSV export.
c_Set=[]; c_Elec=[]; c_Amp=[]; c_Chan=[]; c_Trial=[]; c_Count=[];
t_Set=[]; t_Elec=[]; t_Amp=[]; t_Chan=[]; t_Trial=[]; t_Time=[];

for si = 1:numel(SingleStruct)
    setEntry = SingleStruct(si);

    for ai = 1:numel(setEntry.Amp)
        ampEntry = setEntry.Amp(ai);
        for ei = 1:numel(ampEntry.Electrode)
            elecEntry = ampEntry.Electrode(ei);
            for ci = 1:numel(elecEntry.Channel)
                chanEntry = elecEntry.Channel(ci);
                nTr = numel(chanEntry.SpikeCount);
                trialIdx = (1:nTr)';

                c_Set=[c_Set; repmat(setEntry.SetNumber,nTr,1)]; %#ok<AGROW>
                c_Elec=[c_Elec; repmat(elecEntry.Electrode,nTr,1)]; %#ok<AGROW>
                c_Amp=[c_Amp; repmat(ampEntry.Amplitude_uA,nTr,1)]; %#ok<AGROW>
                c_Chan=[c_Chan; repmat(chanEntry.Channel,nTr,1)]; %#ok<AGROW>
                c_Trial=[c_Trial; trialIdx]; %#ok<AGROW>
                c_Count=[c_Count; chanEntry.SpikeCount]; %#ok<AGROW>

                for k = 1:nTr
                    nSpk = numel(chanEntry.SpikeTimes_ms{k});
                    if nSpk==0, continue; end
                    t_Set=[t_Set; repmat(setEntry.SetNumber,nSpk,1)]; %#ok<AGROW>
                    t_Elec=[t_Elec; repmat(elecEntry.Electrode,nSpk,1)]; %#ok<AGROW>
                    t_Amp=[t_Amp; repmat(ampEntry.Amplitude_uA,nSpk,1)]; %#ok<AGROW>
                    t_Chan=[t_Chan; repmat(chanEntry.Channel,nSpk,1)]; %#ok<AGROW>
                    t_Trial=[t_Trial; repmat(k,nSpk,1)]; %#ok<AGROW>
                    t_Time=[t_Time; chanEntry.SpikeTimes_ms{k}]; %#ok<AGROW>
                end
            end
        end
    end
end

countsTable = table(c_Set,c_Elec,c_Amp,c_Chan,c_Trial,c_Count, 'VariableNames', ...
    {'Set','Electrode','Amplitude_uA','Channel','TrialIndex','SpikeCount'});
timesTable = table(t_Set,t_Elec,t_Amp,t_Chan,t_Trial,t_Time, 'VariableNames', ...
    {'Set','Electrode','Amplitude_uA','Channel','TrialIndex','SpikeTime_ms'});
end