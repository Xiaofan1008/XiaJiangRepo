%% ============================================================
% QC_SequentialStitch_RasterPSTH_v1.m
% Step 3 of the sequential-stimulation artifact recovery pipeline: visual
% QC of the stitching done by SpikeRecovery_SequentialStitch_v1.m.
%
% For each (channel, amplitude, stim set, PTD), produces one figure with:
%   Row 1 -- raster: single-pulse ch1 (the substitution source)
%   Row 2 -- raster: simultaneous (PTD=0, same stim pair, untouched)
%   Row 3 -- raster: sequential, STITCHED -- single row, but each spike
%            is colored by whether it's SUBSTITUTED (relative time <
%            that trial's blank_end_ms) or ORIGINAL/untouched (>=). A
%            dashed line marks the group's median blank_end_ms.
%   Row 4 -- PSTH overlay: single / simultaneous / sequential-stitched
%            (one combined curve each, not split by substituted/original)
%
% No "raw sequential" row -- the point of this script is to check the
% RECOVERED data against its two reference conditions, not to re-show the
% artifact itself.
%
% The substituted/original split needs no extra per-spike bookkeeping:
% SpikeRecovery_SequentialStitch_v1.m's deletion step always removes
% every original spike below blank_end_ms before substituting, so a
% spike's relative time alone tells you which side of that boundary it's
% on.
% ============================================================
clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/'));

%% ================= USER SETTINGS =================
% -------- Single-pulse (source) dataset --------
single_times_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_filtered.mat';
single_folder     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1';

% -------- Sequential dataset --------
seq_filtered_times_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_filtered.mat';
                  
seq_stitched_times_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_stitched.mat';

seq_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1';

fst_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_FirstSpikeTimes.mat';

Electrode_Type = 3;     % 0:single shank rigid; 1:single shank flex; 2:four shank flex

%% -------- Analysis settings --------
target_channels = [33:40];
plot_amps       = [10];
target_PTDs     = [];   % which sequential PTDs to plot; [] = all > 0
target_sets     = [6]; 

ras_win    = [-10 50];   % ms
bin_ms     = 1;
smooth_ms  = 5;
FS         = 30000;      % overridden below if saved in a loaded file

%% -------- Colors --------
col_single    = [0.30 0.45 0.75];   % steel blue  -- single-pulse ch1
col_sim       = [0.30 0.65 0.30];   % green       -- simultaneous
col_original  = [0.20 0.20 0.20];   % near-black  -- stitched, ORIGINAL spikes
col_substitute = [0.90 0.45 0.10];  % orange      -- stitched, SUBSTITUTED spikes

%% ================= LOAD SINGLE-PULSE (SOURCE) DATASET =================
assert(isfile(single_times_file), 'Cannot find %s.', single_times_file);
S1 = load(single_times_file);
assert(isfield(S1, 'sp_clipped'), 'Variable "sp_clipped" not found in %s.', single_times_file);
sp_single = S1.sp_clipped;
if isfield(S1, 'fs'), FS = S1.fs; end

cd(single_folder);
trig_single = loadTrig(0);
pfile1 = dir(fullfile(single_folder, '*_exp_datafile_*.mat'));
assert(~isempty(pfile1), 'No *_exp_datafile_*.mat found in %s.', single_folder);
Sp1 = load(fullfile(single_folder, pfile1(1).name), 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_single  = Sp1.StimParams;
E_MAP_single       = Sp1.E_MAP;
n_Trials_single    = Sp1.n_Trials;
simultaneous_stim1 = Sp1.simultaneous_stim;

trialAmps_single = cell2mat(StimParams_single(2:end,16));
trialAmps_single = trialAmps_single(1:simultaneous_stim1:end);
stimNames_single = StimParams_single(2:end,1);
[~, idx_all_single] = ismember(stimNames_single, E_MAP_single(2:end));
stimChPerTrial_single = cell(n_Trials_single,1);
for t = 1:n_Trials_single
    rr = (t-1)*simultaneous_stim1 + (1:simultaneous_stim1);
    v  = idx_all_single(rr);
    v  = v(v>0);
    stimChPerTrial_single{t} = v(:).';
end

%% ================= LOAD SEQUENTIAL DATASET (filtered + stitched) =================
assert(isfile(seq_filtered_times_file), 'Cannot find %s.', seq_filtered_times_file);
Sfilt = load(seq_filtered_times_file);
assert(isfield(Sfilt, 'sp_clipped'), 'Variable "sp_clipped" not found in %s.', seq_filtered_times_file);
sp_seq_filtered = Sfilt.sp_clipped;

assert(isfile(seq_stitched_times_file), 'Cannot find %s.', seq_stitched_times_file);
Sstitch = load(seq_stitched_times_file);
assert(isfield(Sstitch, 'sp_clipped'), 'Variable "sp_clipped" not found in %s.', seq_stitched_times_file);
assert(isfield(Sstitch, 'blank_end_ms'), ...
    '"blank_end_ms" not found in %s -- re-run the current SpikeRecovery_SequentialStitch_v1.m to produce it.', ...
    seq_stitched_times_file);
sp_seq_stitched = Sstitch.sp_clipped;
blank_end_ms    = Sstitch.blank_end_ms;

nChn = numel(sp_seq_stitched);
assert(numel(sp_seq_filtered) == nChn, 'Filtered (%d) and stitched (%d) sequential files disagree on channel count.', ...
    numel(sp_seq_filtered), nChn);
assert(numel(sp_single) == nChn, 'Single-pulse (%d) and sequential (%d) files disagree on channel count.', ...
    numel(sp_single), nChn);

cd(seq_folder);
trig_seq = loadTrig(0);
pfile2 = dir(fullfile(seq_folder, '*_exp_datafile_*.mat'));
assert(~isempty(pfile2), 'No *_exp_datafile_*.mat found in %s.', seq_folder);
Sp2 = load(fullfile(seq_folder, pfile2(1).name), 'StimParams','simultaneous_stim','E_MAP','n_Trials');
StimParams_seq     = Sp2.StimParams;
E_MAP_seq          = Sp2.E_MAP;
n_Trials_seq       = Sp2.n_Trials;
simultaneous_stim2 = Sp2.simultaneous_stim;

trialAmps_seq = cell2mat(StimParams_seq(2:end,16));
trialAmps_seq = trialAmps_seq(1:simultaneous_stim2:end);
stimNames_seq = StimParams_seq(2:end,1);
[~, idx_all_seq] = ismember(stimNames_seq, E_MAP_seq(2:end));
stimChPerTrial_seq = cell(n_Trials_seq,1);
for t = 1:n_Trials_seq
    rr = (t-1)*simultaneous_stim2 + (1:simultaneous_stim2);
    v  = idx_all_seq(rr);
    v  = v(v>0);
    stimChPerTrial_seq{t} = v(:).';
end
comb_seq = zeros(n_Trials_seq, simultaneous_stim2);
for t = 1:n_Trials_seq
    v = stimChPerTrial_seq{t};
    comb_seq(t, 1:numel(v)) = v;
end
[uniqueComb_seq, ~, combClass_seq] = unique(comb_seq, 'rows', 'stable');
nSeqSets = size(uniqueComb_seq,1);

% -------- Print detected sets, and resolve which ones to actually plot --------
fprintf('\nDetected stimulation sets:\n');
for si = 1:nSeqSets
    v = uniqueComb_seq(si,:); v = v(v>0);
    fprintf('  [%d] %s\n', si, strjoin(arrayfun(@(x) sprintf('Ch%d',x), v, 'UniformOutput', false), '->'));
end
if isempty(target_sets)
    sets_to_plot = 1:nSeqSets;
else
    sets_to_plot = intersect(1:nSeqSets, target_sets, 'stable');
    if isempty(sets_to_plot)
        warning('None of target_sets (%s) matched a detected set index (1..%d) -- nothing will be plotted.', ...
            mat2str(target_sets), nSeqSets);
    end
end
fprintf('Plotting sets: %s\n', mat2str(sets_to_plot));

%% ================= LOAD PTD / SIMULTANEOUS FLAG =================
assert(isfile(fst_file), 'Cannot find %s.', fst_file);
Sf = load(fst_file);
PTD_ms         = Sf.PTD_ms;
isSimultaneous = Sf.isSimultaneous;
assert(Sf.n_Trials == n_Trials_seq, ...
    'FirstSpikeTimes file (%d trials) does not match the sequential dataset (%d trials).', Sf.n_Trials, n_Trials_seq);

available_PTDs = unique(PTD_ms(~isSimultaneous));
if isempty(target_PTDs)
    use_PTDs = available_PTDs;
else
    use_PTDs = intersect(available_PTDs, target_PTDs);
end
fprintf('Plotting sequential PTDs: %s ms\n', num2str(use_PTDs'));

%% ================= PSTH KERNEL & CHANNEL MAP =================
edges = ras_win(1):bin_ms:ras_win(2);
ctrs  = edges(1:end-1) + diff(edges)/2;
bin_s = bin_ms/1000;
g = exp(-0.5*((0:smooth_ms-1)/(smooth_ms/2)).^2);
g = g/sum(g);
d = ChnMap(Electrode_Type, nChn);   % physical channel per map position

%% ===================== MAIN PLOTTING LOOP =====================
for ich = target_channels
    ch = d(ich);
    if ch > nChn, continue; end

    for amp_val = plot_amps
        mask_single_amp = (trialAmps_single == amp_val);

        for set_id = sets_to_plot
            stimVec = uniqueComb_seq(set_id,:);
            stimVec = stimVec(stimVec > 0);
            if isempty(stimVec), continue; end

            % -------- Identify ch1 (substitution source), matched by NAME --------
            ch_first_idx_seq = stimVec(1);
            ch_name = E_MAP_seq{ch_first_idx_seq + 1};
            ch_idx_in_single_map = find(strcmp(E_MAP_single(2:end), ch_name));
            if isempty(ch_idx_in_single_map)
                continue;
            end
            mask_single_chan = cellfun(@(x) ismember(ch_idx_in_single_map, x), stimChPerTrial_single);
            single_trials_for_group = find(mask_single_amp & mask_single_chan);

            stimLabel = strjoin(arrayfun(@(x) sprintf('Ch%d',x), stimVec, 'UniformOutput', false), '->');

            % -------- Simultaneous trials for this same stim set --------
            sim_trials = find(trialAmps_seq == amp_val & combClass_seq == set_id & PTD_ms == 0);

            for p_idx = 1:numel(use_PTDs)
                ptd_val = use_PTDs(p_idx);
                seq_trials_group = find(trialAmps_seq == amp_val & combClass_seq == set_id & ...
                                         PTD_ms == ptd_val & ~isSimultaneous);
                if isempty(seq_trials_group) && isempty(single_trials_for_group) && isempty(sim_trials)
                    continue;
                end

                figure('Color','w','Position',[200 100 800 900], 'Name', ...
                    sprintf('Ch %d | %s | %g uA | PTD %g ms', ich, stimLabel, amp_val, ptd_val));
                tl = tiledlayout(4,1,'TileSpacing','compact','Padding','compact');
                title(tl, sprintf('Channel %d (phys %d) -- %s -- %g uA -- PTD %g ms', ...
                    ich, ch, stimLabel, amp_val, ptd_val), 'FontSize', 13);
                max_psth_rate = 0;

                % ======== ROW 1: single-pulse ch1 raster ========
                ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'off');
                S_ch = sp_single{ch};
                y = 0;
                for tr = single_trials_for_group'
                    t0 = trig_single(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    for k = 1:numel(tt)
                        plot(ax1,[tt(k) tt(k)],[y y+0.8],'Color',col_single,'LineWidth',1.3);
                    end
                    y = y+1;
                end
                title(ax1, sprintf('Single-pulse Ch%d (source, n=%d trials)', ch_first_idx_seq, y), 'FontSize',10);
                xline(ax1,0,'r--','HandleVisibility','off');
                xlim(ax1,ras_win); ylim(ax1,[0 max(1,y)]);

                % ======== ROW 2: simultaneous raster ========
                ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'off');
                S_ch = sp_seq_filtered{ch};
                y = 0;
                for tr = sim_trials'
                    t0 = trig_seq(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    for k = 1:numel(tt)
                        plot(ax2,[tt(k) tt(k)],[y y+0.8],'Color',col_sim,'LineWidth',1.3);
                    end
                    y = y+1;
                end
                title(ax2, sprintf('Simultaneous (untouched, n=%d trials)', y), 'FontSize',10);
                xline(ax2,0,'r--','HandleVisibility','off');
                xlim(ax2,ras_win); ylim(ax2,[0 max(1,y)]);

                % ======== ROW 3: sequential STITCHED raster, 2-color ========
                ax3 = nexttile(tl); hold(ax3,'on'); box(ax3,'off');
                S_ch = sp_seq_stitched{ch};
                this_blank_end = blank_end_ms{ch};
                group_blank_ends = this_blank_end(seq_trials_group);
                y = 0;
                for tr = seq_trials_group'
                    t0 = trig_seq(tr)/FS*1000;
                    tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                    be = this_blank_end(tr);
                    for k = 1:numel(tt)
                        if isfinite(be) && tt(k) < be
                            col = col_substitute;
                        else
                            col = col_original;
                        end
                        plot(ax3,[tt(k) tt(k)],[y y+0.8],'Color',col,'LineWidth',1.3);
                    end
                    y = y+1;
                end
                med_blank_end = median(group_blank_ends(isfinite(group_blank_ends)));
                if isfinite(med_blank_end)
                    xline(ax3, med_blank_end, '--', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
                end
                title(ax3, sprintf('Sequential, STITCHED (orange=substituted, black=original, n=%d trials)', y), ...
                    'FontSize',10);
                xline(ax3,0,'r--','HandleVisibility','off');
                xline(ax3,ptd_val,'k:','HandleVisibility','off');
                xlim(ax3,ras_win); ylim(ax3,[0 max(1,y)]);

                % ======== ROW 4: PSTH overlay ========
                ax4 = nexttile(tl); hold(ax4,'on'); box(ax4,'off');

                if ~isempty(single_trials_for_group)
                    S_ch = sp_single{ch};
                    counts = zeros(1,numel(edges)-1);
                    for tr = single_trials_for_group'
                        t0 = trig_single(tr)/FS*1000;
                        tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                        counts = counts + histcounts(tt,edges);
                    end
                    rate = filter(g,1, counts/(numel(single_trials_for_group)*bin_s));
                    plot(ax4,ctrs,rate,'Color',col_single,'LineWidth',2,'DisplayName','Single');
                    max_psth_rate = max(max_psth_rate, max(rate));
                end

                if ~isempty(sim_trials)
                    S_ch = sp_seq_filtered{ch};
                    counts = zeros(1,numel(edges)-1);
                    for tr = sim_trials'
                        t0 = trig_seq(tr)/FS*1000;
                        tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                        counts = counts + histcounts(tt,edges);
                    end
                    rate = filter(g,1, counts/(numel(sim_trials)*bin_s));
                    plot(ax4,ctrs,rate,'Color',col_sim,'LineWidth',2,'DisplayName','Simultaneous');
                    max_psth_rate = max(max_psth_rate, max(rate));
                end

                if ~isempty(seq_trials_group)
                    S_ch = sp_seq_stitched{ch};
                    counts = zeros(1,numel(edges)-1);
                    for tr = seq_trials_group'
                        t0 = trig_seq(tr)/FS*1000;
                        tt = S_ch(:,1); tt = tt(tt>=t0+ras_win(1) & tt<=t0+ras_win(2)) - t0;
                        counts = counts + histcounts(tt,edges);
                    end
                    rate = filter(g,1, counts/(numel(seq_trials_group)*bin_s));
                    plot(ax4,ctrs,rate,'Color',col_substitute,'LineWidth',2,'DisplayName','Sequential (stitched)');
                    max_psth_rate = max(max_psth_rate, max(rate));
                end

                xline(ax4,0,'r--','HandleVisibility','off');
                xline(ax4,ptd_val,'k:','HandleVisibility','off');
                xlim(ax4,ras_win); xlabel(ax4,'Time (ms)'); ylabel(ax4,'Rate (sp/s)');
                legend(ax4,'Box','off','Location','northeast');
                if max_psth_rate > 0
                    ylim(ax4,[0 max_psth_rate*1.1]);
                end
            end
        end
    end
end