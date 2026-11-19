%% ========================================================================
% SPIKE COUNT vs DISTANCE - LINEARITY (FIXED WINDOW, NO SHIFT)
%
% PURPOSE
%   For a paired-stimulation (Sim/Seq) dataset, compare the ACTUAL spike
%   count evoked by paired stimulation of two channels against the
%   PREDICTED count (the sum of each channel's own single-pulse spike
%   count), across the different channel-pair DISTANCES present in the
%   dataset. This is the first pass of the linearity analysis: a simple
%   fixed-window spike count, no PSTH time-shifting (that's a later,
%   separate analysis).
%
% MULTIPLE DATASET-PAIRS
%   Each entry in DatasetPairs is one paired (Sim/Seq) recording plus its
%   matching single-pulse recording (same session). Different entries can
%   come from entirely different sessions/animals/days - each one only
%   needs to cover whichever distances it has; the script auto-detects
%   the distances present in every entry and merges them all onto the
%   same figures. All entries share the same Electrode_Type (same probe).
%   Channels between a paired entry and its single entry are matched
%   automatically by channel number - no manual pairing needed.
%
% CHANNEL SCOPE
%   For each paired condition (Set x Amp x PTD), the channel list used for
%   BOTH the actual count and the two single-channel components of the
%   predicted count is that exact condition's own saved responding-channel
%   list (<base_name>_RespondingChannels.mat from the PAIRED dataset).
%   Nothing is unioned across conditions, pairs, or datasets; the single
%   dataset's own responding-channel information (if any) is not used.
%
% METRIC
%   Baseline-corrected spike count in a FIXED analysis_win_ms window
%   (default [0,20] ms, separate from the [2,40] ms window used by the
%   bad-trial QC pipeline), averaged over trials, then averaged over the
%   condition's responding channels:
%     Actual    = mean_channels( mean_trials( corrected count ) )  [paired]
%     Predicted = SingleA + SingleB, each computed the same way from the
%                 matching single-channel condition, over the SAME
%                 channel list as Actual.
%   If a channel's single-pulse match can't be found (channel not tested
%   alone at that amplitude, say), Predicted is left NaN for that point
%   and only Actual is plotted.
%
% BAD TRIALS
%   Each dataset excludes trials using its OWN saved
%   <base_name>_BadTrials.mat (AllBadAbsoluteTrialIDs) - paired bad trials
%   are excluded from the paired dataset only, single bad trials from the
%   single dataset only.
%
% OUTPUT
%   Plot only for now (one figure per PTD value; color = amplitude, solid
%   = actual, dashed = predicted; x = distance, y = mean baseline-
%   corrected spike count). Nothing is saved to disk.
% ========================================================================

clear;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

% One entry per dataset-pair. Add as many as you need - copy the block
% below and increment the index.
DatasetPairs = struct( ...
    'simseq_spike_file',{}, ...
    'simseq_dataset_folder',{}, ...
    'single_spike_file',{}, ...
    'single_dataset_folder',{});

DatasetPairs(1).simseq_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1/Xia_Linearity_S4_100_400um_SimSeq1.sp_xia_stitched.mat';
DatasetPairs(1).simseq_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_SimSeq1';
DatasetPairs(1).single_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1/Xia_Linearity_S4_100_400um_Single1.sp_xia_filtered.mat';
DatasetPairs(1).single_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_100_400um_Single1';

DatasetPairs(2).simseq_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_500_700um_SimSeq1/Xia_Linearity_S4_500_700um_SimSeq1.sp_xia_stitched.mat';
DatasetPairs(2).simseq_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_500_700um_SimSeq1';
DatasetPairs(2).single_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_500_700um_Single1/Xia_Linearity_S4_500_700um_Single1.sp_xia_filtered.mat';
DatasetPairs(2).single_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_500_700um_Single1';

% Electrode type - shared by EVERY dataset-pair above (same probe):
%   0 = rigid single-shank probe
%   1 = flexible single-shank probe
%   2 = four-shank flexible probe
%   3 = 64chn+32chn hybrid
Electrode_Type = 3;

FS_fallback = 30000;   % used only if fs/FS is absent from a spike file

% Counting windows relative to each trial's own stimulation trigger (t=0).
% Upper boundary excluded. analysis_win_ms is fixed at [0,20] per the
% linearity-metric design; baseline_win_ms can be changed if needed.
baseline_win_ms = [-50 -5];
analysis_win_ms = [0 20];

Amp_Tol_uA = 0.001;
PTD_Tol_ms = 0.001;

%% =========================== INITIAL CHECKS ===========================

if isempty(DatasetPairs)
    error('DatasetPairs is empty - add at least one dataset-pair entry.');
end

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('SPIKE COUNT vs DISTANCE - LINEARITY (FIXED WINDOW, NO SHIFT)\n');
fprintf('============================================================\n');
fprintf('Dataset-pairs: %d\n',numel(DatasetPairs));
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Analysis window: [%g,%g) ms\n',analysis_win_ms);

%% =========================== MAIN COMPUTATION ==========================

Result_Set          = [];
Result_StimLabel    = {};
Result_Distance     = [];
Result_Amp          = [];
Result_PTD          = [];
Result_Actual       = [];
Result_Predicted    = [];
Result_NChannels    = [];
Result_NTrials      = [];
Result_DatasetLabel = {};

for dp = 1:numel(DatasetPairs)

    entry = DatasetPairs(dp);

    fprintf('\n============================================================\n');
    fprintf('DATASET-PAIR %d/%d\n',dp,numel(DatasetPairs));
    fprintf('============================================================\n');

    if ~isfile(entry.simseq_spike_file)
        error('Dataset-pair %d: paired spike file does not exist:\n%s',dp,entry.simseq_spike_file);
    end
    if ~isfolder(entry.simseq_dataset_folder)
        error('Dataset-pair %d: paired dataset folder does not exist:\n%s',dp,entry.simseq_dataset_folder);
    end
    if ~isfile(entry.single_spike_file)
        error('Dataset-pair %d: single spike file does not exist:\n%s',dp,entry.single_spike_file);
    end
    if ~isfolder(entry.single_dataset_folder)
        error('Dataset-pair %d: single dataset folder does not exist:\n%s',dp,entry.single_dataset_folder);
    end

    fprintf('Paired dataset: %s\n',entry.simseq_dataset_folder);
    fprintf('Single dataset: %s\n',entry.single_dataset_folder);

    fprintf('\n--- Loading paired dataset ---\n');
    Paired = load_linearity_dataset(entry.simseq_spike_file,entry.simseq_dataset_folder, ...
        Electrode_Type,FS_fallback);

    fprintf('\n--- Loading single dataset ---\n');
    Single = load_linearity_dataset(entry.single_spike_file,entry.single_dataset_folder, ...
        Electrode_Type,FS_fallback);

    dataset_label = Paired.base_name;

    %% ---- load this pair's responding-channel list ----

    responding_path = fullfile(entry.simseq_dataset_folder, ...
        sprintf('%s_RespondingChannels.mat',Paired.base_name));

    if ~isfile(responding_path)
        error('Dataset-pair %d: responding-channel file not found:\n%s',dp,responding_path);
    end

    RespLoad = load(responding_path,'RespondingConfirmed');
    if ~isfield(RespLoad,'RespondingConfirmed')
        error('Dataset-pair %d: responding-channel file does not contain RespondingConfirmed.',dp);
    end
    RespondingConfirmed = RespLoad.RespondingConfirmed;

    fprintf('Responding-channel conditions available: %d\n',numel(RespondingConfirmed));
    fprintf('\nComputing actual vs predicted spike counts...\n');

    for si = 1:Paired.nSets

        stim_channels = Paired.uniqueComb(si,Paired.uniqueComb(si,:) > 0);

        if numel(stim_channels) ~= 2
            continue;   % not a pair - nothing to compare here
        end

        [dist_um,dist_note] = ChnPairDistance(stim_channels(1),stim_channels(2),Electrode_Type);

        if isnan(dist_um)
            warning('Set %d (%s): distance unknown (%s) - excluded from the plot.', ...
                si,format_order(stim_channels),dist_note);
        end

        for ai = 1:numel(Paired.Amps)

            amp_value = Paired.Amps(ai);

            for pi = 1:numel(Paired.PTDs_ms)

                ptd_value = Paired.PTDs_ms(pi);

                trial_ids = find(Paired.combClass == si & Paired.ampIdx == ai & Paired.ptdIdx == pi);
                trial_ids = setdiff(trial_ids,Paired.bad_trial_ids);

                if isempty(trial_ids)
                    continue;
                end

                %% ---- find this condition's saved responding-channel list ----

                match_idx = [];
                for k = 1:numel(RespondingConfirmed)
                    if isequal(RespondingConfirmed(k).stimChannels,stim_channels) && ...
                            abs(RespondingConfirmed(k).amp-amp_value) < Amp_Tol_uA && ...
                            abs(RespondingConfirmed(k).ptd_ms-ptd_value) < PTD_Tol_ms
                        match_idx = k;
                        break;
                    end
                end

                if isempty(match_idx)
                    warning('Set %d (%s) | Amp %g uA | PTD %g ms: no saved responding-channel entry - skipped.', ...
                        si,format_order(stim_channels),amp_value,ptd_value);
                    continue;
                end

                responding_channels = RespondingConfirmed(match_idx).channels;

                if isempty(responding_channels)
                    warning('Set %d (%s) | Amp %g uA | PTD %g ms: responding-channel list is empty - skipped.', ...
                        si,format_order(stim_channels),amp_value,ptd_value);
                    continue;
                end

                %% ---------------------- ACTUAL (paired) -----------------------

                trig_ms_paired = Paired.trig_ms(trial_ids);

                per_channel_actual = nan(numel(responding_channels),1);
                for c = 1:numel(responding_channels)
                    ch = responding_channels(c);
                    if ch < 1 || ch > Paired.nMapPos
                        warning('Set %d: responding channel %d is out of range for the paired dataset (1-%d) - skipped for actual.', ...
                            si,ch,Paired.nMapPos);
                        continue;
                    end
                    per_channel_actual(c) = mean_corrected_count( ...
                        Paired.sp_by_mappos{ch},trig_ms_paired,baseline_win_ms,analysis_win_ms);
                end
                actual_value = mean(per_channel_actual,'omitnan');

                %% --------------------- PREDICTED (single + single) -------------

                [ssA,aiA] = find_single_condition(Single,stim_channels(1),amp_value,Amp_Tol_uA);
                [ssB,aiB] = find_single_condition(Single,stim_channels(2),amp_value,Amp_Tol_uA);

                predicted_value = NaN;

                if isempty(ssA) || isempty(ssB)

                    if isempty(ssA)
                        warning('Set %d (%s) | Amp %g uA: no single-pulse match for Ch%d - predicted left NaN.', ...
                            si,format_order(stim_channels),amp_value,stim_channels(1));
                    end
                    if isempty(ssB)
                        warning('Set %d (%s) | Amp %g uA: no single-pulse match for Ch%d - predicted left NaN.', ...
                            si,format_order(stim_channels),amp_value,stim_channels(2));
                    end

                else
                    trials_A = find(Single.combClass == ssA & Single.ampIdx == aiA);
                    trials_A = setdiff(trials_A,Single.bad_trial_ids);

                    trials_B = find(Single.combClass == ssB & Single.ampIdx == aiB);
                    trials_B = setdiff(trials_B,Single.bad_trial_ids);

                    if isempty(trials_A) || isempty(trials_B)

                        if isempty(trials_A)
                            warning('Set %d (%s) | Amp %g uA: no usable single-pulse trials for Ch%d (all excluded/none found) - predicted left NaN.', ...
                                si,format_order(stim_channels),amp_value,stim_channels(1));
                        end
                        if isempty(trials_B)
                            warning('Set %d (%s) | Amp %g uA: no usable single-pulse trials for Ch%d (all excluded/none found) - predicted left NaN.', ...
                                si,format_order(stim_channels),amp_value,stim_channels(2));
                        end

                    else
                        trig_ms_A = Single.trig_ms(trials_A);
                        trig_ms_B = Single.trig_ms(trials_B);

                        per_channel_A = nan(numel(responding_channels),1);
                        per_channel_B = nan(numel(responding_channels),1);

                        for c = 1:numel(responding_channels)
                            ch = responding_channels(c);
                            if ch < 1 || ch > Single.nMapPos
                                warning('Set %d: responding channel %d is out of range for the single dataset (1-%d) - skipped for predicted.', ...
                                    si,ch,Single.nMapPos);
                                continue;
                            end
                            per_channel_A(c) = mean_corrected_count( ...
                                Single.sp_by_mappos{ch},trig_ms_A,baseline_win_ms,analysis_win_ms);
                            per_channel_B(c) = mean_corrected_count( ...
                                Single.sp_by_mappos{ch},trig_ms_B,baseline_win_ms,analysis_win_ms);
                        end

                        singleA_value = mean(per_channel_A,'omitnan');
                        singleB_value = mean(per_channel_B,'omitnan');
                        predicted_value = singleA_value+singleB_value;
                    end
                end

                %% ------------------------- STORE ROW ---------------------------

                Result_Set(end+1,1)          = si; %#ok<AGROW>
                Result_StimLabel{end+1,1}    = format_order(stim_channels); %#ok<AGROW>
                Result_Distance(end+1,1)     = dist_um; %#ok<AGROW>
                Result_Amp(end+1,1)          = amp_value; %#ok<AGROW>
                Result_PTD(end+1,1)          = ptd_value; %#ok<AGROW>
                Result_Actual(end+1,1)       = actual_value; %#ok<AGROW>
                Result_Predicted(end+1,1)    = predicted_value; %#ok<AGROW>
                Result_NChannels(end+1,1)    = numel(responding_channels); %#ok<AGROW>
                Result_NTrials(end+1,1)      = numel(trial_ids); %#ok<AGROW>
                Result_DatasetLabel{end+1,1} = dataset_label; %#ok<AGROW>

                fprintf('  Set %d (%s) | Amp %g uA | PTD %g ms | Dist %s | Actual %.3g | Predicted %.3g\n', ...
                    si,format_order(stim_channels),amp_value,ptd_value, ...
                    dist_text(dist_um,dist_note),actual_value,predicted_value);
            end
        end
    end

end   % end of dataset-pair loop (dp)

Results = table(Result_Set,Result_StimLabel,Result_Distance,Result_Amp,Result_PTD, ...
    Result_Actual,Result_Predicted,Result_NChannels,Result_NTrials,Result_DatasetLabel, ...
    'VariableNames',{'Set','StimLabel','Distance_um','Amp_uA','PTD_ms', ...
    'Actual','Predicted','N_Channels','N_PairedTrials','Dataset_Label'});

fprintf('\nTotal condition points computed: %d\n',height(Results));

%% =============================== PLOT ==================================

plot_rows = ~isnan(Results.Distance_um);
if any(~plot_rows)
    fprintf('Excluding %d point(s) with unknown distance from the plot.\n',sum(~plot_rows));
end
PlotData = Results(plot_rows,:);

%% ------ average across stim orders AND dataset-pairs (same distance/amp/PTD) ------
% A given distance/amp/PTD point can appear more than once - either
% because both stimulation orders (ChA->ChB and ChB->ChA) were tested, or
% because more than one dataset-pair happened to test the same distance.
% Either way, those rows are averaged together into a single point per
% distance/amp/PTD for the plot.

[group_keys,~,group_idx] = unique( ...
    [PlotData.Distance_um,PlotData.Amp_uA,PlotData.PTD_ms],'rows');

nGroups = size(group_keys,1);

Avg_Distance = group_keys(:,1);
Avg_Amp      = group_keys(:,2);
Avg_PTD      = group_keys(:,3);
Avg_Actual    = nan(nGroups,1);
Avg_Predicted = nan(nGroups,1);
Avg_NSets     = zeros(nGroups,1);

fprintf('\nAveraging across stim orders / dataset-pairs for the plot:\n');

for g = 1:nGroups
    rows_in_group = group_idx == g;
    Avg_Actual(g)    = mean(PlotData.Actual(rows_in_group),'omitnan');
    Avg_Predicted(g) = mean(PlotData.Predicted(rows_in_group),'omitnan');
    Avg_NSets(g)      = sum(rows_in_group);

    if Avg_NSets(g) > 1
        row_indices = find(rows_in_group);
        entry_labels = cell(numel(row_indices),1);
        for j = 1:numel(row_indices)
            r = row_indices(j);
            entry_labels{j} = sprintf('%s [%s]',PlotData.StimLabel{r},PlotData.Dataset_Label{r});
        end
        fprintf('  Dist %g um | Amp %g uA | PTD %g ms: averaged %d points (%s)\n', ...
            Avg_Distance(g),Avg_Amp(g),Avg_PTD(g),Avg_NSets(g),strjoin(entry_labels,', '));
    end
end

PlotData = table(Avg_Distance,Avg_Amp,Avg_PTD,Avg_Actual,Avg_Predicted,Avg_NSets, ...
    'VariableNames',{'Distance_um','Amp_uA','PTD_ms','Actual','Predicted','N_Sets_Averaged'});

unique_PTDs = unique(PlotData.PTD_ms);

for p = 1:numel(unique_PTDs)

    this_ptd = unique_PTDs(p);
    subset = PlotData(abs(PlotData.PTD_ms-this_ptd) < PTD_Tol_ms,:);

    if isempty(subset)
        continue;
    end

    unique_amps = unique(subset.Amp_uA);
    amp_colors = lines(numel(unique_amps));

    figure('Color','w','Position',[100 100 800 600]);
    hold on;

    for ai = 1:numel(unique_amps)

        amp_value = unique_amps(ai);
        sub2 = subset(abs(subset.Amp_uA-amp_value) < Amp_Tol_uA,:);
        sub2 = sortrows(sub2,'Distance_um');

        col = amp_colors(ai,:);

        valid_actual = ~isnan(sub2.Actual);
        if any(valid_actual)
            plot(sub2.Distance_um(valid_actual),sub2.Actual(valid_actual),'-o', ...
                'Color',col,'LineWidth',2,'MarkerFaceColor',col,'MarkerSize',7, ...
                'DisplayName',sprintf('%g uA (actual)',amp_value));
        end

        valid_pred = ~isnan(sub2.Predicted);
        if any(valid_pred)
            plot(sub2.Distance_um(valid_pred),sub2.Predicted(valid_pred),'--s', ...
                'Color',col,'LineWidth',2,'MarkerFaceColor','w','MarkerSize',7, ...
                'DisplayName',sprintf('%g uA (predicted)',amp_value));
        end
    end

    xlabel('Distance (\mum)','FontWeight','bold','FontSize',12);
    ylabel('Mean baseline-corrected spike count (per trial, per channel)','FontWeight','bold','FontSize',12);
    title(sprintf('PTD %g ms',this_ptd),'FontWeight','bold','FontSize',14);
    box off;
    lgd = legend('Location','best','Box','off');
    title(lgd,'Amplitude');
end

fprintf('\n============================================================\n');
fprintf('DONE (plot only - no files saved)\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function validate_window(window_ms,variable_name)
if ~isnumeric(window_ms) || numel(window_ms) ~= 2 || any(~isfinite(window_ms)) || window_ms(2) <= window_ms(1)
    error('%s must be [start end], with end > start.',variable_name);
end
end

function files = remove_metadata_and_backups(files)
if isempty(files), return; end
keep = true(size(files));
for file_index = 1:numel(files)
    file_name = files(file_index).name;
    if startsWith(file_name,'._') || contains(file_name,'BACKUP','IgnoreCase',true)
        keep(file_index) = false;
    end
end
files = files(keep);
end

function label = format_order(stim_channels)
stim_channels = stim_channels(stim_channels > 0);
if isempty(stim_channels)
    label = '(none)';
elseif numel(stim_channels) == 1
    label = sprintf('Ch%d',stim_channels(1));
else
    labels = arrayfun(@(c) sprintf('Ch%d',c),stim_channels,'UniformOutput',false);
    label = strjoin(labels,' -> ');
end
end

function output = dist_text(dist_um,dist_note)
if isnan(dist_um)
    output = sprintf('n/a (%s)',dist_note);
else
    output = sprintf('%g um',dist_um);
end
end

function value = mean_corrected_count(spike_times,trig_ms_list,baseline_win_ms,analysis_win_ms)
% Mean baseline-corrected analysis-window spike count across the given
% trigger times (ms), for one channel's spike times (ms).

nTrials = numel(trig_ms_list);
if nTrials == 0
    value = NaN;
    return;
end

if isempty(spike_times)
    value = 0;
    return;
end

baseline_duration_ms = diff(baseline_win_ms);
analysis_duration_ms = diff(analysis_win_ms);

total_corrected = 0;
for k = 1:nTrials
    t0 = trig_ms_list(k);

    baseline_count = sum(spike_times >= t0+baseline_win_ms(1) & spike_times < t0+baseline_win_ms(2));
    analysis_count = sum(spike_times >= t0+analysis_win_ms(1) & spike_times < t0+analysis_win_ms(2));

    expected_baseline = baseline_count*(analysis_duration_ms/baseline_duration_ms);
    total_corrected = total_corrected+(analysis_count-expected_baseline);
end

value = total_corrected/nTrials;
end

function [matched_set,matched_amp_idx] = find_single_condition(Single,target_channel,target_amp,amp_tol)
% Find the single-channel stimulation condition in Single matching
% target_channel, at an amplitude matching target_amp. Returns empty for
% either output if no match exists.

matched_set = [];
matched_amp_idx = [];

for ss = 1:Single.nSets
    stim_channels = Single.uniqueComb(ss,Single.uniqueComb(ss,:) > 0);
    if numel(stim_channels) == 1 && stim_channels(1) == target_channel
        ai = find(abs(Single.Amps-target_amp) < amp_tol,1);
        if ~isempty(ai)
            matched_set = ss;
            matched_amp_idx = ai;
        end
        return;
    end
end
end

function DS = load_linearity_dataset(spike_file,dataset_folder,Electrode_Type,FS_fallback)
% Load one dataset (paired or single) for the linearity analysis: spike
% times by map position, triggers in ms, decoded amp/PTD/stim-set
% conditions, and that dataset's own saved bad-trial exclusion list (if
% any). Mirrors the loading section of TrialSpikeCounts_v2.m.

%% -------- spike times --------

spike_vars_present = who('-file',spike_file);

if ismember('sp_clipped',spike_vars_present)
    SpikeLoad = load(spike_file,'sp_clipped');
    sp_raw = SpikeLoad.sp_clipped;
elseif ismember('sp_waveforms',spike_vars_present)
    SpikeLoad = load(spike_file,'sp_waveforms');
    sp_raw = SpikeLoad.sp_waveforms;
else
    error('Spike file contains neither sp_clipped nor sp_waveforms:\n%s',spike_file);
end

nCh = numel(sp_raw);

sp = cell(nCh,1);
for ch = 1:nCh
    cellData = sp_raw{ch};
    if isempty(cellData)
        sp{ch} = [];
    else
        spike_times = double(cellData(:,1));
        spike_times = spike_times(isfinite(spike_times));
        if ~issorted(spike_times)
            spike_times = sort(spike_times);
        end
        sp{ch} = spike_times;
    end
end
clear sp_raw SpikeLoad;

if ismember('fs',spike_vars_present)
    FsLoad = load(spike_file,'fs');
    FS = double(FsLoad.fs);
elseif ismember('FS',spike_vars_present)
    FsLoad = load(spike_file,'FS');
    FS = double(FsLoad.FS);
else
    warning('fs/FS not found in spike file - using fallback FS = %g Hz.',FS_fallback);
    FS = FS_fallback;
end

fprintf('Physical channels in spike file: %d\n',nCh);
fprintf('Sampling rate: %g Hz\n',FS);

%% -------- electrode mapping --------

map_nums_plus = double(ChnMap(Electrode_Type,nCh));
map_nums_plus = map_nums_plus(:);

if isempty(map_nums_plus)
    error('ChnMap(%d,%d) returned an empty channel map.',Electrode_Type,nCh);
end

nMapPos = numel(map_nums_plus);
fprintf('Map positions: %d\n',nMapPos);

valid_channel_mapping = isfinite(map_nums_plus) & map_nums_plus >= 1 & ...
    map_nums_plus <= nCh & fix(map_nums_plus) == map_nums_plus;

sp_by_mappos = cell(nMapPos,1);
for channel_position = 1:nMapPos
    if ~valid_channel_mapping(channel_position)
        continue;
    end
    phys_ch = map_nums_plus(channel_position);
    sp_by_mappos{channel_position} = sp{phys_ch};
end
clear sp;

%% -------- triggers --------

cd(dataset_folder);

if isempty(dir('*.trig.dat'))
    current_folder = pwd;
    cleanTrig_sabquick;
    cd(current_folder);
end

trigger_files = dir('*.trig.dat');
if isempty(trigger_files)
    error('No *.trig.dat file could be found or created in:\n%s',dataset_folder);
end

trig = double(loadTrig(0));
trig = trig(:);
nTrig = numel(trig);

%% -------- experiment parameters --------

experiment_files = dir('*_exp_datafile_*.mat');
experiment_files = remove_metadata_and_backups(experiment_files);

if isempty(experiment_files)
    error('No *_exp_datafile_*.mat file was found in:\n%s',dataset_folder);
end
if numel(experiment_files) > 1
    error('Multiple current experiment files were found in:\n%s',dataset_folder);
end

experiment_file = experiment_files(1).name;

ExpLoad = load(experiment_file,'StimParams','simultaneous_stim','E_MAP','n_Trials');
required_variables = {'StimParams','simultaneous_stim','E_MAP','n_Trials'};
for variable_index = 1:numel(required_variables)
    if ~isfield(ExpLoad,required_variables{variable_index})
        error('Experiment file is missing variable: %s',required_variables{variable_index});
    end
end

StimParams = ExpLoad.StimParams;
sim_stim = double(ExpLoad.simultaneous_stim);
E_MAP = ExpLoad.E_MAP;
n_Trials = double(ExpLoad.n_Trials);

if nTrig < n_Trials
    error('Only %d triggers were loaded for %d trials.',nTrig,n_Trials);
elseif nTrig > n_Trials
    trig = trig(1:n_Trials);
end

trig_ms = trig(:)/FS*1000;

fprintf('Experiment file: %s\n',experiment_file);
fprintf('Stimulation events per trial (simultaneous_stim): %d\n',sim_stim);
fprintf('Trials: %d\n',n_Trials);

%% -------- decode amplitudes --------

amplitudes_all = cell2mat(StimParams(2:end,16));
amplitudes_all = double(amplitudes_all(:));

first_event_rows = 1:sim_stim:numel(amplitudes_all);
trialAmps = amplitudes_all(first_event_rows);
trialAmps = trialAmps(1:n_Trials);
trialAmps(trialAmps == -1) = 0;

[Amps,~,ampIdx] = unique(trialAmps(:));

%% -------- decode PTDs --------

if sim_stim >= 2
    PTD_all_us = cell2mat(StimParams(2:end,6));
    PTD_all_us = double(PTD_all_us(:));

    second_event_rows = 2:sim_stim:numel(PTD_all_us);
    trialPTD_us = PTD_all_us(second_event_rows);
    trialPTD_us = trialPTD_us(1:n_Trials);

    [PTDs,~,ptdIdx] = unique(trialPTD_us(:));
    PTDs_ms = PTDs/1000;
else
    PTDs_ms = 0;
    ptdIdx = ones(n_Trials,1);
end

%% -------- decode stimulation channel orders --------

stimNames = StimParams(2:end,1);
stimNames = stimNames(1:n_Trials*sim_stim);

[isMapped,idx_all] = ismember(stimNames,E_MAP(2:end));
if any(~isMapped)
    error('%d stimulation entries could not be mapped through E_MAP.',sum(~isMapped));
end

stimSeq = zeros(n_Trials,sim_stim);
for trial_id = 1:n_Trials
    rows_this_trial = (trial_id-1)*sim_stim+(1:sim_stim);
    mapped_channels = idx_all(rows_this_trial);
    mapped_channels = mapped_channels(mapped_channels > 0);
    stimSeq(trial_id,1:numel(mapped_channels)) = mapped_channels(:).';
end

[uniqueComb,~,combClass] = unique(stimSeq,'rows','stable');
nSets = size(uniqueComb,1);

fprintf('Ordered stimulation sets: %d\n',nSets);

%% -------- bad-trial exclusion list (own file, if any) --------

base_name = regexprep(experiment_file,'_exp_datafile_.*$','');
badtrials_path = fullfile(dataset_folder,sprintf('%s_BadTrials.mat',base_name));

if isfile(badtrials_path)
    BadLoad = load(badtrials_path,'AllBadAbsoluteTrialIDs');
    if isfield(BadLoad,'AllBadAbsoluteTrialIDs')
        bad_trial_ids = double(BadLoad.AllBadAbsoluteTrialIDs(:));
        fprintf('Bad-trial file found: %s (%d trials excluded)\n',badtrials_path,numel(bad_trial_ids));
    else
        warning('Bad-trial file found but missing AllBadAbsoluteTrialIDs - no trials excluded:\n%s',badtrials_path);
        bad_trial_ids = [];
    end
else
    warning('No bad-trial file found - no trials excluded for this dataset:\n%s',dataset_folder);
    bad_trial_ids = [];
end

%% -------- package --------

DS = struct();
DS.dataset_folder = dataset_folder;
DS.experiment_file = experiment_file;
DS.base_name = base_name;
DS.FS = FS;
DS.nMapPos = nMapPos;
DS.sp_by_mappos = sp_by_mappos;
DS.trig_ms = trig_ms;
DS.n_Trials = n_Trials;
DS.sim_stim = sim_stim;
DS.Amps = Amps;
DS.ampIdx = ampIdx;
DS.PTDs_ms = PTDs_ms;
DS.ptdIdx = ptdIdx;
DS.uniqueComb = uniqueComb;
DS.combClass = combClass;
DS.nSets = nSets;
DS.bad_trial_ids = bad_trial_ids;
end