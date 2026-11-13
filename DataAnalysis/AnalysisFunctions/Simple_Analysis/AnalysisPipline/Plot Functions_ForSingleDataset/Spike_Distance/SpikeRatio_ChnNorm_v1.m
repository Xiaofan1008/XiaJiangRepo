%% ========================================================================
% SPIKE COUNT vs DISTANCE - LINEARITY (FIXED WINDOW, NO SHIFT)  [v7: + channel-max normalisation]
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
%   The responding-channel file is NOT used. Every condition (Set x Amp x
%   PTD) counts the SAME fixed channel list, Channel_Range (below), for
%   BOTH the actual count and the two single-channel components of the
%   predicted count. Channel_Range must not be empty.
%
% METRIC
%   Baseline-corrected spike count in a FIXED analysis_win_ms window
%   (default [0,20] ms, separate from the [2,40] ms window used by the
%   bad-trial QC pipeline), averaged over trials, then averaged over the
%   Channel_Range channels:
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
%   single dataset only. Apply_BadTrial_Exclusion (below) can turn this
%   off entirely (both paired and single, every dataset-pair).
%
% LINEARITY INDEX (extra figures + table)
%   LI = (Actual - Predicted)/(Actual + Predicted), from the pooled means
%   used in the population plot. 0 = linear (actual = predicted),
%   < 0 = sub-additive, > 0 = supra-additive, bounded to [-1,1] when both
%   are positive. This removes the electrode-strength differences that make
%   the raw counts move with distance. The table also prints Ratio =
%   Actual/Predicted (NaN where Predicted <= 0).
%
% NORMALISED RATIO (extra figures)
%   Ratio = Actual / Predicted per distance/amp/PTD (pooled means), so every
%   pair is its own control and predicted becomes 1 for all distances.
%   A second figure plots the ratio against predicted size to check whether
%   electrode strength still leaks into the ratio (e.g. saturation): a
%   ratio that falls as predicted grows means the normalisation is not clean.
%
% OUTPUT
%   Plot only for now (one figure per AMPLITUDE; color = PTD, solid
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

DatasetPairs(1).simseq_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_SimSeq1/Xia_Linearity_S1_100_400um_SimSeq1.sp_xia_stitched.mat';
DatasetPairs(1).simseq_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_SimSeq1';
DatasetPairs(1).single_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1/Xia_Linearity_S1_100_400um_Single1.sp_xia_filtered.mat';
DatasetPairs(1).single_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1';

DatasetPairs(2).simseq_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_SimSeq1/Xia_Linearity_S1_500_700um_SimSeq1.sp_xia_stitched.mat';
DatasetPairs(2).simseq_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_SimSeq1';
DatasetPairs(2).single_spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1/Xia_Linearity_S1_500_700um_Single1.sp_xia_filtered.mat';
DatasetPairs(2).single_dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1';

Channel_Range = [1:64]; % Channels to count (REQUIRED - used for every condition;
                        % no responding-channel file is read)
Apply_BadTrial_Exclusion = true; % load bad trial
Channels_To_Plot = [];
Plot_PerChannel_Diagnostics = false; % per chn plot
Plot_LinearityIndex = true;          % extra figures: linearity index per amplitude
LI_Min_Denominator  = 0.1;           % index set to NaN where Actual+Predicted is not above this
Plot_NormalizedRatio = true;         % extra figures: actual / predicted (predicted = 1), and ratio vs predicted size
Ratio_Min_Predicted  = 0.2;          % ratio set to NaN where Predicted is not above this (avoids blow-ups near 0)

% ---- per-channel normalisation by each channel's own maximum SINGLE response ----
Plot_ChannelMaxNorm = true;          % extra figures: channel-normalised ratio + actual-vs-predicted curve
Norm_Max_Scope      = 'all';         % 'all' = max over all single responses of the channel (all amps, all stim electrodes)
                                     % 'amp' = max separately for each amplitude
Norm_Min_Max        = 0.5;           % channels whose max single response (spikes/trial) is below this are dropped
Norm_Min_PredMean   = 0.05;          % ratio set to NaN where mean normalised predicted is not above this
IO_Bin_Edges        = [0 0.1 0.25 0.5 0.75 1 1.5 Inf];  % bins of normalised predicted (1 = channel's max single response)
IO_Min_Points       = 10;            % bins with fewer channel-points are not drawn
Near_Max_um         = 300;           % "near" group: distance <= this
Far_Min_um          = 500;           % "far"  group: distance >= this

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

if isempty(Channel_Range)
    error('Channel_Range is empty - set the channels to count (e.g. 1:64). The responding-channel file is no longer used.');
end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('SPIKE COUNT vs DISTANCE - LINEARITY (FIXED WINDOW, NO SHIFT)\n');
fprintf('============================================================\n');
fprintf('Dataset-pairs: %d\n',numel(DatasetPairs));
fprintf('Baseline window: [%g,%g) ms\n',baseline_win_ms);
fprintf('Analysis window: [%g,%g) ms\n',analysis_win_ms);
fprintf('Channel_Range (used for every condition, no responding-channel file): %s\n', ...
    num2str(unique(double(Channel_Range(:).'))));
if Apply_BadTrial_Exclusion
    fprintf('Apply_BadTrial_Exclusion: true (each dataset''s own saved bad trials excluded)\n');
else
    fprintf('Apply_BadTrial_Exclusion: FALSE - no trials excluded (bad-trial files not loaded)\n');
end

%% =========================== MAIN COMPUTATION ==========================

Result_Set              = [];
Result_StimLabel        = {};
Result_Distance         = [];
Result_Amp              = [];
Result_PTD              = [];
Result_Actual           = [];
Result_Predicted        = [];
Result_NChannels        = [];
Result_NTrials          = [];
Result_DatasetLabel     = {};
Result_ActualPerChn     = {};
Result_PredictedPerChn  = {};
Result_SingleAPerChn    = {};
Result_SingleBPerChn    = {};
Result_ChannelList      = {};

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
        Electrode_Type,FS_fallback,Apply_BadTrial_Exclusion);

    fprintf('\n--- Loading single dataset ---\n');
    Single = load_linearity_dataset(entry.single_spike_file,entry.single_dataset_folder, ...
        Electrode_Type,FS_fallback,Apply_BadTrial_Exclusion);

    dataset_label = Paired.base_name;

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

                responding_channels = unique(double(Channel_Range(:).'));

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
                per_channel_predicted = nan(numel(responding_channels),1);
                per_channel_singleA   = nan(numel(responding_channels),1);
                per_channel_singleB   = nan(numel(responding_channels),1);

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
                        per_channel_predicted = per_channel_A+per_channel_B;
                        per_channel_singleA = per_channel_A;
                        per_channel_singleB = per_channel_B;
                    end
                end

                %% ------------------------- STORE ROW ---------------------------

                Result_Set(end+1,1)             = si; %#ok<AGROW>
                Result_StimLabel{end+1,1}       = format_order(stim_channels); %#ok<AGROW>
                Result_Distance(end+1,1)        = dist_um; %#ok<AGROW>
                Result_Amp(end+1,1)             = amp_value; %#ok<AGROW>
                Result_PTD(end+1,1)             = ptd_value; %#ok<AGROW>
                Result_Actual(end+1,1)          = actual_value; %#ok<AGROW>
                Result_Predicted(end+1,1)       = predicted_value; %#ok<AGROW>
                Result_NChannels(end+1,1)       = numel(responding_channels); %#ok<AGROW>
                Result_NTrials(end+1,1)         = numel(trial_ids); %#ok<AGROW>
                Result_DatasetLabel{end+1,1}    = dataset_label; %#ok<AGROW>
                Result_ActualPerChn{end+1,1}    = per_channel_actual(:); %#ok<AGROW>
                Result_PredictedPerChn{end+1,1} = per_channel_predicted(:); %#ok<AGROW>
                Result_SingleAPerChn{end+1,1}   = per_channel_singleA(:); %#ok<AGROW>
                Result_SingleBPerChn{end+1,1}   = per_channel_singleB(:); %#ok<AGROW>
                Result_ChannelList{end+1,1}     = responding_channels(:); %#ok<AGROW>

                fprintf('  Set %d (%s) | Amp %g uA | PTD %g ms | Dist %s | Actual %.3g | Predicted %.3g\n', ...
                    si,format_order(stim_channels),amp_value,ptd_value, ...
                    dist_text(dist_um,dist_note),actual_value,predicted_value);
            end
        end
    end

end   % end of dataset-pair loop (dp)

Results = table(Result_Set,Result_StimLabel,Result_Distance,Result_Amp,Result_PTD, ...
    Result_Actual,Result_Predicted,Result_NChannels,Result_NTrials,Result_DatasetLabel, ...
    Result_ActualPerChn,Result_PredictedPerChn,Result_ChannelList, ...
    Result_SingleAPerChn,Result_SingleBPerChn, ...
    'VariableNames',{'Set','StimLabel','Distance_um','Amp_uA','PTD_ms', ...
    'Actual','Predicted','N_Channels','N_PairedTrials','Dataset_Label', ...
    'Actual_PerChannel','Predicted_PerChannel','Channel_List', ...
    'SingleA_PerChannel','SingleB_PerChannel'});

fprintf('\nTotal condition points computed: %d\n',height(Results));

%% =============================== PLOT ==================================

plot_rows = ~isnan(Results.Distance_um);
if any(~plot_rows)
    fprintf('Excluding %d point(s) with unknown distance from the plot.\n',sum(~plot_rows));
end
PlotData = Results(plot_rows,:);

%% ------ pool across stim orders AND dataset-pairs (same distance/amp/PTD) ------
% A given distance/amp/PTD point can appear more than once - either
% because both stimulation orders (ChA->ChB and ChB->ChA) were tested, or
% because more than one dataset-pair happened to test the same distance.
% Every per-channel value contributing to that point (across however many
% Sets/datasets land on it) is pooled into one distribution, and the plot
% point is the mean of that pooled distribution with the SD as its error
% bar - this is the across-CHANNEL spread (does one channel dominate?),
% not trial-to-trial noise.

[group_keys,~,group_idx] = unique( ...
    [PlotData.Distance_um,PlotData.Amp_uA,PlotData.PTD_ms],'rows');

nGroups = size(group_keys,1);

Avg_Distance    = group_keys(:,1);
Avg_Amp         = group_keys(:,2);
Avg_PTD         = group_keys(:,3);
Avg_Actual      = nan(nGroups,1);
Avg_ActualSD    = nan(nGroups,1);
Avg_Predicted   = nan(nGroups,1);
Avg_PredictedSD = nan(nGroups,1);
Avg_NSets       = zeros(nGroups,1);
Avg_NChannels   = zeros(nGroups,1);

fprintf('\nPooling per-channel values across stim orders / dataset-pairs for the plot:\n');

for g = 1:nGroups
    rows_in_group = group_idx == g;
    row_indices = find(rows_in_group);

    pooled_actual    = vertcat(PlotData.Actual_PerChannel{row_indices});
    pooled_predicted = vertcat(PlotData.Predicted_PerChannel{row_indices});

    Avg_Actual(g)      = mean(pooled_actual,'omitnan');
    Avg_ActualSD(g)     = std(pooled_actual,'omitnan');
    Avg_Predicted(g)   = mean(pooled_predicted,'omitnan');
    Avg_PredictedSD(g)  = std(pooled_predicted,'omitnan');
    Avg_NSets(g)        = numel(row_indices);
    Avg_NChannels(g)     = sum(~isnan(pooled_actual));

    if Avg_NSets(g) > 1
        entry_labels = cell(numel(row_indices),1);
        for j = 1:numel(row_indices)
            r = row_indices(j);
            entry_labels{j} = sprintf('%s [%s]',PlotData.StimLabel{r},PlotData.Dataset_Label{r});
        end
        fprintf('  Dist %g um | Amp %g uA | PTD %g ms: pooled %d points, %d channels total (%s)\n', ...
            Avg_Distance(g),Avg_Amp(g),Avg_PTD(g),Avg_NSets(g),Avg_NChannels(g),strjoin(entry_labels,', '));
    end
end

PlotData = table(Avg_Distance,Avg_Amp,Avg_PTD,Avg_Actual,Avg_ActualSD, ...
    Avg_Predicted,Avg_PredictedSD,Avg_NSets,Avg_NChannels, ...
    'VariableNames',{'Distance_um','Amp_uA','PTD_ms','Actual','Actual_SD', ...
    'Predicted','Predicted_SD','N_Sets_Pooled','N_Channels_Pooled'});

unique_amps = unique(PlotData.Amp_uA);
unique_PTDs = unique(PlotData.PTD_ms);
ptd_colors  = lines(numel(unique_PTDs));

% One figure per AMPLITUDE. Each PTD is one colored line (PTD 0 = simultaneous,
% PTD > 0 = sequential). Predicted (single A + single B) should be the same
% for every PTD, so it is drawn once as a black dashed line when it is
% identical across PTDs; otherwise it is drawn per PTD in that PTD's color.

for a = 1:numel(unique_amps)

    this_amp = unique_amps(a);
    subset = PlotData(abs(PlotData.Amp_uA-this_amp) < Amp_Tol_uA,:);

    if isempty(subset)
        continue;
    end

    % ---- is predicted identical across PTDs? ----
    ref_sub = [];
    pred_identical = true;
    for p = 1:numel(unique_PTDs)
        s = subset(abs(subset.PTD_ms-unique_PTDs(p)) < PTD_Tol_ms,:);
        s = sortrows(s,'Distance_um');
        if isempty(s), continue; end
        if isempty(ref_sub)
            ref_sub = s;
        else
            if height(s) ~= height(ref_sub) || ...
                    any(abs(s.Distance_um-ref_sub.Distance_um) > 1e-9) || ...
                    any(isnan(s.Predicted) ~= isnan(ref_sub.Predicted)) || ...
                    any(abs(s.Predicted-ref_sub.Predicted) > 1e-9)
                pred_identical = false;
                break;
            end
        end
    end

    figure('Color','w','Position',[100 100 800 600]);
    hold on;

    for p = 1:numel(unique_PTDs)

        this_ptd = unique_PTDs(p);
        sub2 = subset(abs(subset.PTD_ms-this_ptd) < PTD_Tol_ms,:);
        if isempty(sub2), continue; end
        sub2 = sortrows(sub2,'Distance_um');

        col = ptd_colors(p,:);

        if abs(this_ptd) < PTD_Tol_ms
            cond_label = 'Simultaneous (PTD 0 ms)';
        else
            cond_label = sprintf('Sequential (PTD %g ms)',this_ptd);
        end

        valid_actual = ~isnan(sub2.Actual);
        if any(valid_actual)
            plot(sub2.Distance_um(valid_actual),sub2.Actual(valid_actual),'-o', ...
                'Color',col,'LineWidth',2,'MarkerFaceColor',col,'MarkerSize',7, ...
                'DisplayName',[cond_label ' - actual']);
        end

        if ~pred_identical
            valid_pred = ~isnan(sub2.Predicted);
            if any(valid_pred)
                plot(sub2.Distance_um(valid_pred),sub2.Predicted(valid_pred),'--s', ...
                    'Color',col,'LineWidth',2,'MarkerFaceColor','w','MarkerSize',7, ...
                    'DisplayName',[cond_label ' - predicted']);
            end
        end
    end

    if pred_identical && ~isempty(ref_sub)
        valid_pred = ~isnan(ref_sub.Predicted);
        if any(valid_pred)
            plot(ref_sub.Distance_um(valid_pred),ref_sub.Predicted(valid_pred),'--s', ...
                'Color','k','LineWidth',2,'MarkerFaceColor','w','MarkerSize',7, ...
                'DisplayName','Predicted (single A + single B)');
        end
    end

    xlabel('Distance (\mum)','FontWeight','bold','FontSize',12);
    ylabel('Mean baseline-corrected spike count (per trial, per channel)','FontWeight','bold','FontSize',12);
    title(sprintf('Amplitude %g \\muA',this_amp),'FontWeight','bold','FontSize',14);
    box off;
    lgd = legend('Location','best','Box','off');
    title(lgd,'Condition');
end

%% ======================= LINEARITY INDEX PLOTS ==========================
% One figure per amplitude, one line per PTD: (Actual-Predicted)/(Actual+Predicted)
% against distance. Dotted zero line = perfectly linear.

if Plot_LinearityIndex

    denom = PlotData.Actual+PlotData.Predicted;
    LI = (PlotData.Actual-PlotData.Predicted)./denom;
    LI(~(denom > LI_Min_Denominator)) = NaN;

    ratio = PlotData.Actual./PlotData.Predicted;
    ratio(~(PlotData.Predicted > 0)) = NaN;

    PlotData.Ratio = ratio;
    PlotData.LinearityIndex = LI;

    fprintf('\n=========== ACTUAL vs PREDICTED (pooled), RATIO AND LINEARITY INDEX ===========\n');
    disp(sortrows(PlotData(:,{'Amp_uA','PTD_ms','Distance_um','Actual','Predicted','Ratio','LinearityIndex'}), ...
        {'Amp_uA','PTD_ms','Distance_um'}));

    for a = 1:numel(unique_amps)

        this_amp = unique_amps(a);
        subset = PlotData(abs(PlotData.Amp_uA-this_amp) < Amp_Tol_uA,:);
        if isempty(subset), continue; end

        figure('Color','w','Position',[120 120 800 600]);
        hold on;

        for p = 1:numel(unique_PTDs)

            this_ptd = unique_PTDs(p);
            sub2 = subset(abs(subset.PTD_ms-this_ptd) < PTD_Tol_ms,:);
            if isempty(sub2), continue; end
            sub2 = sortrows(sub2,'Distance_um');

            if abs(this_ptd) < PTD_Tol_ms
                cond_label = 'Simultaneous (PTD 0 ms)';
            else
                cond_label = sprintf('Sequential (PTD %g ms)',this_ptd);
            end

            valid = ~isnan(sub2.LinearityIndex);
            if any(valid)
                plot(sub2.Distance_um(valid),sub2.LinearityIndex(valid),'-o', ...
                    'Color',ptd_colors(p,:),'LineWidth',2,'MarkerFaceColor',ptd_colors(p,:), ...
                    'MarkerSize',7,'DisplayName',cond_label);
            end
        end

        yline(0,':k','LineWidth',1.5,'HandleVisibility','off');
        ylim([-1 1]);
        xlabel('Distance (\mum)','FontWeight','bold','FontSize',12);
        ylabel('Linearity index  (actual - predicted)/(actual + predicted)','FontWeight','bold','FontSize',12);
        title(sprintf('Amplitude %g \\muA  (0 = linear, <0 = sub-additive)',this_amp),'FontWeight','bold','FontSize',14);
        box off;
        lgd = legend('Location','best','Box','off');
        title(lgd,'Condition');
    end
end

%% ===================== NORMALISED RATIO PLOTS (actual / predicted) ==============
% Figure per amplitude: ratio vs distance, one line per PTD, dashed black line
% at 1 = the prediction. Then one figure: ratio vs predicted size.

if Plot_NormalizedRatio

    NormRatio = PlotData.Actual./PlotData.Predicted;
    NormRatio(~(PlotData.Predicted > Ratio_Min_Predicted)) = NaN;
    PlotData.NormRatio = NormRatio;

    fprintf('\n=========== NORMALISED RATIO (actual / predicted; NaN if predicted <= %g) ===========\n', ...
        Ratio_Min_Predicted);
    disp(sortrows(PlotData(:,{'Amp_uA','PTD_ms','Distance_um','Actual','Predicted','NormRatio'}), ...
        {'Amp_uA','PTD_ms','Distance_um'}));

    for a = 1:numel(unique_amps)

        this_amp = unique_amps(a);
        subset = PlotData(abs(PlotData.Amp_uA-this_amp) < Amp_Tol_uA,:);
        if isempty(subset), continue; end

        figure('Color','w','Position',[140 140 800 600]);
        hold on;

        x_all = [];
        y_all = [];
        for p = 1:numel(unique_PTDs)

            this_ptd = unique_PTDs(p);
            sub2 = subset(abs(subset.PTD_ms-this_ptd) < PTD_Tol_ms,:);
            if isempty(sub2), continue; end
            sub2 = sortrows(sub2,'Distance_um');

            if abs(this_ptd) < PTD_Tol_ms
                cond_label = 'Simultaneous (PTD 0 ms)';
            else
                cond_label = sprintf('Sequential (PTD %g ms)',this_ptd);
            end

            valid = ~isnan(sub2.NormRatio);
            if any(valid)
                plot(sub2.Distance_um(valid),sub2.NormRatio(valid),'-o', ...
                    'Color',ptd_colors(p,:),'LineWidth',2,'MarkerFaceColor',ptd_colors(p,:), ...
                    'MarkerSize',7,'DisplayName',cond_label);
                x_all = [x_all; sub2.Distance_um(valid)]; %#ok<AGROW>
                y_all = [y_all; sub2.NormRatio(valid)]; %#ok<AGROW>
            end
        end

        if ~isempty(x_all)
            plot([min(x_all) max(x_all)],[1 1],'--k','LineWidth',2, ...
                'DisplayName','Predicted (= 1)');
            ylim([0 max(1.5,1.1*max(y_all))]);
        end

        xlabel('Distance (\mum)','FontWeight','bold','FontSize',12);
        ylabel('Actual / Predicted','FontWeight','bold','FontSize',12);
        title(sprintf('Normalised to predicted - Amplitude %g \\muA',this_amp),'FontWeight','bold','FontSize',14);
        box off;
        lgd = legend('Location','best','Box','off');
        title(lgd,'Condition');
    end

    % ---- does electrode strength still leak into the ratio? ----
    figure('Color','w','Position',[160 160 800 600]);
    hold on;
    amp_markers = {'o','s','^','d','v','p','h'};

    fprintf('\nCorrelation between ratio and predicted size (all amplitudes and distances pooled, per PTD):\n');

    for p = 1:numel(unique_PTDs)

        this_ptd = unique_PTDs(p);
        mP = abs(PlotData.PTD_ms-this_ptd) < PTD_Tol_ms & ~isnan(PlotData.NormRatio);

        if abs(this_ptd) < PTD_Tol_ms
            cond_label = 'Simultaneous';
        else
            cond_label = sprintf('Sequential %g ms',this_ptd);
        end

        for a = 1:numel(unique_amps)
            m = mP & abs(PlotData.Amp_uA-unique_amps(a)) < Amp_Tol_uA;
            if ~any(m), continue; end
            mk = amp_markers{mod(a-1,numel(amp_markers))+1};
            scatter(PlotData.Predicted(m),PlotData.NormRatio(m),80,ptd_colors(p,:),mk,'filled', ...
                'DisplayName',sprintf('%s, %g \\muA',cond_label,unique_amps(a)));
            for q = find(m).'
                text(PlotData.Predicted(q)+0.05,PlotData.NormRatio(q), ...
                    sprintf('%g',PlotData.Distance_um(q)),'FontSize',8);
            end
        end

        if sum(mP) >= 3
            cc = corrcoef(PlotData.Predicted(mP),PlotData.NormRatio(mP));
            fprintf('  %-16s: r = %.2f  (n = %d)\n',cond_label,cc(1,2),sum(mP));
        end
    end

    plot(xlim,[1 1],'--k','LineWidth',1.5,'HandleVisibility','off');
    xlabel('Predicted (mean per channel)','FontWeight','bold','FontSize',12);
    ylabel('Actual / Predicted','FontWeight','bold','FontSize',12);
    title('Does the ratio depend on predicted size?  (labels = distance, \mum)','FontWeight','bold','FontSize',13);
    box off;
    legend('Location','best','Box','off','FontSize',8);
    fprintf('(a clearly negative r means electrode strength / saturation still leaks into the ratio)\n');
end

%% ============ CHANNEL-MAX NORMALISATION (each channel / its own max single response) ============
% WHY: Actual and Predicted are divided, channel by channel, by that
% channel's largest SINGLE-electrode response (taken from singles only, so the
% paired data never enter their own normaliser). Every channel then counts
% equally and is expressed as a fraction of its own ceiling (1 = max single).
%   Figure A (per amplitude): normalised ratio vs distance + normalised predicted
%            vs distance (shows how much predicted size still changes with distance).
%   Figure B: ACTUAL vs PREDICTED in normalised units, binned by predicted
%            size, near vs far pairs, one panel per PTD. Points on y = x are
%            linear; bending below at high predicted = saturation. Near vs far
%            curves overlapping = no distance effect at matched predicted size.
% CAUTION: predicted is noisy, and binning on a noisy x makes the highest bins
% look more sub-additive (regression dilution). This affects near and far
% alike, so the near-vs-far comparison is still fair, but the absolute bend is
% slightly exaggerated. Channels are reused across pairs, so SEM bars are rough.

if Plot_ChannelMaxNorm

    R = Results(~isnan(Results.Distance_um),:);
    nR = height(R);

    if nR == 0 || all(cellfun(@(x) all(isnan(x)),R.SingleA_PerChannel))
        warning('Channel-max normalisation: no usable single responses - skipped.');
    else
        chan_len = unique(cellfun(@numel,R.Actual_PerChannel));
        if numel(chan_len) ~= 1
            error('Channel-max normalisation needs the same channel list in every row.');
        end

        rowmat = @(c) cell2mat(cellfun(@(x) x(:).',c,'UniformOutput',false));
        ActM = rowmat(R.Actual_PerChannel);
        PrdM = rowmat(R.Predicted_PerChannel);
        SAM  = rowmat(R.SingleA_PerChannel);
        SBM  = rowmat(R.SingleB_PerChannel);
        nCh  = size(ActM,2);

        % ---- channel maximum from SINGLE responses only ----
        ChMax = nan(nR,nCh);
        switch lower(Norm_Max_Scope)
            case 'all'
                m = max([SAM;SBM],[],1,'omitnan');
                ChMax = repmat(m,nR,1);
            case 'amp'
                for ua = unique(R.Amp_uA).'
                    rr = abs(R.Amp_uA-ua) < Amp_Tol_uA;
                    m = max([SAM(rr,:);SBM(rr,:)],[],1,'omitnan');
                    ChMax(rr,:) = repmat(m,nnz(rr),1);
                end
            otherwise
                error('Norm_Max_Scope must be ''all'' or ''amp''.');
        end

        keep = ChMax >= Norm_Min_Max;     % false for NaN as well
        ActN = ActM./ChMax;  ActN(~keep) = NaN;
        PrdN = PrdM./ChMax;  PrdN(~keep) = NaN;

        fprintf('\n=========== CHANNEL-MAX NORMALISATION (scope = %s) ===========\n',Norm_Max_Scope);
        fprintf('Channels kept (max single response >= %g spikes/trial): %d of %d (per row, min..max)\n', ...
            Norm_Min_Max,min(sum(keep,2)),nCh);
        fprintf('Predicted in these units is up to 2 by construction (sum of two singles, each <= 1).\n');

        % ---- per (distance, amp, PTD) summary, channels weighted equally ----
        [gk,~,gi] = unique([R.Distance_um,R.Amp_uA,R.PTD_ms],'rows');
        nG = size(gk,1);
        G_ActN = nan(nG,1); G_PrdN = G_ActN; G_Ratio = G_ActN; G_nCh = G_ActN;
        for g = 1:nG
            rr = gi == g;
            a = ActN(rr,:); p = PrdN(rr,:);
            ok = ~isnan(a) & ~isnan(p);
            G_ActN(g) = mean(a(ok));
            G_PrdN(g) = mean(p(ok));
            G_nCh(g)  = nnz(ok);
            if G_PrdN(g) > Norm_Min_PredMean
                G_Ratio(g) = G_ActN(g)/G_PrdN(g);
            end
        end
        NormTable = table(gk(:,1),gk(:,2),gk(:,3),G_ActN,G_PrdN,G_Ratio,G_nCh, ...
            'VariableNames',{'Distance_um','Amp_uA','PTD_ms','Actual_N','Predicted_N','Ratio_N','N_ChPoints'});
        NormTable = sortrows(NormTable,{'Amp_uA','PTD_ms','Distance_um'});
        fprintf('\nNormalised actual / predicted per condition (ratio = Actual_N/Predicted_N):\n');
        disp(NormTable);

        amps_n = unique(NormTable.Amp_uA);
        ptds_n = unique(NormTable.PTD_ms);
        cols_n = lines(numel(ptds_n));

        % ---------------- Figure A: per amplitude ----------------
        for a = 1:numel(amps_n)
            figure('Color','w','Position',[80+30*a 100 1250 480]);
            sub = NormTable(abs(NormTable.Amp_uA-amps_n(a)) < Amp_Tol_uA,:);

            subplot(1,2,1); hold on;
            for p = 1:numel(ptds_n)
                s_ = sub(abs(sub.PTD_ms-ptds_n(p)) < PTD_Tol_ms,:);
                if isempty(s_), continue; end
                if ptds_n(p) == 0
                    nm = 'Simultaneous (PTD 0 ms)';
                else
                    nm = sprintf('Sequential (PTD %g ms)',ptds_n(p));
                end
                plot(s_.Distance_um,s_.Ratio_N,'-o','Color',cols_n(p,:),'LineWidth',2.5, ...
                    'MarkerFaceColor',cols_n(p,:),'MarkerSize',8,'DisplayName',nm);
            end
            plot(xlim,[1 1],'--k','LineWidth',1.5,'DisplayName','Linear (= 1)');
            xlabel('Distance (\mum)','FontWeight','bold','FontSize',12);
            ylabel('Normalised actual / normalised predicted','FontWeight','bold','FontSize',12);
            title(sprintf('Channel-max normalised ratio - %g \\muA',amps_n(a)),'FontWeight','bold','FontSize',13);
            legend('Location','best','Box','off'); box off;

            subplot(1,2,2); hold on;
            s_ = sub(abs(sub.PTD_ms-ptds_n(1)) < PTD_Tol_ms,:);
            plot(s_.Distance_um,s_.Predicted_N,'-ks','LineWidth',2.5,'MarkerFaceColor','k','MarkerSize',8, ...
                'DisplayName','Predicted (normalised)');
            for p = 1:numel(ptds_n)
                s2 = sub(abs(sub.PTD_ms-ptds_n(p)) < PTD_Tol_ms,:);
                if isempty(s2), continue; end
                plot(s2.Distance_um,s2.Actual_N,'-o','Color',cols_n(p,:),'LineWidth',2, ...
                    'MarkerFaceColor',cols_n(p,:),'MarkerSize',7, ...
                    'DisplayName',sprintf('Actual, PTD %g ms',ptds_n(p)));
            end
            xlabel('Distance (\mum)','FontWeight','bold','FontSize',12);
            ylabel('Mean fraction of channel max','FontWeight','bold','FontSize',12);
            title('Normalised actual and predicted','FontWeight','bold','FontSize',13);
            legend('Location','best','Box','off'); box off;
        end

        % ---------------- Figure B: actual vs predicted, near vs far ----------------
        distM = repmat(R.Distance_um,1,nCh);
        ptdM  = repmat(R.PTD_ms,1,nCh);
        okAll = ~isnan(ActN) & ~isnan(PrdN);

        edges = IO_Bin_Edges(:).';
        nb = numel(edges)-1;
        grp_names = {sprintf('Near (<= %g \\mum)',Near_Max_um),sprintf('Far (>= %g \\mum)',Far_Min_um)};
        grp_cols  = [0.20 0.45 0.75; 0.85 0.40 0.10];

        figure('Color','w','Position',[60 80 450*numel(ptds_n) 460]);
        fprintf('\nActual vs predicted (normalised), binned by predicted, near vs far:\n');
        for p = 1:numel(ptds_n)
            subplot(1,numel(ptds_n),p); hold on;
            inPTD = okAll & abs(ptdM-ptds_n(p)) < PTD_Tol_ms;
            for gname = 1:2
                if gname == 1
                    inG = inPTD & distM <= Near_Max_um;
                else
                    inG = inPTD & distM >= Far_Min_um;
                end
                px = PrdN(inG); ay = ActN(inG);
                bx = nan(nb,1); by = bx; be = bx; bn = zeros(nb,1);
                for b = 1:nb
                    mb = px >= edges(b) & px < edges(b+1);
                    bn(b) = nnz(mb);
                    if bn(b) >= IO_Min_Points
                        bx(b) = mean(px(mb));
                        by(b) = mean(ay(mb));
                        be(b) = std(ay(mb))/sqrt(bn(b));
                    end
                end
                errorbar(bx,by,be,'-o','Color',grp_cols(gname,:),'LineWidth',2, ...
                    'MarkerFaceColor',grp_cols(gname,:),'MarkerSize',7,'DisplayName',grp_names{gname});
                if gname == 1, bn_near = bn; by_near = by; bx_near = bx; end
            end
            % print comparison for this PTD
            fprintf('  PTD %g ms:  bin of predicted | N near | N far | mean pred near/far | actual near/far | ratio near/far\n',ptds_n(p));
            % recompute far values for the printout
            inF = inPTD & distM >= Far_Min_um;
            pxF = PrdN(inF); ayF = ActN(inF);
            for b = 1:nb
                mbF = pxF >= edges(b) & pxF < edges(b+1);
                nF = nnz(mbF);
                if bn_near(b) >= IO_Min_Points && nF >= IO_Min_Points
                    fprintf('    [%4.2f,%4.2f) | %5d | %5d | %.2f / %.2f | %.2f / %.2f | %.2f / %.2f\n', ...
                        edges(b),min(edges(b+1),99),bn_near(b),nF,bx_near(b),mean(pxF(mbF)), ...
                        by_near(b),mean(ayF(mbF)),by_near(b)/bx_near(b),mean(ayF(mbF))/mean(pxF(mbF)));
                end
            end
            lim = max([xlim ylim]);
            plot([0 lim],[0 lim],'--k','LineWidth',1.5,'DisplayName','Linear (y = x)');
            xlabel('Predicted (fraction of channel max)','FontWeight','bold','FontSize',11);
            ylabel('Actual (fraction of channel max)','FontWeight','bold','FontSize',11);
            if ptds_n(p) == 0
                title('Simultaneous (PTD 0 ms)','FontWeight','bold','FontSize',12);
            else
                title(sprintf('Sequential (PTD %g ms)',ptds_n(p)),'FontWeight','bold','FontSize',12);
            end
            legend('Location','northwest','Box','off','FontSize',9); box off;
        end
        fprintf('(near and far overlapping = no distance effect at matched predicted size; curves below y = x at high predicted = saturation)\n');
    end
end

%% ======================= PER-CHANNEL DIAGNOSTIC PLOTS ===================
% One figure per physical channel (map position), amp and PTD both shown
% within that one figure (color = amp, marker = PTD, solid/dashed =
% actual/predicted). Channel numbers are treated as the same physical
% channel across every dataset-pair (same session, same probe, no
% re-insertion between them). A channel only has a point at a distance
% where it was in that condition's saved responding-channel list - gaps
% are expected, not an error, and are printed below for each channel.
% Skipped entirely if Plot_PerChannel_Diagnostics is false.

if Plot_PerChannel_Diagnostics

fprintf('\n============================================================\n');
fprintf('PER-CHANNEL DIAGNOSTIC PLOTS\n');
fprintf('============================================================\n');

%% ---- build long-format table: one row per (condition, channel) ----

Chn_Channel  = [];
Chn_Distance = [];
Chn_Amp      = [];
Chn_PTD      = [];
Chn_Actual   = [];
Chn_Predicted = [];

for row_index = 1:height(Results)

    if isnan(Results.Distance_um(row_index))
        continue;   % same exclusion as the population plot
    end

    channel_list = Results.Channel_List{row_index};
    actual_vals  = Results.Actual_PerChannel{row_index};
    pred_vals    = Results.Predicted_PerChannel{row_index};

    for c = 1:numel(channel_list)
        Chn_Channel(end+1,1)   = channel_list(c); %#ok<AGROW>
        Chn_Distance(end+1,1)  = Results.Distance_um(row_index); %#ok<AGROW>
        Chn_Amp(end+1,1)       = Results.Amp_uA(row_index); %#ok<AGROW>
        Chn_PTD(end+1,1)       = Results.PTD_ms(row_index); %#ok<AGROW>
        Chn_Actual(end+1,1)    = actual_vals(c); %#ok<AGROW>
        Chn_Predicted(end+1,1) = pred_vals(c); %#ok<AGROW>
    end
end

PerChannelLong = table(Chn_Channel,Chn_Distance,Chn_Amp,Chn_PTD,Chn_Actual,Chn_Predicted, ...
    'VariableNames',{'Channel','Distance_um','Amp_uA','PTD_ms','Actual','Predicted'});

%% ---- average duplicates (same channel/distance/amp/PTD, e.g. both stim orders) ----

[chn_group_keys,~,chn_group_idx] = unique( ...
    [PerChannelLong.Channel,PerChannelLong.Distance_um,PerChannelLong.Amp_uA,PerChannelLong.PTD_ms],'rows');

nChnGroups = size(chn_group_keys,1);

ChnG_Channel  = chn_group_keys(:,1);
ChnG_Distance = chn_group_keys(:,2);
ChnG_Amp      = chn_group_keys(:,3);
ChnG_PTD      = chn_group_keys(:,4);
ChnG_Actual   = nan(nChnGroups,1);
ChnG_Predicted = nan(nChnGroups,1);

for g = 1:nChnGroups
    rows_in_group = chn_group_idx == g;
    ChnG_Actual(g)    = mean(PerChannelLong.Actual(rows_in_group),'omitnan');
    ChnG_Predicted(g) = mean(PerChannelLong.Predicted(rows_in_group),'omitnan');
end

PerChannelAvg = table(ChnG_Channel,ChnG_Distance,ChnG_Amp,ChnG_PTD,ChnG_Actual,ChnG_Predicted, ...
    'VariableNames',{'Channel','Distance_um','Amp_uA','PTD_ms','Actual','Predicted'});

%% ---- consistent color (amp) / marker (PTD) mapping across all figures ----

all_amps = unique(PerChannelAvg.Amp_uA);
all_ptds = unique(PerChannelAvg.PTD_ms);

amp_colors_chn = lines(numel(all_amps));
ptd_markers = {'o','s','^','d','v','p','h'};
if numel(all_ptds) > numel(ptd_markers)
    ptd_markers = repmat(ptd_markers,1,ceil(numel(all_ptds)/numel(ptd_markers)));
end

unique_channels = unique(PerChannelAvg.Channel);

fprintf('Channels with at least one plotted point: %d\n',numel(unique_channels));
fprintf('Distances present overall: %s\n',num2str(unique(PerChannelAvg.Distance_um).'));

if isempty(Channels_To_Plot)
    channels_to_plot = unique_channels;
else
    requested_channels = unique(double(Channels_To_Plot(:).'));
    channels_to_plot = intersect(unique_channels,requested_channels,'stable');
    missing_channels = setdiff(requested_channels,unique_channels);
    if ~isempty(missing_channels)
        warning('Channels_To_Plot includes channel(s) with no data at all - skipped: %s', ...
            num2str(missing_channels));
    end
end

fprintf('Channels selected for plotting: %d (%s)\n', ...
    numel(channels_to_plot),num2str(channels_to_plot(:).'));

for ch_i = 1:numel(channels_to_plot)

    this_channel = channels_to_plot(ch_i);
    chn_rows = PerChannelAvg(PerChannelAvg.Channel == this_channel,:);

    distances_covered = unique(chn_rows.Distance_um);
    fprintf('  Channel %d: distances covered = %s\n', ...
        this_channel,num2str(distances_covered(:).'));

    figure('Color','w','Position',[100 100 800 600]);
    hold on;

    for ai = 1:numel(all_amps)
        for pj = 1:numel(all_ptds)

            this_amp = all_amps(ai);
            this_ptd = all_ptds(pj);

            sub = chn_rows(abs(chn_rows.Amp_uA-this_amp) < Amp_Tol_uA & ...
                abs(chn_rows.PTD_ms-this_ptd) < PTD_Tol_ms,:);

            if isempty(sub)
                continue;
            end

            sub = sortrows(sub,'Distance_um');
            col = amp_colors_chn(ai,:);
            mk = ptd_markers{pj};

            valid_actual = ~isnan(sub.Actual);
            if any(valid_actual)
                plot(sub.Distance_um(valid_actual),sub.Actual(valid_actual), ...
                    '-','Marker',mk,'Color',col,'LineWidth',1.5,'MarkerFaceColor',col,'MarkerSize',6, ...
                    'DisplayName',sprintf('%g uA, PTD %g ms (actual)',this_amp,this_ptd));
            end

            valid_pred = ~isnan(sub.Predicted);
            if any(valid_pred)
                plot(sub.Distance_um(valid_pred),sub.Predicted(valid_pred), ...
                    '--','Marker',mk,'Color',col,'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',6, ...
                    'DisplayName',sprintf('%g uA, PTD %g ms (predicted)',this_amp,this_ptd));
            end
        end
    end

    xlabel('Distance (\mum)','FontWeight','bold','FontSize',12);
    ylabel('Baseline-corrected spike count (per trial)','FontWeight','bold','FontSize',12);
    title(sprintf('Channel %d',this_channel),'FontWeight','bold','FontSize',14);
    box off;
    lgd = legend('Location','best','Box','off','FontSize',8);
    title(lgd,'Condition');
end

end   % end of "if Plot_PerChannel_Diagnostics"

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

function DS = load_linearity_dataset(spike_file,dataset_folder,Electrode_Type,FS_fallback,apply_badtrial_exclusion)
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

if ~apply_badtrial_exclusion
    fprintf('Apply_BadTrial_Exclusion is false - bad-trial file not loaded, no trials excluded:\n%s\n',dataset_folder);
    bad_trial_ids = [];
elseif isfile(badtrials_path)
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