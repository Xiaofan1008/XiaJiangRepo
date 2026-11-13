%% ========================================================================
% CHECK_SINGLERESPONSES_V1 - single-electrode responses, bad trials, and an
% independent re-computation of PREDICTED
%
% PURPOSE
%   Predicted = single(anchor) + single(B). This script looks only at the
%   SINGLE-pulse datasets (nothing from the paired data) so we can see:
%     1. single(anchor Ch1) and single(each ChB) SEPARATELY vs distance
%        from the anchor  -> which term carries the rise in Predicted
%     2. what the single bad-trial exclusion does (and whether a
%        <base_name>_BadTrials.mat exists for each single dataset at all)
%     3. trial-by-trial values vs absolute trial ID (drift, outliers, and
%        which trials are marked bad - red x)
%     4. early-window [0,3) ms vs late-window [3,20) ms counts (a residual
%        stimulation artifact would show up as a large early count)
%     5. an independent re-computation of Predicted from the same data, to
%        cross-check the linearity script (should match it when the same
%        channel list and the same bad-trial setting are used)
%
% Plots + printed tables only. Nothing is saved.
% ========================================================================

clear;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

Singles = struct('spike_file',{},'folder',{});

Singles(1).spike_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1/Xia_Linearity_S1_100_400um_Single1.sp_xia_filtered.mat';
Singles(1).folder     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1';

Singles(2).spike_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1/Xia_Linearity_S1_500_700um_Single1.sp_xia_filtered.mat';
Singles(2).folder     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1';

Channel_List   = 1:64;    % fixed recording channels used for every electrode
Anchor_Channel = 1;       % the electrode that is in every pair

Electrode_Type = 3;
FS_fallback    = 30000;

baseline_win_ms = [-50 -5];
analysis_win_ms = [0 20];
early_win_ms    = [0 3];     % early part of the window (artifact check)
late_win_ms     = [3 20];

Amp_Tol_uA = 0.001;

%% =========================== MAIN COMPUTATION ==========================

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

Row_ds = []; Row_dsLabel = {}; Row_ch = []; Row_amp = []; Row_dist = [];
Row_nTr = []; Row_nBad = []; Row_all = []; Row_good = [];
Row_early = []; Row_late = []; Row_ids = {}; Row_vals = {}; Row_isbad = {};

for d = 1:numel(Singles)

    fprintf('\n============================================================\n');
    fprintf('SINGLE DATASET %d/%d\n%s\n',d,numel(Singles),Singles(d).folder);
    fprintf('Spike file: %s\n',Singles(d).spike_file);
    fprintf('============================================================\n');

    DS = load_linearity_dataset(Singles(d).spike_file,Singles(d).folder, ...
        Electrode_Type,FS_fallback,true);

    if isempty(DS.bad_trial_ids)
        fprintf('\n*** NO BAD TRIALS in effect for this single dataset (file missing or empty) -\n');
        fprintf('*** all trials are used. Run the bad-trial QC on this dataset.\n');
    else
        fprintf('\nBad trials loaded for this single dataset: %d\n',numel(DS.bad_trial_ids));
    end

    for ss = 1:DS.nSets

        stim_channels = DS.uniqueComb(ss,DS.uniqueComb(ss,:) > 0);
        if numel(stim_channels) ~= 1, continue; end
        stim_ch = stim_channels(1);

        dist_um = ChnPairDistance(Anchor_Channel,stim_ch,Electrode_Type);

        for ai = 1:numel(DS.Amps)

            if DS.Amps(ai) <= 0, continue; end

            trials_all = find(DS.combClass == ss & DS.ampIdx == ai);
            if isempty(trials_all), continue; end

            is_bad  = ismember(trials_all,DS.bad_trial_ids);
            trig    = DS.trig_ms(trials_all);

            v_full  = trial_values(DS.sp_by_mappos,trig,Channel_List,baseline_win_ms,analysis_win_ms);
            v_early = trial_values(DS.sp_by_mappos,trig,Channel_List,baseline_win_ms,early_win_ms);
            v_late  = trial_values(DS.sp_by_mappos,trig,Channel_List,baseline_win_ms,late_win_ms);

            Row_ds(end+1,1)      = d; %#ok<AGROW>
            Row_dsLabel{end+1,1} = DS.base_name; %#ok<AGROW>
            Row_ch(end+1,1)      = stim_ch; %#ok<AGROW>
            Row_amp(end+1,1)     = DS.Amps(ai); %#ok<AGROW>
            Row_dist(end+1,1)    = dist_um; %#ok<AGROW>
            Row_nTr(end+1,1)     = numel(trials_all); %#ok<AGROW>
            Row_nBad(end+1,1)    = sum(is_bad); %#ok<AGROW>
            Row_all(end+1,1)     = mean(v_full,'omitnan'); %#ok<AGROW>
            Row_good(end+1,1)    = mean(v_full(~is_bad),'omitnan'); %#ok<AGROW>
            Row_early(end+1,1)   = mean(v_early(~is_bad),'omitnan'); %#ok<AGROW>
            Row_late(end+1,1)    = mean(v_late(~is_bad),'omitnan'); %#ok<AGROW>
            Row_ids{end+1,1}     = trials_all(:); %#ok<AGROW>
            Row_vals{end+1,1}    = v_full(:); %#ok<AGROW>
            Row_isbad{end+1,1}   = is_bad(:); %#ok<AGROW>
        end
    end
end

%% ============================ PRINTED TABLES ===========================

T = table(Row_ds,Row_ch,Row_dist,Row_amp,Row_nTr,Row_nBad,Row_all,Row_good,Row_early,Row_late, ...
    'VariableNames',{'Dataset','StimCh','Dist_from_anchor_um','Amp_uA','N_Trials','N_Bad', ...
    'Single_AllTrials','Single_GoodTrials','Early_0to3ms','Late_3to20ms'});
T = sortrows(T,{'Amp_uA','Dataset','StimCh'});

fprintf('\n=========== SINGLE-ELECTRODE RESPONSES (mean over channels %d-%d) ===========\n', ...
    min(Channel_List),max(Channel_List));
fprintf('Early/Late columns use GOOD trials only.\n');
disp(T);

% ---- independent re-computation of Predicted = single(anchor)+single(B) ----
P_ds = []; P_ch = []; P_dist = []; P_amp = []; P_anchor = []; P_B = []; P_pred = []; P_predAll = [];
for r = 1:numel(Row_ch)
    if Row_ch(r) == Anchor_Channel, continue; end
    ra = find(Row_ds == Row_ds(r) & Row_ch == Anchor_Channel & ...
        abs(Row_amp-Row_amp(r)) < Amp_Tol_uA,1);
    if isempty(ra), continue; end
    P_ds(end+1,1) = Row_ds(r); P_ch(end+1,1) = Row_ch(r); P_dist(end+1,1) = Row_dist(r); %#ok<AGROW>
    P_amp(end+1,1) = Row_amp(r); %#ok<AGROW>
    P_anchor(end+1,1) = Row_good(ra); P_B(end+1,1) = Row_good(r); %#ok<AGROW>
    P_pred(end+1,1) = Row_good(ra)+Row_good(r); %#ok<AGROW>
    P_predAll(end+1,1) = Row_all(ra)+Row_all(r); %#ok<AGROW>
end

TP = table(P_ds,P_ch,P_dist,P_amp,P_anchor,P_B,P_pred,P_predAll, ...
    'VariableNames',{'Dataset','ChB','Dist_um','Amp_uA','Single_Anchor','Single_B', ...
    'Predicted_GoodTrials','Predicted_AllTrials'});
TP = sortrows(TP,{'Amp_uA','Dist_um'});

fprintf('\n=========== RE-COMPUTED PREDICTED (compare with the linearity script) ===========\n');
disp(TP);

%% ================================ FIGURES ==============================
% One figure per amplitude. Left: single(B) vs distance from the anchor
% (star = the anchor itself, one per dataset). Right: early vs late window.

amps = unique(Row_amp);

for a = 1:numel(amps)

    m  = abs(Row_amp-amps(a)) < Amp_Tol_uA;
    mb = m & Row_ch ~= Anchor_Channel;
    ma = m & Row_ch == Anchor_Channel;

    [xs,ord] = sort(Row_dist(mb));
    chs = Row_ch(mb);    chs = chs(ord);
    ya  = Row_all(mb);   ya  = ya(ord);
    yg  = Row_good(mb);  yg  = yg(ord);
    ye  = Row_early(mb); ye  = ye(ord);
    yl  = Row_late(mb);  yl  = yl(ord);

    figure('Color','w','Position',[80 100 1250 500]);

    subplot(1,2,1); hold on;
    plot(xs,ya,'--s','Color',[0.55 0.55 0.55],'LineWidth',1.5,'MarkerSize',7, ...
        'DisplayName','single(ChB): all trials');
    plot(xs,yg,'-o','Color',[0.1 0.4 0.8],'LineWidth',2,'MarkerFaceColor',[0.1 0.4 0.8],'MarkerSize',7, ...
        'DisplayName','single(ChB): bad trials removed');
    for k = 1:numel(xs)
        text(xs(k)+8,yg(k),sprintf('Ch%d',chs(k)),'FontSize',8);
    end
    for r = find(ma).'
        plot(0,Row_good(r),'p','MarkerSize',14,'MarkerFaceColor','k','Color','k', ...
            'DisplayName',sprintf('single(Ch%d) anchor: %s',Anchor_Channel,Row_dsLabel{r}));
    end
    xlabel('Distance of ChB from the anchor (\mum)','FontWeight','bold');
    ylabel('Mean baseline-corrected count (per trial, per channel)','FontWeight','bold');
    title(sprintf('Single-electrode responses - %g \\muA',amps(a)),'FontWeight','bold');
    legend('Location','best','Box','off','Interpreter','none','FontSize',8);
    box off;

    subplot(1,2,2); hold on;
    plot(xs,ye,'-o','LineWidth',2,'MarkerSize',7,'DisplayName',sprintf('early [%g,%g) ms',early_win_ms));
    plot(xs,yl,'-s','LineWidth',2,'MarkerSize',7,'DisplayName',sprintf('late [%g,%g) ms',late_win_ms));
    xlabel('Distance of ChB from the anchor (\mum)','FontWeight','bold');
    ylabel('Mean baseline-corrected count','FontWeight','bold');
    title('Early vs late window (good trials)','FontWeight','bold');
    legend('Location','best','Box','off');
    box off;
end

% Trial-by-trial values vs absolute trial ID (drift / outliers / bad trials)
for d = 1:numel(Singles)

    rows = find(Row_ds == d);
    if isempty(rows), continue; end
    sets_here = unique(Row_ch(rows),'stable');

    figure('Color','w','Position',[60 60 1500 360]);
    for s = 1:numel(sets_here)
        subplot(1,numel(sets_here),s); hold on;
        for r = rows(Row_ch(rows) == sets_here(s)).'
            plot(Row_ids{r},Row_vals{r},'-o','MarkerSize',4,'DisplayName',sprintf('%g uA',Row_amp(r)));
            bad = Row_isbad{r};
            if any(bad)
                plot(Row_ids{r}(bad),Row_vals{r}(bad),'rx','MarkerSize',9,'LineWidth',1.5, ...
                    'HandleVisibility','off');
            end
        end
        title(sprintf('Ch%d',sets_here(s)),'FontWeight','bold');
        xlabel('Absolute trial ID');
        if s == 1
            ylabel('Mean over channels (corrected count)');
            legend('Location','best','Box','off');
        end
        box off;
    end
    sgtitle(sprintf('Single dataset %d: per-trial response (red x = bad trial)',d),'Interpreter','none');
end

fprintf('\n============================================================\n');
fprintf('DONE (plots + tables only - nothing saved)\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function vals = trial_values(sp_by_mappos,trig_ms_list,channels,baseline_win_ms,win_ms)
% Per-trial baseline-corrected spike count in win_ms, averaged over the
% given channels (same correction as mean_corrected_count in the linearity
% script). Returns one value per trigger.

nT = numel(trig_ms_list);
base_dur = diff(baseline_win_ms);
win_dur  = diff(win_ms);
per = nan(nT,numel(channels));

for c = 1:numel(channels)
    ch = channels(c);
    if ch < 1 || ch > numel(sp_by_mappos), continue; end
    st = sp_by_mappos{ch};
    if isempty(st)
        per(:,c) = 0;
        continue;
    end
    for k = 1:nT
        t0 = trig_ms_list(k);
        b = sum(st >= t0+baseline_win_ms(1) & st < t0+baseline_win_ms(2));
        w = sum(st >= t0+win_ms(1) & st < t0+win_ms(2));
        per(k,c) = w-b*(win_dur/base_dur);
    end
end

vals = mean(per,2,'omitnan');
end


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