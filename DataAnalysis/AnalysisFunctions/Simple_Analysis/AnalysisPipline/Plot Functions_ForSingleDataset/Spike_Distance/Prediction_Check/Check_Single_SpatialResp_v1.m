%% ========================================================================
% CHECK_SINGLESPATIALPROFILE_V1 - where on the probe does each single
% stimulation electrode evoke spikes?
%
% PURPOSE
%   Single-electrode datasets only. For every stimulation electrode and
%   amplitude, measure the baseline-corrected response of EVERY recording
%   channel, then ask:
%     1. HOT SPOT CHECK - is the strongest response near the stimulated
%        electrode? (tests that electrode labels / mapping are consistent)
%     2. LOCAL vs FAR - the whole-array mean mixes a local response near the
%        electrode with distant channels. Which one carries the rise of
%        single(ChB) with channel number?
%     3. ROW vs COLUMN - is the strength explained by the RECORDING channel
%        (some channels respond to everything) or by the STIMULATED
%        electrode (response moves with the electrode)?
%
% Plots + printed tables only. Nothing is saved.
% ========================================================================

clear;clc;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

Singles = struct('spike_file',{},'folder',{});

Singles(1).spike_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1/Xia_Linearity_S1_100_400um_Single1.sp_xia_filtered.mat';
Singles(1).folder     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1';

Singles(2).spike_file = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1/Xia_Linearity_S1_500_700um_Single1.sp_xia_filtered.mat';
Singles(2).folder     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1';

Recording_Channels = 1:64;   % four-shank probe only (distances undefined for 65-96)
Anchor_Channel     = 1;      % electrode present in every pair

Apply_BadTrial_Exclusion = true;   % each dataset's own saved bad trials

Local_Radius_um     = 150;   % "local" = recording channels within this distance of the stim electrode
Ring_Edges_um       = [0 100 200 300 400 500 600 800];   % distance rings (first ring excludes d = 0)
Exclude_StimChannel = true;  % leave the stimulated electrode's own channel out of every summary

Map_Amp_uA = 5;              % amplitude used for the probe maps / profiles / matrix
Amp_Tol_uA = 0.001;

Electrode_Type = 3;
FS_fallback    = 30000;

baseline_win_ms = [-50 -5];
analysis_win_ms = [0 20];

%% =========================== MAIN COMPUTATION ==========================

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

nRec = numel(Recording_Channels);
xy = nan(nRec,2);
for i = 1:nRec
    [xy(i,1),xy(i,2)] = ChnMap2Coord(Recording_Channels(i),Electrode_Type);
end

Row_ds = []; Row_stim = []; Row_amp = []; Row_nTr = [];
Prof = zeros(0,nRec);

for d = 1:numel(Singles)

    fprintf('\n============================================================\n');
    fprintf('SINGLE DATASET %d/%d\n%s\n',d,numel(Singles),Singles(d).folder);
    fprintf('============================================================\n');

    DS = load_linearity_dataset(Singles(d).spike_file,Singles(d).folder, ...
        Electrode_Type,FS_fallback,Apply_BadTrial_Exclusion);

    for ss = 1:DS.nSets

        stim_channels = DS.uniqueComb(ss,DS.uniqueComb(ss,:) > 0);
        if numel(stim_channels) ~= 1, continue; end

        for ai = 1:numel(DS.Amps)

            if DS.Amps(ai) <= 0, continue; end

            trials_all = find(DS.combClass == ss & DS.ampIdx == ai);
            trials_use = setdiff(trials_all,DS.bad_trial_ids);
            if isempty(trials_use), continue; end

            trig = DS.trig_ms(trials_use);

            prof = nan(1,nRec);
            for c = 1:nRec
                ch = Recording_Channels(c);
                if ch < 1 || ch > DS.nMapPos, continue; end
                prof(c) = mean_corrected_count(DS.sp_by_mappos{ch},trig, ...
                    baseline_win_ms,analysis_win_ms);
            end

            Row_ds(end+1,1)   = d; %#ok<AGROW>
            Row_stim(end+1,1) = stim_channels(1); %#ok<AGROW>
            Row_amp(end+1,1)  = DS.Amps(ai); %#ok<AGROW>
            Row_nTr(end+1,1)  = numel(trials_use); %#ok<AGROW>
            Prof(end+1,:)     = prof; %#ok<AGROW>
        end
    end
end

nRows = numel(Row_stim);
nRing = numel(Ring_Edges_um)-1;

All_ = nan(nRows,1); Local_ = nan(nRows,1); Far_ = nan(nRows,1);
Ring_ = nan(nRows,nRing);
DistMat = nan(nRows,nRec);
Top3Text = cell(nRows,1);
Peak1Ch = nan(nRows,1); Peak1Dist = nan(nRows,1);

for r = 1:nRows

    [sx,sy] = ChnMap2Coord(Row_stim(r),Electrode_Type);
    dvec = sqrt((xy(:,1)-sx).^2+(xy(:,2)-sy).^2);
    DistMat(r,:) = dvec(:).';

    prof = Prof(r,:).';
    valid = ~isnan(dvec) & ~isnan(prof);
    if Exclude_StimChannel
        valid = valid & (Recording_Channels(:) ~= Row_stim(r));
    end

    All_(r)   = mean(prof(valid));
    Local_(r) = mean(prof(valid & dvec <= Local_Radius_um));
    Far_(r)   = mean(prof(valid & dvec >  Local_Radius_um));

    for k = 1:nRing
        m = valid & dvec > Ring_Edges_um(k) & dvec <= Ring_Edges_um(k+1);
        Ring_(r,k) = mean(prof(m));
    end

    tmp = prof; tmp(~valid) = -Inf;
    [~,ord] = sort(tmp,'descend');
    top = ord(1:min(3,numel(ord)));
    Peak1Ch(r)   = Recording_Channels(top(1));
    Peak1Dist(r) = dvec(top(1));
    parts = cell(1,numel(top));
    for t = 1:numel(top)
        parts{t} = sprintf('Ch%d(%.0fum)',Recording_Channels(top(t)),dvec(top(t)));
    end
    Top3Text{r} = strjoin(parts,' ');
end

%% ============================ PRINTED TABLES ===========================

T = table(Row_ds,Row_stim,Row_amp,Row_nTr,All_,Local_,Far_,Peak1Ch,Peak1Dist,Top3Text, ...
    'VariableNames',{'Dataset','StimCh','Amp_uA','N_Trials','Mean_AllChn', ...
    sprintf('Local_le%gum',Local_Radius_um),'Far','Peak_Ch','Peak_Dist_um','Top3_Channels_Dist'});
T = sortrows(T,{'Amp_uA','Dataset','StimCh'});

fprintf('\n=========== SINGLE RESPONSE: ALL vs LOCAL vs FAR, AND HOT SPOT ===========\n');
fprintf('Recording channels %d-%d%s. Local = within %g um of the stim electrode.\n', ...
    min(Recording_Channels),max(Recording_Channels), ...
    ternary(Exclude_StimChannel,' (stim channel itself excluded)',''),Local_Radius_um);
disp(T);

fprintf('\n=========== ROW (RECORDING CHANNEL) vs COLUMN (STIM ELECTRODE) ===========\n');
fprintf('Share of the variance across the channel x electrode matrix:\n');
amps = unique(Row_amp);
for a = 1:numel(amps)
    rows = find(abs(Row_amp-amps(a)) < Amp_Tol_uA);
    M = Prof(rows,:).';
    if Exclude_StimChannel
        for j = 1:numel(rows)
            M(Recording_Channels == Row_stim(rows(j)),j) = NaN;
        end
    end
    [pr,pc,pi_] = decompose_matrix(M);
    fprintf('  %g uA: recording-channel effect %.0f%% | stim-electrode effect %.0f%% | interaction (local) %.0f%%\n', ...
        amps(a),pr,pc,pi_);
end
fprintf('(large "recording-channel" = some channels respond to everything; large\n');
fprintf(' "interaction" = the response moves with the stimulated electrode)\n');

%% ================================ FIGURES ==============================

rows_map = find(abs(Row_amp-Map_Amp_uA) < Amp_Tol_uA);
[~,o] = sortrows([Row_ds(rows_map),Row_stim(rows_map)]);
rows_map = rows_map(o);
nMap = numel(rows_map);
labels = arrayfun(@(r) sprintf('Ch%d (ds%d)',Row_stim(r),Row_ds(r)),rows_map,'UniformOutput',false);

if nMap == 0
    warning('No rows at Map_Amp_uA = %g - maps skipped.',Map_Amp_uA);
else

    % ---------- Figure A: probe maps ----------
    lo = min(min(Prof(rows_map,:))); hi = max(max(Prof(rows_map,:)));
    nCols = min(nMap,5); nR = ceil(nMap/nCols);
    figure('Color','w','Position',[40 60 1500 330*nR]);
    for k = 1:nMap
        r = rows_map(k);
        subplot(nR,nCols,k); hold on;
        scatter(xy(:,1),xy(:,2),110,Prof(r,:).','filled','s');
        [sx,sy] = ChnMap2Coord(Row_stim(r),Electrode_Type);
        plot(sx,sy,'kp','MarkerSize',16,'MarkerFaceColor','y');
        axis equal; xlim([-100 700]); ylim([-60 800]);
        caxis([lo hi]);
        title(labels{k},'FontWeight','bold');
        if mod(k-1,nCols) == 0, ylabel('y along shank (\mum)'); end
        xlabel('x across shanks (\mum)'); box off;
    end
    cb = colorbar; cb.Position(1) = 0.93; ylabel(cb,'baseline-corrected count');
    sgtitle(sprintf('Single-electrode response on the probe layout, %g \\muA (star = stimulated electrode)',Map_Amp_uA));

    % ---------- Figure B: profiles ----------
    cols = parula(nMap+1);
    ring_centres = (Ring_Edges_um(1:end-1)+Ring_Edges_um(2:end))/2;
    shank1 = find(Recording_Channels >= 1 & Recording_Channels <= 16);

    figure('Color','w','Position',[60 80 1250 480]);
    subplot(1,2,1); hold on;
    for k = 1:nMap
        r = rows_map(k);
        plot(Recording_Channels(shank1),Prof(r,shank1),'-o','Color',cols(k,:),'LineWidth',1.5, ...
            'MarkerSize',5,'DisplayName',labels{k});
        if ismember(Row_stim(r),Recording_Channels(shank1))
            plot(Row_stim(r),Prof(r,Recording_Channels == Row_stim(r)),'k^','MarkerSize',9, ...
                'MarkerFaceColor',cols(k,:),'HandleVisibility','off');
        end
    end
    xlabel('Recording channel on shank 1 (map position)','FontWeight','bold');
    ylabel('Baseline-corrected count','FontWeight','bold');
    title(sprintf('Along shank 1, %g \\muA (triangle = stimulated channel)',Map_Amp_uA),'FontWeight','bold');
    legend('Location','best','Box','off','FontSize',8); box off;

    subplot(1,2,2); hold on;
    for k = 1:nMap
        r = rows_map(k);
        plot(ring_centres,Ring_(r,:),'-o','Color',cols(k,:),'LineWidth',1.5,'MarkerSize',5, ...
            'DisplayName',labels{k});
    end
    xlabel('Distance of recording channel from the stim electrode (\mum)','FontWeight','bold');
    ylabel('Mean baseline-corrected count','FontWeight','bold');
    title('Response vs distance from the stimulated electrode','FontWeight','bold');
    legend('Location','best','Box','off','FontSize',8); box off;

    % ---------- Figure D: channel x electrode matrix ----------
    M = Prof(rows_map,:).';
    figure('Color','w','Position',[80 100 1100 600]);
    subplot(1,2,1);
    imagesc(1:nMap,Recording_Channels,M); set(gca,'YDir','normal');
    set(gca,'XTick',1:nMap,'XTickLabel',labels,'XTickLabelRotation',45);
    colorbar; box off;
    ylabel('Recording channel','FontWeight','bold');
    title(sprintf('Response matrix, %g \\muA',Map_Amp_uA),'FontWeight','bold');
    subplot(1,2,2);
    plot(mean(M,2,'omitnan'),Recording_Channels,'-k','LineWidth',1.5);
    ylim([min(Recording_Channels) max(Recording_Channels)]);
    xlabel('Mean over stim electrodes','FontWeight','bold'); ylabel('Recording channel');
    title('Recording-channel effect (row mean)','FontWeight','bold'); box off;
end

% ---------- Figure C: ALL / LOCAL / FAR vs stim channel, one panel per amp ----------
figure('Color','w','Position',[60 60 1500 450]);
for a = 1:numel(amps)
    subplot(1,numel(amps),a); hold on;
    rows = find(abs(Row_amp-amps(a)) < Amp_Tol_uA);
    nonAnchor = rows(Row_stim(rows) ~= Anchor_Channel);
    anchors   = rows(Row_stim(rows) == Anchor_Channel);
    [~,o] = sort(Row_stim(nonAnchor)); nonAnchor = nonAnchor(o);

    plot(Row_stim(nonAnchor),All_(nonAnchor),'-ok','LineWidth',2,'MarkerFaceColor','k','DisplayName','all channels');
    plot(Row_stim(nonAnchor),Local_(nonAnchor),'-ob','LineWidth',2,'MarkerFaceColor','b', ...
        'DisplayName',sprintf('local (<= %g \\mum)',Local_Radius_um));
    plot(Row_stim(nonAnchor),Far_(nonAnchor),'-or','LineWidth',2,'MarkerFaceColor','r','DisplayName','far');
    for r = anchors.'
        plot(Row_stim(r),All_(r),'kp','MarkerSize',13,'MarkerFaceColor','k','HandleVisibility','off');
        plot(Row_stim(r),Local_(r),'bp','MarkerSize',13,'MarkerFaceColor','b','HandleVisibility','off');
        plot(Row_stim(r),Far_(r),'rp','MarkerSize',13,'MarkerFaceColor','r','HandleVisibility','off');
    end
    xlabel('Stimulated electrode (map position)','FontWeight','bold');
    if a == 1, ylabel('Mean baseline-corrected count','FontWeight','bold'); end
    title(sprintf('%g \\muA (stars = anchor, both sessions)',amps(a)),'FontWeight','bold');
    if a == 1, legend('Location','best','Box','off'); end
    box off;
end

fprintf('\n============================================================\n');
fprintf('DONE (plots + tables only - nothing saved)\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end

function [pRow,pCol,pInt] = decompose_matrix(M)
% Additive decomposition of a (recording channel x stim electrode) matrix
% with NaNs allowed. Returns the percentage of total variance explained by
% the row effect, the column effect, and the remainder (interaction).
valid = ~isnan(M);
gm = mean(M(valid));
rowMean = mean(M,2,'omitnan');
colMean = mean(M,1,'omitnan');
R = repmat(rowMean,1,size(M,2));
C = repmat(colMean,size(M,1),1);
ssTot = sum((M(valid)-gm).^2);
ssRow = sum((R(valid)-gm).^2);
ssCol = sum((C(valid)-gm).^2);
if ssTot <= 0
    pRow = NaN; pCol = NaN; pInt = NaN; return;
end
pRow = 100*ssRow/ssTot;
pCol = 100*ssCol/ssTot;
pInt = max(0,100-pRow-pCol);
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