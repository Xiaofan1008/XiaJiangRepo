%% ========================================================================
% CLEAN Ch1 SINGLE-STIMULATION TRIALS (ALL SINGLE DATASETS OF ONE ANIMAL) - v1
%
% PURPOSE
%   Ch1 is the anchor electrode in every pair, so its single-pulse mean goes
%   into EVERY prediction. This script pools the Ch1 single trials of all
%   single datasets of ONE animal (every distance block), computes a per-trial
%   spike-count value, flags suspicious trials with the same rules for every
%   dataset, and prints / plots them together so you can decide the bad
%   trials in one go. SPIKE COUNTS ONLY (no linearity here - a single
%   electrode has nothing to be linear against).
%
% PER-TRIAL VALUES (same counting as SpikeTRatio: baseline-corrected count
% in analysis_win_ms, baseline scaled by window length), over Channel_Range:
%   Metric : MEAN over channels of the baseline-corrected count  (spikes / chn / trial)
%   MaxCh  : the LARGEST single-channel baseline-corrected count in the trial
%            (catches a one-channel artifact that the mean can hide)
%   Base   : MEAN over channels of the raw baseline-window count (spikes / chn)
%
% RULES (each can be switched on/off; a trial can fire several)
%   LOW       Metric < lowThreshold
%   HIGH      Metric > highThreshold
%   RESP-MAD  |Metric - median| > k * MAD        (two-sided, per group)
%   MAXCH     MaxCh  > maxChThreshold
%   BASE-ABS  Base   > baseAbsThreshold
%   BASE-MAD  Base   > median + k * MAD          (per group)
%   A group = one amplitude within one dataset (MAD_Scope = 'dataset'), or one
%   amplitude pooled over all datasets (MAD_Scope = 'pooled'). If a group's
%   MAD is 0 the MAD rules cannot be defined there and are skipped (noted).
%   The thresholds below are PLACEHOLDERS - set them after reading the
%   distribution table this script prints.
%
% EXISTING BAD TRIALS
%   Each dataset's saved <base_name>_BadTrials.mat (if any) is loaded and
%   those trials are shown as grey crosses. They are NOT counted in the
%   statistics when Exclude_ExistingBad_FromStats = true, and only trials
%   that are NEW (not already saved as bad) are listed as flagged. So after
%   you enter the bad trials and re-run, you see what is left to decide.
%
% OUTPUT (nothing is saved to disk)
%   1. distribution table per dataset x amplitude (to choose thresholds)
%   2. BLOCK CHECK: Ch1 mean +- SE per dataset and amplitude, with and
%      without the newly flagged trials (does Ch1 differ between blocks?
%      does the mean move when the flagged trials go?)
%   3. flagged-trial list grouped by dataset > amplitude, giving the set,
%      amplitude, RELATIVE trial number (ConditionTrialIndex) and absolute
%      trial id, the values and which rules fired
%   4. a ready-to-paste Overrides block for BadTrial_ManualConfirm_v1 for
%      each dataset
%   5. figures: one per amplitude (rows: Metric, MaxCh, Base; columns:
%      datasets) and one block-check figure
% ========================================================================

clear; clc;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

% All SINGLE datasets of ONE animal (every distance block). Add / remove
% entries as needed (copy a block and increment the index).
SingleDatasets = struct('spike_file',{},'dataset_folder',{});

SingleDatasets(1).spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1/Xia_Linearity_S1_100_400um_Single1.sp_xia_filtered.mat';
SingleDatasets(1).dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1';

SingleDatasets(2).spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1/Xia_Linearity_S1_500_700um_Single1.sp_xia_filtered.mat';
SingleDatasets(2).dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1';

% SingleDatasets(3).spike_file     = '.../Xia_Linearity_..._Single2/..._Single2.sp_xia_filtered.mat';
% SingleDatasets(3).dataset_folder = '.../Xia_Linearity_..._Single2';

Target_Channel = 1;              % the stimulation channel to clean (anchor electrode)
Channel_Range  = 1:64;           % channels counted (same idea as in SpikeTRatio)
Electrode_Type = 3;              % 0 rigid | 1 flexible | 2 four-shank | 3 64chn+32chn hybrid
FS_fallback    = 30000;

baseline_win_ms = [-50 -5];
analysis_win_ms = [0 20];

Amp_Tol_uA   = 0.001;
Amps_To_Show = [];               % [] = all amplitudes found; otherwise e.g. [5 10]

% ---- rules (thresholds are in spikes per channel per trial, PLACEHOLDERS) ----
Use_LowRule      = true;   lowThreshold     = 0;
Use_HighRule     = true;   highThreshold    = 5;
Use_RespMADRule  = true;   respMAD_k        = 5;
Use_MaxChRule    = true;   maxChThreshold   = 20;
Use_BaseAbsRule  = true;   baseAbsThreshold = 2;
Use_BaseMADRule  = true;   baseMAD_k        = 5;

MAD_Scope = 'dataset';           % 'dataset' (each amp within a dataset) or 'pooled' (each amp over all datasets)
Exclude_ExistingBad_FromStats = true;   % ignore already-saved bad trials when computing medians / MADs / block means

% ---- figure appearance ----
Label_All_Trials = true;         % write the relative trial number above EVERY point (black; flagged = bold)
Measures_To_Plot = [1 3];            % 1 = RESPONSE | 2 = MAX single chn | 3 = BASELINE   (e.g. [1 2 3] for all)
Figure_Layout    = 'combined';   % 'separate': one wide figure per amplitude AND measure (roomiest)
                                 % 'combined': one figure per amplitude, one row per measure
Fig_Width_Per_DS   = 760;        % figure width per dataset column (pixels)
Fig_Height_Per_Row = 560;        % figure height per row (pixels)
Marker_Size_Kept    = 6;         % blue dots
Marker_Size_Flagged = 8;         % red dots / grey crosses
Label_FontSize      = 8;         % trial numbers

%% =========================== INITIAL CHECKS ===========================

if isempty(SingleDatasets)
    error('SingleDatasets is empty.');
end
validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');
if ~ismember(lower(MAD_Scope),{'dataset','pooled'})
    error('MAD_Scope must be ''dataset'' or ''pooled''.');
end
MAD_Scope = lower(MAD_Scope);

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

nDS = numel(SingleDatasets);
scale_factor = diff(analysis_win_ms)/diff(baseline_win_ms);

fprintf('\n============================================================\n');
fprintf('CLEAN Ch%d SINGLE TRIALS (v1)\n',Target_Channel);
fprintf('============================================================\n');
fprintf('Datasets: %d | Channel_Range: %d channels | MAD scope: %s\n', ...
    nDS,numel(unique(Channel_Range)),MAD_Scope);
fprintf('Baseline [%g,%g) ms | Analysis [%g,%g) ms\n',baseline_win_ms,analysis_win_ms);
fprintf('Rules: LOW(<%g) %s | HIGH(>%g) %s | RESP-MAD(k=%g) %s | MAXCH(>%g) %s | BASE-ABS(>%g) %s | BASE-MAD(k=%g) %s\n', ...
    lowThreshold,onoff(Use_LowRule),highThreshold,onoff(Use_HighRule), ...
    respMAD_k,onoff(Use_RespMADRule),maxChThreshold,onoff(Use_MaxChRule), ...
    baseAbsThreshold,onoff(Use_BaseAbsRule),baseMAD_k,onoff(Use_BaseMADRule));

%% ====================== LOAD DATASETS, PER-TRIAL VALUES =================

T_ds = []; T_set = []; T_ai = []; T_amp = []; T_rel = []; T_abs = [];
T_metric = []; T_maxch = []; T_base = []; T_exist = false(0,1);

DSinfo = struct('label',{},'folder',{},'base_name',{},'set',{},'n_existing',{});

for d = 1:nDS

    fprintf('\n--- Loading dataset %d of %d ---\n',d,nDS);
    DS = load_linearity_dataset(SingleDatasets(d).spike_file,SingleDatasets(d).dataset_folder, ...
        Electrode_Type,FS_fallback,true);

    if DS.sim_stim ~= 1
        warning('Dataset %d (%s) has simultaneous_stim = %d (not a single-pulse dataset) - skipped.', ...
            d,DS.base_name,DS.sim_stim);
        continue;
    end

    set_idx = [];
    for ss = 1:DS.nSets
        sc = DS.uniqueComb(ss,DS.uniqueComb(ss,:) > 0);
        if numel(sc) == 1 && sc(1) == Target_Channel
            set_idx = ss;
            break;
        end
    end
    if isempty(set_idx)
        warning('Dataset %d (%s): no single-pulse set for Ch%d - skipped.',d,DS.base_name,Target_Channel);
        continue;
    end

    chans = unique(double(Channel_Range(:).'));
    chans = chans(chans >= 1 & chans <= DS.nMapPos);

    DSinfo(d).label     = DS.base_name;
    DSinfo(d).folder    = DS.dataset_folder;
    DSinfo(d).base_name = DS.base_name;
    DSinfo(d).set       = set_idx;
    DSinfo(d).n_existing = numel(DS.bad_trial_ids);

    for ai = 1:numel(DS.Amps)

        ids = find(DS.combClass == set_idx & DS.ampIdx == ai);   % absolute ids, in trial order
        if isempty(ids)
            continue;
        end
        amp_value = DS.Amps(ai);
        nT = numel(ids);
        t0 = DS.trig_ms(ids);

        corr_counts = zeros(nT,numel(chans));
        base_counts = zeros(nT,numel(chans));
        for c = 1:numel(chans)
            st = DS.sp_by_mappos{chans(c)};
            if isempty(st)
                continue;
            end
            b0 = count_less_than(st,t0+baseline_win_ms(1));
            b1 = count_less_than(st,t0+baseline_win_ms(2));
            a0 = count_less_than(st,t0+analysis_win_ms(1));
            a1 = count_less_than(st,t0+analysis_win_ms(2));
            bc = b1-b0;
            ac = a1-a0;
            corr_counts(:,c) = ac-bc*scale_factor;
            base_counts(:,c) = bc;
        end

        T_ds     = [T_ds;     repmat(d,nT,1)];                 %#ok<AGROW>
        T_set    = [T_set;    repmat(set_idx,nT,1)];           %#ok<AGROW>
        T_ai     = [T_ai;     repmat(ai,nT,1)];                %#ok<AGROW>
        T_amp    = [T_amp;    repmat(amp_value,nT,1)];         %#ok<AGROW>
        T_rel    = [T_rel;    (1:nT).'];                       %#ok<AGROW>
        T_abs    = [T_abs;    ids(:)];                         %#ok<AGROW>
        T_metric = [T_metric; mean(corr_counts,2)];            %#ok<AGROW>
        T_maxch  = [T_maxch;  max(corr_counts,[],2)];          %#ok<AGROW>
        T_base   = [T_base;   mean(base_counts,2)];            %#ok<AGROW>
        T_exist  = [T_exist;  ismember(ids(:),DS.bad_trial_ids(:))]; %#ok<AGROW>
    end

    fprintf('Dataset %d (%s): Ch%d = Set %d | %d trials | %d already saved as bad (whole dataset) | %d channels counted\n', ...
        d,DS.base_name,Target_Channel,set_idx,sum(T_ds == d),numel(DS.bad_trial_ids),numel(chans));

    clear DS;
end

if isempty(T_ds)
    error('No Ch%d single-pulse trials were found in any dataset.',Target_Channel);
end

T_amp = round(T_amp*1000)/1000;
all_amps = unique(T_amp(:).');
if ~isempty(Amps_To_Show)
    keep_amp = arrayfun(@(a) any(abs(Amps_To_Show(:)-a) < Amp_Tol_uA),T_amp);
    T_ds = T_ds(keep_amp); T_set = T_set(keep_amp); T_ai = T_ai(keep_amp); T_amp = T_amp(keep_amp);
    T_rel = T_rel(keep_amp); T_abs = T_abs(keep_amp); T_metric = T_metric(keep_amp);
    T_maxch = T_maxch(keep_amp); T_base = T_base(keep_amp); T_exist = T_exist(keep_amp);
    if isempty(T_ds)
        error('No trials match Amps_To_Show. Available amplitudes (uA): %s',num2str(all_amps));
    end
end
amps = unique(T_amp(:).');
nTrialsTotal = numel(T_ds);
ds_present = unique(T_ds(:).');

%% ============== GROUP STATISTICS (median / MAD, existing bad optional) ==

use_stat = true(nTrialsTotal,1);
if Exclude_ExistingBad_FromStats
    use_stat = ~T_exist;
end

grp_key = T_amp;                              % pooled: one group per amplitude
if strcmp(MAD_Scope,'dataset')
    grp_key = T_ds*1000+T_amp;                % dataset x amplitude (amp < 1000 uA)
end
[ugrp,~,grp_id] = unique(grp_key);
nGrp = numel(ugrp);

resp_med = nan(nGrp,1); resp_mad = nan(nGrp,1);
base_med = nan(nGrp,1); base_mad = nan(nGrp,1);
for g = 1:nGrp
    sel = grp_id == g & use_stat;
    if any(sel)
        resp_med(g) = median(T_metric(sel));
        resp_mad(g) = median(abs(T_metric(sel)-resp_med(g)));
        base_med(g) = median(T_base(sel));
        base_mad(g) = median(abs(T_base(sel)-base_med(g)));
    end
end

trial_resp_med = resp_med(grp_id); trial_resp_mad = resp_mad(grp_id);
trial_base_med = base_med(grp_id); trial_base_mad = base_mad(grp_id);

%% ============================ FLAG TRIALS ==============================

isLow     = Use_LowRule     & (T_metric < lowThreshold);
isHigh    = Use_HighRule    & (T_metric > highThreshold);
isRespMAD = Use_RespMADRule & (trial_resp_mad > 0) & ...
    (abs(T_metric-trial_resp_med) > respMAD_k*trial_resp_mad);
isMaxCh   = Use_MaxChRule   & (T_maxch > maxChThreshold);
isBaseAbs = Use_BaseAbsRule & (T_base > baseAbsThreshold);
isBaseMAD = Use_BaseMADRule & (trial_base_mad > 0) & ...
    (T_base > trial_base_med+baseMAD_k*trial_base_mad);

anyRule = isLow | isHigh | isRespMAD | isMaxCh | isBaseAbs | isBaseMAD;
isNew   = anyRule & ~T_exist;           % flagged and not yet saved as bad

%% ===================== 1. DISTRIBUTION TABLE ===========================

fprintf('\n------------------------------------------------------------\n');
fprintf('1. Distribution per dataset x amplitude (use it to set the thresholds)\n');
fprintf('   (Metric / MaxCh / Base in spikes per channel per trial; stats exclude\n');
fprintf('    already-saved bad trials: %s)\n',onoff(Exclude_ExistingBad_FromStats));
fprintf('------------------------------------------------------------\n');
fprintf('%-3s %-8s %4s %5s | %7s %6s %7s %7s | %7s %7s | %7s %6s %7s\n', ...
    'DS','Amp(uA)','N','Exist','MetMed','MetMAD','MetMin','MetMax','MaxMed','MaxMax','BasMed','BasMAD','BasMax');
for d = ds_present
    for a = amps
        sel = T_ds == d & T_amp == a;
        if ~any(sel), continue; end
        s2 = sel & use_stat;
        if ~any(s2), s2 = sel; end
        m = T_metric(s2);
        fprintf('%-3d %-8g %4d %5d | %7.3g %6.3g %7.3g %7.3g | %7.3g %7.3g | %7.3g %6.3g %7.3g\n', ...
            d,a,sum(sel),sum(T_exist(sel)), ...
            median(m),median(abs(m-median(m))),min(m),max(m), ...
            median(T_maxch(s2)),max(T_maxch(s2)), ...
            median(T_base(s2)),median(abs(T_base(s2)-median(T_base(s2)))),max(T_base(s2)));
    end
end
fprintf('(DS = dataset number in SingleDatasets)\n');
for d = ds_present
    fprintf('  DS %d = %s\n',d,DSinfo(d).label);
end

%% ===================== 2. BLOCK CHECK (print) ==========================

fprintf('\n------------------------------------------------------------\n');
fprintf('2. BLOCK CHECK - Ch%d response (Metric) per dataset, mean +- SE across trials\n',Target_Channel);
fprintf('   "all" = every trial except already-saved bad ones; "clean" = also without the\n');
fprintf('   newly flagged trials. If the datasets differ clearly, trial cleaning cannot fix it.\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-8s %-3s | %4s %8s %7s | %4s %8s %7s | %12s\n', ...
    'Amp(uA)','DS','nAll','meanAll','SE','nCln','meanCln','SE','clean/first DS');

blk_mean_all = nan(numel(ds_present),numel(amps));
blk_se_all   = nan(numel(ds_present),numel(amps));
blk_mean_cln = nan(numel(ds_present),numel(amps));
blk_se_cln   = nan(numel(ds_present),numel(amps));
blk_n_all    = zeros(numel(ds_present),numel(amps));
blk_n_cln    = zeros(numel(ds_present),numel(amps));

for ai_ = 1:numel(amps)
    a = amps(ai_);
    first_mean = NaN;
    for di = 1:numel(ds_present)
        d = ds_present(di);
        sel_all = T_ds == d & T_amp == a & ~T_exist;
        sel_cln = sel_all & ~isNew;
        if ~any(sel_all), continue; end
        blk_n_all(di,ai_) = sum(sel_all);
        blk_n_cln(di,ai_) = sum(sel_cln);
        blk_mean_all(di,ai_) = mean(T_metric(sel_all));
        blk_se_all(di,ai_)   = std(T_metric(sel_all))/sqrt(max(1,sum(sel_all)));
        if any(sel_cln)
            blk_mean_cln(di,ai_) = mean(T_metric(sel_cln));
            blk_se_cln(di,ai_)   = std(T_metric(sel_cln))/sqrt(max(1,sum(sel_cln)));
        end
        if isnan(first_mean) && ~isnan(blk_mean_cln(di,ai_))
            first_mean = blk_mean_cln(di,ai_);
            rel_txt = '1 (ref)';
        elseif ~isnan(first_mean) && first_mean ~= 0
            rel_txt = sprintf('%.2f',blk_mean_cln(di,ai_)/first_mean);
        else
            rel_txt = 'n/a';
        end
        fprintf('%-8g %-3d | %4d %8.3g %7.2g | %4d %8.3g %7.2g | %12s\n', ...
            a,d,blk_n_all(di,ai_),blk_mean_all(di,ai_),blk_se_all(di,ai_), ...
            blk_n_cln(di,ai_),blk_mean_cln(di,ai_),blk_se_cln(di,ai_),rel_txt);
    end
end

%% ===================== 3. FLAGGED TRIAL LIST ===========================

fprintf('\n------------------------------------------------------------\n');
fprintf('3. Flag counts by rule (all trials; "new" = not yet saved as bad)\n');
fprintf('------------------------------------------------------------\n');
fprintf('  LOW      (<%g):               %d\n',lowThreshold,sum(isLow));
fprintf('  HIGH     (>%g):               %d\n',highThreshold,sum(isHigh));
fprintf('  RESP-MAD (k=%g):              %d\n',respMAD_k,sum(isRespMAD));
fprintf('  MAXCH    (>%g):               %d\n',maxChThreshold,sum(isMaxCh));
fprintf('  BASE-ABS (>%g):               %d\n',baseAbsThreshold,sum(isBaseAbs));
fprintf('  BASE-MAD (k=%g):              %d\n',baseMAD_k,sum(isBaseMAD));
fprintf('  Flagged by at least one rule: %d of %d (%.1f%%)\n',sum(anyRule),nTrialsTotal,100*sum(anyRule)/nTrialsTotal);
fprintf('  ... of which NEW (to decide): %d | already saved as bad: %d\n',sum(isNew),sum(anyRule & T_exist));

no_mad = find(~(resp_mad > 0));
if (Use_RespMADRule || Use_BaseMADRule) && ~isempty(no_mad)
    fprintf('  Note: MAD = 0 (MAD rules skipped) in %d group(s).\n',numel(no_mad));
end

fprintf('\n------------------------------------------------------------\n');
fprintf('Flagged trials (NEW only), grouped by dataset > amplitude\n');
fprintf('   rel = relative trial number within the condition (what you type into\n');
fprintf('   BadTrial_ManualConfirm); abs = absolute trial id in the dataset\n');
fprintf('------------------------------------------------------------\n');

if ~any(isNew)
    fprintf('  (no new flagged trials)\n');
end

for d = ds_present
    if ~any(isNew & T_ds == d), continue; end
    fprintf('\n[DS %d] %s  (Set %d, Ch%d)\n',d,DSinfo(d).label,DSinfo(d).set,Target_Channel);
    for a = amps
        sel = find(isNew & T_ds == d & T_amp == a);
        if isempty(sel), continue; end
        n_cond = sum(T_ds == d & T_amp == a);
        fprintf('  Amp %g uA: %d new flagged of %d trials\n',a,numel(sel),n_cond);
        for k = sel(:).'
            rules_txt = {};
            if isLow(k),     rules_txt{end+1} = 'LOW'; end      %#ok<SAGROW>
            if isHigh(k),    rules_txt{end+1} = 'HIGH'; end     %#ok<SAGROW>
            if isRespMAD(k), rules_txt{end+1} = 'RESP-MAD'; end %#ok<SAGROW>
            if isMaxCh(k),   rules_txt{end+1} = 'MAXCH'; end    %#ok<SAGROW>
            if isBaseAbs(k), rules_txt{end+1} = 'BASE-ABS'; end %#ok<SAGROW>
            if isBaseMAD(k), rules_txt{end+1} = 'BASE-MAD'; end %#ok<SAGROW>
            fprintf('    rel %3d (abs %4d) | Metric %7.3g | MaxCh %7.3g | Base %6.3g | %s\n', ...
                T_rel(k),T_abs(k),T_metric(k),T_maxch(k),T_base(k),strjoin(rules_txt,', '));
        end
    end
end

%% ===================== 4. OVERRIDES TO PASTE ===========================

fprintf('\n------------------------------------------------------------\n');
fprintf('4. Overrides for BadTrial_ManualConfirm_v1 (NEW flagged trials, as a patch:\n');
fprintf('   ResetAllFirst = false adds them to whatever is already saved).\n');
fprintf('   Edit the list first - these are only the rule-flagged trials, not final decisions.\n');
fprintf('------------------------------------------------------------\n');

for d = ds_present
    if ~any(isNew & T_ds == d), continue; end
    fprintf('\n%% [DS %d] %s\n',d,DSinfo(d).label);
    fprintf('%% matFilePath = ''%s'';\n', ...
        fullfile(DSinfo(d).folder,sprintf('%s_TrialSpikeCounts.mat',DSinfo(d).base_name)));
    fprintf('Overrides = {\n');
    fprintf('  %% Set  Amp  PTD(ms)  ResetAllFirst  AddBad  RemoveBad\n');
    for a = amps
        sel = find(isNew & T_ds == d & T_amp == a);
        if isempty(sel), continue; end
        fprintf('    %d,   %g,   0,       false,         [%s],  [];\n', ...
            DSinfo(d).set,a,strjoin(arrayfun(@(x) sprintf('%d',x),sort(T_rel(sel)).','UniformOutput',false),','));
    end
    fprintf('};\n');
end

%% ===================== 5. FIGURES ======================================

col_good = [0.20 0.45 0.80];
col_new  = [0.85 0.15 0.15];
col_old  = [0.55 0.55 0.55];

row_title = {'Metric: mean baseline-corrected count per channel', ...
    'MaxCh: largest single-channel count','Base: mean baseline count per channel'};
panel_name = {'RESPONSE (mean per chn)','MAX single chn','BASELINE (pre-stim, per chn)'};

Measures_To_Plot = unique(Measures_To_Plot(:).');
if any(~ismember(Measures_To_Plot,1:3))
    error('Measures_To_Plot may only contain 1, 2, 3.');
end
Figure_Layout = lower(Figure_Layout);
if ~ismember(Figure_Layout,{'separate','combined'})
    error('Figure_Layout must be ''separate'' or ''combined''.');
end

for a = amps

    if strcmp(Figure_Layout,'separate')
        fig_groups = num2cell(Measures_To_Plot);      % one figure per measure
    else
        fig_groups = {Measures_To_Plot};              % one figure, one row per measure
    end

    for fg = 1:numel(fig_groups)

        meas = fig_groups{fg};
        nRows = numel(meas);
        fig_name = sprintf('Ch%d single trials - %g uA - %s',Target_Channel,a, ...
            strjoin(arrayfun(@(m) panel_name{m},meas,'UniformOutput',false),' + '));
        fig_w = min(1900,Fig_Width_Per_DS*numel(ds_present)+110);
        fig_h = min(1000,Fig_Height_Per_Row*nRows+90);
        figure('Color','w','Name',fig_name,'Position',[40 40 fig_w fig_h]);

        for di = 1:numel(ds_present)
            d = ds_present(di);
            sel = T_ds == d & T_amp == a;
            if ~any(sel), continue; end
            x = T_rel(sel);
            vals = {T_metric(sel),T_maxch(sel),T_base(sel)};
            ex = T_exist(sel);
            nw = isNew(sel);
            gi = grp_id(find(sel,1));

            for rr = 1:nRows
                r = meas(rr);
                ax = subplot(nRows,numel(ds_present),(rr-1)*numel(ds_present)+di);
                hold(ax,'on');
                y = vals{r};

                good = ~ex & ~nw;
                plot(ax,x(good),y(good),'o','Color',col_good,'MarkerFaceColor',col_good,'MarkerSize',Marker_Size_Kept);
                if any(ex)
                    plot(ax,x(ex),y(ex),'x','Color',col_old,'MarkerSize',Marker_Size_Flagged,'LineWidth',1.5);
                end
                if any(nw)
                    plot(ax,x(nw),y(nw),'o','Color',col_new,'MarkerFaceColor',col_new,'MarkerSize',Marker_Size_Flagged);
                end

                if Label_All_Trials
                    % black number just above every point; flagged trials in bold
                    for k = 1:numel(x)
                        if nw(k), fw = 'bold'; else, fw = 'normal'; end
                        text(ax,x(k),y(k),sprintf('%d',x(k)),'Color',[0 0 0],'FontSize',Label_FontSize, ...
                            'FontWeight',fw,'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                            'Clipping','on');
                    end
                end

                if r == 1
                    if ~isnan(resp_med(gi))
                        yline(ax,resp_med(gi),'-','Color',[0.4 0.4 0.4]);
                        if Use_RespMADRule && resp_mad(gi) > 0
                            yline(ax,resp_med(gi)+respMAD_k*resp_mad(gi),':','Color',[0.3 0.3 0.3]);
                            yline(ax,resp_med(gi)-respMAD_k*resp_mad(gi),':','Color',[0.3 0.3 0.3]);
                        end
                    end
                    if Use_HighRule, yline(ax,highThreshold,'--','Color',[0.8 0.4 0.0]); end
                    if Use_LowRule,  yline(ax,lowThreshold,'--','Color',[0.8 0.4 0.0]);  end
                elseif r == 2
                    if Use_MaxChRule, yline(ax,maxChThreshold,'--','Color',[0.8 0.4 0.0]); end
                else
                    if ~isnan(base_med(gi)), yline(ax,base_med(gi),'-','Color',[0.4 0.4 0.4]); end
                    if Use_BaseMADRule && base_mad(gi) > 0
                        yline(ax,base_med(gi)+baseMAD_k*base_mad(gi),':','Color',[0.3 0.3 0.3]);
                    end
                    if Use_BaseAbsRule, yline(ax,baseAbsThreshold,'--','Color',[0.8 0.4 0.0]); end
                end

                title(ax,sprintf('%s  |  DS %d (%s)  |  %g uA  |  n = %d', ...
                    panel_name{r},d,block_label(DSinfo(d).label),a,sum(sel)),'FontSize',9,'Interpreter','none');
                if di == 1, ylabel(ax,row_title{r},'FontSize',8); end
                if rr == nRows, xlabel(ax,'Condition trial # (relative)'); end

                % one tick per trial + faint vertical grid: follow any dot down to its number
                nx = max(x);
                set(ax,'XTick',1:nx,'XGrid','on','YGrid','off','GridAlpha',0.12, ...
                    'GridColor',[0 0 0],'FontSize',7,'Layer','bottom');
                xlim(ax,[0 nx+1]);
                ylim(ax,pad_ylim(ax));
                box(ax,'off');
            end
        end
        sgtitle(sprintf('Ch%d single trials - %g uA   (blue = kept, red = NEW flagged, grey x = already saved bad;  grey line = median, dotted = k*MAD, dashed = manual limits)', ...
            Target_Channel,a),'FontSize',9);
    end
end

% ---- block-check figure ----
figure('Color','w','Name',sprintf('Ch%d block check',Target_Channel),'Position',[100 100 760 480]);
axB = axes; hold(axB,'on');
cmap = lines(max(numel(ds_present),3));
off = ((1:numel(ds_present))-(numel(ds_present)+1)/2)*0.22;
leg_h = gobjects(numel(ds_present),1);
for di = 1:numel(ds_present)
    d = ds_present(di);
    for ai_ = 1:numel(amps)
        a = amps(ai_);
        sel = T_ds == d & T_amp == a & ~T_exist;
        if ~any(sel), continue; end
        xj = ai_+off(di)+(rand(sum(sel),1)-0.5)*0.08;
        plot(axB,xj,T_metric(sel),'.','Color',[cmap(di,:) 0.35],'MarkerSize',8);
        h1 = errorbar(axB,ai_+off(di),blk_mean_all(di,ai_),blk_se_all(di,ai_),'s', ...
            'Color',cmap(di,:),'MarkerFaceColor',cmap(di,:),'LineWidth',1.4,'CapSize',8);
        if ~isnan(blk_mean_cln(di,ai_))
            errorbar(axB,ai_+off(di)+0.07,blk_mean_cln(di,ai_),blk_se_cln(di,ai_),'o', ...
                'Color',cmap(di,:),'MarkerFaceColor','w','LineWidth',1.4,'CapSize',8);
        end
        if ai_ == 1 || ~isgraphics(leg_h(di))
            leg_h(di) = h1;
        end
    end
end
set(axB,'XTick',1:numel(amps),'XTickLabel',arrayfun(@(v) sprintf('%g uA',v),amps,'UniformOutput',false));
xlim(axB,[0.4 numel(amps)+0.6]);
ylabel(axB,'Ch1 response (spikes / chn / trial)');
title(axB,sprintf('Ch%d block check: filled = all kept trials, open = also without NEW flagged (mean +- SE)',Target_Channel), ...
    'FontSize',10);
valid_leg = isgraphics(leg_h);
if any(valid_leg)
    legend(axB,leg_h(valid_leg),arrayfun(@(d) sprintf('DS %d',d),ds_present(valid_leg),'UniformOutput',false), ...
        'Location','best','Box','off');
end
box(axB,'off');

fprintf('\nDone. Nothing was saved. Enter the final bad trials with BadTrial_ManualConfirm_v1, then re-run.\n');

%% =========================== LOCAL FUNCTIONS ==========================

function yl = pad_ylim(ax)
% current y-limits with a little headroom so the numbers above the top points are not clipped
yl = ylim(ax);
span = diff(yl);
yl = [yl(1)-0.03*span, yl(2)+0.08*span];
end

function lbl = block_label(base_name)
% short block name from the dataset name, e.g. '100_400um' (falls back to the full name)
tok = regexp(base_name,'(\d+_\d+um)','tokens','once');
if isempty(tok), lbl = base_name; else, lbl = tok{1}; end
end

function s = onoff(flag)
if flag, s = 'ON'; else, s = 'OFF'; end
end

function n = count_less_than(st,t)
% Number of elements of the sorted vector st that are < t, for every
% element of t (vectorised binary search).
n_st = numel(st);
lo = zeros(size(t));
hi = repmat(n_st,size(t));
while any(lo(:) < hi(:))
    active = lo < hi;
    mid = floor((lo+hi)/2);
    idx = min(max(mid+1,1),n_st);
    val = st(idx);
    go_right = active & (val < t);
    go_left  = active & ~(val < t);
    lo(go_right) = mid(go_right)+1;
    hi(go_left)  = mid(go_left);
end
n = lo;
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