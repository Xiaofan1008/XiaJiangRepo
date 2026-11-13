%% ========================================================================
% CHECK THE SIM / SEQ (PAIRED) TRIALS FOR ONE AMPLITUDE - v2
%
% PURPOSE
%   v2 works on ONE amplitude at a time (Amp_To_Clean) across all distances:
%   rule limits are set for that amplitude, and the paired bad trials are read
%   from / (later) saved to a file tagged with the amplitude, e.g.
%   <base_name>_10uA_BadTrials.mat. The singles keep using their normal
%   <base_name>_BadTrials.mat (a 10 uA single trial only feeds 10 uA).
%
%   Run this AFTER the single datasets are cleaned (Ch1 + ChB bad-trial
%   scripts). It flags suspicious paired trials with the same rules as the
%   single scripts and shows, for every distance, the trials, ACTUAL vs
%   PREDICTED, and the whole linearity curve, with and without the flagged
%   trials. Everything is computed on exactly the channels and the cleaned
%   singles that SpikeTRatio_v4 uses, so what you clean here is what gets
%   analysed there.
%
% PER-TRIAL VALUES (over the channel list of that pair + amplitude)
%   Metric : mean baseline-corrected count per channel      (1 = RESPONSE)
%   MaxCh  : largest single-channel baseline-corrected count (2 = MAX single chn)
%   Base   : mean raw baseline-window count per channel      (3 = BASELINE)
%   Channel list = Channel_Range, or (if Use_RespondingChannels and the paired
%   dataset has a _RespondingChannels.mat) the union of the responding channels
%   of ALL PTDs and both stim orders of that pair + amplitude, intersected with
%   Channel_Range - the same rule as SpikeTRatio_v4.
%
% PREDICTED
%   mean(single A trials) + mean(single B trials), same channels, singles read
%   WITHOUT their saved bad trials. One value per pair + amplitude, the same
%   for Sim and Seq and for both stim orders.
%
% RULES (spike counts; same as the single scripts, thresholds are PLACEHOLDERS)
%   LOW, HIGH, RESP-MAD, MAXCH, BASE-ABS, BASE-MAD. MAD groups = one
%   condition (stim order x amplitude x PTD) of one dataset.
%
% LINEARITY LIMIT FLAG (display only by default)
%   Each paired trial gets  ratio_i = Metric_i / Predicted.  A trial outside
%   [Lin_Ratio_Low, Lin_Ratio_High] gets a black RING and is listed. It is NOT
%   a bad-trial flag unless you set Lin_Counts_As_Flag = true (not recommended:
%   removing trials because of their linearity pulls the ratio toward what you
%   expect). Use the ring to decide which trials to look at.
%
% "NOW" vs "CLEAN"
%   now   = every trial except those already saved as bad in the paired file
%   clean = also without the trials flagged (NEW) in this run
%
% OUTPUT (nothing is saved)
%   1. condition table: trials, flags, predicted
%   2. flagged trials with the effect of removing each one on the ratio
%   3. trials outside the linearity limits
%   4. ACTUAL / PREDICTED / RATIO per distance, amplitude and PTD (now vs clean,
%      bootstrap SD, odd/even split of the ratio)
%   5. a ready-made template (Dist, Order, PTD, trials of THIS amplitude) per dataset for the
%      amplitude-specific entry script. Each bracket holds the trials ALREADY saved as bad
%      plus the NEW flagged ones, i.e. the suggested FINAL list for that condition.
%   6. one figure per dataset and distance (both stim orders on the same panels)
%   7. two curve figures: actual + predicted vs distance, ratio vs distance
% ========================================================================

clear; clc;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

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

Channel_Range = 1:64;
Use_RespondingChannels = true;   % as in SpikeTRatio_v4
Electrode_Type = 3;              % 0 rigid | 1 flexible | 2 four-shank | 3 64chn+32chn hybrid
FS_fallback    = 30000;

baseline_win_ms = [-50 -5];
analysis_win_ms = [0 20];

% AMP change
Amp_To_Clean      = 5;          % the ONE amplitude (uA) cleaned in this run
BadTrial_Tag      = 'auto';      % tag of the PAIRED bad-trial file: 'auto' -> '<amp>uA' (e.g. '10uA',
                                 % file <base_name>_10uA_BadTrials.mat); '' = the old untagged file;
                                 % or any text, e.g. '5uA'. Singles always use the untagged file.
Distances_To_Show = [];          % [] = all; otherwise e.g. [100 300 500]
Distance_Tol_um   = 1;
Amp_Tol_uA = 0.001;
PTD_Tol_ms = 0.001;

% ---- spike-count rules (spikes per channel per trial, PLACEHOLDERS) ----
Use_LowRule      = false;  lowThreshold     = 0;
Use_HighRule     = true;   highThreshold    = 5;
Use_RespMADRule  = true;   respMAD_k        = 5;
Use_MaxChRule    = true;   maxChThreshold   = 20;
Use_BaseAbsRule  = true;   baseAbsThreshold = 3;
Use_BaseMADRule  = true;   baseMAD_k        = 5;
Exclude_ExistingBad_FromStats = true;

% ---- linearity limit flag (per trial: Metric / Predicted) ----
Use_LinearityLimit = true;
Lin_Ratio_Low      = 0.2;        % ring if Metric/Predicted < this
Lin_Ratio_High     = 2.0;        % ring if Metric/Predicted > this
Lin_Min_Predicted  = 0.2;        % no per-trial ratio where Predicted is not above this
Lin_Counts_As_Flag = false;      % true: ringed trials also count as NEW flagged (not recommended)

% ---- uncertainty ----
Boot_N    = 300;                 % bootstrap repeats for the error bars (0 = off)
Boot_Seed = 1;

% ---- figure appearance ----
Trial_Measure    = [1 3];      % ROWS of each distance figure: 1 RESPONSE | 2 MAX single chn | 3 BASELINE | 4 TRIAL RATIO (Metric/Predicted)
Label_All_Trials = true;         % black number next to every point (order 1 above, order 2 below)
Show_Median_Lines = true;
Marker_Size_Kept    = 6;
Marker_Size_Flagged = 8;
Label_FontSize      = 8;
Plot_DistanceFigures = true;
Plot_CurveFigures    = true;
Plot_AmpCurve        = false;        % needs several amplitudes - not used in the one-amplitude check
Amp_Curve_Measure    = 'ratio';      % 'ratio' (Actual/Predicted, 1 = linear) | 'LI' (linearity index, 0 = linear) | 'both'
Amp_Curve_ShowNow    = false;        % true: also draw the 'now' curves (dashed) next to the clean ones
Plot_TrialRatioVsDistance = true;    % every trial's ratio_i against distance, with mean / median / pooled ratio
Ratio_YLim = [];                     % [] = automatic; otherwise e.g. [0 3] for the ratio rows / figure

%% =========================== INITIAL CHECKS ===========================

if isempty(DatasetPairs), error('DatasetPairs is empty.'); end
validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');
Trial_Measure = unique(Trial_Measure(:).');
if isempty(Trial_Measure) || any(~ismember(Trial_Measure,1:4)), error('Trial_Measure may only contain 1, 2, 3, 4.'); end
if isempty(Amp_To_Clean) || ~isscalar(Amp_To_Clean) || ~isnumeric(Amp_To_Clean)
    error('Amp_To_Clean must be ONE amplitude in uA, e.g. 10.');
end
Amps_To_Show = Amp_To_Clean;
if strcmpi(BadTrial_Tag,'auto')
    paired_tag = sprintf('%guA',Amp_To_Clean);
else
    paired_tag = BadTrial_Tag;
end
if Lin_Ratio_Low >= Lin_Ratio_High, error('Lin_Ratio_Low must be below Lin_Ratio_High.'); end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>
rng(Boot_Seed);

scale_factor = diff(analysis_win_ms)/diff(baseline_win_ms);

fprintf('\n============================================================\n');
fprintf('CHECK SIM / SEQ TRIALS, ACTUAL vs PREDICTED, LINEARITY (v1)\n');
fprintf('============================================================\n');
fprintf('Dataset-pairs: %d | baseline [%g,%g) ms | analysis [%g,%g) ms\n', ...
    numel(DatasetPairs),baseline_win_ms,analysis_win_ms);
fprintf('Amplitude cleaned: %g uA | paired bad-trial file tag: ''%s'' (singles: untagged)\n',Amp_To_Clean,paired_tag);
fprintf('Rules: LOW(<%g) %s | HIGH(>%g) %s | RESP-MAD(k=%g) %s | MAXCH(>%g) %s | BASE-ABS(>%g) %s | BASE-MAD(k=%g) %s\n', ...
    lowThreshold,onoff(Use_LowRule),highThreshold,onoff(Use_HighRule), ...
    respMAD_k,onoff(Use_RespMADRule),maxChThreshold,onoff(Use_MaxChRule), ...
    baseAbsThreshold,onoff(Use_BaseAbsRule),baseMAD_k,onoff(Use_BaseMADRule));
fprintf('Linearity limit: %s (ratio_i outside [%g, %g] = ring; counts as a flag: %s)\n', ...
    onoff(Use_LinearityLimit),Lin_Ratio_Low,Lin_Ratio_High,onoff(Lin_Counts_As_Flag));

%% ================ LOAD DATASETS, BUILD ONE ROW PER CONDITION ===========

Cond = struct('dp',{},'si',{},'stims',{},'dist',{},'amp',{},'ptd',{},'ord',{}, ...
    'rel',{},'abs',{},'metric',{},'maxch',{},'base',{},'exist',{},'list',{}, ...
    'kA',{},'kB',{},'vA',{},'vB',{},'pred',{});

SingleStore = struct('corr',{},'bad',{});   % single-trial count matrices, by key
single_key_map = containers.Map('KeyType','char','ValueType','double');

DPinfo = struct('label',{},'folder',{},'base_name',{});

for dp = 1:numel(DatasetPairs)

    entry = DatasetPairs(dp);
    fprintf('\n============================================================\n');
    fprintf('DATASET-PAIR %d/%d\n',dp,numel(DatasetPairs));
    fprintf('============================================================\n');

    if ~isfile(entry.simseq_spike_file),   error('Pair %d: paired spike file missing:\n%s',dp,entry.simseq_spike_file); end
    if ~isfolder(entry.simseq_dataset_folder), error('Pair %d: paired folder missing:\n%s',dp,entry.simseq_dataset_folder); end
    if ~isfile(entry.single_spike_file),   error('Pair %d: single spike file missing:\n%s',dp,entry.single_spike_file); end
    if ~isfolder(entry.single_dataset_folder), error('Pair %d: single folder missing:\n%s',dp,entry.single_dataset_folder); end

    fprintf('\n--- Loading paired dataset ---\n');
    Paired = load_linearity_dataset(entry.simseq_spike_file,entry.simseq_dataset_folder, ...
        Electrode_Type,FS_fallback,true,paired_tag);
    fprintf('\n--- Loading single dataset ---\n');
    Single = load_linearity_dataset(entry.single_spike_file,entry.single_dataset_folder, ...
        Electrode_Type,FS_fallback,true,'');

    DPinfo(dp).label     = Paired.base_name;
    DPinfo(dp).folder    = Paired.dataset_folder;
    DPinfo(dp).base_name = Paired.base_name;

    chans = unique(double(Channel_Range(:).'));
    chans = chans(chans >= 1 & chans <= min(Paired.nMapPos,Single.nMapPos));

    % ---- responding-channel file (optional) ----
    RespondingConfirmed = [];
    have_resp_file = false;
    if Use_RespondingChannels
        responding_path = fullfile(entry.simseq_dataset_folder,sprintf('%s_RespondingChannels.mat',Paired.base_name));
        if isfile(responding_path)
            RespLoad = load(responding_path,'RespondingConfirmed');
            if isfield(RespLoad,'RespondingConfirmed') && ~isempty(RespLoad.RespondingConfirmed)
                RespondingConfirmed = RespLoad.RespondingConfirmed;
                have_resp_file = true;
                fprintf('Responding-channel file found (%d conditions).\n',numel(RespondingConfirmed));
            end
        else
            fprintf('No responding-channel file for this pair - Channel_Range only.\n');
        end
    end

    for si = 1:Paired.nSets

        stims = Paired.uniqueComb(si,Paired.uniqueComb(si,:) > 0);
        if numel(stims) ~= 2, continue; end

        dist_um = NaN;
        try
            [dist_um,~] = ChnPairDistance(stims(1),stims(2),Electrode_Type);
        catch
        end
        if ~isempty(Distances_To_Show) && (isnan(dist_um) || ~any(abs(Distances_To_Show(:)-dist_um) <= Distance_Tol_um))
            continue;
        end
        if stims(1) < stims(2), ord = 1; else, ord = 2; end

        for ai = 1:numel(Paired.Amps)

            amp_value = round(Paired.Amps(ai)*1000)/1000;
            if ~isempty(Amps_To_Show) && ~any(abs(Amps_To_Show(:)-amp_value) < Amp_Tol_uA), continue; end

            % ---- channel list of this pair + amplitude (shared by all PTDs and both orders) ----
            list = chans;
            if have_resp_file
                target_pair = sort(double(stims(:).'));
                resp_union = []; n_used = 0;
                for k = 1:numel(RespondingConfirmed)
                    sc = double(RespondingConfirmed(k).stimChannels(:).');
                    sc = sort(sc(sc > 0));
                    if isequal(sc,target_pair) && abs(RespondingConfirmed(k).amp-amp_value) < Amp_Tol_uA
                        resp_union = union(resp_union,double(RespondingConfirmed(k).channels(:).'));
                        n_used = n_used+1;
                    end
                end
                if n_used == 0
                    warning('Pair %d, Set %d (%s), %g uA: no responding-channel entry - Channel_Range only.', ...
                        dp,si,format_order(stims),amp_value);
                else
                    list = intersect(chans,resp_union);
                    if isempty(list)
                        warning('Pair %d, Set %d (%s), %g uA: responding channels and Channel_Range share nothing - skipped.', ...
                            dp,si,format_order(stims),amp_value);
                        continue;
                    end
                end
            end
            cols = ismember(chans,list);

            % ---- the two singles (cached count matrices) ----
            kk = [0 0]; vv = {[],[]};
            for w = 1:2
                [ss,aiS] = find_single_condition(Single,stims(w),amp_value,Amp_Tol_uA);
                if isempty(ss)
                    warning('Pair %d, Set %d (%s), %g uA: no single-pulse match for Ch%d - predicted NaN.', ...
                        dp,si,format_order(stims),amp_value,stims(w));
                    continue;
                end
                key = sprintf('%d_%d_%d',dp,ss,aiS);
                if isKey(single_key_map,key)
                    idx_k = single_key_map(key);
                else
                    ids_s = find(Single.combClass == ss & Single.ampIdx == aiS);
                    corr_s = trial_counts(Single,ids_s,chans,baseline_win_ms,analysis_win_ms,scale_factor);
                    SingleStore(end+1).corr = corr_s;                        %#ok<AGROW>
                    SingleStore(end).bad = ismember(ids_s(:),Single.bad_trial_ids(:));
                    idx_k = numel(SingleStore);
                    single_key_map(key) = idx_k;
                end
                good_s = ~SingleStore(idx_k).bad;
                if ~any(good_s)
                    warning('Pair %d, Set %d, %g uA: all Ch%d single trials are saved as bad - predicted NaN.',dp,si,amp_value,stims(w));
                    continue;
                end
                kk(w) = idx_k;
                vv{w} = mean(SingleStore(idx_k).corr(good_s,cols),2);
            end
            if all(kk > 0), pred = mean(vv{1})+mean(vv{2}); else, pred = NaN; end

            for pi = 1:numel(Paired.PTDs_ms)

                ids = find(Paired.combClass == si & Paired.ampIdx == ai & Paired.ptdIdx == pi);
                if isempty(ids), continue; end

                [corr_p,base_p] = trial_counts(Paired,ids,chans,baseline_win_ms,analysis_win_ms,scale_factor);

                Cond(end+1).dp = dp; %#ok<AGROW>
                n = numel(Cond);
                Cond(n).si    = si;
                Cond(n).stims = stims;
                Cond(n).dist  = dist_um;
                Cond(n).amp   = amp_value;
                Cond(n).ptd   = Paired.PTDs_ms(pi);
                Cond(n).ord   = ord;
                Cond(n).rel   = (1:numel(ids)).';
                Cond(n).abs   = ids(:);
                Cond(n).metric = mean(corr_p(:,cols),2);
                Cond(n).maxch  = max(corr_p(:,cols),[],2);
                Cond(n).base   = mean(base_p(:,cols),2);
                Cond(n).exist  = ismember(ids(:),Paired.bad_trial_ids(:));
                Cond(n).list   = list;
                Cond(n).kA = kk(1); Cond(n).kB = kk(2);
                Cond(n).vA = vv{1}; Cond(n).vB = vv{2};
                Cond(n).pred = pred;
            end
        end
    end

    fprintf('Pair %d: %d conditions so far | %d paired trials saved as bad | %d single trials saved as bad\n', ...
        dp,sum([Cond.dp] == dp),numel(Paired.bad_trial_ids),numel(Single.bad_trial_ids));
    clear Paired Single;
end

nCond = numel(Cond);
if nCond == 0, error('No paired conditions match the selection.'); end

C_dp   = [Cond.dp].';
C_si   = [Cond.si].';
C_amp  = [Cond.amp].';
C_ptd  = [Cond.ptd].';
C_ord  = [Cond.ord].';
C_dist = [Cond.dist].';
C_dkey = C_dist; C_dkey(isnan(C_dkey)) = -1;        % NaN distance -> key -1
C_nchn = arrayfun(@(c) numel(c.list),Cond).';
C_ntr  = arrayfun(@(c) numel(c.rel),Cond).';

amps_all = unique(C_amp).';
ptds_all = unique(C_ptd).';
dists_all = unique(C_dkey).';
dp_all = unique(C_dp).';

%% ============================ FLAG TRIALS ==============================

Flag = struct('low',{},'high',{},'rmad',{},'maxch',{},'babs',{},'bmad',{},'ratio',{},'lin',{}, ...
    'any',{},'new',{},'med',{});

for r = 1:nCond
    m = Cond(r).metric; b = Cond(r).base; mx = Cond(r).maxch; ex = Cond(r).exist;
    use = true(size(m));
    if Exclude_ExistingBad_FromStats, use = ~ex; end
    if ~any(use), use = true(size(m)); end

    rmed = median(m(use)); rmad = median(abs(m(use)-rmed));
    bmed = median(b(use)); bmad = median(abs(b(use)-bmed));

    F.low   = Use_LowRule     & (m < lowThreshold);
    F.high  = Use_HighRule    & (m > highThreshold);
    F.rmad  = Use_RespMADRule & (rmad > 0) & (abs(m-rmed) > respMAD_k*rmad);
    F.maxch = Use_MaxChRule   & (mx > maxChThreshold);
    F.babs  = Use_BaseAbsRule & (b > baseAbsThreshold);
    F.bmad  = Use_BaseMADRule & (bmad > 0) & (b > bmed+baseMAD_k*bmad);

    F.ratio = nan(size(m));
    if Cond(r).pred > Lin_Min_Predicted, F.ratio = m/Cond(r).pred; end
    F.lin = Use_LinearityLimit & ~isnan(F.ratio) & (F.ratio < Lin_Ratio_Low | F.ratio > Lin_Ratio_High);

    F.any = F.low | F.high | F.rmad | F.maxch | F.babs | F.bmad;
    if Lin_Counts_As_Flag, F.any = F.any | F.lin; end
    F.new = F.any & ~ex;
    F.med = rmed;
    Flag(r) = F;
end

% keep masks (cell per row)
keepNow   = cell(nCond,1);
keepClean = cell(nCond,1);
for r = 1:nCond
    keepNow{r}   = ~Cond(r).exist;
    keepClean{r} = ~Cond(r).exist & ~Flag(r).new;
end

nNew = sum(arrayfun(@(f) sum(f.new),Flag));
nAny = sum(arrayfun(@(f) sum(f.any),Flag));
nTr  = sum(C_ntr);
nLin = sum(arrayfun(@(f) sum(f.lin & ~f.any),Flag));

fprintf('\nTrials: %d | flagged by a spike rule%s: %d | NEW (not yet saved bad): %d | ring only (linearity limit, not flagged): %d\n', ...
    nTr,ternary(Lin_Counts_As_Flag,' or linearity limit',''),nAny,nNew,nLin);

%% ===================== 1. CONDITION TABLE ==============================

fprintf('\n------------------------------------------------------------\n');
fprintf('1. Conditions (one row per pair-dataset x stim order x amplitude x PTD)\n');
fprintf('   Pred = predicted (spikes/chn/trial, same for Sim, Seq and both orders)\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-3s %-9s %-14s %-8s %-7s %4s %5s %4s %4s | %7s %7s | %6s %6s\n', ...
    'DP','Dist(um)','Order','Amp(uA)','PTD(ms)','N','Exist','NEW','Ring','Pred','MetMed','MetMax','nChn');
[~,ordR] = sortrows([C_dp,C_dkey,C_ord,C_amp,C_ptd]);
prev_dp = NaN; prev_d = NaN;
for ii = 1:numel(ordR)
    r = ordR(ii);
    if ii > 1 && (C_dp(r) ~= prev_dp || C_dkey(r) ~= prev_d), fprintf('\n'); end
    prev_dp = C_dp(r); prev_d = C_dkey(r);
    fprintf('%-3d %-9s %-14s %-8g %-7g %4d %5d %4d %4d | %7.3g %7.3g | %6.3g %6d\n', ...
        C_dp(r),dist_label(C_dist(r)),format_order(Cond(r).stims),C_amp(r),C_ptd(r),C_ntr(r), ...
        sum(Cond(r).exist),sum(Flag(r).new),sum(Flag(r).lin & ~Flag(r).new), ...
        Cond(r).pred,median(Cond(r).metric),max(Cond(r).metric),C_nchn(r));
end
for dp = dp_all, fprintf('  DP %d = %s\n',dp,DPinfo(dp).label); end

%% ===================== 2. FLAGGED TRIALS + INFLUENCE ===================

fprintf('\n------------------------------------------------------------\n');
fprintf('2. Flagged trials (NEW only), grouped by pair-dataset > distance > order > amplitude > PTD\n');
fprintf('   rel = relative trial number in the condition (what you enter); abs = absolute trial id\n');
fprintf('   ratio_i = Metric/Predicted. "pooled ratio" = ratio of this distance/amp/PTD (both orders of this\n');
fprintf('   dataset) now -> without THIS trial only.\n');
fprintf('------------------------------------------------------------\n');
if nNew == 0, fprintf('  (no new flagged trials)\n'); end

for ii = 1:numel(ordR)
    r = ordR(ii);
    f = Flag(r);
    idxs = find(f.new);
    if isempty(idxs), continue; end
    grp_rows = find(C_dp == C_dp(r) & C_dkey == C_dkey(r) & C_amp == C_amp(r) & C_ptd == C_ptd(r));
    S0 = group_stats(Cond,grp_rows,keepNow,0);
    fprintf('\n[DP %d] %s | %s | %s [%s, Set %d] | %g uA | PTD %g ms  (%d new flagged of %d)\n', ...
        C_dp(r),DPinfo(C_dp(r)).label,dist_label(C_dist(r)),format_order(Cond(r).stims),shape_name(Cond(r).ord),C_si(r), ...
        C_amp(r),C_ptd(r),numel(idxs),C_ntr(r));
    for t = idxs(:).'
        rt = {};
        if f.low(t),   rt{end+1} = 'LOW'; end       %#ok<SAGROW>
        if f.high(t),  rt{end+1} = 'HIGH'; end      %#ok<SAGROW>
        if f.rmad(t),  rt{end+1} = 'RESP-MAD'; end  %#ok<SAGROW>
        if f.maxch(t), rt{end+1} = 'MAXCH'; end     %#ok<SAGROW>
        if f.babs(t),  rt{end+1} = 'BASE-ABS'; end  %#ok<SAGROW>
        if f.bmad(t),  rt{end+1} = 'BASE-MAD'; end  %#ok<SAGROW>
        if f.lin(t),   rt{end+1} = 'LIN-LIMIT'; end %#ok<SAGROW>
        kk_tmp = keepNow; kk_tmp{r}(t) = false;
        S1 = group_stats(Cond,grp_rows,kk_tmp,0);
        fprintf('    rel %3d (abs %4d) | Metric %7.3g | MaxCh %7.3g | Base %6.3g | ratio_i %6.2f | pooled ratio %.2f -> %.2f | %s\n', ...
            Cond(r).rel(t),Cond(r).abs(t),Cond(r).metric(t),Cond(r).maxch(t),Cond(r).base(t),f.ratio(t), ...
            S0.ratio,S1.ratio,strjoin(rt,', '));
    end
end

%% ===================== 3. LINEARITY-LIMIT TRIALS =======================

fprintf('\n------------------------------------------------------------\n');
fprintf('3. Trials outside the linearity limits (ratio_i = Metric/Predicted outside [%g, %g])\n',Lin_Ratio_Low,Lin_Ratio_High);
fprintf('   display only (ring on the figures)%s\n',ternary(Lin_Counts_As_Flag,' - but they ALSO count as flagged (Lin_Counts_As_Flag = true)',''));
fprintf('------------------------------------------------------------\n');
if ~any(arrayfun(@(f) any(f.lin),Flag)), fprintf('  (none)\n'); end
for ii = 1:numel(ordR)
    r = ordR(ii);
    idxs = find(Flag(r).lin & ~Cond(r).exist);
    if isempty(idxs), continue; end
    fprintf('[DP %d] %s | %s [%s, Set %d] | %g uA | PTD %g ms | predicted %.3g: ', ...
        C_dp(r),dist_label(C_dist(r)),format_order(Cond(r).stims),shape_name(Cond(r).ord),C_si(r),C_amp(r),C_ptd(r),Cond(r).pred);
    fprintf('%s\n',strjoin(arrayfun(@(t) sprintf('rel %d (%.2f)',Cond(r).rel(t),Flag(r).ratio(t)),idxs(:).','UniformOutput',false),', '));
end

%% ===================== 4. ACTUAL / PREDICTED / RATIO ===================

fprintf('\n------------------------------------------------------------\n');
fprintf('4. ACTUAL vs PREDICTED per pair-dataset x distance x amplitude x PTD (both stim orders pooled)\n');
fprintf('   Actual / Pred in spikes per channel per trial; r = Actual/Pred; +-SD = bootstrap SD over trials\n');
fprintf('   (singles resampled too); odd / even = ratio from odd / even clean trials only.\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-3s %-9s %-8s %-7s | %4s | %7s | %7s %7s | %6s %6s %6s | %6s %6s\n', ...
    'DP','Dist(um)','Amp(uA)','PTD(ms)','nTr','Pred','ActNow','ActCln','r now','r cln','+-SD','odd','even');

[gkeys,~,~] = unique([C_dp,C_dkey,C_amp,C_ptd],'rows');
nG = size(gkeys,1);
G = struct('dp',{},'dkey',{},'amp',{},'ptd',{},'rows',{},'now',{},'cln',{},'odd',{},'even',{});
for g = 1:nG
    rows = find(C_dp == gkeys(g,1) & C_dkey == gkeys(g,2) & C_amp == gkeys(g,3) & C_ptd == gkeys(g,4));
    G(g).dp = gkeys(g,1); G(g).dkey = gkeys(g,2); G(g).amp = gkeys(g,3); G(g).ptd = gkeys(g,4);
    G(g).rows = rows;
    G(g).now  = group_stats(Cond,rows,keepNow,Boot_N);
    G(g).cln  = group_stats(Cond,rows,keepClean,Boot_N);
    [G(g).odd,G(g).even] = group_split(Cond,rows,keepClean);
end
[~,ordG] = sortrows(gkeys(:,[1 2 3 4]));
prev_dp = NaN; prev_d = NaN;
for ii = 1:nG
    g = ordG(ii);
    if ii > 1 && (G(g).dp ~= prev_dp || G(g).dkey ~= prev_d), fprintf('\n'); end
    prev_dp = G(g).dp; prev_d = G(g).dkey;
    fprintf('%-3d %-9s %-8g %-7g | %4d | %7.3g | %7.3g %7.3g | %6.2f %6.2f %6.2f | %6.2f %6.2f\n', ...
        G(g).dp,dist_label(G(g).dkey),G(g).amp,G(g).ptd,G(g).cln.n_trials,G(g).cln.pred, ...
        G(g).now.actual,G(g).cln.actual,G(g).now.ratio,G(g).cln.ratio,G(g).cln.ratio_sd,G(g).odd,G(g).even);
end

%% ===================== 5. TEMPLATE FOR THE ENTRY SCRIPT ================

fprintf('\n------------------------------------------------------------\n');
fprintf('5. TEMPLATE for the %g uA entry script (paired bad-trial file tag ''%s'')\n',Amp_To_Clean,paired_tag);
fprintf('   One row per condition that EXISTS (a PTD may have only one stim order - e.g. Sim).\n');
fprintf('   Key = Dist (um, as in the figure title), Order (1 = circle, 2 = diamond), PTD (0 = Sim, else Seq).\n');
fprintf('   Each bracket = trials ALREADY saved as bad + NEW flagged = the suggested FINAL list for that\n');
fprintf('   condition (a saved condition is REPLACED by what you enter). Edit inside the brackets only.\n');
fprintf('   Pairs with unknown distance are not listed.\n');
fprintf('------------------------------------------------------------\n');
for dp = dp_all
    rows_dp = find(C_dp == dp & C_dkey > 0);
    if isempty(rows_dp), continue; end
    fprintf('\n%% [DP %d] %s\n',dp,DPinfo(dp).label);
    fprintf('%% matFilePath = ''%s'';\n',fullfile(DPinfo(dp).folder,sprintf('%s_TrialSpikeCounts.mat',DPinfo(dp).base_name)));
    fprintf('Amp = %g;     %% tag ''%s''\n',Amp_To_Clean,paired_tag);
    fprintf('bad = {\n  %% Dist  Order  PTD   trials (%g uA)\n',Amp_To_Clean);
    tkeys = unique([C_dkey(rows_dp),C_ord(rows_dp),C_ptd(rows_dp)],'rows');
    prev_dist = NaN;
    for q = 1:size(tkeys,1)
        rr = rows_dp(C_dkey(rows_dp) == tkeys(q,1) & C_ord(rows_dp) == tkeys(q,2) & abs(C_ptd(rows_dp)-tkeys(q,3)) < PTD_Tol_ms);
        tl_ = Cond(rr(1)).rel(Cond(rr(1)).exist | Flag(rr(1)).new);
        cell_txt = ['[' strjoin(arrayfun(@(x) sprintf('%d',x),tl_(:).','UniformOutput',false),' ') ']'];
        if ~isnan(prev_dist) && tkeys(q,1) ~= prev_dist, fprintf('\n'); end
        prev_dist = tkeys(q,1);
        fprintf('    %g,   %d,   %g,   %s;   %% %s (%s)  %s  [%d saved + %d new]\n',tkeys(q,1),tkeys(q,2),tkeys(q,3), ...
            cell_txt,format_order(Cond(rr(1)).stims),shape_name(tkeys(q,2)),short_ptd(tkeys(q,3)), ...
            sum(Cond(rr(1)).exist),sum(Flag(rr(1)).new));
    end
    fprintf('};\n');
end

%% ===================== 6. FIGURES: ONE PER DATASET x DISTANCE ==========
% Columns = conditions that exist (Sim, Seq order 1, Seq order 2); rows = the
% measures in Trial_Measure; the last column (spanning all rows) = actual vs predicted.

col_good = [0.20 0.45 0.80];
col_new  = [0.85 0.15 0.15];
col_old  = [0.55 0.55 0.55];
col_pred = [0.10 0.10 0.10];
col_ptd  = lines(max(numel(ptds_all),3));

row_title  = {'Metric: mean baseline-corrected count per channel', ...
    'MaxCh: largest single-channel count','Base: mean baseline count per channel', ...
    'ratio_i = Metric / Predicted  (1 = prediction)'};
panel_name = {'RESPONSE','MAX single chn','BASELINE','TRIAL RATIO'};

if Plot_DistanceFigures

    nR = numel(Trial_Measure);

    for dp = dp_all
        for dk = dists_all

            rows_here = find(C_dp == dp & C_dkey == dk);
            if isempty(rows_here), continue; end
            ptds_here = unique(C_ptd(rows_here)).';
            cond_cols = unique([C_ptd(rows_here),C_ord(rows_here)],'rows');
            nCC = size(cond_cols,1);
            nPb = numel(ptds_here);

            fig = figure('Color','w','Units','normalized','Position',[0.01 0.04 0.98 0.86], ...
                'Name',sprintf('DP %d | %s | %g uA',dp,dist_label(dk),Amp_To_Clean));
            tl = tiledlayout(fig,nR,nCC+1,'TileSpacing','compact','Padding','compact');

            for ir = 1:nR
                mm = Trial_Measure(ir);
                for ic = 1:nCC
                    pt = cond_cols(ic,1); oo = cond_cols(ic,2);
                    ax = nexttile(tl,(ir-1)*(nCC+1)+ic);
                    hold(ax,'on');
                    rr = rows_here(abs(C_ptd(rows_here)-pt) < PTD_Tol_ms & C_ord(rows_here) == oo);
                    draw_trial_panel(ax,Cond,Flag,rr,mm,col_good,col_new,col_old,col_pred, ...
                        Marker_Size_Kept,Marker_Size_Flagged,Label_FontSize,Label_All_Trials,Show_Median_Lines, ...
                        Use_HighRule,highThreshold,Use_LowRule,lowThreshold,Use_MaxChRule,maxChThreshold, ...
                        Use_BaseAbsRule,baseAbsThreshold,Use_LinearityLimit,Lin_Ratio_Low,Lin_Ratio_High);
                    if mm == 4 && ~isempty(Ratio_YLim), ylim(ax,Ratio_YLim); end
                    title(ax,{sprintf('%s | %s',panel_name{mm},short_ptd(pt)),cond_text(Cond,rr)}, ...
                        'FontSize',9,'Interpreter','none');
                    if ic == 1, ylabel(ax,row_title{mm},'FontSize',8); end
                    if ir == nR, xlabel(ax,'Condition trial # (relative)'); end
                end
            end

            % ---- actual vs predicted (spans all rows) ----
            ax = nexttile(tl,nCC+1,[nR 1]);
            hold(ax,'on');
            gs = find([G.dp] == dp & [G.dkey] == dk);
            if ~isempty(gs)
                gp = G(gs(1));
                vals = nan(1+nPb,2); sds = nan(1+nPb,2); rtxt = cell(1+nPb,1);
                vals(1,:) = gp.cln.pred; sds(1,:) = gp.cln.pred_sd;
                for ip = 1:nPb
                    gi = gs(abs([G(gs).ptd]-ptds_here(ip)) < PTD_Tol_ms);
                    if isempty(gi), continue; end
                    vals(1+ip,:) = [G(gi(1)).now.actual, G(gi(1)).cln.actual];
                    sds(1+ip,:)  = [G(gi(1)).now.actual_sd, G(gi(1)).cln.actual_sd];
                    rtxt{1+ip} = sprintf('r %.2f -> %.2f',G(gi(1)).now.ratio,G(gi(1)).cln.ratio);
                end
                hb = bar(ax,1:(1+nPb),vals,'grouped','BarWidth',0.8);
                hb(1).FaceColor = [0.35 0.35 0.35]; hb(2).FaceColor = [0.75 0.75 0.75];
                xo = [-0.14 0.14];
                for k = 1:2
                    errorbar(ax,(1:(1+nPb))+xo(k),vals(:,k),sds(:,k),'k.','CapSize',6);
                end
                ytop = max(vals(:)+sds(:),[],'omitnan');
                if isempty(ytop) || isnan(ytop), ytop = 1; end
                for ip = 1:nPb
                    if ~isempty(rtxt{1+ip})
                        text(ax,1+ip,max(vals(1+ip,:)+sds(1+ip,:),[],'omitnan'),rtxt{1+ip},'HorizontalAlignment','center', ...
                            'VerticalAlignment','bottom','FontSize',9,'FontWeight','bold');
                    end
                end
                set(ax,'XTick',1:(1+nPb),'XTickLabel',[{'Predicted'},arrayfun(@(p_) short_ptd(p_),ptds_here,'UniformOutput',false)], ...
                    'FontSize',9,'XTickLabelRotation',0);
                ylim(ax,[min(0,min(vals(:)-sds(:),[],'omitnan'))*1.1, ytop*1.25+eps]);
                legend(ax,hb,{'now','without NEW flagged'},'Location','best','Box','off','FontSize',8);
            end
            title(ax,{'ACTUAL vs PREDICTED',sprintf('%s | %g uA',dist_label(dk),Amp_To_Clean)},'FontSize',9,'Interpreter','none');
            ylabel(ax,'spikes / chn / trial','FontSize',8);
            box(ax,'off');

            title(tl,sprintf(['DP %d (%s) | %s | %g uA   (blue = kept, red = NEW flagged, grey x = already saved bad, black ring = outside linearity limits;  ' ...
                'grey line = median, solid black line = predicted / ratio 1; each stim order has its own panel)'], ...
                dp,DPinfo(dp).label,dist_label(dk),Amp_To_Clean),'FontSize',9,'Interpreter','none');
        end
    end
end

%% ===================== 7. CURVE FIGURES ================================

if Plot_CurveFigures

    [ukeys,~,~] = unique([C_dkey,C_amp,C_ptd],'rows');
    GA = struct('dkey',{},'amp',{},'ptd',{},'now',{},'cln',{});
    for g = 1:size(ukeys,1)
        rows = find(C_dkey == ukeys(g,1) & C_amp == ukeys(g,2) & C_ptd == ukeys(g,3));
        GA(g).dkey = ukeys(g,1); GA(g).amp = ukeys(g,2); GA(g).ptd = ukeys(g,3);
        GA(g).now = group_stats(Cond,rows,keepNow,Boot_N);
        GA(g).cln = group_stats(Cond,rows,keepClean,Boot_N);
    end

    nAc = numel(amps_all);

    % ---- (a) actual + predicted vs distance ----
    figure('Color','w','Name','Curve: actual and predicted vs distance','Position',[40 60 max(900,min(1900,520*nAc+100)) 520]);
    for ia = 1:nAc
        a = amps_all(ia);
        ax = subplot(1,nAc,ia); hold(ax,'on');
        hl = gobjects(0); ll = {};
        for ip = 1:numel(ptds_all)
            sel = find([GA.amp] == a & [GA.ptd] == ptds_all(ip) & [GA.dkey] > 0);
            if isempty(sel), continue; end
            [~,o] = sort([GA(sel).dkey]); sel = sel(o);
            xd = [GA(sel).dkey];
            an = arrayfun(@(g) g.now.actual,GA(sel)); asn = arrayfun(@(g) g.now.actual_sd,GA(sel));
            ac = arrayfun(@(g) g.cln.actual,GA(sel)); asc = arrayfun(@(g) g.cln.actual_sd,GA(sel));
            h1 = errorbar(ax,xd,an,asn,'-o','Color',col_ptd(ip,:),'MarkerFaceColor',col_ptd(ip,:),'LineWidth',1.6,'MarkerSize',7,'CapSize',6);
            errorbar(ax,xd,ac,asc,'--o','Color',col_ptd(ip,:),'MarkerFaceColor','w','LineWidth',1.2,'MarkerSize',7,'CapSize',6);
            hl(end+1) = h1; ll{end+1} = ['Actual ' ptd_label(ptds_all(ip))]; %#ok<SAGROW>
            if ip == 1
                pn = arrayfun(@(g) g.cln.pred,GA(sel)); ps = arrayfun(@(g) g.cln.pred_sd,GA(sel));
                h2 = errorbar(ax,xd,pn,ps,'-s','Color',col_pred,'MarkerFaceColor',col_pred,'LineWidth',1.6,'MarkerSize',6,'CapSize',6);
                hl(end+1) = h2; ll{end+1} = 'Predicted'; %#ok<SAGROW>
            end
        end
        title(ax,sprintf('%g uA',a),'FontSize',10);
        xlabel(ax,'Distance (um)');
        if ia == 1, ylabel(ax,'spikes / chn / trial'); end
        if ~isempty(hl), legend(ax,hl,ll,'Location','best','Box','off','FontSize',8); end
        grid(ax,'on'); ax.GridAlpha = 0.15; box(ax,'off');
    end
    sgtitle('Actual and predicted vs distance (solid = now, dashed/open = also without NEW flagged; error bar = bootstrap SD over trials; pooled over stim orders and dataset-pairs)','FontSize',9);

    % ---- (b) ratio vs distance ----
    figure('Color','w','Name','Curve: ratio actual/predicted vs distance','Position',[60 80 max(900,min(1900,520*nAc+100)) 520]);
    for ia = 1:nAc
        a = amps_all(ia);
        ax = subplot(1,nAc,ia); hold(ax,'on');
        hl = gobjects(0); ll = {};
        for ip = 1:numel(ptds_all)
            sel = find([GA.amp] == a & [GA.ptd] == ptds_all(ip) & [GA.dkey] > 0);
            if isempty(sel), continue; end
            [~,o] = sort([GA(sel).dkey]); sel = sel(o);
            xd = [GA(sel).dkey];
            rn = arrayfun(@(g) g.now.ratio,GA(sel)); rsn = arrayfun(@(g) g.now.ratio_sd,GA(sel));
            rc = arrayfun(@(g) g.cln.ratio,GA(sel)); rsc = arrayfun(@(g) g.cln.ratio_sd,GA(sel));
            h1 = errorbar(ax,xd,rn,rsn,'-o','Color',col_ptd(ip,:),'MarkerFaceColor',col_ptd(ip,:),'LineWidth',1.6,'MarkerSize',7,'CapSize',6);
            errorbar(ax,xd,rc,rsc,'--o','Color',col_ptd(ip,:),'MarkerFaceColor','w','LineWidth',1.2,'MarkerSize',7,'CapSize',6);
            hl(end+1) = h1; ll{end+1} = ptd_label(ptds_all(ip)); %#ok<SAGROW>
        end
        yline(ax,1,':k','LineWidth',1.3);
        title(ax,sprintf('%g uA',a),'FontSize',10);
        xlabel(ax,'Distance (um)');
        if ia == 1, ylabel(ax,'Ratio = Actual / Predicted  (1 = linear)'); end
        if ~isempty(hl), legend(ax,hl,ll,'Location','best','Box','off','FontSize',8); end
        grid(ax,'on'); ax.GridAlpha = 0.15; box(ax,'off');
    end
    sgtitle('Ratio vs distance (solid = now, dashed/open = also without NEW flagged; error bar = bootstrap SD over trials, singles resampled too)','FontSize',9);

    % ---- (c) every trial's ratio vs distance ----
    if Plot_TrialRatioVsDistance
        nP2 = numel(ptds_all);
        figure('Color','w','Name','Curve: per-trial ratio vs distance','Position',[80 100 max(900,min(1900,480*nAc+100)) min(1000,330*nP2+110)]);
        for ip = 1:nP2
            for ia = 1:nAc
                a = amps_all(ia);
                ax = subplot(nP2,nAc,(ip-1)*nAc+ia); hold(ax,'on');
                rows = find(C_amp == a & abs(C_ptd-ptds_all(ip)) < PTD_Tol_ms & C_dkey > 0);
                dvals = unique(C_dkey(rows)).';
                for r = rows(:).'
                    f = Flag(r);
                    if C_ord(r) == 1, mk = 'o'; xo = -8; else, mk = 'd'; xo = 8; end
                    xj = C_dist(r)+xo+(rand(numel(f.ratio),1)-0.5)*6;
                    ex = Cond(r).exist; nw = f.new; good = ~ex & ~nw & ~isnan(f.ratio);
                    plot(ax,xj(good),f.ratio(good),mk,'Color',[col_good 0.55],'MarkerFaceColor',col_good,'MarkerSize',4);
                    sel2 = nw & ~isnan(f.ratio);
                    if any(sel2), plot(ax,xj(sel2),f.ratio(sel2),mk,'Color',col_new,'MarkerFaceColor','none','MarkerSize',6,'LineWidth',1.2); end
                end
                mean_c = nan(size(dvals)); med_c = nan(size(dvals)); pool_c = nan(size(dvals));
                for q = 1:numel(dvals)
                    rq = rows(C_dkey(rows) == dvals(q));
                    vals = [];
                    for r = rq(:).'
                        k = keepClean{r} & ~isnan(Flag(r).ratio);
                        vals = [vals; Flag(r).ratio(k)]; %#ok<AGROW>
                    end
                    if ~isempty(vals), mean_c(q) = mean(vals); med_c(q) = median(vals); end
                    ga = find([GA.dkey] == dvals(q) & [GA.amp] == a & abs([GA.ptd]-ptds_all(ip)) < PTD_Tol_ms,1);
                    if ~isempty(ga), pool_c(q) = GA(ga).cln.ratio; end
                end
                h_mean = plot(ax,dvals,mean_c,'ks','MarkerFaceColor','k','MarkerSize',7);
                h_med  = plot(ax,dvals,med_c,'k_','MarkerSize',16,'LineWidth',2);
                h_pool = plot(ax,dvals,pool_c,'-','Color',[0.85 0.4 0.0],'LineWidth',1.6);
                yline(ax,1,':k','LineWidth',1.2);
                if Use_LinearityLimit
                    yline(ax,Lin_Ratio_Low,':','Color',[0.6 0.1 0.6]);
                    yline(ax,Lin_Ratio_High,':','Color',[0.6 0.1 0.6]);
                end
                if ~isempty(Ratio_YLim), ylim(ax,Ratio_YLim); end
                title(ax,sprintf('%s | %g uA',ptd_label(ptds_all(ip)),a),'FontSize',9,'Interpreter','none');
                xlabel(ax,'Distance (um)');
                if ia == 1, ylabel(ax,'ratio_i = Metric / Predicted'); end
                if ip == 1 && ia == 1
                    legend(ax,[h_mean h_med h_pool],{'mean of trials','median of trials','pooled ratio (actual/predicted, clean)'}, ...
                        'Location','best','Box','off','FontSize',7);
                end
                grid(ax,'on'); ax.GridAlpha = 0.15; box(ax,'off');
            end
        end
        sgtitle('Every trial ratio_i vs distance: blue = kept trial, red ring = NEW flagged (not in clean), circle / diamond = stim order; black = mean and median of the clean trials; orange = pooled ratio','FontSize',9);
    end

    % ---- (d) linearity vs amplitude: one line per distance, one panel per PTD ----
    if Plot_AmpCurve && numel(amps_all) > 1
        switch lower(Amp_Curve_Measure)
            case 'ratio', ac_list = {'ratio'};
            case 'li',    ac_list = {'LI'};
            otherwise,    ac_list = {'ratio','LI'};
        end
        dvals_all = unique([GA.dkey]);
        dvals_all = dvals_all(dvals_all > 0);
        cm_d = parula(numel(dvals_all)+1);

        for im = 1:numel(ac_list)
            meas_name = ac_list{im};
            fig_a = figure('Color','w','Units','normalized','Position',[0.05 0.12 0.9 0.62], ...
                'Name',sprintf('Curve: %s vs amplitude',meas_name));
            tla = tiledlayout(fig_a,1,numel(ptds_all),'TileSpacing','compact','Padding','compact');
            for ip = 1:numel(ptds_all)
                ax = nexttile(tla,ip); hold(ax,'on');
                hl = gobjects(0); ll = {};
                for id = 1:numel(dvals_all)
                    sel = find([GA.dkey] == dvals_all(id) & abs([GA.ptd]-ptds_all(ip)) < PTD_Tol_ms);
                    if isempty(sel), continue; end
                    [~,o] = sort([GA(sel).amp]); sel = sel(o);
                    xa = [GA(sel).amp];
                    rc = arrayfun(@(g) g.cln.ratio,GA(sel));    sc_ = arrayfun(@(g) g.cln.ratio_sd,GA(sel));
                    rn = arrayfun(@(g) g.now.ratio,GA(sel));    sn_ = arrayfun(@(g) g.now.ratio_sd,GA(sel));
                    if strcmp(meas_name,'LI')
                        yc = (rc-1)./(rc+1);  ec = 2*sc_./(rc+1).^2;
                        yn = (rn-1)./(rn+1);  en = 2*sn_./(rn+1).^2;
                    else
                        yc = rc; ec = sc_; yn = rn; en = sn_;
                    end
                    if Amp_Curve_ShowNow
                        errorbar(ax,xa,yn,en,'--o','Color',cm_d(id,:),'MarkerFaceColor','w','LineWidth',1,'MarkerSize',6,'CapSize',5, ...
                            'HandleVisibility','off');
                    end
                    h1 = errorbar(ax,xa,yc,ec,'-o','Color',cm_d(id,:),'MarkerFaceColor',cm_d(id,:),'LineWidth',1.6,'MarkerSize',7,'CapSize',6);
                    hl(end+1) = h1; ll{end+1} = sprintf('%g um',dvals_all(id)); %#ok<SAGROW>
                end
                if strcmp(meas_name,'LI')
                    yline(ax,0,':k','LineWidth',1.3); ylim(ax,[-1 1]);
                    ylab = 'Linearity index (Actual - Predicted)/(Actual + Predicted)   (0 = linear)';
                else
                    yline(ax,1,':k','LineWidth',1.3);
                    ylab = 'Ratio = Actual / Predicted   (1 = linear)';
                end
                set(ax,'XTick',amps_all);
                xlim(ax,[min(amps_all)-0.05*range(amps_all)-0.5, max(amps_all)+0.05*range(amps_all)+0.5]);
                title(ax,ptd_label(ptds_all(ip)),'FontSize',10);
                xlabel(ax,'Amplitude (uA)');
                if ip == 1, ylabel(ax,ylab,'FontSize',8); end
                if ~isempty(hl)
                    lg = legend(ax,hl,ll,'Location','best','Box','off','FontSize',8);
                    title(lg,'Distance');
                end
                grid(ax,'on'); ax.GridAlpha = 0.15; box(ax,'off');
            end
            title(tla,sprintf('%s vs amplitude - one line per distance (clean trials; error bar = bootstrap SD%s)', ...
                meas_name,ternary(Amp_Curve_ShowNow,'; dashed/open = now','')),'FontSize',9);
        end
    end
end

fprintf('\nDone. Nothing was saved.\n');

%% =========================== LOCAL FUNCTIONS ==========================

function txt = cond_text(Cond,rr)
% 'Ch1 -> Ch5 (Set 5)' for the condition(s) in rr
if isempty(rr)
    txt = '(no trials)';
else
    txt = sprintf('%s (Set %d)',format_order(Cond(rr(1)).stims),Cond(rr(1)).si);
end
end

function txt = short_ptd(p)
if abs(p) < 1e-6, txt = 'Sim'; else, txt = sprintf('Seq %g ms',p); end
end

function txt = shape_name(ord)
if ord == 1, txt = 'circle'; else, txt = 'diamond'; end
end

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end

function txt = dist_label(d)
if isnan(d) || d < 0, txt = 'n/a'; else, txt = sprintf('%g um',d); end
end

function txt = ptd_label(p)
if abs(p) < 1e-6, txt = 'Sim (PTD 0)'; else, txt = sprintf('Seq (PTD %g)',p); end
end

function txt = order_name(Cond,rows,ord)
txt = '(none)';
for k = rows(:).'
    if Cond(k).ord == ord, txt = format_order(Cond(k).stims); return; end
end
end

function [corr_counts,base_counts] = trial_counts(DS,ids,chans,baseline_win_ms,analysis_win_ms,scale_factor)
% Baseline-corrected analysis-window count per trial x channel, and raw baseline counts.
nT = numel(ids);
t0 = DS.trig_ms(ids);
corr_counts = zeros(nT,numel(chans));
base_counts = zeros(nT,numel(chans));
for c = 1:numel(chans)
    st = DS.sp_by_mappos{chans(c)};
    if isempty(st), continue; end
    b0 = count_less_than(st,t0+baseline_win_ms(1));
    b1 = count_less_than(st,t0+baseline_win_ms(2));
    a0 = count_less_than(st,t0+analysis_win_ms(1));
    a1 = count_less_than(st,t0+analysis_win_ms(2));
    corr_counts(:,c) = (a1-a0)-(b1-b0)*scale_factor;
    base_counts(:,c) = b1-b0;
end
end

function S = group_stats(Cond,rows,keep,nBoot)
% Pooled actual / predicted / ratio of one group of condition rows (both stim
% orders, and dataset-pairs if several), weighted by number of channels - the
% same pooling SpikeTRatio does. keep{r} = logical trial mask of row r.
% nBoot > 0: bootstrap SD (actual trials of every row resampled; every single
% resampled ONCE per repeat and shared by all rows that use it).
S = struct('actual',NaN,'pred',NaN,'ratio',NaN,'actual_sd',NaN,'pred_sd',NaN,'ratio_sd',NaN,'n_trials',0);
usable = [];
for j = 1:numel(rows)
    r = rows(j);
    if any(keep{r}) && ~isnan(Cond(r).pred), usable(end+1) = r; end %#ok<AGROW>
end
if isempty(usable), return; end

nU = numel(usable);
w  = arrayfun(@(r) numel(Cond(r).list),usable);
vc = cell(1,nU);
av = zeros(1,nU); pv = zeros(1,nU);
for j = 1:nU
    r = usable(j);
    vc{j} = Cond(r).metric(keep{r});
    av(j) = mean(vc{j});
    pv(j) = Cond(r).pred;
end
S.actual = sum(w.*av)/sum(w);
S.pred   = sum(w.*pv)/sum(w);
if S.pred > 0, S.ratio = S.actual/S.pred; end
S.n_trials = sum(cellfun(@numel,vc));

if nBoot > 0
    kA = arrayfun(@(r) Cond(r).kA,usable);
    kB = arrayfun(@(r) Cond(r).kB,usable);
    nk = zeros(max([kA kB]),1);
    for j = 1:nU
        nk(kA(j)) = numel(Cond(usable(j)).vA);
        nk(kB(j)) = numel(Cond(usable(j)).vB);
    end
    keys = unique([kA kB]);
    ba = nan(nBoot,1); bp = nan(nBoot,1); br = nan(nBoot,1);
    for b = 1:nBoot
        idx = cell(numel(nk),1);
        for k = keys
            idx{k} = randi(nk(k),nk(k),1);
        end
        a_b = zeros(1,nU); p_b = zeros(1,nU);
        for j = 1:nU
            v = vc{j};
            a_b(j) = mean(v(randi(numel(v),numel(v),1)));
            r = usable(j);
            p_b(j) = mean(Cond(r).vA(idx{kA(j)}))+mean(Cond(r).vB(idx{kB(j)}));
        end
        ba(b) = sum(w.*a_b)/sum(w);
        bp(b) = sum(w.*p_b)/sum(w);
        if bp(b) > 0, br(b) = ba(b)/bp(b); end
    end
    S.actual_sd = std(ba);
    S.pred_sd   = std(bp);
    S.ratio_sd  = std(br,'omitnan');
end
end

function [r_odd,r_even] = group_split(Cond,rows,keep)
% Ratio from the odd / even kept trials only (position among the kept trials of each row).
r_odd = NaN; r_even = NaN;
for half = 1:2
    num = 0; den = 0; pnum = 0;
    for j = 1:numel(rows)
        r = rows(j);
        if isnan(Cond(r).pred), continue; end
        idx = find(keep{r});
        if half == 1, idx = idx(1:2:end); else, idx = idx(2:2:end); end
        if isempty(idx), continue; end
        wj = numel(Cond(r).list);
        num  = num+wj*mean(Cond(r).metric(idx));
        pnum = pnum+wj*Cond(r).pred;
        den  = den+wj;
    end
    if den > 0 && pnum > 0
        val = (num/den)/(pnum/den);
        if half == 1, r_odd = val; else, r_even = val; end
    end
end
end

function draw_trial_panel(ax,Cond,Flag,rr,mm,col_good,col_new,col_old,col_pred, ...
    ms_kept,ms_flag,fs,label_all,show_med,Use_High,high,Use_Low,low,Use_MaxCh,maxch, ...
    Use_BaseAbs,baseAbs,Use_Lin,lin_lo,lin_hi)
% Trials of one condition (PTD) at one distance and amplitude; both stim orders
% on the same axes: order 1 = circles (numbers above), order 2 = diamonds (numbers below).

if isempty(rr)
    text(ax,0.5,0.5,'no trials','Units','normalized','HorizontalAlignment','center');
    box(ax,'off');
    return;
end

two = numel(rr) > 1;
nx = 0; pred = NaN;
leg_h = gobjects(0); leg_txt = {};
for k = 1:numel(rr)
    r = rr(k);
    ord = Cond(r).ord;
    if ord == 1, mk = 'o'; va = 'bottom'; off = -0.15*two; ls = '-';  fa = 'normal';
    else,        mk = 'd'; va = 'top';    off =  0.15*two; ls = '--'; fa = 'italic';
    end
    if ~two, va = 'bottom'; fa = 'normal'; end      % one order per panel: numbers above, upright
    x = Cond(r).rel+off;
    switch mm
        case 1, y = Cond(r).metric;
        case 2, y = Cond(r).maxch;
        case 4, y = Flag(r).ratio;
        otherwise, y = Cond(r).base;
    end
    ex = Cond(r).exist; nw = Flag(r).new; good = ~ex & ~nw;
    plot(ax,x(good),y(good),mk,'Color',col_good,'MarkerFaceColor',col_good,'MarkerSize',ms_kept);
    if any(ex), plot(ax,x(ex),y(ex),'x','Color',col_old,'MarkerSize',ms_flag,'LineWidth',1.5); end
    if any(nw), plot(ax,x(nw),y(nw),mk,'Color',col_new,'MarkerFaceColor',col_new,'MarkerSize',ms_flag); end
    if (mm == 1 || mm == 4) && Use_Lin
        ring = Flag(r).lin & ~ex;
        if any(ring)
            plot(ax,x(ring),y(ring),'o','Color',[0 0 0],'MarkerFaceColor','none','MarkerSize',ms_flag+7,'LineWidth',1.3);
        end
    end
    if label_all
        for t = 1:numel(x)
            if nw(t), fw = 'bold'; else, fw = 'normal'; end
            text(ax,x(t),y(t),sprintf('%d',Cond(r).rel(t)),'Color',[0 0 0],'FontSize',fs,'FontWeight',fw,'FontAngle',fa, ...
                'HorizontalAlignment','center','VerticalAlignment',va,'Clipping','on');
        end
    end
    if mm == 1 && show_med && ~isnan(Flag(r).med)
        yline(ax,Flag(r).med,ls,'Color',[0.55 0.55 0.55]);
    elseif mm == 4 && show_med && ~isnan(Flag(r).med) && Cond(r).pred > 0
        yline(ax,Flag(r).med/Cond(r).pred,ls,'Color',[0.55 0.55 0.55]);
    end
    nx = max(nx,max(Cond(r).rel));
    pred = Cond(r).pred;
    leg_h(end+1) = plot(ax,NaN,NaN,mk,'Color',[0 0 0],'MarkerFaceColor','none','MarkerSize',7); %#ok<AGROW>
    leg_txt{end+1} = sprintf('%s (Set %d)',format_order(Cond(r).stims),Cond(r).si);              %#ok<AGROW>
end

if mm == 1
    if ~isnan(pred)
        yline(ax,pred,'-','Color',col_pred,'LineWidth',1.4);
        if Use_Lin && pred > 0
            yline(ax,pred*lin_lo,':','Color',[0.6 0.1 0.6]);
            yline(ax,pred*lin_hi,':','Color',[0.6 0.1 0.6]);
        end
    end
    if Use_High, yline(ax,high,'--','Color',[0.8 0.4 0.0]); end
    if Use_Low,  yline(ax,low,'--','Color',[0.8 0.4 0.0]);  end
elseif mm == 4
    yline(ax,1,'-','Color',col_pred,'LineWidth',1.4);
    if Use_Lin
        yline(ax,lin_lo,':','Color',[0.6 0.1 0.6]);
        yline(ax,lin_hi,':','Color',[0.6 0.1 0.6]);
    end
elseif mm == 2
    if Use_MaxCh, yline(ax,maxch,'--','Color',[0.8 0.4 0.0]); end
else
    if Use_BaseAbs, yline(ax,baseAbs,'--','Color',[0.8 0.4 0.0]); end
end

set(ax,'XTick',1:nx,'XGrid','on','YGrid','off','GridAlpha',0.12,'GridColor',[0 0 0],'FontSize',7,'Layer','bottom');
xlim(ax,[0 nx+1]);
ylim(ax,pad_ylim(ax));
if ~isempty(leg_h) && two
    legend(ax,leg_h,leg_txt,'Location','best','Box','off','FontSize',7,'Interpreter','none');
end
box(ax,'off');
end


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


function DS = load_linearity_dataset(spike_file,dataset_folder,Electrode_Type,FS_fallback,apply_badtrial_exclusion,tag)
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
if isempty(tag)
    badtrials_path = fullfile(dataset_folder,sprintf('%s_BadTrials.mat',base_name));
else
    badtrials_path = fullfile(dataset_folder,sprintf('%s_%s_BadTrials.mat',base_name,tag));
end

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
    warning('No bad-trial file found (looked for %s) - no trials excluded for this dataset.',badtrials_path);
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