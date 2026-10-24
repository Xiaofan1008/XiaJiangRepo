%% ========================================================================
% PLOT DISTANCE vs SPIKE COUNT and DISTANCE vs LINEARITY - SPECIFIC AMPLITUDE(S) - v1
%
% PURPOSE
%   For the amplitude(s) in Amps_To_Plot, compute from the raw spikes
%     - Actual    : mean baseline-corrected spike count per channel per trial of the
%                   Sim / Seq paired stimulation
%     - Predicted : mean(single A) + mean(single B), same channels
%     - Ratio     : Actual / Predicted  (1 = linear),  LI = (A-P)/(A+P)  (0 = linear)
%   and plot them against electrode-pair distance. Same counting, channel lists
%   (responding channels, one union list per pair + amplitude), pooling and
%   bootstrap error bars as Clean_Paired_v2 / SpikeTRatio_v4.
%
% BAD TRIALS (the point of this script)
%   Every amplitude loads ITS OWN paired bad-trial file:
%       <base_name>_<tag>_BadTrials.mat      tag = '<amp>uA' (e.g. 10uA)
%   so the 10 uA result never uses the 5 uA trial marks. The single datasets use
%   the untagged <base_name>_BadTrials.mat (only affects Predicted).
%   If an amplitude's file is missing: Missing_BadFile = 'warn' continues with NO
%   exclusion for that amplitude and marks it "NOT CLEANED" in tables and figures;
%   'error' stops.
%
% FIGURES (one panel per amplitude unless noted; pooled over stim orders and
%   dataset-pairs; error bar = bootstrap SD over trials, singles resampled too)
%   1. Distance vs spike count : Actual (one line per PTD) and Predicted
%   2. Distance vs linearity   : ratio and/or LI, one line per PTD
%   3. Amplitude overlay       : linearity vs distance, one panel per PTD, one line
%                                per amplitude (only if several amplitudes)
%
% OUTPUT
%   Printed table + variable Results (struct array) in the workspace; optional CSV.
%   Nothing is written to the data folders.
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

Amp_Tol_uA = 0.001;
PTD_Tol_ms = 0.001;
Distance_Tol_um = 1;
Distances_To_Show = [];          % [] = all; otherwise e.g. [100 300 500]

% ---- amplitudes and their bad-trial files ----
Amps_To_Plot = [5];             % e.g. 10  or  [5 10]  or  [3 5 10]
BadTrial_Tags = 'auto';          % 'auto' -> '<amp>uA' for every amplitude; or a cell, one tag per amplitude,
                                 % e.g. {'3uA','5uA','10uA'} (same order as Amps_To_Plot)
Missing_BadFile = 'warn';        % 'warn' = continue uncleaned (flagged in output) | 'error' = stop

% ---- uncertainty ----
Boot_N    = 300;                 % bootstrap repeats for the error bars (0 = off)
Boot_Seed = 1;

% ---- plots ----
Plot_SpikeVsDistance   = true;   % figure 1
Plot_LinearityVsDistance = true; % figure 2
Linearity_Measure      = 'ratio';% 'ratio' (Actual/Predicted, 1 = linear) | 'LI' (0 = linear) | 'both'
Plot_AmpOverlay        = true;   % figure 3 (needs >1 amplitude)
Ratio_YLim             = [];     % [] = automatic; otherwise e.g. [0 2]
Save_CSV_Path          = '';     % e.g. '/Users/me/linearity_10uA.csv' ('' = do not save)

%% =========================== INITIAL CHECKS ===========================

if isempty(DatasetPairs), error('DatasetPairs is empty.'); end
validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');
if isempty(Amps_To_Plot) || ~isnumeric(Amps_To_Plot), error('Amps_To_Plot must be one or more amplitudes in uA.'); end
Amps_To_Plot = unique(double(Amps_To_Plot(:).'));
nAmpReq = numel(Amps_To_Plot);

if ischar(BadTrial_Tags) && strcmpi(BadTrial_Tags,'auto')
    amp_tags = arrayfun(@(a) sprintf('%guA',a),Amps_To_Plot,'UniformOutput',false);
elseif iscell(BadTrial_Tags) && numel(BadTrial_Tags) == nAmpReq
    amp_tags = BadTrial_Tags(:).';
else
    error('BadTrial_Tags must be ''auto'' or a cell with one tag per amplitude (%d).',nAmpReq);
end
if ~any(strcmpi(Missing_BadFile,{'warn','error'})), error('Missing_BadFile must be ''warn'' or ''error''.'); end
if ~any(strcmpi(Linearity_Measure,{'ratio','li','both'})), error('Linearity_Measure must be ''ratio'', ''LI'' or ''both''.'); end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>
rng(Boot_Seed);
scale_factor = diff(analysis_win_ms)/diff(baseline_win_ms);

fprintf('\n============================================================\n');
fprintf('DISTANCE vs SPIKE COUNT / LINEARITY, SPECIFIC AMPLITUDE(S) (v1)\n');
fprintf('============================================================\n');
fprintf('Dataset-pairs: %d | baseline [%g,%g) ms | analysis [%g,%g) ms\n', ...
    numel(DatasetPairs),baseline_win_ms,analysis_win_ms);
for ia = 1:nAmpReq
    fprintf('  %g uA -> paired bad-trial file tag ''%s''\n',Amps_To_Plot(ia),amp_tags{ia});
end

%% ================ LOAD DATASETS, BUILD ONE ROW PER CONDITION ===========

Cond = struct('dp',{},'si',{},'stims',{},'dist',{},'amp',{},'ptd',{},'ord',{}, ...
    'rel',{},'abs',{},'metric',{},'exist',{},'list',{}, ...
    'kA',{},'kB',{},'vA',{},'vB',{},'pred',{},'cleaned',{});

SingleStore = struct('corr',{},'bad',{});
single_key_map = containers.Map('KeyType','char','ValueType','double');
n_bad_total = zeros(1,nAmpReq);

for dp = 1:numel(DatasetPairs)

    entry = DatasetPairs(dp);
    fprintf('\n============================================================\n');
    fprintf('DATASET-PAIR %d/%d\n',dp,numel(DatasetPairs));
    fprintf('============================================================\n');

    if ~isfile(entry.simseq_spike_file),   error('Pair %d: paired spike file missing:\n%s',dp,entry.simseq_spike_file); end
    if ~isfolder(entry.simseq_dataset_folder), error('Pair %d: paired folder missing:\n%s',dp,entry.simseq_dataset_folder); end
    if ~isfile(entry.single_spike_file),   error('Pair %d: single spike file missing:\n%s',dp,entry.single_spike_file); end
    if ~isfolder(entry.single_dataset_folder), error('Pair %d: single folder missing:\n%s',dp,entry.single_dataset_folder); end

    fprintf('\n--- Loading paired dataset (spikes only; bad trials are read per amplitude below) ---\n');
    Paired = load_linearity_dataset(entry.simseq_spike_file,entry.simseq_dataset_folder, ...
        Electrode_Type,FS_fallback,false,'');
    fprintf('\n--- Loading single dataset (untagged bad-trial file) ---\n');
    Single = load_linearity_dataset(entry.single_spike_file,entry.single_dataset_folder, ...
        Electrode_Type,FS_fallback,true,'');

    % ---- this pair's bad-trial file for every requested amplitude ----
    bad_ids = cell(1,nAmpReq); has_file = false(1,nAmpReq);
    fprintf('\n--- Paired bad-trial files (one per amplitude) ---\n');
    for ia = 1:nAmpReq
        bf = fullfile(Paired.dataset_folder,sprintf('%s_%s_BadTrials.mat',Paired.base_name,amp_tags{ia}));
        if isfile(bf)
            BL = load(bf,'AllBadAbsoluteTrialIDs');
            if isfield(BL,'AllBadAbsoluteTrialIDs')
                bad_ids{ia} = double(BL.AllBadAbsoluteTrialIDs(:));
                has_file(ia) = true;
                n_bad_total(ia) = n_bad_total(ia)+numel(bad_ids{ia});
                fprintf('  %g uA: %s (%d trials excluded)\n',Amps_To_Plot(ia),bf,numel(bad_ids{ia}));
            else
                msg = sprintf('%g uA: file has no AllBadAbsoluteTrialIDs:\n%s',Amps_To_Plot(ia),bf);
                if strcmpi(Missing_BadFile,'error'), error('%s',msg); else, warning('%s',msg); end
            end
        else
            msg = sprintf('Pair %d, %g uA: bad-trial file not found (looked for %s)',dp,Amps_To_Plot(ia),bf);
            if strcmpi(Missing_BadFile,'error'), error('%s',msg); else, warning('%s - this amplitude is NOT CLEANED.',msg); end
        end
    end

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
        if isnan(dist_um), continue; end
        if ~isempty(Distances_To_Show) && ~any(abs(Distances_To_Show(:)-dist_um) <= Distance_Tol_um)
            continue;
        end
        if stims(1) < stims(2), ord = 1; else, ord = 2; end

        for ai = 1:numel(Paired.Amps)

            amp_value = round(Paired.Amps(ai)*1000)/1000;
            ia = find(abs(Amps_To_Plot-amp_value) < Amp_Tol_uA,1);
            if isempty(ia), continue; end

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

                [corr_p,~] = trial_counts(Paired,ids,chans,baseline_win_ms,analysis_win_ms,scale_factor);

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
                Cond(n).exist  = ismember(ids(:),bad_ids{ia});
                Cond(n).list   = list;
                Cond(n).kA = kk(1); Cond(n).kB = kk(2);
                Cond(n).vA = vv{1}; Cond(n).vB = vv{2};
                Cond(n).pred = pred;
                Cond(n).cleaned = has_file(ia);
            end
        end
    end

    fprintf('Pair %d: %d conditions so far | %d single trials saved as bad\n', ...
        dp,sum([Cond.dp] == dp),numel(Single.bad_trial_ids));
    clear Paired Single;
end

nCond = numel(Cond);
if nCond == 0, error('No paired conditions match the selection (check Amps_To_Plot / Distances_To_Show).'); end

C_amp  = [Cond.amp].';
C_ptd  = [Cond.ptd].';
C_dkey = [Cond.dist].';
amps_all = unique(C_amp).';
ptds_all = unique(C_ptd).';
keep = cell(nCond,1);
for r = 1:nCond, keep{r} = ~Cond(r).exist; end

%% ===================== ACTUAL / PREDICTED / RATIO ======================

[ukeys,~] = unique([C_dkey,C_amp,C_ptd],'rows');
nG = size(ukeys,1);
Results = struct('dist_um',{},'amp_uA',{},'ptd_ms',{},'n_trials',{},'n_excluded',{},'cleaned',{}, ...
    'predicted',{},'predicted_sd',{},'actual',{},'actual_sd',{},'ratio',{},'ratio_sd',{},'LI',{},'LI_sd',{});
for g = 1:nG
    rows = find(C_dkey == ukeys(g,1) & C_amp == ukeys(g,2) & C_ptd == ukeys(g,3));
    S = group_stats(Cond,rows,keep,Boot_N);
    R.dist_um = ukeys(g,1); R.amp_uA = ukeys(g,2); R.ptd_ms = ukeys(g,3);
    R.n_trials = S.n_trials;
    R.n_excluded = sum(arrayfun(@(r) sum(Cond(r).exist),rows));
    R.cleaned = all([Cond(rows).cleaned]);
    R.predicted = S.pred; R.predicted_sd = S.pred_sd;
    R.actual = S.actual;  R.actual_sd = S.actual_sd;
    R.ratio = S.ratio;    R.ratio_sd = S.ratio_sd;
    R.LI = (S.ratio-1)/(S.ratio+1);
    R.LI_sd = 2*S.ratio_sd/(S.ratio+1)^2;
    Results(g) = R;
end

fprintf('\n------------------------------------------------------------\n');
fprintf('RESULTS (pooled over stim orders and dataset-pairs; +-SD = bootstrap SD)\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-9s %-7s %-12s | %5s %5s | %8s %8s | %6s %6s | %6s | %s\n', ...
    'Dist(um)','Amp','PTD','nTr','nExcl','Pred','Actual','Ratio','+-SD','LI','bad file');
[~,oo] = sortrows([[Results.amp_uA].',[Results.ptd_ms].',[Results.dist_um].']);
prev_amp = NaN;
for g = oo(:).'
    R = Results(g);
    if ~isnan(prev_amp) && R.amp_uA ~= prev_amp, fprintf('\n'); end
    prev_amp = R.amp_uA;
    fprintf('%-9g %-7g %-12s | %5d %5d | %8.3g %8.3g | %6.2f %6.2f | %6.2f | %s\n', ...
        R.dist_um,R.amp_uA,short_ptd(R.ptd_ms),R.n_trials,R.n_excluded,R.predicted,R.actual, ...
        R.ratio,R.ratio_sd,R.LI,ternary(R.cleaned,'ok','NOT CLEANED'));
end

if any(~[Results.cleaned])
    bad_amps = unique([Results(~[Results.cleaned]).amp_uA]);
    fprintf('\nWARNING: no paired bad-trial file for %s uA - those results use ALL trials.\n',num2str(bad_amps));
end
for ia = 1:nAmpReq
    if ~any(abs(amps_all-Amps_To_Plot(ia)) < Amp_Tol_uA)
        fprintf('WARNING: %g uA was requested but no paired condition was found at that amplitude.\n',Amps_To_Plot(ia));
    end
end

if ~isempty(Save_CSV_Path)
    T = struct2table(Results);
    writetable(T,Save_CSV_Path);
    fprintf('\nResults table saved: %s\n',Save_CSV_Path);
end

%% ============================= FIGURES =================================

col_pred = [0.10 0.10 0.10];
col_ptd  = lines(max(numel(ptds_all),3));
nAc = numel(amps_all);
fig_w = max(900,min(1900,520*nAc+100));

amp_title = @(a) sprintf('%g uA%s',a,ternary(all([Results([Results.amp_uA] == a).cleaned]),'',' - NOT CLEANED (all trials)'));

% ---- 1. distance vs spike count ----
if Plot_SpikeVsDistance
    figure('Color','w','Name','Distance vs spike count','Position',[40 60 fig_w 520]);
    for ia = 1:nAc
        a = amps_all(ia);
        ax = subplot(1,nAc,ia); hold(ax,'on');
        hl = gobjects(0); ll = {};
        for ip = 1:numel(ptds_all)
            sel = find([Results.amp_uA] == a & abs([Results.ptd_ms]-ptds_all(ip)) < PTD_Tol_ms);
            if isempty(sel), continue; end
            [~,o] = sort([Results(sel).dist_um]); sel = sel(o);
            xd = [Results(sel).dist_um];
            h1 = errorbar(ax,xd,[Results(sel).actual],[Results(sel).actual_sd]*1.5,'-o','Color',col_ptd(ip,:), ...
                'MarkerFaceColor',col_ptd(ip,:),'LineWidth',1.6,'MarkerSize',7,'CapSize',6);
            hl(end+1) = h1; ll{end+1} = ['Actual ' ptd_label(ptds_all(ip))]; %#ok<SAGROW>
            if ip == 1
                h2 = errorbar(ax,xd,[Results(sel).predicted],[Results(sel).predicted_sd]*1.5,'-s','Color',col_pred, ...
                    'MarkerFaceColor',col_pred,'LineWidth',1.6,'MarkerSize',6,'CapSize',6);
                hl(end+1) = h2; ll{end+1} = 'Predicted'; %#ok<SAGROW>
            end
        end
        title(ax,amp_title(a),'FontSize',10,'Interpreter','none');
        xlabel(ax,'Distance (um)');
        if ia == 1, ylabel(ax,'spikes / chn / trial (baseline-corrected)'); end
        if ~isempty(hl), legend(ax,hl,ll,'Location','best','Box','off','FontSize',8); end
        grid(ax,'off'); ax.GridAlpha = 0.15; box(ax,'off');
    end
    sgtitle('Spike count vs distance (Actual of Sim / Seq vs Predicted = single A + single B)','FontSize',9);
end

% ---- 2. distance vs linearity ----
if Plot_LinearityVsDistance
    switch lower(Linearity_Measure)
        case 'ratio', meas_list = {'ratio'};
        case 'li',    meas_list = {'LI'};
        otherwise,    meas_list = {'ratio','LI'};
    end
    for im = 1:numel(meas_list)
        mname = meas_list{im};
        figure('Color','w','Name',sprintf('Distance vs linearity (%s)',mname),'Position',[60 80 fig_w 520]);
        for ia = 1:nAc
            a = amps_all(ia);
            ax = subplot(1,nAc,ia); hold(ax,'on');
            hl = gobjects(0); ll = {};
            for ip = 1:numel(ptds_all)
                sel = find([Results.amp_uA] == a & abs([Results.ptd_ms]-ptds_all(ip)) < PTD_Tol_ms);
                if isempty(sel), continue; end
                [~,o] = sort([Results(sel).dist_um]); sel = sel(o);
                if strcmp(mname,'LI'), yv = [Results(sel).LI]; ev = [Results(sel).LI_sd];
                else,                  yv = [Results(sel).ratio]; ev = [Results(sel).ratio_sd]; end
                h1 = errorbar(ax,[Results(sel).dist_um],yv,ev*1.5,'-o','Color',col_ptd(ip,:), ...
                    'MarkerFaceColor',col_ptd(ip,:),'LineWidth',1.6,'MarkerSize',7,'CapSize',6);
                hl(end+1) = h1; ll{end+1} = ptd_label(ptds_all(ip)); %#ok<SAGROW>
            end
            if strcmp(mname,'LI')
                yline(ax,0,':k','LineWidth',1.3); ylim(ax,[-1 1]);
                ylab = 'LI = (Actual - Predicted)/(Actual + Predicted)   (0 = linear)';
            else
                yline(ax,1,':k','LineWidth',1.3);
                if ~isempty(Ratio_YLim), ylim(ax,Ratio_YLim); end
                ylab = 'Ratio = Actual / Predicted   (1 = linear)';
            end
            title(ax,amp_title(a),'FontSize',10,'Interpreter','none');
            xlabel(ax,'Distance (um)');
            if ia == 1, ylabel(ax,ylab,'FontSize',8); end
            if ~isempty(hl), legend(ax,hl,ll,'Location','best','Box','off','FontSize',8); end
            grid(ax,'off'); ax.GridAlpha = 0.15; box(ax,'off');
        end
        sgtitle(sprintf('Linearity (%s) vs distance - clean trials of each amplitude',mname),'FontSize',9);
    end
end

% ---- 3. amplitude overlay: one panel per PTD, one line per amplitude ----
if Plot_AmpOverlay && nAc > 1
    nP = numel(ptds_all);
    cm_a = parula(nAc+1);
    meas_overlay = meas_list_or_default(Linearity_Measure);
    for im = 1:numel(meas_overlay)
        mname = meas_overlay{im};
        figure('Color','w','Name',sprintf('Amplitude overlay (%s)',mname),'Position',[80 100 max(900,min(1900,520*nP+100)) 520]);
        for ip = 1:nP
            ax = subplot(1,nP,ip); hold(ax,'on');
            hl = gobjects(0); ll = {};
            for ia = 1:nAc
                a = amps_all(ia);
                sel = find([Results.amp_uA] == a & abs([Results.ptd_ms]-ptds_all(ip)) < PTD_Tol_ms);
                if isempty(sel), continue; end
                [~,o] = sort([Results(sel).dist_um]); sel = sel(o);
                if strcmp(mname,'LI'), yv = [Results(sel).LI]; ev = [Results(sel).LI_sd];
                else,                  yv = [Results(sel).ratio]; ev = [Results(sel).ratio_sd]; end
                h1 = errorbar(ax,[Results(sel).dist_um],yv,ev*1.5,'-o','Color',cm_a(ia,:), ...
                    'MarkerFaceColor',cm_a(ia,:),'LineWidth',1.6,'MarkerSize',7,'CapSize',6);
                hl(end+1) = h1; ll{end+1} = sprintf('%g uA',a); %#ok<SAGROW>
            end
            if strcmp(mname,'LI')
                yline(ax,0,':k','LineWidth',1.3); ylim(ax,[-1 1]); ylab = 'LI   (0 = linear)';
            else
                yline(ax,1,':k','LineWidth',1.3);
                if ~isempty(Ratio_YLim), ylim(ax,Ratio_YLim); end
                ylab = 'Ratio = Actual / Predicted   (1 = linear)';
            end
            title(ax,ptd_label(ptds_all(ip)),'FontSize',10);
            xlabel(ax,'Distance (um)');
            if ip == 1, ylabel(ax,ylab,'FontSize',8); end
            if ~isempty(hl), legend(ax,hl,ll,'Location','best','Box','off','FontSize',8); end
            grid(ax,'on'); ax.GridAlpha = 0.15; box(ax,'off');
        end
        sgtitle(sprintf('Linearity (%s) vs distance - one line per amplitude (each with its own bad-trial file)',mname),'FontSize',9);
    end
end

fprintf('\nDone. Nothing was saved to the data folders.\n');

%% =========================== LOCAL FUNCTIONS ==========================

function lst = meas_list_or_default(m)
switch lower(m)
    case 'ratio', lst = {'ratio'};
    case 'li',    lst = {'LI'};
    otherwise,    lst = {'ratio','LI'};
end
end

function txt = short_ptd(p)
if abs(p) < 1e-6, txt = 'Sim'; else, txt = sprintf('Seq %g ms',p); end
end

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end

function txt = ptd_label(p)
if abs(p) < 1e-6, txt = 'Sim (PTD 0)'; else, txt = sprintf('Seq (PTD %g)',p); end
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