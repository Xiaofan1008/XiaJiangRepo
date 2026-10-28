%% ========================================================================
% CHECK THE OTHER (ChB) SINGLE-STIMULATION TRIALS + THE PREDICTION - v1
%
% PURPOSE
%   Run this AFTER Ch1 is cleaned (Clean_Ch1_Single_v1 + BadTrial_Ch1_Single_v1).
%   It pools every single-pulse electrode except the anchor (Ch1) of all
%   single datasets of ONE animal, flags suspicious trials with the same
%   rules for every dataset, and shows, for every pair Ch1 + ChB, what the
%   PREDICTION (single Ch1 + single ChB) looks like and how much it moves when
%   the flagged trials are taken out. SPIKE COUNTS ONLY.
%
% PER-TRIAL VALUES (same counting as SpikeTRatio, over Channel_Range)
%   Metric : mean baseline-corrected count per channel          (1 = RESPONSE)
%   MaxCh  : largest single-channel baseline-corrected count    (2 = MAX single chn)
%   Base   : mean raw baseline-window count per channel         (3 = BASELINE)
%
% RULES  (switchable; thresholds are PLACEHOLDERS - set them from table 1)
%   LOW, HIGH, RESP-MAD (two-sided k*MAD), MAXCH, BASE-ABS, BASE-MAD.
%   Group for the MAD rules: one electrode x amplitude within one dataset
%   (MAD_Scope = 'dataset') or pooled over datasets (MAD_Scope = 'pooled').
%
% PREDICTION (per dataset, as SpikeTRatio does it)
%   Ch1 and ChB are taken from the SAME single dataset:
%     predicted = mean(Ch1 trials) + mean(ChB trials)   (Metric, spikes/chn/trial)
%   "kept now"        : every trial except the ones already saved as bad
%   "w/o NEW flagged" : also without the trials flagged in this run
%   Error bar = sqrt(SE_Ch1^2 + SE_ChB^2). Ch1 is shown in the first column of
%   each pair figure so you can see it is the same clean Ch1 every time.
%
% EXISTING BAD TRIALS
%   Each dataset's saved <base_name>_BadTrials.mat is loaded; those trials are
%   grey crosses and are left out of the statistics (Exclude_ExistingBad_...).
%   Only NEW flags are listed.
%
% OUTPUT (nothing is saved)
%   1. distribution table per dataset x electrode x amplitude
%   2. flag counts + flagged trials grouped by dataset > electrode > amplitude
%      (relative trial number = what you type in the entry script)
%   3. prediction table per pair (now vs without NEW flagged)
%   4. a ready-to-fill block of {StimCh, Amp, [trials]} rows per dataset
%   5. one figure per pair: rows = amplitudes; columns = Ch1 trials | ChB trials
%      | prediction (stacked bars: Ch1 part + ChB part)
% ========================================================================

clear; clc;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ============================ USER SETTINGS ===========================

% All SINGLE datasets of ONE animal (every distance block).
SingleDatasets = struct('spike_file',{},'dataset_folder',{});

SingleDatasets(1).spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1/Xia_Linearity_S1_100_400um_Single1.sp_xia_filtered.mat';
SingleDatasets(1).dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_Single1';

SingleDatasets(2).spike_file     = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1/Xia_Linearity_S1_500_700um_Single1.sp_xia_filtered.mat';
SingleDatasets(2).dataset_folder = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_500_700um_Single1';

Anchor_Channel = 1;              % the cleaned anchor electrode (Ch1)
Channel_Range  = 1:64;
Electrode_Type = 3;              % 0 rigid | 1 flexible | 2 four-shank | 3 64chn+32chn hybrid
FS_fallback    = 30000;

baseline_win_ms = [-50 -5];
analysis_win_ms = [0 20];

Amp_Tol_uA     = 0.001;
Amps_To_Show   = [5];             % [] = all amplitudes; otherwise e.g. [5 10]
Electrodes_To_Show = [];         % [] = every ChB; otherwise e.g. [3 5 7] (stim channel numbers)

% ---- rules (spikes per channel per trial, PLACEHOLDERS) ----
Use_LowRule      = false;  lowThreshold     = 0;     % LOW is OFF by default: a weak response is not an artifact
Use_HighRule     = true;   highThreshold    = 5;
Use_RespMADRule  = true;   respMAD_k        = 4;
Use_MaxChRule    = true;   maxChThreshold   = 20;
Use_BaseAbsRule  = true;   baseAbsThreshold = 3;
Use_BaseMADRule  = true;   baseMAD_k        = 4;

MAD_Scope = 'dataset';           % 'dataset' or 'pooled'
Exclude_ExistingBad_FromStats = true;

% ---- figure appearance ----
Trial_Measure    = [1 3];        % what the trial columns show: 1 RESPONSE | 2 MAX single chn | 3 BASELINE
                                 % one or several, e.g. 1, [1 3] or [1 2 3]; each measure gets a Ch1 + ChB panel pair
Label_All_Trials = true;         % black trial number above every point (flagged = bold)
Fig_Width_Per_Col = 470;        % pixels per panel column (width grows with the number of measures)
Fig_Height_Per_Amp = 330;        % pixels per amplitude row
Marker_Size_Kept    = 6;
Marker_Size_Flagged = 8;
Label_FontSize      = 8;
Plot_PredictionCurve = true;     % one extra figure: the whole prediction curve, distance vs spike count

%% =========================== INITIAL CHECKS ===========================

if isempty(SingleDatasets), error('SingleDatasets is empty.'); end
validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(analysis_win_ms,'analysis_win_ms');
MAD_Scope = lower(MAD_Scope);
if ~ismember(MAD_Scope,{'dataset','pooled'}), error('MAD_Scope must be ''dataset'' or ''pooled''.'); end
Trial_Measure = unique(Trial_Measure(:).');
if isempty(Trial_Measure) || any(~ismember(Trial_Measure,1:3)), error('Trial_Measure may only contain 1, 2, 3.'); end

starting_folder = pwd;
cleanup_object = onCleanup(@() cd(starting_folder)); %#ok<NASGU>

nDS = numel(SingleDatasets);
scale_factor = diff(analysis_win_ms)/diff(baseline_win_ms);

fprintf('\n============================================================\n');
fprintf('CHECK ChB SINGLE TRIALS + PREDICTION (v1)\n');
fprintf('============================================================\n');
fprintf('Datasets: %d | anchor Ch%d | MAD scope: %s\n',nDS,Anchor_Channel,MAD_Scope);
fprintf('Rules: LOW(<%g) %s | HIGH(>%g) %s | RESP-MAD(k=%g) %s | MAXCH(>%g) %s | BASE-ABS(>%g) %s | BASE-MAD(k=%g) %s\n', ...
    lowThreshold,onoff(Use_LowRule),highThreshold,onoff(Use_HighRule), ...
    respMAD_k,onoff(Use_RespMADRule),maxChThreshold,onoff(Use_MaxChRule), ...
    baseAbsThreshold,onoff(Use_BaseAbsRule),baseMAD_k,onoff(Use_BaseMADRule));

%% ====================== LOAD DATASETS, PER-TRIAL VALUES =================

T_ds = []; T_set = []; T_stim = []; T_amp = []; T_rel = []; T_abs = [];
T_metric = []; T_maxch = []; T_base = []; T_exist = false(0,1);

DSinfo = struct('label',{},'folder',{},'base_name',{},'n_existing',{});

for d = 1:nDS

    fprintf('\n--- Loading dataset %d of %d ---\n',d,nDS);
    DS = load_linearity_dataset(SingleDatasets(d).spike_file,SingleDatasets(d).dataset_folder, ...
        Electrode_Type,FS_fallback,true);

    if DS.sim_stim ~= 1
        warning('Dataset %d (%s) has simultaneous_stim = %d (not single-pulse) - skipped.',d,DS.base_name,DS.sim_stim);
        continue;
    end

    chans = unique(double(Channel_Range(:).'));
    chans = chans(chans >= 1 & chans <= DS.nMapPos);

    DSinfo(d).label      = DS.base_name;
    DSinfo(d).folder     = DS.dataset_folder;
    DSinfo(d).base_name  = DS.base_name;
    DSinfo(d).n_existing = numel(DS.bad_trial_ids);

    for ss = 1:DS.nSets

        sc = DS.uniqueComb(ss,DS.uniqueComb(ss,:) > 0);
        if numel(sc) ~= 1, continue; end
        stim_ch = sc(1);

        for ai = 1:numel(DS.Amps)

            ids = find(DS.combClass == ss & DS.ampIdx == ai);
            if isempty(ids), continue; end
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

            T_ds     = [T_ds;     repmat(d,nT,1)];                       %#ok<AGROW>
            T_set    = [T_set;    repmat(ss,nT,1)];                      %#ok<AGROW>
            T_stim   = [T_stim;   repmat(stim_ch,nT,1)];                 %#ok<AGROW>
            T_amp    = [T_amp;    repmat(DS.Amps(ai),nT,1)];             %#ok<AGROW>
            T_rel    = [T_rel;    (1:nT).'];                             %#ok<AGROW>
            T_abs    = [T_abs;    ids(:)];                               %#ok<AGROW>
            T_metric = [T_metric; mean(corr_counts,2)];                  %#ok<AGROW>
            T_maxch  = [T_maxch;  max(corr_counts,[],2)];                %#ok<AGROW>
            T_base   = [T_base;   mean(base_counts,2)];                  %#ok<AGROW>
            T_exist  = [T_exist;  ismember(ids(:),DS.bad_trial_ids(:))]; %#ok<AGROW>
        end
    end

    fprintf('Dataset %d (%s): %d single-pulse trials | %d already saved as bad | %d channels counted\n', ...
        d,DS.base_name,sum(T_ds == d),numel(DS.bad_trial_ids),numel(chans));
    clear DS;
end

if isempty(T_ds), error('No single-pulse trials were found.'); end

T_amp = round(T_amp*1000)/1000;
all_amps = unique(T_amp(:).');

keep = true(size(T_ds));
if ~isempty(Amps_To_Show)
    keep = keep & arrayfun(@(a) any(abs(Amps_To_Show(:)-a) < Amp_Tol_uA),T_amp);
end
if ~isempty(Electrodes_To_Show)
    % always keep the anchor (needed for the prediction)
    keep = keep & (ismember(T_stim,Electrodes_To_Show) | T_stim == Anchor_Channel);
end
T_ds = T_ds(keep); T_set = T_set(keep); T_stim = T_stim(keep); T_amp = T_amp(keep);
T_rel = T_rel(keep); T_abs = T_abs(keep); T_metric = T_metric(keep);
T_maxch = T_maxch(keep); T_base = T_base(keep); T_exist = T_exist(keep);
if isempty(T_ds)
    error('Nothing left after the amplitude / electrode selection. Amplitudes available (uA): %s',num2str(all_amps));
end

amps = unique(T_amp(:).');
nTrialsTotal = numel(T_ds);
ds_present = unique(T_ds(:).');
if ~any(T_stim == Anchor_Channel)
    warning('No Ch%d (anchor) single trials found - the prediction cannot be computed.',Anchor_Channel);
end

%% ============== GROUP STATISTICS (median / MAD) ========================

use_stat = true(nTrialsTotal,1);
if Exclude_ExistingBad_FromStats, use_stat = ~T_exist; end

if strcmp(MAD_Scope,'dataset')
    [~,~,grp_id] = unique([T_ds,T_stim,T_amp],'rows');
else
    [~,~,grp_id] = unique([T_stim,T_amp],'rows');
end
nGrp = max(grp_id);

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
isRespMAD = Use_RespMADRule & (trial_resp_mad > 0) & (abs(T_metric-trial_resp_med) > respMAD_k*trial_resp_mad);
isMaxCh   = Use_MaxChRule   & (T_maxch > maxChThreshold);
isBaseAbs = Use_BaseAbsRule & (T_base > baseAbsThreshold);
isBaseMAD = Use_BaseMADRule & (trial_base_mad > 0) & (T_base > trial_base_med+baseMAD_k*trial_base_mad);

anyRule = isLow | isHigh | isRespMAD | isMaxCh | isBaseAbs | isBaseMAD;
isNew   = anyRule & ~T_exist;

%% ===================== 1. DISTRIBUTION TABLE ===========================

fprintf('\n------------------------------------------------------------\n');
fprintf('1. Distribution per dataset x electrode x amplitude (spikes per channel per trial)\n');
fprintf('   stats exclude already-saved bad trials: %s\n',onoff(Exclude_ExistingBad_FromStats));
fprintf('------------------------------------------------------------\n');
fprintf('%-3s %-5s %-9s %-8s %4s %5s | %7s %6s %7s %7s | %7s %7s | %7s %6s %7s\n', ...
    'DS','Stim','Dist(um)','Amp(uA)','N','Exist','MetMed','MetMAD','MetMin','MetMax','MaxMed','MaxMax','BasMed','BasMAD','BasMax');
for d = ds_present
    stims = unique(T_stim(T_ds == d)).';
    for e = stims
        dist_txt = pair_distance_text(Anchor_Channel,e,Electrode_Type);
        for a = amps
            sel = T_ds == d & T_stim == e & T_amp == a;
            if ~any(sel), continue; end
            s2 = sel & use_stat;
            if ~any(s2), s2 = sel; end
            m = T_metric(s2);
            fprintf('%-3d Ch%-3d %-9s %-8g %4d %5d | %7.3g %6.3g %7.3g %7.3g | %7.3g %7.3g | %7.3g %6.3g %7.3g\n', ...
                d,e,dist_txt,a,sum(sel),sum(T_exist(sel)), ...
                median(m),median(abs(m-median(m))),min(m),max(m), ...
                median(T_maxch(s2)),max(T_maxch(s2)), ...
                median(T_base(s2)),median(abs(T_base(s2)-median(T_base(s2)))),max(T_base(s2)));
        end
    end
end
for d = ds_present
    fprintf('  DS %d = %s\n',d,DSinfo(d).label);
end

%% ===================== 2. FLAGS ========================================

fprintf('\n------------------------------------------------------------\n');
fprintf('2. Flag counts by rule (all trials; NEW = not yet saved as bad)\n');
fprintf('------------------------------------------------------------\n');
fprintf('  LOW (<%g): %d | HIGH (>%g): %d | RESP-MAD (k=%g): %d | MAXCH (>%g): %d | BASE-ABS (>%g): %d | BASE-MAD (k=%g): %d\n', ...
    lowThreshold,sum(isLow),highThreshold,sum(isHigh),respMAD_k,sum(isRespMAD), ...
    maxChThreshold,sum(isMaxCh),baseAbsThreshold,sum(isBaseAbs),baseMAD_k,sum(isBaseMAD));
fprintf('  Flagged by at least one rule: %d of %d (%.1f%%) | NEW: %d | already saved bad: %d\n', ...
    sum(anyRule),nTrialsTotal,100*sum(anyRule)/nTrialsTotal,sum(isNew),sum(anyRule & T_exist));

fprintf('\nFlagged trials (NEW only), grouped by dataset > electrode > amplitude\n');
fprintf('   rel = relative trial number within the condition; abs = absolute trial id\n');
if ~any(isNew), fprintf('  (no new flagged trials)\n'); end

for d = ds_present
    if ~any(isNew & T_ds == d), continue; end
    fprintf('\n[DS %d] %s\n',d,DSinfo(d).label);
    for e = unique(T_stim(isNew & T_ds == d)).'
        fprintf('  Ch%d (distance from Ch%d: %s)\n',e,Anchor_Channel,pair_distance_text(Anchor_Channel,e,Electrode_Type));
        for a = amps
            sel = find(isNew & T_ds == d & T_stim == e & T_amp == a);
            if isempty(sel), continue; end
            fprintf('    Amp %g uA: %d new flagged of %d trials\n',a,numel(sel),sum(T_ds == d & T_stim == e & T_amp == a));
            for k = sel(:).'
                rt = {};
                if isLow(k),     rt{end+1} = 'LOW'; end      %#ok<SAGROW>
                if isHigh(k),    rt{end+1} = 'HIGH'; end     %#ok<SAGROW>
                if isRespMAD(k), rt{end+1} = 'RESP-MAD'; end %#ok<SAGROW>
                if isMaxCh(k),   rt{end+1} = 'MAXCH'; end    %#ok<SAGROW>
                if isBaseAbs(k), rt{end+1} = 'BASE-ABS'; end %#ok<SAGROW>
                if isBaseMAD(k), rt{end+1} = 'BASE-MAD'; end %#ok<SAGROW>
                fprintf('      rel %3d (abs %4d) | Metric %7.3g | MaxCh %7.3g | Base %6.3g | %s\n', ...
                    T_rel(k),T_abs(k),T_metric(k),T_maxch(k),T_base(k),strjoin(rt,', '));
            end
        end
    end
end

%% ===================== 3. PREDICTION TABLE =============================

fprintf('\n------------------------------------------------------------\n');
fprintf('3. PREDICTION = mean Ch%d + mean ChB (same dataset), spikes per channel per trial\n',Anchor_Channel);
fprintf('   now = all kept trials | clean = also without NEW flagged | err = sqrt(SE_A^2+SE_B^2)\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-3s %-5s %-9s %-8s | %4s %4s | %7s %7s | %8s %7s | %8s %7s | %8s\n', ...
    'DS','ChB','Dist(um)','Amp(uA)','nA','nB','meanA','meanB','pred now','err','pred cln','err','change');

Pred = struct('ds',{},'stim',{},'amp',{},'dist',{},'dist_um',{}, ...
    'A_now',{},'B_now',{},'A_cln',{},'B_cln',{},'err_now',{},'err_cln',{});

for d = ds_present
    stims = unique(T_stim(T_ds == d)).';
    stims = stims(stims ~= Anchor_Channel);
    for e = stims
        dist_txt = pair_distance_text(Anchor_Channel,e,Electrode_Type);
        for a = amps
            selA = T_ds == d & T_stim == Anchor_Channel & T_amp == a & ~T_exist;
            selB = T_ds == d & T_stim == e              & T_amp == a & ~T_exist;
            if ~any(selA) || ~any(selB), continue; end
            selAc = selA & ~isNew;  selBc = selB & ~isNew;

            [A_now,seA_now] = mean_se(T_metric(selA));
            [B_now,seB_now] = mean_se(T_metric(selB));
            [A_cln,seA_cln] = mean_se(T_metric(selAc));
            [B_cln,seB_cln] = mean_se(T_metric(selBc));

            P = struct('ds',d,'stim',e,'amp',a,'dist',dist_txt,'dist_um',pair_distance_value(Anchor_Channel,e,Electrode_Type), ...
                'A_now',A_now,'B_now',B_now,'A_cln',A_cln,'B_cln',B_cln, ...
                'err_now',sqrt(seA_now^2+seB_now^2),'err_cln',sqrt(seA_cln^2+seB_cln^2));
            Pred(end+1) = P; %#ok<SAGROW>

            pn = A_now+B_now; pc = A_cln+B_cln;
            if pn ~= 0 && ~isnan(pc), chg = sprintf('%+.0f%%',100*(pc-pn)/abs(pn)); else, chg = 'n/a'; end
            fprintf('%-3d Ch%-3d %-9s %-8g | %4d %4d | %7.3g %7.3g | %8.3g %7.2g | %8.3g %7.2g | %8s\n', ...
                d,e,dist_txt,a,sum(selA),sum(selB),A_now,B_now,pn,P.err_now,pc,P.err_cln,chg);
        end
    end
end

%% ===================== 4. BLOCK TO FILL IN =============================

fprintf('\n------------------------------------------------------------\n');
fprintf('4. NEW flagged trials as {StimCh, Amp, [relative trials]} rows, per dataset.\n');
fprintf('   These are only the rule-flagged trials - edit them before you use them.\n');
fprintf('------------------------------------------------------------\n');

for d = ds_present
    if ~any(isNew & T_ds == d), continue; end
    fprintf('\n%% [DS %d] %s\n',d,DSinfo(d).label);
    fprintf('%% matFilePath = ''%s'';\n', ...
        fullfile(DSinfo(d).folder,sprintf('%s_TrialSpikeCounts.mat',DSinfo(d).base_name)));
    fprintf('bad = {\n');
    for e = unique(T_stim(isNew & T_ds == d)).'
        for a = amps
            sel = find(isNew & T_ds == d & T_stim == e & T_amp == a);
            if isempty(sel), continue; end
            fprintf('    %d,  %g,  [%s];\n',e,a, ...
                strjoin(arrayfun(@(x) sprintf('%d',x),sort(T_rel(sel)).','UniformOutput',false),','));
        end
    end
    fprintf('};\n');
end

%% ===================== 5. FIGURES (one per pair) =======================

col_good = [0.20 0.45 0.80];
col_new  = [0.85 0.15 0.15];
col_old  = [0.55 0.55 0.55];
col_B    = [0.90 0.60 0.15];

row_title = {'Metric: mean baseline-corrected count per channel', ...
    'MaxCh: largest single-channel count','Base: mean baseline count per channel'};
panel_name = {'RESPONSE','MAX single chn','BASELINE'};

C = struct('col_good',col_good,'col_new',col_new,'col_old',col_old, ...
    'ms_kept',Marker_Size_Kept,'ms_flag',Marker_Size_Flagged,'fs',Label_FontSize, ...
    'label_all',Label_All_Trials, ...
    'Use_High',Use_HighRule,'high',highThreshold,'Use_Low',Use_LowRule,'low',lowThreshold, ...
    'Use_RespMAD',Use_RespMADRule,'respK',respMAD_k,'Use_MaxCh',Use_MaxChRule,'maxch',maxChThreshold, ...
    'Use_BaseAbs',Use_BaseAbsRule,'baseAbs',baseAbsThreshold,'Use_BaseMAD',Use_BaseMADRule,'baseK',baseMAD_k);

for d = ds_present

    stims = unique(T_stim(T_ds == d)).';
    stims = stims(stims ~= Anchor_Channel);

    for e = stims

        amps_here = [];
        for a = amps
            if any(T_ds == d & T_stim == e & T_amp == a), amps_here(end+1) = a; end %#ok<SAGROW>
        end
        if isempty(amps_here), continue; end
        nA = numel(amps_here);

        dist_txt = pair_distance_text(Anchor_Channel,e,Electrode_Type);
        nM = numel(Trial_Measure);
        nCols = 2*nM+1;                 % [Ch1, ChB] per measure, then the prediction
        figure('Color','w','Name',sprintf('DS %d | Ch%d + Ch%d | %s',d,Anchor_Channel,e,dist_txt), ...
            'Position',[20 40 min(2000,Fig_Width_Per_Col*nCols+80) min(1000,Fig_Height_Per_Amp*nA+90)]);

        for ia = 1:nA
            a = amps_here(ia);

            % ---- trial panels: for every measure, Ch1 then ChB ----
            for im = 1:nM
                m = Trial_Measure(im);
                for col = 1:2
                    if col == 1, stim_here = Anchor_Channel; else, stim_here = e; end
                    sel = T_ds == d & T_stim == stim_here & T_amp == a;
                    ax = subplot(nA,nCols,(ia-1)*nCols+2*(im-1)+col);
                    hold(ax,'on');
                    if any(sel)
                        plot_trials(ax,T_rel(sel),pick_measure(m,T_metric(sel),T_maxch(sel),T_base(sel)), ...
                            T_exist(sel),isNew(sel),C);
                        gi = grp_id(find(sel,1));
                        add_rule_lines(ax,m,resp_med(gi),resp_mad(gi),base_med(gi),base_mad(gi),C);
                    end
                    title(ax,sprintf('%s | Ch%d | DS %d (%s) | %g uA | n = %d', ...
                        panel_name{m},stim_here,d,block_label(DSinfo(d).label),a,sum(sel)), ...
                        'FontSize',9,'Interpreter','none');
                    if col == 1, ylabel(ax,row_title{m},'FontSize',8); end
                    if ia == nA, xlabel(ax,'Condition trial # (relative)'); end
                    if any(sel)
                        nx = max(T_rel(sel));
                        set(ax,'XTick',1:nx,'XGrid','on','YGrid','off','GridAlpha',0.12,'GridColor',[0 0 0], ...
                            'FontSize',7,'Layer','bottom');
                        xlim(ax,[0 nx+1]);
                        ylim(ax,pad_ylim(ax));
                    end
                    box(ax,'off');
                end
            end

            % ---- last column: prediction ----
            ax = subplot(nA,nCols,(ia-1)*nCols+nCols);
            hold(ax,'on');
            pidx = find([Pred.ds] == d & [Pred.stim] == e & abs([Pred.amp]-a) < Amp_Tol_uA,1);
            if ~isempty(pidx)
                P = Pred(pidx);
                parts = [P.A_now P.B_now; P.A_cln P.B_cln];
                parts(isnan(parts)) = 0;
                hb = bar(ax,1:2,parts,'stacked','BarWidth',0.55);
                hb(1).FaceColor = col_good; hb(2).FaceColor = col_B;
                tot = [P.A_now+P.B_now, P.A_cln+P.B_cln];
                errs = [P.err_now P.err_cln];
                errorbar(ax,1:2,tot,errs,'k.','LineWidth',1.3,'CapSize',10);
                for k = 1:2
                    text(ax,k,tot(k)+errs(k),sprintf('%.3g',tot(k)),'HorizontalAlignment','center', ...
                        'VerticalAlignment','bottom','FontSize',9,'FontWeight','bold');
                end
                set(ax,'XTick',1:2,'XTickLabel',{'kept now','w/o NEW flagged'},'FontSize',8);
                xlim(ax,[0.4 2.6]);
                ylim(ax,[min(0,min(tot-errs))*1.1, max(tot+errs)*1.2+eps]);
                if ia == 1
                    legend(ax,hb,{sprintf('Ch%d part',Anchor_Channel),sprintf('Ch%d part',e)},'Location','best','Box','off','FontSize',7);
                end
            else
                text(ax,0.5,0.5,'no prediction (missing Ch1 or ChB trials)','Units','normalized', ...
                    'HorizontalAlignment','center');
            end
            title(ax,sprintf('PREDICTION Ch%d + Ch%d | %g uA',Anchor_Channel,e,a),'FontSize',9,'Interpreter','none');
            if ia == 1, ylabel(ax,'spikes / chn / trial','FontSize',8); end
            box(ax,'off');
        end

        sgtitle(sprintf('DS %d (%s) | Ch%d + Ch%d | %s   (blue = kept, red = NEW flagged, grey x = already saved bad;  grey line = median, dotted = k*MAD, dashed = manual limits)', ...
            d,block_label(DSinfo(d).label),Anchor_Channel,e,dist_txt),'FontSize',9,'Interpreter','none');
    end
end

%% ===================== 6. PREDICTION CURVE (distance vs spike count) ===

if Plot_PredictionCurve && ~isempty(Pred)

    pd = [Pred.dist_um];
    if all(isnan(pd))
        warning('No numeric distances available - prediction curve not drawn.');
    else
        amps_curve = unique([Pred.amp]);
        nAc = numel(amps_curve);
        ds_cmap = lines(max(numel(ds_present),3));

        figure('Color','w','Name','Prediction curve: distance vs spike count', ...
            'Position',[60 60 min(1900,520*nAc+100) 520]);

        for ia = 1:nAc
            a = amps_curve(ia);
            ax = subplot(1,nAc,ia);
            hold(ax,'on');
            h_now = gobjects(0); h_cln = gobjects(0); lbl_now = {}; lbl_cln = {};

            for di = 1:numel(ds_present)
                d = ds_present(di);
                idx = find([Pred.ds] == d & abs([Pred.amp]-a) < Amp_Tol_uA & ~isnan([Pred.dist_um]));
                if isempty(idx), continue; end
                [~,ord] = sort([Pred(idx).dist_um]);
                idx = idx(ord);
                xd  = [Pred(idx).dist_um];
                yn  = [Pred(idx).A_now]+[Pred(idx).B_now];
                yc  = [Pred(idx).A_cln]+[Pred(idx).B_cln];
                en  = [Pred(idx).err_now];
                ec  = [Pred(idx).err_cln];

                h1 = errorbar(ax,xd,yn,en,'-o','Color',ds_cmap(di,:),'MarkerFaceColor',ds_cmap(di,:), ...
                    'LineWidth',1.5,'MarkerSize',7,'CapSize',6);
                h2 = errorbar(ax,xd,yc,ec,'--o','Color',ds_cmap(di,:),'MarkerFaceColor','w', ...
                    'LineWidth',1.2,'MarkerSize',7,'CapSize',6);
                h_now(end+1) = h1; lbl_now{end+1} = sprintf('DS %d (%s): kept now',d,block_label(DSinfo(d).label)); %#ok<SAGROW>
                h_cln(end+1) = h2; lbl_cln{end+1} = sprintf('DS %d: w/o NEW flagged',d); %#ok<SAGROW>
            end

            title(ax,sprintf('PREDICTION (Ch%d + ChB) | %g uA',Anchor_Channel,a),'FontSize',10);
            xlabel(ax,sprintf('Distance of ChB from Ch%d (um)',Anchor_Channel));
            if ia == 1, ylabel(ax,'Predicted spikes / chn / trial (Ch1 + ChB)'); end
            if ~isempty(h_now)
                legend(ax,[h_now h_cln],[lbl_now lbl_cln],'Location','best','Box','off','FontSize',7,'Interpreter','none');
            end
            grid(ax,'on'); ax.GridAlpha = 0.15;
            box(ax,'off');
        end
        sgtitle('Prediction curve: solid = kept now, dashed/open = also without NEW flagged trials, colour = dataset (distance block); error bar = sqrt(SE_Ch1^2 + SE_ChB^2)', ...
            'FontSize',9);
    end
end

fprintf('\nDone. Nothing was saved.\n');

%% =========================== LOCAL FUNCTIONS ==========================

function v = pair_distance_value(ch_a,ch_b,Electrode_Type)
v = NaN;
try
    [dist_um,~] = ChnPairDistance(ch_a,ch_b,Electrode_Type);
    v = double(dist_um);
catch
end
end

function txt = pair_distance_text(ch_a,ch_b,Electrode_Type)
txt = 'n/a';
try
    [dist_um,~] = ChnPairDistance(ch_a,ch_b,Electrode_Type);
    if ~isnan(dist_um), txt = sprintf('%g um',dist_um); end
catch
end
end

function [m,se] = mean_se(v)
v = v(~isnan(v));
if isempty(v)
    m = NaN; se = NaN;
else
    m = mean(v);
    se = std(v)/sqrt(max(1,numel(v)));
    if numel(v) < 2, se = 0; end
end
end

function y = pick_measure(m,metric,maxch,base)
switch m
    case 1, y = metric;
    case 2, y = maxch;
    otherwise, y = base;
end
end

function plot_trials(ax,x,y,ex,nw,C)
good = ~ex & ~nw;
plot(ax,x(good),y(good),'o','Color',C.col_good,'MarkerFaceColor',C.col_good,'MarkerSize',C.ms_kept);
if any(ex), plot(ax,x(ex),y(ex),'x','Color',C.col_old,'MarkerSize',C.ms_flag,'LineWidth',1.5); end
if any(nw), plot(ax,x(nw),y(nw),'o','Color',C.col_new,'MarkerFaceColor',C.col_new,'MarkerSize',C.ms_flag); end
if C.label_all
    for k = 1:numel(x)
        if nw(k), fw = 'bold'; else, fw = 'normal'; end
        text(ax,x(k),y(k),sprintf('%d',x(k)),'Color',[0 0 0],'FontSize',C.fs,'FontWeight',fw, ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','Clipping','on');
    end
end
end

function add_rule_lines(ax,m,rmed,rmad,bmed,bmad,C)
if m == 1
    if ~isnan(rmed)
        yline(ax,rmed,'-','Color',[0.4 0.4 0.4]);
        if C.Use_RespMAD && rmad > 0
            yline(ax,rmed+C.respK*rmad,':','Color',[0.3 0.3 0.3]);
            yline(ax,rmed-C.respK*rmad,':','Color',[0.3 0.3 0.3]);
        end
    end
    if C.Use_High, yline(ax,C.high,'--','Color',[0.8 0.4 0.0]); end
    if C.Use_Low,  yline(ax,C.low,'--','Color',[0.8 0.4 0.0]);  end
elseif m == 2
    if C.Use_MaxCh, yline(ax,C.maxch,'--','Color',[0.8 0.4 0.0]); end
else
    if ~isnan(bmed), yline(ax,bmed,'-','Color',[0.4 0.4 0.4]); end
    if C.Use_BaseMAD && bmad > 0, yline(ax,bmed+C.baseK*bmad,':','Color',[0.3 0.3 0.3]); end
    if C.Use_BaseAbs, yline(ax,C.baseAbs,'--','Color',[0.8 0.4 0.0]); end
end
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