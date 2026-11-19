%% ========================================================================
% BAD-TRIAL REFERENCE CRITERIA (PRINT ONLY) - v3
%
% PURPOSE
%   Load a saved <base_name>_TrialSpikeCounts.mat file (from
%   TrialSpikeCounts_v2.m / Trial_Info) and flag trials by up to FOUR rules,
%   then print the flagged trials - grouped by condition (Set x Amp x PTD)
%   - using the RELATIVE trial ID (ConditionTrialIndex, i.e. the trial's
%   position within that specific condition).
%
%   RESPONSE-WINDOW RULES (original, unchanged)
%     LOW   : baseline-corrected analysis-window total  < lowThreshold
%     HIGH  : baseline-corrected analysis-window total  > highThreshold
%
%   BASELINE RULES (new) - use the pre-stimulus window only, so they do not
%   depend on the evoked response. Both are printed, side by side:
%     BASE-ABS : total baseline spike count > baselineAbsThreshold
%                (one fixed number for the whole dataset)
%     BASE-MAD : total baseline spike count > median + k*MAD, with the
%                median and MAD computed WITHIN EACH CONDITION
%                (Set x Amp x PTD). If a condition's MAD is 0 the rule
%                cannot be defined there and is skipped (noted in the table).
%
%   PRINT LAYOUT (v3): Print_GroupBy chooses how the flagged-trial list and
%   the per-condition summary table are organised - by 'set', 'distance',
%   'amplitude' or 'ptd', or a combination (outer to inner), or {} for the
%   old flat list. Each group header shows how many trials were flagged in
%   that group. The rules themselves are unchanged from v2.
%
% METRICS
%   Corrected analysis-window total: the analysis window baked into the
%   saved file (TotalBaselineCorrectedAnalysisWindowSpikes).
%   Baseline total: spikes in the saved baseline window
%   (TotalBaselineSpikes), RAW count (not scaled), summed over the same
%   channels as the corrected total.
% ========================================================================

clear; clc;

%% ============================ USER SETTINGS ===========================

% Full path to the saved <base_name>_TrialSpikeCounts.mat file
% matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S4_500_700um_Single1/Linearity_S4_500_700um_Single1_TrialSpikeCounts.mat';
matFilePath = '/Volumes/MACData/Data/Data_Xia/DX037/Xia_Linearity_S1_100_400um_SimSeq1/Linearity_S1_100_400um_SimSeq1_TrialSpikeCounts.mat';

% ---- response-window rules (flag if total < low OR > high) ----
Use_LowRule   = true;
Use_HighRule  = true;
lowThreshold  = 50;
highThreshold = 100;

% ---- baseline rules (new) ----
Use_BaselineAbsRule = true;
baselineAbsThreshold = 250;     % flag if total baseline spikes > this
                               % (PLACEHOLDER - set it after looking at the
                               % baseline distribution printed below)

Use_BaselineMADRule = true;
baselineMAD_k = 4;             % flag if baseline total > median + k*MAD
                               % (per condition; raw MAD, not scaled)

% Which channels to sum for the per-trial totals - double choice:
%   Leave empty [] to use the SAME channels the file was originally saved
%   with (count_channels), and reuse the already-computed TrialSummary
%   totals directly (fastest).
%   Or specify a custom channel list, e.g. InspectChannels = 1:32, to
%   recompute totals (corrected AND baseline) over a different channel
%   subset - this works because Channel_Results (saved per trial) covers
%   every map-position channel, not just the channels originally used.
InspectChannels = [1:64];

% ---- print layout ----
% How the flagged trials (and the summary table) are grouped. Allowed keys:
%   'set'        one group per ordered stim set (e.g. Ch1 -> Ch9)
%   'distance'   one group per electrode-pair distance (both stim orders together)
%   'amplitude'  one group per amplitude
%   'ptd'        one group per PTD
% Examples:
%   Print_GroupBy = {'set'};                  % group by set
%   Print_GroupBy = {'distance'};             % group by distance
%   Print_GroupBy = {'amplitude'};            % group by amplitude
%   Print_GroupBy = {'distance','amplitude'}; % distance, then amplitude inside it
%   Print_GroupBy = {};                       % flat list, no group headers (v2 style)
% (distance is n/a for single-electrode datasets - use 'set' or 'amplitude' there)
Print_GroupBy = {'amplitude'};

%% =========================== LOAD THE FILE ============================

if ~isfile(matFilePath)
    error('File does not exist:\n%s',matFilePath);
end

fprintf('\n============================================================\n');
fprintf('BAD-TRIAL REFERENCE CRITERIA (v3: response + baseline rules, grouped print)\n');
fprintf('============================================================\n');
fprintf('File: %s\n',matFilePath);

loadedData = load(matFilePath,'SpikeCounts','TrialSummary','count_channels', ...
    'selected_channels','baseline_win_ms','analysis_win_ms','Electrode_Type','sim_stim');

SpikeCounts       = loadedData.SpikeCounts;
TrialSummary      = loadedData.TrialSummary;
saved_count_chns  = loadedData.count_channels;
selected_channels = loadedData.selected_channels;
baseline_win_ms   = loadedData.baseline_win_ms;
analysis_win_ms   = loadedData.analysis_win_ms;

nTrials = height(TrialSummary);

fprintf('Baseline window used in this file: [%g,%g) ms\n',baseline_win_ms);
fprintf('Analysis window used in this file: [%g,%g) ms\n',analysis_win_ms);
fprintf('Trials in file: %d\n',nTrials);

fprintf('\nRules (ON/OFF):\n');
fprintf('  LOW      (corrected total < %g):                 %s\n',lowThreshold,onoff(Use_LowRule));
fprintf('  HIGH     (corrected total > %g):                 %s\n',highThreshold,onoff(Use_HighRule));
fprintf('  BASE-ABS (baseline total > %g):                  %s\n',baselineAbsThreshold,onoff(Use_BaselineAbsRule));
fprintf('  BASE-MAD (baseline total > median + %g*MAD, per condition): %s\n',baselineMAD_k,onoff(Use_BaselineMADRule));

%% ===================== DETERMINE CHANNELS TO USE ======================

useCustomChannels = ~isempty(InspectChannels);

if useCustomChannels

    inspectChannels = unique(double(InspectChannels(:).'),'stable');

    validChannels = isfinite(inspectChannels) & ...
        inspectChannels >= 1 & ...
        inspectChannels <= numel(selected_channels) & ...
        fix(inspectChannels) == inspectChannels;

    if any(~validChannels)
        warning('Ignoring invalid InspectChannels: %s', ...
            num2str(inspectChannels(~validChannels)));
        inspectChannels = inspectChannels(validChannels);
    end

    if isempty(inspectChannels)
        error('No valid InspectChannels remain.');
    end

    fprintf('Using CUSTOM channel subset (%d channels): %s\n', ...
        numel(inspectChannels),num2str(inspectChannels));
else
    inspectChannels = saved_count_chns;
    fprintf('Using the ORIGINAL saved count_channels (%d channels): %s\n', ...
        numel(inspectChannels),num2str(inspectChannels));
end

%% ===================== COMPUTE PER-TRIAL METRICS ======================

metricValue   = nan(nTrials,1);   % baseline-corrected analysis-window total
baselineValue = nan(nTrials,1);   % raw baseline-window total

hasBaselineVar = ismember('TotalBaselineSpikes',TrialSummary.Properties.VariableNames);
needLoop = useCustomChannels || ~hasBaselineVar;

if ~needLoop
    % Fast path: reuse the already-computed trial totals directly.
    metricValue   = TrialSummary.TotalBaselineCorrectedAnalysisWindowSpikes;
    baselineValue = TrialSummary.TotalBaselineSpikes;
else
    if ~useCustomChannels
        metricValue = TrialSummary.TotalBaselineCorrectedAnalysisWindowSpikes;
        fprintf('TotalBaselineSpikes not in TrialSummary - recomputing baseline totals from Channel_Results.\n');
    end

    for trial_id = 1:nTrials

        si  = TrialSummary.SetIndex(trial_id);
        ai  = TrialSummary.AmplitudeIndex(trial_id);
        pi  = TrialSummary.PTDIndex(trial_id);
        cti = TrialSummary.ConditionTrialIndex(trial_id);

        T = SpikeCounts.set(si).amp(ai).ptd(pi).trial(cti);
        CR = T.Channel_Results;

        [foundChn,chnPos] = ismember(inspectChannels,CR.Channel_Index);

        if any(~foundChn)
            warning(['Trial %d: some InspectChannels not found in ' ...
                'Channel_Results, skipping those for this trial.'],trial_id);
        end

        chnPos = chnPos(foundChn);
        if useCustomChannels
            metricValue(trial_id) = sum(CR.BaselineCorrected_AnalysisWindow_Spikes(chnPos));
        end
        baselineValue(trial_id) = sum(CR.Baseline_Spikes(chnPos));
    end
end

%% ================ BASELINE DISTRIBUTION (to help set thresholds) =======

fprintf('\n------------------------------------------------------------\n');
fprintf('Baseline total (raw spikes in [%g,%g) ms), all trials in file:\n',baseline_win_ms);
fprintf('------------------------------------------------------------\n');
sortedBase = sort(baselineValue(~isnan(baselineValue)));
pc = sortedBase(max(1,ceil([0.50 0.90 0.95 0.99]*numel(sortedBase))));
fprintf('  min %g | median %g | 90th %g | 95th %g | 99th %g | max %g\n', ...
    sortedBase(1),pc(1),pc(2),pc(3),pc(4),sortedBase(end));
fprintf('  (use this to choose baselineAbsThreshold)\n');

%% ===================== FLAG TRIALS BY EACH RULE ========================

isLow     = Use_LowRule  & (metricValue < lowThreshold);
isHigh    = Use_HighRule & (metricValue > highThreshold);
isBaseAbs = Use_BaselineAbsRule & (baselineValue > baselineAbsThreshold);

% BASE-MAD: median and MAD computed within each condition
[uCond,~,condId] = unique([TrialSummary.SetIndex,TrialSummary.AmplitudeIndex, ...
    TrialSummary.PTDIndex],'rows');
nCond = size(uCond,1);

condMedian = nan(nCond,1);
condMAD    = nan(nCond,1);
condThr    = nan(nCond,1);
isBaseMAD  = false(nTrials,1);

for c = 1:nCond
    idx = find(condId == c);
    b = baselineValue(idx);
    condMedian(c) = median(b);
    condMAD(c)    = median(abs(b-condMedian(c)));

    if condMAD(c) > 0
        condThr(c) = condMedian(c)+baselineMAD_k*condMAD(c);
        if Use_BaselineMADRule
            isBaseMAD(idx) = b > condThr(c);
        end
    end
end

anyFlag = isLow | isHigh | isBaseAbs | isBaseMAD;

fprintf('\n------------------------------------------------------------\n');
fprintf('Total trials flagged, by rule:\n');
fprintf('------------------------------------------------------------\n');
fprintf('  LOW      (<%g):                     %d\n',lowThreshold,sum(isLow));
fprintf('  HIGH     (>%g):                     %d\n',highThreshold,sum(isHigh));
fprintf('  BASE-ABS (>%g):                     %d\n',baselineAbsThreshold,sum(isBaseAbs));
fprintf('  BASE-MAD (>median+%g*MAD):          %d\n',baselineMAD_k,sum(isBaseMAD));
fprintf('  Flagged by at least one rule:       %d of %d (%.1f%%)\n', ...
    sum(anyFlag),nTrials,100*sum(anyFlag)/nTrials);
fprintf('  Flagged by BASE-ABS only:           %d\n',sum(isBaseAbs & ~isBaseMAD & ~isLow & ~isHigh));
fprintf('  Flagged by BASE-MAD only:           %d\n',sum(isBaseMAD & ~isBaseAbs & ~isLow & ~isHigh));
fprintf('  Flagged by BOTH baseline rules:     %d\n',sum(isBaseAbs & isBaseMAD));

%% ============== CONDITION LABELS AND PRINT ORDER (Print_GroupBy) ==============

if ischar(Print_GroupBy) || isstring(Print_GroupBy)
    Print_GroupBy = cellstr(Print_GroupBy);
end
Print_GroupBy = lower(Print_GroupBy(:).');
allKeys = {'set','distance','amplitude','ptd'};
badKeys = ~ismember(Print_GroupBy,allKeys);
if any(badKeys)
    error('Print_GroupBy has unknown key(s): %s. Allowed: set, distance, amplitude, ptd.', ...
        strjoin(Print_GroupBy(badKeys),', '));
end
Print_GroupBy = reshape(unique(Print_GroupBy,'stable'),1,[]);
nKeys = numel(Print_GroupBy);

condSetIdx  = uCond(:,1);
condAmp     = nan(nCond,1);
condPTD     = nan(nCond,1);
condDist    = nan(nCond,1);
condSetName = cell(nCond,1);
for c = 1:nCond
    si = uCond(c,1); ai = uCond(c,2); pi = uCond(c,3);
    condAmp(c)     = SpikeCounts.set(si).amp(ai).amp_value;
    condPTD(c)     = SpikeCounts.set(si).amp(ai).ptd(pi).PTD_ms;
    condDist(c)    = SpikeCounts.set(si).distance_um;
    condSetName{c} = SpikeCounts.set(si).set_name;
end

LBL = struct('setIdx',condSetIdx,'setName',{condSetName},'dist',condDist, ...
    'amp',condAmp,'ptd',condPTD);

keyVal = struct('set',condSetIdx,'distance',condDist,'amplitude',condAmp,'ptd',condPTD);
keyVal.distance(isnan(keyVal.distance)) = Inf;     % unknown distance (singles) sorts last

restKeys  = setdiff(allKeys,Print_GroupBy,'stable');   % default order: set, distance, amplitude, ptd
orderKeys = [Print_GroupBy,restKeys];
sortMat = zeros(nCond,numel(orderKeys));
for k = 1:numel(orderKeys)
    sortMat(:,k) = keyVal.(orderKeys{k});
end
[~,ordC] = sortrows(sortMat);

if nKeys == 0
    fprintf('\nPrint layout: flat list (no grouping).\n');
else
    fprintf('\nPrint layout: grouped by %s.\n',strjoin(Print_GroupBy,' > '));
end

%% ================ PER-CONDITION BASELINE / FLAG SUMMARY ================
% Check here that exclusion rates are balanced across conditions.
% Rows follow the same order as Print_GroupBy; a blank line separates the
% outermost groups.

fprintf('\n------------------------------------------------------------\n');
fprintf('Per-condition summary (baseline stats + number flagged per rule):\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-4s %-14s %-8s %-8s %-8s %4s | %8s %6s %8s | %4s %4s %5s %5s\n', ...
    'Set','Order','Dist(um)','Amp(uA)','PTD(ms)','N','BaseMed','MAD','MADthr','LOW','HIGH','B-ABS','B-MAD');

prevOuter = [];
for ii = 1:nCond
    c = ordC(ii);
    si = uCond(c,1);
    idx = condId == c;

    if nKeys >= 1
        curOuter = keyVal.(Print_GroupBy{1})(c);
        if ~isempty(prevOuter) && curOuter ~= prevOuter
            fprintf('\n');
        end
        prevOuter = curOuter;
    end

    thr_text = 'n/a(MAD=0)';
    if ~isnan(condThr(c)), thr_text = sprintf('%.1f',condThr(c)); end

    dist_str = 'n/a';
    if ~isnan(condDist(c)), dist_str = sprintf('%g',condDist(c)); end

    fprintf('%-4d %-14s %-8s %-8g %-8g %4d | %8.1f %6.1f %8s | %4d %4d %5d %5d\n', ...
        si,condSetName{c},dist_str,condAmp(c),condPTD(c), ...
        sum(idx),condMedian(c),condMAD(c),thr_text, ...
        sum(isLow & idx),sum(isHigh & idx),sum(isBaseAbs & idx),sum(isBaseMAD & idx));
end

if Use_BaselineMADRule && any(isnan(condThr))
    fprintf('\nNote: BASE-MAD could not be defined (MAD = 0) in %d condition(s) - no trials flagged by it there.\n', ...
        sum(isnan(condThr)));
end

%% ===================== PRINT FLAGGED TRIALS (GROUPED) ==================

fprintf('\n------------------------------------------------------------\n');
fprintf('Flagged trials (relative trial ID = position within that\n');
fprintf('condition, e.g. 1-30). Group headers show trials flagged in the group.\n');
fprintf('------------------------------------------------------------\n');

prevKey = [];
for ii = 1:nCond

    c = ordC(ii);
    si = uCond(c,1); ai = uCond(c,2); pi = uCond(c,3);

    % ---- group headers (only the levels that changed) ----
    curKey = zeros(1,nKeys);
    for k = 1:nKeys
        curKey(k) = keyVal.(Print_GroupBy{k})(c);
    end

    if nKeys > 0
        if isempty(prevKey)
            chg = 1;
        else
            chg = find(curKey ~= prevKey,1);
        end
        if ~isempty(chg)
            for L = chg:nKeys
                groupMask = true(nCond,1);
                for k = 1:L
                    groupMask = groupMask & (keyVal.(Print_GroupBy{k}) == curKey(k));
                end
                trialMask = ismember(condId,find(groupMask));
                nGrp  = sum(trialMask);
                nFlag = sum(anyFlag & trialMask);

                pad = repmat(' ',1,2*(L-1));
                labels = cell(1,L);
                for k = 1:L
                    labels{k} = key_label(Print_GroupBy{k},c,LBL);
                end
                fprintf('\n%s=== %s ===  %d trials | flagged %d (%.1f%%) | LOW %d  HIGH %d  B-ABS %d  B-MAD %d\n', ...
                    pad,strjoin(labels,' | '),nGrp,nFlag,100*nFlag/max(1,nGrp), ...
                    sum(isLow & trialMask),sum(isHigh & trialMask), ...
                    sum(isBaseAbs & trialMask),sum(isBaseMAD & trialMask));
                if L == nKeys && nFlag == 0
                    fprintf('%s  (no flagged trials in this group)\n',pad);
                end
            end
        end
    end
    prevKey = curKey;

    % ---- this condition's flagged trials ----
    conditionMask = ...
        TrialSummary.SetIndex == si & ...
        TrialSummary.AmplitudeIndex == ai & ...
        TrialSummary.PTDIndex == pi;

    if ~any(conditionMask)
        continue;
    end

    ids_low  = TrialSummary.ConditionTrialIndex(conditionMask & isLow);
    ids_high = TrialSummary.ConditionTrialIndex(conditionMask & isHigh);
    ids_babs = TrialSummary.ConditionTrialIndex(conditionMask & isBaseAbs);
    ids_bmad = TrialSummary.ConditionTrialIndex(conditionMask & isBaseMAD);

    if isempty(ids_low) && isempty(ids_high) && isempty(ids_babs) && isempty(ids_bmad)
        continue;
    end

    % label: every key that is NOT already shown in a group header
    lineKeys = restKeys;
    lineLabels = cell(1,numel(lineKeys));
    for k = 1:numel(lineKeys)
        lineLabels{k} = key_label(lineKeys{k},c,LBL);
    end
    pad = repmat(' ',1,2*nKeys);

    if isempty(lineLabels)
        fprintf('%s(condition):\n',pad);
    else
        fprintf('%s%s:\n',pad,strjoin(lineLabels,' | '));
    end

    if ~isempty(ids_low)
        fprintf('%s  LOW      (<%g):                trials %s\n', ...
            pad,lowThreshold,num2str(sort(ids_low(:).')));
    end

    if ~isempty(ids_high)
        fprintf('%s  HIGH     (>%g):              trials %s\n', ...
            pad,highThreshold,num2str(sort(ids_high(:).')));
    end

    if ~isempty(ids_babs)
        fprintf('%s  BASE-ABS (>%g spikes):         trials %s\n', ...
            pad,baselineAbsThreshold,num2str(sort(ids_babs(:).')));
    end

    if ~isempty(ids_bmad)
        fprintf('%s  BASE-MAD (>%.1f = med %.1f + %g*MAD %.1f): trials %s\n', ...
            pad,condThr(c),condMedian(c),baselineMAD_k,condMAD(c), ...
            num2str(sort(ids_bmad(:).')));
    end
end

fprintf('\n============================================================\n');
fprintf('INSPECTION COMPLETE (no files modified or saved)\n');
fprintf('============================================================\n');

%% =========================== LOCAL FUNCTIONS ==========================

function s = key_label(key,c,L)
% Text label of one grouping key for condition c.
switch key
    case 'set'
        s = sprintf('Set %d (%s)',L.setIdx(c),L.setName{c});
    case 'distance'
        if isnan(L.dist(c))
            s = 'Dist n/a';
        else
            s = sprintf('Dist %g um',L.dist(c));
        end
    case 'amplitude'
        s = sprintf('%g uA',L.amp(c));
    case 'ptd'
        s = sprintf('PTD %g ms',L.ptd(c));
    otherwise
        s = key;
end
end

function s = onoff(flag)
if flag, s = 'ON'; else, s = 'OFF'; end
end