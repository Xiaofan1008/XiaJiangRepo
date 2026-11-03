%% ============================================================
% NN_RecoveryCountCheck_v1.m
% Numeric check of NN_SpikeRecovery_v1: is exactly what the scheme says being added?
% READ-ONLY: only prints tables, saves nothing.
%
% TABLE 1 (all electrical target conditions, all channels):
%   For every target trial/channel the spikes inside the window [0, blank_end) of the
%   RECOVERED target trial must equal (count AND times) the spikes of its paired SOURCE
%   trial in the same relative window (recovered data). Prints how many pass, and the
%   mean spikes per trial/channel: original target (should be ~0, blanked), recovered
%   target, source.
% TABLE 2 (conditions in condList): spikes right after each pulse, [pulse, pulse+gap_ms],
%   mean per trial/channel, in the ORIGINAL data, the RECOVERED data and the source level.
%   This is the "2 to 3 ms after each pulse" stretch discussed: it shows how many spikes
%   those stretches get from the lower level compared with the rest of the train
%   ([pulse+gap_ms, next pulse) is also printed for comparison).
% ============================================================
clear all; clc;

%% ================= USER SETTINGS =================
analysisPath = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/';
data_folder  = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses';
times_file   = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses/Xia_Train_Freq4_100Hz_20Pulses.sp_xia_filtered.mat';
condList = [];      % conditions for TABLE 2 ([] = skip table 2)
gap_ms   = 3;       % length of the stretch right after each pulse (ms)
FS = 30000;

%% ================= LOAD =================
if ~isempty(analysisPath), addpath(genpath(analysisPath)); end
assert(isfolder(data_folder), 'Folder not found: %s', data_folder);
cd(data_folder);
rec_file = strrep(times_file, '.sp_xia_filtered.mat', '.sp_xia_TrainRecovery.mat');
assert(isfile(times_file) && isfile(rec_file), 'Need %s and %s.', times_file, rec_file);
O = load(times_file);  sp_org = O.sp_clipped;  if isfield(O,'fs'), FS = O.fs; end
R = load(rec_file);    sp_rec = R.sp_clipped;  be = R.blank_end_ms;  srcTrial = R.srcTrial(:);
trialCond = R.Recovery_params.trialCond(:);  condTable = R.Recovery_params.condTable;
nChn = numel(sp_org);  n_Trials = numel(trialCond);  nCond = numel(condTable);
clear O R

f = dir('*.trig.dat'); f = f(~startsWith({f.name}, '.'));
fid = fopen(f(1).name, 'r'); trig = fread(fid, [1 inf], 'double'); fclose(fid);
t0 = nan(n_Trials,1);  nUse = min(numel(trig), n_Trials);  t0(1:nUse) = trig(1:nUse)/FS*1000;

%% ================= TABLE 1: WINDOW = SOURCE? =================
nTot = zeros(nCond,1); nPass = zeros(nCond,1);
sumOrg = zeros(nCond,1); sumRec = zeros(nCond,1); sumSrc = zeros(nCond,1);
for ch = 1:nChn
    so = sp_org{ch}(:,1);  sr = sp_rec{ch}(:,1);
    for tr = find(isfinite(be{ch}) & isfinite(t0) & isfinite(srcTrial))'
        k = trialCond(tr);  s = srcTrial(tr);  w = be{ch}(tr);
        rec = sr(lb_(sr, t0(tr)) : lb_(sr, t0(tr)+w) - 1) - t0(tr);
        src = sr(lb_(sr, t0(s))  : lb_(sr, t0(s)+w)  - 1) - t0(s);
        org = so(lb_(so, t0(tr)) : lb_(so, t0(tr)+w) - 1);
        nTot(k) = nTot(k) + 1;
        if numel(rec) == numel(src) && (isempty(rec) || max(abs(rec - src)) < 1e-6), nPass(k) = nPass(k) + 1; end
        sumOrg(k) = sumOrg(k) + numel(org);  sumRec(k) = sumRec(k) + numel(rec);  sumSrc(k) = sumSrc(k) + numel(src);
    end
end
fprintf('TABLE 1: recovered window == source window?  (mean spikes per trial and channel inside [0, blank_end))\n');
fprintf(' cond  level  source  trial-ch  identical(%%)  original  recovered  source   label\n');
for k = find(nTot' > 0)
    fprintf(' %3d   %2d/%-3d  %4g   %7d    %6.1f       %6.2f    %6.2f   %6.2f   %s\n', k, condTable(k).level, condTable(k).nLevels, ...
        condTable(k).sourceCond, nTot(k), 100*nPass(k)/nTot(k), sumOrg(k)/nTot(k), sumRec(k)/nTot(k), sumSrc(k)/nTot(k), condTable(k).label);
end
if all(nPass(nTot>0) == nTot(nTot>0)), fprintf('=> every window equals its source window exactly.\n');
else, fprintf('=> SOME windows differ from their source: investigate (see identical(%%) column).\n'); end

%% ================= TABLE 2: STRETCH AFTER EACH PULSE =================
if ~isempty(condList)
    fe = dir('*_exp_datafile_*.mat'); fe = fe(~startsWith({fe.name}, '.'));
    E = load(fe(1).name, 'StimParams', 'simultaneous_stim');  nEl = E.simultaneous_stim;
    for k = condList(:)'
        t1 = find(trialCond == k, 1);
        if isempty(t1), fprintf('\ncond %d has no trials.\n', k); continue; end
        pt = [];
        for e = 1:nEl
            row = 1 + (t1-1)*nEl + e;
            if double(E.StimParams{row,16}) > 0
                n = double(E.StimParams{row,8});  per = double(E.StimParams{row,9});  p0 = double(E.StimParams{row,6});
                pt = [pt, (p0 + (0:n-1)*per)/1000]; %#ok<AGROW>
            end
        end
        pt = unique(pt);
        trK = find(trialCond == k & isfinite(t0));
        srcK = condTable(k).sourceCond;
        trS = []; if isfinite(srcK), trS = find(trialCond == srcK & isfinite(t0)); end
        fprintf('\nTABLE 2: cond %d (%s), mean spikes per trial and channel\n', k, condTable(k).label);
        fprintf(' pulse  time(ms)   [p, p+%g): original  recovered  source  |  [p+%g, next): original  recovered  source\n', gap_ms, gap_ms);
        for j = 1:numel(pt)
            a1 = pt(j);  a2 = pt(j) + gap_ms;
            if j < numel(pt), a3 = pt(j+1); else, a3 = pt(j) + 3*gap_ms; end
            v = [cnt_(sp_org, t0(trK), a1, a2), cnt_(sp_rec, t0(trK), a1, a2), cnt_(sp_rec, t0(trS), a1, a2), ...
                 cnt_(sp_org, t0(trK), a2, a3), cnt_(sp_rec, t0(trK), a2, a3), cnt_(sp_rec, t0(trS), a2, a3)];
            fprintf(' %4d  %7.1f           %7.3f   %7.3f   %7.3f  |            %7.3f   %7.3f   %7.3f\n', j, pt(j), v);
        end
    end
end

%% ================= HELPERS =================
function i = lb_(s, x)
lo = 1; hi = numel(s) + 1;
while lo < hi
    m = floor((lo + hi)/2);
    if s(m) < x, lo = m + 1; else, hi = m; end
end
i = lo;
end

function m = cnt_(sp, t0s, a, b)
% mean spikes per trial and channel in [t0+a, t0+b)
if isempty(t0s) || b <= a, m = NaN; return; end
tot = 0;
for ch = 1:numel(sp)
    s = sp{ch}(:,1);
    for i = 1:numel(t0s)
        tot = tot + lb_(s, t0s(i)+b) - lb_(s, t0s(i)+a);
    end
end
m = tot / (numel(t0s) * numel(sp));
end
