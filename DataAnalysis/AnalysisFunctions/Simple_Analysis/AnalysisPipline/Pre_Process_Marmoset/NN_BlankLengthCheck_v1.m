%% NN_BlankLengthCheck_v1  (script; READ-ONLY, prints a table and one figure, saves nothing)
% How long is the artifact blank after the LAST pulse, for each condition?
%
% The blanked file (amplifier_dn_sab.dat) has the blanked stretch replaced by a straight
% line (linear interpolation from the first event -1 ms to the recovery point after the last
% event). This script finds that straight stretch in a few trials per condition and channel
% and measures where it ends relative to the START of the last pulse:
%     blank after last pulse (ms) = end of straight line - start of last pulse
% (the normal pulse width is part of it). It also reports
%     - whether the straight line covers the WHOLE train (first event to last event),
%     - how often it hits the cap used by the blanking code (180 samples = 6 ms after the last event),
%     - how often no straight stretch was found (NaN).
% Run it AFTER run_preprocess_marmoset has finished on that folder.

%% ============================ SETTINGS ============================
folderPath = '/Volumes/MACData/Data/Data_Xia/CJ286/Xia_Train_Freq4_100Hz_20Pulses';   % <-- EDIT
channels   = [];        % amplifier channel numbers to measure; [] = all
plotCond   = [];        % condition IDs; [] = all electrical conditions
nTrialsPerCond = 8;     % random trials measured per condition
tolLSB     = 2;         % |second difference| <= this (int16 steps) counts as "straight line"
minRunSamp = 6;         % shortest straight run accepted (samples)
capSamp    = 180;       % cap used by the blanking code (samples after the last event)

%% ============================ LOAD ============================
origDir = pwd;
cleanupObj = onCleanup(@() cd(origDir)); %#ok<NASGU>
cd(folderPath);
[amplifier_channels, fp] = read_Intan_RHS2000_file; %#ok<ASGLU>
nChn = size(amplifier_channels,2);
FS   = fp.amplifier_sample_rate;
if isempty(channels), channels = 1:nChn; end

dnName = 'amplifier_dn_sab.dat';
assert(isfile(dnName), 'No %s - run run_preprocess_marmoset first.', dnName);
fDn = fopen(dnName,'r');
closeFiles = onCleanup(@() fclose('all')); %#ok<NASGU>

tf = dir('*.trig.dat'); tf = tf(~startsWith({tf.name}, '.'));
ft = fopen(tf(1).name,'r'); trig = fread(ft,[1 inf],'double'); fclose(ft);

ef = dir('*_exp_datafile_*.mat'); ef = ef(~startsWith({ef.name}, '.'));
S = load(ef(1).name,'TrialParams','StimParams');
condAll = cell2mat(S.TrialParams(2:end,2))';
nElect  = min(diff(find(diff(condAll)~=0)));
trialCond = condAll(1:nElect:end);
nTr = min(numel(trialCond), numel(trig));
trialCond = trialCond(1:nTr);  trig = trig(1:nTr);

condIDs = unique(trialCond);
if ~isempty(plotCond), condIDs = intersect(condIDs, plotCond); end

%% ============================ MEASURE ============================
rng(1);
res = struct('cond',{},'lastMs',{},'period',{},'blank',{},'coversTrain',{},'nNaN',{},'nCap',{},'n',{});
for c = condIDs
    t0 = find(trialCond == c, 1);
    rows = S.StimParams(1 + (t0-1)*nElect + (1:nElect), :);
    ev = []; per = NaN;
    for r = 1:size(rows,1)
        amp = rows{r,16};
        if isnumeric(amp) && ~isempty(amp) && amp > 0
            ptd = rows{r,6}/1000;  nP = rows{r,8};  per = rows{r,9}/1000;
            ev = [ev, ptd + (0:nP-1)*per]; %#ok<AGROW>
        end
    end
    if isempty(ev), continue; end                       % blank / visual condition
    firstMs = min(ev);  lastMs = max(ev);
    idx = find(trialCond == c & trig > 0);
    if isempty(idx), continue; end
    idx = idx(randperm(numel(idx), min(nTrialsPerCond, numel(idx))));

    preS  = round(2*FS/1000);                           % 2 ms before the trigger
    nSeg  = preS + round((lastMs + 20)*FS/1000) + 1;
    iFirst = preS + round(firstMs*FS/1000) + 1;
    iLast  = preS + round(lastMs*FS/1000) + 1;
    blank = nan(numel(channels), numel(idx));
    cover = false(numel(channels), numel(idx));
    for k = 1:numel(idx)
        s0 = round(trig(idx(k))) - preS;
        X = readSeg_(fDn, nChn, s0, nSeg);
        for ci = 1:numel(channels)
            x = X(channels(ci), :);
            if any(isnan(x)), continue; end
            lin = abs(diff(x,2)) <= tolLSB;             % lin(j): samples j..j+2 are straight
            [rs, re] = runs_(lin);                      % run j-range
            re = re + 2;                                % in sample units: start .. end sample
            keep = (re - rs + 1) >= minRunSamp;
            rs = rs(keep); re = re(keep);
            hit = find(rs <= iLast + 1 & re >= iLast - 1, 1, 'first');
            if isempty(hit), continue; end
            blank(ci,k) = (re(hit) - iLast) / FS * 1000;      % ms after the START of the last pulse
            cover(ci,k) = rs(hit) <= iFirst;
        end
    end
    v = blank(~isnan(blank));
    res(end+1) = struct('cond',c,'lastMs',lastMs,'period',per,'blank',v(:), ...
        'coversTrain', mean(cover(~isnan(blank))), 'nNaN', sum(isnan(blank(:))), ...
        'nCap', sum(v >= capSamp/FS*1000 - 0.05), 'n', numel(blank)); %#ok<SAGROW>
end

%% ============================ REPORT ============================
fprintf('\nBlank after the START of the last pulse (ms), measured on %d channels x %d trials per condition\n', numel(channels), nTrialsPerCond);
fprintf(' cond  lastPulse  period  | median   p5    p95   min   max | whole-train(%%)  at cap(%%)  not found(%%)\n');
med = nan(1,numel(res)); lo = med; hi = med;
for i = 1:numel(res)
    r = res(i);  v = r.blank;
    if isempty(v), fprintf(' %3d   %7.1f   %5.1f  |  no straight stretch found in any trial\n', r.cond, r.lastMs, r.period); continue; end
    med(i) = median(v);  lo(i) = prctile(v,5);  hi(i) = prctile(v,95);
    fprintf(' %3d   %7.1f   %5.1f  | %5.2f  %5.2f  %5.2f  %5.2f %5.2f |     %5.1f        %5.1f        %5.1f\n', ...
        r.cond, r.lastMs, r.period, med(i), lo(i), hi(i), min(v), max(v), 100*r.coversTrain, ...
        100*r.nCap/numel(v), 100*r.nNaN/r.n);
end
fprintf('\nRead: "median" is the blank length after the last pulse start. About 2 ms is the typical value you expect;\n');
fprintf('values near %.0f ms mean the cap was hit (blank may be too long). "whole-train" should be ~100%% for trains.\n', capSamp/FS*1000);

figure('Color','w','Name','blank length');
x = 1:numel(res);
errorbar(x, med, med-lo, hi-med, 'ko-', 'LineWidth', 1); hold on
yline(capSamp/FS*1000, 'r--', 'cap');
set(gca, 'XTick', x, 'XTickLabel', arrayfun(@(r) num2str(r.cond), res, 'UniformOutput', false));
xlabel('condition'); ylabel('blank after last pulse start (ms)  median, 5-95%'); box off

%% ============================ HELPERS ============================
function [rs, re] = runs_(m)
m = m(:)';
d = diff([0 m 0]);
rs = find(d == 1);
re = find(d == -1) - 1;
end

function X = readSeg_(fid, nChn, s0, n)
X = nan(nChn, n);
first = max(s0,1);
skip  = first - s0;
if first - 1 < 0, return; end
fseek(fid, (first-1)*nChn*2, 'bof');
d = fread(fid, [nChn, n-skip], 'int16');
if isempty(d), return; end
X(:, skip+(1:size(d,2))) = double(d);
end
