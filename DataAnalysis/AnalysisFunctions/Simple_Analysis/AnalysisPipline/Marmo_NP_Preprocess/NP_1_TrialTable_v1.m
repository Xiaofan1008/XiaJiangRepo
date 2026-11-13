%% ============================================================
% NP_1_TrialTable_v1.m   (step 1 of the NP pipeline)
% NP probe (Neuropixels, Open Ephys) | marmoset pulse-TRAIN sessions | TRIAL TABLE | v1
%
% PURPOSE
%   Builds ONE file per session that every later NP script reads, so the trial bookkeeping is done once:
%   which trials exist, their condition, the trigger time on the probe clock, the stimulation pulse times,
%   the photodiode (visual) timing, and where every channel sits on the probe.
%   The recording is only READ. The only file written is <session>.NP_TrialTable.mat
%   (plus a dated backup of an older one when overwrite = true).
%
% INPUTS (found automatically from the session folder)
%   recording1/structure.oebin, continuous/ (AP, NI), events/ (trial messages, 1 Hz sync pulses)
%   recording1/<stim>.conds.mat    stim-PC design (cond, info, args, chnList). The neurostim .mat log is NOT needed.
%   Record Node xxx/settings.xml   probe channel positions
%
% STEPS
%   1. Streams from structure.oebin, channel positions from settings.xml.
%   2. Clock sync: the 1 Hz TTL pulses seen by the NI board and by the probe -> linear fit NI time -> probe time.
%   3. Triggers: rising edges of NI AI3 (100 ms pulse) -> probe sample.
%   4. Trial messages 'TR:n' and 'TRIALID c-n': each trigger is paired with the latest TR message before it.
%      Messages that never got a trigger (session stopped by hand) are dropped and listed.
%   5. Design from .conds.mat: pulse times of every condition (electrode A = chnList{1}, B = chnList{2}).
%   6. Photodiode (NI AI0): frame pulses after each trigger -> visual onset / offset per trial.
%   7. Artifact check (a sanity check, deliberately NOT strict): in a few trials per condition the stimulation
%      pulses are found in the AP signal. A condition is OK when the number of pulses is right and every pulse is
%      within artTolMs (1 ms) of the design. The delay trigger -> FIRST pulse is measured and added to the design times.
%      (The stimulator stretches the pulse spacing by ~0.3 % - 10.033 ms instead of 10 ms - so the last pulse of a
%       15-pulse train sits ~0.5 ms later than the design. That is accepted: blanking ends adaptively after the
%       last pulse, and later analyses use bins of 2 ms or more.)
%   8. Prints a summary, saves TT, draws check figures.
%
% OUTPUT  <outDir>/<short name>.NP_TrialTable.mat  with the struct TT
%         short name = the stim file name up to the first '.', e.g. CJ287_15Pulses3_100Hz.NP_TrialTable.mat
%         (or outName if you type one). The struct TT:
%   TT.info     session name, folders, files, sampling rates, channel count, uV per bit, stim names, settings, version
%   TT.sync     NI -> probe clock fit (P: probe_s = polyval(P, ni_s)), pulses used, residual, check against Open Ephys timestamps
%   TT.chan     ch (1..384), x_um, y_um (distance from the probe TIP), bank, orderTipToTop (channel numbers sorted by y),
%               preset, filterCut.  NOTE: channel number is NOT depth order for the 'Single column' preset.
%   TT.cond     one entry per condition of the conds file: id, label, kind ('sim','AB','BA','blank','visual'), level,
%               nLevels, isFinal, amp [A B], ptd_us [A B], nsp [A B], designMs (merged pulse times, ms after the
%               first possible pulse), designMsA / designMsB, nTrials (kept), contrast, phase
%   TT.pattern  one entry per train pattern (sim / AB / BA): kind, conds in level order, levels, finalCond
%   TT.trial    one entry per KEPT trial (columns):
%       trialNum   trial number from the stim PC (TR:n)
%       cond       condition number (index into TT.cond)
%       rep        repetition number of that condition (1, 2, ... in recording order)
%       trigNI     trigger sample in the NI file (1-based)
%       trigAP     trigger sample in the AP file (1-based)      <- use this to cut AP windows
%       trigAP_s   trigger time in s from the start of the AP file
%       nextTrigAP AP sample of the next trigger (end of file + 1 for the last trial)
%       trigWidthMs, msgLagS (TR message -> trigger, s)
%       pulseMs    {cell} pulse times in ms after the trigger (design + measured delay); empty for blank / visual
%       nPulses, firstPulseMs, lastPulseMs
%       isVisualDesign (condition kind is 'visual'), diodeN (frame pulses in the window), diodePreN (before the trigger),
%       isVisualDiode, diodeOnMs / diodeOffMs (ms after the trigger), visOnAP (AP sample of the visual onset)
%   TT.dropped  trials removed: trialNum, cond, reason
%   TT.artCheck per condition: trials checked, expected / detected pulses, missing, extra, worst timing error;
%               delayMs (median trigger -> pulse delay used for pulseMs)
%   TT.args, TT.chnList, TT.condsFile   copied from the conds file
%
% Written without access to MATLAB: first run it on one session and check the printout.
% Needs MATLAB R2018b or newer (jsondecode, contains, startsWith, sgtitle, xline).
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
% one session folder (the folder that contains 'Record Node xxx'), or a cell array of folders processed in turn
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
condsFile  = '';      % '' = the only *.conds.mat in the recording1 folder (a full path can only be given for ONE session)
outDir     = '';      % '' = save in the recording1 folder of each session
outName    = '';      % '' = short name from the stim file (e.g. CJ287_15Pulses3_100Hz); a name can only be typed for ONE session
overwrite  = false;   % true = replace an existing TrialTable (the old file is first copied to a dated backup)
do_save    = true;    % false = only check, print and plot; nothing is written

% ---------------- streams and lines ----------------
S.apMatch   = '-AP';      % text that identifies the AP (30 kHz spike band) stream folder
S.niMatch   = 'NI-DAQ';   % text that identifies the NI-DAQ stream folder (trigger, photodiode)
S.trigAI    = 3;          % NI analog input (0-based) carrying the trial trigger
S.diodeAI   = 0;          % NI analog input (0-based) carrying the photodiode
S.syncState = 1;          % TTL line of the 1 Hz sync pulse (the same line on the NI and probe streams)

% ---------------- trigger ----------------
S.trigThrV     = 2;        % V: AI3 above this = trigger high (the pulse is ~4.1 V)
S.trigMinGapMs = 500;      % ms: a rising edge closer than this to the previous one is ignored (glitch protection)
S.trigWidthMs  = [50 150]; % ms: expected trigger pulse width; trials outside are flagged in the printout (not removed)
S.maxMsgLagS   = 30;       % s: a trigger must come at most this long after its TR message, else it is flagged

% ---------------- photodiode ----------------
S.diodeHiV   = 0.10;   % V: the signal rising above this counts as one frame pulse (baseline ~0.015 V, pulses 0.17-0.7 V)
S.diodeLoV   = 0.05;   % V: it must fall below this before the next pulse is counted (hysteresis against noise)
S.diodeWinMs = 600;    % ms after the trigger searched for frame pulses
S.diodePreMs = 200;    % ms before the trigger checked for pulses that are NOT from the stimulus (screen flicker)
S.diodeMinN  = 5;      % at least this many frame pulses after the trigger = the photodiode says "visual trial"

% ---------------- artifact check (AP signal) ----------------
S.artTrialsPerCond = 2;     % trials per condition checked (the first ones); 0 = skip the check
S.artPreMs   = 10;          % ms before the trigger read to measure the quiet level
S.artK       = 15;          % threshold = artK x (median quiet level of the common step signal) + artAbsUV
S.artAbsUV   = 20;          % uV, see artK
S.artGroupMs = 0.8;         % ms: threshold crossings closer than this belong to the same pulse
S.artTolMs   = 1.0;         % ms: a detected pulse matches the design when it is within this of (design time + delay); loose on purpose
S.pulseDelay = 'measured';  % 'measured' = add the measured trigger->pulse delay to the design times; or a number (ms)

% ---------------- figures ----------------
S.showFigures = true;       % false = no figures
%% =================================================

if ischar(sessionDir) || isstring(sessionDir), sessionList = {char(sessionDir)}; else, sessionList = sessionDir; end
assert(isempty(condsFile) || numel(sessionList) == 1, 'condsFile can only be set when ONE session is processed.');
assert(isempty(outName) || numel(sessionList) == 1, 'outName can only be set when ONE session is processed.');
for iS = 1:numel(sessionList)
    fprintf('\n################ session %d of %d ################\n', iS, numel(sessionList));
    TT = buildTrialTable_(sessionList{iS}, condsFile, S);
    saveTT_(TT, outDir, outName, overwrite, do_save);
    if S.showFigures, plotTT_(TT); end
end

%% ======================================================================
%% ============================ MAIN FUNCTION ===========================
%% ======================================================================
function TT = buildTrialTable_(sessionDir, condsFile, S)
codeVersion = 'NP_1_TrialTable_v1 (2026-11-12, loose artifact check)';
tStart = tic;

%% ---- 1. folders, streams, probe positions --------------------------------------------
recDir = resolveRec_(sessionDir);
[sessName, nodeDir, expNum] = sessionParts_(recDir);
oe = jsondecode(fileread(fullfile(recDir, 'structure.oebin')));
cont = oe.continuous;  if isstruct(cont), cont = num2cell(cont); end
fol = cellfun(@(c) stripSep_(c.folder_name), cont, 'UniformOutput', false);
iAP = find(contains(fol, S.apMatch), 1);   iNI = find(contains(fol, S.niMatch), 1);
assert(~isempty(iAP) && ~isempty(iNI), 'AP or NI stream not found in structure.oebin: %s', strjoin(fol, ', '));
apDir = fol{iAP};  niDir = fol{iNI};
fsAP = double(cont{iAP}.sample_rate);  nCh = double(cont{iAP}.num_channels);
fsNI = double(cont{iNI}.sample_rate);  nNI = double(cont{iNI}.num_channels);
apCh = chanList_(cont{iAP});  niCh = chanList_(cont{iNI});
uVperBit  = double(apCh{1}.bit_volts);                       % AP: microvolts per bit
niVperBit = cellfun(@(c) double(c.bit_volts), niCh);          % NI: volts per bit, per input
apFile = fullfile(recDir, 'continuous', apDir, 'continuous.dat');
niFile = fullfile(recDir, 'continuous', niDir, 'continuous.dat');
dA = dir(apFile);  nS  = floor(dA.bytes/(2*nCh));             % AP samples in the file
dN = dir(niFile);  nSN = floor(dN.bytes/(2*nNI));             % NI samples in the file
apSN0 = readNPY_(fullfile(recDir, 'continuous', apDir, 'sample_numbers.npy'), 1);   % sample number of the first AP sample
niSN0 = readNPY_(fullfile(recDir, 'continuous', niDir, 'sample_numbers.npy'), 1);   % sample number of the first NI sample

PR = readProbe_(nodeDir, expNum, nCh);

%% ---- 2. clock sync NI -> probe ----------------------------------------------------------
SY = syncFit_(recDir, niDir, apDir, fsNI, fsAP, S.syncState, niSN0, apSN0);

%% ---- 3. triggers on AI3 ------------------------------------------------------------------
fprintf('reading NI inputs AI%d (trigger) and AI%d (photodiode) ...\n', S.trigAI, S.diodeAI);
R = readNIrows_(niFile, nNI, nSN, [S.trigAI S.diodeAI] + 1);
ai3 = single(R(1,:)) * single(niVperBit(S.trigAI+1));
ai0 = single(R(2,:)) * single(niVperBit(S.diodeAI+1));
clear R
hi = ai3 > S.trigThrV;
up = find(~hi(1:end-1) & hi(2:end)) + 1;              % rising edges (NI file sample, 1-based)
dn = find(hi(1:end-1) & ~hi(2:end)) + 1;              % falling edges
keep = true(size(up));  lastUp = -Inf;
for k = 1:numel(up)
    if (up(k) - lastUp)/fsNI*1e3 < S.trigMinGapMs, keep(k) = false; else, lastUp = up(k); end
end
up = up(keep);
nTrig = numel(up);
assert(nTrig > 0, 'No trigger found on AI%d.', S.trigAI);
trigWidthMs = nan(1, nTrig);
for k = 1:nTrig
    j = find(dn > up(k), 1);
    if ~isempty(j), trigWidthMs(k) = (dn(j) - up(k))/fsNI*1e3; end
end
trigAP_s = polyval(SY.P, (up - 1)/fsNI);              % s from the start of the AP file
trigAP   = round(trigAP_s*fsAP) + 1;                  % AP file sample (1-based)

% check against Open Ephys' own synchronised timestamps
tsNI = readNPYat_(fullfile(recDir, 'continuous', niDir, 'timestamps.npy'), up);
okI  = tsNI ~= -1 & trigAP >= 1 & trigAP <= nS;
tsAP = nan(size(tsNI));
tsAP(okI) = readNPYat_(fullfile(recDir, 'continuous', apDir, 'timestamps.npy'), trigAP(okI));
dTs = (tsNI - tsAP)*1e3;
SY.oeCheckN = sum(okI);
SY.oeCheckMedAbsMs = median(abs(dTs(okI)));
SY.oeCheckMaxAbsMs = max(abs(dTs(okI)));

%% ---- 4. trial messages and pairing -------------------------------------------------------
md  = fullfile(recDir, 'events', 'MessageCenter');
txt = readNPYtext_(fullfile(md, 'text.npy'));
msn = readNPY_(fullfile(md, 'sample_numbers.npy'));
trN = [];  trSN = [];  idN = [];  idC = [];  stimName = '';
for k = 1:numel(txt)
    t = txt{k};
    tok = regexp(t, '^TR:(\d+)$', 'tokens', 'once');
    if ~isempty(tok), trN(end+1) = str2double(tok{1}); trSN(end+1) = msn(k); continue; end %#ok<AGROW>
    tok = regexp(t, '^TRIALID\s+(\d+)-(\d+)', 'tokens', 'once');
    if ~isempty(tok), idC(end+1) = str2double(tok{1}); idN(end+1) = str2double(tok{2}); continue; end %#ok<AGROW>
    if startsWith(t, 'NEUROSTIM FILE') && isempty(stimName), stimName = strtrim(strrep(t, 'NEUROSTIM FILE', '')); end
end
assert(~isempty(trN), 'No TR:n trial messages in %s', md);
[trN, ia] = unique(trN);  trMsgAP = trSN(ia) - apSN0 + 1;      % each message is logged twice: keep the first
[idN, ib] = unique(idN);  idC = idC(ib);
[~, srt] = sort(trMsgAP);  trN = trN(srt);  trMsgAP = trMsgAP(srt);
msgCond = nan(size(trN));
[tf, loc] = ismember(trN, idN);  msgCond(tf) = idC(loc(tf));

pairMsg = nan(1, nTrig);                                       % index into trN for every trigger
for k = 1:nTrig
    j = find(trMsgAP <= trigAP(k), 1, 'last');
    if ~isempty(j), pairMsg(k) = j; end
end
msgLagS = nan(1, nTrig);
okP = ~isnan(pairMsg);
msgLagS(okP) = (trigAP(okP) - trMsgAP(pairMsg(okP)))/fsAP;
[cnt, ~] = histcounts(pairMsg(okP), 0.5:1:numel(trN)+0.5);
if any(cnt > 1)
    bad = find(cnt > 1);
    error(['%d trial message(s) got more than one trigger (e.g. trial %d). A trial message is probably missing; ' ...
           'check the recording before using it.'], numel(bad), trN(bad(1)));
end
if any(diff(pairMsg(okP)) <= 0)
    error('Triggers and trial messages are not in the same order (a message came before the previous trigger). Check this session.');
end
noMsgTrig = find(~okP);                                        % triggers before any message (should not happen)

% trials dropped: messages without trigger, and triggers without message/condition
dropped = struct('trialNum', {}, 'cond', {}, 'reason', {});
for j = setdiff(1:numel(trN), pairMsg(okP))
    dropped(end+1) = struct('trialNum', trN(j), 'cond', msgCond(j), 'reason', 'trial message without trigger'); %#ok<AGROW>
end
for k = noMsgTrig
    dropped(end+1) = struct('trialNum', NaN, 'cond', NaN, 'reason', sprintf('trigger %d (%.1f s) before any trial message', k, trigAP_s(k))); %#ok<AGROW>
end
useTrig = find(okP);
useTrig = useTrig(~isnan(msgCond(pairMsg(useTrig))));
for k = setdiff(find(okP), useTrig)
    dropped(end+1) = struct('trialNum', trN(pairMsg(k)), 'cond', NaN, 'reason', 'no TRIALID (condition unknown)'); %#ok<AGROW>
end

%% ---- 5. design from the conds file ------------------------------------------------------
if isempty(condsFile)
    f = dir(fullfile(recDir, '*.conds.mat'));  f = f(~startsWith({f.name}, '._'));
    assert(numel(f) == 1, 'Found %d *.conds.mat in %s: copy the stim-PC conds file there, or set condsFile.', numel(f), recDir);
    condsFile = fullfile(f(1).folder, f(1).name);
end
C = load(condsFile, 'cond', 'info', 'args', 'chnList');
nCond = numel(C.cond);
assert(numel(C.info) == nCond, 'conds file: cond and info have different lengths.');
assert(max(msgCond(pairMsg(useTrig))) <= nCond, 'Trial messages use condition %d but the conds file has only %d.', ...
       max(msgCond(pairMsg(useTrig))), nCond);
freq = double(C.args.freq);  periodMs = 1000/freq;
cond = struct('id', {}, 'label', {}, 'kind', {}, 'level', {}, 'nLevels', {}, 'isFinal', {}, 'amp', {}, 'ptd_us', {}, ...
              'nsp', {}, 'designMs', {}, 'designMsA', {}, 'designMsB', {}, 'contrast', {}, 'phase', {}, 'nTrials', {});
for c = 1:nCond
    cc = C.cond(c);  ii = C.info(c);
    amp = double(cc.amp(:)');  ptd = double(cc.ptd(:)');  nsp = double(cc.nsp(:)');
    pe = {zeros(1,0), zeros(1,0)};
    for e = 1:numel(amp)
        if amp(e) > 0 && nsp(e) > 0, pe{e} = ptd(e)/1000 + (0:nsp(e)-1)*periodMs; end
    end
    allP = sort([pe{:}]);
    if ~isempty(allP), allP = allP([true, diff(allP) > 1e-6]); end     % simultaneous pulses counted once
    cond(c) = struct('id', c, 'label', char(cc.label), 'kind', char(ii.kind), 'level', double(ii.level), ...
        'nLevels', double(ii.nLevels), 'isFinal', logical(ii.isFinal), 'amp', amp, 'ptd_us', ptd, 'nsp', nsp, ...
        'designMs', allP, 'designMsA', pe{1}, 'designMsB', pe{2}, 'contrast', getf_(ii, 'contrast', NaN), ...
        'phase', getf_(ii, 'phase', NaN), 'nTrials', 0);
end
kinds = {cond.kind};
pattern = struct('kind', {}, 'conds', {}, 'levels', {}, 'finalCond', {});
for pk = {'sim', 'AB', 'BA'}
    cs = find(strcmp(kinds, pk{1}));
    if isempty(cs), continue; end
    [lv, o] = sort([cond(cs).level]);  cs = cs(o);
    fc = cs([cond(cs).isFinal]);  if isempty(fc), fc = cs(end); end
    pattern(end+1) = struct('kind', pk{1}, 'conds', cs, 'levels', lv, 'finalCond', fc(1)); %#ok<AGROW>
end

%% ---- 6. photodiode per trigger -----------------------------------------------------------
nWin = round(S.diodeWinMs*fsNI/1e3);  nPre = round(S.diodePreMs*fsNI/1e3);
diodeN = zeros(1, nTrig);  diodePreN = zeros(1, nTrig);  diodeOnMs = nan(1, nTrig);  diodeLastMs = nan(1, nTrig);
allGaps = [];
for k = 1:nTrig
    a = ai0(up(k) : min(nSN, up(k) + nWin));
    e = edgesHyst_(a, S.diodeLoV, S.diodeHiV);
    b = ai0(max(1, up(k) - nPre) : up(k) - 1);
    diodePreN(k) = numel(edgesHyst_(b, S.diodeLoV, S.diodeHiV));
    diodeN(k) = numel(e);
    if ~isempty(e)
        diodeOnMs(k) = (e(1) - 1)/fsNI*1e3;  diodeLastMs(k) = (e(end) - 1)/fsNI*1e3;
        if numel(e) >= S.diodeMinN, allGaps = [allGaps, diff(e)/fsNI*1e3]; end %#ok<AGROW>
    end
end
framePeriodMs = median(allGaps);  if isempty(framePeriodMs) || isnan(framePeriodMs), framePeriodMs = 1000/120; end
diodeOffMs = diodeLastMs + framePeriodMs;
isVisualDiode = diodeN >= S.diodeMinN;
visOnAP = nan(1, nTrig);
iv = find(isVisualDiode);
visOnAP(iv) = round(polyval(SY.P, (up(iv) - 1 + round(diodeOnMs(iv)*fsNI/1e3))/fsNI)*fsAP) + 1;
clear ai0 ai3

%% ---- 7. artifact check in the AP signal --------------------------------------------------
trialCondAll = nan(1, nTrig);  trialCondAll(useTrig) = msgCond(pairMsg(useTrig));
maxLast = max([0, cellfun(@(v) max([v 0]), {cond.designMs})]);
AC = struct('cond', num2cell(1:nCond), 'trialsChecked', [], 'nExpected', NaN, 'nDetected', [], 'nMissing', NaN, ...
            'nExtra', NaN, 'maxErrMs', NaN, 'ok', false, 'detMs', []);
delayMs = 0;  delayN = 0;
if S.artTrialsPerCond > 0
    mm = memmapfile(apFile, 'Format', {'int16', [nCh nS], 'x'});
    preS = round(S.artPreMs*fsAP/1e3);  postS = round((maxLast + 15)*fsAP/1e3);
    gapS = max(1, round(S.artGroupMs*fsAP/1e3));
    fprintf('artifact check: %d trial(s) per condition ...\n', S.artTrialsPerCond);
    det = cell(nCond, S.artTrialsPerCond);
    for c = 1:nCond
        tk = find(trialCondAll == c, S.artTrialsPerCond);
        AC(c).trialsChecked = tk;
        for q = 1:numel(tk)
            i0 = trigAP(tk(q));  i1 = i0 - preS;  i2 = min(nS, i0 + postS);
            if i1 < 1, continue; end
            X  = double(mm.Data.x(:, i1:i2));
            cm = median(abs(diff(X, 1, 2)), 1) * uVperBit;      % common step signal: the artifact steps on all channels at once
            thr = S.artK*median(cm(1:preS-5)) + S.artAbsUV;
            idx = find(cm > thr);
            if isempty(idx), det{c,q} = zeros(1,0); continue; end
            g = idx([true, diff(idx) > gapS]);
            det{c,q} = (g - preS)/fsAP*1e3;                     % ms after the trigger
        end
    end
    % pass 1: delay trigger -> FIRST pulse (the later pulses drift slightly, so only the first one is used)
    dd = [];
    for c = 1:nCond
        for q = 1:size(det, 2)
            if isempty(det{c,q}) || isempty(cond(c).designMs), continue; end
            e = cond(c).designMs(1);
            [m, j] = min(abs(det{c,q} - e));
            if m <= 1.5, dd(end+1) = det{c,q}(j) - e; end %#ok<AGROW>
        end
    end
    if ~isempty(dd), delayMs = median(dd);  delayN = numel(dd); end
    % pass 2: per condition, missing / extra pulses with the measured delay
    for c = 1:nCond
        ex = cond(c).designMs + delayMs;
        AC(c).nExpected = numel(ex);  AC(c).nMissing = 0;  AC(c).nExtra = 0;  AC(c).maxErrMs = 0;
        nd = [];
        for q = 1:numel(AC(c).trialsChecked)
            d = det{c,q};  if isempty(d), d = zeros(1,0); end
            nd(end+1) = numel(d); %#ok<AGROW>
            used = false(size(d));
            for e = ex
                [m, j] = min(abs(d - e));
                if ~isempty(m) && m <= S.artTolMs && ~used(j)
                    used(j) = true;  AC(c).maxErrMs = max(AC(c).maxErrMs, m);
                else
                    AC(c).nMissing = AC(c).nMissing + 1;
                end
            end
            AC(c).nExtra = AC(c).nExtra + sum(~used);
        end
        AC(c).nDetected = nd;
        if ~isempty(AC(c).trialsChecked), AC(c).detMs = det{c,1}; end
        AC(c).ok = ~isempty(AC(c).trialsChecked) && AC(c).nMissing == 0 && AC(c).nExtra == 0;
    end
end
if ischar(S.pulseDelay) && strcmpi(S.pulseDelay, 'measured')
    if delayN == 0, warning('NPTT:delay', 'No delay measured (artifact check off or no pulses found): delay = 0 ms.'); end
    useDelay = delayMs;
else
    useDelay = double(S.pulseDelay);
end

%% ---- 8. the trial table ------------------------------------------------------------------
nT = numel(useTrig);
tr = struct();
tr.trialNum  = trN(pairMsg(useTrig))';
tr.cond      = msgCond(pairMsg(useTrig))';
tr.rep       = zeros(nT, 1);
for c = 1:nCond
    w = find(tr.cond == c);  tr.rep(w) = (1:numel(w))';  cond(c).nTrials = numel(w);
end
tr.trigNI      = up(useTrig)';
tr.trigAP      = trigAP(useTrig)';
tr.trigAP_s    = trigAP_s(useTrig)';
nxt = [trigAP(2:end), nS + 1];
tr.nextTrigAP  = nxt(useTrig)';
tr.trigWidthMs = trigWidthMs(useTrig)';
tr.msgLagS     = msgLagS(useTrig)';
tr.pulseMs     = arrayfun(@(c) cond(c).designMs + useDelay, tr.cond, 'UniformOutput', false);
tr.nPulses     = cellfun(@numel, tr.pulseMs);
tr.firstPulseMs = cellfun(@(v) min([v NaN]), tr.pulseMs);
tr.lastPulseMs  = cellfun(@(v) max([v NaN]), tr.pulseMs);
tr.isVisualDesign = strcmp(kinds(tr.cond), 'visual')';
tr.diodeN      = diodeN(useTrig)';
tr.diodePreN   = diodePreN(useTrig)';
tr.isVisualDiode = isVisualDiode(useTrig)';
tr.diodeOnMs   = diodeOnMs(useTrig)';   tr.diodeOnMs(~tr.isVisualDiode) = NaN;
tr.diodeOffMs  = diodeOffMs(useTrig)';  tr.diodeOffMs(~tr.isVisualDiode) = NaN;
tr.visOnAP     = visOnAP(useTrig)';

%% ---- assemble ---------------------------------------------------------------------------------
TT = struct();
[~, cn] = fileparts(condsFile);  shortName = strtok(cn, '.');          % e.g. CJ287_15Pulses3_100Hz
if isempty(shortName), shortName = sessName; end
TT.info = struct('session', sessName, 'shortName', shortName, 'recDir', recDir, 'nodeDir', nodeDir, 'apFile', apFile, 'niFile', niFile, ...
    'fsAP', fsAP, 'fsNI', fsNI, 'nChAP', nCh, 'nSampAP', nS, 'nSampNI', nSN, 'uVperBit', uVperBit, ...
    'apSN0', apSN0, 'niSN0', niSN0, 'stimName', stimName, 'condsFile', condsFile, 'framePeriodMs', framePeriodMs, ...
    'delayMs', useDelay, 'settings', S, 'version', codeVersion, 'date', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
TT.sync = SY;
TT.chan = PR;
TT.cond = cond;
TT.pattern = pattern;
TT.trial = tr;
TT.dropped = dropped;
TT.artCheck = AC;
TT.artDelay = struct('delayMs', delayMs, 'nPulses', delayN);
TT.args = C.args;
TT.chnList = C.chnList;
TT.condsFile = condsFile;

%% ---- printout ---------------------------------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('NP TRIAL TABLE (%s)\n', codeVersion);
fprintf('============================================================\n');
fprintf('Session   : %s  (short name %s)\nRecording : %s\n', sessName, shortName, recDir);
fprintf('Streams   : AP %d ch @ %g Hz, %.1f s (%.4f uV/bit) | NI %d ch @ %g Hz, %.1f s\n', nCh, fsAP, nS/fsAP, uVperBit, nNI, fsNI, nSN/fsNI);
fprintf('Probe     : preset "%s", hardware filterCut "%s", y = %d..%d um from the tip\n', PR.preset, PR.filterCut, min(PR.y_um), max(PR.y_um));
fprintf('            channel number is NOT depth order: ch 1-4 sit at y = %s um. Use TT.chan.orderTipToTop.\n', mat2str(PR.y_um(1:min(4,nCh))'));
fprintf('Sync      : %d of %d NI pulses matched, slope %+.2f ppm, residual SD %.4f ms\n', SY.nMatched, SY.nNI, SY.slopePpm, SY.residSDms);
fprintf('            vs Open Ephys timestamps: median |diff| %.3f ms, max %.3f ms (%d triggers)\n', SY.oeCheckMedAbsMs, SY.oeCheckMaxAbsMs, SY.oeCheckN);
fprintf('Triggers  : %d on AI%d, width median %.1f ms [%.1f %.1f]\n', nTrig, S.trigAI, median(trigWidthMs, 'omitnan'), min(trigWidthMs), max(trigWidthMs));
wBad = find(trigWidthMs < S.trigWidthMs(1) | trigWidthMs > S.trigWidthMs(2) | isnan(trigWidthMs));
if ~isempty(wBad), fprintf('            %d trigger(s) with unusual width: %s\n', numel(wBad), mat2str(wBad(1:min(20,end)))); end
fprintf('Messages  : %d trials (TR %d..%d); %d paired with a trigger; message -> trigger lag median %.2f s [%.2f %.2f]\n', ...
        numel(trN), min(trN), max(trN), sum(okP), median(msgLagS, 'omitnan'), min(msgLagS), max(msgLagS));
lagBad = find(msgLagS > S.maxMsgLagS);
if ~isempty(lagBad), fprintf('            %d trigger(s) more than %g s after their message: %s\n', numel(lagBad), S.maxMsgLagS, mat2str(lagBad(1:min(20,end)))); end
if isempty(dropped)
    fprintf('Dropped   : none\n');
else
    fprintf('Dropped   : %d\n', numel(dropped));
    for q = 1:numel(dropped), fprintf('            trial %g (cond %g): %s\n', dropped(q).trialNum, dropped(q).cond, dropped(q).reason); end
end
kk = unique(kinds, 'stable');
fprintf('Design    : %s\n            %d conditions (%s), %g Hz, electrodes A = %s, B = %s\n', condsFile, nCond, ...
        strjoin(cellfun(@(k) sprintf('%s %d', k, sum(strcmp(kinds, k))), kk, 'UniformOutput', false), ', '), freq, ...
        char(C.chnList{1}), char(C.chnList{min(2, end)}));
nTr = [cond.nTrials];
fprintf('Kept      : %d trials; per condition %d..%d', nT, min(nTr), max(nTr));
if min(nTr) ~= max(nTr)
    fprintf(' (conditions with fewer than %d: %s)', max(nTr), mat2str(find(nTr < max(nTr))));
end
fprintf('\n');
for p = 1:numel(pattern)
    fprintf('            pattern %-4s: conds %d..%d, %d levels, final cond %d\n', pattern(p).kind, min(pattern(p).conds), ...
            max(pattern(p).conds), numel(pattern(p).conds), pattern(p).finalCond);
end
vd = tr.isVisualDesign;  vp = tr.isVisualDiode;
fprintf('Photodiode: visual by design %d, by photodiode %d, agree %d of %d; frame period %.2f ms\n', ...
        sum(vd), sum(vp), sum(vd == vp), nT, framePeriodMs);
if any(vd), fprintf('            visual onset %.1f ms median [%.1f %.1f] after the trigger, offset %.1f ms median\n', ...
        median(tr.diodeOnMs(vd & vp)), min(tr.diodeOnMs(vd & vp)), max(tr.diodeOnMs(vd & vp)), median(tr.diodeOffMs(vd & vp))); end
dis = find(vd ~= vp);
for q = dis(:)'
    fprintf('            DISAGREE: trial %d cond %d (%s): %d frame pulses\n', tr.trialNum(q), tr.cond(q), cond(tr.cond(q)).label, tr.diodeN(q));
end
nPreFl = sum(tr.diodePreN >= S.diodeMinN);
if nPreFl > 0
    w = find(tr.diodePreN >= S.diodeMinN);
    fprintf('            %d trial(s) have photodiode pulses in the %g ms BEFORE the trigger (harmless if few): trial %s\n', ...
            nPreFl, S.diodePreMs, mat2str(tr.trialNum(w(1:min(10,end)))'));
end
if S.artTrialsPerCond > 0
    okC = [AC.ok];  chk = arrayfun(@(a) ~isempty(a.trialsChecked), AC);
    fprintf('Artifacts : delay trigger -> first pulse %.3f ms (median of %d trials)\n', delayMs, delayN);
    fprintf('            %d of %d conditions OK (right number of pulses, each within %.1f ms of the design)\n', ...
            sum(okC & chk), sum(chk), S.artTolMs);
    lateMs = arrayfun(@(a) a.maxErrMs, AC(okC & chk));
    if ~isempty(lateMs), fprintf('            largest timing difference in the OK conditions: %.2f ms (usually the last pulse of the longest train)\n', max(lateMs)); end
    for c = find(chk & ~okC)
        fprintf('            cond %2d %-28s expected %2d, detected %s, missing %d, extra %d\n', c, cond(c).label, ...
                AC(c).nExpected, mat2str(AC(c).nDetected), AC(c).nMissing, AC(c).nExtra);
    end
end
fprintf('Pulse times in TT.trial.pulseMs = design + %.3f ms\n', useDelay);
fprintf('Done in %.0f s.\n', toc(tStart));
end

%% ======================================================================
%% ============================ SAVE / PLOT =============================
%% ======================================================================
function saveTT_(TT, outDir, outName, overwrite, do_save)
if isempty(outDir), od = TT.info.recDir; else, od = outDir; end
if isempty(outName), outName = TT.info.shortName; end
outName = regexprep(outName, '(\.NP_TrialTable)?\.mat$', '');
outFile = fullfile(od, [outName '.NP_TrialTable.mat']);
fprintf('Output    : %s\n', outFile);
if ~do_save, fprintf('do_save = false: nothing saved.\n'); return; end
if ~isfolder(od), mkdir(od); end
if isfile(outFile)
    if ~overwrite
        error('NPTT:exists', 'The file already exists and overwrite = false. Nothing was changed.');
    end
    bak = strrep(outFile, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
    copyfile(outFile, bak);
    fprintf('Old file copied to %s\n', bak);
end
save(outFile, 'TT');
fprintf('Saved.\n');
end

function plotTT_(TT)
tr = TT.trial;  nm = strrep(TT.info.session, '_', '\_');
% ---- figure 1: trials ----
figure('Name', ['Trials: ' TT.info.session], 'Color', 'w');
subplot(3,1,1); plot(tr.trialNum, tr.msgLagS, 'k.'); ylabel('message -> trigger (s)'); grid on;
title(sprintf('%d kept trials, %d dropped', numel(tr.trialNum), numel(TT.dropped)));
subplot(3,1,2); plot(tr.trialNum, tr.trigWidthMs, 'k.'); ylabel('trigger width (ms)'); grid on;
subplot(3,1,3); hold on;
plot(tr.trialNum(~tr.isVisualDesign), tr.diodeN(~tr.isVisualDesign), 'k.');
plot(tr.trialNum(tr.isVisualDesign), tr.diodeN(tr.isVisualDesign), 'r.');
ylabel('photodiode pulses'); xlabel('trial number'); legend({'electrical / blank', 'visual (design)'}, 'Location', 'best'); grid on;
sgtitle(nm);
% ---- figure 2: artifact check ----
AC = TT.artCheck;
if any(arrayfun(@(a) ~isempty(a.trialsChecked), AC))
    figure('Name', ['Artifact check: ' TT.info.session], 'Color', 'w'); hold on;
    nC = numel(TT.cond);
    for c = 1:nC
        ex = TT.cond(c).designMs + TT.info.delayMs;
        plot(ex, c*ones(size(ex)), 'o', 'Color', [0.85 0.2 0.2], 'MarkerSize', 5);
        d = AC(c).detMs;  plot(d, c*ones(size(d)), 'k+', 'MarkerSize', 6);
        if ~AC(c).ok && ~isempty(AC(c).trialsChecked), plot(-3, c, 'r<', 'MarkerFaceColor', 'r'); end
    end
    set(gca, 'YDir', 'reverse', 'YTick', 1:nC, 'YTickLabel', {TT.cond.label}, 'FontSize', 7);
    xlabel('ms after the trigger'); ylim([0 nC+1]); grid on;
    title(sprintf('%s: red o = design + %.2f ms, black + = detected (first checked trial); red < = wrong pulse count or > %.1f ms off', nm, TT.info.delayMs, TT.info.settings.artTolMs));
end
% ---- figure 3: probe channels ----
figure('Name', ['Probe channels: ' TT.info.session], 'Color', 'w'); hold on;
ch = TT.chan;  b = ch.bank;  ub = unique(b);  col = lines(max(numel(ub), 1));
for q = 1:numel(ub)
    w = b == ub(q);  plot(ch.x_um(w), ch.y_um(w), 's', 'Color', col(q,:), 'MarkerFaceColor', col(q,:), 'MarkerSize', 4);
end
for c = 1:16:numel(ch.ch), text(ch.x_um(c) + 8, ch.y_um(c), sprintf('ch %d', ch.ch(c)), 'FontSize', 7); end
xlabel('x (um)'); ylabel('y from the tip (um)'); grid on;
legend(arrayfun(@(v) sprintf('bank %d', v), ub, 'UniformOutput', false), 'Location', 'bestoutside');
title(sprintf('%s: preset "%s" (labels every 16th channel)', nm, ch.preset));
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function rec = resolveRec_(p)
% accepts the session folder, the Record Node / experiment folder, or the recording folder itself
p = regexprep(char(p), '[\\/]+$', '');
assert(isfolder(p), 'Folder not found: %s', p);
if isfile(fullfile(p, 'structure.oebin')), rec = p; return; end
c = dir(fullfile(p, '**', 'structure.oebin'));  c = c(~startsWith({c.name}, '._'));
assert(~isempty(c), 'No structure.oebin below %s', p);
assert(numel(c) == 1, 'Several recordings below %s - give the recording folder:\n%s', p, strjoin(fullfile({c.folder}, {c.name}), newline));
rec = c(1).folder;
end

function [sessName, nodeDir, expNum] = sessionParts_(recDir)
parts = strsplit(recDir, filesep);
k = find(startsWith(parts, 'Record Node'), 1, 'last');
assert(~isempty(k) && k > 1, 'Cannot find the "Record Node" folder in %s', recDir);
sessName = parts{k-1};
nodeDir  = strjoin(parts(1:k), filesep);
expNum = 1;
if numel(parts) > k
    t = regexp(parts{k+1}, 'experiment(\d+)', 'tokens', 'once');
    if ~isempty(t), expNum = str2double(t{1}); end
end
end

function s = stripSep_(s)
s = regexprep(char(s), '[\\/]+$', '');
end

function L = chanList_(stream)
L = stream.channels;
if isstruct(L), L = num2cell(L); end
end

function v = getf_(s, f, def)
if isfield(s, f) && ~isempty(s.(f)), v = double(s.(f)); else, v = def; end
end

function PR = readProbe_(nodeDir, expNum, nCh)
% channel positions (um) and banks from settings.xml of the Record Node
f = fullfile(nodeDir, 'settings.xml');
if expNum > 1 && isfile(fullfile(nodeDir, sprintf('settings_%d.xml', expNum))), f = fullfile(nodeDir, sprintf('settings_%d.xml', expNum)); end
assert(isfile(f), 'No settings.xml in %s', nodeDir);
txt = fileread(f);
yT = regexp(txt, '<ELECTRODE_YPOS\s([^/>]*)/?>', 'tokens');
xT = regexp(txt, '<ELECTRODE_XPOS\s([^/>]*)/?>', 'tokens');
bT = regexp(txt, '<CHANNELS\s([^/>]*)/?>', 'tokens');
assert(~isempty(yT), 'No ELECTRODE_YPOS in %s', f);
if numel(yT) > 1, warning('NPTT:probes', '%d probes in settings.xml: the first one is used.', numel(yT)); end
y = parseCH_(yT{1}{1}, nCh);
x = zeros(nCh, 1);  if ~isempty(xT), x = parseCH_(xT{1}{1}, nCh); end
b = zeros(nCh, 1);  if ~isempty(bT), b = parseCH_(bT{1}{1}, nCh); end
assert(~any(isnan(y)), 'settings.xml has no y position for some of the %d channels.', nCh);
[~, o] = sortrows([y(:) (1:nCh)']);
PR = struct('ch', (1:nCh)', 'x_um', x, 'y_um', y, 'bank', b, 'orderTipToTop', o(:), ...
    'preset', tokOrQ_(txt, 'electrodeConfigurationPreset'), 'filterCut', tokOrQ_(txt, 'filterCut'), 'settingsFile', f);
end

function v = parseCH_(s, nCh)
tk = regexp(s, 'CH(\d+)="([^"]*)"', 'tokens');
v = nan(nCh, 1);
for q = 1:numel(tk)
    i = str2double(tk{q}{1}) + 1;
    val = tk{q}{2};  c = strfind(val, ':');  if ~isempty(c), val = val(c(end)+1:end); end
    if i >= 1 && i <= nCh, v(i) = str2double(val); end
end
end

function s = tokOrQ_(txt, name)
t = regexp(txt, [name '="([^"]*)"'], 'tokens', 'once');
if isempty(t), s = '?'; else, s = t{1}; end
end

function SY = syncFit_(rec, niDir, apDir, fsNI, fsAP, state, niSN0, apSN0)
% linear fit: probe time (s from AP file start) = polyval(P, NI time (s from NI file start)), from the 1 Hz sync pulses
sdir = fullfile(rec, 'events', niDir, 'TTL');  pdir = fullfile(rec, 'events', apDir, 'TTL');
dsn = readNPY_(fullfile(sdir, 'sample_numbers.npy'));  dst = readNPY_(fullfile(sdir, 'states.npy'));
psn = readNPY_(fullfile(pdir, 'sample_numbers.npy'));  pst = readNPY_(fullfile(pdir, 'states.npy'));
sd = (dsn(dst == state) - niSN0)/fsNI;
sp = (psn(pst == state) - apSN0)/fsAP;
assert(numel(sd) >= 20 && numel(sp) >= 20, 'Too few sync pulses (NI %d, probe %d).', numel(sd), numel(sp));
j  = interp1(sp, (1:numel(sp))', sd, 'nearest', 'extrap');
dt = sp(j) - sd;
ok = abs(dt - median(dt)) < 0.05;
[~, uj] = unique(j(ok), 'stable');  ia = find(ok);  ia = ia(uj);
assert(numel(ia) >= 20, 'Only %d sync pulses matched.', numel(ia));
xS = sd(ia);  yS = sp(j(ia));
keep = true(size(xS));
for it = 1:5
    Pk = polyfit(xS(keep), yS(keep), 1);
    rk = (yS - polyval(Pk, xS))*1e3;
    s = 1.4826*median(abs(rk - median(rk)));
    nk = abs(rk - median(rk)) <= 5*max(s, 1e-3);
    if isequal(nk, keep), break; end
    keep = nk;
end
P = polyfit(xS(keep), yS(keep), 1);
res = (yS - polyval(P, xS))*1e3;
SY = struct('P', P, 'nNI', numel(sd), 'nProbe', numel(sp), 'nMatched', sum(keep), 'slopePpm', (P(1)-1)*1e6, ...
            'residSDms', std(res(keep)), 'residMaxMs', max(abs(res(keep))));
end

function R = readNIrows_(file, nRow, nSamp, rows)
% reads only the wanted rows of an interleaved int16 file, in chunks
R = zeros(numel(rows), nSamp, 'int16');
fid = fopen(file, 'r');  cl = onCleanup(@() fclose(fid));
chunk = 5e6;  pos = 0;
while pos < nSamp
    n = min(chunk, nSamp - pos);
    X = fread(fid, [nRow n], '*int16');
    if isempty(X), break; end
    R(:, pos+1:pos+size(X,2)) = X(rows, :);
    pos = pos + size(X, 2);
end
end

function e = edgesHyst_(seg, lo, hi)
% rising crossings of hi that follow a dip below lo
st = zeros(size(seg));  st(seg > hi) = 1;  st(seg < lo) = -1;
nz = find(st ~= 0);  s2 = st(nz);
e = nz(find(diff(s2) == 2) + 1);
end

function [fid, prec, nb, endian, nEl, off, descr] = openNPY_(fpath)
assert(isfile(fpath), 'Missing: %s', fpath);
fid = fopen(fpath, 'r');
magic = fread(fid, 6, '*uint8')';
assert(isequal(magic, uint8([147 78 85 77 80 89])), 'Not a .npy file: %s', fpath);
major = fread(fid, 1, 'uint8');  fread(fid, 1, 'uint8');
if major == 1, hlen = fread(fid, 1, 'uint16', 'l'); else, hlen = fread(fid, 1, 'uint32', 'l'); end
header = fread(fid, hlen, '*char')';
off = ftell(fid);
descr = regexp(header, '''descr''\s*:\s*''([^'']+)''', 'tokens', 'once');  descr = descr{1};
shpT  = regexp(header, '''shape''\s*:\s*\(([^)]*)\)', 'tokens', 'once');
shp   = str2num(['[' strrep(shpT{1}, ',', ' ') ']']); %#ok<ST2NM>
if isempty(shp), shp = 1; end
nEl = prod(shp);
endian = 'l';  if startsWith(descr, '>'), endian = 'b'; end
switch descr(2:end)
    case 'i8', prec = 'int64=>double';  nb = 8;   case 'i4', prec = 'int32=>double';  nb = 4;
    case 'i2', prec = 'int16=>double';  nb = 2;   case 'i1', prec = 'int8=>double';   nb = 1;
    case 'u8', prec = 'uint64=>double'; nb = 8;   case 'u4', prec = 'uint32=>double'; nb = 4;
    case 'u2', prec = 'uint16=>double'; nb = 2;   case 'u1', prec = 'uint8=>double';  nb = 1;
    case 'f8', prec = 'double=>double'; nb = 8;   case 'f4', prec = 'single=>double'; nb = 4;
    otherwise, prec = '';  nb = NaN;               % text types are handled by readNPYtext_
end
end

function A = readNPY_(fpath, maxCount)
if nargin < 2, maxCount = Inf; end
[fid, prec, ~, endian, nEl, ~, descr] = openNPY_(fpath);  cl = onCleanup(@() fclose(fid));
assert(~isempty(prec), 'Unsupported dtype "%s" in %s', descr, fpath);
A = fread(fid, min(nEl, maxCount), prec, 0, endian);
end

function A = readNPYat_(fpath, idx)
% reads only the elements idx (1-based) of a numeric .npy file
[fid, prec, nb, endian, nEl, off, descr] = openNPY_(fpath);  cl = onCleanup(@() fclose(fid));
assert(~isempty(prec), 'Unsupported dtype "%s" in %s', descr, fpath);
A = nan(size(idx));
for q = 1:numel(idx)
    if idx(q) < 1 || idx(q) > nEl, continue; end
    fseek(fid, off + (idx(q)-1)*nb, 'bof');
    A(q) = fread(fid, 1, prec, 0, endian);
end
end

function T = readNPYtext_(fpath)
% fixed-width byte strings ('|S513') -> cell array of char
[fid, ~, ~, ~, nEl, ~, descr] = openNPY_(fpath);  cl = onCleanup(@() fclose(fid));
w = regexp(descr, '^\|S(\d+)$', 'tokens', 'once');
assert(~isempty(w), 'Expected a byte-string .npy (|S...), got "%s" in %s', descr, fpath);
w = str2double(w{1});
B = fread(fid, [w nEl], '*uint8');
T = cell(nEl, 1);
for q = 1:nEl
    b = B(:, q)';  z = find(b == 0, 1);  if ~isempty(z), b = b(1:z-1); end
    T{q} = strtrim(char(b));
end
end
