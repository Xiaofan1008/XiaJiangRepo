%% ============================================================
% NP_2_SpikeDetect_v1.m   (step 2 of the NP pipeline)
% NP probe (Neuropixels, Open Ephys) | marmoset pulse-TRAIN sessions | ARTIFACT BLANKING + SPIKE DETECTION | v1
%
% PURPOSE
%   Detects multi-unit spikes on all AP channels in a window around every trial, after blanking the
%   stimulation artifact the NN way. Reads the trial table of step 1 (NP_1_TrialTable_v1); the AP file is
%   read one trial window at a time, nothing large is written, the recording is never changed.
%   No figures: rasters / PSTHs are drawn by later scripts, after the train recovery.
%
% PER TRIAL (window -preMs ... +postMs around the trigger, plus padMs on both sides for the filter)
%   1. Raw AP signal in uV, minus each channel's median before the trigger.
%   2. Blanking (only trials with pulses; blank and visual trials are untouched):
%        blankMode = 'nn'       (DEFAULT, = simpleBlank_PTD_PulseTrain of the NN pipeline): ONE blank per train,
%                               from blankPreMs before the FIRST pulse to the recovery point after the LAST pulse.
%        blankMode = 'perPulse' the same rule around EACH pulse on its own (blanks that touch are merged).
%      Recovery point, per channel: the first 46-sample window (recWin) after the pulse whose max-min range is
%      <= rangeMaxUV, + bufferSamp, but never more than capSamp after the pulse.
%      Fill: straight line across the blank; the DC step at the blank end is removed from everything after it.
%   3. Band-pass bpHz (Butterworth, forwards + backwards = zero phase), every channel on its own (no common reference).
%   4. Threshold = -thrCoeff x noise; noise = median(|x|)/0.6745 of the band-passed signal BEFORE the trigger,
%      per channel, median over nNoiseTrials trials spread over the session (measured once, before detection).
%   5. Spike time = most negative point within peakLookMs after the crossing; refractoryMs; snippet checks
%      (peak-to-peak and absolute limits) remove events too small or too large to be spikes.
%      Spikes inside a blank (+/- guardSpkMs) are removed.
%
% OUTPUT FILES (in the folder of the trial table, i.e. recording1, unless outDir is set)
%   <short name>.NP_Spikes.mat            (small; the file the next steps read)
%     sp_clipped{ch}      spike times in ms from the START OF THE AP FILE (column), as sp_xia.mat in the NN pipeline.
%                         Relative to a trial: t - (TT.trial.trigAP(k)-1)/FS*1000
%     sp_trial{ch}        row of TT.trial each spike belongs to (same length as sp_clipped{ch})
%     nSpk(ch, k)         spikes per channel and trial
%     blank_ms{ch, k}     blanked stretches of channel ch in trial k: [start end] in ms after the trigger, one row
%                         per blank ([] = no blank). 'nn' mode: one row.
%     blankStart_ms, blankEnd_ms (ch x trial)   start of the first / end of the last blank (NaN = no blank)
%     blankCapped (ch x trial)  true if a blank of that channel and trial hit capSamp (artifact never settled)
%     noise_uV, thr_uV (ch)     noise and threshold per channel;  goodCh  channels that were processed
%     FS, preMs, postMs, ttFile, settings S, version, date
%   <short name>.NP_Spikes_waveforms.mat  ONLY if saveWaveforms = true (can be large; delete it whenever you like):
%     sp_waveforms{ch}    int16 matrix, one row per spike (same order as sp_clipped{ch}), uV = value * wfScale_uV
%     wfScale_uV, wfTime_ms (time of each sample relative to the spike peak)
%
% Written without access to MATLAB: run it on one session first (nTrialsMax = 20 for a quick test).
% Needs the Signal Processing Toolbox (butter, filtfilt) and MATLAB R2018b or newer (movmax / movmin).
% ============================================================
clear; clc;

%% ================= USER SETTINGS =================
% session folder (the one containing 'Record Node xxx'), or a cell array of session folders processed in turn.
% The trial table of step 1 (*.NP_TrialTable.mat) must be in its recording1 folder.
sessionDir = '/Volumes/LaCie/Marmoset_Data/CJ287/Xia/CJ287_15Pulses3_100Hz_expTrainFreq_100Hz_080922_2026-11-10_08-09-29';
ttFile     = '';      % full path of the trial table; '' = the only *.NP_TrialTable.mat in recording1 (ONE session only)
outDir     = '';      % '' = save next to the trial table (recording1)
overwrite  = false;   % true = replace existing output (the old files are first copied to dated backups)
nTrialsMax = Inf;     % process only the first N trials (e.g. 20 for a quick test; then nothing is saved unless do_save)
do_save    = true;    % false = run and print only

% ---------------- trial window ----------------
S.preMs   = 300;      % ms BEFORE the trigger kept (baseline, noise)
S.postMs  = 700;      % ms AFTER the trigger kept (longest train ends ~146 ms, visual spot off ~190 ms)
S.padMs   = 100;      % extra ms read on both sides only so the filter settles; spikes there are discarded
S.guardMs = 5;        % the window never reaches closer than this to the NEXT trigger

% ---------------- channels ----------------
S.badChannels = [];   % AP channels (1-based) to skip entirely (they get no spikes), e.g. [83 376]

% ---------------- artifact blanking (NN values) ----------------
S.blankMode   = 'nn'; % 'nn' = one blank per train (default, as the NN pipeline) | 'perPulse' = one blank per pulse
S.blankPreMs  = 1;    % ms the blank starts BEFORE the (first) pulse
S.recWin      = 46;   % samples (~1.5 ms): sliding window used to test whether the artifact has settled
S.rangeMaxUV  = 150;  % uV: settled when max - min inside that window is at most this
S.bufferSamp  = 10;   % samples added after the settling point
S.capSamp     = 180;  % samples (6 ms): a blank never ends later than this after its (last) pulse
S.guardSpkMs  = 0;    % ms: spikes closer than this to a blank edge are also removed (0 = only inside the blank)

% ---------------- spike detection ----------------
S.bpHz         = [300 6000]; % band-pass corners (Hz)
S.bpOrder      = 3;          % Butterworth order (run forwards and backwards, so the effective order is double)
S.thrCoeff     = 3.5;        % threshold = -thrCoeff x noise (as the old NP script)
S.nNoiseTrials = 100;        % trials (spread over the session) used to measure the noise of each channel
S.noiseEndMs   = 2;          % the noise stretch ends this many ms before the trigger
S.peakLookMs   = 0.4;        % after a crossing, the most negative point within this many ms is the spike time
S.refractoryMs = 0.2;        % two spikes on one channel closer than this: only the larger one is kept
S.snippetMs    = 1.6;        % ms of waveform around the peak used for the checks below (and saved; NN: 1.6 ms)
S.ptpMinUV     = 20;         % uV: reject if the snippet's peak-to-peak is smaller than this (noise)
S.ptpMaxUV     = 500;        % uV: ... or larger than this (artifact)
S.absMaxUV     = 500;        % uV: ... or if any sample of the snippet is beyond +/- this

% ---------------- waveforms ----------------
S.saveWaveforms = true;     % true = also write <short name>.NP_Spikes_waveforms.mat (separate file, can be large)
S.wfScale_uV    = 0.2;       % uV per stored int16 count (0.2 uV resolution, range +/- 6.5 mV)
%% =================================================

if ischar(sessionDir) || isstring(sessionDir), sessionList = {char(sessionDir)}; else, sessionList = sessionDir; end
assert(isempty(ttFile) || numel(sessionList) == 1, 'ttFile can only be set when ONE session is processed.');
assert(any(strcmp(S.blankMode, {'nn', 'perPulse'})), 'blankMode must be ''nn'' or ''perPulse''.');
for iS = 1:numel(sessionList)
    fprintf('\n################ session %d of %d ################\n', iS, numel(sessionList));
    tf = ttFile;
    if isempty(tf), tf = findTT_(sessionList{iS}); end
    detectSession_(tf, outDir, overwrite, nTrialsMax, do_save, S);
end

%% ======================================================================
%% ============================ MAIN FUNCTION ===========================
%% ======================================================================
function detectSession_(ttFile, outDir, overwrite, nTrialsMax, do_save, S)
codeVersion = 'NP_2_SpikeDetect_v1 (2026-11-12)';
tStart = tic;

%% ---- trial table, AP file -------------------------------------------------------------------------
L  = load(ttFile, 'TT');  TT = L.TT;  clear L
fs = TT.info.fsAP;  nCh = TT.info.nChAP;  bv = TT.info.uVperBit;
apFile = TT.info.apFile;
if ~isfile(apFile)                                    % drive renamed / moved: look next to the trial table
    rel = strrep(apFile, TT.info.recDir, '');
    apFile = [fileparts(ttFile) rel];
end
assert(isfile(apFile), 'AP file not found: %s', apFile);
dA = dir(apFile);  nS = floor(dA.bytes/(2*nCh));
assert(nS == TT.info.nSampAP, 'The AP file has %d samples but the trial table says %d: wrong file?', nS, TT.info.nSampAP);
mm = memmapfile(apFile, 'Format', {'int16', [nCh nS], 'x'});

tr = TT.trial;
nTr = min(numel(tr.trigAP), nTrialsMax);
goodCh = true(nCh, 1);  goodCh(S.badChannels(S.badChannels >= 1 & S.badChannels <= nCh)) = false;
gc = find(goodCh);

[bF, aF] = butter(S.bpOrder, S.bpHz/(fs/2), 'bandpass');
cfg = S;  cfg.fs = fs;  cfg.nCh = nCh;  cfg.bv = bv;  cfg.bF = bF;  cfg.aF = aF;
preS  = round(S.preMs*fs/1e3);  postS = round(S.postMs*fs/1e3);
Lk    = round(S.peakLookMs*fs/1e3);  refr = round(S.refractoryMs*fs/1e3);
half  = round(S.snippetMs/2*fs/1e3);  gS = round(S.guardSpkMs*fs/1e3);

fprintf('\n%s\nsession  : %s (%s)\ntrials   : %d of %d | channels: %d (%d skipped) | blank mode: %s\n', ...
        codeVersion, TT.info.shortName, TT.info.session, nTr, numel(tr.trigAP), nCh, nCh - numel(gc), S.blankMode);

%% ---- pass 1: noise per channel (pre-trigger, band-passed) ----------------------------------------
nn = unique(round(linspace(1, nTr, min(S.nNoiseTrials, nTr))));
noiseTr = nan(nCh, numel(nn));
fprintf('pass 1: noise from %d trials ...\n', numel(nn));
for q = 1:numel(nn)
    k = nn(q);
    R = prepTrial_(mm, nS, tr.trigAP(k), tr.nextTrigAP(k), [], cfg);
    c1 = max(1, R.iT - preS);  c2 = R.iT - round(S.noiseEndMs*fs/1e3);
    if c2 > c1 + 100, noiseTr(:, q) = median(abs(R.Xf(:, c1:c2)), 2)/0.6745; end
end
noise_uV = median(noiseTr, 2, 'omitnan');
thr_uV = -S.thrCoeff*noise_uV;
fprintf('noise    : median %.1f uV (range %.1f - %.1f over the processed channels) -> threshold median %.1f uV\n', ...
        median(noise_uV(gc)), min(noise_uV(gc)), max(noise_uV(gc)), median(thr_uV(gc)));

%% ---- pass 2: blank, filter, detect --------------------------------------------------------------------
spkRel   = cell(nCh, nTr);                 % spike samples relative to the trigger (int32)
wfAll    = cell(nCh, nTr);                 % waveforms (int16), only if saveWaveforms
blank_ms = cell(nCh, nTr);
blankStart_ms = nan(nCh, nTr);  blankEnd_ms = nan(nCh, nTr);  blankCapped = false(nCh, nTr);
fprintf('pass 2: %d trials ...\n', nTr);
for k = 1:nTr
    pS = round(tr.pulseMs{k}(:)'*fs/1e3);              % pulse onsets, samples after the trigger
    R  = prepTrial_(mm, nS, tr.trigAP(k), tr.nextTrigAP(k), pS, cfg);
    iT = R.iT;  Xf = R.Xf;  nW = size(Xf, 2);
    w1 = max(1, iT - preS);  w2 = min(nW, iT + postS - 1);
    for c = 1:nCh
        iv = R.blank{c};                                   % [start end] in window samples
        if ~isempty(iv)
            blank_ms{c,k} = (iv - iT)/fs*1e3;
            blankStart_ms(c,k) = blank_ms{c,k}(1,1);  blankEnd_ms(c,k) = blank_ms{c,k}(end,2);
        end
    end
    blankCapped(:,k) = R.capped;
    for c = gc'
        pk = detectChannel_(Xf(c,:), thr_uV(c), Lk, refr);
        pk = pk(pk >= w1 & pk <= w2 & pk > half & pk <= nW - half);
        iv = R.blank{c};
        for b = 1:size(iv, 1)
            pk = pk(~(pk >= iv(b,1) - gS & pk <= iv(b,2) + gS));
        end
        if isempty(pk), continue; end
        xc = Xf(c, :);
        W = xc(bsxfun(@plus, pk(:), -half:half));           % one row per spike
        ptp = max(W, [], 2) - min(W, [], 2);  am = max(abs(W), [], 2);
        ok = ptp >= S.ptpMinUV & ptp <= S.ptpMaxUV & am <= S.absMaxUV;
        pk = pk(ok);
        spkRel{c,k} = int32(pk(:) - iT);
        if S.saveWaveforms && any(ok), wfAll{c,k} = int16(round(W(ok,:)/S.wfScale_uV)); end
    end
    if mod(k, 50) == 0 || k == nTr
        fprintf('  trial %d/%d (%.0f s)\n', k, nTr, toc(tStart));
    end
end

%% ---- output variables (NN layout) -----------------------------------------------------------------
nSpk = cellfun(@numel, spkRel);
sp_clipped = cell(1, nCh);  sp_trial = cell(1, nCh);  sp_waveforms = cell(1, nCh);
t0 = double(tr.trigAP(1:nTr)) - 1;                    % trigger, samples from the start of the AP file (0-based)
for c = 1:nCh
    a = cell(nTr, 1);  b = cell(nTr, 1);
    for k = 1:nTr
        a{k} = (t0(k) + double(spkRel{c,k}))/fs*1e3;
        b{k} = k*ones(numel(spkRel{c,k}), 1);
    end
    sp_clipped{c} = vertcat(a{:});  sp_trial{c} = vertcat(b{:});
    if S.saveWaveforms, sp_waveforms{c} = vertcat(wfAll{c,:}); end
end
FS = fs;  preMs = S.preMs;  postMs = S.postMs;  settings = S;  version = codeVersion;  date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
nTrialsProcessed = nTr;

%% ---- printout -------------------------------------------------------------------------------------------
durS = (S.preMs + S.postMs)/1e3;
rate = sum(nSpk, 2)/max(1, nTr)/durS;                % spikes/s over the whole window, per channel
fprintf('\n============================================================\n');
fprintf('NP SPIKE DETECTION (%s)\n', codeVersion);
fprintf('============================================================\n');
fprintf('Session   : %s\nTrial tab.: %s\n', TT.info.shortName, ttFile);
fprintf('Window    : -%g ... +%g ms around the trigger, %d trials\n', S.preMs, S.postMs, nTr);
fprintf('Detection : band-pass %g-%g Hz, threshold -%.1f x noise, %d channels\n', S.bpHz, S.thrCoeff, numel(gc));
fprintf('Spikes    : %d in total; rate per channel (whole window) median %.1f /s, 90th pct %.1f /s, max %.1f /s (ch %d)\n', ...
        sum(nSpk(:)), median(rate(gc)), pct_(rate(gc), 90), max(rate(gc)), gc(find(rate(gc) == max(rate(gc)), 1)));
o = TT.chan.orderTipToTop;  o = o(ismember(o, gc));
nb = 8;  edges = round(linspace(0, numel(o), nb+1));
fprintf('            rate along the probe (tip -> top, %d blocks of channels):\n', nb);
for q = 1:nb
    ch_ = o(edges(q)+1:edges(q+1));
    fprintf('              y %4d-%4d um : median %.1f /s\n', min(TT.chan.y_um(ch_)), max(TT.chan.y_um(ch_)), median(rate(ch_)));
end
fprintf('Blanking  : mode "%s"; blank length after the LAST pulse per condition (median over channels and trials):\n', S.blankMode);
condK = tr.cond(1:nTr);
for c = unique(condK(:))'
    kk = find(condK == c);
    if isempty(tr.pulseMs{kk(1)}), continue; end
    lastP = tr.lastPulseMs(kk)';
    after = bsxfun(@minus, blankEnd_ms(gc, kk), lastP);
    capP  = 100*mean(reshape(blankCapped(gc, kk), [], 1));
    fprintf('            cond %2d %-26s %5.2f ms after the last pulse (90th pct %.2f), capped %4.1f %%\n', ...
            c, TT.cond(c).label, median(after(:), 'omitnan'), pct_(after(:), 90), capP);
end
fprintf('Done in %.0f s.\n', toc(tStart));

%% ---- save -------------------------------------------------------------------------------------------------
if isempty(outDir), od = fileparts(ttFile); else, od = outDir; end
spFile = fullfile(od, [TT.info.shortName '.NP_Spikes.mat']);
wfFile = fullfile(od, [TT.info.shortName '.NP_Spikes_waveforms.mat']);
fprintf('Output    : %s\n', spFile);
if S.saveWaveforms, fprintf('            %s\n', wfFile); end
if ~do_save, fprintf('do_save = false: nothing saved.\n'); return; end
if nTr < numel(tr.trigAP), warning('NPSD:partial', 'Only %d of %d trials were processed (nTrialsMax).', nTr, numel(tr.trigAP)); end
if ~isfolder(od), mkdir(od); end
outs = {spFile};  if S.saveWaveforms, outs{end+1} = wfFile; end
for q = 1:numel(outs)
    if isfile(outs{q})
        if ~overwrite, error('NPSD:exists', '%s already exists and overwrite = false. Nothing was saved.', outs{q}); end
        bak = strrep(outs{q}, '.mat', ['_backup_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']);
        copyfile(outs{q}, bak);  fprintf('Old file copied to %s\n', bak);
    end
end
save(spFile, 'sp_clipped', 'sp_trial', 'nSpk', 'blank_ms', 'blankStart_ms', 'blankEnd_ms', 'blankCapped', ...
     'noise_uV', 'thr_uV', 'goodCh', 'FS', 'preMs', 'postMs', 'nTrialsProcessed', 'ttFile', 'settings', 'version', 'date', '-v7');
if S.saveWaveforms
    wfScale_uV = S.wfScale_uV;  wfTime_ms = (-half:half)/fs*1e3; %#ok<NASGU>
    save(wfFile, 'sp_waveforms', 'wfScale_uV', 'wfTime_ms', 'FS', 'ttFile', 'version', 'date', '-v7.3');
end
fprintf('Saved.\n');
end

%% ======================================================================
%% ============================ HELPERS =================================
%% ======================================================================
function tf = findTT_(sessionDir)
p = regexprep(char(sessionDir), '[\\/]+$', '');
assert(isfolder(p), 'Folder not found: %s', p);
c = dir(fullfile(p, '**', '*.NP_TrialTable.mat'));
c = c(~startsWith({c.name}, '._') & ~contains({c.name}, '_backup_'));
assert(~isempty(c), 'No *.NP_TrialTable.mat below %s: run NP_1_TrialTable first.', p);
assert(numel(c) == 1, 'Several trial tables below %s - set ttFile:\n%s', p, strjoin(fullfile({c.folder}, {c.name}), newline));
tf = fullfile(c(1).folder, c(1).name);
end

function R = prepTrial_(mm, nS, trigAP, nextAP, pS, cfg)
% reads one trial window, blanks the artifact (NN rule), band-passes every channel
fs = cfg.fs;  nCh = cfg.nCh;
pre  = round((cfg.preMs + cfg.padMs)*fs/1e3);
post = round((cfg.postMs + cfg.padMs)*fs/1e3);
s1 = max(1, trigAP - pre);  s2 = min(nS, trigAP + post);
s2 = min(s2, nextAP - round(cfg.guardMs*fs/1e3));
X  = double(mm.Data.x(:, s1:s2)) * cfg.bv;
iT = trigAP - s1 + 1;  nW = size(X, 2);
if iT > 20, X = bsxfun(@minus, X, median(X(:, 1:iT-1), 2)); else, X = bsxfun(@minus, X, median(X, 2)); end
R.iT = iT;  R.blank = cell(nCh, 1);  R.capped = false(nCh, 1);
if ~isempty(pS)
    pS = sort(pS);
    if strcmp(cfg.blankMode, 'nn'), grp = {pS};  else, grp = num2cell(pS); end
    preB = round(cfg.blankPreMs*fs/1e3);
    ivAll = cell(nCh, 1);
    for g = 1:numel(grp)
        bS = iT + grp{g}(1) - preB;                    % blank start (same for all channels)
        o0 = iT + grp{g}(end);                         % onset of the (last) pulse of the group
        nK = cfg.capSamp;
        if o0 + nK + cfg.recWin - 1 > nW, nK = nW - o0 - cfg.recWin + 1; end
        if bS <= 6 || nK < 1, continue; end
        seg = X(:, o0 + (1:nK + cfg.recWin - 1));
        rg  = movmax(seg, [0 cfg.recWin-1], 2, 'Endpoints', 'discard') - movmin(seg, [0 cfg.recWin-1], 2, 'Endpoints', 'discard');
        ok  = rg <= cfg.rangeMaxUV;
        [anyOK, kk] = max(double(ok), [], 2);
        anyOK = anyOK > 0;  kk(~anyOK) = nK;
        endOff = min(kk + cfg.bufferSamp, cfg.capSamp);
        bE = min(o0 + endOff, nW - 6);
        R.capped = R.capped | ~anyOK | (kk + cfg.bufferSamp >= cfg.capSamp);
        for c = 1:nCh, ivAll{c}(end+1, :) = [bS, bE(c)]; end
    end
    for c = 1:nCh
        iv = mergeIv_(ivAll{c});
        for b = 1:size(iv, 1)                          % fill: as simpleBlank_PTD_PulseTrain (NN)
            bS = iv(b,1);  e = iv(b,2);
            shift = X(c, e) - X(c, bS);
            X(c, e:end) = X(c, e:end) - shift;
            X(c, bS:e+1) = linspace(X(c, bS), X(c, e+1), e - bS + 2);
        end
        R.blank{c} = iv;
    end
end
R.Xf = filtfilt(cfg.bF, cfg.aF, X')';
end

function v = pct_(x, p)
% percentile without the Statistics Toolbox (linear interpolation, NaN ignored)
x = sort(x(~isnan(x(:))));
if isempty(x), v = NaN; return; end
r = 1 + (numel(x) - 1)*p/100;  lo = floor(r);  hi = ceil(r);
v = x(lo) + (r - lo)*(x(hi) - x(lo));
end

function iv = mergeIv_(iv)
% merges overlapping / touching [start end] rows
if size(iv, 1) <= 1, return; end
iv = sortrows(iv);  out = iv(1, :);
for b = 2:size(iv, 1)
    if iv(b,1) <= out(end,2) + 1, out(end,2) = max(out(end,2), iv(b,2)); else, out(end+1, :) = iv(b, :); end %#ok<AGROW>
end
iv = out;
end

function pk = detectChannel_(x, thr, L, refr)
% negative threshold crossings -> most negative point within L samples; refractory: keep the larger spike
cr = find(x(1:end-1) > thr & x(2:end) <= thr) + 1;
if isempty(cr), pk = zeros(0, 1); return, end
J = min(bsxfun(@plus, cr(:)', (0:L)'), numel(x));
[~, m] = min(x(J), [], 1);
pk = unique(cr(:)' + m - 1);
keep = true(size(pk));  last = 0;
for q = 1:numel(pk)
    if last > 0 && pk(q) - pk(last) < refr
        if x(pk(q)) < x(pk(last)), keep(last) = false; last = q; else, keep(q) = false; end
    else
        last = q;
    end
end
pk = pk(keep)';
end
