%% NP STAGE 1: spike detection on Neuropixels AP data, aligned to the trigger (script)
% Standalone: needs only the Open Ephys recording folder and the *_exp_datafile_*.mat of the
% session. Uses NO result of the RF-mapping scripts. Nothing big is written to disk: the AP file is
% read trial by trial (a window around each trigger), blanked / filtered in memory, and only the
% detected spikes are kept. Nothing in the recording folder is modified.
%
% WHAT IT DOES (follows the NN pipeline: denoiseIntan_sab -> simpleBlank_PTD_PulseTrain_Xia_v1 -> spike detection):
%   1. Triggers: NI analog input AI3 rising edges, converted to the probe clock with the 1 Hz sync TTL (linear fit).
%   2. Experiment file (StimParams / TrialParams): for every trial the stimulation pulse times
%      (PTD + k * period, only rows with amplitude > 0). Zero-current trials are not blanked.
%   3. Artifact blanking per trial and channel (blankMode 'nn' = same as simpleBlank_PTD_PulseTrain_Xia_v1):
%        start = first pulse - 1 ms; end = after the LAST pulse, a 46-sample window slides until its
%        peak-to-peak range is <= 150 uV, + 10 samples, capped at 180 samples (6 ms).
%        The blank is filled with a straight line and the DC step after it is removed.
%   4. Zero-phase band-pass 300-6000 Hz of each channel on its own (no re-referencing across channels).
%   5. Threshold = -thrCoeff x noise (median|x|/0.6745 of the pre-trigger signal, median over trials);
%      the negative peak within 0.4 ms of each crossing is the spike time; refractory period; snippet checks
%      (peak-to-peak and absolute limits). Spikes inside a trial's blank (plus a guard) are removed.
%   6. Photodiode: a trial is "visual" if regular flashes start after the trigger and were NOT already
%      there before it (some sessions flicker at 120 Hz all the time -> then nobody is flagged).
%   7. NO trials are rejected and no BadTrials file is written: every trial is kept. Per-trial spike
%      count / noise are only plotted (diagnostic figure).
%   8. Saves <name>.sp_np.mat (same variables as .sp_xia.mat + NP extras) and <name>.sp_np_check.mat
%      (raw / blanked / band-passed traces of a few trials and channels: look at them!).
%      Figures: blank diagnostics, per-trial diagnostics, depth x time map, raster + PSTH per channel,
%      blank check per condition.
%
% Written without access to MATLAB. First run: set nTrialsMax = 12 and check the printout. Then run all trials.
% Needs the Signal Processing Toolbox (butter, filtfilt) and MATLAB R2018b or newer (sgtitle, xline, yline).

clear; clc; close all;
codeVersion = 'np_spike_stage1 v1.6 (2026-11-01)';

%% ===================== SETTINGS: EDIT HERE ==========================================
% ---------------- what to process and where to save --------------------------------------
sessionDir = '/Volumes/Shared/Xia/DATA/CJ282/Xia/CJ282_Xia_Exp1_SimSeq4';  % Open Ephys session folder (or its recording1 folder)
expFile    = '';          % experiment .mat with StimParams/TrialParams; '' = the only *exp_datafile*.mat in the recording folder
outDir     = '';          % where the output files are saved; '' = the recording folder
outName    = '';          % start of the output file names; '' = name of the session folder
nTrialsMax = Inf;         % process only the first N trials (e.g. 12 for a quick test); Inf = all trials

% ---------------- streams and trigger lines --------------------------------------------
apMatch   = '-AP';        % text that identifies the AP (30 kHz spike band) stream folder
niMatch   = 'NI-DAQmx';   % text that identifies the NI-DAQ stream folder (trigger and photodiode)
trigAI    = 3;            % NI analog input (0-based) that carries the trial trigger
diodeAI   = 0;            % NI analog input (0-based) that carries the photodiode
syncState = 1;            % TTL line (state number) of the 1 Hz sync pulse, present on both NI and probe streams
allowMismatch = false;    % false = stop with an error if number of triggers ~= number of trials in the experiment file
badChannels   = [];       % AP channels (1-based) to skip entirely, e.g. [83 376]; they get no spikes

% ---------------- trial window (ms, relative to the trigger) -----------------------------
preMs   = 200;  % ms BEFORE the trigger kept per trial (raster, noise estimate, baseline)
postMs  = 300;  % ms AFTER the trigger kept per trial (spikes are only detected up to here)
padMs   = 100;  % extra ms read before and after the window only so the filter has settled; spikes in the pad are discarded
guardMs = 5;    % the data read for a trial stops this many ms before the next trigger (never includes the next artifact)

% ---------------- artifact blanking (NN values) ------------------------------------------
blankMode     = 'nn';  % 'nn' = exactly the NN rule (fixed range limit, straight-line fill between the two blank edge samples)
                       % 'adaptive' = my variant (range limit scales with the channel's noise, flat fill at the pre-blank mean)
blankPreMs    = 1;     % ms the blank starts BEFORE the first pulse (NN: 1 ms)
recWin        = 46;    % samples (~1.5 ms): length of the sliding window used to test whether the artifact has recovered
rangeMaxUV    = 150;   % uV: the artifact counts as recovered when max-min inside that window is below this (NN: 150)
rangeNoiseMult= 8;     % only for blankMode 'adaptive': limit = max(rangeMaxUV, this x the channel's own high-frequency noise)
bufferSamp    = 10;    % extra samples added after the recovery point (NN: 10)
capSamp       = 180;   % samples: blank never lasts longer than this after the LAST pulse (NN: 180 = 6 ms)
guardSpkMs    = 0;   % spikes closer than this (ms) to the blank edges are dropped (filter ringing, residual artifact)

% ---------------- spike detection -----------------------------------------------------
bpHz      = [300 6000];  % band-pass corner frequencies (Hz) applied to every channel
bpOrder   = 3;           % Butterworth order of the band-pass; filtering runs forwards and backwards (zero phase), so effectively double
thrCoeff  = 3.5;         % threshold = -thrCoeff x noise, per channel (noise = median|x|/0.6745); 3.5 leaves ~2-3 Hz of noise crossings per channel
nNoiseTrials = 100;       % number of trials (spread over the recording) used to measure each channel's noise in the pre-trigger period
peakLookMs   = 0.4;      % after the signal crosses the threshold, the most negative point within this many ms is taken as the spike time
refractoryMs = 0.2;        % two spikes on the same channel closer than this (ms): only the larger one is kept
snippetMs    = 3;        % length (ms) of the waveform cut around each spike (centred on the peak) for the checks below and for saving
ptpMinUV = 20;           % uV: reject a spike if its snippet's peak-to-peak (max - min) is smaller than this (too small = noise)
ptpMaxUV = 500;          % uV: ... or larger than this (too large = artifact)
absMaxUV = 500;          % uV: ... or if any sample of the snippet goes beyond +/- this
saveWaveforms = false;   % true = also save every spike snippet in <name>.sp_np_waveforms.mat (can be very large)

% ---------------- photodiode: which trials are visual -----------------------------------
diodeMinFlashes   = 5;        % at least this many diode flashes in the 500 ms after the trigger ...
diodeIntervalMs   = [6 18];   % ... with a typical gap between flashes in this range (ms): 120 Hz screen = 8.3 ms
diodeCV           = 0.25;     % ... and a regular rhythm (std of the gaps / mean gap below this)
diodePreMaxEdges  = 3;        % but fewer than this many flashes in the 200 ms BEFORE the trigger (else it is continuous flicker, not a stimulus)
diodeMinVisTrials = 5;        % if fewer trials than this are flagged, the session is treated as having no visual trials

% ---------------- check file and figures ------------------------------------------------
nCheckCh    = 10;       % number of channels (spread evenly over the probe) whose traces are saved in the check file
nChkPerCond = 5;       % number of trials per condition whose traces are saved in the check file
showBlankCheck = true; % draw the blank-check figures: one panel per condition with raw / blanked / band-passed traces, pulse times and blank window
blankCheckCh   = [100:105];   % channel(s) (1-based) for the blank-check figures; [] = the middle one of the check channels
binMs        = 1;      % bin width (ms) of the PSTH and the depth x time map
smoothSigmaMs = 3;    % PSTH only: Gaussian smoothing, standard deviation in ms (0 = no smoothing; kernel spans +/-3 sigma; use a value >= binMs)
psthKeepBlankGap = false; % true = leave the PSTH empty where (almost) no data exist because of the blank; false = let the smoothing bridge those bins
plotChannels = [];     % channels for the raster + PSTH figures (one figure each); [] = the 4 channels with the largest change 10-100 ms after the trigger
showFigures  = true;   % false = no figures at all (only the saved files)
%% ====================================================================================

tStart = tic;
%% ---- recording folder, streams -----------------------------------------------------
recDir = resolveRec_(sessionDir);
oe = jsondecode(fileread(fullfile(recDir,'structure.oebin')));
cont = oe.continuous;  if isstruct(cont), cont = num2cell(cont); end
fol  = cellfun(@(c) c.folder_name, cont, 'UniformOutput', false);
iAP  = find(contains(fol, apMatch), 1);   iNI = find(contains(fol, niMatch), 1);
assert(~isempty(iAP) && ~isempty(iNI), 'AP or NI stream not found in structure.oebin: %s', strjoin(fol, ', '));
apDir = stripSep_(fol{iAP});  niDir = stripSep_(fol{iNI});
fsAP = double(cont{iAP}.sample_rate);  nCh = double(cont{iAP}.num_channels);
fsNI = double(cont{iNI}.sample_rate);  nNI = double(cont{iNI}.num_channels);
chI  = cont{iAP}.channels;  if iscell(chI), chI = chI{1}; else, chI = chI(1); end
bv   = double(chI.bit_volts);
fprintf('\n%s\nrecording: %s\nAP: %s, %g Hz, %d ch, %.4f uV/bit | NI: %s, %g Hz, %d ch\n', ...
        codeVersion, recDir, apDir, fsAP, nCh, bv, niDir, fsNI, nNI);

apFile = fullfile(recDir,'continuous',apDir,'continuous.dat');
niFile = fullfile(recDir,'continuous',niDir,'continuous.dat');
dA = dir(apFile);  nS  = floor(dA.bytes/(2*nCh));
dN = dir(niFile);  nSN = floor(dN.bytes/(2*nNI));
mm   = memmapfile(apFile, 'Format', {'int16',[nCh nS],'x'});
mmNI = memmapfile(niFile, 'Format', {'int16',[nNI nSN],'x'});
fprintf('AP %.1f s, NI %.1f s\n', nS/fsAP, nSN/fsNI);

if isempty(outDir), outDir = recDir; end
if isempty(outName), [~,outName] = fileparts(regexprep(sessionDir,'[\\/]+$','')); end
if ~isfolder(outDir), mkdir(outDir); end

%% ---- experiment file ---------------------------------------------------------------
if isempty(expFile)
    f = dir(fullfile(recDir,'*exp_datafile*.mat'));  f = f(~startsWith({f.name},'._'));
    assert(numel(f) == 1, 'Found %d *exp_datafile*.mat in %s: set expFile.', numel(f), recDir);
    expFile = fullfile(f(1).folder, f(1).name);
end
Texp = loadExp_(expFile);
fprintf('experiment file: %s\n  %d trials, %d electrode rows per trial, %d conditions, %d zero-current trials\n', ...
        expFile, Texp.nTrials, Texp.nElec, numel(unique(Texp.condID)), sum(Texp.isZero));

%% ---- triggers (NI AI3) -> probe clock ---------------------------------------------------
fprintf('\nreading NI analog inputs ...\n');
aiRows = double(mmNI.Data.x([trigAI+1, diodeAI+1], :));
aiTrig = aiRows(1,:);  aiDiode = aiRows(2,:);  clear aiRows
lim = prctile(aiTrig(1:5:end), [1 99]);
hiT = aiTrig > mean(lim);
upN = find(diff(hiT) == 1) + 1;                          % NI sample index of each rising edge
fprintf('triggers on AI%d: %d (width %.0f ms median)\n', trigAI, numel(upN), ...
        1e3*median(local_widths_(hiT, upN))/fsNI);
P = syncFit_(recDir, niDir, apDir, fsNI, fsAP, syncState);
trigAPs = polyval(P, (upN-1)/fsNI);                       % seconds on the AP clock
trigS   = round(trigAPs*fsAP) + 1;                        % AP sample index (1-based)
nTrig = numel(trigS);
if nTrig ~= Texp.nTrials
    msg = sprintf('%d triggers but %d trials in the experiment file.', nTrig, Texp.nTrials);
    if allowMismatch, warning('%s Using the first %d.', msg, min(nTrig,Texp.nTrials)); else, error('%s', msg); end
end
nTr = min([nTrig, Texp.nTrials, nTrialsMax]);
trigS = trigS(1:nTr);  upN = upN(1:nTr);
fields = fieldnames(Texp);
for q = 1:numel(fields)
    v = Texp.(fields{q});
    if (iscell(v) || isnumeric(v) || islogical(v)) && numel(v) == Texp.nTrials && ~strcmp(fields{q},'nTrials')
        Texp.(fields{q}) = v(1:nTr);
    end
end
sp = diff(double(trigS))/fsAP*1e3;
if any(sp < postMs + padMs), warning('%d trial(s) start closer than %g ms after the previous one: their window is cut at the next trigger.', sum(sp < postMs+padMs), postMs+padMs); end

%% ---- photodiode flag ---------------------------------------------------------------------
[isVis, nFlashPost, nFlashPre] = diodeFlags_(aiDiode, upN, fsNI, diodeMinFlashes, diodeIntervalMs, diodeCV, diodePreMaxEdges);
clear aiDiode aiTrig
if mean(nFlashPre >= diodePreMaxEdges) > 0.5
    fprintf('photodiode: flashes were already present BEFORE the trigger in %.0f%% of the trials (continuous flicker): visual flag = 0 for all.\n', 100*mean(nFlashPre >= diodePreMaxEdges));
end
if sum(isVis) < diodeMinVisTrials && any(isVis)
    fprintf('photodiode: only %d trial(s) flagged (< %d): treated as noise, visual flag = 0 for all.\n', sum(isVis), diodeMinVisTrials);
    isVis(:) = false;
end
fprintf('photodiode: %d of %d trials flagged visual\n', sum(isVis), nTr);

%% ---- probe geometry ------------------------------------------------------------------------
[ypos, xpos] = readProbePositions_(sessionDir);
if isempty(ypos) || numel(ypos) < nCh
    warning('No electrode positions in settings.xml: depth = channel number.');
    ypos = 1:nCh;  xpos = zeros(1,nCh);
end
depth_um = ypos(1:nCh);  x_um = xpos(1:nCh);
goodMask = true(nCh,1);  goodMask(badChannels) = false;

%% ---- configuration struct for the helper functions -------------------------------------------
cfg = struct('fs',fsAP,'nCh',nCh,'bv',bv,'preMs',preMs,'postMs',postMs,'padMs',padMs,'guardMs',guardMs, ...
    'blankPreMs',blankPreMs,'recWin',recWin,'rangeMaxUV',rangeMaxUV,'rangeNoiseMult',rangeNoiseMult,'blankMode',blankMode, ...
    'bufferSamp',bufferSamp,'capSamp',capSamp,'goodMask',goodMask,'chkCh',[]);
[bF, aF] = butter(bpOrder, bpHz/(fsAP/2), 'bandpass');
cfg.bF = bF;  cfg.aF = aF;
preSamp  = round(preMs*fsAP/1e3);  postSamp = round(postMs*fsAP/1e3);

%% ---- pass 1: noise per channel -----------------------------------------------------------------
nn = unique(round(linspace(1, nTr, min(nNoiseTrials, nTr))));
noiseTr = nan(nCh, numel(nn));
fprintf('\npass 1: noise from %d trials ...\n', numel(nn));
for k = 1:numel(nn)
    t = nn(k);
    R = prepTrial_(mm, nS, trigS(t), nextTrig_(trigS,t), Texp.events{t}, cfg, false);
    c2 = R.iT - round(2*fsAP/1e3);  c1 = max(1, R.iT - preSamp);
    if c2 > c1 + 50, noiseTr(:,k) = median(abs(R.Xf(:, c1:c2)), 2)/0.6745; end
end
noise_uV = median(noiseTr, 2, 'omitnan');
thr_uV = -thrCoeff*noise_uV;
fprintf('noise (MAD/0.6745) of the band-passed signal: median %.1f uV, range %.1f to %.1f; threshold median %.1f uV\n', ...
        median(noise_uV(goodMask)), min(noise_uV(goodMask)), max(noise_uV(goodMask)), median(thr_uV(goodMask)));

%% ---- check-file selection -------------------------------------------------------------------------
gc = find(goodMask);  chkCh = gc(unique(round(linspace(1, numel(gc), nCheckCh))));
if isempty(blankCheckCh), blankCheckCh = chkCh(ceil(numel(chkCh)/2)); end
chkCh = unique([chkCh(:)' blankCheckCh(:)']);
cfg.chkCh = chkCh;
ids = unique(Texp.condID);  chkTrials = [];
for q = 1:numel(ids)
    tt = find(Texp.condID == ids(q));
    chkTrials = [chkTrials tt(1:min(nChkPerCond, numel(tt)))]; %#ok<AGROW>
end
chkTrials = sort(chkTrials);
chk = struct('trial',{},'cond',{},'t_ms',{},'raw',{},'blanked',{},'filtered',{},'spike_rel_ms',{});

%% ---- pass 2: all trials ---------------------------------------------------------------------------
spk = cell(nCh, nTr);                       % relative sample (to the trigger), int32
blkS = nan(nCh, nTr, 'single');  blkE = nan(nCh, nTr, 'single');   % blank start/end, samples after trigger
capped = false(nCh, nTr);
trialMax = nan(nCh, nTr);  trialNoise = nan(nCh, nTr);
snips = cell(nCh, 1);  snipT = cell(nCh, 1);
half = round(snippetMs/2*fsAP/1e3);  L = round(peakLookMs*fsAP/1e3);  refr = round(refractoryMs*fsAP/1e3);
gSpk = round(guardSpkMs*fsAP/1e3);
fprintf('\npass 2: %d trials ...\n', nTr);
for t = 1:nTr
    isChk = any(chkTrials == t);
    R = prepTrial_(mm, nS, trigS(t), nextTrig_(trigS,t), Texp.events{t}, cfg, isChk);
    iT = R.iT;  Xf = R.Xf;  nW = size(Xf,2);
    w1 = max(1, iT - preSamp);  w2 = min(nW, iT + postSamp - 1);
    blkS(:,t) = R.blankS - iT;  blkE(:,t) = R.blankE - iT;  capped(:,t) = R.capped;
    c2 = iT - round(2*fsAP/1e3);
    if c2 > w1 + 50, trialNoise(:,t) = median(abs(Xf(:, w1:c2)), 2)/0.6745; end
    if any(~isnan(R.blankE)), s0 = max(R.blankE) + round(5e-3*fsAP); else, s0 = iT + round(5e-3*fsAP); end
    if s0 < w2, trialMax(:,t) = max(abs(Xf(:, s0:w2)), [], 2); end
    for c = gc'
        pk = detectChannel_(Xf(c,:), thr_uV(c), L, refr);
        pk = pk(pk >= w1 & pk <= w2);
        if ~isnan(R.blankS(c)), pk = pk(~(pk >= R.blankS(c) - gSpk & pk <= R.blankE(c) + gSpk)); end
        pk = pk(pk > half & pk <= nW - half);
        if isempty(pk), continue; end
        S = Xf(c, pk + (-half:half));  if isvector(S) && numel(pk) == 1, S = S(:)'; end
        ptp = max(S,[],2) - min(S,[],2);  am = max(abs(S),[],2);
        ok = ptp >= ptpMinUV & ptp <= ptpMaxUV & am <= absMaxUV;
        pk = pk(ok);  S = S(ok,:);
        spk{c,t} = int32(pk(:) - iT);
        if saveWaveforms && ~isempty(pk)
            snips{c} = [snips{c}; int16(round(S*5))];  snipT{c} = [snipT{c}; double(pk(:)-iT)/fsAP*1e3 + (trigS(t)-1)/fsAP*1e3];
        end
    end
    if isChk
        sp_rel = cell(numel(chkCh),1);
        for q = 1:numel(chkCh), sp_rel{q} = double(spk{chkCh(q),t})/fsAP*1e3; end
        cols = R.chkCols;
        chk(end+1) = struct('trial',t,'cond',Texp.condID(t),'t_ms',((cols - iT)/fsAP*1e3), ...
            'raw',single(R.rawChk),'blanked',single(R.blkChk),'filtered',single(Xf(chkCh, cols)),'spike_rel_ms',{{sp_rel}}); %#ok<SAGROW>
    end
    if mod(t,25) == 0 || t == nTr
        fprintf('  trial %d/%d (%.0f s elapsed)\n', t, nTr, toc(tStart));
    end
end

%% ---- spike counts (no trial is rejected) ---------------------------------------------------------------
nSpk = cellfun(@numel, spk);                  % nCh x nTr
trial_bad = false(1, nTr);                    % kept for format compatibility: always all false
cnt = sum(nSpk(goodMask,:), 1);
tr_idx = 1:nTr;  good = tr_idx;

%% ---- build the output variables (same names as sp_xia.mat) -------------------------------------------
sp_clipped = cell(1, nCh);
AllSpikes = struct('rel_spike_times',[],'trial_numbers',[],'channel',[],'probe',[],'coeff',[], ...
                   'threshold_uV',[],'noise_uV',[],'good_trials',[],'n_good_trials',[]);
for c = 1:nCh
    relA = cell(nTr,1);  trl = cell(nTr,1);
    for t = 1:nTr
        r = double(spk{c,t});
        relA{t} = r/fsAP*1e3;  trl{t} = t*ones(numel(r),1);
    end
    absMs = cell(nTr,1);
    for t = 1:nTr, absMs{t} = (double(trigS(t)-1) + double(spk{c,t}))/fsAP*1e3; end
    sp_clipped{c} = vertcat(absMs{:});
    rg = relA(good);  tg = trl(good);
    AllSpikes(c).rel_spike_times = vertcat(rg{:});
    AllSpikes(c).trial_numbers   = vertcat(tg{:});
    AllSpikes(c).channel = c;  AllSpikes(c).probe = 1;  AllSpikes(c).coeff = thrCoeff;
    AllSpikes(c).threshold_uV = thr_uV(c);  AllSpikes(c).noise_uV = noise_uV(c);
    AllSpikes(c).good_trials = good;  AllSpikes(c).n_good_trials = numel(good);
end
fs = fsAP;  pre_ms = preMs;  post_ms = postMs;  stim_dur_ms = NaN;
good_trials_probe = {good};  trial_bad_probe = trial_bad;
probe_names = {'Neuropixels'};  probe_elecs = {1:nCh};  hs_names = {'Neuropixels'};  hs_elecs = {1:nCh};
coeff_ch = thrCoeff*ones(1,nCh);  thr_ch = thr_uV(:)';  threshold_method = 'mad';
np = struct();
np.code = codeVersion;  np.sessionDir = sessionDir;  np.recDir = recDir;  np.expFile = expFile;
np.depth_um = depth_um;  np.x_um = x_um;  np.badChannels = badChannels;
np.trig_sample = trigS(:)';  np.trig_ms = (double(trigS(:)')-1)/fsAP*1e3;
np.cond_id = Texp.condID;  np.cond_label = Texp.label;  np.event_times_us = Texp.events;
np.is_zero_current = Texp.isZero;  np.amp_max_uA = Texp.ampMax;  np.elec_names = Texp.elecNames;
np.is_visual = isVis(:)';  np.n_flash_post = nFlashPost(:)';  np.n_flash_pre = nFlashPre(:)';
np.blank_start_rel_samples = blkS;  np.blank_end_rel_samples = blkE;  np.blank_capped = capped;
np.n_spikes = nSpk;  np.trial_noise_uV = trialNoise;  np.trial_max_uV = trialMax;
np.sync_fit = P;
np.settings = struct('preMs',preMs,'postMs',postMs,'padMs',padMs,'blankPreMs',blankPreMs,'recWin',recWin, ...
    'rangeMaxUV',rangeMaxUV,'rangeNoiseMult',rangeNoiseMult,'blankMode',blankMode,'bufferSamp',bufferSamp,'capSamp',capSamp, ...
    'bpHz',bpHz,'bpOrder',bpOrder,'thrCoeff',thrCoeff,'peakLookMs',peakLookMs, ...
    'refractoryMs',refractoryMs,'snippetMs',snippetMs,'ptpMinUV',ptpMinUV,'ptpMaxUV',ptpMaxUV,'absMaxUV',absMaxUV);

spFile = fullfile(outDir, [outName '.sp_np.mat']);
save(spFile, 'sp_clipped','AllSpikes','fs','pre_ms','post_ms','stim_dur_ms','good_trials_probe','trial_bad_probe', ...
     'tr_idx','probe_names','probe_elecs','hs_names','hs_elecs','coeff_ch','thr_ch', ...
     'threshold_method','np','-v7.3');
save(fullfile(outDir,[outName '.sp_np_check.mat']), 'chk','chkCh','chkTrials','depth_um','fsAP','-v7.3');
if saveWaveforms
    sp_waveforms = cell(1,nCh);
    for c = 1:nCh, sp_waveforms{c} = {snipT{c}, snips{c}}; end %#ok<NASGU>
    snippet_ms = snippetMs;  waveform_scale_uV_per_count = 1/5; %#ok<NASGU>
    save(fullfile(outDir,[outName '.sp_np_waveforms.mat']), 'sp_waveforms','fs','snippet_ms','waveform_scale_uV_per_count','-v7.3');
end
fprintf('\nsaved:\n  %s\n  %s\n', spFile, fullfile(outDir,[outName '.sp_np_check.mat']));

%% ---- printout ------------------------------------------------------------------------------------------
rate = sum(nSpk(:, ~trial_bad), 2) / max(1, sum(~trial_bad)) / ((preMs+postMs)/1e3);
fprintf('\n--- summary ---\n');
fprintf('spikes: %d in total; rate per channel (good trials) median %.1f Hz, range %.1f to %.1f Hz\n', ...
        sum(nSpk(:)), median(rate(goodMask)), min(rate(goodMask)), max(rate(goodMask)));
isEl = ~Texp.isZero(:)';
bl = double(blkE - blkS)/fsAP*1e3;
fprintf('\nblank length per condition (good channels x trials), mode "%s":\n', blankMode);
for id = unique(Texp.condID)
    tt = find(Texp.condID == id);
    if all(Texp.isZero(tt)), fprintf('  cond %2d: %3d trials, no stimulation (no blank)\n', id, numel(tt)); continue; end
    bb = bl(goodMask, tt);  bb = bb(~isnan(bb));
    fprintf('  cond %2d: %3d trials, blank median %.1f ms, min %.1f, max %.1f, at cap %.1f%%   [%s]\n', id, numel(tt), ...
            median(bb), min(bb), max(bb), 100*mean(capped(goodMask, tt), 'all'), Texp.label{tt(1)});
end
if any(isEl)
    b = bl(goodMask, isEl);  b = b(~isnan(b));
    fprintf('blank length (electrical trials): median %.1f ms, 90%% %.1f ms, max %.1f ms; hit the cap in %.1f%% of channel-trials\n', ...
            median(b), prctile(b,90), max(b), 100*mean(capped(goodMask, isEl), 'all'));
end
fprintf('done in %.0f s\n', toc(tStart));

%% ---- figures ------------------------------------------------------------------------------------------
if showFigures
    tms = -preMs:binMs:postMs;  nB = numel(tms)-1;  tc = tms(1:end-1) + binMs/2;
    [~, ord] = sortrows([depth_um(:) (1:nCh)']);   ord = ord(:)';
    % ---- Fig 1: blanking
    figure('Color','w','Position',[40 60 1250 420],'Name','blanking','NumberTitle','off');
    subplot(1,3,1);
    if any(isEl), histogram(b, 40); end
    xlabel('blank length per channel and trial (ms)'); ylabel('count'); grid on;
    title(sprintf('electrical trials, cap %.1f ms reached in %.1f%%', capSamp/fsAP*1e3, 100*mean(capped(goodMask,isEl),'all')));
    subplot(1,3,2);
    imagesc(1:nTr, 1:nCh, double(capped(ord,:))); colormap(gca, [1 1 1; 0.8 0.1 0.1]);
    xlabel('trial'); ylabel('channel (sorted by depth)'); title('red = blank hit the cap');
    subplot(1,3,3);
    if ~isempty(chk)
        q = find(~Texp.isZero(chkTrials), 1);  if isempty(q), q = 1; end
        cc = 1;
        plot(chk(q).t_ms, chk(q).raw(cc,:), 'Color', [0.7 0.7 0.7]); hold on;
        plot(chk(q).t_ms, chk(q).blanked(cc,:), 'b'); plot(chk(q).t_ms, chk(q).filtered(cc,:), 'r');
        xlabel('ms from trigger'); ylabel('\muV'); grid on; legend({'raw','blanked','band-passed'});
        title(sprintf('trial %d (cond %d), channel %d', chk(q).trial, chk(q).cond, chkCh(cc)));
    end
    % ---- Fig 2: per-trial diagnostics (nothing is rejected)
    figure('Color','w','Position',[60 80 900 420],'Name','per-trial diagnostics','NumberTitle','off');
    yyaxis left;  plot(tr_idx, cnt, 'k.'); ylabel('BLACK dots: spikes per trial (all good channels)');
    yyaxis right; plot(tr_idx, median(trialNoise(goodMask,:),1,'omitnan'), '.', 'Color', [0.2 0.4 0.9]); ylabel('BLUE dots: pre-trigger noise (\muV)');
    ax = gca;  ax.YAxis(1).Color = 'k';  ax.YAxis(2).Color = [0.2 0.4 0.9];
    grid on; xlabel('trial'); title(sprintf('%d trials, all kept', nTr));
    % ---- depth x time maps (rate, blank time excluded)
    groups = {find(~trial_bad), 'good trials'};
    if any(isVis & ~trial_bad) && any(~isVis & ~trial_bad)
        groups = {find(isVis & ~trial_bad), 'visual (diode flashes)'; find(~isVis & ~trial_bad), 'no flashes (electrical / blank)'};
    end
    nG = size(groups,1);
    Z = cell(nG,1);  Rt = cell(nG,1);  Rs = cell(nG,1);
    for g = 1:nG
        tt = groups{g,1};
        cntB = zeros(nCh, nB);  valid = zeros(nCh, nB);
        for t = tt(:)'
            for c = gc'
                r = double(spk{c,t})/fsAP*1e3;
                if ~isempty(r), cntB(c,:) = cntB(c,:) + histcounts(r, tms); end
            end
            bs = double(blkS(:,t))/fsAP*1e3;  be = double(blkE(:,t))/fsAP*1e3;
            ov = max(0, min(repmat(tms(2:end),nCh,1), repmat(be,1,nB)) - max(repmat(tms(1:end-1),nCh,1), repmat(bs,1,nB)));
            ov(isnan(ov)) = 0;
            valid = valid + (binMs - ov)/1e3;
        end
        rt = cntB ./ max(valid, eps);  rt(valid < 0.2*numel(tt)*binMs/1e3) = NaN;
        % smoothed PSTH: smooth the spike counts and the valid (non-blanked) time separately, then divide
        if smoothSigmaMs > 0
            sb = smoothSigmaMs/binMs;  hw = ceil(3*sb);  kx = -hw:hw;
            ker = exp(-0.5*(kx/sb).^2);  ker = ker/sum(ker);
            cS = conv2(cntB, ker, 'same');  vS = conv2(valid, ker, 'same');
            rs = cS ./ max(vS, eps);
        else
            rs = cntB ./ max(valid, eps);  vS = valid;
        end
        if psthKeepBlankGap
            rs(valid < 0.2*numel(tt)*binMs/1e3) = NaN;
        else
            rs(vS < 0.05*numel(tt)*binMs/1e3) = NaN;
        end
        Rs{g} = rs;  Rs{g}(~goodMask,:) = NaN;
        bb = tc < -5;  mu = mean(rt(:,bb), 2, 'omitnan');  sd = std(rt(:,bb), 0, 2, 'omitnan');
        Z{g} = (rt - mu) ./ max(sd, eps);  Rt{g} = rt;  Rt{g}(~goodMask,:) = NaN;  Z{g}(~goodMask,:) = NaN;
    end
    figure('Color','w','Position',[80 100 520*nG 620],'Name','depth x time','NumberTitle','off');
    for g = 1:nG
        subplot(1,nG,g);
        imagesc(tc, 1:nCh, Z{g}(ord,:), [-4 4]);  set(gca,'YDir','normal');
        colormap(gca, redblue_());  colorbar;
        tk = 1:40:nCh;  set(gca,'YTick',tk,'YTickLabel',arrayfun(@(k) sprintf('%.0f',depth_um(ord(k))), tk,'UniformOutput',false));
        xlabel('ms from trigger'); ylabel('depth from tip (\mum)');
        title(sprintf('%s (n = %d): rate, z vs baseline', groups{g,2}, numel(groups{g,1})));
        xline(0,'k:');
    end
    % ---- raster + PSTH: one figure per channel
    [~, gBig] = max(cellfun(@numel, groups(:,1)));
    if isempty(plotChannels)
        win = tc >= 10 & tc <= 100;
        sc = max(abs(Z{gBig}(:,win)), [], 2, 'omitnan');  sc(~goodMask | isnan(sc)) = -Inf;
        [~, o2] = sort(sc, 'descend');  plotChannels = o2(1:min(4, numel(o2)))';
    end
    [~, tord] = sortrows([double(isVis(:)) double(Texp.condID(:))]);  tord = tord(:)';
    tord = tord(~trial_bad(tord));
    for q = 1:numel(plotChannels)
        c = plotChannels(q);
        figure('Color','w','Position',[100+30*q 60 760 700],'Name',sprintf('raster + PSTH ch %d',c),'NumberTitle','off');
        subplot(5,1,1:3); hold on;
        row = 0;  prevKey = [];  edgesRow = [];  lab = {};  labRow = [];
        for t = tord
            row = row + 1;
            r = double(spk{c,t})/fsAP*1e3;
            plot(r, row*ones(size(r)), 'k.', 'MarkerSize', 4);
            key = [isVis(t) Texp.condID(t)];
            if isempty(prevKey) || ~isequal(key, prevKey)
                if ~isempty(prevKey), yline(row-0.5, 'Color', [0.8 0.3 0.3]); end
                edgesRow(end+1) = row; %#ok<AGROW>
                lab{end+1} = sprintf('cond %d', Texp.condID(t)); %#ok<AGROW>
            end
            prevKey = key;
        end
        set(gca,'YTick',edgesRow,'YTickLabel',lab);
        xlim([-preMs postMs]); ylim([0 row+1]); xline(0,'b:'); set(gca,'YDir','reverse'); box on;
        title(sprintf('ch %d (%.0f \\mum depth), threshold %.1f \\muV', c, depth_um(c), thr_uV(c)));
        ylabel('trial (sorted by condition)');
        subplot(5,1,4:5); hold on;
        for g = 1:nG
            if numel(groups{g,1}) >= 3, plot(tc, Rs{g}(c,:), 'LineWidth', 1.2); end
        end
        xlim([-preMs postMs]); xline(0,'b:'); grid on; xlabel('ms from trigger'); ylabel(sprintf('rate (Hz), Gaussian smooth \\sigma = %g ms', smoothSigmaMs));
        legend(groups(cellfun(@numel,groups(:,1)) >= 3, 2), 'Location','northeast');
    end
end
if showBlankCheck
    for ch = blankCheckCh(:)'
        blankCheckFig_(chk, chkCh, ch, Texp, blkS, blkE, fsAP, thr_uV(ch));
    end
end
fprintf('\nStage 1 finished.\n');

%% ============================ local functions =====================================================
function rec = resolveRec_(p)
p = regexprep(p,'[\\/]+$','');
if exist(fullfile(p,'structure.oebin'),'file'), rec = p; return, end
c = dir(fullfile(p,'**','structure.oebin'));  c = c(~startsWith({c.name},'._'));
assert(~isempty(c), 'No structure.oebin below %s', p);
if numel(c) > 1
    [~, i] = sort({c.folder});  c = c(i);
    warning('Several recordings below %s; using the first: %s', p, c(1).folder);
end
rec = c(1).folder;
end

function s = stripSep_(s)
s = regexprep(char(s),'[\\/]+$','');
end

function n = nextTrig_(trigS, t)
if t < numel(trigS), n = trigS(t+1); else, n = []; end
end

function w = local_widths_(hi, up)
w = nan(size(up));
for k = 1:numel(up)
    j = find(~hi(up(k):min(numel(hi), up(k)+1e5)), 1);
    if ~isempty(j), w(k) = j-1; end
end
w = w(~isnan(w));  if isempty(w), w = NaN; end
end

function T = loadExp_(expFile)
S = load(expFile);
assert(isfield(S,'StimParams') && isfield(S,'TrialParams'), '%s has no StimParams / TrialParams.', expFile);
SP = S.StimParams;  TP = S.TrialParams;
nRows = size(SP,1) - 1;
nTr = [];  if isfield(S,'n_Trials'), nTr = double(S.n_Trials(1)); end
if isempty(nTr) || nTr < 1 || mod(nRows, nTr) ~= 0
    tn = cell2mat(TP(2:end,1));  nElec = sum(tn == tn(1));  nTr = nRows/nElec;
else
    nElec = nRows/nTr;
end
T = struct('nTrials',nTr,'nElec',nElec,'condID',zeros(1,nTr),'label',{cell(1,nTr)},'events',{cell(1,nTr)}, ...
           'isZero',false(1,nTr),'ampMax',zeros(1,nTr),'elecNames',{cell(1,nTr)});
for k = 1:nTr
    rows = 1 + (k-1)*nElec + (1:nElec);
    amp = cell2mat(SP(rows,16));  act = amp > 0;
    ptd = cell2mat(SP(rows,6));  npl = cell2mat(SP(rows,8));  per = cell2mat(SP(rows,9));
    ev = [];
    for r = find(act(:))'
        if npl(r) > 0, ev = [ev, ptd(r) + (0:npl(r)-1)*per(r)]; end %#ok<AGROW>
    end
    T.events{k} = unique(sort(ev));
    T.isZero(k) = ~any(act);
    T.ampMax(k) = max(amp);
    T.condID(k) = double(cell2mat(TP(rows(1),2)));
    T.elecNames{k} = SP(rows(act),1)';
    % label comes from StimParams (what was really delivered); StimMeta can be stale in "changed_electrodes" files
    ia = find(act(:))';  parts = cell(1, numel(ia));
    for ii = 1:numel(ia)
        r = ia(ii);
        nm = regexp(char(SP{rows(r),1}), '\d+', 'match', 'once');   % Intan name 'A-016' -> channel number 16
        parts{ii} = sprintf('Ch%s %gms', nm, ptd(r)/1e3);
    end
    if isempty(parts)
        T.label{k} = sprintf('cond %d: no stim', T.condID(k));
    elseif all(amp(ia) == amp(ia(1)))
        T.label{k} = sprintf('cond %d: %s, %g uA', T.condID(k), strjoin(parts, ' + '), amp(ia(1)));
    else
        for ii = 1:numel(ia), parts{ii} = sprintf('%s %guA', parts{ii}, amp(ia(ii))); end
        T.label{k} = sprintf('cond %d: %s', T.condID(k), strjoin(parts, ' + '));
    end
end
end

function P = syncFit_(rec, niDir, apDir, fsNI, fsAP, state)
d0 = readNPY_(fullfile(rec,'continuous',niDir,'sample_numbers.npy'), 1);
p0 = readNPY_(fullfile(rec,'continuous',apDir,'sample_numbers.npy'), 1);
sdir = fullfile(rec,'events',niDir,'TTL');  pdir = fullfile(rec,'events',apDir,'TTL');
dsn = readNPY_(fullfile(sdir,'sample_numbers.npy'));  dst = readNPY_(fullfile(sdir,'states.npy'));
psn = readNPY_(fullfile(pdir,'sample_numbers.npy'));  pst = readNPY_(fullfile(pdir,'states.npy'));
sd = (dsn(dst == state) - double(d0))/fsNI;
sn = (psn(pst == state) - double(p0))/fsAP;
assert(numel(sd) >= 20 && numel(sn) >= 20, 'Too few sync pulses (NI %d, probe %d).', numel(sd), numel(sn));
j  = interp1(sn, (1:numel(sn))', sd, 'nearest', 'extrap');
dt = sn(j) - sd;
ok = abs(dt - median(dt)) < 0.05;
[~, uj] = unique(j(ok), 'stable');  ia = find(ok);  ia = ia(uj);
assert(numel(ia) >= 20, 'Only %d sync pulses matched.', numel(ia));
xS = sd(ia);  yS = sn(j(ia));
keep = true(size(xS));
for it = 1:5
    Pk = polyfit(xS(keep), yS(keep), 1);
    rk = (yS - polyval(Pk, xS))*1e3;
    s = 1.4826*median(abs(rk - median(rk)));
    nk = abs(rk - median(rk)) <= 5*max(s, eps);
    if isequal(nk, keep), break; end
    keep = nk;
end
P = polyfit(xS(keep), yS(keep), 1);
res = (yS - polyval(P, xS))*1e3;
fprintf('sync: %d/%d pulses matched, slope %+.2f ppm, residual SD %.4f ms\n', numel(ia), numel(sd), (P(1)-1)*1e6, std(res(keep)));
end

function [isVis, nPost, nPre] = diodeFlags_(ai, upN, fs, minFl, ivMs, cvMax, preMax)
lim = prctile(ai(1:7:end), [1 99]);
lo = lim(1) + 0.3*diff(lim);  hi = lim(1) + 0.7*diff(lim);
n = numel(upN);  isVis = false(1,n);  nPost = zeros(1,n);  nPre = zeros(1,n);
for k = 1:n
    a = ai(max(1, upN(k) - round(0.2*fs)) : upN(k)-1);
    b = ai(upN(k) : min(numel(ai), upN(k) + round(0.5*fs)));
    ea = edges_(a, lo, hi);  eb = edges_(b, lo, hi);
    nPre(k) = numel(ea);  nPost(k) = numel(eb);
    if numel(eb) >= minFl
        d = diff(eb)/fs*1e3;
        isVis(k) = median(d) >= ivMs(1) && median(d) <= ivMs(2) && std(d)/mean(d) < cvMax && numel(ea) < preMax;
    end
end
end

function e = edges_(seg, lo, hi)
st = zeros(size(seg));  st(seg > hi) = 1;  st(seg < lo) = -1;
nz = find(st ~= 0);  s2 = st(nz);
e = nz(find(diff(s2) == 2) + 1);
end

function [ypos, xpos] = readProbePositions_(root)
ypos = [];  xpos = [];
f = dir(fullfile(root,'**','settings.xml'));  f = f(~startsWith({f.name},'._'));
for q = 1:numel(f)
    try, txt = fileread(fullfile(f(q).folder, f(q).name)); catch, continue; end
    yTok = regexp(txt,'<ELECTRODE_YPOS\s([^/>]*)/?>','tokens','once');
    xTok = regexp(txt,'<ELECTRODE_XPOS\s([^/>]*)/?>','tokens','once');
    if isempty(yTok), continue; end
    ypos = parseCH_(yTok{1});
    if ~isempty(xTok), xpos = parseCH_(xTok{1}); else, xpos = zeros(size(ypos)); end
    if ~isempty(ypos), return; end
end
end

function v = parseCH_(s)
t = regexp(s,'CH(\d+)\s*=\s*"(-?[\d.]+)"','tokens');
if isempty(t), v = []; return; end
idx = cellfun(@(c) str2double(c{1}), t) + 1;  val = cellfun(@(c) str2double(c{2}), t);
v = nan(1, max(idx));  v(idx) = val;
end

function R = prepTrial_(mm, nS, trigS, nextS, ev_us, cfg, keep)
fs = cfg.fs;  nCh = cfg.nCh;
pre  = round((cfg.preMs + cfg.padMs)*fs/1e3);
post = round((cfg.postMs + cfg.padMs)*fs/1e3);
s1 = max(1, trigS - pre);  s2 = min(nS, trigS + post);
if ~isempty(nextS), s2 = min(s2, nextS - round(cfg.guardMs*fs/1e3)); end
X = double(mm.Data.x(:, s1:s2)) * cfg.bv;
iT = trigS - s1 + 1;
nW = size(X,2);
if iT > 20, X = X - median(X(:, 1:iT-1), 2); else, X = X - median(X, 2); end
R.iT = iT;  R.blankS = nan(nCh,1);  R.blankE = nan(nCh,1);  R.capped = false(nCh,1);
if keep
    cols = max(1, iT - round(5e-3*fs)) : min(nW, iT + round(30e-3*fs));
    R.chkCols = cols;  R.rawChk = X(cfg.chkCh, cols);
end
% ---- artifact blanking (adaptive end, as simpleBlank_PTD_PulseTrain_Xia_v1) -----------------------
if ~isempty(ev_us)
    evS = round(ev_us*fs/1e6);
    bS = iT + min(evS) - round(cfg.blankPreMs*fs/1e3);
    o0 = iT + max(evS);
    nK = cfg.capSamp;
    if o0 + nK + cfg.recWin - 1 > nW, nK = max(1, nW - o0 - cfg.recWin + 1); end
    if bS > 6 && nK >= 1 && o0 + nK + cfg.recWin - 1 <= nW
        seg = X(:, o0 + (1 : nK + cfg.recWin - 1));
        rg = movmax(seg, [0 cfg.recWin-1], 2, 'Endpoints', 'discard') - movmin(seg, [0 cfg.recWin-1], 2, 'Endpoints', 'discard');
        p0 = max(1, iT - round(cfg.preMs*fs/1e3));
        xp = X(:, p0 : max(p0+20, iT-5));
        sig = median(abs(diff(xp, 1, 2)), 2)/0.6745/sqrt(2);
        if strcmp(cfg.blankMode, 'nn')
            rthr = cfg.rangeMaxUV*ones(nCh,1);               % NN: fixed limit
        else
            rthr = max(cfg.rangeMaxUV, cfg.rangeNoiseMult*sig);
        end
        ok = bsxfun(@le, rg, rthr);
        [anyOK, k] = max(double(ok), [], 2);
        anyOK = anyOK > 0;  k(~anyOK) = nK;
        endOff = min(k + cfg.bufferSamp, cfg.capSamp);
        bE = min(o0 + endOff, nW - 6);
        for c = 1:nCh
            e = bE(c);
            if strcmp(cfg.blankMode, 'nn')
                % as simpleBlank_PTD_PulseTrain_Xia_v1: straight line between the two edge samples,
                % subtract the step at the blank end from everything after, re-draw the line to e+1
                shift = X(c, e) - X(c, bS);
                X(c, e:end) = X(c, e:end) - shift;
                X(c, bS:e+1) = linspace(X(c, bS), X(c, e+1), e - bS + 2);
            else
                a = bS - 1;
                preM = mean(X(c, max(1,a-4):a));  postM = mean(X(c, e+1:min(nW,e+5)));
                X(c, e+1:end) = X(c, e+1:end) - (postM - preM);
                X(c, bS:e) = preM;
            end
        end
        R.blankS(:) = bS;  R.blankE = double(bE(:));  R.capped = ~anyOK(:) | (k(:) + cfg.bufferSamp >= cfg.capSamp);
    end
end
if keep, R.blkChk = X(cfg.chkCh, cols); end
% ---- filter -------------------------------------------------------------------------
R.Xf = filtfilt(cfg.bF, cfg.aF, X')';
end

function blankCheckFig_(chk, chkCh, ch, Texp, blkS, blkE, fs, thr)
% One panel per condition: all kept check trials of that condition overlaid, for one channel.
q = find(chkCh == ch, 1);
conds = unique([chk.cond]);  nC = numel(conds);
nCol = min(5, nC);  nRow = ceil(nC/nCol);
for zoom = 0:1
    if zoom
        nm = sprintf('blank check ch %d: ZOOM band-passed', ch);
    else
        nm = sprintf('blank check ch %d: raw / blanked / band-passed', ch);
    end
    figure('Color','w','Position',[20 40 330*nCol 290*nRow+30],'Name',nm,'NumberTitle','off');
    for k = 1:nC
        subplot(nRow, nCol, k); hold on;
        idx = find([chk.cond] == conds(k));
        t1 = chk(idx(1)).trial;
        if zoom, yl = [-200 200]; else, yl = [-1300 1300]; end
        bl = nan(1, numel(idx));
        for ii = 1:numel(idx)
            d = chk(idx(ii));  t = d.trial;
            if ~isnan(blkS(ch,t))
                bl(ii) = double(blkE(ch,t) - blkS(ch,t))/fs*1e3;
                patch([blkS(ch,t) blkE(ch,t) blkE(ch,t) blkS(ch,t)]/fs*1e3, [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
            end
        end
        for ii = 1:numel(idx)
            d = chk(idx(ii));
            if ~zoom
                plot(d.t_ms, d.raw(q,:), 'Color', [0.8 0.8 0.8]);
                plot(d.t_ms, d.blanked(q,:), 'Color', [0.3 0.5 1]);
            end
            plot(d.t_ms, d.filtered(q,:), 'Color', [0.85 0.1 0.1]);
        end
        ev = Texp.events{t1};
        for e = ev(:)'/1e3, xline(e, 'm-'); end
        if zoom, yline(thr, 'k--'); yline(-thr, 'k:'); end
        xlim([-3 25]);  ylim(yl);  grid on;
        if isempty(ev)
            ttl2 = 'no stimulation: no blank';
        else
            ttl2 = sprintf('blank %.1f-%.1f ms (%d trials)', min(bl), max(bl), numel(idx));
        end
        title({Texp.label{t1}, ttl2}, 'FontSize', 7, 'Interpreter', 'none');
        if k > (nRow-1)*nCol, xlabel('ms from trigger'); end
        if mod(k-1, nCol) == 0, ylabel('\muV'); end
    end
    if zoom
        sgtitle(sprintf('ch %d: band-passed (red), threshold dashed, magenta lines = pulse times, grey band = blank', ch));
    else
        sgtitle(sprintf('ch %d: raw (grey), blanked (blue), band-passed (red); magenta lines = pulse times; grey band = blank', ch));
    end
end
end

function pk = detectChannel_(x, thr, L, refr)
cr = find(x(1:end-1) > thr & x(2:end) <= thr) + 1;
if isempty(cr), pk = zeros(0,1); return, end
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

function Z = robustZ_(M, grp)
% distance from the median in robust SDs, across trials (columns), separately for each group of trials
Z = nan(size(M));
for g = unique(grp(:))'
    ix = find(grp == g);
    A = M(:, ix);
    md = median(A, 2, 'omitnan');
    s = 1.4826*median(abs(bsxfun(@minus, A, md)), 2, 'omitnan');
    s(s == 0 | isnan(s)) = NaN;
    Z(:, ix) = bsxfun(@rdivide, bsxfun(@minus, A, md), s);
end
end

function A = readNPY_(fpath, maxCount)
if nargin < 2, maxCount = Inf; end
assert(isfile(fpath), 'Missing: %s', fpath);
fid = fopen(fpath,'r');  cl = onCleanup(@() fclose(fid)); %#ok<NASGU>
magic = fread(fid,6,'*uint8')';
assert(isequal(magic, uint8([147 78 85 77 80 89])), 'Not a .npy file: %s', fpath);
major = fread(fid,1,'uint8');  fread(fid,1,'uint8');
if major == 1, hlen = fread(fid,1,'uint16','l'); else, hlen = fread(fid,1,'uint32','l'); end
header = fread(fid,hlen,'*char')';
descr = regexp(header,'''descr''\s*:\s*''([^'']+)''','tokens','once');  descr = descr{1};
shpT  = regexp(header,'''shape''\s*:\s*\(([^)]*)\)','tokens','once');
shp   = str2num(['[' strrep(shpT{1},',',' ') ']']); %#ok<ST2NM>
if isempty(shp), shp = 1; end
endian = 'l';  if startsWith(descr,'>'), endian = 'b'; end
switch descr(2:end)
    case 'i8', prec = 'int64=>double';  case 'i4', prec = 'int32=>double';
    case 'i2', prec = 'int16=>double';  case 'i1', prec = 'int8=>double';
    case 'u8', prec = 'uint64=>double'; case 'u4', prec = 'uint32=>double';
    case 'u2', prec = 'uint16=>double'; case 'u1', prec = 'uint8=>double';
    case 'f8', prec = 'double=>double'; case 'f4', prec = 'single=>double';
    otherwise, error('Unsupported dtype "%s"', descr);
end
A = fread(fid, min(prod(shp), maxCount), prec, 0, endian);
end

function cm = redblue_()
n = 128;  a = linspace(0,1,n)';
cm = [[a, a, ones(n,1)]; [ones(n,1), flipud(a), flipud(a)]];
cm(n,:) = [1 1 1];
end