%% ====================================================================
%  s1_noisegrid_events.m  —  CJ280 noisegrid, STEP 1 of 4
%
%  Photodiode (AI0, NI-DAQ) -> frame times -> noise-update onsets, carried
%  onto the probe clock with the 1 Hz sync pulses. Writes everything the
%  later steps need into the recording folder.
%
%  Established by the v4 diagnostic on this session:
%    - AI0 strobes once per video frame, 2.4 ms high, 28.8 % duty
%    - frame period 8.3333 ms  ->  120 Hz display
%    - ZERO dropped frames across 90 trials
%    - updates every 4th frame, starting at frame 1 of each trial
%      -> floor((2401-1)/4)+1 = 601 per trial, 90*601 = 54090 updates,
%         which is exactly nUpdates in the parameter file
%
%  Pipeline:
%    s1_noisegrid_events   <- you are here
%    s2_noisegrid_stimulus    which randels were white/black per update
%    s3_noisegrid_mua         threshold crossings per channel (slow, cached)
%    s4_noisegrid_rf          STA maps
% ====================================================================

clear; clc; close all;

%% ==================== CONFIG =======================================
REC_ROOT = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/CJ283/CJ283_noisegrid_110713_2026-09-16_11-07-18';
OUT_NAME = 'noisegrid_events.mat';   % written into the recording folder

AP_DIR   = '';          % empty = auto-detect from structure.oebin
DAQ_DIR  = '';
DIODE_CH = 1;           % AI0

% detection
THR_HI_FRAC  = 0.70;    % Schmitt trigger, fractions of the p1..p99 range
THR_LO_FRAC  = 0.30;
THR_MANUAL   = [];      % [hi lo] raw counts, overrides the fractions
MIN_STATE_US = 100;

% frame -> update mapping.
FRAMES_PER_UPDATE  = [];   % expected 4. Empty = derive from the frame period.
                           % Used here only for reporting; s2 confirms it
                           % against the logged frameInterval.
NOMINAL_UPDATE_MS  = 33.33;

% sync
SYNC_STATE  = 1;
MATCH_TOL_S = 0.25;
OUTLIER_MAD = 5;        % pulses this many robust SDs off the fit are dropped

% tolerances
MAX_SYNC_RESID_MS = 0.5;
MAX_DROPPED_FRAC  = 0.001;
MIN_TRIAL_FRAMES  = 0.5;   % Reject "trials" with fewer frames than this.
                           % Below 1 it is a fraction of the median frames
                           % per trial; 1 or more is an absolute count.
                           % A stray flash separated from the train by a
                           % long gap otherwise becomes a one-frame trial:
                           % CJ283 gave 52 trials against 50 in the
                           % datafile, and s2 then refused to run.

SCAN_DECIM = 25;
DO_PLOTS   = true;

%% ==================== LOCATE THE RECORDING =========================
ob = dir(fullfile(REC_ROOT,'**','structure.oebin'));
assert(~isempty(ob),'No structure.oebin under %s',REC_ROOT);
rec_dir = ob(1).folder;
info    = jsondecode(fileread(fullfile(rec_dir,'structure.oebin')));
cont    = info.continuous;
if isstruct(cont), cont = num2cell(cont(:)'); end
names   = cellfun(@(s) string(s.folder_name), cont);

if isempty(AP_DIR)
    k = find(contains(names,'ProbeA-AP','IgnoreCase',true),1);
    assert(~isempty(k),'No ProbeA-AP stream'); AP_DIR = char(names(k));
end
if isempty(DAQ_DIR)
    k = find(contains(names,'NI-DAQ','IgnoreCase',true),1);
    assert(~isempty(k),'No NI-DAQ stream'); DAQ_DIR = char(names(k));
end
iAP  = find(names==string(AP_DIR),1);
iDAQ = find(names==string(DAQ_DIR),1);
FS_AP  = double(cont{iAP}.sample_rate);
N_AP   = double(cont{iAP}.num_channels);
BV_AP  = cont{iAP}.channels(1).bit_volts;
FS_DAQ = double(cont{iDAQ}.sample_rate);
N_DAQ  = double(cont{iDAQ}.num_channels);
BV     = cont{iDAQ}.channels(1).bit_volts;
chName = arrayfun(@(c) string(c.channel_name), cont{iDAQ}.channels);

fprintf('recording : %s\n', rec_dir);
fprintf('AP        : %s  %g Hz  %d ch  bit_volts %g\n', AP_DIR, FS_AP, N_AP, BV_AP);
fprintf('DAQ       : %s  %g Hz  %d ch\n', DAQ_DIR, FS_DAQ, N_DAQ);
fprintf('diode     : ch %d = %s\n\n', DIODE_CH, chName(DIODE_CH));

daqDat = fullfile(rec_dir,'continuous',stripSep(DAQ_DIR),'continuous.dat');
f  = dir(daqDat);  nD = f.bytes/(N_DAQ*2);
mm = memmapfile(daqDat,'Format',{'int16',[N_DAQ nD],'x'});
d  = double(mm.Data.x(DIODE_CH,:))';
t  = (0:nD-1)'/FS_DAQ;
fprintf('%s: %d samples, %.1f s, %d .. %d counts (%.3f .. %.3f V)\n', ...
        chName(DIODE_CH), nD, t(end), min(d), max(d), min(d)*BV, max(d)*BV);

if DO_PLOTS
    nDec = floor(nD/SCAN_DECIM);
    fid  = fopen(daqDat,'r');
    X    = fread(fid,[N_DAQ nDec],sprintf('%d*int16=>double',N_DAQ),(SCAN_DECIM-1)*N_DAQ*2);
    fclose(fid);
    tDec = (0:size(X,2)-1)'*SCAN_DECIM/FS_DAQ;
    xd   = X(DIODE_CH,:);
end

%% ==================== FRAME DETECTION ==============================
lo = prctile(d,1);  hi = prctile(d,99);
if isempty(THR_MANUAL)
    thrHi = lo + THR_HI_FRAC*(hi-lo);
    thrLo = lo + THR_LO_FRAC*(hi-lo);
else
    thrHi = max(THR_MANUAL);  thrLo = min(THR_MANUAL);
end
fprintf('p1 %.0f | p99 %.0f | Schmitt hi %.0f / lo %.0f\n', lo, hi, thrHi, thrLo);

upC = find(d(1:end-1) <= thrHi & d(2:end) >  thrHi) + 1;
dnC = find(d(1:end-1) >= thrLo & d(2:end) <  thrLo) + 1;
ev  = sortrows([upC(:) ones(numel(upC),1); dnC(:) -ones(numel(dnC),1)], 1);

up = zeros(numel(upC),1);  dn = zeros(numel(upC),1);
k = 0;  armed = true;
for e = 1:size(ev,1)
    if armed && ev(e,2) > 0
        k = k+1;  up(k) = ev(e,1);  armed = false;
    elseif ~armed && ev(e,2) < 0
        dn(k) = ev(e,1);  armed = true;
    end
end
if k > 0 && dn(k) == 0, k = k-1; end
up = up(1:k);  dn = dn(1:k);

minSamp = max(1, round(MIN_STATE_US/1e6*FS_DAQ));
keep = (dn-up) >= minSamp;
up = up(keep);  dn = dn(keep);
fprintf('frames detected: %d\n', numel(up));

tFrame = t(up);
durMs  = (t(dn)-t(up))*1000;
perMs  = diff(tFrame)*1000;

%% ==================== TRIALS =======================================
% Segment, reject trials too short to be a real presentation, then segment
% again on the surviving frames. A stray flash separated from the train by
% a long gap becomes a one-frame "trial", and the spurious trials then fail
% against callbackCounter in s2. Re-segmenting is necessary because
% removing the stray frame also removes the two artificial gaps it created.
dnT = dn;                       % falling edges, kept alongside tFrame
for pass = 1:2
    perMs  = diff(tFrame)*1000;
    gapThr = max(6*median(perMs), 100);
    brk    = find(perMs > gapThr);
    trS    = [1; brk+1];   trE = [brk; numel(tFrame)];
    nFr    = trE - trS + 1;

    if MIN_TRIAL_FRAMES < 1
        minFr = MIN_TRIAL_FRAMES*median(nFr);
    else
        minFr = MIN_TRIAL_FRAMES;
    end
    bad = find(nFr < minFr);
    if isempty(bad) || pass == 2, break; end

    fprintf('\nrejecting %d spurious trial(s) with fewer than %.0f frames:\n', ...
            numel(bad), minFr);
    for i = bad(:)'
        fprintf('   %d frame(s) at %.2f s\n', nFr(i), tFrame(trS(i)));
    end
    kf = true(numel(tFrame),1);
    for i = bad(:)', kf(trS(i):trE(i)) = false; end
    tFrame = tFrame(kf);  durMs = durMs(kf);  dnT = dnT(kf);
    fprintf('%d frames remain\n', numel(tFrame));
end

trT0   = tFrame(trS);  trT1 = t(dnT(trE));
nTrial = numel(trS);

inTrial = true(numel(perMs),1);  inTrial(brk) = false;
pIn     = perMs(inTrial);
frameMs = median(pIn);
nSkip   = round(pIn/frameMs);
nDrop   = sum(nSkip(nSkip>=2) - 1);

fprintf('\n--- frames -------------------------------------------------\n');
fprintf('high duration  : median %.3f ms (duty %.1f %%)\n', ...
        median(durMs), 100*median(durMs)/frameMs);
fprintf('frame period   : median %.4f ms -> %.2f Hz\n', frameMs, 1000/frameMs);
fprintf('dropped frames : %d (%.4f %%)\n', nDrop, 100*nDrop/numel(pIn));
fprintf('trials         : %d, %.4f s median, %d..%d frames each\n', ...
        nTrial, median(trT1-trT0), min(nFr), max(nFr));
assert(nDrop/numel(pIn) <= MAX_DROPPED_FRAC, ...
    ['%d dropped frames (%.3f %%) exceeds the tolerance. Counting updates ' ...
     'in fours from each trial start is no longer safe — the frame-drop log ' ...
     'from the parameter file is needed to reindex.'], nDrop, 100*nDrop/numel(pIn));

%% ==================== FRAME -> UPDATE ==============================
if isempty(FRAMES_PER_UPDATE)
    FRAMES_PER_UPDATE = round(NOMINAL_UPDATE_MS/frameMs);
end

% NOTE: s1 deliberately does NOT decide how many updates a trial contained.
% A trial's frame count does not determine it. After the last update the
% display keeps showing that same pattern until the trial ends, so 2401 and
% 2406 frames are both 601 updates, and floor((nFr-1)/F)+1 invents an update
% at frame 2405 that never occurred. The authority is callbackCounter in the
% datafile; s2 reads it and builds update_t from the frame times saved here.
fprintf('\n--- frames per trial ---------------------------------------\n');
fprintf('frames per trial   : min %d, median %.0f, max %d\n', ...
        min(nFr), median(nFr), max(nFr));
fprintf('frames per update  : %d nominal (%.4f ms)\n', ...
        FRAMES_PER_UPDATE, FRAMES_PER_UPDATE*frameMs);
fprintf('implied updates    : ~%.0f per trial, ~%.0f total (INDICATIVE ONLY)\n', ...
        median(nFr)/FRAMES_PER_UPDATE, sum(nFr)/FRAMES_PER_UPDATE);
fprintf('s2 sets the exact count from callbackCounter and builds update_t.\n');


%% ==================== DAQ -> PROBE CLOCK ===========================
d0 = readNPYlocal(fullfile(rec_dir,'continuous',stripSep(DAQ_DIR),'sample_numbers.npy'),1);
p0 = readNPYlocal(fullfile(rec_dir,'continuous',stripSep(AP_DIR), 'sample_numbers.npy'),1);
sdir = fullfile(rec_dir,'events',stripSep(DAQ_DIR),'TTL');
pdir = fullfile(rec_dir,'events',stripSep(AP_DIR), 'TTL');
dsn = readNPYlocal(fullfile(sdir,'sample_numbers.npy'));
dst = readNPYlocal(fullfile(sdir,'states.npy'));
psn = readNPYlocal(fullfile(pdir,'sample_numbers.npy'));
pst = readNPYlocal(fullfile(pdir,'states.npy'));

sync_daq = (dsn(dst==SYNC_STATE) - double(d0))/FS_DAQ;
sync_npx = (psn(pst==SYNC_STATE) - double(p0))/FS_AP;

fprintf('\n--- sync ---------------------------------------------------\n');
% Which TTL lines are present. If SYNC_STATE picks up more than the 1 Hz
% line, the extra events are matched as if they were sync pulses and the
% residual inflates for a reason that has nothing to do with clock drift.
uD = unique(dst);  uP = unique(pst);
fprintf('DAQ TTL states  :');
for s = uD(:)', fprintf(' %d(n=%d)', s, sum(dst==s)); end
fprintf('\nprobe TTL states:');
for s = uP(:)', fprintf(' %d(n=%d)', s, sum(pst==s)); end
fprintf('\nusing state %d: DAQ %d pulses, probe %d pulses\n', ...
        SYNC_STATE, numel(sync_daq), numel(sync_npx));
fprintf('ITI: DAQ %.6f s (SD %.6f) | probe %.6f s (SD %.6f)\n', ...
        mean(diff(sync_daq)), std(diff(sync_daq)), ...
        mean(diff(sync_npx)), std(diff(sync_npx)));

j  = interp1(sync_npx,(1:numel(sync_npx))',sync_daq,'nearest','extrap');
dt = sync_npx(j) - sync_daq;
ok = abs(dt-median(dt)) < MATCH_TOL_S;
[~,uj] = unique(j(ok),'stable');  ia = find(ok);  ia = ia(uj);
assert(numel(ia) >= 20,'Only %d sync pulses matched.',numel(ia));

xS = sync_daq(ia);  yS = sync_npx(j(ia));
P0 = polyfit(xS, yS, 1);
r0 = (yS - polyval(P0,xS))*1000;

% Robust refit: drop pulses more than OUTLIER_MAD robust SDs from the fit
% and refit on the rest. A few mis-timed pulses inflate an ordinary
% least-squares residual badly while leaving the transform itself fine.
keepR = true(size(xS));
for it = 1:5
    Pk = polyfit(xS(keepR), yS(keepR), 1);
    rk = (yS - polyval(Pk,xS))*1000;
    s  = 1.4826*median(abs(rk - median(rk)));
    nk = abs(rk - median(rk)) <= OUTLIER_MAD*max(s,eps);
    if isequal(nk,keepR), break; end
    keepR = nk;
end
P     = polyfit(xS(keepR), yS(keepR), 1);
resid = (yS - polyval(P,xS))*1000;
residK= resid(keepR);

fprintf('matched %d/%d | slope %.9f (%+.2f ppm) | offset %+.6f s\n', ...
        numel(ia), numel(sync_daq), P(1), (P(1)-1)*1e6, P(2));
fprintf('residual: all pulses %.4f ms SD | robust (%d kept, %d dropped) %.4f ms SD\n', ...
        std(r0), sum(keepR), sum(~keepR), std(residK));
fprintf('residual median |dev| %.4f ms, max %.4f ms\n', ...
        median(abs(resid-median(resid))), max(abs(resid)));

% Shape of the residual decides whether this is a few bad pulses or a
% genuinely wrong pairing.
pr = polyfit(xS, resid, 1);
fprintf('residual trend %+.6f ms/s over the session\n', pr(1));
nBig = sum(abs(resid) > 3*std(residK));
fprintf('%d pulse(s) beyond 3x the robust SD (%.1f%%)\n', nBig, 100*nBig/numel(resid));

figure('Color','w','Position',[80 80 1150 420],'Name','sync diagnostic','NumberTitle','off');
subplot(1,3,1);
plot(xS, resid, '.', 'MarkerSize', 8); hold on;
plot(xS(~keepR), resid(~keepR), 'ro','MarkerSize',6);
yline(0,'k:'); grid on; box on;
xlabel('DAQ time (s)'); ylabel('residual (ms)');
title(sprintf('all %.3f ms SD, robust %.3f ms', std(r0), std(residK)));
subplot(1,3,2);
histogram(resid, 100); grid on; box on;
xlabel('residual (ms)'); ylabel('count'); title('residual distribution');
subplot(1,3,3);
plot(diff(sync_daq),'.'); hold on; plot(diff(sync_npx),'.');
grid on; box on; xlabel('pulse'); ylabel('ITI (s)');
legend('DAQ','probe','Location','best'); title('pulse intervals');

if std(residK) > MAX_SYNC_RESID_MS
    error(['Sync residual %.4f ms (robust) exceeds %.2f ms, and %d of %d ' ...
           'pulses were dropped as outliers. The residual trend is %+.4f ms/s. ' ...
           'A ramp or step means the pulse correspondence is wrong over part ' ...
           'of the recording; scattered outliers with a good robust residual ' ...
           'would be harmless. Check the sync diagnostic figure.'], ...
           std(residK), MAX_SYNC_RESID_MS, sum(~keepR), numel(keepR), pr(1));
elseif std(r0) > MAX_SYNC_RESID_MS
    warning(['Raw residual %.4f ms exceeds the tolerance but the robust fit ' ...
             'is %.4f ms after dropping %d outlying pulse(s). Proceeding on ' ...
             'the robust fit -- inspect the sync figure to confirm the ' ...
             'outliers are isolated rather than a trend.'], ...
             std(r0), std(residK), sum(~keepR));
end

frame_t  = polyval(P, tFrame);          % every frame, probe clock
trial_t0 = polyval(P, trT0);
trial_t1 = polyval(P, trT1);

fprintf('frames on probe clock: %.4f .. %.4f s\n', frame_t(1), frame_t(end));

%% ==================== SAVE =========================================
ev = struct();
ev.rec_dir      = rec_dir;
ev.ap_dir       = AP_DIR;
ev.daq_dir      = DAQ_DIR;
ev.fs_ap        = FS_AP;
ev.n_ap         = N_AP;
ev.bitvolts_ap  = BV_AP;
ev.fs_daq       = FS_DAQ;
ev.diode_ch     = DIODE_CH;
ev.diode_name   = char(chName(DIODE_CH));

ev.frame_t      = frame_t;              % [nFrames x 1] probe clock, seconds
ev.frame_dur_ms = durMs;
ev.frame_ms     = frameMs;
ev.frame_hz     = 1000/frameMs;
ev.frames_per_update = FRAMES_PER_UPDATE;   % nominal; s2 confirms it
ev.n_dropped_frames  = nDrop;

ev.trial_frame0 = trS(:);               % index into frame_t of each trial's
ev.trial_frame1 = trE(:);               % first and last detected frame
ev.trial_t0     = trial_t0;
ev.trial_t1     = trial_t1;
ev.trial_nframes= nFr(:);
ev.n_trials     = nTrial;
ev.min_trial_frames = MIN_TRIAL_FRAMES;

ev.sync_P       = P;
ev.sync_resid_ms= std(residK);
ev.sync_resid_all_ms = std(r0);
ev.sync_ppm     = (P(1)-1)*1e6;
ev.sync_n       = sum(keepR);
ev.sync_n_dropped = sum(~keepR);

ev.thr_hi       = thrHi;
ev.thr_lo       = thrLo;
ev.min_state_us = MIN_STATE_US;
ev.created      = datetime('now');
ev.script       = mfilename('fullpath');

outFile = fullfile(rec_dir, OUT_NAME);
save(outFile, '-struct', 'ev', '-v7.3');
fprintf('\nsaved -> %s\n', outFile);

%% ==================== PLOTS ========================================
if ~DO_PLOTS, return; end
figure('Color','w','Position',[40 40 1260 900],'Name','s1 events');

subplot(4,2,1:2);
plot(tDec,xd,'k-'); hold on; yline(thrHi,'r--'); yline(thrLo,'b--');
axis tight; box on; xlabel('DAQ time (s)'); ylabel('AI0 (counts)');
title(sprintf('%s, whole session (decimated 1:%d), %d frames, %d trials', ...
      chName(DIODE_CH), SCAN_DECIM, numel(tFrame), nTrial));

subplot(4,2,3);
histogram(xd,200); set(gca,'YScale','log'); hold on;
xline(thrHi,'r--'); xline(thrLo,'b--');
xlabel('counts'); ylabel('count (log)'); title('amplitude histogram');

subplot(4,2,4);
w = tFrame(1) + [0.5 0.6];  sel = t>=w(1) & t<=w(2);
plot(t(sel),d(sel),'k.-','MarkerSize',3); hold on;
yline(thrHi,'r--'); yline(thrLo,'b--');
fw = tFrame(tFrame>=w(1) & tFrame<=w(2));
uw = tFrame(trS(1) + (0:FRAMES_PER_UPDATE:(nFr(1)-1))');
uw = uw(uw>=w(1) & uw<=w(2));
plot(fw, thrHi*ones(size(fw)),'g^','MarkerFaceColor','g','MarkerSize',4);
plot(uw, thrHi*ones(size(uw)),'mv','MarkerFaceColor','m','MarkerSize',7);
axis tight; box on; xlabel('DAQ time (s)');
title('100 ms: frames (green) and update onsets (magenta)');

subplot(4,2,5);
w = tFrame(1) + [0.5 0.53];  sel = t>=w(1) & t<=w(2);
plot(t(sel),d(sel),'k.-','MarkerSize',6); hold on;
yline(thrHi,'r--'); yline(thrLo,'b--');
axis tight; box on; xlabel('DAQ time (s)'); title('30 ms: single flash');

subplot(4,2,6);
histogram(pIn(pIn < 4*frameMs), 300);
xlabel('frame period (ms)'); ylabel('count');
title(sprintf('%.4f ms = %.2f Hz, %d dropped', frameMs, 1000/frameMs, nDrop));

subplot(4,2,7);
plot(trial_t0, nFr, 'o-','MarkerSize',4); grid on;
xlabel('trial onset, probe clock (s)'); ylabel('frames in trial');
title(sprintf('%d trials, %d..%d frames each', nTrial, min(nFr), max(nFr)));

subplot(4,2,8);
plot(xS, resid,'.','MarkerSize',8); hold on;
plot(xS(~keepR), resid(~keepR),'ro','MarkerSize',5);
yline(0,'k:'); grid on;
xlabel('DAQ time (s)'); ylabel('residual (ms)');
title(sprintf('sync %.4f ms SD robust, %+.2f ppm', std(residK), (P(1)-1)*1e6));

%% ==================== LOCAL FUNCTIONS ==============================
function s = stripSep(s)
    s = regexprep(char(s),'[\\/]+$','');
end

function data = readNPYlocal(fpath, maxCount)
    if nargin < 2, maxCount = Inf; end
    assert(isfile(fpath),'Missing: %s',fpath);
    fid = fopen(fpath,'r');
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    magic = fread(fid,6,'*uint8')';
    assert(isequal(magic,uint8([147 78 85 77 80 89])),'Not a .npy file: %s',fpath);
    major = fread(fid,1,'uint8'); fread(fid,1,'uint8');
    if major==1, hlen = fread(fid,1,'uint16','l');
    else,        hlen = fread(fid,1,'uint32','l'); end
    header = fread(fid,hlen,'*char')';
    descr = regexp(header,'''descr''\s*:\s*''([^'']+)''','tokens','once'); descr = descr{1};
    fo    = regexp(header,'''fortran_order''\s*:\s*(\w+)','tokens','once');
    shpT  = regexp(header,'''shape''\s*:\s*\(([^)]*)\)','tokens','once');
    shp   = str2num(['[' strrep(shpT{1},',',' ') ']']); %#ok<ST2NM>
    if isempty(shp), shp = 1; end
    endian = 'l'; if startsWith(descr,'>'), endian='b'; end
    switch descr(2:end)
        case 'i8', prec='int64=>double';  case 'i4', prec='int32=>double';
        case 'i2', prec='int16=>double';  case 'i1', prec='int8=>double';
        case 'u8', prec='uint64=>double'; case 'u4', prec='uint32=>double';
        case 'u2', prec='uint16=>double'; case 'u1', prec='uint8=>double';
        case 'f8', prec='double=>double'; case 'f4', prec='single=>double';
        case 'b1', prec='uint8=>double';
        otherwise, error('Unsupported dtype "%s"',descr);
    end
    nRead = min(prod(shp),maxCount);
    data  = fread(fid,nRead,prec,0,endian);
    if numel(shp)>1 && isinf(maxCount)
        if strcmpi(fo{1},'True'), data = reshape(data,shp);
        else,                     data = reshape(data,fliplr(shp))'; end
    end
end