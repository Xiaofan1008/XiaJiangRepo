%% ====================================================================
%  s2_noisegrid_stimulus.m  —  CJ280 noisegrid, STEP 2 of 4
%
%  Reconstructs which randels were white and which were black on each of
%  the 54090 noise updates, by replaying the stimulus RNG stream.
%
%  Neurostim must be on the path so that load() can instantiate the
%  neurostim.cic object. Nothing from Neurostim is CALLED except
%  parameter.get() to read the logged values -- the replay below is our
%  own, and reproduces mapping.noisegridPassive/ternary exactly:
%
%      vals   = bgColor*ones(1,nRands);
%      nOn    = round(sparsity*nRands);
%      inds   = randperm(nRands,nOn);
%      factor = 1 + contrast*sign(rand(1,nOn)-0.5);   % contrast forced to 1
%      vals(inds) = bgColor*factor;                   % -> 0 or 2*bgColor
%
%  so factor==2 is white (ON) and factor==0 is black (OFF).
%
%  Three things from noiseclut.m that the replay depends on:
%
%  1. THE DRY RUN. setupCallbacks() logs rngState (line 414), sets
%     callbackCtr=-1 (line 422), then calls runCallbacks() once (line 429)
%     to size the CLUT. The logged state is therefore PRE-dry-run, and the
%     replay is: restore state -> one discarded call -> callbackCounter
%     real calls. Omit it and every update is off by one draw.
%
%  2. RANDEL INDEXING. noisegrid.beforeTrial builds
%     im = reshape(1:prod(size), size_h, size_w), and MATLAB is
%     column-major, so randel r sits at row mod(r-1,size_h)+1 and column
%     floor((r-1)/size_h)+1. Index runs DOWN columns.
%
%  3. VALIDATION. logInfo() stores clut(:,1:nRandelsToLog) at the end of
%     each trial, nRandelsToLog = round(0.1*6400) = 640. The last replayed
%     frame of each trial must match those 640 values exactly. Every trial
%     is checked below; one bad draw fails all of them.
%
%  Writes noisegrid_stimulus.mat into the recording folder.
% ====================================================================

clear; clc; close all;

%% ==================== CONFIG =======================================
NEUROSTIM_DIR = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Marmo_Visual_Code/neurostim';
REC_ROOT      = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/CJ283/CJ283_noisegrid_110713_2026-09-16_11-07-18';

% The Neurostim datafile sits in the recording folder alongside the raw
% data, e.g. ...\recording1\CJ280.noisegrid.162113.mat. Leave empty to
% find it there; set a full path to override.
PARAM_FILE    = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/CJ283/CJ283_noisegrid_110713_2026-09-16_11-07-18/Record Node 101/experiment1/recording1/CJ283.noisegrid.110713.mat';

STIM_NAME = 'noise';                    % c.noise
OUT_NAME  = 'noisegrid_stimulus.mat';
EVENTS_IN = 'noisegrid_events.mat';     % from s1, for the cross-check

STRICT_VALIDATE = true;   % abort if any trial fails the logged-CLUT check
DO_PLOTS        = true;
ALLOW_SHORT_DATAFILE = true;  % If Neurostim was interrupted (c.stage == 1)
                              % the datafile holds fewer trials than the
                              % diode saw. Use the first N diode trials, on
                              % the assumption that trials are logged in
                              % order so the saved ones are the first
                              % presented. Set false to refuse instead.
DROP_BAD_TRIALS = true;       % Drop trials where the diode saw fewer frames
                              % than callbackCounter requires. Those trials
                              % lost frames to a photodiode dropout, and the
                              % missing update onsets cannot be recovered --
                              % the alternative is to discard the session.

%% ==================== LOAD THE PARAMETER FILE ======================
assert(isfolder(NEUROSTIM_DIR), 'Neurostim folder not found: %s', NEUROSTIM_DIR);
addpath(genpath(NEUROSTIM_DIR));
assert(exist('neurostim.cic','class') == 8, ...
    ['neurostim.cic is not visible after addpath. The +neurostim package ' ...
     'folder must sit directly inside %s.'], NEUROSTIM_DIR);

rec_dir = findRecDir(REC_ROOT);
fprintf('recording: %s\n', rec_dir);

if isempty(PARAM_FILE)
    cand = dir(fullfile(rec_dir,'*.mat'));
    mine = [{OUT_NAME, EVENTS_IN}, {'liveacute-noisegrid-rf-cache.mat'}];
    cand = cand(~ismember({cand.name}, mine) & ~startsWith({cand.name},'noisegrid_'));
    assert(~isempty(cand), ...
        ['No Neurostim datafile found in %s. Set PARAM_FILE explicitly.'], rec_dir);
    if numel(cand) > 1
        fprintf('%d candidate datafiles, using the first:\n', numel(cand));
        for i = 1:numel(cand), fprintf('   %s\n', cand(i).name); end
    end
    PARAM_FILE = fullfile(cand(1).folder, cand(1).name);
end
assert(isfile(PARAM_FILE), 'Parameter file not found: %s', PARAM_FILE);
fprintf('loading %s ...\n', PARAM_FILE);
S  = load(PARAM_FILE);
fn = fieldnames(S);
ix = find(cellfun(@(f) isa(S.(f),'neurostim.cic'), fn), 1);
assert(~isempty(ix), ['No neurostim.cic in the file (variables: %s). If ' ...
    'MATLAB warned about an undefined class, the Neurostim version on the ' ...
    'path does not match the one that wrote this file.'], strjoin(fn,', '));
c = S.(fn{ix});
fprintf('loaded cic "%s", paradigm "%s", %d trials\n', c.name, c.paradigm, c.trial);

o = c.(STIM_NAME);
assert(isa(o,'neurostim.stimuli.noisegrid'), ...
       'c.%s is a %s, expected neurostim.stimuli.noisegrid', STIM_NAME, class(o));

trials  = 1:c.trial;
nTrials = numel(trials);

%% ==================== READ THE LOGGED PARAMETERS ===================
gv = @(p) get(o.prms.(p), 'trial', trials, 'atTrialTime', Inf);

rngSt   = gv('rngState');            % [nTrials x nStateElements]
cbCtr   = gv('callbackCounter');     % updates per trial
sparsity= gv('sparsity');
colorL  = gv('color');
size_w  = gv('size_w');
size_h  = gv('size_h');
loggedC = get(o.prms.clutVals,'trial',trials,'atTrialTime',Inf,'matrixIfPossible',false);
ixImage = get(o.prms.ixImage, 'trial',trials,'atTrialTime',Inf,'matrixIfPossible',false);

if ~iscell(loggedC), loggedC = num2cell(loggedC,2); end
if ~iscell(ixImage), ixImage = num2cell(ixImage,2); end

W = unique(size_w);  H = unique(size_h);
assert(isscalar(W) && isscalar(H), 'Grid size changed across trials.');
nRandels = W*H;
fprintf('grid %d x %d (size_h x size_w) = %d randels\n', H, W, nRandels);

sp = unique(sparsity);
assert(isscalar(sp), 'sparsity changed across trials: %s', mat2str(sp'));
nOn = round(sp*nRandels);
fprintf('sparsity %g -> %d active randels per update\n', sp, nOn);

% background luminance: ternary uses o.color(3)
if size(colorL,2) >= 3, bg = colorL(:,3); else, bg = colorL(:,end); end
assert(numel(unique(bg))==1, 'Background colour changed across trials.');
fprintf('background %g -> white %g, black %g\n', bg(1), 2*bg(1), 0);

fprintf('updates per trial: min %d, median %d, max %d (total %d)\n', ...
        min(cbCtr), median(cbCtr), max(cbCtr), sum(cbCtr));

% how many randels were logged for validation
propLog  = get(o.prms.propOfRandelsToLog,'trial',trials,'atTrialTime',Inf);
nLog     = max(min(nRandels,10), round(propLog(1)*nRandels));
fprintf('logged CLUT values per trial: %d randels\n', nLog);

%% ==================== REPLAY =======================================
rs = o.rng;
if iscell(rs), rs = rs{1}; end
assert(isa(rs,'RandStream'), 'o.rng is not a RandStream (%s)', class(rs));
fprintf('RNG stream: %s, seed %d, stream %d of %d\n', ...
        rs.Type, rs.Seed, rs.StreamIndex, rs.NumStreams);

nUpd = sum(cbCtr);
act_idx   = zeros(nOn, nUpd, 'uint16');   % which randels were active
act_pol   = zeros(nOn, nUpd, 'int8');     % +1 white, -1 black
upd_trial = zeros(nUpd, 1, 'uint16');
upd_k     = zeros(nUpd, 1, 'uint16');

prevStream = RandStream.setGlobalStream(rs);
restoreStream = onCleanup(@() RandStream.setGlobalStream(prevStream));

fprintf('\nreplaying %d trials ...\n', nTrials);
tic; p = 0; nBad = 0; badTrials = [];
for i = 1:nTrials
    rs.State = rngSt(i,:);

    % --- the dry run in setupCallbacks (noiseclut.m line 429) ---------
    ternaryDraw(nRandels, nOn);

    lastVals = [];
    for j = 1:cbCtr(i)
        [inds, isWhite] = ternaryDraw(nRandels, nOn);
        p = p + 1;
        act_idx(:,p)   = uint16(inds);
        act_pol(:,p)   = int8(2*isWhite - 1);
        upd_trial(p)   = i;
        upd_k(p)       = j;
        if j == cbCtr(i)
            lastVals = bg(i)*ones(1,nRandels);
            lastVals(inds) = bg(i)*(1 + (2*isWhite-1));
        end
    end

    % --- validate against the logged CLUT ----------------------------
    lg = loggedC{i};
    if ~isempty(lg)
        lg = lg(end,:);                       % luminance row (XYL: row 3)
        lg = lg(1:min(nLog,numel(lg)));
        if ~isequal(lastVals(1:numel(lg)), lg)
            nBad = nBad + 1;  badTrials(end+1) = i; %#ok<SAGROW>
        end
    end

    if mod(i,10)==0 || i==nTrials
        fprintf('  trial %3d/%d  (%d updates, %.1f s)\n', i, nTrials, p, toc);
    end
end
clear restoreStream

fprintf('\nvalidation: %d/%d trials matched the logged CLUT values\n', ...
        nTrials-nBad, nTrials);
if nBad > 0
    msg = sprintf(['%d trial(s) failed: %s. The replay is not reproducing ' ...
        'what was shown -- most likely the dry-run call or the RNG stream ' ...
        'is wrong. Do NOT use these maps.'], nBad, mat2str(badTrials));
    if STRICT_VALIDATE, error(msg); else, warning(msg); end %#ok<SPERR>
end

%% ==================== CROSS-CHECK AGAINST s1 =======================
%% ==================== UPDATE ONSETS ================================
% callbackCounter is the number of times the noise callback actually ran in
% a trial, so it is the update count. The frame times come from s1; update j
% of trial i is that trial's frame 1+(j-1)*F.
evFile = fullfile(rec_dir, EVENTS_IN);
assert(isfile(evFile), 'Run s1 first: %s not found.', evFile);
E = load(evFile,'frame_t','trial_frame0','trial_frame1','trial_nframes', ...
                'n_trials','frame_ms','frames_per_update','sync_resid_ms');
if E.n_trials ~= nTrials
    if E.n_trials > nTrials && ALLOW_SHORT_DATAFILE
        warning(['The diode saw %d trials but the datafile has %d. Using the ' ...
                 'FIRST %d diode trials and discarding the last %d, which ' ...
                 'were presented but never logged. This assumes Neurostim ' ...
                 'logs trials in order.'], ...
                 E.n_trials, nTrials, nTrials, E.n_trials - nTrials);
    else
        error(['Diode found %d trials, the datafile has %d. The datafile ' ...
               'cannot have MORE trials than were presented, so this is not ' ...
               'an interrupted run.'], E.n_trials, nTrials);
    end
end
nUse = nTrials;                      % trials with both diode and stimulus data

% frames per update: prefer the logged frameInterval over s1's nominal guess
F = E.frames_per_update;
try
    fi = gv('frameInterval');
    fi = median(fi(isfinite(fi)));
    Ffi = round(fi/E.frame_ms);
    fprintf('\nframeInterval logged as %.4f ms / frame %.4f ms -> %d frames\n', ...
            fi, E.frame_ms, Ffi);
    if Ffi ~= F
        warning('s1 assumed %d frames per update, frameInterval says %d. Using %d.', ...
                F, Ffi, Ffi);
    end
    F = Ffi;
catch
    fprintf('\nframeInterval not logged; using s1''s %d frames per update\n', F);
end

need  = (cbCtr(:)-1)*F + 1;                 % frames the updates require
nfr   = E.trial_nframes(1:nUse);
slack = nfr(:) - need;          % frames left over after the last
fprintf('trailing frames after the last update: min %d, median %.0f, max %d\n', ...
        min(slack), median(slack), max(slack));
keepTr = true(nUse,1);
if any(slack < 0)
    b = find(slack < 0);
    fprintf('\n%-7s %10s %10s %10s %10s\n', ...
            'trial','frames','needed','callbacks','short by');
    for i = b(1:min(20,numel(b)))'
        fprintf('%-7d %10d %10d %10d %10d\n', i, nfr(i), need(i), cbCtr(i), -slack(i));
    end
    if DROP_BAD_TRIALS
        keepTr(b) = false;
        warning(['%d trial(s) lost frames to a photodiode dropout and are ' ...
                 'being DROPPED. Their update onsets cannot be placed: the ' ...
                 'frame times inside the gap were never recorded. %d of %d ' ...
                 'trials remain.'], numel(b), sum(keepTr), nUse);
    else
        error(['%d trial(s) report more updates than the diode saw frames ' ...
               'for. Either flashes were missed at the end of those trials ' ...
               'or the frames-per-update is wrong. Set DROP_BAD_TRIALS to ' ...
               'discard them instead.'], numel(b));
    end
end
if any(slack(keepTr) >= F)
    n = sum(slack(keepTr) >= F);
    fprintf(['%d trial(s) have %d or more trailing frames. That is expected: ' ...
             'the last pattern stays on screen until the trial ends.\n'], n, F);
end

% Update onsets for the surviving trials only, and the replayed noise
% subset to match. The two must stay column-for-column aligned.
use = find(keepTr(:))';
nUpdEv   = sum(cbCtr(use));
update_t = zeros(nUpdEv,1);
q = 0;
for i = use
    n = cbCtr(i);
    update_t(q+(1:n)) = E.frame_t(E.trial_frame0(i) + (0:n-1)'*F);
    q = q + n;
end
colKeep = ismember(double(upd_trial), use);
assert(q == sum(colKeep), 'Built %d update times for %d replayed updates.', ...
       q, sum(colKeep));
if ~all(colKeep)
    act_idx   = act_idx(:,colKeep);
    act_pol   = act_pol(:,colKeep);
    upd_trial = upd_trial(colKeep);
    upd_k     = upd_k(colKeep);
    nUpd      = numel(upd_trial);
    fprintf('dropped %d update(s) belonging to excluded trials\n', ...
            sum(~colKeep));
end
assert(issorted(update_t), 'Update times are not monotonic.');
fprintf('update onsets: %d over %d trial(s), %.4f .. %.4f s on the probe clock\n', ...
        nUpdEv, numel(use), update_t(1), update_t(end));
fprintf('sync residual from s1: %.4f ms\n', E.sync_resid_ms);

%% ==================== GEOMETRY =====================================
% randel r -> row mod(r-1,H)+1, col floor((r-1)/H)+1  (column-major)
randel_row = uint16(mod((1:nRandels)'-1, H) + 1);
randel_col = uint16(floor(((1:nRandels)'-1)/H) + 1);

% verify against the logged ixImage rather than trusting the derivation
ix1 = ixImage{1};
if ~isempty(ix1)
    ix1 = reshape(ix1, H, W);
    ok  = isequal(double(ix1), reshape(1:nRandels, H, W));
    fprintf('\nixImage matches reshape(1:%d,%d,%d): %d\n', nRandels, H, W, ok);
    if ~ok
        warning(['Logged ixImage is not the plain column-major reshape. ' ...
                 'Using it directly for the randel->pixel map.']);
        [randel_row, randel_col] = deal(zeros(nRandels,1,'uint16'));
        for r = 1:nRandels
            [rr,cc] = find(ix1==r,1);
            randel_row(r) = rr;  randel_col(r) = cc;
        end
    end
end

gp = @(p) get(o.prms.(p),'trial',trials,'atTrialTime',Inf);
stim = struct();
try
    stim.width  = unique(gp('width'));
    stim.height = unique(gp('height'));
    stim.X      = unique(gp('X'));
    stim.Y      = unique(gp('Y'));
catch ME
    warning('Could not read stimulus geometry: %s', ME.message);
end

% Screen parameters. Nothing downstream hardcodes a viewing distance -- the
% axes in s4/s5 are whatever units Neurostim logged. Those are degrees only
% if c.screen.distance and the physical screen size were configured at
% runtime; otherwise they are pixels, and s4 would still label the axis
% "deg". Recording them here makes the units checkable after the fact.
scr = struct();
for fld = {'distance','width','height','xpixels','ypixels','frameRate'}
    try, scr.(fld{1}) = c.screen.(fld{1}); catch, scr.(fld{1}) = []; end
end
if ~isempty(scr.distance) && ~isempty(scr.width)
    fprintf('\nscreen: distance %g, physical %g x %g, %g x %g px, %g Hz\n', ...
            scr.distance, scr.width, scr.height, scr.xpixels, scr.ypixels, scr.frameRate);
    fprintf('full screen subtends %.1f deg horizontally\n', ...
            2*atand(scr.width/2/scr.distance));
    if isfield(stim,'width') && ~isempty(stim.width)
        fprintf('noise grid logged as %g x %g at (%g, %g) -> ', ...
                stim.width, stim.height, stim.X, stim.Y);
        if stim.width < 200
            fprintf('plausibly DEGREES\n');
        else
            fprintf(['looks like PIXELS. s4/s5 label the axes "deg" ' ...
                     'regardless, so RF sizes and any deg/mm gradient would ' ...
                     'be unscaled.\n']);
        end
    end
else
    warning(['c.screen.distance or width not logged, so the stimulus units ' ...
             'cannot be confirmed as degrees.']);
end
stim.screen = scr;

%% ==================== SAVE =========================================
st = struct();
st.act_idx    = act_idx;      % [nOn x nUpd] uint16, randel indices
st.act_pol    = act_pol;      % [nOn x nUpd] int8, +1 white / -1 black
st.update_t   = update_t;     % [nUpd x 1] onset on the probe clock, seconds
st.upd_trial  = upd_trial;
st.upd_k      = upd_k;
st.n_updates  = nUpd;
st.n_trials   = nTrials;
st.updates_per_trial = cbCtr(:);
st.trials_used       = use(:);
st.trials_dropped    = find(~keepTr(:));
st.frames_per_update = F;
st.n_trials_diode    = E.n_trials;
st.n_trials_unused   = E.n_trials - nUse;
st.frame_ms   = E.frame_ms;
st.sync_resid_ms = E.sync_resid_ms;

st.n_randels  = nRandels;
st.size_w     = W;
st.size_h     = H;
st.randel_row = randel_row;
st.randel_col = randel_col;
st.sparsity   = sp;
st.n_active   = nOn;
st.bg         = bg(1);
st.lum_white  = 2*bg(1);
st.lum_black  = 0;
st.geometry   = stim;

st.param_file = PARAM_FILE;
st.rng_type   = rs.Type;
st.validated  = (nBad == 0);
st.created    = datetime('now');
st.script     = mfilename('fullpath');

outFile = fullfile(rec_dir, OUT_NAME);
save(outFile, '-struct', 'st', '-v7.3');
fprintf('\nsaved -> %s\n', outFile);

%% ==================== DIAGNOSTIC PLOTS =============================
if ~DO_PLOTS, return; end

% how often was each randel white / black? should be flat at ~nUpd*sp/2
cntOn  = accumarray(double(act_idx(act_pol>0)), 1, [nRandels 1]);
cntOff = accumarray(double(act_idx(act_pol<0)), 1, [nRandels 1]);
mapOn  = reshape(cntOn,  H, W);
mapOff = reshape(cntOff, H, W);
expct  = nUpd*sp/2;

figure('Color','w','Position',[60 60 1150 760],'Name','s2 stimulus');

subplot(2,3,1);
imagesc(mapOn); axis image; colorbar;
title(sprintf('times white (expect %.0f)', expct)); xlabel('col'); ylabel('row');

subplot(2,3,2);
imagesc(mapOff); axis image; colorbar;
title(sprintf('times black (expect %.0f)', expct)); xlabel('col'); ylabel('row');

subplot(2,3,3);
histogram(cntOn, 40); hold on; histogram(cntOff, 40);
xline(expct,'r--','LineWidth',1.2);
legend('white','black','expected','Location','best');
xlabel('presentations per randel'); ylabel('count');
title(sprintf('uniformity: SD %.1f vs Poisson %.1f', std(cntOn), sqrt(expct)));

% three example updates rendered as images
for k = 1:3
    subplot(2,3,3+k);
    u  = round(k*nUpd/4);
    im = st.bg*ones(nRandels,1);
    im(double(act_idx(:,u))) = st.bg*(1 + double(act_pol(:,u)));
    imagesc(reshape(im,H,W), [0 2*st.bg]); axis image; colormap(gca,gray);
    title(sprintf('update %d (trial %d, k %d)', u, upd_trial(u), upd_k(u)));
end

fprintf('\nrandel occupancy: white %.1f +/- %.1f, black %.1f +/- %.1f (Poisson %.1f)\n', ...
        mean(cntOn), std(cntOn), mean(cntOff), std(cntOff), sqrt(expct));

%% ==================== LOCAL FUNCTIONS ==============================
function [inds, isWhite] = ternaryDraw(nRands, nOn)
% One invocation of mapping.noisegridPassive/ternary. The two RNG calls
% must happen in this order, with these shapes, or the stream desynchronises.
    inds    = randperm(nRands, nOn);
    isWhite = sign(rand(1,nOn) - 0.5) > 0;
end

function rd = findRecDir(root)
    ob = dir(fullfile(root,'**','structure.oebin'));
    assert(~isempty(ob),'No structure.oebin under %s',root);
    rd = ob(1).folder;
end