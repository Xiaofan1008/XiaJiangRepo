%% ====================================================================
%  s4_noisegrid_rf.m  —  CJ280 noisegrid, STEP 4 of 4
%
%  ON and OFF receptive fields per channel, by reverse correlation of MUA
%  threshold crossings against the reconstructed ternary noise.
%
%  Method
%  ------
%  For every noise update the spikes are counted in BIN_MS bins across
%  LAG_MS. For each randel p and each lag b:
%
%      mOn(p,b)  = mean spike count over the ~2704 updates where p was white
%      mAll(b)   = mean spike count over ALL 54090 updates
%      z(p,b)    = (mOn(p,b) - mAll(b)) / (sd(b)/sqrt(nOn(p)))
%
%  Subtracting mAll matters: 640 randels flash on every single update, so
%  every channel has a stimulus-locked transient that has nothing to do
%  with its receptive field. mAll removes exactly that.
%
%  The z normalisation makes channels comparable. A channel with 3116
%  crossings and one with 35021 both end up on a scale where |z| ~ 1 is
%  noise, so the same colour limits are meaningful everywhere and a weak
%  channel is not mistaken for a strong one just because it fired more.
%
%  No lag window is fixed in advance. The full kernel is computed and each
%  channel's own peak lag is found from the data, because latency varies
%  along the probe and a window chosen up front biases that.
%
%  Requires: noisegrid_events.mat (s1), noisegrid_stimulus.mat (s2),
%            noisegrid_mua.mat (s3)
% ====================================================================

clear; clc; close all;

%% ==================== CONFIG =======================================
REC_ROOT = '/Volumes/Shared/Xia/MarmosetData/Saline_Test/CJ283/CJ283_noisegrid_110713_2026-09-16_11-07-18';

CHANNELS = [1:384];            % PROCESSED channels. [] = every channel present in
                          % the MUA file, or e.g. 1:10. These are the ones
                          % computed and written to noisegrid_rf.mat.

% reverse correlation
LAG_MS   = [-50 250];     % relative to update onset
BIN_MS   = 10;
PEAK_WIN = 3;             % bins averaged around the peak for the summary map

% display
PLOT_CHANNELS   = [5:8:384];     % DISPLAYED channels, a subset of CHANNELS.
                          % [] = all of them. Anything left out is still
                          % processed, still saved, and still listed in the
                          % summary table -- it is only not drawn.
MOSAIC_CHANNELS = [1:384];     % same idea for the mosaic. [] = all processed.
MAX_FIGURES     = 40;     % guard: opening one window per channel across 384
                          % channels will exhaust the display
SMOOTH_SD   = 1.0;        % Gaussian SD in randels, DISPLAY ONLY (0 = none)
CLIM_MODE   = 'robust';   % 'robust' = 99.5th pct of |z|, or 'fixed'
CLIM_FIXED  = 6;
FLIP_Y      = true;       % see the note on orientation below
IN_DEGREES  = true;       % axes in visual degrees if geometry is available

% output
PLOT_PER_CHANNEL = true;
PLOT_MOSAIC      = true;
SAVE_MAPS        = true;   % write noisegrid_rf.mat (needed by s5)
SAVE_RATES       = false;  % also store rOn/rOff. Adds ~300 MB per array at
                           % 384 channels and s5 does not use them; the z
                           % maps carry the same structure.
SAVE_FIGS        = false;
FIG_DIR          = '';

CH_BLOCK = 16;            % channels per matrix multiply

%% ==================== LOAD =========================================
ob = dir(fullfile(REC_ROOT,'**','structure.oebin'));
assert(~isempty(ob),'No structure.oebin under %s',REC_ROOT);
rec_dir = ob(1).folder;

E = load(fullfile(rec_dir,'noisegrid_events.mat'));
S = load(fullfile(rec_dir,'noisegrid_stimulus.mat'));
M = load(fullfile(rec_dir,'noisegrid_mua.mat'));

fprintf('events   : %d frames, %d trials, %.4f ms/frame\n', ...
        numel(E.frame_t), E.n_trials, E.frame_ms);
fprintf('stimulus : %d updates, %d randels (%dx%d), %d active/update, validated %d\n', ...
        S.n_updates, S.n_randels, S.size_h, S.size_w, S.n_active, S.validated);
fprintf('mua      : %d channels, %d spikes total\n', ...
        numel(M.channels), sum(M.n_spikes));

assert(isfield(S,'update_t'), ...
    ['This noisegrid_stimulus.mat predates the change that moved update ' ...
     'timing into s2. Re-run s1 then s2.']);
if ~S.validated
    warning('The stimulus reconstruction did not validate. Maps are unreliable.');
end

update_t = S.update_t(:);
nUpd     = numel(update_t);
assert(nUpd == S.n_updates, 'update_t has %d entries, n_updates says %d.', ...
       nUpd, S.n_updates);
assert(size(S.act_idx,2) == nUpd, 'act_idx has %d columns, expected %d.', ...
       size(S.act_idx,2), nUpd);
nPix     = double(S.n_randels);
H        = double(S.size_h);
W        = double(S.size_w);

if isempty(CHANNELS)
    chSel = M.channels(:)';
else
    chSel = CHANNELS(:)';
    missing = setdiff(chSel, M.channels);
    assert(isempty(missing), ...
        'Channel(s) %s are not in the MUA file. Re-run s3 with a wider CHANNELS.', ...
        mat2str(missing));
end
[~, chIx] = ismember(chSel, M.channels);
nCh = numel(chSel);

% which of the processed channels get drawn
plotIx = resolveSubset(PLOT_CHANNELS,   chSel, 'PLOT_CHANNELS');
mosIx  = resolveSubset(MOSAIC_CHANNELS, chSel, 'MOSAIC_CHANNELS');

fprintf('processing %d channel(s); %d drawn individually, %d in the mosaic\n\n', ...
        nCh, numel(plotIx), numel(mosIx));

%% ==================== LAG BINS =====================================
edges  = LAG_MS(1):BIN_MS:LAG_MS(2);
ctrs   = edges(1:end-1) + BIN_MS/2;
nBin   = numel(ctrs);
binSec = BIN_MS/1000;
fprintf('lags: %g to %g ms in %g ms bins (%d bins)\n', LAG_MS(1), LAG_MS(2), BIN_MS, nBin);

%% ==================== STIMULUS INDICATORS ==========================
% Sparse [nPix x nUpd] indicators. Built once and reused for all channels;
% each has ~320*54090 = 17.3M nonzeros.
fprintf('building stimulus indicators ...\n');
tic;
cols   = repmat(uint32(1:nUpd), S.n_active, 1);
isOn   = S.act_pol > 0;
Mon    = sparse(double(S.act_idx(isOn)),  double(cols(isOn)),  1, nPix, nUpd);
Moff   = sparse(double(S.act_idx(~isOn)), double(cols(~isOn)), 1, nPix, nUpd);
clear cols isOn
nOn    = full(sum(Mon, 2));
nOff   = full(sum(Moff,2));
fprintf('  ON %.0f +/- %.1f, OFF %.0f +/- %.1f presentations per randel (%.1f s)\n', ...
        mean(nOn), std(nOn), mean(nOff), std(nOff), toc);
assert(min(nOn) > 50 && min(nOff) > 50, 'Some randels were almost never shown.');

%% ==================== REVERSE CORRELATION ==========================
zOn  = zeros(nPix, nBin, nCh, 'single');
zOff = zeros(nPix, nBin, nCh, 'single');
rOn  = zeros(nPix, nBin, nCh, 'single');
rOff = zeros(nPix, nBin, nCh, 'single');
psth = zeros(nBin, nCh);

fprintf('\nreverse correlation ...\n');
tic;
for b0 = 1:CH_BLOCK:nCh
    b1 = min(b0+CH_BLOCK-1, nCh);
    nb = b1-b0+1;

    % spike counts per update per lag bin, for this block of channels
    C = zeros(nUpd, nBin*nb);
    mAll = zeros(nBin, nb);
    sdAll= zeros(nBin, nb);
    for ci = 1:nb
        st = M.spikes{chIx(b0+ci-1)};
        cc = countPerUpdate(st, update_t, edges);   % [nUpd x nBin]
        C(:, (ci-1)*nBin + (1:nBin)) = cc;
        mAll(:,ci)  = mean(cc,1)';
        sdAll(:,ci) = std(cc,0,1)';
    end

    Son  = reshape(full(Mon  * C), nPix, nBin, nb);
    Soff = reshape(full(Moff * C), nPix, nBin, nb);

    for ci = 1:nb
        k = b0+ci-1;
        mOnP  = Son(:,:,ci)  ./ nOn;      % mean count per update, per randel
        mOffP = Soff(:,:,ci) ./ nOff;

        seOn  = sdAll(:,ci)' ./ sqrt(nOn);
        seOff = sdAll(:,ci)' ./ sqrt(nOff);

        zOn(:,:,k)  = single((mOnP  - mAll(:,ci)') ./ max(seOn, eps));
        zOff(:,:,k) = single((mOffP - mAll(:,ci)') ./ max(seOff,eps));
        rOn(:,:,k)  = single(mOnP  / binSec);
        rOff(:,:,k) = single(mOffP / binSec);
        psth(:,k)   = mAll(:,ci) / binSec;
    end
    fprintf('  channels %d-%d of %d (%.1f s)\n', b0, b1, nCh, toc);
end
clear C Son Soff

%% ==================== PEAK LAG PER CHANNEL =========================
% Under the null, z has unit SD across randels. A real RF shows up as an
% inflated spread at the responsive lag, so the spread is the statistic.
sprOn  = squeeze(std(double(zOn), 0, 1));    % [nBin x nCh]
sprOff = squeeze(std(double(zOff),0, 1));
if nCh == 1, sprOn = sprOn(:); sprOff = sprOff(:); end
spr    = max(sprOn, sprOff);

% only consider positive lags for the peak; negative lags are the control
posBin = find(ctrs > 0);
[pkVal, pkRel] = max(spr(posBin,:), [], 1);
pkBin  = posBin(pkRel);
pkLag  = ctrs(pkBin);

% noise level from the negative lags, where there can be no response
negBin = find(ctrs < -BIN_MS);
if isempty(negBin), nullLvl = ones(1,nCh);
else,               nullLvl = mean(spr(negBin,:), 1); end
snr = pkVal ./ max(nullLvl, eps);

fprintf('\n%-6s %8s %8s %8s %8s %8s\n','ch','spikes','peak ms','spread','null','ratio');
for k = 1:nCh
    fprintf('%-6d %8d %8.0f %8.2f %8.2f %8.2f\n', ...
            chSel(k), M.n_spikes(chIx(k)), pkLag(k), pkVal(k), nullLvl(k), snr(k));
end
fprintf(['\nratio is peak spread over the spread at negative lags. Near 1 ' ...
         'means no receptive field on that channel.\n']);

%% ==================== SPATIAL AXES =================================
row = double(S.randel_row);  col = double(S.randel_col);
haveGeom = isfield(S,'geometry') && isstruct(S.geometry) && ...
           isfield(S.geometry,'width') && ~isempty(S.geometry.width);
if IN_DEGREES && haveGeom
    g  = S.geometry;
    xC = g.X - g.width/2  + ((1:W)-0.5)*g.width/W;
    yC = g.Y + g.height/2 - ((1:H)-0.5)*g.height/H;
    xLab = 'azimuth (deg)'; yLab = 'elevation (deg)';
    fprintf('\ngrid spans %.2f to %.2f deg horizontally, %.2f to %.2f vertically\n', ...
            min(xC), max(xC), min(yC), max(yC));
else
    xC = 1:W;  yC = 1:H;
    xLab = 'randel column'; yLab = 'randel row';
end
if FLIP_Y, yC = fliplr(yC); end

cmap = divergingMap(256);

%% ==================== PER-CHANNEL FIGURES ==========================
if SAVE_FIGS
    if isempty(FIG_DIR), FIG_DIR = fullfile(rec_dir,'rf_figs'); end
    if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
end

if PLOT_PER_CHANNEL && ~isempty(plotIx)
    if numel(plotIx) > MAX_FIGURES && ~SAVE_FIGS
        warning(['%d channels selected for individual figures. Drawing the ' ...
                 'first %d; set PLOT_CHANNELS to choose, or SAVE_FIGS = true ' ...
                 'to write them all to disk instead of opening windows.'], ...
                 numel(plotIx), MAX_FIGURES);
        plotIx = plotIx(1:MAX_FIGURES);
    end
    for k = plotIx(:)'
        b   = pkBin(k);
        wIx = max(1,b-floor(PEAK_WIN/2)) : min(nBin,b+floor(PEAK_WIN/2));

        mapOn  = prep(mean(double(zOn(:,wIx,k)),2),  H, W, SMOOTH_SD);
        mapOff = prep(mean(double(zOff(:,wIx,k)),2), H, W, SMOOTH_SD);
        mapDif = mapOn - mapOff;

        switch lower(CLIM_MODE)
            case 'robust'
                cl = prctile(abs([mapOn(:); mapOff(:)]), 99.5);
                cl = max(cl, 1);
            otherwise
                cl = CLIM_FIXED;
        end

        fig = figure('Color','w','Position',[50 50 1250 720], ...
                     'Name',sprintf('channel %d',chSel(k)),'NumberTitle','off');

        subplot(2,3,1);
        imagesc(xC,yC,mapOn,[-cl cl]); axis image xy; colormap(gca,cmap); colorbar;
        title(sprintf('ON  (white), %g ms', pkLag(k))); xlabel(xLab); ylabel(yLab);

        subplot(2,3,2);
        imagesc(xC,yC,mapOff,[-cl cl]); axis image xy; colormap(gca,cmap); colorbar;
        title(sprintf('OFF (black), %g ms', pkLag(k))); xlabel(xLab);

        subplot(2,3,3);
        imagesc(xC,yC,mapDif,[-cl cl]); axis image xy; colormap(gca,cmap); colorbar;
        title('ON - OFF'); xlabel(xLab);

        subplot(2,3,4);
        plot(ctrs, sprOn(:,k), '-o','MarkerSize',3); hold on;
        plot(ctrs, sprOff(:,k),'-s','MarkerSize',3);
        yline(nullLvl(k),'k:'); xline(pkLag(k),'r--');
        grid on; xlabel('lag (ms)'); ylabel('SD of z across randels');
        legend('ON','OFF','null','Location','best');
        title(sprintf('peak %g ms, ratio %.2f', pkLag(k), snr(k)));

        subplot(2,3,5);
        [~,pOn]  = max(abs(double(zOn(:,b,k))));
        [~,pOff] = max(abs(double(zOff(:,b,k))));
        plot(ctrs, double(zOn(pOn,:,k)),  '-o','MarkerSize',3); hold on;
        plot(ctrs, double(zOff(pOff,:,k)),'-s','MarkerSize',3);
        yline(0,'k:'); xline(0,'r--'); grid on;
        xlabel('lag (ms)'); ylabel('z');
        legend(sprintf('ON randel %d',pOn), sprintf('OFF randel %d',pOff), ...
               'Location','best');
        title('time course at the strongest randel');

        subplot(2,3,6);
        plot(ctrs, psth(:,k), '-k','LineWidth',1.2); grid on;
        xline(0,'r--'); xlabel('lag (ms)'); ylabel('rate (sp/s)');
        title(sprintf('ch %d: %d spikes, %.2f sp/s mean', ...
              chSel(k), M.n_spikes(chIx(k)), mean(psth(:,k))));

        if SAVE_FIGS
            exportgraphics(fig, fullfile(FIG_DIR,sprintf('rf_ch%03d.png',chSel(k))), ...
                           'Resolution',150);
            close(fig);
        end
    end
end

%% ==================== MOSAIC =======================================
if PLOT_MOSAIC && numel(mosIx) > 1
    nM = numel(mosIx);
    nc = ceil(sqrt(nM*1.6));  nr = ceil(nM/nc);
    for pass = 1:2
        ttl = {'ON','OFF'};
        figure('Color','w','Position',[40 40 1350 850], ...
               'Name',sprintf('%s mosaic',ttl{pass}),'NumberTitle','off');
        for kk = 1:nM
            k   = mosIx(kk);
            b   = pkBin(k);
            wIx = max(1,b-floor(PEAK_WIN/2)) : min(nBin,b+floor(PEAK_WIN/2));
            if pass==1, mp = mean(double(zOn(:,wIx,k)),2);
            else,       mp = mean(double(zOff(:,wIx,k)),2); end
            mp = prep(mp, H, W, SMOOTH_SD);
            cl = max(prctile(abs(mp(:)),99.5), 1);
            subplot(nr,nc,kk);
            imagesc(mp,[-cl cl]); axis image off xy; colormap(gca,cmap);
            title(sprintf('%d | %gms | %.1f', chSel(k), pkLag(k), snr(k)), ...
                  'FontSize',8);
        end
        sgtitle(sprintf('%s receptive fields (channel | peak lag | ratio)', ttl{pass}));
    end
end

%% ==================== SAVE =========================================
if SAVE_MAPS
    rf = struct();
    rf.channels   = chSel(:);
    rf.zOn = zOn;  rf.zOff = zOff;
    if SAVE_RATES, rf.rOn = rOn;  rf.rOff = rOff; end
    rf.psth = psth;
    rf.lag_ms = ctrs(:);
    rf.bin_ms = BIN_MS;
    rf.peak_lag_ms = pkLag(:);
    rf.peak_bin = pkBin(:);
    rf.snr = snr(:);
    rf.size_h = H;  rf.size_w = W;
    rf.randel_row = S.randel_row;  rf.randel_col = S.randel_col;
    rf.n_on = nOn;  rf.n_off = nOff;
    rf.x_deg = xC;  rf.y_deg = yC;
    rf.created = datetime('now');
    rf.script  = mfilename('fullpath');
    save(fullfile(rec_dir,'noisegrid_rf.mat'), '-struct','rf','-v7.3');
    fprintf('\nsaved -> %s\n', fullfile(rec_dir,'noisegrid_rf.mat'));
end

%% ==================== LOCAL FUNCTIONS ==============================
function ix = resolveSubset(want, have, name)
% Map a list of channel numbers onto positions in the processed list.
% Requesting a channel that was not processed is a warning, not an error:
% the rest of the selection is still drawn.
    if isempty(want), ix = 1:numel(have); return; end
    [tf, loc] = ismember(want(:)', have);
    if any(~tf)
        warning('%s lists channel(s) %s that are not being processed; skipping them.', ...
                name, mat2str(want(~tf)));
    end
    ix = sort(loc(tf));
    if isempty(ix)
        warning('%s matched no processed channel; nothing will be drawn.', name);
    end
end

function cnt = countPerUpdate(st, u, edges)
% Spike counts in each lag bin after every update.
% st and u are both sorted, so the count of spikes below each query time
% comes from one histcounts per edge rather than a search per update.
    nU = numel(u);  nE = numel(edges);
    cum = zeros(nU, nE);
    if isempty(st)
        cnt = zeros(nU, nE-1); return;
    end
    st = st(:);
    for e = 1:nE
        q = u + edges(e)/1000;                 % ascending, since u is
        h = histcounts(st, [-Inf; q; Inf]);
        cum(:,e) = cumsum(h(1:end-1))';        % spikes strictly before q_k
    end
    cnt = diff(cum, 1, 2);
end

function m = prep(v, H, W, sd)
% randel vector -> H x W image, optionally smoothed for display only
    m = reshape(v, H, W);
    if sd > 0
        r = ceil(3*sd);
        [gx,gy] = meshgrid(-r:r,-r:r);
        g = exp(-(gx.^2+gy.^2)/(2*sd^2));  g = g/sum(g(:));
        m = conv2(m, g, 'same');
    end
end

function cm = divergingMap(n)
    if nargin < 1, n = 256; end
    m = floor(n/2);
    lo = [linspace(0.02,1,m)' linspace(0.19,1,m)' linspace(0.38,1,m)'];
    hi = [linspace(1,0.70,n-m)' linspace(1,0.09,n-m)' linspace(1,0.17,n-m)'];
    cm = [lo; hi];
end