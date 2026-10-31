%% ====================================================================
%  s5_noisegrid_rf_progression.m  —  CJ280 noisegrid, STEP 5
%
%  Every channel's receptive field drawn in the same visual coordinates,
%  coloured by position along the probe, so that retinotopic progression
%  is visible: whether the RFs move as the probe goes from channel 1
%  (deepest) to channel 384 (shallowest), and in which direction.
%
%  Requires noisegrid_rf.mat, i.e. s3 run over CHANNELS = 1:384 and s4 run
%  with SAVE_MAPS = true.
%
%  How each RF is reduced to an outline
%  ------------------------------------
%  The z map at the channel's own peak lag is thresholded at
%  max(Z_ABS, Z_FRAC x peak), the connected region containing the peak is
%  kept, and a weighted centroid and covariance give an ellipse. Keeping
%  only the connected region matters: with 6400 randels a plain threshold
%  picks up scattered noise pixels far from the RF, which drags the
%  centroid toward the middle of the grid and manufactures a spurious
%  progression toward the centre.
%
%  Channels whose peak-to-null ratio falls below SNR_MIN are excluded
%  rather than drawn faintly. A channel with no receptive field still has
%  a brightest pixel, and plotting it adds a random point that looks like
%  scatter around the trend when it is really nothing at all.
% ====================================================================

clear; clc; close all;

%% ==================== CONFIG =======================================
REC_ROOT = 'D:\Data\CJ283\CJ283_noisegrid_110713_2026-09-16_11-07-18';

CHANNELS  = [1:384];        % [] = all channels in the file
SNR_MIN   = 1;         % peak/null spread ratio required to include a channel

% RF outline
PEAK_WIN   = 3;        % lag bins averaged around the peak
SMOOTH_SD  = 1.0;      % Gaussian SD in randels, applied before the fit
POLARITY   = 'positive';  % 'positive' (facilitation) or 'abs'
Z_ABS      = 1;        % absolute z floor for the threshold. This is applied to
                       % the SMOOTHED map, whose noise SD is well below 1 --
                       % the script reports the measured null level so this can
                       % be set from the data rather than guessed.
Z_FRAC     = 0.5;      % and this fraction of the channel's peak
ELLIPSE_SD = 1.5;      % ellipse drawn at this many SDs of the weighted spread
MIN_PIXELS = 1;        % smaller regions are treated as noise

% probe geometry
USE_PROBE_XML = true;  % read the real electrode positions from the Open Ephys
                       % settings.xml. This is the only reliable source: a
                       % single-column bank does not necessarily step one row
                       % per channel, and assuming the wrong pitch rescales
                       % every deg/mm gradient below.
ROW_PITCH_UM = 20;     % fallback only, used if the XML cannot be read
N_COLS       = 1;      % fallback only

% vertical orientation
FLIP_Y = true;         % reverses the row -> elevation mapping, i.e. treats
                       % randel row 1 as the BOTTOM of the screen rather than
                       % the top. This is applied on top of whatever s4 saved:
                       % if s4 already ran with its own FLIP_Y = true, setting
                       % it true here flips it back.

% display
BIN_CH     = 0;        % >0 averages centres over blocks of this many channels
                       % for the trajectory panel; 0 = every channel
SHOW_ELLIPSES = true;
SAVE_FIGS  = false;
FIG_DIR    = '';

%% ==================== LOAD =========================================
ob = dir(fullfile(REC_ROOT,'**','structure.oebin'));
assert(~isempty(ob),'No structure.oebin under %s',REC_ROOT);
rec_dir = ob(1).folder;
rfFile  = fullfile(rec_dir,'noisegrid_rf.mat');
assert(isfile(rfFile), ...
    ['%s not found. Run s4 with SAVE_MAPS = true (and s3 over ' ...
     'CHANNELS = 1:384) first.'], rfFile);

R = load(rfFile, 'channels','zOn','zOff','lag_ms','bin_ms','peak_bin', ...
                 'peak_lag_ms','snr','size_h','size_w','x_deg','y_deg');

H = double(R.size_h);  W = double(R.size_w);
xC = R.x_deg(:)';  yC = R.y_deg(:)';
if FLIP_Y, yC = fliplr(yC); end
inDeg = ~isequal(xC, 1:W);
unitStr = 'deg'; if ~inDeg, unitStr = 'randels'; end

allCh = R.channels(:)';
if isempty(CHANNELS)
    sel = 1:numel(allCh);
else
    [tf, sel] = ismember(CHANNELS, allCh);
    assert(all(tf), 'Channel(s) %s not in the file.', mat2str(CHANNELS(~tf)));
    sel = sel(:)';
end
nCh = numel(sel);

fprintf('loaded %d channels, grid %dx%d, axes in %s\n', nCh, H, W, unitStr);
fprintf('lags %g to %g ms\n', min(R.lag_ms), max(R.lag_ms));

depth_um = [];
if USE_PROBE_XML
    % Electrode positions from the Open Ephys settings.xml, which records
    % ELECTRODE_XPOS / ELECTRODE_YPOS as attributes CH0..CHn for the sites
    % the bank actually selected. Parsed inline so this block needs no
    % helper function.
    fx = dir(fullfile(REC_ROOT,'**','settings.xml'));
    ypos = [];  xpos = [];  xmlFile = '';
    for q = 1:numel(fx)
        cand = fullfile(fx(q).folder, fx(q).name);
        try, txt = fileread(cand); catch, continue; end
        yTok = regexp(txt,'<ELECTRODE_YPOS\s([^/>]*)/?>','tokens','once');
        if isempty(yTok), continue; end
        tY = regexp(yTok{1},'CH(\d+)\s*=\s*"(-?[\d.]+)"','tokens');
        if isempty(tY), continue; end
        iy = cellfun(@(c) str2double(c{1}), tY) + 1;
        vy = cellfun(@(c) str2double(c{2}), tY);
        ypos = nan(1, max(iy));  ypos(iy) = vy;

        xpos = zeros(size(ypos));
        xTok = regexp(txt,'<ELECTRODE_XPOS\s([^/>]*)/?>','tokens','once');
        if ~isempty(xTok)
            tX = regexp(xTok{1},'CH(\d+)\s*=\s*"(-?[\d.]+)"','tokens');
            if ~isempty(tX)
                jx = cellfun(@(c) str2double(c{1}), tX) + 1;
                vx = cellfun(@(c) str2double(c{2}), tX);
                xpos = nan(1, max(jx));  xpos(jx) = vx;
            end
        end
        xmlFile = cand;  break;
    end

    if ~isempty(ypos) && numel(ypos) >= max(allCh)
        depth_um = ypos(allCh);
        ux = unique(xpos(allCh));
        dy = unique(round(diff(sort(unique(ypos(allCh)))), 3));
        fprintf('\nprobe geometry from %s\n', xmlFile);
        fprintf('  %d distinct x position(s): %s um\n', numel(ux), mat2str(ux(:)',5));
        fprintf('  y pitch between occupied rows: %s um\n', mat2str(dy(:)',5));
        fprintf('  channel 1 at y = %.0f, channel %d at y = %.0f um\n', ...
                depth_um(1), allCh(end), depth_um(end));
        if numel(ux) == 1
            fprintf('  -> single column confirmed\n');
        else
            fprintf('  -> %d columns; channels at the same y are a stagger pair\n', numel(ux));
        end
    else
        warning(['Could not read electrode positions from settings.xml; ' ...
                 'falling back to ROW_PITCH_UM = %g with N_COLS = %d.'], ...
                 ROW_PITCH_UM, N_COLS);
    end
end
if isempty(depth_um)
    depth_um = floor((double(allCh)-1)/N_COLS) * ROW_PITCH_UM;
end
depth_um = double(depth_um(:))';
fprintf('probe span: %.0f to %.0f um from the tip\n', ...
        min(depth_um(sel)), max(depth_um(sel)));

%% ==================== NULL LEVEL OF THE SMOOTHED MAP ===============
% Smoothing correlates neighbouring randels and shrinks the noise, so the
% right Z_ABS depends on SMOOTH_SD and is not 4 just because raw z is unit
% variance. Measure it: negative lags carry no response by construction, so
% the largest smoothed value there is what noise alone produces.
nBinAll = numel(R.lag_ms);
negIx = find(R.lag_ms < 0);
if numel(negIx) >= 2
    nsamp = min(nCh, 60);
    pick  = round(linspace(1, nCh, nsamp));
    nullPk = zeros(nsamp,2);
    for q = 1:nsamp
        k = sel(pick(q));
        for pol = 1:2
            Zp = R.(ternary(pol==1,'zOn','zOff'));
            mm = smoothMap(mean(double(Zp(:,negIx,k)),2), H, W, SMOOTH_SD);
            nullPk(q,pol) = max(mm(:));
        end
    end
    fprintf('\nsmoothed-map noise (negative lags, %d channels sampled):\n', nsamp);
    fprintf('  peak value: median %.2f, 95th pct %.2f, max %.2f\n', ...
            median(nullPk(:)), prctile(nullPk(:),95), max(nullPk(:)));
    fprintf('  Z_ABS = %.2f  (currently %.2f)\n', ...
            max(nullPk(:)), Z_ABS);
    if Z_ABS < prctile(nullPk(:),95)
        warning(['Z_ABS = %.2f is below the 95th percentile of pure noise ' ...
                 '(%.2f). Most channels will produce a blob whether or not ' ...
                 'they have a receptive field.'], Z_ABS, prctile(nullPk(:),95));
    end
end

%% ==================== FIT EACH RF ==================================
nBin = numel(R.lag_ms);
F = struct('ch',{},'pol',{},'cx',{},'cy',{},'sx',{},'sy',{},'theta',{}, ...
           'area',{},'peakz',{},'npix',{},'snr',{},'lag',{},'depth',{});

fprintf('\nfitting ...\n');
for pol = 1:2
    Z = R.(ternary(pol==1,'zOn','zOff'));
    for ii = 1:nCh
        k  = sel(ii);
        if R.snr(k) < SNR_MIN, continue; end
        b   = double(R.peak_bin(k));
        wIx = max(1,b-floor(PEAK_WIN/2)) : min(nBin, b+floor(PEAK_WIN/2));
        m   = smoothMap(mean(double(Z(:,wIx,k)),2), H, W, SMOOTH_SD);

        if strcmpi(POLARITY,'abs'), s = abs(m); else, s = m; end
        [pk, pkIx] = max(s(:));
        thr = max(Z_ABS, Z_FRAC*pk);
        if pk < Z_ABS, continue; end

        mask = growFrom(s >= thr, pkIx, H, W);
        if nnz(mask) < MIN_PIXELS, continue; end

        [rr, cc] = find(mask);
        w  = s(mask) - thr;   w = w(:);
        if sum(w) <= 0, continue; end
        x  = xC(cc)';  y = yC(rr)';
        cx = sum(w.*x)/sum(w);
        cy = sum(w.*y)/sum(w);
        C  = [sum(w.*(x-cx).^2) sum(w.*(x-cx).*(y-cy)); 0 sum(w.*(y-cy).^2)]/sum(w);
        C(2,1) = C(1,2);
        [V,D] = eig(C);
        d = sqrt(max(diag(D),0));
        [d, o] = sort(d,'descend');  V = V(:,o);

        F(end+1) = struct('ch',allCh(k),'pol',pol,'cx',cx,'cy',cy, ...
            'sx',d(1),'sy',d(2),'theta',atan2(V(2,1),V(1,1)), ...
            'area',pi*d(1)*d(2),'peakz',pk,'npix',nnz(mask), ...
            'snr',R.snr(k),'lag',R.peak_lag_ms(k),'depth',depth_um(k)); %#ok<SAGROW>
    end
end

onF  = F([F.pol]==1);
offF = F([F.pol]==2);
fprintf('included: %d ON, %d OFF (of %d channels, SNR_MIN %.2f)\n', ...
        numel(onF), numel(offF), nCh, SNR_MIN);
if isempty(onF) && isempty(offF)
    error(['No channel passed the criteria. Either there are no receptive ' ...
           'fields in this data, or SNR_MIN / Z_ABS are too strict.']);
end

%% ==================== PROGRESSION STATISTICS =======================
fprintf('\n%-5s %10s %8s %8s %8s %8s %8s %7s\n', ...
        'pol','depth(um)','ch','x','y','sx','sy','lag');
for p = 1:2
    G = F([F.pol]==p);
    if isempty(G), continue; end
    [~,o] = sort([G.depth]);  G = G(o);
    step = max(1, round(numel(G)/12));
    for i = 1:step:numel(G)
        fprintf('%-5s %10.0f %8d %8.2f %8.2f %8.2f %8.2f %7.0f\n', ...
            ternary(p==1,'ON','OFF'), G(i).depth, G(i).ch, ...
            G(i).cx, G(i).cy, G(i).sx, G(i).sy, G(i).lag);
    end
end

fprintf('\n--- retinotopic gradient (linear fit vs distance from tip) ---\n');
grad = struct();
for p = 1:2
    G = F([F.pol]==p);  nm = ternary(p==1,'ON','OFF');
    if numel(G) < 5
        fprintf('%-4s too few channels to fit\n', nm); continue;
    end
    dmm = [G.depth]'/1000;
    for ax = 1:2
        v = ternary(ax==1, [G.cx]', [G.cy]');
        [pf,Sf] = polyfit(dmm, v, 1);
        rr = v - polyval(pf,dmm);
        r2 = 1 - sum(rr.^2)/max(sum((v-mean(v)).^2),eps);
        fprintf('%-4s %s: %+7.3f %s/mm  (R2 %.3f, n %d)\n', ...
                nm, ternary(ax==1,'azimuth  ','elevation'), pf(1), unitStr, r2, numel(G));
        grad.(sprintf('%s_%s',lower(nm),ternary(ax==1,'x','y'))) = pf;
    end
end
fprintf(['\nA gradient much larger than its scatter means the RFs shift ' ...
         'systematically along the probe.\n']);

%% ==================== FIGURES ======================================
cmap = parula(256);
dmin = min(depth_um(sel));  dmax = max(depth_um(sel));
col  = @(d) cmap(1+round(255*(d-dmin)/max(dmax-dmin,eps)), :);

if SAVE_FIGS
    if isempty(FIG_DIR), FIG_DIR = fullfile(rec_dir,'rf_figs'); end
    if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end
end

for p = 1:2
    G = F([F.pol]==p);  nm = ternary(p==1,'ON','OFF');
    if isempty(G), continue; end
    [~,o] = sort([G.depth]);  G = G(o);

    fig = figure('Color','w','Position',[40 40 1300 820], ...
                 'Name',sprintf('%s progression',nm),'NumberTitle','off');

    % --- all RFs overlaid in visual space -------------------------
    subplot(2,3,[1 2 4 5]); hold on;
    for i = 1:numel(G)
        c = col(G(i).depth);
        if SHOW_ELLIPSES
            [ex,ey] = ellipsePts(G(i).cx,G(i).cy,G(i).sx,G(i).sy,G(i).theta,ELLIPSE_SD);
            plot(ex, ey, '-', 'Color', [c 0.55], 'LineWidth', 1.0);
        end
        plot(G(i).cx, G(i).cy, '.', 'Color', c, 'MarkerSize', 9);
    end
    % centre trajectory, optionally averaged in blocks of channels
    if BIN_CH > 0
        e = dmin:(BIN_CH*ROW_PITCH_UM/N_COLS):dmax+1;
        bx = nan(numel(e)-1,1); by = bx;
        for i = 1:numel(e)-1
            in = [G.depth] >= e(i) & [G.depth] < e(i+1);
            if any(in), bx(i) = mean([G(in).cx]); by(i) = mean([G(in).cy]); end
        end
        ok = ~isnan(bx);
        plot(bx(ok), by(ok), 'k-', 'LineWidth', 2);
    else
        plot([G.cx], [G.cy], 'k-', 'LineWidth', 0.8);
    end
    plot(G(1).cx,   G(1).cy,   'ks','MarkerFaceColor','w','MarkerSize',10);
    plot(G(end).cx, G(end).cy, 'kp','MarkerFaceColor','w','MarkerSize',14);
    axis equal; grid on; box on;
    xlim([min(xC) max(xC)]); ylim([min(yC) max(yC)]);
    xlabel(sprintf('azimuth (%s)',unitStr)); ylabel(sprintf('elevation (%s)',unitStr));
    title(sprintf(['%s receptive fields, %d channels\n' ...
                   'square = nearest the tip (ch %d), star = furthest (ch %d)'], ...
                   nm, numel(G), G(1).ch, G(end).ch));
    cb = colorbar; colormap(gca,cmap); caxis([dmin dmax]);
    cb.Label.String = 'distance from probe tip (\mum)';

    % --- azimuth vs depth -----------------------------------------
    subplot(2,3,3); hold on;
    scatter([G.depth],[G.cx],18,[G.depth],'filled'); colormap(gca,cmap); caxis([dmin dmax]);
    if numel(G) >= 5
        dd = linspace(min([G.depth]),max([G.depth]),50);
        pf = polyfit([G.depth]'/1000,[G.cx]',1);
        plot(dd, polyval(pf,dd/1000), 'r-','LineWidth',1.4);
        title(sprintf('azimuth: %+.3f %s/mm', pf(1), unitStr));
    end
    grid on; box on; xlabel('from tip (\mum)'); ylabel(sprintf('azimuth (%s)',unitStr));

    % --- elevation vs depth ---------------------------------------
    subplot(2,3,6); hold on;
    scatter([G.depth],[G.cy],18,[G.depth],'filled'); colormap(gca,cmap); caxis([dmin dmax]);
    if numel(G) >= 5
        pf = polyfit([G.depth]'/1000,[G.cy]',1);
        plot(dd, polyval(pf,dd/1000), 'r-','LineWidth',1.4);
        title(sprintf('elevation: %+.3f %s/mm', pf(1), unitStr));
    end
    grid on; box on; xlabel('from tip (\mum)'); ylabel(sprintf('elevation (%s)',unitStr));

    if SAVE_FIGS
        exportgraphics(fig, fullfile(FIG_DIR,sprintf('progression_%s.png',nm)), ...
                       'Resolution',150);
    end
end

%% ==================== ON vs OFF AND RF SIZE ========================
figure('Color','w','Position',[60 60 1250 420],'Name','ON vs OFF','NumberTitle','off');

subplot(1,3,1); hold on;
if ~isempty(onF),  plot([onF.cx], [onF.cy], 'o','Color',[0.85 0.15 0.15],'MarkerSize',4); end
if ~isempty(offF), plot([offF.cx],[offF.cy],'s','Color',[0.10 0.30 0.80],'MarkerSize',4); end
axis equal; grid on; box on;
xlim([min(xC) max(xC)]); ylim([min(yC) max(yC)]);
xlabel(sprintf('azimuth (%s)',unitStr)); ylabel(sprintf('elevation (%s)',unitStr));
legend('ON','OFF','Location','best'); title('ON and OFF centres');

subplot(1,3,2); hold on;
if ~isempty(onF),  plot([onF.depth], 2*ELLIPSE_SD*sqrt([onF.sx].*[onF.sy]), 'o','MarkerSize',4); end
if ~isempty(offF), plot([offF.depth],2*ELLIPSE_SD*sqrt([offF.sx].*[offF.sy]),'s','MarkerSize',4); end
grid on; box on; xlabel('from tip (\mum)');
ylabel(sprintf('RF diameter (%s)',unitStr));
legend('ON','OFF','Location','best'); title('RF size along the probe');

subplot(1,3,3); hold on;
if ~isempty(onF),  plot([onF.depth], [onF.lag], 'o','MarkerSize',4); end
if ~isempty(offF), plot([offF.depth],[offF.lag],'s','MarkerSize',4); end
grid on; box on; xlabel('from tip (\mum)'); ylabel('peak lag (ms)');
legend('ON','OFF','Location','best'); title('response latency along the probe');

%% ==================== LOCAL FUNCTIONS ==============================
function [ypos, xpos, xmlFile] = readProbePositions(root)
% Electrode positions from the Open Ephys settings.xml, which records
% ELECTRODE_XPOS / ELECTRODE_YPOS attributes CH0..CHn for the sites actually
% selected. Returns [] if it cannot be read, so the caller can fall back.
    ypos = []; xpos = []; xmlFile = '';
    f = dir(fullfile(root,'**','settings.xml'));
    if isempty(f), return; end
    for q = 1:numel(f)
        xmlFile = fullfile(f(q).folder, f(q).name);
        try, txt = fileread(xmlFile); catch, continue; end
        yTok = regexp(txt,'<ELECTRODE_YPOS\s([^/>]*)/?>','tokens','once');
        xTok = regexp(txt,'<ELECTRODE_XPOS\s([^/>]*)/?>','tokens','once');
        if isempty(yTok), continue; end
        ypos = parseCH(yTok{1});
        if ~isempty(xTok), xpos = parseCH(xTok{1}); else, xpos = zeros(size(ypos)); end
        if ~isempty(ypos), return; end
    end
    ypos = []; xpos = []; xmlFile = '';
end

function v = parseCH(attrStr)
% 'CH0="0" CH1="0" CH2="20" ...' -> vector indexed from channel 1
    t = regexp(attrStr,'CH(\d+)\s*=\s*"(-?[\d.]+)"','tokens');
    if isempty(t), v = []; return; end
    idx = cellfun(@(c) str2double(c{1}), t) + 1;
    val = cellfun(@(c) str2double(c{2}), t);
    v = nan(1, max(idx));
    v(idx) = val;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function m = smoothMap(v, H, W, sd)
    m = reshape(v, H, W);
    if sd > 0
        r = ceil(3*sd);
        [gx,gy] = meshgrid(-r:r,-r:r);
        g = exp(-(gx.^2+gy.^2)/(2*sd^2));  g = g/sum(g(:));
        m = conv2(m, g, 'same');
    end
end

function mask = growFrom(bin, seed, H, W)
% 8-connected flood fill from the peak pixel. Written out rather than using
% bwconncomp so the script does not need the Image Processing Toolbox.
    mask = false(H,W);
    if ~bin(seed), return; end
    mask(seed) = true;
    stack = seed;
    while ~isempty(stack)
        p = stack(end);  stack(end) = [];
        r = mod(p-1,H)+1;  c = floor((p-1)/H)+1;
        for dr = -1:1
            for dc = -1:1
                rr = r+dr;  cc = c+dc;
                if rr>=1 && rr<=H && cc>=1 && cc<=W && bin(rr,cc) && ~mask(rr,cc)
                    mask(rr,cc) = true;
                    stack(end+1) = (cc-1)*H + rr; %#ok<AGROW>
                end
            end
        end
    end
end

function [x,y] = ellipsePts(cx,cy,sx,sy,th,k)
    t = linspace(0,2*pi,64);
    x = cx + k*(sx*cos(t)*cos(th) - sy*sin(t)*sin(th));
    y = cy + k*(sx*cos(t)*sin(th) + sy*sin(t)*cos(th));
end