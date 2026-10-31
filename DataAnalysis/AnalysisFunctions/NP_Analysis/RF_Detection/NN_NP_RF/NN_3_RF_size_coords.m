%% rf_fit.m — precise RF centres and sizes from the NeuroNexus map
%
% Reads nnx_rf.mat (n1) and fits a 2D Gaussian to each channel's receptive
% field. Outputs rfX, rfY and diameter with confidence intervals.
%
% WHY A FIT RATHER THAN A THRESHOLDED CENTROID
% The usual approach - smooth the map, threshold at some fraction of the
% peak, take second moments of what survives - has three biases that cannot
% be removed afterwards:
%   * the answer depends on the threshold, and there is no principled value
%   * smoothing inflates every size, and the correction is only approximate
%   * the mask is a hard edge on a noisy map, so noise pixels touching the
%     blob get included whole
% A model fit uses every pixel at its own weight, needs no threshold and no
% smoothing, and gives a real error bar.
%
% WHAT IS FITTED
%   z(x,y) = base + A * exp(-0.5*[(xr/sx)^2 + (yr/sy)^2])
% with xr,yr the coordinates rotated by theta. Seven free parameters,
% fitted by weighted least squares on the UNSMOOTHED z map. (n1 saves
% zOn/zOff unsmoothed - it only smooths inside its own figure loop.)
%
% THE ERROR BARS ARE REAL
% Noise is measured from the map's own NEGATIVE lags, averaged exactly the
% way the signal is, so it carries the same averaging. The parameters then
% get a parametric bootstrap: resample the map from the fit plus that noise,
% refit, repeat. The spread of the refits is the confidence interval. No
% assumption that the likelihood is quadratic.
%
% WHAT IS STILL SYSTEMATIC, and no amount of fitting fixes:
%   * FLIP_Y. Sets the SIGN of elevation. Verify it against a stimulus at a
%     known position before trusting rfY.
%   * the stimulus geometry (g.X, g.Y, g.width) from the reconstruction. If
%     the grid's own centre is wrong, every coordinate shifts with it.
%   * eye position. Anaesthetised, but not nailed down.
clear; clc; close all;

%% ---- CONFIG -------------------------------------------------------
NNX_FILE  = 'E:\Data\CJ282\Penetration 3\noisegrid_260908_000959\nnx_rf.mat';

% ---- COORDINATE CONVENTION -----------------------------------------
% n1 writes x_deg/y_deg under its own conventions, and on 301026 those put
% a right-hemisphere dorsal penetration in the UPPER RIGHT quadrant when
% the anatomy demands LOWER LEFT -- a 180 deg rotation, i.e. both axes
% mirrored. Eyes were paralysed, so that is the pipeline, not physiology.
%
% Correct it HERE, explicitly, and let the script check the answer against
% the anatomy rather than you eyeballing it. n2 cannot catch this: it
% applies FLIP_Y to both probes together, so NN and NPX agreeing says
% nothing about the absolute sign.
FLIP_X = false;         % mirror azimuth read from the file
FLIP_Y = true;         % mirror elevation read from the file (ON TOP of n1's)

% What the anatomy predicts. Used only to check the fitted centre.
HEMISPHERE = 'right';  % 'right' | 'left'  -> contralateral hemifield
SURFACE    = 'dorsal'; % 'dorsal' | 'ventral' -> lower | upper field

POLARITY  = 'both';    % 'both' | 'on' | 'off' -- headline fit
PEAK_WIN  = 3;         % lag bins averaged around n1's peak bin
NBOOT     = 200;       % bootstrap draws per channel
MIN_R2    = 0.25;      % fits below this are not receptive fields
MIN_SNR   = 2.0;       % peak amplitude / noise
MAX_ITER  = 4000;

DIAM_RULE  = 2.0;      % protocol: stimulus diameter = 2 x aggregate RF
DIAM_FLOOR = 4.0;      % protocol: never below 4 deg
SIG_TO_DIAM = 2.0;     % RF "diameter" = this many sigma. 2 sigma is the
                       % usual convention and what rf_params used.
SF_CPD      = 1.0;     % the sf the protocols run. Sets the SMALLEST usable
                       % aperture, independently of the RF -- see below.

PLOT_N = 6;            % channels to show fit/data/residual for

%% ---- load ---------------------------------------------------------
R = load(NNX_FILE);
H = double(R.size_h);  W = double(R.size_w);
xC = R.x_deg(:)';      yC = R.y_deg(:)';
nCh = numel(R.channels);
if isequal(xC,1:W)
    error(['Axes are in randels, not degrees - n1 had no stimulus geometry. ' ...
           'Cannot produce rfX/rfY/diameter.']);
end
dx = median(abs(diff(xC)));  dy = median(abs(diff(yC)));
randelDeg = mean([dx dy]);

% Apply the correction to the AXES, so every coordinate downstream -- fits,
% CIs, aggregate, figures -- is already in the corrected frame.
if FLIP_X, xC = -xC; end
if FLIP_Y, yC = -yC; end
[XX,YY] = meshgrid(xC,yC);          % H x W, matching reshape(v,H,W)

fprintf('=== PRECISE RF FIT =========================================\n');
fprintf('%s\n', NNX_FILE);
fprintf('grid %d x %d randels of %.4f deg; azimuth %.1f..%.1f, elevation %.1f..%.1f\n', ...
        H, W, randelDeg, min(xC), max(xC), min(yC), max(yC));
if isfield(R,'flip_y') && R.flip_y
    fprintf('FLIP_Y was TRUE in n1 -- the SIGN of rfY depends on it.\n');
end

nB = numel(R.lag_ms);
lag = R.lag_ms(:)';

%% ---- per-channel fit ----------------------------------------------
P = struct('ch',[],'depth',[],'A',[],'cx',[],'cy',[],'sx',[],'sy',[],'th',[], ...
           'r2',[],'snr',[],'noise',[],'ciX',[],'ciY',[],'ciSx',[],'ciSy',[]);
nTried = 0;

fprintf('\n%-5s %7s %8s %8s %9s %9s %8s %8s %7s\n', ...
        'elec','depth','cx','cy','sx(deg)','sy(deg)','R2','pk/sd','status');
for e = 1:nCh
    b = double(R.peak_bin(e));
    w = max(1,b-floor(PEAK_WIN/2)):min(nB,b+floor(PEAK_WIN/2));
    z = peakMap(R, e, w, POLARITY);                 % nPix x 1, unsmoothed

    % --- noise, measured the same way the signal is -----------------
    % Average the SAME number of consecutive bins, but at negative lags
    % where there can be no response. The SD across those maps is the
    % noise on a peak-window-averaged map, averaging included.
    negIx = find(lag < -R.bin_ms);
    nullSD = NaN;
    if numel(negIx) >= PEAK_WIN
        nulls = [];
        for s0 = 1:(numel(negIx)-PEAK_WIN+1)
            nulls = [nulls, peakMap(R, e, negIx(s0:s0+PEAK_WIN-1), POLARITY)]; %#ok<AGROW>
        end
        nullSD = std(nulls(:));
    end
    if ~isfinite(nullSD) || nullSD <= 0, nullSD = 1; end

    M = reshape(z,H,W);
    [pk,pkIx] = max(M(:));
    amp = pk - median(M(:));
    snrE = amp/nullSD;
    if snrE < MIN_SNR
        fprintf('%-5d %7.0f %8s %8s %9s %9s %8s %8.2f %7s\n', ...
                e, dep(R,e), '-','-','-','-','-', snrE, 'low SNR');
        continue;
    end
    nTried = nTried + 1;

    [fit, r2, win] = fitOne(M, XX, YY, pkIx, nullSD, randelDeg, MAX_ITER);
    if isempty(fit) || r2 < MIN_R2
        fprintf('%-5d %7.0f %8s %8s %9s %9s %8.2f %8.2f %7s\n', ...
                e, dep(R,e), '-','-','-','-', r2, snrE, 'bad fit');
        continue;
    end

    % --- parametric bootstrap --------------------------------------
    mdl = gauss2d(fit, XX(win.r,win.c), YY(win.r,win.c));
    bs = nan(NBOOT,4);
    for q = 1:NBOOT
        sim = mdl + nullSD*randn(size(mdl));
        f2 = fitCore(sim, XX(win.r,win.c), YY(win.r,win.c), fit, nullSD, MAX_ITER);
        if ~isempty(f2)
            bs(q,:) = [f2(2) f2(3) exp(f2(4)) exp(f2(5))];
        end
    end
    ci = @(v) prctile(v(~isnan(v)),[2.5 97.5]);

    P.ch(end+1)=R.channels(e);  P.depth(end+1)=dep(R,e);
    P.A(end+1)=fit(1);  P.cx(end+1)=fit(2);  P.cy(end+1)=fit(3);
    P.sx(end+1)=exp(fit(4)); P.sy(end+1)=exp(fit(5)); P.th(end+1)=fit(6);
    P.r2(end+1)=r2;  P.snr(end+1)=snrE;  P.noise(end+1)=nullSD;
    P.ciX(end+1,:)=ci(bs(:,1));  P.ciY(end+1,:)=ci(bs(:,2));
    P.ciSx(end+1,:)=ci(bs(:,3)); P.ciSy(end+1,:)=ci(bs(:,4));

    fprintf('%-5d %7.0f %8.3f %8.3f %9.3f %9.3f %8.2f %8.2f %7s\n', ...
            e, dep(R,e), fit(2), fit(3), exp(fit(4)), exp(fit(5)), r2, snrE, 'ok');
end
nOK = numel(P.ch);
fprintf('\n%d of %d channels fitted (%d passed SNR)\n', nOK, nCh, nTried);
if nOK < 3
    error('Only %d usable fit(s). Not enough to define an aggregate RF.', nOK);
end

%% ---- remove the randel box ----------------------------------------
% The map is the RF convolved with one randel, a square of width w. That
% convolution adds exactly w^2/12 to the second moment, so it comes out in
% quadrature. Nothing else is subtracted: no smoothing was applied.
sBox = randelDeg/sqrt(12);
sxT = sqrt(max(P.sx(:).^2 - sBox^2, 0));
syT = sqrt(max(P.sy(:).^2 - sBox^2, 0));
dFit = SIG_TO_DIAM*sqrt(sxT.*syT);

% [301026] THE RESOLUTION LIMIT, and why a fitted number is not a size.
% Sparse noise cannot resolve an RF smaller than one randel. Below that the
% map is a single bright randel, and a Gaussian fitted to one pixel returns
% a width that describes the RANDEL, not the receptive field. Earlier this
% script carried those numbers into the aggregate and produced a confident
% 0.357 deg for RFs it had not actually measured.
%
% So: anything whose fitted sigma is not comfortably above the randel box
% is reported as an UPPER BOUND and contributes the resolution limit to the
% aggregate, never the fitted value.
% [301026] The criterion was wrong and far too harsh. It demanded the RF
% DIAMETER exceed two randels (0.8 deg here), which confuses "resolvable"
% with "big". What actually has to be true is that the fitted sigma stands
% clear of the randel's own blur: sigma_meas > sqrt(2)*box puts sigma_true
% above one box, which is the point where deconvolution carries real
% information. With that, 11 of these 15 channels resolve at ~0.3-0.5 deg
% instead of all 15 being thrown away.
dLimit   = SIG_TO_DIAM*sBox;                       % diameter of a point source
resolved = P.sx(:) > 1.45*sBox & P.sy(:) > 1.45*sBox;
dCh      = max(dFit, dLimit);                      % never claim below a point

fprintf('\n--- SIZE ---------------------------------------------------\n');
fprintf('randel %.3f deg; box %.4f removed in quadrature\n', randelDeg, sBox);
fprintf('a point source would fit %.3f deg; sigma must exceed %.3f to resolve\n', ...
        dLimit, 1.45*sBox);
fprintf('channels with a RESOLVED size: %d of %d\n', sum(resolved), nOK);
if any(resolved)
    fprintf('  RF DIAMETER: median %.3f deg, range %.3f..%.3f  (n=%d)\n', ...
            median(dFit(resolved)), min(dFit(resolved)), max(dFit(resolved)), sum(resolved));
    r = median(dFit(resolved))/dLimit;
    fprintf('  that is %.1fx the point-source floor', r);
    if r < 2
        fprintf('  -- within a factor of 2, so read it as\n');
        fprintf('     "about %.1f deg, possibly smaller". Randels of %.2f deg\n', ...
                median(dFit(resolved)), randelDeg/3);
        fprintf('     would turn it into a proper measurement.\n');
    else
        fprintf('  -- comfortably measured.\n');
    end
end
if ~all(resolved)
    fprintf(['  %d channel(s) fitted BELOW the randel blur (sigma %.3f..%.3f).\n' ...
             '  Those are single-randel blobs; they contribute %.3f deg to the\n' ...
             '  aggregate rather than their fitted value.\n'], ...
             sum(~resolved), min(P.sx(~resolved)), max(P.sx(~resolved)), dLimit);
end

%% ---- precision actually achieved -----------------------------------
wX = diff(P.ciX,1,2);  wY = diff(P.ciY,1,2);
fprintf('\ncentre 95%% CI width: median %.3f deg azimuth, %.3f deg elevation\n', ...
        median(wX), median(wY));
fprintf('  = %.2f / %.2f randels\n', median(wX)/randelDeg, median(wY)/randelDeg);
fprintf('sigma  95%% CI width: median %.3f / %.3f deg\n', ...
        median(diff(P.ciSx,1,2)), median(diff(P.ciSy,1,2)));

%% ---- aggregate -----------------------------------------------------
% Centre: inverse-variance weighted, so a sloppy channel cannot drag it.
vx = max((wX/3.92).^2, eps);  vy = max((wY/3.92).^2, eps);   % CI -> variance
cx0 = sum(P.cx(:)./vx)/sum(1./vx);
cy0 = sum(P.cy(:)./vy)/sum(1./vy);
seX = sqrt(1/sum(1./vx));  seY = sqrt(1/sum(1./vy));

off   = hypot(P.cx(:)-cx0, P.cy(:)-cy0);
reach = off + dCh/2;
aggR  = prctile(reach,90);
aggD  = 2*aggR;
stimD = max(DIAM_RULE*aggD, DIAM_FLOOR);

fprintf('\n=== AGGREGATE ==============================================\n');
fprintf('centre            %+.3f +/- %.3f , %+.3f +/- %.3f deg\n', cx0, seX, cy0, seY);
fprintf('eccentricity      %.2f deg\n', hypot(cx0,cy0));
% Split the aggregate so you can see WHICH term sets it. If the scatter of
% centres dominates, the number is about the penetration's spread through
% cortex. If the RF radius dominates and the RFs were unresolved, the
% number is really the resolution limit wearing a disguise.
fprintf('\naggregate %.2f deg is built from:\n', aggD);
fprintf('  centre spread (90th pct)  %.3f deg\n', prctile(off,90));
fprintf('  RF radius used            %.3f deg  %s\n', median(dCh)/2, ...
        tern(all(resolved),'(measured)','(RESOLUTION LIMIT for some/all)'));
fprintf('  max over channels         %.2f deg\n', 2*max(reach));
fprintf('\nstimulus diameter %.1f deg  (= max(%.1f x aggregate, %.1f floor))\n', ...
        stimD, DIAM_RULE, DIAM_FLOOR);
% [301026] WHY THE FLOOR IS NOT JUST CAUTION.
% An aperture cannot be shrunk to the receptive field, because the stimulus
% has to remain a GRATING. At SF_CPD cycles/deg a patch of D degrees holds
% D*SF cycles, and below about 3 cycles there is no spatial frequency to
% speak of: orientation and SF selectivity stop being engaged and you are
% showing a blob. That, not RF size, is what sets the smallest usable
% aperture here.
cyc  = stimD*SF_CPD;
dMin = 3/SF_CPD;
fprintf('\naperture vs spatial frequency, at %.2f cpd:\n', SF_CPD);
fprintf('  %.1f deg holds %.1f cycles; 3 cycles needs %.1f deg\n', stimD, cyc, dMin);
if DIAM_RULE*aggD < DIAM_FLOOR
    fprintf(['  -> the FLOOR sets the diameter, not the RF measurement.\n' ...
             '     It is justified: %.1f deg from the RF rule would hold only\n' ...
             '     %.1f cycles. You cannot match the aperture to a %.1f deg\n' ...
             '     RF and still show a grating.\n'], ...
             DIAM_RULE*aggD, DIAM_RULE*aggD*SF_CPD, median(dCh));
    fprintf(['     CONSEQUENCE for exp1C: the disc necessarily contains the\n' ...
             '     RF plus near surround. The annulus still starts at %.1f deg\n' ...
             '     (innerMult 1.5), well outside the RF, so the test survives --\n' ...
             '     but it compares "RF + near surround" against "far surround",\n' ...
             '     not "RF" against "everything else". Say so in the methods.\n'], ...
             0.75*stimD);
end
if cyc < 3
    fprintf('*** %.1f cycles is below 3. Raise diameter to %.1f or sf to %.2f.\n', ...
            cyc, dMin, 3/stimD);
end

%% ---- does the centre agree with the anatomy? -----------------------
switch lower(HEMISPHERE), case 'right', wantX=-1; otherwise, wantX=+1; end
switch lower(SURFACE),    case 'dorsal', wantY=-1; otherwise, wantY=+1; end
qs = @(s) tern(s>0,'upper','lower');
fprintf('\n--- ANATOMY CHECK ------------------------------------------\n');
fprintf('%s hemisphere, %s surface -> expect %s %s quadrant\n', ...
        HEMISPHERE, SURFACE, qs(wantY), tern(wantX<0,'left','right'));
fprintf('fitted centre is %s %s  (FLIP_X %d, FLIP_Y %d applied)\n', ...
        qs(cy0), tern(cx0<0,'left','right'), FLIP_X, FLIP_Y);
if sign(cx0)==wantX && sign(cy0)==wantY
    fprintf('AGREES.\n');
else
    fprintf(['*** DISAGREES. Do not paste this into a protocol until a\n' ...
             '    grating at (%+.2f,%+.2f) has been shown to drive the\n' ...
             '    column and one at (%+.2f,%+.2f) has been shown not to.\n'], ...
             cx0, cy0, -cx0, -cy0);
end

% the stimulus must sit inside what was actually mapped
margin = min([cx0-min(xC), max(xC)-cx0, cy0-min(yC), max(yC)-cy0]);
if stimD/2 > margin
    fprintf(['*** stimulus radius %.2f exceeds the %.2f deg to the nearest\n' ...
             '    grid edge. The RF may continue past the mapped area, so\n' ...
             '    this diameter is a LOWER bound.\n'], stimD/2, margin);
end

%% ---- ON vs OFF, as an independent check ----------------------------
if strcmpi(POLARITY,'both')
    cOn  = quickCentres(R,H,W,XX,YY,PEAK_WIN,'on', MIN_SNR,MIN_R2,randelDeg,MAX_ITER);
    cOff = quickCentres(R,H,W,XX,YY,PEAK_WIN,'off',MIN_SNR,MIN_R2,randelDeg,MAX_ITER);
    [~,ia,ib] = intersect(cOn.ch, cOff.ch);
    if numel(ia) >= 3
        sep = hypot(cOn.cx(ia)-cOff.cx(ib), cOn.cy(ia)-cOff.cy(ib));
        fprintf('\nON vs OFF centres: median %.3f deg apart over %d channels\n', ...
                median(sep), numel(ia));
        fprintf('  %s\n', tern(median(sep) < median(dCh), ...
            'inside one RF -- the two polarities agree', ...
            '*** further apart than an RF. One polarity may be noise.'));
    end
end

%% ---- what to paste -------------------------------------------------
fprintf('\n=== PASTE INTO THE PROTOCOLS ===============================\n');
fprintf('''rfX'',%.2f,''rfY'',%.2f,''diameter'',%.1f\n', ...
        round(cx0*100)/100, round(cy0*100)/100, round(stimD*2)/2);

%% ---- figures -------------------------------------------------------
[~,ord] = sort(P.snr,'descend');
show = ord(1:min(PLOT_N,nOK));
figure('Color','w','Position',[40 40 1400 240*numel(show)/2],'Name','fits','NumberTitle','off');
for i = 1:numel(show)
    k = show(i);  e = find(R.channels==P.ch(k),1);
    bb = double(R.peak_bin(e));
    ww = max(1,bb-floor(PEAK_WIN/2)):min(nB,bb+floor(PEAK_WIN/2));
    M  = reshape(peakMap(R,e,ww,POLARITY),H,W);
    p  = [P.A(k) P.cx(k) P.cy(k) log(P.sx(k)) log(P.sy(k)) P.th(k) 0];
    F  = gauss2d(p, XX, YY);  F = F + (median(M(:))-median(F(:)));
    cl = max(abs(M(:)));
    subplot(numel(show),3,3*i-2);
    imagesc(xC,yC,M,[-cl cl]); axis image xy; title(sprintf('E%d data',P.ch(k)));
    subplot(numel(show),3,3*i-1);
    imagesc(xC,yC,F,[-cl cl]); axis image xy;
    title(sprintf('fit R^2=%.2f',P.r2(k)));
    subplot(numel(show),3,3*i);
    imagesc(xC,yC,M-F,[-cl cl]); axis image xy; title('residual');
end

figure('Color','w','Position',[60 60 760 640],'Name','aggregate','NumberTitle','off');
hold on;
for i = 1:nOK
    [ex,ey] = ell(P.cx(i),P.cy(i),sxT(i),syT(i),P.th(i),SIG_TO_DIAM/2);
    plot(ex,ey,'-','Color',[0.85 0.2 0.15 0.35]);
end
errorbar(P.cx,P.cy,P.cy(:)'-P.ciY(:,1)',P.ciY(:,2)'-P.cy(:)', ...
         P.cx(:)'-P.ciX(:,1)',P.ciX(:,2)'-P.cx(:)','.','Color',[0.4 0.4 0.4]);
scatter(P.cx,P.cy,40,P.depth,'filled');
[ex,ey] = ell(cx0,cy0,aggR,aggR,0,1);      plot(ex,ey,'k-','LineWidth',2);
[ex,ey] = ell(cx0,cy0,stimD/2,stimD/2,0,1); plot(ex,ey,'b--','LineWidth',2);
plot(cx0,cy0,'kp','MarkerFaceColor','y','MarkerSize',18);
axis equal; grid on; box on;
xlim([min(xC) max(xC)]); ylim([min(yC) max(yC)]);
xlabel('azimuth (deg)'); ylabel('elevation (deg)');
cb = colorbar; cb.Label.String = 'depth (\mum)';
title(sprintf('aggregate %.2f deg (black), stimulus %.1f deg (blue)', aggD, stimD));

%% ===================== LOCAL FUNCTIONS =============================
function v = peakMap(R, e, w, pol)
    switch lower(pol)
        case 'on',  v = mean(double(R.zOn (:,w,e)),2);
        case 'off', v = mean(double(R.zOff(:,w,e)),2);
        otherwise,  v = (mean(double(R.zOn(:,w,e)),2) + ...
                         mean(double(R.zOff(:,w,e)),2))/2;
    end
end

function d = dep(R,e)
    if isfield(R,'depth_um'), d = R.depth_um(e); else, d = e; end
end

function m = gauss2d(p, X, Y)
% p = [A cx cy log(sx) log(sy) theta base]. sigma through a log so the
% optimiser cannot walk it negative.
    A=p(1); cx=p(2); cy=p(3); sx=exp(p(4)); sy=exp(p(5)); th=p(6); b=p(7);
    xr =  (X-cx)*cos(th) + (Y-cy)*sin(th);
    yr = -(X-cx)*sin(th) + (Y-cy)*cos(th);
    m = b + A*exp(-0.5*((xr/sx).^2 + (yr/sy).^2));
end

function [p, r2, win] = fitOne(M, XX, YY, pkIx, sd, randelDeg, maxIter)
% Fit on a window around the peak, then widen once if the RF turned out
% bigger than the window could hold.
    [H,W] = size(M);
    [pr,pc] = ind2sub([H W], pkIx);
    Rw = 12;
    for pass = 1:2
        rr = max(1,pr-Rw):min(H,pr+Rw);
        cc = max(1,pc-Rw):min(W,pc+Rw);
        Z = M(rr,cc);  X = XX(rr,cc);  Y = YY(rr,cc);
        p0 = [max(Z(:))-median(Z(:)), XX(pkIx), YY(pkIx), ...
              log(2*randelDeg), log(2*randelDeg), 0, median(Z(:))];
        p = fitCore(Z, X, Y, p0, sd, maxIter);
        if isempty(p), r2 = NaN; win = struct('r',rr,'c',cc); return; end
        need = ceil(4*max(exp(p(4)),exp(p(5)))/randelDeg);
        if pass==1 && need > Rw, Rw = min(need, floor(min(H,W)/2)); continue; end
        break;
    end
    F  = gauss2d(p, X, Y);
    ss = sum((Z(:)-F(:)).^2);  st = sum((Z(:)-mean(Z(:))).^2);
    r2 = 1 - ss/max(st,eps);
    win = struct('r',rr,'c',cc);
end

function p = fitCore(Z, X, Y, p0, sd, maxIter)
% Weighted least squares. Nelder-Mead so nothing outside base MATLAB is
% needed; the initial guess is close enough that it is reliable here.
    cost = @(p) sum(((gauss2d(p,X,Y) - Z)./sd).^2, 'all');
    o = optimset('MaxFunEvals',maxIter,'MaxIter',maxIter, ...
                 'TolX',1e-6,'TolFun',1e-6,'Display','off');
    try
        p = fminsearch(cost, p0, o);
    catch
        p = []; return;
    end
    % reject runaways: a centre outside the window or an absurd sigma
    if p(2) < min(X(:)) || p(2) > max(X(:)) || ...
       p(3) < min(Y(:)) || p(3) > max(Y(:)) || ...
       exp(p(4)) > range(X(:)) || exp(p(5)) > range(Y(:))
        p = [];
    end
end

function C = quickCentres(R,H,W,XX,YY,pw,pol,minSNR,minR2,randelDeg,maxIter)
    C = struct('ch',[],'cx',[],'cy',[]);
    nB = numel(R.lag_ms);
    for e = 1:numel(R.channels)
        b = double(R.peak_bin(e));
        w = max(1,b-floor(pw/2)):min(nB,b+floor(pw/2));
        M = reshape(peakMap(R,e,w,pol),H,W);
        [pk,pkIx] = max(M(:));
        sd = 1;
        if (pk-median(M(:)))/sd < minSNR, continue; end
        [p,r2] = fitOne(M,XX,YY,pkIx,sd,randelDeg,maxIter);
        if isempty(p) || r2 < minR2, continue; end
        C.ch(end+1)=R.channels(e); C.cx(end+1)=p(2); C.cy(end+1)=p(3);
    end
end

function [x,y] = ell(cx,cy,sx,sy,th,k)
    t = linspace(0,2*pi,128);
    x = cx + k*(sx*cos(t)*cos(th) - sy*sin(t)*sin(th));
    y = cy + k*(sx*cos(t)*sin(th) + sy*sin(t)*cos(th));
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end