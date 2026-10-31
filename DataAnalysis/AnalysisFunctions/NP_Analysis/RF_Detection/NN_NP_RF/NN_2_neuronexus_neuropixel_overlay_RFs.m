%% ====================================================================
%  n2_rf_overlay.m  —  DO THE TWO PROBES SEE THE SAME VISUAL SPACE?
%
%  Overlays the receptive fields from the Neuronexus probe (n1) and the
%  Neuropixels probe (s4) in the same visual coordinates.
%
%  The comparison is only meaningful because both probes recorded the SAME
%  presentation, through the same photodiode, with the same noise grid.
%  The randel grid is therefore physically identical, so the two z maps
%  are directly comparable pixel by pixel -- no registration, no assumed
%  offset, nothing fitted between them. If they had been separate sessions
%  any apparent shift could just as easily be eye position or a different
%  stimulus centre, and none of what follows would mean anything.
%
%  Three ways of asking whether they overlap, from coarse to fine:
%
%   1. CENTRES. Where each probe's receptive fields sit, and how far apart
%      the two clouds are relative to their own scatter.
%   2. COVERAGE. For each probe, the fraction of its responsive channels
%      covering each randel. Overlaid as contours, this shows the region
%      of visual space each probe samples, not just its centre.
%   3. OVERLAP INDEX. Jaccard (intersection over union) of the two
%      coverage regions, with the amount expected by chance if the same
%      areas were placed independently, so the number has a reference.
%
%  Requires nnx_rf.mat (n1) and noisegrid_rf.mat (s4 with SAVE_MAPS).
% ====================================================================

clear; clc; close all;

%% ==================== CONFIG =======================================
NNX_FILE = 'D:\Data\CJ283\noisegrid_260916_110717\nnx_rf.mat';
NPX_FILE = ['D:\Data\CJ283\CJ283_noisegrid_110713_2026-09-16_11-07-18\' ...
             'Record Node 101\experiment1\recording1\noisegrid_rf.mat'];

% Restrict either probe to a span of the shank, by DEPTH RANK: 1 = the site
% nearest the tip, ascending toward the shallowest. Rank rather than
% channel number, because consecutive Neuropixels channels sit 3840 um
% apart on this bank -- 170:384 in channel numbers is not a span of tissue.
NPX_DEPTH_SHOW = [];   % [] = all, or e.g. 170:384
NNX_DEPTH_SHOW = [];   % [] = all

SNR_MIN   = 1.;      % channels below this are not plotted
Z_ABS     = 1;         % threshold on the smoothed map
Z_FRAC    = 0.5;
SMOOTH_SD = 1.0;
PEAK_WIN  = 3;
MIN_PIX   = 1;
POLARITY  = 'both';    % 'both' | 'on' | 'off'
COVER_FRAC = 0.2;      % a randel is "covered" if this fraction of a
                       % probe's responsive channels have it above threshold
FLIP_Y = true;         % reverses the row -> elevation mapping for BOTH
                       % probes, applied on top of whatever n1 and s4
                       % saved. Since it flips both together it mirrors the
                       % figure without changing their relative positions,
                       % so the separation and overlap numbers are unchanged.

%% ==================== LOAD =========================================
A = load(NNX_FILE);
B = load(NPX_FILE);
fprintf('neuronexus: %d channels\n', numel(A.channels));
fprintf('neuropixels: %d channels\n', numel(B.channels));

assert(isequal(double(A.size_h),double(B.size_h)) && ...
       isequal(double(A.size_w),double(B.size_w)), ...
    'The two files use different grids (%dx%d vs %dx%d).', ...
    A.size_h, A.size_w, B.size_h, B.size_w);
H = double(A.size_h);  W = double(A.size_w);

xC = A.x_deg(:)';  yC = A.y_deg(:)';
if ~isequal(round(xC,4), round(B.x_deg(:)',4)) || ...
   ~isequal(round(yC,4), round(B.y_deg(:)',4))
    warning(['The two files carry different axis vectors. If FLIP_Y was set ' ...
             'differently in n1 and s4 the overlay will be mirrored. Check ' ...
             'before believing any vertical offset.']);
end
unitStr = 'deg';
if isequal(xC, 1:W), unitStr = 'randels'; end
if FLIP_Y, yC = fliplr(yC); end

%% ==================== DEPTH RANKS AND SELECTION ====================
% n1 saves depth_um; s4 does not, so the Neuropixels depths are read from
% ELECTRODE_YPOS in the recording folder that NPX_FILE sits in. Without it
% a "depth index" would silently mean channel number, which on this probe
% is not depth order at all.
if isfield(A,'depth_um'), depA = A.depth_um(:);
else, depA = (1:numel(A.channels))';  end

if isfield(B,'depth_um')
    depB = B.depth_um(:);
else
    yp = readYposO(fileparts(NPX_FILE));
    if ~isempty(yp) && numel(yp) >= max(B.channels)
        depB = yp(B.channels)';
        fprintf('neuropixels depths from settings.xml: %.0f to %.0f um\n', ...
                min(depB), max(depB));
    else
        depB = (1:numel(B.channels))';
        warning(['Could not read ELECTRODE_YPOS; using channel order as depth. ' ...
                 'NPX_DEPTH_SHOW will select channel numbers, not a span of ' ...
                 'tissue.']);
    end
end
A.depth_um = depA;  B.depth_um = depB;

rankA = zeros(numel(depA),1);  [~,oA] = sort(depA);  rankA(oA) = 1:numel(depA);
rankB = zeros(numel(depB),1);  [~,oB] = sort(depB);  rankB(oB) = 1:numel(depB);

keepA = true(numel(depA),1);
if ~isempty(NNX_DEPTH_SHOW), keepA = ismember(rankA, NNX_DEPTH_SHOW); end
keepB = true(numel(depB),1);
if ~isempty(NPX_DEPTH_SHOW), keepB = ismember(rankB, NPX_DEPTH_SHOW); end

if ~isempty(NNX_DEPTH_SHOW)
    fprintf('neuronexus restricted to depth rank %d-%d: %d of %d channels, %.0f to %.0f um\n', ...
            min(NNX_DEPTH_SHOW), max(NNX_DEPTH_SHOW), sum(keepA), numel(keepA), ...
            min(depA(keepA)), max(depA(keepA)));
end
if ~isempty(NPX_DEPTH_SHOW)
    fprintf('neuropixels restricted to depth rank %d-%d: %d of %d channels, %.0f to %.0f um\n', ...
            min(NPX_DEPTH_SHOW), max(NPX_DEPTH_SHOW), sum(keepB), numel(keepB), ...
            min(depB(keepB)), max(depB(keepB)));
end
assert(any(keepA) && any(keepB), 'A depth filter selected no channels.');

%% ==================== FIT BOTH PROBES ==============================
[Pa, covA] = fitProbe(A, H, W, xC, yC, SNR_MIN, Z_ABS, Z_FRAC, SMOOTH_SD, ...
                      PEAK_WIN, MIN_PIX, POLARITY, keepA);
[Pb, covB] = fitProbe(B, H, W, xC, yC, SNR_MIN, Z_ABS, Z_FRAC, SMOOTH_SD, ...
                      PEAK_WIN, MIN_PIX, POLARITY, keepB);
fprintf('\nfitted: %d neuronexus, %d neuropixels channels\n', ...
        numel(Pa.ch), numel(Pb.ch));
assert(~isempty(Pa.ch) && ~isempty(Pb.ch), ...
    'One probe has no fittable receptive field; nothing to overlay.');

%% ==================== CENTRES ======================================
ma = [mean(Pa.cx) mean(Pa.cy)];  mb = [mean(Pb.cx) mean(Pb.cy)];
sa = [std(Pa.cx)  std(Pa.cy)];   sb = [std(Pb.cx)  std(Pb.cy)];
sep = hypot(ma(1)-mb(1), ma(2)-mb(2));
pooled = hypot(mean([sa(1) sb(1)]), mean([sa(2) sb(2)]));

fprintf('\n--- centres (%s) ---\n', unitStr);
fprintf('neuronexus  %+7.2f, %+7.2f  (SD %.2f, %.2f)\n', ma(1), ma(2), sa(1), sa(2));
fprintf('neuropixels %+7.2f, %+7.2f  (SD %.2f, %.2f)\n', mb(1), mb(2), sb(1), sb(2));
fprintf('separation %.2f %s = %.2f x the within-probe scatter\n', ...
        sep, unitStr, sep/max(pooled,eps));
if sep < pooled
    fprintf('-> the two clouds sit inside each other''s spread: same region\n');
elseif sep < 2*pooled
    fprintf('-> partially separated: neighbouring but offset regions\n');
else
    fprintf('-> clearly separated: the probes are sampling different visual space\n');
end

%% ==================== COVERAGE AND OVERLAP =========================
ca = covA >= COVER_FRAC;
cb2 = covB >= COVER_FRAC;
inter = nnz(ca & cb2);  uni = nnz(ca | cb2);
jac = inter/max(uni,1);
% chance: same areas, independent placement across the grid
pa = nnz(ca)/(H*W);  pb = nnz(cb2)/(H*W);
jacChance = (pa*pb)/max(pa+pb-pa*pb, eps);

fprintf('\n--- coverage (>= %.0f%% of a probe''s responsive channels) ---\n', 100*COVER_FRAC);
fprintf('neuronexus covers %d randels (%.1f%% of the grid)\n', nnz(ca), 100*pa);
fprintf('neuropixels covers %d randels (%.1f%% of the grid)\n', nnz(cb2), 100*pb);
fprintf('intersection %d, union %d\n', inter, uni);
fprintf('Jaccard overlap %.3f (chance for these areas %.3f)\n', jac, jacChance);
if inter == 0
    fprintf('-> NO overlap at all\n');
elseif jac > 3*jacChance
    fprintf('-> substantially more overlap than chance: shared visual territory\n');
else
    fprintf('-> overlap is near what independent placement would give\n');
end
if nnz(ca) > 0
    fprintf('%.0f%% of the neuronexus region is also covered by neuropixels\n', ...
            100*inter/nnz(ca));
end
if nnz(cb2) > 0
    fprintf('%.0f%% of the neuropixels region is also covered by neuronexus\n', ...
            100*inter/nnz(cb2));
end

%% ==================== FIGURES ======================================
colA = [0.85 0.20 0.15];   colB = [0.10 0.35 0.80];

figure('Color','w','Position',[40 40 1350 760],'Name','RF overlay','NumberTitle','off');

subplot(2,3,[1 2 4 5]); hold on;
for i = 1:numel(Pa.ch)
    [ex,ey] = ellipseN(Pa.cx(i),Pa.cy(i),Pa.sx(i),Pa.sy(i),Pa.th(i),1.5);
    plot(ex,ey,'-','Color',[colA 0.45],'LineWidth',1.0);
end
for i = 1:numel(Pb.ch)
    [ex,ey] = ellipseN(Pb.cx(i),Pb.cy(i),Pb.sx(i),Pb.sy(i),Pb.th(i),1.5);
    plot(ex,ey,'-','Color',[colB 0.30],'LineWidth',0.8);
end
hA = plot(Pa.cx,Pa.cy,'o','MarkerFaceColor',colA,'MarkerEdgeColor','none','MarkerSize',6);
hB = plot(Pb.cx,Pb.cy,'s','MarkerFaceColor',colB,'MarkerEdgeColor','none','MarkerSize',4);
plot(ma(1),ma(2),'p','MarkerFaceColor',colA,'MarkerEdgeColor','k','MarkerSize',18);
plot(mb(1),mb(2),'p','MarkerFaceColor',colB,'MarkerEdgeColor','k','MarkerSize',18);
plot([ma(1) mb(1)],[ma(2) mb(2)],'k-','LineWidth',1.5);
axis equal; grid on; box on;
xlim([min(xC) max(xC)]); ylim([min(yC) max(yC)]);
xlabel(sprintf('azimuth (%s)',unitStr)); ylabel(sprintf('elevation (%s)',unitStr));
legend([hA hB],{sprintf('neuronexus (n=%d)',numel(Pa.ch)), ...
                sprintf('neuropixels (n=%d)',numel(Pb.ch))},'Location','best');
title(sprintf('receptive fields, both probes | centres %.2f %s apart (%.2fx scatter)', ...
      sep, unitStr, sep/max(pooled,eps)));

subplot(2,3,3); hold on;
imagesc(xC,yC,zeros(H,W)); axis image xy;
contour(xC,yC,covA,[COVER_FRAC COVER_FRAC],'LineColor',colA,'LineWidth',2);
contour(xC,yC,covB,[COVER_FRAC COVER_FRAC],'LineColor',colB,'LineWidth',2);
xlim([min(xC) max(xC)]); ylim([min(yC) max(yC)]);
grid on; box on; xlabel(sprintf('azimuth (%s)',unitStr));
title(sprintf('coverage outlines, Jaccard %.3f', jac));

subplot(2,3,6);
rgb = zeros(H,W,3);
rgb(:,:,1) = covA/max(covA(:)+eps);
rgb(:,:,3) = covB/max(covB(:)+eps);
image(xC,yC,rgb); axis image xy;
xlabel(sprintf('azimuth (%s)',unitStr)); ylabel(sprintf('elevation (%s)',unitStr));
title('red = neuronexus, blue = neuropixels, magenta = both');

figure('Color','w','Position',[70 70 1100 380],'Name','depth','NumberTitle','off');
subplot(1,2,1); hold on;
plot(Pa.depth, Pa.cx, 'o-','Color',colA); 
plot(Pb.depth, Pb.cx, 's','Color',colB,'MarkerSize',3);
grid on; box on; xlabel('distance along probe (\mum)');
ylabel(sprintf('azimuth (%s)',unitStr));
legend('neuronexus','neuropixels','Location','best'); title('azimuth vs position');
subplot(1,2,2); hold on;
plot(Pa.depth, Pa.cy, 'o-','Color',colA);
plot(Pb.depth, Pb.cy, 's','Color',colB,'MarkerSize',3);
grid on; box on; xlabel('distance along probe (\mum)');
ylabel(sprintf('elevation (%s)',unitStr)); title('elevation vs position');

%% ==================== LOCAL FUNCTIONS ==============================
function [P, cover] = fitProbe(R, H, W, xC, yC, snrMin, zAbs, zFrac, sd, pw, minPix, pol, keep)
    n = numel(R.channels);
    if nargin < 13 || isempty(keep), keep = true(n,1); end
    P = struct('ch',[],'depth',[],'cx',[],'cy',[],'sx',[],'sy',[],'th',[]);
    cover = zeros(H,W);  nResp = 0;
    if isfield(R,'depth_um'), dpt = R.depth_um(:); else, dpt = (1:n)'; end
    for k = 1:n
        if ~keep(k), continue; end
        if R.snr(k) < snrMin, continue; end
        b = double(R.peak_bin(k));  nB = numel(R.lag_ms);
        w = max(1,b-floor(pw/2)):min(nB,b+floor(pw/2));
        switch lower(pol)
            case 'on',  v = mean(double(R.zOn(:,w,k)),2);
            case 'off', v = mean(double(R.zOff(:,w,k)),2);
            otherwise,  v = (mean(double(R.zOn(:,w,k)),2)+mean(double(R.zOff(:,w,k)),2))/2;
        end
        m = smoothN2(v,H,W,sd);
        [pk,pkIx] = max(m(:));
        if pk < zAbs, continue; end
        thr = max(zAbs, zFrac*pk);
        mask = growN2(m>=thr, pkIx, H, W);
        if nnz(mask) < minPix, continue; end
        cover = cover + double(mask);  nResp = nResp + 1;
        [rr,cc] = find(mask);  ww = m(mask)-thr;
        if sum(ww) <= 0, continue; end
        x = xC(cc)'; y = yC(rr)';
        cx = sum(ww(:).*x)/sum(ww);  cy = sum(ww(:).*y)/sum(ww);
        C = [sum(ww(:).*(x-cx).^2) sum(ww(:).*(x-cx).*(y-cy)); ...
             0                     sum(ww(:).*(y-cy).^2)]/sum(ww);
        C(2,1) = C(1,2);
        [V,D] = eig(C);  d = sqrt(max(diag(D),0));
        [d,o] = sort(d,'descend');  V = V(:,o);
        P.ch(end+1)    = R.channels(k);
        P.depth(end+1) = dpt(k);
        P.cx(end+1)=cx; P.cy(end+1)=cy;
        P.sx(end+1)=d(1); P.sy(end+1)=d(2);
        P.th(end+1)=atan2(V(2,1),V(1,1));
    end
    if nResp > 0, cover = cover/nResp; end
end

function m = smoothN2(v,H,W,sd)
    m = reshape(v,H,W);
    if sd > 0
        r = ceil(3*sd);
        [gx,gy] = meshgrid(-r:r,-r:r);
        g = exp(-(gx.^2+gy.^2)/(2*sd^2)); g = g/sum(g(:));
        m = conv2(m,g,'same');
    end
end

function mask = growN2(bin,seed,H,W)
    mask = false(H,W);
    if ~bin(seed), return; end
    mask(seed) = true; stack = seed;
    while ~isempty(stack)
        p = stack(end); stack(end) = [];
        r = mod(p-1,H)+1; c = floor((p-1)/H)+1;
        for dr = -1:1
            for dc = -1:1
                rr = r+dr; cc = c+dc;
                if rr>=1&&rr<=H&&cc>=1&&cc<=W&&bin(rr,cc)&&~mask(rr,cc)
                    mask(rr,cc)=true; stack(end+1)=(cc-1)*H+rr; %#ok<AGROW>
                end
            end
        end
    end
end

function ypos = readYposO(startDir)
% ELECTRODE_YPOS from the Open Ephys settings.xml, searched for by walking
% up from the folder holding noisegrid_rf.mat.
    ypos = [];
    d = startDir;
    for up = 1:6
        fx = dir(fullfile(d,'**','settings.xml'));
        for q = 1:numel(fx)
            try, txt = fileread(fullfile(fx(q).folder,fx(q).name)); catch, continue; end
            yT = regexp(txt,'<ELECTRODE_YPOS\s([^/>]*)/?>','tokens','once');
            if isempty(yT), continue; end
            t = regexp(yT{1},'CH(\d+)\s*=\s*"(-?[\d.]+)"','tokens');
            if isempty(t), continue; end
            iy = cellfun(@(a) str2double(a{1}),t)+1;
            vy = cellfun(@(a) str2double(a{2}),t);
            ypos = nan(1,max(iy));  ypos(iy) = vy;  return;
        end
        pd = fileparts(d);
        if isempty(pd) || strcmp(pd,d), break; end
        d = pd;
    end
end

function [x,y] = ellipseN(cx,cy,sx,sy,th,k)
    t = linspace(0,2*pi,64);
    x = cx + k*(sx*cos(t)*cos(th) - sy*sin(t)*sin(th));
    y = cy + k*(sx*cos(t)*sin(th) + sy*sin(t)*cos(th));
end