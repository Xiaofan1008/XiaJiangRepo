%% CHECK 4: raw-trace viewer around the trigger (script)
% Self-contained (no cd, no loadTrig, no Depth_s). Reads amplifier.dat
% (int16 x 0.195 uV), the trigger from digitalin.dat and the conditions from
% the exp_datafile in the folder. Only plots; nothing is saved.
%
% One FIGURE per amplitude level x type (sim / seq / ctrl) x channel, with one
% subplot per condition. The two orders of a sequential pair (A->B k/10 and
% B->A k/10) can be overlaid in one subplot (mergeOrders): red = A->B,
% blue = B->A. Grey/black = single condition.
% Small dark triangles at the top = the pulse times StimParams says should occur.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
intanFolder = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_27/Train_Amp_Test2_261027_192530';
fs      = 30000;
nCh     = 32;
trigBit = 0;

% ---- what to plot ----
plotMode    = 'by_channel';  % 'by_channel': figure per (amp, type, channel), subplot per condition
                             % 'by_condition': figure per condition (pair), subplot per channel
chans       = [32];        % amplifier rows, or 'all'  (row r = A-(r-1) if A-000..A-031 order)
types       = {'sim','seq'}; % any of 'sim', 'seq', 'ctrl' (ctrl = visual / blank, no pulse)
ampOnly     = [];            % amplitude levels, e.g. [50]. [] = all
conds       = [];            % condition IDs from the printed table. [] = all
labelHas    = '';            % regexp on the label. '' = no filter
mergeOrders = true;          % true: overlay A->B k/n and B->A k/n in one subplot
nTrialsPlot = 8;             % trials per condition
maxFigures  = 12;            % safety cap

% ---- window and display ----
winMs        = [-10 80];     % ms around the trigger
ylimUV       = [];           % [] = auto, or e.g. [-6000 6000]
showExpected = true;
pulseStartMs = 0;            % pulses start at the trigger in this setup (checked on Train_Amp_test2)
figSize      = [80 60 1500 900];
%% ================================

%% ---- load time, trigger, conditions ----
fid = fopen(fullfile(intanFolder,'time.dat'),'r');      t = fread(fid,Inf,'int32');  fclose(fid);
fid = fopen(fullfile(intanFolder,'digitalin.dat'),'r'); D = fread(fid,Inf,'uint16'); fclose(fid);
N = numel(t);
brk = find(diff(t) ~= 1);  segStart = [1; brk+1];  segEnd = [brk; N];
up = find(diff(bitget(D,trigBit+1))==1) + 1;  clear D

f = dir(fullfile(intanFolder,'*_exp_datafile_*.mat'));
assert(~isempty(f), 'No *_exp_datafile_*.mat in the folder.');
S = load(fullfile(intanFolder,f(1).name),'StimParams','TrialTable');
nT = min(numel(up), numel(S.TrialTable));
cnd = double([S.TrialTable(1:nT).cond])';
lab = cellstr({S.TrialTable(1:nT).label}');

amp = zeros(nT,2);  ptd = zeros(nT,2);  nsp = zeros(nT,2);  per = zeros(nT,1);
for n = 1:nT
    for e = 1:2
        r = 1 + (n-1)*2 + e;
        a = double(S.StimParams{r,16});
        if a > 0, amp(n,e) = a; ptd(n,e) = double(S.StimParams{r,6}); nsp(n,e) = double(S.StimParams{r,8}); end
        per(n) = double(S.StimParams{r,9});
    end
end
ampMax = max(amp,[],2);

%% ---- condition table ----
ids = unique(cnd)';
fprintf('\nConditions in this session:\n%5s %-26s %9s %13s %9s %6s\n', 'ID', 'label', 'amp E1/E2', 'ptd E1/E2 us', 'nsp E1/E2', 'trials');
for k = ids
    i = find(cnd==k, 1);
    fprintf('%5d %-26s %4g/%-4g %6g/%-6g %4g/%-4g %6d\n', k, lab{i}, amp(i,1), amp(i,2), ptd(i,1), ptd(i,2), nsp(i,1), nsp(i,2), sum(cnd==k));
end

%% ---- build plotting units (a unit = one subplot) ----
sel = ids;
if ~isempty(conds), sel = intersect(sel, conds); end
U = struct('type',{},'amp',{},'key',{},'num',{},'ids',{},'tags',{});
for k = sel
    i = find(cnd==k,1);  L = lab{i};
    if startsWith(L,'sim'),       ty = 'sim'; tag = 'sim';
    elseif contains(L,'A->B'),    ty = 'seq'; tag = 'A->B';
    elseif contains(L,'B->A'),    ty = 'seq'; tag = 'B->A';
    else,                         ty = 'ctrl'; tag = 'ctrl'; end
    if ~isempty(labelHas) && isempty(regexp(L, labelHas, 'once')), continue; end
    tk = regexp(L,'(\d+)/(\d+)','tokens','once');
    if isempty(tk), key = L; num = k; else, key = [tk{1} '/' tk{2}]; num = str2double(tk{1}); end
    a = ampMax(i);
    merged = false;
    if strcmp(ty,'seq') && mergeOrders
        for u = 1:numel(U)
            if strcmp(U(u).type,'seq') && U(u).amp==a && strcmp(U(u).key,key)
                U(u).ids(end+1) = k;  U(u).tags{end+1} = tag;  merged = true;  break
            end
        end
    end
    if ~merged, U(end+1) = struct('type',ty,'amp',a,'key',key,'num',num,'ids',k,'tags',{{tag}}); end %#ok<SAGROW>
end
keep = ismember({U.type}, types);
if ~isempty(ampOnly), keep = keep & ismember([U.amp], ampOnly); end
U = U(keep);
assert(~isempty(U), 'No condition matches your filters.');
tnum = zeros(numel(U),1);
for u = 1:numel(U), tnum(u) = find(strcmp(U(u).type, {'sim','seq','ctrl'})); end
[~, ord] = sortrows([[U.amp]', tnum, [U.num]']);
U = U(ord);  tnum = tnum(ord);
[~,~,gidx] = unique([[U.amp]', tnum], 'rows', 'stable');

if ischar(chans) || isstring(chans), chans = 1:nCh; end
chans = chans(:)';
fprintf('\n%d subplot unit(s) in %d group(s) (amp x type); channels %s; mode %s\n', numel(U), max(gidx), mat2str(chans), plotMode);

%% ---- read the windows once ----
d = dir(fullfile(intanFolder,'amplifier.dat'));
assert(d.bytes == nCh*N*2, 'amplifier.dat size does not match %d ch x %d samples x 2 bytes', nCh, N);
m = memmapfile(fullfile(intanFolder,'amplifier.dat'),'Format',{'int16',[nCh N],'a'});
w = round(winMs/1e3*fs);  idxRel = w(1):w(2);  tms = idxRel/fs*1e3;
X = cell(numel(U),1);  firstTrial = zeros(numel(U),1);
for u = 1:numel(U)
    X{u} = cell(1, numel(U(u).ids));
    for j = 1:numel(U(u).ids)
        blk = zeros(numel(chans), numel(idxRel), 0, 'single');
        for tr = find(cnd == U(u).ids(j))'
            if size(blk,3) >= nTrialsPlot, break; end
            uu = up(tr);  s = find(segStart <= uu, 1, 'last');
            if uu+w(1) < segStart(s) || uu+w(2) > segEnd(s), continue; end
            blk(:,:,end+1) = single(double(m.Data.a(chans, uu+idxRel)) * 0.195); %#ok<SAGROW>
            if j==1 && firstTrial(u)==0, firstTrial(u) = tr; end
        end
        X{u}{j} = blk;
    end
end

%% ---- plot ----
tn = {'sim','seq','ctrl'};
nFig = 0;
if strcmp(plotMode,'by_channel')
    for ci = 1:numel(chans)
        for g = 1:max(gidx)
            us = find(gidx==g);
            if nFig >= maxFigures, warning('maxFigures (%d) reached.', maxFigures); break, end
            nFig = nFig + 1;
            figure('Color','w','Position',figSize);
            nc = min(5, numel(us));  nr = ceil(numel(us)/nc);
            for s = 1:numel(us)
                u = us(s);  subplot(nr, nc, s);
                Ys = cell(1,numel(U(u).ids));
                for j = 1:numel(Ys), Ys{j} = reshape(X{u}{j}(ci,:,:), numel(tms), []); end
                drawPanel_(Ys, U(u).tags, tms, tickTimes_(firstTrial(u), amp, ptd, nsp, per, pulseStartMs, showExpected), ylimUV, ...
                    panelTitle_(U(u), lab, cnd, nsp, firstTrial(u)));
                if s > (nr-1)*nc, xlabel('ms from trigger'); end
                if mod(s-1,nc)==0, ylabel('\muV'); end
            end
            ttl = sprintf('%s | %g \\muA | amplifier row %d', tn{tnum(us(1))}, U(us(1)).amp, chans(ci));
            if any(strcmp(U(us(1)).type,'seq')) && mergeOrders, ttl = [ttl ' | red = A->B, blue = B->A']; end
            sgtitle(ttl);
        end
    end
else   % by_condition: one figure per unit, subplot per channel
    for u = 1:numel(U)
        if nFig >= maxFigures, warning('maxFigures (%d) reached.', maxFigures); break, end
        nFig = nFig + 1;
        figure('Color','w','Position',figSize);
        nc = min(6, numel(chans));  nr = ceil(numel(chans)/nc);
        for ci = 1:numel(chans)
            subplot(nr, nc, ci);
            Ys = cell(1,numel(U(u).ids));
            for j = 1:numel(Ys), Ys{j} = reshape(X{u}{j}(ci,:,:), numel(tms), []); end
            drawPanel_(Ys, U(u).tags, tms, tickTimes_(firstTrial(u), amp, ptd, nsp, per, pulseStartMs, showExpected), ylimUV, ...
                sprintf('row %d', chans(ci)));
            if ci > (nr-1)*nc, xlabel('ms from trigger'); end
            if mod(ci-1,nc)==0, ylabel('\muV'); end
        end
        sgtitle([tn{tnum(u)} ' | ' panelTitle_(U(u), lab, cnd, nsp, firstTrial(u)) sprintf(' | %g \\muA', U(u).amp)]);
    end
end

%% ======================== local functions ==============================
function T = tickTimes_(i, amp, ptd, nsp, per, startMs, on)
T = [];
if ~on || i==0, return, end
for e = 1:2
    if amp(i,e) > 0, T = [T, startMs + ptd(i,e)/1e3 + (0:nsp(i,e)-1)*per(i)/1e3]; end %#ok<AGROW>
end
T = unique(T);
end

function s = panelTitle_(Ui, lab, cnd, nsp, i)
if numel(Ui.ids) > 1
    j2 = find(cnd==Ui.ids(2),1);
    s = sprintf('seq %s  A->B [%g %g] | B->A [%g %g]', Ui.key, nsp(i,1), nsp(i,2), nsp(j2,1), nsp(j2,2));
else
    s = lab{find(cnd==Ui.ids(1),1)};
end
end

function drawPanel_(Ys, tags, tms, tickT, ylimUV, ttl)
hold on; merged = numel(Ys) > 1;
cols = zeros(numel(Ys),3);
for j = 1:numel(Ys)
    Y = Ys{j};  if isempty(Y), continue, end
    if merged
        if strcmp(tags{j},'A->B'), c = [0.85 0.2 0.2]; else, c = [0.15 0.35 0.9]; end
        plot(tms, Y, 'Color', [c 0.25]);
    else
        c = [0 0 0];  plot(tms, Y, 'Color', [0.72 0.72 0.72]);
    end
    cols(j,:) = c;
end
for j = 1:numel(Ys)
    if ~isempty(Ys{j}), plot(tms, mean(Ys{j},2), 'Color', cols(j,:)*0.8, 'LineWidth', 1.6); end
end
if ~isempty(ylimUV), ylim(ylimUV); end
yl = ylim;  rg = diff(yl);  ylim([yl(1) yl(2)+0.1*rg]);
plot([0 0], [yl(1) yl(2)], ':', 'Color', [0.4 0.4 0.4]);
if ~isempty(tickT), plot(tickT, (yl(2)+0.05*rg)*ones(size(tickT)), 'v', 'MarkerSize', 6, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none'); end
xlim(tms([1 end]));  title(ttl, 'Interpreter', 'none', 'FontSize', 10);
end