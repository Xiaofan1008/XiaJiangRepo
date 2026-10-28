%% CHECK 6: visual-only / blank trials must have NO stimulation artifact (script)
% For every trial and every channel: biggest fast step in 0..150 ms after the
% trigger, in units of that channel's own baseline noise (z). Electrical
% trials should be huge, visual / blank trials should look like baseline.
% Prints a summary per control label and plots (nothing saved):
%   Fig 1: per-trial max z over channels (log), controls vs electrical
%   Fig 2: control trials overlaid on rows 9 / 10 / 32 at full scale and zoomed in

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
intanFolder = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_27/Train_Amp_Test2_261027_192530';
fs      = 30000;
nCh     = 32;
trigBit = 0;
preMs   = 600;  postMs = 550;
thrZ    = 15;                % a step bigger than this many baseline noise units counts as "pulse-like"
showRows = [9 10 32];        % channels for Fig 2
nOverlay = 10;               % control trials overlaid in Fig 2
%% ================================

%% ---- load ----
fid = fopen(fullfile(intanFolder,'time.dat'),'r');      t = fread(fid,Inf,'int32');  fclose(fid);
fid = fopen(fullfile(intanFolder,'digitalin.dat'),'r'); D = fread(fid,Inf,'uint16'); fclose(fid);
N = numel(t);
brk = find(diff(t) ~= 1);  segStart = [1; brk+1];  segEnd = [brk; N];
up = find(diff(bitget(D,trigBit+1))==1) + 1;  clear D
f = dir(fullfile(intanFolder,'*_exp_datafile_*.mat'));
assert(~isempty(f), 'No *_exp_datafile_*.mat in the folder.');
S = load(fullfile(intanFolder,f(1).name),'StimParams','TrialTable');
nT = min(numel(up), numel(S.TrialTable));
lab = cellstr({S.TrialTable(1:nT).label}');
ampMax = zeros(nT,1);
for n = 1:nT, ampMax(n) = max([double(S.StimParams{2+(n-1)*2,16}), double(S.StimParams{3+(n-1)*2,16})]); end
isCtrl = ampMax <= 0;
fprintf('\nTrials: %d electrical, %d visual/blank (amp 0)\n', sum(~isCtrl), sum(isCtrl));

d = dir(fullfile(intanFolder,'amplifier.dat'));
assert(d.bytes == nCh*N*2, 'amplifier.dat size does not match %d ch x %d samples x 2 bytes', nCh, N);
m = memmapfile(fullfile(intanFolder,'amplifier.dat'),'Format',{'int16',[nCh N],'a'});
pre = round(preMs/1e3*fs);  post = round(postMs/1e3*fs);

%% ---- artifact size per trial and channel ----
ok = false(nT,1);  Z = nan(nT,nCh);  stepUV = nan(nT,nCh);
for i = 1:nT
    u = up(i);  s = find(segStart <= u, 1, 'last');
    idx = (u-pre):(u+post);
    if idx(1) < segStart(s) || idx(end) > segEnd(s), continue, end
    ok(i) = true;
    dW = abs(diff(double(m.Data.a(:, idx)) * 0.195, 1, 2));
    rel = (idx(2:end)-u)/fs*1e3;
    bm = rel >= -600 & rel <= -150;  am = rel >= 0 & rel <= 150;
    sg = max(1.4826*median(abs(dW(:,bm) - median(dW(:,bm),2)), 2), eps);
    stepUV(i,:) = max(dW(:,am), [], 2)';
    Z(i,:) = (max(dW(:,am), [], 2) ./ sg)';
end
zMax = max(Z, [], 2);
E = ok & ~isCtrl;  C = ok & isCtrl;

%% ---- summary ----
fprintf('\nMax z over channels (a step of z > %g is pulse-like):\n', thrZ);
fprintf('  electrical trials: median %.0f | min %.0f\n', median(zMax(E)), min(zMax(E)));
fprintf('  visual/blank:      median %.1f | max %.1f\n', median(zMax(C)), max(zMax(C)));
fprintf('  visual/blank trials with ANY channel above threshold: %d of %d\n', sum(zMax(C) > thrZ), sum(C));
fprintf('  electrical trials BELOW threshold on every channel:   %d of %d\n', sum(zMax(E) <= thrZ), sum(E));
fprintf('  largest single step on a control trial: %.1f uV (electrical median of the same: %.0f uV)\n', ...
    max(stepUV(C,:),[],'all'), median(max(stepUV(E,:),[],2)));
fprintf('\nPer control label: [n, median max-z, max max-z, trials above threshold]\n');
cl = unique(lab(C));
for k = 1:numel(cl)
    r = C & strcmp(lab, cl{k});
    fprintf('  %-28s n %3d  median %5.1f  max %6.1f  above %d\n', cl{k}, sum(r), median(zMax(r)), max(zMax(r)), sum(zMax(r) > thrZ));
end
badC = find(C & zMax > thrZ);
if ~isempty(badC)
    fprintf('\nControl trials that look like they have a pulse (trial, label, max z, channel):\n');
    for q = badC(1:min(20,end))'
        [~, c] = max(Z(q,:));  fprintf('  %3d  %-26s z %7.1f  row %d\n', q, lab{q}, zMax(q), c);
    end
end

%% ---- Fig 1: per-trial max z ----
figure('Color','w','Position',[80 80 1100 450]); hold on;
semilogy(find(E), zMax(E), 'o', 'MarkerSize', 4, 'DisplayName', 'electrical');
semilogy(find(C), zMax(C), 's', 'MarkerSize', 5, 'MarkerFaceColor', [0.9 0.4 0.1], 'DisplayName', 'visual / blank');
plot([1 nT], [thrZ thrZ], 'k--', 'DisplayName', 'threshold'); set(gca, 'YScale', 'log');
xlabel('trial'); ylabel('max over channels of step / baseline noise (z)'); legend('Location', 'east');
title('Controls should sit near baseline noise, far below the electrical trials');

%% ---- Fig 2: overlay of control trials ----
cIdx = find(C);  cIdx = cIdx(1:min(nOverlay, numel(cIdx)));
eFirst = find(E, 1);
w = round([-20 150]/1e3*fs);  idxRel = w(1):w(2);  tms = idxRel/fs*1e3;
figure('Color','w','Position',[100 60 1500 800]);
for p = 1:numel(showRows)
    for rr = 1:2
        subplot(2, numel(showRows), (rr-1)*numel(showRows)+p); hold on;
        for i = cIdx(:)'
            x = double(m.Data.a(showRows(p), up(i)+idxRel)) * 0.195;
            plot(tms, x - mean(x(tms<0)), 'Color', [0.9 0.4 0.1 0.5]);
        end
        if rr == 1
            x = double(m.Data.a(showRows(p), up(eFirst)+idxRel)) * 0.195;
            plot(tms, x - mean(x(tms<0)), 'Color', [0.2 0.2 0.2]);
            ylim([-6500 6500]);
            title(sprintf('row %d: controls (orange) vs one electrical trial (black), full scale', showRows(p)), 'FontSize', 8);
        else
            ylim([-300 300]);
            title(sprintf('row %d: controls zoomed (+-300 uV)', showRows(p)), 'FontSize', 8);
            xlabel('ms from trigger');
        end
        xlim(tms([1 end])); ylabel('\muV');
    end
end