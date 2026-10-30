%% CHECK 9: raw vs denoised overlay, does the blanking work? (script)
% One figure per channel. One column per condition you list.
%   top row   : raw (grey) and denoised (blue), full scale
%   bottom row: denoised only, zoomed, so leftover artifact / tails are visible
% Black triangles = pulse times StimParams says should occur (pulses start at the trigger).
% Prints, per condition, how much the denoised data differs from raw BEFORE the
% first pulse (should be ~0 if only the artifact is touched).
% Nothing is saved.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
intanFolder = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_27/Train_Amp_Test2';
fs      = 30000;
nCh     = 32;
trigBit = 0;

chans   = [32 9];                 % amplifier rows (one figure each)
condSubstr = {'sim 3/5 @50', 'sim 3/5 @75', 'A->B 10/10 [5 4] @50', 'A->B 10/10 [5 4] @75'};   % label substrings (see check4's table)
nTrials = 8;                      % trials overlaid per condition
winMs   = [-10 100];              % ms around the trigger
zoomUV  = 300;                    % bottom row y-limit (+-)
%% ================================

%% ---- load ----
fid = fopen(fullfile(intanFolder,'time.dat'),'r');      t = fread(fid,Inf,'int32');  fclose(fid);
fid = fopen(fullfile(intanFolder,'digitalin.dat'),'r'); D = fread(fid,Inf,'uint16'); fclose(fid);
N = numel(t);
brk = find(diff(t) ~= 1);  segStart = [1; brk+1];  segEnd = [brk; N];
up = find(diff(bitget(D,trigBit+1))==1) + 1;  clear D
f = dir(fullfile(intanFolder,'*_exp_datafile_*.mat'));
S = load(fullfile(intanFolder,f(1).name),'StimParams','TrialTable');
nT = min(numel(up), numel(S.TrialTable));
lab = cellstr({S.TrialTable(1:nT).label}');
amp = zeros(nT,2); ptd = zeros(nT,2); nsp = zeros(nT,2); per = zeros(nT,1);
for n = 1:nT
    for e = 1:2
        r = 1 + (n-1)*2 + e;  a = double(S.StimParams{r,16});
        if a > 0, amp(n,e) = a; ptd(n,e) = double(S.StimParams{r,6}); nsp(n,e) = double(S.StimParams{r,8}); end
        per(n) = double(S.StimParams{r,9});
    end
end

rawF = fullfile(intanFolder,'amplifier.dat');  dnF = fullfile(intanFolder,'amplifier_dn_sab.dat');
assert(exist(dnF,'file')==2, 'amplifier_dn_sab.dat not found in %s', intanFolder);
mR = memmapfile(rawF,'Format',{'int16',[nCh N],'a'});
mD = memmapfile(dnF, 'Format',{'int16',[nCh N],'a'});
w = round(winMs/1e3*fs);  idxRel = w(1):w(2);  tms = idxRel/fs*1e3;

%% ---- trials per condition ----
nC = numel(condSubstr);
tr = cell(1,nC);
for c = 1:nC
    cand = find(contains(lab, condSubstr{c}));
    keep = [];
    for i = cand'
        u = up(i);  s = find(segStart <= u, 1, 'last');
        if u+w(1) >= segStart(s) && u+w(2) <= segEnd(s), keep(end+1) = i; end %#ok<SAGROW>
        if numel(keep) >= nTrials, break, end
    end
    tr{c} = keep;
    if isempty(keep), fprintf('No usable trial for "%s"\n', condSubstr{c}); end
end

%% ---- figures ----
fprintf('\nMax |denoised - raw| in the 5 ms BEFORE the first pulse (uV, median over trials):\n');
for ci = 1:numel(chans)
    figure('Color','w','Position',[60+30*ci 60+20*ci 1500 800]);
    for c = 1:nC
        i1 = tr{c};  if isempty(i1), continue, end
        subplot(2, nC, c); hold on;
        subplot(2, nC, nC+c); hold on;
        pre = zeros(1,numel(i1));  R = zeros(numel(tms), numel(i1));  Dn = R;
        for k = 1:numel(i1)
            u = up(i1(k));
            R(:,k)  = double(mR.Data.a(chans(ci), u+idxRel)) * 0.195;
            Dn(:,k) = double(mD.Data.a(chans(ci), u+idxRel)) * 0.195;
            pre(k) = max(abs(Dn(tms >= -5 & tms < -0.1, k) - R(tms >= -5 & tms < -0.1, k)));
        end
        subplot(2, nC, c);
        plot(tms, R, 'Color', [0.75 0.75 0.75]);  plot(tms, Dn, 'Color', [0.2 0.45 0.95 0.5]);
        plot(tms, mean(Dn,2), 'Color', [0.05 0.25 0.8], 'LineWidth', 1.4);
        ylim([-7000 7600]);  xlim(tms([1 end]));
        % expected pulse times
        j = i1(1);  T = [];
        for e = 1:2, if amp(j,e) > 0, T = [T, ptd(j,e)/1e3 + (0:nsp(j,e)-1)*per(j)/1e3]; end, end %#ok<AGROW>
        T = unique(T);  plot(T, 7200*ones(size(T)), 'v', 'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none');
        title(sprintf('%s\nrow %d | raw grey, denoised blue', lab{j}, chans(ci)), 'Interpreter', 'none', 'FontSize', 8);
        subplot(2, nC, nC+c);
        plot(tms, Dn, 'Color', [0.2 0.45 0.95 0.5]);  plot(tms, mean(Dn,2), 'Color', [0.05 0.25 0.8], 'LineWidth', 1.4);
        ylim([-zoomUV zoomUV]);  xlim(tms([1 end]));  xlabel('ms from trigger');
        plot(T, zoomUV*0.95*ones(size(T)), 'v', 'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', 'none');
        title(sprintf('denoised zoomed (+-%g uV)', zoomUV), 'FontSize', 8);  if c == 1, ylabel('\muV'); end
        fprintf('  row %2d | %-26s | %.2f uV (n=%d)\n', chans(ci), lab{j}, median(pre), numel(i1));
    end
end