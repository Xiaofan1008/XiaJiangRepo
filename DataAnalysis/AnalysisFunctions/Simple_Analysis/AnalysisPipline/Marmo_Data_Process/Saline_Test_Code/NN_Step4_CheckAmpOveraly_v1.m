%% CHECK 5: does the artifact scale with the programmed amplitude? (script)
% Takes ONE condition (e.g. sim 5/5, or A->B 3/10) at every amplitude level,
% averages the trials, and
%   Figure 1: overlays the amplitudes on every channel (colour = amplitude),
%   Figure 2: artifact peak-to-peak vs amplitude per channel (log-log). If the
%             current scales, far channels follow the dashed slope-1 line;
%             channels where the amplifier clips flatten out.
% Prints a table: peak-to-peak ratio vs the lowest amplitude (expected = amp ratio).
% Nothing is saved.

clear; clc; close all;

%% ===== SETTINGS: EDIT HERE =====
intanFolder = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_27/Train_Amp_Test2_261027_192530';
fs      = 30000;
nCh     = 32;
trigBit = 0;

condType = 'sim';     % 'sim' or 'seq'
condKey  = '5/5';     % the k/n part of the label: sim 1/5..5/5 ; seq 1/10..10/10
order    = 'A->B';    % only for 'seq': 'A->B' or 'B->A'
chansShow = 'all';    % rows for Figure 1, or e.g. [9 10 12 16 24 32]
nTrials  = 20;        % trials averaged per amplitude
winMs    = [-2 10];   % ms around the trigger (covers the first pulse)
stimRows = [9 10];    % highlighted in Figure 2 (the two stim sites)
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
cnd = double([S.TrialTable(1:nT).cond])';  lab = cellstr({S.TrialTable(1:nT).label}');
ampMax = zeros(nT,1);
for n = 1:nT
    ampMax(n) = max([double(S.StimParams{2+(n-1)*2,16}), double(S.StimParams{3+(n-1)*2,16})]);
end

%% ---- find this condition at every amplitude ----
cand = [];
for k = unique(cnd)'
    L = lab{find(cnd==k,1)};
    if strcmp(condType,'sim'), isType = startsWith(L,'sim'); else, isType = startsWith(L,order); end
    tk = regexp(L,'(\d+/\d+)','tokens','once');
    if isType && ~isempty(tk) && strcmp(tk{1},condKey) && ampMax(find(cnd==k,1)) > 0, cand(end+1) = k; end %#ok<SAGROW>
end
assert(~isempty(cand), 'No condition matches type %s, key %s.', condType, condKey);
a = arrayfun(@(k) ampMax(find(cnd==k,1)), cand);
[a, o] = sort(a);  cand = cand(o);
nA = numel(a);
fprintf('\nCondition %s %s%s at amplitudes: %s uA (cond IDs %s)\n', condType, ...
    ternary_(strcmp(condType,'seq'), [order ' '], ''), condKey, mat2str(a), mat2str(cand));

%% ---- mean trace per channel and amplitude ----
d = dir(fullfile(intanFolder,'amplifier.dat'));
assert(d.bytes == nCh*N*2, 'amplifier.dat size does not match %d ch x %d samples x 2 bytes', nCh, N);
m = memmapfile(fullfile(intanFolder,'amplifier.dat'),'Format',{'int16',[nCh N],'a'});
w = round(winMs/1e3*fs);  idxRel = w(1):w(2);  tms = idxRel/fs*1e3;
M = zeros(nCh, numel(idxRel), nA);  nUsed = zeros(1,nA);
for q = 1:nA
    acc = zeros(nCh, numel(idxRel));
    for tr = find(cnd == cand(q))'
        if nUsed(q) >= nTrials, break, end
        u = up(tr);  s = find(segStart <= u, 1, 'last');
        if u+w(1) < segStart(s) || u+w(2) > segEnd(s), continue, end
        x = double(m.Data.a(:, u+idxRel)) * 0.195;
        acc = acc + (x - mean(x(:, idxRel < 0), 2));      % baseline-subtract each trial
        nUsed(q) = nUsed(q) + 1;
    end
    M(:,:,q) = acc / max(nUsed(q),1);
end
fprintf('Trials averaged per amplitude: %s\n', mat2str(nUsed));
p2p = squeeze(max(M,[],2) - min(M,[],2));       % nCh x nA   (uV)

%% ---- table: ratio vs lowest amplitude ----
fprintf('\nPeak-to-peak ratio vs %g uA (expected: %s)\n', a(1), mat2str(round(a/a(1),2)));
fprintf('%4s %10s', 'row', sprintf('p2p@%guA', a(1)));
for q = 2:nA, fprintf(' %9s', sprintf('x@%g', a(q))); end, fprintf('\n');
for c = 1:nCh
    fprintf('%4d %10.0f', c, p2p(c,1));
    for q = 2:nA, fprintf(' %9.2f', p2p(c,q)/p2p(c,1)); end
    fprintf('\n');
end

%% ---- Figure 1: overlay of amplitudes per channel ----
if ischar(chansShow) || isstring(chansShow), chansShow = 1:nCh; end
cm = parula(nA+1);
figure('Color','w','Position',[60 60 1500 900]);
nc = min(6, numel(chansShow));  nr = ceil(numel(chansShow)/nc);
for ci = 1:numel(chansShow)
    subplot(nr, nc, ci); hold on;
    for q = 1:nA, plot(tms, squeeze(M(chansShow(ci),:,q)), 'Color', cm(q,:), 'LineWidth', 1.2); end
    xlim(tms([1 end]));  title(sprintf('row %d', chansShow(ci)), 'FontSize', 9);
    if ci == 1, legend(arrayfun(@(x) sprintf('%g uA',x), a, 'UniformOutput', false), 'Location', 'best', 'FontSize', 7); end
    if ci > (nr-1)*nc, xlabel('ms from trigger'); end
end
sgtitle(sprintf('%s %s: mean trace per amplitude (baseline-subtracted, uV)', condType, condKey));

%% ---- Figure 2: peak-to-peak vs amplitude ----
figure('Color','w','Position',[100 100 800 600]); hold on;
for c = 1:nCh, plot(a, p2p(c,:), '-o', 'Color', [0.75 0.75 0.75], 'MarkerSize', 3); end
for c = stimRows(stimRows <= nCh), plot(a, p2p(c,:), '-o', 'LineWidth', 1.8, 'DisplayName', sprintf('row %d', c)); end
ref = median(p2p(:,1)) * a/a(1);
plot(a, ref, 'k--', 'LineWidth', 1.2, 'DisplayName', 'slope 1 (proportional)');
set(gca, 'XScale', 'log', 'YScale', 'log'); xlabel('programmed amplitude (\muA)'); ylabel('peak-to-peak (\muV)');
legend('Location', 'northwest'); title('Grey = all channels. Proportional scaling follows the dashed line.');

function o = ternary_(c,a,b)
if c, o = a; else, o = b; end
end