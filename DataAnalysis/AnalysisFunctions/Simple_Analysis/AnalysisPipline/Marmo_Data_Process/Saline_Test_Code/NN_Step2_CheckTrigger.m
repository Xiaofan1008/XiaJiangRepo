%% CHECK 2: trigger edges in digitalin.dat vs. time.dat segments (script)
% Read-only. Finds which digital input carries the trigger, counts the edges
% per recording segment, and shows where the trigger sits inside each segment
% (so we know whether the stimulus + response are inside the recording or hit
% the restart gap). Paste me the whole printed output.

clear; clc;

%% ===== SETTINGS: EDIT HERE =====
intanFolder = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_27/Train_Amp_Test2_261027_192530';
fs      = 30000;
nTrials = 275;     % trials that ran (Train_Amp_test2 = 275)
trigBit = [];      % [] = pick automatically; or a number 0-15 (digital input index)
needMs  = 600;     % recording we want after the trigger (trial 500/550 ms + margin)
%% ================================

%% ---- load ----
fid = fopen(fullfile(intanFolder,'time.dat'),'r');    t = fread(fid,Inf,'int32');  fclose(fid);
fid = fopen(fullfile(intanFolder,'digitalin.dat'),'r'); D = fread(fid,Inf,'uint16'); fclose(fid);
assert(numel(D)==numel(t), 'digitalin (%d) and time (%d) differ in length', numel(D), numel(t));
brk = find(diff(t) ~= 1);
segStart = [1; brk+1];  segEnd = [brk; numel(t)];  nSeg = numel(segStart);
segOf = zeros(numel(t),1);
for s = 1:nSeg, segOf(segStart(s):segEnd(s)) = s; end

%% ---- which bits are active ----
fprintf('\nDigital inputs with activity (bit = 0..15):\n');
cnt = zeros(1,16);  fall = zeros(1,16);
for b = 0:15
    x = bitget(D,b+1);
    if any(x)
        cnt(b+1)  = sum(diff(x)==1);
        fall(b+1) = sum(diff(x)==-1);
        fprintf('  bit %2d: %6d rising, %6d falling, high %.2f%% of samples\n', ...
            b, cnt(b+1), fall(b+1), 100*mean(x));
    end
end
if isempty(trigBit)
    [~,i] = min(abs(cnt - nTrials) + 1e9*(cnt==0));
    trigBit = i-1;
    fprintf('Auto-picked trigBit = %d (rising edges closest to %d trials)\n', trigBit, nTrials);
end
x = bitget(D,trigBit+1);
up = find(diff(x)==1)+1;                   % rising edge sample index
dn = find(diff(x)==-1)+1;
fprintf('\nTrigger bit %d: %d rising edges (trials ran: %d)\n', trigBit, numel(up), nTrials);
if ~isempty(dn)
    w = []; for k = 1:numel(up), j = find(dn>up(k),1); if ~isempty(j), w(end+1) = dn(j)-up(k); end, end %#ok<AGROW>
    fprintf('Pulse width (ms): min %.2f median %.2f max %.2f\n', min(w)/fs*1e3, median(w)/fs*1e3, max(w)/fs*1e3);
end

%% ---- edges per segment ----
nEdge = accumarray(segOf(up), 1, [nSeg 1]);
fprintf('\nEdges per segment: 0 edges in %d segments | 1 edge in %d | 2+ edges in %d\n', ...
    sum(nEdge==0), sum(nEdge==1), sum(nEdge>=2));
fprintf('Segments with 0 edges: %s\n', mat2str(find(nEdge==0)'));
fprintf('Segments with 2+ edges: %s\n', mat2str(find(nEdge>=2)'));

% first edge in each segment: time from segment start, and recording left after it
first = nan(nSeg,1); left = nan(nSeg,1);
for s = 1:nSeg
    k = up(segOf(up)==s);
    if ~isempty(k)
        first(s) = (k(1)-segStart(s))/fs*1e3;
        left(s)  = (segEnd(s)-k(1))/fs*1e3;
    end
end
segMs = (segEnd-segStart+1)/fs*1e3;
mid = 2:nSeg-1;     % trial-sized segments (skip the long first and last)
fprintf('\nTrial-sized segments (2..%d):\n', nSeg-1);
fprintf('  trigger time from segment start (ms): min %.0f | median %.0f | max %.0f\n', ...
    min(first(mid)), median(first(mid),'omitnan'), max(first(mid)));
fprintf('  recording left AFTER trigger (ms):    min %.0f | median %.0f | max %.0f\n', ...
    min(left(mid)), median(left(mid),'omitnan'), max(left(mid)));
fprintf('  segments with < %d ms after trigger: %d of %d\n', needMs, sum(left(mid)<needMs), numel(mid));
fprintf('  first (long) segment: edges %d, last (long) segment: edges %d\n', nEdge(1), nEdge(end));

fprintf('\nFirst 15 segments: [seg, length ms, edges, trigger at ms, ms left after trigger]\n');
for s = 1:min(15,nSeg)
    fprintf('  %3d  len %8.1f  edges %d  trig %8.1f  left %8.1f\n', s, segMs(s), nEdge(s), first(s), left(s));
end
fprintf('\nAny other bit with activity inside the segments is listed at the top; tell me what each one is wired to.\n');