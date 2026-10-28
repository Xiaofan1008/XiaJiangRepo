%% CHECK 1: Intan folder inventory + time.dat segments (script)
% Read-only. Lists the files in the Intan folder and, if time.dat exists,
% how many separate recordings (segments) it contains and how long each is.
% Run it, then paste me the whole printed output.

clear; clc;

%% ===== SETTINGS: EDIT HERE =====
intanFolder = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_27/Train_Amp_Test2_261027_192530';   % folder with the Intan files
fs          = 30000;     % Intan sample rate (Hz); info.rhs says if unsure
%% ================================

assert(exist(intanFolder,'dir')==7, 'Folder not found: %s', intanFolder);

%% ---- 1. everything in the folder ----
L = dir(intanFolder);
L = L(~[L.isdir]);
fprintf('\nFolder: %s\n%d files\n\n', intanFolder, numel(L));
[~,ord] = sort({L.name});
for i = ord
    fprintf('  %-55s %12.2f MB\n', L(i).name, L(i).bytes/1e6);
end
exts = cellfun(@(n) lower(n(find(n=='.',1,'last'):end)), {L.name}, 'UniformOutput', false);
[u,~,ix] = unique(exts);
fprintf('\nBy extension:\n');
for k = 1:numel(u), fprintf('  %-8s x %d\n', u{k}, sum(ix==k)); end
fprintf('\nexp_datafile present: %d\n', any(contains({L.name},'_exp_datafile_')));
fprintf('info.rhs present: %d | time.dat: %d | amplifier.dat: %d | board-DIN: %d | stim.dat: %d\n', ...
    exist(fullfile(intanFolder,'info.rhs'),'file')==2, ...
    exist(fullfile(intanFolder,'time.dat'),'file')==2, ...
    exist(fullfile(intanFolder,'amplifier.dat'),'file')==2, ...
    any(contains({L.name},'board-DIN')), exist(fullfile(intanFolder,'stim.dat'),'file')==2);
fprintf('per-channel amp files (amp-*.dat): %d\n', sum(startsWith({L.name},'amp-')));

%% ---- 2. time.dat segments ----
tf = fullfile(intanFolder,'time.dat');
if exist(tf,'file')~=2
    fprintf('\nNo time.dat here: tell me what the recording files are called.\n');
    return
end
fid = fopen(tf,'r'); t = fread(fid,Inf,'int32'); fclose(fid);
fprintf('\ntime.dat: %d samples = %.1f s total\n', numel(t), numel(t)/fs);
d = diff(t);
brk = find(d ~= 1);                    % any reset or jump
fprintf('Breaks (diff ~= 1): %d\n', numel(brk));
segStart = [1; brk+1];  segEnd = [brk; numel(t)];
segLen   = segEnd - segStart + 1;
fprintf('Segments: %d\n', numel(segLen));
fprintf('Segment length (ms): min %.0f | median %.0f | max %.0f\n', ...
    min(segLen)/fs*1e3, median(segLen)/fs*1e3, max(segLen)/fs*1e3);
fprintf('First sample value of each segment: min %d, max %d (0 = fresh restart)\n', ...
    min(t(segStart)), max(t(segStart)));
nShow = min(15, numel(segLen));
fprintf('\nFirst %d segments: [start sample, length ms, first t value, step before it]\n', nShow);
for s = 1:nShow
    if s==1, pre = NaN; else, pre = d(brk(s-1)); end
    fprintf('  %3d  start %9d  len %8.1f ms  t0 %9d  jump %d\n', ...
        s, segStart(s), segLen(s)/fs*1e3, t(segStart(s)), pre);
end
fprintf('\nLast 5 segments (ms): %s\n', mat2str(round(segLen(max(1,end-4):end)'/fs*1e3)));