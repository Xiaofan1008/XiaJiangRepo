%% NS INSPECT STIM FILE: what is inside a neurostim .mat log? (script, read-only)
% Loads one neurostim stimulus log and prints a plain report:
%   1. the variables in the file
%   2. the experiment ('cic') basics: subject, paradigm, date, screen, number of trials, durations
%   3. every plugin / stimulus (name + class)
%   4. for each one: every logged parameter with the number of log entries and a preview of its value
%   5. detail on the noise grid stimulus (settings, RNG state, logged clut values, frame/trial timing)
% Nothing is modified or saved except an optional text copy of the report (saveTxt).
% The file contains neurostim objects, so the neurostim folder must be on the MATLAB path
% (and should be the version that wrote the file).

clear; clc;

%% ===== SETTINGS: EDIT HERE =====
matFile  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/NoiseGrid_Files/26_09_15/CJ283.noisegrid.143121.mat';   % stim log to inspect
nsPath   = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Marmo_Visual_Code/neurostim';          % neurostim folder ('' = already on the path)
mstPath  = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Marmo_Visual_Code/MATLAB/marmolab/marmolab-stimuli';    % marmolab-stimuli folder ('' = skip)
maxPrev  = 6;       % how many numbers of a long parameter to show
saveTxt  = '';      % e.g. '/Users/xiaofan/Desktop/stimfile_report.txt' to keep a copy of the report; '' = screen only
%% ================================

if ~isempty(saveTxt), if exist(saveTxt,'file'), delete(saveTxt); end, diary(saveTxt); end
if ~isempty(nsPath),  addpath(genpath(nsPath));  end
if ~isempty(mstPath), addpath(genpath(mstPath)); end

fprintf('\n===== FILE =====\n%s\n', matFile);
d = dir(matFile);  assert(~isempty(d), 'File not found: %s', matFile);
fprintf('size %.2f MB, saved %s\n', d.bytes/1e6, d.date);

%% 1. variables
fprintf('\n===== 1. VARIABLES IN THE FILE =====\n');
w = whos('-file', matFile);
for i = 1:numel(w), fprintf('  %-14s %-18s %s\n', w(i).name, w(i).class, mat2str(w(i).size)); end

S = load(matFile);
if isfield(S,'c'), c = S.c;
else
    k = find(structfun(@(x) isobject(x) && isprop(x,'stimuli'), S), 1);
    fn = fieldnames(S);  assert(~isempty(k), 'No neurostim cic object found in this file');
    c = S.(fn{k});
end

%% 2. experiment basics
fprintf('\n===== 2. EXPERIMENT (cic) =====\n');
fprintf('class: %s\n', class(c));
tryShow_(c, {'subject','paradigm','file','fullFile','date','runNr','expScript','nrTrialsTotal','nrConditions','trial','blockCount'});
try, fprintf('screen: %s\n', evalc('disp(c.screen)')); catch, end
try, fprintf('frame rate (screen.frameRate) = %g Hz\n', c.screen.frameRate); catch, end
try, fprintf('screen pixels = %g x %g, physical = %g x %g, viewing distance = %g\n', c.screen.xpixels, c.screen.ypixels, c.screen.width, c.screen.height, c.screen.viewDist); catch, end
try, fprintf('trials in the log: %d\n', c.nrTrialsTotal); catch, end

%% 3. plugins / stimuli
fprintf('\n===== 3. PLUGINS AND STIMULI =====\n');
names = {};
for src = {'pluginOrder','plugins','stimuli'}
    if ~isempty(names), break, end
    try, names = namesOf_(c.(src{1})); catch, end
end
names = names(~strcmp(names, 'cic'));
fprintf('(%d plugins found)\n', numel(names));
for i = 1:numel(names)
    try, p = c.(names{i}); fprintf('  %-18s %s\n', names{i}, class(p)); catch, fprintf('  %-18s (cannot access)\n', names{i}); end
end

%% 4. parameters of each plugin
fprintf('\n===== 4. LOGGED PARAMETERS (name : number of log entries : value preview) =====\n');
for i = 1:numel(names)
    try, p = c.(names{i}); catch, continue, end
    fprintf('\n--- %s (%s) ---\n', names{i}, class(p));
    prms = [];
    try, prms = p.prms; catch, end
    if isempty(prms) || ~isstruct(prms), fprintf('  (no parameter list)\n'); continue, end
    f = fieldnames(prms);
    for j = 1:numel(f)
        pr = prms.(f{j});
        nLog = NaN;  prev = '';
        try, nLog = numel(pr.log); catch, end
        try
            v = pr.log;  if iscell(v) && ~isempty(v), v = v{1}; end
            prev = preview_(v, maxPrev);
        catch
            try, prev = preview_(p.(f{j}), maxPrev); catch, end
        end
        fprintf('  %-22s : %5g : %s\n', f{j}, nLog, prev);
    end
end

%% 5. noise grid detail
fprintf('\n===== 5. NOISE GRID DETAIL =====\n');
ng = '';
for i = 1:numel(names)
    try
        if contains(lower(class(c.(names{i}))), 'noisegrid') || contains(lower(class(c.(names{i}))), 'noiseclut'), ng = names{i}; break, end
    catch, end
end
if isempty(ng)
    fprintf('No noisegrid / noiseclut stimulus found by class name; look at section 3.\n');
else
    p = c.(ng);
    fprintf('stimulus object: %s (%s)\n', ng, class(p));
    props = {'width','height','size','sparsity','lifetime','sampleFun','rngState','callbackCounter','clutVals','nrRandels','nrChannels','X','Y','position','alphaMask'};
    for j = 1:numel(props)
        try, v = p.(props{j}); fprintf('  %-16s %s\n', props{j}, preview_(v, maxPrev)); catch, end
    end
    % logged series: first, last, how many
    for q = {'rngState','clutVals','callbackCounter','startTime','stopTime','state','startFrame'}
        try
            pr = p.prms.(q{1});
            [vals, tr, tm] = deal(pr.log, pr.trial, pr.time);
            fprintf('  log %-14s entries %5d  | first trial %g, last trial %g | first time %.3f s, last time %.3f s\n', ...
                q{1}, numel(vals), tr(1), tr(end), tm(1)/1000, tm(end)/1000);
        catch, end
    end
end

%% 7. update timing from the callbackCounter log
fprintf('\n===== 7. NOISE UPDATE TIMING (callbackCounter log) =====\n');
if ~isempty(ng)
    try
        pr = c.(ng).prms.callbackCounter;
        [v, tr, trt] = get(pr);                 % value, trial number, time within the trial (ms)
        if iscell(v), v = cellfun(@(x) double(x(1)), v(:)); end
        v = double(v(:));  tr = tr(:);  tm = trt(:);
        fprintf('  get() returned %d entries\n', numel(v));
        trs = unique(tr(tr > 0 & isfinite(tr)))';
        span = nan(size(trs));  cnt = nan(size(trs));
        for k = 1:numel(trs)
            ix = find(tr == trs(k));
            span(k) = (max(tm(ix)) - min(tm(ix)))/1000;  cnt(k) = max(v(ix));
            if k <= 3 || k == numel(trs)
                fprintf('  trial %2d: %d entries, values %s, trial time (s) %s\n', trs(k), numel(ix), mat2str(v(ix)', 6), mat2str(tm(ix)'/1000, 5));
            end
        end
        fprintf('  per trial: logged span median %.3f s (min %.3f, max %.3f); final counter median %g (min %g, max %g)\n', median(span), min(span), max(span), median(cnt), min(cnt), max(cnt));
        fprintf('  implied update interval = span / counter: median %.1f ms\n', 1000*median(span./cnt));
        fprintf('  short trials (span < 0.9*median): %s\n', mat2str(trs(span < 0.9*median(span))));
    catch err
        fprintf('  could not read: %s\n', err.message);
    end
end

%% trial timing (any plugin that logs 'firstFrame' at trial start is cic itself)
fprintf('\n===== 6. TRIAL TIMING FROM cic =====\n');
for q = {'trialStartTime','trialStopTime','firstFrame','frame'}
    try
        pr = c.prms.(q{1});
        t = pr.log;  tm = pr.time;
        fprintf('  %-16s entries %6d | first %s | last %s\n', q{1}, numel(t), preview_(t(1), 2), preview_(t(end), 2));
        if strcmp(q{1}, 'trialStartTime') || strcmp(q{1}, 'firstFrame')
            tt = pr.time(:)/1000;
            fprintf('    trial-to-trial interval: median %.3f s, min %.3f s, max %.3f s\n', median(diff(tt)), min(diff(tt)), max(diff(tt)));
        end
    catch, end
end

fprintf('\nDone. Copy this whole printout back to Claude.\n');
if ~isempty(saveTxt), diary off; fprintf('Report copy saved to %s\n', saveTxt); end

%% ======================== local functions ==============================
function tryShow_(obj, fields)
for i = 1:numel(fields)
    try, v = obj.(fields{i}); fprintf('  %-16s %s\n', fields{i}, preview_(v, 6)); catch, end
end
end

function nm = namesOf_(x)
% plugin names from a cell of plugins, an object array, or a cell of char
nm = {};
if ischar(x), nm = cellstr(x); return, end
if iscell(x)
    for k = 1:numel(x)
        if ischar(x{k}), nm{end+1} = x{k}; else, nm{end+1} = x{k}.name; end %#ok<AGROW>
    end
else
    for k = 1:numel(x), nm{end+1} = x(k).name; end %#ok<AGROW>
end
end

function s = preview_(v, n)
try
    if isnumeric(v) || islogical(v)
        if numel(v) <= n, s = mat2str(double(v), 5);
        else, s = sprintf('%s [%s %s, first %d: %s ...]', class(v), mat2str(size(v)), 'array', n, mat2str(double(v(1:n)), 5)); end
    elseif ischar(v)
        s = ['''' v(1:min(end,80)) ''''];
    elseif isstring(v)
        s = ['"' char(v(1)) '"'];
    elseif iscell(v)
        s = sprintf('cell %s', mat2str(size(v)));
        if ~isempty(v), s = [s ' first: ' preview_(v{1}, n)]; end
    elseif isstruct(v)
        s = sprintf('struct with fields: %s', strjoin(fieldnames(v)', ', '));
    elseif isa(v, 'function_handle')
        s = func2str(v);
    elseif isobject(v)
        s = sprintf('object of class %s', class(v));
    else
        s = class(v);
    end
catch
    s = '(cannot preview)';
end
end