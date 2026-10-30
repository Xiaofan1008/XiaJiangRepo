%% BUILD_EXP_DATAFILE  (script)
% Builds the old-style "*_exp_datafile_001.mat" for the NN (Intan)
% pre-processing from the files the stimulus PC saved.
%
% Standalone: no ProbeMAP, no E_MAP script, no working-folder requirement.
% The probe maps are built in. Neurostim is only needed to read the .mat
% (it is added to the path automatically if it is found).
%
% The old pipeline (run_sab_preprocess / denoiseIntan_sab / loadTrialParams /
% loadStimParams) reads StimParams, TrialParams, simultaneous_stim, E_MAP and
% n_Trials from ONE *_exp_datafile_*.mat in the Intan folder.
%
% HOW TO USE: edit the settings below, then run the script.
% Needs MATLAB R2016b or newer. Written without access to MATLAB: compare the
% printed check at the end with the stim-file summary (trials, reps/cond).

clear; clc;

addpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Marmo_Exp');   % folder that holds estimIntanMulti.m

%% ======================== SETTINGS: EDIT HERE ==========================
stimFile   = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_ISI_Test2_261029_145947/Train_ISI_test2.expTrainISI.145939.mat';        % Neurostim file
condFile   = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_ISI_Test2_261029_145947/Train_ISI_test2.expTrainISI.145939.conds.mat';  % conds file
outFolder  = '/Volumes/MACData/Data/Data_Xia/Marmo_Saline_Test/26_10_29/Train_ISI_Test2_261029_145947';      % Intan folder; the file is written here

probe      = 'flex1';   % 'flex1' 1-shank flexible (32 ch) | 'flex4' 4-shank flexible (128 ch)
                        % | 'hybrid' 64 ch 4-shank + 32 ch single shank | 'file' (see probeFile)
probeFile  = '';        % only for probe = 'file': an old *_exp_datafile_*.mat to copy E_MAP from

visualTrials = 'keep';  % 'keep': visual / blank trials stay in the file with amplitude -1 on
                        %         both rows (a trial with NO pulse)
                        % 'drop': remove them. Only if the trigger line does NOT fire on them.
maxTrials  = Inf;       % only use the first N trials that ran
neurostimPath = '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Marmo_Visual_Code/neurostim';     % folder containing +neurostim. '' = find it automatically
trialCsv   = '';        % '' = read Neurostim. Or a *_trials.csv (trial,cond,label,startMs)
                        %      to skip Neurostim completely
%% ======================================================================

%% ---- 0. check the paths ----------------------------------------------
[~, stimName] = fileparts(stimFile);
if endsWith(stimName,'.conds'), error('stimFile must be the Neurostim .mat, not the .conds.mat'); end
assert(exist(stimFile,'file')==2, 'Stim file not found: %s', stimFile);
assert(exist(condFile,'file')==2, 'Conds file not found: %s', condFile);
if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% ---- 1. what the stimulus PC says RAN --------------------------------
if ~isempty(trialCsv)
    Tt = readtable(trialCsv);
    cnd = double(Tt.cond(:));  tStart = double(Tt.startMs(:));
else
    cleanupObj = ensureNeurostim_(neurostimPath); %#ok<NASGU>
    S = load(stimFile,'c');
    c = S.c;
    cnd = get(c.prms.condition,'atTrialTime',Inf);   % one value per trial
    if iscell(cnd), cnd = cell2mat(cnd); end
    cnd = double(cnd(:));
    tStart = nan(size(cnd));                         % trial start (Neurostim clock, ms)
    try
        ts = trialStartTime(c.prms.trial);
        m = min(numel(ts),numel(cnd));  tStart(1:m) = ts(1:m);
    catch
        warning('Could not read trial start times; TrialTable.startMs will be NaN.');
    end
end
nRan = numel(cnd);
bad = find(isnan(cnd), 1);
if ~isempty(bad)
    warning('Trial %d has no condition (NaN). Using only trials 1..%d.', bad, bad-1);
    nRan = bad-1;
end
nRan = min(nRan, maxTrials);
cnd = cnd(1:nRan);  tStart = tStart(1:nRan);

%% ---- 2. condition table and arguments --------------------------------
Q = load(condFile);
cond = Q.cond;  args = Q.args;
nCond = numel(cond);
assert(all(cnd>=1 & cnd<=nCond & cnd==round(cnd)), 'Condition index outside 1..%d', nCond);

if isfield(Q,'chnList') && ~isempty(Q.chnList), chn = Q.chnList; else, chn = args.chn; end
chn = cellstr(chn);  assert(numel(chn)==2, 'Expected 2 stimulation channels.');

if isfield(args,'freq') && ~isempty(args.freq)
    period_us = round(1e6/double(args.freq));      % train period, microseconds
else
    period_us = 10000;                              % single pulse files: unused
end
fpd = argOr_(args,'fpd',200);  spd = argOr_(args,'spd',200);  ipi = argOr_(args,'ipi',100);

isVis = false(nCond,1);
for k = 1:nCond, isVis(k) = all(cond(k).amp==0); end   % visual-only and blank: no pulse
if strcmpi(visualTrials,'drop'), keep = ~isVis(cnd); else, keep = true(size(cnd)); end
trialIdx = find(keep);
nT = numel(trialIdx);

%% ---- 3. E_MAP (probe) ------------------------------------------------
E_MAP = probeMap_(probe, probeFile);
names = E_MAP(2:end);
[tf, eIdx] = ismember(chn, names);
assert(all(tf), 'Channel(s) %s not in the %s map.', strjoin(chn(~tf),', '), E_MAP{1});

%% ---- 4. StimParams (25 columns) and TrialParams ----------------------
hdr = {'CHANNEL','ENABLE STIMULATION','TRIGGER SOURCE','TRIGGER EDGE/LEVEL', ...
    'TRIGGER HIGH/LOW','POST-TRIGGER DELAY','PULSE OR TRAIN','NUMBER OF PULSES IN TRAIN', ...
    'PULSE TRAIN PERIOD','REFRACTORY PERIOD','STIMULUS SHAPE','FIRST PHASE POLARITY', ...
    'FIRST PHASE DURATION','SECOND PHASE DURATION','INTER-PHASE DELAY','FIRST PHASE AMPLITUDE', ...
    'SECOND PHASE AMPLITUDE','ENABLE AMP SETTLE','PRE-STIM AMP SETTLE','POST-STIM AMP SETTLE', ...
    'MAINTAIN AMP SETTLE DURING TRAIN','ENABLE CHARGE RECOVERY','POST STIM CHARGE RECOVERY ON', ...
    'POST STIM CHARGE RECOVERY OFF',[]};
nrow = 2*nT;
StimParams  = cell(1+nrow, 25);
StimParams(1,:) = hdr;
TrialParams = cell(1+nrow, 3);
TrialParams(1,:) = {'Trial Number','Trial ID','CHANNEL'};
TrialTable = struct('trial',{},'cond',{},'label',{},'isVisual',{},'startMs',{});

for n = 1:nT
    t = trialIdx(n);  k = cnd(t);  cc = cond(k);
    for e = 1:2
        r = 1 + (n-1)*2 + e;
        amp = double(cc.amp(e));
        if amp <= 0, amp = -1; end                 % -1 = no pulse on this electrode
        StimParams(r,:) = {chn{e}, 1, 0, 0, 0, double(cc.ptd(e)), 0, double(cc.nsp(e)), ...
            period_us, 1000, 1, 0, fpd, spd, ipi, amp, amp, 1, 200, 160000, 1, 1, 100000, 150000, 500};
        TrialParams(r,:) = {double(n), double(k), double(eIdx(e))};
    end
    TrialTable(n).trial = n;  TrialTable(n).cond = k;  TrialTable(n).label = cc.label;
    TrialTable(n).isVisual = isVis(k);  TrialTable(n).startMs = tStart(t);
end
simultaneous_stim = 2;
n_Trials = nT;

%% ---- 5. save ---------------------------------------------------------
parts = strsplit(stimName,'.');
outFile = fullfile(outFolder,[parts{1} '_exp_datafile_001.mat']);
save(outFile,'StimParams','TrialParams','simultaneous_stim','n_Trials','E_MAP','TrialTable','-v7');

%% ---- 6. printed check ------------------------------------------------
fprintf('\nWrote %s\n', outFile);
fprintf('Trials that ran: %d   written: %d   (visual/blank trials: %d, %s)\n', ...
    nRan, nT, sum(isVis(cnd)), visualTrials);
cnt = accumarray(cnd, 1, [nCond 1])';
fprintf('Repeats per condition (min %d, max %d):\n', min(cnt), max(cnt)); disp(cnt);
fprintf('Immediate repeats of the previous condition: %d\n', sum(diff(cnd)==0));
fprintf('Median start-to-start gap: %.0f ms\n', median(diff(tStart),'omitnan'));
fprintf('Probe: %s, %d entries. %s -> E%d,  %s -> E%d\n', E_MAP{1}, numel(names), ...
    chn{1}, eIdx(1), chn{2}, eIdx(2));
fprintf('Put the file in the Intan folder (only one *_exp_datafile_*.mat there).\n');

%% ======================== local functions ==============================
function v = argOr_(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = double(s.(f)); else, v = d; end
end

function E_MAP = probeMap_(probe, probeFile)
% Probe maps copied from ProbeMAP.m (checked against old experiment files).
lab = @(port,ch) arrayfun(@(x) sprintf('%c-%03d',port,x), ch(:), 'UniformOutput', false);
flex  = [8,7,9,6,10,5,12,3,13,2,14,1,23,24,22,25,21,26,19,28,18,29,17,30,16,31,20,27,15,0,11,4];
shA   = [31,7,0,24,30,6,1,25,29,5,2,26,28,4,3,27,  8,16,23,15,9,17,22,14,10,18,21,13,11,19,20,12];
shB   = [8,11,9,15,10,19,12,23,13,22,14,21,16,20,17,18,  4,7,0,6,28,5,24,3,25,2,26,1,27,31,29,30];
switch lower(probe)
    case 'flex1'
        E_MAP = [{'FLEXIBLE I-Series NN'}; lab('A',flex)];
    case 'flex4'
        E_MAP = [{'4SHANK FLEXIBLE I-Series NN'}; lab('A',shA); lab('B',shB); lab('C',shA); lab('D',shB)];
    case 'hybrid'
        E_MAP = [{'64chn+32chn HYBRID I-Series NN'}; lab('A',shA); lab('B',shB); lab('C',flex)];
    case 'file'
        assert(exist(probeFile,'file')==2, 'probeFile not found.');
        L = load(probeFile,'E_MAP');  E_MAP = L.E_MAP;
    otherwise
        error('Unknown probe "%s" (flex1, flex4, hybrid or file).', probe);
end
end

function cleanupObj = ensureNeurostim_(userPath)
% Make sure +neurostim can be found; undo any path change when finished.
cleanupObj = [];
if ~isempty(userPath)
    cand = {userPath};
elseif ~isempty(which('neurostim.cic'))
    return                                         % already on the path
else
    here = fileparts(mfilename('fullpath'));
    home = char(java.lang.System.getProperty('user.home'));
    cand = {here, fullfile(here,'neurostim'), fullfile(here,'neurostim-master'), ...
        fullfile(home,'Desktop','PhD Study','Code','Marmo_Code','neurostim-master'), ...
        fullfile(home,'Desktop','PhD Study','Code','Marmo_Code','Marmo_Visual_Code','neurostim')};
end
for i = 1:numel(cand)
    if exist(fullfile(cand{i},'+neurostim'),'dir')
        old = path;  addpath(cand{i});
        cleanupObj = onCleanup(@() path(old));
        fprintf('Using Neurostim from %s\n', cand{i});
        return
    end
end
error(['Neurostim was not found. Set neurostimPath to its folder, e.g. ' ...
    '.../neurostim-master, or set trialCsv to the *_trials.csv.']);
end