function run_sab_preprocess(folderPath)
% RUN_SAB_PREPROCESS  Run artifact blanking, filtering, and spike
% detection (Sections 1-2 of sab_quickanalysis.m) on a given data folder,
% regardless of MATLAB's current working directory.
%
% Usage:
%   run_sab_preprocess('D:\data\marmoset1\session3')
%   run_sab_preprocess(pwd)   % same as running the old script in place
%
% This wrapper cd's into folderPath, runs the same processing that used
% to require you to manually cd there first, and always restores your
% original working directory afterwards -- even if an error occurs.

if nargin < 1 || isempty(folderPath)
    error('run_sab_preprocess:noPath', ...
        'You must provide the data folder path, e.g. run_sab_preprocess(''D:\\data\\session3'')');
end

if ~isfolder(folderPath)
    error('run_sab_preprocess:badPath', 'Folder does not exist: %s', folderPath);
end

% Remember where we started, and guarantee we return here no matter what
% (error, Ctrl-C, or normal completion).
origDir = pwd;
cleanupObj = onCleanup(@() cd(origDir)); %#ok<NASGU>

cd(folderPath);
fprintf('Running preprocessing in: %s\n', folderPath);

%% ===================== Parameters to alter =====================
Startpoint_analyse = 0;            % set to 0 for no input
Overall_time_to_analyse = 0;       % time from beginning of Startpoint_analyse; set to 0 for no input
artefact = -500;                   % removes spikes below this threshold
artefact_high = 500;               % removes spikes above this threshold
startpointseconds = 2;             % how long after the trigger to skip spike analysis (ms)
secondstoanalyse = 8;              % how long after the trigger to analyse spikes for (ms)
printspiking = 0;
par = 0;

%% ===================== 1. Blank stimulus =====================
FS = 30000;
filepath = pwd;   % now guaranteed to equal folderPath

fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn = 32;
    E_Mapnumber = 0;
else
    E_Mapnumber = loadMapNum;
    if E_Mapnumber > 0
        nChn = 64;
    else
        nChn = 32;
    end
end

dName = 'amplifier';
vFID = fopen([filepath filesep dName '.dat'], 'r');
mem_check = dir([filepath filesep 'amplifier.dat']);
T = mem_check.bytes ./ (2 * nChn * 30000);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes / (nChn * 2 * 30000);

if t_len < (Overall_time_to_analyse + Startpoint_analyse)
    error('Time to analyse exceeds the data recording period')
end
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
elseif Overall_time_to_analyse ~= 0
    T = Overall_time_to_analyse;
end

info = fileinfo.bytes / 2;
nL = (ceil(info / (nChn * FS * double(T))) + 1); %#ok<NASGU>
vblank = []; %#ok<NASGU>
BREAK = 1;   %#ok<NASGU>
N = 1;       %#ok<NASGU>

denoiseIntan_sab(filepath, dName, T, par, Startpoint_analyse, Overall_time_to_analyse);
trig = loadTrig(0);
theseTrig = trig; %#ok<NASGU>

if ~isempty(vFID) && vFID ~= -1
    fclose(vFID);
end

%% ===================== 2. Thresholds & Mu =====================
allExtract_sab_1(dName, filepath, T, par, artefact, artefact_high);
% alternate: allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');

fprintf('Done. Working directory restored to: %s\n', origDir);

end
