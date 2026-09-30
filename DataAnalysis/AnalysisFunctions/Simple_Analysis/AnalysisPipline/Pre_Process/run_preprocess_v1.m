%% RUN_SAB_PREPROCESS (script version)
% Runs artifact blanking, threshold estimation, and MU filtering
% (Sections 1-2 of sab_quickanalysis.m) on a given data folder,
% regardless of MATLAB's current working directory. Stops after writing
% <name>.mu_sab.dat -- does NOT run spike detection / save .sp.mat
% (that's handled by a separate script now).
%
% Usage: edit folderPath below, then press Run.

folderPath = '/Volumes/MACData/Data/Data_Xia/DX036/Xia_Linearity_600_700um_Single1_260929_155142';   % <-- EDIT THIS, then Run

%% ===================== Setup =====================
if isempty(folderPath)
    error('run_sab_preprocess:noPath', ...
        'Set folderPath above to your data folder before running.');
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
artefact = -500;                   % kept for call-signature compatibility (unused now spike detection is skipped)
artefact_high = 500;               % kept for call-signature compatibility (unused now spike detection is skipped)
par = 0;

%% ===================== 1. Blank stimulus =====================
FS = 30000;
filepath = pwd;   % now guaranteed to equal folderPath
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn = 32;
    E_Mapnumber = 0; %#ok<NASGU>
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

%% ===================== 2. Thresholds & Mu (stops before spike detection) =====================
allExtract_sab_1_stopAtMu(dName, filepath, T, par, artefact, artefact_high);

fclose('all');
fprintf('Done. Working directory restored to: %s\n', origDir);