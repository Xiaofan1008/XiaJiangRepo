% ---------- Strobe Analysis Code ----------

% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-500; %removes spikes below this threshold
artefact_high=500; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=8; %How long after the trigger do you want to analyse spikes for(ms)? 
printspiking=0;
par=0;

FS=30000;
filepath = pwd;
Electrode_Type = 2;
if Electrode_Type == 2
    nChn=64;
else 
    nChn=32;
end
dName='amplifier';
vFID = fopen([filepath filesep dName '.dat'],'r');
mem_check=dir('amplifier.dat');
T = mem_check.bytes ./ (2 * nChn * 30000);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * 30000);
if t_len < (Overall_time_to_analyse+Startpoint_analyse)
    error('Time to analyse exceeds the data recording period')
end
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
elseif Overall_time_to_analyse~=0
    T=Overall_time_to_analyse;
end
info = fileinfo.bytes/2;
nL = (ceil(info / (nChn*FS*double(T)))+1);
vblank=[];
BREAK = 1;
N=1;

%% Load Trigger
if isempty(dir('*.trig.dat'))
cleanTrig_Strobequick;
end
trig = loadTrig(0);
theseTrig = trig;

%% 2. Thresholds & Mu
allExtract_sab_1(dName,filepath,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');