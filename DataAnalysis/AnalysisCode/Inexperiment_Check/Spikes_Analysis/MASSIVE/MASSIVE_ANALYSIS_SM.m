%%for analysing data on massive and saving denoised files

% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-150; %removes spikes below this threshold
artefact_high=100; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=12; %How long after the trigger do you want to analyse spikes for(ms)? 
printspiking=0;
par=0;

%% 1. Blank stimulus

FS=30000;
folder = fileparts(which('MASSIVE_ANALYSIS_SM')); % Determines filepath to folder containing your .m file.
addpath(genpath(folder)); % Add that folder plus all subfolders to the path.
fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd([fileptochange filesep 'Data' filesep 'dual_vary_noshift_pen4_002_200707_205943'])%change working directory to where data is stored - currently manually input
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=loadMapNum;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
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
denoiseIntan_sab(filepath, dName, T, par, Startpoint_analyse, Overall_time_to_analyse);
trig = loadTrig(0);
theseTrig = trig;

%% 2. Thresholds & Mu
allExtract_sab_1(dName,filepath,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');

%% 4. Calculate Structure of sorted trials according to IDs
starttrial=1;
trialjump=1;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
[IDstruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
save('IDstruct.mat', 'IDstruct')

%% 5. Calculates template of trials and spiking responses (Output in true electrode order)
[avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct);
save('Averagetrialresponse.mat','avgnospT','stderrspktrial')

fprintf('endanalysis')