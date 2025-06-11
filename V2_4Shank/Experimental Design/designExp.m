%% Script for generating experimental protocols
function designExp(DIR,AMP,DUR,IPD,E_DIAM,E_MAT,n_REP,CHN,PARAM,PULSE)
% Arguments:
% AMP: A range of stimulation amplitudes. Expects e.g. [1:2:19], in units
% of uA.
% DUR: A range of stimulation durations. Expects e.g. [100:100:400], in
% units of us.
% IPD: Inter-phase delay. Expects e.g. 100, in units of us.
% E_DIAM: Electrode area diameter. Expects e.g. 25, in units of um.
% E_MAT: Electrode material. Expects e.g.'IrOx', in string format. List of
% currently available materials: 'IrOx'.
% E_IMP: Electrode impedance. Expects e.g. 'E_IMP.mat', as a file handle to
% a struct containing impedance values.
% REPEAT: The number of repeats per trial condition. Expects e.g. 75.
% CHN: A cell array range of electrodes to be used for stimulation. Expects e.g.
% [1,4,7]. Number corresponds with designation in E_MAP, e.g. CHN 1 =
% 'A-019' in Intan. Note that this relies on the MAP
% file - be very careful to select the right electrodes.
% ID: A unique identifier for this stimulation set. Expects e.g. '001'
% PARAM: How long to wait after a parameter update. Expects e.g. 500
% (msec)
% PULSE: How long to wait after a stimulation pulse. Expects e.g. 1000
% (msec)
HAPPY = 0;
while (~HAPPY)
    if nargin == 0                      % nargin refers to the num of input variables. 
        DIR = uigetdir('Z:\Tim');
        CASE = 0;
    elseif nargin == 1
        CASE = 0;
    else
        CASE = 1;
    end
    if ~(CASE)
        E_DIAM = 25;
        E_MAT = 'IrOx';
        AMP = input('Please enter stimulus amplitude in uA like so: [x:y:z]\n');
        DUR = input('Please enter stimulus duration in us like so: [x:y:z]\n');
        MAX_CHARGE = AMP(end) * DUR(end) * 1e-12;
        if (strcmp(E_MAT,'IrOx'))
            SAFE_CHARGE = IrOx_safe(E_DIAM);
        end
        if (MAX_CHARGE > SAFE_CHARGE)
            disp('The amount of charge indicated exceeds the safe charge threshold.');
            disp('Please try again');
            continue;
        end
        IPD = input('Please enter interphase delay like so: x\n');
        n_REP = input('Please enter number of repeats like so: x\n');
        isTRAIN = input('Please enter a "1" for stimulus train or "0" for single pulse\n');
        if isTRAIN
            nTRAIN = input('Please enter the maximum number of stimulus pulses\n');
            FREQ = input('Please enter the frequencies of the train like so: [x,y,z]\n');
        else
            nTRAIN = 1;
            FREQ = 100;
        end
        CHN = input('Please enter channels like so: [x,y,z]\n');
        % How long do we need to recover from a parameter update: about 200
        % msec seems normal
        checkPARAM = input('Does this experiment require LFP recovery before stimulation?\nType "1" for YES\n');
        if checkPARAM
            PARAM = 520;
        else
            PARAM = 150;
        end
        % How long do we need to recover from a stimulation pulse: about
        % 150 msec seems normal
        %PULSE = input('Please enter delay after stimulation in ms like so: x\n');
        PULSE = 150;
        settings = newSettings; % This function generates a stimulation settings list of default values
        if max([settings{20},settings{24}])/1000 > PULSE
            PULSE = max([settings{20},settings{24}])/1000 + 10; % Forcibly push the parameter update to occur after amp settle has a chance to turn off, just in case
        end
    end
    MAX_CHARGE = AMP(end) * DUR(end) * 1e-12;
    if (strcmp(E_MAT,'IrOx'))
        SAFE_CHARGE = IrOx_safe(E_DIAM);
    end
    if (MAX_CHARGE > SAFE_CHARGE)
        disp('The amount of charge indicated exceeds the safe charge threshold.');
        disp('Please try again');
        continue;
    end
    delay = 325; % This magic number comes from the mean jitter (125 msec/trial) plus the parameter update delay inherent to Intan (~200 msec/trial)
    DURATION = length(CHN)*length(DUR)*length(AMP)*n_REP*length(FREQ)*nTRAIN*(PARAM+PULSE+delay);
    minutes = floor((DURATION)/60000);
    seconds = rem((DURATION),60000)/1000;
    disp(['The duration of this experiment set is approximately ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
    HAPPY = input('Is this duration acceptable? Type "1" for YES.\n');
end
%% Initialise
E_MAP = ProbeMAP;
E_MAP = E_MAP(:,3); % USE THIS FIELD TO ADJUST THE ELECTRODE MAPPING
%% Design the settings for each trial
settings = newSettings; % This function generates a stimulation settings list of default values
if max([settings{20},settings{24}])/1000 > PULSE
    PULSE = max([settings{20},settings{24}])/1000 + 10; % Forcibly push the parameter update to occur after amp settle has a chance to turn off, just in case
end
n_Trials = length(CHN)*length(DUR)*length(AMP)*n_REP*length(FREQ)*nTRAIN;
temp = cell(n_Trials,27);
TRAIN = 1:nTRAIN;
FREQ = (1e6)./FREQ;
for C = 1:length(CHN)
    for D = 1:length(DUR)
        for A = 1:length(AMP)
            for N = 1:n_REP
                for P = 1:nTRAIN
                    for F = 1:length(FREQ)
                        trial_number = ((C-1)*length(DUR)*length(AMP)*n_REP*nTRAIN*length(FREQ)) + ((D-1)*length(AMP)*n_REP*nTRAIN*length(FREQ)) + ((A-1)*n_REP*nTRAIN*length(FREQ)) + ((N-1)*nTRAIN*length(FREQ)) + (P-1)*length(FREQ) + F;
                        settings{1} = E_MAP{CHN(C)+1,1};  % Channel number
                        settings{7} = isTRAIN; % Pulses or Train
                        settings{8} = TRAIN(P); % Number of pulses
                        settings{9} = FREQ(F); % Pulse Train Period
                        settings{13} = DUR(D); % First phase duration
                        settings{14} = DUR(D); % Second phase duration
                        settings{15} = IPD;    % Inter-phase delay
                        settings{16} = AMP(A); % First phase amplitude
                        settings{17} = AMP(A); % Second phase amplitude
                        for n = 1:24
                            temp{trial_number,n} = settings{n};
                        end
                        temp{trial_number,25} = ((C-1)*length(DUR)*length(AMP)*nTRAIN*length(FREQ)) + ((D-1)*length(AMP)*nTRAIN*length(FREQ)) + ((A-1)*nTRAIN*length(FREQ)) + ((P-1)*length(FREQ)) + F;
                        temp{trial_number,26} = CHN(C);
                        temp{trial_number,27} = trial_number;
                    end
                end
            end
        end
    end
end

%% Randomize the order of trials
rand_order = randperm(n_Trials);
StimParams = newStimParams(n_Trials);
TrialParams = newTrialParams(n_Trials);
for n = 1:24
    StimParams(2:end,n) = temp(rand_order,n);
end
TrialParams(2:end,2) = temp(rand_order,25);
TrialParams(2:end,3) = temp(rand_order,26);
TrialParams(2:end,1) = temp(:,27);

% Save used datafiles to a new folder to preserve what happened.
here = pwd;
fileID = '001';
C = strsplit(DIR,'\');
NAME = C{size(C,2)};
NEWFILE = strcat(NAME,'_exp_datafile_',fileID,'.mat');
id = 1;
cd(DIR);
dataPlace = pwd;
while exist(strcat(dataPlace,'\',NEWFILE),'file')
    id = id + 1;
    if (id < 10)
        fileID = strcat('00',num2str(id));
    else
        fileID = strcat('0',num2str(id));
    end
    NEWFILE = strcat(NAME,'_exp_datafile_',fileID,'.mat');
end
FREQ = FREQ./(1e6);
save(NEWFILE,'TrialParams','StimParams','n_Trials','E_MAP','E_MAT','E_DIAM','n_REP','rand_order','AMP','DUR','CHN','PARAM','PULSE','FREQ','TRAIN');
disp(['Experimental datafile: ' NAME '_exp_datafile_' fileID '.mat has been saved']);
cd (here);
end