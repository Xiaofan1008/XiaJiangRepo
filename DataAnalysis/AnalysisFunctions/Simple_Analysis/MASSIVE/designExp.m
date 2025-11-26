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
HAPPY = 1;
while (HAPPY)
    if nargin == 0 %determines number of input functions
        DIR = uigetdir('\\ad.monash.edu\home\User042\smei0006\Documents\2020 Research\Experimental Design\TESTING CODE'); %starts new save directory
        CASE = 0;
    elseif nargin == 1 %if you already specified the directory when calling the function
        CASE = 0;
    else
        CASE = 1; %catch
    end
    if ~(CASE)
        E_DIAM = 25; %electrode diameter
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
        
        
       DUALSTIM= input('Please enter whether you would like to stimulate on two electrodes YES=1 NO=0: \n');
       if DUALSTIM==1
           NORECORDELECT= input('Please enter the number of recording electrodes you would like in between the stimulating electrodes like so: x\n');
           
           checklength=CHN+NORECORDELECT+1; %if electrode spacing causes it to exceed probe length
           if any(checklength(:) > 32)%probe length
               disp('The requested electrode spacing is not possible.');
               error('Please try again');
           end
           singanddual= input('Please enter whether you would like to stimulate with both single electrodes and multiple electrodes YES=1 NO=0: \n');
           if singanddual==1
               percentsingdual=input('Please enter number of dual electrode stimulation trials out of 10 i.e. input of "7"=7dual:3single: \n');
           end
        end
        % How long do we need to recover from a parameter update: about 200
        % msec seems normal
        %PARAM = input('Please enter delay after parameter update in ms like so: x\n');
        PARAM = 520;
        % How long do we need to recover from a stimulation pulse: about
        % 150 msec seems normal
        %PULSE = input('Please enter delay after stimulation in ms like so: x\n');
        PULSE = 150;
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
    DURATION = length(CHN)*length(DUR)*length(AMP)*n_REP*length(FREQ)*nTRAIN*(PARAM+PULSE+125);
    minutes = floor((DURATION)/60000);
    seconds = rem((DURATION),60000)/1000;
    disp(['The duration of this experiment set is approximately ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
    HAPPY = input('Is this duration acceptable? Type "0" for YES.\n');
end
%% Initialise
E_MAP = ProbeMAP;
E_MAP = E_MAP(:,3); % USE THIS FIELD TO ADJUST THE ELECTRODE MAPPING
%% Design the settings for each trial
settings = newSettings; % This function generates a stimulation settings list of default values
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
                        temp{trial_number,25} = ((C-1)*length(DUR)*length(AMP)) + ((D-1)*length(AMP)) + A;
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
maxID=max(cell2mat(TrialParams(2:end,2)));
if DUALSTIM==1
    StimParams=repelem(StimParams,2,1);
    StimParams(1,:)=[];
    TrialParams=repelem(TrialParams,2,1);
    TrialParams(1,:)=[];
    if singanddual==1 %creates randomised vector of ones and zeros according to ratio to enable both single and dual stim in one experiemnt
        NumberOfElements = (length(StimParams(:,1))-1)/2;
        numberOfOnes = round(NumberOfElements * percentsingdual*10 / 100);
        % Make initial signal with proper number of 0's and 1's.
        singledualpulsesignal = [ones(1, numberOfOnes), zeros(1, NumberOfElements - numberOfOnes)];
        % Scramble them up with randperm
        singledualpulsesignal = singledualpulsesignal(randperm(length(singledualpulsesignal)));
    end
    for i=3:2:length(StimParams(:,1))
        chnstr=StimParams(i,1);
        z=0;
        count=0;
        while ~z
            count=count+1;
            z=strcmp(chnstr,E_MAP(count,1));
        end
        if singanddual==0
            StimParams(i,1)=E_MAP(count+NORECORDELECT+1,1);
            TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT+1};
        elseif singanddual==1
            if singledualpulsesignal(i/2-0.5)==0
                StimParams(i,1)={F};
                TrialParams(i,3)={0};
            else
                StimParams(i,1)=E_MAP(count+NORECORDELECT+1,1);
                TrialParams(i-1,2)={cell2mat(TrialParams(i,2))+maxID};
                TrialParams(i,2)={cell2mat(TrialParams(i,2))+maxID};%assign new parameters trial iDs based on original trial numbers 
                TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT+1};
            end
        else
            StimParams(i,1)=E_MAP(count+NORECORDELECT+1,1);
            TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT+1};
        end
    end

end




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