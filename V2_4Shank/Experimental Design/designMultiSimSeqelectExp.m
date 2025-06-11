%% Script for generating experimental protocols
function designMultiSimSeqelectExp(DIR,AMP,DUR,IPD,E_DIAM,E_MAT,n_REP,CHN,PARAM,PULSE,PTD)
% Arguments: AMP: A range of stimulation amplitudes. Expects e.g. [1:2:19],in units of uA. 
% DUR: A range of stimulation durations. Expects e.g.[100:100:400], in units of us. 
% IPD: Inter-phase delay. Expects e.g. 100,in units of us. 
% E_DIAM: Electrode area diameter. Expects e.g. 25, in units of um. 
% E_MAT: Electrode material. Expects e.g.'IrOx', in string
% format. List of currently available materials: 'IrOx'. 
% E_IMP: Electrode impedance. Expects e.g. 'E_IMP.mat', as a file handle to a struct containing impedance values. 
% REPEAT: The number of repeats per trial condition. Expects e.g. 75. 
% CHN: A cell array range of electrodes to be used for stimulation. Expects e.g. [1,4,7]. Number corresponds with designation in 
% E_MAP, e.g. CHN 1 = 'A-019' in Intan. Note that this relies on the MAP file - be very careful to select the right electrodes.
% ID: A unique identifier for this stimulation set. Expects e.g. '001'
% PARAM: How long to wait after a parameter update. Expects e.g. 500 (msec)
% PULSE: How long to wait after a stimulation pulse. Expects e.g. 1000 (msec)
% PTD: post trigger delay: Time delay for the pulse after the trigger. 

% Dual stim order of variation
% 1. Variation of amplitude - 0/100 25/75 50/50 75/25 100/0
% 2. Pulse train period
% 3. Number of pulses
% 4. Amplitude - your chosen levels
% 5. Duration
% 6. Channel


HAPPY = 1;
while (HAPPY)
    if nargin == 0 %determines number of input functions
        DIR = uigetdir('Z:\Shared\Sabrina\Experimental Design\TESTING CODE'); %starts new save directory
        CASE = 0;
    elseif nargin == 1 %if you already specified the directory when calling the function
        CASE = 0;
    else
        CASE = 1; %catch
    end
    if ~(CASE)
        %% Channel map handling
        E_DIAM = 15; %electrode diameter
        E_MAT = 'IrOx';
        E_Mapnumber = input('Please enter "1" for four shank or "0" for single shank: \n');
        if E_Mapnumber>0
            nChn=64;
        else
            nChn=32;
        end

        %% Amplitude
        AMP = input('Please enter stimulus amplitude in uA like so: [x:y:z]\n');
        if any(AMP==-1) || any(AMP==0)
            zeroampflag=1;
            AMP(AMP==-1)=[];
            AMP(AMP==0)=[];
        else
            zeroampflag=0;
        end

        %% Duration
        DUR = 200; 
        % DUR = input('Please enter stimulus duration in us like so: [x:y:z]\n');
        MAX_CHARGE = AMP(end) * DUR(end) * 1e-6;
        if (strcmp(E_MAT,'IrOx'))
            SAFE_CHARGE = 2124 * 1e-6; %IrOx_safe(E_DIAM);
        end
        if (MAX_CHARGE > SAFE_CHARGE)
            disp('The amount of charge indicated exceeds the safe charge threshold.');
            fprintf('Your charge: %.6f, Safe charge: %.6f\n', MAX_CHARGE, SAFE_CHARGE);
            disp('Please try again');
            continue;
        end

        %% Interphase Interval
        IPD = 100;
        % IPD = input('Please enter interphase delay like so: x\n');
        n_REP = input('Please enter number of repeats like so: x\n');

        %% Pulse trian or Signle pulse
        isTRAIN = input('Please enter a "1" for stimulus train or "0" for single pulse\n');
        if isTRAIN
            % nTRAIN = 1;
            nTRAIN = input('Please enter the maximum number of stimulus pulses\n');
            FREQ = input('Please enter the frequencies of the train in Hz like so: [x,y,z]\n');
            numpulses = input('Please enter the number of stimulus pulses\n');
        else
            numpulses =1;
            nTRAIN = 1;
            FREQ = 100;
        end
        simultaneous_stim= input('How many channels would you like to stimulate on for each tiral ("1" or "2" or "3"...)? \n');
        if simultaneous_stim==1 % Signle channel stimulation
            if E_Mapnumber==0 % single shank
                CHN = input('Please enter channels like so: [x,y,z]\n');
            else % four shank 
                %%
                CHN = input('Please enter channels like so: PAS1E5,PBS2E10,...,PDS3E13 \nWhere P is port letter, S is #shank, E is #electrode\n','s'); % channels on port B are +32 shanks are 16 electrodes each
                port_shank_electrode=split(CHN,["P","S","E",","]);
                shank_electrode=(cellfun(@(x) str2double(x),port_shank_electrode,'UniformOutput',false));
                Portnum=(cellfun(@(x) x(isletter(x)),port_shank_electrode,'UniformOutput',false));
                Portnum=Portnum(~cellfun('isempty',Portnum));
                shank_electrode=cell2mat(shank_electrode);
                shank_electrode=shank_electrode(~isnan(shank_electrode));
                loopcounter=0;
                CHN=zeros(1,length(shank_electrode)/2);
                for shank=1:2:length(shank_electrode)
                    loopcounter=loopcounter+1;
                    CHN(loopcounter)=shank_electrode(shank+1);
                    if shank_electrode(shank)==4 || shank_electrode(shank)==3
                        CHN(loopcounter)=CHN(loopcounter)+16;
                    elseif shank_electrode(shank)>4
                        error('Shank does not exist')
                    end
                    if strcmp(Portnum(loopcounter),'B')
                        CHN(loopcounter)=CHN(loopcounter)+32;
                    elseif strcmp(Portnum(loopcounter),'C')
                        CHN(loopcounter)=CHN(loopcounter)+64;
                    elseif strcmp(Portnum(loopcounter),'D')
                        CHN(loopcounter)=CHN(loopcounter)+96;
                    end
                end
            end
            CHN=CHN';
        else % Multi-channel stimulation

            sequential_stim= input('Are you testing Sequential Stimulation, Yes = 1, No = 0 \n');

            checkEndpairs='1';
            loopiterate=1;
            if E_Mapnumber==0 % single shank electrode array
                while (checkEndpairs)~=0
                    if loopiterate==1
                        CHN(loopiterate,:)= input(['Input first ' num2str(simultaneous_stim) ' electrodes for simultaneous stim like so: [x,y,z] \n']);
                    else
                        CHN(loopiterate,:)= input(['Input next ' num2str(simultaneous_stim) ' electrodes for simultaneous stim or 0 if there are no more sets \n']);
                    end
                    checkEndpairs=CHN(loopiterate,1);
                    loopiterate=loopiterate+1;
                end
                CHN(loopiterate-1,:)=[];
                if size(CHN,2)~=simultaneous_stim
                    error("Too many or too few electrodes input")
                end
            else % four shank electrode array


                while str2double(checkEndpairs)~=0
                    if loopiterate==1
                        CHN_temp = input(['Input first ' num2str(simultaneous_stim) ' electrodes for simultaneous stim like so: PAS1E5,PBS2E10,...,PDS3E13 \nWhere P is port letter, S is #shank, E is #electrode\n'],'s'); % channels on port B are +32 shanks are 16 electrodes each
                    else
                        CHN_temp= input(['Input next ' num2str(simultaneous_stim) ' electrodes for simultaneous stim or 0 if there are no more sets \n'],'s');
                    end
                    checkEndpairs=CHN_temp;
                    if str2double(checkEndpairs)==0
                        break;
                    end
                    %convert to electrode numbers
                    port_shank_electrode=split(CHN_temp,["P","S","E",","]);
                    shank_electrode=(cellfun(@(x) str2double(x),port_shank_electrode,'UniformOutput',false));
                    Portnum=(cellfun(@(x) x(isletter(x)),port_shank_electrode,'UniformOutput',false));
                    Portnum=Portnum(~cellfun('isempty',Portnum));
                    shank_electrode=cell2mat(shank_electrode);
                    shank_electrode=shank_electrode(~isnan(shank_electrode));
                    loopcounter=0;
                    CHN_num=zeros(1,length(shank_electrode)/2);
                    for shank=1:2:length(shank_electrode)
                        loopcounter=loopcounter+1;
                        CHN_num(loopcounter)=shank_electrode(shank+1);
                        if shank_electrode(shank)==4 || shank_electrode(shank)==3
                            CHN_num(loopcounter)=CHN_num(loopcounter)+16;
                        elseif shank_electrode(shank)>4
                            error('Shank does not exist')
                        end
                        if strcmp(Portnum(loopcounter),'B')
                            CHN_num(loopcounter)=CHN_num(loopcounter)+32;
                        elseif strcmp(Portnum(loopcounter),'C')
                            CHN_num(loopcounter)=CHN_num(loopcounter)+64;
                        elseif strcmp(Portnum(loopcounter),'D')
                            CHN_num(loopcounter)=CHN_num(loopcounter)+96;
                        end
                    end
                    if length(CHN_num)~=simultaneous_stim
                        error("Too many or too few electrodes input")
                    end
                    CHN(loopiterate,:)=CHN_num;
                    loopiterate=loopiterate+1;
                end

            end
        end

        %% Multi-channel stimualtion setting
        if simultaneous_stim~=1
            if simultaneous_stim==2
                varamplitude= input('Are you testing current steering? YES=1 NO=0: \n'); %default 0/100 25/75 50/50 75/25 100/0
                singanddual=1;
                if varamplitude==0
                    singanddual= input('Do you want whitenoise stim (i.e. single, dual... simultaneous)? YES=1 NO=0: \n');
                end
            else
                varamplitude=0;
                singanddual= input('Do you want whitenoise stim (i.e. single, dual... simultaneous)? YES=1 NO=0: \n');
            end
            if singanddual==1
                if varamplitude==1 % current steering
                    AMPtemp=repmat(AMP',1,size(CHN, 2));
                    AMPvar=[1,0; 0.75,0.25; 0.5,0.5; 0.25,0.75; 0,1];
                    row=1;
                    AMP=[];
                    for i=1:size(AMPtemp,1)
                        AMP(row:row+size(AMPvar,1)-1,:)=AMPtemp(i,:).*AMPvar;
                        row=row+size(AMPvar,1);
                    end
                    missamp=setdiff(AMP,[0; AMPtemp(:,1)]);
                    %fill in missing single amplitudes for comparison
                    AMP(row:row+length(missamp)-1,1:2)=[missamp, zeros(length(missamp),1)];
                    AMP(row+length(missamp):row+length(missamp)+length(missamp)-1,1:2)=[zeros(length(missamp),1), missamp];
                    varamplitude=4;
                else
                    n=1:size(CHN,2);
                    row=0;
                    AMPcombo=zeros(((2^length(AMP))),simultaneous_stim);
                    for k=1:size(CHN,2)
                        permutes=nchoosek(n,k);
                        allcurrents = repelem(AMP, size(permutes,2));
                        for i=1:size(permutes,1)
                            currentcombo=uniqueperms(allcurrents);
                            currentcombo=unique(currentcombo(:,1:k),'rows');
                            for currentcycle=1:size(currentcombo,1)
                                row=row+1;
                                AMPcombo(row,permutes(i,:))=currentcombo(currentcycle,:);
                            end
                        end
                    end
                    AMP=AMPcombo;
                end
            else
                AMP=repmat(AMP',1,size(CHN, 2));
                % scales = [0.25, 0.5, 0.25];
                % AMP = AMP' * scales;
                varamplitude=0;
            end
            %% Post Trigger Delay
            PTD = zeros(1,simultaneous_stim); % Post Trigger Delay
            if sequential_stim == 1 % Sequential stimulation 
                for j = 1:simultaneous_stim
                    PTD(j) = input('Please input the Post Trigger Delay for each electrode in one trail (in µs); \n'); 
                end
            else % Simultaneous stimulation                
            end

        else
            AMP=repmat(AMP',1,size(CHN, 2));
            singanddual=0;
            varamplitude=0;
        end


        % How long do we need to recover from a parameter update: about 200
        % msec seems normal
        %PARAM = input('Please enter delay after parameter update in ms like so: x\n');
        PARAM = 200; %320%520
        % How long do we need to recover from a stimulation pulse: about
        % 150 msec seems normal
        PULSE = 600; %500; %input('Please enter delay after stimulation in ms like so: [x y z]\n');
        %PULSE = 180;%150
    end
    AMP(AMP==0)=-1;
    MAX_CHARGE = AMP(end) * DUR(end) * 1e-12;
    if (strcmp(E_MAT,'IrOx'))
        SAFE_CHARGE = 2124 * 1e-12; %IrOx_safe(E_DIAM);
    end
    if (MAX_CHARGE > SAFE_CHARGE)
        disp('The amount of charge indicated exceeds the safe charge threshold.');
        disp('Please try again');
        continue;
    end
    %AMP = AMP(1:35, :);
    n_Trials = size(CHN,1)*length(DUR)*size(AMP,1)*n_REP*length(FREQ)*nTRAIN*length(PULSE);
    DURATION = n_Trials*(PARAM+mean(PULSE)+125+200)+40000;
    minutes = floor((DURATION)/60000);
    seconds = rem((DURATION),60000)/1000;
    disp(['The duration of this experiment set is approximately ' num2str(minutes) ' minutes and ' num2str(seconds) ' seconds.']);
    HAPPY = input('Is this duration acceptable? Type "0" for YES.\n');
end
%% Initialise
E_MAP = ProbeMAP;
E_MAP = E_MAP(:,E_Mapnumber+5);

if isempty(E_MAP{35})
    E_MAP=E_MAP(1:33);
end

%% Design the settings for each trial
settings = newSettings; % This function generates a stimulation settings list of default values
temp = cell(n_Trials*simultaneous_stim,28);
TRAIN = 1:nTRAIN;
FREQ = (1e6)./FREQ;
%Polarity = {0,1,0};
for T=1:length(PULSE)
    for C = 1:size(CHN,1)
        for D = 1:length(DUR)
            for A = 1:size(AMP,1)
                for P = 1:nTRAIN
                    for F = 1:length(FREQ)
                        for N = 1:n_REP
                            %%%%assign all variation in same loop?
                            %trial_number = ((C-1)*length(DUR)*length(AMP)*n_REP*nTRAIN*length(FREQ)) + ((D-1)*length(AMP)*n_REP*nTRAIN*length(FREQ)) + ((A-1)*n_REP*nTRAIN*length(FREQ)) + ((N-1)*nTRAIN*length(FREQ)) + (P-1)*length(FREQ) + F;
                            trial_number = (((T-1)*size(CHN,1)*length(DUR)*size(AMP,1)*nTRAIN*length(FREQ)*n_REP)+((C-1)*length(DUR)*size(AMP,1)*nTRAIN*length(FREQ)*n_REP)+((D-1)*size(AMP,1)*nTRAIN*length(FREQ)*n_REP)+((A-1)*nTRAIN*length(FREQ)*n_REP)+((P-1)*length(FREQ)*n_REP)+((F-1)*n_REP))+N;
                            settings{1} = E_MAP{CHN(C,1)+1,1};  % Channel number
                            settings{7} = isTRAIN; % Pulses or Train
                            settings{8} = numpulses; %TRAIN(P); % Number of pulses
                            settings{9} = FREQ(F); % Pulse Train Period
                            settings{13} = DUR(D); % First phase duration
                            settings{14} = DUR(D); % Second phase duration
                            settings{15} = IPD;    % Inter-phase delay
                            settings{16} = AMP(A,1); % First phase amplitude
                            settings{17} = AMP(A,1); % Second phase amplitude
                            settings{25} = PULSE(T);%length of recovery time 
                            settings{26} = (((C-1)*length(DUR)*size(AMP,1)*nTRAIN*length(FREQ))+(((D-1)*size(AMP,1)*nTRAIN*length(FREQ)))+((A-1)*nTRAIN*length(FREQ))+(((P-1)*length(FREQ))+((F-1))))+1;
                            settings{27} = CHN(C,1);
                            settings{28} = trial_number;
                            for simstimchns=1:simultaneous_stim % multielectrode stim
                                if CHN(C,simstimchns)~=0
                                    settings{1} = E_MAP{CHN(C,simstimchns)+1,1};  % Channel number
                                else %set flag to not update row to intan
                                    settings{1} = '0-000';  % Channel number
                                end
                                settings{27} = CHN(C,simstimchns);
                                settings{16} = AMP(A,simstimchns); % First phase amplitude
                                settings{17} = AMP(A,simstimchns); % Second phase amplitude
                                %settings{12} = Polarity{simstimchns}; %make polarity of other channels than the first one anodic
                                settings{6} = PTD(simstimchns); % Post-trigger delay
                                for n = 1:28
                                    temp{(trial_number-1)*simultaneous_stim+simstimchns,n} = settings{n};
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%% one set of zero trials
if zeroampflag==1 || singanddual~=0
    temp(trial_number*simultaneous_stim+1:trial_number*simultaneous_stim+n_REP*simultaneous_stim,:)=temp(trial_number*simultaneous_stim-n_REP*simultaneous_stim+1:trial_number*simultaneous_stim,:);
    for i=trial_number*simultaneous_stim+1:simultaneous_stim:(trial_number*simultaneous_stim+n_REP*simultaneous_stim)
        temp(i:i+simultaneous_stim-1,28)={temp{i-1,28}+1};
        temp(i:i+simultaneous_stim-1,26)={temp{trial_number*simultaneous_stim,26}+1};
        temp(i:i+simultaneous_stim-1,16)={-1};
        temp(i:i+simultaneous_stim-1,17)={-1};
    end
    n_Trials=n_Trials+n_REP;
end

%% Assign stim params and trialparams
StimParams = newStimParams(n_Trials.*simultaneous_stim);
TrialParams = newTrialParams(n_Trials.*simultaneous_stim);
for n = 1:25
    StimParams(2:end,n) = temp(1:(n_Trials*simultaneous_stim),n);
end

TrialParams(2:end,2) = temp(1:(n_Trials*simultaneous_stim),26);
TrialParams(2:end,3) = temp(1:(n_Trials*simultaneous_stim),27);
TrialParams(2:end,1) = temp(:,28);
TrialParams([false; cell2mat(StimParams(2:end,16))==0],3)={0};
maxID=max(cell2mat(TrialParams(2:end,2)));

%% Randomize the order of trials
rand_order = randperm(n_Trials).*simultaneous_stim;
rand_order = repelem(rand_order,simultaneous_stim);
for i=1:simultaneous_stim-1
    rand_order(i:simultaneous_stim:end)=rand_order(i:simultaneous_stim:end)-simultaneous_stim+i;
end
temp1=StimParams(1,:);
StimParams(1,:)=[];
StimParams = StimParams(rand_order,:);
StimParams=[temp1(1,:);StimParams(:,:)];

temp1=TrialParams(1,:);
TrialParams(1,:)=[];
TrialParams = TrialParams(rand_order,:);
TrialParams=[temp1(1,:); TrialParams(:,:)];


%% Save used datafiles to a new folder to preserve what happened.
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
save(NEWFILE,'TrialParams','StimParams','n_Trials','E_MAP','E_MAT','E_DIAM','n_REP','rand_order','AMP','DUR','CHN','PARAM','PULSE','FREQ','TRAIN','varamplitude','singanddual','E_Mapnumber','simultaneous_stim');
disp(['Experimental datafile: ' NAME '_exp_datafile_' fileID '.mat has been saved']);
cd (here);
end