%% Script for generating experimental protocols
function designExp(DIR,AMP,DUR,IPD,E_DIAM,E_MAT,n_REP,CHN,PARAM,PULSE)
% Arguments: AMP: A range of stimulation amplitudes. Expects e.g. [1:2:19],
% in units of uA. DUR: A range of stimulation durations. Expects e.g.
% [100:100:400], in units of us. IPD: Inter-phase delay. Expects e.g. 100,
% in units of us. E_DIAM: Electrode area diameter. Expects e.g. 25, in
% units of um. E_MAT: Electrode material. Expects e.g.'IrOx', in string
% format. List of currently available materials: 'IrOx'. E_IMP: Electrode
% impedance. Expects e.g. 'E_IMP.mat', as a file handle to a struct
% containing impedance values. REPEAT: The number of repeats per trial
% condition. Expects e.g. 75. CHN: A cell array range of electrodes to be
% used for stimulation. Expects e.g. [1,4,7]. Number corresponds with
% designation in E_MAP, e.g. CHN 1 = 'A-019' in Intan. Note that this
% relies on the MAP file - be very careful to select the right electrodes.
% ID: A unique identifier for this stimulation set. Expects e.g. '001'
% PARAM: How long to wait after a parameter update. Expects e.g. 500 (msec)
% PULSE: How long to wait after a stimulation pulse. Expects e.g. 1000
% (msec)

%Dual stim order of variation
% 1. Varation of amplitude - 0/100 25/75 50/50 75/25 100/0
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
        E_DIAM = 15; %electrode diameter
        E_MAT = 'IrOx';
        E_Mapnumber = input('Please enter "1" for four shank or "0" for single shank: \n');
        if E_Mapnumber>0
            nChn=64;
        else
            nChn=32;
        end

        %%
        AMP = input('Please enter stimulus amplitude in uA like so: [x:y:z]\n');
        DUR = input('Please enter stimulus duration in us like so: [x:y:z]\n');
        MAX_CHARGE = AMP(end) * DUR(end) * 1e-6;
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
        n_REP_true=n_REP;
        isTRAIN = input('Please enter a "1" for stimulus train or "0" for single pulse\n');
        if isTRAIN
            nTRAIN = input('Please enter the maximum number of stimulus pulses\n');
            FREQ = input('Please enter the frequencies of the train like so: [x,y,z]\n');
        else
            nTRAIN = 1;
            FREQ = 100;
        end
        if E_Mapnumber==0
            CHN = input('Please enter channels like so: [x,y,z]\n'); 
        else
            CHN = input('Please enter channels like so: S1E5,S2E10,...,S3E13\n','s'); % channels on port B are +32 shanks are 16 electrodes each
            shank_electrode=split(CHN,["S","E",","]);
            shank_electrode=(cellfun(@(x) str2double(x),shank_electrode,'UniformOutput',false));
            shank_electrode=cell2mat(shank_electrode);
            shank_electrode=shank_electrode(~isnan(shank_electrode));
            loopcounter=0;
            CHN=zeros(1,length(shank_electrode)/2);
            for shank=1:2:length(shank_electrode)
                loopcounter=loopcounter+1;
                if shank_electrode(shank)==1
                    CHN(loopcounter)=shank_electrode(shank+1);
                elseif shank_electrode(shank)==2
                    CHN(loopcounter)=shank_electrode(shank+1)+32;
                elseif shank_electrode(shank)==3
                    CHN(loopcounter)=shank_electrode(shank+1)+32+16;
                elseif shank_electrode(shank)==4
                    CHN(loopcounter)=shank_electrode(shank+1)+16;
                elseif shank_electrode(shank)>4
                    error('Shank does not exist')
                end                
            end
        end
        
       DUALSTIM= input('Please enter whether you would like to stimulate on two electrodes YES=1 NO=0: \n');
       if DUALSTIM==1
           if E_Mapnumber==1
               lam_cross = input('Please enter "1" for laminar or "0" for across shanks: \n');
           end
           if lam_cross==0
               NORECORDELECT=CHN(2);
               CHN=CHN(1);
           else
               NORECORDELECT= input('Please enter the number of recording electrodes you would like in between the stimulating electrodes like so: [x,y,z]\n');
               lam_cross=-1;
               checklength=CHN+NORECORDELECT(end)+1; %if electrode spacing causes it to exceed probe length
               n_REP=n_REP*length(NORECORDELECT);
               if any(checklength(:) > nChn)%probe length
                   disp('The requested electrode spacing is not possible.');
                   error('Please try again');
               end
           end
            singanddual= input('Please enter whether you would like to stimulate with both single electrodes and multiple electrodes YES=1 NO=0: \n');
            
            if singanddual==1
                %percentsingdual=5;
                
                %             percentsingdual=input('Please enter number of dual electrode stimulation trials out of 10 i.e. input of "7"=7dual:3single: \n');
                %             n_REP=ceil((n_REP/(10-percentsingdual))+((n_REP/(10-percentsingdual))*percentsingdual));
                varamplitude= input('Please enter whether you would like to vary the amplitude between the stimulating electrodes during the dual stimulation experiment YES=1 NO=0: \n'); %default 0/100 25/75 50/50 75/25 100/0
                if varamplitude==1
                    n_REP=n_REP*5-(n_REP_true*(length(NORECORDELECT)-1));
                    MissAmpFill= input('Please enter whether you would like to fill in missing single stim amplitudes YES=1 NO=0: \n'); %i.e. the 75% amplitude may not be accounted for

                    %             percentsingdual=input('Please enter number of dual electrode stimulation trials out of 10 i.e. input of "7"=7dual:3single: \n');
                    %             n_REP=ceil((n_REP/(10-percentsingdual))+((n_REP/(10-percentsingdual))*percentsingdual));
                else
                    MissAmpFill=0;
                    n_REP=n_REP+n_REP_true*length(NORECORDELECT)+n_REP_true;
                    NORECORDELECT=[NORECORDELECT, 0, -1.*NORECORDELECT];
                end
            else 
                MissAmpFill=0;
                varamplitude=0;
                lam_cross=-1;
            end

        else
            MissAmpFill=0;
            singanddual=0;
            varamplitude=0;
            NORECORDELECT=-1;
            DUALSTIM=1;
            lam_cross=-1;
        end
        % How long do we need to recover from a parameter update: about 200
        % msec seems normal
        %PARAM = input('Please enter delay after parameter update in ms like so: x\n');
        PARAM = 520;%520
        % How long do we need to recover from a stimulation pulse: about
        % 150 msec seems normal
        %PULSE = input('Please enter delay after stimulation in ms like so: x\n');
        PULSE = 150;%150
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
        E_MAP = E_MAP(:,E_Mapnumber+5);
        
        if isempty(E_MAP{35})
            E_MAP=E_MAP(1:33);
        end

%% Design the settings for each trial
settings = newSettings; % This function generates a stimulation settings list of default values
n_Trials = length(CHN)*length(DUR)*length(AMP)*n_REP*length(FREQ)*nTRAIN;
temp = cell(n_Trials,27);
TRAIN = 1:nTRAIN;
FREQ = (1e6)./FREQ;
for C = 1:length(CHN)
    for D = 1:length(DUR)
        for A = 1:length(AMP)
            for P = 1:nTRAIN
                for F = 1:length(FREQ)
                    for N = 1:n_REP
                        %trial_number = ((C-1)*length(DUR)*length(AMP)*n_REP*nTRAIN*length(FREQ)) + ((D-1)*length(AMP)*n_REP*nTRAIN*length(FREQ)) + ((A-1)*n_REP*nTRAIN*length(FREQ)) + ((N-1)*nTRAIN*length(FREQ)) + (P-1)*length(FREQ) + F;
                        trial_number = (((C-1)*length(DUR)*length(AMP)*nTRAIN*length(FREQ)*n_REP)+(((D-1)*length(AMP)*nTRAIN*length(FREQ)*n_REP))+((A-1)*nTRAIN*length(FREQ)*n_REP)+(((P-1)*length(FREQ)*n_REP)+((F-1)*n_REP)))+N;
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
                        %temp{trial_number,25} = ((C-1)*length(DUR)*length(AMP)) + ((D-1)*length(AMP)) + A;
                        temp{trial_number,25} = (((C-1)*length(DUR)*length(AMP)*nTRAIN*length(FREQ))+(((D-1)*length(AMP)*nTRAIN*length(FREQ)))+((A-1)*nTRAIN*length(FREQ))+(((P-1)*length(FREQ))+((F-1))))+1;
                        %temp{trial_number,25} = (((F-1)*length(DUR)*length(AMP)*nTRAIN*length(FREQ))+(((D-1)*length(AMP)*nTRAIN*length(FREQ)))+((A-1)*nTRAIN*length(FREQ))+(((P-1)*length(FREQ))+((F-1))))+1;

                        temp{trial_number,26} = CHN(C);
                        temp{trial_number,27} = trial_number;
                    end
                end
            end
        end
    end
end


%% Randomize the order of trials
StimParams = newStimParams(n_Trials);
TrialParams = newTrialParams(n_Trials);
for n = 1:24
    StimParams(2:end,n) = temp(1:n_Trials,n);
end

TrialParams(2:end,2) = temp(1:n_Trials,25);
TrialParams(2:end,3) = temp(1:n_Trials,26);
TrialParams(2:end,1) = temp(:,27);
maxID=max(cell2mat(TrialParams(2:end,2)));

if DUALSTIM==1
    StimParams=repelem(StimParams,2,1);
    StimParams(1,:)=[];
    TrialParams=repelem(TrialParams,2,1);
    TrialParams(1,:)=[];
    if singanddual==1 %creates randomised vector of ones and zeros according to ratio to enable both single and dual stim in one experiemnt
        if varamplitude==1 %for varing the amplitude of current between stimulating electrodes
%             NumberOfElements = (length(StimParams(:,1))-1)/2;
%             numberOfOnes = round((NumberOfElements*length(NORECORDELECT)/(length(NORECORDELECT)+1)) * 20 / 100); % number of dual electrode params+number or varying amplitude
%             
            % Make initial signal with proper number of 0's and 1's.
            temp1 = [4.*ones(1, n_REP_true),3.*ones(1, n_REP_true),2.*ones(1, n_REP_true),ones(1, n_REP_true), zeros(1, n_REP_true)];
            temp2 = [4.*ones(1, n_REP_true),3.*ones(1, n_REP_true),2.*ones(1, n_REP_true),ones(1, n_REP_true)]; %to remove repetitions
            % Scramble them up with randperm
            %singledualpulsesignal = singledualpulsesignal(randperm(length(singledualpulsesignal)));
            
            singledualpulsesignal=temp1;
            if (length(NORECORDELECT)>1)
                for j=1:(length(StimParams)-1)/(length(temp1)*length(NORECORDELECT)-n_REP_true*(length(NORECORDELECT)-1))
                    for i=1:(length(NORECORDELECT)-1)
                        singledualpulsesignal =[singledualpulsesignal, temp2];
                    end
                    singledualpulsesignal =[singledualpulsesignal, temp1];
                end
            else
                for j=1:((length(StimParams)-1)/(length(temp1))-1)
%                     for i=1:(length(NORECORDELECT)-1)
%                         singledualpulsesignal =[singledualpulsesignal, temp2];
%                     end
                    singledualpulsesignal =[singledualpulsesignal, temp1];
                end
            end

%         else
%             singledualpulsesignal = [4.*ones(1, n_REP_true),3.*ones(1, n_REP_true),2.*ones(1, n_REP_true),ones(1, n_REP_true), zeros(1, n_REP_true)];
%              if (length(NORECORDELECT)-1)>1
%                  temp1=singledualpulsesignal;
%                  for i=1:(length(StimParams)-1)/length(temp1)
%                      singledualpulsesignal =[singledualpulsesignal, temp1];
%                  end
%              end
%             %numberOfOnes = round(length(StimParams(:,1))/(length(NORECORDELECT)*maxID))/2;%gives even number of trials per number of spacing plus single trials
%             % Make initial signal with proper number of 0's and 1's.
%             singledualpulsesignal = [ones(1, n_REP-n_REP_true), zeros(1, n_REP_true)];
%             % Scramble them up with randperm
%             %singledualpulsesignal = singledualpulsesignal(randperm(length(singledualpulsesignal)));
%             if length(NORECORDELECT)~=1
%                 temp1=singledualpulsesignal;
%                 for i=1:((length(StimParams)-1)/2/length(temp1))-1
%                     singledualpulsesignal =[singledualpulsesignal, temp1];
%                 end
%             end
        end
    end
        counter=1;
        electrecord=1;
        ampcount=1;
        numberofrecordelec=round(length(StimParams(:,1))/(length(NORECORDELECT)*maxID));
        for i=3:2:(length(StimParams(:,1)))
            if (i==numberofrecordelec*counter+3)
                if electrecord==length(NORECORDELECT)
                    ampcount=ampcount+1;
                    electrecord=0;
                end
                electrecord=electrecord+1;
                counter=counter+1;
            end
            
%             electrecord=electrecord+1;
%             if electrecord>length(NORECORDELECT)
%                 electrecord=1;
%             end
            chnstr=StimParams(i,1);
            z=0;
            count=0;
            while ~z
                count=count+1;
                z=strcmp(chnstr,E_MAP(count,1));
            end
            if singanddual==0

                StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
                TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
                if (exist('NORECORDELECT(2)'))
                    TrialParams(i-1,2)={(counter)};
                    TrialParams(i,2)={(counter)};%assign new parameters trial iDs based on original trial numbers
                end
             elseif singanddual==1
                if varamplitude==0
                    if NORECORDELECT(electrecord)>0
                        StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
                        TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
                    elseif NORECORDELECT(electrecord)==0
                         TrialParams(i,3)={0};
                    elseif NORECORDELECT(electrecord)<0
                        StimParams(i,1)=E_MAP(count+(-1*NORECORDELECT(electrecord))+1,1);
                        TrialParams(i,3)={TrialParams{i,3}+(-1*NORECORDELECT(electrecord))+1};
                        StimParams(i-1,1)=StimParams(i,1);%assigns parameters of second renamed electrode to the original channel for 100%
                        TrialParams(i-1,3)=TrialParams(i,3);
                        TrialParams(i,3)={0};
                    end

                    TrialParams(i-1,2)={floor((i-3)/(n_REP_true*2))+1};
                    TrialParams(i,2)={floor((i-3)/(n_REP_true*2))+1};%assign new parameters trial iDs based on original trial numbers

                else
                    if singledualpulsesignal(i/2-0.5)==0
                        TrialParams(i,3)={0};
%                     elseif singledualpulsesignal(i/2-0.5)==5
%                         if StimParams{i,16}~=-1
%                             StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
%                             TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
%                         end
%                         TrialParams(i-1,2)={cell2mat(TrialParams(i,2))+maxID*electrecord};
%                         TrialParams(i,2)={cell2mat(TrialParams(i,2))+maxID*electrecord};%assign new parameters trial iDs based on original trial numbers
                    elseif singledualpulsesignal(i/2-0.5)==1
                       % if StimParams{i,16}~=-1
                            if lam_cross==0
                                StimParams(i,1)=E_MAP(NORECORDELECT+1,1);
                                TrialParams(i,3)={NORECORDELECT};
                            else
                            StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
                            TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
                            end
                            
                            if StimParams{i,16}~=-1
                                StimParams{i,16}=StimParams{i,16}*0.25;%assigns amplitude of 25%
                                StimParams{i,17}=StimParams{i,16};%ensures both pos and neg phase amplitude balanced
                                StimParams{i-1,16}=StimParams{i-1,16}*0.75;%assigns amplitudes of 75%
                                StimParams{i-1,17}=StimParams{i-1,16};%ensures both pos and neg phase amplitude balanced
                            end
                            
%                         TrialParams(i-1,2)={cell2mat(TrialParams(i,2))+maxID*electrecord};
%                         TrialParams(i,2)={cell2mat(TrialParams(i,2))+maxID*electrecord};%assign new parameters trial iDs based on original trial numbers
                    elseif singledualpulsesignal(i/2-0.5)==2
                        % if StimParams{i,16}~=-1
                        if lam_cross==0
                            StimParams(i,1)=E_MAP(NORECORDELECT+1,1);
                            TrialParams(i,3)={NORECORDELECT};
                        else
                            StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
                            TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
                        end

                            if StimParams{i,16}~=-1
                                StimParams{i,16}=StimParams{i,16}*0.5;%assigns amplitude
                                StimParams{i,17}=StimParams{i,16};%ensures both pos and neg phase amplitude balanced
                                StimParams{i-1,16}=StimParams{i-1,16}*0.5;%assigns amplitudes
                                StimParams{i-1,17}=StimParams{i-1,16};%ensures both pos and neg phase amplitude balanced
                            end
%                         TrialParams(i-1,2)={cell2mat(TrialParams(i,2))+maxID*2*electrecord};
%                         TrialParams(i,2)={cell2mat(TrialParams(i,2))+maxID*2*electrecord};%assign new parameters trial iDs based on original trial numbers
                    elseif singledualpulsesignal(i/2-0.5)==3
                      %  if StimParams{i,16}~=-1
                        if lam_cross==0
                            StimParams(i,1)=E_MAP(NORECORDELECT+1,1);
                            TrialParams(i,3)={NORECORDELECT};
                        else
                            StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
                            TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
                        end
                            if StimParams{i,16}~=-1
                                StimParams{i,16}=StimParams{i,16}*0.75;%assigns amplitude
                                StimParams{i,17}=StimParams{i,16};%ensures both pos and neg phase amplitude balanced
                                StimParams{i-1,16}=StimParams{i-1,16}*0.25;%assigns amplitudes
                                StimParams{i-1,17}=StimParams{i-1,16};%ensures both pos and neg phase amplitude balanced
                            end
%                         TrialParams(i-1,2)={cell2mat(TrialParams(i,2))+maxID*3*electrecord};
%                         TrialParams(i,2)={cell2mat(TrialParams(i,2))+maxID*3*electrecord};%assign new parameters trial iDs based on original trial numbers
                    elseif singledualpulsesignal(i/2-0.5)==4
%                        if StimParams{i,16}~=-1
                        if lam_cross==0
                            StimParams(i,1)=E_MAP(NORECORDELECT+1,1);
                            TrialParams(i,3)={NORECORDELECT};
                        else
                            StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
                            TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
                        end
                            StimParams(i-1,1)=StimParams(i,1);%assigns parameters of second renamed electrode to the original channel for 100%
                            TrialParams(i-1,3)=TrialParams(i,3);
                            TrialParams(i,3)={0};
%                        end
%                         TrialParams(i-1,2)={cell2mat(TrialParams(i,2))+maxID*4*electrecord};
%                         TrialParams(i,2)={cell2mat(TrialParams(i,2))+maxID*4*electrecord};%assign new parameters trial iDs based on original trial numbers
                    end
                    TrialParams(i-1,2)={floor((i-3)/(n_REP_true*2))+1};
                    TrialParams(i,2)={floor((i-3)/(n_REP_true*2))+1};%assign new parameters trial iDs based on original trial numbers
                end
            else
                if lam_cross==0
                    StimParams(i,1)=E_MAP(NORECORDELECT+1,1);
                    TrialParams(i,3)={NORECORDELECT};
                else
                    StimParams(i,1)=E_MAP(count+NORECORDELECT(electrecord)+1,1);
                    TrialParams(i,3)={TrialParams{i,3}+NORECORDELECT(electrecord)+1};
                end
            end
        end

end


%% For missing single trials enabling a prediction to be made
if MissAmpFill==1
    AllstimAMP = unique(cell2mat(StimParams(2:end,16)));
    MissingSingleAMP = setdiff(AllstimAMP,AMP);
    stimCHN = unique(cell2mat(TrialParams(2:end,3)));
    stimCHN(stimCHN==0)=[];
    originalEND=size(StimParams,1);
    StimParams=[StimParams;repmat(StimParams(end,:),n_REP_true*length(MissingSingleAMP)*length(stimCHN)*2,1)];
    TrialParams=[TrialParams;repmat(TrialParams(end,:),n_REP_true*length(MissingSingleAMP)*length(stimCHN)*2,1)];
    ar=(originalEND-1)/2+1:1:(size(StimParams,1)-1)/2;
    ar = repelem(ar,2)';
    TrialParams(originalEND+1:end,1)=num2cell(ar);
    TrialParams(originalEND+1:end,3)={0};
    for chosenchn=1:length(stimCHN)
        for chosenamp=1:length(MissingSingleAMP)
            StimParams(originalEND+n_REP_true*(chosenamp-1)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP)+1:originalEND+n_REP_true*(chosenamp)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP),16)={MissingSingleAMP(chosenamp)};
            StimParams(originalEND+n_REP_true*(chosenamp-1)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP)+1:originalEND+n_REP_true*(chosenamp)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP),17)={MissingSingleAMP(chosenamp)};
            StimParams(originalEND+n_REP_true*(chosenamp-1)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP)+1:originalEND+n_REP_true*(chosenamp)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP),1)={E_MAP{stimCHN(chosenchn)+1,1}};
            TrialParams(originalEND+n_REP_true*(chosenamp-1)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP)+1:originalEND+n_REP_true*(chosenamp)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP),2)={cell2mat(TrialParams(originalEND,2))+(chosenamp)+(chosenchn-1)*length(MissingSingleAMP)};
            TrialParams(originalEND+n_REP_true*(chosenamp-1)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP)+1:2:originalEND+n_REP_true*(chosenamp)*2+n_REP_true*(chosenchn-1)*2*length(MissingSingleAMP),3)={stimCHN(chosenchn)};
        end
    end
    n_Trials=n_Trials+n_REP_true*length(MissingSingleAMP)*length(stimCHN);
else
    originalEND=size(StimParams,1);
end
%%

if DUALSTIM==1
    rand_order = randperm(n_Trials).*2;
    rand_order_1 = rand_order-1;
    rand_order=[rand_order_1;rand_order];
    rand_order=rand_order(:)';
    temp=StimParams(1,:);
    StimParams(1,:)=[];

    StimParams = StimParams(rand_order,:);
   
    StimParams=[temp(1,:);StimParams(:,:)];

temp=TrialParams(1,:);
TrialParams(1,:)=[];

TrialParams = TrialParams(rand_order,:);
TrialParams=[temp(1,:); TrialParams(:,:)];
    if NORECORDELECT==-1
        TrialParams(3:2:end,3)={0};
    end
else
    rand_order = randperm(n_Trials);
    for n=1:24
        StimParams(2:end,n) = temp(rand_order,n);
    end
    
    
    TrialParams(2:end,2) = temp(rand_order,25);
    TrialParams(2:end,3) = temp(rand_order,26);
    TrialParams(2:end,1) = temp(:,27);
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
if ~(exist('NORECORDELECT(2)')) &&  (NORECORDELECT(1)==-1)
    DUALSTIM=0;
end
if DUALSTIM==1
    save(NEWFILE,'TrialParams','StimParams','n_Trials','E_MAP','NORECORDELECT','E_MAT','E_DIAM','n_REP','n_REP_true','rand_order','AMP','DUR','CHN','PARAM','PULSE','FREQ','TRAIN','DUALSTIM','varamplitude','singanddual','E_Mapnumber','originalEND');
else
    save(NEWFILE,'TrialParams','StimParams','n_Trials','E_MAP','E_MAT','E_DIAM','n_REP','n_REP_true','rand_order','AMP','DUR','CHN','PARAM','PULSE','FREQ','TRAIN','DUALSTIM','varamplitude','singanddual','E_Mapnumber','originalEND');

end
disp(['Experimental datafile: ' NAME '_exp_datafile_' fileID '.mat has been saved']);
cd (here);
end