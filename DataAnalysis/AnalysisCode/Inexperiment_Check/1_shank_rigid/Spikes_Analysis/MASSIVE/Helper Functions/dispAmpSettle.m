%% Script that demonstrates Amp Settle on each channel in this experiment
E_MAP = ProbeMAP;
E_MAP = E_MAP(:,E_Mapnumber+5);

if isempty(E_MAP{35})
    E_MAP=E_MAP(1:33);
end
CHNlist = cell2mat(unique(TRIAL.StimParams(2:end,1)));

% For each channel . . .

for i = 1:size(CHNlist,1)
        % Configure Intan:
        fprintf(t,['channel=' CHNlist(i,:)]);
        fprintf(t,'enabled=1');
        fprintf(t,'firstPhaseAmplitude=0');
        fprintf(t,'secondPhaseAmplitude=0');
        fprintf(t,'stimShape=1');
        fprintf(t,'preStimAmpSettle=200');
        fprintf(t,'postStimAmpSettle=500000');
        fprintf(t,'enableAmpSettle=1');
        fprintf(t,'SET');
        % Pause
        pause(0.5);%0.1
        FRAME = CHNlist(i,:);
        for j = 1:64
            if (strcmp(E_MAP{j+1},FRAME))
                break;
            end
        end
        if j>(32)
            j=j-32;
        end
        if strncmpi(CHNlist(i,:), 'A',1)
            fprintf(t,['changePort=' num2str(0)]);
        elseif strncmpi(CHNlist(i,:), 'B',1)
            fprintf(t,['changePort=' num2str(1)]);%Port A,B,C,D correspond with 0,1,2,3 while Analog are 4,5 Dig 6,7
        end
        fprintf(t,['changeFrame=' num2str(j-1)]);
        pause(0.5);%0.1
        fprintf(t,'OPENSTIM');
        pause(0.5);%0.1
        fprintf(t,'CLOSESTIM');
        pause(0.5);%0.1
        % Reset the TCP line
        OK = 0;
        while(~OK)
            OK = checkTCPOK(t);
            if OK ~= 1
                disp('WARNING');
                recoverTCP(t);
                pause(0.50)
            end
        end
        % Display the settings
        fprintf(t,'RUN');
        pause(1);
        DaqAOut(DAQ,0,1);
        pause(0.1);%0.01
        DaqAOut(DAQ,0,0);
        pause(1);
        fprintf(t,'STOP');
end
