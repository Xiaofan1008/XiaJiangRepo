%% Script that demonstrates Amp Settle on each channel in this experiment
CHNlist = cell2mat(unique(TRIAL.StimParams(2:end,1)));
% For each channel . . .
E_MAP = ProbeMAP;
E_MAP = E_MAP(:,3);
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
    pause(0.1);
    FRAME = CHNlist(i,:);
    for j = 1:32
        if (strcmp(E_MAP{j+1},FRAME))
            break;
        end        
    end
    fprintf(t,['changeFrame=' num2str(j-1)]);
    pause(0.1);
    fprintf(t,'OPENSTIM');
    pause(0.1);
    fprintf(t,'CLOSESTIM');
    pause(0.1);
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
    pause(0.01);
    DaqAOut(DAQ,0,0);
    pause(1);
    fprintf(t,'STOP');
end