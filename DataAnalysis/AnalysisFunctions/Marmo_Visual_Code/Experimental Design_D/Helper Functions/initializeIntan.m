function initializeIntan(t,DAQ,nChn,port)
%% Function initializes Intan settings for experiment
while(1)
    E = input('Is Headstage Global Amp Settle enabled?\n','s');
    if (strcmpi(E,'Y'))
        D = input('Is the Digital Line connected and enabled?\n','s');
        if (strcmpi(D,'Y'))
            break;
        end
    end    
end

fprintf(t,'MCS');

for n = 1:nChn
    if (n-1) < 10
        chn = [port '-00' num2str(n-1)];
    else
        chn = [port '-0' num2str(n-1)];
    end
    fprintf(t,['channel=' chn]);
    fprintf(t,'enabled=1');
    fprintf(t,'postStimAmpSettle=500000');
    fprintf(t,'SET');
    pause(0.1);    
end

for n = 1:nChn
    fprintf(t,['changeFrame=' num2str(n-1)]);
    pause(0.1);    
    fprintf(t,'OPENSTIM');
    pause(0.1);
    fprintf(t,'CLOSESTIM');
    pause(0.1);
end

%% Check that amp settle is being applied
disp('Please observe Intan.');
disp('Amp Settle should be applied.');
input('Are you ready?\n','s');

%% Send a pulse
fprintf(t,'RUN');
pause(1);
DaqAOut(DAQ,0,1);
pause(0.01);
DaqAOut(DAQ,0,0);
pause(1);
fprintf(t,'STOP');

%% Clear the stimulator settings
OK = 0;
while(~OK)    
    OK = checkTCPOK(t);
    if OK ~= 1
        disp('WARNING');
        recoverTCP(t);
        pause(0.50)
    end
end
pause(0.01);
fprintf(t,'DALL');
pause(3);

%% Check that all channels are off
disp('Please observe Intan.');
disp('All channels should be off.');
input('Are you ready?\n','s');

fprintf(t,'RUN');
pause(1);
DaqAOut(DAQ,0,1);
pause(0.01);
DaqAOut(DAQ,0,0);
pause(1);
DaqAOut(DAQ,0,1);
pause(0.01);
DaqAOut(DAQ,0,0);
pause(1);
fprintf(t,'STOP');
pause(1)

fprintf(t,'NOMCS');

disp('Initialization complete.');

end