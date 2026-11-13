function lesion(chn)
mu = char(hex2dec('03BC'));
test = 0;
%% Selecting a channel
if nargin < 1
    chn = input('Please enter a channel to lesion\n');
end
fprintf('Performing a lesion on channel %1.0f\n',chn);
%% Initial settings
amp = 5;
dur = 5;
%% Verifying settings
fprintf(['Lesioning will use %1.0f ' mu 'A of current delivered monophasically for %1.0f seconds.\n'],amp,dur);
chk = input('Are these settings appropriate? Please enter "1" for YES.');
while (~chk)
    amp = input(['Please enter a lesion amplitude in ' mu 'A.\n']);
    dur = input('Please enter a lesion duration in seconds.\n');    
    fprintf(['Lesioning will use %1.0f ' mu 'A of current delivered monophasically for %1.0f seconds.\n'],amp,dur);
    chk = input('Are these settings appropriate? Please enter "1" for YES.');
end
%% Set up for lesioning
% TCP/IP connection
if ~test
    fprintf('Please connect Intan to the TCP/IP server\n');
    t = local_TCPIP_server;
end
% DAQ connection
if ~test
    DAQ = DaqDeviceIndex;
    DaqAOut(DAQ,0,0);
end
% Electrode mapping
E_MAP = ProbeMAP;
E_MAP = E_MAP(:,3); % USE THIS FIELD TO ADJUST THE ELECTRODE MAPPING
if ~test
    input('Have you verified the Intan filename and directory?');
    % Prompt to begin
    input('Are you ready to begin the experiment?','s');
    fprintf(t,'NOMCS');
end
%% Applying settings
lesionStimSet = newSettings;
lesionStimSet{1} = E_MAP{chn+1,1};      % Channel
lesionStimSet{4} = 1;                   % Trigger set to LEVEL
lesionStimSet{7} = 0;                   % Single pulse or train
lesionStimSet{10} = 1000;               % Post-stim refractory period
lesionStimSet{11} = 0;                  % Stimulation shape
lesionStimSet{13} = 5e3;                % First phase duration
lesionStimSet{14} = 0;                  % Second phase duration
lesionStimSet{16} = amp;                % First phase amplitude
lesionStimSet{17} = 0;                  % Second phase amplitude
lesionStimSet{18} = 0;                  % Disable amp settle
lesionStimSet{19} = 0;                  % Amp settle pre-stim
lesionStimSet{20} = 0;                  % Amp settle post-stim
lesionStimSet{22} = 0;                  % Disable charge recovery
lesionStimSet{23} = 0;                  % Charge recovery pre-stim
lesionStimSet{24} = 0;                  % Charge recovery post-stim
%% Run the experiment
if ~test
    DaqAOut(DAQ,0,0);
    fprintf(t,'RECORD'); % Start the recording
end
tStart = tic;
pause(5); % Allow Intan to generate a small header
OK = 0;
while(~OK)
    if ~test
        parameter_update(t,lesionStimSet,1);
        % Send the parameter update
        OK = checkTCPOK(t);
        if OK ~= 1
            disp('WARNING');
            recoverTCP(t);
            pause(0.01)
        end
    else
        OK = 1;
    end
end
if ~test
    fprintf(t,'SET');
    % Verify successful send of parameters
    S = verifyParameters(t,1);
    while ~(S)
        fprintf('WARNING: COMMUNICATIONS INTERRUPTED. RESENDING PARAMETERS\n');
        % Resend the parameters
        while(~OK)
            if ~test
                parameter_update(t,lesionStimSet,1);
                % Send the parameter update
                OK = checkTCPOK(t);
                if OK ~= 1
                    disp('WARNING');
                    recoverTCP(t);
                    pause(0.01)
                end
            else
                OK = 1;
            end
        end
        fprintf(t,'SET');
        % Verify successful send of parameters
        S = verifyParameters(t,n);
    end
else
    disp('Waiting 0.2');
    pause(0.2);
end
% Send a pulse
chk = input('Last chance - do you want to lesion here? Enter "1" for YES.');
if ~test && chk
    DaqAOut(DAQ,0,1);
    pause(dur);
    DaqAOut(DAQ,0,0);
end
pause(5);
TIME = toc(tStart);
minutes = floor(TIME/60);
seconds = floor(rem(TIME,60));
ms = rem(TIME,1)*1000;
fprintf('The lesion took: %03.0f:%02.0f:%03.0f.\n',minutes,seconds,ms);
%% End of experiment
fprintf('All protocols completed successfully\n');
if ~test
    fprintf(t,'STOP'); % End the recording
    % Make sure all channels are switched off
    fprintf(t,'DALL');
end
% Close the TCP port
if ~test
    fclose(t);
    clear t;
end
end