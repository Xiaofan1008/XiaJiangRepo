%% Script for running a stimulation experiment with Intan
function run_exp
home = pwd; %returns the current directory in the string "home"
[DATAFILE,FOLDER] = uigetfile; % The experimental protocol file
cd(FOLDER);
test = 1; % Test mode, disables communication with Intan and assumes everything will work. For duration assessments.
%% Local variables
exp_datafile = dir(DATAFILE);

disp(['Now running experimental datafile "' DATAFILE '"']);

%% Set up the experiment
% First, we need to create a local TCPIP server for communication
if ~test
disp('Please connect Intan to the TCP/IP server');
t = local_TCPIP_server;
end

%% Grab the DAQ
if ~test
DAQ = DaqDeviceIndex;
end

% Load in the experimental datafiles for quick accessdes
TRIAL = load(exp_datafile.name);
n_Trials = TRIAL.n_Trials;
param_to_pulse = TRIAL.PARAM; % Time (msec) between parameter update packet and stimulus pulse
pulse_to_param = TRIAL.PULSE; % Time (msec) between stimulation pulse and parameter update packet

% Create timer objects
t_pulse_param = timer('TimerFcn','disp(''PARAMETER UPDATE'')','StartDelay',(pulse_to_param/1000));

if ~test
% Prompt the user for anything outstanding they may have forgotten
input('Have you set the Intan recording directory and file format?','s');
input('Have you loaded the correct experimental datafile?','s');
input('Are all communication protocols connected?','s');
input('Have you checked the stimulus timing inputs?','s');

% Prompt to begin
input('Are you ready to begin the experiment?\n','s');
fprintf(t,'NOMCS');
end

exp_len = TRIAL.n_Trials*(TRIAL.PARAM+TRIAL.PULSE+125);
minutes = floor((exp_len)/60000);
seconds = rem((exp_len),60000)/1000;

if ~test
disp('We will now demonstrate the Amp Settle function on each channel to be used in this experiment.');
disp('Please monitor Intan. If you do not see a large yellow window after stimulation, abort and investigate.');

dispAmpSettle;
end
DUALSTIM= input('Please enter whether you would like to stimulate on two electrodes YES=1 NO=0: \n');

input('Have you verified the charge recovery settings?');

fprintf('The experiment should run for a little more than: %04.0f minutes and %02.0f seconds.\n',(minutes),(seconds));

%% Run the experiment
if ~test
fprintf(t,'RECORD'); % Start the recording
end
tStart = tic;
pause(20); % Allow Intan to generate a small header
count = 0;
for n = 1:n_Trials
    count = count + 1;
    OK = 0;
    disp(['Running trial: ' num2str(n) ' of: ' num2str(n_Trials) '.']);
    % Start a timer
    start(t_pulse_param);
    % Initialise the parameter update
    disp(['Waiting ' num2str(pulse_to_param/1000)]);
    wait(t_pulse_param);
    % Reset the line for the next pulse
    if ~test
    DaqAOut(DAQ,0,0);
    end
    if (count == 500)
        count = 0;
        fprintf('Pausing 5 seconds for indexing...\n');
        pause(5);
    end
    while(~OK)
        if ~test
        parameter_update(t,TRIAL.StimParams(n+1,:),n);      
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
        S = verifyParameters(t,n);
        while ~(S)
            fprintf('WARNING: COMMUNICATIONS INTERRUPTED. RESENDING PARAMETERS\n');
            % Resend the parameters
            while(~OK)
                if ~test
                    parameter_update(t,TRIAL.StimParams(n+1,:),n);
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
        
        
        if DUALSTIM==1
            if TRIAL.StimParams{n+2,1}~=1
            while(~OK)
                if ~test
                    parameter_update(t,TRIAL.StimParams(n+2,:),n);
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
                S = verifyParameters(t,n);
                while ~(S)
                    fprintf('WARNING: COMMUNICATIONS INTERRUPTED. RESENDING PARAMETERS\n');
                    % Resend the parameters
                    while(~OK)
                        if ~test
                            parameter_update(t,TRIAL.StimParams(n2,:),n);
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
            
            end
            end
            n=n+1;
        end
        
    end   
    disp(['Sent stimulus condition ID: ' num2str(TRIAL.TrialParams{n+1,2})]);
    % Jitter the stimulus to the 4 Hz band, with a minimum duration
    r = randi([param_to_pulse 125+param_to_pulse])/1e3;
    disp(['Waiting ' num2str(r)]);
    pause(r);    
    % Send a pulse
    if ~test
    DaqAOut(DAQ,0,1);
    end
    exp_len = exp_len-TRIAL.PARAM-TRIAL.PULSE-125;
    minutes = floor((exp_len-TRIAL.PARAM-TRIAL.PULSE)/60000);
    seconds = rem((exp_len-TRIAL.PARAM-TRIAL.PULSE),60000)/1000;
    if ~(n == n_Trials)
        fprintf('The experiment should last a little less than: %04.0f minutes and %02.0f seconds.\n',(minutes),(seconds));
    end
    TIME = toc(tStart);
    minutes = floor(TIME/60);
    seconds = floor(rem(TIME,60));
    ms = rem(TIME,1)*1000;
    fprintf('The experiment has been running for: %03.0f:%02.0f:%03.0f.\n',minutes,seconds,ms);
end
DaqAOut(DAQ,0,0);
pause(20); % Allow Intan to generate a closing buffer

%% End of experiment
fprintf('All protocols completed successfully');

if ~test
fprintf(t,'STOP'); % End the recording

% Make sure all channels are switched off
fprintf(t,'DALL');    
end
delete(t_pulse_param);

% Close the TCP port
if ~test
fclose(t);
clear t;
end

cd(home);
end