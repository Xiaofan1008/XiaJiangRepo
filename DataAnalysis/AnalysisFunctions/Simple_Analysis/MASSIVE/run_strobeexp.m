%% Script for running a stimulation experiment with Intan
function run_strobeexp
home = pwd; %returns the current directory in the string "home"

test = 0; % Test mode, disables communication with Intan and assumes everything will work. For duration assessments.
overalltime = input('How long would you like to run for (s)? \n');
flashnum = input('How many flashes per second? \n');
jitter = input('Add random 100ms jitter to timing (1=YES, 0=NO)? \n');
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
n_Trials = overalltime*flashnum;

if ~test
% Prompt the user for anything outstanding they may have forgotten
input('Have you set the Intan recording directory and file format?','s');
input('Are all communication protocols connected?','s');


% Prompt to begin
input('Are you ready to begin the experiment?\n','s');
fprintf(t,'NOMCS');
end

%% Run the experiment
if ~test
fprintf(t,'RECORD'); % Start the recording
end
tStart = tic;
% Set DAQ high
if ~test
    DaqAOut(DAQ,0,1);
end
pause(20); % Allow Intan to generate a small header
count = 0;

for n = 1:n_Trials
    count = count + 1;
    disp(['Running trial: ' num2str(n) ' of: ' num2str(n_Trials) '.']);

    if (count == 500)
        count = 0;
        fprintf('Pausing 5 seconds for indexing...\n');
        pause(5);
    end
    if ~test
        DaqAOut(DAQ,0,0);
    end
    pause(1/(flashnum*2));
    % Send a pulse
    if ~test
        DaqAOut(DAQ,0,1);
    end
    
    pause(1/(flashnum*2));
    TIME = toc(tStart);
    minutes = floor(TIME/60);
    seconds = floor(rem(TIME,60));
    ms = rem(TIME,1)*1000;
    fprintf('The experiment has been running for: %03.0f:%02.0f:%03.0f.\n',minutes,seconds,ms);
    if jitter==1
        jittime=0.1*rand(1,1);
        pause(jittime)
    end
end
DaqAOut(DAQ,0,1);
pause(20); % Allow Intan to generate a closing buffer

%% End of experiment
fprintf('All protocols completed successfully');

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

cd(home);
end