home = pwd;
[DATAFILE,FOLDER] = uigetfile; % The experimental protocol file
cd(FOLDER);
exp_datafile = dir(DATAFILE);
t = local_TCPIP_server;

% Load in the experimental datafiles for quick access
TRIAL = load(exp_datafile.name);

timeoutSeconds = 10;
responseReceived = false;

%% Run the handshake
%write(t, uint8('set runmode record')); % Start the recording
n = 4;
write(t, uint8('execute clearalldataoutputs'));
%write(t, uint8('get type'), "char");
tic;
parameter_update(t,TRIAL.StimParams(n+1,:),n);

while toc < timeoutSeconds
    if t.NumBytesAvailable > 0
        elapsedTime = toc;
        responseReceived = true;
        break;
    end
end

if responseReceived
    fprintf('Response received in %.10f seconds.\n', elapsedTime);
else
    warning('Timeout: Expected response not received within %.1f seconds.', timeoutSeconds);
end

%write(t, uint8('set runmode stop')); % End the recording
