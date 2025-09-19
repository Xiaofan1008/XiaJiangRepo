%% Opens a local TCP/IP port for network-based communications
function t = local_TCPIP_server
CONNECTION = '0.0.0.0'; % This will only allow connections from this machine
% Use '0.0.0.0' to allow any connection, from anywhere
PORT = 9004; % Can be anything, but must be consistent. Don't use a common port.
TYPE = 'server'; % Use 'client' to connect to the other side of this connection.
%myIP=regexp(webread('http://checkip.dyndns.org','Timeout',30),'(\d+)(\.\d+){3}','match');
%disp(['The IP address is: ' myIP{1}]);
disp(['The port is: ' num2str(PORT)]);
t = tcpip(CONNECTION,PORT,'NetworkRole',TYPE);
fopen(t);
data = fscanf(t);
data = string(data(1:end-1));
if ~strcmp(data,"READY")
    disp('ERROR. BAD CONNECTION.')
end
flushinput(t)
end