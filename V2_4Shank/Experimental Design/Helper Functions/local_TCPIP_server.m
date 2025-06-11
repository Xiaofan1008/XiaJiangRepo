%% Opens a local TCP/IP port for network-based communications
function t = local_TCPIP_server
CONNECTION = '49.127.50.247'; % This will only allow connections from this machine
% Use '0.0.0.0' to allow any connection, from anywhere
PORT = 9001; % Can be anything, but must be consistent. Don't use a common port.
% PORT = 5001; % Can be anything, but must be consistent. Don't use a common port.
TYPE = 'server'; % Use 'client' to connect to the other side of this connection.
% [~,myIP]=system('ifconfig'); 
% myIP = strsplit(myIP,'IPv4'); 
% % myIP = myIP{3}(33:43);
% % myIP = '49.127.50.247';
% disp(['The IP address is: ' myIP]);
% disp(['The port is: ' num2str(PORT)]);
t = tcpip(CONNECTION,PORT,'NetworkRole',TYPE);
% t = tcpip(myIP,PORT,'NetworkRole',TYPE);
fopen(t);
data = fscanf(t);
data = string(data(1:end-1));
if ~strcmp(data,"READY")
    disp('ERROR. BAD CONNECTION.')
end
flushinput(t)
end