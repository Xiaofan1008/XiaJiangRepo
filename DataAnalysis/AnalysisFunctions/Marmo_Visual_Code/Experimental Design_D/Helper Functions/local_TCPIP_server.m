function t = local_TCPIP_server
PORT = 9000;
disp(['The port is: ' num2str(PORT)]);
t = TcpLink(PORT);
data = fscanf(t);
data = string(data(1:end-1));
if ~strcmp(data,"READY")
    disp('ERROR. BAD CONNECTION.')
end
flushinput(t)
end