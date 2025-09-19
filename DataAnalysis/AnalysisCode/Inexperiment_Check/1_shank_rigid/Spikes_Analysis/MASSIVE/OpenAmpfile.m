
nChn=32;
FS=30000;%sampling rate Hz
T=5; %amount of time you want to read in seconds
filepath = pwd;


vFID = fopen([filepath filesep 'amplifier.dat'],'r');
%%%% if you want to offset from beginning of file uncomment these

% offset_seconds=5; %time from beginning of file to offset
% offset=offset_seconds*FS*nChn*shortbytes;%offset from beginning of file
% ftell(vFID)
% fseek(vFID,offset,'bof');
% ftell(vFID)

%%%

v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195; % reads data 32 channels

figure
plot(v(1,:)) %plots channel one(based on headstage) for 5 seconds