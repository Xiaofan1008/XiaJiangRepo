function shortIntan(dir)
%% This function takes an experimental Intan directory and creates a 60 second Intan datafile from the amplifier.dat file
% This is to be used for quickly testing blanking functions and spike
% extraction
%% Variables
len = 60 * 1e3; % Length in seconds
FS = 30000;     % Sampling rate
nChn = 32;      % Number of channels
%% Logic
iFile = [dir '\amplifier.dat'];
kFile = [dir '\short.dat'];
iFID = fopen(iFile,'r');
kFID = fopen(kFile,'w');
v = fread(iFID,[nChn,(len * (FS/1e3))],'int16');
fwrite(kFID,v,'int16');
end