function stitchDataFiles
% This function requests two Intan datafiles, and stitches together the
% amplifier.dat and digitalin.dat files
nChn = 32;
FS = 30000;
T = 256; % Amount to read/write per cycle in seconds
% File Selection
[~, path1, ~] = ...
        uigetfile('*.rhs', 'Select the first RHS2000 Data File', 'MultiSelect', 'off');
[~, path2, ~] = ...
        uigetfile('*.rhs', 'Select the second RHS2000 Data File', 'MultiSelect', 'off');
% First File
cd(path1);
new_fid = fopen('amplifier_stitched.dat','w');
dig_fid = fopen('digitalin_stitched.dat','w');
f_fid = fopen('amplifier.dat','r')';
fd_fid = fopen('digitalin.dat','r');
while ~feof(f_fid)
    f = fread(f_fid,[nChn FS*T],'int16');
    fwrite(new_fid,f,'int16');    
end
fclose(f_fid);
nSam = dir('digitalin.dat');
nSam = nSam.bytes/2;
fd = fread(fd_fid, nSam, 'uint16');
fwrite(dig_fid,fd,'uint16');
fclose(fd_fid);
% Second File
cd(path2);
s_fid = fopen('amplifier.dat','r')';
sd_fid = fopen('digitalin.dat','r');
while ~feof(s_fid)
    s = fread(s_fid,[nChn FS*T],'int16');
    fwrite(new_fid,s,'int16');
end
nSam = dir('digitalin.dat');
nSam = nSam.bytes/2;
sd = fread(sd_fid, nSam, 'uint16');
fwrite(dig_fid,sd,'uint16');
fclose(sd_fid);
% Close Down
fclose(new_fid);
fclose(dig_fid);
% Rename files
movefile([path1 'amplifier.dat'],[path1 'amplifier_old.dat']);
movefile([path1 'amplifier_stitched.dat'],[path1 'amplifier.dat']);
movefile([path1 'digitalin.dat'],[path1 'digitalin_old.dat']);
movefile([path1 'digitalin_stitched.dat'],[path1 'digitalin.dat']);
end
