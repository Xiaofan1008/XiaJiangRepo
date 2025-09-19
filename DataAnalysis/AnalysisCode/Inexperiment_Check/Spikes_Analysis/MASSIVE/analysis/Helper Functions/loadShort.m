here = pwd;
file = [here '\short.dat'];
FID = fopen(file,'r');
v = fread(FID,[32, 1800000],'int16') .* 0.195;
clear here file FID