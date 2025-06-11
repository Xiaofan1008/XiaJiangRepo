function extractMUA(nChn,FS)
%% This function attempts to high-pass filter the blanked data
% You must be in the right directory, and have a blanked amplifier data file or raw amplifier data file present
if nargin<2
    FS = 30000;
end
if nargin<1
    nChn = 32;
end
SCALEFACTOR = 10;
chk = 1; N = 0;
Ninput = 1e6;
name = strsplit(pwd,'\');
name = name{end};
name = name(1:end-14);
%% Load in the amplifier waveform data
try   
    v_fid = fopen('amplifier_dn.dat', 'r');
catch
    try
        v_fid = fopen('amplifier.dat', 'r');
    catch
        error('No amplifier data found\n');
    end
end
%% Generate filters
[Mufilt,~] = generate_Filters;
MuNf = length(Mufilt);
%% Loop through the data
mu_fid = fopen([name '.mu.dat'],'W');
while (chk && N < Ninput)
    N = N + 1;
    data = fread(v_fid, [nChn, FS*256], 'int16') * 0.195;
    if (size(data,2))
        Ndata = size(data,2);
        MuOut = cell(1,nChn);        
        for iChn = 1:nChn
            flip_data = fliplr(data(iChn,:));
            tmp = conv(flip_data,Mufilt);
            mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
            MuOut{iChn} = mu;
        end
        mu2 = zeros(nChn,Ndata);
        for iChn = 1:nChn
            mu2(iChn,:) = MuOut{iChn};
        end        
        fwrite(mu_fid,SCALEFACTOR*mu2,'short');
    end
    if (size(data,2) < FS * 256)
        chk = 0;
    end
end
fclose(mu_fid);
end

