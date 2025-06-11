function extractLFP(nChn,FS)
%% This function attempts to low-pass filter the raw data
% You must be in the right directory, and have a raw amplifier data file present
if nargin<2
    FS = 30000;
end
if nargin<1
    nChn = 32;
end
name = strsplit(pwd,'\');
name = name{end};
name = name(1:end-14);
Ninput = 1e6;
lv_fid = fopen('amplifier.dat','r');
lfp_fid = fopen([name '.lfp.dat'],'W');
%% Generate filters
[~,Lfpfilt] = generate_Filters;
LfpNf = length(Lfpfilt);
%% Calculate LFP
chk = 1; N = 0;
while (chk && N < Ninput)
    N = N + 1;
    datalfp = fread(lv_fid, [nChn, FS*256], 'int16') * 0.195;
    if (size(datalfp,2))
        Ndata = size(datalfp,2);
        LfpOut = cell(1,nChn);
        for iChn = 1:nChn
            tmp = conv(datalfp(iChn,:),Lfpfilt);
            lfp = tmp(1,LfpNf/2:FS/1e3:Ndata+LfpNf/2-1);
            LfpOut{iChn} = lfp;
        end
        lfp2 = zeros(nChn,size(LfpOut{1},2));
        for iChn = 1:nChn
            lfp2(iChn,:) = LfpOut{iChn};
        end
        fwrite(lfp_fid,lfp2,'float');
    end
    if (size(datalfp,2) < FS * 256)
        chk = 0;
    end
end
fclose(lfp_fid);
end

