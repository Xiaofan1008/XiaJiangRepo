function LFPextract
%% Extracts a filtered lfp.dat file from an amplifier_lfp.dat file
%% Variables
FS = 30000;
nChn = 32;
T = 63;                         % Number of seconds per data loop
chk = 1; N = 0; time = 0;
Ninput = 1e6;
name = pwd;
name = strsplit(name,'\');
name = name{end};
name = name(1:end-14);
%% Setup
lfp_fid = fopen([name '.lfp2.dat'],'w');
%% Load in the amplifier waveform data
try
    fileinfo = dir('amplifier_lfp.dat');
    v_fid = fopen('amplifier_lfp.dat', 'r');
catch
    return;
end
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);
%% Generate filters
[~,Lfpfilt] = generate_Filters;
LfpNf = length(Lfpfilt);
%% Filtering
while (chk && N < Ninput)
    tic;
    N = N + 1;
    disp(['Processing LFP: Loop ' num2str(N) ' of ' num2str(ntimes)]);
    data = fread(v_fid,[nChn,FS*T],'int16') .* 0.195;
    if (size(data,2))
        Ndata = size(data,2);
        LfpOut = cell(1,nChn);        
        for iChn = 1:nChn                        
            tmp = conv(data(iChn,:),Lfpfilt);
            lfp = tmp(1,LfpNf/2:FS/1e3:Ndata+LfpNf/2-1);            
            LfpOut{iChn} = lfp;
        end
        lfp2 = zeros(nChn,size(LfpOut{1},2));
        for iChn = 1:nChn
            lfp2(iChn,:) = LfpOut{iChn};
        end
        fwrite(lfp_fid,lfp2,'float');
        time = time + T*1e3;
    end
    if (size(data,2) < FS * T)
        chk = 0; 
    end
    Elapsed = toc; disp([num2str(Elapsed) 's to process ' num2str(T) 's']);
end
fclose(lfp_fid);
%% Closing

end