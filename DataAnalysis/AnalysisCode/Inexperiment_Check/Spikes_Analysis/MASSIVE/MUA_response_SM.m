function MUA_response_SM
%Used to save a MUA response file
par=0;
nChn=32;
FS=30000;
fileinfo = dir('*.mu_sab.dat');
filename= fileinfo.name;
v_fid = fopen(filename, 'r');
mem = memory;
T = mem.MaxPossibleArrayBytes ./ (2 * 32 * 30000);
t_len = fileinfo.bytes/(32 * 2 * 30000);
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
end
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);
ARTMAX = 0.5e3;
threshfac = 2;
Ninput = 1e6;
[Mufilt,Lfpfilt] = generate_Filters;
thresh = cell(1, nChn);
artchk = zeros(1,nChn);
munoise = cell(1,nChn);
LfpNf = length(Lfpfilt);
%% 1. Clipping extreme values +-2SD   2. RMS   3.Low pass filter
if isempty(dir('MUAdata.dat'))

LP_fid = fopen('MUAdata.dat','W');

    dispstat('','init');
    dispstat(sprintf('Processing SD threshold . . .'),'keepthis','n');
    chk = 1; N = 1; 
    while sum(artchk) < nChn
        N = N + 1;
        mu = fread(v_fid,[nChn, FS*15],'short') ./ 10;
        
        if ~isempty(mu)
            % Checks for artifact
            for iChn = 1:nChn
                if isempty(thresh{iChn})
                    munoise{iChn} = [];
                    if (max(abs(mu(1,:))) < ARTMAX)||(max(abs(mu(15,:))) < ARTMAX)
                        artchk(iChn) = 1;
                        sd = median(abs(mu(iChn)))./0.6745;
                        thresh{iChn} = threshfac*sd;
                    else
                        munoise{iChn} = [munoise{iChn} mu(MuNf+1:end-MuNf)];
                    end
                end
            end
        else
            for iChn = 1:nChn
                if isempty(thresh{iChn}) % If threshold is still 0 - noisy channel
                    artchk(iChn) = 1;
                    sd = median(abs(mu(iChn)))./0.6745;
                    thresh{iChn} = threshfac*sd;
                end
            end
        end
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/(T*ntimes/20))),'timestamp');
    end
    fclose(v_fid);
    dispstat('','init');
    dispstat(sprintf('Processing MUA response . . .'),'keepthis','n');
    N=0;
    v_fid = fopen(filename, 'r');
    while (chk && N < Ninput)
        N = N + 1;
        data = fread(v_fid,[nChn, FS*T],'short') ./ 10;
        if (size(data,2))
            Ndata = size(data,2);
            MuOut = cell(1,nChn);
            if (par==1)
                parfor iChn = 1:nChn
                    data(iChn,:)=min(data(iChn,:),-thresh{iChn});
                    data(iChn,:)=max(data(iChn,:),thresh{iChn});
                    data(iChn,:)=sqrt(data(iChn,:).^2);
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data ,Lfpfilt);
                    LP_data = fliplr(tmp(1,LfpNf/2:Ndata+LfpNf/2-1));
                    LPOut{iChn} = LP_data;
                end
            else
                for iChn = 1:nChn
                    data(iChn,:)=min(data(iChn,:),thresh{iChn});
                    data(iChn,:)=max(data(iChn,:),-thresh{iChn});
                    data(iChn,:)=sqrt(data(iChn,:).^2);
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data ,Lfpfilt);
                    LP_data = fliplr(tmp(1,LfpNf/2:Ndata+LfpNf/2-1));
                    LPOut{iChn} = LP_data;
                end
            end
            LP2 = zeros(nChn,Ndata);
            for iChn = 1:nChn
                LP2(iChn,:) = LPOut{iChn};
            end
            fwrite(LP_fid,10*LP2,'short');            
            if (size(data,2) < FS * T)
                chk = 0;
            end
            dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
        end
    end
    fclose(LP_fid);
    
end







%% 4.downsample?

end

