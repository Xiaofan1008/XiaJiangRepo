function allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue)
%% This function takes a blanked Intan datafile (amplifier_dn)
% and constructs mu
%% Variables
FS = 30000;
if strcmp(dName,'analogin')
    nChn = 1;
elseif strcmp(dName,'amplifier')
    nChn = 32;
end
threshfac = -3;
shortbytes=2;
sp = cell(1, nChn);
thresh = cell(1, nChn);
SCALEFACTOR = 10;
chk = 1; N = 0; time = 0;
NSp = zeros(1,nChn);
Ninput = 1e6;
ARTMAX = 0.5e3;
name = pwd;
name = strsplit(name,'\');
name = name{end};
name = name(1:end-14);
if isempty(dir('*.mu_sab3.dat'))
    justMu = 1;
else
    justMu = 0;
end
%% Load in the amplifier waveform data
try
    fileinfo = dir([dName '_dn_sab.dat']);
    nSam = fileinfo.bytes/(nChn * 2);
    v_fid = fopen([dName '_dn_sab.dat'], 'r');
catch
    try
        fileinfo = dir([dName '.dat']);
        nSam = fileinfo.bytes/(nChn * 2);
        v_fid = fopen([dName '.dat'], 'r');
    catch
        fprintf('No amplifier data found. Aborting . . .\n');
        return;
    end
end
lv_fid = fopen([dName '.dat'],'r');
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);
%% Generate filters
[Mufilt,Lfpfilt] = generate_Filters;
LfpNf = length(Lfpfilt);
MuNf = length(Mufilt);
%% Calculate thresholds
dispstat('','init');
dispstat(sprintf('Processing thresholds . . .'),'keepthis','n');
for iChn = 1:nChn
    sp{iChn} = zeros(ceil(ntimes * double(T) * 20), FS * 1.6 / 1e3 + 1 + 1);
end
artchk = zeros(1,nChn);
munoise = cell(1,nChn);
nRun = 0;
tRun = ceil(nSam / FS / 10) + 1; % Number of 10 second chunks of time in the data
while sum(artchk) < nChn
    if strcmp(dName,'analogin')        
        v = fread(v_fid, [nChn, FS*20], 'uint16');
        v = (v - 32768) .* 0.0003125;
    elseif strcmp(dName,'amplifier')  
        v = fread(v_fid, [nChn, FS*20], 'int16') * 0.195; % reading 20s 
        if size(v,2)<trig(1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            v = fread(v_fid, [nChn, FS*20], 'int16') * 0.195; % reading 20s 
            trig=trig-FS*20;
        end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    nRun = nRun + 1;
    dispstat(sprintf('Progress %03.2f%%',(100*(nRun/tRun))),'timestamp');
    if ~isempty(v)
        % Checks for artifact
        for iChn = 1:nChn
            if isempty(thresh{iChn})
                munoise{iChn} = [];
                mu = conv(v(iChn,:),Mufilt); %why conv not flipped data
                mu = detrend(mu(trig(3)+275:trig(3)+100*FS/1000));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if max(abs(mu(MuNf+1:end-MuNf))) < ARTMAX
                    artchk(iChn) = 1;
                    sd = median(abs(mu(MuNf+1:end-MuNf)))./0.6745;
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
                sd = median(abs(munoise{iChn}))./0.6745;
                thresh{iChn} = threshfac*sd;
            end
        end
    end
end
disp(['Total recording time: ' num2str(nSam / FS) ' seconds.']);
disp(['Time analysed per loop: ' num2str(T) ' seconds.']);
fseek(v_fid,0,'bof'); % Returns the data pointer to the beginning of the file
%% Loop through the data
if (justMu)
    mu_fid = fopen([name '.mu_sab3.dat'],'W');
    dispstat('','init');
    dispstat(sprintf('Processing MU . . .'),'keepthis','n');
    while (chk && N < Ninput)
        N = N + 1;        
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/ntimes)),'timestamp');
        if strcmp(dName,'analogin')
            data = fread(v_fid, [nChn, FS*T], 'uint16');
            data = (data - 32768) .* 0.0003125;
        elseif strcmp(dName,'amplifier')
            data = fread(v_fid, [nChn, FS*T], 'int16') * 0.195;
        end
        if (size(data,2))
            Ndata = size(data,2);
            MuOut = cell(1,nChn);
            if (par)
                parfor iChn = 1:nChn
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data,Mufilt);
                    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
                    MuOut{iChn} = mu;
                end
            else
                for iChn = 1:nChn
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data,Mufilt);
                    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
                    MuOut{iChn} = mu;
                end
            end
            mu2 = zeros(nChn,Ndata);
            for iChn = 1:nChn
                mu2(iChn,:) = MuOut{iChn};
            end
            fwrite(mu_fid,SCALEFACTOR*mu2,'short');
            time = time + T*1e3;
        end
        if (size(data,2) < FS * T)
            chk = 0;
        end
        dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
    end
    fclose(mu_fid);
    clear data mu mu2 tmp flip_data
end
%% fixing amp issue
% if (amp_issue==1)&& isempty(dir('*.muamp_sab.mat'))
%     m_fid = fopen([name '.mu_sab3.dat'],'r');
%     dispstat('','init');
%     dispstat(sprintf('Amplifying . . .'),'keepthis','n');
%     chk = 1; N = 0; time = 0;
%     while (chk && N < Ninput)
%          N = N + 1;
%         dispstat(sprintf('Progress %03.2f%%',100*((N-1)/ntimes)),'timestamp');
%         mu = fread(m_fid,[nChn, FS*T],'short') ./ 10;
%         if (size(mu,2))
%             tnum=1;
%             shortbytes=2;
%             offset_test=(trig(tnum)+2)*shortbytes;%offset from beginning of file to trigger
%             while offset_test<size(mu,2)
%                 ftell(m_fid)
%                 offset=(trig(tnum)+2)*nChn*shortbytes;%offset from beginning of file to trigger
%                 fseek(m_fid,offset,'bof');
%                 ftell(m_fid)
%                 offsetmu = fread(m_fid,[nChn, 0.1*FS],'short')./10; %plots o5ne second from trigger
%                 sd = median(abs(offsetmu(iChn,:)))./0.6745;
%                 if thresh{iChn}<(sd*threshfac-2)
%                     fseek(m_fid,1,'bof');
%                     ftell(m_fid)
%                     mu_temp=mu;
%                     mu_temp(:,trig(tnum):trig(tnum)+100*FS/1000)=detrend(mu(:,trig(tnum):trig(tnum)+100*FS/1000).*thresh{iChn}/(sd*threshfac));
%                     
%                         figure
%                         plot (mu_temp(1,trig(tnum):trig(tnum)+30000))
%                         title(['Channel ' num2str(1)])
%                         xlabel('Time (ms)')
%                         ylabel('Voltage (uV)')
%                         
%                 end
%                 tnum=tnum+1;
%                 offset_test=(trig(tnum)+2)*shortbytes;%offset from beginning of file to trigger
%             end
%         end
%     end
%     
% end

%%detrend after trigger for 100ms
trig = loadTrig(0);
theseTrig = trig;
lengthTrig=length(trig);
%% Calculate SP
if isempty(dir('*.sp.mat'))
    m_fid = fopen([name '.mu_sab3.dat'],'r');
    dispstat('','init');
    dispstat(sprintf('Processing SP . . .'),'keepthis','n');
    chk = 1; N = 0; time = 0;
    TrialParams=loadTrialParams;
%     for tID=1:5
%         TrialParamstID = find(cell2mat(TrialParams(1:end,2)) == tID); %identifies trial row matching trial ID
%         TrialParamstID(1:2:end)=[];
%         TrialParamstID=TrialParamstID./2;
%         tmptrig=trig;
%         trig(TrialParamstID)=0;
%     end
    while (chk && N < Ninput)
        N = N + 1;
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/ntimes)),'timestamp');
        mu = fread(m_fid,[nChn, FS*T],'short') ./ 10;
        if (size(mu,2))
            if (par)
                parfor iChn = 1:nChn
                    [Sp_tmp, Spktimes_tmp] = spikeextract(mu(iChn,:), thresh{iChn},FS,artefact,artefact_high);
                    NSp_tmp = length(Spktimes_tmp);
                    if NSp_tmp > 1
                        SpMat = [Spktimes_tmp+double(time) Sp_tmp];
                        NSp_tmp = size(SpMat,1);
                        sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
                        NSp(iChn) = NSp(iChn) + NSp_tmp;
                    
                    end
                end
            else
                for iChn = 1:nChn
                    for tnum=1:length(trig)
                        %if trig(tnum)==0

                        if(trig(tnum)-(N-1)*T*FS)<size(mu,2)
                            if ((trig(tnum)+100*FS/1000-(N-1)*T*FS)>size(mu,2))&&((size(mu,2)-(trig(tnum)+100*FS/1000-(N-1)*T*FS))<1000)
                                mu_temp=(mu(:,(trig(tnum)-(N-1)*T*FS):end));
                            else
                                mu_temp=(mu(:,(trig(tnum)-(N-1)*T*FS):trig(tnum)+100*FS/1000-(N-1)*T*FS));
                            end
                            sd = median(abs(mu_temp(iChn,MuNf+1+50:end-MuNf-50)))./0.6745;
                            if sd>5
                                thresh = -4*sd;
                            else
                                thresh = -4*sd;
                            end
                            [Sp_tmp, Spktimes_tmp] = spikeextract(mu_temp(iChn,:), thresh, FS,artefact,artefact_high);
                            if (~isempty(Sp_tmp) && any(Sp_tmp(:)<(-200)))
                                stop=0;
                            end
                            NSp_tmp = length(Spktimes_tmp);
                            if NSp_tmp > 1
                                %SpMat = [Spktimes_tmp+double(time) Sp_tmp];
                                SpMat = [Spktimes_tmp+trig(tnum)*1000/FS Sp_tmp];
                                NSp_tmp = size(SpMat,1);
                                sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
                                NSp(iChn) = NSp(iChn) + NSp_tmp;
                            end
                            num_trig=tnum;
                        end
                    end
                end
                if (size(trig,2) > 0)
                    %check=lengthTrig-length(trig);
                    trig(1:num_trig)=[];
                end
            end
            %time = time + T*1e3;
        end
        if (size(mu,2) < FS * T)
            chk = 0;
        end
        dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
    end
    
    %% Calculate a zero-condition spike template and r2t
    if ~isempty(dir('*_exp_datafile_*'))
        trig = loadTrig(0); loadNREP; loadNTrials; loadDUALSTIM;
        TrialParams = loadTrialParams;
        TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == 1));
        if DUALSTIM
            TrialParams(2:2:length(TrialParams),:) = [];
        end
        TrialParams = TrialParams(1:n_REP);
        trig = trig(TrialParams);
        nTrig = length(trig);
        spWN = [-300 300]; template = cell(1,nChn); r2t = cell(1,nChn); spt = cell(1,nChn);
        for c = 1:nChn
            for n = 1:nTrig
                spt{c} = [spt{c}; sp{c}(sp{c}(:,1) > trig(n)/30+spWN(1) & sp{c}(:,1) < trig(n)/30+spWN(2),2:end)];
            end
            template{c} = mean(cell2mat(spt(c)),1);
            r2t{c} = mean((cell2mat(spt(c)) - template{c}).^2);
        end
        %% Save the datafiles
        disp('Saving spikes');
        for iChn = 1:nChn
            sp{iChn} = sp{iChn}(1:NSp(iChn),:);
        end
        save([name '.sp.mat'],'sp','thresh','threshfac','template','r2t','-v7.3');
    else
        disp('Saving spikes');
        for iChn = 1:nChn
            sp{iChn} = sp{iChn}(1:NSp(iChn),:);
        end
        save([name '.sp.mat'],'sp','thresh','threshfac','-v7.3');
    end
end