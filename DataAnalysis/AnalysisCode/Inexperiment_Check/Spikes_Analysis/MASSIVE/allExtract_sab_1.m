function allExtract_sab_1(dName,filepath,T,par,artefact,artefact_high)
%% This function takes a blanked Intan datafile (amplifier_dn)
% and constructs mu
%% Variables
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
FS=frequency_parameters.amplifier_sample_rate;
if strcmp(dName,'analogin')
    nChn = 1;
elseif strcmp(dName,'amplifier')
    nChn=size(amplifier_channels,2);  
end
FS=frequency_parameters.amplifier_sample_rate;
threshfac = -4.5;
sp = cell(1, nChn);
thresh = cell(1, nChn);
SCALEFACTOR = 10;
chk = 1; N = 0; time = 0;
NSp = zeros(1,nChn);
Ninput = 1e6;
ARTMAX = 0.5e3; % artifact threshold 
name = filepath;
name = strsplit(name, filesep);
name = name{end};
name = name(1:end-14);
if isempty(dir('*.mu_sab.dat'))
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
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T); % number of blocks
%% Generate filters
[Mufilt,Lfpfilt] = generate_Filters;
LfpNf = length(Lfpfilt);
MuNf = length(Mufilt);
%% Calculate thresholds
%if isempty(dir('*.sp.mat'))
    dispstat('','init');
    dispstat(sprintf('Processing thresholds . . .'),'keepthis','n');
    for iChn = 1:nChn
        sp{iChn} = zeros(ceil(ntimes * double(T) * 20), FS * 1.6 / 1e3 + 1 + 1); 
        % Rows: Numebr of spike events; Column: waveform length (number of samples in 1.6 ms)
    end
    artchk = zeros(1,nChn); % flages to indicate threshold calculated 
    munoise = cell(1,nChn);  % accumulate noise for noisy chunks 
    nRun = 0;
    tRun = ceil(nSam / FS / 10) + 1; % Number of 10 second chunks of time in the data
    while sum(artchk) < nChn
        if strcmp(dName,'analogin')
            v = fread(v_fid, [nChn, FS*20], 'uint16'); % 20s of data
            v = (v - 32768) .* 0.0003125;
        elseif strcmp(dName,'amplifier')
            %v = vblank(1:nChn, 1:FS*20);
            v = fread(v_fid, [nChn, FS*20], 'int16') * 0.195; % reading 20s
        end
        nRun = nRun + 1;
        dispstat(sprintf('Progress %03.2f%%',(100*(nRun/tRun))),'timestamp');
        if ~isempty(v)
            % Checks for artifact
            for iChn = 1:nChn
                if isempty(thresh{iChn})
                    munoise{iChn} = [];
                    mu = conv(v(iChn,:),Mufilt); %why conv not flipped data
                    if max(abs(mu(MuNf+1:end-MuNf))) < ARTMAX % (MuNf+1:end-MuNf): ignore the convolution edges (MuNf is the half filter lenght)
                        artchk(iChn) = 1;
                        sd = median(abs(mu(MuNf+1:end-MuNf)))./0.6745; % estimate the standard deviation using median absolutie deviation(MAD) sd ~median(|x|)/0.6745
                        thresh{iChn} = threshfac*sd;
                    else
                        munoise{iChn} = [munoise{iChn} mu(MuNf+1:end-MuNf)]; % if noisy, save the filtered data to try again later
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
%end
%% Loop through the data
if (justMu)
    mu_fid = fopen([name '.mu_sab.dat'],'W');
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
                parfor iChn = 1:nChn % zero-phase filtering
                    flip_data = fliplr(data(iChn,:));
                    tmp = conv(flip_data,Mufilt);
                    mu = fliplr(tmp(1,MuNf/2:Ndata+MuNf/2-1));
                    MuOut{iChn} = mu;
                end
            else
                for iChn = 1:nChn  % zero-phase filtering
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
            dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
        end
        if (size(data,2) < FS * T)
            chk = 0;
        end
    end
    fclose(mu_fid);
    clear data mu mu2 tmp flip_data
end
%% Calculate SP
if isempty(dir('*.sp.mat'))
    m_fid = fopen([name '.mu_sab.dat'],'r');
    dispstat('','init');
    dispstat(sprintf('Processing SP . . .'),'keepthis','n');
    chk = 1; N = 0; time = 0;
    while (chk && N < Ninput)
        N = N + 1;
        dispstat(sprintf('Progress %03.2f%%',100*((N-1)/ntimes)),'timestamp');
        mu = fread(m_fid,[nChn, FS*T],'short') ./ 10; % read data and scale down to real µV values
        if (size(mu,2))
            if (par)
                parfor iChn = 1:nChn
                    [Sp_tmp, Spktimes_tmp] = spikeextract(mu(iChn,:), thresh{iChn}, FS,artefact,artefact_high);
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
                    [Sp_tmp, Spktimes_tmp] = spikeextract(mu(iChn,:), thresh{iChn}, FS,artefact,artefact_high);
                    if (~isempty(Sp_tmp) && any(Sp_tmp(:)<(-200)))
                        stop=0;
                    end
                    NSp_tmp = length(Spktimes_tmp);
                    if NSp_tmp > 1
                        SpMat = [Spktimes_tmp+double(time) Sp_tmp];
                        NSp_tmp = size(SpMat,1);
                        sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
                        NSp(iChn) = NSp(iChn) + NSp_tmp;
                    end
                end
            end
            time = time + T*1e3;
            dispstat(sprintf('Progress %03.2f%%',100*((N)/ntimes)),'timestamp');
        end
        if (size(mu,2) < FS * T)
            chk = 0;
        end
    end
    
    %% Calculate a zero-condition spike template and r2t
    expfile=dir('*_exp_datafile_*');
    if ~isempty(expfile)
        trig = loadTrig(0); loadNREP; loadNTrials; loadDUALSTIM;
        if ~exist('DUALSTIM','var')
            simstim=load(expfile.name,'simultaneous_stim');
            DUALSTIM=simstim.simultaneous_stim;
        end
        TrialParams = loadTrialParams;
        TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == 1)); %% find the baseline trial (TrialParams(:2) = 1)'s trial number (TrialParams(:,1)). 
%         if DUALSTIM
%             TrialParams(DUALSTIM:DUALSTIM:length(TrialParams),:) = [];
%         end
        % DUALSTIM = 2;
        TrialParams = TrialParams(DUALSTIM:DUALSTIM:length(TrialParams));
        trig = trig(TrialParams);
        nTrig = length(trig);
        spWN = [-300 300]; template = cell(1,nChn); r2t = cell(1,nChn); spt = cell(1,nChn);
        for c = 1:nChn
            for n = 1:nTrig
                % trig(n)/30 convert trig from samples to ms
                % spt{c} = spike waveforms around each baseline trigger
                spt{c} = [spt{c}; sp{c}(sp{c}(:,1) > trig(n)/30+spWN(1) & sp{c}(:,1) < trig(n)/30+spWN(2),2:end)];
            end
            template{c} = mean(cell2mat(spt(c)),1); % average basedline spike shape 
            r2t{c} = mean((cell2mat(spt(c)) - template{c}).^2); 
        end
        %% Save the datafiles
        disp('Saving spikes');
        for iChn = 1:nChn
            sp{iChn} = sp{iChn}(1:NSp(iChn),:); % first column: spike time; following columns: spike waveform values
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
end