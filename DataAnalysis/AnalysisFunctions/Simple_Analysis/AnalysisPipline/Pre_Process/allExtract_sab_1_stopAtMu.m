function allExtract_sab_1_stopAtMu(dName,filepath,T,par,artefact,artefact_high)
%% This is allExtract_sab_1.m, unchanged, EXCEPT that it returns right
% after writing <name>.mu_sab.dat and does not run the "Calculate SP"
% section (threshold-based spike detection / .sp.mat saving). Spike
% detection is handled by a separate script now.
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
sp = cell(1, nChn); %#ok<NASGU>
thresh = cell(1, nChn);
SCALEFACTOR = 10;
chk = 1; N = 0; time = 0;
NSp = zeros(1,nChn); %#ok<NASGU>
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
lv_fid = fopen([dName '.dat'],'r'); %#ok<NASGU>
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T); % number of blocks
%% Generate filters
[Mufilt,Lfpfilt] = generate_Filters; %#ok<NASGU>
LfpNf = length(Lfpfilt); %#ok<NASGU>
MuNf = length(Mufilt);
%% Calculate thresholds
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

%% ===== STOP HERE: spike detection / .sp.mat saving intentionally removed =====
% The original function continued into a "Calculate SP" section here
% (threshold-based spike detection + saving <name>.sp.mat). That section
% is handled by a separate spike-detection script now, so this function
% returns immediately after writing the MU file.
if ~isempty(v_fid) && v_fid ~= -1
    fclose(v_fid);
end
fprintf('Done: wrote %s.mu_sab.dat (thresholds computed but not saved; spike detection skipped).\n', name);
end