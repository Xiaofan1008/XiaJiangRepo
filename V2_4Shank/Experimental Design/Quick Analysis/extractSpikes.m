function extractSpikes(nChn,FS,threshfac)
%% This function attempts to extract spikes from the high-pass filtered data
% You must be in the right directory, and have a high-pass filtered data
% file present
if nargin<3
    threshfac = -4.5;
end
if nargin<2
    FS = 30000;
end
if nargin<1
    nChn = 32;
end
sp = cell(1, nChn);
thresh = cell(1, nChn);
NSp = zeros(1,nChn);
Ninput = 1e6;
ARTMAX = 1.5e3;
name = strsplit(pwd,'\');
name = name{end};
name = name(1:end-14);
%% Load in the amplifier waveform data
try
    fileinfo = dir('amplifier_dn.dat');
    v_fid = fopen('amplifier_dn.dat', 'r');
catch
    try
        fileinfo = dir('amplifier.dat');
        v_fid = fopen('amplifier.dat', 'r');
    catch
        error('No amplifier data found\n');
    end
end
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / 256);
%% Generate filters
[Mufilt,~] = generate_Filters;
MuNf = length(Mufilt);
%% Calculate thresholds
for iChn = 1:nChn
    sp{iChn} = zeros(ntimes * ceil(double(256) * 20), FS * 1.6 / 1e3 + 1 + 1);    
end
artchk = zeros(1,nChn);
munoise = cell(1,nChn);
while sum(artchk) < nChn       
    v = fread(v_fid, [nChn, FS*20], 'int16') * 0.195;       
    if ~isempty(v)
        % Checks for artifact
        for iChn = 1:nChn
            if isempty(thresh{iChn})
                munoise{iChn} = [];
                mu = conv(v(iChn,:),Mufilt);
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
fclose(v_fid);
%% Extract spiking
m_fid = fopen([name '.mu.dat'],'r');
chk = 1; N = 0; time = 0;
while (chk && N < Ninput)
    N = N + 1;
    mu = fread(m_fid,[nChn, FS*256],'short') ./ 10;
    for iChn = 1:nChn
        [Sp_tmp, Spktimes_tmp] = spikeextract(mu(iChn,:), thresh{iChn}, FS);
        NSp_tmp = length(Spktimes_tmp);
        if NSp_tmp > 1
            SpMat = [Spktimes_tmp+double(time) Sp_tmp];
            NSp_tmp = size(SpMat,1);
            sp{iChn}(NSp(iChn)+1:NSp(iChn)+NSp_tmp,:) = SpMat;
            NSp(iChn) = NSp(iChn) + NSp_tmp;
        end
    end
    time = time + 256*1e3;
    if (size(mu,2) < FS * 256)
        chk = 0;
    end
end
%% Calculate a zero-condition spike template and r2t
if ~isempty(dir('*_exp_datafile_*'))
    trig = loadTrig(0); n_REP = loadNREP;
    TrialParams = loadTrialParams;
    TrialParams = cell2mat(TrialParams(cell2mat(TrialParams(:,2)) == 1));
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
        % Measure Spike Width at 1/3rd spike threshold
        meanSp = mean(sp{c}(:,2:50));
        T = meanSp(13); % Definitionally, this is the trough
        P = max(meanSp(13:end)); % The peak of the spike.
        H = (P - T)/2 + T; % The value we will calculate width at. 50% the difference between peak and trough
        [x,~] = polyxpoly(1:49,meanSp,1:49,H.*ones(1,49)); % Find the intersections of meanSp with H
        Spkwidth(c) = diff(x(1:2)) ./ (FS/1000); % Convert to msec
    end
    %% Save the datafiles 
    for iChn = 1:nChn
        sp{iChn} = sp{iChn}(1:NSp(iChn),:);
    end
    save([name '.sp.mat'],'sp','thresh','threshfac','template','r2t','Spkwidth','-v7.3');
else
    trig = loadTrig(0); nTrig = length(trig);
    spWN = [-300 300]; template = cell(1,nChn); r2t = cell(1,nChn); spt = cell(1,nChn);
    for c = 1:nChn
        for n = 1:nTrig
            spt{c} = [spt{c}; sp{c}(sp{c}(:,1) > trig(n)/30+spWN(1) & sp{c}(:,1) < trig(n)/30+spWN(2),2:end)];
        end
        template{c} = mean(cell2mat(spt(c)),1);
        r2t{c} = mean((cell2mat(spt(c)) - template{c}).^2);              
%       % Measure Spike Width at 1/3rd spike threshold
        T = template{c}(13); % Definitionally, this is the trough
        P = max(template{c}(13:end)); % The peak of the spike.
        H = (P - T)/2 + T; % The value we will calculate width at. 50% the difference between peak and trough
        [x,~] = polyxpoly(1:49,template{c},1:49,H.*ones(1,49)); % Find the intersections of meanSp with H
        Spkwidth(c) = diff(x(1:2)) ./ (FS/1000); % Convert to msec
    end
    for iChn = 1:nChn
        sp{iChn} = sp{iChn}(1:NSp(iChn),:);
    end    
    save([name '.sp.mat'],'sp','thresh','threshfac','template','r2t','Spkwidth','-v7.3');
end
end

