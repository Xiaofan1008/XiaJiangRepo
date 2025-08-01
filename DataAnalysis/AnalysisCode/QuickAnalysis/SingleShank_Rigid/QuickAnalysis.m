%% Quick Analysis Code
clear all
close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Sabrina/Experimental_Design'));
% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-500; %removes spikes below this threshold
artefact_high=500; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)?
secondstoanalyse=8; %How long after the trigger do you want to analyse spikes for(ms)?
printspiking=0;
par=0;

%% Initial
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
dName='amplifier';
mem_check=dir('amplifier.dat');
T = mem_check.bytes ./ (2 * nChn * FS);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * FS);
if t_len < (Overall_time_to_analyse+Startpoint_analyse)
    error('Time to analyse exceeds the data recording period')
end
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
elseif Overall_time_to_analyse~=0
    T=Overall_time_to_analyse;
end

%% Load Trigger
cleanTrig_sabquick;
% cleanTrig_xiaquick;
trig = loadTrig(0);

%% MUA analysis
% Blank artifact

% -- Check trigger -- 
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end
trig = loadTrig(0);
theseTrig = trig;
shortbytes=2; % 2 bytes per data
if isempty(dir([dName '_dn_xia.dat'])) % denoise file does not exit
    if ~empty(trig) % trigger exist
        % -- File initialization -- 
        info = dir([filepath filesep dName '.dat']);
        info = info.bytes/2;
        nL = (ceil(info / (nChn*FS*double(T)))+1);
        vFID = fopen([filepath filesep dName '.dat'],'r'); % read data file
        vdnFID = fopen([filepath filesep dName '_dn_xia.dat'],'W'); % write data file
        N = 1;
        CHK = split(filepath,'_');
        blank = 1; 
        for i = 1:size(CHK,1)
            if strcmp(CHK{i},'CSD') || strcmp(dName,'analogin') || isempty(dir('*exp_datafile*.mat'))
                blank = 0;
            end
        end 

        dispstat('','init');
        if (blank)
            dispstat(sprintf('Blanking artefact . . .'),'keepthis','n');
        else
            dispstat(sprintf('Skipping blanking . . .'),'keepthis','n');
        end

        BREAK = 1;

        while (N) && (BREAK)
            % -- load data --
            if strcmp(dName,'analogin') % load 'analogin.dat' 
                v = fread(vFID,[nChn, (FS * T)],'uint16');
                v = (v - 32768) .* 0.0003125;
            elseif strcmp(dName,'amplifier') % load 'amplifier.dat'
                if Startpoint_analyse~=0 % shift start point 
                    offset=Startpoint_analyse*FS*nChn*shortbytes;%offset from beginning of file
                    ftell(vFID)
                    fseek(vFID,offset,'bof');
                    ftell(vFID)
                end
                v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195; % one chunk data read
            end
            if ~size(v,2) % if no data
                BREAK = 0;
            end

            % -- Append saved data to the start to the next loop --
            if (N~=1)
                v = [hv v];
            end

            if (blank)
                d = Depth_s(0); 
                % 0 - single shank rigid;
                % 1 - single shank flex;
                % 2 - four shank flex
                for iChn = 1:nChn
                    % -- load trial ID --
                    TrialParams = loadTrialParams; 
                    TrialParams = cell2mat(TrialParams(:,2));
                    num_elect=min(diff(find(diff(TrialParams)~=0))); % number of electrodes per trial
                    TrialParams=TrialParams(1:num_elect:end,:);
                    % -- load 'num_pulse per train' and 'Pulse train period' 
                    loadStimParams;
                    StimParams_pulseinfo=cell2mat(StimParams(2:num_elect:end,8:9)); 
                    data = v(iChn,:);
                    nChn = size(data,1);

                    % -- trials cross two chunks --
                    % search the triggers around [-0.1s, 0.1s] of the
                    % boundary of prev and current chunk
                    missed_trials = TrialParams(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
                    missed_trig = trig(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
                    % remove missing triggers
                    missed_trials((missed_trig==-500))=[];
                    missed_trig((missed_trig==-500))=[];

                    % -- triggers login --
                    % find trials within current data chunk
                    index=trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000;
                    trials = TrialParams(index);
                    pulsetrain_info=StimParams_pulseinfo(index,:);
                    trig = trig(index);

                    if ~isempty(missed_trig) % trials cross two chunks 
                        missed_pulseinfo=cell2mat(StimParams(temp+1,8:9));    missed_pulseinfo((missed_trig==-500))=[];
                        trials = [missed_trials trials];
                        trig = [missed_trig trig];
                        pulsetrain_info = [missed_pulseinfo; pulsetrain_info];
                    end

                    % -- realign time of triggers --
                    if (N>1) % compensate 1s data overlap in the following chunks
                        trig = trig - ((N-1)*T*FS) + FS;
                    end
                    wd=dir;

                    % -- safe range after artifact blanking -- 
                    shifttime=min(diff(trig))-0.101*FS;

                    % -- Dynamic artifact
                    for t = 1:length(trig) % Loop through trials
                        


                    end

                end

            end


        end



    end
end