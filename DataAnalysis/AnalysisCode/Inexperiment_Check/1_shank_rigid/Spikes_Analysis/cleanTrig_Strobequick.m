%% Initial Setup
filepath = pwd;
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
if nChn>32
    E_Mapnumber=1;
else
    E_Mapnumber=0;
end

name = pwd;
[FP,name,ext] = fileparts(name);
name = name(1:end-14);
%% Deal with the digital lines
disp('Cleaning the digital lines');
fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);
stimDig = flip(find(digital_in == 1)); % Fix for finding the falling line instead of the rising line
visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
dt = diff(visDig);
kill = dt == -1;
visDig(kill) = [];
nStamps = max([length(stimDig); length(visDig)]);
time_stamps = nan(1,nStamps);
time_stamps(1,1:length(stimDig)) = flip(stimDig);
%time_stamps(2,1:length(visDig)) = flip(visDig);
% if ~isempty(time_stamps)
%     if isnan(time_stamps(2,2))
%         time_stamps = time_stamps(1,:);
%     else
%         time_stamps = time_stamps(2,:);
%     end
% end
% Correct for 200 us delay
loadDelay;
time_stamps = time_stamps + delay*(1e-6)*FS;
% At this point, we need to correct for jitter
% d = Depth(E_Mapnumber); 
d = Depth_s(0); % 0-Single Shank Rigid, 1-Single Shank Flex, 2-Four Shanks Flex
SKIP = 1;
try 
    loadStimChn;
catch 
    stimChn = 16;
end
trig = time_stamps;
if (SKIP == 0)
warning('off','signal:findpeaks:largeMinPeakHeight');
    for t = 1:length(trig)
        v = loadRAW(trig(t));
        [~,adj] = findpeaks(abs(diff(v(d(stimChn(1)),:))),'MinPeakHeight',500);
        if isempty(adj)
            continue;
        end
        % adj = adj(1) - (500*FS/1e3);
        adj = adj(1);
        trig(t) = trig(t) + adj;
        v = loadRAW(trig(t));
    end
end
% Adjust trigger lines to match datafile
% loadNTrials; loadNREP;
% if length(trig) > n_Trials
%     % Delete trigger lines
%     TrialParams = loadTrialParams; AMP = loadAMP; DUR = loadDUR;
%     tParams = cell2mat(TrialParams(1:end,2));
%     for nS = 1:length(CHN)
%         for nD = 1:length(DUR)
%             for nA = 1:length(AMP)
%                 id = nA + (nD-1)*length(AMP) + (nS-1)*length(DUR)*length(AMP);
%                 tmp = find(tParams == id);
%                 mk = tmp(n_REP);
%                 trig(tmp(tmp > mk)) = NaN;
%             end
%         end
%     end
%     trig(isnan(trig)) = [];
% end
if trig(1)<30000
    trig(1)=[];
end
fprintf('There are %.0f digital lines here.\n',length(trig));
trig_fid = fopen([name '.trig.dat'],'w');
fwrite(trig_fid,trig,'double');
fclose(trig_fid);