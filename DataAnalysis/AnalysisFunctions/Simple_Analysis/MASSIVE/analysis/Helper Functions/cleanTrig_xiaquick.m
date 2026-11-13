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
d = Depth(E_Mapnumber); SKIP = 1;
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
TP=loadTrialParams;
numtrial=find(diff(cell2mat(TP(:,1)))~=0,1,'first');
MissingTrig = 0;
if length(trig)<(length(TP)/numtrial)
    %Dealing with missing triggers- find first ‘no stim’ trial with a pulse
    %and then find next ‘no stim’ trial along and assign -500 to all trigs
    %in between, then shift by one value to add the missing trigger. The
    %-500 is used as a flag later on in the code.
    file = pwd;
    mDIR = dir('amplifier.dat');
    mNAME = mDIR.name;
    mFID = fopen([file filesep mNAME], 'r');
    fprintf('Correcting for missing triggers.\n');
    loadStimParams;
    StimParams(1,:) = [];
    StimParams = StimParams(numtrial:numtrial:end,:);
    
    MissingTrig = 1;
    trig_diff = diff(trig)/FS; % time difference between adjacent triggers
    trig_std = std(trig_diff);
    trig_ave = mean(trig_diff);
    peak_cutoff = trig_ave + 2 * trig_std;
    [val,locs] = findpeaks(trig_diff, 'MinPeakProminence', 2*trig_std, 'MinPeakHeight', peak_cutoff); % find locations for missing triggers
    trig_to_check = trig_diff(locs);
    % in the runExp code, there is 5s pause set for each 500 trial,
    % resulting in large trigger time gaps. 
    num_block = (length(TP)/numtrial)/500; % 500 trials with 5s pause = 1 block, calculate num of block
    num_peak = fix(num_block); % how many high peaks related to the 5s pause. 
    num_block_to_process = ceil(num_block); % number of block to process
    [~,trig_pause_idx] = maxk(val,num_peak); % find the location of the peaks of 5s pause
    trig_to_check(trig_pause_idx) = [];
    trig_to_check_idx = locs;
    trig_to_check_idx(trig_pause_idx) = [];
    thresh1 = trig_ave + trig_std;
    thresh2 = 2*trig_ave + trig_std;
    new_trig = trig;
    shift = 0; 
    for i=1:length(trig_to_check)
        idx_org = trig_to_check_idx(i); % original index in trig_diff
        idx = idx_org + shift;      % shift to the current new_trig index
        Insert_Trig=-500;
        if trig_to_check(i) > thresh1 && trig_to_check(i) <= thresh2
            % one missing trigger
                new_trig = [new_trig(1:idx) Insert_Trig new_trig(idx+1:end)];
                shift = shift + 1;
        elseif trig_to_check(i) > thresh2
            % two missing trigger
                new_trig = [new_trig(1:idx) Insert_Trig Insert_Trig new_trig(idx+1:end)];
        end
    end
else 
    new_trig = trig;
end

if length(new_trig)~=size(TP,1)/numtrial
    new_trig=[new_trig -500*ones(size(TP,1)/numtrial-length(trig),1)'];
end
trig_fid = fopen([name '.trig.dat'],'w');
fwrite(trig_fid,new_trig,'double');
fclose(trig_fid);

%    % error('missing triggers')
%     %%need to fix code below to deal with whitenoise
%     NoStimTrials=find(cell2mat(StimParams(:,16))==-1);
%     for i=length(NoStimTrials):-1:1
%         if NoStimTrials(i)>length(trig)
%             NoStimTrials(i)=[];
%         end
%     end
%     Trigs_to_Check = trig(NoStimTrials);
%     BIN = [-50 50];
%     for t=1:length(Trigs_to_Check)
%         OFFSET = cast(nChn*2*(Trigs_to_Check(t)+(BIN(1)*FS/1e3)),'int64');
%         fseek(mFID,OFFSET,'bof'); % used to point to specific point in the data indicated by offset from begginning of file
%         v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
%         if any(v(5,:)>2000) || any(v(15,:)>2000)
%             if t~=1
%             Start_trig=find(trig==Trigs_to_Check(t-1));
%             End_trig=find(trig==Trigs_to_Check(t));
%             % trig(Start_trig+1:End_trig)=-500;
%             Insert_Trig=-500;
%             trig = [trig(1,1:Start_trig) Insert_Trig trig(1,Start_trig+1:end)];
%             Trigs_to_Check = trig(NoStimTrials);
%             else
%                 End_trig=find(trig==Trigs_to_Check(t));
%                 Insert_Trig=-500;
%                 % trig(1:End_trig)=-500;
%                 % trig = [Insert_Trig trig];
%                 trig = [trig(1:End_trig-1), Insert_Trig, trig(End_trig:end)];
%                 Trigs_to_Check = trig(NoStimTrials);
% 
%             end
%         end
%     end
% 
% end

% 
% if length(trig)~=size(TP,1)/numtrial
%     trig=[trig -500*ones(size(TP,1)/numtrial-length(trig),1)'];
% end
% trig_fid = fopen([name '.trig.dat'],'w');
% fwrite(trig_fid,trig,'double');
% fclose(trig_fid);