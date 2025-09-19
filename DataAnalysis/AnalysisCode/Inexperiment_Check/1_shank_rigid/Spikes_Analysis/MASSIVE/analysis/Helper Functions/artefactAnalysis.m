%% Script for analysing artefact removal strategies
close all
if ~exist('filepath','var')
    filepath = uigetdir;
end
datafile = dir([filepath, '\*_exp_datafile*.mat']);
datafile = ['\', datafile.name];
info = read_Intan_RHS2000_file(filepath,'\info.rhs',0);
FS = info.frequency_parameters.amplifier_sample_rate;       % Sampling frequency
nChn = length(info.amplifier_channels);
% Time-Stamps
generateTimeStamps;
% Data-File
load([filepath, datafile],'AMP','DUR','n_REP','CHN','StimParams','n_Trials','TrialParams');
if (n_Trials ~= nStamps)
    disp('Warning: Digital Lines do not match expected number of trials');
    %return;
end
uniqueTrials = length(AMP)*length(DUR)*length(CHN);
% Filters
[Mufilt,Lfpfilt] = generate_Filters;
MuNf = length(Mufilt);
LfpNf = length(Lfpfilt);
BIN = [-100, 200]; % Extraction window for reading in data
nData = diff(BIN)*30+1; % Number of datapoints
X = 1/30:1/30:nData/30; % For plotting in msec
% Open the data
v_fid = fopen([filepath, '\amplifier.dat'],'r');
% Load in a sample series of data
PULSES = 10; % Number of pulses to stack into the data
CHN = 23; % Which channel to look at
v = zeros(PULSES,nData);
for r = 1:PULSES
    OFFSET = cast(nChn*2*(FS/1e3)*(time_stamps(r)+BIN(1)),'int64');
    fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
    tmp = fread(v_fid,[nChn,(FS/1e3)*diff(BIN)+1],'int16') * 0.195;
    v(r,:) = tmp(CHN,:);
end
%% Start trying to remove artefact
S = struct;
% High=Frequency Area
% Strategy 1. Blanking, then taper filters
S.s1 = v; % Retains some trend. Long ringing.
for r = 1:PULSES
    S.s1(r,:) = BlankArte(v(r,:),99,103);
    tmp = conv(S.s1(r,:),Mufilt);
    S.s1(r,:) = tmp(1,MuNf/2:nData+MuNf/2-1);
end
% Strategy 2. Blanking, detrending, then taper filters
S.s2 = v; % Removes trend. Long ringing.
for r = 1:PULSES
    S.s2(r,:) = BlankArte(v(r,:),99,103);
    S.s2(r,:) = locdetrend(S.s2(r,:),FS,[10/1e3 1/1e3]);
    tmp = conv(S.s2(r,:),Mufilt);
    S.s2(r,:) = tmp(1,MuNf/2:nData+MuNf/2-1);
end
% Strategy 3. Blanking, detrending, re_blanking.
S.s3 = v;
for r = 1:PULSES
    S.s3(r,:) = BlankArte(v(r,:),99,103);
    S.s3(r,:) = locdetrend(S.s3(r,:),FS,[10/1e3 1/1e3]);
    tmp = conv(S.s3(r,:),Mufilt);
    S.s3(r,:) = tmp(1,MuNf/2:nData+MuNf/2-1);
    S.s3(r,98*30:104*30) = 0;
end
% Strategy 4. Blanking, Butterworth filter
S.s4 = v; % Less ringing. Butterworth seems to do a better job of removing trend.
[b,a] = butter(3,300/(FS/2),'high');
for r = 1:PULSES
    S.s4(r,:) = BlankArte(v(r,:),99,103);    
    S.s4(r,:) = filtfilt(b,a,S.s4(r,:));
end
% Strategy 5. Tailored blanking window based on artefact features.
S.s5 = v;
for r = 1:PULSES
    % Detrend first    
    S.s5(r,:) = locdetrend(S.s5(r,:),FS,[10/1e3 1/1e3]);
    % Calculate features of artefact
    MAX = find(S.s5(r,:) >= floor(max(S.s5(r,:))),1);
    MIN = find(S.s5(r,:) <= ceil(min(S.s5(r,:))),1);
    if (MAX > MIN)
        % Negative facing artefact
        XSTART = find(S.s5(r,1:MIN) >= -50,1,'last');
        XSTOP = find(S.s5(r,MAX:end) <= 50,1,'first') + MAX;
    else
        % Positive facing artefact
        XSTART = find(S.s5(r,1:MAX) <= 50,1,'last');
        XSTOP = find(S.s5(r,MIN:end) >= -50,1,'first') + MIN;
    end
    % Prevent estimates wildly exceeding the window
    if isempty(XSTART) || (XSTART) < 98*30
        XSTART = 98*30;
    end    
    if isempty(XSTOP) || (XSTOP) > 104*30
        XSTOP = 104*30;
    end
    S.s5(r,:) = BlankArte(S.s5(r,:),floor(XSTART/30),ceil(XSTOP/30));
    S.s5(r,:) = filtfilt(b,a,S.s5(r,:));
end
%% Conclusions:
   % Extremely difficult to deal with artefact
   % in situations where the amplifier has saturated.
   % Amp Settle is absolutely necessary.
   % Must investigate better strategies - at a minimum,
   % might get away with robust PCA spike sorting.
   % Need to talk to Yan/Nic/Mau for help with PCA though.
% Output
% for sol = 1:5 % Number of solutions
%     figure; hold on
%     for r = 1:PULSES
%         plot(X,S.(['s' num2str(sol)])(r,:))
%     end
%     xlim([80 120])
%     ylim([-1000 1000])
% end
%% Spike Extraction and Sorting
% Generate a threshold
threshfac = -4.2;
v_info = dir([filepath, '\amplifier.dat']);
threshBYTES = min([(FS/1e3)*(100*1000),v_info.bytes/2]);
fseek(v_fid,(nChn*2*(FS/1e3)*threshBYTES),'bof');
v_thresh = fread(v_fid,[nChn, threshBYTES],'int16') * 0.195;
v_thresh = v_thresh(CHN,:);
tmp = conv(v_thresh,Mufilt);
filt_v = tmp(1,MuNf/2:threshBYTES+MuNf/2-1);
thresh = threshfac*median(abs(filt_v))./0.6795;
sample = filt_v(7.96e5:8.08e5); % Magic numbers drawn from looking at the data
[WAVES,Sp_tmp] = spikeextract(sample,thresh,FS);
SAMPLESPIKES = [Sp_tmp WAVES];
[~,SAMPLEPCS,~] = spikepcs(SAMPLESPIKES);
% From the above, I see three spike templates
TEMPLATE(1,:) = [0.4, -0.25];
TEMPLATE(2,:) = [0.05, -0.65];
TEMPLATE(3,:) = [0.4, 0.05];
TEMPLATE(4,:) = [0.1, 0.1];
TEMPLATE(5,:) = [0.7, -0.03];
TEMPLATE(6,:) = [0.7, -2];
% Extract spikes
for sol = 1:5
    for r = 1:PULSES
        [WAVES,Sp_tmp] = spikeextract(S.(['s' num2str(sol)])(r,:),thresh,FS);
        if (Sp_tmp ~= 0)
            S.(['sp' num2str(sol)]).(['rep' num2str(r)]) = [Sp_tmp WAVES];
        else
            S.(['sp' num2str(sol)]).(['rep' num2str(r)]) = 0;
        end
    end
end
% Spike Sorting
% At present this is poorly calibrated. It tends to allow non-spikes
% through, and reject spikes at the same time.
T = 1.4;
for sol = 1:5
    for r = 1:PULSES
        [~,pcs,~] = spikepcs(S.(['sp' num2str(sol)]).(['rep' num2str(r)]));
        for sp = 1:size(pcs,2)
            chk = 0;
            for t = 1:size(TEMPLATE,1)
                if size(pcs,2) < 2
                    continue;
                end
                if TEMPLATE(t,1) > 0
                    if pcs(sp,1) > T*TEMPLATE(t,1) || pcs(sp,1) < (1/T)*TEMPLATE(t,1)
                        continue;
                    end
                    if TEMPLATE(t,2) > 0
                        if pcs(sp,2) > T*TEMPLATE(t,2) || pcs(sp,2) < (1/T)*TEMPLATE(t,2)
                            continue;
                        end
                    else
                        if pcs(sp,2) < T*TEMPLATE(t,2) || pcs(sp,2) > (1/T)*TEMPLATE(t,2)
                            continue;
                        end
                    end
                else
                    if pcs(sp,1) < T*TEMPLATE(t,1) || pcs(sp,1) > (1/T)*TEMPLATE(t,1)
                        continue;
                    end
                    if TEMPLATE(t,2) > 0
                        if pcs(sp,2) > T*TEMPLATE(t,2) || pcs(sp,2) < (1/T)*TEMPLATE(t,2)
                            continue;
                        end
                    else
                        if pcs(sp,2) < T*TEMPLATE(t,2) || pcs(sp,2) > (1/T)*TEMPLATE(t,2)
                            continue;
                        end
                    end
                end                
                chk = 1;
            end
            if (chk == 0)
                S.(['sp' num2str(sol)]).(['rep' num2str(r)])(sp,:) = NaN;
            end
        end
        S.(['sp' num2str(sol)]).(['rep' num2str(r)])(isnan(S.(['sp' num2str(sol)]).(['rep' num2str(r)])(:,1)),:) = [];
    end    
end