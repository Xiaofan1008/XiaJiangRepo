% This helper function loads in post-mortem LFP data in order to help
% detrend
function lfp = detrendLFP(lfp,chn,filepath, time_stamp)
%% Load in Intan data
generateIntanData;
%% Generate time stamps
generateTimeStamps;
%% Load in the datafile and select the Trial ID for analysis
load([filepath, datafile],'AMP','DUR','n_REP','CHN','StimParams','n_Trials','TrialParams');
if (n_Trials ~= nStamps)
    disp('Warning: Digital Lines do not match expected number of trials');
    %return;
end
%% Generate filters and thresholds
generateFiltersandThresholds;
X = BIN(1):BIN(2);
%% Load the desired data
OFFSET = cast(nChn*2*(FS/1e3)*(time_stamps(time_stamp)+BIN(1)),'int64');
fseek(v_fid,OFFSET,'bof'); % Find the relevant timestamp
v = fread(v_fid,[nChn,(FS/1e3)*diff(BIN)],'int16') * 0.195;
b_v = BlankArte(v(chn,:));
tmp = conv(b_v,Lfpfilt);
pm_LFP = tmp(1,LfpNf/2:FS/1e3:nData+LfpNf/2-1);
end