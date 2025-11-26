% Generate a MUA filter
[Mufilt,Lfpfilt] = generate_Filters;
MuNf = length(Mufilt);
LfpNf = length(Lfpfilt);
[b,a] = butter(3,300/(FS/2),'high');
if exist([filepath, '\amplifier_dn.dat'],'file')
    v_fid = fopen([filepath, '\amplifier_dn.dat'],'r');
else
    v_fid = fopen([filepath, '\amplifier.dat'],'r');
end
nData = diff(BIN)*30+1;  % Number of datapoints
% Load in a little data for thresholding
if exist([filepath, '\amplifier_dn.dat'],'file')
    v_info = dir([filepath, '\amplifier_dn.dat']);
else
    v_info = dir([filepath, '\amplifier.dat']);
end
threshBYTES = min([(FS/1e3)*(100*1000),v_info.bytes/2]);
v_thresh = fread(v_fid,[nChn, threshBYTES],'int16') * 0.195;
% Generate a threshold
thresh = zeros(2,nChn);
if ~exist('threshfac','var')
    threshfac = -3.5;
end
for c = 1:nChn
    chan = depth(c);
    tmp = conv(v_thresh(chan,:),Mufilt);
    filt_v = tmp(1,MuNf/2:threshBYTES+MuNf/2-1);
    thresh(1,c) = threshfac*median(abs(filt_v))./0.6795;
    thresh(2,c) = -threshfac*median(abs(filt_v))./0.6795;
end