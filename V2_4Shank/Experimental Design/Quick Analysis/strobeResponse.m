function strobeResponse(nChn,BIN,SPACING)
%% This function attempts to produce a strobe response plot from the LFP
% You must be in the right directory, and have an LFP and trigger data file
% present
if nargin<3
    SPACING = 20;
end
if nargin<2
    BIN = [-50 400];
end
if nargin<1
    nChn = 32;
end
%% Load Data
tmp = dir([pwd filesep '*.lfp.dat']);
if ~isempty(tmp)
    sz = tmp.bytes/2;
    tmp = [pwd filesep tmp.name];
    l_fid = fopen(tmp,'r');
    lfp = fread(l_fid,[nChn, sz/nChn],'float');
    fclose(l_fid);
else
    error('No LFP data file present\n');
end
tmp = dir([pwd filesep '*.trig.dat']);
if ~isempty(tmp)
    sz = tmp.bytes/2;
    tmp = [pwd filesep tmp.name];
    t_fid = fopen(tmp,'r');
    trig = fread(t_fid,[1, sz],'double');
    fclose(t_fid);
    trig = trig ./ 30;
    trig = cast(trig,'int64');
    trig(trig(1)+BIN(1) <= 0) = [];
    trig(trig(end)+BIN(2) >= size(lfp,2)) = [];
    nTrig = size(trig,2);
else
    error('No trig data file present\n');
end
%% Process Data
mLFP = zeros(size(lfp,1),diff(BIN)+1);
x = BIN(1):1:BIN(2);
count = 0;
for t = 2:nTrig-1
    mLFP = mLFP + lfp(:,trig(t)+BIN(1):trig(t)+BIN(2));
    count = count + 1;
end
mLFP = mLFP ./ count;
%% Plot Data
figure; hold on;
d = Depth;
for c = 1:nChn
    plot(x,mLFP(d(c),:)+(c-1)*SPACING,'Color','k','LineWidth',2');
end
xlabel('Time (ms)'); ylabel('Voltage (\muV)'); title('Strobe Flash Response');
YLIM = ylim; line([0 0],[YLIM(1) YLIM(2)],'Color','r','LineWidth',2);
beautifyPlot(24);
end