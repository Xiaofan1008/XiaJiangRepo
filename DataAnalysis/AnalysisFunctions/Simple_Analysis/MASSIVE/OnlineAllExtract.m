function [Spike_trialstruct,baslinespike_trialstruct]=OnlineAllExtract(chnnum,nChn,varargin)

%get small amplifier file
if nargin>2
    timestart=varargin{1};
    timeend=varargin{2};
end
filepath = pwd;
dName='amplifier';
FS=30000;
vFID = fopen([filepath filesep dName '.dat'],'r');
BREAK=1;
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * 30000);
T=256;
if t_len < T
   T=t_len;
end
temp=dir('*.trig.dat');
if isempty(temp)
    cleanTrig_sabquick;
end
trig=loadTrig(0);
temp=dir('*.sp.mat');
temp2=dir('Spk_array.mat');
if isempty(temp) && isempty(temp2)
Spk_array=cell(size(chnnum,1),1,1);
N=0;
thresh=ones(nChn,1);
while BREAK
    v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
    N=N+1;
    
    if ~size(v,2)
        BREAK = 0;
    else
        if (N ~= 1)
            v = [hv v];% Append saved data to the start of the next loop
        end
        for chnit=1:length(chnnum)
            v(chnit,:) = simpleBlank(v(chnit,:),N,T,trig,1,FS);
            if N>1
                [Sp,Spktimes,thresh(chnit)]=OnlineSpkExtract(v((chnit),FS+1:end),thresh(chnit));%don't double up on spikes
            else
                [Sp,Spktimes,thresh(chnit)]=OnlineSpkExtract(v((chnit),:),thresh(chnit));
            end
            if ~isempty(Sp)
                Spk_array{chnit}=[Spk_array{chnit}; Spktimes+(N-1)*T*1000, Sp];
            end
        end
        hv = v(:,end-FS+1:end);
        dispstat(sprintf('Progress %03.2f%%',(100*(N*T/t_len))),'timestamp');
    end
end
fclose all;
save('Spk_array.mat','Spk_array')
elseif ~isempty(temp)
    Spk_array=loadSpikes;
else
    load('Spk_array.mat','Spk_array');
end
%process max four current levels - find max/min and select closest to 0.25
%and 0.75
%min 2 current levels
[Spike_trialstruct,baslinespike_trialstruct] = OnlinesortTrials(trig,Spk_array,chnnum,timestart,timeend);

save('IDstruct_online.mat','Spike_trialstruct','baslinespike_trialstruct')
end