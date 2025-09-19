function Spk_array=OnlineFRresponse(varargin)
%nChn is the total number of channels recorded simulataneously
%ID is the stimulus ID w/ stim params. Check loadTrialInfo to see IDs

tic
filepath = pwd;
dName='amplifier';

%user input channels of choice?? 
if nargin==0
    Chns=input('Please input the channels to average i.e. 1:64:\n');
    nChn=input('Please input the total number of channels recorded from all headstages:\n');
    ID=input('Please input the trial ID you want to process:\n');
elseif nargin==1
    fprintf('Not enough input parameters\n')
    Chns=input('Please input the channels to average i.e. 1:64:\n');
    nChn=input('Please input the total number of channels recorded from all headstages:\n');
    ID=input('Please input the trial ID you want to process:\n');
else
    Chns=varargin{1};
    nChn=varargin{2};
    ID=varargin{3};
end
chnsplit=split(Chns,",");
chnnum=zeros(size(chnsplit,1),1);
for i=1:size(chnsplit,1)
    if strcmp(chnsplit{i}(1),'A')
        chnnum(i)=str2double(chnsplit{i}(end-2:end))+1;
    elseif strcmp(chnsplit{i}(1),'B')
        chnnum(i)=str2double(chnsplit{i}(end-2:end))+33;
    elseif strcmp(chnsplit{i}(1),'C')
         chnnum(i)=str2double(chnsplit{i}(end-2:end))+65;
    elseif strcmp(chnsplit{i}(1),'D')
        chnnum(i)=str2double(chnsplit{i}(end-2:end))+97;
    end
end
Chns=chnnum;
%get small amplifier file
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
    cleanTrig;
end
trig=loadTrig(0);
Spk_array=cell(nChn,1,1);
N=0;
while BREAK
    v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
    N=N+1;
    if ~size(v,2)
        BREAK = 0;
    else
        trigN=trig(trig>T*FS*(N-1) & trig<=T*FS*N)-T*FS*(N-1);
        for iterate=1:length(Chns)
            chnit=Chns(iterate);
            for trigit=1:length(trigN)
                if size(v,2)>=trigN(trigit)+FS*0.100
                    v(chnit,trigN(trigit):trigN(trigit)+FS*0.100)=v(chnit,trigN(trigit):trigN(trigit)+FS*0.100)-(mean(v(chnit,trigN(trigit)+FS*0.005:trigN(trigit)+FS*0.025))-mean(v(chnit,trigN(trigit)-FS*0.020:trigN(trigit)-FS*0.001)));
                else
                    v(chnit,trigN(trigit):size(v,2))=v(chnit,trigN(trigit):size(v,2))-(mean(v(chnit,trigN(trigit)+FS*0.005:trigN(trigit)+FS*0.025))-mean(v(chnit,trigN(trigit)-FS*0.020:trigN(trigit)-FS*0.001)));
                end
                v(chnit,trigN(trigit)-FS*0.001:trigN(trigit)+FS*0.002)=interpolate(v(chnit,trigN(trigit)-FS*0.001:trigN(trigit)+FS*0.002),1);%%%blank artefact
            end
            if N==1
                thresh=1;
            end
            [Sp,Spktimes,thresh]=OnlineSpkExtract(v(chnit,:),thresh);
            Spk_array{iterate}=[Spk_array{iterate}; Spktimes+(N-1)*T*1000, Sp];
        end
    end
end
fclose all;
toc


%% group FR response
tic
trig = loadTrig(0);
TP = loadTrialParams;
tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
theseTrig = trig(tID)./30;
nT=length(theseTrig);
BIN = [-200 200]; 
SMOOTHING = 2; MAX = 400;
xdata = [];
ydata = [];

for iterate=1:length(Chns)
    recordchn=Chns(iterate);
    sp = Spk_array{iterate};
    for tr = 1:nT
        theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
        for i = 1:length(theseSp)
            xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
            ydata = [ydata, tr*(MAX/nT)];
        end
    end
end

Z = hist(xdata,0:400); %#ok<HIST>
    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
    rate = (1000/nT)*conv(Z,window);
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING)./length(Chns);
    figure; hold on; axM = gca;
    plot(axM,rate,'b','LineWidth',2);
    ylim(axM,[0 MAX]); xlim(axM,[150 275]);
    xlabel(axM,'Time (ms)');
    ylabel(axM,'Firing Rate (Sp/s)');
    xticks(axM,150:50:275);
    xticklabels(axM,-50:50:75);
    xline(200,'r')
    title('All channels collapsed')
    
toc
%% raster
%ID=34; % there will be an increase at 100ms due to offset
for recordchn=1:length(chnnum)
Online_raster(ID,Spk_array{recordchn})
end

end