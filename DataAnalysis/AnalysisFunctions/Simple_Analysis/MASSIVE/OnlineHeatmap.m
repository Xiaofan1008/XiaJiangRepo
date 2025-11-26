function [array,arraycolmod,xvalues,yvalues]=OnlineHeatmap(varargin)
%nChn is the total number of channels recorded simulataneously
%ID is the stimulus ID w/ stim params. Check loadTrialInfo to see IDs
%% Online sigmoid generation
tic
filepath = pwd;
dName='amplifier';

if nargin==0
    nChn=input('Please input the total number of channels recorded from all headstages:\n');
    ID=input('Please input the trial ID you want to process:\n');
elseif nargin==1
    fprintf('Not enough input parameters\n')
    nChn=input('Please input the total number of channels recorded from all headstages:\n');
    ID=input('Please input the trial ID you want to process:\n');
elseif nargin==3
    nChn=varargin{1};
    ID=varargin{2};
    Spk_array=varargin{3};%to skip 
else
    nChn=varargin{1};
    ID=varargin{2};
end
%% get small amplifier file
[~,frequency_parameters]=read_Intan_RHS2000_file;
FS=frequency_parameters.amplifier_sample_rate;
if exist('Spk_array','var')==0

vFID = fopen([filepath filesep dName '.dat'],'r');
BREAK=1;
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * FS);
T=256;
if t_len < T
   T=t_len;
end
temp=dir('*.trig.dat');
if isempty(temp)
    cleanTrig_sabquick;
end
trig=loadTrig(0);
Spk_array=cell(nChn,1,1);
N=0;
if nargin~=3
while BREAK
    v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
    N=N+1;
    if ~size(v,2)
        BREAK = 0;
    else
        trigN=trig(trig>T*FS*(N-1) & trig<=T*FS*N)-T*FS*(N-1);
        for chnit=1:nChn
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
            Spk_array{chnit}=[Spk_array{chnit}; Spktimes+(N-1)*T*1000, Sp];
        end
    end
end
fclose all;
end
toc
end

%% group FR response
tic
trig = loadTrig(0);
TP = loadTrialParams;
tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,1:2));
stimchn=trialinfo(trialinfo(:,1)==ID,2);
theseTrig = trig(tID)./30;
theseTrig=unique(theseTrig);
nT=length(theseTrig);
BIN = [2 8]; 
%%

chnpos=Depth(1);
array=zeros(16,nChn/16);
for recordchn=1:nChn
    sp = Spk_array{recordchn};
    for tr = 1:nT
        theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
        array(chnpos==recordchn)=array(chnpos==recordchn)+length(theseSp);
    end
end
%% this will need to be modified depending on the array configuration
arraycolmod=zeros(16,nChn/16);
arraytemp=array;
TI=loadTrialInfo;
TI=cell2mat(TI(2:end,[1 2 4:end]));
TIstimelect=TI((TI(:,1)==ID) & (TI(:,17)~=-1),2);
array(TIstimelect)=nan;
for i=1:nChn/64
    arraycolmod(1:16,(1:4)+((i-1)*4))=[array(:,1+(i-1)*4),array(:,3+(i-1)*4),array(:,4+(i-1)*4),array(:,2+(i-1)*4)];
end
arraycolmod=flipud(arraycolmod);
yvalues=fliplr({'Deep','2','3','4','5','6','7','8','9','10','11','12','13','14','15','Shallow'});
xvalues=[];
for i=1:nChn/16
xvalues=[xvalues {num2str(i)}];
end
arraycolmod=(arraycolmod./nT)./(BIN(2)-BIN(1)).*1000;
array=arraytemp;

%plotting
%  figure
%  heatmap(xvalues,yvalues,arraycolmod)




% stimchn(stimchn==0)=[];
% shank=zeros(length(stimchn),1);
% shank(stimchn<17)=1;
% shank(stimchn>16)=4;
% shank(stimchn>32)=2;
% shank(stimchn>48)=3;
% stimchn(stimchn>48)=stimchn(stimchn>48)-48;
% stimchn(stimchn>32)=stimchn(stimchn>32)-32;
% stimchn(stimchn>16)=stimchn(stimchn>16)-16;
% fprintf('Stimulating channels: \n')
% for i=1:length(stimchn)
%     fprintf(['S' num2str(shank(i)), 'E', num2str(stimchn(i)) '\n'])
% end

end