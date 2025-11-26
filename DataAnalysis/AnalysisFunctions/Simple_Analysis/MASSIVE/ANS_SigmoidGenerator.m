function [spk_all]=ANS_SigmoidGenerator(varargin)
%chns input 'A-001,B-013'
%nChn is the total number of channels recorded simulataneously
%% Online sigmoid generation
tic
%user input channels of choice?? 
if nargin==0
    chns=input('Please input channels to analyse like so "A-001,B-016":\n','s');
    nChn=input('Please input the total number of channels recorded from all headstages\n');
elseif nargin==1
    fprintf('Not enough input parameters\n')
    chns=input('Please input channels to analyse like so "A-001,B-016":\n','s');
    nChn=input('Please input the total number of channels recorded from all headstages\n');
else
    chnnum=varargin{1};
    nChn=varargin{2};
    timestart=varargin{3};
    timeend=varargin{4};
    IDs=varargin{5}-4:varargin{5};
    Spk_array=varargin{6};
end
ti=loadTrialInfo;
IDs=[ti{end,1} IDs];

trig=loadTrig(0);
 [Spike_trialstruct,baslinespike_trialstruct] = OnlinesortTrials(trig,Spk_array,chnnum,timestart,timeend);

%% plot sigmoid -  this can only do one channel at the moment
%rejig so record electrode just loops through the outside
spk_all=zeros(1,5);
      
        for recordchn=1:length(chnnum)
            
            for ID=1:size(IDs,2)
                t=['ID_' num2str(IDs(ID))];
                spkcount=zeros(40,1);
                spstim=zeros(40,1);
                bspk=zeros(40,1);
                for j=1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t),'AsArray',true),2)
                    t2=['Trial_' num2str(j)];
                    spkcount(j)=size(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1)-size(baslinespike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1);
                    spstim(j)=size(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1);
                    bspk(j)=size(baslinespike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1);
                end
                spk_all(ID)=mean(spkcount(1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t),'AsArray',true),2)))/((timeend-timestart)/1000);
                
            end
        end

toc



end