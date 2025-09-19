function [electsig,electnonsig,electall, electsig75,electnonsig75,electall75,electsig25,electnonsig25,electall25,p,stimshank,othershank, meansig50, meansig75,meansig25,stdersig50, stdersig75,stdersig25,stimChn,currentavg50,currentavg25,currentavg75,electfitratio]=RATIOLOOP_MASSIVE(SubDir_Path)
%%for analysing data on massive and saving denoised files
folder = fileparts(which('RATIOLOOP_MASSIVE')); % Determines filepath to folder containing your .m file.
addpath(genpath(folder)); % Add that folder plus all subfolders to the path.
%fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd(SubDir_Path)%change working directory to where data is stored - currently manually input
ratioALLnormalise;
loadStimChn;
load('Significant.mat','ratiosingduals')
sigratio=ratiosingduals;
load('NonSignificant.mat','ratiosingdualn')
nonsigratio=ratiosingdualn;
nonsigratio(isnan(nonsigratio))=0;
sigratio(isnan(sigratio))=0;
nonsigratio(isinf(nonsigratio))=0;
sigratio(isinf(sigratio))=0;
electall=sum(sigratio+nonsigratio)./sum((sigratio+nonsigratio)~=0);
electsig=mean(sigratio(sigratio~=0));
electnonsig=mean(nonsigratio(nonsigratio~=0));

load('Significant.mat','ratiosingduals75')
sigratio=ratiosingduals75;
load('NonSignificant.mat','ratiosingdualn75')
nonsigratio=ratiosingdualn75;
nonsigratio(isnan(nonsigratio))=0;
sigratio(isnan(sigratio))=0;
nonsigratio(isinf(nonsigratio))=0;
sigratio(isinf(sigratio))=0;
electall75=sum(sigratio+nonsigratio)./sum((sigratio+nonsigratio)~=0);
electsig75=mean(sigratio(sigratio~=0));
electnonsig75=mean(nonsigratio(nonsigratio~=0));

load('Significant.mat','ratiosingduals25')
sigratio=ratiosingduals25;
load('NonSignificant.mat','ratiosingdualn25')
nonsigratio=ratiosingdualn25;
nonsigratio(isnan(nonsigratio))=0;
sigratio(isnan(sigratio))=0;
nonsigratio(isinf(nonsigratio))=0;
sigratio(isinf(sigratio))=0;
electall25=sum(sigratio+nonsigratio)./sum((sigratio+nonsigratio)~=0);
electsig25=mean(sigratio(sigratio~=0));
electnonsig25=mean(nonsigratio(nonsigratio~=0));

load('Significantnm.mat','meansig50')
load('Significantnm.mat','meansig75')
load('Significantnm.mat','meansig25')

load('Significantnm.mat','stdersig50')
load('Significantnm.mat','stdersig75')
load('Significantnm.mat','stdersig25')
load('Significantnm.mat','currentavg50')
load('Significantnm.mat','currentavg75')
load('Significantnm.mat','currentavg25')
load('Significantnm.mat','electfitratio')


% 
% load('Significant.mat','sigtest','counter')

% sigtest1=sigtest(1:16,:);
% sigtest = sigtest(any(sigtest,2),:);
% [h,p] = ttest(sigtest(1:32,1),sigtest(1:32,2));
% 
ap=ratiosingduals;
 norm=jbtest(ap(ap~=0));
 if norm==1
    [~,p]=ttest(ap(ap~=0),0);
 else
     p=signrank(ap(ap~=0));
 end
 filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=loadMapNum;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
if nChn==64
s1sig=ratiosingduals(1:16);
s1sig=mean(s1sig(s1sig~=0));
s2sig=ratiosingduals(17:32);
s2sig=mean(s2sig(s2sig~=0));
s3sig=ratiosingduals(33:48);
s3sig=mean(s3sig(s3sig~=0));
s4sig=ratiosingduals(49:64);
s4sig=mean(s4sig(s4sig~=0));
loadStimChn;
if stimChn(1)<17
    i=1;
elseif stimChn(1)<33
    i=2;
elseif stimChn(1)<49
    i=3;
elseif stimChn(1)<65
    i=4;
end

if stimChn(2)<17
    if i==1
        stimshank=s1sig;
        othershank=(s2sig+s3sig+s4sig)/3;
    elseif i==2
        stimshank=(s1sig+s2sig)/2;
        othershank=(s3sig+s4sig)/2;
    elseif i==3
        stimshank=(s1sig+s3sig)/2;
        othershank=(s2sig+s4sig)/2;
    elseif i==4
        stimshank=(s1sig+s4sig)/2;
        othershank=(s2sig+s3sig)/2;
    end
elseif stimChn(2)<33
    if i==1
        stimshank=(s1sig+s2sig)/2;
        othershank=(s3sig+s4sig)/2;
    elseif i==2
        stimshank=(s2sig);
        othershank=(s1sig+s3sig+s4sig)/3;
    elseif i==3
        stimshank=(s2sig+s3sig)/2;
        othershank=(s1sig+s4sig)/2;
    elseif i==4
        stimshank=(s2sig+s4sig)/2;
        othershank=(s1sig+s3sig)/2;
    end
elseif stimChn(2)<49
    if i==1
        stimshank=(s1sig+s3sig)/2;
        othershank=(s2sig+s4sig)/2;
    elseif i==2
        stimshank=(s2sig+s3sig)/2;
        othershank=(s1sig+s4sig)/2;
    elseif i==3
        stimshank=(s3sig);
        othershank=(s1sig+s4sig+s2sig)/3;
    elseif i==4
        stimshank=(s3sig+s4sig)/2;
        othershank=(s1sig+s2sig)/2;
    end
elseif stimChn(2)<65
    if i==1
        stimshank=(s1sig+s4sig)/2;
        othershank=(s2sig+s3sig)/2;
    elseif i==2
        stimshank=(s2sig+s4sig)/2;
        othershank=(s1sig+s3sig)/2;
    elseif i==3
        stimshank=(s3sig+s4sig)/2;
        othershank=(s1sig+s2sig)/2;
    elseif i==4
        stimshank=(s4sig);
        othershank=(s1sig+s2sig+s3sig)/3;
    end
end
    
save('sigshank.mat','s1sig','s2sig','s3sig','s4sig')
else
    stimshank=0;
    othershank=0;
end

end