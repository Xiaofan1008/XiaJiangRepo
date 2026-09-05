%ratioLoop
Ratnum='Rat_015';
D = dir; % D is a struct ... first elements are '.' and '..' used for navigation.
folder = fileparts(which('MASSIVE_ANALYSIS_LOOP')); % Determines filepath to folder containing your .m file even if your is not stored in same dir as long as path is added
addpath(genpath(folder)); % adds folder and all subfolders to path for all functions
fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'

cd([fileptochange filesep 'Data' filesep Ratnum])%change working directory to where data is stored - currently manually input
D_data=dir;

electsigall=0;
electnonsigall=0;
electallall=0;
electsigall75=0;
electnonsigall75=0;
electallall75=0;
electsigall25=0;
electnonsigall25=0;
electallall25=0;
stimshankall=0;
othershankall=0;
pall=0;
num=0;
meansig50array=[];
meansig7525array=[];
stdersig50array=[];
stdersig7525array=[];
Stimchnall=[];
Ratiosingelect=[];
currentvariation50=[];
currentvariation75=[];
currentvariation25=[];
MSEalltrials=0;
clear electarray75 electarray25 electarray stimshankarray othershankarray realelectarray75 realelectarray25 
for k = 3:length(D_data) % avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    if strcmp(currD,'pen1_dualvary_15_002_200904_130226')
        continue
    end
    try
    
        [electsig,electnonsig,electall, electsig75,electnonsig75,electall75,electsig25,electnonsig25,electall25,p,stimshank,othershank,meansig50, meansig75,meansig25, stdersig50, stdersig75,stdersig25,Stimchn, currentavg50,currentavg25,currentavg75, electfitratio]=RATIOLOOP_MASSIVE([D_data(k).folder filesep currD])

        electsigall=electsig+electsigall;
        electarray(k-2)=electsig;
        
        electnonsigall=electnonsig+electnonsigall;
        electallall=electallall+electall;
        realelectarray75(k-2)=electsig75;
        realelectarray25(k-2)=electsig25;
        if electsig75>electsig25
            electsigall75=electsig75+electsigall75;
        electarray75(k-2)=electsig75;
        electnonsigall75=electnonsig75+electnonsigall75;
        electallall75=electallall75+electall75;
        electsigall25=electsig25+electsigall25;
        electarray25(k-2)=electsig25;
        electnonsigall25=electnonsig25+electnonsigall25;
        electallall25=electallall25+electall25;
            
        else
            electsigall75=electsig25+electsigall75;
            electarray75(k-2)=electsig25;
            electnonsigall75=electnonsig25+electnonsigall75;
            electallall75=electallall75+electall25;
            electsigall25=electsig75+electsigall25;
            electarray25(k-2)=electsig75;
            electnonsigall25=electnonsig75+electnonsigall25;
            electallall25=electallall25+electall75;
            
        end
        
        pall=p+pall;
        stimshankall=stimshank+stimshankall;
        stimshankarray(k-2)=stimshank;
        othershankall=othershank+othershankall;
        othershankarray(k-2)=othershank;
       %[electfit,MSEtotal]=SigmoidfunctionNoDual([D_data(k).folder filesep currD]);
       %fprintf([num2str(MSEtotal) newline])
       %MSEalltrials=MSEtotal+MSEalltrials;
       num=num+1;
       meansig50array=[meansig50array meansig50];
       meansig7525array=[meansig7525array meansig75 meansig25];
         
       stdersig50array=[stdersig50array stdersig50];
       stdersig7525array=[stdersig7525array stdersig75 stdersig25];
       SubDir_Path=[D_data(k).folder filesep currD];
       fprintf(['End of Analysis for: %s ' newline], SubDir_Path)
       subdirectpath.(['N' num2str(k)])=SubDir_Path;
       Stimchnall=[Stimchnall; Stimchn'];
       currentvariation50=[currentvariation50; currentavg50];
       currentvariation75=[currentvariation75; currentavg75];
       currentvariation25=[currentvariation25; currentavg25];

       Ratiosingelect=[Ratiosingelect;electfitratio];
    catch
    end

    cd(folder)

end
electsigerror=std(electarray(electarray~=0));
electsigerror75=std(electarray75(electarray75~=0));
electsigerror25=std(electarray25(electarray25~=0));

stimshankerror=std(stimshankarray(stimshankarray~=0));
othershankerror=std(othershankarray(othershankarray~=0));
cd([fileptochange filesep 'Data' filesep Ratnum])%change working directory to where data is stored - currently manually input
save('Ratio_all.mat', 'electsigall', 'electnonsigall', 'electallall', ...
    'electsigall75', 'electnonsigall75', 'electallall75', ...
    'electsigall25', 'electnonsigall25', 'electallall25', 'MSEalltrials','pall','num',...
    'meansig50array', 'meansig7525array', 'stdersig50array', 'stdersig7525array',...
    'realelectarray75','realelectarray25','Stimchnall', 'subdirectpath', 'Ratiosingelect',...
    'currentvariation50','currentvariation75','currentvariation25')

fprintf(['Number of successfully completed trials: ' num2str(num) newline])
