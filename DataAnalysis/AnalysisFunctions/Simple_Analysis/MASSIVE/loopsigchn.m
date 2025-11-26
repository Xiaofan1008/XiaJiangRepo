D = dir; % D is a struct ... first elements are '.' and '..' used for navigation.
folder = fileparts(which('loopsigchn')); % Determines filepath to folder containing your .m file even if your is not stored in same dir as long as path is added
addpath(genpath(folder)); % adds folder and all subfolders to path for all functions
fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd([fileptochange filesep 'Data' filesep 'Rat_004'])%change working directory to where data is stored - currently manually input
D_data=dir;
filepath = pwd;
nChn=32;
stimelectrodes=zeros(nChn,1);
numsigchn=zeros(nChn,1);
for k = 3:length(D_data) % avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    cd ([D_data(k).folder filesep currD])
    load('sigchn.mat')
    stimelectrodes=stimelectrodes+whichchn;
    numsigchn=numsigchn+numofsigperchn;
    cd (D_data(k).folder)
end
save('overallsigstim.mat','numsigchn','stimelectrodes')