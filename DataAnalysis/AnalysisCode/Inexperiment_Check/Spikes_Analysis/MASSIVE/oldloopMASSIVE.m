
D = dir; % D is a struct ... first elements are '.' and '..' used for navigation.
folder = fileparts(which('MASSIVE_ANALYSIS_LOOP')); % Determines filepath to folder containing your .m file even if your is not stored in same dir as long as path is added
addpath(genpath(folder)); % adds folder and all subfolders to path for all functions
fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd([fileptochange filesep 'Data' filesep 'Rat_006'])%change working directory to where data is stored - currently manually input
D_data=dir;


for k = 3:length(D_data) % avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    if strcmp(currD(1:12),'Pen2_accross')||strcmp(currD(1:11),'Pen2_across')
        continue
    end
    fprintf([currD,'\n']);
    try
    MASSIVE_ANALYSIS_LOOP([D_data(k).folder filesep currD]); %Perform data analysis and save files in appropriate folders
    catch
    end 
    
end

%% old loop

% for k = 3:length(D_data) % avoid using the first ones
%     currD = D_data(k).name; % Get the current subdirectory name
%     MASSIVE_ANALYSIS_LOOP([D_data(k).folder filesep currD]); %Perform data analysis and save files in appropriate folders
% end