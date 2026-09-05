%% Loop through subdirectories function

D = dir; % D is a struct ... first elements are '.' and '..' used for navigation.
folder = fileparts(which('MASSIVE_ANALYSISFLASH_LOOP')); % Determines filepath to folder containing your .m file even if your is not stored in same dir as long as path is added
addpath(genpath(folder)); % adds folder and all subfolders to path for all functions
fileptochange=folder(1:end-8);%Path name minus '/MASSIVE'
cd([fileptochange filesep 'Data' filesep 'Rat_021_FLASH'])%change working directory to where data is stored - currently manually input
D_data=dir;
idStr = getenv('SLURM_ARRAY_TASK_ID'); % get task ID string
arrayTaskID = str2double(idStr); % conv string to number
currD = D_data(arrayTaskID+2).name; % Get the current subdirectory name  % avoid using the first ones '.' '..'
fprintf([currD,'\n']);
MASSIVE_ANALYSISFLASH_LOOP([D_data(arrayTaskID+2).folder filesep currD]); %Perform data analysis and save files in appropriate folders

%% old loop

% for k = 3:length(D_data) % avoid using the first ones
%     currD = D_data(k).name; % Get the current subdirectory name
%     MASSIVE_ANALYSIS_LOOP([D_data(k).folder filesep currD]); %Perform data analysis and save files in appropriate folders
% end