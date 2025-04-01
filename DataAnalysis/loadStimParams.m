function StimParams =  loadStimParams(data_folder)
global DATA_FOLDER;
data_files = dir(fullfile(DATA_FOLDER, '*_exp_datafile_*.mat'));
if isempty(data_files)
    error('No data file found in %s', DATA_FOLDER);
end
data_file = fullfile(DATA_FOLDER, data_files(1).name);
load(data_file, 'StimParams');
end

% data = dir('*_exp_datafile_*.mat');
% if isempty(data)
%     return
% end
% data = data.name;
% load(data,'StimParams');
% '/Users/xiaofan/Desktop/PhD Study/Data/Sabrina/S3E2_9elect_001_210511_104603'
