function TrialParams = loadTrialParams

global DATA_FOLDER;
data_files = dir(fullfile(DATA_FOLDER, '*_exp_datafile_*.mat'));
if isempty(data_files)
   error('No data file found in %s', DATA_FOLDER);
end
data_file = fullfile(DATA_FOLDER, data_files(1).name);
load(data_file, 'TrialParams');
TrialParams = TrialParams(2:end, :);
end

% tparams = dir('*_exp_datafile_*.mat');
% if isempty(tparams)
%     TrialParams = [];
%     return
% end
% tparams = tparams.name;
% load(tparams,'TrialParams');
% TrialParams = TrialParams(2:end,:); % 1: Trial number 2: Trial ID 3: Channel
% end