%% ============================================================
%   QC_SaveBadChannels.m
%   Manually enter bad channels for each stimulation set.
%   Saves results in order-sensitive format:
%       BadCh_perSet{setID} = [ch1 ch2 ch3 ...];
% ============================================================

clear; clc;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ================= USER SETTINGS =================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX012/Xia_Exp1_Sim6_251125_181554';
Electrode_Type = 1;    % 0 rigid, 1 single-shank flex, 2 four-shank flex

%% ================= CHECK FOLDER =================
if ~isfolder(data_folder)
    error('Folder not found: %s', data_folder);
end
cd(data_folder);

fprintf('QC Bad Channel Selection\nFolder: %s\n\n', data_folder);

%% ================= BASE NAME =================
parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u)>=4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

%% =============== LOAD StimParams (order-sensitive) ===============
fileDIR = dir('*_exp_datafile_*.mat');
assert(~isempty(fileDIR), 'No exp data file found.');

S = load(fileDIR(1).name, 'StimParams', 'simultaneous_stim', 'E_MAP', 'n_Trials');

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;
E_MAP             = S.E_MAP;
n_Trials          = S.n_Trials;

%% === Decode stimulation sets (ORDER-SENSITIVE) ===
stimNames = StimParams(2:end,1);
E_NAME    = E_MAP(2:end);
[~, idx_all] = ismember(stimNames, E_NAME);

stimChPerTrial_all = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);

    % KEEP ORDER (do NOT use unique)
    v = idx_all(rr);
    v = v(v>0);
    stimChPerTrial_all{t} = v(:)';     % row vector
end

comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial_all{t};
    comb(t,1:numel(v)) = v;
end

[uniqueComb, ~, combClass] = unique(comb, 'rows', 'stable');
nSets = size(uniqueComb,1);

%% =============== Create storage structure ===============
BadCh_perSet = cell(nSets,1);

fprintf('Detected %d stimulation sets (order-sensitive).\n\n', nSets);

%% ============= User Input Loop =============
for s = 1:nSets
    stimSeq = uniqueComb(s, uniqueComb(s,:) > 0);
    stimLabel = strjoin(string(stimSeq), '-');

    fprintf('-----------------------------------------------------\n');
    fprintf(' Stimulation Set %d :  [%s]\n', s, stimLabel);
    fprintf(' Enter bad channels for this set (using Depth_s index).\n');
    fprintf(' Example:  [3 14 22]   or   [] if no bad channels.\n');

    user_input = input(sprintf(' Bad channels for Set %d = ', s));

    % ensure numeric row vector
    if isempty(user_input)
        BadCh_perSet{s} = [];
    else
        BadCh_perSet{s} = unique(user_input(:)');
    end

    fprintf(' Saved for Set %d:  [%s]\n\n', s, num2str(BadCh_perSet{s}));
end

%% ============= SAVE RESULTS =============
save_name = sprintf('%s.BadChannels.mat', base_name);

save(save_name, 'BadCh_perSet', 'uniqueComb', 'combClass', 'stimChPerTrial_all');

fprintf('\n=============================================\n');
fprintf(' Manual bad-channel decisions saved to:\n   %s\n', save_name);
fprintf('=============================================\n');