%% Artifact blank-time estimate from processed spikes
% Estimate the post-stimulation recovery time using the earliest detected
% spike for each recording channel and trial.

clear;
clc;

addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));
%% Settings
data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
save_folder = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Artifact_Time_Est';

Electrode_Type = 2;      % 0: 1-rigid, 1: 1-flex, 2: 4-shank flex
max_blank_ms   = 4;

save_result = true;
make_plot   = true;

if ~isfolder(save_folder)
    mkdir(save_folder);
end
%% Load data

cd(data_folder);
parts = split(data_folder, filesep);
last_folder = parts{end};
underscores = strfind(last_folder, '_');
if numel(underscores) > 4
    base_name = last_folder(1:underscores(end-1)-1);
else
    base_name = last_folder;
end
fprintf('\nDataset: %s\n', last_folder);

% Load filtered spikes
ssd_file = [base_name '.sp_xia_SSD.mat'];
assert(isfile(ssd_file), ...
    'Cannot find %s', ssd_file);
S = load(ssd_file);
sp_corr = S.sp_corr;

% Load triggers
if isempty(dir('*.trig.dat'))
    cur_dir = pwd;
    cleanTrig_sabquick;
    cd(cur_dir);
end
trig = loadTrig(0);
FS = 30000;


% Load stimulation parameters
%% Load stimulation parameters
fileDIR = dir('*exp_datafile*.mat');
assert(~isempty(fileDIR), 'No exp_datafile found.');
S = load(fileDIR(1).name, ...
    'StimParams', ...
    'simultaneous_stim', ...
    'E_MAP', ...
    'n_Trials');

StimParams         = S.StimParams;
simultaneous_stim  = S.simultaneous_stim;
E_MAP              = S.E_MAP;
n_Trials           = S.n_Trials;

%% Decode amplitude

trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
trialAmps(trialAmps == -1) = 0;

%% Decode stimulation sets

stimNames = StimParams(2:end,1);
E_NAME = E_MAP(2:end);
[~, idx_all] = ismember(stimNames, E_NAME);
stimChPerTrial = cell(n_Trials,1);
for t = 1:n_Trials
    rr = (t-1)*simultaneous_stim + (1:simultaneous_stim);
    v = idx_all(rr);
    v = v(v > 0);
    stimChPerTrial{t} = v(:)';
end
comb = zeros(n_Trials, simultaneous_stim);
for t = 1:n_Trials
    v = stimChPerTrial{t};
    comb(t,1:numel(v)) = v;
end

[uniqueComb, ~, combClass] = unique(comb, 'rows', 'stable');
nSets = size(uniqueComb,1);


%% Decode PTD

PTD_all_us = cell2mat(StimParams(2:end,6));
PTD_matrix = reshape(PTD_all_us, ...
    simultaneous_stim, [])';
trialPTD_ms = max(PTD_matrix, [], 2) / 1000;

%% Recording channel order
d = Depth_s(Electrode_Type);
nRecCh = numel(d);
%% Estimate earliest spike for each channel x trial

fprintf('Estimating recovery times...\n');

Trial = [];
StimSet = [];
Amplitude = [];
PTD = [];
Condition = strings(0,1);

ChannelIndex = [];
Channel = [];

SearchEnd_ms = [];
EarliestSpike_ms = [];


for tr = 1:n_Trials
    amp_val = trialAmps(tr);
    % Exclude 0 uA
    if amp_val == 0
        continue;
    end

    ptd_val = trialPTD_ms(tr);
    if ptd_val == 0

        condition = "Sim";
        search_end = max_blank_ms;
    else
        condition = "Seq";
        % Do not allow the second sequential pulse to enter the search
        search_end = min(max_blank_ms, ptd_val);

    end


    t0 = trig(tr) / FS * 1000;


    for ich = 1:nRecCh

        ch = d(ich);

        earliest = NaN;

        if ch <= numel(sp_corr) && ~isempty(sp_corr{ch})

            spike_times = sp_corr{ch}(:,1) - t0;

            early_spikes = spike_times( ...
                spike_times > 0 & ...
                spike_times < search_end);

            if ~isempty(early_spikes) % Estimate time detection
                earliest = min(early_spikes);
                earliest = earliest -1;
                % earliest = earliest -1.4;
                earliest = max(earliest,0);
            end
        end
        Trial(end+1,1) = tr;
        StimSet(end+1,1) = combClass(tr);
        Amplitude(end+1,1) = amp_val;
        PTD(end+1,1) = ptd_val;
        Condition(end+1,1) = condition;
        ChannelIndex(end+1,1) = ich;
        Channel(end+1,1) = ch;
        SearchEnd_ms(end+1,1) = search_end;
        EarliestSpike_ms(end+1,1) = earliest;
    end
end

RawTable = table( ...
    Trial, ...
    StimSet, ...
    Amplitude, ...
    PTD, ...
    Condition, ...
    ChannelIndex, ...
    Channel, ...
    SearchEnd_ms, ...
    EarliestSpike_ms);

%% Channel summary
% Average trials first for each recording channel.
groupVars = { ...
    'StimSet', ...
    'Amplitude', ...
    'Condition', ...
    'ChannelIndex', ...
    'Channel'};
[G, ChannelGroups] = findgroups(RawTable(:,groupVars));
NTrials = splitapply(@numel, ...
    RawTable.EarliestSpike_ms, G);
NValid = splitapply(@(x) sum(~isnan(x)), ...
    RawTable.EarliestSpike_ms, G);
Mean_ms = splitapply(@mean_omitnan, ...
    RawTable.EarliestSpike_ms, G);
Median_ms = splitapply(@median_omitnan, ...
    RawTable.EarliestSpike_ms, G);
SD_ms = splitapply(@std_omitnan, ...
    RawTable.EarliestSpike_ms, G);
ValidPercent = NValid ./ NTrials * 100;

ChannelSummary = ChannelGroups;

ChannelSummary.NTrials = NTrials;
ChannelSummary.NValid = NValid;
ChannelSummary.ValidPercent = ValidPercent;

ChannelSummary.Mean_ms = Mean_ms;
ChannelSummary.Median_ms = Median_ms;
ChannelSummary.SD_ms = SD_ms;


%% Stimulation-set summary
% Average channel means within each stimulation set.
validChannel = ...
    ChannelSummary.NValid > 0 & ...
    ~isnan(ChannelSummary.Mean_ms);
CS = ChannelSummary(validChannel,:);
groupVars = { ...
    'StimSet', ...
    'Amplitude', ...
    'Condition'};
[G, SetGroups] = findgroups(CS(:,groupVars));
NChannels = splitapply(@numel, ...
    CS.Mean_ms, G);
Mean_ms = splitapply(@mean, ...
    CS.Mean_ms, G);
SD_ms = splitapply(@std, ...
   CS.Mean_ms, G);
SEM_ms = SD_ms ./ sqrt(NChannels);

SetSummary = SetGroups;

SetSummary.NChannels = NChannels;
SetSummary.Mean_ms = Mean_ms;
SetSummary.SD_ms = SD_ms;
SetSummary.SEM_ms = SEM_ms;

%% Dataset summary
% Average stimulation-set means.
groupVars = { ...
    'Amplitude', ...
    'Condition'};
[G, DatasetGroups] = findgroups(SetSummary(:,groupVars));
NSets = splitapply(@numel, ...
    SetSummary.Mean_ms, G);
Mean_ms = splitapply(@mean, ...
    SetSummary.Mean_ms, G);
SD_ms = splitapply(@std, ...
    SetSummary.Mean_ms, G);
SEM_ms = SD_ms ./ sqrt(NSets);

DatasetSummary = DatasetGroups;

DatasetSummary.NSets = NSets;
DatasetSummary.Mean_ms = Mean_ms;
DatasetSummary.SD_ms = SD_ms;
DatasetSummary.SEM_ms = SEM_ms;
%% Print summary
fprintf('\nDataset summary\n');
fprintf('-------------------------------------------\n');
fprintf('Amp     Condition     Mean      SEM    Sets\n');
fprintf('-------------------------------------------\n');
for i = 1:height(DatasetSummary)
    fprintf('%2g uA    %-4s        %5.3f     %5.3f     %d\n', ...
        DatasetSummary.Amplitude(i), ...
        DatasetSummary.Condition(i), ...
        DatasetSummary.Mean_ms(i), ...
        DatasetSummary.SEM_ms(i), ...
        DatasetSummary.NSets(i));
end
fprintf('-------------------------------------------\n');
%% Check valid spike percentage
fprintf('\nChannel validity\n');
fprintf('Mean valid trials: %.1f %%\n', ...
    mean(ChannelSummary.ValidPercent, 'omitnan'));
fprintf('Median valid trials: %.1f %%\n', ...
    median(ChannelSummary.ValidPercent, 'omitnan'));
%% Plot
if make_plot
    figure('Color','w');
    hold on;
    conditions = ["Sim","Seq"];
    for ci = 1:numel(conditions)
        idx = DatasetSummary.Condition == conditions(ci);
        T = DatasetSummary(idx,:);
        T = sortrows(T, 'Amplitude');
        errorbar( ...
            T.Amplitude, ...
            T.Mean_ms, ...
            T.SEM_ms, ...
            '-o', ...
            'LineWidth',1.5, ...
            'MarkerSize',6);

    end
    xlabel('Stimulation amplitude (\muA)');
    ylabel('Estimated recovery time (ms)');
    legend({'Simultaneous','Sequential'}, ...
        'Location','best', ...
        'Box','off');
    xlim([min(DatasetSummary.Amplitude)-0.5, ...
          max(DatasetSummary.Amplitude)+0.5]);
    ylim([0 max_blank_ms]);
    box off;
end

%% Save
if save_result
    result_file = fullfile(save_folder, ...
        [base_name '_DX010_ArtifactEst.mat']);  
    save(result_file, ...
        'RawTable', ...
        'ChannelSummary', ...
        'SetSummary', ...
        'DatasetSummary', ...
        'uniqueComb', ...
        'max_blank_ms');
    
    fprintf('Saved: %s\n', result_file);
end

%% Functions
function y = mean_omitnan(x)
x = x(~isnan(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end


function y = median_omitnan(x)
x = x(~isnan(x));
if isempty(x)
    y = NaN;
else
    y = median(x);
end
end

function y = std_omitnan(x)
x = x(~isnan(x));
if numel(x) < 2
    y = NaN;
else
    y = std(x);
end
end