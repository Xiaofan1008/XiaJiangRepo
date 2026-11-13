%% Cross-dataset artifact recovery-time analysis
clear;
clc;

%% Settings
result_folder = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Artifact_Time_Est';
save_result = false;
make_plot = true;

%% Find files
files = dir(fullfile(result_folder, '*_ArtifactEst.mat'));
assert(~isempty(files), 'No artifact estimate files found.');

fprintf('Found %d files.\n', numel(files));

%% Combine file-level results
AllData = table();

for i = 1:numel(files)
    file_path = fullfile(files(i).folder, files(i).name);
    S = load(file_path, 'DatasetSummary');

    if ~isfield(S, 'DatasetSummary') || isempty(S.DatasetSummary)
        fprintf('Skipped: %s\n', files(i).name);
        continue;
    end

    T = S.DatasetSummary;
    filename = string(files(i).name);

    % Extract animal ID, e.g. DX010
    animal_match = regexp(filename, 'DX\d+', 'match', 'once');

    if isempty(animal_match)
        fprintf('No animal ID: %s\n', files(i).name);
        continue;
    end

    T.Animal = repmat(string(animal_match), height(T), 1);
    T.File = repmat(filename, height(T), 1);

    AllData = [AllData; T];
end

AllData = movevars(AllData, {'Animal','File'}, 'Before', 1);
AllData = AllData(~isnan(AllData.Mean_ms), :);

fprintf('Loaded %d valid file-level rows.\n', height(AllData));

%% Animal-level summary
% Average all files from the same animal first.

groupVars = {'Animal','Amplitude','Condition'};
[G, AnimalGroups] = findgroups(AllData(:,groupVars));

NFiles = splitapply(@numel, AllData.Mean_ms, G);
Mean_ms = splitapply(@mean, AllData.Mean_ms, G);
Median_ms = splitapply(@median, AllData.Mean_ms, G);
SD_ms = splitapply(@std, AllData.Mean_ms, G);

AnimalSummary = AnimalGroups;
AnimalSummary.NFiles = NFiles;
AnimalSummary.Mean_ms = Mean_ms;
AnimalSummary.Median_ms = Median_ms;
AnimalSummary.SD_ms = SD_ms;

AnimalSummary = sortrows(AnimalSummary, ...
    {'Amplitude','Condition','Animal'});

%% Population summary
% Average animal means across animals.

groupVars = {'Amplitude','Condition'};
[G, PopulationGroups] = findgroups(AnimalSummary(:,groupVars));

NAnimals = splitapply(@numel, AnimalSummary.Mean_ms, G);
Mean_ms = splitapply(@mean, AnimalSummary.Mean_ms, G);
Median_ms = splitapply(@median, AnimalSummary.Mean_ms, G);
SD_ms = splitapply(@std, AnimalSummary.Mean_ms, G);
SEM_ms = SD_ms ./ sqrt(NAnimals);

PopulationSummary = PopulationGroups;
PopulationSummary.NAnimals = NAnimals;
PopulationSummary.Mean_ms = Mean_ms;
PopulationSummary.Median_ms = Median_ms;
PopulationSummary.SD_ms = SD_ms;
PopulationSummary.SEM_ms = SEM_ms;

PopulationSummary = sortrows(PopulationSummary, ...
    {'Amplitude','Condition'});

%% Print population summary
fprintf('\nPopulation summary\n');
fprintf('---------------------------------------------------\n');
fprintf('Amp     Condition     Mean      SEM      Animals\n');
fprintf('---------------------------------------------------\n');

for i = 1:height(PopulationSummary)
    fprintf('%2g uA    %-4s        %5.3f     %5.3f       %d\n', ...
        PopulationSummary.Amplitude(i), ...
        PopulationSummary.Condition(i), ...
        PopulationSummary.Mean_ms(i), ...
        PopulationSummary.SEM_ms(i), ...
        PopulationSummary.NAnimals(i));
end

fprintf('---------------------------------------------------\n');

%% Overall recovery time across animals
% First average all amplitude/condition values within each animal.

[G, AnimalOverallGroups] = findgroups(AnimalSummary.Animal);

AnimalMean_ms = splitapply(@mean, ...
    AnimalSummary.Mean_ms, G);

AnimalOverall = table( ...
    AnimalOverallGroups, ...
    AnimalMean_ms, ...
    'VariableNames', {'Animal','Mean_ms'});

% Population mean and SEM across animals
NAnimalsOverall = height(AnimalOverall);
OverallMean_ms = mean(AnimalOverall.Mean_ms);
OverallSD_ms = std(AnimalOverall.Mean_ms);
OverallSEM_ms = OverallSD_ms / sqrt(NAnimalsOverall);

% Range of amplitude x condition population means
ConditionMeanMin_ms = min(PopulationSummary.Mean_ms);
ConditionMeanMax_ms = max(PopulationSummary.Mean_ms);

fprintf('\nOverall recovery time\n');
fprintf('------------------------------------\n');
fprintf('Mean +/- SEM: %.3f +/- %.3f ms\n', ...
    OverallMean_ms, OverallSEM_ms);
fprintf('Animals: %d\n', NAnimalsOverall);
fprintf('Condition mean range: %.3f - %.3f ms\n', ...
    ConditionMeanMin_ms, ConditionMeanMax_ms);
fprintf('------------------------------------\n');

%% Paired Sim - Seq difference at animal level
amps = unique(AnimalSummary.Amplitude);
DifferenceSummary = table();

for ai = 1:numel(amps)
    amp = amps(ai);

    Tsim = AnimalSummary( ...
        AnimalSummary.Amplitude == amp & ...
        AnimalSummary.Condition == "Sim", ...
        {'Animal','Mean_ms'});

    Tseq = AnimalSummary( ...
        AnimalSummary.Amplitude == amp & ...
        AnimalSummary.Condition == "Seq", ...
        {'Animal','Mean_ms'});

    Tsim.Properties.VariableNames{'Mean_ms'} = 'Sim_ms';
    Tseq.Properties.VariableNames{'Mean_ms'} = 'Seq_ms';

    Tpair = innerjoin(Tsim, Tseq, 'Keys','Animal');

    if isempty(Tpair)
        continue;
    end

    diff_ms = Tpair.Sim_ms - Tpair.Seq_ms;

    newRow = table( ...
        amp, ...
        height(Tpair), ...
        mean(Tpair.Sim_ms), ...
        mean(Tpair.Seq_ms), ...
        mean(diff_ms), ...
        std(diff_ms), ...
        std(diff_ms)/sqrt(height(Tpair)), ...
        'VariableNames', ...
        {'Amplitude','NAnimals','SimMean_ms','SeqMean_ms', ...
        'SimMinusSeq_ms','DifferenceSD_ms','DifferenceSEM_ms'});

    DifferenceSummary = [DifferenceSummary; newRow];
end

fprintf('\nPaired Sim - Seq difference\n');
fprintf('---------------------------------------------------------\n');
fprintf('Amp      Sim       Seq      Sim-Seq     SEM      Animals\n');
fprintf('---------------------------------------------------------\n');

for i = 1:height(DifferenceSummary)
    fprintf('%2g uA    %5.3f     %5.3f     %+5.3f     %5.3f       %d\n', ...
        DifferenceSummary.Amplitude(i), ...
        DifferenceSummary.SimMean_ms(i), ...
        DifferenceSummary.SeqMean_ms(i), ...
        DifferenceSummary.SimMinusSeq_ms(i), ...
        DifferenceSummary.DifferenceSEM_ms(i), ...
        DifferenceSummary.NAnimals(i));
end

fprintf('---------------------------------------------------------\n');

%% Overall range
overall_min = min(PopulationSummary.Mean_ms);
overall_max = max(PopulationSummary.Mean_ms);

fprintf('\nPopulation mean recovery range: %.3f - %.3f ms\n', ...
    overall_min, overall_max);

%% Plot
if make_plot
    figure('Color','w');
    hold on;

    conditions = ["Sim","Seq"];

    for ci = 1:numel(conditions)
        idx = PopulationSummary.Condition == conditions(ci);
        T = PopulationSummary(idx,:);
        T = sortrows(T, 'Amplitude');

        errorbar(T.Amplitude, T.Mean_ms, T.SEM_ms, ...
            '-o', 'LineWidth',1.5, 'MarkerSize',7);
    end

    xlabel('Stimulation amplitude (\muA)');
    ylabel('Estimated recovery time (ms)');
    legend({'Simultaneous','Sequential'}, ...
        'Location','best','Box','off');

    ylim([0 4]);
    box off;
end

%% Save
if save_result
    save_file = fullfile(result_folder, ...
        'CrossAnimal_ArtifactEstimate.mat');

    save(save_file, ...
        'AllData', ...
        'AnimalSummary', ...
        'AnimalOverall', ...
        'PopulationSummary', ...
        'DifferenceSummary', ...
        'OverallMean_ms', ...
        'OverallSEM_ms', ...
        'OverallSD_ms', ...
        'NAnimalsOverall', ...
        'ConditionMeanMin_ms', ...
        'ConditionMeanMax_ms');

    fprintf('\nSaved: %s\n', save_file);
end