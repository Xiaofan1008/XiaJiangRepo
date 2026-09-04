%% Cross-dataset artifact recovery-time analysis
clear;
clc;

%% Settings
result_folder = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Artifact_Time_Est';
save_result = false;
make_plot = true;

%% Find result files
files = dir(fullfile(result_folder, '*_ArtifactEst.mat'));
assert(~isempty(files), 'No artifact blank estimate files found.');

fprintf('Found %d datasets.\n', numel(files));

%% Combine dataset-level results
AllData = table();

for i = 1:numel(files)
    file_path = fullfile(files(i).folder, files(i).name);
    S = load(file_path, 'DatasetSummary');

    if ~isfield(S, 'DatasetSummary') || isempty(S.DatasetSummary)
        fprintf('Skipped: %s\n', files(i).name);
        continue;
    end

    T = S.DatasetSummary;
    dataset_name = erase(files(i).name, '_ArtifactBlankEstimate.mat');
    T.Dataset = repmat(string(dataset_name), height(T), 1);

    AllData = [AllData; T];
    fprintf('Loaded: %s\n', dataset_name);
end

AllData = movevars(AllData, 'Dataset', 'Before', 1);

%% Remove invalid rows
valid = ~isnan(AllData.Mean_ms);
AllData = AllData(valid,:);

%% Population summary: mean across datasets
groupVars = {'Amplitude','Condition'};
[G, PopulationGroups] = findgroups(AllData(:,groupVars));

NDatasets = splitapply(@numel, AllData.Mean_ms, G);
Mean_ms = splitapply(@mean, AllData.Mean_ms, G);
Median_ms = splitapply(@median, AllData.Mean_ms, G);
SD_ms = splitapply(@std, AllData.Mean_ms, G);
SEM_ms = SD_ms ./ sqrt(NDatasets);

PopulationSummary = PopulationGroups;
PopulationSummary.NDatasets = NDatasets;
PopulationSummary.Mean_ms = Mean_ms;
PopulationSummary.Median_ms = Median_ms;
PopulationSummary.SD_ms = SD_ms;
PopulationSummary.SEM_ms = SEM_ms;

PopulationSummary = sortrows(PopulationSummary, {'Amplitude','Condition'});

%% Print population summary
fprintf('\nCross-dataset summary\n');
fprintf('--------------------------------------------------\n');
fprintf('Amp     Condition     Mean      SEM       N\n');
fprintf('--------------------------------------------------\n');

for i = 1:height(PopulationSummary)
    fprintf('%2g uA    %-4s        %5.3f     %5.3f      %d\n', ...
        PopulationSummary.Amplitude(i), ...
        PopulationSummary.Condition(i), ...
        PopulationSummary.Mean_ms(i), ...
        PopulationSummary.SEM_ms(i), ...
        PopulationSummary.NDatasets(i));
end

fprintf('--------------------------------------------------\n');

%% Sim - Seq difference for each amplitude
amps = unique(AllData.Amplitude);
DifferenceSummary = table();

for ai = 1:numel(amps)
    amp = amps(ai);

    Tsim = AllData(AllData.Amplitude == amp & AllData.Condition == "Sim", ...
        {'Dataset','Mean_ms'});

    Tseq = AllData(AllData.Amplitude == amp & AllData.Condition == "Seq", ...
        {'Dataset','Mean_ms'});

    Tsim.Properties.VariableNames{'Mean_ms'} = 'Sim_ms';
    Tseq.Properties.VariableNames{'Mean_ms'} = 'Seq_ms';

    Tpair = innerjoin(Tsim, Tseq, 'Keys','Dataset');

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
        std(diff_ms) / sqrt(height(Tpair)), ...
        'VariableNames', ...
        {'Amplitude','NDatasets','SimMean_ms','SeqMean_ms', ...
        'SimMinusSeq_ms','DifferenceSEM_ms'});

    DifferenceSummary = [DifferenceSummary; newRow];
end

fprintf('\nSim - Seq difference\n');
fprintf('------------------------------------------------------\n');
fprintf('Amp      Sim       Seq      Sim-Seq     SEM      N\n');
fprintf('------------------------------------------------------\n');

for i = 1:height(DifferenceSummary)
    fprintf('%2g uA    %5.3f     %5.3f     %+5.3f     %5.3f    %d\n', ...
        DifferenceSummary.Amplitude(i), ...
        DifferenceSummary.SimMean_ms(i), ...
        DifferenceSummary.SeqMean_ms(i), ...
        DifferenceSummary.SimMinusSeq_ms(i), ...
        DifferenceSummary.DifferenceSEM_ms(i), ...
        DifferenceSummary.NDatasets(i));
end

fprintf('------------------------------------------------------\n');

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
        'CrossDataset_ArtifactBlankEstimate.mat');

    save(save_file, ...
        'AllData', ...
        'PopulationSummary', ...
        'DifferenceSummary', ...
        'overall_min', ...
        'overall_max');

    fprintf('\nSaved: %s\n', save_file);
end