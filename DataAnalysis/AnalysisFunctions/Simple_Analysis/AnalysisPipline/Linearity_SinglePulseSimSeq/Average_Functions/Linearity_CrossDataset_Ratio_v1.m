%% ========================================================================
%  CROSS-DATASET E1 LINEARITY SUMMARY
%
%  Creates two figures from saved Linearity1Saved MAT files:
%    1. Simultaneous and sequential linearity ratio versus amplitude
%    2. Sequential benefit versus amplitude
%
%  Definitions
%    Sim ratio       = matched simultaneous / linear reference
%    Seq ratio       = sequential / linear reference
%    Seq benefit     = Seq ratio - Sim ratio
%
%  Averaging
%    - Ratios are calculated separately for every saved pair/order branch.
%    - Valid branches are averaged within each dataset first.
%    - Dataset values are then averaged with equal dataset weight.
%    - Sim and Seq always use the same valid branches.
%    - Missing amplitudes remain NaN; no interpolation is performed.
% ========================================================================

clear;

%% ============================= USER SETTINGS ============================

Results_Folder = ...
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq';

% Empty loads every file matching Result_File_Pattern in Results_Folder.
% Alternatively, enter full MAT-file paths, for example:
% Result_Files = {
%     '/path/DX001_dataset_Linearity1_v1.mat'
%     '/path/DX002_dataset_Linearity1_v1.mat'
% };
Result_Files = {
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX005_Xia_Exp1_Sim_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX009_Xia_Exp1_Sim5_Linearity1_v1.mat';
    
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim7_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX010_Xia_Exp1_Sim8_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim7_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim8_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX011_Xia_Exp1_Sim9_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX012_Xia_Exp1_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX012_Xia_Exp1_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX012_Xia_Exp1_Sim6_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim7_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX013_Xia_Exp1_Seq_Sim8_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX014_Xia_Seq_Sim6_Linearity1_v1.mat';

    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim1_Linearity1_v1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim3_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim4_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim5_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim6_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX015_Xia_Seq_Sim7_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX016_Xia_Exp1_Seq_Full_1_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX016_Xia_Exp1_Seq_Full_2_Linearity1_v1.mat';
    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX016_Xia_Exp1_Seq_Full_3_Linearity1_v1.mat';

    '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim1_Linearity1_v1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim2_Linearity1_v1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim3_Linearity1_v1.mat';
    % '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Linearity_SinglePulse_SimSeq/DX018_Xia_Exp1_Sim4_Linearity1_v1.mat';
    
};
% Result_File_Pattern = '*_Linearity1_v1.mat';

% Empty uses the union of all positive amplitudes found in the saved files.
% Missing amplitudes remain NaN and do not contribute to that amplitude.
Amplitudes_To_Plot = [];
Include_Zero_Amplitude = false;
Amplitude_Tolerance_uA = 1e-6;

% Ratios with a small or non-positive linear prediction are unstable.
Minimum_Linear_Response = 0.1;  % baseline-corrected spikes/trial

Show_Individual_Datasets = true;
Show_Error_Bars = true;

%% ============================== FIND FILES ===============================

if isempty(Result_Files)
    if ~isfolder(Results_Folder)
        error('CrossDataset:FolderNotFound', ...
            'Results folder not found: %s',Results_Folder);
    end
    file_info = dir(fullfile(Results_Folder,Result_File_Pattern));
    file_info = file_info(~startsWith({file_info.name},'._'));
    Result_Files = arrayfun(@(x) fullfile(x.folder,x.name), ...
        file_info,'UniformOutput',false);
else
    Result_Files = cellstr(string(Result_Files(:)));
end

Result_Files = unique(Result_Files,'stable');
if isempty(Result_Files)
    error('CrossDataset:NoResultFiles', ...
        'No result files matched %s in %s.', ...
        Result_File_Pattern,Results_Folder);
end

%% ============================== LOAD FILES ===============================

nDatasets = numel(Result_Files);
SavedData = cell(nDatasets,1);
DatasetLabels = strings(nDatasets,1);
all_amplitudes = [];

for id = 1:nDatasets
    file_path = Result_Files{id};
    if ~isfile(file_path)
        error('CrossDataset:FileNotFound','Result file not found: %s',file_path);
    end

    L = load(file_path,'Linearity1Saved');
    if ~isfield(L,'Linearity1Saved')
        error('CrossDataset:MissingSavedVariable', ...
            'Linearity1Saved is missing from %s.',file_path);
    end
    Saved = L.Linearity1Saved;
    if ~isfield(Saved,'branch_results') || isempty(Saved.branch_results)
        error('CrossDataset:MissingBranchResults', ...
            'No branch_results were found in %s.',file_path);
    end

    for ib = 1:numel(Saved.branch_results)
        P = Saved.branch_results(ib);
        if ~isfield(P,'average_summary') || ~istable(P.average_summary)
            error('CrossDataset:MissingAverageSummary', ...
                'A branch in %s does not contain average_summary.',file_path);
        end
        required_columns = {'Amplitude_uA','NUsed','Linear', ...
            'SimMatched','Sequential'};
        if ~all(ismember(required_columns, ...
                P.average_summary.Properties.VariableNames))
            error('CrossDataset:InvalidAverageSummary', ...
                'average_summary has missing columns in %s.',file_path);
        end
        all_amplitudes = [all_amplitudes; ...
            double(P.average_summary.Amplitude_uA(:))]; %#ok<AGROW>
    end

    SavedData{id} = Saved;
    DatasetLabels(id) = make_dataset_label(Saved,file_path);
end

all_amplitudes = unique(all_amplitudes(isfinite(all_amplitudes))).';
if ~Include_Zero_Amplitude
    all_amplitudes(abs(all_amplitudes) < Amplitude_Tolerance_uA) = [];
end

if isempty(Amplitudes_To_Plot)
    Amplitudes = sort(all_amplitudes);
else
    Amplitudes = unique(double(Amplitudes_To_Plot(:).'),'stable');
    if ~Include_Zero_Amplitude
        Amplitudes(abs(Amplitudes) < Amplitude_Tolerance_uA) = [];
    end
end
if isempty(Amplitudes)
    error('CrossDataset:NoAmplitudes','No amplitudes remain for analysis.');
end

%% ======================= DATASET-LEVEL RATIOS ============================

nAmp = numel(Amplitudes);
DatasetSimRatio = nan(nDatasets,nAmp);
DatasetSeqRatio = nan(nDatasets,nAmp);
DatasetSeqBenefit = nan(nDatasets,nAmp);
NBranchesUsed = zeros(nDatasets,nAmp);

for id = 1:nDatasets
    BranchResults = SavedData{id}.branch_results;

    for ia = 1:nAmp
        amp = Amplitudes(ia);
        branch_sim_ratio = [];
        branch_seq_ratio = [];

        for ib = 1:numel(BranchResults)
            T = BranchResults(ib).average_summary;
            row = find(abs(double(T.Amplitude_uA)-amp) < ...
                Amplitude_Tolerance_uA,1);
            if isempty(row)
                continue;
            end

            linear_value = double(T.Linear(row));
            sim_value = double(T.SimMatched(row));
            seq_value = double(T.Sequential(row));
            n_used = double(T.NUsed(row));

            % Require a complete paired comparison from this branch.
            if n_used <= 0 || ~isfinite(linear_value) || ...
                    ~isfinite(sim_value) || ~isfinite(seq_value) || ...
                    linear_value <= Minimum_Linear_Response
                continue;
            end

            branch_sim_ratio(end+1) = sim_value/linear_value; %#ok<SAGROW>
            branch_seq_ratio(end+1) = seq_value/linear_value; %#ok<SAGROW>
        end

        if isempty(branch_sim_ratio)
            continue;
        end

        % Equal branch/pair weight within this dataset.
        DatasetSimRatio(id,ia) = mean(branch_sim_ratio);
        DatasetSeqRatio(id,ia) = mean(branch_seq_ratio);
        DatasetSeqBenefit(id,ia) = mean(branch_seq_ratio-branch_sim_ratio);
        NBranchesUsed(id,ia) = numel(branch_sim_ratio);
    end
end

%% ============================ SUMMARIZE =================================

MeanSimRatio = nan(nAmp,1);
SEMSimRatio = nan(nAmp,1);
MeanSeqRatio = nan(nAmp,1);
SEMSeqRatio = nan(nAmp,1);
MeanSeqBenefit = nan(nAmp,1);
SEMSeqBenefit = nan(nAmp,1);
NDatasetsUsed = zeros(nAmp,1);

for ia = 1:nAmp
    paired_valid = isfinite(DatasetSimRatio(:,ia)) & ...
        isfinite(DatasetSeqRatio(:,ia));

    [MeanSimRatio(ia),SEMSimRatio(ia),n_sim] = ...
        mean_sem(DatasetSimRatio(paired_valid,ia));
    [MeanSeqRatio(ia),SEMSeqRatio(ia),n_seq] = ...
        mean_sem(DatasetSeqRatio(paired_valid,ia));
    [MeanSeqBenefit(ia),SEMSeqBenefit(ia),n_benefit] = ...
        mean_sem(DatasetSeqBenefit(paired_valid,ia));

    if ~(n_sim == n_seq && n_seq == n_benefit)
        error('CrossDataset:InternalPairingError', ...
            'Sim/Seq dataset counts became inconsistent at %.3g uA.', ...
            Amplitudes(ia));
    end
    NDatasetsUsed(ia) = n_sim;
end

SummaryTable = table(Amplitudes(:),NDatasetsUsed,MeanSimRatio, ...
    SEMSimRatio,MeanSeqRatio,SEMSeqRatio,MeanSeqBenefit,SEMSeqBenefit, ...
    'VariableNames',{'Amplitude_uA','NDatasets','MeanSimRatio', ...
    'SEMSimRatio','MeanSeqRatio','SEMSeqRatio', ...
    'MeanSeqBenefit','SEMSeqBenefit'});

fprintf('\nCross-dataset E1 linearity ratios\n');
fprintf('Datasets loaded: %d\n',nDatasets);
fprintf('Valid ratio rule: Linear > %.3g spikes/trial\n\n', ...
    Minimum_Linear_Response);
disp(SummaryTable);

%% ============================== FIGURE 1 ================================

sim_color = [0 0.35 0.85];
seq_color = [0.85 0.33 0.10];
light_sim = 0.72*[1 1 1] + 0.28*sim_color;
light_seq = 0.72*[1 1 1] + 0.28*seq_color;

figure('Color','w','Position',[100 100 780 560]);
ax = axes();
hold(ax,'on');
yline(ax,1,'--','Color',[0.25 0.25 0.25], ...
    'LineWidth',1.2,'HandleVisibility','off');

if Show_Individual_Datasets
    for id = 1:nDatasets
        plot(ax,Amplitudes,DatasetSimRatio(id,:),'-o', ...
            'Color',light_sim,'MarkerFaceColor',light_sim, ...
            'LineWidth',0.8,'MarkerSize',4,'HandleVisibility','off');
        plot(ax,Amplitudes,DatasetSeqRatio(id,:),'-o', ...
            'Color',light_seq,'MarkerFaceColor',light_seq, ...
            'LineWidth',0.8,'MarkerSize',4,'HandleVisibility','off');
    end
end

plot_summary_curve(ax,Amplitudes,MeanSimRatio,SEMSimRatio, ...
    sim_color,'Simultaneous',Show_Error_Bars);
plot_summary_curve(ax,Amplitudes,MeanSeqRatio,SEMSeqRatio, ...
    seq_color,'Sequential',Show_Error_Bars);

xlabel(ax,'Amplitude (uA)');
ylabel(ax,'Observed / linear prediction');
title(ax,'Cross-dataset linearity ratio','FontWeight','bold');
legend(ax,'Location','best','Box','off');
box(ax,'off');
axis square;
set(ax,'XTick',Amplitudes);
add_n_labels(ax,Amplitudes,NDatasetsUsed);

%% ============================== FIGURE 2 ================================

benefit_color = [0.45 0.15 0.65];
light_benefit = 0.75*[1 1 1] + 0.25*benefit_color;

figure('Color','w','Position',[150 130 780 560]);
ax = axes();
hold(ax,'on');
yline(ax,0,'--','Color',[0.25 0.25 0.25], ...
    'LineWidth',1.2,'HandleVisibility','off');

if Show_Individual_Datasets
    for id = 1:nDatasets
        plot(ax,Amplitudes,DatasetSeqBenefit(id,:),'-o', ...
            'Color',light_benefit,'MarkerFaceColor',light_benefit, ...
            'LineWidth',0.8,'MarkerSize',4,'HandleVisibility','off');
    end
end

plot_summary_curve(ax,Amplitudes,MeanSeqBenefit,SEMSeqBenefit, ...
    benefit_color,'Sequential - simultaneous',Show_Error_Bars);

xlabel(ax,'Amplitude (uA)');
ylabel(ax,'Sequential ratio - simultaneous ratio');
title(ax,'Sequential benefit across datasets','FontWeight','bold');
legend(ax,'Location','best','Box','off');
box(ax,'off');
set(ax,'XTick',Amplitudes);
axis square;
add_n_labels(ax,Amplitudes,NDatasetsUsed);

%% ============================== FUNCTIONS ================================

function label = make_dataset_label(Saved,file_path)
[~,file_name] = fileparts(file_path);
animal = "";
dataset = "";
if isfield(Saved,'animal_id') && ~isempty(Saved.animal_id)
    animal = string(Saved.animal_id);
end
if isfield(Saved,'dataset_id') && ~isempty(Saved.dataset_id)
    dataset = string(Saved.dataset_id);
end
if strlength(animal) > 0 && strlength(dataset) > 0
    label = animal+" | "+dataset;
elseif strlength(dataset) > 0
    label = dataset;
else
    label = string(file_name);
end
end

function [avg,se,n] = mean_sem(values)
values = double(values(:));
values = values(isfinite(values));
n = numel(values);
if n == 0
    avg = NaN;
    se = NaN;
    return;
end
avg = mean(values);
if n > 1
    se = std(values,0)/sqrt(n);
else
    se = NaN;
end
end

function plot_summary_curve(ax,x,y,se,color,label,show_errors)
x = double(x(:).');
y = double(y(:).');
se = double(se(:).');
valid = isfinite(x) & isfinite(y);
if ~any(valid)
    return;
end
if show_errors
    se_plot = se;
    se_plot(~isfinite(se_plot)) = 0;
    errorbar(ax,x(valid),y(valid),se_plot(valid),'-o', ...
        'Color',color,'MarkerFaceColor',color,'LineWidth',2, ...
        'MarkerSize',6,'CapSize',7,'DisplayName',label);
else
    plot(ax,x(valid),y(valid),'-o','Color',color, ...
        'MarkerFaceColor',color,'LineWidth',2,'MarkerSize',6, ...
        'DisplayName',label);
end
end

function add_n_labels(ax,x,n_values)
valid_n = n_values > 0;
if ~any(valid_n)
    return;
end
yl = ylim(ax);
span = yl(2)-yl(1);
if ~isfinite(span) || span <= 0
    return;
end
for k = find(valid_n(:).')
    text(ax,x(k),yl(1)+0.025*span,sprintf('n=%d',n_values(k)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',8,'Color',[0.3 0.3 0.3]);
end
end
