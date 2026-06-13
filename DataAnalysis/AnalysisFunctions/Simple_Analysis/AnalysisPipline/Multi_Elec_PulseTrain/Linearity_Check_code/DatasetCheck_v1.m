%% ============================================================
%  Step 1: Condition Inspector for Single- and Multi-electrode
%  Pulse-train Datasets
%
%  Purpose:
%    Before doing linearity / prediction analysis, inspect conditions in
%    the single-electrode and two-electrode pulse-train datasets.
%
%  Why:
%    The single-electrode pulse-train dataset is stored separately from the
%    two-electrode pulse-train dataset.
%
%    Therefore, before predicting:
%
%       Predicted Seq = Single A + shifted Single B
%       Predicted AutoSim = Single A + Single B
%
%    we need to confirm which condition corresponds to each stimulation
%    channel, train level, amplitude, and event timing.
%
%  Output:
%    ConditionSummaryTable printed in MATLAB command window.
%
%  Important:
%    This code does NOT load spikes.
%    This code does NOT calculate spike count.
%    This code only inspects StimMeta / StimParams / TrialParams.
% ============================================================

clear all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ====================== USER SETTINGS ======================

% ----- Change these two paths -----
single_data_folder = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Single2';
multi_data_folder  = '/Volumes/MACData/Data/Data_Xia/DX028/Xia_Elec2_Train2';

% Whether to save the summary table.
save_summary_table = true;

% Name for saved summary table.
out_table_name = 'PulseTrain_ConditionInspection_SingleAndMulti.mat';

%% ====================== INSPECT BOTH DATASETS ======================

DatasetFolders = {single_data_folder, multi_data_folder};
DatasetNames   = {'SingleDataset', 'MultiDataset'};

AllRows = {};
row_i = 0;

for ds = 1:numel(DatasetFolders)

    data_folder = DatasetFolders{ds};
    dataset_name = DatasetNames{ds};

    fprintf('\n============================================================\n');
    fprintf('Inspecting dataset: %s\n', dataset_name);
    fprintf('Folder:\n%s\n', data_folder);
    fprintf('============================================================\n');

    if ~isfolder(data_folder)
        error('Dataset folder does not exist:\n%s', data_folder);
    end

    %% ----- Load exp datafile -----
    fileDIR = dir(fullfile(data_folder, '*_exp_datafile_*.mat'));
    assert(~isempty(fileDIR), 'No *_exp_datafile_*.mat found in:\n%s', data_folder);

    exp_file = fullfile(data_folder, fileDIR(1).name);
    S_exp = load(exp_file);

    fprintf('Loaded exp datafile:\n%s\n', exp_file);

    %% ----- Required variables -----
    StimParams        = S_exp.StimParams;
    TrialParams       = S_exp.TrialParams;
    simultaneous_stim = S_exp.simultaneous_stim;
    n_Trials          = S_exp.n_Trials;
    E_MAP             = S_exp.E_MAP;

    if isfield(S_exp, 'StimMeta')
        StimMeta = S_exp.StimMeta;
    else
        error('StimMeta was not found in:\n%s', exp_file);
    end

    if isfield(S_exp, 'trainMode')
        trainMode = S_exp.trainMode;
    else
        trainMode = NaN;
    end

    fprintf('n_Trials          : %d\n', n_Trials);
    fprintf('simultaneous_stim : %d\n', simultaneous_stim);
    fprintf('trainMode         : %g\n', trainMode);
    fprintf('Number of StimMeta conditions: %d\n', numel(StimMeta));

    %% ----- Remove header rows -----
    StimParams_data  = StimParams(2:end,:);
    TrialParams_data = TrialParams(2:end,:);

    expected_rows = n_Trials * simultaneous_stim;

    if size(StimParams_data,1) ~= expected_rows
        warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
            size(StimParams_data,1), expected_rows);
    end

    if size(TrialParams_data,1) ~= expected_rows
        warning('TrialParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d).', ...
            size(TrialParams_data,1), expected_rows);
    end

    %% ----- Trial-level condition IDs -----
    % TrialParams columns:
    %   col 1 = trial number
    %   col 2 = condition ID
    %   col 3 = internal electrode ID

    firstRow_eachTrial = 1:simultaneous_stim:size(TrialParams_data,1);

    conditionID_trial = cell2mat(TrialParams_data(firstRow_eachTrial,2));
    conditionID_trial = conditionID_trial(:);

    %% ----- Count trials per condition -----
    condIDs_all = unique(conditionID_trial(~isnan(conditionID_trial)));
    condIDs_all = condIDs_all(:)';

    %% ----- Loop through all conditions -----
    for cc = 1:numel(condIDs_all)

        condID = condIDs_all(cc);

        if condID < 1 || condID > numel(StimMeta)
            warning('Skipping invalid condition ID %d.', condID);
            continue;
        end

        meta = StimMeta(condID);

        trial_list = find(conditionID_trial == condID);
        nTrials_thisCond = numel(trial_list);

        if isempty(trial_list)
            continue;
        end

        % Use first trial of this condition to inspect randomized StimParams.
        tr_example = trial_list(1);
        rr = (tr_example-1)*simultaneous_stim + (1:simultaneous_stim);

        if max(rr) > size(StimParams_data,1)
            warning('Condition %d example trial rows exceed StimParams size. Skipped.', condID);
            continue;
        end

        %% ----- Basic StimMeta fields -----
        stimSet = getFieldOrNaN(meta, 'StimSetIndex');
        trainLevel = getFieldOrNaN(meta, 'TrainLevel');
        totalLevels = getFieldOrNaN(meta, 'TotalTrainLevels');

        isAutoSim = false;
        if isfield(meta, 'IsAutoSimultaneous')
            isAutoSim = logical(meta.IsAutoSimultaneous);
        end

        isZero = false;
        if isfield(meta, 'IsZeroCurrentControl')
            isZero = logical(meta.IsZeroCurrentControl);
        end

        if isfield(meta, 'EventTimesIncluded_ms')
            eventTimes_ms = meta.EventTimesIncluded_ms;
        else
            eventTimes_ms = [];
        end

        if isfield(meta, 'EventEndTime_ms')
            eventEnd_ms = meta.EventEndTime_ms;
        elseif ~isempty(eventTimes_ms)
            eventEnd_ms = max(eventTimes_ms);
        else
            eventEnd_ms = NaN;
        end

        if isfield(meta, 'PulseCountPerElectrode')
            pulseCountPerElectrode = meta.PulseCountPerElectrode;
        else
            pulseCountPerElectrode = [];
        end

        %% ----- StimParams information from example trial -----
        stimNames_all = StimParams_data(rr,1);
        ptd_us_all    = cell2mat(StimParams_data(rr,6));
        pulseNum_all  = cell2mat(StimParams_data(rr,8));
        period_us_all = cell2mat(StimParams_data(rr,9));
        amp_all       = cell2mat(StimParams_data(rr,16));

        amp_all = amp_all(:)';
        ptd_us_all = ptd_us_all(:)';
        pulseNum_all = pulseNum_all(:)';
        period_us_all = period_us_all(:)';

        activeRows = amp_all > 0;

        activeStimNames = stimNames_all(activeRows);
        activeAmpValues = amp_all(activeRows);
        activePTD_us    = ptd_us_all(activeRows);
        activePulseNum  = pulseNum_all(activeRows);
        activePeriod_us = period_us_all(activeRows);

        if isempty(activeStimNames)
            activeStimNamesText = 'None';
            mappedChannelText = 'None';
            nActiveStimChannels = 0;
            ampValue = 0;
        else
            activeStimNames = unique(activeStimNames, 'stable');
            activeStimNamesText = strjoin(activeStimNames(:)', ',');

            mappedChannelText = buildStimLabelFromStimNames(activeStimNames, E_MAP, isAutoSim);

            nActiveStimChannels = numel(activeStimNames);

            % Usually all active rows have same amplitude.
            % Use max as condition amplitude.
            ampValue = max(activeAmpValues);
        end

        %% ----- Convert numeric arrays to readable text -----
        eventTimesText = numericVectorToText(eventTimes_ms);
        pulseCountText = numericVectorToText(pulseCountPerElectrode);
        activePTDText  = numericVectorToText(activePTD_us ./ 1000); % ms
        activePulseNumText = numericVectorToText(activePulseNum);
        activePeriodText = numericVectorToText(activePeriod_us ./ 1000); % ms

        %% ----- Condition class label -----
        if isZero
            conditionClass = 'ZeroControl';
        elseif nActiveStimChannels == 1
            conditionClass = 'SingleElectrode';
        elseif nActiveStimChannels >= 2 && isAutoSim
            conditionClass = 'MultiElectrode_AutoSim';
        elseif nActiveStimChannels >= 2 && ~isAutoSim
            conditionClass = 'MultiElectrode_Seq';
        else
            conditionClass = 'Unknown';
        end

        %% ----- Add row -----
        row_i = row_i + 1;

        AllRows(row_i,1:20) = { ...
            dataset_name, ...
            data_folder, ...
            condID, ...
            conditionClass, ...
            stimSet, ...
            trainLevel, ...
            totalLevels, ...
            isAutoSim, ...
            isZero, ...
            nActiveStimChannels, ...
            activeStimNamesText, ...
            mappedChannelText, ...
            ampValue, ...
            eventTimesText, ...
            eventEnd_ms, ...
            pulseCountText, ...
            activePTDText, ...
            activePulseNumText, ...
            activePeriodText, ...
            nTrials_thisCond};
    end
end

%% ====================== CREATE TABLE ======================

ConditionSummaryTable = cell2table(AllRows, ...
    'VariableNames', { ...
    'Dataset', ...
    'Folder', ...
    'ConditionID', ...
    'ConditionClass', ...
    'StimSetIndex', ...
    'TrainLevel', ...
    'TotalTrainLevels', ...
    'IsAutoSimultaneous', ...
    'IsZeroCurrentControl', ...
    'NActiveStimChannels', ...
    'ActiveStimChannelNames', ...
    'MappedChannelLabel', ...
    'Amplitude_uA', ...
    'EventTimesIncluded_ms', ...
    'EventEndTime_ms', ...
    'PulseCountPerElectrode', ...
    'ActivePTD_ms', ...
    'ActivePulseNum', ...
    'ActivePulsePeriod_ms', ...
    'NTrials'});

%% ====================== PRINT TABLE ======================

fprintf('\n\n================ FULL CONDITION SUMMARY TABLE ================\n');
disp(ConditionSummaryTable);
fprintf('==============================================================\n');

%% ====================== PRINT SIMPLIFIED MATCHING VIEW ======================
% This view is easier for matching single and multi conditions.

SimpleTable = ConditionSummaryTable(:, { ...
    'Dataset', ...
    'ConditionID', ...
    'ConditionClass', ...
    'StimSetIndex', ...
    'TrainLevel', ...
    'IsAutoSimultaneous', ...
    'NActiveStimChannels', ...
    'MappedChannelLabel', ...
    'Amplitude_uA', ...
    'EventTimesIncluded_ms', ...
    'PulseCountPerElectrode', ...
    'NTrials'});

fprintf('\n\n================ SIMPLIFIED MATCHING VIEW ================\n');
disp(SimpleTable);
fprintf('==========================================================\n');

%% ====================== SAVE ======================

if save_summary_table

    % Save to the multi dataset folder by default.
    out_path = fullfile(multi_data_folder, out_table_name);

    save(out_path, ...
        'ConditionSummaryTable', ...
        'SimpleTable', ...
        'single_data_folder', ...
        'multi_data_folder', ...
        '-v7.3');

    fprintf('\nSaved condition inspection table to:\n%s\n', out_path);
end

fprintf('\nFinished condition inspection.\n');

%% ====================== LOCAL HELPER FUNCTIONS ======================

function val = getFieldOrNaN(S, fieldName)
    if isfield(S, fieldName)
        val = S.(fieldName);
        if isempty(val)
            val = NaN;
        end
    else
        val = NaN;
    end
end

function txt = numericVectorToText(v)
    % Convert a numeric vector to compact text for tables.

    if isempty(v)
        txt = '[]';
        return;
    end

    if iscell(v)
        try
            v = cell2mat(v);
        catch
            txt = '[cell]';
            return;
        end
    end

    if ~isnumeric(v)
        txt = '[non-numeric]';
        return;
    end

    v = v(:)';

    if all(isnan(v))
        txt = '[NaN]';
        return;
    end

    txt = ['[' strtrim(sprintf('%.3g ', v)) ']'];
end

function setLabel = buildStimLabelFromStimNames(stimNames_active, E_MAP, isAutoSim)
    % Build stimulation label from active Intan stimulation channel names.
    %
    % Output examples:
    %   Single:
    %       Ch18
    %
    %   Sequential:
    %       Ch18→Ch22
    %
    %   AutoSim:
    %       Ch18+Ch22

    if isempty(stimNames_active)
        setLabel = 'NoActiveStim';
        return;
    end

    stimNames_active = unique(stimNames_active, 'stable');

    labelParts = cell(1, numel(stimNames_active));

    for i = 1:numel(stimNames_active)

        stimName = stimNames_active{i};

        chNum = convertStimNameUsingEMap(stimName, E_MAP);

        if isnan(chNum)
            labelParts{i} = sprintf('%s', stimName);
        else
            labelParts{i} = sprintf('Ch%d', chNum);
        end
    end

    if numel(labelParts) == 1
        setLabel = labelParts{1};
    elseif isAutoSim
        setLabel = strjoin(labelParts, '+');
    else
        setLabel = strjoin(labelParts, '→');
    end
end

function chNum = convertStimNameUsingEMap(stimName, E_MAP)
    % Convert Intan stimulation label such as 'A-017' to the corresponding
    % channel number using E_MAP.
    %
    % In your E_MAP format:
    %   E_MAP{1,1} is the array/map name.
    %   E_MAP{2,1} is channel 1.
    %   E_MAP{3,1} is channel 2.
    %   ...
    %
    % Therefore:
    %   if E_MAP{row,1} matches stimName,
    %   channel number = row - 1.

    chNum = NaN;

    if isempty(stimName)
        return;
    end

    if isstring(stimName)
        stimName = char(stimName);
    end

    stimName = strtrim(stimName);

    %% ----- Main case: E_MAP is a cell array -----
    if iscell(E_MAP)

        for r = 1:size(E_MAP,1)

            thisName = E_MAP{r,1};

            if isstring(thisName)
                thisName = char(thisName);
            end

            if ischar(thisName)

                thisName = strtrim(thisName);

                if strcmp(thisName, stimName)
                    chNum = r - 1;
                    return;
                end
            end
        end
    end

    %% ----- If E_MAP is a string array -----
    if isstring(E_MAP)

        hit = find(strcmp(strtrim(E_MAP), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    %% ----- If E_MAP is a char matrix -----
    if ischar(E_MAP)

        nameList = cellstr(E_MAP);
        hit = find(strcmp(strtrim(nameList), stimName), 1, 'first');

        if ~isempty(hit)
            chNum = hit - 1;
            return;
        end
    end

    %% ----- Fallback only if E_MAP matching fails -----
    tok = regexp(stimName, '(\d+)', 'tokens', 'once');

    if ~isempty(tok)
        chNum = str2double(tok{1});
        warning('Stim channel %s was not found in E_MAP. Falling back to parsed number %d.', ...
            stimName, chNum);
    else
        warning('Stim channel %s was not found in E_MAP and could not be parsed.', stimName);
    end
end