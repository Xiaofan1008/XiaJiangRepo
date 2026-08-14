%% Stage D3: plot spike count/trial versus stimulation amplitude
%
% PURPOSE
%   This script loads only the Stage D1 ModelData MAT file and plots
%   amplitude-response curves for the five stimulation conditions:
%       A alone, B alone, A+B simultaneous, A->B, and B->A.
%
%   A simple linear A+B prediction is also calculated from the independent
%   A-alone and B-alone trial groups.
%
% IMPORTANT CALCULATION DETAILS
%   1. Spike counts are recalculated directly from relative spike times.
%      Therefore, baseline_win_ms and count_win_ms can be changed without
%      rerunning Stage C1 or Stage D1.
%
%   2. With baseline correction enabled:
%          evoked count = post-window count
%                         - expected baseline count in the post window
%
%   3. The linear prediction is:
%          predicted A+B = mean evoked A + mean evoked B
%
%      A and B trials are independent; they are not artificially paired.
%      The individual-channel linear SEM is propagated from the two
%      independent single-condition SEM values.
%
%   4. The final population figure averages channel means. Its error bars
%      are SEM across responding channels, not SEM across trials or animals.
%
% TRIAL MODES
%   'all'      uses every trial retained in the modelling package.
%   'balanced' uses the fixed balanced subset created by Stage D1.
%
% OUTPUT SAFETY
%   Figures are displayed only. This script does not save or modify files.

clear;
clc;

%% ========================= USER SETTINGS ==============================
model_data_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/stage_a_check_v5/check_results/Luke_ModelPackage/DX014_D4_Pair_A030_A031_ModelData.mat'];

% Choose 'all' or 'balanced'.
Trial_Mode = 'all';

% Select broad curve groups. Any combination can be entered:
%   Curve_Groups_To_Plot = 'single';
%   Curve_Groups_To_Plot = 'simultaneous';
%   Curve_Groups_To_Plot = 'sequential';
%   Curve_Groups_To_Plot = {'single','simultaneous'};
%   Curve_Groups_To_Plot = {'simultaneous','sequential'};
%   Curve_Groups_To_Plot = 'all';
%
% Group definitions:
%   single       = A alone and B alone
%   simultaneous = A+B simultaneous
%   sequential   = A->B and B->A, or their combined curve when
%                  Combine_Sequential_Orders=true
Curve_Groups_To_Plot = 'all';

% OPTIONAL ADVANCED OVERRIDE
% Leave empty to use Curve_Groups_To_Plot. To choose exact conditions,
% enter condition codes or indices here, for example:
%   Conditions_To_Plot = {'A','AB','A_to_B'};
%   Conditions_To_Plot = [1 3 4];
%
% Condition codes:
%   A, B, AB, A_to_B, B_to_A
Conditions_To_Plot = [];

% Empty means every amplitude stored in ModelData.
% Example: Amplitudes_To_Plot = [5 10];
Amplitudes_To_Plot = [];

% Select anatomical/depth-channel numbers. Use 'all' for every exported
% responding channel, or a numeric list such as [1 3 5 8].
Channels_To_Plot = 'all';

% Analysis windows in milliseconds relative to the first pulse.
% Both must lie inside ModelData.metadata.stored_spike_window_ms.
% The interval convention is [start,end): start included, end excluded.
baseline_win_ms = [-50 -5];
count_win_ms = [2 20];

% true  = plot baseline-corrected evoked spike count/trial
% false = plot raw post-window spike count/trial
Use_Baseline_Correction = true;

% Add the prediction calculated from A-alone and B-alone responses.
Plot_Linear_Prediction = true;

% Choose which figure type to display:
%   'individual' = recording-channel panels only
%   'average'    = responding-channel average only
%   'both'       = individual panels and the average
%   'all'        = accepted as an alias for 'both'
Plot_Figure_Mode = 'average';

% false = keep A->B and B->A as two separate curves.
% true  = replace them with one equally weighted, order-averaged
%         Sequential curve. The two means are averaged and their
%         independent SEMs are propagated; trials are not artificially
%         paired and neither stimulation order receives extra weight.
Combine_Sequential_Orders = true;

% Individual recording-channel layout.
channels_per_figure = 4;
Use_Shared_YAxis_Per_Figure = true;

% Keeping zero visible prevents small differences from appearing larger
% merely because the y-axis starts near the smallest response.
Include_Zero_On_YAxis = true;

% Presentation only. No figures are saved.
Close_Existing_Figures = true;
figure_position = [80 70 1300 820];

%% ========================== LOAD DATA =================================
if Close_Existing_Figures
    close all;
end

if ~isfile(model_data_file)
    error('StageD3:MissingModelData', ...
        'ModelData file not found:\n%s',model_data_file);
end

S = load(model_data_file,'ModelData');
if ~isfield(S,'ModelData')
    error('StageD3:InvalidModelData', ...
        'The selected MAT file does not contain ModelData.');
end
D = S.ModelData;

required_top_fields = {'format_name','metadata','amplitudes_uA', ...
    'condition_definitions','channels','conditions','pair_key', ...
    'trial_selection'};
require_fields(D,required_top_fields,'ModelData');
if ~strcmp(D.format_name,'SinglePairStimulationModelData')
    error('StageD3:UnexpectedFormat', ...
        'The selected MAT file is not a Stage D1 single-pair package.');
end
if numel(D.conditions) ~= 5
    error('StageD3:ConditionCount','The package must contain five conditions.');
end

required_channel_variables = {'ChannelIndex','DepthChannel','IsStimA','IsStimB'};
for k = 1:numel(required_channel_variables)
    if ~ismember(required_channel_variables{k},D.channels.Properties.VariableNames)
        error('StageD3:MissingChannelVariable', ...
            'ModelData.channels is missing %s.',required_channel_variables{k});
    end
end

%% ====================== RESOLVE USER CHOICES ==========================
trial_mode = lower(strtrim(char(string(Trial_Mode))));
if ~ismember(trial_mode,{'all','balanced'})
    error('StageD3:InvalidTrialMode', ...
        'Trial_Mode must be either all or balanced.');
end
[plot_individual,plot_average,plot_figure_mode] = ...
    resolve_plot_figure_mode(Plot_Figure_Mode);

if isempty(Conditions_To_Plot)
    condition_indices = resolve_curve_groups(Curve_Groups_To_Plot,D);
    condition_selection_source = 'Curve_Groups_To_Plot';
else
    condition_indices = resolve_conditions(Conditions_To_Plot,D);
    condition_selection_source = 'Conditions_To_Plot override';
end
[amplitude_indices,selected_amplitudes] = ...
    resolve_amplitudes(Amplitudes_To_Plot,D.amplitudes_uA);
[selected_channel_columns,selected_depth_channels] = ...
    resolve_channels(Channels_To_Plot,D.channels);

validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(count_win_ms,'count_win_ms');
stored_window = double(D.metadata.stored_spike_window_ms(:).');
if baseline_win_ms(1) < stored_window(1) || ...
        baseline_win_ms(2) > stored_window(2) || ...
        count_win_ms(1) < stored_window(1) || ...
        count_win_ms(2) > stored_window(2)
    error('StageD3:WindowOutsideStoredData', ...
        ['Baseline and count windows must lie inside the stored spike ' ...
         'window [%g %g] ms.'],stored_window);
end
if ~isscalar(channels_per_figure) || channels_per_figure < 1 || ...
        mod(channels_per_figure,1) ~= 0
    error('StageD3:InvalidPageSize', ...
        'channels_per_figure must be a positive integer.');
end

%% ====================== DISPLAY DEFINITIONS ===========================
condition_labels = {D.conditions.label};
condition_colors = [0.10 0.45 0.80; ... % A
                    0.85 0.33 0.10; ... % B
                    0.15 0.65 0.30; ... % A+B
                    0.60 0.25 0.70; ... % A->B
                    0.15 0.65 0.70];    % B->A
linear_color = [0.10 0.10 0.10];

if Use_Baseline_Correction
    metric_name = 'Baseline-corrected spike count / trial';
    metric_short = 'baseline-corrected';
else
    metric_name = 'Raw spike count / trial';
    metric_short = 'raw';
end

%% =================== RECALCULATE TRIAL METRICS ========================
% These arrays store one mean and SEM for every exported channel,
% condition, and amplitude. Calculating all exported channels also allows
% the population figure to remain independent of Channels_To_Plot.
nChannelsAll = height(D.channels);
nConditionsAll = numel(D.conditions);
nAmplitudesAll = numel(D.amplitudes_uA);

channel_mean = nan(nChannelsAll,nConditionsAll,nAmplitudesAll);
channel_sem = nan(nChannelsAll,nConditionsAll,nAmplitudesAll);
trial_count_used = zeros(nConditionsAll,nAmplitudesAll);

% Store A and B trial values separately because they are required for the
% independent linear prediction even if A or B is hidden from the plot.
A_values = cell(1,nAmplitudesAll);
B_values = cell(1,nAmplitudesAll);
A_raw_components = cell(1,nAmplitudesAll);
B_raw_components = cell(1,nAmplitudesAll);

for ic = 1:nConditionsAll
    for ia = 1:nAmplitudesAll
        block = D.conditions(ic).amplitude(ia);
        validate_amplitude_block(block,nChannelsAll);
        trial_indices = select_trial_indices(block,trial_mode);
        trial_count_used(ic,ia) = numel(trial_indices);

        metrics = calculate_trial_metrics( ...
            block.spike_times_ms(trial_indices,:), ...
            baseline_win_ms,count_win_ms);

        if Use_Baseline_Correction
            values = metrics.evoked_count;
        else
            values = metrics.post_count;
        end

        channel_mean(:,ic,ia) = mean(values,1).';
        channel_sem(:,ic,ia) = ...
            (std(values,0,1)/sqrt(size(values,1))).';

        if ic == 1
            A_values{ia} = metrics.evoked_count;
            % For a raw-count linear prediction, A and B each contribute
            % half of one baseline copy. Their sum therefore contains one
            % expected baseline rather than two.
            A_raw_components{ia} = metrics.post_count- ...
                0.5*metrics.expected_baseline_post_count;
        elseif ic == 2
            B_values{ia} = metrics.evoked_count;
            B_raw_components{ia} = metrics.post_count- ...
                0.5*metrics.expected_baseline_post_count;
        end
    end
end

%% ===================== LINEAR A+B PREDICTION ==========================
linear_mean = nan(nChannelsAll,nAmplitudesAll);
linear_sem = nan(nChannelsAll,nAmplitudesAll);

for ia = 1:nAmplitudesAll
    if isempty(A_values{ia}) || isempty(B_values{ia})
        error('StageD3:MissingSingleConditions', ...
            'A-alone and B-alone data are required for the linear prediction.');
    end

    if Use_Baseline_Correction
        values_A = A_values{ia};
        values_B = B_values{ia};
    else
        values_A = A_raw_components{ia};
        values_B = B_raw_components{ia};
    end

    mean_A = mean(values_A,1);
    mean_B = mean(values_B,1);
    sem_A = std(values_A,0,1)/sqrt(size(values_A,1));
    sem_B = std(values_B,0,1)/sqrt(size(values_B,1));

    linear_mean(:,ia) = (mean_A+mean_B).';
    linear_sem(:,ia) = sqrt(sem_A.^2+sem_B.^2).';
end

%% ================== PREPARE DISPLAYED CONDITIONS ======================
% channel_mean and channel_sem above always retain the original five
% conditions. The arrays below optionally replace A->B and B->A with one
% equally weighted order-averaged Sequential condition for visualization.
[plot_channel_mean,plot_channel_sem,plot_condition_codes, ...
    plot_condition_labels,plot_condition_colors,report_condition_indices] = ...
    prepare_display_conditions(channel_mean,channel_sem,condition_indices, ...
    condition_labels,condition_colors,Combine_Sequential_Orders);
nPlotConditions = numel(plot_condition_codes);

%% ======================= COMMAND-WINDOW SUMMARY =======================
nChannelFigures = ceil(numel(selected_channel_columns)/channels_per_figure);
nFiguresExpected = double(plot_individual)*nChannelFigures+double(plot_average);

fprintf('\n============================================================\n');
fprintf('STAGE D3: SPIKE COUNT/TRIAL VS AMPLITUDE\n');
fprintf('============================================================\n');
fprintf('ModelData file: %s\n',model_data_file);
fprintf('Pair: %s\n',D.pair_key);
fprintf('Trial mode: %s\n',upper(trial_mode));
fprintf('Metric: %s\n',metric_name);
fprintf('Baseline window: [%g %g) ms\n',baseline_win_ms);
fprintf('Count window: [%g %g) ms\n',count_win_ms);
fprintf('Conditions plotted: %s\n', ...
    strjoin(plot_condition_codes,', '));
fprintf('Condition selection source: %s\n',condition_selection_source);
fprintf('Sequential orders combined: %s\n', ...
    logical_text(Combine_Sequential_Orders));
fprintf('Linear A+B prediction: %s\n',logical_text(Plot_Linear_Prediction));
fprintf('Amplitudes: %s uA\n',num2str(selected_amplitudes));
fprintf('Figure mode: %s\n',upper(plot_figure_mode));
if plot_individual
    fprintf('Individual depth channels: %s\n',num2str(selected_depth_channels));
end
fprintf('Responding-channel average: %s (%d exported channels)\n', ...
    logical_text(plot_average),nChannelsAll);
fprintf('Expected figures: %d\n',nFiguresExpected);
fprintf('Figures will be displayed only. Files written: NONE\n');

if isfield(D,'qc') && isfield(D.qc,'single_trial_review_complete') && ...
        ~D.qc.single_trial_review_complete && ...
        (Plot_Linear_Prediction || any(ismember(report_condition_indices,[1 2])))
    fprintf(['NOTICE: The displayed single-condition data and/or linear ' ...
        'prediction use trials exported before single-trial review was completed.\n']);
end

% Print the number of selected trials because the all-trial mode can have
% different counts across conditions.
SelectedTrialCounts = array2table( ...
    trial_count_used(report_condition_indices,amplitude_indices).', ...
    'VariableNames',{D.conditions(report_condition_indices).code}, ...
    'RowNames',cellstr(compose('%g_uA',selected_amplitudes(:))));
fprintf('\n[TRIAL COUNTS USED]\n%s\n',evalc('disp(SelectedTrialCounts)'));

%% ================= INDIVIDUAL-CHANNEL FIGURES =========================
figure_counter = 0;
nPages = ceil(numel(selected_channel_columns)/channels_per_figure);
nTileColumns = ceil(sqrt(channels_per_figure));
nTileRows = ceil(channels_per_figure/nTileColumns);

if plot_individual
  for page = 1:nPages
    first_local = (page-1)*channels_per_figure+1;
    last_local = min(page*channels_per_figure,numel(selected_channel_columns));
    page_local_indices = first_local:last_local;
    page_channel_columns = selected_channel_columns(page_local_indices);

    figure_counter = figure_counter+1;
    figure('Color','w', ...
        'Name',sprintf('D3 amplitude curves page %d',page), ...
        'Position',figure_position);
    tl = tiledlayout(nTileRows,nTileColumns, ...
        'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s | %s | %s trials | [%g, %g) ms', ...
        D.pair_key,metric_short,upper(trial_mode),count_win_ms), ...
        'Interpreter','none','FontWeight','bold');

    page_ymin = inf;
    page_ymax = -inf;
    if Use_Shared_YAxis_Per_Figure
        values_for_limits = [];
        for jc = page_channel_columns
            observed = reshape(plot_channel_mean(jc,:,amplitude_indices), ...
                nPlotConditions,numel(amplitude_indices));
            observed_sem = reshape(plot_channel_sem(jc,:,amplitude_indices), ...
                nPlotConditions,numel(amplitude_indices));
            values_for_limits = [values_for_limits; ... %#ok<AGROW>
                observed-observed_sem; observed+observed_sem];
            if Plot_Linear_Prediction
                lm = linear_mean(jc,amplitude_indices);
                ls = linear_sem(jc,amplitude_indices);
                values_for_limits = [values_for_limits; lm-ls; lm+ls]; %#ok<AGROW>
            end
        end
        [page_ymin,page_ymax] = padded_limits(values_for_limits, ...
            Include_Zero_On_YAxis);
    end

    for q = 1:numel(page_channel_columns)
        jc = page_channel_columns(q);
        ax = nexttile(tl); hold(ax,'on');
        curve_handles = gobjects(0);
        legend_labels = {};

        observed = reshape(plot_channel_mean(jc,:,amplitude_indices), ...
            nPlotConditions,numel(amplitude_indices));
        observed_sem = reshape(plot_channel_sem(jc,:,amplitude_indices), ...
            nPlotConditions,numel(amplitude_indices));

        for icLocal = 1:nPlotConditions
            h = errorbar(ax,selected_amplitudes,observed(icLocal,:), ...
                observed_sem(icLocal,:),'-o', ...
                'Color',plot_condition_colors(icLocal,:), ...
                'MarkerFaceColor',plot_condition_colors(icLocal,:), ...
                'LineWidth',1.5,'MarkerSize',5,'CapSize',6);
            curve_handles(end+1) = h; %#ok<SAGROW>
            legend_labels{end+1} = plot_condition_labels{icLocal}; %#ok<SAGROW>
        end

        if Plot_Linear_Prediction
            h = errorbar(ax,selected_amplitudes, ...
                linear_mean(jc,amplitude_indices), ...
                linear_sem(jc,amplitude_indices),'--d', ...
                'Color',linear_color,'MarkerFaceColor','w', ...
                'LineWidth',1.8,'MarkerSize',6,'CapSize',7);
            curve_handles(end+1) = h;
            legend_labels{end+1} = 'Linear A+B';
        end

        if Use_Baseline_Correction
            yline(ax,0,':','Color',[0.45 0.45 0.45], ...
                'HandleVisibility','off');
        end
        if Use_Shared_YAxis_Per_Figure
            ylim(ax,[page_ymin page_ymax]);
        else
            local_limit_values = [observed-observed_sem; ...
                observed+observed_sem];
            if Plot_Linear_Prediction
                local_limit_values = [local_limit_values; ...
                    linear_mean(jc,amplitude_indices)- ...
                        linear_sem(jc,amplitude_indices); ...
                    linear_mean(jc,amplitude_indices)+ ...
                        linear_sem(jc,amplitude_indices)];
            end
            [local_min,local_max] = padded_limits(local_limit_values, ...
                Include_Zero_On_YAxis);
            ylim(ax,[local_min local_max]);
        end

        depth_ch = D.channels.DepthChannel(jc);
        channel_tag = stimulation_channel_tag(D.channels(jc,:));
        title(ax,sprintf('Recording Ch %d%s',depth_ch,channel_tag), ...
            'Interpreter','none','FontWeight','bold');
        xlabel(ax,'Amplitude (uA)');
        ylabel(ax,metric_name);
        xticks(ax,selected_amplitudes);
        grid(ax,'on'); box(ax,'off');

        if q == 1
            legend(ax,curve_handles,legend_labels, ...
                'Location','best','Box','off');
        end
    end
    fprintf('Figure %d/%d: individual channels %s\n', ...
        figure_counter,nFiguresExpected, ...
        num2str(D.channels.DepthChannel(page_channel_columns).'));
    drawnow;
  end
end

%% =============== RESPONDING-CHANNEL AVERAGE FIGURE ===================
if plot_average
    % Each channel contributes one condition mean. This prevents channels
    % with more trials from receiving more weight in the channel average.
    population_mean = squeeze(mean(plot_channel_mean(:,:, ...
        amplitude_indices),1,'omitnan'));
    population_n = squeeze(sum(~isnan(plot_channel_mean(:,:, ...
        amplitude_indices)),1));
    population_sem = squeeze(std(plot_channel_mean(:,:, ...
        amplitude_indices),0,1,'omitnan'))./sqrt(population_n);

    population_linear_mean = mean(linear_mean(:,amplitude_indices),1,'omitnan');
    population_linear_n = sum(~isnan(linear_mean(:,amplitude_indices)),1);
    population_linear_sem = std(linear_mean(:,amplitude_indices),0,1,'omitnan') ./ ...
        sqrt(population_linear_n);

    % Preserve condition x amplitude orientation when only one condition
    % or one amplitude is selected.
    population_mean = reshape(population_mean, ...
        nPlotConditions,numel(amplitude_indices));
    population_sem = reshape(population_sem, ...
        nPlotConditions,numel(amplitude_indices));

    figure_counter = figure_counter+1;
    fig_population = figure('Color','w', ...
        'Name','D3 responding-channel average','Position',figure_position);
    ax = axes(fig_population); hold(ax,'on');
    curve_handles = gobjects(0);
    legend_labels = {};

    for icLocal = 1:nPlotConditions
        h = errorbar(ax,selected_amplitudes,population_mean(icLocal,:), ...
            population_sem(icLocal,:),'-o', ...
            'Color',plot_condition_colors(icLocal,:), ...
            'MarkerFaceColor',plot_condition_colors(icLocal,:), ...
            'LineWidth',1.8,'MarkerSize',6,'CapSize',7);
        curve_handles(end+1) = h; %#ok<SAGROW>
        legend_labels{end+1} = plot_condition_labels{icLocal}; %#ok<SAGROW>
    end

    if Plot_Linear_Prediction
        h = errorbar(ax,selected_amplitudes,population_linear_mean, ...
            population_linear_sem,'--d','Color',linear_color, ...
            'MarkerFaceColor','w','LineWidth',2.1,'MarkerSize',7,'CapSize',8);
        curve_handles(end+1) = h;
        legend_labels{end+1} = 'Linear A+B';
    end

    if Use_Baseline_Correction
        yline(ax,0,':','Color',[0.45 0.45 0.45], ...
            'HandleVisibility','off');
    end
    population_limit_values = [population_mean-population_sem; ...
        population_mean+population_sem];
    if Plot_Linear_Prediction
        population_limit_values = [population_limit_values; ...
            population_linear_mean-population_linear_sem; ...
            population_linear_mean+population_linear_sem];
    end
    [population_ymin,population_ymax] = padded_limits( ...
        population_limit_values,Include_Zero_On_YAxis);
    ylim(ax,[population_ymin population_ymax]);
    xlabel(ax,'Amplitude (uA)','FontWeight','bold');
    ylabel(ax,['Mean ' lower(metric_name)],'FontWeight','bold');
    title(ax,sprintf('%s | Average across %d responding channels | %s trials', ...
        D.pair_key,nChannelsAll,upper(trial_mode)), ...
        'Interpreter','none','FontWeight','bold');
    subtitle(ax,'Error bars are SEM across channels, not across animals.');
    xticks(ax,selected_amplitudes);
    legend(ax,curve_handles,legend_labels,'Location','bestoutside','Box','off');
    grid(ax,'on'); box(ax,'off');

    fprintf('Figure %d/%d: average across %d responding channels\n', ...
        figure_counter,nFiguresExpected,nChannelsAll);
    drawnow;
end

fprintf('\n============================================================\n');
fprintf('STAGE D3 COMPLETE\n');
fprintf('Figures displayed: %d\n',figure_counter);
fprintf('Files written: NONE\n');
fprintf('============================================================\n\n');

%% =========================== FUNCTIONS ================================
function [plot_individual,plot_average,mode] = ...
        resolve_plot_figure_mode(request)
mode = lower(strtrim(char(string(request))));
switch mode
    case 'individual'
        plot_individual = true;
        plot_average = false;
    case 'average'
        plot_individual = false;
        plot_average = true;
    case {'both','all'}
        plot_individual = true;
        plot_average = true;
        mode = 'both';
    otherwise
        error('StageD3:InvalidFigureMode', ...
            'Plot_Figure_Mode must be individual, average, both, or all.');
end
end

function [means_out,sems_out,codes_out,labels_out,colors_out,report_indices] = ...
        prepare_display_conditions(means_in,sems_in,requested_indices, ...
        original_labels,original_colors,combine_sequential)
% Build the curves shown to the user without changing the original data.
% When requested, conditions 4 and 5 are replaced by one equally weighted
% order average. SEM propagation assumes their condition means are
% independent:
%       SEM_combined = 0.5*sqrt(SEM_AtoB^2 + SEM_BtoA^2)

if ~isscalar(combine_sequential) || ...
        ~(islogical(combine_sequential) || isnumeric(combine_sequential))
    error('StageD3:InvalidSequentialOption', ...
        'Combine_Sequential_Orders must be one true/false value.');
end
combine_sequential = logical(combine_sequential);

nChannels = size(means_in,1);
nAmplitudes = size(means_in,3);
means_out = nan(nChannels,0,nAmplitudes);
sems_out = nan(nChannels,0,nAmplitudes);
codes_out = {};
labels_out = {};
colors_out = zeros(0,3);
report_indices = [];
sequential_added = false;

original_codes = {'A','B','AB','A_to_B','B_to_A'};
combined_color = [0.15 0.55 0.65];

for k = 1:numel(requested_indices)
    ic = requested_indices(k);
    is_sequential_order = ic == 4 || ic == 5;

    if combine_sequential && is_sequential_order
        if sequential_added
            continue;
        end
        means_out(:,end+1,:) = 0.5*(means_in(:,4,:)+means_in(:,5,:));
        sems_out(:,end+1,:) = 0.5*sqrt(sems_in(:,4,:).^2+sems_in(:,5,:).^2);
        codes_out{end+1} = 'Sequential'; %#ok<AGROW>
        labels_out{end+1} = 'Sequential (order average)'; %#ok<AGROW>
        colors_out(end+1,:) = combined_color; %#ok<AGROW>
        report_indices = [report_indices 4 5]; %#ok<AGROW>
        sequential_added = true;
    else
        means_out(:,end+1,:) = means_in(:,ic,:);
        sems_out(:,end+1,:) = sems_in(:,ic,:);
        codes_out{end+1} = original_codes{ic}; %#ok<AGROW>
        labels_out{end+1} = original_labels{ic}; %#ok<AGROW>
        colors_out(end+1,:) = original_colors(ic,:); %#ok<AGROW>
        report_indices(end+1) = ic; %#ok<AGROW>
    end
end

report_indices = unique(report_indices,'stable');
if isempty(codes_out)
    error('StageD3:NoDisplayConditions','No conditions remain to plot.');
end
end

function M = calculate_trial_metrics(spike_cells,baseline_window,count_window)
% Return trial x channel matrices calculated from relative spike vectors.
[nTrials,nChannels] = size(spike_cells);
baseline_count = zeros(nTrials,nChannels);
post_count = zeros(nTrials,nChannels);

for it = 1:nTrials
    for jc = 1:nChannels
        tt = double(spike_cells{it,jc}(:));
        baseline_count(it,jc) = sum(tt >= baseline_window(1) & ...
            tt < baseline_window(2));
        post_count(it,jc) = sum(tt >= count_window(1) & ...
            tt < count_window(2));
    end
end

baseline_duration_ms = diff(baseline_window);
count_duration_ms = diff(count_window);
expected_baseline_post_count = baseline_count* ...
    (count_duration_ms/baseline_duration_ms);

M = struct('baseline_count',baseline_count, ...
    'post_count',post_count, ...
    'expected_baseline_post_count',expected_baseline_post_count, ...
    'evoked_count',post_count-expected_baseline_post_count);
end

function indices = resolve_curve_groups(request,D)
% Translate user-friendly curve groups into the five stored conditions.
% The package order is expected to be A, B, AB, A_to_B, B_to_A.
if numel(D.conditions) ~= 5
    error('StageD3:CurveGroupConditionCount', ...
        'Curve-group selection requires the standard five conditions.');
end

if ischar(request) || (isstring(request) && isscalar(request))
    groups = {char(string(request))};
elseif isstring(request) || iscell(request)
    groups = cellstr(string(request));
else
    error('StageD3:InvalidCurveGroups', ...
        'Curve_Groups_To_Plot must be text or a cell array of group names.');
end
if isempty(groups)
    error('StageD3:EmptyCurveGroups','No curve groups were selected.');
end

indices = [];
for k = 1:numel(groups)
    group = lower(strtrim(groups{k}));
    switch group
        case 'all'
            indices = 1:5;
            return;
        case {'single','singles','single_electrode'}
            indices = [indices 1 2]; %#ok<AGROW>
        case {'simultaneous','sim'}
            indices = [indices 3]; %#ok<AGROW>
        case {'sequential','seq'}
            indices = [indices 4 5]; %#ok<AGROW>
        otherwise
            error('StageD3:UnknownCurveGroup', ...
                ['Unknown curve group: %s. Use single, simultaneous, ' ...
                 'sequential, or all.'],groups{k});
    end
end
indices = unique(indices,'stable');
if isempty(indices)
    error('StageD3:EmptyCurveGroups','No curve groups were selected.');
end
end

function indices = resolve_conditions(request,D)
codes = {D.conditions.code};
if ischar(request) || (isstring(request) && isscalar(request))
    value = char(string(request));
    if strcmpi(strtrim(value),'all')
        indices = 1:numel(codes);
    else
        indices = match_condition_codes({value},codes);
    end
elseif isstring(request) || iscell(request)
    requested_codes = cellstr(string(request));
    indices = match_condition_codes(requested_codes,codes);
elseif isnumeric(request)
    indices = unique(double(request(:).'),'stable');
    if isempty(indices) || any(indices < 1) || any(indices > numel(codes)) || ...
            any(mod(indices,1) ~= 0)
        error('StageD3:InvalidConditionIndices', ...
            'Numeric condition indices must be integers from 1 to %d.',numel(codes));
    end
else
    error('StageD3:InvalidConditions','Invalid Conditions_To_Plot setting.');
end
end

function indices = match_condition_codes(requested,codes)
indices = zeros(1,numel(requested));
for k = 1:numel(requested)
    match = find(strcmpi(strtrim(requested{k}),codes),1);
    if isempty(match)
        error('StageD3:UnknownCondition','Unknown condition code: %s',requested{k});
    end
    indices(k) = match;
end
indices = unique(indices,'stable');
end

function [indices,amps] = resolve_amplitudes(request,saved_amps)
saved_amps = double(saved_amps(:).');
if isempty(request)
    indices = 1:numel(saved_amps);
else
    requested = unique(double(request(:).'),'stable');
    indices = zeros(size(requested));
    for k = 1:numel(requested)
        match = find(abs(saved_amps-requested(k)) <= 1e-6,1);
        if isempty(match)
            error('StageD3:UnavailableAmplitude', ...
                'Amplitude %g uA is not present in ModelData.',requested(k));
        end
        indices(k) = match;
    end
end
amps = saved_amps(indices);
end

function [columns,depth_channels] = resolve_channels(request,ChannelTable)
saved_depth = double(ChannelTable.DepthChannel(:).');
if ischar(request) || (isstring(request) && isscalar(request))
    if strcmpi(strtrim(char(string(request))),'all')
        columns = 1:numel(saved_depth);
    else
        error('StageD3:InvalidChannelText', ...
            'Text Channels_To_Plot must be all. Otherwise use depth-channel numbers.');
    end
elseif isnumeric(request)
    requested = unique(double(request(:).'),'stable');
    if isempty(requested) || any(mod(requested,1) ~= 0)
        error('StageD3:InvalidChannels','Depth-channel selection is invalid.');
    end
    [found,columns] = ismember(requested,saved_depth);
    if any(~found)
        error('StageD3:UnavailableChannel', ...
            'Depth channel(s) not present in ModelData: %s', ...
            num2str(requested(~found)));
    end
else
    error('StageD3:InvalidChannels','Invalid Channels_To_Plot setting.');
end
depth_channels = saved_depth(columns);
end

function validate_amplitude_block(B,nChannels)
required = {'n_trials_all','n_trials_balanced','all_trial_indices', ...
    'balanced_trial_indices','spike_times_ms'};
require_fields(B,required,'condition amplitude block');
if ~iscell(B.spike_times_ms) || size(B.spike_times_ms,1) ~= B.n_trials_all || ...
        size(B.spike_times_ms,2) ~= nChannels
    error('StageD3:SpikeShape','The exported spike_times_ms shape is invalid.');
end
end

function indices = select_trial_indices(B,mode)
switch mode
    case 'all'
        indices = double(B.all_trial_indices(:));
    case 'balanced'
        indices = double(B.balanced_trial_indices(:));
end
if isempty(indices) || any(indices < 1) || any(indices > B.n_trials_all) || ...
        any(mod(indices,1) ~= 0) || numel(unique(indices)) ~= numel(indices)
    error('StageD3:InvalidTrialIndices', ...
        'The selected trial indices in ModelData are invalid.');
end
end

function [lower,upper] = padded_limits(values,include_zero)
values = values(isfinite(values));
if isempty(values)
    lower = -1;
    upper = 1;
    return;
end
lower = min(values);
upper = max(values);
if include_zero
    lower = min(lower,0);
    upper = max(upper,0);
end
span = upper-lower;
if span <= 0
    span = max(1,abs(upper));
end
padding = 0.10*span;
lower = lower-padding;
upper = upper+padding;
end

function tag = stimulation_channel_tag(channel_row)
if channel_row.IsStimA && channel_row.IsStimB
    tag = ' [STIM A/B]';
elseif channel_row.IsStimA
    tag = ' [STIM A]';
elseif channel_row.IsStimB
    tag = ' [STIM B]';
else
    tag = '';
end
end

function validate_window(value,name)
if ~isnumeric(value) || numel(value) ~= 2 || any(~isfinite(value)) || ...
        value(2) <= value(1)
    error('StageD3:InvalidWindow','%s must be [start end], with end>start.',name);
end
end

function require_fields(S,names,label)
for k = 1:numel(names)
    if ~isfield(S,names{k})
        error('StageD3:MissingField','%s.%s is missing.',label,names{k});
    end
end
end

function value = logical_text(tf)
if tf, value = 'true'; else, value = 'false'; end
end
