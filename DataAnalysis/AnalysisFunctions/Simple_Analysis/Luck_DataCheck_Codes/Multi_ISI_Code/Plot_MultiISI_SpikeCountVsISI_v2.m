%% MULTI-ISI SPIKE COUNT VERSUS ISI
%
% PURPOSE
%   Plot mean spike count per trial as a function of inter-stimulation
%   interval (ISI/PTD) from one exported MultiISI ModelData MAT file.
%
% DATA SOURCE
%   This script reads only the file created by Run_MultiISI_Export.m.
%   It does not need the original experiment folders, trigger files,
%   sp_corr files, bad-trial files, or external analysis functions.
%
% DEFAULT PRESENTATION
%   - A+B simultaneous is plotted at ISI = 0 ms.
%   - The simultaneous point is connected to each sequential curve.
%   - A->B and B->A are kept separate unless combining is requested.
%   - The linear prediction shifts the second single-pulse response by
%     each PTD before adding it to the first response.
%   - Population error bars are SEM across recording channels.
%
% IMPORTANT STATISTICAL RULE
%   Each recording channel is first averaged across its own clean trials.
%   Population statistics are then calculated across channel means. Thus,
%   a channel with more retained trials does not receive extra weight.
%
% This script only opens figures. It does not save results or change data.

clear;
close all;

%% ========================= USER SETTINGS ==============================

% Full path to one exported *_MultiISI_ModelData.mat file.
% No folder-selection window is used.
model_data_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/Multi_ISI_Output/DX020_XiaISISimSeq2_MultiISI_ModelData.mat'];

% Needed only to translate hardware names such as A-021 into anatomical
% channel numbers such as Ch 27. This calculation is self-contained.
%   0 = single-shank rigid
%   1 = single-shank flexible
%   2 = four-shank flexible
Electrode_Type = 2;

% Which recording-channel results should be displayed?
%   'average'    = one population-average figure per amplitude
%   'individual' = individual-channel figures only
%   'all'        = both average and individual figures
Plot_Mode = 'average';

% Empty means use every amplitude represented in sequential conditions.
% Conditions that do not contain a selected amplitude are simply skipped.
Amplitudes_To_Plot = [];

% Empty means use every observed positive ISI/PTD.
PTDs_To_Plot = [];

% These are exported ChannelIndex values in ModelData.channels.ChannelIndex.
% Empty means use every exported responding channel.
Channels_To_Plot = [];

% Available line names:
%   'A'       : response to stimulation electrode A alone
%   'B'       : response to stimulation electrode B alone
%   'AB'      : simultaneous A+B response at ISI 0 ms
%   'A_to_B'  : A followed by B
%   'B_to_A'  : B followed by A
%   'Linear'  : PTD-dependent shifted single-A + single-B prediction
%
% A and B alone are PTD-independent and appear as horizontal lines. The
% shifted linear prediction can vary with PTD because the same paired-data
% response window is applied after shifting the second response.
Lines_To_Plot = {'AB','A_to_B','B_to_A','Linear'};

% false: retain A->B and B->A as separate curves.
% true:  average the available order means within each channel. If only one
%        order exists at a PTD, that available order is retained.
Combine_Sequential_Orders = false;

% When true, prepend the simultaneous A+B point at 0 ms to every displayed
% sequential curve. Change to false later if a separate point is preferred.
Connect_Simultaneous_To_Sequential = true;

% Spike-count windows in ms relative to the first stimulation pulse.
% Recommended convention: exclude the immediate 0-2 ms period and retain
% enough time for the response to a second pulse delivered at 20 ms.
% To reproduce the older analysis windows, use [0 40] and [-50 -10].
Response_Window_ms = [2 45];
Baseline_Window_ms = [-50 -5];

% true  = subtract the expected baseline count from every trial
% false = plot the raw response-window spike count
Subtract_Baseline = true;

% Applied after baseline subtraction, separately for every observed trial
% AND every single-electrode trial used by the shifted prediction:
% true  = max(0, corrected count), matching the previous analysis
% false = preserve positive and negative corrected values (recommended)
Force_Nonnegative = true;

% true uses a separate manually corrected response population for A->B and
% B->A. This requires a ModelData v1.1 file from the revised exporter.
% false uses the complete exported responding-channel union for both orders.
Use_Order_Specific_Response_Channels = true;

% Show SEM across channels on population-average figures.
Show_Population_SEM = true;

% Maximum number of individual recording channels placed in one figure.
% Additional channels are continued in a new figure.
Channels_Per_Figure = 9;

% Figure appearance.
figure_position_average    = [100 100 760 560];
figure_position_individual = [60 60 1500 900];
line_width = 2.0;
marker_size = 7;

%% ======================= LOAD AND VALIDATE ============================

if ~isfile(model_data_file)
    error('MultiISIPlot:FileNotFound', ...
        'ModelData file was not found:\n%s',model_data_file);
end

Loaded = load(model_data_file,'ModelData');
if ~isfield(Loaded,'ModelData')
    error('MultiISIPlot:MissingModelData', ...
        'The selected MAT file does not contain ModelData.');
end
ModelData = Loaded.ModelData;

required_top_fields = {'format_name','conditions','channels','stimulation'};
for k = 1:numel(required_top_fields)
    if ~isfield(ModelData,required_top_fields{k})
        error('MultiISIPlot:InvalidModelData', ...
            'ModelData.%s is missing.',required_top_fields{k});
    end
end
if ~strcmp(ModelData.format_name,'MultiISIStimulationModelData')
    error('MultiISIPlot:WrongFormat', ...
        'The selected file is not a multi-ISI modelling export.');
end
if isempty(ModelData.conditions) || height(ModelData.channels) == 0
    error('MultiISIPlot:EmptyData','The exported ModelData is empty.');
end
if ~ismember(Electrode_Type,[0 1 2])
    error('MultiISIPlot:ElectrodeType','Electrode_Type must be 0, 1, or 2.');
end
[stim_channel_A,stim_channel_B] = stimulation_depth_channels( ...
    ModelData,Electrode_Type);
ModelData.display_stimulation_A = sprintf('Ch %d',stim_channel_A);
ModelData.display_stimulation_B = sprintf('Ch %d',stim_channel_B);

valid_modes = {'average','individual','all'};
Plot_Mode = lower(strtrim(Plot_Mode));
if ~any(strcmp(Plot_Mode,valid_modes))
    error('MultiISIPlot:PlotMode', ...
        'Plot_Mode must be ''average'', ''individual'', or ''all''.');
end
if ~isscalar(Channels_Per_Figure) || Channels_Per_Figure < 1 || ...
        fix(Channels_Per_Figure) ~= Channels_Per_Figure
    error('MultiISIPlot:PageSize', ...
        'Channels_Per_Figure must be a positive integer.');
end

allowed_lines = {'A','B','AB','A_to_B','B_to_A','Linear'};
Lines_To_Plot = cellstr(string(Lines_To_Plot));
for k = 1:numel(Lines_To_Plot)
    match = find(strcmpi(Lines_To_Plot{k},allowed_lines),1);
    if isempty(match)
        error('MultiISIPlot:LineSelection', ...
            'Unknown Lines_To_Plot entry: %s',Lines_To_Plot{k});
    end
    Lines_To_Plot{k} = allowed_lines{match};
end
Lines_To_Plot = unique(Lines_To_Plot,'stable');

validate_count_window(Response_Window_ms,'Response_Window_ms');
validate_count_window(Baseline_Window_ms,'Baseline_Window_ms');
if ~islogical(Subtract_Baseline) || ~isscalar(Subtract_Baseline) || ...
        ~islogical(Force_Nonnegative) || ~isscalar(Force_Nonnegative)
    error('MultiISIPlot:CountOptions', ...
        'All baseline and nonnegative settings must be true or false.');
end
if ~isfield(ModelData,'metadata') || ...
        ~isfield(ModelData.metadata,'stored_spike_window_ms')
    error('MultiISIPlot:StoredWindow', ...
        'ModelData.metadata.stored_spike_window_ms is missing.');
end
stored_window = double(ModelData.metadata.stored_spike_window_ms(:).');
if Response_Window_ms(1)<stored_window(1) || ...
        Response_Window_ms(2)>stored_window(2) || ...
        Baseline_Window_ms(1)<stored_window(1) || ...
        Baseline_Window_ms(2)>stored_window(2)
    error('MultiISIPlot:WindowOutsideExport', ...
        ['The response and baseline windows must remain inside the stored ' ...
         'spike window [%g,%g] ms.'],stored_window(1),stored_window(2));
end

CountSettings = struct('response_window_ms',double(Response_Window_ms(:).'), ...
    'baseline_window_ms',double(Baseline_Window_ms(:).'), ...
    'subtract_baseline',Subtract_Baseline, ...
    'force_nonnegative',Force_Nonnegative);

all_channel_ids = double(ModelData.channels.ChannelIndex(:).');
if isempty(Channels_To_Plot)
    selected_channel_ids = all_channel_ids;
else
    Channels_To_Plot = unique(double(Channels_To_Plot(:).'),'stable');
    selected_channel_ids = Channels_To_Plot(ismember(Channels_To_Plot,all_channel_ids));
end
if isempty(selected_channel_ids)
    error('MultiISIPlot:NoChannels', ...
        'None of Channels_To_Plot exists in ModelData.channels.ChannelIndex.');
end
[~,selected_channel_rows] = ismember(selected_channel_ids,all_channel_ids);

if Use_Order_Specific_Response_Channels
    required_columns = {'Responsive_A_to_B','Responsive_B_to_A'};
    for k = 1:numel(required_columns)
        if ~ismember(required_columns{k},ModelData.channels.Properties.VariableNames)
            error('MultiISIPlot:ResponseMasksMissing', ...
                ['Order-specific response masks are missing. Rerun the revised ' ...
                 'Run_MultiISI_Export.m, then select the new ModelData file.']);
        end
    end
    response_mask_A_to_B = logical( ...
        ModelData.channels.Responsive_A_to_B(selected_channel_rows));
    response_mask_B_to_A = logical( ...
        ModelData.channels.Responsive_B_to_A(selected_channel_rows));
else
    response_mask_A_to_B = true(numel(selected_channel_rows),1);
    response_mask_B_to_A = true(numel(selected_channel_rows),1);
end

% Determine available sequential PTDs and amplitudes from the actual data.
is_seq = strcmp({ModelData.conditions.stimulation_type},'sequential');
available_ptds = unique([ModelData.conditions(is_seq).PTD_ms]);
available_ptds = sort(available_ptds(isfinite(available_ptds) & available_ptds>0));
if isempty(PTDs_To_Plot)
    selected_ptds = available_ptds;
else
    selected_ptds = requested_values_present(PTDs_To_Plot,available_ptds);
end
if isempty(selected_ptds)
    error('MultiISIPlot:NoPTDs','No selected positive PTDs are available.');
end

% The second single-pulse response is shifted by PTD. Equivalently, its
% unshifted spikes are counted inside Response_Window_ms-PTD. Confirm that
% every required shifted window was retained in the exported spike data.
largest_ptd = max(selected_ptds);
largest_shifted_window = CountSettings.response_window_ms-largest_ptd;
if largest_shifted_window(1)<stored_window(1) || ...
        largest_shifted_window(2)>stored_window(2)
    error('MultiISIPlot:ShiftedWindowOutsideExport', ...
        ['The largest PTD requires single-pulse spikes in [%g,%g] ms, ' ...
         'outside the stored range [%g,%g] ms.'],largest_shifted_window,stored_window);
end

available_amps = sequential_amplitudes(ModelData.conditions);
if isempty(Amplitudes_To_Plot)
    selected_amps = available_amps;
else
    selected_amps = requested_values_present(Amplitudes_To_Plot,available_amps);
end
if isempty(selected_amps)
    error('MultiISIPlot:NoAmplitudes', ...
        'No selected amplitudes are available in sequential conditions.');
end

plot_average = any(strcmp(Plot_Mode,{'average','all'}));
plot_individual = any(strcmp(Plot_Mode,{'individual','all'}));

%% =========================== PLOTTING ================================

for amp = selected_amps

    % Each matrix contains one row per selected recording channel and one
    % column per PTD. Entries are channel-level trial means.
    A_values  = condition_channel_means(ModelData.conditions,'A',0,amp, ...
        selected_channel_rows,CountSettings);
    B_values  = condition_channel_means(ModelData.conditions,'B',0,amp, ...
        selected_channel_rows,CountSettings);
    AB_values = condition_channel_means(ModelData.conditions,'AB',0,amp, ...
        selected_channel_rows,CountSettings);

    nCh = numel(selected_channel_rows);
    AtoB = nan(nCh,numel(selected_ptds));
    BtoA = nan(nCh,numel(selected_ptds));
    for ip = 1:numel(selected_ptds)
        ptd = selected_ptds(ip);
        AtoB(:,ip) = condition_channel_means( ...
            ModelData.conditions,'A_to_B',ptd,amp,selected_channel_rows,CountSettings);
        BtoA(:,ip) = condition_channel_means( ...
            ModelData.conditions,'B_to_A',ptd,amp,selected_channel_rows,CountSettings);
    end

    % Apply each order's manually corrected response population. The same
    % order-specific population is used for its sequential curve, its
    % simultaneous anchor, and its linear prediction.
    AtoB(~response_mask_A_to_B,:) = NaN;
    BtoA(~response_mask_B_to_A,:) = NaN;
    AB_A_to_B = apply_population_mask(AB_values,response_mask_A_to_B);
    AB_B_to_A = apply_population_mask(AB_values,response_mask_B_to_A);

    % Shift the second single-pulse response by every PTD before addition.
    % PTD=0 is included so the predicted simultaneous value is visible.
    prediction_ptds = [0 selected_ptds];
    Linear_A_to_B = shifted_linear_prediction(ModelData.conditions, ...
        'A','B',amp,prediction_ptds,selected_channel_rows,CountSettings);
    Linear_B_to_A = shifted_linear_prediction(ModelData.conditions, ...
        'B','A',amp,prediction_ptds,selected_channel_rows,CountSettings);
    Linear_A_to_B = apply_population_mask(Linear_A_to_B,response_mask_A_to_B);
    Linear_B_to_A = apply_population_mask(Linear_B_to_A,response_mask_B_to_A);

    % Combining orders gives equal weight to the available order means.
    % It does not pool all trials and therefore does not require balancing.
    Combined = mean_available_orders(AtoB,BtoA);
    Combined_AB = mean_available_orders(AB_A_to_B,AB_B_to_A);
    Combined_Linear = mean_available_orders(Linear_A_to_B,Linear_B_to_A);

    if plot_average
        plot_population_figure(amp,selected_ptds,prediction_ptds,A_values,B_values, ...
            AB_A_to_B,AB_B_to_A,Combined_AB,AtoB,BtoA,Combined, ...
            Linear_A_to_B,Linear_B_to_A,Combined_Linear,Lines_To_Plot, ...
            Combine_Sequential_Orders,Connect_Simultaneous_To_Sequential, ...
            Show_Population_SEM,CountSettings,ModelData,figure_position_average, ...
            line_width,marker_size);
    end

    if plot_individual
        nPages = ceil(nCh/Channels_Per_Figure);
        for page = 1:nPages
            first_ch = (page-1)*Channels_Per_Figure+1;
            last_ch = min(page*Channels_Per_Figure,nCh);
            page_rows = first_ch:last_ch;
            plot_individual_page(amp,selected_ptds,prediction_ptds,A_values,B_values, ...
                AB_A_to_B,AB_B_to_A,Combined_AB,AtoB,BtoA,Combined, ...
                Linear_A_to_B,Linear_B_to_A,Combined_Linear, ...
                page_rows,selected_channel_rows, ...
                Lines_To_Plot,Combine_Sequential_Orders, ...
                Connect_Simultaneous_To_Sequential,CountSettings,ModelData, ...
                figure_position_individual,line_width,marker_size,page,nPages);
        end
    end
end

%% ========================= LOCAL FUNCTIONS ============================

function validate_count_window(window,name)
if ~isnumeric(window) || numel(window)~=2 || ...
        any(~isfinite(window)) || window(2)<=window(1)
    error('MultiISIPlot:InvalidWindow','%s must be [start end], with end > start.',name);
end
end

function values = requested_values_present(requested,available)
% Retain requested numeric values that exist, using a small tolerance.
requested = unique(double(requested(:).'),'stable');
available = double(available(:).');
values = [];
for value = requested
    idx = find(abs(available-value)<1e-6,1);
    if ~isempty(idx), values(end+1) = available(idx); end %#ok<AGROW>
end
end

function amps = sequential_amplitudes(Conditions)
% Return the union of amplitudes represented by sequential conditions.
amps = [];
for ic = 1:numel(Conditions)
    if ~strcmp(Conditions(ic).stimulation_type,'sequential'), continue; end
    if ~isempty(Conditions(ic).amplitude)
        amps = [amps [Conditions(ic).amplitude.amplitude_uA]]; %#ok<AGROW>
    end
end
amps = sort(unique(amps(isfinite(amps))));
end

function means = condition_channel_means(Conditions,code,ptd,amp,channel_rows,settings)
% Recalculate per-trial counts from saved relative spikes, then average each
% channel. This permits user-selected response and baseline windows.
means = nan(numel(channel_rows),1);
condition_index = [];
for ic = 1:numel(Conditions)
    if strcmp(Conditions(ic).code,code) && ...
            abs(Conditions(ic).PTD_ms-ptd)<1e-6
        condition_index = ic;
        break;
    end
end
if isempty(condition_index), return; end

C = Conditions(condition_index);
amp_index = [];
for ia = 1:numel(C.amplitude)
    if abs(C.amplitude(ia).amplitude_uA-amp)<1e-6
        amp_index = ia;
        break;
    end
end
if isempty(amp_index), return; end

channels = C.amplitude(amp_index).channel;
for k = 1:numel(channel_rows)
    row = channel_rows(k);
    if row>numel(channels) || ~isfield(channels(row),'spike_times_ms')
        continue;
    end
    means(k) = mean_count_in_window(channels(row).spike_times_ms, ...
        settings.response_window_ms,settings.baseline_window_ms, ...
        settings.subtract_baseline,settings.force_nonnegative);
end
end

function prediction = shifted_linear_prediction(Conditions,first_code, ...
        second_code,amp,ptds,channel_rows,settings)
% Integrating a shifted single-pulse PSTH over response window W is
% equivalent to counting the unshifted second response inside W-PTD.
% The same Force_Nonnegative choice used for observed trials is also applied
% to every single-electrode trial before the two contributions are added.
first_channels = condition_amplitude_channels(Conditions,first_code,amp);
second_channels = condition_amplitude_channels(Conditions,second_code,amp);
prediction = nan(numel(channel_rows),numel(ptds));
if isempty(first_channels) || isempty(second_channels), return; end

for k = 1:numel(channel_rows)
    row = channel_rows(k);
    if row>numel(first_channels) || row>numel(second_channels), continue; end
    first_value = mean_count_in_window(first_channels(row).spike_times_ms, ...
        settings.response_window_ms,settings.baseline_window_ms, ...
        settings.subtract_baseline,settings.force_nonnegative);
    for ip = 1:numel(ptds)
        shifted_second_window = settings.response_window_ms-ptds(ip);
        second_value = mean_count_in_window(second_channels(row).spike_times_ms, ...
            shifted_second_window,settings.baseline_window_ms, ...
            settings.subtract_baseline,settings.force_nonnegative);
        value = first_value+second_value;
        prediction(k,ip) = value;
    end
end
end

function channels = condition_amplitude_channels(Conditions,code,amp)
channels = [];
for ic = 1:numel(Conditions)
    if ~strcmp(Conditions(ic).code,code), continue; end
    for ia = 1:numel(Conditions(ic).amplitude)
        if abs(Conditions(ic).amplitude(ia).amplitude_uA-amp)<1e-6
            channels = Conditions(ic).amplitude(ia).channel;
            return;
        end
    end
end
end

function value = mean_count_in_window(spike_trials,response_window, ...
        baseline_window,subtract_baseline,force_nonnegative)
trial_values = nan(numel(spike_trials),1);
response_duration = diff(response_window);
baseline_duration = diff(baseline_window);
for it = 1:numel(spike_trials)
    spikes = double(spike_trials{it}(:));
    response_count = sum(spikes>=response_window(1) & spikes<response_window(2));
    if subtract_baseline
        baseline_count = sum(spikes>=baseline_window(1) & spikes<baseline_window(2));
        trial_values(it) = response_count- ...
            baseline_count*(response_duration/baseline_duration);
        if force_nonnegative
            trial_values(it) = max(0,trial_values(it));
        end
    else
        trial_values(it) = response_count;
    end
end
trial_values = trial_values(isfinite(trial_values));
if isempty(trial_values), value=NaN;
else, value=mean(trial_values); end
end

function values = apply_population_mask(values,mask)
values(~mask,:) = NaN;
end

function combined = mean_available_orders(order1,order2)
% Average order means without requiring both orders to be present.
combined = nan(size(order1));
only1 = isfinite(order1) & ~isfinite(order2);
only2 = ~isfinite(order1) & isfinite(order2);
both = isfinite(order1) & isfinite(order2);
combined(only1) = order1(only1);
combined(only2) = order2(only2);
combined(both) = (order1(both)+order2(both))/2;
end

function plot_population_figure(amp,ptds,prediction_ptds,A,B, ...
        AB_AtoB,AB_BtoA,AB_Combined, ...
        AtoB,BtoA,Combined,Linear_AtoB,Linear_BtoA,Linear_Combined, ...
        lines,combine_orders,connect_ab,show_sem,settings,M,position,lw,ms)
% Plot the population mean. SEM is calculated across channel means.
name = sprintf('Average spike count vs ISI | %.4g uA',amp);
figure('Color','w','Name',name,'NumberTitle','off','Position',position);
ax = axes; hold(ax,'on');

cols.A = [0.20 0.55 0.25];
cols.B = [0.55 0.25 0.70];
cols.AtoB = [0.85 0.25 0.12];
cols.BtoA = [0.95 0.55 0.05];
cols.Combined = [0.75 0.15 0.55];
use_ab_connection = connect_ab && selected_line(lines,'AB');
[electrode_A,electrode_B] = stimulation_names(M);
label_A_to_B = sprintf('%s -> %s',electrode_A,electrode_B);
label_B_to_A = sprintf('%s -> %s',electrode_B,electrode_A);

x_reference = [0 ptds];
if selected_line(lines,'A')
    plot_reference(ax,x_reference,A,cols.A,'--', ...
        sprintf('%s alone',electrode_A),show_sem,lw,ms);
end
if selected_line(lines,'B')
    plot_reference(ax,x_reference,B,cols.B,'--', ...
        sprintf('%s alone',electrode_B),show_sem,lw,ms);
end
if combine_orders
    show_seq = selected_line(lines,'A_to_B') || selected_line(lines,'B_to_A');
    if selected_line(lines,'Linear')
        plot_prediction_population(ax,prediction_ptds,Linear_Combined, ...
            [0.15 0.15 0.15],'Shifted linear prediction',lw);
    end
    if show_seq
        plot_sequence_population(ax,ptds,Combined,AB_Combined,use_ab_connection,cols.Combined, ...
            'Sequential (orders combined)',show_sem,lw,ms);
    end
    if selected_line(lines,'AB')
        plot_ab_population_marker(ax,AB_Combined,cols.Combined, ...
            'A+B simultaneous (combined population)',show_sem,lw,ms);
    end
else
    if selected_line(lines,'Linear') && selected_line(lines,'A_to_B')
        plot_prediction_population(ax,prediction_ptds,Linear_AtoB,cols.AtoB, ...
            sprintf('Shifted linear prediction (%s)',label_A_to_B),lw);
    end
    if selected_line(lines,'Linear') && selected_line(lines,'B_to_A')
        plot_prediction_population(ax,prediction_ptds,Linear_BtoA,cols.BtoA, ...
            sprintf('Shifted linear prediction (%s)',label_B_to_A),lw);
    end
    if selected_line(lines,'A_to_B')
        plot_sequence_population(ax,ptds,AtoB,AB_AtoB,use_ab_connection,cols.AtoB, ...
            label_A_to_B,show_sem,lw,ms);
    end
    if selected_line(lines,'B_to_A')
        plot_sequence_population(ax,ptds,BtoA,AB_BtoA,use_ab_connection,cols.BtoA, ...
            label_B_to_A,show_sem,lw,ms);
    end
    if selected_line(lines,'AB') && selected_line(lines,'A_to_B')
        plot_ab_population_marker(ax,AB_AtoB,cols.AtoB, ...
            sprintf('A+B simultaneous (%s population)',label_A_to_B),show_sem,lw,ms);
    end
    if selected_line(lines,'AB') && selected_line(lines,'B_to_A')
        plot_ab_population_marker(ax,AB_BtoA,cols.BtoA, ...
            sprintf('A+B simultaneous (%s population)',label_B_to_A),show_sem,lw,ms);
    end
end

decorate_axes(ax,amp,settings,M);
end

function plot_individual_page(amp,ptds,prediction_ptds,A,B, ...
        AB_AtoB,AB_BtoA,AB_Combined, ...
        AtoB,BtoA,Combined,Linear_AtoB,Linear_BtoA,Linear_Combined, ...
        page_rows,channel_rows,lines,combine_orders,connect_ab,settings,M, ...
        position,lw,ms,page,nPages)
% Plot channel-level trial means. Trial counts are intentionally omitted.
name = sprintf('Individual channels | %.4g uA | page %d of %d',amp,page,nPages);
figure('Color','w','Name',name,'NumberTitle','off','Position',position);
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
title(tl,sprintf('%s | %.4g uA',stimulation_pair_label(M),amp), ...
    'FontWeight','bold','Interpreter','none');

for local_row = page_rows
    ax = nexttile(tl); hold(ax,'on');
    plot_one_channel(ax,ptds,prediction_ptds,A(local_row),B(local_row), ...
        AB_AtoB(local_row),AB_BtoA(local_row),AB_Combined(local_row), ...
        AtoB(local_row,:),BtoA(local_row,:),Combined(local_row,:), ...
        Linear_AtoB(local_row,:),Linear_BtoA(local_row,:),Linear_Combined(local_row,:), ...
        lines,combine_orders,connect_ab,M,lw,ms);

    table_row = channel_rows(local_row);
    channel_id = M.channels.ChannelIndex(table_row);
    depth_id = M.channels.DepthChannel(table_row);
    title(ax,sprintf('Recording channel %d | depth channel %d', ...
        channel_id,depth_id),'Interpreter','none');
    decorate_axes(ax,amp,settings,M);
end
end

function plot_one_channel(ax,ptds,prediction_ptds,A,B, ...
        AB_AtoB,AB_BtoA,AB_Combined, ...
        AtoB,BtoA,Combined,Linear_AtoB,Linear_BtoA,Linear_Combined, ...
        lines,combine_orders,connect_ab,M,lw,ms)
% Draw one recording channel without trial-count annotations.
cols.A = [0.20 0.55 0.25]; cols.B = [0.55 0.25 0.70];
cols.AtoB = [0.85 0.25 0.12];
cols.BtoA = [0.95 0.55 0.05]; cols.Combined = [0.75 0.15 0.55];
x_reference = [0 ptds];
use_ab_connection = connect_ab && selected_line(lines,'AB');
[electrode_A,electrode_B] = stimulation_names(M);
label_A_to_B = sprintf('%s -> %s',electrode_A,electrode_B);
label_B_to_A = sprintf('%s -> %s',electrode_B,electrode_A);

if selected_line(lines,'A') && isfinite(A)
    plot(ax,x_reference,repmat(A,size(x_reference)),'--', ...
        'Color',cols.A,'LineWidth',lw,'DisplayName',sprintf('%s alone',electrode_A));
end
if selected_line(lines,'B') && isfinite(B)
    plot(ax,x_reference,repmat(B,size(x_reference)),'--', ...
        'Color',cols.B,'LineWidth',lw,'DisplayName',sprintf('%s alone',electrode_B));
end
if combine_orders
    if selected_line(lines,'Linear') && any(isfinite(Linear_Combined))
        valid = isfinite(Linear_Combined);
        plot(ax,prediction_ptds(valid),Linear_Combined(valid),'--', ...
            'Color',[0.15 0.15 0.15],'LineWidth',lw, ...
            'DisplayName','Shifted linear prediction');
    end
    if selected_line(lines,'A_to_B') || selected_line(lines,'B_to_A')
        plot_sequence_values(ax,ptds,Combined,AB_Combined,use_ab_connection,cols.Combined, ...
            'Sequential combined',lw,ms);
    end
    if selected_line(lines,'AB') && isfinite(AB_Combined)
        plot(ax,0,AB_Combined,'o','Color',cols.Combined, ...
            'MarkerFaceColor',cols.Combined,'MarkerSize',ms,'LineWidth',lw, ...
            'DisplayName','A+B simultaneous');
    end
else
    if selected_line(lines,'Linear') && selected_line(lines,'A_to_B') && ...
            any(isfinite(Linear_AtoB))
        valid = isfinite(Linear_AtoB);
        plot(ax,prediction_ptds(valid),Linear_AtoB(valid),'--', ...
            'Color',cols.AtoB,'LineWidth',lw, ...
            'DisplayName',sprintf('Shifted linear (%s)',label_A_to_B));
    end
    if selected_line(lines,'Linear') && selected_line(lines,'B_to_A') && ...
            any(isfinite(Linear_BtoA))
        valid = isfinite(Linear_BtoA);
        plot(ax,prediction_ptds(valid),Linear_BtoA(valid),'--', ...
            'Color',cols.BtoA,'LineWidth',lw, ...
            'DisplayName',sprintf('Shifted linear (%s)',label_B_to_A));
    end
    if selected_line(lines,'A_to_B')
        plot_sequence_values(ax,ptds,AtoB,AB_AtoB,use_ab_connection,cols.AtoB,label_A_to_B,lw,ms);
    end
    if selected_line(lines,'B_to_A')
        plot_sequence_values(ax,ptds,BtoA,AB_BtoA,use_ab_connection,cols.BtoA,label_B_to_A,lw,ms);
    end
    if selected_line(lines,'AB') && selected_line(lines,'A_to_B') && isfinite(AB_AtoB)
        plot(ax,0,AB_AtoB,'o','Color',cols.AtoB,'MarkerFaceColor',cols.AtoB, ...
            'MarkerSize',ms,'LineWidth',lw, ...
            'DisplayName',sprintf('A+B (%s population)',label_A_to_B));
    end
    if selected_line(lines,'AB') && selected_line(lines,'B_to_A') && isfinite(AB_BtoA)
        plot(ax,0,AB_BtoA,'o','Color',cols.BtoA,'MarkerFaceColor',cols.BtoA, ...
            'MarkerSize',ms,'LineWidth',lw, ...
            'DisplayName',sprintf('A+B (%s population)',label_B_to_A));
    end
end
end

function plot_ab_population_marker(ax,values,col,label,show_sem,lw,ms)
[mu,se] = population_stats(values);
if ~isfinite(mu), return; end
if show_sem && isfinite(se)
    errorbar(ax,0,mu,se,'o','Color',col,'MarkerFaceColor',col, ...
        'MarkerSize',ms,'LineWidth',lw,'CapSize',8,'DisplayName',label);
else
    plot(ax,0,mu,'o','Color',col,'MarkerFaceColor',col, ...
        'MarkerSize',ms,'LineWidth',lw,'DisplayName',label);
end
end

function plot_sequence_population(ax, ptds, data, AB, connect_ab, ...
    col, label, show_sem, lw, ms)
% Plot the observed sequential response averaged across recording channels.
%
% INPUTS
%   data       - channels × PTDs response matrix
%   AB         - simultaneous A+B response for each channel
%   connect_ab - if true, add simultaneous A+B at PTD = 0
%   show_sem   - if true, show SEM error bars

% Calculate population mean and SEM at every PTD
mu = nan(1, numel(ptds));
se = nan(1, numel(ptds));

for k = 1:numel(ptds)
    [mu(k), se(k)] = population_stats(data(:,k));
end

% Calculate simultaneous A+B population mean and SEM
[ab_mu, ab_se] = population_stats(AB);

% Optionally connect simultaneous A+B at PTD = 0 to the sequence
if connect_ab && isfinite(ab_mu)
    x   = [0, ptds];
    y   = [ab_mu, mu];
    err = [ab_se, se];
else
    x   = ptds;
    y   = mu;
    err = se;
end

% Remove unavailable conditions
valid = isfinite(x) & isfinite(y);

if ~any(valid)
    return;
end

% Plot observed results
if show_sem
    errorbar(ax, x(valid), y(valid), err(valid), '-o', ...
        'Color', col, ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', ms, ...
        'LineWidth', lw, ...
        'CapSize', 7, ...
        'DisplayName', label);
else
    plot(ax, x(valid), y(valid), '-o', ...
        'Color', col, ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', ms, ...
        'LineWidth', lw, ...
        'DisplayName', label);
end
end
function plot_sequence_values(ax,ptds,values,AB,connect_ab,col,label,lw,ms)
% Plot one channel's sequential condition means.
if connect_ab && isfinite(AB)
    x = [0 ptds]; y = [AB values];
else
    x = ptds; y = values;
end

valid = isfinite(y);
if any(valid)
    plot(ax,x(valid),y(valid),'-o','Color',col,'MarkerFaceColor','w', ...
        'MarkerSize',ms,'LineWidth',lw,'DisplayName',label);
end
end

function plot_prediction_population(ax, ptds, data, col, label, lw)
% Plot the mean PTD-dependent linear prediction.
% Error bars and markers are intentionally omitted.

mu = nan(1, numel(ptds));

for k = 1:numel(ptds)
    mu(k) = population_stats(data(:,k));
end

valid = isfinite(ptds) & isfinite(mu);

if ~any(valid)
    return;
end

plot(ax, ptds(valid), mu(valid), '--', ...
    'Color', col, ...
    'LineWidth', lw, ...
    'DisplayName', label);
end

function plot_reference(ax,x,values,col,style,label,show_sem,lw,ms)
% Plot a PTD-independent population reference and its SEM band.
[mu,se] = population_stats(values);
if ~isfinite(mu), return; end
if show_sem && isfinite(se)
    fill(ax,[x fliplr(x)],[repmat(mu+se,size(x)) repmat(mu-se,size(x))], ...
        col,'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
end
plot(ax,x,repmat(mu,size(x)),style,'Color',col,'LineWidth',lw, ...
    'MarkerSize',ms,'DisplayName',label);
end

function [mu,se] = population_stats(values)
% Mean and SEM across finite recording-channel means.
values = double(values(:)); values = values(isfinite(values));
if isempty(values), mu=NaN; se=NaN; return; end
mu = mean(values);
if numel(values)>1, se=std(values,0)/sqrt(numel(values));
else, se=NaN; end
end

function tf = selected_line(lines,name)
tf = any(strcmp(lines,name));
end

function [channel_A,channel_B] = stimulation_depth_channels(M,electrode_type)
% Reproduce the relevant ProbeMAP lookup locally, so the plotter remains
% independent of Depth_s, ProbeMAP, and the original experiment folder.
if ~isfield(M,'stimulation') || ~isfield(M.stimulation,'electrode_A') || ...
        ~isfield(M.stimulation,'electrode_B')
    error('MultiISIPlot:StimulationNames', ...
        'ModelData stimulation electrode names are missing.');
end
channel_A = native_name_to_depth_channel( ...
    char(string(M.stimulation.electrode_A)),electrode_type);
channel_B = native_name_to_depth_channel( ...
    char(string(M.stimulation.electrode_B)),electrode_type);
end

function depth_channel = native_name_to_depth_channel(name,electrode_type)
token = regexp(name,'^([A-D])-([0-9]+)$','tokens','once');
if isempty(token)
    error('MultiISIPlot:ElectrodeName','Cannot interpret electrode name %s.',name);
end
bank = token{1}; native_number = str2double(token{2});

map_rigid = [22 27 18 28 13 4 9 29 19 3 14 2 12 6 17 30 ...
    20 25 15 31 11 1 16 0 21 5 10 26 23 7 8 24];
map_flexible = [8 7 9 6 10 5 12 3 13 2 14 1 23 24 22 25 ...
    21 26 19 28 18 29 17 30 16 31 20 27 15 0 11 4];
map_four_A = [31 7 0 24 30 6 1 25 29 5 2 26 28 4 3 27 ...
    8 16 23 15 9 17 22 14 10 18 21 13 11 19 20 12];
map_four_B = [8 11 9 15 10 19 12 23 13 22 14 21 16 20 17 18 ...
    4 7 0 6 28 5 24 3 25 2 26 1 27 31 29 30];

switch electrode_type
    case 0
        if ~strcmp(bank,'A'), error('MultiISIPlot:ElectrodeBank', ...
                'Rigid single-shank mapping expects bank A.'); end
        position = find(map_rigid==native_number,1); offset = 0;
    case 1
        if ~strcmp(bank,'A'), error('MultiISIPlot:ElectrodeBank', ...
                'Flexible single-shank mapping expects bank A.'); end
        position = find(map_flexible==native_number,1); offset = 0;
    case 2
        switch bank
            case 'A', map=map_four_A; offset=0;
            case 'B', map=map_four_B; offset=32;
            case 'C', map=map_four_A; offset=64;
            case 'D', map=map_four_B; offset=96;
        end
        position = find(map==native_number,1);
end
if isempty(position)
    error('MultiISIPlot:ElectrodeMapping', ...
        'Electrode %s was not found in the selected probe map.',name);
end
depth_channel = offset+position;
end

function [A,B] = stimulation_names(M)
% Prefer anatomical channel numbers derived in the main script.
if isfield(M,'display_stimulation_A') && isfield(M,'display_stimulation_B')
    A = M.display_stimulation_A;
    B = M.display_stimulation_B;
elseif isfield(M,'stimulation') && isfield(M.stimulation,'electrode_A') && ...
        isfield(M.stimulation,'electrode_B')
    A = char(string(M.stimulation.electrode_A));
    B = char(string(M.stimulation.electrode_B));
else
    A = 'A'; B = 'B';
end
end

function label = stimulation_pair_label(M)
[A,B] = stimulation_names(M);
label = sprintf('%s + %s',A,B);
end

function decorate_axes(ax,amp,settings,M)
% Apply consistent labels without displaying trial counts.
xline(ax,0,':','Color',[0.35 0.35 0.35],'HandleVisibility','off');
xlabel(ax,'ISI / PTD (ms)');
ylabel(ax,count_axis_label(settings));
title_text = sprintf('%s | %.4g uA',stimulation_pair_label(M),amp);
if isempty(ax.Title.String)
    title(ax,title_text,'Interpreter','none');
end
grid(ax,'on'); box(ax,'off');
legend(ax,'Location','best','Box','off');
end

function label = count_axis_label(settings)
w = settings.response_window_ms;
if ~settings.subtract_baseline
    label = sprintf('Raw spike count / trial [%g,%g) ms',w(1),w(2));
elseif settings.force_nonnegative
    label = sprintf('Nonnegative baseline-corrected spikes / trial [%g,%g) ms',w(1),w(2));
else
    label = sprintf('Baseline-corrected spike count / trial [%g,%g) ms',w(1),w(2));
end
end
