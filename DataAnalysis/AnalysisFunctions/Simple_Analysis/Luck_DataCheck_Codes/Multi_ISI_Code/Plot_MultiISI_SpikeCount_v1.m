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
%   - The linear prediction is response(A) + response(B).
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
model_data_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/Multi_ISI_Output/DX020_XiaISISimSeq1_Pair_A017_A018_MultiISI_ModelData.mat'];

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
%   'Linear'  : single-A response + single-B response
%
% A and B are PTD-independent and therefore appear as horizontal lines.
Lines_To_Plot = {'AB','A_to_B','B_to_A','Linear'};

% false: retain A->B and B->A as separate curves.
% true:  average the available order means within each channel. If only one
%        order exists at a PTD, that available order is retained.
Combine_Sequential_Orders = false;

% When true, prepend the simultaneous A+B point at 0 ms to every displayed
% sequential curve. Change to false later if a separate point is preferred.
Connect_Simultaneous_To_Sequential = true;

% Metric stored for every clean trial in the exported ModelData file.
% Recommended default:
%   'baseline_corrected_whole_count' = count in [2,45) ms minus the
%                                      expected baseline count
% Other useful choices:
%   'whole_response_count'           = raw count in [2,45) ms
%   'early_window_count'             = raw count in [2,20) ms
%   'second_aligned_count'           = raw count in [PTD+2,PTD+20) ms
Metric = 'baseline_corrected_whole_count';

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

valid_metrics = {'baseline_corrected_whole_count', ...
    'whole_response_count','early_window_count','second_aligned_count'};
if ~ismember(Metric,valid_metrics)
    error('MultiISIPlot:Metric', ...
        'Metric must be one of: %s',strjoin(valid_metrics,', '));
end

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
        selected_channel_rows,Metric);
    B_values  = condition_channel_means(ModelData.conditions,'B',0,amp, ...
        selected_channel_rows,Metric);
    AB_values = condition_channel_means(ModelData.conditions,'AB',0,amp, ...
        selected_channel_rows,Metric);

    nCh = numel(selected_channel_rows);
    AtoB = nan(nCh,numel(selected_ptds));
    BtoA = nan(nCh,numel(selected_ptds));
    for ip = 1:numel(selected_ptds)
        ptd = selected_ptds(ip);
        AtoB(:,ip) = condition_channel_means( ...
            ModelData.conditions,'A_to_B',ptd,amp,selected_channel_rows,Metric);
        BtoA(:,ip) = condition_channel_means( ...
            ModelData.conditions,'B_to_A',ptd,amp,selected_channel_rows,Metric);
    end

    % The linear reference is calculated separately for every recording
    % channel before any population averaging is performed.
    Linear_values = A_values + B_values;

    % Combining orders gives equal weight to the available order means.
    % It does not pool all trials and therefore does not require balancing.
    Combined = mean_available_orders(AtoB,BtoA);

    if plot_average
        plot_population_figure(amp,selected_ptds,A_values,B_values,AB_values, ...
            AtoB,BtoA,Combined,Linear_values,Lines_To_Plot, ...
            Combine_Sequential_Orders,Connect_Simultaneous_To_Sequential, ...
            Show_Population_SEM,Metric,ModelData,figure_position_average, ...
            line_width,marker_size);
    end

    if plot_individual
        nPages = ceil(nCh/Channels_Per_Figure);
        for page = 1:nPages
            first_ch = (page-1)*Channels_Per_Figure+1;
            last_ch = min(page*Channels_Per_Figure,nCh);
            page_rows = first_ch:last_ch;
            plot_individual_page(amp,selected_ptds,A_values,B_values,AB_values, ...
                AtoB,BtoA,Combined,Linear_values,page_rows,selected_channel_rows, ...
                Lines_To_Plot,Combine_Sequential_Orders, ...
                Connect_Simultaneous_To_Sequential,Metric,ModelData, ...
                figure_position_individual,line_width,marker_size,page,nPages);
        end
    end
end

%% ========================= LOCAL FUNCTIONS ============================

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

function means = condition_channel_means(Conditions,code,ptd,amp,channel_rows,metric)
% Find one condition/amplitude and calculate its per-channel trial means.
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
    if row>numel(channels) || ~isfield(channels(row).trial_metrics,metric)
        continue;
    end
    trial_values = double(channels(row).trial_metrics.(metric)(:));
    trial_values = trial_values(isfinite(trial_values));
    if ~isempty(trial_values)
        means(k) = mean(trial_values);
    end
end
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

function plot_population_figure(amp,ptds,A,B,AB,AtoB,BtoA,Combined,Linear, ...
        lines,combine_orders,connect_ab,show_sem,metric,M,position,lw,ms)
% Plot the population mean. SEM is calculated across channel means.
name = sprintf('Average spike count vs ISI | %.4g uA',amp);
figure('Color','w','Name',name,'NumberTitle','off','Position',position);
ax = axes; hold(ax,'on');

cols.A = [0.20 0.55 0.25];
cols.B = [0.55 0.25 0.70];
cols.AB = [0.10 0.25 0.75];
cols.AtoB = [0.85 0.25 0.12];
cols.BtoA = [0.95 0.55 0.05];
cols.Combined = [0.75 0.15 0.55];
cols.Linear = [0.15 0.15 0.15];
use_ab_connection = connect_ab && selected_line(lines,'AB');

x_reference = [0 ptds];
if selected_line(lines,'A')
    plot_reference(ax,x_reference,A,cols.A,'--','A alone',show_sem,lw,ms);
end
if selected_line(lines,'B')
    plot_reference(ax,x_reference,B,cols.B,'--','B alone',show_sem,lw,ms);
end
if selected_line(lines,'Linear')
    plot_reference(ax,x_reference,Linear,cols.Linear,'--', ...
        'Linear prediction: A + B',show_sem,lw,ms);
end

if combine_orders
    show_seq = selected_line(lines,'A_to_B') || selected_line(lines,'B_to_A');
    if show_seq
        plot_sequence_population(ax,ptds,Combined,AB,use_ab_connection,cols.Combined, ...
            'Sequential (orders combined)',show_sem,lw,ms);
    end
else
    if selected_line(lines,'A_to_B')
        plot_sequence_population(ax,ptds,AtoB,AB,use_ab_connection,cols.AtoB, ...
            'A -> B',show_sem,lw,ms);
    end
    if selected_line(lines,'B_to_A')
        plot_sequence_population(ax,ptds,BtoA,AB,use_ab_connection,cols.BtoA, ...
            'B -> A',show_sem,lw,ms);
    end
end

if selected_line(lines,'AB')
    [mu,se] = population_stats(AB);
    if isfinite(mu)
        if show_sem && isfinite(se)
            errorbar(ax,0,mu,se,'o','Color',cols.AB,'MarkerFaceColor',cols.AB, ...
                'MarkerSize',ms,'LineWidth',lw,'CapSize',8, ...
                'DisplayName','A+B simultaneous');
        else
            plot(ax,0,mu,'o','Color',cols.AB,'MarkerFaceColor',cols.AB, ...
                'MarkerSize',ms,'LineWidth',lw,'DisplayName','A+B simultaneous');
        end
    end
end

decorate_axes(ax,amp,metric,M);
end

function plot_individual_page(amp,ptds,A,B,AB,AtoB,BtoA,Combined,Linear, ...
        page_rows,channel_rows,lines,combine_orders,connect_ab,metric,M, ...
        position,lw,ms,page,nPages)
% Plot channel-level trial means. Trial counts are intentionally omitted.
name = sprintf('Individual channels | %.4g uA | page %d of %d',amp,page,nPages);
figure('Color','w','Name',name,'NumberTitle','off','Position',position);
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
title(tl,sprintf('%s | %.4g uA',char(M.pair_key),amp), ...
    'FontWeight','bold','Interpreter','none');

for local_row = page_rows
    ax = nexttile(tl); hold(ax,'on');
    plot_one_channel(ax,ptds,A(local_row),B(local_row),AB(local_row), ...
        AtoB(local_row,:),BtoA(local_row,:),Combined(local_row,:), ...
        Linear(local_row),lines,combine_orders,connect_ab,lw,ms);

    table_row = channel_rows(local_row);
    channel_id = M.channels.ChannelIndex(table_row);
    depth_id = M.channels.DepthChannel(table_row);
    title(ax,sprintf('Recording channel %d | depth channel %d', ...
        channel_id,depth_id),'Interpreter','none');
    decorate_axes(ax,amp,metric,M);
end
end

function plot_one_channel(ax,ptds,A,B,AB,AtoB,BtoA,Combined,Linear, ...
        lines,combine_orders,connect_ab,lw,ms)
% Draw one recording channel without trial-count annotations.
cols.A = [0.20 0.55 0.25]; cols.B = [0.55 0.25 0.70];
cols.AB = [0.10 0.25 0.75]; cols.AtoB = [0.85 0.25 0.12];
cols.BtoA = [0.95 0.55 0.05]; cols.Combined = [0.75 0.15 0.55];
cols.Linear = [0.15 0.15 0.15];
x_reference = [0 ptds];
use_ab_connection = connect_ab && selected_line(lines,'AB');

if selected_line(lines,'A') && isfinite(A)
    plot(ax,x_reference,repmat(A,size(x_reference)),'--', ...
        'Color',cols.A,'LineWidth',lw,'DisplayName','A alone');
end
if selected_line(lines,'B') && isfinite(B)
    plot(ax,x_reference,repmat(B,size(x_reference)),'--', ...
        'Color',cols.B,'LineWidth',lw,'DisplayName','B alone');
end
if selected_line(lines,'Linear') && isfinite(Linear)
    plot(ax,x_reference,repmat(Linear,size(x_reference)),'--', ...
        'Color',cols.Linear,'LineWidth',lw,'DisplayName','Linear: A+B');
end

if combine_orders
    if selected_line(lines,'A_to_B') || selected_line(lines,'B_to_A')
        plot_sequence_values(ax,ptds,Combined,AB,use_ab_connection,cols.Combined, ...
            'Sequential combined',lw,ms);
    end
else
    if selected_line(lines,'A_to_B')
        plot_sequence_values(ax,ptds,AtoB,AB,use_ab_connection,cols.AtoB,'A -> B',lw,ms);
    end
    if selected_line(lines,'B_to_A')
        plot_sequence_values(ax,ptds,BtoA,AB,use_ab_connection,cols.BtoA,'B -> A',lw,ms);
    end
end
if selected_line(lines,'AB') && isfinite(AB)
    plot(ax,0,AB,'o','Color',cols.AB,'MarkerFaceColor',cols.AB, ...
        'MarkerSize',ms,'LineWidth',lw,'DisplayName','A+B simultaneous');
end
end

function plot_sequence_population(ax,ptds,data,AB,connect_ab,col,label,show_sem,lw,ms)
% Plot population sequential values, optionally beginning at A+B/ISI 0.
mu = nan(1,numel(ptds)); se = nan(size(mu));
for k = 1:numel(ptds), [mu(k),se(k)] = population_stats(data(:,k)); end
[ab_mu,ab_se] = population_stats(AB);
if connect_ab && isfinite(ab_mu)
    x = [0 ptds]; y = [ab_mu mu]; err = [ab_se se];
else
    x = ptds; y = mu; err = se;
end
valid = isfinite(y);
if ~any(valid), return; end
if show_sem
    errorbar(ax,x(valid),y(valid),err(valid),'-o','Color',col, ...
        'MarkerFaceColor','w','MarkerSize',ms,'LineWidth',lw,'CapSize',7, ...
        'DisplayName',label);
else
    plot(ax,x(valid),y(valid),'-o','Color',col,'MarkerFaceColor','w', ...
        'MarkerSize',ms,'LineWidth',lw,'DisplayName',label);
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

function decorate_axes(ax,amp,metric,M)
% Apply consistent labels without displaying trial counts.
xline(ax,0,':','Color',[0.35 0.35 0.35],'HandleVisibility','off');
xlabel(ax,'ISI / PTD (ms)');
ylabel(ax,metric_axis_label(metric,M));
title_text = sprintf('%s | %.4g uA',char(M.pair_key),amp);
if isempty(ax.Title.String)
    title(ax,title_text,'Interpreter','none');
end
grid(ax,'on'); box(ax,'off');
legend(ax,'Location','best','Box','off');
end

function label = metric_axis_label(metric,M)
switch metric
    case 'baseline_corrected_whole_count'
        if isfield(M,'metadata') && isfield(M.metadata,'whole_response_window_ms')
            w = M.metadata.whole_response_window_ms;
            label = sprintf('Baseline-corrected spike count / trial [%g,%g) ms',w(1),w(2));
        else
            label = 'Baseline-corrected spike count / trial';
        end
    case 'whole_response_count'
        label = 'Spike count / trial in whole-response window';
    case 'early_window_count'
        label = 'Spike count / trial in early-response window';
    case 'second_aligned_count'
        label = 'Spike count / trial relative to second pulse';
end
end
