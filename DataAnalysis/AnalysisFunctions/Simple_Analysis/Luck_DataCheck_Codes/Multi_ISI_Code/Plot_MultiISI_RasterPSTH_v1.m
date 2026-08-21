%% MULTI-ISI RASTER + PSTH PLOTTER
%
% PURPOSE
%   Loads only one exported MultiISI ModelData MAT file and plots clean,
%   trigger-aligned spikes. Original experiment files, trigger files,
%   sp_corr files, and external analysis functions are not required.
%
% FIGURE ORGANIZATION
%   - One figure per available condition x amplitude.
%   - Every selected responding recording channel occupies one tile.
%   - Left y-axis: PSTH firing rate (spikes/second).
%   - Right y-axis: clean trial-by-trial raster.
%
% IMPORTANT
%   Different recording channels may contain different clean trial counts
%   because channel-specific bad trials were removed before export. Each
%   PSTH is normalized by that channel's own number of clean trials.
%
% OUTPUT
%   Figures are displayed only. No files are saved or modified.
%   Routine messages are not printed in the Command Window.

clear;
clc;

%% ========================= USER SETTINGS ==============================

model_data_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/Multi_ISI_Output/DX020_XiaISISimSeq1_Pair_A017_A018_MultiISI_ModelData.mat'];

% Select conditions using their general codes:
%   A, B, AB, A_to_B, B_to_A
%
% Examples:
%   Condition_Codes_To_Plot = 'all';
%   Condition_Codes_To_Plot = {'AB','A_to_B','B_to_A'};
%   Condition_Codes_To_Plot = {'A_to_B'};
Condition_Codes_To_Plot = 'all';

% PTD selection applies only to sequential conditions. Empty means every
% positive PTD stored in ModelData. A small default avoids opening dozens
% of figures during the first check.
PTDs_To_Plot = [];

% When true, requested A, B, and AB reference conditions remain visible
% even when PTDs_To_Plot contains only positive sequential PTDs.
Include_Reference_Conditions = true;

% Empty means every amplitude stored under each selected condition.
% Conditions do not need to contain identical amplitude sets.
Amplitudes_To_Plot = [10];

% Use 'all' for all exported responding channels, or enter anatomical/depth
% channel numbers such as [35 36 40 44].
Channels_To_Plot = 'all';

% true: A->B and B->A figures use their own manually corrected response
% populations saved by the revised MultiISI exporter (ModelData v1.1).
% false: every condition uses the complete exported channel union.
Use_Order_Specific_Response_Channels = true;

% A, B, and simultaneous AB are reference conditions without one unique
% sequential direction. Choose which response population they display:
%   'union'  = all channels responsive to either order (recommended)
%   'A_to_B' = only the A->B response population
%   'B_to_A' = only the B->A response population
Reference_Response_Population = 'union';

% Display window relative to the first stimulation pulse. It must remain
% inside ModelData.metadata.stored_spike_window_ms.
raster_win_ms = [-50 80];

% PSTH settings. Smoothing uses a centred moving mean and no toolbox.
psth_bin_ms = 1;
psth_smooth_ms = 5;

% When true, all channel tiles within one figure use the same PSTH y-limit.
Use_Shared_PSTH_YAxis = false;

% Safety limit. Increase this deliberately, or use Inf, when plotting many
% PTDs and amplitudes in one run.
Maximum_Figures = 35;

% Figure presentation only.
Close_Existing_Figures = true;
figure_position = [40 40 1600 900];

%% ========================== LOAD DATA =================================
if Close_Existing_Figures
    close all;
end
if ~isfile(model_data_file)
    error('MultiISIRaster:MissingFile', ...
        'ModelData file not found:\n%s',model_data_file);
end
S = load(model_data_file,'ModelData');
if ~isfield(S,'ModelData')
    error('MultiISIRaster:InvalidFile', ...
        'The selected MAT file does not contain ModelData.');
end
D = S.ModelData;

required_top = {'format_name','metadata','pair_key','channels','conditions'};
require_fields(D,required_top,'ModelData');
if ~strcmp(D.format_name,'MultiISIStimulationModelData')
    error('MultiISIRaster:WrongFormat', ...
        'The selected file is not a multi-ISI ModelData package.');
end
if isempty(D.conditions) || height(D.channels) < 1
    error('MultiISIRaster:EmptyPackage','ModelData contains no plottable data.');
end

required_channel_variables = {'ChannelIndex','DepthChannel'};
for k = 1:numel(required_channel_variables)
    if ~ismember(required_channel_variables{k},D.channels.Properties.VariableNames)
        error('MultiISIRaster:MissingChannelField', ...
            'ModelData.channels is missing %s.',required_channel_variables{k});
    end
end

if Use_Order_Specific_Response_Channels
    response_variables = {'Responsive_A_to_B','Responsive_B_to_A'};
    for k = 1:numel(response_variables)
        if ~ismember(response_variables{k},D.channels.Properties.VariableNames)
            error('MultiISIRaster:ResponseMasksMissing', ...
                ['Order-specific response masks are missing. Rerun the revised ' ...
                 'Run_MultiISI_Export.m and use the new ModelData file.']);
        end
    end
end
valid_reference_populations = {'union','A_to_B','B_to_A'};
if ~any(strcmpi(Reference_Response_Population,valid_reference_populations))
    error('MultiISIRaster:ReferencePopulation', ...
        'Reference_Response_Population must be union, A_to_B, or B_to_A.');
end

%% ====================== RESOLVE USER CHOICES ==========================
condition_indices = resolve_conditions(Condition_Codes_To_Plot, ...
    PTDs_To_Plot,Include_Reference_Conditions,D.conditions);
[channel_columns,~] = ...
    resolve_channels(Channels_To_Plot,D.channels);

validate_window(raster_win_ms,'raster_win_ms');
stored_window = double(D.metadata.stored_spike_window_ms(:).');
if raster_win_ms(1) < stored_window(1) || raster_win_ms(2) > stored_window(2)
    error('MultiISIRaster:WindowOutsideData', ...
        ['raster_win_ms [%g %g] lies outside the stored range ' ...
         '[%g %g] ms.'],raster_win_ms,stored_window);
end
if ~isscalar(psth_bin_ms) || ~isfinite(psth_bin_ms) || psth_bin_ms <= 0
    error('MultiISIRaster:Bin','psth_bin_ms must be one positive number.');
end
if ~isscalar(psth_smooth_ms) || ~isfinite(psth_smooth_ms) || ...
        psth_smooth_ms < 0
    error('MultiISIRaster:Smoothing', ...
        'psth_smooth_ms must be one nonnegative number.');
end
if ~isscalar(Maximum_Figures) || isnan(Maximum_Figures) || Maximum_Figures < 1
    error('MultiISIRaster:FigureLimit', ...
        'Maximum_Figures must be a positive number or Inf.');
end

plot_plan = build_plot_plan(condition_indices,Amplitudes_To_Plot,D.conditions);
if isempty(plot_plan)
    error('MultiISIRaster:NoConditions', ...
        'No stored condition/amplitude combination matches the selections.');
end
if numel(plot_plan) > Maximum_Figures
    error('MultiISIRaster:TooManyFigures', ...
        ['The current settings would create %d figures, exceeding ' ...
         'Maximum_Figures=%g. Select fewer PTDs/amplitudes or increase ' ...
         'Maximum_Figures.'],numel(plot_plan),Maximum_Figures);
end

%% ====================== PSTH BIN DEFINITIONS ==========================
edges = raster_win_ms(1):psth_bin_ms:raster_win_ms(2);
if isempty(edges) || edges(1) ~= raster_win_ms(1)
    edges = raster_win_ms(1);
end
if edges(end) < raster_win_ms(2)
    edges(end+1) = raster_win_ms(2);
end
centres = edges(1:end-1)+diff(edges)/2;
smooth_bins = max(1,round(psth_smooth_ms/psth_bin_ms));

%% ========================== MAIN PLOTTING =============================
for iPlan = 1:numel(plot_plan)
    ic = plot_plan(iPlan).condition_index;
    ia = plot_plan(iPlan).amplitude_index;
    condition = D.conditions(ic);
    block = condition.amplitude(ia);
    validate_amplitude_block(block,height(D.channels));

    [condition_channel_columns,population_label] = ...
        channels_for_condition(channel_columns,condition.code,D.channels, ...
        Use_Order_Specific_Response_Channels,Reference_Response_Population);
    if isempty(condition_channel_columns)
        continue;
    end
    condition_depth_channels = double( ...
        D.channels.DepthChannel(condition_channel_columns).');
    nSelectedChannels = numel(condition_channel_columns);
    nTileColumns = ceil(sqrt(nSelectedChannels));
    nTileRows = ceil(nSelectedChannels/nTileColumns);

    condition_color = color_for_condition(condition.code);
    channel_spikes = cell(1,nSelectedChannels);
    channel_rates = cell(1,nSelectedChannels);
    individual_ymax = zeros(1,nSelectedChannels);
    shared_rate_max = 0;

    % Precalculate all selected channels so a common y-limit can be used.
    for jcLocal = 1:nSelectedChannels
        jc = condition_channel_columns(jcLocal);
        channel_data = block.channel(jc);
        validate_channel_data(channel_data);
        spikes = restrict_spikes(channel_data.spike_times_ms,raster_win_ms);
        rate = calculate_psth(spikes,edges,smooth_bins);
        channel_spikes{jcLocal} = spikes;
        channel_rates{jcLocal} = rate;
        if isempty(rate), this_max = 0; else, this_max = max(rate); end
        individual_ymax(jcLocal) = nice_rate_limit(this_max);
        shared_rate_max = max(shared_rate_max,this_max);
    end
    shared_ymax = nice_rate_limit(shared_rate_max);

    amp = block.amplitude_uA;
    figure_name = sprintf('%s | %s | %g uA', ...
        D.pair_key,condition.label,amp);
    figure('Color','w','Name',figure_name,'Position',figure_position);
    tl = tiledlayout(nTileRows,nTileColumns, ...
        'TileSpacing','compact','Padding','compact');

    title(tl,sprintf('%s | %s | %g uA | %s | %s', ...
        D.pair_key,condition.label,amp,pulse_description(condition), ...
        population_label), ...
        'Interpreter','none','FontWeight','bold');

    for jcLocal = 1:nSelectedChannels
        jc = condition_channel_columns(jcLocal);
        depth_ch = condition_depth_channels(jcLocal);
        spikes = channel_spikes{jcLocal};
        rate = channel_rates{jcLocal};
        nTrials = numel(spikes);

        ax = nexttile(tl); hold(ax,'on');

        % Left axis: PSTH rate normalized by this channel's clean trials.
        yyaxis(ax,'left');
        plot(ax,centres,rate,'Color',condition_color,'LineWidth',1.4);
        if Use_Shared_PSTH_YAxis
            ylim(ax,[0 shared_ymax]);
        else
            ylim(ax,[0 individual_ymax(jcLocal)]);
        end
        ax.YAxis(1).Color = condition_color;
        if mod(jcLocal-1,nTileColumns) == 0
            ylabel(ax,'Rate (sp/s)');
        end

        % Right axis: one row for every clean trial retained for this channel.
        yyaxis(ax,'right');
        for it = 1:nTrials
            tt = spikes{it};
            if ~isempty(tt)
                plot(ax,tt,it*ones(size(tt)),'.k','MarkerSize',4);
            end
        end
        ylim(ax,[0 max(1,nTrials+1)]);
        ax.YAxis(2).Color = [0 0 0];
        if mod(jcLocal,nTileColumns) == 0 || jcLocal == nSelectedChannels
            ylabel(ax,'Trial');
            if nTrials == 1
                yticks(ax,1);
            elseif nTrials > 1
                yticks(ax,[1 nTrials]);
            else
                set(ax,'YTick',[]);
            end
        else
            set(ax,'YTick',[]);
        end

        % The first pulse is at 0 ms. Sequential stimulation also receives
        % a marker at its positive PTD.
        xline(ax,0,'r--','LineWidth',1);
        later_pulses = condition.pulse_times_ms( ...
            condition.pulse_times_ms > 0);
        for ip = 1:numel(later_pulses)
            xline(ax,later_pulses(ip),'k:','LineWidth',1.1);
        end
        xlim(ax,raster_win_ms);

        title(ax,sprintf('Rec Ch %d%s',depth_ch, ...
            stimulation_channel_tag(depth_ch,D)), ...
            'Interpreter','none','FontSize',10,'FontWeight','bold');
        if ceil(jcLocal/nTileColumns) == nTileRows
            xlabel(ax,'Time (ms)');
        end
        box(ax,'off');
    end
    drawnow;
end

%% =========================== FUNCTIONS ================================

function indices = resolve_conditions(request,ptd_request,include_refs,conditions)
codes = {conditions.code};
if ischar(request) || (isstring(request) && isscalar(request))
    if strcmpi(strtrim(char(string(request))),'all')
        requested_codes = unique(codes,'stable');
    else
        requested_codes = {strtrim(char(string(request)))};
    end
elseif isstring(request) || iscell(request)
    requested_codes = cellstr(string(request));
else
    error('MultiISIRaster:ConditionSelection', ...
        'Condition_Codes_To_Plot must be all, text, or a cell array of codes.');
end
valid_codes = {'A','B','AB','A_to_B','B_to_A'};
for k = 1:numel(requested_codes)
    if ~any(strcmpi(requested_codes{k},valid_codes))
        error('MultiISIRaster:UnknownCode', ...
            'Unknown condition code: %s',requested_codes{k});
    end
end

ptd_request = unique(double(ptd_request(:).'));
if any(ptd_request <= 0)
    error('MultiISIRaster:PTDSelection', ...
        'PTDs_To_Plot must contain positive sequential PTDs only.');
end
keep = false(1,numel(conditions));
for ic = 1:numel(conditions)
    code_requested = any(strcmpi(codes{ic},requested_codes));
    if ~code_requested, continue; end
    is_sequential = strcmp(conditions(ic).stimulation_type,'sequential');
    if is_sequential
        keep(ic) = isempty(ptd_request) || ...
            any(abs(conditions(ic).PTD_ms-ptd_request) <= 1e-6);
    else
        keep(ic) = include_refs;
    end
end
indices = find(keep);
if isempty(indices)
    error('MultiISIRaster:NoSelectedConditions', ...
        'No stored conditions match the requested codes and PTDs.');
end
end

function plan = build_plot_plan(condition_indices,amp_request,conditions)
amp_request = unique(double(amp_request(:).'));
plan = repmat(struct('condition_index',[],'amplitude_index',[]),1,0);
for ic = condition_indices(:).'
    if ~isfield(conditions(ic),'amplitude') || isempty(conditions(ic).amplitude)
        continue;
    end
    stored_amps = [conditions(ic).amplitude.amplitude_uA];
    if isempty(amp_request)
        amp_indices = 1:numel(stored_amps);
    else
        amp_indices = [];
        for requested = amp_request
            match = find(abs(stored_amps-requested) <= 1e-6,1);
            if ~isempty(match), amp_indices(end+1) = match; end %#ok<AGROW>
        end
        amp_indices = unique(amp_indices,'stable');
    end
    for ia = amp_indices
        item = struct('condition_index',ic,'amplitude_index',ia);
        plan(end+1) = item; %#ok<AGROW>
    end
end
end

function [columns,depth_channels] = resolve_channels(request,T)
saved = double(T.DepthChannel(:).');
if ischar(request) || (isstring(request) && isscalar(request))
    if strcmpi(strtrim(char(string(request))),'all')
        columns = 1:numel(saved);
    else
        error('MultiISIRaster:ChannelText', ...
            'Text Channels_To_Plot must be all; otherwise enter depth numbers.');
    end
elseif isnumeric(request)
    requested = unique(double(request(:).'),'stable');
    if isempty(requested) || any(mod(requested,1)~=0)
        error('MultiISIRaster:ChannelSelection','Invalid depth-channel list.');
    end
    [found,columns] = ismember(requested,saved);
    if any(~found)
        error('MultiISIRaster:UnavailableChannel', ...
            'Depth channel(s) not in ModelData: %s',num2str(requested(~found)));
    end
else
    error('MultiISIRaster:ChannelSelection','Invalid Channels_To_Plot setting.');
end
depth_channels = saved(columns);
end

function [columns,label] = channels_for_condition(base_columns,code,T, ...
        use_order_specific,reference_population)
% Intersect the user's channel selection with the response population that
% belongs to the plotted stimulation order.
if ~use_order_specific
    columns = base_columns;
    label = 'Response population: union';
    return;
end

mask_A_to_B = logical(T.Responsive_A_to_B(:).');
mask_B_to_A = logical(T.Responsive_B_to_A(:).');
switch code
    case 'A_to_B'
        response_mask = mask_A_to_B;
        label = 'Response population: A->B';
    case 'B_to_A'
        response_mask = mask_B_to_A;
        label = 'Response population: B->A';
    otherwise
        switch lower(reference_population)
            case 'union'
                response_mask = mask_A_to_B | mask_B_to_A;
                label = 'Response population: order union';
            case 'a_to_b'
                response_mask = mask_A_to_B;
                label = 'Response population: A->B';
            case 'b_to_a'
                response_mask = mask_B_to_A;
                label = 'Response population: B->A';
        end
end
columns = base_columns(response_mask(base_columns));
end

function validate_amplitude_block(B,nChannels)
required = {'amplitude_uA','channel'};
require_fields(B,required,'condition amplitude block');
if numel(B.channel) ~= nChannels
    error('MultiISIRaster:ChannelCount', ...
        'Amplitude block channel count does not match ModelData.channels.');
end
end

function validate_channel_data(C)
required = {'n_trials','spike_times_ms'};
require_fields(C,required,'condition amplitude channel block');
if ~iscell(C.spike_times_ms) || numel(C.spike_times_ms) ~= C.n_trials
    error('MultiISIRaster:SpikeShape', ...
        'A channel spike_times_ms cell array has an invalid size.');
end
end

function spikes_out = restrict_spikes(spikes_in,window)
spikes_out = cell(size(spikes_in));
for k = 1:numel(spikes_in)
    tt = double(spikes_in{k}(:));
    spikes_out{k} = tt(tt >= window(1) & tt <= window(2));
end
end

function rate = calculate_psth(spikes,edges,smooth_bins)
counts = zeros(1,numel(edges)-1);
for k = 1:numel(spikes)
    counts = counts+histcounts(spikes{k},edges);
end
if isempty(spikes)
    rate = zeros(size(counts));
else
    bin_duration_s = diff(edges)/1000;
    rate = counts./(numel(spikes)*bin_duration_s);
end
if smooth_bins > 1
    rate = movmean(rate,smooth_bins,'Endpoints','shrink');
end
end

function limit = nice_rate_limit(max_rate)
if ~isfinite(max_rate) || max_rate <= 0
    limit = 10;
else
    limit = max(10,ceil(max_rate*1.10/10)*10);
end
end

function color = color_for_condition(code)
switch code
    case 'A',      color = [0.10 0.45 0.80];
    case 'B',      color = [0.85 0.33 0.10];
    case 'AB',     color = [0.15 0.65 0.30];
    case 'A_to_B', color = [0.60 0.25 0.70];
    case 'B_to_A', color = [0.15 0.65 0.70];
    otherwise
        error('MultiISIRaster:UnknownCondition','Unknown condition code: %s',code);
end
end

function text = pulse_description(condition)
parts = cell(1,numel(condition.electrode_order));
for k = 1:numel(parts)
    parts{k} = sprintf('%s at %g ms',condition.electrode_order{k}, ...
        condition.pulse_times_ms(k));
end
text = strjoin(parts,', ');
end

function tag = stimulation_channel_tag(depth_ch,D)
tag = '';
row = find(double(D.channels.DepthChannel)==depth_ch,1);
if isempty(row), return; end
isA = ismember('DistanceToA_um',D.channels.Properties.VariableNames) && ...
    abs(D.channels.DistanceToA_um(row)) <= 1e-9;
isB = ismember('DistanceToB_um',D.channels.Properties.VariableNames) && ...
    abs(D.channels.DistanceToB_um(row)) <= 1e-9;
if isA && isB
    tag = ' [STIM A/B]';
elseif isA
    tag = ' [STIM A]';
elseif isB
    tag = ' [STIM B]';
end
end

function validate_window(value,name)
if ~isnumeric(value) || numel(value)~=2 || any(~isfinite(value)) || ...
        value(2)<=value(1)
    error('MultiISIRaster:Window','%s must be [start end], with end>start.',name);
end
end

function require_fields(S,names,label)
for k = 1:numel(names)
    if ~isfield(S,names{k})
        error('MultiISIRaster:MissingField','%s.%s is missing.',label,names{k});
    end
end
end
