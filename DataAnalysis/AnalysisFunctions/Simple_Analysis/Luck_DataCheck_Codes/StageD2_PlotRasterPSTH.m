%% Stage D2: plot raster and PSTH for the exported modelling data
%
% PURPOSE
%   This script loads only the Stage D1 ModelData MAT file. It does not
%   require the original experiment folders, trigger files, sp_corr files,
%   Stage A/B/C results, or any external analysis functions.
%
% FIGURE ORGANIZATION
%   One figure is created for each selected condition x amplitude.
%   Every selected responding recording channel is shown as one tile.
%   Within each tile:
%       left y-axis  = PSTH firing rate (spikes/second)
%       right y-axis = trial-by-trial spike raster
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
model_data_file = ['/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/auto_pair_pipeline/auto_pair_results/Luke_ModelPackage/DX014_D3_Pair_A029_A030_ModelData.mat'];

% Choose 'all' or 'balanced'.
Trial_Mode = 'all';

% Conditions can be selected by code or index.
% Use 'all' for all five conditions, or examples such as:
%   Conditions_To_Plot = {'A','AB','A_to_B'};
%   Conditions_To_Plot = [1 3 4];
%
% Condition codes:
%   A, B, AB, A_to_B, B_to_A
Conditions_To_Plot = 'all';

% Empty means all amplitudes contained in the exported package.
% Example: Amplitudes_To_Plot = [5 10];
Amplitudes_To_Plot = [];

% Channels are selected using anatomical/depth-channel numbers shown in
% ModelData.channels.DepthChannel. Use 'all' for every exported responding
% channel, or enter a numeric list such as [1 3 5 8].
Channels_To_Plot = 'all';

% Display window in milliseconds relative to the first stimulation pulse.
% It must remain inside ModelData.metadata.stored_spike_window_ms.
raster_win_ms = [-50 80];

% PSTH settings. The firing rate is calculated from the selected trials.
% Smoothing uses a centred moving average and requires no toolbox.
psth_bin_ms = 1;
psth_smooth_ms = 5;

% Use one common PSTH y-axis maximum across all channels within each
% condition x amplitude figure. This makes channel responses comparable.
Use_Shared_PSTH_YAxis = true;

% Figure presentation only. No figures are saved by this script.
Close_Existing_Figures = true;
figure_position = [40 40 1600 900];

%% ========================== LOAD DATA =================================
if Close_Existing_Figures
    close all;
end

if ~isfile(model_data_file)
    error('StageD2:MissingModelData', ...
        'ModelData file not found:\n%s',model_data_file);
end

S = load(model_data_file,'ModelData');
if ~isfield(S,'ModelData')
    error('StageD2:InvalidModelData', ...
        'The selected MAT file does not contain ModelData.');
end
D = S.ModelData;

required_top_fields = {'format_name','metadata','amplitudes_uA', ...
    'condition_definitions','channels','conditions','pair_key', ...
    'trial_selection'};
require_fields(D,required_top_fields,'ModelData');
if ~strcmp(D.format_name,'SinglePairStimulationModelData')
    error('StageD2:UnexpectedFormat', ...
        'The selected MAT file is not a Stage D1 single-pair package.');
end
if numel(D.conditions) ~= 5 || numel(D.condition_definitions) ~= 5
    error('StageD2:ConditionCount','The package must contain five conditions.');
end

required_channel_variables = {'ChannelIndex','DepthChannel','IsStimA','IsStimB'};
for k = 1:numel(required_channel_variables)
    if ~ismember(required_channel_variables{k},D.channels.Properties.VariableNames)
        error('StageD2:MissingChannelVariable', ...
            'ModelData.channels is missing %s.',required_channel_variables{k});
    end
end

%% ====================== RESOLVE USER CHOICES ==========================
trial_mode = lower(strtrim(char(string(Trial_Mode))));
if ~ismember(trial_mode,{'all','balanced'})
    error('StageD2:InvalidTrialMode', ...
        'Trial_Mode must be either all or balanced.');
end

condition_indices = resolve_conditions(Conditions_To_Plot,D);
[amplitude_indices,selected_amplitudes] = ...
    resolve_amplitudes(Amplitudes_To_Plot,D.amplitudes_uA);
[channel_columns,selected_depth_channels] = ...
    resolve_channels(Channels_To_Plot,D.channels);

validate_window(raster_win_ms,'raster_win_ms');
stored_window = double(D.metadata.stored_spike_window_ms(:).');
if raster_win_ms(1) < stored_window(1) || raster_win_ms(2) > stored_window(2)
    error('StageD2:WindowOutsideStoredData', ...
        ['Requested raster window [%g %g] ms lies outside the stored ' ...
         'window [%g %g] ms.'],raster_win_ms,stored_window);
end
if ~isscalar(psth_bin_ms) || ~isfinite(psth_bin_ms) || psth_bin_ms <= 0
    error('StageD2:InvalidBin','psth_bin_ms must be one positive number.');
end
if ~isscalar(psth_smooth_ms) || ~isfinite(psth_smooth_ms) || ...
        psth_smooth_ms < 0
    error('StageD2:InvalidSmoothing', ...
        'psth_smooth_ms must be one nonnegative number.');
end

%% ====================== PSTH BIN DEFINITIONS ==========================
% The final bin is allowed to be shorter if the requested window width is
% not an exact multiple of psth_bin_ms. Its own duration is used when the
% rate is calculated.
edges = raster_win_ms(1):psth_bin_ms:raster_win_ms(2);
if isempty(edges) || edges(1) ~= raster_win_ms(1)
    edges = raster_win_ms(1);
end
if edges(end) < raster_win_ms(2)
    edges(end+1) = raster_win_ms(2);
end
centres = edges(1:end-1)+diff(edges)/2;
smooth_bins = max(1,round(psth_smooth_ms/psth_bin_ms));

condition_colors = [0.10 0.45 0.80; ... % A
                    0.85 0.33 0.10; ... % B
                    0.15 0.65 0.30; ... % A+B
                    0.60 0.25 0.70; ... % A->B
                    0.15 0.65 0.70];    % B->A

%% ======================= COMMAND-WINDOW SUMMARY =======================
nFigures = numel(condition_indices)*numel(amplitude_indices);
fprintf('\n============================================================\n');
fprintf('STAGE D2: RASTER + PSTH\n');
fprintf('============================================================\n');
fprintf('ModelData file: %s\n',model_data_file);
fprintf('Pair: %s\n',D.pair_key);
fprintf('Trial mode: %s\n',upper(trial_mode));
fprintf('Conditions: %s\n',strjoin({D.conditions(condition_indices).code},', '));
fprintf('Amplitudes: %s uA\n',num2str(selected_amplitudes));
fprintf('Depth channels: %s\n',num2str(selected_depth_channels));
fprintf('Raster window: [%g %g] ms\n',raster_win_ms);
fprintf('PSTH bin: %g ms; smoothing: %g ms (%d bins)\n', ...
    psth_bin_ms,psth_smooth_ms,smooth_bins);
fprintf('Expected figures: %d\n',nFigures);
fprintf('Figures will be displayed only. Files written: NONE\n');

if isfield(D,'qc') && isfield(D.qc,'single_trial_review_complete') && ...
        ~D.qc.single_trial_review_complete && ...
        any(ismember(condition_indices,[1 2]))
    fprintf(['NOTICE: A-alone and B-alone trials were exported before the ' ...
        'single-trial manual review was completed.\n']);
end

%% ========================== MAIN PLOTTING =============================
figure_counter = 0;
nSelectedChannels = numel(channel_columns);
nTileColumns = ceil(sqrt(nSelectedChannels));
nTileRows = ceil(nSelectedChannels/nTileColumns);

for icLocal = 1:numel(condition_indices)
    ic = condition_indices(icLocal);
    condition = D.conditions(ic);
    definition = D.condition_definitions(ic);

    for iaLocal = 1:numel(amplitude_indices)
        ia = amplitude_indices(iaLocal);
        block = condition.amplitude(ia);
        validate_amplitude_block(block,height(D.channels));

        trial_indices = select_trial_indices(block,trial_mode);
        nTrials = numel(trial_indices);
        if nTrials < 1
            error('StageD2:NoSelectedTrials', ...
                'No %s trials are available for %s at %g uA.', ...
                trial_mode,condition.code,selected_amplitudes(iaLocal));
        end

        % Precalculate every channel PSTH before plotting. This allows a
        % shared y-axis to be selected without drawing each curve twice.
        channel_spikes = cell(1,nSelectedChannels);
        channel_rates = cell(1,nSelectedChannels);
        individual_ymax = zeros(1,nSelectedChannels);
        shared_ymax = 0;

        for jcLocal = 1:nSelectedChannels
            jc = channel_columns(jcLocal);
            spikes = block.spike_times_ms(trial_indices,jc);
            spikes = restrict_spikes(spikes,raster_win_ms);
            rate = calculate_psth(spikes,edges,smooth_bins);
            channel_spikes{jcLocal} = spikes;
            channel_rates{jcLocal} = rate;
            if isempty(rate)
                this_max = 0;
            else
                this_max = max(rate);
            end
            individual_ymax(jcLocal) = nice_rate_limit(this_max);
            shared_ymax = max(shared_ymax,this_max);
        end
        shared_ymax = nice_rate_limit(shared_ymax);

        figure_counter = figure_counter+1;
        figure_name = sprintf('D2 %s %.6g uA %s trials', ...
            condition.code,selected_amplitudes(iaLocal),trial_mode);
        figure('Color','w','Name',figure_name,'Position',figure_position);
        tl = tiledlayout(nTileRows,nTileColumns, ...
            'TileSpacing','compact','Padding','compact');

        pulse_text = pulse_description(definition);
        title(tl,sprintf('%s | %s | %.6g uA | %s trials (n=%d) | %s', ...
            D.pair_key,condition.label,selected_amplitudes(iaLocal), ...
            upper(trial_mode),nTrials,pulse_text), ...
            'Interpreter','none','FontWeight','bold');

        for jcLocal = 1:nSelectedChannels
            jc = channel_columns(jcLocal);
            depth_ch = selected_depth_channels(jcLocal);
            spikes = channel_spikes{jcLocal};
            rate = channel_rates{jcLocal};

            ax = nexttile(tl); hold(ax,'on');

            % ----- Left axis: population firing rate across trials -----
            yyaxis(ax,'left');
            plot(ax,centres,rate,'Color',condition_colors(ic,:), ...
                'LineWidth',1.4);
            if Use_Shared_PSTH_YAxis
                ylim(ax,[0 shared_ymax]);
            else
                ylim(ax,[0 individual_ymax(jcLocal)]);
            end
            ax.YAxis(1).Color = condition_colors(ic,:);
            if mod(jcLocal-1,nTileColumns) == 0
                ylabel(ax,'Rate (sp/s)');
            end

            % ----- Right axis: one raster row per selected trial --------
            yyaxis(ax,'right');
            for it = 1:nTrials
                tt = spikes{it};
                if ~isempty(tt)
                    plot(ax,tt,it*ones(size(tt)),'.k','MarkerSize',4);
                end
            end
            ylim(ax,[0 nTrials+1]);
            ax.YAxis(2).Color = [0 0 0];
            if mod(jcLocal,nTileColumns) == 0 || jcLocal == nSelectedChannels
                ylabel(ax,'Trial');
                if nTrials == 1
                    yticks(ax,1);
                else
                    yticks(ax,[1 nTrials]);
                end
            else
                set(ax,'YTick',[]);
            end

            % Pulse 1 is always aligned to zero. Sequential conditions
            % additionally receive a marker at the second-pulse time.
            xline(ax,0,'r--','LineWidth',1);
            later_pulses = definition.pulse_times_ms( ...
                definition.pulse_times_ms > 0);
            for ip = 1:numel(later_pulses)
                xline(ax,later_pulses(ip),'k:','LineWidth',1.1);
            end
            xlim(ax,raster_win_ms);

            channel_tag = stimulation_channel_tag(D.channels(jc,:));
            title(ax,sprintf('Rec Ch %d%s',depth_ch,channel_tag), ...
                'Interpreter','none','FontSize',10,'FontWeight','bold');
            if ceil(jcLocal/nTileColumns) == nTileRows
                xlabel(ax,'Time (ms)');
            end
            box(ax,'off');
        end

        fprintf('Figure %d/%d: %s | %g uA | %d %s trials\n', ...
            figure_counter,nFigures,condition.code, ...
            selected_amplitudes(iaLocal),nTrials,trial_mode);

        fprintf('Figure %d/%d: %s | %g uA \n', ...
            figure_counter,nFigures,condition.code, ...
            selected_amplitudes(iaLocal));
        drawnow;
    end
end

fprintf('\n============================================================\n');
fprintf('STAGE D2 COMPLETE\n');
fprintf('Figures displayed: %d\n',figure_counter);
fprintf('Files written: NONE\n');
fprintf('============================================================\n\n');

%% =========================== FUNCTIONS ================================
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
        error('StageD2:InvalidConditionIndices', ...
            'Numeric condition indices must be integers from 1 to %d.',numel(codes));
    end
else
    error('StageD2:InvalidConditions','Invalid Conditions_To_Plot setting.');
end
end

function indices = match_condition_codes(requested,codes)
indices = zeros(1,numel(requested));
for k = 1:numel(requested)
    match = find(strcmpi(strtrim(requested{k}),codes),1);
    if isempty(match)
        error('StageD2:UnknownCondition', ...
            'Unknown condition code: %s',requested{k});
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
            error('StageD2:UnavailableAmplitude', ...
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
        error('StageD2:InvalidChannelText', ...
            'Text Channels_To_Plot must be all. Otherwise use depth-channel numbers.');
    end
elseif isnumeric(request)
    requested = unique(double(request(:).'),'stable');
    if isempty(requested) || any(mod(requested,1) ~= 0)
        error('StageD2:InvalidChannels','Depth-channel selection is invalid.');
    end
    [found,columns] = ismember(requested,saved_depth);
    if any(~found)
        error('StageD2:UnavailableChannel', ...
            'Depth channel(s) not present in ModelData: %s', ...
            num2str(requested(~found)));
    end
else
    error('StageD2:InvalidChannels','Invalid Channels_To_Plot setting.');
end
depth_channels = saved_depth(columns);
end

function validate_amplitude_block(B,nChannels)
required = {'n_trials_all','n_trials_balanced','all_trial_indices', ...
    'balanced_trial_indices','spike_times_ms'};
require_fields(B,required,'condition amplitude block');
if ~iscell(B.spike_times_ms) || size(B.spike_times_ms,1) ~= B.n_trials_all || ...
        size(B.spike_times_ms,2) ~= nChannels
    error('StageD2:SpikeShape','The exported spike_times_ms shape is invalid.');
end
if numel(B.all_trial_indices) ~= B.n_trials_all || ...
        numel(B.balanced_trial_indices) ~= B.n_trials_balanced
    error('StageD2:TrialIndexCount','Exported trial-index counts are inconsistent.');
end
end

function indices = select_trial_indices(B,mode)
switch mode
    case 'all'
        indices = double(B.all_trial_indices(:));
    case 'balanced'
        indices = double(B.balanced_trial_indices(:));
end
if any(indices < 1) || any(indices > B.n_trials_all) || ...
        any(mod(indices,1) ~= 0) || numel(unique(indices)) ~= numel(indices)
    error('StageD2:InvalidTrialIndices', ...
        'The selected trial indices in ModelData are invalid.');
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

function text_value = pulse_description(definition)
electrodes = definition.electrode_order;
times = definition.pulse_times_ms;
parts = cell(1,numel(electrodes));
for k = 1:numel(electrodes)
    parts{k} = sprintf('%s at %g ms',electrodes{k},times(k));
end
text_value = strjoin(parts,', ');
end

function validate_window(value,name)
if ~isnumeric(value) || numel(value) ~= 2 || any(~isfinite(value)) || ...
        value(2) <= value(1)
    error('StageD2:InvalidWindow','%s must be [start end], with end>start.',name);
end
end

function require_fields(S,names,label)
for k = 1:numel(names)
    if ~isfield(S,names{k})
        error('StageD2:MissingField','%s.%s is missing.',label,names{k});
    end
end
end
