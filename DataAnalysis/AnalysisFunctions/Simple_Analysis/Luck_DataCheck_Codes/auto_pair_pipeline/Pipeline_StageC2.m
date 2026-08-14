%% Stage C2: visualize one prepared stimulation pair
% Loads only a Stage C1 MAT result. Displays figures but does not save files
% and does not modify the Stage C1 result or any experiment data.

clear;
clc;

%% ====================== PIPELINE SETTINGS =============================
config_file = getenv('AUTO_PAIR_QC_CONFIG');
if isempty(config_file) || ~isfile(config_file)
    error('AutoPairQC:MissingConfig', ...
        'Run Pipeline_StageC2 through Run_AutoPairQC.m.');
end
tmp_cfg = load(config_file,'PipelineConfig');
PipelineConfig = tmp_cfg.PipelineConfig;
stage_c1_file = PipelineConfig.current_stage_c1_file;
Plot_Figures = PipelineConfig.Plot_Figures;
Channels_To_Plot = PipelineConfig.Channels_To_Plot;
Amplitudes_To_Plot = PipelineConfig.Amplitudes_To_Plot;
raster_win_ms = PipelineConfig.raster_win_ms;
baseline_win_ms = PipelineConfig.baseline_win_ms;
count_win_ms = PipelineConfig.count_win_ms;
psth_bin_ms = PipelineConfig.psth_bin_ms;
psth_smooth_ms = PipelineConfig.psth_smooth_ms;
baseline_correct_amplitude_curves = ...
    PipelineConfig.baseline_correct_amplitude_curves;
channels_per_curve_figure = PipelineConfig.channels_per_curve_figure;
figure_position_large = PipelineConfig.figure_position_large;
figure_position_medium = PipelineConfig.figure_position_medium;

%% ========================= LOAD AND CHECK =============================
if ~isfile(stage_c1_file)
    error('StageC2:MissingFile','Stage C1 file not found:\n%s',stage_c1_file);
end

S = load(stage_c1_file,'StageC1');
if ~isfield(S,'StageC1')
    error('StageC2:InvalidFile','The selected MAT file does not contain StageC1.');
end
C = S.StageC1;

required_fields = {'Data','ChannelInfo','ConditionTrialCounts', ...
    'amplitudes_uA','condition_order','selected_depth_channels','windows'};
for k = 1:numel(required_fields)
    if ~isfield(C,required_fields{k})
        error('StageC2:MissingField','StageC1.%s is missing.',required_fields{k});
    end
end
condition_order = C.condition_order;
nCondition = numel(condition_order);
if size(C.Data,3) ~= nCondition || nCondition < 4 || ...
        ~all(ismember({'A','B','AB'},condition_order)) || ...
        ~any(ismember({'A_to_B','B_to_A'},condition_order))
    error('StageC2:ConditionCount', ...
        'C2 requires A, B, AB, and at least one sequential condition.');
end

validate_window(raster_win_ms,'raster_win_ms');
validate_window(baseline_win_ms,'baseline_win_ms');
validate_window(count_win_ms,'count_win_ms');
saved_raster = C.windows.raster_ms;
if raster_win_ms(1) < saved_raster(1) || raster_win_ms(2) > saved_raster(2)
    error('StageC2:RasterOutsideSavedRange', ...
        'raster_win_ms must remain inside the C1 range [%g %g] ms.',saved_raster);
end
if baseline_win_ms(1) < saved_raster(1) || baseline_win_ms(2) > saved_raster(2) || ...
        count_win_ms(1) < saved_raster(1) || count_win_ms(2) > saved_raster(2)
    error('StageC2:AnalysisWindowOutsideSavedRange', ...
        'The baseline and count windows must remain inside the C1 raster range.');
end
if ~isscalar(psth_bin_ms) || psth_bin_ms <= 0 || ...
        ~isscalar(psth_smooth_ms) || psth_smooth_ms < 0
    error('StageC2:InvalidPSTHSettings', ...
        'psth_bin_ms must be >0 and psth_smooth_ms must be >=0.');
end

valid_figure_ids = 1:5;
Plot_Figures = unique(double(Plot_Figures(:).'),'stable');
if any(~ismember(Plot_Figures,valid_figure_ids))
    error('StageC2:InvalidFigureID', ...
        'Plot_Figures can contain only 1, 2, 3, 4, and 5.');
end

[selected_rows,selected_channels] = resolve_channels(Channels_To_Plot,C);
[selected_amp_indices,selected_amps] = resolve_amplitudes(Amplitudes_To_Plot,C);

[condition_names,condition_short,condition_colors] = ...
    condition_plot_style(condition_order,C);
idx_A = find(strcmp(condition_order,'A'),1);
idx_B = find(strcmp(condition_order,'B'),1);
linear_color = [0.10 0.10 0.10];
dot_labels = [condition_short {'Linear'}];

fprintf('\n============================================================\n');
fprintf('STAGE C2: PAIR VISUALIZATION\n');
fprintf('============================================================\n');
fprintf('C1 source: %s\n',stage_c1_file);
fprintf('Pair: %s\n',C.pair_key);
fprintf('Conditions: %s\n',strjoin(condition_short,', '));
fprintf('Amplitudes: %s uA\n',num2str(selected_amps));
fprintf('Depth channels: %s\n',num2str(selected_channels));
fprintf('Figures requested: %s\n',num2str(Plot_Figures));
fprintf('Raster window: [%g %g] ms\n',raster_win_ms);
fprintf('Baseline window: [%g %g) ms\n',baseline_win_ms);
fprintf('Spike-count window: [%g %g) ms\n',count_win_ms);
fprintf('Figures will be displayed only. No files will be written.\n');

nFigExpected = estimate_figure_count(Plot_Figures,numel(selected_rows), ...
    channels_per_curve_figure);
fprintf('Expected number of figures: %d\n',nFigExpected);
if isfield(C,'single_qc') && ~C.single_qc.review_complete
    fprintf(['NOTICE: Single-stimulation trials are not yet marked as fully ' ...
        'reviewed. Inspect A-alone and B-alone rasters carefully.\n']);
end

%% ================= FIGURE 1: TRIAL OVERVIEW ==========================
if ismember(1,Plot_Figures)
    counts = table2array(C.ConditionTrialCounts(selected_amp_indices,:));
    figure('Color','w','Name','C2 Figure 1 - Trial overview', ...
        'Position',figure_position_medium);
    ax = axes; hold(ax,'on');
    b = bar(ax,selected_amps,counts,'grouped');
    for ic = 1:nCondition
        b(ic).FaceColor = condition_colors(ic,:);
    end
    xlabel(ax,'Stimulation amplitude (uA)','FontWeight','bold');
    ylabel(ax,'Number of retained trials','FontWeight','bold');
    title(ax,sprintf('%s: retained trials (no balancing)',C.pair_key), ...
        'Interpreter','none','FontWeight','bold');
    legend(ax,condition_names,'Location','bestoutside','Box','off');
    xticks(ax,selected_amps);
    grid(ax,'on'); box(ax,'off');

    info_text = sprintf(['Selected channels: %d    Pair spacing: %.0f um    ' ...
        'Single review complete: %s'],numel(selected_channels), ...
        C.stimulation.pair_distance_um,logical_text(C.single_qc.review_complete));
    subtitle(ax,info_text,'Interpreter','none');
end

%% ================= FIGURE 2: RASTER + PSTH ==========================
if ismember(2,Plot_Figures)
    edges = raster_win_ms(1):psth_bin_ms:raster_win_ms(2);
    if edges(end) < raster_win_ms(2)
        edges(end+1) = raster_win_ms(2);
    end
    centres = edges(1:end-1)+diff(edges)/2;
    smooth_bins = max(1,round(psth_smooth_ms/psth_bin_ms));

    for ir = 1:numel(selected_rows)
        row = selected_rows(ir);
        depth_ch = selected_channels(ir);
        stim_tag = stimulation_channel_tag(depth_ch,C);
        figure('Color','w', ...
            'Name',sprintf('C2 Figure 2 - Ch %d raster PSTH',depth_ch), ...
            'Position',figure_position_large);
        tl = tiledlayout(numel(selected_amp_indices),nCondition, ...
            'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('%s | Recording Ch %d%s | Raster + PSTH', ...
            C.pair_key,depth_ch,stim_tag),'Interpreter','none','FontWeight','bold');

        for iaLocal = 1:numel(selected_amp_indices)
            ia = selected_amp_indices(iaLocal);
            max_rate_this_amp = 0;
            cached_rates = cell(1,nCondition);
            cached_spikes = cell(1,nCondition);
            for ic = 1:nCondition
                D = C.Data(row,ia,ic);
                spikes = restrict_spikes(D.spike_times_ms,raster_win_ms);
                cached_spikes{ic} = spikes;
                rate = calculate_psth(spikes,edges,smooth_bins);
                cached_rates{ic} = rate;
                if ~isempty(rate)
                    max_rate_this_amp = max(max_rate_this_amp,max(rate));
                end
            end
            psth_ymax = max(10,ceil(max_rate_this_amp*1.10/10)*10);

            for ic = 1:nCondition
                ax = nexttile(tl); hold(ax,'on');
                spikes = cached_spikes{ic};
                rate = cached_rates{ic};
                nTrials = numel(spikes);

                yyaxis(ax,'left');
                plot(ax,centres,rate,'Color',condition_colors(ic,:), ...
                    'LineWidth',1.3);
                ylim(ax,[0 psth_ymax]);
                if ic == 1
                    ylabel(ax,sprintf('%g uA\nRate (sp/s)',selected_amps(iaLocal)), ...
                        'FontWeight','bold');
                end
                ax.YColor = condition_colors(ic,:);

                yyaxis(ax,'right');
                for it = 1:nTrials
                    tt = spikes{it};
                    if ~isempty(tt)
                        plot(ax,tt,it*ones(size(tt)),'.k','MarkerSize',4);
                    end
                end
                ylim(ax,[0 nTrials+1]);
                set(ax,'YTick',[]);
                ax.YColor = [0 0 0];
                xline(ax,0,'r--','LineWidth',1);
                if ismember(condition_order{ic},{'A_to_B','B_to_A'})
                    xline(ax,C.Data(row,ia,ic).second_pulse_ms,'k:', ...
                        'LineWidth',1);
                end
                xlim(ax,raster_win_ms);
                title(ax,sprintf('%s | n=%d',condition_short{ic},nTrials), ...
                    'Interpreter','none','FontSize',10);
                if iaLocal == numel(selected_amp_indices)
                    xlabel(ax,'Time (ms)');
                end
                box(ax,'off');
            end
        end
        drawnow;
    end
end

%% =============== FIGURE 3: TRIAL COUNT DISTRIBUTIONS ================
if ismember(3,Plot_Figures)
    for ir = 1:numel(selected_rows)
        row = selected_rows(ir);
        depth_ch = selected_channels(ir);
        stim_tag = stimulation_channel_tag(depth_ch,C);
        figure('Color','w', ...
            'Name',sprintf('C2 Figure 3 - Ch %d trial counts',depth_ch), ...
            'Position',figure_position_medium);
        tl = tiledlayout(1,numel(selected_amp_indices), ...
            'TileSpacing','compact','Padding','compact');
        title(tl,sprintf('%s | Recording Ch %d%s | Per-trial counts [%g, %g) ms', ...
            C.pair_key,depth_ch,stim_tag,count_win_ms), ...
            'Interpreter','none','FontWeight','bold');

        yMax = 1;
        all_counts = cell(numel(selected_amp_indices),nCondition);
        dot_linear_mean = nan(numel(selected_amp_indices),1);
        dot_linear_sem = nan(numel(selected_amp_indices),1);
        for iaLocal = 1:numel(selected_amp_indices)
            ia = selected_amp_indices(iaLocal);
            for ic = 1:nCondition
                metrics = recalculate_metrics(C.Data(row,ia,ic).spike_times_ms, ...
                    baseline_win_ms,count_win_ms);
                all_counts{iaLocal,ic} = metrics.post_count;
                if ~isempty(metrics.post_count)
                    yMax = max(yMax,max(metrics.post_count));
                end
            end
            mA = recalculate_metrics(C.Data(row,ia,idx_A).spike_times_ms, ...
                baseline_win_ms,count_win_ms);
            mB = recalculate_metrics(C.Data(row,ia,idx_B).spike_times_ms, ...
                baseline_win_ms,count_win_ms);
            % Figure 3 displays raw post-window counts. On that same raw
            % scale, subtract one shared baseline copy from A+B.
            contribution_A = mA.post_count-0.5*mA.expected_baseline_post_count;
            contribution_B = mB.post_count-0.5*mB.expected_baseline_post_count;
            dot_linear_mean(iaLocal) = mean(contribution_A)+mean(contribution_B);
            dot_linear_sem(iaLocal) = sqrt( ...
                (std(contribution_A,0)/sqrt(numel(contribution_A)))^2 + ...
                (std(contribution_B,0)/sqrt(numel(contribution_B)))^2);
            yMax = max(yMax,dot_linear_mean(iaLocal)+dot_linear_sem(iaLocal));
        end
        yMax = max(1,ceil(yMax+0.5));

        for iaLocal = 1:numel(selected_amp_indices)
            ax = nexttile(tl); hold(ax,'on');
            for ic = 1:nCondition
                values = all_counts{iaLocal,ic};
                jitter = deterministic_jitter(numel(values),0.24);
                scatter(ax,ic+jitter,values,18,condition_colors(ic,:), ...
                    'filled','MarkerFaceAlpha',0.35,'MarkerEdgeAlpha',0.35);
                mu = mean(values);
                sem = std(values,0)/sqrt(numel(values));
                errorbar(ax,ic,mu,sem,'Color',[0 0 0], ...
                    'LineStyle','none','LineWidth',1.3,'CapSize',8);
                plot(ax,ic,mu,'o','MarkerSize',6,'MarkerFaceColor', ...
                    condition_colors(ic,:),'MarkerEdgeColor','k');
            end
            % A and B are independent trial groups, so there is no genuine
            % trial-by-trial linear prediction to scatter. Show its mean
            % and propagated SEM without inventing paired observations.
            linear_x = nCondition+1;
            errorbar(ax,linear_x,dot_linear_mean(iaLocal),dot_linear_sem(iaLocal), ...
                'Color',linear_color,'LineStyle','none','LineWidth',1.6, ...
                'CapSize',9);
            plot(ax,linear_x,dot_linear_mean(iaLocal),'d','MarkerSize',7, ...
                'MarkerFaceColor','w','MarkerEdgeColor',linear_color, ...
                'LineWidth',1.4);
            xlim(ax,[0.5 nCondition+1.5]); ylim(ax,[-0.2 yMax]);
            xticks(ax,1:nCondition+1); xticklabels(ax,dot_labels);
            xtickangle(ax,25);
            ylabel(ax,'Spike count / trial');
            title(ax,sprintf('%g uA',selected_amps(iaLocal)));
            grid(ax,'on'); box(ax,'off');
        end
        drawnow;
    end
end

%% ===== FIGURE 5: AVERAGE CURVES ACROSS RESPONDING CHANNELS ===========
if ismember(5,Plot_Figures)
    fprintf('Creating Figure 5: average across responding channels...\n');
    if ~ismember('IsRespondingUnion',C.ChannelInfo.Properties.VariableNames)
        error('StageC2:MissingRespondingFlag', ...
            'ChannelInfo does not contain IsRespondingUnion.');
    end
    population_rows = find(C.ChannelInfo.IsRespondingUnion(:)).';
    if isempty(population_rows)
        error('StageC2:NoRespondingPopulation', ...
            'No responding-union channels are stored in the C1 file.');
    end

    nPopulation = numel(population_rows);
    channel_condition_means = nan(nPopulation,nCondition,numel(selected_amp_indices));
    channel_linear_means = nan(nPopulation,numel(selected_amp_indices));

    for ip = 1:nPopulation
        row = population_rows(ip);
        for iaLocal = 1:numel(selected_amp_indices)
            ia = selected_amp_indices(iaLocal);
            for ic = 1:nCondition
                metrics = recalculate_metrics(C.Data(row,ia,ic).spike_times_ms, ...
                    baseline_win_ms,count_win_ms);
                if baseline_correct_amplitude_curves
                    values = metrics.evoked_count;
                else
                    values = metrics.post_count;
                end
                channel_condition_means(ip,ic,iaLocal) = mean(values);
            end

            mA = recalculate_metrics(C.Data(row,ia,idx_A).spike_times_ms, ...
                baseline_win_ms,count_win_ms);
            mB = recalculate_metrics(C.Data(row,ia,idx_B).spike_times_ms, ...
                baseline_win_ms,count_win_ms);
            if baseline_correct_amplitude_curves
                channel_linear_means(ip,iaLocal) = ...
                    mean(mA.evoked_count)+mean(mB.evoked_count);
            else
                contribution_A = mA.post_count- ...
                    0.5*mA.expected_baseline_post_count;
                contribution_B = mB.post_count- ...
                    0.5*mB.expected_baseline_post_count;
                channel_linear_means(ip,iaLocal) = ...
                    mean(contribution_A)+mean(contribution_B);
            end
        end
    end

    population_mean = squeeze(mean(channel_condition_means,1,'omitnan'));
    population_n = squeeze(sum(~isnan(channel_condition_means),1));
    population_sem = squeeze(std(channel_condition_means,0,1,'omitnan')) ./ ...
        sqrt(population_n);
    population_linear_mean = mean(channel_linear_means,1,'omitnan');
    population_linear_sem = std(channel_linear_means,0,1,'omitnan') ./ ...
        sqrt(sum(~isnan(channel_linear_means),1));

    fig_population = figure('Color','w', ...
        'Name','C2 Figure 5 - responding-channel average', ...
        'Position',figure_position_medium);
    ax = axes; hold(ax,'on');
    curve_handles = gobjects(1,nCondition+1);
    for ic = 1:nCondition
        curve_handles(ic) = errorbar(ax,selected_amps, ...
            population_mean(ic,:),population_sem(ic,:),'-o', ...
            'Color',condition_colors(ic,:), ...
            'MarkerFaceColor',condition_colors(ic,:), ...
            'LineWidth',1.7,'MarkerSize',6,'CapSize',7);
    end
    curve_handles(nCondition+1) = errorbar(ax,selected_amps,population_linear_mean, ...
        population_linear_sem,'--d','Color',linear_color, ...
        'MarkerFaceColor','w','LineWidth',2,'MarkerSize',7,'CapSize',7);
    if baseline_correct_amplitude_curves
        yline(ax,0,':','Color',[0.45 0.45 0.45], ...
            'HandleVisibility','off');
        y_label = 'Mean baseline-corrected spike count / trial';
    else
        y_label = 'Mean raw spike count / trial';
    end
    xlabel(ax,'Amplitude (uA)','FontWeight','bold');
    ylabel(ax,y_label,'FontWeight','bold');
    title(ax,sprintf('%s | Average across %d responding channels', ...
        C.pair_key,nPopulation),'Interpreter','none','FontWeight','bold');
    subtitle(ax,'Error bars are SEM across channels, not across animals.');
    xticks(ax,selected_amps);
    legend(ax,curve_handles,[condition_names {'Linear A+B'}], ...
        'Location','bestoutside','Box','off');
    grid(ax,'on'); box(ax,'off');
    drawnow;
    fprintf('Figure 5 created using %d responding channels.\n',nPopulation);
end

%% ============ FIGURE 4: AMPLITUDE-RESPONSE CURVES ===================
if ismember(4,Plot_Figures)
    nPages = ceil(numel(selected_rows)/channels_per_curve_figure);
    for page = 1:nPages
        first_local = (page-1)*channels_per_curve_figure+1;
        last_local = min(page*channels_per_curve_figure,numel(selected_rows));
        page_rows = first_local:last_local;

        figure('Color','w', ...
            'Name',sprintf('C2 Figure 4 - amplitude curves %d',page), ...
            'Position',figure_position_medium);
        tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        if baseline_correct_amplitude_curves
            metric_title = 'baseline-corrected spike count / trial';
        else
            metric_title = 'raw spike count / trial';
        end
        title(tl,sprintf('%s | Amplitude-response curves | %s', ...
            C.pair_key,metric_title),'Interpreter','none','FontWeight','bold');

        for q = 1:numel(page_rows)
            ir = page_rows(q);
            row = selected_rows(ir);
            depth_ch = selected_channels(ir);
            ax = nexttile(tl); hold(ax,'on');

            means = nan(nCondition,numel(selected_amp_indices));
            sems = nan(nCondition,numel(selected_amp_indices));
            linear_mean = nan(1,numel(selected_amp_indices));
            linear_sem = nan(1,numel(selected_amp_indices));

            for iaLocal = 1:numel(selected_amp_indices)
                ia = selected_amp_indices(iaLocal);
                for ic = 1:nCondition
                    metrics = recalculate_metrics( ...
                        C.Data(row,ia,ic).spike_times_ms, ...
                        baseline_win_ms,count_win_ms);
                    if baseline_correct_amplitude_curves
                        values = metrics.evoked_count;
                    else
                        values = metrics.post_count;
                    end
                    means(ic,iaLocal) = mean(values);
                    sems(ic,iaLocal) = std(values,0)/sqrt(numel(values));
                end

                % Linear prediction is defined from independent A-alone and
                % B-alone trials. Use evoked responses to avoid adding the
                % baseline twice.
                mA = recalculate_metrics(C.Data(row,ia,idx_A).spike_times_ms, ...
                    baseline_win_ms,count_win_ms);
                mB = recalculate_metrics(C.Data(row,ia,idx_B).spike_times_ms, ...
                    baseline_win_ms,count_win_ms);
                evA = mA.evoked_count;
                evB = mB.evoked_count;
                linear_evoked_mean = mean(evA)+mean(evB);
                linear_evoked_sem = sqrt((std(evA,0)/sqrt(numel(evA)))^2 + ...
                    (std(evB,0)/sqrt(numel(evB)))^2);

                if baseline_correct_amplitude_curves
                    linear_mean(iaLocal) = linear_evoked_mean;
                else
                    % For a raw-count display, add back one pooled estimate
                    % of the expected spontaneous count in the post window.
                    pooled_base = [mA.expected_baseline_post_count; ...
                        mB.expected_baseline_post_count];
                    linear_mean(iaLocal) = linear_evoked_mean+mean(pooled_base);
                end
                linear_sem(iaLocal) = linear_evoked_sem;
            end

            curve_handles = gobjects(1,nCondition+1);
            for ic = 1:nCondition
                curve_handles(ic) = errorbar(ax,selected_amps,means(ic,:), ...
                    sems(ic,:),'-o','Color',condition_colors(ic,:), ...
                    'MarkerFaceColor',condition_colors(ic,:), ...
                    'LineWidth',1.4,'MarkerSize',5,'CapSize',6);
            end
            curve_handles(nCondition+1) = errorbar(ax,selected_amps,linear_mean,linear_sem, ...
                '--d','Color',linear_color,'MarkerFaceColor','w', ...
                'LineWidth',1.7,'MarkerSize',6,'CapSize',6);
            if baseline_correct_amplitude_curves
                yline(ax,0,':','Color',[0.45 0.45 0.45], ...
                    'HandleVisibility','off');
            end
            xlabel(ax,'Amplitude (uA)');
            ylabel(ax,metric_title);
            title(ax,sprintf('Recording Ch %d%s',depth_ch, ...
                stimulation_channel_tag(depth_ch,C)),'Interpreter','none');
            xticks(ax,selected_amps);
            grid(ax,'on'); box(ax,'off');
            if q == 1
                legend(ax,curve_handles,[condition_names {'Linear A+B'}], ...
                    'Location','best','Box','off');
            end
        end
        drawnow;
    end
end

% Figure 5 is calculated before the Figure 4 pages. Bring the population
% average back to the front so it is easy to find after a full run.
if ismember(5,Plot_Figures) && exist('fig_population','var') && ...
        isgraphics(fig_population)
    figure(fig_population);
    drawnow;
end

fprintf('\n============================================================\n');
fprintf('STAGE C2 COMPLETE\n');
fprintf('Figures displayed: %d expected\n',nFigExpected);
fprintf('Files written: NONE\n');
fprintf('Experiment files modified: NO\n');
fprintf('============================================================\n\n');

%% =========================== FUNCTIONS ================================
function validate_window(value,name)
if ~isnumeric(value) || numel(value) ~= 2 || any(~isfinite(value)) || ...
        value(2) <= value(1)
    error('StageC2:InvalidWindow','%s must be [start end], with end>start.',name);
end
end

function [rows,channels] = resolve_channels(request,C)
saved = double(C.selected_depth_channels(:).');
if ischar(request) || (isstring(request) && isscalar(request))
    mode = lower(strtrim(char(request)));
    switch mode
        case 'all'
            rows = 1:numel(saved);
        case 'responding'
            if ~ismember('IsRespondingUnion',C.ChannelInfo.Properties.VariableNames)
                error('StageC2:MissingRespondingFlag', ...
                    'ChannelInfo does not contain IsRespondingUnion.');
            end
            rows = find(C.ChannelInfo.IsRespondingUnion(:)).';
        otherwise
            error('StageC2:InvalidChannelMode', ...
                'Channels_To_Plot must be all, responding, or a numeric list.');
    end
elseif isnumeric(request)
    requested = unique(double(request(:).'),'stable');
    if isempty(requested) || any(mod(requested,1) ~= 0)
        error('StageC2:InvalidChannels','The numeric channel list is invalid.');
    end
    [tf,rows] = ismember(requested,saved);
    if any(~tf)
        error('StageC2:UnavailableChannel', ...
            'Requested depth channel(s) not stored by C1: %s', ...
            num2str(requested(~tf)));
    end
else
    error('StageC2:InvalidChannels','Invalid Channels_To_Plot setting.');
end
if isempty(rows)
    error('StageC2:NoChannels','No channels were selected.');
end
channels = saved(rows);
end

function [indices,amps] = resolve_amplitudes(request,C)
saved = double(C.amplitudes_uA(:).');
if isempty(request)
    indices = 1:numel(saved);
else
    requested = unique(double(request(:).'),'stable');
    indices = zeros(size(requested));
    for k = 1:numel(requested)
        match = find(abs(saved-requested(k)) <= 1e-6,1);
        if isempty(match)
            error('StageC2:UnavailableAmplitude', ...
                'Amplitude %g uA was not stored by C1.',requested(k));
        end
        indices(k) = match;
    end
end
amps = saved(indices);
end

function tag = stimulation_channel_tag(depth_ch,C)
isA = depth_ch == C.stimulation.A.depth_channel;
isB = depth_ch == C.stimulation.B.depth_channel;
if isA && isB
    tag = ' [STIM A/B]';
elseif isA
    tag = ' [STIM A]';
elseif isB
    tag = ' [STIM B]';
else
    tag = '';
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
    bin_seconds = diff(edges)/1000;
    rate = counts./(numel(spikes)*bin_seconds);
end
if smooth_bins > 1
    rate = movmean(rate,smooth_bins,'Endpoints','shrink');
end
end

function M = recalculate_metrics(spikes,baseline_window,count_window)
n = numel(spikes);
baseline_count = zeros(n,1);
post_count = zeros(n,1);
base_duration = diff(baseline_window);
post_duration = diff(count_window);
for k = 1:n
    tt = double(spikes{k}(:));
    baseline_count(k) = sum(tt >= baseline_window(1) & tt < baseline_window(2));
    post_count(k) = sum(tt >= count_window(1) & tt < count_window(2));
end
expected_baseline_post_count = baseline_count*(post_duration/base_duration);
M = struct('baseline_count',baseline_count,'post_count',post_count, ...
    'expected_baseline_post_count',expected_baseline_post_count, ...
    'evoked_count',post_count-expected_baseline_post_count);
end

function jitter = deterministic_jitter(n,width)
if n <= 1
    jitter = zeros(n,1);
else
    jitter = linspace(-width,width,n).';
end
end

function n = estimate_figure_count(ids,nChannels,channelsPerCurveFigure)
n = 0;
if ismember(1,ids), n = n+1; end
if ismember(2,ids), n = n+nChannels; end
if ismember(3,ids), n = n+nChannels; end
if ismember(4,ids), n = n+ceil(nChannels/channelsPerCurveFigure); end
if ismember(5,ids), n = n+1; end
end

function [names,short_names,colors] = condition_plot_style(codes,C)
names = cell(size(codes));
short_names = cell(size(codes));
colors = zeros(numel(codes),3);
for k = 1:numel(codes)
    switch codes{k}
        case 'A'
            names{k} = sprintf('%s alone',C.electrode_A);
            short_names{k} = 'A';
            colors(k,:) = [0.10 0.45 0.80];
        case 'B'
            names{k} = sprintf('%s alone',C.electrode_B);
            short_names{k} = 'B';
            colors(k,:) = [0.85 0.33 0.10];
        case 'AB'
            names{k} = 'A+B simultaneous';
            short_names{k} = 'A+B';
            colors(k,:) = [0.15 0.65 0.30];
        case 'A_to_B'
            names{k} = 'A->B sequential';
            short_names{k} = 'A->B';
            colors(k,:) = [0.60 0.25 0.70];
        case 'B_to_A'
            names{k} = 'B->A sequential';
            short_names{k} = 'B->A';
            colors(k,:) = [0.15 0.65 0.70];
        otherwise
            error('StageC2:UnknownCondition','Unknown condition code %s.',codes{k});
    end
end
end

function value = logical_text(tf)
if tf
    value = 'true';
else
    value = 'false';
end
end
