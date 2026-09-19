%% ============================================================
%   PLOT: RASTER + PSTH, ALL CHANNELS PER CONDITION
%   Reads ONLY the exported .mat file (no access to the lab's raw data,
%   internal helper functions, or shared drive needed).
%
%   One figure per condition (Set x Amp x ISI for Paired, or
%   Set x Amp x Electrode for Single) -- every channel for that
%   condition is a tile inside that one figure, not a separate figure
%   per channel.
%
%   Just set mat_path below and run. Leave the Filter* settings empty to
%   plot every matching condition found in the file; narrow them down to
%   get just the ones you want. MaxFigures stops with a warning instead
%   of silently opening a huge number of windows.
% ============================================================
clear;

%% ================= USER SETTINGS ============================
mat_path = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Luke_Data/DX023_10uA_v1/DX023_10uA_Export.mat';

DatasetType = 'Paired';   % 'Paired' or 'Single'

FilterSets       = [];   % Set number(s), e.g. [2]; [] = all Sets in the file
FilterAmps       = [];   % Amplitude(s) in uA; [] = all
FilterISIs       = [0 5];   % ISI(s) in ms -- Paired only, ignored for Single; [] = all
FilterElectrodes = [];   % Electrode(s) -- Single only, ignored for Paired; [] = all

PSTH_bin_ms    = 1;
PSTH_smooth_ms = 5;

Figure_Window_Style = 'docked';   % 'docked' (tabs in one window) or 'normal'

MaxFigures = 12;   % safety cap -- narrow the filters above if you hit this

%% ================= LOAD DATA =====================
Loaded = load(mat_path, DatasetType, 'ExportMetadata');
Data = Loaded.(DatasetType);
Meta = Loaded.ExportMetadata;

fprintf('Loaded export for AnimalID = %s\n', Meta.AnimalID);
fprintf('Spike-time window: [%g, %g] ms\n', Meta.SpikeTimeWindow_ms);
fprintf('Found %d Set(s) of type "%s" in this file.\n', numel(Data), DatasetType);

%% ================= COLLECT MATCHING CONDITIONS =====================
% Each entry: {titleStr, xline2_ms (or []), ChannelStructArray}
condJobs = {};

for si = 1:numel(Data)
    setEntry = Data(si);

    if ~isempty(FilterSets) && ~ismember(setEntry.SetNumber, FilterSets)
        continue;
    end

    for ai = 1:numel(setEntry.Amp)
        ampEntry = setEntry.Amp(ai);

        if ~isempty(FilterAmps) && ~any(abs(FilterAmps - ampEntry.Amplitude_uA) < 1e-4)
            continue;
        end

        if strcmp(DatasetType, 'Paired')
            for pi = 1:numel(ampEntry.ISI)
                isiEntry = ampEntry.ISI(pi);

                if ~isempty(FilterISIs) && ~any(abs(FilterISIs - isiEntry.ISI_ms) < 1e-4)
                    continue;
                end
                if isempty(isiEntry.Channel), continue; end

                titleStr = sprintf('Set %d (Elec %s) | %.1f uA | ISI %.1f ms', ...
                    setEntry.SetNumber, num2str(setEntry.Electrodes), ...
                    ampEntry.Amplitude_uA, isiEntry.ISI_ms);

                secondPulseMark = isiEntry.ISI_ms; % 0 = no separate mark needed
                condJobs{end+1} = {titleStr, secondPulseMark, isiEntry.Channel}; %#ok<AGROW>
            end

        else % 'Single'
            for ei = 1:numel(ampEntry.Electrode)
                elecEntry = ampEntry.Electrode(ei);

                if ~isempty(FilterElectrodes) && ~ismember(elecEntry.Electrode, FilterElectrodes)
                    continue;
                end
                if isempty(elecEntry.Channel), continue; end

                titleStr = sprintf('Set %d | Elec %d alone | %.1f uA', ...
                    setEntry.SetNumber, elecEntry.Electrode, ampEntry.Amplitude_uA);

                condJobs{end+1} = {titleStr, 0, elecEntry.Channel}; %#ok<AGROW>
            end
        end
    end
end

fprintf('%d matching condition(s) found.\n', numel(condJobs));

if numel(condJobs) > MaxFigures
    error(['%d conditions match your filters, which is more than MaxFigures (%d). ' ...
        'Narrow FilterSets/FilterAmps/FilterISIs/FilterElectrodes, ' ...
        'or raise MaxFigures if you really want that many figures.'], ...
        numel(condJobs), MaxFigures);
end

if isempty(condJobs)
    warning('No data matched your filters -- nothing to plot.');
end

%% ================= PLOT ONE FIGURE PER CONDITION =====================
win = Meta.SpikeTimeWindow_ms;
edges = win(1):PSTH_bin_ms:win(2);
ctrs = edges(1:end-1) + diff(edges)/2;
bin_s = PSTH_bin_ms/1000;
kernel_size = 2*ceil(2*PSTH_smooth_ms)+1;
g = gausswin(kernel_size); g = g/sum(g);

for j = 1:numel(condJobs)
    titleStr = condJobs{j}{1};
    secondPulseMark = condJobs{j}{2};
    ChannelArray = condJobs{j}{3};
    nCh = numel(ChannelArray);

    if strcmpi(Figure_Window_Style, 'docked')
        fig = figure('Color','w', 'Name', titleStr, 'NumberTitle','off', 'WindowStyle','docked');
    else
        fig = figure('Color','w', 'Name', titleStr, 'NumberTitle','off', 'Position',[50 50 1600 900]);
    end

    layout = tiledlayout(fig, 'flow', 'TileSpacing','compact', 'Padding','compact');
    title(layout, titleStr, 'FontSize', 14, 'FontWeight','bold', 'Interpreter','none');

    for ci = 1:nCh
        chanEntry = ChannelArray(ci);
        spikeTimesPerTrial = chanEntry.SpikeTimes_ms;
        nTrials = numel(spikeTimesPerTrial);

        allSpikes = vertcat(spikeTimesPerTrial{:});
        counts = histcounts(allSpikes, edges);
        rate = counts / (max(nTrials,1) * bin_s);
        rate_smooth = conv(rate, g, 'same');
        yMaxPSTH = max(50, ceil(max(rate_smooth)*1.1/10)*10);

        ax = nexttile(layout); hold(ax,'on');

        % ---- Left axis: PSTH ----
        yyaxis(ax,'left');
        plot(ax, ctrs, rate_smooth, 'LineWidth', 1.4);
        ylim(ax, [0 yMaxPSTH]); ylabel(ax, 'Rate (sp/s)');

        % ---- Right axis: raster ----
        yyaxis(ax,'right');
        for ti = 1:nTrials
            tt = spikeTimesPerTrial{ti};
            if isempty(tt), continue; end
            plot(ax, tt, ti*ones(size(tt)), '.', 'Color','k', 'MarkerSize', 4);
        end
        ylim(ax, [0 max(nTrials,1)+1]); set(ax,'YTick',[]);

        xline(ax, 0, 'r--', 'LineWidth', 1);
        if secondPulseMark > 0
            xline(ax, secondPulseMark, 'k:', 'LineWidth', 1);
        end
        xlim(ax, win);

        title(ax, sprintf('Ch %d', chanEntry.Channel), 'FontSize', 11, 'FontWeight','bold');
        if ci > nCh - ceil(sqrt(nCh))
            xlabel(ax, 'Time (ms)');
        end
    end
end