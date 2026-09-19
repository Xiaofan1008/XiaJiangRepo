%% ============================================================
%   PLOT: ISI SPIKE-COUNT TUNING CURVE
%   Reads ONLY the exported .mat file (no access to the lab's raw data,
%   internal helper functions, or shared drive needed).
%
%   Just set mat_path below and run. Every electrode pair (Set) found in
%   the file gets its own figure automatically -- no need to know or
%   enter electrode numbers. Leave FilterAmps/FilterISIs empty to plot
%   everything found for that Set; narrow them down to specific values
%   if you only want a subset.
% ============================================================
clear;

%% ================= USER SETTINGS ============================
mat_path = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Luke_Data/DX023_10uA_v2/DX023_10uA_Export.mat';

FilterAmps = [];   % e.g. [5 10] to restrict to those amplitudes; [] = all
FilterISIs = [];   % e.g. [0 5 10] to restrict to those ISIs; [] = all

%% ================= LOAD DATA =====================
Loaded = load(mat_path, 'Paired', 'ExportMetadata');
Paired = Loaded.Paired;
Meta   = Loaded.ExportMetadata;

fprintf('Loaded export for AnimalID = %s\n', Meta.AnimalID);
fprintf('Spike-count window: [%g, %g] ms; baseline: [%g, %g] ms\n', ...
    Meta.FixedMacroWindow_ms, Meta.BaselineWindow_ms);
fprintf('Found %d electrode pair(s) (Sets) in this file.\n', numel(Paired));

%% ================= PLOT: ONE FIGURE PER SET =====================
for si = 1:numel(Paired)

    setEntry = Paired(si);
    elecs = setEntry.Electrodes;

    ampsAvail = [setEntry.Amp.Amplitude_uA];
    if isempty(FilterAmps)
        ampsToPlot = ampsAvail;
    else
        ampsToPlot = ampsAvail(ismember(round(ampsAvail,4), round(FilterAmps,4)));
    end
    ampsToPlot = sort(unique(ampsToPlot));

    if isempty(ampsToPlot)
        fprintf('Set %d (Electrodes %s): no matching amplitudes, skipping.\n', ...
            setEntry.SetNumber, num2str(elecs));
        continue;
    end

    colors = lines(numel(ampsToPlot));
    figure('Color','w','Position',[100 100 800 600]); hold on;

    for a = 1:numel(ampsToPlot)
        target_amp = ampsToPlot(a);
        ai = find(abs(ampsAvail - target_amp) < 1e-4, 1);
        ampEntry = setEntry.Amp(ai);

        isisAvail = [ampEntry.ISI.ISI_ms];
        if isempty(FilterISIs)
            isisToPlot = isisAvail;
        else
            isisToPlot = isisAvail(ismember(round(isisAvail,4), round(FilterISIs,4)));
        end
        isisToPlot = sort(unique(isisToPlot));

        y_mean = nan(numel(isisToPlot),1);
        y_sem  = nan(numel(isisToPlot),1);

        for p = 1:numel(isisToPlot)
            pi = find(abs(isisAvail - isisToPlot(p)) < 1e-4, 1);
            isiEntry = ampEntry.ISI(pi);

            % Pool every trial from every channel at this (Set, Amp, ISI)
            % into one combined set of values for the mean/SEM.
            allVals = vertcat(isiEntry.Channel.SpikeCount);

            if ~isempty(allVals)
                y_mean(p) = mean(allVals);
                y_sem(p)  = std(allVals) / sqrt(numel(allVals));
            end
        end

        valid = ~isnan(y_mean);
        if ~any(valid), continue; end

        errorbar(isisToPlot(valid), y_mean(valid), y_sem(valid), '-o', ...
            'Color', colors(a,:), 'LineWidth', 2, 'MarkerFaceColor','w', ...
            'MarkerSize', 8, 'DisplayName', sprintf('%.1f uA', target_amp));
    end

    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold');
    ylabel('Mean Net Spike Count / Trial', 'FontWeight','bold');
    title(sprintf('Set %d (Electrodes %s) -- Spike Count vs ISI', ...
        setEntry.SetNumber, num2str(elecs)), 'FontWeight','bold');
    box off; legend('Location','best','Box','off');
end