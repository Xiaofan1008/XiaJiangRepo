%% ============================================================
%   PLOT: LINEARITY CURVE (Actual Paired Response vs Predicted Sum)
%   Reads ONLY the exported .mat file (no access to the lab's raw data,
%   internal helper functions, or shared drive needed).
%
%   For each electrode pair (Set) and amplitude:
%     R_A, R_B      = single-electrode-alone responses (from the Single
%                     data in this same file)
%     R_predicted   = R_A + R_B                (constant across ISI)
%     R_actual(ISI) = the measured paired response at each ISI (from the
%                     Paired data in this same file)
%     Ratio(ISI)    = R_actual(ISI) / R_predicted
%                     ~1 = linear summation; <1 = sub-additive;
%                     >1 = supra-additive
%
%   Just set mat_path below and run. Every electrode pair (Set) found in
%   the file gets its own pair of figures automatically. Leave
%   FilterSets/FilterAmps empty to plot everything found; narrow them
%   down to specific values if you only want a subset.
% ============================================================
clear;

%% ================= USER SETTINGS ============================
mat_path = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Luke_Data/DX023_10uA_v2/DX023_10uA_Export.mat';

FilterSets = [];   % Set number(s), e.g. [2]; [] = all Sets in the file
FilterAmps = [];   % Amplitude(s) in uA; [] = all

%% ================= LOAD DATA =====================
Loaded = load(mat_path, 'Paired', 'Single', 'ExportMetadata');
Paired = Loaded.Paired;
Single = Loaded.Single;
Meta   = Loaded.ExportMetadata;

fprintf('Loaded export for AnimalID = %s\n', Meta.AnimalID);
fprintf('Found %d electrode pair(s) (Sets) in this file.\n', numel(Paired));

%% ================= PLOT: ONE PAIR OF FIGURES PER SET =====================
amp_colors_base = lines(8); % reused/cycled per Set as needed

for si = 1:numel(Paired)

    setEntry = Paired(si);
    if ~isempty(FilterSets) && ~ismember(setEntry.SetNumber, FilterSets)
        continue;
    end

    % Find this Set's matching entry in Single (same SetNumber)
    singleIdx = find(arrayfun(@(s) s.SetNumber == setEntry.SetNumber, Single), 1);
    if isempty(singleIdx)
        fprintf('Set %d: no matching Single-electrode data found, skipping.\n', setEntry.SetNumber);
        continue;
    end
    singleEntry = Single(singleIdx);

    ampsAvail = [setEntry.Amp.Amplitude_uA];
    if isempty(FilterAmps)
        ampsToPlot = ampsAvail;
    else
        ampsToPlot = ampsAvail(ismember(round(ampsAvail,4), round(FilterAmps,4)));
    end
    ampsToPlot = sort(unique(ampsToPlot));

    if isempty(ampsToPlot)
        fprintf('Set %d: no matching amplitudes, skipping.\n', setEntry.SetNumber);
        continue;
    end

    colors = amp_colors_base(mod((1:numel(ampsToPlot))-1, size(amp_colors_base,1)) + 1, :);

    % ---- Compute R_predicted, R_actual, Ratio for each amplitude ----
    isisAll = [];
    R_actual_byAmp = cell(numel(ampsToPlot),1);
    R_predicted_byAmp = nan(numel(ampsToPlot),1);

    for a = 1:numel(ampsToPlot)
        target_amp = ampsToPlot(a);
        ai = find(abs(ampsAvail - target_amp) < 1e-4, 1);
        ampEntry = setEntry.Amp(ai);

        isisAvail = [ampEntry.ISI.ISI_ms];
        isisAll = union(isisAll, isisAvail);

        y = nan(numel(isisAvail),1);
        for p = 1:numel(isisAvail)
            allVals = vertcat(ampEntry.ISI(p).Channel.SpikeCount);
            if ~isempty(allVals), y(p) = mean(allVals); end
        end
        R_actual_byAmp{a} = [isisAvail(:), y(:)];

        % Matching Single-electrode amplitude entry
        ai_se = find(abs([singleEntry.Amp.Amplitude_uA] - target_amp) < 1e-4, 1);
        if isempty(ai_se)
            fprintf('Set %d, %.1f uA: no matching Single-electrode amplitude, skipping predicted line.\n', ...
                setEntry.SetNumber, target_amp);
            continue;
        end
        singleAmpEntry = singleEntry.Amp(ai_se);

        if numel(singleAmpEntry.Electrode) < 2
            fprintf('Set %d, %.1f uA: Single-electrode data missing one electrode, skipping predicted line.\n', ...
                setEntry.SetNumber, target_amp);
            continue;
        end

        R_A = mean(vertcat(singleAmpEntry.Electrode(1).Channel.SpikeCount));
        R_B = mean(vertcat(singleAmpEntry.Electrode(2).Channel.SpikeCount));
        R_predicted_byAmp(a) = R_A + R_B;
    end

    isisAll = sort(isisAll);

    %% ---- Plot A: Ratio vs ISI ----
    figure('Color','w','Position',[100 100 800 600]); hold on;
    yline(1.0, 'k--', 'LineWidth', 1);

    for a = 1:numel(ampsToPlot)
        if isnan(R_predicted_byAmp(a)), continue; end
        vals = R_actual_byAmp{a};
        ratio = vals(:,2) / R_predicted_byAmp(a);
        valid = ~isnan(ratio);
        if ~any(valid), continue; end

        plot(vals(valid,1), ratio(valid), '-o', 'Color', colors(a,:), ...
            'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerSize', 8, ...
            'DisplayName', sprintf('%.1f uA', ampsToPlot(a)));
    end

    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold');
    ylabel('Ratio (Actual / Predicted Linear Sum)', 'FontWeight','bold');
    ylim([0.5 1.4]);
    title(sprintf('Set %d -- Linearity Ratio', ...
        setEntry.SetNumber), 'FontWeight','bold');
    box off; legend('Location','best','Box','off');

    %% ---- Plot B: Actual vs ISI, with predicted lines overlaid ----
    figure('Color','w','Position',[950 100 800 600]); hold on;

    for a = 1:numel(ampsToPlot)
        vals = R_actual_byAmp{a};
        valid = ~isnan(vals(:,2));
        if any(valid)
            plot(vals(valid,1), vals(valid,2), '-o', 'Color', colors(a,:), ...
                'LineWidth', 2, 'MarkerFaceColor','w', 'MarkerSize', 8, ...
                'DisplayName', sprintf('%.1f uA', ampsToPlot(a)));
        end
        if ~isnan(R_predicted_byAmp(a))
            yline(R_predicted_byAmp(a), '--', 'Color', colors(a,:), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%.1f uA predicted (A+B)', ampsToPlot(a)));
        end
    end

    xlabel('Inter-Stimulus Interval (ms)', 'FontWeight','bold');
    ylabel('Baseline corrected spike count / Trial', 'FontWeight','bold');
    ylim([0.5 1.4]);
    title(sprintf('Set %d -- Actual vs Predicted', setEntry.SetNumber), 'FontWeight','bold');
    box off; legend('Location','best','Box','off');
end