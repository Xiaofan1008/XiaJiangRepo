%% CHECK MULTI-ISI LINEAR-PREDICTION BIAS
%
% PURPOSE
%   Diagnose why the shifted single-pulse linear prediction may be much
%   higher than the observed simultaneous/sequential responses.
%
% THREE CALCULATION METHODS ARE COMPARED
%
%   1. Signed
%      Baseline-correct every trial and preserve negative values.
%
%   2. Trial-clipped
%      Force every baseline-corrected trial to be >= 0 before averaging.
%      This is the method suspected of inflating the prediction.
%
%   3. Mean-clipped
%      Average signed trials first, then force the final channel mean to
%      be >= 0. For the prediction, A and shifted B are added before the
%      final clipping operation.
%
% OUTPUT
%   Figure 1: prediction-method comparison
%   Figure 2: separate prediction components
%   Figure 3: response-window sensitivity
%   Figure 4: bias caused by trial-level clipping
%   Numerical summaries in the MATLAB Command Window
%
% This script DOES NOT save or modify any files.

clear;
close all;

%% ======================== USER SETTINGS ===============================

% Exported Multi-ISI ModelData file
model_data_file = ...
    '/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/Luck_DataCheck_Codes/Multi_ISI_Output/DX020_XiaISISimSeq3_MultiISI_ModelData.mat';

% Empty means use all sequential amplitudes
Amplitudes_To_Check = [];

% Empty means use all positive PTDs
PTDs_To_Check = [];

% ModelData.channels.ChannelIndex values.
% Empty means use all exported responding channels.
Channels_To_Check = [];

% Primary windows used for Figures 1, 2, and 4
Response_Window_ms = [2 45];
Baseline_Window_ms = [-50 -5];

% Windows compared in Figure 3
Response_Windows_To_Compare = {
    [2 20]
    [2 30]
    [2 45]
};

% Use manually corrected response populations for each sequence order
Use_Order_Specific_Response_Channels = true;

% Choose diagnostic figures
Figures_To_Plot = [1 2 3 4];

figure_position = [80 80 1250 520];
line_width = 2;
marker_size = 7;

%% ======================= LOAD MODELDATA ================================

if ~isfile(model_data_file)
    error('Diagnostic:FileNotFound', ...
        'ModelData file not found:\n%s', model_data_file);
end

Loaded = load(model_data_file, 'ModelData');

if ~isfield(Loaded, 'ModelData')
    error('Diagnostic:MissingModelData', ...
        'The selected file does not contain ModelData.');
end

ModelData = Loaded.ModelData;

if ~isfield(ModelData, 'conditions') || ...
        ~isfield(ModelData, 'channels')
    error('Diagnostic:InvalidModelData', ...
        'ModelData.conditions or ModelData.channels is missing.');
end

Conditions = ModelData.conditions;

%% ======================== CHANNEL SELECTION ===========================

all_channel_ids = double(ModelData.channels.ChannelIndex(:).');

if isempty(Channels_To_Check)
    selected_channel_ids = all_channel_ids;
else
    requested = unique(double(Channels_To_Check(:).'), 'stable');
    selected_channel_ids = requested(ismember(requested, all_channel_ids));
end

if isempty(selected_channel_ids)
    error('Diagnostic:NoChannels', ...
        'None of Channels_To_Check is available.');
end

[~, channel_rows] = ismember(selected_channel_ids, all_channel_ids);
nChannels = numel(channel_rows);

%% ================= ORDER-SPECIFIC RESPONSE MASKS ======================

mask_A_to_B = true(nChannels, 1);
mask_B_to_A = true(nChannels, 1);

if Use_Order_Specific_Response_Channels

    variables = ModelData.channels.Properties.VariableNames;

    if ismember('Responsive_A_to_B', variables)
        mask_A_to_B = logical( ...
            ModelData.channels.Responsive_A_to_B(channel_rows));
    else
        warning(['Responsive_A_to_B is missing. All exported channels ' ...
                 'will be used for A->B.']);
    end

    if ismember('Responsive_B_to_A', variables)
        mask_B_to_A = logical( ...
            ModelData.channels.Responsive_B_to_A(channel_rows));
    else
        warning(['Responsive_B_to_A is missing. All exported channels ' ...
                 'will be used for B->A.']);
    end
end

%% ====================== AVAILABLE CONDITIONS =========================

available_ptds = [];
available_amps = [];

for ic = 1:numel(Conditions)

    if strcmp(Conditions(ic).stimulation_type, 'sequential')

        if isfinite(Conditions(ic).PTD_ms) && Conditions(ic).PTD_ms > 0
            available_ptds(end+1) = Conditions(ic).PTD_ms; %#ok<SAGROW>
        end

        if isfield(Conditions(ic), 'amplitude')
            available_amps = [available_amps, ...
                [Conditions(ic).amplitude.amplitude_uA]]; %#ok<AGROW>
        end
    end
end

available_ptds = sort(unique(available_ptds));
available_amps = sort(unique(available_amps));

if isempty(PTDs_To_Check)
    selected_ptds = available_ptds;
else
    selected_ptds = keep_available_values( ...
        PTDs_To_Check, available_ptds);
end

if isempty(Amplitudes_To_Check)
    selected_amps = available_amps;
else
    selected_amps = keep_available_values( ...
        Amplitudes_To_Check, available_amps);
end

if isempty(selected_ptds)
    error('Diagnostic:NoPTDs', 'No selected PTDs are available.');
end

if isempty(selected_amps)
    error('Diagnostic:NoAmplitudes', ...
        'No selected amplitudes are available.');
end

prediction_ptds = [0, selected_ptds];

%% ======================= WINDOW VALIDATION ============================

stored_window = double( ...
    ModelData.metadata.stored_spike_window_ms(:).');

all_test_windows = [{Response_Window_ms}; ...
                    Response_Windows_To_Compare(:)];

for iw = 1:numel(all_test_windows)

    this_window = all_test_windows{iw};

    if this_window(1) < stored_window(1) || ...
            this_window(2) > stored_window(2)
        error('Diagnostic:WindowOutsideData', ...
            'Response window [%g,%g) is outside stored window [%g,%g].', ...
            this_window, stored_window);
    end

    shifted_window = this_window - max(selected_ptds);

    if shifted_window(1) < stored_window(1) || ...
            shifted_window(2) > stored_window(2)
        error('Diagnostic:ShiftedWindowOutsideData', ...
            ['Response window [%g,%g) at PTD %g ms requires stored ' ...
             'single-pulse spikes in [%g,%g) ms.'], ...
            this_window, max(selected_ptds), shifted_window);
    end
end

%% ========================= GENERAL INFORMATION ========================

[electrode_A, electrode_B] = stimulation_names(ModelData);

fprintf('\n');
fprintf('============================================================\n');
fprintf('MULTI-ISI LINEAR-PREDICTION DIAGNOSTIC\n');
fprintf('============================================================\n');
fprintf('Pair: %s + %s\n', electrode_A, electrode_B);
fprintf('Channels selected: %d\n', nChannels);
fprintf('A->B responsive channels: %d\n', sum(mask_A_to_B));
fprintf('B->A responsive channels: %d\n', sum(mask_B_to_A));
fprintf('Primary response window: [%g,%g) ms\n', Response_Window_ms);
fprintf('Baseline window: [%g,%g) ms\n', Baseline_Window_ms);
fprintf('PTDs: %s ms\n', num2str(selected_ptds));
fprintf('Amplitudes: %s uA\n', num2str(selected_amps));
fprintf('No files will be modified or saved.\n');

%% ======================== MAIN DIAGNOSTIC =============================

for amp = selected_amps

    fprintf('\n------------------------------------------------------------\n');
    fprintf('AMPLITUDE: %g uA\n', amp);
    fprintf('------------------------------------------------------------\n');

    % Simultaneous A+B response at PTD 0
    AB = calculate_condition_methods( ...
        Conditions, 'AB', 0, amp, channel_rows, ...
        Response_Window_ms, Baseline_Window_ms);

    AB_AtoB = apply_mask_to_methods(AB, mask_A_to_B);
    AB_BtoA = apply_mask_to_methods(AB, mask_B_to_A);

    % Sequential observed responses
    Seq_AtoB = initialise_method_matrix(nChannels, numel(selected_ptds));
    Seq_BtoA = initialise_method_matrix(nChannels, numel(selected_ptds));

    for ip = 1:numel(selected_ptds)

        ptd = selected_ptds(ip);

        this_AtoB = calculate_condition_methods( ...
            Conditions, 'A_to_B', ptd, amp, channel_rows, ...
            Response_Window_ms, Baseline_Window_ms);

        this_BtoA = calculate_condition_methods( ...
            Conditions, 'B_to_A', ptd, amp, channel_rows, ...
            Response_Window_ms, Baseline_Window_ms);

        Seq_AtoB = insert_method_column(Seq_AtoB, this_AtoB, ip);
        Seq_BtoA = insert_method_column(Seq_BtoA, this_BtoA, ip);
    end

    Seq_AtoB = apply_mask_to_methods(Seq_AtoB, mask_A_to_B);
    Seq_BtoA = apply_mask_to_methods(Seq_BtoA, mask_B_to_A);

    % Shifted single-pulse prediction
    Pred_AtoB = calculate_prediction_methods( ...
        Conditions, 'A', 'B', amp, prediction_ptds, ...
        channel_rows, Response_Window_ms, Baseline_Window_ms);

    Pred_BtoA = calculate_prediction_methods( ...
        Conditions, 'B', 'A', amp, prediction_ptds, ...
        channel_rows, Response_Window_ms, Baseline_Window_ms);

    Pred_AtoB = apply_mask_to_methods(Pred_AtoB, mask_A_to_B);
    Pred_BtoA = apply_mask_to_methods(Pred_BtoA, mask_B_to_A);

    %% ---------------- FIGURE 1: METHOD COMPARISON ---------------------

    if ismember(1, Figures_To_Plot)

        figure('Color', 'w', ...
            'Name', sprintf('Method comparison | %g uA', amp), ...
            'NumberTitle', 'off', ...
            'Position', figure_position);

        tiledlayout(1, 2, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        ax = nexttile;
        plot_method_comparison(ax, selected_ptds, prediction_ptds, ...
            Seq_AtoB, AB_AtoB, Pred_AtoB, ...
            sprintf('%s -> %s', electrode_A, electrode_B), ...
            line_width, marker_size);

        ax = nexttile;
        plot_method_comparison(ax, selected_ptds, prediction_ptds, ...
            Seq_BtoA, AB_BtoA, Pred_BtoA, ...
            sprintf('%s -> %s', electrode_B, electrode_A), ...
            line_width, marker_size);

        sgtitle(sprintf('%s + %s | %g uA | clipping-method comparison', ...
            electrode_A, electrode_B, amp), ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end

    %% ---------------- FIGURE 2: PREDICTION COMPONENTS ----------------

    if ismember(2, Figures_To_Plot)

        figure('Color', 'w', ...
            'Name', sprintf('Prediction components | %g uA', amp), ...
            'NumberTitle', 'off', ...
            'Position', figure_position);

        tiledlayout(1, 2, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        ax = nexttile;
        plot_prediction_components(ax, prediction_ptds, Pred_AtoB, ...
            electrode_A, electrode_B, line_width, marker_size);

        ax = nexttile;
        plot_prediction_components(ax, prediction_ptds, Pred_BtoA, ...
            electrode_B, electrode_A, line_width, marker_size);

        sgtitle(sprintf('%s + %s | %g uA | signed prediction components', ...
            electrode_A, electrode_B, amp), ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end

    %% ---------------- FIGURE 3: WINDOW SENSITIVITY -------------------

    if ismember(3, Figures_To_Plot)

        figure('Color', 'w', ...
            'Name', sprintf('Window sensitivity | %g uA', amp), ...
            'NumberTitle', 'off', ...
            'Position', figure_position);

        tiledlayout(1, 2, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        ax_AtoB = nexttile;
        hold(ax_AtoB, 'on');

        ax_BtoA = nexttile;
        hold(ax_BtoA, 'on');

        colours = lines(numel(Response_Windows_To_Compare));

        for iw = 1:numel(Response_Windows_To_Compare)

            response_window = Response_Windows_To_Compare{iw};

            [W_AB_AtoB, W_Seq_AtoB, W_Pred_AtoB] = ...
                calculate_window_result(Conditions, amp, ...
                'A_to_B', 'A', 'B', selected_ptds, ...
                channel_rows, mask_A_to_B, ...
                response_window, Baseline_Window_ms);

            [W_AB_BtoA, W_Seq_BtoA, W_Pred_BtoA] = ...
                calculate_window_result(Conditions, amp, ...
                'B_to_A', 'B', 'A', selected_ptds, ...
                channel_rows, mask_B_to_A, ...
                response_window, Baseline_Window_ms);

            window_label = sprintf('[%g,%g) ms', response_window);

            plot_window_result(ax_AtoB, selected_ptds, prediction_ptds, ...
                W_AB_AtoB, W_Seq_AtoB, W_Pred_AtoB, ...
                colours(iw,:), window_label, line_width, marker_size);

            plot_window_result(ax_BtoA, selected_ptds, prediction_ptds, ...
                W_AB_BtoA, W_Seq_BtoA, W_Pred_BtoA, ...
                colours(iw,:), window_label, line_width, marker_size);
        end

        format_diagnostic_axis(ax_AtoB, ...
            sprintf('%s -> %s', electrode_A, electrode_B));

        format_diagnostic_axis(ax_BtoA, ...
            sprintf('%s -> %s', electrode_B, electrode_A));

        sgtitle(sprintf('%s + %s | %g uA | response-window sensitivity', ...
            electrode_A, electrode_B, amp), ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end

    %% ---------------- FIGURE 4: CLIPPING BIAS ------------------------

    if ismember(4, Figures_To_Plot)

        figure('Color', 'w', ...
            'Name', sprintf('Clipping bias | %g uA', amp), ...
            'NumberTitle', 'off', ...
            'Position', figure_position);

        tiledlayout(1, 2, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        ax = nexttile;
        plot_clipping_bias(ax, selected_ptds, prediction_ptds, ...
            Seq_AtoB, AB_AtoB, Pred_AtoB, ...
            sprintf('%s -> %s', electrode_A, electrode_B), ...
            line_width, marker_size);

        ax = nexttile;
        plot_clipping_bias(ax, selected_ptds, prediction_ptds, ...
            Seq_BtoA, AB_BtoA, Pred_BtoA, ...
            sprintf('%s -> %s', electrode_B, electrode_A), ...
            line_width, marker_size);

        sgtitle(sprintf('%s + %s | %g uA | trial-clipping inflation', ...
            electrode_A, electrode_B, amp), ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end

    %% ---------------- COMMAND-WINDOW SUMMARY -------------------------

    print_order_summary( ...
        sprintf('%s -> %s', electrode_A, electrode_B), ...
        selected_ptds, AB_AtoB, Seq_AtoB, Pred_AtoB);

    print_order_summary( ...
        sprintf('%s -> %s', electrode_B, electrode_A), ...
        selected_ptds, AB_BtoA, Seq_BtoA, Pred_BtoA);
end

fprintf('\n============================================================\n');
fprintf('DIAGNOSTIC COMPLETE\n');
fprintf('No data or analysis files were modified.\n');
fprintf('============================================================\n');

%% ========================== FUNCTIONS =================================

function Methods = initialise_method_matrix(nChannels, nConditions)

Methods.signed = nan(nChannels, nConditions);
Methods.trial_clipped = nan(nChannels, nConditions);
Methods.mean_clipped = nan(nChannels, nConditions);

end

function Output = insert_method_column(Output, Input, column)

Output.signed(:,column) = Input.signed;
Output.trial_clipped(:,column) = Input.trial_clipped;
Output.mean_clipped(:,column) = Input.mean_clipped;

end

function Methods = calculate_condition_methods(Conditions, code, ptd, amp, ...
        channel_rows, response_window, baseline_window)

Methods = initialise_method_matrix(numel(channel_rows), 1);

channels = find_condition_channels(Conditions, code, ptd, amp);

if isempty(channels)
    return;
end

for k = 1:numel(channel_rows)

    row = channel_rows(k);

    if row > numel(channels) || ...
            ~isfield(channels(row), 'spike_times_ms')
        continue;
    end

    trial_values = baseline_corrected_trial_counts( ...
        channels(row).spike_times_ms, ...
        response_window, baseline_window);

    Methods.signed(k) = mean_finite(trial_values);

    Methods.trial_clipped(k) = ...
        mean_finite(max(0, trial_values));

    if isfinite(Methods.signed(k))
        Methods.mean_clipped(k) = max(0, Methods.signed(k));
    end
end

end

function Prediction = calculate_prediction_methods(Conditions, ...
        first_code, second_code, amp, ptds, channel_rows, ...
        response_window, baseline_window)

nChannels = numel(channel_rows);
nPTDs = numel(ptds);

Prediction = initialise_method_matrix(nChannels, nPTDs);

Prediction.first_signed = nan(nChannels, nPTDs);
Prediction.second_signed = nan(nChannels, nPTDs);
Prediction.first_trial_clipped = nan(nChannels, nPTDs);
Prediction.second_trial_clipped = nan(nChannels, nPTDs);

first_channels = find_condition_channels( ...
    Conditions, first_code, 0, amp);

second_channels = find_condition_channels( ...
    Conditions, second_code, 0, amp);

if isempty(first_channels) || isempty(second_channels)
    return;
end

for k = 1:nChannels

    row = channel_rows(k);

    if row > numel(first_channels) || ...
            row > numel(second_channels)
        continue;
    end

    first_trials = baseline_corrected_trial_counts( ...
        first_channels(row).spike_times_ms, ...
        response_window, baseline_window);

    first_signed = mean_finite(first_trials);
    first_trial_clipped = mean_finite(max(0, first_trials));

    for ip = 1:nPTDs

        shifted_second_window = response_window - ptds(ip);

        second_trials = baseline_corrected_trial_counts( ...
            second_channels(row).spike_times_ms, ...
            shifted_second_window, baseline_window);

        second_signed = mean_finite(second_trials);
        second_trial_clipped = mean_finite(max(0, second_trials));

        Prediction.first_signed(k,ip) = first_signed;
        Prediction.second_signed(k,ip) = second_signed;

        Prediction.first_trial_clipped(k,ip) = ...
            first_trial_clipped;

        Prediction.second_trial_clipped(k,ip) = ...
            second_trial_clipped;

        % Method 1: retain signed values
        Prediction.signed(k,ip) = ...
            first_signed + second_signed;

        % Method 2: clip both single responses before addition
        Prediction.trial_clipped(k,ip) = ...
            first_trial_clipped + second_trial_clipped;

        % Method 3: add signed means, then clip once
        if isfinite(Prediction.signed(k,ip))
            Prediction.mean_clipped(k,ip) = ...
                max(0, Prediction.signed(k,ip));
        end
    end
end

end

function trial_values = baseline_corrected_trial_counts( ...
        spike_trials, response_window, baseline_window)

nTrials = numel(spike_trials);
trial_values = nan(nTrials, 1);

response_duration = diff(response_window);
baseline_duration = diff(baseline_window);

for it = 1:nTrials

    spikes = double(spike_trials{it}(:));

    response_count = sum( ...
        spikes >= response_window(1) & ...
        spikes < response_window(2));

    baseline_count = sum( ...
        spikes >= baseline_window(1) & ...
        spikes < baseline_window(2));

    expected_baseline_count = baseline_count * ...
        (response_duration / baseline_duration);

    trial_values(it) = response_count - expected_baseline_count;
end

end

function channels = find_condition_channels( ...
        Conditions, code, ptd, amp)

channels = [];

for ic = 1:numel(Conditions)

    if ~strcmp(Conditions(ic).code, code)
        continue;
    end

    if abs(Conditions(ic).PTD_ms - ptd) >= 1e-6
        continue;
    end

    for ia = 1:numel(Conditions(ic).amplitude)

        if abs(Conditions(ic).amplitude(ia).amplitude_uA - amp) < 1e-6
            channels = Conditions(ic).amplitude(ia).channel;
            return;
        end
    end
end

end

function Output = apply_mask_to_methods(Input, mask)

Output = Input;
field_names = fieldnames(Output);

for k = 1:numel(field_names)

    field_name = field_names{k};
    values = Output.(field_name);

    if isnumeric(values) && size(values,1) == numel(mask)
        values(~mask,:) = NaN;
        Output.(field_name) = values;
    end
end

end

function [AB, Seq, Pred] = calculate_window_result(Conditions, amp, ...
        sequence_code, first_code, second_code, ptds, ...
        channel_rows, response_mask, response_window, baseline_window)

AB = calculate_condition_methods( ...
    Conditions, 'AB', 0, amp, channel_rows, ...
    response_window, baseline_window);

Seq = initialise_method_matrix(numel(channel_rows), numel(ptds));

for ip = 1:numel(ptds)

    this_result = calculate_condition_methods( ...
        Conditions, sequence_code, ptds(ip), amp, channel_rows, ...
        response_window, baseline_window);

    Seq = insert_method_column(Seq, this_result, ip);
end

Pred = calculate_prediction_methods( ...
    Conditions, first_code, second_code, amp, [0, ptds], ...
    channel_rows, response_window, baseline_window);

AB = apply_mask_to_methods(AB, response_mask);
Seq = apply_mask_to_methods(Seq, response_mask);
Pred = apply_mask_to_methods(Pred, response_mask);

end

function plot_method_comparison(ax, ptds, prediction_ptds, ...
        Seq, AB, Pred, order_label, lw, ms)

hold(ax, 'on');

method_names = {'signed', 'trial_clipped', 'mean_clipped'};
labels = {'Signed', 'Trial-clipped', 'Mean-clipped'};
colours = [
    0.35 0.35 0.35
    0.85 0.25 0.12
    0.10 0.40 0.80
];

for im = 1:numel(method_names)

    field_name = method_names{im};

    observed_y = [ ...
        population_mean(AB.(field_name)), ...
        population_curve(Seq.(field_name))];

    observed_x = [0, ptds];

    prediction_y = population_curve(Pred.(field_name));

    valid_observed = isfinite(observed_y);
    valid_prediction = isfinite(prediction_y);

    plot(ax, observed_x(valid_observed), ...
        observed_y(valid_observed), '-o', ...
        'Color', colours(im,:), ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', ms, ...
        'LineWidth', lw, ...
        'DisplayName', sprintf('Observed: %s', labels{im}));

    plot(ax, prediction_ptds(valid_prediction), ...
        prediction_y(valid_prediction), '--', ...
        'Color', colours(im,:), ...
        'LineWidth', lw, ...
        'DisplayName', sprintf('Prediction: %s', labels{im}));
end

format_diagnostic_axis(ax, order_label);

end

function plot_prediction_components(ax, ptds, Prediction, ...
        first_name, second_name, lw, ms)

hold(ax, 'on');

first_component = population_curve(Prediction.first_signed);
second_component = population_curve(Prediction.second_signed);
linear_sum = population_curve(Prediction.signed);

valid = isfinite(first_component);
plot(ax, ptds(valid), first_component(valid), '-o', ...
    'Color', [0.20 0.55 0.25], ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', ms, ...
    'LineWidth', lw, ...
    'DisplayName', sprintf('%s contribution', first_name));

valid = isfinite(second_component);
plot(ax, ptds(valid), second_component(valid), '-o', ...
    'Color', [0.55 0.25 0.70], ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', ms, ...
    'LineWidth', lw, ...
    'DisplayName', sprintf('Shifted %s contribution', second_name));

valid = isfinite(linear_sum);
plot(ax, ptds(valid), linear_sum(valid), '--', ...
    'Color', [0.10 0.10 0.10], ...
    'LineWidth', lw, ...
    'DisplayName', 'Signed linear sum');

format_diagnostic_axis(ax, ...
    sprintf('%s -> %s', first_name, second_name));

end

function plot_window_result(ax, ptds, prediction_ptds, ...
        AB, Seq, Pred, colour, window_label, lw, ms)

observed_y = [ ...
    population_mean(AB.mean_clipped), ...
    population_curve(Seq.mean_clipped)];

observed_x = [0, ptds];

prediction_y = population_curve(Pred.mean_clipped);

valid = isfinite(observed_y);
plot(ax, observed_x(valid), observed_y(valid), '-o', ...
    'Color', colour, ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', ms, ...
    'LineWidth', lw, ...
    'DisplayName', sprintf('Observed %s', window_label));

valid = isfinite(prediction_y);
plot(ax, prediction_ptds(valid), prediction_y(valid), '--', ...
    'Color', colour, ...
    'LineWidth', lw, ...
    'DisplayName', sprintf('Prediction %s', window_label));

end

function plot_clipping_bias(ax, ptds, prediction_ptds, ...
        Seq, AB, Pred, order_label, lw, ms)

hold(ax, 'on');

observed_trial = [ ...
    population_mean(AB.trial_clipped), ...
    population_curve(Seq.trial_clipped)];

observed_mean = [ ...
    population_mean(AB.mean_clipped), ...
    population_curve(Seq.mean_clipped)];

prediction_trial = population_curve(Pred.trial_clipped);
prediction_mean = population_curve(Pred.mean_clipped);

observed_bias = observed_trial - observed_mean;
prediction_bias = prediction_trial - prediction_mean;

observed_x = [0, ptds];

valid = isfinite(observed_bias);
plot(ax, observed_x(valid), observed_bias(valid), '-o', ...
    'Color', [0.10 0.45 0.80], ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', ms, ...
    'LineWidth', lw, ...
    'DisplayName', 'Observed clipping bias');

valid = isfinite(prediction_bias);
plot(ax, prediction_ptds(valid), prediction_bias(valid), '--o', ...
    'Color', [0.85 0.25 0.12], ...
    'MarkerFaceColor', 'w', ...
    'MarkerSize', ms, ...
    'LineWidth', lw, ...
    'DisplayName', 'Prediction clipping bias');

yline(ax, 0, ':k', 'HandleVisibility', 'off');

ylabel(ax, ...
    'Trial-clipped minus mean-clipped spike count');

format_diagnostic_axis(ax, order_label);

end

function print_order_summary(order_label, ptds, AB, Seq, Pred)

last_index = numel(ptds);
longest_ptd = ptds(last_index);
prediction_index = last_index + 1;

fprintf('\n%s\n', order_label);
fprintf('Longest PTD checked: %g ms\n', longest_ptd);

fprintf('  Simultaneous A+B, mean-clipped: %.4f\n', ...
    population_mean(AB.mean_clipped));

fprintf('  Sequential signed:             %.4f\n', ...
    population_mean(Seq.signed(:,last_index)));

fprintf('  Sequential trial-clipped:      %.4f\n', ...
    population_mean(Seq.trial_clipped(:,last_index)));

fprintf('  Sequential mean-clipped:       %.4f\n', ...
    population_mean(Seq.mean_clipped(:,last_index)));

fprintf('  First single contribution:     %.4f\n', ...
    population_mean(Pred.first_signed(:,prediction_index)));

fprintf('  Shifted second contribution:   %.4f\n', ...
    population_mean(Pred.second_signed(:,prediction_index)));

fprintf('  Prediction signed:             %.4f\n', ...
    population_mean(Pred.signed(:,prediction_index)));

fprintf('  Prediction trial-clipped:      %.4f\n', ...
    population_mean(Pred.trial_clipped(:,prediction_index)));

fprintf('  Prediction mean-clipped:       %.4f\n', ...
    population_mean(Pred.mean_clipped(:,prediction_index)));

prediction_bias = ...
    population_mean(Pred.trial_clipped(:,prediction_index)) - ...
    population_mean(Pred.mean_clipped(:,prediction_index));

observed_bias = ...
    population_mean(Seq.trial_clipped(:,last_index)) - ...
    population_mean(Seq.mean_clipped(:,last_index));

fprintf('  Trial-clipping bias, observed: %.4f\n', observed_bias);
fprintf('  Trial-clipping bias, predicted: %.4f\n', prediction_bias);

end

function curve = population_curve(values)

curve = nan(1, size(values,2));

for k = 1:size(values,2)
    curve(k) = population_mean(values(:,k));
end

end

function value = population_mean(values)

values = double(values(:));
values = values(isfinite(values));

if isempty(values)
    value = NaN;
else
    value = mean(values);
end

end

function value = mean_finite(values)

values = double(values(:));
values = values(isfinite(values));

if isempty(values)
    value = NaN;
else
    value = mean(values);
end

end

function values = keep_available_values(requested, available)

requested = unique(double(requested(:).'), 'stable');
available = double(available(:).');
values = [];

for value = requested

    index = find(abs(available - value) < 1e-6, 1);

    if ~isempty(index)
        values(end+1) = available(index); %#ok<AGROW>
    end
end

end

function [A, B] = stimulation_names(ModelData)

if isfield(ModelData, 'stimulation') && ...
        isfield(ModelData.stimulation, 'electrode_A') && ...
        isfield(ModelData.stimulation, 'electrode_B')

    A = char(string(ModelData.stimulation.electrode_A));
    B = char(string(ModelData.stimulation.electrode_B));
else
    A = 'A';
    B = 'B';
end

end

function format_diagnostic_axis(ax, plot_title)

xlabel(ax, 'ISI / PTD (ms)');
ylabel(ax, 'Baseline-corrected spike count / trial');
title(ax, plot_title, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'none');

grid(ax, 'on');
box(ax, 'off');
legend(ax, 'Location', 'best', 'Box', 'off');

end