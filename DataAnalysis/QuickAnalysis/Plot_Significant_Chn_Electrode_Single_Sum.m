%% ============================================================
% Plot Combined Significant Channels (Sum of Two Single-Stim Sets)
% ============================================================
%  - Merges significance from both single-stimulation sets
%  - Red = significant in either set (union)
%  - Gray = not significant in any set
%  - Uses per-set amplitudes (FR_summary(s).amps)
% ============================================================

clear
%% === User Inputs ===
data_folder = '/Volumes/MACData/Data/Data_Xia/DX010/Xia_Exp1_Single1';
amp_to_plot = 5; % µA
nShank = 1;
chPerShank = 32;

if ~isfolder(data_folder)
    error('The specified folder does not exist.');
end

cd(data_folder);
mat_files = dir(fullfile(data_folder, '*FR_SigCh_ByAmp.mat'));
assert(~isempty(mat_files), 'No FR_SigCh_ByAmp.mat file found.');
load(fullfile(data_folder, mat_files(1).name), ...
    'FR_summary', 'E_MAP', 'uniqueComb');

fprintf('Loaded file: %s\n', mat_files(1).name);

%% === Detect stimulation type from folder name ===
if contains(lower(data_folder), 'single')
    stim_type = 'Separate Stimulation';
elseif contains(lower(data_folder), 'sim')
    stim_type = 'Simultaneous Stimulation';
elseif contains(lower(data_folder), 'seq')
    stim_type = 'Sequential Stimulation';
else
    stim_type = 'Stimulation';
end

fprintf('Combining significant channels for %.2f µA (%s)...\n', amp_to_plot, stim_type);

%% === Detect number of total channels ===
nCh_total = numel(E_MAP) - 1;  % exclude header
fprintf('Detected %d total channels in dataset.\n', nCh_total);

%% === Detect number of stimulation sets ===
nSets = numel(FR_summary);
fprintf('Detected %d stimulation set(s).\n', nSets);

if nSets < 2
    error('This script is designed for at least two stimulation sets.');
end

%% === Collect and combine significant channels ===
sig_all_sets = false(nCh_total, nSets);

for setIdx = 1:nSets
    amps_thisSet = FR_summary(setIdx).amps;
    amp_idx = find(abs(amps_thisSet - amp_to_plot) < 1e-6);  % tolerance to match

    if isempty(amp_idx)
        fprintf('⚠️  Amplitude %.2f µA not found in Set %d — skipping.\n', amp_to_plot, setIdx);
        continue;
    end

    if amp_idx <= size(FR_summary(setIdx).sig, 2)
        sig_idx = find(FR_summary(setIdx).sig(:, amp_idx));
        sig_all_sets(sig_idx, setIdx) = true;
        fprintf('Set %d: found %d significant channels at %.2f µA.\n', ...
            setIdx, numel(sig_idx), amp_to_plot);
    end
end

sig_combined = any(sig_all_sets, 2); % logical OR (union)
n_sig_total = sum(sig_combined);

fprintf('\nCombined significant channels across sets: %d total\n', n_sig_total);
fprintf('Channels: %s\n', num2str(find(sig_combined)'));

%% === Assign channels to shanks ===
fprintf('\nAssigning channels into %d shanks, %d channels per shank...\n', nShank, chPerShank);
channel_groups = cell(1, nShank);
ch_counter = 1;
for s = 1:nShank
    start_idx = ch_counter;
    end_idx = min(start_idx + chPerShank - 1, nCh_total);
    channel_groups{s} = start_idx:end_idx;
    ch_counter = end_idx + 1;
end

%% === Create coordinate grid for plotting ===
spacing_y = 1;   % vertical spacing (arbitrary)
spacing_x = 4;   % distance between shanks
x_pos = [];
y_pos = [];
for s = 1:nShank
    nCh_shank = numel(channel_groups{s});
    y_vals = (1:nCh_shank) * spacing_y;
    x_vals = ones(1, nCh_shank) * (s-1) * spacing_x;
    x_pos = [x_pos, x_vals];
    y_pos = [y_pos, y_vals];
end

%% === Plot combined figure ===
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.4, 0.3, 0.22, 0.55]);
hold on;

stimCh_all = unique(uniqueComb(uniqueComb > 0));
title(sprintf('%s\n%.1f µA | Combined Stim Ch %s', ...
    stim_type, amp_to_plot, num2str(stimCh_all)), ...
    'FontSize', 12, 'FontWeight', 'bold');

axis equal; axis off;
label_offset = 1.2;
label_font = 9;

for s = 1:nShank
    for ch = channel_groups{s}
        idx = find(ch == [channel_groups{:}]);
        if isempty(idx), continue; end
        x = x_pos(idx);
        y = y_pos(idx);

        if sig_combined(ch)
            markerColor = 'r';
        else
            markerColor = [0.75 0.75 0.75];
        end

        plot(x, y, 'o', 'MarkerSize', 9, ...
            'MarkerFaceColor', markerColor, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        text(x + label_offset, y, sprintf('Ch%d', ch), ...
            'FontSize', label_font, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end
end

% --- Draw tapered triangular border for the shank ---
for s = 1:nShank
    nCh_shank = numel(channel_groups{s});
    tip_y = -2;
    base_y = (nCh_shank + 1) * spacing_y;
    half_width_base = 1.4;
    half_width_tip = 0.2;
    x_center = (s-1) * spacing_x;
    x_border = [x_center-half_width_base, x_center+half_width_base, ...
                x_center+half_width_tip, x_center-half_width_tip];
    y_border = [base_y, base_y, tip_y, tip_y];
    patch(x_border, y_border, [0.85 0.85 0.85], ...
        'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 0.25, 'LineWidth', 1);
end

% --- Adjust axis limits for compact view ---
margin = 1;
xlim([min(x_pos)-2, max(x_pos)+2]);
ylim([-2, max(y_pos)+margin]);
hold off;

fprintf('\nCombined figure generated successfully.\n');