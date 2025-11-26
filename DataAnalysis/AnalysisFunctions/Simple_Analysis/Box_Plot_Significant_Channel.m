%% ============================================================
%  Plot Box Plot of Responding Channel Percentages
% =============================================================
%  This script visualizes the percentage of responding (significant)
%  channels for each stimulation type: Separate Ch1, Separate Ch2,
%  Simultaneous, and Sequential.
%
%  You can manually enter data from your summary tables (DX010, DX011)
%  as shown below.
% ============================================================

clear;

%% === Enter Data from Tables (example from DX010 & DX011) ===
% Each row = one stimulation set
% Each column = one stimulation type
% Order: [Separate Ch1, Separate Ch2, Simultaneous, Sequential]

DX010 = [
0.5625 0.5313 0.6563 0.5313;
0.5625 0.5938 0.4688 0.9063;
0.8125 0.6875 0.5313 0.7500;
0.4688 0.5000 0.5000 0.5625;
0.4063 0.6250 0.3438 0.7188;
0.6875 0.7188 0.3125 0.6250;
0.5000 0.4688 0.2813 0.3438
];

DX011 = [
0.7813 0.5938 0.6563 0.7500;
0.7188 0.6250 0.7500 0.8125;
0.8438 0.8125 0.8125 0.9063;
0.8125 0.8438 0.9688 0.9375;
0.8125 0.8750 0.8438 0.8438;
0.8125 0.7500 0.8438 0.9375;
0.6250 0.8438 0.8125 0.9375;
0.6875 0.7813 0.9063 0.8750
];

%% === Combine All Data ===
% allData = [DX010; DX011];  % combine datasets
allData = [DX011];
stimLabels = {'Separate Ch1', 'Separate Ch2', 'Simultaneous', 'Sequential'};

%% === Convert to percentage (optional) ===
allData = allData * 100;   % convert fractions (0–1) to percentages

%% === Plot Box Plot ===
figure; 
boxplot(allData, 'Labels', stimLabels, 'Whisker', 1.5, ...
    'Colors', lines(4), 'Symbol', 'o');

title('Responding Channel Percentage', 'FontWeight', 'bold');
ylabel('Responding Channels (%)', 'FontWeight', 'bold');
ylim([0 100]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

%% === Optional: Overlay Mean Values ===
hold on;
means = mean(allData, 'omitnan');
plot(1:4, means, 'kd-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Mean');

legend('Mean', 'Location', 'northwest');
box off;

fprintf('Box plot generated successfully.\n');