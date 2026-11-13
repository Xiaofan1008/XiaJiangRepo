%% ============================================================
%   IDEAL PSTH ILLUSTRATION
%   Generates a smooth bell curve to explain Peak and Duration
%% ============================================================
clear; 

% 1. Define the Gaussian (Bell Curve) parameters
x = linspace(-50, 50, 500); % Time axis (e.g., milliseconds)
mu = 0;                     % Center of the peak
sigma = 10;                 % Width/Spread of the curve
A = 100;                    % Peak amplitude (e.g., Spikes/s)

save_figure = 1;
% Calculate the curve
y = A * exp(-(x - mu).^2 / (2 * sigma^2));

% 2. Setup the Figure
% figure('Color', 'w', 'Position', [200, 200, 550, 400]);
fig1 = figure('Units', 'centimeters', 'Position', [2, 2, 8.8, 8], 'Color', 'w');
hold on;

% 3. Plot the Baseline (Gray solid line)
plot([-50 50], [0 0], 'Color', [0.6 0.6 0.6], 'LineWidth', 2);

% 4. Plot the Ideal PSTH Curve (Thick Yellow line)
% plot(x, y, 'Color', [0.95 0.8 0.2], 'LineWidth', 4);
plot(x, y, 'Color', [0.6 0.6 0.6], 'LineWidth', 3);

% 5. Add "Peak" Annotation (Dashed vertical line)
% plot([mu mu], [0 A], '--', 'Color', [0.95 0.8 0.2], 'LineWidth', 2);
% scatter(mu, A, 50, [0.95 0.8 0.2], 'filled'); % Dot at the top

% plot([mu mu], [0 A], '--', 'Color', [0 0 0], 'LineWidth', 1);
% scatter(mu, A, 50, [0 0 0], 'filled'); % Dot at the top
% text(mu + 2, A, 'Peak Spike Rate', 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');

% 6. Add "Duration" Annotation (Horizontal double-arrow)
% Let's use the Full-Width at Half-Maximum (FWHM) as the "Duration" example
half_max = A / 2;
x_left = mu - sqrt(-2 * sigma^2 * log(0.5));
x_right = mu + sqrt(-2 * sigma^2 * log(0.5));

% Draw the line
plot([x_left+3 x_right-3], [half_max half_max], '-k', 'LineWidth', 0.75);

% Draw the left and right arrowheads
plot(x_left+3, half_max, '<k', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
plot(x_right-3, half_max, '>k', 'MarkerSize', 7, 'MarkerFaceColor', 'k');

% Add the text
% text(mu, half_max + 6, 'Response Duration', 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontName', 'Arial');

% 7. Clean up the aesthetics (Remove standard axes to look like a diagram)
axis off; 
axis square;
ylim([-10 110]); % Give it some breathing room


% 8. SAVE THE DIAGRAM
if save_figure
    save_dir = '/Users/xiaofan/Desktop/PhD Study/Paper/IEEE_TBME/Figures/Figure3/Sample';
    
    % Ensure the directory exists
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    % Set the filename
    filename = fullfile(save_dir, 'Ideal_PSTH_Diagram.tiff');
    
    % Export at 600 DPI for IEEE Publication
    exportgraphics(gcf, filename, 'Resolution', 600, 'BackgroundColor', 'w');
    
    fprintf('>>> Diagram successfully saved to:\n%s\n', filename);
end