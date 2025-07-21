function RasterPlot_amp(FS,nChn,pre_ms,post_ms,chn)
%% Raster Plot over all amplitudes
% FS - Sampling Rate (30kHz normally)
% nChn - number of channels (assume 64 for now, update more later)
% pre_ms & post_ms: Timing window
% chn - specific channel
% Raster plot across all amplitude


TrialParams = loadTrialParams;
StimParams = loadStimParams;
trig = loadTrig(0);
files = dir('*.sp.mat'); 
if ~isempty(files)
    load(files(1).name);  % Load the first matching file
else
    error('No .sp.mat files found in the current folder.');
end 
d = Depth;

stim_amp = cell2mat(StimParams(2:end,16));
TrialParams(:,4) = num2cell(stim_amp);
unique_amp = unique(stim_amp);

thisChn = d(chn); % channel file map 
sp_chn = sp{thisChn};
sp_times_ms = sp_chn(:,1); % spike time in ms

trig_ms = trig / (FS/1000);
nTrials = length(trig);
numAmp = length(unique_amp);

figure('Color','w');

% Loop through each amp 
for ii = 1:numAmp
    amp = unique_amp(ii);
    % find trials with this amplitude
    trial_idx = find(stim_amp == amp);
    trial_idx = ceil(trial_idx/2);
    trial_idx = unique(trial_idx);

    raster_x = [];
    raster_y = [];

    % loop over trials
    for t = 1:length(trial_idx)
        t0 = trig_ms(trial_idx(t));
        sp_idx = find(sp_times_ms >= (t0 - pre_ms) & sp_times_ms <= (t0 + post_ms));
        sp_in_trial = sp_times_ms(sp_idx) - t0;
        raster_x = [raster_x; sp_in_trial];
        raster_y = [raster_y; t * ones(size(sp_in_trial))];
    end

    % Plot raster for this amplitude
    subplot(numAmp, 1, ii); hold on; yyaxis right;
    set(gca, 'YColor', 'none');
    xvals = [raster_x(:)'; raster_x(:)'];
    yvals = [raster_y(:)'-0.6; raster_y(:)'+0.6];
    plot(xvals, yvals, 'k', 'LineStyle', '-', 'Marker', 'none') 
    ylim([0 length(trial_idx)+1]);

    smoothing = 5;
    Z = hist(raster_x, (-pre_ms): post_ms);
    window = normpdf((-3*smoothing:3*smoothing),0,smoothing);
    % rate = (1000/n_REP_true)*conv(Z,window);
    rate = (1000/length(trial_idx))*conv(Z,window);
    rate = rate(3*smoothing:end-3*smoothing-1);
    yyaxis left

    xlim([-pre_ms post_ms]);
    if (max(rate) >=150)
        plot((-pre_ms:post_ms),rate, 'r', 'LineWidth',2)
    else 
        ylim([0 150])
        plot((-pre_ms:post_ms),rate, 'k', 'LineWidth',2)
    end    
    title(sprintf('Amplitude = %.1f', amp));
    hold off
end

han = axes('Position',[0 0 1 1],'Visible','off','Units','normalized');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
han.Title.Visible = 'on';
xlabel(han, 'Time from Stimulation (ms)', ...
    'FontSize', 12, 'Units','normalized', ...
    'Position',[0.5, 0.02, 0]);
ylabel(han, 'Firing Rate (spikes/s)', ...
    'FontSize', 12, 'Units','normalized', ...
    'Position',[0.02, 0.5, 0]);  % Left middle
title(han, sprintf('Raster Plot over Amplitudes - Channel %d', chn), ...
      'FontSize',12, ...
      'FontWeight','bold', ...
      'Units','normalized', ...
      'Position',[0.5, 0.97, 0]); 
uistack(han, 'top');