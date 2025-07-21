function RasterPlot_trials(FS,nChn,pre_ms,post_ms)
%% Raster Plot over all trials
% FS - Sampling Rate (30kHz normally)
% nChn - number of channels (assume 64 for now, update more later)
% mode - 0:Horizontal version (rows:shanks columns:electrodes) 1: Vertical
% version (rows:electrodes, columns: shanks)
% Spike plot across all channels (Horizonal version or vertical version)

% load files
TrialParams = loadTrialParams; d = Depth; AMP = loadAMP;
trig = loadTrig(0);
loadnREPtrue; loadDUALSTIM;

files = dir('*.sp.mat'); 
if ~isempty(files)
    load(files(1).name);  % Load the first matching file
else
    error('No .sp.mat files found in the current folder.');
end

% units in ms
trig_ms = trig/(FS/1000); % convert trigger to ms
pre_samp = pre_ms * FS / 1000;
post_samp = post_ms * FS / 1000;

% find stimulation channels
stim_chns = cell2mat(TrialParams(:,3));
unique_stim_chns = unique(stim_chns);
unique_stim_chns(unique_stim_chns == 0) = [];

% Spike timestaps and waveforms
% sp_chn = sp{thischn};
% sp_times_ms = sp_chn(:,1); % timestampes in ms
% waveforms = sp_chn(:,2:end);

% time axis for waveform
% nSamples = size(waveforms,2);
% timeaxis = (0:nSamples-1)*(1000/FS); % in ms

%Initiialize raster and MUA
nTrials = length(trig_ms);
raster_x = []; 
raster_y = [];

% Plot setting
col_order = [1 4 2 3];
row = 16:-1:1;
nRows = 16; nCols = 4;
left_margin = 0.03; right_margin = 0.01; top_margin = 0.04; bottom_margin = 0.04;
h_spacing = 0.02; v_spacing = 0.02;
plot_width = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
plot_height = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
figure('Units','normalized','Position',[0 0 1 1],'Color','w'); % fullscreen

% ALl channels
for i= 1:nChn
    thischn = d(i);
    sp_chn = sp{thischn};
    sp_times_ms = sp_chn(:,1); % timestampes in ms
    waveforms = sp_chn(:,2:end);
    raster_x = []; 
    raster_y = [];

    % Loop over tials
    for t = 1:nTrials
        t0 = trig_ms(t);
        sp_idx = find(sp_times_ms >= (t0-pre_ms) & sp_times_ms <= (t0+post_ms));
        sp_in_trial = sp_times_ms(sp_idx) - (t0);
        raster_x = [raster_x; sp_in_trial];
        raster_y = [raster_y; t*ones(size(sp_in_trial))];
    end

    row_idx = row(mod(i-1, 16) + 1);
    col_idx = col_order(ceil(i / 16));
    left = left_margin + (col_idx-1)*(plot_width + h_spacing);
    bottom = 1 - top_margin - row_idx*plot_height - (row_idx-1)*v_spacing;
    axes('Position', [left bottom plot_width plot_height]);

    yyaxis right
    set(gca, 'YColor', 'none');
    hold on
    % ax = gca;
    % ax.YAxis(2).Color = 'k';
    xvals = [raster_x(:)'; raster_x(:)'];
    yvals = [raster_y(:)'-0.6; raster_y(:)'+0.6];
    plot(xvals, yvals, 'k', 'LineStyle', '-', 'Marker', 'none') 
  

    smoothing = 5;
    Z = hist(raster_x, (-pre_ms): post_ms);
    window = normpdf((-3*smoothing:3*smoothing),0,smoothing);
    % rate = (1000/n_REP_true)*conv(Z,window);
    rate = (1000/nTrials)*conv(Z,window);
    rate = rate(3*smoothing:end-3*smoothing-1);
    yyaxis left
    % plot((-pre_ms:post_ms),rate, 'k', 'LineWidth',2)
    if ismember(i, unique_stim_chns)
        plot((-pre_ms:post_ms), rate, 'b', 'LineWidth', 2)
    elseif (max(rate)<= 40)
        ylim([0 50]);
        plot((-pre_ms:post_ms),rate, 'k', 'LineWidth',2)
    else
        plot((-pre_ms:post_ms),rate, 'r', 'LineWidth',2)
    end
    hold off

    text(-45, 50, sprintf('Ch %02d', i), 'FontSize', 10);
end

% Create invisible overall axes for global labels
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
title(han, 'Raster Plot over All trials', ...
      'FontSize',12, ...
      'FontWeight','bold', ...
      'Units','normalized', ...
      'Position',[0.5, 0.97, 0]); 
uistack(han, 'top');