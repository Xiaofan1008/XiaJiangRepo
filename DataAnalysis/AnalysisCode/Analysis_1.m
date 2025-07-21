%% Analysis Code
clear all
close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
addpath(genpath('/Volumes/MACData/Data/Data_Sabrina/Experimental_Design'));
% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-500; %removes spikes below this threshold
artefact_high=500; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)?
secondstoanalyse=8; %How long after the trigger do you want to analyse spikes for(ms)?
printspiking=0;
par=0;

%% Initial
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
dName='amplifier';
mem_check=dir('amplifier.dat');
T = mem_check.bytes ./ (2 * nChn * FS);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * FS);
if t_len < (Overall_time_to_analyse+Startpoint_analyse)
    error('Time to analyse exceeds the data recording period')
end
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
elseif Overall_time_to_analyse~=0
    T=Overall_time_to_analyse;
end

%% Load Trigger
% cleanTrig_sabquick;
cleanTrig_xiaquick;
trig = loadTrig(0);

%% Remove Artifact
denoiseIntan_sab(filepath, dName, T, par, Startpoint_analyse, Overall_time_to_analyse,FS);

% %% Load Raw Data 
% fileinfo = dir([dName '.dat']);
% nSam = fileinfo.bytes/(nChn * 2);
% v_fid = fopen([dName '.dat'], 'r');
% lv_fid = fopen([dName '.dat'],'r');
% ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);
% 
% Trial_to_extract = 274;
% 
% if MissingTrig == 1
%     amp = StimParams{Trial_to_extract,16};
% else 
%     amp = StimParams{Trial_to_extract*2,16};
% end
% 
% % shift to the first trigger 
% T = 20; % s
% samples_before = 10 * FS; % in s
% samples_after = T*FS; % in s
% % start_index = trig(6) - samples_before; 
% start_index = trig(Trial_to_extract) - samples_before;
% start_byte = start_index * nChn *2; 
% fseek(v_fid, start_byte, 'bof');
% total_sample = samples_before + samples_after;
% data = fread(v_fid, [nChn, total_sample], 'int16') * 0.195;
% 
% % export_chn = 55;
% % figure;
% % time_axis = (-samples_before : samples_after-1) / FS;
% % plot (time_axis, data(export_chn,:))
% % title(sprintf('Trial %d - data in CHN %d', Trial_to_extract, export_chn));
% % xlabel("Time (ms)")
% % ylabel("Amplitude (µA)")
% 
% % plot Trigger data
% export_chn = 11;
% figure;
% start_time = -1; % in ms
% end_time = 10; % in ms
% 
% start_offset = round(start_time * FS/1000);
% end_offset = round(end_time * FS/1000);
% trigger_sample = samples_before;
% window_range = (trigger_sample + start_offset) : (trigger_sample + end_offset);
% time_axis = (start_offset : end_offset)/FS*1000;
% plot(time_axis,data(export_chn,window_range))
% title(sprintf('Trial %d - Trigger data in CHN %d, Amp= %d µA', Trial_to_extract, export_chn, amp));
% xlabel('Time (ms)')
% ylabel("Amplitude (µA)")
% box off
% 
% %% FFT plot
% figure;
% data_fft = data(export_chn,window_range);
% L = length(data_fft);
% % data_fft_noDC= data_fft - mean(data_fft);
% N = length(data_fft);
% Y=fft(data_fft);  
% Y_abs = abs(Y)/N; % calculate the power
% Y_abs_s = fftshift(Y_abs);
% f = (-L/2 : L/2-1) * (FS / L);  % frequency axis
% plot(f,Y_abs_s)
% xlim([-250 250])
% xlabel('Frequency (Hz)')
% title(sprintf('FFT Plot, Trial %d - Trigger data in CHN %d, Amp= %d µA', Trial_to_extract, export_chn, amp));
% box off

%% Threshold & MU & Spikes
allExtract_sab_1(dName,filepath,T,par,artefact,artefact_high);
fclose('all');

%% Spike look through all channels (horizational version)
starttrial=1;
trialjump=1;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
[IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
save('IDstruct.mat', 'IDstruct','baslinespikestruct')

% find max and min values
files = dir('*.sp.mat'); 
if ~isempty(files)
    load(files(1).name);  % Load the first matching file
else
    error('No .sp.mat files found in the current folder.');
end
d = Depth;
nSamples = size(sp{find(~cellfun(@isempty,sp),1)},2)-1;
time_axis = (0:nSamples-1)*(1000/FS); % time axis in ms

nRows = 4; nCols = 16;
left_margin = 0.02; right_margin = 0.01; top_margin = 0.04; bottom_margin = 0.04;
h_spacing = 0.01; v_spacing = 0.01;
plot_width = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
plot_height = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
figure('Units','normalized','Position',[0 0 1 1],'Color','w'); % fullscreen
global_max = -inf;
global_min = inf;
% find average spikes
for i = 1:nChn
    if ~isempty(sp{i})
       ave_spike(i,:) = mean(sp{i}(:,2:end),1); % average spikes 
       ymin = min(ave_spike(i,:));
       ymax = max(ave_spike(i,:));
       global_max = max(global_max,ymax);
       global_min = min(global_min,ymin);        
    end
end
column = 16:-1:1;
row = [4,1,3,2]; 
for ch = 1:nChn
    row_idx = ceil(ch/16);
    col_idx = mod(ch-1,16)+1;

    shank = column(col_idx); % colomn position 
    elc = row(row_idx);% row position 

    left = left_margin + (shank-1)*(plot_width + h_spacing);
    bottom = 1 - top_margin - elc*plot_height - (elc-1)*v_spacing;

    axes('Position', [left bottom plot_width plot_height]);
    
    thisChn = d(ch);
     if ~isempty(sp{thisChn})
       plot(time_axis,ave_spike(thisChn,:),'k','LineWidth',2);
       xlim([time_axis(1),time_axis(end)]);
       xlim([0 1.5]);
       ylim([global_min,global_max]);
     end

    if elc == nRows    
        xlabel('Time(ms)','FontSize',7);
    else 
        set(gca,'XTickLabel',[]);
    end

    if shank == 1
        ylabel('µV','FontSize',7);
    else 
        % set(gca,'YTickLabel',[]);
    end
    text(time_axis(2), ymin, num2str(ch), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
    set(gca,'FontSize',6,'Box','off');
end

%% Spike look through all channels (Vertical version)
starttrial=1;
trialjump=1;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
[IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
save('IDstruct.mat', 'IDstruct','baslinespikestruct')

% find max and min values
files = dir('*.sp.mat'); 
if ~isempty(files)
    load(files(1).name);  % Load the first matching file
else
    error('No .sp.mat files found in the current folder.');
end
d = Depth;
nSamples = size(sp{find(~cellfun(@isempty,sp),1)},2)-1;
time_axis = (0:nSamples-1)*(1000/FS); % time axis in ms

nRows = 16;
nCols = 4;
left_margin = 0.02;right_margin = 0.01;top_margin = 0.04;bottom_margin = 0.04;
h_spacing = 0.01;v_spacing = 0.01;
plot_width = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
plot_height = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
figure('Units','normalized','Position',[0 0 1 1],'Color','w'); % fullscreen
global_max = -inf; global_min = inf;
% find average spikes
for i = 1:nChn
    if ~isempty(sp{i})
       ave_spike(i,:) = mean(sp{i}(:,2:end),1); % average spikes 
       ymin = min(ave_spike(i,:));
       ymax = max(ave_spike(i,:));
       global_max = max(global_max,ymax);
       global_min = min(global_min,ymin);        
    end
end
row = 16:-1:1;
column = [1,4,2,3]; 
for ch = 1:nChn    
    temp = mod(ch,16);
    if temp == 0
        loc_col = fix(ch/16);
        loc_row = 16; 
    else 
        loc_col = fix(ch/16)+1;
        loc_row = temp;
    end
    
    shank = column(loc_col); % colomn position 
    elc = row(loc_row);% row position 

    left = left_margin + (shank-1)*(plot_width + h_spacing);
    bottom = 1 - top_margin - elc*plot_height - (elc-1)*v_spacing;

    axes('Position', [left bottom plot_width plot_height]);
    
    thisChn = d(ch);
     if ~isempty(sp{thisChn})
       plot(time_axis,ave_spike(thisChn,:),'k','LineWidth',2);
       xlim([time_axis(1),time_axis(end)]);
       xlim([0 1.5]);
       ylim([global_min,global_max]);
     end

    % if row == nRows
    if elc == nRows    
        xlabel('Time(ms)','FontSize',7);
    else 
        set(gca,'XTickLabel',[]);
    end

    if shank == 1
        ylabel('µV','FontSize',7);
    else 
        % set(gca,'YTickLabel',[]);
    end
    text(time_axis(2), ymin, num2str(ch), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
    set(gca,'FontSize',6,'Box','off');
end

%% Raster Plot 
TrialParams = loadTrialParams; d = Depth; AMP = loadAMP;
loadnREPtrue; loadDUALSTIM;
n_REP = n_REP_true;
% Initialisation of data strucutre
if ~isempty(TrialParams)
    nTr = size(TrialParams,1)/n_REP/(DUALSTIM+1);
else
    nTr = length(trig);
end 
nT = zeros(n_REP*(DUALSTIM+1),nTr);
TrialParams = cell2mat(TrialParams);

RATE = 1;
trig_ms = trig./(FS/1000); 
for i = 1:nTr
    if ~isempty(TrialParams)
        nT(:,i) = TrialParams(TrialParams(:,2) == i,1);        
    else
        nT(i,:) = i;
    end
end
if (DUALSTIM == 1)
    nT = nT(1:2:end,:);
end

trial = 5;
r = trial; 
MR = 500;
XLIM = [-20 150];

col_order = [1 4 2 3];
row = 16:-1:1;
nRows = 16; nCols = 4;
left_margin = 0.02; right_margin = 0.01; top_margin = 0.04; bottom_margin = 0.04;
h_spacing = 0.02; v_spacing = 0.02;
plot_width = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
plot_height = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
figure('Units','normalized','Position',[0 0 1 1],'Color','w'); % fullscreen
for i = 1:nChn
    chn = d(i);
    sp_chn = sp{chn};
    SPIKES = cell(1,nTr);
       
    for t = 1:n_REP
        thisSp = sp_chn(sp_chn(:,1) >= (trig_ms(nT(t,r)) - 500) & sp_chn(:,1) <= (trig_ms(nT(t,r)) + 500),:);
        SPIKES{t} = thisSp(:,1) - (trig_ms(nT(t,r)) - 500); % time 0 in SPIKES = trig_ms - 500ms
    end    
    row_idx = row(mod(i-1, 16) + 1);
    col_idx = col_order(ceil(i / 16));
    left = left_margin + (col_idx-1)*(plot_width + h_spacing);
    bottom = 1 - top_margin - row_idx*plot_height - (row_idx-1)*v_spacing;
    axes('Position', [left bottom plot_width plot_height]);
    if r > 1
        % psth(SPIKES(:),[-500 500],1,MR,[],RATE,MR*(r-2));
        psth(SPIKES(:),[-500 500],1,MR,[],RATE,0);
        line(XLIM,[MR*(r-2) MR*(r-2)],'Color','b');
    else
        % psth(SPIKES(:),[-500 500],1,MR,[],RATE,MR*(r-1));     
        psth(SPIKES(:),[-500 500],1,MR,[],RATE,0); 
        line(XLIM,[MR*(r-1) MR*(r-1)],'Color','b');
    end

    xlim(XLIM); 
    text(1, MR*0.85, sprintf('Ch %02d', i), 'FontSize', 7);
    line([3 3],[0 MR*(r-1)],'Color','r');
    line([8 8],[0 MR*(r-1)],'Color','r');
end

% han = axes(gcf, 'visible', 'off'); 
% han.XLabel.Visible = 'on';
% han.YLabel.Visible = 'on';
% xlabel(han, 'Time from Stimulation (ms)');
% ylabel(han, 'Firing Rate (spikes/s)');

xlabel('Time from Stimulation (ms)');
ylabel('Firing Rate (spikes/s)');



%% Raster Plot (Hand coded)
TrialParams = loadTrialParams; d = Depth; AMP = loadAMP;

trig = loadTrig(0);
loadnREPtrue; loadDUALSTIM;

chn = 44;
thischn = d(chn);
FS = 30000;
pre_ms = 50; post_ms = 150; % in ms
trig_ms = trig/(FS/1000); % convert trigger to ms
% trig_ms = trig_ms - 20*1000;
% convert to samples
pre_samp = pre_ms * FS / 1000;
post_samp = post_ms * FS / 1000;

% Spike timestaps and waveforms
sp_chn = sp{thischn};
sp_times_ms = sp_chn(:,1); % timestampes in ms
waveforms = sp_chn(:,2:end);

% time axis for waveform
nSamples = size(waveforms,2);
timeaxis = (0:nSamples-1)*(1000/FS); % in ms

%Initiialize raster and MUA
nTrials = length(trig_ms);
raster_x = []; 
raster_y = [];
% MUA_trials = zeros(nTrials,nSamples);

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
    if (max(rate)<= 40)
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
xlabel(han, 'Time from Stimulation (ms)', ...
    'FontSize', 12, 'Units','normalized', ...
    'Position',[0.5, 0.02, 0]);
ylabel(han, 'Firing Rate (spikes/s)', ...
    'FontSize', 12, 'Units','normalized', ...
    'Position',[0.02, 0.5, 0]);  % Left middle
uistack(han, 'top');

% Loop over tials
for t = 1:nTrials
    t0 = trig_ms(t);
    sp_idx = find(sp_times_ms >= (t0-pre_ms) & sp_times_ms <= (t0+post_ms));
    sp_in_trial = sp_times_ms(sp_idx) - (t0);
    raster_x = [raster_x; sp_in_trial];
    raster_y = [raster_y; t*ones(size(sp_in_trial))];

    % % Compute average waveform per trial
    % if ~isempty(sp_idx)
    %     MUA_trials(t,:) = mean(waveforms(sp_idx,:),1);
    % end
end

% compute avergae MUA across trials
% avg_MUA = mean(MUA_trials,1);

% Plotting
% figure('Color','w');

% % Raster plot
% subplot(2,1,1); hold on;
% for i = 1:length(raster_x)
%     line([raster_x(i) raster_x(i)], [raster_y(i)-0.4 raster_y(i)+0.4], 'Color', 'k');
% end
% xlim([-pre_ms post_ms]);
% ylim([0 nTrials+1]);
% ylabel('Trial');
% title(['Raster Plot - Channel ' num2str(thischn)]);
% 
% % MUA plot (average spike shape)
% subplot(2,1,2); hold on;
% plot(time_axis, avg_MUA, 'k', 'LineWidth', 1.5);
% xlabel('Time relative to spike (ms)');
% ylabel('Average µV');
% title('Average MUA Waveform');

% Raster plot
% hold on;
% 
% % yyaxis left
% % plot(time_axis, avg_MUA, 'k', 'LineWidth', 1.5);
% yyaxis right
% 
% 
% for i = 1:length(raster_x)
%     line([raster_x(i) raster_x(i)], [raster_y(i)-0.4 raster_y(i)+0.4], 'Color', 'k');
% end
% % plot(time_axis, avg_MUA, 'k', 'LineWidth', 1.5);
% xlabel('Time relative to spike (ms)');
% ylabel('Average µV');
% title('Average MUA Waveform');

figure;
yyaxis right
hold on
ax = gca;
ax.YAxis(2).Color = 'k';
xvals = [raster_x(:)'; raster_x(:)'];
yvals = [raster_y(:)'-0.4; raster_y(:)'+0.4];
plot(xvals, yvals, 'k', 'LineStyle', '-', 'Marker', 'none')  % <-- enforce no markers

% plot(time_axis, avg_MUA, 'k', 'LineWidth', 1.5);

smoothing = 1;
Z = hist(raster_x, (-pre_ms): post_ms);
window = normpdf((-3*smoothing:3*smoothing),0,smoothing);
rate = (1000/n_REP_true)*conv(Z,window);
rate = rate(3*smoothing:end-3*smoothing-1);
% figure(3);
yyaxis left
plot((-pre_ms:post_ms),rate, 'k', 'LineWidth',2)
xlabel('Time relative to stimulation (ms)');