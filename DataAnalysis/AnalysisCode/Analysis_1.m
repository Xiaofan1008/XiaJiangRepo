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

%% Spike look through all channels 
starttrial=1;
trialjump=1;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
[IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
save('IDstruct.mat', 'IDstruct','baslinespikestruct')

% find max and min values

nSamples = size(sp{find(~cellfun(@isempty,sp),1)},2)-1;
time_axis = (0:nSamples-1)*(1000/FS); % time axis in ms

nRows = 16;
nCols = 4;
left_margin = 0.02;
right_margin = 0.01;
top_margin = 0.04;
bottom_margin = 0.04;
h_spacing = 0.005;
v_spacing = 0.005;
plot_width = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
plot_height = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
figure('Units','normalized','Position',[0 0 1 1],'Color','w'); % fullscreen
global_max = -inf;
global_min = inf;
% layout = reshape(1:64,16,4)';
% find average spikes
for i = 1:nChn
    if ~isempty(sp{ch})
       ave_spike(i,:) = mean(sp{ch}(:,2:end),1); % average spikes 
       ymin = min(ave_spike(i,:));
       ymax = max(ave_spike(i,:));
       global_max = max(global_max,ymax);
       global_min = min(global_min,ymin);
        
    end
end
row = 16:-1:1;
column = [1,4,2,3]; 
for ch = 1:nChn
    % subplot(16,4,ch);
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
    

    % row = ceil(ch / nCols);     % row index
    % col = mod(ch-1, nCols) + 1; % column index

    % left = left_margin + (col-1)*(plot_width + h_spacing);
    % bottom = 1 - top_margin - row*plot_height - (row-1)*v_spacing;

    left = left_margin + (shank-1)*(plot_width + h_spacing);
    bottom = 1 - top_margin - elc*plot_height - (elc-1)*v_spacing;

    axes('Position', [left bottom plot_width plot_height]);

    if ~isempty(sp{ch})
       % ave_spike = mean(sp{ch}(:,2:end),1); % average spikes 

       plot(time_axis,ave_spike(ch,:),'k','LineWidth',2);
       xlim([time_axis(1),time_axis(end)]);
       xlim([0 1]);
       ylim([global_min,global_max]);
    end

    % if row == nRows
    if elc == nRows    
        xlabel('Time(ms)','FontSize',7);
    else 
        set(gca,'XTickLabel',[]);
    end

    % if col == 1
    if shank == 1
        ylabel('µV','FontSize',7);
    else 
        % set(gca,'YTickLabel',[]);
    end
    text(time_axis(2), ymin, num2str(ch), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
    set(gca,'FontSize',6,'Box','off');
    
 
end
