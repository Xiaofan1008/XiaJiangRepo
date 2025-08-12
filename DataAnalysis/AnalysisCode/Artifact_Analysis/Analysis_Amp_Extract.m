%% Analysis Code
clear all
close all
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
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

%% Load Raw Data & Plot
fileinfo = dir([dName '.dat']);
nSam = fileinfo.bytes/(nChn * 2);
v_fid = fopen([dName '.dat'], 'r');
lv_fid = fopen([dName '.dat'],'r');
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);


Trial_to_extract = 40;
if MissingTrig == 1
    amp = StimParams{Trial_to_extract,16};
else 
    amp = StimParams{Trial_to_extract*2,16};
end
% shift to the first trigger 
T = 20; % s
samples_before = 10 * FS; % in s
samples_after = T*FS; % in s
start_index = trig(Trial_to_extract) - samples_before; 
start_byte = start_index * nChn *2; 
fseek(v_fid, start_byte, 'bof');

total_sample = samples_before + samples_after;
data = fread(v_fid, [nChn, total_sample], 'int16') * 0.195;

% plot 
export_chn = 53;
figure;
start_time = -2; % in ms
end_time = 1000; % in ms

start_offset = round(start_time * FS/1000);
end_offset = round(end_time * FS/1000);
trigger_sample = samples_before;
window_range = (trigger_sample + start_offset) : (trigger_sample + end_offset);
time_axis = (start_offset : end_offset)/FS*1000;
plot(time_axis,data(export_chn,window_range))
title(sprintf('Trial %d - Trigger data in CHN %d, Amp= %d µA', Trial_to_extract, export_chn, amp));
xlabel('Time (ms)')
ylabel("Amplitude (µA)")
box off


%% FFT plot
figure;
data_fft = data(export_chn,window_range);
L = length(data_fft);
% data_fft_noDC= data_fft - mean(data_fft);
N = length(data_fft);
Y=fft(data_fft);  
Y_abs = abs(Y)/N; % calculate the power
Y_abs_s = fftshift(Y_abs);
f = (-L/2 : L/2-1) * (FS / L);  % frequency axis
plot(f,Y_abs_s)
xlim([-500 500])
xlabel('Frequency (Hz)')
title(sprintf('FFT Plot, Trial %d - Trigger data in CHN %d, Amp= %d µA', Trial_to_extract, export_chn, amp));
box off

%% Amp trials extract
AMP = [-1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; % Actual AMP setting for the experiment
AMP_ref = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; % Actual AMP setting for the experiment

T = 5; % in s
samples_before = 2 * FS; % seconds before trigger
samples_after = T*FS; % seconds after trigger
total_sample = samples_after + samples_before + 1;

% Extract the amp setting for each trial
amps = cell2mat(StimParams(2:end,16));

% Extract channel for each trial
channels = cell2mat(TrialParams(2:end,3));

% Extract data 
nAmps = length(AMP);
% amp_data = nan(nTrials*2, total_sample);
amp_ave_all = nan(nAmps,total_sample);

for a = 1:nAmps
    i_amp = AMP(a); 

    % Find the corresponding trials and channels
    idx_amp = find(amps == i_amp); % trials
    idx_trig = idx_amp(2:2:end)/2; % trigger index (2 rows for 1 trial)
    idx_chn = channels(idx_amp); % channels

    nTrials = length(idx_trig);
    amp_data = nan(nTrials*2, total_sample);

    for i = 1:nTrials
        t0 = idx_trig(i); % trigger

        ch_1 = channels((i-1)*2+1); % first channel
        ch_2 = channels(2*i); % second channel

        % read data
        start_index = new_trig(t0) - samples_before; 
        start_byte = start_index * nChn *2;
        fseek(v_fid, start_byte, 'bof'); % move pointer
        data = fread(v_fid, [nChn, total_sample], 'int16') * 0.195; % all channel data

        % extract corresponding channel
        data_ch1 = data(ch_1,:); % channel 1 data
        data_ch2 = data(ch_2,:); % channel 2 data

        % store data
        amp_data((i-1)*2+1,:) = data_ch1;
        amp_data(i*2,:) = data_ch2;
    end 

    % average data 
    amp_ave_all(a,:) = mean(amp_data,1);
end

%% Plot Average Artifacts
figure;
start_time = -2; % in ms
end_time = 300; % in ms
start_offset = round(start_time * FS/1000);
end_offset = round(end_time * FS/1000);
trigger_sample = samples_before;
window_range = (trigger_sample + start_offset) : (trigger_sample + end_offset);
time_axis = (start_offset : end_offset)/FS*1000;
hold on 
colors = lines(length(AMP));
for i = 1:length(AMP)
    plot(time_axis,amp_ave_all(i,window_range),'Color',colors(i,:), 'LineWidth', 2.0);
end
xlabel('Time (ms)');
ylabel('Amplitude (µV)')
title('Average Artifacts')
legend_strings = arrayfun(@(x) sprintf('%d µA', x), AMP_ref, 'UniformOutput', false);
legend(legend_strings,'Location','northeast');
hold off
box off

%% Plot average artifact across all trials 
amp_all = mean(amp_ave_all,1);
figure;
start_time = -2; % in ms
end_time = 300; % in ms
start_offset = round(start_time * FS/1000);
end_offset = round(end_time * FS/1000);
trigger_sample = samples_before;
window_range = (trigger_sample + start_offset) : (trigger_sample + end_offset);
time_axis = (start_offset : end_offset)/FS*1000;
plot(time_axis,amp_all(window_range), 'LineWidth', 2.0);
xlabel('Time (ms)');
ylabel('Amplitude (µV)')
title('Average Artifacts (All Trials)')
box off

%% Plot trigger time difference
trig_diff = diff(trig)/FS; % time difference between adjacent triggers
figure;
plot(trig_diff);
xlabel('Triger time difference index')
ylabel('Time(s)')
title("Trigger Time Differences")