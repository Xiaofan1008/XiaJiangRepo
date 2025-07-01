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

%% Load Raw Data 
fileinfo = dir([dName '.dat']);
nSam = fileinfo.bytes/(nChn * 2);
v_fid = fopen([dName '.dat'], 'r');
lv_fid = fopen([dName '.dat'],'r');
ntimes = ceil(fileinfo.bytes / 2 / nChn / FS / T);

Trial_to_extract = 13;

% shift to the first trigger 
T = 20; % s
samples_before = 10 * FS; % in s
samples_after = T*FS; % in s
% start_index = trig(6) - samples_before; 
start_index = trig(Trial_to_extract) - samples_before;
start_byte = start_index * nChn *2; 
fseek(v_fid, start_byte, 'bof');
total_sample = samples_before + samples_after;
data = fread(v_fid, [nChn, total_sample], 'int16') * 0.195;

% export_chn = 55;
% figure;
% time_axis = (-samples_before : samples_after-1) / FS;
% plot (time_axis, data(export_chn,:))
% title(sprintf('Trial %d - data in CHN %d', Trial_to_extract, export_chn));
% xlabel("Time (ms)")
% ylabel("Amplitude (µA)")

% plot Trigger data
export_chn = 53;
figure;
start_time = -1; % in ms
end_time = 5; % in ms

start_offset = round(start_time * FS/1000);
end_offset = round(end_time * FS/1000);
trigger_sample = samples_before;
window_range = (trigger_sample + start_offset) : (trigger_sample + end_offset);
time_axis = (start_offset : end_offset)/FS*1000;
plot(time_axis,data(export_chn,window_range))
title(sprintf('Trial %d - Trigger data in CHN %d', Trial_to_extract, export_chn));
xlabel("Time (ms)")
ylabel("Amplitude (µA)")