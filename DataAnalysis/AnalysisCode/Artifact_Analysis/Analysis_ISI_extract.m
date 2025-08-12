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

Trial_to_extract = 340;

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

% plot Trigger data
export_chn = 37;
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

%% ISI plots
loadStimParams;
post_trig = cell2mat(StimParams(2:end, 6));
isi_values = post_trig(2:2:end);
unique_isi = unique(isi_values);
nISI = length(unique_isi);

% Extract all channels
isi_channels = cell2mat(TrialParams(2:end, 3)); 

isi_trial_idx = cell(nISI,1);
isi_trial_chns = cell(nISI,1);
isi_data_ch1 = cell(nISI,1); % first channel data
isi_data_ch2 = cell(nISI,1); % second channel data

samples_before = 2 * FS;
samples_after = 5 * FS;
total_sample = samples_before + samples_after + 1;

% loop over
for i = 1:nISI
   isi = unique_isi(i);
   
   % find matched isi
   match_isi = find(isi_values == isi);

   isi_trial_idx{i} = match_isi; 
   chn1 = isi_channels((match_isi-1)*2+1);
   chn2 = isi_channels((match_isi-1)*2+2);
   isi_trial_chns{i} = [chn1(:) chn2(:)];

   nTrials = length(match_isi); % numebr of trails per condition

   data_temp_ch1 = nan(nTrials, total_sample);
   data_temp_ch2 = nan(nTrials, total_sample);

   for j = 1:nTrials
        t_idx = match_isi(j);
        ch1 = chn1(j);
        ch2 = chn2(j);

        if any(new_trig(t_idx) == -500 | (new_trig(t_idx) - samples_before) < 1)
            continue
        end

        % Calculate file read position
        start_index = new_trig(t_idx) - samples_before;
        start_byte = start_index * nChn * 2;
        fseek(v_fid, start_byte, 'bof');

        % Read data
        raw = fread(v_fid, [nChn, total_sample], 'int16') * 0.195;

        % Store data
        data_temp_ch1(j, :) = raw(ch1, :);
        data_temp_ch2(j, :) = raw(ch2, :);
   end
    isi_data_ch1{i} = mean(data_temp_ch1,1,'omitnan');
    isi_data_ch2{i} = mean(data_temp_ch2,1,'omitnan');
end

%% Plot Average artifact of different distance
figure;
start_time = -1; % in ms
end_time = 3; % in ms
start_offset = round(start_time * FS/1000);
end_offset = round(end_time * FS/1000);
trigger_sample = samples_before;
window_range = (trigger_sample + start_offset) : (trigger_sample + end_offset);
time_axis = (start_offset : end_offset)/FS*1000;
hold on 
colors = turbo(length(DIST));

for i = 1:length(DIST)
    plot(time_axis,dist_all_chn1(i,window_range),'Color',colors(i,:), 'LineWidth', 2.0);
end
xlabel('Time (ms)');
ylabel('Amplitude (µV)')
title('Average Artifacts in Channel 1')
legend_strings = arrayfun(@(x) sprintf('%d µm', x), DIST, 'UniformOutput', false);
legend(legend_strings,'Location','northeast');
hold off
box off

figure;
hold on
for i = 1:length(DIST)
    plot(time_axis,dist_all_chn2(i,window_range),'Color',colors(i,:), 'LineWidth', 2.0);
end
xlabel('Time (ms)');
ylabel('Amplitude (µV)')
title('Average Artifacts in Channel 2')
legend_strings = arrayfun(@(x) sprintf('%d mm', x), DIST, 'UniformOutput', false);
legend(legend_strings,'Location','northeast');
hold off
box off
