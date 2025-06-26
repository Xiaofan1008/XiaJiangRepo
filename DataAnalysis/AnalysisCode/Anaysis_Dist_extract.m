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

%% Electrode seperation plots
% setting
DIST = [50,100,150,200,300,400];    % differnet distances
% Chn_50 =    [40,39; 56,55; 38,37];    % 50µm  channels
% Chn_100 =   [41,39; 55,53; 58,56];    % 100µm channels
% Chn_150 =   [40,37; 56,53; 44,41];    % 150µm channels
% Chn_200 =   [43,39; 57,53; 42,58];    % 200µm channels
% Chn_300 =   [43,37; 59,53; 57,51];    % 300µm channels
% Chn_400 =   [44,36; 60,52; 9,57];     % 400µm channels

Chn_50 =    [40,39];    % 50µm  channels
Chn_100 =   [41,39];    % 100µm channels
Chn_150 =   [40,37];    % 150µm channels
Chn_200 =   [43,39];    % 200µm channels
Chn_300 =   [43,37];    % 300µm channels
Chn_400 =   [44,36];     % 400µm channels
Chn_map = {Chn_50,Chn_100,Chn_150,Chn_200,Chn_300,Chn_400};

samples_before = 2 * FS; % seconds before trigger
samples_after = 5*FS; % seconds after trigger
total_sample = samples_after + samples_before + 1;
loadStimParams;

% extract channels
channels = cell2mat(TrialParams(2:end,3));
channels = reshape(channels,2,[]).';
nTrials = length(channels)/2; % number of trials
ch1_all = channels(1:2:end); % channel 1 
ch2_all = channels(2:2:end); % channel 2

dist_all_chn1 = nan(length(DIST),total_sample);
dist_all_chn2 = nan(length(DIST),total_sample);
% Loop through each distance conditions
for d = 1:length(DIST)
    pairs = Chn_map{d}; 

    % Find the matched trial index
    valid_idx = find(ismember(channels,pairs,'rows'));
    
    % Preallocate data
    nValid = length(valid_idx);
    data_dist_chn1 = nan(nValid,total_sample);
    data_dist_chn2 = nan(nValid,total_sample);

    % Extract data from valid trials
    for i = 1:nValid
        trial_idx = valid_idx(i);

        % Check trigger
        if any((new_trig(trial_idx) == -500) | ((new_trig(trial_idx) - samples_before) < 1))
            continue;
        end

        start_index = new_trig(trial_idx) - samples_before;
        start_byte = start_index * nChn *2; 
        fseek(v_fid, start_byte, 'bof');
        raw = fread(v_fid, [nChn, total_sample], 'int16') * 0.195;

        chA = ch1_all(trial_idx);
        chB = ch2_all(trial_idx);
        data_dist_chn1(i,:) = raw(chA,:);
        data_dist_chn2(i,:) = raw(chB,:);
        
    end
    % artifacts_by_dist{d} = data_dist;
    dist_all_chn1(d,:) = mean(data_dist_chn1,1,'omitnan');
    dist_all_chn2(d,:) = mean(data_dist_chn2,1,'omitnan');
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