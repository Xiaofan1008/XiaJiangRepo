% Quick check of the Triggers and file size
clear 
close all;
clc;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
% ---- Paths / filenames ---- % 
% filepath = pwd;
filepath ='/Volumes/MACData/Data/Data_Xia/DX002/Exp1_Sim_Xia_250808_132343';
%% Trigger check 
clear_trig_check;

%% File size check 
File_size_check;

% Simple consistency check: last trigger time should be <= file duration 
if trig_ms(end) > duration_ms + 0.5   % small tolerance
    warning('Last trigger (%.1f ms) exceeds file duration (%.1f ms). Check sync.', trig_ms(end), duration_ms);
end

%% Artifact plot
Artifact_check(25,[-2,10]);

%% Bandwidth information 
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
% Print the values
fprintf('Lower Bandwidth: %.2f Hz\n', frequency_parameters.actual_lower_bandwidth);
fprintf('Upper Bandwidth: %.2f Hz\n', frequency_parameters.actual_upper_bandwidth);
fprintf('Notch Filter Frequency: %.2f Hz\n', frequency_parameters.notch_filter_frequency);
