% Quick check of the Triggers and file size
clear 
close all;
clc;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/Functions/MASSIVE'));
% ---- Paths / filenames ---- % 
filepath = pwd;
%% Trigger check 
clear_trig_check;

%% File size check 
File_size_check;

% Simple consistency check: last trigger time should be <= file duration 
if trig_ms(end) > duration_ms + 0.5   % small tolerance
    warning('Last trigger (%.1f ms) exceeds file duration (%.1f ms). Check sync.', trig_ms(end), duration_ms);
end

%% Artifact plot
Artifact_check(25);
