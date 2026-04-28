%% ========================================================================
%  Batch Manual Override for Responding Channels
%  Safely forces specific noisy channels to SILENT, or rescued channels
%  to RESPONDING, across multiple experimental conditions at once.
%  Automatically creates a safe backup before overwriting.
% ========================================================================
clear;

%% ========================================================================
%  1. FILE PATH SETTINGS
% ========================================================================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX019/Xia_ISI_SimSeq1';

%% ========================================================================
%  2. OVERRIDE CONTROL PANEL
%  Add as many blocks as you need below to fix your noisy channels.
% ========================================================================
overrides = {};

%% ================ Condition 1 =============== 
% --- Oms
idx = 1;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;          % 10 µA
overrides{idx}.target_ptd_ms  = 0;           % 0 ms (Simultaneous)
overrides{idx}.force_silent   = [1:64]; % Channels to turn OFF (noise)
overrides{idx}.force_respond  = [58,63];          % Channels to turn ON (rescue)

% --- 3ms
idx = 2;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 3;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [11:15,19:21,23:26];        % Turn ON

% ---4ms 
idx = 3;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 4;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [11:15,19:21,23:26];        % Turn ON

% --- 5ms
idx = 4;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 5;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [11:15,19:21,23:26];        % Turn ON

% --- 6ms
idx = 5;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 6;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,6,10:15,19:21,23:26];        % Turn ON

% --- 7ms
idx = 6;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 7;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,3,9:16,19:21,23:26];

% --- 8ms
idx = 7;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 8;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,6,10:15,19:21,23:26];

% --- 9ms
idx = 8;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 9;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,6,10:15,19:21,23:26];

% --- 10ms
idx = 9;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 10;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,3,6,11:15,19:21,23:26];

% --- 11ms
idx = 10;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 11;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,6,10:15,19:21,23:26];

% --- 12ms
idx = 11;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 12;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,3,10:15,19:21,23:26];

% --- 13ms
idx = 12;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 13;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,10:15,19:21,23:26];

% --- 14ms
idx = 13;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           %  µA
overrides{idx}.target_ptd_ms  = 14;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,10:15,19:21,23:26];    % Turn ON

% --- 15ms
idx = 14;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 15;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,10:15,19:21,23:26];

% --- 17ms
idx = 15;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           %  µA
overrides{idx}.target_ptd_ms  = 17;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,4,10:15,19:21,23:26];    % Turn ON

% --- 20ms
idx = 16;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 5;           %  µA
overrides{idx}.target_ptd_ms  = 20;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [1,10:15,19:21,23:26];    % Turn ON

%% ============== Condition 2 ==============
% --- Oms
idx = 17;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;          % 10 µA
overrides{idx}.target_ptd_ms  = 0;           % 0 ms (Simultaneous)
overrides{idx}.force_silent   = [1:64]; % Channels to turn OFF (noise)
overrides{idx}.force_respond  = [18:22,24:26];          % Channels to turn ON (rescue)

% --- 3ms
idx = 18;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 3;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];        % Turn ON

% ---4ms 
idx = 19;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 4;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:16,18:24,26:29];        % Turn ON

% --- 5ms
idx = 20;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 5;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];        % Turn ON

% --- 6ms
idx = 21;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 6;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,18:26];        % Turn ON

% --- 7ms
idx = 22;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 7;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:16,18:26];

% --- 8ms
idx = 23;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 8;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:16,18:26];

% --- 9ms
idx = 24;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 9;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:16,18:26];

% --- 10ms
idx = 25;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 10;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:16,18:26];

% --- 11ms
idx = 26;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 11;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:16,19:26];

% --- 12ms
idx = 27;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 12;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];

% --- 13ms
idx = 28;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 13;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];

% --- 14ms
idx = 29;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           %  µA
overrides{idx}.target_ptd_ms  = 14;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];    % Turn ON

% --- 15ms
idx = 30;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 15;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];

% --- 17ms
idx = 31;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           %  µA
overrides{idx}.target_ptd_ms  = 17;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];    % Turn ON

% --- 20ms
idx = 32;
overrides{idx}.target_set_idx = 1;           % Set 1
overrides{idx}.target_amp_uA  = 10;           %  µA
overrides{idx}.target_ptd_ms  = 20;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];    % Turn ON


%% =========== Condition 3 ============
% --- Oms
idx = 33;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;          % 10 µA
overrides{idx}.target_ptd_ms  = 0;           % 0 ms (Simultaneous)
overrides{idx}.force_silent   = [1:64]; % Channels to turn OFF (noise)
overrides{idx}.force_respond  = [13,14,19,20,21,23:25];          % Channels to turn ON (rescue)

% --- 3ms
idx = 34;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 3;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [12:15,19:20,23:26];        % Turn ON

% ---4ms 
idx = 35;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 4;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [12:15,19:20,23:26];        % Turn ON

% --- 5ms
idx = 36;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 5;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:14,19:21,23:26];        % Turn ON

% --- 6ms
idx = 37;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 6;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:16,19:21,23:26];        % Turn ON

% --- 7ms
idx = 38;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 7;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:16,19:21,23:26];

% --- 8ms
idx = 39;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 8;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:17,19:21,23:27];

% --- 9ms
idx = 40;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 9;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:17,19:21,23:27];

% --- 10ms
idx = 41;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 10;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [9:17,19:21,23:27];

% --- 11ms
idx = 42;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 11;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:21,23:26];

% --- 12ms
idx = 43;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 12;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:26];

% --- 13ms
idx = 44;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 13;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:21,23:26];

% --- 14ms
idx = 45;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           %  µA
overrides{idx}.target_ptd_ms  = 14;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:21,23:26];    % Turn ON

% --- 15ms
idx = 46;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           % 4 µA
overrides{idx}.target_ptd_ms  = 15;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:21,23:26];

% --- 17ms
idx = 47;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           %  µA
overrides{idx}.target_ptd_ms  = 17;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:21,23:26];    % Turn ON

% --- 20ms
idx = 48;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 5;           %  µA
overrides{idx}.target_ptd_ms  = 20;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:21,23:26];    % Turn ON

%% ================== Condition 4 =================
% --- Oms
idx = 49;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;          % 10 µA
overrides{idx}.target_ptd_ms  = 0;           % 0 ms (Simultaneous)
overrides{idx}.force_silent   = [1:64]; % Channels to turn OFF (noise)
overrides{idx}.force_respond  = [11,13,14,15,19:22,24,25];          % Channels to turn ON (rescue)

% --- 3ms
idx = 50;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 3;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [10:15,19:21,23:26];        % Turn ON

% ---4ms 
idx = 51;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 4;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [7,8,11:15,18:26];        % Turn ON

% --- 5ms
idx = 52;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 5;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:15,19:26];        % Turn ON

% --- 6ms
idx = 53;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 6;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:15,18:26];        % Turn ON

% --- 7ms
idx = 54;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 7;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,18:26];

% --- 8ms
idx = 55;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 8;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:15,19:26];

% --- 9ms
idx = 56;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 9;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:15,19:26];

% --- 10ms
idx = 57;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 10;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,18:26];

% --- 11ms
idx = 58;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 11;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,18:26];

% --- 12ms
idx = 59;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 12;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,18:26];

% --- 13ms
idx = 60;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 13;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,18:26];

% --- 14ms
idx = 61;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           %  µA
overrides{idx}.target_ptd_ms  = 14;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,18:26];    % Turn ON

% --- 15ms
idx = 62;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           % 4 µA
overrides{idx}.target_ptd_ms  = 15;           % 5 ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,19:26];

% --- 17ms
idx = 63;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           %  µA
overrides{idx}.target_ptd_ms  = 17;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,19:26];    % Turn ON

% --- 20ms
idx = 64;
overrides{idx}.target_set_idx = 2;           % Set 1
overrides{idx}.target_amp_uA  = 10;           %  µA
overrides{idx}.target_ptd_ms  = 20;           %  ms (Sequential)
overrides{idx}.force_silent   = [1:64];        % Turn OFF
overrides{idx}.force_respond  = [8,10:16,18:26];    % Turn ON

%% ========================================================================
%  3. INITIALIZATION & BACKUP
% ========================================================================
if ~isfolder(data_folder)
    error('Data folder does not exist: %s', data_folder);
end
cd(data_folder);

% Find the RespondingChannels file
file_list = dir('*_MultiISIRespondingChannels.mat');
if isempty(file_list)
    error('Could not find a *_MultiISIRespondingChannels.mat file in this folder.');
end

target_file = file_list(1).name;
full_path = fullfile(data_folder, target_file);

fprintf('\nLoading File: %s\n', target_file);
load(full_path, 'Responding'); % Load the main structure
% Also load the other variables to save them back safely
S_all = load(full_path); 

% Create a Backup (Only if one doesn't exist yet to prevent overwriting the original backup)
backup_name = strrep(target_file, '.mat', '_BACKUP.mat');
if ~isfile(backup_name)
    copyfile(target_file, backup_name);
    fprintf('Created Safe Backup: %s\n', backup_name);
else
    fprintf('Backup already exists. Safe to proceed.\n');
end

fprintf('\nStarting Batch Overrides...\n');
fprintf('--------------------------------------------------\n');

%% ========================================================================
%  4. EXECUTE OVERRIDES
% ========================================================================
total_changed = 0;

for k = 1:length(overrides)
    ov = overrides{k};
    
    si = ov.target_set_idx;
    target_amp = ov.target_amp_uA;
    target_ptd = ov.target_ptd_ms;
    
    % Failsafe: Check if the Set exists
    if si > numel(Responding.set)
        fprintf('WARNING [Override %d]: Set %d does not exist. Skipping.\n', k, si);
        continue;
    end
    
    % Search for the correct Amplitude index
    ai_found = 0;
    for ai = 1:numel(Responding.set(si).amp)
        % Using tolerance for floating point safety
        if abs(Responding.set(si).amp(ai).amp_value - target_amp) < 1e-4
            ai_found = ai;
            break;
        end
    end
    
    if ai_found == 0
        fprintf('WARNING [Override %d]: Amp %.1f uA not found in Set %d. Skipping.\n', k, target_amp, si);
        continue;
    end
    
    % Search for the correct PTD index
    pi_found = 0;
    for pi = 1:numel(Responding.set(si).amp(ai_found).ptd)
        if abs(Responding.set(si).amp(ai_found).ptd(pi).PTD_ms - target_ptd) < 1e-4
            pi_found = pi;
            break;
        end
    end
    
    if pi_found == 0
        fprintf('WARNING [Override %d]: PTD %.1f ms not found. Skipping.\n', k, target_ptd);
        continue;
    end
    
    % --- APPLY THE FIXES ---
    ai = ai_found;
    pi = pi_found;
    
    fprintf('Override %d -> Set %d | %.1f uA | %.1f ms:\n', k, si, target_amp, target_ptd);
    
    % Force Silent (Turn OFF)
    for ch = ov.force_silent
        if ch <= numel(Responding.set(si).amp(ai).ptd(pi).channel)
            Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = false;
            fprintf('   [-] Ch %02d forced SILENT\n', ch);
            total_changed = total_changed + 1;
        else
            fprintf('   [!] Ch %02d does not exist in data.\n', ch);
        end
    end
    
    % Force Respond (Turn ON)
    for ch = ov.force_respond
        if ch <= numel(Responding.set(si).amp(ai).ptd(pi).channel)
            Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = true;
            fprintf('   [+] Ch %02d forced RESPONDING\n', ch);
            total_changed = total_changed + 1;
        else
            fprintf('   [!] Ch %02d does not exist in data.\n', ch);
        end
    end
    
    if isempty(ov.force_silent) && isempty(ov.force_respond)
        fprintf('   (No channels modified for this condition)\n');
    end
end

fprintf('--------------------------------------------------\n');

%% ========================================================================
%  5. SAVE FINAL DATA
% ========================================================================
if total_changed > 0
    % Overwrite the 'Responding' struct in the loaded data
    S_all.Responding = Responding;
    
    % Save back to the original file
    save(full_path, '-struct', 'S_all');
    fprintf('COMPLETE: %d channel states updated and saved to main file.\n\n', total_changed);
else
    fprintf('No changes were made. File not overwritten.\n\n');
end