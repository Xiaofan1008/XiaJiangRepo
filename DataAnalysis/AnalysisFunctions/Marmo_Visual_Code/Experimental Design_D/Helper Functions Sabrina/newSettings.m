%% Preset some xml settings for quick testing
function settings = newSettings
settings{1} = 'A-001'; % Channel name
settings{2} = 1; % Stimulation enabled
settings{3} = 0; % Trigger Source
settings{4} = 0; % Trigger Edge or Level
settings{5} = 0; % Trigger High or Low
settings{6} = 0; % Post-trigger delay
settings{7} = 0; % Pulses or Train
settings{8} = 2; % Number of pulses
settings{9} = 10000; % Pulse Train Period
settings{10} = 1000; % Refractory Period (us)
settings{11} = 1; % Stim Shape
settings{12} = 0; % Stim Polarity
settings{13} = 100; % First Phase Duration (us)
settings{14} = 100; % Second Phase Duration (us)
settings{15} = 100; % Inter-phase Delay (us)
settings{16} = 20; % First phase Amplitude (uA)
settings{17} = 20; % Second phase Amplitude (uA)
settings{18} = 1; % Enable Amp Settle
settings{19} = 200; % Pre-Stim Amp Settle (us)
settings{20} = 160000; % Post-Stim Amp Settle (us) %was 160000
settings{21} = 1; % Maintain Amp Settle during inter-phase delay
settings{22} = 1; % Enable Charge Recovery
settings{23} = 100000; % Post Stim Charge Recovery On (us)
settings{24} = 150000; % Post Stim Charge Recovery Off (us)
end