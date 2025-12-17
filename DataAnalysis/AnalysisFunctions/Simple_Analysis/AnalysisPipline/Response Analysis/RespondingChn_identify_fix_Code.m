% 1. Define the path to your file (Copy-paste the full path from your previous code)
file_path = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq2_5ms/Xia_Exp1_Seq2_5ms_RespondingChannels.mat';

% 2. Save the modified 'Responding' structure back to the file
save(file_path,'Responding', 'Detection_Mode','baseline_win_ms','post_win_ms','post_win_singlepulse_ms', ...
    'k_SD','min_FR_post','min_frac_trials','Amps','PTDs');
fprintf('Successfully saved manual changes to:\n%s\n', file_path);