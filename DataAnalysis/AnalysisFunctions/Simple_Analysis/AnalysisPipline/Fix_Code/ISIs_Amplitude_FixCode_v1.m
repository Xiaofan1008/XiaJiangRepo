%% ========================================================================
%  Data Patch: Swap Reversed Amplitude Data in Saved Results
%  This script swaps the data slices for two specific amplitudes
%  without needing to re-run the spike counting logic.
% ========================================================================
clear;

% --- 1. USER SETTINGS ---
data_folder  = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_ISI_SimSeq2';

% Paste the exact path to your saved results file here:
results_path = '/Volumes/MACData/Data/Data_Xia/Analyzed_Results/Multi_ISIs_SpikeCount/DX018/Result_SpikeCount_FixWin_5_10uA_Xia_ISI_SimSeq2.mat';

amp_A = 5;  % The first reversed amplitude
amp_B = 10; % The second reversed amplitude

% --- 2. LOAD HARDWARE AMPS TO FIND INDICES ---
% We load the exp_datafile just to figure out which matrix index corresponds to 5 and 10
cd(data_folder);
fileDIR = dir('*_exp_datafile_*.mat');
if isempty(fileDIR)
    error('Could not find the exp_datafile in the data folder.');
end
S_exp = load(fileDIR(1).name, 'StimParams', 'simultaneous_stim');
simN  = S_exp.simultaneous_stim;
amps_all  = cell2mat(S_exp.StimParams(2:end,16)); 
trialAmps = amps_all(1:simN:end);
[Amps, ~, ~] = unique(trialAmps); 
Amps(Amps == -1) = 0;

% Find the exact dimension indices for 5 and 10
idx_A = find(abs(Amps - amp_A) < 0.001);
idx_B = find(abs(Amps - amp_B) < 0.001);

if isempty(idx_A) || isempty(idx_B)
    error('Could not find the specified amplitudes in the hardware file.');
end

% --- 3. LOAD AND SWAP THE DATA ---
fprintf('Loading results file...\n');
load(results_path, 'ResultFR');

fprintf('Swapping data between %.1f uA (Index %d) and %.1f uA (Index %d)...\n', amp_A, idx_A, amp_B, idx_B);

% 1. Swap SpikeCounts.Sim [Channels x Amps x Sets]
temp = ResultFR.SpikeCounts.Sim(:, idx_A, :);
ResultFR.SpikeCounts.Sim(:, idx_A, :) = ResultFR.SpikeCounts.Sim(:, idx_B, :);
ResultFR.SpikeCounts.Sim(:, idx_B, :) = temp;

% 2. Swap SpikeCounts.Seq [Channels x Amps x Sets x PTDs]
temp = ResultFR.SpikeCounts.Seq(:, idx_A, :, :);
ResultFR.SpikeCounts.Seq(:, idx_A, :, :) = ResultFR.SpikeCounts.Seq(:, idx_B, :, :);
ResultFR.SpikeCounts.Seq(:, idx_B, :, :) = temp;

% 3. Swap PSTH.Sim [Channels x Bins x Amps x Sets]
temp = ResultFR.PSTH.Sim(:, :, idx_A, :);
ResultFR.PSTH.Sim(:, :, idx_A, :) = ResultFR.PSTH.Sim(:, :, idx_B, :);
ResultFR.PSTH.Sim(:, :, idx_B, :) = temp;

% 4. Swap PSTH.Seq [Channels x Bins x Amps x Sets x PTDs]
temp = ResultFR.PSTH.Seq(:, :, idx_A, :, :);
ResultFR.PSTH.Seq(:, :, idx_A, :, :) = ResultFR.PSTH.Seq(:, :, idx_B, :, :);
ResultFR.PSTH.Seq(:, :, idx_B, :, :) = temp;

% --- 4. RE-SAVE THE FILE ---
% Create a backup just in case
backup_path = strrep(results_path, '.mat', '_PRE_SWAP_BACKUP.mat');
if ~isfile(backup_path)
    copyfile(results_path, backup_path);
    fprintf('Created a safety backup at:\n%s\n', backup_path);
end

save(results_path, 'ResultFR');
fprintf('\nSUCCESS! Amplitudes swapped and file updated.\n');