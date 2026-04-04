%% ============================================
%   MANUAL BAD TRIAL ENTRY (GLOBAL REMOVAL)
%   - Purpose: Manually flag trials as "Bad" for ALL CHANNELS.
%   - Input Format: [Stim_Set_ID, Amplitude, PTD, Relative_Index]
%   - Logic: Maps (Set + Amp + PTD) -> Absolute Trial ID -> Applies to ALL channels
%   - Safety: Appends to existing file.
%   - FIXED: Uses fullfile() for robust saving
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX020/Xia_ISI_SimSeq1';
Electrode_Type = 2;       

% --- MANUAL BAD TRIALS LIST ---
% Format: [ Stim_Set_ID,  Amplitude_uA,  PTD_ms,  Relative_Trial_Index ]
% NOTE: PTD = 0 for Simultaneous.
manual_additions = [
    % Example: Set 2, 6uA, 5ms PTD, 6th trial shown
    % 2, 6, 5, 6;
    
    % ENTER DATA HERE:
    % 1, 3, 0, 30;  % Set 1, 3uA, Sim(0ms), Trial 30
    % 1, 2, 5, 37;  % Set 1, 2uA, Seq(5ms), Trial 37
    


   1,10,3,1;
   1,10,3,2;

   1,10,5,1;
   1,10,5,2;
   1,10,5,3;
   1,10,5,4;
   1,10,5,14;
   1,10,5,24;

   1,10,6,2;
   1,10,6,3;
   1,10,6,6;
   1,10,6,19;
   1,10,6,24;

   1,10,7,1;
   1,10,7,2;
   1,10,7,16;
   1,10,7,18;
   1,10,7,20;
   1,10,7,22;
   1,10,7,24;
   1,10,7,28;
   1,10,7,29;
   1,10,7,30;

   1,10,8,1;
   1,10,8,2;
   1,10,8,3;
   1,10,8,7;
   1,10,8,12;
   1,10,8,13;
   1,10,8,18;
   1,10,8,19;
   1,10,8,23;
   1,10,8,24;
   1,10,8,25;
   1,10,8,27;
   1,10,8,29;

   1,10,9,1;
   1,10,9,13;
   1,10,9,14;
   1,10,9,16;
   1,10,9,18;
   1,10,9,20;
   1,10,9,21;
   1,10,9,23;
   1,10,9,24;
   1,10,9,25;
   1,10,9,30;

   1,10,10,1;
   1,10,10,3;
   1,10,10,4;
   1,10,10,5;
   1,10,10,6;
   1,10,10,7;
   1,10,10,8;

   1,10,10,15;
   % 1,10,10,19;
   % 1,10,10,23;

   1,10,11,1;
   1,10,11,2;
   1,10,11,3;
   1,10,11,5;
   1,10,11,19;
   1,10,11,23;
   1,10,11,25;
   1,10,11,27;
   1,10,11,28;

   1,10,11,6;
   1,10,11,21;
   1,10,11,28;
   1,10,11,30;

   1,10,12,1;
   1,10,12,2;
   1,10,12,4;
   1,10,12,5;
   1,10,12,10;
   1,10,12,20;
   1,10,12,21;
   1,10,12,25;

   1,10,12,14;
   1,10,12,19;
   1,10,12,23;
   1,10,12,26;
   1,10,12,30;

   1,10,13,1;
   1,10,13,2;
   1,10,13,3;
   1,10,13,4;
   1,10,13,6;
   1,10,13,25;
   1,10,13,26;
   1,10,13,28;

   1,10,13,13;
   1,10,13,14;
   1,10,13,16;
   1,10,13,19;
   1,10,13,11;
   1,10,13,10;

   1,10,14,2;
   1,10,14,3;
   1,10,14,5;
   1,10,14,7;
   1,10,14,30;

   1,10,14,20;
   1,10,14,14;
   1,10,14,16;
   1,10,14,21;
   1,10,14,25;

   1,10,15,1;
   1,10,15,2;
   1,10,15,6; 
   % 1,10,15,7;
   1,10,15,9;
   % 1,10,15,17;
   1,10,15,21;
   1,10,15,24; 
   % 1,10,15,25;
   1,10,15,26;
   % 1,10,15,28;

   1,10,15,20;
   1,10,15,21;
   1,10,15,26;
   1,10,15,30;
   
   % 1,10,17,1;
   % 1,10,17,2;
   1,10,17,13;
   % 1,10,17,17;
   1,10,17,18;
   1,10,17,19;
   1,10,17,21;
   % 1,10,17,24;

   1,10,17,7;
   1,10,17,8;
   1,10,17,9;
   1,10,17,5;
   1,10,17,22;
   1,10,17,24;

   % 1,10,20,8;
   % 1,10,20,14;
   % 1,10,20,17;
   % 1,10,20,20;

   1,10,20,6;
   1,10,20,13;
   1,10,20,19;
   1,10,20,21; 
   1,10,20,23;
   1,10,20,24;
   1,10,20,28;
   1,10,20,30;



   2,10,3,1;
   2,10,3,2;
   2,10,3,3;
   2,10,3,20;

   2,10,4,1;
   2,10,4,2;
   2,10,4,3;
   2,10,4,4;
   2,10,4,24;


   2,10,5,1;
   2,10,5,2;
   2,10,5,3;
   2,10,5,4;
   2,10,5,5;
   2,10,5,6;
   2,10,5,11;
   2,10,5,24;
   2,10,5,28;
   2,10,5,29;

   2,10,6,1;
   2,10,6,4;
   2,10,6,5;
   2,10,6,8;
   2,10,6,17;
   2,10,6,19;
   2,10,6,22;

   2,10,7,1;
   2,10,7,2;
   2,10,7,3;
   2,10,7,4;
   2,10,7,5;
   2,10,7,29;

   2,10,8,1;
   2,10,8,2;
   2,10,8,3;
   2,10,8,9;
   2,10,8,17;
   2,10,8,27;
   2,10,8,28;

   2,10,9,1;
   2,10,9,2;
   2,10,9,9;
   2,10,9,25;
   2,10,9,26;


   2,10,10,2;
   2,10,10,3;
   2,10,10,4;
   2,10,10,5;
   2,10,10,7;
   2,10,10,28;

   2,10,10,11;
   2,10,10,13;
   2,10,10,16;

   2,10,11,1;
   2,10,11,7;
   2,10,11,8;
   2,10,11,15;

   % 2,10,11,12;
   2,10,11,18;
   % 2,10,11,21;


   2,10,12,1;
   2,10,12,3;
   2,10,12,20;
   2,10,12,27;
   2,10,12,29;

   2,10,12,18;
   2,10,12,19;
   % 2,10,12,22;
   % 2,10,12,24;

   2,10,13,1;
   2,10,13,2;
   2,10,13,28;

   2,10,13,12;
   2,10,13,14;
   2,10,13,17;
   2,10,13,21;
   2,10,13,25;
   2,10,13,26;
   2,10,13,27;
   2,10,13,30;
   

   % 2,10,14,1;
   % 2,10,14,2;
   
   2,10,14,5;
   2,10,14,10;
   2,10,14,13;
   2,10,14,14;
   2,10,14,19;
   2,10,14,21;
   2,10,14,22;
   2,10,14,23;
   2,10,14,30;


   % 2,10,15,1;
   % 2,10,15,2;
   % 2,10,15,4;
   % 2,10,15,5;

   2,10,15,7;
   2,10,15,13;
   2,10,15,18;
   2,10,15,25;
   2,10,15,26;
   2,10,15,27;

   
    2,10,17,7;
    2,10,17,12;
    2,10,17,14;
    2,10,17,18;
    2,10,17,21;
    2,10,17,26;

    2,10,20,8;
    2,10,20,10;
    2,10,20,13;
    2,10,20,16;
    2,10,20,17;
    2,10,20,18;
    2,10,20,24;
    2,10,20,25;
    2,10,20,26;
    2,10,20,29;
    2,10,20,30;


];

% ---------------- SETUP & CHECKS ----------------
if ~isfolder(data_folder), error('Folder not found'); end
cd(data_folder);
fprintf('\nRunning Global Manual Bad Trial Entry in:\n%s\n\n', data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u) > 4, base_name = last_folder(1 : u(end-1)-1);
else, base_name = last_folder; end

% ---------------- LOAD METADATA & PARSE SETS ----------------
d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

f_exp = dir('*_exp_datafile_*.mat');
if isempty(f_exp), error('No exp_datafile found.'); end
S_exp = load(f_exp(1).name);

StimParams = S_exp.StimParams;
simN       = S_exp.simultaneous_stim;

% 1. Extract Amplitudes
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:simN:end); 
nTrials       = numel(trialAmps);

% 2. Extract PTDs (Pulse Train Delay)
if simN > 1
    PTD_all = cell2mat(StimParams(3:simN:end,6)); % µs from 2nd pulse row
else
    PTD_all = zeros(nTrials,1);
end
trialPTDs = PTD_all / 1000; % Convert to ms

% 3. Extract Stimulation Sets (CombClass)
E_MAP = S_exp.E_MAP;
stimNames = StimParams(2:end,1); 
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

% Reconstruct trial combinations
comb_seq = zeros(nTrials, simN);
for t = 1:nTrials
    rr = (t-1)*simN + (1:simN); 
    v = idx_all(rr); v = v(v>0); 
    comb_seq(t,1:numel(v)) = v(:).'; 
end
[~,~,combClass] = unique(comb_seq,'rows','stable');

% ---------------- LOAD OR CREATE BAD TRIALS FILE ----------------
% [FIX] Use fullfile to ensure correct path even if PWD changes
bad_filename = sprintf('%s.MultiISIsBadTrials.mat', base_name);
save_file_path = fullfile(data_folder, bad_filename);

if isfile(save_file_path)
    fprintf('Found existing file: %s\nLoading previous bad trials...\n', bad_filename);
    tmp = load(save_file_path);
    if isfield(tmp, 'MultiISIsBadTrials'), BadTrials = tmp.BadTrials;
    else, BadTrials = cell(nCh_Total, 1); end
else
    fprintf('Creating new BadTrials file: %s\n', bad_filename);
    BadTrials = cell(nCh_Total, 1);
end
if numel(BadTrials) < nCh_Total, BadTrials{nCh_Total} = []; end

%% ---------------- PROCESS MANUAL ADDITIONS ----------------
if isempty(manual_additions)
    fprintf('\nNo entries. Done.\n');
else
    fprintf('\nProcessing %d manual entries (Applying to ALL channels)...\n', size(manual_additions, 1));
    count_added = 0;
    
    for i = 1:size(manual_additions, 1)
        tgt_set  = manual_additions(i, 1); % Set ID
        tgt_amp  = manual_additions(i, 2); % Amplitude
        tgt_ptd  = manual_additions(i, 3); % PTD (ms)
        tgt_idx  = manual_additions(i, 4); % Relative Index
        
        % Find Absolute Trials matching Set + Amp + PTD
        matching_trials = find(combClass == tgt_set & ...
                               abs(trialAmps - tgt_amp) < 0.01 & ...
                               abs(trialPTDs - tgt_ptd) < 0.01);
        
        if isempty(matching_trials)
            fprintf('  [SKIP] Set %d, Amp %.1f, PTD %.1f not found.\n', tgt_set, tgt_amp, tgt_ptd);
            continue;
        end
        
        if tgt_idx > numel(matching_trials)
            fprintf('  [SKIP] Set %d, Amp %.1f, PTD %.1f: Req #%d, but only %d trials exist.\n', ...
                tgt_set, tgt_amp, tgt_ptd, tgt_idx, numel(matching_trials));
            continue;
        end
        
        % Convert to Absolute ID
        abs_trial = matching_trials(tgt_idx);
        
        % --- APPLY TO ALL CHANNELS ---
        is_new = false;
        for ch = 1:nCh_Total
            current_list = BadTrials{ch};
            new_list = unique([current_list(:); abs_trial(:)]);
            
            if numel(new_list) > numel(current_list)
                BadTrials{ch} = new_list;
                is_new = true;
            end
        end
        
        if is_new
            fprintf('  -> Added GLOBAL Bad Trial: Set %d | %.1f uA | PTD %.1f | Trial #%d (Abs: %d)\n', ...
                tgt_set, tgt_amp, tgt_ptd, tgt_idx, abs_trial);
            count_added = count_added + 1;
        end
    end
    fprintf('\nDone. Added %d new global bad trials.\n', count_added);
end

%% ---------------- SAVE ----------------
GoodTrials = cell(nCh_Total, 1);
for ich = 1:nCh_Total
    GoodTrials{ich} = setdiff(1:nTrials, BadTrials{ich});
end

fprintf('\nSaving to: %s ...\n', save_file_path);

% Save to full_file_path
if isfile(save_file_path)
    save(save_file_path, 'BadTrials', 'GoodTrials', '-append');
else
    detect_win_ms=[0 0]; slide_win_size_ms=0; max_spikes_in_slide=0;
    save(save_file_path, 'BadTrials', 'GoodTrials', 'detect_win_ms', ...
        'slide_win_size_ms', 'max_spikes_in_slide');
end
fprintf('Save Complete.\n');