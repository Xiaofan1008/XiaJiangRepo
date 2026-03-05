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
data_folder = '/Volumes/MACData/Data/Data_Xia/DX016/Xia_Exp1_Seq_Full_2';
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
    

   1,10,0,1;
   1,10,0,2;
   1,10,0,5;
   1,10,0,9;
   1,10,0,10;
   1,10,0,11;
   1,10,0,14;
   1,10,0,25;
   1,10,0,27;

   1,10,5,11;
   1,10,5,18;

   1,10,8,6;
   1,10,8,7;
   1,10,8,15;
   1,10,8,16;
   1,10,8,23;
   1,10,8,26;

   % 1,10,10,2;
   % 1,10,10,11;
   % 1,10,10,12;
   % 1,10,10,30;

   1,10,12,13;
   1,10,12,21;

   1,10,15,1;
   1,10,15,2;
   1,10,15,6;
   1,10,15,10;
   1,10,15,11;
   1,10,15,12;
   1,10,15,16;
   1,10,15,17;
   1,10,15,19;
   1,10,15,20;
   1,10,15,21;
   1,10,15,22;
   1,10,15,23;
   1,10,15,24;
   1,10,15,25;

   1,10,17,8;
   1,10,17,9;
   1,10,17,10;
   1,10,17,11;
   1,10,17,13;
   1,10,17,14;
   1,10,17,15;
   1,10,17,19;
   1,10,17,20;
   1,10,17,23;
   1,10,17,1;
   1,10,17,2;

   1,10,20,1;
   1,10,20,2;
   1,10,20,5;
   1,10,20,6;
   1,10,20,7;
   1,10,20,8;
   1,10,20,10;
   1,10,20,11;
   1,10,20,14;
   1,10,20,15;
   1,10,20,19;
   1,10,20,20;
   1,10,20,25;
   1,10,20,26;
   1,10,20,27;
   1,10,20,28;
   1,10,20,30;

   1,10,25,1;
   1,10,25,2;
   1,10,25,3;
   1,10,25,5;
   1,10,25,10;
   1,10,25,12;
   1,10,25,13;
   1,10,25,14;
   1,10,25,15;
   1,10,25,18;
   1,10,25,20;
   1,10,25,23;
   1,10,25,30;


   2,10,0,3;
   2,10,0,5;
   2,10,0,7;
   2,10,0,8;
   2,10,0,11;
   2,10,0,12;
   2,10,0,14;
   2,10,0,16;
   2,10,0,23;

   2,10,5,22;

   2,10,8,11;
   2,10,8,23;

   2,10,10,9;
   2,10,10,18;
   2,10,10,21;

   2,10,12,19;
   2,10,12,20;
   2,10,12,30;

   2,10,15,2;
   2,10,15,4;
   2,10,15,15;
   2,10,15,17;
   2,10,15,19;
   2,10,15,28;

   2,10,17,1;
   2,10,17,2;
   2,10,17,6;
   2,10,17,8;
   2,10,17,9;
   2,10,17,11;
   2,10,17,13;
   2,10,17,14;
   2,10,17,15;
   2,10,17,17;
   2,10,17,19;
   2,10,17,22;
   2,10,17,26;

   2,10,20,3;
   2,10,20,4;
   2,10,20,6;
   2,10,20,9;
   2,10,20,8;
   2,10,20,14;
   2,10,20,16;
   2,10,20,19;
   2,10,20,20;
   2,10,20,22;
   2,10,20,27;
   2,10,20,28;
   2,10,20,30;

   2,10,25,2;
   2,10,25,3;
   2,10,25,4;
   2,10,25,6;
   2,10,25,8;
   2,10,25,13;
   2,10,25,16;
   2,10,25,17;
   2,10,25,18;
   2,10,25,19;
   2,10,25,23;
   2,10,25,28;
   2,10,25,29;
   2,10,25,30;

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
bad_filename = sprintf('%s.MultiISIs_BadTrials.mat', base_name);
save_file_path = fullfile(data_folder, bad_filename);

if isfile(save_file_path)
    fprintf('Found existing file: %s\nLoading previous bad trials...\n', bad_filename);
    tmp = load(save_file_path);
    if isfield(tmp, 'MultiISIs_BadTrials'), BadTrials = tmp.BadTrials;
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