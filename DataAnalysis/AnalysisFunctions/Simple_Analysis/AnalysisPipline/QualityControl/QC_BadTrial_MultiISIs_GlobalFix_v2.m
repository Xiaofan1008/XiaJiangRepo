%% ============================================
%   MANUAL BAD TRIAL ENTRY (GLOBAL REMOVAL) - GROUPED VERSION
%   - Purpose: Manually flag trials as "Bad" for ALL CHANNELS.
%   - Input Format: Grouped by [Set, Amp], then lists {PTD, [Indices]}
%   - Logic: Maps (Set + Amp + PTD) -> Absolute Trial ID -> Applies to all channels
%   - Safety: Appends to existing file.
%   - FIXED: Uses fullfile() for robust saving
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX020/Xia_ISI_SimSeq2';
Electrode_Type = 2;       

% --- MANUAL BAD TRIALS LIST (GROUPED FORMAT) ---
Bad_Groups = [];

% ================ Condition 1: Set 1, 5 uA ================
c = 1;
Bad_Groups(c).Set = 1;
Bad_Groups(c).Amp = 5;
Bad_Groups(c).Trials = {
    % PTD (ms),  [List of Relative Bad Trial Indices]
      0, [];
      3, [6,7,8]; % Good
      4, [10,11,13]; % Bad
      5, [3,5,6,7]; % Good
      6, [12,17,18,20,21]; % Bad
      7, [8,9,11,14,15,16]; % Bad
      8, [];
      9, [1,6,8,13,15,17,18]; % Bad ,20,21,28
      10,[4,14,15,19,20,22,23,25]; % Bad ,26,27,28
      11,[]; % Good 1,2,11,14,15,21,22,23,25
      12,[]; 
      13,[];
      14,[1,3,5,7,8,9,10]; % Good ,11,14,16
      15,[1,3,4,6,7]; % Good ,8,10,12,14,17,19
      17,[1:7]; % Good ,10,17,21,24,26
      20,[1:3,5,6,8,10,12,13,14,16,17,19,21,22,23,28,30]; % Good
};

% ================ Condition 2: Set 1, 10 uA ================
c = 2;
Bad_Groups(c).Set = 1;
Bad_Groups(c).Amp = 10;
Bad_Groups(c).Trials = {
      0, [1,3,5,8,9]; % Bad 10,14,16
      3, [4,12,13,19,21,22]; % Bad ,15,16,18
      4, [1,17,18,24,25,27:30]; % Bad ,2,6,9
      5, [7,9,11:13,14,16:20,24,30]; % Bad
      6, [1,10,12,18,25,28]; % Bad ,19,23,29,30
      7, [3,12,13,18,27,29,30]; % Bad  ,7,9,23,24
      8, [6,10,14,18,19,29]; % Bad ,17,21,22,25,27
      9, [6,10,14,17:19,29]; % Bad 25,27,28,
      10,[1,3,6,8,12,17,28]; % Bad ,21,27
      11,[11,16,20,21,27,30]; % Bad 10,18,22,23
      12,[5,17]; % Bad ,16,19,26,18,22,27,29,30
      13,[1:6,8,16,20]; % Good ,11:14,21,22,23
      14,[2:5,8,10:14]; % Good ,16:18,20,24,
      15,[1,2,4,6,8,9,10,11,12,14,27]; % Good 17,19,21
      17,[1:3,5:7,9:11,13]; % Good ,16,18,22,26:27,
      20,[4,5,7,8,10:14,18:19]; % Good 21,22,28,29
};
% ================ Condition 3: Set 2, 5 uA ================
c = 3;
Bad_Groups(c).Set = 2;
Bad_Groups(c).Amp = 5;
Bad_Groups(c).Trials = {
      0, [];
      3, [];
      4, [];
      5, [8,13,16,17,21]; % Bad ,27,29 ,10,5,28
      6, [9,10,13,14,25,26,27,29]; % Bad 16,20,23,
      7, [10,16,18,25,]; % Bad 7,26,27,30
      8, [];
      9, [5,22,28]; % Bad 11,17,29
      10,[]; % Bad
      11,[];
      12,[];
      13,[]; % Good
      14,[]; % Good
      15,[]; % Good
      17,[]; % Good
      20,[]; % Good
};

% ================ Condition 4: Set 2, 10 uA ================
c = 4;
Bad_Groups(c).Set = 2;
Bad_Groups(c).Amp = 10;
Bad_Groups(c).Trials = {
      0, [];
      3, [];
      4, [];
      5, [];
      6, [];
      7, [];
      8, [];
      9, [];
      10,[];
      11,[];
      12,[];
      13,[];
      14,[];
      15,[];
      17,[];
      20,[];
};


%% ---------------- SETUP & CHECKS ----------------
if ~isfolder(data_folder), error('Folder not found'); end
cd(data_folder);
fprintf('\nRunning Global Manual Bad Trial Entry in:\n%s\n\n', data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u) > 4, base_name = last_folder(1 : u(end-1)-1);
else, base_name = last_folder; end

%% ---------------- LOAD METADATA & PARSE SETS ----------------
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

%% ---------------- LOAD OR CREATE BAD TRIALS FILE ----------------
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

%% ---------------- PROCESS MANUAL ADDITIONS (GROUPED LOGIC) ----------------
if isempty(Bad_Groups)
    fprintf('\nNo entries. Done.\n');
else
    fprintf('\nProcessing manual entries (Applying to ALL channels)...\n');
    count_added = 0;
    
    % Outer Loop: Iterate through each Condition Group
    for c = 1:length(Bad_Groups)
        tgt_set = Bad_Groups(c).Set;
        tgt_amp = Bad_Groups(c).Amp;
        
        % Middle Loop: Iterate through each PTD row in the Trials table
        for f = 1:size(Bad_Groups(c).Trials, 1)
            tgt_ptd = Bad_Groups(c).Trials{f, 1};
            bad_indices_list = Bad_Groups(c).Trials{f, 2};
            
            % Find Absolute Trials matching Set + Amp + PTD
            matching_trials = find(combClass == tgt_set & ...
                                   abs(trialAmps - tgt_amp) < 0.01 & ...
                                   abs(trialPTDs - tgt_ptd) < 0.01);
            
            if isempty(matching_trials)
                fprintf('  [SKIP] Set %d, Amp %.1f, PTD %.1f not found.\n', tgt_set, tgt_amp, tgt_ptd);
                continue;
            end
            
            % Inner Loop: Process each bad index in the array
            for b = 1:length(bad_indices_list)
                tgt_idx = bad_indices_list(b);
                
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
            end % End Inner Loop (Bad Indices)
        end % End Middle Loop (PTD rows)
    end % End Outer Loop (Condition Groups)
    
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