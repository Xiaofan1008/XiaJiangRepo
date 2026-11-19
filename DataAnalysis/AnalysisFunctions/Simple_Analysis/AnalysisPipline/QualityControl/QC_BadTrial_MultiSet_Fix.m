%% ============================================
%   MANUAL BAD TRIAL ENTRY (Sequential Support)
%   - Purpose: Manually flag trials as "Bad".
%   - Input Format: [Stim_Set_ID, Amplitude, Channel, Relative_Index]
%   - Logic: Maps (Set + Amp) -> Absolute Trial ID
%   - Safety: Appends to existing file.
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq3';
Electrode_Type = 1;       

% --- MANUAL BAD TRIALS LIST ---
% NEW Format: [ Stim_Set_ID,  Amplitude_uA,  Channel_Number,  Relative_Trial_Index ]
manual_additions = [
    % Example: Set 1, 5uA, Ch 26, 26th trial
    % 1, 5, 26, 26;
    
    % ENTER DATA HERE:
    1,5,21,29;1,5,21,20;
    1,5,26,29;1,5,26,20;
    1,5,27,29;1,5,27,20;
    1,5,28,29;1,5,28,20;
    1,5,29,20;

    1,4,21,25;1,4,21,24;1,4,21,19;1,4,21,10;1,4,21,16;
    1,4,19,21;
    1,4,20,21;
    1,4,24,8;
    1,4,25,24;
    1,4,26,24;1,4,26,21;
    1,4,27,24;
    1,4,28,24;
    1,4,29,20;

    1,3,8,10;
    1,3,19,4;
    1,3,20,4;1,3,20,19;
    1,3,21,3;1,3,21,24;
    1,3,25,17;
    1,3,26,17;1,3,26,19;1,3,26,3;1,3,26,5;
    1,3,27,3;1,3,27,5;
    1,3,28,3;1,3,28,5;
    1,3,29,3;

    1,6,21,9;1,6,21,12;
    1,6,24,9;1,6,25,9;1,6,26,9;1,6,27,9;1,6,28,9;1,6,29,9;

    1,2,20,16;
    1,2,21,26;1,2,21,8;
    1,2,24,26;
    1,2,26,26;1,2,26,16;1,2,26,7;1,2,26,13;1,2,26,11;
    1,2,27,11;1,2,27,26;1,2,27,16;1,2,27,8;
    1,2,28,11;1,2,28,26;
    1,2,29,11;1,2,29,26;
    1,2,30,11;1,2,30,26;1,2,30,16;
    1,2,31,26;
    1,2,32,26;

    1,1,10,18;
    1,1,19,28;1,1,19,29;1,1,19,26;
    1,1,20,28;1,1,20,29;1,1,20,23;1,1,20,19;1,1,20,18;1,1,20,15;
    1,1,21,28;1,1,21,29;1,1,21,23;1,1,21,19;1,1,21,18;1,1,21,15;1,1,21,16;1,1,21,17;
    1,1,22,28;1,1,22,18;1,1,22,23;
    1,1,23,21;
    1,1,24,7;1,1,24,10;1,1,24,21;1,1,24,1;1,1,24,3;1,1,24,13;1,1,24,16;1,1,24,29;
    1,1,25,13;1,1,25,16;1,1,25,15;1,1,25,29;1,1,25,28;1,1,25,21;1,1,25,1;1,1,25,2;1,1,25,3;1,1,25,10;1,1,25,7;
    1,1,26,13;1,1,26,16;1,1,26,15;1,1,26,29;1,1,26,28;1,1,26,21;1,1,26,1;1,1,26,3;1,1,26,10;1,1,26,7;1,1,26,6;
    1,1,27,1;1,1,27,3;1,1,27,6;1,1,27,7;1,1,27,10;1,1,27,13;1,1,27,15;1,1,27,21;1,1,27,28;1,1,27,29;
    1,1,28,3;1,1,28,6;1,1,28,10;1,1,28,15;1,1,28,16;1,1,28,28;1,1,28,29;
    1,1,29,14;1,1,29,21;1,1,29,29;1,1,29,28;
    1,1,30,10;1,1,30,28;


    % ... Add your other entries below ...
];

% ---------------- SETUP & CHECKS ----------------
if ~isfolder(data_folder), error('Folder not found'); end
cd(data_folder);
fprintf('\nRunning Manual Bad Trial Entry in:\n%s\n\n', data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u) >= 4, base_name = last_folder(1 : u(end-1)-1);
else, base_name = last_folder; end

% ---------------- LOAD METADATA & PARSE SETS ----------------
d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

f_exp = dir('*_exp_datafile_*.mat');
if isempty(f_exp), error('No exp_datafile found.'); end
S_exp = load(f_exp(1).name);

% 1. Extract Amplitudes
StimParams = S_exp.StimParams;
simN       = S_exp.simultaneous_stim;
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:simN:end); 
nTrials       = numel(trialAmps);

% 2. Extract Stimulation Sets (CombClass)
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
% Identify Set IDs (1, 2, 3...) matching the analysis order
[~,~,combClass] = unique(comb_seq,'rows','stable');

% ---------------- LOAD OR CREATE BAD TRIALS FILE ----------------
save_file = sprintf('%s.BadTrials.mat', base_name);

if isfile(save_file)
    fprintf('Found existing file: %s\nLoading previous bad trials...\n', save_file);
    tmp = load(save_file);
    if isfield(tmp, 'BadTrials'), BadTrials = tmp.BadTrials;
    else, BadTrials = cell(nCh_Total, 1); end
else
    fprintf('Creating new BadTrials file.\n');
    BadTrials = cell(nCh_Total, 1);
end
if numel(BadTrials) < nCh_Total, BadTrials{nCh_Total} = []; end

%% ---------------- PROCESS MANUAL ADDITIONS ----------------
if isempty(manual_additions)
    fprintf('\nNo entries. Done.\n');
else
    fprintf('\nProcessing %d manual entries...\n', size(manual_additions, 1));
    count_added = 0;
    
    for i = 1:size(manual_additions, 1)
        tgt_set  = manual_additions(i, 1); % Set ID
        tgt_amp  = manual_additions(i, 2); % Amplitude
        tgt_ch   = manual_additions(i, 3); % Channel
        tgt_idx  = manual_additions(i, 4); % Relative Index
        
        if tgt_ch > nCh_Total || tgt_ch < 1
            fprintf('  [SKIP] Invalid Ch %d\n', tgt_ch); continue; 
        end
        
        % Find Absolute Trials matching BOTH Set ID AND Amplitude
        matching_trials = find(combClass == tgt_set & abs(trialAmps - tgt_amp) < 0.01);
        
        if isempty(matching_trials)
            fprintf('  [SKIP] Set %d, Amp %.1f not found.\n', tgt_set, tgt_amp);
            continue;
        end
        
        if tgt_idx > numel(matching_trials)
            fprintf('  [SKIP] Set %d, Amp %.1f (Ch %d): Req #%d, but only %d trials exist.\n', ...
                tgt_set, tgt_amp, tgt_ch, tgt_idx, numel(matching_trials));
            continue;
        end
        
        % Convert to Absolute ID
        abs_trial = matching_trials(tgt_idx);
        
        % Add to list (Column Vector Force)
        current_list = BadTrials{tgt_ch};
        BadTrials{tgt_ch} = unique([current_list(:); abs_trial(:)]);
        
        if numel(BadTrials{tgt_ch}) > numel(current_list)
            fprintf('  -> Added: Set %d | %.1f uA | Ch %d | Trial #%d (Abs: %d)\n', ...
                tgt_set, tgt_amp, tgt_ch, tgt_idx, abs_trial);
            count_added = count_added + 1;
        end
    end
    fprintf('\nDone. Added %d new bad trials.\n', count_added);
end

%% ---------------- SAVE ----------------
GoodTrials = cell(nCh_Total, 1);
for ich = 1:nCh_Total
    GoodTrials{ich} = setdiff(1:nTrials, BadTrials{ich});
end

fprintf('\nSaving to: %s ...\n', save_file);
if isfile(save_file)
    save(save_file, 'BadTrials', 'GoodTrials', '-append');
else
    detect_win_ms=[0 0]; slide_win_size_ms=0; max_spikes_in_slide=0;
    save(save_file, 'BadTrials', 'GoodTrials', 'detect_win_ms', ...
        'slide_win_size_ms', 'max_spikes_in_slide');
end
fprintf('Save Complete.\n');