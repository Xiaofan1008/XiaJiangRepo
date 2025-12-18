%% ============================================
%   MANUAL BAD TRIAL ENTRY (GLOBAL REMOVAL)
%   - Purpose: Manually flag trials as "Bad" for ALL CHANNELS.
%   - Input Format: [Stim_Set_ID, Amplitude, Relative_Index]
%   - Logic: Maps (Set + Amp) -> Absolute Trial ID -> Applies to ALL channels
%   - Safety: Appends to existing file.
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq3';
Electrode_Type = 1;       

% --- MANUAL BAD TRIALS LIST ---
% NEW Format: [ Stim_Set_ID,  Amplitude_uA,  Relative_Trial_Index ]
% NOTE: Channel is no longer needed because the bad trial applies to ALL channels.
manual_additions = [
    % Example: Set 2, 6uA, 6th trial shown
    % 2, 6, 6;
    
    % ENTER DATA HERE:
    2, 6, 6;
    2, 6, 8;
    2, 6, 10;
    2, 6, 24;
    
    % ... Add your other entries below ...
];

% ---------------- SETUP & CHECKS ----------------
if ~isfolder(data_folder), error('Folder not found'); end
cd(data_folder);
fprintf('\nRunning Global Manual Bad Trial Entry in:\n%s\n\n', data_folder);

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
    fprintf('\nProcessing %d manual entries (Applying to ALL channels)...\n', size(manual_additions, 1));
    count_added = 0;
    
    for i = 1:size(manual_additions, 1)
        tgt_set  = manual_additions(i, 1); % Set ID
        tgt_amp  = manual_additions(i, 2); % Amplitude
        tgt_idx  = manual_additions(i, 3); % Relative Index
        
        % Find Absolute Trials matching BOTH Set ID AND Amplitude
        matching_trials = find(combClass == tgt_set & abs(trialAmps - tgt_amp) < 0.01);
        
        if isempty(matching_trials)
            fprintf('  [SKIP] Set %d, Amp %.1f not found.\n', tgt_set, tgt_amp);
            continue;
        end
        
        if tgt_idx > numel(matching_trials)
            fprintf('  [SKIP] Set %d, Amp %.1f: Req #%d, but only %d trials exist.\n', ...
                tgt_set, tgt_amp, tgt_idx, numel(matching_trials));
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
            fprintf('  -> Added GLOBAL Bad Trial: Set %d | %.1f uA | Trial #%d (Abs: %d)\n', ...
                tgt_set, tgt_amp, tgt_idx, abs_trial);
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

fprintf('\nSaving to: %s ...\n', save_file);
if isfile(save_file)
    save(save_file, 'BadTrials', 'GoodTrials', '-append');
else
    detect_win_ms=[0 0]; slide_win_size_ms=0; max_spikes_in_slide=0;
    save(save_file, 'BadTrials', 'GoodTrials', 'detect_win_ms', ...
        'slide_win_size_ms', 'max_spikes_in_slide');
end
fprintf('Save Complete.\n');