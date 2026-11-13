%% ============================================
%   MANUAL BAD TRIAL ENTRY
%   - Purpose: Manually flag specific trials as "Bad" based on visual inspection.
%   - Input: [Channel, Amplitude, Relative_Index] (e.g., "Ch 22, 10uA, 5th trial")
%   - Logic: Converts input -> Absolute Trial ID -> Updates .BadTrials.mat
%   - Safety: Loads existing file first (appends, does not overwrite previous work).
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Seq3';
Electrode_Type = 1;       

% --- MANUAL BAD TRIALS LIST ---
% Format: [ Channel_Number,  Amplitude_uA,  Relative_Trial_Index ]
% Relative_Index = The n-th trial for that amplitude as seen in the raster/plot.
manual_additions = [
    % Example: [22, 10, 5]; -> Channel 22, 10uA, the 5th trial displayed
    
    % ENTER DATA HERE:
    

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

% ---------------- LOAD METADATA ----------------
% 1. Map (Depth -> Channel)
d = Depth_s(Electrode_Type); 
nCh_Total = length(d);

% 2. StimParams (To map Amplitude -> Absolute Trial ID)
f_exp = dir('*_exp_datafile_*.mat');
if isempty(f_exp), error('No *_exp_datafile_*.mat found. Cannot map trials.'); end
S_exp = load(f_exp(1).name);
StimParams = S_exp.StimParams;
sim_stim   = S_exp.simultaneous_stim;

% Extract Amps for every absolute trial
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps     = trialAmps_all(1:sim_stim:end); % Downsample to 1 per trial
nTrials       = numel(trialAmps);

% ---------------- LOAD OR CREATE BAD TRIALS FILE ----------------
save_file = sprintf('%s.BadTrials.mat', base_name);

if isfile(save_file)
    fprintf('Found existing file: %s\nLoading previous bad trials...\n', save_file);
    tmp = load(save_file);
    if isfield(tmp, 'BadTrials')
        BadTrials = tmp.BadTrials;
    else
        fprintf('Warning: File exists but variable missing. Creating new.\n');
        BadTrials = cell(nCh_Total, 1);
    end
else
    fprintf('No existing file found. Creating new BadTrials file.\n');
    BadTrials = cell(nCh_Total, 1);
end

% Ensure cell array size is correct
if numel(BadTrials) < nCh_Total
    BadTrials{nCh_Total} = []; 
end

%% ---------------- PROCESS MANUAL ADDITIONS ----------------
if isempty(manual_additions)
    fprintf('\nNo manual additions entered in list. Nothing to do.\n');
else
    fprintf('\nProcessing %d manual entries...\n', size(manual_additions, 1));
    count_added = 0;
    
    for i = 1:size(manual_additions, 1)
        tgt_ch_depth = manual_additions(i, 1); % This is likely the Depth Index (Plot Ch)
        tgt_amp      = manual_additions(i, 2);
        tgt_idx      = manual_additions(i, 3);
        
        % Safety Check
        if tgt_ch_depth > nCh_Total || tgt_ch_depth < 1
            fprintf('  [SKIP] Invalid Channel %d (Max %d)\n', tgt_ch_depth, nCh_Total);
            continue;
        end
        
        % Find all Absolute Trial IDs matching this Amplitude
        matching_trials = find(abs(trialAmps - tgt_amp) < 0.01);
        
        if isempty(matching_trials)
            fprintf('  [SKIP] Amp %.1f uA not found in experiment.\n', tgt_amp);
            continue;
        end
        
        if tgt_idx > numel(matching_trials)
            fprintf('  [SKIP] Ch %d Amp %.1f: Requesting trial #%d, but only %d exist.\n', ...
                tgt_ch_depth, tgt_amp, tgt_idx, numel(matching_trials));
            continue;
        end
        
        % Convert to Absolute Trial Number
        abs_trial = matching_trials(tgt_idx);
        
        % --- FIX: Force column vector using (:) to prevent vertcat error ---
        prev_count = numel(BadTrials{tgt_ch_depth});
        
        current_list = BadTrials{tgt_ch_depth};
        % Combine existing list (forced column) with new trial (forced column)
        BadTrials{tgt_ch_depth} = unique([current_list(:); abs_trial(:)]);
        
        new_count  = numel(BadTrials{tgt_ch_depth});
        
        if new_count > prev_count
            fprintf('  -> Ch %d | %.1f uA | Trial #%d (Abs ID: %d) MARKED BAD.\n', ...
                tgt_ch_depth, tgt_amp, tgt_idx, abs_trial);
            count_added = count_added + 1;
        else
            fprintf('  -> Ch %d | %.1f uA | Trial #%d (Abs ID: %d) ALREADY BAD.\n', ...
                tgt_ch_depth, tgt_amp, tgt_idx, abs_trial);
        end
    end
    fprintf('\nDone. Added %d new bad trials.\n', count_added);
end

%% ---------------- UPDATE GOOD TRIALS & SAVE ----------------
GoodTrials = cell(nCh_Total, 1);
for ich = 1:nCh_Total
    GoodTrials{ich} = setdiff(1:nTrials, BadTrials{ich});
end

fprintf('\nSaving updated file: %s ...\n', save_file);
% Preserve existing settings if they exist, otherwise save defaults
if isfile(save_file)
    save(save_file, 'BadTrials', 'GoodTrials', '-append');
else
    % New file needs basic params to avoid loading errors in other scripts
    detect_win_ms = [0 0]; slide_win_size_ms = 0; max_spikes_in_slide = 0;
    save(save_file, 'BadTrials', 'GoodTrials', 'detect_win_ms', ...
        'slide_win_size_ms', 'max_spikes_in_slide');
end
fprintf('Save Complete.\n');