%% ============================================
%   MANUAL BAD TRIAL ENTRY (Channel-Specific) - FIXED
%   - Purpose: Flag trials as "Bad" for SPECIFIC CHANNELS.
%   - Input Format: Cell Array {Set, Amp, Rel_Idx, [Channels]}
%   - Logic: Maps (Set + Amp) -> Absolute Trial ID -> Applies to Listed Channels
%   - Safety: Appends to existing file.
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions/Simple_Analysis/MASSIVE'));

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX011/Xia_Exp1_Sim9';
Electrode_Type = 1;       

% --- MANUAL BAD TRIALS LIST (CELL ARRAY) ---
% Format: { Stim_Set_ID,  Amplitude_uA,  Relative_Trial_Index,  [List_of_Bad_Channels] }
% Note: Use curly braces {} and brackets [] for the channel list.

manual_additions = {
    % Example: Set 2, 6uA, 6th trial is bad on Ch 25 and 26
    % {2, 6, 6, [25, 26]};
    
    % Example: Set 2, 6uA, 8th trial is bad on ALL channels (1 to 32)
    % {2, 6, 8, [1:32]};
    
    % --- ENTER DATA HERE ---
    % {2, 6, 6,  [25 26]};      % Specific channels
    % {2, 6, 8,  [1:32]};       % All channels
    % {2, 6, 10, [19]};         % Single channel
    % {2, 6, 24, [20 21 22]};   % Multiple specific channels
    
    {1,8,2,[21,22,26,27,28,29,31]};
    {1,8,4,[21,24,26,27,31,32]};
    {1,8,5,[22,24,25,26,27]};
    {1,8,10,[22,26,25]};
    {1,8,15,[21,22,23,24,25,26,27,28,29,30,31,32]};
    {1,8,18,[26,28,29]};
    {1,8,19,[21,22]};
    {1,8,20,[21,22,23,26]};
    {1,8,21,[21,22,23,24,25,27,28,29]};
    {1,8,23,[21,25]};
    {1,8,24,[21,22,27,28,29]};
    {1,8,27,[26,27,28,29]};
    {1,8,30,[21,22,23,24,25,26,27,28,29,30,31,32]};
    {1,8,3,[22,26,27,28]};
    {1,8,11,[26,27,28,29,30,31,32]};
    {1,8,6,[22,26,27,28,29,30,31]};
    {1,8,16,[22,26,27,28,29,30,31]};
    {1,8,17,[26,27,28,29]};
    {1,8,1,[26,27,29]};
    
   

};

% ---------------- SETUP & CHECKS ----------------
if ~isfolder(data_folder), error('Folder not found'); end
cd(data_folder);
fprintf('\nRunning Channel-Specific Manual Bad Trial Entry in:\n%s\n\n', data_folder);

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
    fprintf('\nProcessing %d manual entries...\n', size(manual_additions, 1));
    count_added = 0;
    
    for i = 1:size(manual_additions, 1)
        % --- FIX: Unpack the inner cell array correctly ---
        current_row = manual_additions{i}; 
        
        tgt_set  = current_row{1}; 
        tgt_amp  = current_row{2}; 
        tgt_idx  = current_row{3}; 
        tgt_chs  = current_row{4}; % This is a vector of channels
        
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
        
        % --- APPLY TO SPECIFIED CHANNELS ---
        is_new_entry = false;
        
        for k = 1:length(tgt_chs)
            ch = tgt_chs(k);
            
            % Safety check for channel range
            if ch < 1 || ch > nCh_Total
                fprintf('    Warning: Channel %d out of range (1-%d). Skipped.\n', ch, nCh_Total);
                continue;
            end
            
            current_list = BadTrials{ch};
            new_list = unique([current_list(:); abs_trial(:)]);
            
            if numel(new_list) > numel(current_list)
                BadTrials{ch} = new_list;
                is_new_entry = true;
            end
        end
        
        if is_new_entry
            fprintf('  -> Added Bad Trial: Set %d | %.1f uA | Trial #%d (Abs: %d) -> applied to %d channels\n', ...
                tgt_set, tgt_amp, tgt_idx, abs_trial, length(tgt_chs));
            count_added = count_added + 1;
        end
    end
    fprintf('\nDone. Processed entries.\n');
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