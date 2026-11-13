%% ========================================================================
%  Batch Manual Override for Responding Channels (Compact Table Version)
%
%  One flat table, one row per (Set, Amp, PTD) you need to fix. Each row
%  lists ONLY the channels to force RESPONDING and/or force SILENT for
%  that condition -- any channel not listed is left exactly as it was
%  (whatever the automatic detection produced).
%
%  Same safety behavior as before: creates a backup before overwriting,
%  and only saves if something actually changed.
% ========================================================================
clear;

%% ========================================================================
%  1. FILE PATH SETTINGS
% ========================================================================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX023/Xia_ISI_SimSeq1';

%% ========================================================================
%  2. OVERRIDES
%     Only list rows where you actually need to change something.
%     ResetAllFirst: if true, ALL channels present for this condition are
%       force-silenced BEFORE ForceRespond/ForceSilent are applied -- use
%       this when a condition's auto-detection is bad enough that you
%       want to wipe it and hand-pick the real responders from scratch.
%     ForceRespond : channels forced to RESPONDING for this PTD.
%     ForceSilent  : channels forced to SILENT for this PTD (redundant
%       with ResetAllFirst=true, but useful on its own when
%       ResetAllFirst=false and you just want to silence a few).
%     Leave ForceRespond/ForceSilent empty ([]) if not needed for that
%     row. A channel should never appear in both lists for the same row
%     -- that's flagged as a conflict and skipped for safety.
% ========================================================================
Overrides = {
  % Set  Amp   PTD(ms) ResetAllFirst ForceRespond ForceSilent
    1,   5,    0,        true,           [3,4,6:11,13,15:16,33:36,39,41,44:46],     [];
    1,   5,    3,        true,           [1:11,13:15,32:37,39:42,44:50,54:55,57:58,61:62],     [];
    1,   5,    4,        true,           [1:11,13,15,33:48,50:52,54:55,57:58,61:63],     [];
    1,   5,    5,        true,           [1:15,33:50,52,54:55,57:58,61:63],     [];
    1,   5,    6,        true,           [1:15,33:52,54:55,57:58,61:63],     [];
    1,   5,    7,        true,           [1:14,15,33:52,53,57:58,61:63],     [];    
    1,   5,    8,        true,           [1:15,33:52,54:55,57:58,61:63],     [];
    1,   5,    9,        true,           [1:15,33:52,54:58,61:63],     [];
    1,   5,    10,       true,           [1:15,33:37,39:52,54:55,57:58,61:63],     [];
    1,   5,    11,       true,           [1:15,33:42,44:52,54:55,57:58,61:63],     [];
    1,   5,    12,       true,           [1:15,33:48,50:52,54:55,57:58,61:63],     [];
    1,   5,    13,       true,           [1:15,33:37,39:42,44:48,51,54:55.57:58,61:62],     [];
    1,   5,    14,       true,           [1:15,33:52,54:55,57:58,61:63],     [];
    1,   5,    15,       true,           [1:13,15,33:37,39:42,44:48,50:51,54:55,57:58,61:63],     [];
    1,   5,    17,       true,           [1:15,33:48,50:52,54:55,57:58,61:62],     [];
    1,   5,    20,       true,           [1:15,33:37,39:42,44:48,50:52,54:55,57:58,61:63],     [];

    1,   10,    0,        true,           [2:11,13,15,43,45],     [];
    1,   10,    3,        true,           [1:16,33:48,50:52,54:55,57,61:63],     [];
    1,   10,    4,        true,           [1:16,33:52,54:55,57:58,62:63],     [];
    1,   10,    5,        true,           [1:15,33:52,54:58,61:63],     [];
    1,   10,    6,        true,           [1:16,33:52,54:55,57:58,61:63],     [];
    1,   10,    7,        true,           [1:16,33:48,50:52,54:58,61:63],     [];    
    1,   10,    8,        true,           [1:16,33:52,54:55,57:58,61:63],     [];
    1,   10,    9,        true,           [1:16,33:52,54:58,61:62],     [];
    1,   10,    10,       true,           [1:16,33:52,54:58,61:63],     [];
    1,   10,    11,       true,           [1:16,33:51,54:55,57:58,61:63],     [];
    1,   10,    12,       true,           [1:16,33:51,54:55,57:58,61:63],     [];
    1,   10,    13,       true,           [1:16,33:51,54:55,57:58,61:63],     [];
    1,   10,    14,       true,           [1:16,33:51,54:55,57:58,61:63],     [];
    1,   10,    15,       true,           [1:16,33:51,54:55,57:58,61:63],     [];
    1,   10,    17,       true,           [1:16,33:51,54:55,57:58,61:63],     [];
    1,   10,    20,       true,           [1:16,33:51,54:55,57:58,61:63],     [];

    2,   5,    0,        true,           [2:11,13,15,33:34,44:45,47],     [];
    2,   5,    3,        true,           [1:15,33:37,39:42,44:52,54:55,57:58,61:62],     [];
    2,   5,    4,        true,           [1:15,33:42,44:52,54:55,57:58,61:62],     [];
    2,   5,    5,        true,           [1:15,33:48,50:52,54:55,57:57,61:62],     [];
    2,   5,    6,        true,           [1:15,33:52,54,57:58,61:62],     [];
    2,   5,    7,        true,           [1:15,33:51,54:55,57:58,61:62],     [];    
    2,   5,    8,        true,           [1:15,33:52,55,57:58,61:63],     [];
    2,   5,    9,        true,           [1:15,33:52,54:55,57:58,61:63],     [];
    2,   5,    10,       true,           [1:15,33:37,39:52,54:55,58,61:63],     [];
    2,   5,    11,       true,           [1:15,33:37,39:42,44:52,54:55,57:58,61:63],     [];
    2,   5,    12,       true,           [1:15,33:37,39:42,44:52,54:55,57:58,61:63],     [];
    2,   5,    13,       true,           [1:15,33:37,39:42,44:52,54:55,57:58,61:63],     [];
    2,   5,    14,       true,           [1:15,33:37,39:42,44:52,54,58,62:63],     [];
    2,   5,    15,       true,           [1:15,33:37,39:42,44:52,54:55,58,62],     [];
    2,   5,    17,       true,           [1:15,33:37,39:42,44:52,54,57:58,62:63],     [];
    2,   5,    20,       true,           [1:15,33:37,39:42,44:51,58,61:62],     [];

    2,   10,    0,        true,           [1:11,13,15,34:37,40,42:43,45:48],     [];
    2,   10,    3,        true,           [1:14,15,33:52,54:55,57:58,62:63],     [];
    2,   10,    4,        true,           [1:15,33:42,44:52,54:55,57:58,60:63],     [];
    2,   10,    5,        true,           [1:16,33:52,54:55,57:58,61:63],     [];
    2,   10,    6,        true,           [1:16,33:52,54:55,57:58,60:63],     [];
    2,   10,    7,        true,           [1:16,33:52,54:55,57:58,60:63],     [];    
    2,   10,    8,        true,           [1:16,33:52,54:55,57:58,60:63],     [];
    2,   10,    9,        true,           [1:16,33:52,54:55,57:58,60:63],     [];
    2,   10,    10,       true,           [1:16,33:52,54:55,57:58,61:63],     [];
    2,   10,    11,       true,           [1:16,33:52,54:55,57:58,61:63],     [];
    2,   10,    12,       true,           [1:16,33:52,54:55,57:58,60:63],     [];
    2,   10,    13,       true,           [1:16,33:48,51:52,54:55,57:58,61:63],     [];
    2,   10,    14,       true,           [1:16,33:48,51:52,54:55,58,61:63],     [];
    2,   10,    15,       true,           [1:16,33:48,51:52,54:55,57:58,61:63],     [];
    2,   10,    17,       true,           [1:16,33:48,51:52,54:55,57:58,61:63],     [];
    2,   10,    20,       true,           [1:16,33:48,51:52,54:55,57:58,61:63],     [];





};

%% ========================================================================
%  3. INITIALIZATION & BACKUP
% ========================================================================
if ~isfolder(data_folder)
    error('Data folder does not exist: %s', data_folder);
end
cd(data_folder);

file_list = dir('*_MultiISI_RespondingChannels.mat');
if isempty(file_list)
    error('Could not find a *_MultiISI_RespondingChannels.mat file in this folder.');
end
target_file = file_list(1).name;
full_path = fullfile(data_folder, target_file);

fprintf('\nLoading File: %s\n', target_file);
load(full_path, 'Responding');
S_all = load(full_path);

backup_name = strrep(target_file, '.mat', '_BACKUP.mat');
if ~isfile(backup_name)
    copyfile(target_file, backup_name);
    fprintf('Created Safe Backup: %s\n', backup_name);
else
    fprintf('Backup already exists. Safe to proceed.\n');
end

fprintf('\nStarting Batch Overrides...\n');
fprintf('--------------------------------------------------\n');

total_changed = 0;

%% ========================================================================
%  4. APPLY OVERRIDES
% ========================================================================
for row = 1:size(Overrides,1)

    si            = Overrides{row,1};
    target_amp    = Overrides{row,2};
    target_ptd    = Overrides{row,3};
    reset_all_first = Overrides{row,4};
    force_resp    = Overrides{row,5};
    force_sil     = Overrides{row,6};

    if si > numel(Responding.set)
        fprintf('WARNING [Row %d]: Set %d does not exist. Skipping.\n', row, si);
        continue;
    end

    ai = find_index_by_value( ...
        [Responding.set(si).amp.amp_value], target_amp);
    if isempty(ai)
        fprintf('WARNING [Row %d]: Amp %.1f uA not found in Set %d. Skipping.\n', ...
            row, target_amp, si);
        continue;
    end

    pi = find_index_by_value( ...
        [Responding.set(si).amp(ai).ptd.PTD_ms], target_ptd);
    if isempty(pi)
        fprintf('WARNING [Row %d]: PTD %.1f ms not found. Skipping.\n', row, target_ptd);
        continue;
    end

    % Conflict check: a channel cannot be forced both ways in one row.
    conflict = intersect(force_resp, force_sil);
    if ~isempty(conflict)
        fprintf(['WARNING [Row %d]: channel(s) %s listed in BOTH ' ...
            'ForceRespond and ForceSilent -- skipping those channels ' ...
            'for this row.\n'], row, num2str(conflict));
        force_resp = setdiff(force_resp, conflict);
        force_sil  = setdiff(force_sil, conflict);
    end

    fprintf('Row %d -> Set %d | %.1f uA | %.1f ms:\n', row, si, target_amp, target_ptd);

    if reset_all_first
        nChannelsThisCondition = numel(Responding.set(si).amp(ai).ptd(pi).channel);
        for ch = 1:nChannelsThisCondition
            Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = false;
            total_changed = total_changed + 1;
        end
        fprintf('   [reset] All %d channels forced SILENT first\n', nChannelsThisCondition);
    end

    for ch = force_sil
        if ch <= numel(Responding.set(si).amp(ai).ptd(pi).channel)
            Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = false;
            fprintf('   [-] Ch %02d forced SILENT\n', ch);
            total_changed = total_changed + 1;
        else
            fprintf('   [!] Ch %02d does not exist in data.\n', ch);
        end
    end

    for ch = force_resp
        if ch <= numel(Responding.set(si).amp(ai).ptd(pi).channel)
            Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = true;
            fprintf('   [+] Ch %02d forced RESPONDING\n', ch);
            total_changed = total_changed + 1;
        else
            fprintf('   [!] Ch %02d does not exist in data.\n', ch);
        end
    end

    if isempty(force_resp) && isempty(force_sil)
        fprintf('   (No channels modified for this row)\n');
    end
end

fprintf('--------------------------------------------------\n');

%% ========================================================================
%  5. SAVE FINAL DATA
% ========================================================================
if total_changed > 0
    S_all.Responding = Responding;
    save(full_path, '-struct', 'S_all');
    fprintf('COMPLETE: %d channel states updated and saved to main file.\n\n', total_changed);
else
    fprintf('No changes were made. File not overwritten.\n\n');
end

%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function idx = find_index_by_value(valueArray, targetValue)
% Find the index of the closest match within tolerance, or [] if none.
idx = find(abs(valueArray - targetValue) < 1e-4, 1);
end