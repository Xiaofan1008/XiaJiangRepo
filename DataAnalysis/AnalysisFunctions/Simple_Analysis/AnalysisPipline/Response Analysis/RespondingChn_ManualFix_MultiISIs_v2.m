%% ========================================================================
%  Batch Manual Override for Responding Channels (Grouped Version)
%  Safely forces specific noisy channels to SILENT, or rescued channels
%  to RESPONDING, across multiple experimental conditions at once.
%  Automatically creates a safe backup before overwriting.
% ========================================================================
clear;

%% ========================================================================
%  1. FILE PATH SETTINGS
% ========================================================================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX020/Xia_ISI_10uA_SimSeq1';

%% ========================================================================
%  2. OVERRIDE CONTROL PANEL (Grouped by Set and Amplitude)
%  Define the base parameters for a condition, then list the ISIs and 
%  their corresponding responding channels in the "Fixes" table.
% ========================================================================
Conditions = []; % Initialize

% ================ Condition 1 ================
% c = 1;
% Conditions(c).Set    = 1;
% Conditions(c).Amp    = 5;
% Conditions(c).Silent = 1:64; % Default silent channels for this entire block
% Conditions(c).Fixes  = {
%     % PTD (ms),  [Force_Respond]
%       0,         [61];
%       3,         [41,54,61];
%       4,         [35,40,41,42,50,54,56,61];
%       5,         [35,40,41,42,45,47,50,53,54,56,61,63];
%       6,         [35,40,41,42,45,47,50,53,54,56,61,62,63];
%       7,         [35,37,39:42,44:48,50,53:57,59,61:63];
%       8,         [35,37,39:42,44:48,50,53:57,59,61:63];
%       9,         [35:37,39:42,44:48,50,53:59,61:64];
%       10,        [35:37,39:42,44:48,50,53:59,61:64];
%       11,        [35:37,39:42,44:48,50,54:57,59,61:64];
%       12,        [35:37,39:42,44:48,50,56:59,62:64];
%       13,        [35:37,39:42,44:46,50,56:59,62:63];
%       14,        [35:37,39:41,44:47,50,57,59,62:64];
%       15,        [35:37,40:42,44:47,50,62:63];
%       17,        [35:37,39:41,45:47,50,57,59,62,64];
%       20,        [36:37,39:42,45:47,53,59];
% };

% ================ Condition 2 ================
c = 1;
Conditions(c).Set    = 1;
Conditions(c).Amp    = 10;
Conditions(c).Silent = 1:64;
Conditions(c).Fixes  = {
      0,         [33:40,42:53,55:64];
      3,         [33:40,42:64];
      4,         [33:40,42:64];
      5,         [33:40,42:64];
      6,         [33:40,42:64];
      7,         [33:40,42:64];
      8,         [33:40,42:64];
      9,         [33:40,42:64];
      10,        [33:40,42:64];
      11,        [33:40,42:64];
      12,        [33:40,42:64];
      13,        [33:40,42:64];
      14,        [33:40,42:64];
      15,        [33:40,42:64];
      17,        [33:40,42:64];
      20,        [33:40,42:64];
};

% ================ Condition 3 ================
% c = 3;
% Conditions(c).Set    = 2;
% Conditions(c).Amp    = 5;
% Conditions(c).Silent = 1:64;
% Conditions(c).Fixes  = {
%       0,         [50,55];
%       3,         [40,50,52,57,63];
%       4,         [40:42,45,48,50,52,53,61,62];
%       5,         [37,41:42,44,45,48,50,52,53,55,57,61,62];
%       6,         [35:37,41:42,44:47,48,50:53,55:63];
%       7,         [35:37,39:41,44:48,50:53,55:63];
%       8,         [35:37,39:42,44:48,50:63];
%       9,         [35:37,39:42,44:48,50:63];
%       10,        [35:37,39:42,44:48,50:61,63];
%       11,        [35:37,39:42,44:48,50:63];
%       12,        [35:37,39:41,44:46,48,50:51,55:56,59:61,63];
%       13,        [35:37,39:41,44:46,48,50:51,55,59:61,63];
%       14,        [35:37,39:41,44:46,48,50:51,55,59,61,63];
%       15,        [35:37,39:41,45:46,48,50:51,55,59,61,63];
%       17,        [35:37,39:41,45:46,48,50:51,55,59,61,63];
%       20,        [35:37,39:41,44:46,48,50:51,54,55,59,61,63];
% };

% ================ Condition 4 ================
c = 2;
Conditions(c).Set    = 2;
Conditions(c).Amp    = 10;
Conditions(c).Silent = 1:64;
Conditions(c).Fixes  = {
      0,         [33:40,42:64];
      3,         [33:40,42:64];
      4,         [33:40,42:64];
      5,         [33:40,42:64];
      6,         [33:40,42:64];
      7,         [33:40,42:64];
      8,         [33:40,42:64];
      9,         [33:40,42:64];
      10,        [33:40,42:64];
      11,        [33:40,42:64];
      12,        [33:40,42:64];
      13,        [33:40,42:64];
      14,        [33:40,42:64];
      15,        [33:40,42:64];
      17,        [33:40,42:64];
      20,        [33:40,42:64];
};

%% ========================================================================
%  3. INITIALIZATION & BACKUP
% ========================================================================
if ~isfolder(data_folder)
    error('Data folder does not exist: %s', data_folder);
end
cd(data_folder);

% Find the RespondingChannels file
file_list = dir('*_MultiISIRespondingChannels.mat');
if isempty(file_list)
    error('Could not find a *_MultiISIRespondingChannels.mat file in this folder.');
end
target_file = file_list(1).name;
full_path = fullfile(data_folder, target_file);

fprintf('\nLoading File: %s\n', target_file);
load(full_path, 'Responding'); % Load the main structure

% Also load the other variables to save them back safely
S_all = load(full_path); 

% Create a Backup (Only if one doesn't exist yet to prevent overwriting the original backup)
backup_name = strrep(target_file, '.mat', '_BACKUP.mat');
if ~isfile(backup_name)
    copyfile(target_file, backup_name);
    fprintf('Created Safe Backup: %s\n', backup_name);
else
    fprintf('Backup already exists. Safe to proceed.\n');
end

fprintf('\nStarting Batch Overrides...\n');
fprintf('--------------------------------------------------\n');

%% ========================================================================
%  4. EXECUTE OVERRIDES
% ========================================================================
total_changed = 0;

% Outer Loop: Iterate through each major Condition (Set & Amp)
for c = 1:length(Conditions)
    cond = Conditions(c);
    si = cond.Set;
    target_amp = cond.Amp;
    force_silent = cond.Silent;
    
    % Failsafe: Check if the Set exists
    if si > numel(Responding.set)
        fprintf('WARNING [Condition %d]: Set %d does not exist. Skipping.\n', c, si);
        continue;
    end
    
    % Search for the correct Amplitude index
    ai_found = 0;
    for ai = 1:numel(Responding.set(si).amp)
        if abs(Responding.set(si).amp(ai).amp_value - target_amp) < 1e-4
            ai_found = ai;
            break;
        end
    end
    
    if ai_found == 0
        fprintf('WARNING [Condition %d]: Amp %.1f uA not found in Set %d. Skipping.\n', c, target_amp, si);
        continue;
    end
    
    ai = ai_found;
    
    % Inner Loop: Iterate through the specific Fixes (ISIs) for this condition
    for f = 1:size(cond.Fixes, 1)
        target_ptd = cond.Fixes{f, 1};
        force_respond = cond.Fixes{f, 2};
        
        % Search for the correct PTD index
        pi_found = 0;
        for pi = 1:numel(Responding.set(si).amp(ai).ptd)
            if abs(Responding.set(si).amp(ai).ptd(pi).PTD_ms - target_ptd) < 1e-4
                pi_found = pi;
                break;
            end
        end
        
        if pi_found == 0
            fprintf('WARNING [Cond %d]: PTD %.1f ms not found. Skipping.\n', c, target_ptd);
            continue;
        end
        
        pi = pi_found;
        fprintf('Cond %d -> Set %d | %.1f uA | %2.1f ms:\n', c, si, target_amp, target_ptd);
        
        % Force Silent (Turn OFF)
        for ch = force_silent
            if ch <= numel(Responding.set(si).amp(ai).ptd(pi).channel)
                Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = false;
                fprintf('   [-] Ch %02d forced SILENT\n', ch);
                total_changed = total_changed + 1;
            else
                fprintf('   [!] Ch %02d does not exist in data.\n', ch);
            end
        end
        
        % Force Respond (Turn ON)
        for ch = force_respond
            if ch <= numel(Responding.set(si).amp(ai).ptd(pi).channel)
                Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = true;
                fprintf('   [+] Ch %02d forced RESPONDING\n', ch);
                total_changed = total_changed + 1;
            else
                fprintf('   [!] Ch %02d does not exist in data.\n', ch);
            end
        end
        
        if isempty(force_silent) && isempty(force_respond)
            fprintf('   (No channels modified for this condition)\n');
        end
    end
end
fprintf('--------------------------------------------------\n');

%% ========================================================================
%  5. SAVE FINAL DATA
% ========================================================================
if total_changed > 0
    % Overwrite the 'Responding' struct in the loaded data
    S_all.Responding = Responding;
    
    % Save back to the original file
    save(full_path, '-struct', 'S_all');
    fprintf('COMPLETE: %d channel states updated and saved to main file.\n\n', total_changed);
else
    fprintf('No changes were made. File not overwritten.\n\n');
end