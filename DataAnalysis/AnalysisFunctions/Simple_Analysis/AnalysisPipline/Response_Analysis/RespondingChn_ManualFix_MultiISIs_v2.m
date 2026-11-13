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
data_folder = '/Volumes/MACData/Data/Data_Xia/DX014/Xia_Seq_Sim3';

%% ========================================================================
%  2. OVERRIDE CONTROL PANEL (Grouped by Set and Amplitude)
%  Define the base parameters for a condition, then list the ISIs and 
%  their corresponding responding channels in the "Fixes" table.
% ========================================================================
Conditions = []; % Initialize

% ================ Condition 1 ================
c = 1;
Conditions(c).Set    = 1;
Conditions(c).Amp    = 5;
Conditions(c).Silent = 1:64; % Default silent channels for this entire block
Conditions(c).Fixes  = {
    % PTD (ms),  [Force_Respond]
      0, [1:3,5,9:12,14,16];
      3, [];
      4, [];
      5, [2:4,7:16];
      6, [];
      7, [];
      8, [1:4,6:16];
      9, [];
      10,[1:16];
      11,[];
      12,[1:16];
      13,[];
      14,[];
      15,[2:16];
      17,[2:16];
      20,[2:16];
      25,[1:14];
};

% ================ Condition 2 ================
c = 2;
Conditions(c).Set    = 1;
Conditions(c).Amp    = 10;
Conditions(c).Silent = 1:64;
Conditions(c).Fixes  = {
      0, [1:6,10,12:16];
      3, [];
      4, [];
      5, [1:16];
      6, [];
      7, [];
      8, [1:16];
      9, [];
      10,[1:16];
      11,[];
      12,[1:16];
      13,[];
      14,[];
      15,[1:16];
      17,[1:16];
      20,[1:16];
      25,[1:16];
};

% ================ Condition 3 ================
c = 3;
Conditions(c).Set    = 2;
Conditions(c).Amp    = 5;
Conditions(c).Silent = 1:64;
Conditions(c).Fixes  = {
      0, [1,3,6,8,13:16];
      3, [];
      4, [];
      5, [1:16];
      6, [];
      7, [];
      8, [1:16];
      9, [];
      10,[1:16];
      11,[];
      12,[1:16];
      13,[];
      14,[];
      15,[1:16];
      17,[1:16];
      20,[1:16];
      25,[1:16];
};

% ================ Condition 4 ================
c = 4;
Conditions(c).Set    = 2;
Conditions(c).Amp    = 10;
Conditions(c).Silent = 1:64;
Conditions(c).Fixes  = {
      0, [1:3,6,8,12,14:16];
      3, [];
      4, [];
      5, [1:16];
      6, [];
      7, [];
      8, [1:16];
      9, [];
      10,[1:16];
      11,[];
      12,[1:16];
      13,[];
      14,[];
      15,[1:10,12:16];
      17,[1:10,12:16];
      20,[1:16];
      25,[1:16];
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