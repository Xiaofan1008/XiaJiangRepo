%% ========================================================================
%  Batch Manual Override for Responding Channels
%  Fixed-ISI Simultaneous vs Sequential Version
%
%  Purpose:
%  Safely rebuilds responding-channel labels manually for each:
%       stimulation set × amplitude × selected PTD
%
%  Data structure expected:
%       Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive
%
%  Manual rebuild logic:
%       1) First force selected channels to SILENT
%       2) Then force selected channels to RESPONDING
%
%  If Silent = 1:64, then the final responding channels for that condition
%  will be exactly the channels listed in Respond.
%
%  Automatically creates a safe backup before overwriting.
% ========================================================================
clear;

%% ========================================================================
%  1. FILE PATH SETTINGS
% ========================================================================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_Exp1_Sim4';

% ========================================================================
%  [MODIFIED SECTION 1]
%  PTD SELECTION
%
%  Your Responding file still has a .ptd layer:
%       Responding.set(si).amp(ai).ptd(pi).channel
%
%  For current sequential 5 ms analysis:
%       Target_PTD_ms = 5;
%
%  For simultaneous 0 ms analysis:
%       Target_PTD_ms = 0;
%
%  This script will search for the matching PTD inside each Set/Amp block.
% ========================================================================
Target_PTD_ms = 0;
PTD_Tol_ms    = 0.01;

%% ========================================================================
%  2. OVERRIDE CONTROL PANEL
%  Grouped by stimulation set.
%
%  Each row in AmpFixes is:
%       Amp (uA),  Force_Silent,  Force_Respond
%
%  For full manual rebuild:
%       Force_Silent  = 1:64
%       Force_Respond = final manually selected responding channels
% ========================================================================
Conditions = []; % Initialize

% ================ Stimulation Set 1 ================
c = 1;
Conditions(c).Set = 1;
Conditions(c).AmpFixes = {
    % Amp (uA),  Force_Silent,  Force_Respond
      % 0,         1:64,          [];
      1,         1:64,          [36:40,42:48,52:53,56:58,62,64];
      2,         1:64,          [36:40,42:48,50:64];
      3,         1:64,          [36:40,42:48,50:64];
      4,         1:64,          [36:40,42:48,50:64];
      5,         1:64,          [36:40,42:48,50:64];
      % 6,         1:64,        [];
      8,         1:64,          [36:40,42:48,50:64];
      10,        1:64,          [36:40,42:48,50:64];
};

% ================ Stimulation Set 2 ================
c = 2;
Conditions(c).Set = 2;
Conditions(c).AmpFixes = {
    % Amp (uA),  Force_Silent,  Force_Respond
      % 0,         1:64,          [];
      1,         1:64,          [];
      2,         1:64,          [];
      3,         1:64,          [];
      4,         1:64,          [];
      5,         1:64,          [];
      % 6,         1:64,        [];
      8,         1:64,          [];
      10,        1:64,          [];
};

%% ========================================================================
%  3. INITIALIZATION & BACKUP
% ========================================================================
if ~isfolder(data_folder)
    error('Data folder does not exist: %s', data_folder);
end
cd(data_folder);

% Find the RespondingChannels file
file_list = dir('*_RespondingChannels.mat');

if isempty(file_list)
    error('Could not find a *_RespondingChannels.mat file in this folder.');
end

target_file = file_list(1).name;
full_path = fullfile(data_folder, target_file);

fprintf('\nLoading File: %s\n', target_file);
load(full_path, 'Responding'); % Load the main structure

% Also load the other variables to save them back safely
S_all = load(full_path);

% Create a Backup
% Only create one backup if it does not already exist, to avoid overwriting
% the original backup.
backup_name = strrep(target_file, '.mat', '_BACKUP.mat');
if ~isfile(backup_name)
    copyfile(target_file, backup_name);
    fprintf('Created Safe Backup: %s\n', backup_name);
else
    fprintf('Backup already exists. Safe to proceed.\n');
end

fprintf('\nStarting Batch Overrides...\n');
fprintf('Target PTD: %.3f ms\n', Target_PTD_ms);
fprintf('--------------------------------------------------\n');

%% ========================================================================
%  4. EXECUTE OVERRIDES
% ========================================================================
total_changed = 0;

% Outer Loop:
% Iterate through each stimulation set.
for c = 1:length(Conditions)

    cond = Conditions(c);
    si = cond.Set;

    % Failsafe:
    % Check whether the stimulation set exists in Responding.
    if si > numel(Responding.set)
        fprintf('WARNING [Condition %d]: Set %d does not exist. Skipping.\n', c, si);
        continue;
    end

    % Inner Loop:
    % Iterate through each amplitude listed under this stimulation set.
    for f = 1:size(cond.AmpFixes, 1)

        target_amp    = cond.AmpFixes{f, 1};
        force_silent  = cond.AmpFixes{f, 2};
        force_respond = cond.AmpFixes{f, 3};

        % Search for the correct amplitude index.
        ai_found = 0;
        for ai = 1:numel(Responding.set(si).amp)
            if abs(Responding.set(si).amp(ai).amp_value - target_amp) < 1e-4
                ai_found = ai;
                break;
            end
        end

        if ai_found == 0
            fprintf('WARNING [Condition %d]: Amp %.1f uA not found in Set %d. Skipping.\n', ...
                c, target_amp, si);
            continue;
        end

        ai = ai_found;

        % =================================================================
        % [MODIFIED SECTION 2]
        % Search for the correct PTD index inside this Set/Amp condition.
        %
        % Your file has:
        %     Responding.set(si).amp(ai).ptd(pi).channel
        %
        % Therefore we need to find pi before modifying channels.
        % =================================================================
        if ~isfield(Responding.set(si).amp(ai), 'ptd')
            fprintf('WARNING [Condition %d]: Set %d | Amp %.1f uA has no .ptd field. Skipping.\n', ...
                c, si, target_amp);
            continue;
        end

        pi_found = 0;
        for pi = 1:numel(Responding.set(si).amp(ai).ptd)

            % Prefer PTD_ms if it exists
            if isfield(Responding.set(si).amp(ai).ptd(pi), 'PTD_ms')
                this_ptd_ms = Responding.set(si).amp(ai).ptd(pi).PTD_ms;

            % Fallback: use PTD_us if only PTD_us exists
            elseif isfield(Responding.set(si).amp(ai).ptd(pi), 'PTD_us')
                this_ptd_ms = Responding.set(si).amp(ai).ptd(pi).PTD_us / 1000;

            else
                fprintf('WARNING: Set %d | Amp %.1f uA | PTD index %d has no PTD_ms or PTD_us field.\n', ...
                    si, target_amp, pi);
                continue;
            end

            if abs(this_ptd_ms - Target_PTD_ms) < PTD_Tol_ms
                pi_found = pi;
                break;
            end
        end

        if pi_found == 0
            fprintf('WARNING [Condition %d]: Set %d | Amp %.1f uA | PTD %.3f ms not found. Skipping.\n', ...
                c, si, target_amp, Target_PTD_ms);
            continue;
        end

        pi = pi_found;

        fprintf('Cond %d -> Set %d | %.1f uA | PTD %.3f ms:\n', ...
            c, si, target_amp, Target_PTD_ms);

        % =================================================================
        % [MODIFIED SECTION 3]
        % Channel access now uses:
        %     Responding.set(si).amp(ai).ptd(pi).channel
        %
        % instead of:
        %     Responding.set(si).amp(ai).channel
        % =================================================================

        % Force Silent
        % For manual rebuild, this is usually 1:64, meaning all channels
        % are first reset to non-responsive.
        for ch = force_silent
            if ch <= numel(Responding.set(si).amp(ai).ptd(pi).channel)
                Responding.set(si).amp(ai).ptd(pi).channel(ch).is_responsive = false;
                fprintf('   [-] Ch %02d forced SILENT\n', ch);
                total_changed = total_changed + 1;
            else
                fprintf('   [!] Ch %02d does not exist in data.\n', ch);
            end
        end

        % Force Responding
        % These channels are then turned back on as the final manually
        % selected responding channels for this set, amplitude, and PTD.
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

    % Overwrite the Responding struct in the loaded data.
    S_all.Responding = Responding;

    % Save back to the original file.
    save(full_path, '-struct', 'S_all');

    fprintf('COMPLETE: %d channel states updated and saved to main file.\n\n', total_changed);

else
    fprintf('No changes were made. File not overwritten.\n\n');
end