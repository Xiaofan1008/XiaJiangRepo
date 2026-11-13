%% ========================================================================
%  Batch Manual Override for Responding Channels
%  Mixed-Prefix Multi-Electrode Stimulation Version
%  Input format:
%    Conditions(c).Set = setID;
%    Conditions(c).Fixes = {
%        % Type,   Amp,  PrefixOrActiveCount,  ISI_ms,  ForceSilent,  ForceRespond
%          'seq',  5,    5,                    5,       1:64,         [1:10];
%          'sim',  5,    5,                    NaN,     1:64,         [1:8];
%    };
%  Meaning:
%    For 'seq':
%       PrefixOrActiveCount = PrefixLength
%       ISI_ms              = ISI
%       Modifies:
%           Responding.set(si).amp(ai).prefix(pi).channel(ch).is_responsive
%
%    For 'sim':
%       PrefixOrActiveCount = ActiveElectrodeCount
%       ISI_ms              = ignored
%       Modifies:
%           Responding.set(si).amp(ai).sim.channel(ch).is_responsive
% ========================================================================

clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));
%% ========================================================================
%  1. FILE PATH SETTINGS
% ========================================================================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';
% Use 'FullSeqAndSim' if you are only fixing Prefix 5 vs Sim.
% Use 'AllPrefixesAndSim' if you are fixing Prefix 1-5 and Sim.
RespondingFileMode = 'FullSeqAndSim';

% Matching tolerance.
Amp_Tol_uA = 1e-4;
ISI_Tol_ms = 1e-4;

%% ========================================================================
%  2. OVERRIDE CONTROL PANEL
%
%  Each row:
%       Type, Amp, PrefixOrActiveCount, ISI_ms, ForceSilent, ForceRespond
%  Type:
%       'seq' = sequential prefix condition
%       'sim' = simultaneous condition
%  For 'seq':
%       PrefixOrActiveCount = PrefixLength
%       ISI_ms = ISI
%  For 'sim':
%       PrefixOrActiveCount = ActiveElectrodeCount
%       ISI_ms = NaN
%  Example:
%       'seq', 5, 5, 5, 1:64, [1:12]
%       means:
%           Set X, 5 uA, Prefix 5, ISI 5 ms:
%           first turn off channels 1:64,
%           then force channels 1:12 to responsive.
%
%       'sim', 5, 5, NaN, 1:64, [1:8]
%       means:
%           Set X, 5 uA, simultaneous 5-electrode condition:
%           first turn off channels 1:64,
%           then force channels 1:8 to responsive.
% ========================================================================

Conditions = [];

% ================= Set 1 =================
c = 1;
Conditions(c).Set = 1;
Conditions(c).Fixes = {
    % Type,   Amp,  PrefixOrActiveCount,  ISI_ms,  ForceSilent,  ForceRespond
      'seq',  5,    5,                    5,       1:64,         [1:4,6:10,17:18,21,25];
      'seq',  10,   5,                    5,       1:64,         [1:4,6:10,17:19,21,23,25];
      'sim',  5,    5,                    NaN,     1:64,         [1,2,4,6,7,9:10,19,21,23];
      'sim',  10,   5,                    NaN,     1:64,         [1:4,6:10,17:19,21,23,25];
};

% ================= Set 2 =================
c = 2;
Conditions(c).Set = 2;
Conditions(c).Fixes = {
    % Type,   Amp,  PrefixOrActiveCount,  ISI_ms,  ForceSilent,  ForceRespond
      'seq',  5,    5,                    5,       1:64,         [1:4,6:10,17:19,21,23,25];
      'seq',  10,   5,                    5,       1:64,         [1:4,6:10,17:19,21,23,25];
};

% ================= Set 3 =================
c = 3;
Conditions(c).Set = 3;
Conditions(c).Fixes = {
    % Type,   Amp,  PrefixOrActiveCount,  ISI_ms,  ForceSilent,  ForceRespond
      'seq',  5,    5,                    5,       1:64,         [1:4,6:10,17:19,21,23,25];
      'seq',  10,   5,                    5,       1:64,         [1:4,6:10,17:19,21,23,25];
};

%% ========================================================================
%  3. INITIALIZATION
% ========================================================================

if ~isfolder(data_folder)
    error('Data folder does not exist: %s', data_folder);
end

cd(data_folder);

fprintf('\nRunning manual responding-channel override in:\n%s\n\n', data_folder);

%% ----- Build base_name -----
parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');

if numel(u) >= 4
    base_name = last_folder(1:u(end-1)-1);
else
    base_name = last_folder;
end

fprintf('Base name: %s\n', base_name);

%% ----- Target responding-channel file -----
target_file = sprintf('%s_RespondingChannels_%s.mat', base_name, RespondingFileMode);
full_path = fullfile(data_folder, target_file);

if ~isfile(full_path)
    error('Could not find responding-channel file:\n%s', full_path);
end

fprintf('Loading responding-channel file:\n%s\n', target_file);

S_all = load(full_path);

if ~isfield(S_all, 'Responding')
    error('Responding variable not found in %s', target_file);
end

Responding = S_all.Responding;

%% ========================================================================
%  4. CREATE BACKUP
% ========================================================================

backup_name = strrep(target_file, '.mat', '_BACKUP.mat');
backup_path = fullfile(data_folder, backup_name);

if ~isfile(backup_path)
    copyfile(full_path, backup_path);
    fprintf('Created safe backup:\n%s\n', backup_name);
else
    fprintf('Backup already exists:\n%s\n', backup_name);
end

%% ========================================================================
%  5. EXECUTE OVERRIDES
% ========================================================================

fprintf('\nStarting batch overrides...\n');
fprintf('--------------------------------------------------\n');

total_changed = 0;
total_skipped = 0;

OverrideRows = {};
rowCounter = 0;

for c = 1:numel(Conditions)

    cond = Conditions(c);

    if ~isfield(cond, 'Set')
        fprintf('WARNING [Condition block %d]: Missing Set field. Skipping.\n', c);
        total_skipped = total_skipped + 1;
        continue;
    end

    if ~isfield(cond, 'Fixes') || isempty(cond.Fixes)
        fprintf('WARNING [Condition block %d]: Missing or empty Fixes table. Skipping.\n', c);
        total_skipped = total_skipped + 1;
        continue;
    end

    target_set = cond.Set;

    %% ----- Find set index safely by setID -----
    si = find_set_index(Responding, target_set);

    if si == 0
        fprintf('WARNING [Block %d]: SetID %d not found in Responding. Skipping.\n', ...
            c, target_set);
        total_skipped = total_skipped + 1;
        continue;
    end

    FixesTable = cond.Fixes;

    for f = 1:size(FixesTable, 1)

        target_type = lower(string(FixesTable{f, 1}));
        target_amp  = FixesTable{f, 2};
        target_count_or_prefix = FixesTable{f, 3};
        target_isi  = FixesTable{f, 4};
        force_silent  = FixesTable{f, 5};
        force_respond = FixesTable{f, 6};

        %% ----- Find amplitude index -----
        ai = find_amp_index(Responding, si, target_amp, Amp_Tol_uA);

        if ai == 0
            fprintf('WARNING [Block %d Row %d]: Set %d | Amp %.3f uA not found. Skipping.\n', ...
                c, f, target_set, target_amp);
            total_skipped = total_skipped + 1;
            continue;
        end

        %% ----- Find condition object -----
        condition_label = '';
        prefix_length = NaN;
        active_count  = NaN;
        pi = NaN;

        if target_type == "seq" || target_type == "prefix"

            condition_label = 'seq';
            prefix_length = target_count_or_prefix;

            pi = find_prefix_index(Responding, si, ai, prefix_length, target_isi, ISI_Tol_ms);

            if pi == 0
                fprintf('WARNING [Block %d Row %d]: Set %d | Seq | Amp %.1f | Prefix %.1f | ISI %.1f not found. Skipping.\n', ...
                    c, f, target_set, target_amp, prefix_length, target_isi);
                total_skipped = total_skipped + 1;
                continue;
            end

            if ~isfield(Responding.set(si).amp(ai).prefix(pi), 'channel')
                fprintf('WARNING [Block %d Row %d]: Prefix condition found but channel field missing. Skipping.\n', c, f);
                total_skipped = total_skipped + 1;
                continue;
            end

            nChannelStruct = numel(Responding.set(si).amp(ai).prefix(pi).channel);

            fprintf('\nSet %d | Seq | Amp %.1f uA | Prefix %.1f | ISI %.1f ms:\n', ...
                target_set, target_amp, prefix_length, target_isi);

            %% ----- Force silent first -----
            [Responding, changed1, rows1] = apply_channel_override_prefix( ...
                Responding, si, ai, pi, target_set, target_amp, ...
                prefix_length, active_count, target_isi, ...
                force_silent, false, 'ForceSilent', nChannelStruct);

            %% ----- Force responding second -----
            [Responding, changed2, rows2] = apply_channel_override_prefix( ...
                Responding, si, ai, pi, target_set, target_amp, ...
                prefix_length, active_count, target_isi, ...
                force_respond, true, 'ForceRespond', nChannelStruct);

        elseif target_type == "sim" || target_type == "simultaneous"

            condition_label = 'sim';
            active_count = target_count_or_prefix;

            if ~isfield(Responding.set(si).amp(ai), 'sim') || ...
               isempty(Responding.set(si).amp(ai).sim)

                fprintf('WARNING [Block %d Row %d]: Set %d | Sim | Amp %.1f has no sim field. Skipping.\n', ...
                    c, f, target_set, target_amp);
                total_skipped = total_skipped + 1;
                continue;
            end

            if ~isfield(Responding.set(si).amp(ai).sim, 'channel')
                fprintf('WARNING [Block %d Row %d]: Sim condition found but channel field missing. Skipping.\n', c, f);
                total_skipped = total_skipped + 1;
                continue;
            end

            nChannelStruct = numel(Responding.set(si).amp(ai).sim.channel);

            fprintf('\nSet %d | Sim | Amp %.1f uA | ActiveCount %.1f:\n', ...
                target_set, target_amp, active_count);

            %% ----- Force silent first -----
            [Responding, changed1, rows1] = apply_channel_override_sim( ...
                Responding, si, ai, target_set, target_amp, ...
                prefix_length, active_count, target_isi, ...
                force_silent, false, 'ForceSilent', nChannelStruct);

            %% ----- Force responding second -----
            [Responding, changed2, rows2] = apply_channel_override_sim( ...
                Responding, si, ai, target_set, target_amp, ...
                prefix_length, active_count, target_isi, ...
                force_respond, true, 'ForceRespond', nChannelStruct);

        else

            fprintf('WARNING [Block %d Row %d]: Unknown type "%s". Use ''seq'' or ''sim''. Skipping.\n', ...
                c, f, target_type);
            total_skipped = total_skipped + 1;
            continue;
        end

        %% ----- Store override rows -----
        allRows = [rows1; rows2];

        for rr = 1:size(allRows, 1)
            rowCounter = rowCounter + 1;
            OverrideRows(rowCounter,:) = allRows(rr,:); %#ok<SAGROW>
        end

        changed_this_row = changed1 + changed2;
        total_changed = total_changed + changed_this_row;

        if isempty(force_silent) && isempty(force_respond)
            fprintf('   (No channels modified for this condition)\n');
        elseif changed_this_row == 0
            fprintf('   (Requested channels were already in requested states)\n');
        end
    end
end

fprintf('--------------------------------------------------\n');
fprintf('Override finished.\n');
fprintf('Total channel states changed: %d\n', total_changed);
fprintf('Total skipped rows/blocks: %d\n', total_skipped);

%% ========================================================================
%  6. CREATE MANUAL OVERRIDE TABLE
% ========================================================================

if isempty(OverrideRows)

    NewManualOverrideTable = table();

else

    NewManualOverrideTable = cell2table(OverrideRows, ...
        'VariableNames', { ...
        'SetID', ...
        'ConditionType', ...
        'Amp_uA', ...
        'PrefixLength', ...
        'ActiveCount', ...
        'ISI_ms', ...
        'Channel', ...
        'Action', ...
        'NewValue', ...
        'OldValue'});

end

if isfield(S_all, 'ManualOverrideTable') && ~isempty(S_all.ManualOverrideTable)

    if ~isempty(NewManualOverrideTable)
        ManualOverrideTable = [S_all.ManualOverrideTable; NewManualOverrideTable];
    else
        ManualOverrideTable = S_all.ManualOverrideTable;
    end

else

    ManualOverrideTable = NewManualOverrideTable;
end

ManualOverrideInfo = struct();
ManualOverrideInfo.data_folder = data_folder;
ManualOverrideInfo.base_name = base_name;
ManualOverrideInfo.RespondingFileMode = RespondingFileMode;
ManualOverrideInfo.target_file = target_file;
ManualOverrideInfo.backup_file = backup_name;
ManualOverrideInfo.Amp_Tol_uA = Amp_Tol_uA;
ManualOverrideInfo.ISI_Tol_ms = ISI_Tol_ms;
ManualOverrideInfo.format = ...
    'Conditions(c).Fixes = {Type, Amp, PrefixOrActiveCount, ISI_ms, ForceSilent, ForceRespond}';
ManualOverrideInfo.last_modified_on = datestr(now);
ManualOverrideInfo.total_changed_this_run = total_changed;
ManualOverrideInfo.total_skipped_this_run = total_skipped;

%% ========================================================================
%  7. SAVE FINAL DATA
% ========================================================================

if total_changed > 0

    S_all.Responding = Responding;
    S_all.ManualOverrideTable = ManualOverrideTable;
    S_all.ManualOverrideInfo = ManualOverrideInfo;
    S_all.Conditions_last_override = Conditions;

    save(full_path, '-struct', 'S_all');

    fprintf('\nCOMPLETE: %d channel states updated and saved to:\n%s\n\n', ...
        total_changed, target_file);

else

    fprintf('\nNo channel states changed. Responding file was not overwritten.\n\n');

end

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function si = find_set_index(Responding, target_set)

    si = 0;

    if ~isfield(Responding, 'set') || isempty(Responding.set)
        return;
    end

    for k = 1:numel(Responding.set)

        if isfield(Responding.set(k), 'setID')

            if Responding.set(k).setID == target_set
                si = k;
                return;
            end

        else

            % Fallback only if setID field does not exist.
            if k == target_set
                si = k;
                return;
            end
        end
    end
end

function ai = find_amp_index(Responding, si, target_amp, Amp_Tol_uA)

    ai = 0;

    if si > numel(Responding.set)
        return;
    end

    if ~isfield(Responding.set(si), 'amp') || isempty(Responding.set(si).amp)
        return;
    end

    for k = 1:numel(Responding.set(si).amp)

        if isfield(Responding.set(si).amp(k), 'amp_value')

            if abs(Responding.set(si).amp(k).amp_value - target_amp) < Amp_Tol_uA
                ai = k;
                return;
            end

        end
    end
end

function pi = find_prefix_index(Responding, si, ai, target_prefix, target_isi, ISI_Tol_ms)

    pi = 0;

    if ~isfield(Responding.set(si).amp(ai), 'prefix') || ...
       isempty(Responding.set(si).amp(ai).prefix)
        return;
    end

    for k = 1:numel(Responding.set(si).amp(ai).prefix)

        P = Responding.set(si).amp(ai).prefix(k);

        prefix_match = false;
        isi_match = true;

        %% ----- Prefix match -----
        if isfield(P, 'prefix_length')
            prefix_match = abs(P.prefix_length - target_prefix) < 1e-4;
        elseif isfield(P, 'prefixLength')
            prefix_match = abs(P.prefixLength - target_prefix) < 1e-4;
        else
            % Fallback: assume prefix array index equals prefix length.
            prefix_match = abs(k - target_prefix) < 1e-4;
        end

        %% ----- ISI match, only if available in structure -----
        if isfield(P, 'isi_ms')
            isi_match = abs(P.isi_ms - target_isi) < ISI_Tol_ms;
        elseif isfield(P, 'ISI_ms')
            isi_match = abs(P.ISI_ms - target_isi) < ISI_Tol_ms;
        else
            % If the responding file only contains one ISI, do not block matching.
            isi_match = true;
        end

        if prefix_match && isi_match
            pi = k;
            return;
        end
    end
end

function [Responding, changed_count, Rows] = apply_channel_override_prefix( ...
    Responding, si, ai, pi, target_set, target_amp, ...
    prefix_length, active_count, target_isi, ...
    ch_list, new_value, action_label, nChannelStruct)

    changed_count = 0;
    Rows = {};

    if isempty(ch_list)
        return;
    end

    for ch = ch_list

        if ch < 1 || ch > nChannelStruct
            fprintf('   [!] Ch %02d does not exist in prefix channel struct.\n', ch);
            continue;
        end

        old_value = Responding.set(si).amp(ai).prefix(pi).channel(ch).is_responsive;

        Responding.set(si).amp(ai).prefix(pi).channel(ch).is_responsive = new_value;

        if old_value ~= new_value
            changed_count = changed_count + 1;

            if new_value
                fprintf('   [+] Ch %02d forced RESPONDING\n', ch);
            else
                fprintf('   [-] Ch %02d forced SILENT\n', ch);
            end

            Rows(end+1,:) = { ... %#ok<AGROW>
                target_set, ...
                'seq', ...
                target_amp, ...
                prefix_length, ...
                active_count, ...
                target_isi, ...
                ch, ...
                action_label, ...
                new_value, ...
                old_value};
        end
    end
end

function [Responding, changed_count, Rows] = apply_channel_override_sim( ...
    Responding, si, ai, target_set, target_amp, ...
    prefix_length, active_count, target_isi, ...
    ch_list, new_value, action_label, nChannelStruct)

    changed_count = 0;
    Rows = {};

    if isempty(ch_list)
        return;
    end

    for ch = ch_list

        if ch < 1 || ch > nChannelStruct
            fprintf('   [!] Ch %02d does not exist in sim channel struct.\n', ch);
            continue;
        end

        old_value = Responding.set(si).amp(ai).sim.channel(ch).is_responsive;

        Responding.set(si).amp(ai).sim.channel(ch).is_responsive = new_value;

        if old_value ~= new_value
            changed_count = changed_count + 1;

            if new_value
                fprintf('   [+] Ch %02d forced RESPONDING\n', ch);
            else
                fprintf('   [-] Ch %02d forced SILENT\n', ch);
            end

            Rows(end+1,:) = { ... %#ok<AGROW>
                target_set, ...
                'sim', ...
                target_amp, ...
                prefix_length, ...
                active_count, ...
                target_isi, ...
                ch, ...
                action_label, ...
                new_value, ...
                old_value};
        end
    end
end