%% ============================================
%   MANUAL BAD TRIAL ENTRY (GLOBAL REMOVAL)
%   Fixed-ISI Simultaneous vs Sequential Version
%
%   Purpose:
%   Manually flag selected trials as "Bad" for ALL CHANNELS.
%
%   Manual input format:
%       Grouped by stimulation Set.
%       Each row under AmpTrials is:
%           Amp (uA), [List of Relative Bad Trial Indices]
%
%   Matching logic:
%       If Use_PTD_Filter = false:
%           Set + Amp -> Absolute Trial ID
%
%       If Use_PTD_Filter = true:
%           Set + Amp + Target_PTD_ms -> Absolute Trial ID
%
%   Safety:
%       - Appends to existing BadTrials file if it already exists.
%       - Keeps PTD extraction for checking fixed-ISI or mixed datasets.
%       - Uses fullfile() for robust saving.
% ============================================
clear;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

% ========================================================================
% [MODIFIED SECTION 1]
% Optional PTD filtering control
%
% Use_PTD_Filter = false:
%   Use this for current mode-specific folders, e.g. only sequential 5 ms
%   or only simultaneous 0 ms. Trials are grouped by Set + Amp only.
%
% Use_PTD_Filter = true:
%   Use this if one dataset contains mixed timing conditions, e.g. both
%   simultaneous 0 ms and sequential 5 ms trials. Trials are grouped by
%   Set + Amp + Target_PTD_ms.
%
% IMPORTANT:
%   The relative bad trial indices must be collected using the same grouping
%   logic as this script.
% ========================================================================
Use_PTD_Filter = false;   % false = group by Set + Amp only
Target_PTD_ms  = 5;       % only used when Use_PTD_Filter = true
PTD_Tol_ms     = 0.01;    % tolerance for matching PTD in ms

Amp_Tol_uA     = 0.01;    % tolerance for matching amplitude in uA

%% ---------------- USER SETTINGS ----------------
data_folder = '/Volumes/MACData/Data/Data_Xia/DX018/Xia_Exp1_Seq4';
Electrode_Type = 2;

%% ---------------- MANUAL BAD TRIALS LIST ----------------
% ========================================================================
% [MODIFIED SECTION 2]
% Manual bad trial list is now grouped as:
%
%   Bad_Groups(c).Set = stimulation set number;
%   Bad_Groups(c).AmpTrials = {
%       % Amp (uA), [relative bad trial indices]
%         0,        [];
%         1,        [];
%         2,        [1 4 8];
%   };
%
% If Use_PTD_Filter = false:
%   Relative bad trial index = trial number within Set + Amp.
%
% If Use_PTD_Filter = true:
%   Relative bad trial index = trial number within Set + Amp + Target_PTD_ms.
% ========================================================================
Bad_Groups = [];

% ================ Stimulation Set 1 ================
c = 1;
Bad_Groups(c).Set = 1;
Bad_Groups(c).AmpTrials = {
    % Amp (uA),  [List of Relative Bad Trial Indices]
      % 0,         [];
      1,         [1,2,4:6,9,10,12,15,19,22:25,29];
      2,         [2,3,5:8,11:12,15:16,18:19,23,28,30];
      3,         [2:4,8:11,13,18,27,29];
      4,         [1,2,4:6,9,11,13,15,19,21,26,28];
      5,         [1:4,6,8:13,15,17,18,20,21,24,26:28];
      % 6,       [];
      8,         [1,4,15,17,22,26,27,29];
      10,        [2,3,7,11,12,16,29];
};

% ================ Stimulation Set 2 ================
c = 2;
Bad_Groups(c).Set = 2;
Bad_Groups(c).AmpTrials = {
    % Amp (uA),  [List of Relative Bad Trial Indices]
      % 0,         [];
      1,         [2:9,11:13,16:19,24,26]
      2,         [1,2,6,10,15,21:24,27]
      3,         [1:7,9,11,14:15,19,20,26:28]
      4,         [2:6,9,11,14,16:18,20,22,27,28,30]
      5,         [1:4,9,14,24:27]
      % 6,       []
      8,         [4,7,9:11,13,15,17:18,20,23,24,26,29,30]
      10,        [4,6:9,12,14,15,17,21,23,25,27,29]
};

% You can add more stimulation sets using the same format:
%
% c = 3;
% Bad_Groups(c).Set = 3;
% Bad_Groups(c).AmpTrials = {
%     % Amp (uA),  [List of Relative Bad Trial Indices]
%       0,         [];
%       1,         [];
%       2,         [];
%       3,         [];
%       4,         [];
%       5,         [];
%       6,         [];
%       8,         [];
%       10,        [];
% };

%% ---------------- SETUP & CHECKS ----------------
if ~isfolder(data_folder), error('Folder not found'); end
cd(data_folder);
fprintf('\nRunning Global Manual Bad Trial Entry in:\n%s\n\n', data_folder);

parts = split(data_folder, filesep);
last_folder = parts{end};
u = strfind(last_folder, '_');
if numel(u) > 4
    base_name = last_folder(1 : u(end-1)-1);
else
    base_name = last_folder;
end

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

% ========================================================================
% [MODIFIED SECTION 3]
% Keep PTD extraction for safety checking and optional filtering.
%
% For sequential stimulation with simN > 1:
%   PTD is read from column 6 of the second stimulation row.
%
% For single/simultaneous-only style data with simN == 1:
%   PTD is set to 0 for all trials.
%
% trialPTDs is in ms.
% ========================================================================
if simN > 1
    PTD_all = cell2mat(StimParams(3:simN:end,6)); % µs from 2nd pulse row
else
    PTD_all = zeros(nTrials,1);
end
trialPTDs = PTD_all / 1000; % Convert to ms

% Print detected PTD values as a sanity check
unique_ptds_ms = unique(trialPTDs);
fprintf('Detected PTD values in this dataset: ');
fprintf('%.3f ', unique_ptds_ms);
fprintf('ms\n');

% If PTD filtering is enabled, stop if the target PTD does not exist.
if Use_PTD_Filter
    if ~any(abs(unique_ptds_ms - Target_PTD_ms) < PTD_Tol_ms)
        error('Target_PTD_ms %.3f ms was not found in this dataset.', Target_PTD_ms);
    else
        fprintf('PTD filtering ENABLED: using Target_PTD_ms = %.3f ms\n', Target_PTD_ms);
    end
else
    fprintf('PTD filtering DISABLED: grouping trials by Set + Amp only.\n');
end

% 3. Extract Stimulation Sets (CombClass)
E_MAP = S_exp.E_MAP;
stimNames = StimParams(2:end,1);
[~, idx_all] = ismember(stimNames, E_MAP(2:end));

% Reconstruct trial combinations
comb_seq = zeros(nTrials, simN);
for t = 1:nTrials
    rr = (t-1)*simN + (1:simN);
    v = idx_all(rr);
    v = v(v>0);
    comb_seq(t,1:numel(v)) = v(:).';
end
[~,~,combClass] = unique(comb_seq,'rows','stable');

%% ---------------- LOAD OR CREATE BAD TRIALS FILE ----------------
% ========================================================================
% [MODIFIED SECTION 4]
% Changed the output filename from MultiISIsBadTrials to SimSeqBadTrials.
% This avoids mixing this fixed-ISI sim-vs-seq bad-trial file with your
% previous multi-ISI bad-trial file.
% ========================================================================
bad_filename = sprintf('%s.BadTrials.mat', base_name);
save_file_path = fullfile(data_folder, bad_filename);

if isfile(save_file_path)
    fprintf('Found existing file: %s\nLoading previous bad trials...\n', bad_filename);
    tmp = load(save_file_path);
    if isfield(tmp, 'BadTrials')
        BadTrials = tmp.BadTrials;
    else
        BadTrials = cell(nCh_Total, 1);
    end
else
    fprintf('Creating new BadTrials file: %s\n', bad_filename);
    BadTrials = cell(nCh_Total, 1);
end

if numel(BadTrials) < nCh_Total
    BadTrials{nCh_Total} = [];
end

%% ---------------- PROCESS MANUAL ADDITIONS ----------------
% ========================================================================
% [MODIFIED SECTION 5]
% Processing logic now loops through:
%       Bad_Groups(c).Set
%       Bad_Groups(c).AmpTrials rows
%
% Instead of:
%       Bad_Groups(c).Set
%       Bad_Groups(c).Amp
%       Bad_Groups(c).Trials rows by PTD
%
% Matching is controlled by Use_PTD_Filter:
%   false -> Set + Amp
%   true  -> Set + Amp + Target_PTD_ms
% ========================================================================
if isempty(Bad_Groups)
    fprintf('\nNo entries. Done.\n');
else
    fprintf('\nProcessing manual entries (Applying to ALL channels)...\n');
    count_added = 0;

    % Outer Loop:
    % Iterate through each stimulation set group.
    for c = 1:length(Bad_Groups)

        tgt_set = Bad_Groups(c).Set;

        % Inner Loop:
        % Iterate through each amplitude row in the AmpTrials table.
        for f = 1:size(Bad_Groups(c).AmpTrials, 1)

            tgt_amp = Bad_Groups(c).AmpTrials{f, 1};
            bad_indices_list = Bad_Groups(c).AmpTrials{f, 2};

            % ------------------------------------------------------------
            % Find absolute trials matching this condition.
            %
            % If Use_PTD_Filter is false:
            %   matching_trials = all trials matching Set + Amp
            %
            % If Use_PTD_Filter is true:
            %   matching_trials = all trials matching Set + Amp + Target_PTD_ms
            % ------------------------------------------------------------
            if Use_PTD_Filter
                matching_trials = find(combClass == tgt_set & ...
                                       abs(trialAmps - tgt_amp) < Amp_Tol_uA & ...
                                       abs(trialPTDs - Target_PTD_ms) < PTD_Tol_ms);
            else
                matching_trials = find(combClass == tgt_set & ...
                                       abs(trialAmps - tgt_amp) < Amp_Tol_uA);
            end

            if isempty(matching_trials)
                if Use_PTD_Filter
                    fprintf('  [SKIP] Set %d, Amp %.1f, PTD %.3f ms not found.\n', ...
                        tgt_set, tgt_amp, Target_PTD_ms);
                else
                    fprintf('  [SKIP] Set %d, Amp %.1f not found.\n', ...
                        tgt_set, tgt_amp);
                end
                continue;
            end

            % ------------------------------------------------------------
            % Optional warning:
            % If PTD filtering is OFF but this Set + Amp contains multiple
            % PTD values, relative trial labels may mix timing conditions.
            % This does not stop the code, but warns you.
            % ------------------------------------------------------------
            if ~Use_PTD_Filter
                ptds_this_group = unique(trialPTDs(matching_trials));
                if numel(ptds_this_group) > 1
                    fprintf('  [WARNING] Set %d, Amp %.1f contains multiple PTDs: ', ...
                        tgt_set, tgt_amp);
                    fprintf('%.3f ', ptds_this_group);
                    fprintf('ms. Consider Use_PTD_Filter = true.\n');
                end
            end

            % Process each bad relative trial index.
            for b = 1:length(bad_indices_list)

                tgt_idx = bad_indices_list(b);

                if tgt_idx > numel(matching_trials)
                    if Use_PTD_Filter
                        fprintf('  [SKIP] Set %d, Amp %.1f, PTD %.3f ms: Req #%d, but only %d trials exist.\n', ...
                            tgt_set, tgt_amp, Target_PTD_ms, tgt_idx, numel(matching_trials));
                    else
                        fprintf('  [SKIP] Set %d, Amp %.1f: Req #%d, but only %d trials exist.\n', ...
                            tgt_set, tgt_amp, tgt_idx, numel(matching_trials));
                    end
                    continue;
                end

                % Convert relative trial index to absolute trial ID.
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
                    if Use_PTD_Filter
                        fprintf('  -> Added GLOBAL Bad Trial: Set %d | %.1f uA | PTD %.3f ms | Trial #%d (Abs: %d)\n', ...
                            tgt_set, tgt_amp, Target_PTD_ms, tgt_idx, abs_trial);
                    else
                        fprintf('  -> Added GLOBAL Bad Trial: Set %d | %.1f uA | Trial #%d (Abs: %d)\n', ...
                            tgt_set, tgt_amp, tgt_idx, abs_trial);
                    end
                    count_added = count_added + 1;
                end
            end % End bad index loop

        end % End amplitude row loop
    end % End stimulation set group loop

    fprintf('\nDone. Added %d new global bad trials.\n', count_added);
end

%% ---------------- SAVE ----------------
GoodTrials = cell(nCh_Total, 1);
for ich = 1:nCh_Total
    GoodTrials{ich} = setdiff(1:nTrials, BadTrials{ich});
end

fprintf('\nSaving to: %s ...\n', save_file_path);

if isfile(save_file_path)
    save(save_file_path, 'BadTrials', 'GoodTrials', '-append');
else
    detect_win_ms = [0 0];
    slide_win_size_ms = 0;
    max_spikes_in_slide = 0;

    save(save_file_path, 'BadTrials', 'GoodTrials', 'detect_win_ms', ...
        'slide_win_size_ms', 'max_spikes_in_slide');
end

fprintf('Save Complete.\n');