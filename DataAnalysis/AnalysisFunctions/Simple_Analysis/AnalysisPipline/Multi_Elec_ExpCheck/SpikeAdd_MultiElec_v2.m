% ============================================================
%  Prefix-Based Spike Recovery for Mixed-Prefix Stimulation
%  Purpose:
%    Recover spikes removed by artifact blanking in longer prefix trials.
%  Recovery rule:
%    Target Prefix N is recovered using Source Prefix N-1.
%  Important:
%    This version uses CASCADING recovery:
%      Prefix 2 uses original Prefix 1
%      Prefix 3 uses recovered Prefix 2
%      Prefix 4 uses recovered Prefix 3
%      Prefix 5 uses recovered Prefix 4
%  Output:
%    File:
%      *.sp_xia_PrefixRecovery.mat
%    Variable:
%      sp_seq
% ============================================================

clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER PARAMETERS ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX027/Xia_FixISI_5ms1';
% data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

FS = 30000;

% Recovery window starts relative to the first trigger.
% Usually 0 ms because blanking starts around the first artifact.
win_start_ms = 0;

% If no first post-blank spike is detected, use:
%   win_end = lastActivePTD_ms + fallback_ms
use_fallback = true;
fallback_ms  = 2;

% Safety cap:
% Do not recover beyond:
%   lastActivePTD_ms + max_extra_ms
%
% This avoids replacing too much real target response.
max_extra_ms = 10;

% Prefixes to recover.
% [] means recover all possible target prefixes from 2 to max prefix.
target_prefixes_to_recover = [];

% Optional filters.
% [] means all detected values.
amps_to_recover = [];
sets_to_recover = [];
isis_to_recover_ms = [];

% Source trial selection.
% Deterministic cycling means source trials are reused in order.
use_deterministic_cycling = true;

% If true, print one example trial per prefix to check alignment.
debug_print_trial_content = true;

% If true, print recovery window information for first few target trials.
debug_recovery_window = true;
debug_n_trials_per_group = 3;
debug_channel = 1;

%% ====================== LOAD FOLDER ======================
if ~isfolder(data_folder)
    error('Data folder does not exist.');
end
cd(data_folder);

fprintf('\nData folder:\n%s\n', data_folder);

%% ====================== LOAD TARGET SPIKES ======================
% Load the original spike file that will be recovered.
seq_sp_file = dir('*sp_xia.mat');
assert(~isempty(seq_sp_file), 'No *sp_xia.mat file found.');

load(seq_sp_file(1).name, 'sp_clipped');
sp_seq = sp_clipped;   % output variable name required by later analysis

nChn = numel(sp_seq);

fprintf('Loaded spike file: %s\n', seq_sp_file(1).name);

%% ====================== LOAD TRIGGERS ======================
if isempty(dir('*.trig.dat'))
    cleanTrig_sabquick;
end

trig = loadTrig(0);
nTrig = numel(trig);

fprintf('Number of triggers: %d\n', nTrig);

%% ====================== LOAD EXPERIMENT PARAMETERS ======================
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat file found.');

S = load(param_file(1).name);

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;   % rows/slots per trial
n_Trials          = S.n_Trials;

fprintf('n_Trials from exp file: %d\n', n_Trials);
fprintf('Rows/slots per trial: %d\n', simultaneous_stim);

if n_Trials ~= nTrig
    warning('n_Trials (%d) does not match number of triggers (%d). Using min of both.', ...
        n_Trials, nTrig);
end

nTrials_use = min(n_Trials, nTrig);

%% ====================== REMOVE HEADER ROW ======================
StimParams_data = StimParams(2:end,:);

expected_rows = n_Trials * simultaneous_stim;
if size(StimParams_data,1) ~= expected_rows
    warning('StimParams data rows (%d) do not equal n_Trials*simultaneous_stim (%d). Check file.', ...
        size(StimParams_data,1), expected_rows);
end

%% ====================== TRIAL METADATA FROM STIMPARAMS ======================
% Important:
%   Do NOT use separate metadata arrays directly here.
%   They may not be randomized in the same order as StimParams.
%
% Instead, derive metadata directly from randomized StimParams columns:
%   26 = ActiveElectrodeCount
%   27 = PrefixLength
%   28 = ISI_ms
%   29 = ConditionType
%   30 = ConditionSetID

if size(StimParams_data,2) < 30
    error('StimParams does not contain columns 26–30. Cannot use mixed-prefix metadata.');
end

firstRow_eachTrial = 1:simultaneous_stim:size(StimParams_data,1);

activeCount_trial    = cell2mat(StimParams_data(firstRow_eachTrial,26));
prefixLength_trial   = cell2mat(StimParams_data(firstRow_eachTrial,27));
isi_ms_trial         = cell2mat(StimParams_data(firstRow_eachTrial,28));
conditionType_trial  = cell2mat(StimParams_data(firstRow_eachTrial,29));
conditionSetID_trial = cell2mat(StimParams_data(firstRow_eachTrial,30));

% Force column vectors.
activeCount_trial    = activeCount_trial(:);
prefixLength_trial   = prefixLength_trial(:);
isi_ms_trial         = isi_ms_trial(:);
conditionType_trial  = conditionType_trial(:);
conditionSetID_trial = conditionSetID_trial(:);

fprintf('\nUsing trial metadata from StimParams columns 26–30.\n');

%% ====================== EXTRACT AMPLITUDE PER TRIAL ======================
trialAmps_all = cell2mat(StimParams_data(:,16));
trialAmps = trialAmps_all(firstRow_eachTrial);

% Convert inactive/zero-control amplitude from -1 to 0.
trialAmps(trialAmps == -1) = 0;
trialAmps = trialAmps(:);

%% ====================== CALCULATE FINAL ACTIVE ARTIFACT TIME ======================
% lastActivePTD_ms(tr) = final active artifact time for trial tr.
%
% For each active row:
%   final artifact time =
%       PTD_us + (PulseNum - 1) * PulsePeriod_us
%
% Then for the trial:
%   lastActivePTD_us(tr) = max(final artifact time across active rows)

lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials

    if conditionType_trial(tr) == 0
        lastActivePTD_us(tr) = 0;
        continue;
    end

    if conditionType_trial(tr) == 1
        nActive_this = prefixLength_trial(tr);
    elseif conditionType_trial(tr) == 2
        nActive_this = activeCount_trial(tr);
    else
        nActive_this = activeCount_trial(tr);
    end

    if isnan(nActive_this) || nActive_this < 1
        lastActivePTD_us(tr) = 0;
        continue;
    end

    nActive_this = min(round(nActive_this), simultaneous_stim);

    rr = (tr-1)*simultaneous_stim + (1:simultaneous_stim);
    activeRows = rr(1:nActive_this);

    % Column 6: post-trigger delay, us.
    ptd_us = cell2mat(StimParams_data(activeRows,6));

    % Column 8: pulse train number.
    pulseNum = cell2mat(StimParams_data(activeRows,8));

    % Column 9: pulse train period, us.
    pulsePeriod_us = cell2mat(StimParams_data(activeRows,9));

    pulseNum(isnan(pulseNum) | pulseNum < 1) = 1;
    pulsePeriod_us(isnan(pulsePeriod_us)) = 0;

    rowFinalArtifact_us = ptd_us + (pulseNum - 1) .* pulsePeriod_us;

    lastActivePTD_us(tr) = max(rowFinalArtifact_us);
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

%% ====================== DEBUG TRIAL CONTENT CHECK ======================
if debug_print_trial_content

    isPrefixTrial_tmp = conditionType_trial == 1;

    prefix_all = unique(prefixLength_trial(isPrefixTrial_tmp & prefixLength_trial > 0));
    set_all = unique(conditionSetID_trial(conditionSetID_trial > 0));
    amp_all = unique(trialAmps(isPrefixTrial_tmp));
    amp_all = amp_all(amp_all > 0);
    isi_all = unique(isi_ms_trial(isPrefixTrial_tmp));

    if ~isempty(prefix_all) && ~isempty(set_all) && ~isempty(amp_all) && ~isempty(isi_all)

        debug_set = set_all(1);
        debug_amp = amp_all(end);
        debug_isi = isi_all(1);

        fprintf('\n================ DEBUG TRIAL CONTENT CHECK ================\n');
        fprintf('Debug Set = %d | Amp = %.1f uA | ISI = %.1f ms\n', ...
            debug_set, debug_amp, debug_isi);

        for ip = 1:numel(prefix_all)

            prefix_val = prefix_all(ip);

            tr_debug = find(conditionSetID_trial == debug_set & ...
                            conditionType_trial == 1 & ...
                            prefixLength_trial == prefix_val & ...
                            isi_ms_trial == debug_isi & ...
                            trialAmps == debug_amp, ...
                            1, 'first');

            if isempty(tr_debug)
                fprintf('\nPrefix %d: no matching trial found.\n', prefix_val);
                continue;
            end

            rr = (tr_debug-1)*simultaneous_stim + (1:simultaneous_stim);

            stimNames_debug = StimParams_data(rr,1);
            ptd_debug = cell2mat(StimParams_data(rr,6));
            pulseNum_debug = cell2mat(StimParams_data(rr,8));
            pulsePeriod_debug = cell2mat(StimParams_data(rr,9));
            amp_debug = cell2mat(StimParams_data(rr,16));
            activeCount_debug = cell2mat(StimParams_data(rr,26));
            prefix_debug = cell2mat(StimParams_data(rr,27));
            isi_debug = cell2mat(StimParams_data(rr,28));
            condType_debug = cell2mat(StimParams_data(rr,29));
            setID_debug = cell2mat(StimParams_data(rr,30));

            fprintf('\nPrefix %d | Trial %d\n', prefix_val, tr_debug);
            fprintf('  conditionType = %d, setID = %d, ISI = %.1f ms, amp = %.1f uA\n', ...
                conditionType_trial(tr_debug), conditionSetID_trial(tr_debug), ...
                isi_ms_trial(tr_debug), trialAmps(tr_debug));

            fprintf('  StimNames:       ');
            disp(stimNames_debug');

            fprintf('  PTD_us:          ');
            disp(ptd_debug');

            fprintf('  PulseNum:        ');
            disp(pulseNum_debug');

            fprintf('  PulsePeriod_us:  ');
            disp(pulsePeriod_debug');

            fprintf('  Amp_col16:       ');
            disp(amp_debug');

            fprintf('  ActiveCount_col26: ');
            disp(activeCount_debug');

            fprintf('  Prefix_col27:      ');
            disp(prefix_debug');

            fprintf('  ISI_col28:         ');
            disp(isi_debug');

            fprintf('  CondType_col29:    ');
            disp(condType_debug');

            fprintf('  SetID_col30:       ');
            disp(setID_debug');

            fprintf('  FinalArtifact_us = %.1f\n', lastActivePTD_us(tr_debug));
        end

        fprintf('===========================================================\n');
    end
end

%% ====================== LOAD FIRST-SPIKE TIMES ======================
fst_file = dir('FirstSpikeTimes_Prefix.mat');
assert(~isempty(fst_file), 'Missing FirstSpikeTimes_Prefix.mat. Run the first-spike extraction code first.');

FST = load(fst_file(1).name, ...
    'firstSpikeTimes', ...
    'hasSpike', ...
    'lastActivePTD_ms', ...
    'activeCount_trial', ...
    'prefixLength_trial', ...
    'isi_ms_trial', ...
    'conditionType_trial', ...
    'conditionSetID_trial', ...
    'trialAmps');

firstSpikeTimes = FST.firstSpikeTimes;

fprintf('Loaded first-spike file: %s\n', fst_file(1).name);

%% ====================== FIRST-SPIKE ALIGNMENT CHECK ======================
assert(numel(FST.prefixLength_trial) == n_Trials, ...
    'FirstSpikeTimes file does not match n_Trials.');

if any(FST.prefixLength_trial(:) ~= prefixLength_trial(:))
    error('FirstSpikeTimes prefixLength_trial does not match current StimParams metadata.');
end

if any(FST.conditionType_trial(:) ~= conditionType_trial(:))
    error('FirstSpikeTimes conditionType_trial does not match current StimParams metadata.');
end

if any(FST.conditionSetID_trial(:) ~= conditionSetID_trial(:))
    error('FirstSpikeTimes conditionSetID_trial does not match current StimParams metadata.');
end

if any(FST.trialAmps(:) ~= trialAmps(:))
    error('FirstSpikeTimes trialAmps does not match current StimParams metadata.');
end

if isfield(FST, 'lastActivePTD_ms')
    max_diff_ptd = max(abs(FST.lastActivePTD_ms(:) - lastActivePTD_ms(:)));
    if max_diff_ptd > 1e-6
        warning('FirstSpikeTimes lastActivePTD_ms differs from current calculation. Using current calculation.');
    end
end

%% ====================== DEFINE GROUPS TO RECOVER ======================
isPrefixTrial = conditionType_trial == 1;

all_prefixes = unique(prefixLength_trial(isPrefixTrial & prefixLength_trial > 0));
all_prefixes = sort(all_prefixes(:))';

if isempty(target_prefixes_to_recover)
    % Prefix 1 has no shorter source prefix, so start from Prefix 2.
    target_prefixes = all_prefixes(all_prefixes >= 2);
else
    target_prefixes = intersect(all_prefixes, target_prefixes_to_recover);
    target_prefixes = target_prefixes(target_prefixes >= 2);
end

% ================= XIA MODIFICATION =================
% Important for cascading recovery.
% Prefix 2 must be recovered before Prefix 3,
% Prefix 3 before Prefix 4, and so on.
target_prefixes = sort(target_prefixes(:))';
% =====================================================

all_amps = unique(trialAmps(isPrefixTrial));
all_amps = all_amps(all_amps > 0);   % exclude zero-control

if isempty(amps_to_recover)
    amps_sel = all_amps;
else
    amps_sel = intersect(all_amps, amps_to_recover);
end

all_sets = unique(conditionSetID_trial(isPrefixTrial & conditionSetID_trial > 0));

if isempty(sets_to_recover)
    sets_sel = all_sets;
else
    sets_sel = intersect(all_sets, sets_to_recover);
end

all_isis = unique(isi_ms_trial(isPrefixTrial));

if isempty(isis_to_recover_ms)
    isis_sel = all_isis;
else
    isis_sel = intersect(all_isis, isis_to_recover_ms);
end

fprintf('\nTarget prefixes to recover: ');
disp(target_prefixes);

fprintf('Amplitudes to recover: ');
disp(amps_sel');

fprintf('Set IDs to recover: ');
disp(sets_sel');

fprintf('ISIs to recover (ms): ');
disp(isis_sel');

fprintf('Detected final active artifact times (ms): ');
disp(unique(lastActivePTD_ms)');

%% ====================== RECOVERY MAIN LOOP ======================
total_spikes_added_all = 0;
total_spikes_deleted_all = 0;

RecoveryLog = struct();
log_i = 0;

fprintf('\nStarting prefix-based CASCADING recovery...\n');

for ai = 1:numel(amps_sel)

    amp_val = amps_sel(ai);

    total_added_amp = 0;
    total_deleted_amp = 0;

    for si = 1:numel(sets_sel)

        set_id = sets_sel(si);

        for ii = 1:numel(isis_sel)

            isi_val = isis_sel(ii);

            % IMPORTANT:
            % target_prefixes loop must stay inside amp/set/ISI loop and
            % must be sorted ascending for cascading recovery.
            for pi = 1:numel(target_prefixes)

                target_prefix = target_prefixes(pi);
                source_prefix = target_prefix - 1;

                %% ----- Find source and target trials -----
                source_trials = find(conditionType_trial == 1 & ...
                                     prefixLength_trial == source_prefix & ...
                                     conditionSetID_trial == set_id & ...
                                     isi_ms_trial == isi_val & ...
                                     trialAmps == amp_val);

                target_trials = find(conditionType_trial == 1 & ...
                                     prefixLength_trial == target_prefix & ...
                                     conditionSetID_trial == set_id & ...
                                     isi_ms_trial == isi_val & ...
                                     trialAmps == amp_val);

                % Keep only trials with triggers.
                source_trials = source_trials(source_trials <= nTrials_use);
                target_trials = target_trials(target_trials <= nTrials_use);

                if isempty(source_trials) || isempty(target_trials)
                    continue;
                end

                fprintf('\nAmp %g uA | Set %d | ISI %.1f ms | Recover Prefix %d using Prefix %d\n', ...
                    amp_val, set_id, isi_val, target_prefix, source_prefix);

                fprintf('  Source trials: %d | Target trials: %d\n', ...
                    numel(source_trials), numel(target_trials));

                group_added = 0;
                group_deleted = 0;

                %% ----- Loop target trials -----
                for gi = 1:numel(target_trials)

                    tr_target = target_trials(gi);

                    % Select matching source trial.
                    if use_deterministic_cycling
                        idx_source = mod(gi-1, numel(source_trials)) + 1;
                    else
                        idx_source = randi(numel(source_trials));
                    end

                    tr_source = source_trials(idx_source);

                    t0_target_ms = trig(tr_target) / FS * 1000;
                    t0_source_ms = trig(tr_source) / FS * 1000;

                    %% ----- Per recording channel recovery -----
                    for rec_ch = 1:nChn

                        % ================= XIA MODIFICATION =================
                        % CRITICAL FIX:
                        % Use sp_seq as the source, not sp_clipped.
                        %
                        % This allows cascading recovery:
                        %   Prefix 3 can use already-recovered Prefix 2.
                        %   Prefix 4 can use already-recovered Prefix 3.
                        %   Prefix 5 can use already-recovered Prefix 4.
                        %
                        % If sp_clipped is used here, only Prefix 2 can be
                        % properly recovered. Prefix 3/4/5 will not inherit
                        % the early spikes added in previous prefix steps.
                        spikes_source = sp_seq{rec_ch};
                        % =====================================================

                        if isempty(spikes_source)
                            continue;
                        end

                        %% ----- Determine recovery window end -----
                        first_ms = NaN;

                        if rec_ch <= numel(firstSpikeTimes) && ...
                           tr_target <= numel(firstSpikeTimes{rec_ch})
                            first_ms = firstSpikeTimes{rec_ch}(tr_target);
                        end

                        lastPTD_ms = lastActivePTD_ms(tr_target);

                        if isfinite(first_ms) && first_ms > win_start_ms
                            win_end = first_ms;
                        else
                            if ~use_fallback
                                continue;
                            end
                            win_end = lastPTD_ms + fallback_ms;
                        end

                        % Safety cap to prevent over-recovery.
                        win_end = min(win_end, lastPTD_ms + max_extra_ms);

                        if win_end <= win_start_ms
                            continue;
                        end

                        if debug_recovery_window && rec_ch == debug_channel && gi <= debug_n_trials_per_group
                            fprintf(['    Debug Ch %d | Target trial %d | Source trial %d | ' ...
                                     'Prefix %d<-Prefix %d | first_ms %.2f | lastPTD %.2f | win_end %.2f\n'], ...
                                    rec_ch, tr_target, tr_source, target_prefix, source_prefix, ...
                                    first_ms, lastPTD_ms, win_end);
                        end

                        %% ----- Delete target spikes in recovery window -----
                        seq_spikes_ch = sp_seq{rec_ch};
                        n_deleted = 0;

                        if ~isempty(seq_spikes_ch)
                            rel_target = seq_spikes_ch(:,1) - t0_target_ms;

                            del_mask = rel_target >= win_start_ms & ...
                                       rel_target <  win_end;

                            n_deleted = sum(del_mask);

                            if n_deleted > 0
                                seq_spikes_ch(del_mask,:) = [];
                                sp_seq{rec_ch} = seq_spikes_ch;
                            end
                        end

                        %% ----- Copy source spikes from same relative window -----
                        % Important:
                        % spikes_source was captured before target deletion.
                        % This is okay because source and target trials are
                        % different prefix groups.
                        rel_source = spikes_source(:,1) - t0_source_ms;

                        in_win = rel_source >= win_start_ms & ...
                                 rel_source <  win_end;

                        spikes_add = spikes_source(in_win,:);

                        if isempty(spikes_add)
                            total_spikes_deleted_all = total_spikes_deleted_all + n_deleted;
                            total_deleted_amp = total_deleted_amp + n_deleted;
                            group_deleted = group_deleted + n_deleted;
                            continue;
                        end

                        % Shift source spikes into target trial time.
                        spikes_add(:,1) = rel_source(in_win) + t0_target_ms;

                        % Append and sort.
                        sp_seq{rec_ch} = sortrows([sp_seq{rec_ch}; spikes_add], 1);

                        %% ----- Update counters -----
                        n_added = size(spikes_add,1);

                        total_spikes_added_all = total_spikes_added_all + n_added;
                        total_spikes_deleted_all = total_spikes_deleted_all + n_deleted;

                        total_added_amp = total_added_amp + n_added;
                        total_deleted_amp = total_deleted_amp + n_deleted;

                        group_added = group_added + n_added;
                        group_deleted = group_deleted + n_deleted;

                    end % rec_ch
                end % target trial

                %% ----- Save group-level log -----
                log_i = log_i + 1;
                RecoveryLog(log_i).amp = amp_val;
                RecoveryLog(log_i).set_id = set_id;
                RecoveryLog(log_i).isi_ms = isi_val;
                RecoveryLog(log_i).target_prefix = target_prefix;
                RecoveryLog(log_i).source_prefix = source_prefix;
                RecoveryLog(log_i).n_source_trials = numel(source_trials);
                RecoveryLog(log_i).n_target_trials = numel(target_trials);
                RecoveryLog(log_i).n_added = group_added;
                RecoveryLog(log_i).n_deleted = group_deleted;

                fprintf('  Group added %d spikes, deleted %d spikes.\n', ...
                    group_added, group_deleted);

            end % target prefix
        end % ISI
    end % set

    fprintf('\nAmplitude %g uA: added %d spikes, deleted %d spikes.\n', ...
        amp_val, total_added_amp, total_deleted_amp);
end

fprintf('\n==============================================================\n');
fprintf('Total spikes added:   %d\n', total_spikes_added_all);
fprintf('Total spikes deleted: %d\n', total_spikes_deleted_all);
fprintf('==============================================================\n');

%% ====================== SAVE RECOVERED SPIKES ======================
base_sp_name = seq_sp_file(1).name;

% Make output name.
% Example:
%   Xia_xxx.sp_xia.mat
% becomes:
%   Xia_xxx.sp_xia_PrefixRecovery.mat
if contains(base_sp_name, '.sp_xia.mat')
    new_name = strrep(base_sp_name, '.sp_xia.mat', '.sp_xia_PrefixRecovery.mat');
else
    new_name = [erase(base_sp_name, '.mat') '_PrefixRecovery.mat'];
end

save(fullfile(data_folder, new_name), ...
     'sp_seq', ...
     'RecoveryLog', ...
     'total_spikes_added_all', ...
     'total_spikes_deleted_all', ...
     'win_start_ms', ...
     'fallback_ms', ...
     'max_extra_ms', ...
     'use_fallback', ...
     'target_prefixes', ...
     'amps_sel', ...
     'sets_sel', ...
     'isis_sel', ...
     'FS', ...
     'lastActivePTD_ms', ...
     'lastActivePTD_us', ...
     'activeCount_trial', ...
     'prefixLength_trial', ...
     'isi_ms_trial', ...
     'conditionType_trial', ...
     'conditionSetID_trial', ...
     'trialAmps', ...
     'n_Trials', ...
     'nTrials_use', ...
     '-v7.3');

fprintf('Saved prefix-recovered spike file:\n%s\n', fullfile(data_folder, new_name));