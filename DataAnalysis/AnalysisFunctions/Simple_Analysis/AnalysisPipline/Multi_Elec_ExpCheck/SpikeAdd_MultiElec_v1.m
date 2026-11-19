% ============================================================
%  Prefix-Based Spike Recovery for Mixed-Prefix Stimulation
%
%  Purpose:
%    Recover spikes removed by artifact blanking in longer prefix trials.
%
%  Recovery rule:
%    Target Prefix N is recovered using Source Prefix N-1.
%
%  Example:
%    Prefix 2 target uses Prefix 1 source
%    Prefix 3 target uses Prefix 2 source
%    Prefix 4 target uses Prefix 3 source
%    Prefix 5 target uses Prefix 4 source
%
%  Matching:
%    source and target trials must have:
%      - same amplitude
%      - same condition set ID
%      - same ISI
%      - conditionType == 1
%
%  Output:
%    File:
%      *.sp_xia_PrefixRecovery.mat
%
%    Variable:
%      sp_seq
% ============================================================

clear all;
addpath(genpath('/Volumes/MACData/Data/Data_Xia/AnalysisFunctions'));

%% ====================== USER PARAMETERS ======================
data_folder = '/Volumes/MACData/Data/Data_Xia/DX026/Xia_Ele5_SimSeq5Pulse1_260602_182126';

FS = 30000;

% Recovery window starts relative to the first trigger.
% Usually 0 ms because the blanking starts around the first artifact.
win_start_ms = 0;

% If no first post-blank spike is detected, use:
%   win_end = lastActivePTD_ms + fallback_ms
use_fallback = true;
fallback_ms  = 2;

% Safety cap:
% Even if the first detected spike is very late, do not recover beyond:
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

% Source trial selection:
% deterministic cycling means source trials are reused in order.
% This is similar to your previous method.
use_deterministic_cycling = true;

%% ====================== LOAD FOLDER ======================
cd(data_folder);

%% ====================== LOAD TARGET SPIKES ======================
% Load the original spike file that will be recovered.
seq_sp_file = dir('*sp_xia.mat');
assert(~isempty(seq_sp_file), 'No *sp_xia.mat file found.');

load(seq_sp_file(1).name, 'sp_clipped');
sp_seq = sp_clipped;   % output variable name required by later analysis

nChn = numel(sp_seq);

fprintf('Loaded spike file: %s\n', seq_sp_file(1).name);

%% ====================== LOAD TRIGGERS ======================
trig = loadTrig(0);

%% ====================== LOAD EXPERIMENT PARAMETERS ======================
param_file = dir('*_exp_datafile_*.mat');
assert(~isempty(param_file), 'No *_exp_datafile_*.mat file found.');

S = load(param_file(1).name);

StimParams        = S.StimParams;
simultaneous_stim = S.simultaneous_stim;   % rows/slots per trial
n_Trials          = S.n_Trials;

%% ====================== CHECK REQUIRED METADATA ======================
requiredVars = {'active_electrode_count_by_trial', ...
                'prefix_length_by_trial', ...
                'isi_ms_by_trial', ...
                'condition_type_by_trial', ...
                'condition_set_id_by_trial'};

for i = 1:numel(requiredVars)
    assert(isfield(S, requiredVars{i}), ...
        'Missing metadata variable "%s". This code requires the new mixed-prefix file.', ...
        requiredVars{i});
end

activeCount_trial    = S.active_electrode_count_by_trial(:);
prefixLength_trial   = S.prefix_length_by_trial(:);
isi_ms_trial         = S.isi_ms_by_trial(:);
conditionType_trial  = S.condition_type_by_trial(:);
conditionSetID_trial = S.condition_set_id_by_trial(:);

%% ====================== EXTRACT AMPLITUDE PER TRIAL ======================
trialAmps_all = cell2mat(StimParams(2:end,16));
trialAmps = trialAmps_all(1:simultaneous_stim:end);
trialAmps(trialAmps == -1) = 0;

%% ====================== CALCULATE LAST ACTIVE PTD ======================
% This should match what was saved by Code 1.
lastActivePTD_us = zeros(n_Trials,1);

for tr = 1:n_Trials
    activeCount_this = activeCount_trial(tr);

    if isnan(activeCount_this) || activeCount_this < 1
        lastActivePTD_us(tr) = 0;
    else
        activeCount_this = min(round(activeCount_this), simultaneous_stim);
        stimRow = 1 + (tr-1)*simultaneous_stim + activeCount_this; % +1 for header
        lastActivePTD_us(tr) = StimParams{stimRow,6};
    end
end

lastActivePTD_ms = lastActivePTD_us ./ 1000;

%% ====================== LOAD FIRST-SPIKE TIMES ======================
fst_file = dir('FirstSpikeTimes_Prefix.mat');
assert(~isempty(fst_file), 'Missing FirstSpikeTimes_Prefix.mat. Run Code 1 first.');

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

% Basic alignment check.
assert(numel(FST.prefixLength_trial) == n_Trials, ...
    'FirstSpikeTimes file does not match n_Trials.');

fprintf('Loaded first-spike file: %s\n', fst_file(1).name);

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

%% ====================== RECOVERY MAIN LOOP ======================
total_spikes_added_all = 0;
total_spikes_deleted_all = 0;

% Keep a log so you can check what happened.
RecoveryLog = struct();
log_i = 0;

fprintf('\nStarting prefix-based recovery...\n');

for ai = 1:numel(amps_sel)
    amp_val = amps_sel(ai);

    total_added_amp = 0;
    total_deleted_amp = 0;

    for si = 1:numel(sets_sel)
        set_id = sets_sel(si);

        for ii = 1:numel(isis_sel)
            isi_val = isis_sel(ii);

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

                if isempty(source_trials) || isempty(target_trials)
                    continue;
                end

                fprintf('\nAmp %g µA | Set %d | ISI %.1f ms | Recover Prefix %d using Prefix %d\n', ...
                    amp_val, set_id, isi_val, target_prefix, source_prefix);
                fprintf('  Source trials: %d | Target trials: %d\n', ...
                    numel(source_trials), numel(target_trials));

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

                        spikes_source = sp_clipped{rec_ch};
                        if isempty(spikes_source)
                            continue;
                        end

                        %% ----- Determine recovery window end -----
                        % First choice:
                        %   first detected post-blank spike in target trial/channel.
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

                        %% ----- Delete target spikes in recovery window -----
                        % This prevents duplicate spikes after injection.
                        seq_spikes_ch = sp_seq{rec_ch};
                        n_deleted = 0;

                        if ~isempty(seq_spikes_ch)
                            rel_target = seq_spikes_ch(:,1) - t0_target_ms;
                            del_mask = rel_target >= win_start_ms & rel_target < win_end;

                            n_deleted = sum(del_mask);

                            if n_deleted > 0
                                seq_spikes_ch(del_mask,:) = [];
                                sp_seq{rec_ch} = seq_spikes_ch;
                            end
                        end

                        %% ----- Copy source spikes from same relative window -----
                        rel_source = spikes_source(:,1) - t0_source_ms;
                        in_win = rel_source >= win_start_ms & rel_source < win_end;

                        spikes_add = spikes_source(in_win,:);

                        if isempty(spikes_add)
                            total_spikes_deleted_all = total_spikes_deleted_all + n_deleted;
                            total_deleted_amp = total_deleted_amp + n_deleted;
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

            end % target prefix
        end % ISI
    end % set

    fprintf('\nAmplitude %g µA: added %d spikes, deleted %d spikes.\n', ...
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
     '-v7.3');

fprintf('Saved prefix-recovered spike file:\n%s\n', fullfile(data_folder, new_name));