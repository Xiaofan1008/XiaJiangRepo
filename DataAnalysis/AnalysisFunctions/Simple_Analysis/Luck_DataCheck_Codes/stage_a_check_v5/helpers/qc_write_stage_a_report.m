function qc_write_stage_a_report(filename, StageA, print_to_command_window)
%QC_WRITE_STAGE_A_REPORT Write and optionally print the validation report.

if nargin < 3
    print_to_command_window = false;
end

fid = fopen(filename, 'w');
if fid < 0
    error('StageA:ReportWriteFailed', 'Could not write report: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

write_report(fid, StageA);

if print_to_command_window
    fprintf('\n\n');
    fprintf('################ COMPLETE STAGE A REPORT ################\n');
    write_report(1, StageA); % File identifier 1 is the MATLAB Command Window.
    fprintf('############## END COMPLETE STAGE A REPORT ##############\n\n');
end
end

function write_report(fid, StageA)
fprintf(fid, 'Stage A Dataset Validation Report\n');
fprintf(fid, 'Created: %s\n', char(string(StageA.created_on)));
fprintf(fid, 'Overall passed: %s\n\n', ...
    char(string(StageA.validation_passed)));

roles = {'single','sim','seq'};
for k = 1:numel(roles)
    role = roles{k};
    fprintf(fid, '[%s]\n', upper(role));
    if ~isfield(StageA.datasets, role)
        fprintf(fid, 'Dataset result missing.\n\n');
        continue;
    end
    D = StageA.datasets.(role);
    fprintf(fid, 'Folder: %s\n', D.folder);
    fprintf(fid, 'Passed: %s\n', char(string(D.validation_passed)));
    if isfield(D, 'fatal_error')
        fprintf(fid, 'Fatal error: %s\n\n', D.fatal_error);
        continue;
    end
    fprintf(fid, 'Spike file: %s\n', D.files.spike.path);
    fprintf(fid, 'Experiment file: %s\n', D.files.experiment.path);
    fprintf(fid, 'Trigger file: %s\n', D.files.trigger.path);
    fprintf(fid, 'Bad-trial file: %s\n', value_or_none(D.files.bad_trials));
    fprintf(fid, 'Bad-channel file: %s\n', value_or_none(D.files.bad_channels));
    fprintf(fid, 'Responding file: %s\n', value_or_none(D.files.responding));
    fprintf(fid, 'Trials: %d\n', D.experiment.n_trials);
    fprintf(fid, 'Triggers: %d\n', D.trigger.n_triggers);
    fprintf(fid, 'sp_corr channels: %d\n', D.spike.n_channels);
    fprintf(fid, 'Stimulation events/trial: %d\n', ...
        D.experiment.simultaneous_stim);
    fprintf(fid, 'Amplitudes (uA): %s\n', num2str(D.experiment.amplitudes(:).'));
    fprintf(fid, 'PTDs (ms): %s\n', num2str(D.experiment.ptds_ms(:).'));

    if D.bad_trials.file_found
        fprintf(fid, 'Global bad-trial consistency: %s\n', ...
            char(string(D.bad_trials.is_global)));
        fprintf(fid, 'Global bad-trial count: %d\n', ...
            numel(D.bad_trials.global_trials));
        fprintf(fid, 'Global bad-trial IDs: %s\n', ...
            vector_or_none(D.bad_trials.global_trials));
    else
        fprintf(fid, 'Global bad trials: (not available)\n');
    end

    if D.bad_channels.file_found
        fprintf(fid, 'Bad-channel sets: %d\n', ...
            numel(D.bad_channels.per_set));
    else
        fprintf(fid, 'Bad-channel sets: (not available)\n');
    end

    fprintf(fid, 'Responding results available: %s\n\n', ...
        char(string(D.responding.file_found)));
end

if isfield(StageA, 'cross_dataset') && ...
        isfield(StageA.cross_dataset, 'validation_passed')
    C = StageA.cross_dataset;
    fprintf(fid, '[CROSS-DATASET CHECKS]\n');
    fprintf(fid, 'Passed: %s\n', char(string(C.validation_passed)));
    fprintf(fid, 'Sampling rates consistent: %s\n', ...
        char(string(C.fs_consistent)));
    fprintf(fid, 'Channel counts consistent: %s\n', ...
        char(string(C.channel_count_consistent)));
    fprintf(fid, 'Single-to-Sim electrode-name coverage: %.1f%%\n', ...
        100*C.single_to_sim_name_coverage);
    fprintf(fid, 'Single-to-Seq electrode-name coverage: %.1f%%\n', ...
        100*C.single_to_seq_name_coverage);
    fprintf(fid, 'Target Sim PTD available: %s\n', ...
        char(string(C.sim_target_ptd_available)));
    fprintf(fid, 'Target Seq PTD available: %s\n\n', ...
        char(string(C.seq_target_ptd_available)));
end

fprintf(fid, '[MESSAGES]\n');
if isempty(StageA.messages)
    fprintf(fid, 'None\n');
else
    for k = 1:numel(StageA.messages)
        fprintf(fid, '%s\n', StageA.messages{k});
    end
end
end

function value = value_or_none(record)
if record.found
    value = record.path;
else
    value = '(not found)';
end
end

function value = vector_or_none(values)
if isempty(values)
    value = '(none)';
else
    value = num2str(values(:).');
end
end
