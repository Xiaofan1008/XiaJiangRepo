function cross = qc_validate_cross_dataset(DatasetInfo, CFG)
%QC_VALIDATE_CROSS_DATASET Check mappings required by later stages.

cross = struct();
cross.messages = {};
cross.validation_passed = true;

roles = {'single','sim','seq'};
channel_counts = cellfun(@(r) DatasetInfo.(r).spike.n_channels, roles);
cross.channel_count_consistent = numel(unique(channel_counts)) == 1;
if ~cross.channel_count_consistent
    cross = add_problem(cross, sprintf( ...
        'sp_corr channel counts differ: Single=%d, Sim=%d, Seq=%d.', ...
        channel_counts(1), channel_counts(2), channel_counts(3)));
end

fs_values = [DatasetInfo.single.trigger.FS, DatasetInfo.sim.trigger.FS, ...
    DatasetInfo.seq.trigger.FS];
cross.fs_consistent = numel(unique(fs_values)) == 1;
if ~cross.fs_consistent
    cross = add_problem(cross, 'Sampling rates differ across datasets.');
end

single_names = nonempty_names(DatasetInfo.single.experiment.electrode_names);
sim_names = nonempty_names(DatasetInfo.sim.experiment.electrode_names);
seq_names = nonempty_names(DatasetInfo.seq.experiment.electrode_names);

cross.single_to_sim_name_coverage = coverage(single_names, sim_names);
cross.single_to_seq_name_coverage = coverage(single_names, seq_names);

if cross.single_to_sim_name_coverage < 1
    cross.messages{end+1,1} = sprintf( ...
        '[WARNING][CROSS] Only %.1f%% of single electrode names occur in Sim.', ...
        100*cross.single_to_sim_name_coverage);
end
if cross.single_to_seq_name_coverage < 1
    cross.messages{end+1,1} = sprintf( ...
        '[WARNING][CROSS] Only %.1f%% of single electrode names occur in Seq.', ...
        100*cross.single_to_seq_name_coverage);
end

cross.sim_target_ptd_available = any(abs( ...
    DatasetInfo.sim.experiment.ptds_ms - CFG.sim_PTD_ms) ...
    <= CFG.PTD_tolerance_ms);
cross.seq_target_ptd_available = any(abs( ...
    DatasetInfo.seq.experiment.ptds_ms - CFG.seq_PTD_ms) ...
    <= CFG.PTD_tolerance_ms);

if ~cross.sim_target_ptd_available
    cross = add_problem(cross, sprintf( ...
        'Sim target PTD %.3f ms is unavailable.', CFG.sim_PTD_ms));
end
if ~cross.seq_target_ptd_available
    cross = add_problem(cross, sprintf( ...
        'Seq target PTD %.3f ms is unavailable.', CFG.seq_PTD_ms));
end
end

function names = nonempty_names(names)
names = names(~cellfun(@isempty, names));
names = unique(names);
end

function value = coverage(source, target)
if isempty(source)
    value = 0;
else
    value = mean(ismember(source, target));
end
end

function cross = add_problem(cross, message)
cross.messages{end+1,1} = ['[ERROR][CROSS] ' message];
cross.validation_passed = false;
end

