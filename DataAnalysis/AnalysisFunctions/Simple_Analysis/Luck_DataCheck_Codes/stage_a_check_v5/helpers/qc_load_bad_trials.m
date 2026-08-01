function result = qc_load_bad_trials(file_record)
%QC_LOAD_BAD_TRIALS Load and verify a globally repeated BadTrials cell array.

result = struct('file_found',false, 'source_file','', 'raw',{{}}, ...
    'is_global',true, 'global_trials',[], 'nonempty_entries',0);

if ~file_record.found
    return;
end

S = load(file_record.path);
if ~isfield(S, 'BadTrials')
    error('StageA:MissingBadTrialsVariable', ...
        'BadTrials variable missing from %s', file_record.path);
end
if ~iscell(S.BadTrials)
    error('StageA:InvalidBadTrials', ...
        'BadTrials must be a cell array in %s', file_record.path);
end

result.file_found = true;
result.source_file = file_record.path;
result.raw = S.BadTrials;

normalized = cell(size(S.BadTrials));
for k = 1:numel(S.BadTrials)
    value = S.BadTrials{k};
    if isempty(value)
        normalized{k} = [];
    elseif isnumeric(value)
        normalized{k} = unique(double(value(:).'));
    else
        error('StageA:InvalidBadTrialsEntry', ...
            'BadTrials{%d} is not numeric in %s', k, file_record.path);
    end
end

result.nonempty_entries = sum(~cellfun(@isempty, normalized));
if isempty(normalized)
    result.global_trials = [];
    return;
end

reference = normalized{1};
result.is_global = all(cellfun(@(x) isequal(x, reference), normalized));
if result.is_global
    result.global_trials = reference;
else
    result.global_trials = [];
end
end

