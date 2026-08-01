function result = qc_load_bad_channels(file_record)
%QC_LOAD_BAD_CHANNELS Load BadCh_perSet or BadCh without changing indices.

result = struct('file_found',false, 'source_file','', ...
    'variable','', 'per_set',{{}});
if ~file_record.found
    return;
end

S = load(file_record.path);
if isfield(S, 'BadCh_perSet')
    value = S.BadCh_perSet;
    variable = 'BadCh_perSet';
elseif isfield(S, 'BadCh')
    value = S.BadCh;
    variable = 'BadCh';
else
    error('StageA:MissingBadChannelVariable', ...
        'No BadCh_perSet or BadCh variable in %s', file_record.path);
end

if iscell(value)
    per_set = value;
else
    per_set = {value};
end

for k = 1:numel(per_set)
    if isempty(per_set{k})
        per_set{k} = [];
    elseif isnumeric(per_set{k})
        per_set{k} = unique(double(per_set{k}(:).'));
    else
        error('StageA:InvalidBadChannelEntry', ...
            'Bad-channel entry %d is not numeric in %s', k, file_record.path);
    end
end

result.file_found = true;
result.source_file = file_record.path;
result.variable = variable;
result.per_set = per_set;
end

