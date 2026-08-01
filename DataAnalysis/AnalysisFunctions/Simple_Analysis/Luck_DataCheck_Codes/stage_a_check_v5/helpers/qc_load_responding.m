function result = qc_load_responding(file_record)
%QC_LOAD_RESPONDING Load the final saved Responding structure.

result = struct('file_found',false, 'source_file','', 'Responding',[]);
if ~file_record.found
    return;
end

S = load(file_record.path);
if ~isfield(S, 'Responding')
    error('StageA:MissingRespondingVariable', ...
        'Responding variable missing from %s', file_record.path);
end

result.file_found = true;
result.source_file = file_record.path;
result.Responding = S.Responding;

optional_names = {'Detection_Mode','baseline_win_ms','post_win_ms', ...
    'post_win_singlepulse_ms','k_SD','Amps','PTDs'};
for k = 1:numel(optional_names)
    name = optional_names{k};
    if isfield(S, name)
        result.(name) = S.(name);
    end
end
end

