function experiment = qc_decode_stim_params(E)
%QC_DECODE_STIM_PARAMS Decode trial stimulation identities and parameters.

StimParams = E.StimParams;
simN = double(E.simultaneous_stim);
nTrials = double(E.n_Trials);

if ~iscell(StimParams) || size(StimParams,1) < 2
    error('StageA:InvalidStimParams', 'StimParams must be a nonempty cell array.');
end
if ~isscalar(simN) || simN < 1 || mod(simN,1) ~= 0
    error('StageA:InvalidSimultaneousStim', ...
        'simultaneous_stim must be a positive integer.');
end
if ~isscalar(nTrials) || nTrials < 1 || mod(nTrials,1) ~= 0
    error('StageA:InvalidTrialCount', 'n_Trials must be a positive integer.');
end
if size(StimParams,2) < 16
    error('StageA:InvalidStimParams', ...
        'StimParams has fewer than the required 16 columns.');
end

required_data_rows = nTrials * simN;
if size(StimParams,1) - 1 < required_data_rows
    error('StageA:InvalidStimParams', ...
        'StimParams does not contain enough rows for all trials.');
end

stim_names = cell(nTrials, simN);
stim_indices = zeros(nTrials, simN);
amp_per_event = nan(nTrials, simN);
ptd_us_per_event = nan(nTrials, simN);
trial_rows = cell(nTrials,1);

map_names = normalize_names(E.E_MAP(2:end));

for tr = 1:nTrials
    rows = (tr-1)*simN + (2:(simN+1));
    trial_rows{tr} = rows;

    names = normalize_names(StimParams(rows,1));
    stim_names(tr,:) = names(:).';
    [tf, idx] = ismember(names, map_names);
    idx(~tf) = 0;
    stim_indices(tr,:) = idx(:).';

    amp_per_event(tr,:) = numeric_cells(StimParams(rows,16), 'amplitude', tr);
    ptd_us_per_event(tr,:) = numeric_cells(StimParams(rows,6), 'PTD', tr);
end

trial_amplitudes = amp_per_event(:,1);
trial_amplitudes(trial_amplitudes == -1) = 0;

if simN > 1
    trial_ptd_us = ptd_us_per_event(:,2);
else
    trial_ptd_us = zeros(nTrials,1);
end
trial_ptd_ms = trial_ptd_us / 1000;

[unique_sets,~,set_index] = unique(stim_indices, 'rows', 'stable');
[amplitudes,~,amplitude_index] = unique(trial_amplitudes);
[ptds_ms,~,ptd_index] = unique(trial_ptd_ms);

experiment = struct();
experiment.n_trials = nTrials;
experiment.simultaneous_stim = simN;
experiment.E_MAP = E.E_MAP;
experiment.electrode_names = map_names;
experiment.StimParams = StimParams;
experiment.trial_stimparam_rows = trial_rows;
experiment.stim_names_per_trial = stim_names;
experiment.stim_indices_per_trial = stim_indices;
experiment.amp_per_event = amp_per_event;
experiment.ptd_us_per_event = ptd_us_per_event;
experiment.trial_amplitude = trial_amplitudes;
experiment.trial_ptd_ms = trial_ptd_ms;
experiment.unique_sets = unique_sets;
experiment.set_index = set_index;
experiment.amplitudes = amplitudes;
experiment.amplitude_index = amplitude_index;
experiment.ptds_ms = ptds_ms;
experiment.ptd_index = ptd_index;
end

function names = normalize_names(values)
names = cell(size(values));
for k = 1:numel(values)
    value = values{k};
    if isstring(value) || ischar(value)
        names{k} = strtrim(char(string(value)));
    elseif isempty(value)
        names{k} = '';
    else
        names{k} = strtrim(char(string(value)));
    end
end
end

function values = numeric_cells(cells_in, field_name, trial_id)
values = nan(1,numel(cells_in));
for k = 1:numel(cells_in)
    value = cells_in{k};
    if isnumeric(value) && isscalar(value)
        values(k) = double(value);
    else
        error('StageA:NonNumericStimParam', ...
            'Non-numeric %s value in trial %d, event %d.', ...
            field_name, trial_id, k);
    end
end
end

