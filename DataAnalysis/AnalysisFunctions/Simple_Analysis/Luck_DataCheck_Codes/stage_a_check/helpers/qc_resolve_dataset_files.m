function files = qc_resolve_dataset_files(folder)
%QC_RESOLVE_DATASET_FILES Resolve known experiment filenames deterministically.

files = struct();
files.spike = unique_required(folder, '*sp_xia_SSD.mat', 'spike');
files.experiment = unique_required(folder, '*_exp_datafile_*.mat', 'experiment');
files.trigger = unique_required(folder, '*.trig.dat', 'trigger');

files.responding = first_pattern(folder, {
    '*_RespondingChannels.mat'
    }, {'BACKUP'});

files.bad_trials = first_pattern(folder, {
    '*.SimSeqBadTrials.mat'
    '*.MultiSIsBadTrials.mat'
    '*.MultiISIsBadTrials.mat'
    '*._MultiISIsBadTrials.mat'
    '*.BadTrials.mat'
    }, {'BACKUP'});

files.bad_channels = first_pattern(folder, {
    '*.MultiSIsBadChannels.mat'
    '*.MultiISIsBadChannels.mat'
    '*.BadChannels.mat'
    }, {'BACKUP'});
end

function out = unique_required(folder, pattern, label)
listing = filter_listing(dir(fullfile(folder, pattern)), {});
out = empty_file_record();
out.pattern = pattern;
out.candidates = {listing.name};

if isempty(listing)
    return;
end
if numel(listing) > 1
    names = strjoin({listing.name}, ', ');
    error('StageA:AmbiguousFile', ...
        'Multiple %s files match %s in %s: %s', label, pattern, folder, names);
end
out = fill_record(out, listing(1), folder);
end

function out = first_pattern(folder, patterns, exclusions)
out = empty_file_record();
out.patterns_checked = patterns;
out.candidates = {};

for k = 1:numel(patterns)
    listing = filter_listing(dir(fullfile(folder, patterns{k})), exclusions);
    out.candidates = [out.candidates, {listing.name}]; %#ok<AGROW>
    if isempty(listing)
        continue;
    end
    if numel(listing) > 1
        names = strjoin({listing.name}, ', ');
        error('StageA:AmbiguousFile', ...
            'Multiple files match priority pattern %s in %s: %s', ...
            patterns{k}, folder, names);
    end
    out.pattern = patterns{k};
    out = fill_record(out, listing(1), folder);
    return;
end
end

function listing = filter_listing(listing, exclusions)
if isempty(listing), return; end
keep = ~[listing.isdir];
for k = 1:numel(exclusions)
    keep = keep & ~contains({listing.name}, exclusions{k}, 'IgnoreCase', true);
end
listing = listing(keep);
end

function out = empty_file_record()
out = struct('found',false, 'name','', 'path','', 'bytes',NaN, ...
    'pattern','', 'patterns_checked',{{}}, 'candidates',{{}});
end

function out = fill_record(out, item, folder)
out.found = true;
out.name = item.name;
out.path = fullfile(folder, item.name);
out.bytes = item.bytes;
end

